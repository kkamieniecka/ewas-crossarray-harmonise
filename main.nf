#!/usr/bin/env nextflow
/*
 * ewas-harmonise -- reproducible cross-array harmonisation and region
 * detection for longitudinal EWAS (Illumina 450K + EPIC).
 *
 * Extends the EWASGalaxy tool suite (Murat et al. 2019) with:
 *   - a harmonisation stage that keeps each array in its native probe space
 *     until QC is done, then merges only through minfi::combineArrays()
 *   - replacements for sections 2.6 and 2.7 of Aryee et al. 2014, whose
 *     fixed-gap clustering, fixed-span loess and free permutation are not
 *     valid on a harmonised, repeated-measures design
 *   - a legacy baseline run on the identical matrix, so the change is measured
 *     rather than asserted
 */

nextflow.enable.dsl = 2

params.sheet        = null           // minfi-format sample sheet (required)
params.idat_dir     = null           // directory of *_Grn/_Red.idat.gz (required)
params.probe_map    = null           // crossarray_probe_map.csv.gz (required)
// params.outdir default lives in nextflow.config, which references it for the
// timeline/report/trace/dag paths before this script is evaluated

params.exposure     = 'days_on_clozapine'
params.exposure_scale = 100
params.time_scale   = 'per100days'   // 02_probe_model.R accepts a named scale
params.subject      = 'Subject_ID'
params.covars       = 'smoking_score,cd8t,cd4t,bcell,gran,mono,nk'

params.detp         = 0.01
params.detp_frac    = 0.05
params.drop_sex     = true
params.bmiq         = false          // off by default: shared probes are
                                     // chemically identical across arrays
params.max_gap      = 1000
params.rho_min      = 0.30
params.decay_bp     = 1000
params.lam_grid     = '0.05,0.1,0.25,0.5,1.0,2.0,4.0'
params.n_perm       = 200
params.n_boot       = 50
params.min_probes   = 3
params.min_effect   = 0.05
params.cv_folds     = 5
// Stage 04 exists in R and in Python. They are proved equivalent function by
// function (tests/test_equivalence.R) and produce identical regions on the
// same input; the R implementation is the default because the rest of the
// suite is R, and costs about twice the runtime and twice the memory. Set
// --dmr_impl python to run the Python one, e.g. to audit the agreement.
params.dmr_impl     = 'r'            // 'r' | 'python'
params.var_method   = 'limma'        // R only: 'limma' (squeezeVar) | 'mom'

params.block_max_gap     = 1500
params.block_rho_min     = 0.20
params.block_length_scale = 250000
params.block_min_post    = 0.80
// Stage 05 likewise exists in both languages. Baum-Welch uses no RNG, so the
// two produce byte-identical blocks rather than merely agreeing in
// distribution; R is the default for the same reason as stage 04. Set
// --blocks_impl python to audit the agreement on real data.
params.blocks_impl       = 'r'       // 'r' | 'python'

params.run_baseline = true
params.naive_perm   = true           // also run the legacy permutation, to
                                     // quantify its mis-calibration
params.seed         = 1

def required(name, value) {
    if (value == null) exit 1, "ewas-harmonise: --${name} is required"
    return value
}

process HARMONISE {
    tag 'harmonise'
    publishDir params.outdir, mode: 'copy'
    label 'r_heavy'

    input:
    path sheet
    path idat_dir

    output:
    path 'harmonised/harmonised.rds',         emit: rds
    path 'harmonised/mval.f64',               emit: mval
    path 'harmonised/mval_dims.json',         emit: dims
    path 'harmonised/pheno_used.csv',         emit: pheno
    path 'harmonised/probe_annotation.csv',   emit: anno
    path 'harmonised/qc_samples.csv',         emit: qc
    path 'harmonised/probe_filter.csv',       emit: filt

    script:
    """
    Rscript ${projectDir}/bin/01_harmonise.R \\
        --sheet ${sheet} --idat_dir ${idat_dir} --out_dir harmonised \\
        --detp ${params.detp} --detp_frac ${params.detp_frac} \\
        --drop_sex ${params.drop_sex} --bmiq ${params.bmiq} \\
        --keep_gset TRUE --seed ${params.seed}
    """
}

process PROBE_MODEL {
    tag 'probe-model'
    publishDir params.outdir, mode: 'copy'
    label 'r_heavy'

    input:
    path rds

    output:
    path 'probe_model/*', emit: all

    script:
    """
    Rscript ${projectDir}/bin/02_probe_model.R \\
        --harmonised ${rds} --out_dir probe_model \\
        --exposure ${params.exposure} --covars ${params.covars} \\
        --time_scale ${params.time_scale} --seed ${params.seed}
    """
}

process BASELINE {
    tag 'legacy-bumphunter'
    publishDir params.outdir, mode: 'copy'
    label 'r_heavy'

    input:
    path rds

    output:
    path 'baseline/*', emit: all

    script:
    """
    Rscript ${projectDir}/bin/03_baseline_bumphunter.R \\
        --in-dir . --out-dir baseline \\
        --exposure ${params.exposure} --exposure-scale ${params.exposure_scale} \\
        --n-perm ${params.n_perm} --seed ${params.seed}
    """
}

// Stage 04 exists twice, in R and in Python, with identical command-line
// flags. The argument list is built once here so the two processes cannot
// drift apart.
def dmrArgs() {
    def args = [
        '--in-dir .', '--out-dir dmr',
        "--exposure ${params.exposure}", "--exposure-scale ${params.exposure_scale}",
        "--subject ${params.subject}", "--covars ${params.covars}",
        "--max-gap ${params.max_gap}", "--rho-min ${params.rho_min}",
        "--decay-bp ${params.decay_bp}", "--lam-grid ${params.lam_grid}",
        "--cv-folds ${params.cv_folds}",
        "--min-probes ${params.min_probes}", "--min-effect ${params.min_effect}",
        "--n-perm ${params.n_perm}", "--n-boot ${params.n_boot}",
        "--seed ${params.seed}",
    ]
    if( params.naive_perm ) args << '--also-naive-perm'
    return args.join(' ')
}

process DMR_ML {
    tag 'dmr-ml'
    publishDir params.outdir, mode: 'copy'
    label 'r_heavy'

    input:
    path mval
    path dims
    path pheno
    path anno

    output:
    path 'dmr/dmr_ml.csv',            emit: regions
    path 'dmr/lambda_cv.csv',         emit: cv
    path 'dmr/permutation_null.csv',  emit: null_dist
    path 'dmr/run_config.json',       emit: cfg

    script:
    """
    Rscript ${projectDir}/bin/04_dmr_ml.R ${dmrArgs()} \\
        --var-method ${params.var_method}
    """
}

// The Python implementation of the same stage. Kept in the tree so that the
// agreement between the two can be re-checked on real data at any time
// (--dmr_impl python); tests/test_equivalence.R proves it function by function
// on generated inputs. Its outputs carry the same names, so COMPARE and the
// report do not care which one ran.
process DMR_ML_PY {
    tag 'dmr-ml-py'
    publishDir params.outdir, mode: 'copy'
    label 'py_heavy'

    input:
    path mval
    path dims
    path pheno
    path anno

    output:
    path 'dmr/dmr_ml.csv',            emit: regions
    path 'dmr/lambda_cv.csv',         emit: cv
    path 'dmr/permutation_null.csv',  emit: null_dist
    path 'dmr/run_config.json',       emit: cfg

    script:
    """
    python ${projectDir}/bin/04_dmr_ml.py ${dmrArgs()}
    """
}

// Stage 05 arguments, shared by the two implementations so that an audit run
// differs only in which interpreter is invoked.
def blocksArgs() {
    return [
        "--in-dir .", "--out-dir blocks",
        "--exposure ${params.exposure}",
        "--exposure-scale ${params.exposure_scale}",
        "--subject ${params.subject}", "--covars ${params.covars}",
        "--max-gap ${params.block_max_gap}",
        "--rho-min ${params.block_rho_min}",
        "--length-scale ${params.block_length_scale}",
        "--min-post ${params.block_min_post}",
        "--fixed-collapse",
        "--seed ${params.seed}",
    ].join(' ')
}

process BLOCKS_HSMM {
    tag 'blocks-hsmm'
    publishDir params.outdir, mode: 'copy'
    label 'r_heavy'

    input:
    path mval
    path dims
    path pheno
    path anno
    path probe_map

    output:
    path 'blocks/blocks_hsmm.csv',                emit: blocks
    path 'blocks/openSea_cluster_effects.csv.gz', emit: clusters
    path 'blocks/hsmm_params.json',               emit: cfg

    script:
    """
    Rscript ${projectDir}/bin/05_blocks_hsmm.R ${blocksArgs()} \\
        --probe-map ${probe_map} --var-method ${params.var_method}
    """
}

// The Python implementation of stage 05, kept so the agreement can be
// re-checked on real data (--blocks_impl python). Output names match, so
// COMPARE and the report do not care which one ran.
process BLOCKS_HSMM_PY {
    tag 'blocks-hsmm-py'
    publishDir params.outdir, mode: 'copy'
    label 'py_heavy'

    input:
    path mval
    path dims
    path pheno
    path anno
    path probe_map

    output:
    path 'blocks/blocks_hsmm.csv',                emit: blocks
    path 'blocks/openSea_cluster_effects.csv.gz', emit: clusters
    path 'blocks/hsmm_params.json',               emit: cfg

    script:
    """
    python ${projectDir}/bin/05_blocks_hsmm.py ${blocksArgs()} \\
        --probe-map ${probe_map}
    """
}

process COMPARE {
    tag 'compare'
    publishDir params.outdir, mode: 'copy'
    label 'py_light'

    input:
    path 'dmr/*'
    path 'blocks/*'
    path 'probe_model/*'
    path 'baseline/*'

    output:
    path 'comparison/*', emit: all

    script:
    // The baseline directory is empty when --run_baseline false; 06_compare.py
    // then reports the replacement methods only. --probe-model-dir is what
    // produces the *_common columns (every method's per-array region effect
    // re-estimated from the same within-subject per-probe fits), which are the
    // only cross-array columns comparable between methods.
    """
    python ${projectDir}/bin/06_compare.py \\
        --dmr-dir dmr --blocks-dir blocks --baseline-dir baseline \\
        --probe-model-dir probe_model \\
        --out-dir comparison
    """
}

workflow {
    sheet = file(required('sheet', params.sheet))
    idats = file(required('idat_dir', params.idat_dir))
    pmap  = file(required('probe_map', params.probe_map))

    def dmr_impl = "${params.dmr_impl}".toLowerCase()
    if( !(dmr_impl in ['r', 'python']) )
        exit 1, "ewas-harmonise: --dmr_impl must be 'r' or 'python', got '${params.dmr_impl}'"

    def blocks_impl = "${params.blocks_impl}".toLowerCase()
    if( !(blocks_impl in ['r', 'python']) )
        exit 1, "ewas-harmonise: --blocks_impl must be 'r' or 'python', got '${params.blocks_impl}'"

    HARMONISE(sheet, idats)
    PROBE_MODEL(HARMONISE.out.rds)
    if( dmr_impl == 'r' )
        DMR_ML(HARMONISE.out.mval, HARMONISE.out.dims,
               HARMONISE.out.pheno, HARMONISE.out.anno)
    else
        DMR_ML_PY(HARMONISE.out.mval, HARMONISE.out.dims,
                  HARMONISE.out.pheno, HARMONISE.out.anno)
    dmr_out = dmr_impl == 'r' ? DMR_ML.out : DMR_ML_PY.out
    if( blocks_impl == 'r' )
        BLOCKS_HSMM(HARMONISE.out.mval, HARMONISE.out.dims,
                    HARMONISE.out.pheno, HARMONISE.out.anno, pmap)
    else
        BLOCKS_HSMM_PY(HARMONISE.out.mval, HARMONISE.out.dims,
                       HARMONISE.out.pheno, HARMONISE.out.anno, pmap)
    blocks_out = blocks_impl == 'r' ? BLOCKS_HSMM.out : BLOCKS_HSMM_PY.out

    // COMPARE runs whether or not the legacy baseline did: with it, the table
    // is the four-method comparison; without it (--run_baseline false, as in
    // -profile test) it still reports the replacement methods, so the smoke
    // test exercises stage 06 instead of stopping one stage short.
    if (params.run_baseline) {
        BASELINE(HARMONISE.out.rds)
        baseline_files = BASELINE.out.all.collect()
    }
    else {
        baseline_files = Channel.value([])
    }

    COMPARE(dmr_out.regions.mix(dmr_out.cfg).collect(),
            blocks_out.blocks.mix(blocks_out.cfg).collect(),
            PROBE_MODEL.out.all.collect(),
            baseline_files)
}

// The run summary handler lives in nextflow.config: a top-level
// `workflow.onComplete { }` statement here is rejected by the strict script
// parser ("statements cannot be mixed with script declarations").
