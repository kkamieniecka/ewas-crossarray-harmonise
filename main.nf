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
params.outdir       = 'results'

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

params.block_max_gap     = 1500
params.block_rho_min     = 0.20
params.block_length_scale = 250000
params.block_min_post    = 0.80

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
    publishDir "${params.outdir}/harmonised", mode: 'copy'
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
    publishDir "${params.outdir}/probe_model", mode: 'copy'
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
    publishDir "${params.outdir}/baseline", mode: 'copy'
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

process DMR_ML {
    tag 'dmr-ml'
    publishDir "${params.outdir}/dmr", mode: 'copy'
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
    def naive = params.naive_perm ? '--also-naive-perm' : ''
    """
    python ${projectDir}/bin/04_dmr_ml.py \\
        --in-dir . --out-dir dmr \\
        --exposure ${params.exposure} --exposure-scale ${params.exposure_scale} \\
        --subject ${params.subject} --covars ${params.covars} \\
        --max-gap ${params.max_gap} --rho-min ${params.rho_min} \\
        --decay-bp ${params.decay_bp} --lam-grid ${params.lam_grid} \\
        --min-probes ${params.min_probes} --min-effect ${params.min_effect} \\
        --n-perm ${params.n_perm} --n-boot ${params.n_boot} \\
        --seed ${params.seed} ${naive}
    """
}

process BLOCKS_HSMM {
    tag 'blocks-hsmm'
    publishDir "${params.outdir}/blocks", mode: 'copy'
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
    python ${projectDir}/bin/05_blocks_hsmm.py \\
        --in-dir . --probe-map ${probe_map} --out-dir blocks \\
        --exposure ${params.exposure} --exposure-scale ${params.exposure_scale} \\
        --subject ${params.subject} --covars ${params.covars} \\
        --max-gap ${params.block_max_gap} --rho-min ${params.block_rho_min} \\
        --length-scale ${params.block_length_scale} \\
        --min-post ${params.block_min_post} --fixed-collapse \\
        --seed ${params.seed}
    """
}

process COMPARE {
    tag 'compare'
    publishDir "${params.outdir}/comparison", mode: 'copy'
    label 'py_light'

    input:
    path 'dmr/*'
    path 'blocks/*'
    path 'baseline/*'

    output:
    path 'comparison/*', emit: all

    script:
    """
    python ${projectDir}/bin/06_compare.py \\
        --dmr-dir dmr --blocks-dir blocks --baseline-dir baseline \\
        --out-dir comparison
    """
}

workflow {
    sheet = file(required('sheet', params.sheet))
    idats = file(required('idat_dir', params.idat_dir))
    pmap  = file(required('probe_map', params.probe_map))

    HARMONISE(sheet, idats)
    PROBE_MODEL(HARMONISE.out.rds)
    DMR_ML(HARMONISE.out.mval, HARMONISE.out.dims,
           HARMONISE.out.pheno, HARMONISE.out.anno)
    BLOCKS_HSMM(HARMONISE.out.mval, HARMONISE.out.dims,
                HARMONISE.out.pheno, HARMONISE.out.anno, pmap)

    if (params.run_baseline) {
        BASELINE(HARMONISE.out.rds)
        COMPARE(DMR_ML.out.regions.mix(DMR_ML.out.cfg).collect(),
                BLOCKS_HSMM.out.blocks.mix(BLOCKS_HSMM.out.cfg).collect(),
                BASELINE.out.all.collect())
    }
}

workflow.onComplete {
    log.info """
    ewas-harmonise finished
      status    : ${workflow.success ? 'OK' : 'FAILED'}
      duration  : ${workflow.duration}
      outdir    : ${params.outdir}
      revision  : ${workflow.revision ?: 'n/a'}  commit: ${workflow.commitId ?: 'n/a'}
    """.stripIndent()
}
