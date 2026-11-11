# Why the audit does not attempt star-allele calling

Written as an answer to the query on the manuscript. The claim was previously
an argument assembled from two sources; it is now a measurement made on the
harmonised panel itself, so the paper can report it as a result. Supporting
table: `star_allele_probe_coverage.csv`.

---

## Short reply, for the comment thread

Deliberately, and the paper now says so with numbers rather than by assertion.
A star allele is defined by its defining variants, and genotype-to-phenotype
translation requires those positions to be observed rather than inferred
(Gaedigk et al., 2018, 2021; Caudle et al., 2020). The harmonised panel does
not observe any of them: across 32 defining variants for the commonly reported
alleles of the seven CPIC-actionable CYP genes, none coincides with the
interrogated CpG of any of the 454,181 probes shared by the two designs. The
median distance from a defining variant to the nearest shared probe is 3.1 kb,
the largest 18.3 kb (CYP2C19\*2, rs4244285), and 18 of the 32 have no shared
probe within 1 kb. The 59 identity probes, which are the only genotype-bearing
content common to both arrays, all lie elsewhere: none falls inside any of the
seven gene windows, and the nearest sits 1.37 Mb from any defining variant.
Three defining variants do fall within 50 bp of a probe, which is a reason to
distrust those probes rather than to read them as genotypes. A new paragraph in
the Discussion states this with the table in the supplement.

---

## For the Methods or Discussion

> Star-allele calling was not attempted, for three reasons that the audit
> itself makes quantitative. First, a star allele is defined by the presence of
> its defining variants, and the translation of genotype to phenotype requires
> those positions to be observed rather than imputed (Gaedigk et al., 2018,
> 2021; Caudle et al., 2020). None of them is interrogated here: of 32 defining
> variants for the commonly reported alleles of CYP1A2, CYP2B6, CYP2C9,
> CYP2C19, CYP2D6, CYP3A4 and CYP3A5 (Supplementary Table S9), not one
> coincides with the interrogated CpG dinucleotide of any of the 454,181 probes
> shared by the 450K and EPIC designs; the median distance to the nearest
> shared probe is 3,118 bp, the maximum 18,269 bp at CYP2C19\*2, and 18 of the
> 32 variants have no shared probe within 1 kb. The 59 identity probes carried
> by both designs, the only genotype-bearing content on the harmonised panel,
> lie outside all seven gene windows, the nearest 1.37 Mb from any defining
> variant.
>
> Second, where a defining variant does fall within a probe's hybridisation
> footprint, the probe reports an artefact rather than a genotype. Three
> variants lie within 50 bp of a shared probe's CpG: rs28371706 (CYP2D6\*17) 24
> bp from cg07016288, rs28399504 (CYP2C19\*4) 30 bp from cg18564458, and
> rs3892097 (CYP2D6\*4) 37 bp from cg25452165. A polymorphism under the probe
> body perturbs hybridisation in an allele-dependent way and so shifts the beta
> value (Zhou et al., 2017), which is the mechanism that produces a gapped
> distribution at a CpG probe and which is documented here at cg08561071 and
> rs113512828. Reading such a probe as an allele call would invert the
> inference, reporting a measurement artefact as a genotype, with no second
> measurement available to check it and no way to establish which group
> corresponds to which allele, since the direction depends on probe design and
> on the methylation state at the site.
>
> Third, the pharmacologically decisive variation at CYP2D6 is structural:
> whole-gene deletion, duplication and multiplication, and CYP2D6-CYP2D7
> hybrid alleles (Nofziger et al., 2020; Kim et al., 2012). A methylation array
> reports an intensity ratio at a CpG, not the copy number of a locus, so even
> complete SNP coverage of the gene body would not distinguish one functional
> copy from three.
>
> The audit therefore answers the question it is posed, which is what
> pharmacogenetic information harmonised 450K and EPIC data contain, and the
> answer is that they contain sample identity but not clinically actionable
> metaboliser genotype. Obtaining the latter requires a genotyping assay: a
> targeted pharmacogenomic panel, exome or genome sequencing, or long-read
> sequencing where CYP2D6 structure matters.

## Method note, for the supplement

Coordinates are GRCh37, matching the harmonised manifest. Defining-variant
positions were resolved from dbSNP rather than typed, and the allele sets are
the commonly reported alleles for each gene, which is not exhaustive: rarer
alleles add defining positions without adding probes, so an exhaustive list can
only widen the gaps reported here. Gene windows are RefSeq transcript spans
from the UCSC browser plus 1,500 bp, taken deliberately generously (the CYP3A4
and CYP3A5 windows overlap), which only increases the number of probes tested
against each variant. Probe footprints are tested as the interrogated CpG
dinucleotide and, separately, as a symmetric 50-base window around it, because
the manifest records no strand for the probe body; the symmetric window is the
more permissive of the two possible orientations.

## Numbers quoted above, with their source

| Quantity | Value | Where it comes from |
|---|---|---|
| Probes shared by the two designs | 454,181 | `crossarray_probe_map.csv.gz`, `status == shared_450K_EPIC` |
| Shared probes with manifest coordinates | 454,122 | the same file; the 59 identity probes carry no MAPINFO |
| Defining variants tested | 32 | `star_allele_probe_coverage.csv` |
| Variants at an interrogated CpG | 0 | the same table, `at_cpg_site` |
| Variants within 50 bp of a probe | 3 | the same table, `within_50bp_footprint` |
| Median distance to the nearest shared probe | 3,118 bp | the same table, `distance_bp` |
| Largest distance | 18,269 bp (CYP2C19\*2) | the same table |
| Variants with no shared probe within 1 kb | 18 of 32 | the same table, `probes_within_1kb` |
| Identity probes inside the seven gene windows | 0 of 59 | `identity_probe_positions.csv` |
| Nearest identity probe to any defining variant | 1.37 Mb | the same file |
| Harmonised probes in the seven gene windows | 130 | the generous windows described above |

The 130 figure is coordinate-based and so is not produced the same way as the
47 probes reported in the gene-annotation audit, which uses the manifest's
`UCSC_RefGene_Name`; the two should not be quoted side by side without saying
which definition is in use.
