/*
 * subworkflows/local/vcf_postprocess.nf
 *
 * Uniform per-caller postprocess:
 *   BCFTOOLS_VIEW (FILTER=PASS) → BCFTOOLS_NORM (left-align, split multiallelics)
 *   → ENSEMBLVEP_VEP (annotate)  → HAPPY_HAPPY (benchmark vs GIAB truth)
 *
 * Input: tagged VCF channel where meta carries id and caller. Truth assets
 * joined by meta.id; benchmark only runs when both truth_vcf and truth_bed are
 * present for that sample.
 */

include { BCFTOOLS_VIEW   } from '../../modules/local/bcftools_view.nf'
include { BCFTOOLS_NORM   } from '../../modules/local/bcftools_norm.nf'
include { ENSEMBLVEP_VEP  } from '../../modules/local/ensemblvep_vep.nf'
include { HAPPY_HAPPY     } from '../../modules/local/happy_happy.nf'

workflow VCF_POSTPROCESS {
    take:
    tagged_vcfs    // (meta+caller, vcf)
    truth_ch       // (sample_id, truth_vcf, truth_tbi, truth_bed)
    reference
    fai

    main:
    BCFTOOLS_VIEW(tagged_vcfs)
    BCFTOOLS_NORM(BCFTOOLS_VIEW.out.vcf, reference)
    ENSEMBLVEP_VEP(BCFTOOLS_NORM.out.vcf, reference, fai)

    // Re-key the normalised VCFs on meta.id so combine(by:0) joins with truth_ch
    BCFTOOLS_NORM.out.vcf
        .map { meta, vcf -> tuple(meta.id, meta, vcf) }
        .combine(truth_ch, by: 0)
        .map { sample_id, meta, vcf, truth_vcf, truth_tbi, truth_bed ->
            tuple(meta, vcf, truth_vcf, truth_tbi, truth_bed) }
        .filter { meta, vcf, truth_vcf, truth_tbi, truth_bed ->
            truth_vcf != [] && truth_bed != [] }
        .set { benchmark_input }

    HAPPY_HAPPY(benchmark_input, reference, fai)

    emit:
    annotated_vcf     = ENSEMBLVEP_VEP.out.vcf
    benchmark_summary = HAPPY_HAPPY.out.summary
}
