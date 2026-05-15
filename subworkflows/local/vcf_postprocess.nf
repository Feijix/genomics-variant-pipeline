/*
 * subworkflows/local/vcf_postprocess.nf
 *
 * Postprocessing applied uniformly to every caller's output:
 * QC filter → normalise → annotate (VEP) → benchmark (hap.py vs GIAB truth).
 *
 * Input is a channel of tagged VCFs where meta carries both the sample id
 * (meta.id) and the caller name (meta.caller). Benchmark joins each VCF with
 * its per-sample truth assets by meta.id and only runs when both are present.
 */

include { vcfQc        } from '../../modules/local/vcf_qc.nf'
include { normaliseVcf } from '../../modules/local/normalise_vcf.nf'
include { annotateVcf  } from '../../modules/local/annotate_vcf.nf'
include { benchmark    } from '../../modules/local/benchmark.nf'

workflow VCF_POSTPROCESS {
    take:
    tagged_vcfs    // (meta+caller, vcf)
    truth_ch       // (sample_id, truth_vcf, truth_tbi, truth_bed)
    reference
    fai

    main:
    vcfQc(tagged_vcfs)
    normaliseVcf(vcfQc.out.vcf, reference)
    annotateVcf(normaliseVcf.out.vcf, reference, fai)

    // ── Join normalised VCFs with their truth assets by sample id ──
    // combine(by:0) requires a join key in position 0, so first re-key on meta.id
    normaliseVcf.out.vcf
        .map { meta, vcf -> tuple(meta.id, meta, vcf) }
        .combine(truth_ch, by: 0)
        .map { sample_id, meta, vcf, truth_vcf, truth_tbi, truth_bed ->
            tuple(meta, vcf, truth_vcf, truth_tbi, truth_bed) }
        .filter { meta, vcf, truth_vcf, truth_tbi, truth_bed ->
            truth_vcf != [] && truth_bed != [] }
        .set { benchmark_input }

    benchmark(benchmark_input, reference, fai)

    emit:
    annotated_vcf     = annotateVcf.out.vcf
    benchmark_summary = benchmark.out.summary
}
