/*
 * subworkflows/local/vcf_postprocess.nf
 *
 * Postprocessing path applied uniformly to every caller's output:
 * QC filter → normalise → annotate (VEP) → benchmark (hap.py vs GIAB truth).
 *
 * Input is a channel of tagged VCFs (sample, caller, vcf), plus a per-sample
 * truth-asset channel (sample, truth_vcf, truth_tbi, truth_bed). Benchmark
 * only runs for samples that have a truth_vcf and truth_bed declared.
 */

include { vcfQc        } from '../../modules/local/vcf_qc.nf'
include { normaliseVcf } from '../../modules/local/normalise_vcf.nf'
include { annotateVcf  } from '../../modules/local/annotate_vcf.nf'
include { benchmark    } from '../../modules/local/benchmark.nf'

workflow VCF_POSTPROCESS {
    take:
    tagged_vcfs    // (sample, caller, vcf)
    truth_ch       // (sample, truth_vcf, truth_tbi, truth_bed)
    reference
    fai

    main:
    vcfQc(tagged_vcfs)
    normaliseVcf(vcfQc.out.vcf, reference)
    annotateVcf(normaliseVcf.out.vcf, reference, fai)

    // ── Benchmark: join normalised (sample, caller, vcf) with truth_ch (sample, ...)
    //    by sample using combine-by. Only samples with truth assets are kept.
    normaliseVcf.out.vcf
        .combine(truth_ch, by: 0)
        .filter { sample, caller, vcf, truth_vcf, truth_tbi, truth_bed ->
            truth_vcf != [] && truth_bed != [] }
        .set { benchmark_input }

    benchmark(benchmark_input, reference, fai)

    emit:
    annotated_vcf    = annotateVcf.out.vcf
    benchmark_summary = benchmark.out.summary
}
