/*
 * subworkflows/local/quality_control.nf
 *
 * Collects QC from three points in the pipeline and aggregates into one report:
 *   aligned BAM  -> SAMTOOLS_STATS, MOSDEPTH
 *   called VCF   -> BCFTOOLS_STATS
 *   markdup      -> metrics file already produced by GATK_MARKDUPLICATES
 */

include { SAMTOOLS_STATS } from '../../modules/local/samtools_stats.nf'
include { MOSDEPTH       } from '../../modules/local/mosdepth.nf'
include { BCFTOOLS_STATS } from '../../modules/local/bcftools_stats.nf'
include { MULTIQC        } from '../../modules/local/multiqc.nf'

workflow QUALITY_CONTROL {
    take:
    bams            // (meta, bam, bai)
    vcfs            // (meta+caller, vcf)
    extra_metrics   // path channel: MarkDuplicates metrics, VEP summaries, etc.

    main:
    SAMTOOLS_STATS(bams)
    MOSDEPTH(bams)
    BCFTOOLS_STATS(vcfs)

    // Strip meta — MultiQC identifies samples from filenames.
    SAMTOOLS_STATS.out.stats.map { _meta, f -> f }
        .mix(MOSDEPTH.out.summary.map        { _meta, f -> f })
        .mix(BCFTOOLS_STATS.out.stats.map    { _meta, f -> f })
        .mix(extra_metrics)
        .collect()
        .set { all_qc }

    MULTIQC(all_qc)

    emit:
    report = MULTIQC.out.report
}
