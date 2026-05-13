/*
 * subworkflows/local/short_read_calling.nf
 *
 * Short-read variant calling path: align → sort → DeepVariant.
 * Input is a channel of branched.shortread samples (full 9-element tuple).
 * Output is a channel of (sample_id, vcf) — caller tagging happens at the
 * workflow level so the downstream postprocess subworkflow can treat all
 * callers uniformly.
 */

include { align        } from '../../modules/local/align.nf'
include { sortBam      } from '../../modules/local/sort_bam.nf'
include { callVariants } from '../../modules/local/call_variants.nf'

workflow SHORT_READ_CALLING {
    take:
    shortread_samples       // (sample, type, reads, aligned_bam, truth_vcf, truth_bed, clair3_platform, clair3_model, dv_model)
    reference
    reference_index
    fai

    main:
    align(
        shortread_samples.map { it -> tuple(it[0], it[2]) },
        reference,
        reference_index
    )
    sortBam(align.out)

    // Pair sorted BAM with the dv_model for the calling step
    sortBam.out
        .join(shortread_samples.map { it -> tuple(it[0], it[8]) })
        .set { dv_input }

    callVariants(dv_input, reference, fai)

    emit:
    vcf = callVariants.out
}
