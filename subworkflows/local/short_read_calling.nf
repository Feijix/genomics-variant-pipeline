/*
 * subworkflows/local/short_read_calling.nf
 *
 * Short-read variant calling path: align → sort → DeepVariant.
 * Input is the shortread branch of the parsed samplesheet, with shape
 * (meta, [fastq_1, fastq_2], aligned_bam, truth_vcf, truth_bed).
 *
 * Output is (meta, vcf) with meta.caller set to 'deepvariant' for downstream
 * postprocessing.
 */

include { align       } from '../../modules/local/align.nf'
include { sortBam     } from '../../modules/local/sort_bam.nf'
include { deepvariant } from '../../modules/local/deepvariant.nf'

workflow SHORT_READ_CALLING {
    take:
    shortread_samples       // (meta, reads, aligned_bam, truth_vcf, truth_bed)
    reference
    reference_index
    fai

    main:
    align(
        shortread_samples.map { meta, reads, _bam, _vcf, _bed -> tuple(meta, reads) },
        reference,
        reference_index
    )
    sortBam(align.out)

    deepvariant(sortBam.out, reference, fai)

    // Tag with caller for downstream uniformity
    deepvariant.out
        .map { meta, vcf -> tuple(meta + [caller: 'deepvariant'], vcf) }
        .set { tagged_vcf }

    emit:
    vcf = tagged_vcf
}
