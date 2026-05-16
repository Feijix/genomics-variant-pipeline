/*
 * subworkflows/local/short_read_calling.nf
 *
 * Short-read variant calling path: align → sort → DeepVariant.
 * Input shape: (meta, [fastq_1, fastq_2], aligned_bam, truth_vcf, truth_bed).
 * Output: (meta, vcf) with meta.caller = 'deepvariant'.
 */

include { BWA_MEM       } from '../../modules/local/bwa_mem.nf'
include { SAMTOOLS_SORT } from '../../modules/local/samtools_sort.nf'
include { DEEPVARIANT   } from '../../modules/local/deepvariant.nf'

workflow SHORT_READ_CALLING {
    take:
    shortread_samples       // (meta, reads, aligned_bam, truth_vcf, truth_bed)
    reference
    reference_index
    fai

    main:
    BWA_MEM(
        shortread_samples.map { meta, reads, _bam, _vcf, _bed -> tuple(meta, reads) },
        reference,
        reference_index
    )
    SAMTOOLS_SORT(BWA_MEM.out.out)

    DEEPVARIANT(SAMTOOLS_SORT.out.out, reference, fai)

    DEEPVARIANT.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'deepvariant'], vcf) }
        .set { tagged_vcf }

    emit:
    vcf = tagged_vcf
}
