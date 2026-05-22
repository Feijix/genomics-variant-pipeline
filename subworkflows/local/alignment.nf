/*
 * subworkflows/local/alignment.nf
 *
 * Unified alignment for all sample types. Branches by sequencing_type and
 * whether a pre-aligned BAM is provided, converges all paths to a single
 * (meta, bam, bai) channel.
 *
 *   shortread + FASTQ      -> BWA_MEM -> SAMTOOLS_SORT
 *   longread  + FASTQ      -> MINIMAP2_ALIGN -> SAMTOOLS_SORT
 *   longread  + aligned BAM -> SAMTOOLS_INDEX (no realignment)
 *
 * The minimap2 preset is selected from meta.platform (see MINIMAP2_ALIGN module).
 *
 * shortread with a pre-aligned BAM input isn't currently expected but would
 * route through SAMTOOLS_INDEX naturally if it appeared.
 */

include { BWA_MEM         } from '../../modules/local/bwa_mem.nf'
include { MINIMAP2_ALIGN  } from '../../modules/local/minimap2_align.nf'
include { SAMTOOLS_SORT   } from '../../modules/local/samtools_sort.nf'
include { SAMTOOLS_INDEX  } from '../../modules/local/samtools_index.nf'

workflow ALIGNMENT {
    take:
    samples             // (meta, [fastq_1, fastq_2], aligned_bam, truth_vcf, truth_bed)
    reference
    reference_index     // BWA_MEM_INDEX output, only used for shortread
    fai                 // unused here but kept for API symmetry with calling stages

    main:
    // ── Branch by sequencing_type + input flavour ──
    samples
        .branch { meta, _reads, aligned_bam, _vcf, _bed ->
            shortread_fq: meta.sequencing_type == 'shortread' && aligned_bam == []
            longread_fq:  meta.sequencing_type == 'longread'  && aligned_bam == []
            prealigned:   aligned_bam != []
        }
        .set { branched }

    // ── Short-read FASTQ -> BWA_MEM -> SAMTOOLS_SORT ──
    BWA_MEM(
        branched.shortread_fq.map { meta, reads, _bam, _vcf, _bed -> tuple(meta, reads) },
        reference,
        reference_index
    )

    // ── Long-read FASTQ -> MINIMAP2_ALIGN -> SAMTOOLS_SORT ──
    MINIMAP2_ALIGN(
        branched.longread_fq.map { meta, reads, _bam, _vcf, _bed -> tuple(meta, reads[0]) },
        reference
    )

    // ── Both alignment paths produce SAMs that need sorting + indexing ──
    BWA_MEM.out.out
        .mix(MINIMAP2_ALIGN.out.out)
        .set { sams_to_sort }

    SAMTOOLS_SORT(sams_to_sort)

    // ── Pre-aligned BAMs only need indexing ──
    SAMTOOLS_INDEX(
        branched.prealigned.map { meta, _reads, bam, _vcf, _bed -> tuple(meta, bam) }
    )

    // ── Converge all paths to a single (meta, bam, bai) channel ──
    SAMTOOLS_SORT.out.out
        .mix(SAMTOOLS_INDEX.out.out)
        .set { aligned_bams }

    emit:
    bams = aligned_bams
}
