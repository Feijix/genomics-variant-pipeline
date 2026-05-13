/*
 * subworkflows/local/long_read_calling.nf
 *
 * Long-read variant calling path. Accepts both FASTQ-input and pre-aligned-BAM
 * input samples; converges them via mix() after alignment/indexing; then runs:
 *   - Clair3 on all long-read samples
 *   - DeepVariant on samples that declare a dv_model
 *   - Ensemble merge (intersection) and union of Clair3 + DV outputs
 *
 * Emits one channel per caller-flavour. The workflow-level orchestration is
 * responsible for tagging and mixing for downstream postprocess.
 */

include { alignLongRead       } from '../../modules/local/align_long_read.nf'
include { sortBamLong         } from '../../modules/local/sort_bam_long.nf'
include { indexAlignedBam     } from '../../modules/local/index_aligned_bam.nf'
include { callVariantsLong    } from '../../modules/local/call_variants_long.nf'
include { callVariantsLongDV  } from '../../modules/local/call_variants_long_dv.nf'
include { ensembleMergeLong   } from '../../modules/local/ensemble_merge_long.nf'
include { ensembleUnionLong   } from '../../modules/local/ensemble_union_long.nf'

workflow LONG_READ_CALLING {
    take:
    longread_fq_samples     // (sample, type, reads, aligned_bam, truth_vcf, truth_bed, clair3_platform, clair3_model, dv_model)
    longread_bam_samples    // same shape
    reference
    fai

    main:
    // ── Bring both input flavours up to (sample, bam, bai) ──
    alignLongRead(
        longread_fq_samples.map { it -> tuple(it[0], it[2][0]) },
        reference
    )
    sortBamLong(alignLongRead.out)

    indexAlignedBam(
        longread_bam_samples.map { it -> tuple(it[0], it[3]) }
    )

    sortBamLong.out
        .mix(indexAlignedBam.out)
        .set { longread_bams }

    // ── Combined source of per-sample metadata for both input flavours ──
    longread_fq_samples
        .mix(longread_bam_samples)
        .set { all_longread_samples }

    // ── Clair3 on every long-read sample ──
    longread_bams
        .join(
            all_longread_samples.map { it -> tuple(it[0], it[6], it[7]) }
        )
        .set { clair3_input }

    callVariantsLong(clair3_input, reference, fai)

    // ── DeepVariant only for samples that declare a dv_model ──
    longread_bams
        .join(
            all_longread_samples.map { it -> tuple(it[0], it[8]) }
        )
        .filter { sample, bam, bai, dv_model -> dv_model != '' }
        .set { dv_long_input }

    callVariantsLongDV(dv_long_input, reference, fai)

    // ── Ensemble merges: only samples that have both Clair3 and DV outputs ──
    callVariantsLong.out
        .join(callVariantsLongDV.out)
        .set { ensemble_input }

    ensembleMergeLong(ensemble_input)
    ensembleUnionLong(ensemble_input, reference)

    emit:
    clair3_vcf         = callVariantsLong.out
    dv_long_vcf        = callVariantsLongDV.out
    ensemble_vcf       = ensembleMergeLong.out.vcf
    ensemble_union_vcf = ensembleUnionLong.out.vcf
}
