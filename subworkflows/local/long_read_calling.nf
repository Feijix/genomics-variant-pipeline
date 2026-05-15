/*
 * subworkflows/local/long_read_calling.nf
 *
 * Long-read variant calling path. Accepts both FASTQ-input and pre-aligned-BAM
 * input samples; converges them via mix() after alignment/indexing; then runs:
 *   - Clair3 on samples that declare clair3_platform
 *   - DeepVariant on samples that declare dv_model
 *   - Ensemble merge (intersect) and union on samples that ran both
 *
 * Caller tagging happens at the subworkflow boundary — each emitted channel
 * has meta.caller set so downstream postprocess is uniform.
 */

include { alignLongRead     } from '../../modules/local/align_long_read.nf'
include { sortBam           } from '../../modules/local/sort_bam.nf'
include { indexAlignedBam   } from '../../modules/local/index_aligned_bam.nf'
include { callVariantsLong  } from '../../modules/local/call_variants_long.nf'
include { deepvariant       } from '../../modules/local/deepvariant.nf'
include { ensembleMergeLong } from '../../modules/local/ensemble_merge_long.nf'
include { ensembleUnionLong } from '../../modules/local/ensemble_union_long.nf'

workflow LONG_READ_CALLING {
    take:
    longread_fq_samples     // (meta, reads, aligned_bam, truth_vcf, truth_bed)
    longread_bam_samples    // same shape
    reference
    fai

    main:
    // ── Bring both input flavours up to (meta, bam, bai) ──
    alignLongRead(
        longread_fq_samples.map { meta, reads, _bam, _vcf, _bed -> tuple(meta, reads[0]) },
        reference
    )
    sortBam(alignLongRead.out)

    indexAlignedBam(
        longread_bam_samples.map { meta, _reads, bam, _vcf, _bed -> tuple(meta, bam) }
    )

    sortBam.out
        .mix(indexAlignedBam.out)
        .set { longread_bams }

    // ── Clair3 on samples that declare a platform/model ──
    longread_bams
        .filter { meta, _bam, _bai -> meta.clair3_platform != '' && meta.clair3_model != '' }
        .set { clair3_input }

    callVariantsLong(clair3_input, reference, fai)

    // ── DeepVariant on samples that declare a dv_model ──
    longread_bams
        .filter { meta, _bam, _bai -> meta.dv_model != '' }
        .set { dv_input }

    deepvariant(dv_input, reference, fai)

    // ── Ensemble: only samples that produced both Clair3 and DV outputs ──
    callVariantsLong.out
        .join(deepvariant.out)
        .set { ensemble_input }

    ensembleMergeLong(ensemble_input)
    ensembleUnionLong(ensemble_input, reference)

    // ── Tag each caller's output for downstream uniformity ──
    callVariantsLong.out
        .map { meta, vcf -> tuple(meta + [caller: 'clair3'], vcf) }
        .set { tagged_clair3 }

    deepvariant.out
        .map { meta, vcf -> tuple(meta + [caller: 'deepvariant_long'], vcf) }
        .set { tagged_dv_long }

    ensembleMergeLong.out.vcf
        .map { meta, vcf -> tuple(meta + [caller: 'ensemble'], vcf) }
        .set { tagged_ensemble }

    ensembleUnionLong.out.vcf
        .map { meta, vcf -> tuple(meta + [caller: 'ensemble_union'], vcf) }
        .set { tagged_ensemble_union }

    emit:
    clair3_vcf         = tagged_clair3
    dv_long_vcf        = tagged_dv_long
    ensemble_vcf       = tagged_ensemble
    ensemble_union_vcf = tagged_ensemble_union
}
