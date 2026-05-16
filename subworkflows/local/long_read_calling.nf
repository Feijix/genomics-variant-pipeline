/*
 * subworkflows/local/long_read_calling.nf
 *
 * Long-read variant calling path. Accepts FASTQ input (aligned via minimap2)
 * and pre-aligned BAM input (indexed). Converges via mix() then runs:
 *   - Clair3 on samples that declare clair3_platform + clair3_model
 *   - DeepVariant on samples that declare dv_model
 *   - Ensemble intersect (BCFTOOLS_ISEC) and union (BCFTOOLS_CONCAT) on samples
 *     that produced both Clair3 and DeepVariant outputs
 *
 * Caller tagging happens at the subworkflow boundary so downstream postprocess
 * sees a uniform (meta+caller, vcf) channel shape.
 */

include { MINIMAP2_ALIGN  } from '../../modules/local/minimap2_align.nf'
include { SAMTOOLS_SORT   } from '../../modules/local/samtools_sort.nf'
include { SAMTOOLS_INDEX  } from '../../modules/local/samtools_index.nf'
include { CLAIR3          } from '../../modules/local/clair3.nf'
include { DEEPVARIANT     } from '../../modules/local/deepvariant.nf'
include { BCFTOOLS_ISEC   } from '../../modules/local/bcftools_isec.nf'
include { BCFTOOLS_CONCAT } from '../../modules/local/bcftools_concat.nf'

workflow LONG_READ_CALLING {
    take:
    longread_fq_samples     // (meta, reads, aligned_bam, truth_vcf, truth_bed)
    longread_bam_samples    // same shape
    reference
    fai

    main:
    // ── Bring both input flavours up to (meta, bam, bai) ──
    MINIMAP2_ALIGN(
        longread_fq_samples.map { meta, reads, _bam, _vcf, _bed -> tuple(meta, reads[0]) },
        reference
    )
    SAMTOOLS_SORT(MINIMAP2_ALIGN.out.out)

    SAMTOOLS_INDEX(
        longread_bam_samples.map { meta, _reads, bam, _vcf, _bed -> tuple(meta, bam) }
    )

    SAMTOOLS_SORT.out.out
        .mix(SAMTOOLS_INDEX.out.out)
        .set { longread_bams }

    // ── Clair3 on samples that declare platform/model ──
    longread_bams
        .filter { meta, _bam, _bai -> meta.clair3_platform != '' && meta.clair3_model != '' }
        .set { clair3_input }

    CLAIR3(clair3_input, reference, fai)

    // ── DeepVariant on samples that declare a dv_model ──
    longread_bams
        .filter { meta, _bam, _bai -> meta.dv_model != '' }
        .set { dv_input }

    DEEPVARIANT(dv_input, reference, fai)

    // ── Ensemble: samples that have both Clair3 and DV outputs ──
    CLAIR3.out.out
        .join(DEEPVARIANT.out.out)
        .set { ensemble_input }

    BCFTOOLS_ISEC(ensemble_input)
    BCFTOOLS_CONCAT(ensemble_input, reference)

    // ── Tag each caller's output for downstream uniformity ──
    CLAIR3.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'clair3'], vcf) }
        .set { tagged_clair3 }

    DEEPVARIANT.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'deepvariant_long'], vcf) }
        .set { tagged_dv_long }

    BCFTOOLS_ISEC.out.vcf
        .map { meta, vcf -> tuple(meta + [caller: 'ensemble'], vcf) }
        .set { tagged_ensemble }

    BCFTOOLS_CONCAT.out.vcf
        .map { meta, vcf -> tuple(meta + [caller: 'ensemble_union'], vcf) }
        .set { tagged_ensemble_union }

    emit:
    clair3_vcf         = tagged_clair3
    dv_long_vcf        = tagged_dv_long
    ensemble_vcf       = tagged_ensemble
    ensemble_union_vcf = tagged_ensemble_union
}
