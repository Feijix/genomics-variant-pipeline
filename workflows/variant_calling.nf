/*
 * workflows/variant_calling.nf
 *
 * Main orchestration. Loads the samplesheet, branches by sequencing technology
 * and input flavour, dispatches each branch to the appropriate subworkflow,
 * tags VCFs with caller names, and feeds the unified tagged stream into
 * postprocessing.
 */

include { indexReference     } from '../modules/local/index_reference.nf'
include { faiReference       } from '../modules/local/fai_reference.nf'

include { SHORT_READ_CALLING } from '../subworkflows/local/short_read_calling.nf'
include { LONG_READ_CALLING  } from '../subworkflows/local/long_read_calling.nf'
include { VCF_POSTPROCESS    } from '../subworkflows/local/vcf_postprocess.nf'

workflow VARIANT_CALLING {

    // ── Reference indexing ──
    indexReference(file(params.reference))
    faiReference(file(params.reference))

    // ── Load samplesheet, validate, branch ──
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def fastq_2          = row.fastq_2          ? file(row.fastq_2)     : []
            def aligned_bam      = row.aligned_bam      ? file(row.aligned_bam) : []
            def truth_vcf        = row.truth_vcf        ? file(row.truth_vcf)   : []
            def truth_bed        = row.truth_bed        ? file(row.truth_bed)   : []
            def clair3_platform  = row.clair3_platform  ?: ''
            def clair3_model     = row.clair3_model     ?: ''
            def dv_model         = row.dv_model         ?: ''
            def fastq_1          = row.fastq_1          ? file(row.fastq_1)     : []

            // Per-sample validation
            // if (row.sequencing_type == 'longread' && !clair3_model) {
            //     error "Sample '${row.sample}' is longread but has no clair3_model declared"
            // }
            // if (row.sequencing_type == 'longread' && !clair3_platform) {
            //     error "Sample '${row.sample}' is longread but has no clair3_platform declared"
            // }
            if (!fastq_1 && !aligned_bam) {
                error "Sample '${row.sample}' has neither fastq_1 nor aligned_bam declared"
            }

            tuple(row.sample, row.sequencing_type, [fastq_1, fastq_2], aligned_bam,
                  truth_vcf, truth_bed, clair3_platform, clair3_model, dv_model)
        }
        .branch {
            shortread:    it[1] == 'shortread'
            longread_fq:  it[1] == 'longread' && it[3] == []
            longread_bam: it[1] == 'longread' && it[3] != []
        }
        .set { branched }

    // ── Per-technology calling ──
    SHORT_READ_CALLING(
        branched.shortread,
        file(params.reference),
        indexReference.out.index,
        faiReference.out.fai
    )

    LONG_READ_CALLING(
        branched.longread_fq,
        branched.longread_bam,
        file(params.reference),
        faiReference.out.fai
    )

    // ── Tag VCFs with caller name so postprocess can fan out per caller ──
    SHORT_READ_CALLING.out.vcf
        .map { sample, vcf -> tuple(sample, 'deepvariant', vcf) }
        .set { tagged_short }

    LONG_READ_CALLING.out.clair3_vcf
        .map { sample, vcf -> tuple(sample, 'clair3', vcf) }
        .set { tagged_clair3 }

    LONG_READ_CALLING.out.dv_long_vcf
        .map { sample, vcf -> tuple(sample, 'deepvariant_long', vcf) }
        .set { tagged_dv_long }

    LONG_READ_CALLING.out.ensemble_vcf
        .map { sample, vcf -> tuple(sample, 'ensemble', vcf) }
        .set { tagged_ensemble }

    LONG_READ_CALLING.out.ensemble_union_vcf
        .map { sample, vcf -> tuple(sample, 'ensemble_union', vcf) }
        .set { tagged_ensemble_union }

    tagged_short
        .mix(tagged_clair3)
        .mix(tagged_dv_long)
        .mix(tagged_ensemble)
        .mix(tagged_ensemble_union)
        .set { all_vcfs }

    // ── Per-sample truth assets (sample, truth_vcf, truth_tbi, truth_bed) ──
    branched.shortread
        .mix(branched.longread_fq)
        .mix(branched.longread_bam)
        .map { it -> tuple(it[0], it[4], file(it[4].toString() + ".tbi"), it[5]) }
        .set { truth_ch }

    // ── Shared postprocessing ──
    VCF_POSTPROCESS(
        all_vcfs,
        truth_ch,
        file(params.reference),
        faiReference.out.fai
    )
}
