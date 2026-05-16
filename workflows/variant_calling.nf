/*
 * workflows/variant_calling.nf
 *
 * Main orchestration. Loads the samplesheet, builds a meta map per row,
 * branches by sequencing technology and input flavour, dispatches to the
 * appropriate calling subworkflow, and feeds the unified tagged stream into
 * postprocessing.
 */

include { BWA_MEM_INDEX      } from '../modules/local/bwa_mem_index.nf'
include { SAMTOOLS_FAIDX     } from '../modules/local/samtools_faidx.nf'

include { SHORT_READ_CALLING } from '../subworkflows/local/short_read_calling.nf'
include { LONG_READ_CALLING  } from '../subworkflows/local/long_read_calling.nf'
include { VCF_POSTPROCESS    } from '../subworkflows/local/vcf_postprocess.nf'

workflow VARIANT_CALLING {

    // ── Reference indexing ──
    BWA_MEM_INDEX(file(params.reference))
    SAMTOOLS_FAIDX(file(params.reference))

    // ── Load samplesheet, build meta + paths, branch ──
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id              : row.sample,
                sequencing_type : row.sequencing_type,
                clair3_platform : row.clair3_platform ?: '',
                clair3_model    : row.clair3_model    ?: '',
                dv_model        : row.dv_model        ?: '',
            ]
            def fastq_1     = row.fastq_1     ? file(row.fastq_1)     : []
            def fastq_2     = row.fastq_2     ? file(row.fastq_2)     : []
            def aligned_bam = row.aligned_bam ? file(row.aligned_bam) : []
            def truth_vcf   = row.truth_vcf   ? file(row.truth_vcf)   : []
            def truth_bed   = row.truth_bed   ? file(row.truth_bed)   : []

            if (!fastq_1 && !aligned_bam) {
                error "Sample '${meta.id}' has neither fastq_1 nor aligned_bam declared"
            }

            tuple(meta, [fastq_1, fastq_2], aligned_bam, truth_vcf, truth_bed)
        }
        .branch { meta, _reads, aligned_bam, _vcf, _bed ->
            shortread:    meta.sequencing_type == 'shortread'
            longread_fq:  meta.sequencing_type == 'longread' && aligned_bam == []
            longread_bam: meta.sequencing_type == 'longread' && aligned_bam != []
        }
        .set { branched }

    // ── Per-technology calling ──
    SHORT_READ_CALLING(
        branched.shortread,
        file(params.reference),
        BWA_MEM_INDEX.out.index,
        SAMTOOLS_FAIDX.out.fai
    )

    LONG_READ_CALLING(
        branched.longread_fq,
        branched.longread_bam,
        file(params.reference),
        SAMTOOLS_FAIDX.out.fai
    )

    // ── Mix all caller outputs (each subworkflow tags with meta.caller) ──
    SHORT_READ_CALLING.out.vcf
        .mix(LONG_READ_CALLING.out.clair3_vcf)
        .mix(LONG_READ_CALLING.out.dv_long_vcf)
        .mix(LONG_READ_CALLING.out.ensemble_vcf)
        .mix(LONG_READ_CALLING.out.ensemble_union_vcf)
        .set { all_vcfs }

    // ── Per-sample truth assets (sample_id, truth_vcf, truth_tbi, truth_bed) ──
    branched.shortread
        .mix(branched.longread_fq)
        .mix(branched.longread_bam)
        .map { meta, _reads, _bam, truth_vcf, truth_bed ->
            tuple(meta.id, truth_vcf, file(truth_vcf.toString() + ".tbi"), truth_bed)
        }
        .set { truth_ch }

    // ── Shared postprocessing ──
    VCF_POSTPROCESS(
        all_vcfs,
        truth_ch,
        file(params.reference),
        SAMTOOLS_FAIDX.out.fai
    )
}
