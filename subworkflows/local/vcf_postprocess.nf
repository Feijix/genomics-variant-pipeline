/*
 * subworkflows/local/vcf_postprocess.nf
 *
 * Uniform per-caller postprocess applied to every caller's output:
 *   BCFTOOLS_VIEW (FILTER=PASS) -> BCFTOOLS_NORM (left-align, split multiallelics)
 *   -> ENSEMBLVEP_VEP (annotate) -> HAPPY_HAPPY (benchmark vs GIAB truth)
 *   -> BEACON_VCF2BFF (annotated VCF -> Beacon v2 genomicVariations BFF)
 *
 * Benchmark joins by meta.id and only runs when truth_vcf + truth_bed present.
 * Beacon conversion runs on every annotated VCF regardless of truth presence.
 */

include { BCFTOOLS_VIEW   } from '../../modules/local/bcftools_view.nf'
include { BCFTOOLS_NORM   } from '../../modules/local/bcftools_norm.nf'
include { ENSEMBLVEP_VEP  } from '../../modules/local/ensemblvep_vep.nf'
include { HAPPY_HAPPY     } from '../../modules/local/happy_happy.nf'
include { SOMPY_SOMPY     } from '../../modules/local/sompy_sompy.nf'
include { BEACON_VCF2BFF  } from '../../modules/local/beacon_vcf2bff.nf'

workflow VCF_POSTPROCESS {
    take:
    tagged_vcfs    // (meta+caller, vcf)
    truth_ch       // (sample_id, truth_vcf, truth_tbi, truth_bed)
    reference
    fai

    main:
    BCFTOOLS_VIEW(tagged_vcfs)
    BCFTOOLS_NORM(BCFTOOLS_VIEW.out.vcf, reference)
    ENSEMBLVEP_VEP(BCFTOOLS_NORM.out.vcf, reference, fai)

    // ── Benchmark: join annotated VCFs with truth assets by sample id ──
    BCFTOOLS_NORM.out.vcf
        .map { meta, vcf -> tuple(meta.id, meta, vcf) }
        .combine(truth_ch, by: 0)
        .map { sample_id, meta, vcf, truth_vcf, truth_tbi, truth_bed ->
            tuple(meta, vcf, truth_vcf, truth_tbi, truth_bed) }
        .filter { meta, vcf, truth_vcf, truth_tbi, truth_bed ->
            truth_vcf && truth_vcf.toString() != '[]' && truth_vcf.exists() &&
            truth_bed && truth_bed.exists()
        }
        .branch { row ->
            def meta = row[0]
            somatic:  meta.analysis in ['somatic_paired','somatic_tumoronly']
            germline: true
        }
        .set { bench }

    HAPPY_HAPPY(bench.germline, reference, fai)
    SOMPY_SOMPY(bench.somatic,  reference, fai)


    // ── Beacon v2 BFF conversion: every annotated VCF ──
    BEACON_VCF2BFF(ENSEMBLVEP_VEP.out.vcf)

    emit:
    annotated_vcf     = ENSEMBLVEP_VEP.out.vcf
    benchmark_summary = HAPPY_HAPPY.out.summary
    somatic_summary   = SOMPY_SOMPY.out.stats
    beacon_bff        = BEACON_VCF2BFF.out.bff
}
