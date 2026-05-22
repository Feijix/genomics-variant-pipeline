/*
 * workflows/variant_calling.nf
 *
 * Phase E: data-attribute-driven routing.
 *
 * The samplesheet describes the data (platform, chemistry, basecaller).
 * The routing config (conf/routing.config, nested under params.routing) maps
 * data attributes to caller dispatch. The workflow body reads the routing
 * config, picks per-sample callers, and dispatches one channel per caller
 * (Option A).
 *
 * Adding a new caller:
 *   1. Add a module in modules/local/
 *   2. Add a routing entry in conf/routing.config
 *   3. Add a per-caller branch label and process call in this file
 */

include { BWA_MEM_INDEX  } from '../modules/local/bwa_mem_index.nf'
include { SAMTOOLS_FAIDX } from '../modules/local/samtools_faidx.nf'

include { ALIGNMENT       } from '../subworkflows/local/alignment.nf'
include { VCF_POSTPROCESS } from '../subworkflows/local/vcf_postprocess.nf'

include { DEEPVARIANT     } from '../modules/local/deepvariant.nf'
include { CLAIR3          } from '../modules/local/clair3.nf'
include { BCFTOOLS_ISEC   } from '../modules/local/bcftools_isec.nf'
include { BCFTOOLS_CONCAT } from '../modules/local/bcftools_concat.nf'


// ─── Routing helpers (inline) ────────────────────────────────────────────────

/**
 * Build a routing key from samplesheet attributes. Empty trailing segments
 * are stripped so:
 *   ('illumina', '', '')              -> 'illumina'
 *   ('ont', 'r10.4.1', '')            -> 'ont/r10.4.1'
 *   ('ont', 'r10.4.1', 'dorado_sup')  -> 'ont/r10.4.1/dorado_sup'
 */
def buildRoutingKey(meta) {
    def parts = [meta.platform, meta.chemistry, meta.basecaller].findAll { it }
    return parts.join('/')
}

/**
 * Resolve a routing profile by key, with fallback semantics:
 *   1. Exact match
 *   2. On miss: warn, strip last '/'-segment, retry
 *   3. On match with falls_back_to: recurse to that profile
 *   4. On no match at any level: hard error
 */
def resolveProfile(key, profiles, originalKey = null) {
    def origKey = originalKey ?: key
    def profile = profiles[key]

    if (profile == null) {
        if (!key.contains('/')) {
            error "No routing profile found for '${origKey}' (and no fallback at any level)"
        }
        def parent = key.substring(0, key.lastIndexOf('/'))
        log.warn "No exact routing profile for '${origKey}', falling back to '${parent}'"
        return resolveProfile(parent, profiles, origKey)
    }

    if (profile.containsKey('falls_back_to')) {
        return resolveProfile(profile.falls_back_to, profiles, origKey)
    }
    return profile
}

/**
 * Pick which callers to run for a sample.
 *   1. Per-sample samplesheet `callers` (intersect with profile's tools)
 *   2. CLI --callers all  -> all defaults + alternatives
 *   3. CLI --callers <list>  -> intersect with profile's tools
 *   4. Default -> profile.defaults only
 */
def selectCallers(profile, sampleOverride, cliCallers) {
    def available = (profile.defaults ?: []) + (profile.alternatives ?: [])

    if (sampleOverride) {
        def names = sampleOverride.split(',').collect { it.trim() } as Set
        return available.findAll { it.tool in names }
    }
    if (cliCallers == 'all') {
        return available
    }
    if (cliCallers) {
        def names = cliCallers.split(',').collect { it.trim() } as Set
        return available.findAll { it.tool in names }
    }
    return profile.defaults ?: []
}


// ─── Workflow ────────────────────────────────────────────────────────────────

workflow VARIANT_CALLING {

    // ── Reference indexing ──
    BWA_MEM_INDEX(file(params.reference))
    SAMTOOLS_FAIDX(file(params.reference))

    // ── Load samplesheet, build meta, resolve routing, attach selected callers ──
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id              : row.sample,
                sequencing_type : row.sequencing_type,
                platform        : row.platform        ?: '',
                chemistry       : row.chemistry       ?: '',
                basecaller      : row.basecaller      ?: '',
            ]
            def fastq_1     = row.fastq_1     ? file(row.fastq_1)     : []
            def fastq_2     = row.fastq_2     ? file(row.fastq_2)     : []
            def aligned_bam = row.aligned_bam ? file(row.aligned_bam) : []
            def truth_vcf   = row.truth_vcf   ? file(row.truth_vcf)   : []
            def truth_bed   = row.truth_bed   ? file(row.truth_bed)   : []
            def callers_override = row.callers ?: ''

            if (!fastq_1 && !aligned_bam) {
                error "Sample '${meta.id}' has neither fastq_1 nor aligned_bam declared"
            }

            meta.routing_key = buildRoutingKey(meta)
            def profile = resolveProfile(meta.routing_key, params.routing.profiles)
            meta.selected_callers = selectCallers(profile, callers_override, params.callers)

            if (meta.selected_callers.isEmpty()) {
                error "Sample '${meta.id}' (routing '${meta.routing_key}') resolved to zero callers"
            }

            tuple(meta, [fastq_1, fastq_2], aligned_bam, truth_vcf, truth_bed)
        }
        .set { samples }

    // ── Alignment (all samples through one subworkflow) ──
    ALIGNMENT(
        samples,
        file(params.reference),
        BWA_MEM_INDEX.out.index,
        SAMTOOLS_FAIDX.out.fai
    )

    // ── Fan each aligned sample to one item per selected caller ──
    ALIGNMENT.out.bams
        .flatMap { meta, bam, bai ->
            meta.selected_callers.collect { caller_spec ->
                def enriched = meta + caller_spec.params + [caller: caller_spec.tool]
                tuple(enriched, bam, bai)
            }
        }
        .branch { meta, _bam, _bai ->
            deepvariant: meta.caller == 'deepvariant'
            clair3:      meta.caller == 'clair3'
        }
        .set { per_caller }

    // ── Per-caller invocations (Option A dispatch) ──
    DEEPVARIANT(per_caller.deepvariant, file(params.reference), SAMTOOLS_FAIDX.out.fai)
    CLAIR3     (per_caller.clair3,      file(params.reference), SAMTOOLS_FAIDX.out.fai)

    // ── Tag each caller's output for postprocess uniformity ──
    DEEPVARIANT.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'deepvariant'], vcf) }
        .set { tagged_dv }

    CLAIR3.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'clair3'], vcf) }
        .set { tagged_clair3 }

    tagged_dv.mix(tagged_clair3).set { primary_vcfs }

    // ── Ensemble formation: samples with 2+ callers (assumes exactly 2 today) ──
    primary_vcfs
        .map { meta, vcf -> tuple(meta.id, meta, vcf) }
        .groupTuple()
        .filter { _id, _metas, vcfs -> vcfs.size() >= 2 }
        .map { id, metas, vcfs ->
            def base_meta = metas[0].findAll { k, _v ->
                k in ['id','sequencing_type','platform','chemistry','basecaller','routing_key']
            }
            tuple(base_meta, vcfs[0], vcfs[1])
        }
        .set { ensemble_pair }

    if (params.routing.ensembles.intersect) {
        BCFTOOLS_ISEC(ensemble_pair)
        BCFTOOLS_ISEC.out.vcf
            .map { meta, vcf -> tuple(meta + [caller: 'ensemble'], vcf) }
            .set { tagged_isec }
    } else {
        Channel.empty().set { tagged_isec }
    }

    if (params.routing.ensembles.union) {
        BCFTOOLS_CONCAT(ensemble_pair, file(params.reference))
        BCFTOOLS_CONCAT.out.vcf
            .map { meta, vcf -> tuple(meta + [caller: 'ensemble_union'], vcf) }
            .set { tagged_concat }
    } else {
        Channel.empty().set { tagged_concat }
    }

    // ── All VCFs (primary + ensemble) into postprocess ──
    primary_vcfs
        .mix(tagged_isec)
        .mix(tagged_concat)
        .set { all_vcfs }

    // ── Per-sample truth assets ──
    samples
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
