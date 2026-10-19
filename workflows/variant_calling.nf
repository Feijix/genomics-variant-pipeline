/*
 * workflows/variant_calling.nf
 * 
 * The workflow body reads the routing
 * config, picks per-sample callers, and dispatches one channel per caller
 *
 * To add a new caller:
 *   1. Add a module in modules/local/
 *   2. Add a routing entry in conf/routing.config
 *   3. Add a per-caller branch label and process call in this file
 */

include { BWA_MEM_INDEX  } from '../modules/local/bwa_mem_index.nf'
include { SAMTOOLS_FAIDX } from '../modules/local/samtools_faidx.nf'

include { ALIGNMENT       } from '../subworkflows/local/alignment.nf'
include { VCF_POSTPROCESS } from '../subworkflows/local/vcf_postprocess.nf'
include { QUALITY_CONTROL } from '../subworkflows/local/quality_control.nf'            

include { DEEPVARIANT     } from '../modules/local/deepvariant.nf'
include { CLAIR3          } from '../modules/local/clair3.nf'
include { CLAIRS          } from '../modules/local/clairs.nf'
include { SNIFFLES2       } from '../modules/local/sniffles2.nf'
include { CUTESV          } from '../modules/local/cutesv.nf'
include { MANTA           } from '../modules/local/manta.nf'
include { DELLY           } from '../modules/local/delly.nf'
include { CONVERT_DELLY   } from '../modules/local/convert_delly.nf'
include { SMOOVE          } from '../modules/local/smoove.nf'
include { GRIDSS          } from '../modules/local/gridss.nf'
include { CONVERT_GRIDSS  } from '../modules/local/convert_gridss.nf'
include { SURVIVOR        } from '../modules/local/survivor.nf'
include { BCFTOOLS_ISEC   } from '../modules/local/bcftools_isec.nf'
include { BCFTOOLS_CONCAT } from '../modules/local/bcftools_concat.nf'
include { BCFTOOLS_VIEW_SV } from '../modules/local/bcftools_view_sv.nf'
include { TRUVARI_BENCH   } from '../modules/local/truvari_bench.nf'
include { BEACON_SV2BFF   } from '../modules/local/beacon_sv2bff.nf'
include { GATK_CREATESEQUENCEDICTIONARY } from '../modules/local/gatk_createsequencedictionary.nf'
include { GATK_MARKDUPLICATES           } from '../modules/local/gatk_markduplicates.nf'
include { GATK_MUTECT2                  } from '../modules/local/gatk_mutect2.nf'
include { GATK_FILTERMUTECTCALLS        } from '../modules/local/gatk_filtermutectcalls.nf'
include { MANTA_SOMATIC                 } from '../modules/local/manta_somatic.nf'
include { STRELKA2_SOMATIC              } from '../modules/local/strelka2_somatic.nf'
include { BCFTOOLS_CONCAT_SOMATIC       } from '../modules/local/bcftools_concat_somatic.nf'
include { BCFTOOLS_CONCAT_SOMATIC as BCFTOOLS_CONCAT_CLAIRS } from '../modules/local/bcftools_concat_somatic.nf'
include { INFER_READ_LENGTH		} from '../modules/local/infer_read_length.nf'
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

def selectCallers(profile, sampleOverride, cliCallers, analysis, toolRegistry) {
    def isGermline      = (analysis == 'germline')
    def allDefaults     = (isGermline ? profile.defaults     : profile.somatic_defaults)     ?: []
    def allAlternatives = (isGermline ? profile.alternatives : profile.somatic_alternatives) ?: []

    // Everything the profile declares, before preconditions are applied. A name
    // absent from this is an override; a name present here but absent from
    // `available` was withdrawn by a precondition.
    def declared = ((profile.defaults ?: []) + (profile.alternatives ?: []) +
                    (profile.somatic_defaults ?: []) + (profile.somatic_alternatives ?: []))

    def defaults     = allDefaults
    def alternatives = allAlternatives
    if (analysis == 'somatic_tumoronly') {
        defaults     = defaults.findAll     { !it.requires_normal }
        alternatives = alternatives.findAll { !it.requires_normal }
    }
    def available = defaults + alternatives

    def names  = null
    def source = null
    if (sampleOverride) {
        names  = sampleOverride.split(/[;,]/).collect { it.trim() } as Set
        source = 'samplesheet'
    } else if (cliCallers == 'all') {
        return available
    } else if (cliCallers) {
        names  = cliCallers.split(',').collect { it.trim() } as Set
        source = '--callers'
    } else {
        return defaults
    }

    def fromProfile = available.findAll { it.tool in names }

    // Named tools the profile declares but a precondition withdrew (requires_normal on a
    // tumour-only sample) are dropped; report them instead of dropping them silently.
    def withdrawn = names.findAll { n -> (n in declared*.tool) && !(n in available*.tool) }
    if (withdrawn) {
        log.warn "Callers ${withdrawn.join(', ')} (${source}) withdrawn: they require a matched " +
                 "normal, which this sample lacks."
    }

    // Named tools the profile does not declare: an explicit override that departs from
    // attribute-derived selection.
    def forced = names - declared*.tool

    // A name absent from the tool registry is not a known tool (often a typo); without this
    // check it would be dispatched to no process and dropped silently.
    def unknown = forced.findAll { !(it in toolRegistry) }
    if (unknown) {
        error "Callers ${unknown.join(', ')} (${source}) are not defined in routing.tools."
    }

    // Tools that need a profile-supplied model cannot run without one, so forcing them
    // outside a profile is refused before they fail in their container.
    def unrunnable = forced.findAll { toolRegistry[it].requires_model }
    if (unrunnable) {
        error "Callers ${unrunnable.join(', ')} require a model that only a matching profile " +
              "supplies; add them to the profile, or route the sample to a profile that " +
              "includes them."
    }

    // The remaining forced tools run with default parameters.
    def runnable = forced - unrunnable
    if (runnable) {
        log.warn "Callers ${runnable.join(', ')} (${source}) are not declared in this profile and " +
                 "were forced; they run with default parameters, which may be inadequate for " +
                 "this data."
    }
    def fromForced = runnable.collect { t ->
        [tool: t, params: [:], variant_class: toolRegistry[t].variant_class]
    }

    return fromProfile + fromForced
}

workflow VARIANT_CALLING {

    // ── Reference indexing ──
    SAMTOOLS_FAIDX(file(params.reference))

    // ── Load samplesheet, build meta, resolve routing, attach selected callers ──
    // ── 1. Parse rows into meta. No caller selection yet. ──
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id              : row.sample,
                subject         : row.subject ?: row.sample,
                status          : row.status  ?: 'normal',
                sequencing_type : row.sequencing_type,
                platform        : row.platform   ?: '',
                chemistry       : row.chemistry  ?: '',
                basecaller      : row.basecaller ?: '',
                callers_override: row.callers    ?: '',
            ]
            def fastq_1     = row.fastq_1     ? file(row.fastq_1)     : []
            def fastq_2     = row.fastq_2     ? file(row.fastq_2)     : []
            def aligned_bam = row.aligned_bam ? file(row.aligned_bam) : []
            def truth_vcf   = row.truth_vcf   ? file(row.truth_vcf)   : []
            def truth_bed   = row.truth_bed   ? file(row.truth_bed)   : []

            if (!fastq_1 && !aligned_bam)
                error "Sample '${meta.id}' has neither fastq_1 nor aligned_bam declared"

            meta.routing_key = buildRoutingKey(meta)
            tuple(meta, [fastq_1, fastq_2], aligned_bam, truth_vcf, truth_bed)
        }
        .set { rows }

        // ── 1b. Read-length inference for samples with no declared platform. ──
    rows
        .branch { meta, _r, _b, _tv, _tb ->
            declared:   meta.routing_key != ''
            undeclared: true
        }
        .set { by_decl }

    // run inference on the undeclared branch
    INFER_READ_LENGTH(
        by_decl.undeclared.map { meta, r, b, _tv, _tb -> tuple(meta, r.findAll { it }, b) }
    )

    // set platform from the inferred family, rebuild the key, restore the full tuple
    INFER_READ_LENGTH.out.inferred
        .map { meta, family, mean ->
            meta.platform = family
            meta.sequencing_type = family
            meta.platform_inferred = true
            meta.routing_key = family
            log.info "Inferred: ${meta.id} -> family '${family}' from mean read length ${mean} bases"
            tuple(meta.id, meta)
        }
        .set { inferred_meta }

    // join the inferred meta back onto the original undeclared rows to recover reads/truth
    by_decl.undeclared
        .map { meta, r, b, tv, tb -> tuple(meta.id, r, b, tv, tb) }
        .combine(inferred_meta, by: 0)
        .map { _id, r, b, tv, tb, meta -> tuple(meta, r, b, tv, tb) }
        .set { undeclared_routed }

    // recombine declared + inferred-routed samples
    by_decl.declared
        .mix(undeclared_routed)
        .set { rows_routed }

    // ── 2. Metadata-only pre-pass: validate each subject, emit pairing facts. ──
    rows_routed
        .map { meta, _r, _b, _tv, _tb ->
            tuple(meta.subject, [id: meta.id, status: meta.status, seqtype: meta.sequencing_type])
        }
        .groupTuple()
        .map { subj, entries ->
            def bad = entries.findAll { !(it.status in ['tumor', 'normal']) }
            if (bad)
                error "Subject '${subj}': invalid status '${bad[0].status}' on sample '${bad[0].id}'."

            def normals = entries.findAll { it.status == 'normal' }
            def tumors  = entries.findAll { it.status == 'tumor'  }

            if (tumors && normals.size() > 1)
                error "Subject '${subj}': ${normals.size()} normals with ${tumors.size()} tumor(s). Paired somatic requires exactly one normal."

            if (tumors && normals.size() == 1) {
                def n = normals[0]
                tumors.findAll { it.seqtype != n.seqtype }.each { t ->
                    error "Subject '${subj}': cross-technology T/N pairing not supported (tumor '${t.id}' ${t.seqtype} vs normal '${n.id}' ${n.seqtype})."
                }
            }
            tuple(subj, normals.size() == 1, tumors.size() > 0)
        }
        .set { subject_facts }

    // ── 3. Join facts back, resolve analysis, select callers. ──
    rows_routed
        .map { meta, r, b, tv, tb -> tuple(meta.subject, meta, r, b, tv, tb) }
        .combine(subject_facts, by: 0)
        .map { _subj, meta, r, b, tv, tb, hasNormal, hasTumor ->
            meta.analysis = meta.status == 'normal' ? 'germline'
                        : hasNormal               ? 'somatic_paired'
                        :                           'somatic_tumoronly'
            meta.is_matched_normal = (meta.status == 'normal' && hasTumor)

            def profile = resolveProfile(meta.routing_key, params.routing.profiles)
            meta.selected_callers = selectCallers(profile, meta.callers_override, params.callers, meta.analysis, params.routing.tools)
	    log.info "Routing: ${meta.id} -> profile '${meta.routing_key}', analysis '${meta.analysis}', callers ${meta.selected_callers*.tool}"
	    
            if (meta.selected_callers.isEmpty()) {
                if (meta.status == 'normal' && meta.is_matched_normal) {
                    log.info "Sample '${meta.id}' is a matched normal; callers run on the tumour."
                } else {
                    error "Sample '${meta.id}' (routing '${meta.routing_key}', analysis '${meta.analysis}') resolved to zero callers"
                }
            }
                
            def tools = meta.selected_callers*.tool
	    meta.selected_callers.each { e ->
                (params.routing.tools[e.tool]?.requires ?: []).each { dep ->
                    if (!(dep in tools))
                        error "Sample '${meta.id}': ${e.tool} requires ${dep}, which is not selected."
                }
            }

            tuple(meta, r, b, tv, tb)
        }
        .set { samples }

    if (params.skip_bwa_index.toString() == 'true') {
        bwa_index = Channel.empty()
    } else {
        BWA_MEM_INDEX(file(params.reference))
        bwa_index = BWA_MEM_INDEX.out.index
    }

    // ── Alignment (all samples through one subworkflow) ──
    ALIGNMENT(
        samples,
        file(params.reference),
        bwa_index,
        SAMTOOLS_FAIDX.out.fai
    )

    // ── Fan each aligned sample to one item per selected caller ──
    ALIGNMENT.out.bams
        .branch { meta, _bam, _bai ->
            normal: meta.status == 'normal'
            tumor : true
        }
        .set { by_status }

    // Germline fanout: UNCHANGED body, now fed only by normals.
    by_status.normal
        .flatMap { meta, bam, bai ->
            meta.selected_callers.collect { caller_spec ->
                def enriched = meta + caller_spec.params + [caller: caller_spec.tool, variant_class: caller_spec.variant_class ?: 'small']
                tuple(enriched, bam, bai)
            }
        }
        .branch { meta, _bam, _bai ->
            deepvariant: meta.caller == 'deepvariant'
            clair3:      meta.caller == 'clair3'
            sniffles2:   meta.caller == 'sniffles2'
            cutesv:      meta.caller == 'cutesv'
            manta:       meta.caller == 'manta'
            delly:       meta.caller == 'delly'
            smoove:      meta.caller == 'smoove'
            gridss:      meta.caller == 'gridss'
            unmatched:   true
        }
        .set { per_caller }
 
    per_caller.unmatched.map { meta, _b, _i ->
        error "Caller '${meta.caller}' selected for '${meta.id}' has no dispatch label."
    }

    // ── Per-caller invocations (Option A dispatch) ──
    DEEPVARIANT(per_caller.deepvariant, file(params.reference), SAMTOOLS_FAIDX.out.fai)
    CLAIR3     (per_caller.clair3,      file(params.reference), SAMTOOLS_FAIDX.out.fai)
    SNIFFLES2  (per_caller.sniffles2,   file(params.reference), SAMTOOLS_FAIDX.out.fai)
    CUTESV     (per_caller.cutesv,      file(params.reference), SAMTOOLS_FAIDX.out.fai)
    MANTA      (per_caller.manta,       file(params.reference), SAMTOOLS_FAIDX.out.fai)
    DELLY      (per_caller.delly,       file(params.reference), SAMTOOLS_FAIDX.out.fai)
    CONVERT_DELLY(DELLY.out.bcf.join(DELLY.out.csi))
    SMOOVE     (per_caller.smoove, file(params.reference), SAMTOOLS_FAIDX.out.fai)
    GRIDSS(per_caller.gridss, file(params.reference), SAMTOOLS_FAIDX.out.fai)
    CONVERT_GRIDSS(GRIDSS.out.bnd, SAMTOOLS_FAIDX.out.fai)
    DEEPVARIANT.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'deepvariant'], vcf) }
        .set { tagged_dv }

    CLAIR3.out.out
        .map { meta, vcf -> tuple(meta + [caller: 'clair3'], vcf) }
        .set { tagged_clair3 }
    
    SNIFFLES2.out.vcf
        .map { meta, vcf -> tuple(meta + [caller: 'sniffles2'], vcf) }
        .mix(CUTESV.out.vcf.map { meta, vcf -> tuple(meta + [caller: 'cutesv'], vcf) })
        .mix(MANTA.out.vcf.map  { meta, vcf -> tuple(meta + [caller: 'manta'],  vcf) })
        .mix(CONVERT_DELLY.out.vcf.map { meta, vcf -> tuple(meta + [caller: 'delly'], vcf) })
        .mix(SMOOVE.out.vcf.map         { meta, vcf -> tuple(meta + [caller: 'smoove'], vcf) })
        .mix(CONVERT_GRIDSS.out.vcf.map { meta, vcf -> tuple(meta + [caller: 'gridss'], vcf) })
        .set { sv_vcfs }

    

    tagged_dv.mix(tagged_clair3).set { primary_vcfs }

    // ── Ensemble formation: samples with 2+ callers (assumes exactly 2 today) ──
    primary_vcfs
        .map { meta, vcf -> tuple(meta.id, meta, vcf) }
        .groupTuple()
        .filter { _id, _metas, vcfs -> vcfs.size() >= 2 }
        .map { id, metas, vcfs ->
            def priority = ['deepvariant', 'clair3']
            def pairs = [metas, vcfs].transpose().sort { pair -> priority.indexOf(pair[0].caller) }
            def base_meta = pairs[0][0].findAll { k, _v ->
                k in ['id','sequencing_type','platform','chemistry','basecaller','routing_key']
            }
            tuple(base_meta, pairs[0][1], pairs[1][1])
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

    // ── MarkDuplicates on tumors + matched normals only. ──
    GATK_MARKDUPLICATES(
        by_status.tumor.mix(
            by_status.normal.filter { meta, _bam, _bai -> meta.is_matched_normal }
        )
    )

    GATK_MARKDUPLICATES.out.bam
        .branch { meta, _bam, _bai ->
            tumor : meta.status == 'tumor'
            normal: true
        }
        .set { md }

    // ── Pair by subject (Rule 2: N tumors × 1 normal → N pairs, falls out of combine). ──
    md.tumor
        .map { meta, bam, bai -> tuple(meta.subject, meta, bam, bai) }
        .combine(md.normal.map { meta, bam, bai -> tuple(meta.subject, meta.id, bam, bai) }, by: 0)
        .map { _subj, tmeta, tbam, tbai, nid, nbam, nbai ->
            tuple(tmeta + [normal_id: nid], tbam, tbai, nbam, nbai)
        }
        .set { tn_pairs }

    // ── Tumor-only: empty normal paths trigger the '[]' idiom in Mutect2. ──
    md.tumor
        .filter { meta, _bam, _bai -> meta.analysis == 'somatic_tumoronly' }
        .map    { meta, bam, bai -> tuple(meta, bam, bai, [], []) }
        .set { tumor_only }

    // ── Somatic per-caller fanout. ──
    tn_pairs.mix(tumor_only)
        .flatMap { meta, tbam, tbai, nbam, nbai ->
            meta.selected_callers.collect { spec ->
                def enriched = meta + spec.params + [caller: spec.tool, variant_class: spec.variant_class ?: 'small']
                tuple(enriched, tbam, tbai, nbam, nbai)
            }
        }
        .branch { meta, _tb, _ti, _nb, _ni ->
            mutect2      : meta.caller == 'mutect2'
            strelka2     : meta.caller == 'strelka2'
            manta_somatic: meta.caller == 'manta_somatic'
            clairs       : meta.caller == 'clairs'
            unmatched    : true
        }
        .set { per_somatic }

    per_somatic.unmatched.map { meta, _tb, _ti, _nb, _ni ->
        error "Caller '${meta.caller}' selected for '${meta.id}' has no dispatch label."
    }

    // ── Reference dict + calling-reference handles. ──
    GATK_CREATESEQUENCEDICTIONARY(file(params.reference))
    def ref_calling = file(params.reference)
    def fai_calling = SAMTOOLS_FAIDX.out.fai

    def germline_resource = params.germline_resource
        ? [ file(params.germline_resource), file("${params.germline_resource}.tbi") ]
        : [ [], [] ]

    def pon = params.pon
        ? [ file(params.pon), file("${params.pon}.tbi") ]
        : [ [], [] ]

    // ── Mutect2 → FilterMutectCalls. ──
    GATK_MUTECT2(
        per_somatic.mutect2,
        ref_calling, fai_calling, GATK_CREATESEQUENCEDICTIONARY.out.dict,
        germline_resource,
        pon
    )
    GATK_FILTERMUTECTCALLS(
        GATK_MUTECT2.out.vcf.join(GATK_MUTECT2.out.tbi).join(GATK_MUTECT2.out.stats),
        ref_calling, fai_calling, GATK_CREATESEQUENCEDICTIONARY.out.dict
    )

    def call_regions = params.call_regions
        ? [ file(params.call_regions), file("${params.call_regions}.tbi") ]
        : [ [], [] ]
        
    // ── Manta somatic (paired only; per_somatic.manta_somatic is empty for tumor-only). ──
    MANTA_SOMATIC(per_somatic.manta_somatic, file(params.reference), SAMTOOLS_FAIDX.out.fai, call_regions)

    // ── Strelka2: serial on Manta's indel candidates. Join by sample id. ──
    per_somatic.strelka2
        .map { meta, tb, ti, nb, ni -> tuple(meta.id, meta, tb, ti, nb, ni) }
        .join(
            MANTA_SOMATIC.out.indel_candidates
                .join(MANTA_SOMATIC.out.indel_candidates_tbi)
                .map { meta, vcf, tbi -> tuple(meta.id, vcf, tbi) }
        )
        .map { _id, meta, tb, ti, nb, ni, cvcf, ctbi -> tuple(meta, tb, ti, nb, ni, cvcf, ctbi) }
        .set { strelka_in }

    STRELKA2_SOMATIC(strelka_in, ref_calling, fai_calling, call_regions)

    // ── Merge Strelka2 SNVs + indels (aliased CONCAT). ──
    BCFTOOLS_CONCAT_SOMATIC(
        STRELKA2_SOMATIC.out.snvs.join(STRELKA2_SOMATIC.out.indels)
            .map { meta, snvs, indels ->
                tuple(meta + [caller: 'strelka2'], snvs, indels) },
        ref_calling
    )

    CLAIRS(
        per_somatic.clairs,
        file(params.reference),
        SAMTOOLS_FAIDX.out.fai
    )

    BCFTOOLS_CONCAT_CLAIRS(
        CLAIRS.out.snvs.join(CLAIRS.out.indels)
            .map { meta, snvs, indels -> tuple(meta + [caller: 'clairs'], snvs, indels) },
        ref_calling
    )

    // ── Collect somatic small variants (kept OUT of the germline ensemble). ──
    GATK_FILTERMUTECTCALLS.out.vcf
        .map { meta, vcf -> tuple(meta + [caller: 'mutect2'], vcf) }
        .mix(BCFTOOLS_CONCAT_SOMATIC.out.vcf.map { meta, vcf -> tuple(meta + [caller: 'strelka2'], vcf) })
        .mix(BCFTOOLS_CONCAT_CLAIRS.out.vcf)
        .set { somatic_small_vcfs }

    sv_vcfs
        .map { meta, vcf -> tuple(meta.id, vcf) }
        .groupTuple()                                  // [id, [vcf1, vcf2, ...]]
        .map { id, vcfs -> tuple([id: id], vcfs) }
        .set { sv_grouped }
    
    sv_vcfs
        .mix(MANTA_SOMATIC.out.sv.map { meta, vcf -> tuple(meta + [caller: 'manta_somatic'], vcf) })
        .set { sv_vcfs_raw }
	
    BCFTOOLS_VIEW_SV(sv_vcfs_raw)
    BCFTOOLS_VIEW_SV.out.vcf.set { sv_vcfs }

    // SURVIVOR guard: a tumor's only SV caller is manta_somatic → group of one.
    sv_vcfs
        .map { meta, vcf -> tuple(meta.id, vcf) }
        .groupTuple()
        .filter { _id, vcfs -> vcfs.size() >= 2 }     // ← NEW LINE
        .map { id, vcfs -> tuple([id: id], vcfs) }
        .set { sv_grouped }
        
    SURVIVOR(sv_grouped)

    // ── All VCFs (primary + ensemble) into postprocess ──
    primary_vcfs
        .mix(tagged_isec)
        .mix(tagged_concat)
        .mix(somatic_small_vcfs)
        .set { all_vcfs }

    // ── Per-sample truth assets ──
    samples
        .map { meta, _reads, _bam, truth_vcf, truth_bed ->
            def tbi = truth_vcf ? file(truth_vcf.toString() + ".tbi") : []
            tuple(meta.id, truth_vcf ?: [], tbi, truth_bed ?: [])
        }
        .set { truth_ch }

    // ── SV benchmarking ──
    // Runs outside VCF_POSTPROCESS: SVs need Truvari rather than hap.py, and the
    // normalise/annotate chain does not apply to breakend records.
    BCFTOOLS_VIEW_SV.out.vcf
        .join(BCFTOOLS_VIEW_SV.out.tbi)
        .mix(
            SURVIVOR.out.union.join(SURVIVOR.out.union_tbi)
                .map { meta, vcf, tbi -> tuple(meta + [caller: 'survivor_union'], vcf, tbi) },
            SURVIVOR.out.consensus.join(SURVIVOR.out.consensus_tbi)
                .map { meta, vcf, tbi -> tuple(meta + [caller: 'survivor_consensus'], vcf, tbi) }
        )
        .map { meta, vcf, tbi -> tuple(meta.id, meta, vcf, tbi) }
        .combine(truth_ch, by: 0)
        .map { id, meta, vcf, tbi, tv, tt, tb -> tuple(meta, vcf, tbi, tv, tt, tb) }
        .filter { meta, vcf, tbi, tv, tt, tb -> tv && tv.exists() && tb && tb.exists() }
        .set { sv_bench_input }

    TRUVARI_BENCH(
        sv_bench_input,
        file(params.reference),
        SAMTOOLS_FAIDX.out.fai
    )

    // ── Shared postprocessing ──
    VCF_POSTPROCESS(
        all_vcfs,
        truth_ch,
        file(params.reference),
        SAMTOOLS_FAIDX.out.fai
    )

    BEACON_SV2BFF(sv_vcfs, file("${projectDir}/scripts/sv_to_bff.py"))

    QUALITY_CONTROL(
        ALIGNMENT.out.bams,
        all_vcfs,
        GATK_MARKDUPLICATES.out.metrics
    )
}
