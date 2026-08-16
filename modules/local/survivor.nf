/*
 * modules/local/survivor.nf
 *
 * SV-SV ensemble merging with SURVIVOR. Takes all SV-caller VCFs for one sample
 * and produces a union (>=1 caller) and an intersection/consensus (>=N callers)
 * merged call set.
 *
 * Input: (meta, [vcf1, vcf2, ...]) — the list of per-caller SV VCFs for a sample,
 * collected upstream by grouping sv_vcfs on meta.id.
 * Emits union and consensus VCFs (bgzipped, indexed).
 *
 * Notes baked in from manual validation:
 *   - SURVIVOR reads PLAIN (uncompressed) VCF — decompress inputs first.
 *   - Inputs should be PASS-filtered before merging (SURVIVOR doesn't set
 *     FILTER=PASS on its output, so downstream benchmarking must NOT use
 *     --passonly on the merged result).
 *   - merge params: <max_dist=1000> <min_callers> <type_match=1> <strand_match=1>
 *     <dyn_dist=0> <min_size=50> <out.vcf>
 *   - Ensemble benefit requires complementary TRUE positives; in this project's
 *     benchmarks the union did not beat the best single caller's recall, while
 *     the intersection maximised precision. The capability is provided; the
 *     interpretation is the precision/recall trade-off, not a strict gain.
 *
 * params.survivor_consensus = minimum number of callers for the consensus set
 *   (default 2). Union is always min_callers=1.
 *
 * Container: locally built survivor:1.0.7 (ubuntu:22.04 + build-essential +
 * zlib1g-dev, SURVIVOR 1.0.7 from source). bcftools image used for sort/index.
 */

process SURVIVOR {
    tag "${meta.id}"
    container 'survivor:1.0.7'

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("${meta.id}.sv_union.vcf.gz"),         emit: union
    tuple val(meta), path("${meta.id}.sv_union.vcf.gz.tbi"),     emit: union_tbi
    tuple val(meta), path("${meta.id}.sv_consensus.vcf.gz"),     emit: consensus
    tuple val(meta), path("${meta.id}.sv_consensus.vcf.gz.tbi"), emit: consensus_tbi

    script:
    def consensus_n = params.survivor_consensus ?: 2
    """
    # SURVIVOR reads plain VCF — decompress any gzipped inputs.
    : > vcf_list.txt
    for v in \$(ls ${vcfs} | sort); do
      case "\$v" in
        *.gz) zcat "\$v" > "\${v%.gz}.plain.vcf"; echo "\${v%.gz}.plain.vcf" >> vcf_list.txt ;;
        *)    echo "\$v" >> vcf_list.txt ;;
      esac
    done

    # Union (>=1 caller) and consensus (>=N callers).
    SURVIVOR merge vcf_list.txt ${params.survivor_max_dist} 1 1 1 0 50 ${meta.id}.sv_union.vcf
    SURVIVOR merge vcf_list.txt ${params.survivor_max_dist} ${consensus_n} 1 1 0 50 ${meta.id}.sv_consensus.vcf
    
    for f in ${meta.id}.sv_union ${meta.id}.sv_consensus; do
        bcftools sort \$f.vcf -Oz -o \$f.vcf.gz
        bcftools index -t \$f.vcf.gz
    done
    """
}
