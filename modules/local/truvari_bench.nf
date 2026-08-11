// modules/local/truvari_bench.nf
//
// SV benchmarking against a truth callset. Truvari matches on position, size
// and sequence similarity rather than exact coordinates, since two callers
// rarely place the same breakpoint identically.
//
// Parameters match the manual prototype run (data/sv_benchmark_test) so the
// new numbers are comparable to it:
//   --passonly            score confident calls only
//   --sizemin 50          below this is a small variant, not an SV
//   --sizemax 50000       above this, truth-set representation is unreliable
//   --pctseq/--pctsize 0.7  sequence and size similarity thresholds
//   --refdist 500         maximum breakpoint separation for a match
//
// Truvari refuses to write into an existing directory, hence the unique -o.
process TRUVARI_BENCH {
    tag "${meta.id}.${meta.caller}"
    container 'truvari:5.4.0'
    cpus 2
    memory '8 GB'

    input:
    tuple val(meta), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi), path(truth_bed)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.truvari/summary.json"), emit: summary
    path("${meta.id}.${meta.caller}.truvari"),                              emit: outdir

    script:
    def dup_arg = (params.truvari_dup_to_ins.toString() == 'true') ? '--dup-to-ins' : ''

    """
    truvari bench \\
      -b ${truth_vcf} \\
      -c ${vcf} \\
      -o ${meta.id}.${meta.caller}.truvari \\
      --includebed ${truth_bed} \\
      -f ${reference} \\
      --passonly \\
      --sizemin 50 \\
      --sizemax 50000 \\
      --pctseq ${params.truvari_pctseq} \\
      ${dup_arg} \\
      --pctsize 0.7 \\
      --refdist 500 \\
      --pick single
    """
}
