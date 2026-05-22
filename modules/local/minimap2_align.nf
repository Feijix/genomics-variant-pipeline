process MINIMAP2_ALIGN {
    tag "${meta.id}"
    container 'staphb/minimap2:latest'

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: out

    script:
    // ── Preset selection by platform ──
    //   ont    -> map-ont   (default ONT parameters)
    //   pacbio -> map-hifi  (tuned for HiFi's lower error rate)
    //   others -> map-ont   (defensive default; warn at workflow level if reached)
    def preset = meta.platform == 'pacbio' ? 'map-hifi' : 'map-ont'
    def pl_tag = meta.platform == 'pacbio' ? 'PACBIO'  : 'ONT'
    """
    minimap2 -ax ${preset} -t ${task.cpus} \\
      -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:${pl_tag}" \\
      ${reference} ${reads} > ${meta.id}.sam
    """
}
