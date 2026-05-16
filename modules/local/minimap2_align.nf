process MINIMAP2_ALIGN {
    tag "${meta.id}"
    container 'staphb/minimap2:latest'

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: out

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} \\
      -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT" \\
      ${reference} ${reads} > ${meta.id}.sam
    """
}
