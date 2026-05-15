process align {
    tag "${meta.id}"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    memory '32 GB'

    input:
    tuple val(meta), path(reads)
    path(reference)
    path(index)

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: out

    script:
    """
    bwa-mem2 mem -t ${task.cpus} \\
      -R "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA" \\
      ${reference} ${reads[0]} ${reads[1]} > ${meta.id}.sam
    """
}
