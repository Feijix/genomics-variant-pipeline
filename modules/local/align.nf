process align {
    tag "${sample_id}"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    publishDir "${params.outdir}/bam", mode: 'copy'
    memory '32 GB'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    path(index)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: out

    script:
    """
    bwa-mem2 mem -t ${task.cpus} \\
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
      ${reference} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}
