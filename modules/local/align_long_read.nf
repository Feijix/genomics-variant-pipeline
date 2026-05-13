process alignLongRead {
    tag "${sample_id}"
    container 'staphb/minimap2:latest'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: out

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} \\
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ONT" \\
      ${reference} ${reads} > ${sample_id}.sam
    """
}
