process SAMTOOLS_SORT {
    tag "${meta.id}"
    container 'staphb/samtools:latest'

    input:
    tuple val(meta), path(sam_file)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: out

    script:
    """
    samtools sort -o ${meta.id}.sorted.bam ${sam_file}
    samtools index ${meta.id}.sorted.bam
    """
}
