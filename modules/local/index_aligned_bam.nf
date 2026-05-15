process indexAlignedBam {
    tag "${meta.id}"
    container 'staphb/samtools:latest'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}.aligned.bam"), path("${meta.id}.aligned.bam.bai"), emit: out

    script:
    """
    # Rename for consistency with downstream naming, then index
    cp ${bam_file} ${meta.id}.aligned.bam
    samtools index ${meta.id}.aligned.bam
    """
}
