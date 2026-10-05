process SAMTOOLS_INDEX {
    tag "${meta.id}"
    container 'staphb/samtools:1.24'

    input:
    tuple val(meta), path(bam_file)

    output:
    tuple val(meta), path("${meta.id}.aligned.bam"), path("${meta.id}.aligned.bam.bai"), emit: out

    script:
    """
    # Rename for consistency with downstream naming, then index
    if [ "${bam_file}" != "${meta.id}.aligned.bam" ]; then
        ln -s ${bam_file} ${meta.id}.aligned.bam
    fi
    samtools index ${meta.id}.aligned.bam
    """
}
