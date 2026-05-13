process sortBam {
    tag "${sample_id}"
    container 'staphb/samtools:latest'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: out

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam ${sam_file}
    samtools index ${sample_id}.sorted.bam
    """
}
