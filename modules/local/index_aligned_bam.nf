process indexAlignedBam {
    tag "${sample_id}"
    container 'staphb/samtools:latest'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: out

    script:
    """
    # Rename for consistency with downstream naming, then index
    cp ${bam_file} ${sample_id}.aligned.bam
    samtools index ${sample_id}.aligned.bam
    """
}
