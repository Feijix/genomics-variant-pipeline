// Alignment-level QC. MultiQC parses samtools stats natively.
process SAMTOOLS_STATS {
    tag "${meta.id}"
    container 'staphb/samtools:1.24'
    publishDir "${params.outdir}/qc/samtools", mode: 'copy'
    cpus 2

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.samtools_stats.txt"), emit: stats

    script:
    """
    samtools stats -@ ${task.cpus} ${bam} > ${meta.id}.samtools_stats.txt
    """
}
