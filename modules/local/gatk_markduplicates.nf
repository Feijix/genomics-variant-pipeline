process GATK_MARKDUPLICATES {
    tag "${meta.id}"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/bam", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bam.bai"), emit: bam
    path("${meta.id}.md.metrics.txt"),                                          emit: metrics

    script:
    def avail_mem = task.memory ? task.memory.toGiga() - 2 : 6
    """
    gatk --java-options "-Xmx${avail_mem}g" MarkDuplicates \\
      -I ${bam} \\
      -O ${meta.id}.md.bam \\
      -M ${meta.id}.md.metrics.txt \\
      --ASSUME_SORT_ORDER coordinate \\
      --VALIDATION_STRINGENCY LENIENT \\
      --CREATE_INDEX true \\
      --TMP_DIR .

    # Picard writes <prefix>.bai; rename to the <prefix>.bam.bai convention used elsewhere.
    mv ${meta.id}.md.bai ${meta.id}.md.bam.bai
    """
}
