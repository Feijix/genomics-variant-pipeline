process GATK_FILTERMUTECTCALLS {
    tag "${meta.id}"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/somatic/mutect2", mode: 'copy'
    cpus 2
    memory '8 GB'

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats)
    path(reference)
    path(fai)
    path(dict)

    output:
    tuple val(meta), path("${meta.id}.mutect2.filtered.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.mutect2.filtered.vcf.gz.tbi"), emit: tbi
    path("${meta.id}.filteringStats.tsv"),                            emit: filter_stats

    script:
    def avail_mem = task.memory ? task.memory.toGiga() - 2 : 6
    """
    gatk --java-options "-Xmx${avail_mem}g" FilterMutectCalls \\
      -R ${reference} \\
      -V ${vcf} \\
      --stats ${stats} \\
      --filtering-stats ${meta.id}.filteringStats.tsv \\
      -O ${meta.id}.mutect2.filtered.vcf.gz
    """
}
