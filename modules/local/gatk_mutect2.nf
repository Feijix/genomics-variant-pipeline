process GATK_MUTECT2 {
    tag "${meta.id}"
    container 'broadinstitute/gatk:4.6.2.0'
    publishDir "${params.outdir}/somatic/mutect2", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(reference)
    path(fai)
    path(dict)
    tuple path(germline_resource), path(germline_resource_tbi)
    tuple path(pon), path(pon_tbi)

    output:
    tuple val(meta), path("${meta.id}.mutect2.unfiltered.vcf.gz"),       emit: vcf
    tuple val(meta), path("${meta.id}.mutect2.unfiltered.vcf.gz.tbi"),   emit: tbi
    tuple val(meta), path("${meta.id}.mutect2.unfiltered.vcf.gz.stats"), emit: stats

    script:
    def normal_arg = meta.analysis == 'somatic_paired' ? "-I ${normal_bam} -normal ${meta.normal_id}" : ''   
    def interval_arg = params.target_region    ? "-L ${params.target_region}"                         : ''
    def avail_mem    = task.memory ? task.memory.toGiga() - 2 : 6
    """
    gatk --java-options "-Xmx${avail_mem}g" Mutect2 \\
      -R ${reference} \\
      -I ${tumor_bam} \\
      ${normal_arg} \\
      --germline-resource ${germline_resource} \\
      --panel-of-normals ${pon} \\
      ${interval_arg} \\
      --native-pair-hmm-threads ${task.cpus} \\
      -O ${meta.id}.mutect2.unfiltered.vcf.gz
    """
}
