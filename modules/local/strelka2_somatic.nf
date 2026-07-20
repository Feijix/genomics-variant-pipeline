process STRELKA2_SOMATIC {
    tag "${meta.id}"
    container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'
    publishDir "${params.outdir}/somatic/strelka2", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai), path(indel_candidates), path(indel_candidates_tbi)
    path(reference)
    path(fai)
    tuple path(call_regions), path(call_regions_tbi)

    output:
    tuple val(meta), path("${meta.id}.somatic.snvs.vcf.gz"),       emit: snvs
    tuple val(meta), path("${meta.id}.somatic.snvs.vcf.gz.tbi"),   emit: snvs_tbi
    tuple val(meta), path("${meta.id}.somatic.indels.vcf.gz"),     emit: indels
    tuple val(meta), path("${meta.id}.somatic.indels.vcf.gz.tbi"), emit: indels_tbi

    script:
    def regions_arg = call_regions.name != '[]' ? "--callRegions ${call_regions}" : ''
    """
    configureStrelkaSomaticWorkflow.py \\
      --normalBam ${normal_bam} \\
      --tumorBam ${tumor_bam} \\
      --referenceFasta ${reference} \\
      --indelCandidates ${indel_candidates} \\
      ${regions_arg} \\
      --runDir strelka_run

    strelka_run/runWorkflow.py -m local -j ${task.cpus} -g ${task.memory.toGiga()}

    cp strelka_run/results/variants/somatic.snvs.vcf.gz       ${meta.id}.somatic.snvs.vcf.gz
    cp strelka_run/results/variants/somatic.snvs.vcf.gz.tbi   ${meta.id}.somatic.snvs.vcf.gz.tbi
    cp strelka_run/results/variants/somatic.indels.vcf.gz     ${meta.id}.somatic.indels.vcf.gz
    cp strelka_run/results/variants/somatic.indels.vcf.gz.tbi ${meta.id}.somatic.indels.vcf.gz.tbi
    """
}
