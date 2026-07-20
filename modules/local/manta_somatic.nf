process MANTA_SOMATIC {
    tag "${meta.id}"
    container 'quay.io/biocontainers/manta:1.6.0--py27h9948957_6'
    publishDir "${params.outdir}/somatic/manta", mode: 'copy'
    cpus 8
    memory '16 GB'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(reference)
    path(fai)
    tuple path(call_regions), path(call_regions_tbi)

    output:
    tuple val(meta), path("${meta.id}.somaticSV.vcf.gz"),                emit: sv
    tuple val(meta), path("${meta.id}.somaticSV.vcf.gz.tbi"),            emit: sv_tbi
    tuple val(meta), path("${meta.id}.candidateSmallIndels.vcf.gz"),     emit: indel_candidates
    tuple val(meta), path("${meta.id}.candidateSmallIndels.vcf.gz.tbi"), emit: indel_candidates_tbi

    script:
    def regions_arg = call_regions.name != '[]' ? "--callRegions ${call_regions}" : ''
    """
    # BAMs are chr20-restricted (single-contig header), so no --callRegions needed.
    configManta.py \\
      --normalBam ${normal_bam} \\
      --tumorBam ${tumor_bam} \\
      --referenceFasta ${reference} \\
      ${regions_arg} \\
      --runDir manta_run

    manta_run/runWorkflow.py -m local -j ${task.cpus} -g ${task.memory.toGiga()}

    cp manta_run/results/variants/somaticSV.vcf.gz            ${meta.id}.somaticSV.vcf.gz
    cp manta_run/results/variants/somaticSV.vcf.gz.tbi        ${meta.id}.somaticSV.vcf.gz.tbi
    cp manta_run/results/variants/candidateSmallIndels.vcf.gz     ${meta.id}.candidateSmallIndels.vcf.gz
    cp manta_run/results/variants/candidateSmallIndels.vcf.gz.tbi ${meta.id}.candidateSmallIndels.vcf.gz.tbi
    """
}
