process BCFTOOLS_NORM {
    tag "${meta.id}/${meta.caller}"
    container 'staphb/bcftools:latest'

    input:
    tuple val(meta), path(vcf_file)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.normalised.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.${meta.caller}.normalised.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools norm -f ${reference} -m -any ${vcf_file} -Oz -o ${meta.id}.${meta.caller}.normalised.vcf.gz
    bcftools index -t ${meta.id}.${meta.caller}.normalised.vcf.gz
    """
}
