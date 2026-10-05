process BCFTOOLS_VIEW {
    tag "${meta.id}/${meta.caller}"
    container 'staphb/bcftools:1.24'

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.final.vcf.gz"), emit: vcf
    path("${meta.id}.${meta.caller}.stats.txt"),                     emit: stats

    script:
    """
    bcftools view -f PASS ${vcf_file} -Oz -o ${meta.id}.${meta.caller}.final.vcf.gz
    bcftools stats ${meta.id}.${meta.caller}.final.vcf.gz > ${meta.id}.${meta.caller}.stats.txt
    """
}
