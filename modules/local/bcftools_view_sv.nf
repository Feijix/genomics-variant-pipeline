// modules/local/bcftools_view_sv.nf

process BCFTOOLS_VIEW_SV {
    tag "${meta.id}.${meta.caller}"
    container 'staphb/bcftools:latest'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.pass.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.${meta.caller}.pass.vcf.gz.tbi"), emit: tbi
    
    script:
    """
    bcftools view -f 'PASS,.' ${vcf} -Oz -o ${meta.id}.${meta.caller}.pass.vcf.gz
    bcftools index -t ${meta.id}.${meta.caller}.pass.vcf.gz
    """
}
