process BCFTOOLS_ISEC {
    tag "${meta.id}"
    container 'staphb/bcftools:1.24'

    input:
    tuple val(meta), path(clair3_vcf), path(dv_vcf)

    output:
    tuple val(meta), path("${meta.id}.ensemble.vcf.gz"), emit: vcf

    script:
    """
    bcftools index -t -f ${clair3_vcf}
    bcftools index -t -f ${dv_vcf}
    bcftools isec -n=2 -p isec_out -Oz ${clair3_vcf} ${dv_vcf}
    cp isec_out/0000.vcf.gz ${meta.id}.ensemble.vcf.gz
    bcftools index -t ${meta.id}.ensemble.vcf.gz
    """
}
