process BCFTOOLS_ISEC {
    tag "${meta.id}"
    container 'staphb/bcftools:1.24'

    input:
    tuple val(meta), path(first_vcf), path(second_vcf)

    output:
    tuple val(meta), path("${meta.id}.ensemble.vcf.gz"), emit: vcf

    script:
    """
    bcftools index -t -f ${first_vcf}
    bcftools index -t -f ${second_vcf}
    bcftools isec -n=2 -p isec_out -Oz ${first_vcf} ${second_vcf}
    cp isec_out/0000.vcf.gz ${meta.id}.ensemble.vcf.gz
    bcftools index -t ${meta.id}.ensemble.vcf.gz
    """
}
