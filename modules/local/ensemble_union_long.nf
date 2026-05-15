process ensembleUnionLong {
    tag "${meta.id}"
    container 'staphb/bcftools:latest'

    input:
    tuple val(meta), path(clair3_vcf), path(dv_vcf)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.ensemble_union.vcf.gz"), emit: vcf

    script:
    """
    bcftools norm -f ${reference} -m -any ${clair3_vcf} -Oz -o clair3.norm.vcf.gz
    bcftools norm -f ${reference} -m -any ${dv_vcf}     -Oz -o dv.norm.vcf.gz
    bcftools index -t clair3.norm.vcf.gz
    bcftools index -t dv.norm.vcf.gz

    bcftools concat -a clair3.norm.vcf.gz dv.norm.vcf.gz -Oz -o concat.vcf.gz
    bcftools sort concat.vcf.gz -Oz -o sorted.vcf.gz
    bcftools norm -d exact sorted.vcf.gz -Oz -o ${meta.id}.ensemble_union.vcf.gz
    bcftools index -t ${meta.id}.ensemble_union.vcf.gz
    """
}
