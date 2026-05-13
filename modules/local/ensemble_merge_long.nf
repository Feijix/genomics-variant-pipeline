process ensembleMergeLong {
    tag "${sample_id}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/vcf/${sample_id}_ensemble", mode: 'copy'

    input:
    tuple val(sample_id), path(clair3_vcf), path(dv_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.ensemble.vcf.gz"), emit: vcf

    script:
    """
    bcftools index -t -f ${clair3_vcf}
    bcftools index -t -f ${dv_vcf}
    bcftools isec -n=2 -p isec_out -Oz ${clair3_vcf} ${dv_vcf}
    cp isec_out/0000.vcf.gz ${sample_id}.ensemble.vcf.gz
    bcftools index -t ${sample_id}.ensemble.vcf.gz
    """
}
