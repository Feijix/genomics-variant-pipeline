process vcfQc {
    tag "${sample_id}/${caller}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/qc/${caller}", mode: 'copy'

    input:
    tuple val(sample_id), val(caller), path(vcf_file)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.final.vcf.gz"), emit: vcf
    path("${sample_id}.${caller}.stats.txt"), emit: stats

    script:
    """
    bcftools view -f PASS ${vcf_file} -Oz -o ${sample_id}.${caller}.final.vcf.gz
    bcftools stats ${sample_id}.${caller}.final.vcf.gz > ${sample_id}.${caller}.stats.txt
    """
}
