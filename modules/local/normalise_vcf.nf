process normaliseVcf {
    tag "${sample_id}/${caller}"
    container 'staphb/bcftools:latest'
    publishDir path: { "${params.outdir}/normalised/${caller}" }, mode: 'copy'

    input:
    tuple val(sample_id), val(caller), path(vcf_file)
    path(reference)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.normalised.vcf.gz"),     emit: vcf
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.normalised.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools norm -f ${reference} -m -any ${vcf_file} -Oz -o ${sample_id}.${caller}.normalised.vcf.gz
    bcftools index -t ${sample_id}.${caller}.normalised.vcf.gz
    """
}
