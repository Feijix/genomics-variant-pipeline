process annotateVcf {
    tag "${sample_id}/${caller}"
    container 'ensemblorg/ensembl-vep:latest'
    containerOptions "-v ${params.vep_cache}:/opt/vep/.vep"
    publishDir "${params.outdir}/annotated/${caller}", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), val(caller), path(vcf_file)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.annotated.vcf.gz"),     emit: vcf
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.annotated.vcf.gz.tbi"), emit: tbi
    path("${sample_id}.${caller}.vep_summary.html"),                                        emit: summary

    script:
    """
    vep \\
      --input_file ${vcf_file} \\
      --output_file ${sample_id}.${caller}.annotated.vcf \\
      --format vcf \\
      --vcf \\
      --offline \\
      --fasta ${reference} \\
      --species homo_sapiens \\
      --assembly GRCh38 \\
      --cache \\
      --dir_cache /opt/vep/.vep \\
      --fork ${task.cpus} \\
      --everything \\
      --force_overwrite \\
      --stats_file ${sample_id}.${caller}.vep_summary.html

    bgzip ${sample_id}.${caller}.annotated.vcf
    tabix -p vcf ${sample_id}.${caller}.annotated.vcf.gz
    """
}
