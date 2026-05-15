process annotateVcf {
    tag "${meta.id}/${meta.caller}"
    container 'ensemblorg/ensembl-vep:latest'
    containerOptions "-v ${params.vep_cache}:/opt/vep/.vep"
    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(vcf_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.annotated.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.${meta.caller}.annotated.vcf.gz.tbi"), emit: tbi
    path("${meta.id}.${meta.caller}.vep_summary.html"),                      emit: summary

    script:
    """
    vep \\
      --input_file ${vcf_file} \\
      --output_file ${meta.id}.${meta.caller}.annotated.vcf \\
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
      --stats_file ${meta.id}.${meta.caller}.vep_summary.html

    bgzip ${meta.id}.${meta.caller}.annotated.vcf
    tabix -p vcf ${meta.id}.${meta.caller}.annotated.vcf.gz
    """
}
