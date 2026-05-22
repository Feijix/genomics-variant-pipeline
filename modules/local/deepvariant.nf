process DEEPVARIANT {
    tag "${meta.id} (${meta.model_type})"
    container 'google/deepvariant:1.10.0'
    cpus 4
    memory '32 GB'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.deepvariant.vcf.gz"), emit: out

    script:
    """
    mkdir -p output
    /opt/deepvariant/bin/run_deepvariant \\
      --model_type=${meta.model_type} \\
      --ref=${reference} \\
      --reads=${bam_file} \\
      --output_vcf=output/${meta.id}.deepvariant.vcf.gz \\
      --num_shards=${task.cpus} \\
      --sample_name=${meta.id} \\
      --regions=chr20

    mv output/${meta.id}.deepvariant.vcf.gz ${meta.id}.deepvariant.vcf.gz
    """
}
