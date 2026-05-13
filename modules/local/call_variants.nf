process callVariants {
    tag "${sample_id} (${dv_model})"
    container 'google/deepvariant:1.10.0'
    publishDir "${params.outdir}/vcf/${sample_id}_deepvariant", mode: 'copy'
    cpus 4
    memory '32 GB'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), val(dv_model)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), emit: out

    script:
    """
    mkdir -p output
    /opt/deepvariant/bin/run_deepvariant \\
      --model_type=${dv_model} \\
      --ref=${reference} \\
      --reads=${bam_file} \\
      --output_vcf=output/${sample_id}.deepvariant.vcf.gz \\
      --num_shards=${task.cpus} \\
      --sample_name=${sample_id} \\
      --regions=chr20

    mv output/${sample_id}.deepvariant.vcf.gz ${sample_id}.deepvariant.vcf.gz
    """
}
