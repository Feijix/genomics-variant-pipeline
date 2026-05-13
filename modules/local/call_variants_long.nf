process callVariantsLong {
    tag "${sample_id} (${clair3_platform}/${clair3_model})"
    container 'hkubal/clair3:latest'
    publishDir "${params.outdir}/vcf/${sample_id}_clair3", mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), val(clair3_platform), val(clair3_model)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.clair3.vcf.gz"), emit: out

    script:
    """
    if [ ! -d "/opt/models/${clair3_model}" ]; then
        echo "ERROR: Clair3 model directory /opt/models/${clair3_model} not found in container" >&2
        echo "Available models:" >&2
        ls /opt/models/ >&2
        exit 1
    fi

    /opt/bin/run_clair3.sh \\
      --bam_fn=${bam_file} \\
      --ref_fn=${reference} \\
      --threads=${task.cpus} \\
      --platform=${clair3_platform} \\
      --model_path=/opt/models/${clair3_model} \\
      --output=\$PWD/clair3_output \\
      --sample_name=${sample_id} \\
      --ctg_name=chr20

    mv clair3_output/merge_output.vcf.gz ${sample_id}.clair3.vcf.gz
    """
}
