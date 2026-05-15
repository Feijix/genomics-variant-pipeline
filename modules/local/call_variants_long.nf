process callVariantsLong {
    tag "${meta.id} (${meta.clair3_platform}/${meta.clair3_model})"
    container 'hkubal/clair3:latest'
    cpus 4
    memory '8 GB'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.clair3.vcf.gz"), emit: out

    script:
    """
    if [ ! -d "/opt/models/${meta.clair3_model}" ]; then
        echo "ERROR: Clair3 model directory /opt/models/${meta.clair3_model} not found in container" >&2
        echo "Available models:" >&2
        ls /opt/models/ >&2
        exit 1
    fi

    /opt/bin/run_clair3.sh \\
      --bam_fn=${bam_file} \\
      --ref_fn=${reference} \\
      --threads=${task.cpus} \\
      --platform=${meta.clair3_platform} \\
      --model_path=/opt/models/${meta.clair3_model} \\
      --output=\$PWD/clair3_output \\
      --sample_name=${meta.id} \\
      --ctg_name=chr20

    mv clair3_output/merge_output.vcf.gz ${meta.id}.clair3.vcf.gz
    """
}
