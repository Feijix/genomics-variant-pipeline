process HAPPY_HAPPY {
    tag "${meta.id}/${meta.caller}"
    container 'jmcdani20/hap.py:v0.3.12'
    cpus 4
    memory '16 GB'

    input:
    tuple val(meta), path(vcf_file), path(truth_vcf), path(truth_tbi), path(truth_bed)
    path(reference)
    path(fai)

    output:
    path("${meta.id}.${meta.caller}.benchmark.summary.csv"),     emit: summary
    path("${meta.id}.${meta.caller}.benchmark.metrics.json.gz"), emit: metrics
    path("${meta.id}.${meta.caller}.benchmark.extended.csv"),    emit: extended

    script:
    """
    echo -e "chr20\\t0\\t64444167" > chr20.bed

    /opt/hap.py/bin/hap.py \\
      ${truth_vcf} \\
      ${vcf_file} \\
      -r ${reference} \\
      -f ${truth_bed} \\
      -o ${meta.id}.${meta.caller}.benchmark \\
      -T chr20.bed \\
      --threads ${task.cpus} \\
      --engine=vcfeval
    """
}
