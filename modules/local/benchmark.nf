process benchmark {
    tag "${sample_id}/${caller}"
    container 'jmcdani20/hap.py:v0.3.12'
    publishDir path: { "${params.outdir}/benchmark/${caller}" }, mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), val(caller), path(vcf_file), path(truth_vcf), path(truth_tbi), path(truth_bed)
    path(reference)
    path(fai)

    output:
    path("${sample_id}.${caller}.benchmark.summary.csv"),     emit: summary
    path("${sample_id}.${caller}.benchmark.metrics.json.gz"), emit: metrics
    path("${sample_id}.${caller}.benchmark.extended.csv"),    emit: extended

    script:
    """
    echo -e "chr20\\t0\\t64444167" > chr20.bed

    /opt/hap.py/bin/hap.py \\
      ${truth_vcf} \\
      ${vcf_file} \\
      -r ${reference} \\
      -f ${truth_bed} \\
      -o ${sample_id}.${caller}.benchmark \\
      -T chr20.bed \\
      --threads ${task.cpus} \\
      --engine=vcfeval
    """
}
