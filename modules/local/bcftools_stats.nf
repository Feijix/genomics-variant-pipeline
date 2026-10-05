// Variant-level QC: counts, Ts/Tv, indel length distribution.
process BCFTOOLS_STATS {
    tag "${meta.id}.${meta.caller}"
    container 'staphb/bcftools:1.24'
    publishDir "${params.outdir}/qc/bcftools", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.bcftools_stats.txt"), emit: stats

    script:
    """
    bcftools stats ${vcf} > ${meta.id}.${meta.caller}.bcftools_stats.txt
    """
}
