process BCFTOOLS_CONCAT_SOMATIC {
    tag "${meta.id}"
    container 'staphb/bcftools:1.24'
    publishDir "${params.outdir}/somatic", mode: 'copy'

    input:
    tuple val(meta), path(snvs), path(indels)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.${meta.caller}.somatic.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.${meta.caller}.somatic.vcf.gz.tbi"), emit: tbi

    script:
    """
    # Normalise both, concat (SNVs + indels are disjoint variant types),
    # sort. FILTER column preserved — do NOT drop non-PASS here: som.py
    # needs the full FILTER-annotated callset to compute precision/recall.
    bcftools norm -f ${reference} -m -any ${snvs}   -Oz -o snvs.norm.vcf.gz
    bcftools norm -f ${reference} -m -any ${indels} -Oz -o indels.norm.vcf.gz
    bcftools index -t snvs.norm.vcf.gz
    bcftools index -t indels.norm.vcf.gz

    bcftools concat -a snvs.norm.vcf.gz indels.norm.vcf.gz -Oz -o concat.vcf.gz
    bcftools sort concat.vcf.gz -Oz -o ${meta.id}.${meta.caller}.somatic.vcf.gz
    bcftools index -t ${meta.id}.${meta.caller}.somatic.vcf.gz
    """
}