/*
 * modules/local/convert_delly.nf
 *
 * Converts Delly's native BCF output to a bgzipped, tabix-indexed VCF, so it
 * matches the (meta, vcf.gz) contract of the other SV callers and can be folded
 * into sv_vcfs for ensembling/benchmarking. Separate from the DELLY process
 * because the BioContainers Delly image ships only delly (no bcftools).
 *
 * Container: the same bcftools image used elsewhere in the pipeline.
 */

process CONVERT_DELLY {
    tag "${meta.id}"
    container 'staphb/bcftools:latest'

    input:
    tuple val(meta), path(bcf_file), path(csi_file)

    output:
    tuple val(meta), path("${meta.id}.delly.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.delly.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools view ${bcf_file} -Oz -o ${meta.id}.delly.vcf.gz
    bcftools index -t ${meta.id}.delly.vcf.gz
    """
}
