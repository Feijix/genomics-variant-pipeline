/*
 * modules/local/delly.nf
 *
 * Short-read structural variant calling with Delly (second short-read SV
 * caller, paired with Manta for SV-SV ensembles).
 *
 * Input contract matches the other callers: (meta, bam, bai), reference, fai.
 * Emits (meta, vcf). meta.variant_class = 'sv' routes it past the
 * small-variant postprocess, same as the other SV callers.
 *
 * Delly is a single-command caller (simpler than Manta's two-step) but writes
 * BCF natively. The BioContainers Delly image ships ONLY delly (no bcftools /
 * bgzip / tabix), so the BCF -> bgzipped-VCF conversion is done in a following
 * bcftools step in the workflow, not here. This module emits the BCF + its csi
 * index; CONVERT_DELLY (bcftools container) produces the final .vcf.gz.
 *
 * Requirements (handled upstream in the chr20 BAM prep):
 *   - sorted + indexed BAM (done)
 *   - duplicate-marked BAM (samtools markdup) — Delly uses read-depth signal,
 *     so unmarked duplicates inflate calls
 *   - chr20-only reference (GRCh38.chr20.fa), matching the chr20-only BAM, to
 *     satisfy Delly's split-read reference lookups (same as Manta)
 *
 * Container: BioContainers conda build.
 *   quay.io/biocontainers/delly:1.7.3--hd6466ae_0
 */

process DELLY {
    tag "${meta.id}"
    container 'quay.io/biocontainers/delly:1.7.3--hd6466ae_0'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.delly.bcf"),     emit: bcf
    tuple val(meta), path("${meta.id}.delly.bcf.csi"), emit: csi

    script:
    """
    delly call \\
      -g ${reference} \\
      -o ${meta.id}.delly.bcf \\
      ${bam_file}
    """
}
