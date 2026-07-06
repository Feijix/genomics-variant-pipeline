/*
 * modules/local/sniffles2.nf
 *
 * Long-read structural variant calling with Sniffles2.
 *
 * Input contract matches the other callers: (meta, bam, bai), reference, fai.
 * Emits (meta, vcf). The output is a structural-variant VCF (SVTYPE/SVLEN/END
 * INFO fields), NOT the SNV/indel records the small-variant postprocess chain
 * expects — meta.variant_class = 'sv' (set in the workflow) routes it past the
 * small-variant postprocess.
 *
 * Container: locally-built pinned image (see Dockerfile.sniffles).
 *   docker build -t sniffles:2.8.0 -f Dockerfile.sniffles .
 *
 * --reference is supplied so DEL sequences are emitted. Output filename ends in
 * .gz so Sniffles bgzips + tabix-indexes automatically.
 */

process SNIFFLES2 {
    tag "${meta.id}"
    container 'sniffles:2.8.0'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.sniffles2.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.sniffles2.vcf.gz.tbi"), emit: tbi

    script:
    """
    sniffles \\
      --input ${bam_file} \\
      --reference ${reference} \\
      --threads ${task.cpus} \\
      --vcf ${meta.id}.sniffles2.vcf.gz
    """
}
