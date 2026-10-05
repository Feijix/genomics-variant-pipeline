/*
 * modules/local/smoove.nf
 *
 * Short-read structural variant calling with smoove (brentp/smoove), a wrapper
 * around lumpy + lumpy_filter + svtyper that adds read-filtering and genotyping.
 * Fifth SV caller; lumpy-based PE+split-read evidence model.
 *
 * Input contract matches the other SV callers: (meta, bam, bai), reference, fai.
 * Emits (meta, vcf). meta.variant_class = 'sv' routes it past small-variant
 * postprocess, same as the other SV callers.
 *
 * smoove is a single-command caller but has strict input requirements that the
 * chr20 dev BAM must satisfy (all learned the hard way):
 *   - coordinate-sorted with SO:coordinate in @HD (svtyper checks this)
 *   - an @RG line with SM set (svtyper reads bam.header['RG'][0]['SM'])
 *   - SA tags intact (lumpy_filter extracts split reads by SA tag) — so the BAM
 *     must come from bwa-mem2, NOT an SA-stripped BAM
 *   - single-contig (chr20-only) header + matching chr20-only reference, else
 *     smoove trips on other contigs (e.g. chrM) during region iteration
 * The upstream alignment/subset steps must produce such a BAM; this module
 * assumes it (the same chr20-only bwa-mem2 BAM used for the other SV callers).
 *
 * --exclude is intentionally omitted: smoove's BED parser rejects Delly-style
 * exclude files with whole-contig (no-coordinate) lines. A proper 3-column BED
 * could be supplied via params if needed.
 *
 * Container: brentp/smoove (Docker Hub; bundles lumpy, lumpy_filter, svtyper,
 * gsort, bgzip, tabix, samtools). Memory-capped per the shared-host OOM history.
 */

process SMOOVE {
    tag "${meta.id}"
    container 'brentp/smoove:v0.2.8'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.smoove.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.smoove.vcf.gz.tbi"), emit: tbi

    script:
    """
    mkdir -p smoove_out

    smoove call \\
      --name ${meta.id} \\
      --outdir smoove_out \\
      --fasta ${reference} \\
      -p ${task.cpus} \\
      --genotype \\
      ${bam_file}

    # smoove writes <name>-smoove.genotyped.vcf.gz; rename to the pipeline scheme
    # and (re)create a .tbi ourselves, since smoove may emit .csi instead.
    cp smoove_out/${meta.id}-smoove.genotyped.vcf.gz ${meta.id}.smoove.vcf.gz
    tabix -p vcf ${meta.id}.smoove.vcf.gz
    """
}
