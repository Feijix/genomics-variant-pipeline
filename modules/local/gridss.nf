/*
 * modules/local/gridss.nf
 *
 * Short-read structural variant calling with GRIDSS2 (assembly-based, positional
 * de Bruijn). Genuinely different evidence model from Manta/Delly/smoove. Outputs
 * BREAKEND (BND) records, which CONVERT_GRIDSS turns into typed SVs.
 *
 * Two-process chain: GRIDSS (this) -> CONVERT_GRIDSS (BND -> DEL/DUP/INV/INS).
 * meta.variant_class = 'sv'.
 *
 * Heavy and finicky (lessons baked in):
 *   - Needs a CLASSIC bwa index (.amb/.ann/.bwt/.pac/.sa) of the reference, NOT
 *     bwa-mem2's index. Build once: `bwa index <ref>` (provided as params input
 *     here via the reference dir; GRIDSS finds it alongside the FASTA).
 *   - Memory MUST be capped two ways on the shared host: --jvmheap (GRIDSS) and
 *     the Docker --memory cap (set in modules.config containerOptions) so a
 *     runaway can't OOM the host (confirmed multi-day outage previously).
 *   - A stale GRCh38.fa.tmp.gridsslock from a crashed run blocks setupreference;
 *     the script removes it defensively before running.
 *
 * Container: BioContainers conda build (bundles gridss, bwa, samtools, R +
 * StructuralVariantAnnotation).
 *   quay.io/biocontainers/gridss:2.13.2--h96c455f_6
 */

process GRIDSS {
    tag "${meta.id}"
    container 'quay.io/biocontainers/gridss:2.13.2--h96c455f_6'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)       // must have a classic bwa index alongside it
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.gridss.bnd.vcf.gz"), emit: bnd

    script:
    """
    # Clear any stale lock from a previously-killed run (defensive).
    rm -rf ${reference}.tmp.gridsslock || true

    gridss \\
      --reference ${reference} \\
      --output ${meta.id}.gridss.bnd.vcf.gz \\
      --threads ${task.cpus} \\
      --jvmheap ${task.memory.toGiga() - 4}g \\
      --workingdir gridss_work \\
      --picardoptions VALIDATION_STRINGENCY=LENIENT \\
      ${bam_file}
    """
}
