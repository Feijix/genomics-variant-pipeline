/*
 * modules/local/manta.nf
 *
 * Short-read structural variant calling with Manta (Illumina).
 * First short-read SV caller, paired with Delly for SV-SV ensembles.
 *
 * Input contract matches the other callers: (meta, bam, bai), reference, fai.
 * Emits (meta, vcf). Output is a structural-variant VCF; meta.variant_class =
 * 'sv' routes it past the small-variant postprocess, same as Sniffles2/cuteSV.
 *
 * Manta is a TWO-STEP workflow (unlike the single-command long-read callers):
 *   1. configManta.py  — generates a run directory + execution scripts
 *   2. runWorkflow.py  — executes the assembly/scoring locally
 * The germline diploid SV output is results/variants/diploidSV.vcf.gz.
 *
 * The input BAM is pre-restricted to chr20 (single-contig header from the
 * subsetting step), so Manta calls only on chr20 without --callRegions. This
 * also avoids needing bgzip/tabix, which the Manta image does not provide.
 *
 * Container: BioContainers conda build (Python 2 handled internally).
 *   quay.io/biocontainers/manta:1.6.0-6
 */

process MANTA {
    tag "${meta.id}"
    container 'quay.io/biocontainers/manta:1.6.0--py27h9948957_6'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.manta.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.manta.vcf.gz.tbi"), emit: tbi

    script:
    """
    # The input BAM is already restricted to chr20 (single-contig header), so
    # Manta naturally calls only on chr20 — no --callRegions needed (and the
    # Manta image ships no bgzip/tabix to build one). configManta.py then
    # runWorkflow.py is Manta's required two-step pattern.

    # Step 1: configure the run (generates runWorkflow.py under manta_run/).
    configManta.py \\
      --bam ${bam_file} \\
      --referenceFasta ${reference} \\
      --runDir manta_run

    # Step 2: execute locally with the allotted cpus/memory.
    manta_run/runWorkflow.py \\
      -m local \\
      -j ${task.cpus} \\
      -g ${task.memory.toGiga()}

    # Germline diploid SV calls are the primary output.
    cp manta_run/results/variants/diploidSV.vcf.gz     ${meta.id}.manta.vcf.gz
    cp manta_run/results/variants/diploidSV.vcf.gz.tbi ${meta.id}.manta.vcf.gz.tbi
    """
}
