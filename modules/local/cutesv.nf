/*
 * modules/local/cutesv.nf
 *
 * Long-read structural variant calling with cuteSV (second SV caller,
 * alongside Sniffles2, to enable SV-SV ensembles).
 *
 * Input contract matches the other callers: (meta, bam, bai), reference, fai.
 * Emits (meta, vcf). Like Sniffles2, the output is a structural-variant VCF;
 * meta.variant_class = 'sv' routes it past the small-variant postprocess.
 *
 * Platform-specific clustering parameters are REQUIRED for correct calling —
 * cuteSV's INS/DEL cluster-bias and merging-ratio thresholds differ by
 * technology (HiFi vs ONT vs CLR). These are supplied from the routing config
 * via meta.cutesv_params (the cuteSV analogue of Clair3's model selection).
 * Wrong parameters silently degrade accuracy.
 *
 * cuteSV requires a scratch work directory as its 4th positional argument;
 * it is created in the isolated task dir (C1).
 *
 * Container: locally-built pinned image (see Dockerfile.cutesv).
 *   docker build -t cutesv:2.1.3 -f Dockerfile.cutesv .
 */

process CUTESV {
    tag "${meta.id}"
    container 'cutesv:2.1.3'

    input:
    tuple val(meta), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.cutesv.vcf.gz"),     emit: vcf
    tuple val(meta), path("${meta.id}.cutesv.vcf.gz.tbi"), emit: tbi

    script:
    """
    mkdir -p cutesv_work

    # cuteSV iterates over every contig in the BAM header AND follows split-read
    # pointers in each read's SA (supplementary alignment) tag. GIAB/PacBio
    # pbmm2 BAMs carry extra gene-patch decoy contigs (KMT2C/MAP2K3/KCNJ18)
    # absent from plain GRCh38; chr20 reads have SA tags pointing to those decoy
    # contigs, so even a chr20-only BAM crashes cuteSV when it resolves the SA
    # targets against the reference.
    #
    # Method (verified): extract chr20 reads as headerless SAM (referencing
    # 'chr20' by NAME), strip the SA tag so cuteSV cannot chase split alignments
    # onto absent decoy contigs, then re-encode to BAM against a chr20-only
    # dictionary so the header lists only chr20 and the BAM indexes cleanly.
    # SV calling proceeds from primary chr20 alignments; the discarded SA targets
    # were on decoy contigs outside the chr20 development scope anyway.
    samtools faidx ${reference} chr20 > chr20_ref.fa
    samtools faidx chr20_ref.fa                       # -> chr20_ref.fa.fai (1 contig)

    samtools view --remove-tag SA ${bam_file} chr20 \\
      | samtools view -b -t chr20_ref.fa.fai - \\
      > ${meta.id}.chr20.bam
    samtools index ${meta.id}.chr20.bam

    cuteSV \\
      ${meta.id}.chr20.bam \\
      ${reference} \\
      ${meta.id}.cutesv.vcf \\
      cutesv_work \\
      --threads ${task.cpus} \\
      --sample ${meta.id} \\
      ${meta.cutesv_params}

    bgzip ${meta.id}.cutesv.vcf
    tabix -p vcf ${meta.id}.cutesv.vcf.gz
    """
}
