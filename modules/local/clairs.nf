// modules/local/clairs.nf
//
// Long-read somatic small-variant calling, tumour-normal paired.
//
// Model choice (ont_r10_dorado_sup_5khz_ssrs) is deliberate:
//   - `ssrs` = synthetic + real-cancer augmented. Its real-sample training set
//     is HCC1937/HCC1954/H1437/H2009 — HCC1395 is NOT included, so there is no
//     circularity when benchmarking on HCC1395.
//   - All ClairS models exclude chr20 from training (chr1-19, 21, 22 only),
//     so chr20 is a held-out test chromosome.
//   - The model's validated configuration is GRCh38_no_alt + minimap2, which
//     matches the CASTLE alignments used here.
//
// --enable_indel_calling is opt-in: without it ClairS emits SNVs only and the
// indel benchmark would silently report zero recall. ONT R10 only.
process CLAIRS {
    tag "${meta.id}"
    container 'hkubal/clairs:v0.4.5'
    cpus 8
    memory '32 GB'

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    path(reference)
    path(fai)

    output:
    tuple val(meta), path("${meta.id}.clairs.snv.vcf.gz"),       emit: snvs
    tuple val(meta), path("${meta.id}.clairs.snv.vcf.gz.tbi"),   emit: snvs_tbi
    tuple val(meta), path("${meta.id}.clairs.indel.vcf.gz"),     emit: indels
    tuple val(meta), path("${meta.id}.clairs.indel.vcf.gz.tbi"), emit: indels_tbi

    script:
    def ctg_arg = params.target_region ? "-c ${params.target_region}" : ''
    def platform = meta.clairs_platform ?: params.clairs_platform
    """
    run_clairs \\
      -T ${tumor_bam} \\
      -N ${normal_bam} \\
      -R ${reference} \\
      -o clairs_output \\
      -t ${task.cpus} \\
      -p ${platform} \\
      -s ${meta.id} \\
      ${ctg_arg} \\
      --enable_indel_calling

    cp clairs_output/output.vcf.gz     ${meta.id}.clairs.snv.vcf.gz
    cp clairs_output/output.vcf.gz.tbi ${meta.id}.clairs.snv.vcf.gz.tbi
    cp clairs_output/indel.vcf.gz      ${meta.id}.clairs.indel.vcf.gz
    cp clairs_output/indel.vcf.gz.tbi  ${meta.id}.clairs.indel.vcf.gz.tbi
    """
}
