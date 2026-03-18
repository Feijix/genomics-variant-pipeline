#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ── Parameters ────────────────────────────────────────────────────────────────
params.samplesheet = 'data/samplesheet.csv'
params.reference   = 'data/reference/chr20.fa'
params.outdir      = 'data/outputs'

// ── Load Samples ──────────────────────────────────────────────────────────────
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample, [file(row.fastq_1), file(row.fastq_2)]) }
    .set { samples_ch }

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {
    indexReference(file(params.reference))
    faiReference(file(params.reference))

    align(samples_ch, file(params.reference), indexReference.out.index)
    sortBam(align.out)
    callVariants(sortBam.out, file(params.reference), faiReference.out.fai)
    vcfQc(callVariants.out)
}

// ── Reference Indexing ────────────────────────────────────────────────────────

process indexReference {
    tag "bwa-mem2 index"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    storeDir 'data/reference'

    input:
    path(reference)

    output:
    path("${reference}.*"), emit: index

    script:
    """
    bwa-mem2 index ${reference}
    """
}

process faiReference {
    tag "samtools faidx"
    container 'staphb/samtools:latest'
    storeDir 'data/reference'

    input:
    path(reference)

    output:
    path("${reference}.fai"), emit: fai

    script:
    """
    samtools faidx ${reference}
    """
}

// ── Alignment ─────────────────────────────────────────────────────────────────

process align {
    tag "${sample_id}"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)
    path(index)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: out

    script:
    """
    bwa-mem2 mem -t ${task.cpus} \
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
      ${reference} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}

// ── BAM Processing ────────────────────────────────────────────────────────────

process sortBam {
    tag "${sample_id}"
    container 'staphb/samtools:latest'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(sam_file)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: out

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam ${sam_file}
    samtools index ${sample_id}.sorted.bam
    """
}

// ── Variant Calling ───────────────────────────────────────────────────────────

process callVariants {
    tag "${sample_id}"
    container 'google/deepvariant:1.5.0'
    publishDir "${params.outdir}/vcf", mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: out

    script:
    """
    mkdir -p output

    /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --ref=${reference} \
      --reads=${bam_file} \
      --output_vcf=output/${sample_id}.vcf.gz \
      --num_shards=${task.cpus} \
      --sample_name=${sample_id}

    mv output/${sample_id}.vcf.gz ${sample_id}.vcf.gz
    """
}

// ── VCF QC & Filtering ────────────────────────────────────────────────────────

process vcfQc {
    tag "${sample_id}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file)

    output:
    path("${sample_id}.final.vcf.gz"), emit: vcf
    path("${sample_id}.stats.txt"),    emit: stats

    script:
    """
    bcftools view -i 'QUAL>20' ${vcf_file} -Oz -o ${sample_id}.final.vcf.gz
    bcftools stats ${sample_id}.final.vcf.gz > ${sample_id}.stats.txt
    """
}