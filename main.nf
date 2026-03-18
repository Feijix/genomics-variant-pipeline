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
    .map { row ->
        def fastq_2 = row.fastq_2 ? file(row.fastq_2) : []
        tuple(row.sample, row.sequencing_type, [file(row.fastq_1), fastq_2])
    }
    .branch {
        shortread: it[1] == 'short'
        longread:  it[1] == 'long'
    }
    .set { branched }

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {
    indexReference(file(params.reference))
    faiReference(file(params.reference))

    // ── Short-read path
    align(
        branched.shortread.map { sample, type, reads -> tuple(sample, reads) },
        file(params.reference),
        indexReference.out.index
    )
    sortBam(align.out)
    callVariants(sortBam.out, file(params.reference), faiReference.out.fai)

    // ── Long-read path
    alignLongRead(
        branched.longread.map { sample, type, reads -> tuple(sample, reads[0]) },
        file(params.reference)
    )
    sortBamLong(alignLongRead.out)
    callVariantsLong(sortBamLong.out, file(params.reference), faiReference.out.fai)

    // ── Converge both paths at QC
    vcfQc(
        callVariants.out.mix(callVariantsLong.out)
    )
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

// ── Short-read Alignment ──────────────────────────────────────────────────────

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

// ── Long-read Alignment ───────────────────────────────────────────────────────

process alignLongRead {
    tag "${sample_id}"
    container 'staphb/minimap2:latest'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), emit: out

    script:
    // map-ont preset applies the Oxford Nanopore error model
    // -a outputs SAM format (minimap2 defaults to PAF)
    // @RG tag same convention as short-read for consistency downstream
    """
    minimap2 -ax map-ont -t ${task.cpus} \
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ONT" \
      ${reference} ${reads} > ${sample_id}.sam
    """
}

// ── BAM Processing (shared by both paths) ────────────────────────────────────

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

process sortBamLong {
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

// ── Short-read Variant Calling ────────────────────────────────────────────────

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

// ── Long-read Variant Calling ─────────────────────────────────────────────────

process callVariantsLong {
    tag "${sample_id}"
    container 'hkubal/clair3:latest'
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
    // Clair3 requires absolute paths — $PWD resolves this inside the work dir
    // ont is the model for Oxford Nanopore reads
    // --include_all_ctgs calls variants on all contigs, not just chr1-22/X/Y
    """
    /opt/bin/run_clair3.sh \
      --bam_fn=${bam_file} \
      --ref_fn=${reference} \
      --threads=${task.cpus} \
      --platform=ont \
      --model_path=/opt/models/ont \
      --output=\$PWD/clair3_output \
      --sample_name=${sample_id} \
      --include_all_ctgs

    mv clair3_output/merge_output.vcf.gz ${sample_id}.vcf.gz
    """
}

// ── VCF QC & Filtering (shared by both paths) ─────────────────────────────────

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