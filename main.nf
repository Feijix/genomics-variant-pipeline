#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ── Parameters ────────────────────────────────────────────────────────────────
params.samplesheet = 'data/samplesheet.csv'
params.reference   = 'data/reference/GRCh38.fa'
params.outdir      = 'data/outputs'
params.vep_cache   = "$HOME/vep_data"

// ── Load Samples ──────────────────────────────────────────────────────────────
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row ->
        def fastq_2  = row.fastq_2  ? file(row.fastq_2)  : []
        def truth_vcf = row.truth_vcf ? file(row.truth_vcf) : []
        def truth_bed = row.truth_bed ? file(row.truth_bed) : []
        tuple(row.sample, row.sequencing_type, [file(row.fastq_1), fastq_2], truth_vcf, truth_bed)
    }
    .branch {
        shortread: it[1] == 'shortread'
        longread:  it[1] == 'longread'
    }
    .set { branched }

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {
    indexReference(file(params.reference))
    faiReference(file(params.reference))

    // ── Short-read path
    align(
        branched.shortread.map { sample, type, reads, truth_vcf, truth_bed -> tuple(sample, reads) },
        file(params.reference),
        indexReference.out.index
    )
    sortBam(align.out)
    callVariants(sortBam.out, file(params.reference), faiReference.out.fai)

    // ── Long-read path
    alignLongRead(
        branched.longread.map { sample, type, reads, truth_vcf, truth_bed -> tuple(sample, reads[0]) },
        file(params.reference)
    )
    sortBamLong(alignLongRead.out)
    callVariantsLong(sortBamLong.out, file(params.reference), faiReference.out.fai)

    // ── Shared QC
    vcfQc(
        callVariants.out.mix(callVariantsLong.out)
    )

    // ── Normalise all samples
    normaliseVcf(vcfQc.out.vcf, file(params.reference))

    // ── Annotate all samples
    annotateVcf(normaliseVcf.out.vcf, file(params.reference), faiReference.out.fai)

    // ── Benchmark only samples with a truth VCF defined
    branched.shortread
        .mix(branched.longread)
        .map { sample, type, reads, truth_vcf, truth_bed ->
            tuple(sample, truth_vcf, file(truth_vcf.toString() + ".tbi"), truth_bed)
        }
        .set { truth_ch }

    normaliseVcf.out.vcf
        .join(truth_ch)
        .filter { sample, vcf, truth_vcf, truth_tbi, truth_bed -> 
            truth_vcf != [] && truth_bed != [] }
        .view { "DEBUG benchmark input: $it" }
        .set { benchmark_input }

    benchmark(benchmark_input, file(params.reference), faiReference.out.fai)
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
    memory '32 GB'

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
    """
    minimap2 -ax map-ont -t ${task.cpus} \
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ONT" \
      ${reference} ${reads} > ${sample_id}.sam
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
      --sample_name=${sample_id} \
      --regions=chr20

    mv output/${sample_id}.vcf.gz ${sample_id}.vcf.gz
    """
}

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
    """
    /opt/bin/run_clair3.sh \
      --bam_fn=${bam_file} \
      --ref_fn=${reference} \
      --threads=${task.cpus} \
      --platform=ont \
      --model_path=/opt/models/r1041_e82_400bps_sup_v500 \
      --output=\$PWD/clair3_output \
      --sample_name=${sample_id} \
      --ctg_name=chr20

    mv clair3_output/merge_output.vcf.gz ${sample_id}.vcf.gz
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
    tuple val(sample_id), path("${sample_id}.final.vcf.gz"), emit: vcf
    path("${sample_id}.stats.txt"),                          emit: stats

    script:
    """
    bcftools view -f PASS ${vcf_file} -Oz -o ${sample_id}.final.vcf.gz
    bcftools stats ${sample_id}.final.vcf.gz > ${sample_id}.stats.txt
    """
}

// ── VCF Normalisation ─────────────────────────────────────────────────────────

process normaliseVcf {
    tag "${sample_id}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/normalised", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf_file)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.normalised.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.normalised.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools norm \
      -f ${reference} \
      -m -any \
      ${vcf_file} \
      -Oz -o ${sample_id}.normalised.vcf.gz

    bcftools index -t ${sample_id}.normalised.vcf.gz
    """
}

// ── VCF Annotation — VEP ─────────────────────────────────────────────────────

process annotateVcf {
    tag "${sample_id}"
    container 'ensemblorg/ensembl-vep:latest'
    containerOptions "-v ${params.vep_cache}:/opt/vep/.vep"
    publishDir "${params.outdir}/annotated", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), path(vcf_file)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz"),     emit: vcf
    tuple val(sample_id), path("${sample_id}.annotated.vcf.gz.tbi"), emit: tbi
    path("${sample_id}.vep_summary.html"),                            emit: summary

    script:
    """
    vep \
      --input_file ${vcf_file} \
      --output_file ${sample_id}.annotated.vcf \
      --format vcf \
      --vcf \
      --offline \
      --fasta ${reference} \
      --species homo_sapiens \
      --assembly GRCh38 \
      --cache \
      --dir_cache /opt/vep/.vep \
      --fork ${task.cpus} \
      --everything \
      --force_overwrite \
      --stats_file ${sample_id}.vep_summary.html
    
    bgzip ${sample_id}.annotated.vcf
    tabix -p vcf ${sample_id}.annotated.vcf.gz
    """
}

// ── Benchmarking — hap.py ────────────────────────────────────────────────────

process benchmark {
    tag "${sample_id}"
    container 'jmcdani20/hap.py:v0.3.12'
    publishDir "${params.outdir}/benchmark", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), path(vcf_file), path(truth_vcf), path(truth_tbi), path(truth_bed)
    path(reference)
    path(fai)

    output:
    path("${sample_id}.benchmark.summary.csv"), emit: summary
    path("${sample_id}.benchmark.metrics.json.gz"), emit: metrics
    path("${sample_id}.benchmark.extended.csv"),    emit: extended

    script:
    """
    echo -e "chr20\t0\t64444167" > chr20.bed

    /opt/hap.py/bin/hap.py \
      ${truth_vcf} \
      ${vcf_file} \
      -r ${reference} \
      -f ${truth_bed} \
      -o ${sample_id}.benchmark \
      -T chr20.bed \
      --threads ${task.cpus} \
      --engine=vcfeval
    """
}
