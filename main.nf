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
        def fastq_2          = row.fastq_2          ? file(row.fastq_2)     : []
        def aligned_bam      = row.aligned_bam      ? file(row.aligned_bam) : []
        def truth_vcf        = row.truth_vcf        ? file(row.truth_vcf)   : []
        def truth_bed        = row.truth_bed        ? file(row.truth_bed)   : []
        def clair3_platform  = row.clair3_platform  ?: ''
        def clair3_model     = row.clair3_model     ?: ''
        def dv_model         = row.dv_model         ?: ''
        def fastq_1          = row.fastq_1          ? file(row.fastq_1)     : []

        // Validation
        if (row.sequencing_type == 'longread' && !clair3_model) {
            error "Sample '${row.sample}' is longread but has no clair3_model declared"
        }
        if (row.sequencing_type == 'longread' && !clair3_platform) {
            error "Sample '${row.sample}' is longread but has no clair3_platform declared"
        }
        if (!fastq_1 && !aligned_bam) {
            error "Sample '${row.sample}' has neither fastq_1 nor aligned_bam declared"
        }

        tuple(row.sample, row.sequencing_type, [fastq_1, fastq_2], aligned_bam,
              truth_vcf, truth_bed, clair3_platform, clair3_model, dv_model)
    }
    .branch {
        shortread:    it[1] == 'shortread'
        longread_fq:  it[1] == 'longread' && it[3] == []
        longread_bam: it[1] == 'longread' && it[3] != []
    }
    .set { branched }

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {
    indexReference(file(params.reference))
    faiReference(file(params.reference))

    // ── Short-read path
    align(
        branched.shortread.map { it -> tuple(it[0], it[2]) },
        file(params.reference),
        indexReference.out.index
    )
    sortBam(align.out)

    // Pair sorted BAM with dv_model for short-read DV call
    sortBam.out
        .join(branched.shortread.map { it -> tuple(it[0], it[8]) })
        .set { dv_short_input }

    callVariants(dv_short_input, file(params.reference), faiReference.out.fai)

    // ── Long-read path: FASTQ input
    alignLongRead(
        branched.longread_fq.map { it -> tuple(it[0], it[2][0]) },
        file(params.reference)
    )
    sortBamLong(alignLongRead.out)

    // ── Long-read path: pre-aligned BAM input
    indexAlignedBam(
        branched.longread_bam.map { it -> tuple(it[0], it[3]) }
    )

    // Mix the two long-read paths back together (sample, bam, bai)
    sortBamLong.out
        .mix(indexAlignedBam.out)
        .set { longread_bams }

    // ── Clair3 calling on all long-read samples
    longread_bams
        .join(
            branched.longread_fq
                .mix(branched.longread_bam)
                .map { it -> tuple(it[0], it[6], it[7]) }
        )
        .set { clair3_input }

    callVariantsLong(clair3_input, file(params.reference), faiReference.out.fai)

    // ── DeepVariant calling on long-read samples that declare a dv_model
    longread_bams
        .join(
            branched.longread_fq
                .mix(branched.longread_bam)
                .map { it -> tuple(it[0], it[8]) }
        )
        .filter { sample, bam, bai, dv_model -> dv_model != '' }
        .set { dv_long_input }

    callVariantsLongDV(dv_long_input, file(params.reference), faiReference.out.fai)

    // ── Ensemble merges
    callVariantsLong.out
        .join(callVariantsLongDV.out)
        .set { ensemble_input }

    ensembleMergeLong(ensemble_input)
    ensembleUnionLong(ensemble_input, file(params.reference))

    // ── Tag VCFs with caller name for downstream per-caller routing
    callVariants.out
        .map { sample, vcf -> tuple(sample, 'deepvariant', vcf) }
        .set { tagged_short }

    callVariantsLong.out
        .map { sample, vcf -> tuple(sample, 'clair3', vcf) }
        .set { tagged_clair3 }

    callVariantsLongDV.out
        .map { sample, vcf -> tuple(sample, 'deepvariant_long', vcf) }
        .set { tagged_dv_long }

    ensembleMergeLong.out.vcf
        .map { sample, vcf -> tuple(sample, 'ensemble', vcf) }
        .set { tagged_ensemble }

    ensembleUnionLong.out.vcf
        .map { sample, vcf -> tuple(sample, 'ensemble_union', vcf) }
        .set { tagged_ensemble_union }

    tagged_short
        .mix(tagged_clair3)
        .mix(tagged_dv_long)
        .mix(tagged_ensemble)
        .mix(tagged_ensemble_union)
        .set { all_vcfs }

    vcfQc(all_vcfs)
    normaliseVcf(vcfQc.out.vcf, file(params.reference))
    annotateVcf(normaliseVcf.out.vcf, file(params.reference), faiReference.out.fai)

    // ── Benchmark
    branched.shortread
        .mix(branched.longread_fq)
        .mix(branched.longread_bam)
        .map { it -> tuple(it[0], it[4], file(it[4].toString() + ".tbi"), it[5]) }
        .set { truth_ch }

    normaliseVcf.out.vcf
        .combine(truth_ch, by: 0)
        .filter { sample, caller, vcf, truth_vcf, truth_tbi, truth_bed ->
            truth_vcf != [] && truth_bed != [] }
        .set { benchmark_input }

    benchmark(benchmark_input, file(params.reference), faiReference.out.fai)
}

// ── Reference Indexing ────────────────────────────────────────────────────────

process indexReference {
    tag "bwa-mem2 index"
    container 'dceoy/bwa-mem2:latest'
    containerOptions '--entrypoint='
    storeDir 'data/reference'
    input:  path(reference)
    output: path("${reference}.*"), emit: index
    script: "bwa-mem2 index ${reference}"
}

process faiReference {
    tag "samtools faidx"
    container 'staphb/samtools:latest'
    storeDir 'data/reference'
    input:  path(reference)
    output: path("${reference}.fai"), emit: fai
    script: "samtools faidx ${reference}"
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
    bwa-mem2 mem -t ${task.cpus} \\
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \\
      ${reference} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}

// ── Long-read Alignment (FASTQ input) ─────────────────────────────────────────

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
    minimap2 -ax map-ont -t ${task.cpus} \\
      -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ONT" \\
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

// ── Already-aligned BAM input — index only ────────────────────────────────────

process indexAlignedBam {
    tag "${sample_id}"
    container 'staphb/samtools:latest'
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: out

    script:
    """
    # Rename for consistency with downstream naming, then index
    cp ${bam_file} ${sample_id}.aligned.bam
    samtools index ${sample_id}.aligned.bam
    """
}

// ── Variant Calling ───────────────────────────────────────────────────────────

process callVariants {
    tag "${sample_id} (${dv_model})"
    container 'google/deepvariant:1.10.0'
    publishDir "${params.outdir}/vcf/${sample_id}_deepvariant", mode: 'copy'
    cpus 4
    memory '32 GB'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), val(dv_model)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), emit: out

    script:
    """
    mkdir -p output
    /opt/deepvariant/bin/run_deepvariant \\
      --model_type=${dv_model} \\
      --ref=${reference} \\
      --reads=${bam_file} \\
      --output_vcf=output/${sample_id}.deepvariant.vcf.gz \\
      --num_shards=${task.cpus} \\
      --sample_name=${sample_id} \\
      --regions=chr20

    mv output/${sample_id}.deepvariant.vcf.gz ${sample_id}.deepvariant.vcf.gz
    """
}

process callVariantsLong {
    tag "${sample_id} (${clair3_platform}/${clair3_model})"
    container 'hkubal/clair3:latest'
    publishDir "${params.outdir}/vcf/${sample_id}_clair3", mode: 'copy'
    cpus 4
    memory '8 GB'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), val(clair3_platform), val(clair3_model)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.clair3.vcf.gz"), emit: out

    script:
    """
    if [ ! -d "/opt/models/${clair3_model}" ]; then
        echo "ERROR: Clair3 model directory /opt/models/${clair3_model} not found in container" >&2
        echo "Available models:" >&2
        ls /opt/models/ >&2
        exit 1
    fi

    /opt/bin/run_clair3.sh \\
      --bam_fn=${bam_file} \\
      --ref_fn=${reference} \\
      --threads=${task.cpus} \\
      --platform=${clair3_platform} \\
      --model_path=/opt/models/${clair3_model} \\
      --output=\$PWD/clair3_output \\
      --sample_name=${sample_id} \\
      --ctg_name=chr20

    mv clair3_output/merge_output.vcf.gz ${sample_id}.clair3.vcf.gz
    """
}

process callVariantsLongDV {
    tag "${sample_id} (${dv_model})"
    container 'google/deepvariant:1.10.0'
    publishDir "${params.outdir}/vcf/${sample_id}_deepvariant_long", mode: 'copy'
    cpus 4
    memory '32 GB'

    input:
    tuple val(sample_id), path(bam_file), path(bai_file), val(dv_model)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant_long.vcf.gz"), emit: out

    script:
    """
    mkdir -p output
    /opt/deepvariant/bin/run_deepvariant \\
      --model_type=${dv_model} \\
      --ref=${reference} \\
      --reads=${bam_file} \\
      --output_vcf=output/${sample_id}.deepvariant_long.vcf.gz \\
      --num_shards=${task.cpus} \\
      --sample_name=${sample_id} \\
      --regions=chr20

    mv output/${sample_id}.deepvariant_long.vcf.gz ${sample_id}.deepvariant_long.vcf.gz
    """
}

process ensembleMergeLong {
    tag "${sample_id}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/vcf/${sample_id}_ensemble", mode: 'copy'

    input:
    tuple val(sample_id), path(clair3_vcf), path(dv_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.ensemble.vcf.gz"), emit: vcf

    script:
    """
    bcftools index -t -f ${clair3_vcf}
    bcftools index -t -f ${dv_vcf}
    bcftools isec -n=2 -p isec_out -Oz ${clair3_vcf} ${dv_vcf}
    cp isec_out/0000.vcf.gz ${sample_id}.ensemble.vcf.gz
    bcftools index -t ${sample_id}.ensemble.vcf.gz
    """
}

process ensembleUnionLong {
    tag "${sample_id}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/vcf/${sample_id}_ensemble_union", mode: 'copy'

    input:
    tuple val(sample_id), path(clair3_vcf), path(dv_vcf)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.ensemble_union.vcf.gz"), emit: vcf

    script:
    """
    bcftools norm -f ${reference} -m -any ${clair3_vcf} -Oz -o clair3.norm.vcf.gz
    bcftools norm -f ${reference} -m -any ${dv_vcf}     -Oz -o dv.norm.vcf.gz
    bcftools index -t clair3.norm.vcf.gz
    bcftools index -t dv.norm.vcf.gz

    bcftools concat -a clair3.norm.vcf.gz dv.norm.vcf.gz -Oz -o concat.vcf.gz
    bcftools sort concat.vcf.gz -Oz -o sorted.vcf.gz
    bcftools norm -d exact sorted.vcf.gz -Oz -o ${sample_id}.ensemble_union.vcf.gz
    bcftools index -t ${sample_id}.ensemble_union.vcf.gz
    """
}

// ── VCF QC, Normalisation, Annotation, Benchmarking ──────────────────────────

process vcfQc {
    tag "${sample_id}/${caller}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/qc/${caller}", mode: 'copy'

    input:
    tuple val(sample_id), val(caller), path(vcf_file)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.final.vcf.gz"), emit: vcf
    path("${sample_id}.${caller}.stats.txt"), emit: stats

    script:
    """
    bcftools view -f PASS ${vcf_file} -Oz -o ${sample_id}.${caller}.final.vcf.gz
    bcftools stats ${sample_id}.${caller}.final.vcf.gz > ${sample_id}.${caller}.stats.txt
    """
}

process normaliseVcf {
    tag "${sample_id}/${caller}"
    container 'staphb/bcftools:latest'
    publishDir "${params.outdir}/normalised/${caller}", mode: 'copy'

    input:
    tuple val(sample_id), val(caller), path(vcf_file)
    path(reference)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.normalised.vcf.gz"),     emit: vcf
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.normalised.vcf.gz.tbi"), emit: tbi

    script:
    """
    bcftools norm -f ${reference} -m -any ${vcf_file} -Oz -o ${sample_id}.${caller}.normalised.vcf.gz
    bcftools index -t ${sample_id}.${caller}.normalised.vcf.gz
    """
}

process annotateVcf {
    tag "${sample_id}/${caller}"
    container 'ensemblorg/ensembl-vep:latest'
    containerOptions "-v ${params.vep_cache}:/opt/vep/.vep"
    publishDir "${params.outdir}/annotated/${caller}", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), val(caller), path(vcf_file)
    path(reference)
    path(fai)

    output:
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.annotated.vcf.gz"),     emit: vcf
    tuple val(sample_id), val(caller), path("${sample_id}.${caller}.annotated.vcf.gz.tbi"), emit: tbi
    path("${sample_id}.${caller}.vep_summary.html"),                            emit: summary

    script:
    """
    vep \
      --input_file ${vcf_file} \
      --output_file ${sample_id}.${caller}.annotated.vcf \
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
      --stats_file ${sample_id}.${caller}.vep_summary.html
    
    bgzip ${sample_id}.${caller}.annotated.vcf
    tabix -p vcf ${sample_id}.${caller}.annotated.vcf.gz
    """
}

// ── Benchmarking — hap.py ────────────────────────────────────────────────────

process benchmark {
    tag "${sample_id}/${caller}"
    container 'jmcdani20/hap.py:v0.3.12'
    publishDir "${params.outdir}/benchmark/${caller}", mode: 'copy'
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample_id), val(caller), path(vcf_file), path(truth_vcf), path(truth_tbi), path(truth_bed)
    path(reference)
    path(fai)

    output:
    path("${sample_id}.${caller}.benchmark.summary.csv"), emit: summary
    path("${sample_id}.${caller}.benchmark.metrics.json.gz"), emit: metrics
    path("${sample_id}.${caller}.benchmark.extended.csv"),    emit: extended

    script:
    """
    echo -e "chr20\t0\t64444167" > chr20.bed

    /opt/hap.py/bin/hap.py \
      ${truth_vcf} \
      ${vcf_file} \
      -r ${reference} \
      -f ${truth_bed} \
      -o ${sample_id}.${caller}.benchmark \
      -T chr20.bed \
      --threads ${task.cpus} \
      --engine=xcmp
    """
}
