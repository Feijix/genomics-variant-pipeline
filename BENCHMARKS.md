# Benchmarks

Single source of truth for validated benchmark results. Each row is a complete, reproducible experiment: data source, pipeline configuration, and measured accuracy against a GIAB truth set.

All results below: chr20 only, GIAB v4.2.1 high-confidence regions (`*_noinconsistent.bed`), `hap.py v0.3.12` with vcfeval engine, normalised query VCFs (`bcftools norm -m -any`).

## Validated runs

### Short-read

| Sample | Technology | Source | Coverage | Caller + model | Filter | SNP recall | SNP precision | SNP F1 | Indel recall | Indel precision | Indel F1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HG001 | Illumina 2×150bp | NYGC 30× CRAM (ENA ERR3239334), chr20 streamed via HTTP range request | 36.3× | DeepVariant 1.5.0 (`--model_type=WGS`) | `FILTER=PASS` | 0.9956 | 0.9994 | 0.998 | 0.9950 | 0.9977 | 0.996 |

### Long-read ONT R9.4.1

| Sample | Source | Basecaller | Coverage (chr20) | Caller + model | Filter | SNP F1 | Indel F1 | Notes |
|---|---|---|---|---|---|---|---|---|
| HG002_ONT (UCSC) | UCSC ultralong PromethION (May 2020 release) | Unknown Guppy version | 46.5× | Clair3 `r941_prom_hac_g360+g422` | `FILTER=PASS` | 0.987 | 0.215 | Indel F1 below published expectations; attributed to undocumented basecaller provenance and ultralong-class read length distribution |
| HG002_ONT (UCSC) | as above | as above | as above | as above | `FILTER=PASS + QUAL≥15` (offline sweep) | 0.987 | 0.371 | Operating point selected from QUAL sweep peak; precision-recall trade-off |
| HG002_ONT_giab | GIAB community release (Feb 2020) | Guppy 3.4.5 HAC | ~52× | Clair3 `r941_prom_hac_g360+g422` | `FILTER=PASS` | 0.995 | 0.541 | Within published Clair3-on-R9 expectations |
| HG002_ONT_giab | as above | as above | as above | as above | `FILTER=PASS + QUAL≥15` (offline) | 0.995 | 0.534 | Same F1, shifts to high-precision/low-recall operating point |

### Long-read ONT R10.4.1

| Sample | Source | Basecaller | Coverage (chr20) | Caller + model | Filter | SNP F1 | Indel F1 | Notes |
|---|---|---|---|---|---|---|---|---|
| HG002_ONT_R10 | GIAB 2025.01 release, PAW70337 flowcell, chr20 streamed from `s3://ont-open-data/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam` | Dorado v5.0.0 SUP (`dna_r10.4.1_e8.2_400bps@v5.0.0`) | [VERIFY: run `samtools depth` on the new BAM] | Clair3 `r1041_e82_400bps_sup_v500` | `FILTER=PASS` | 0.999 | 0.903 | Production-grade calling on both variant types; SNP F1 matches short-read ceiling |

## Detailed results

### HG001 short-read (HG001.benchmark.summary.csv)

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | F1 | Het/hom (truth → query) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| INDEL | PASS | 9974 | 9924 | 50 | 19819 | 23 | 9882 | 18 | 1 | 0.9950 | 0.9977 | 0.996 | 1.630 → 1.783 |
| SNP | PASS | 69175 | 68874 | 301 | 84183 | 39 | 15226 | 7 | 1 | 0.9956 | 0.9994 | 0.998 | 1.730 → 1.723 |

### HG002 ONT R9.4.1, GIAB v3.4.5 release (HG002_ONT_giab.benchmark.summary.csv)

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | F1 | Het/hom |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| INDEL | PASS | 11256 | 5957 | 5299 | 15221 | 4808 | 4470 | 318 | 925 | 0.5292 | 0.5528 | 0.541 | 1.562 → 1.087 |
| SNP | PASS | 71333 | 70897 | 436 | 110708 | 222 | 39563 | 20 | 165 | 0.9939 | 0.9969 | 0.995 | 1.716 → 1.727 |

### HG002 ONT R10.4.1 (HG002_ONT_R10.benchmark.summary.csv)

| Type | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | Recall | Precision | F1 | Het/hom |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| INDEL | PASS | 11256 | 9927 | 1329 | 20117 | 803 | 9396 | 329 | 292 | 0.8819 | 0.9251 | 0.903 | 1.562 → 1.795 |
| SNP | PASS | 71333 | 71244 | 89 | 123689 | 38 | 52353 | 13 | 16 | 0.9988 | 0.9995 | 0.999 | 1.716 → 2.104 |

## How to reproduce

For each benchmark above, the inputs (FASTQ source URL, samplesheet entry, Clair3 model, filter applied) are sufficient to reproduce the result on any machine running the pipeline at the documented commit. Pipeline commit hashes for each run can be added once results are committed to the repository.

Somatic (SEQC2 HCC1395/HCC1395BL chr20, truth v1.2.1: 1106 SNV / 42 indel)
  rep1 T/N   → data/outputs/benchmark/somatic.rep1_paired.bak/
  tumour-only→ data/outputs/benchmark/somatic.tumoronly.bak/
  rep2 T/N   → data/outputs/benchmark/somatic/WGS_IL_T_2.*
  cmd: nextflow run main.nf -profile docker --samplesheet <sheet> \
       --callers mutect2,strelka2,manta_somatic,deepvariant --skip_bwa_index true