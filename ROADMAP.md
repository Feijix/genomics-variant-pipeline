# Roadmap

## Stage 1 — Short-read pipeline, chr20 ✓
- End-to-end pipeline: FASTQ → filtered VCF
- bwa-mem2 + DeepVariant WGS + bcftools
- Automatic reference indexing with storeDir caching
- Validated on HG001 chr20 synthetic data

## Stage 2 — Multi-sample parallel processing ✓
- N samples run concurrently via Nextflow channel parallelism
- Read group tags (@RG: SM, ID, PL) on all alignments
- Validated on two synthetic samples (HG001/HG002): ~47,000 variants each, Ts/Tv = 0.50 (expected)

## Stage 3 — Long-read support ✓
- minimap2 (`map-ont`) for ONT alignment; Clair3 v2 for variant calling
- `branch` operator routes samples by `sequencing_type` at load time
- Paths converge at QC via `.mix()`
- Validated: synthetic ONT sample (HG003) produced 7 calls — Clair3 correctly suppresses nanopore noise

## Stage 4 — Real data validation ✓
- VCF normalisation: `bcftools norm -m -any` (left-align indels, split multiallelics)
- VEP annotation: consequence, gene, HGVS, gnomAD frequencies, ClinVar
- Benchmarking: `hap.py v0.3.12` with vcfeval engine against GIAB v4.2.1 truth set
- Real data sourced from NYGC 30× GRCh38 CRAM (ENA ERR3239334), chr20 streamed via HTTP range request

### Benchmark results — HG001 chr20 vs GIAB v4.2.1

| Run | Coverage | SNP recall | SNP precision | Indel recall | Indel precision |
|---|---|---|---|---|---|
| Initial (partial FASTQ) | 3.66× | 51.9% | 89.9% | 39.1% | 86.0% |
| NYGC 30× CRAM | 36.3× | 98.95% | 99.97% | 98.97% | 99.82% |

Het/hom ratio recovered from 0.75 → 1.747 (truth: 1.730). Ts/Tv = 2.175.

## Stage 5 — Beacon v2 output
- [ ] Convert annotated VCF to Beacon Friendly Format (BFF)
- [ ] Produce `genomicVariationsVcf.json.gz` per GA4GH Beacon v2 spec
- [ ] Validate against Beacon v2 reference implementation