# Roadmap

## Stage 1 — Short-read pipeline, chr20 ✓
- End-to-end pipeline: FASTQ → filtered VCF
- bwa-mem2 + DeepVariant WGS + bcftools
- Automatic reference indexing with storeDir caching
- Validated on NA12878 chr20 synthetic data

## Stage 2 — Multi-sample + real data
- Download real NA12878 short-read chr20 FASTQs (GIAB)
- Add second short-read sample (NA19240)
- bcftools merge + comparison report
- Validate Ts/Tv ratios against GIAB truth set

## Stage 3 — Long-read support
- Add Clair3 process (`callVariantsLongRead`)
- Add `sequencing_type` column to samplesheet
- Conditional workflow routing (short vs long)
- Test on publicly available HiFi chr20 data

## Stage 4 — Real data validation + Beacon v2
- [x] Download real GIAB HG001 reads (NYGC 30× GRCh38 CRAM, ENA ERR3239334)
- [x] VCF normalisation with `bcftools norm`
- [x] VEP annotation (gene + consequence + gnomAD frequencies)
- [x] Benchmarking with `hap.py` v0.3.12 (vcfeval engine)
- [x] Low-coverage run (3.66×): Ts/Tv 2.17 ✓, SNP recall 51.9% (coverage-limited, documented)
- [ ] 30× run in progress — SNP recall target > 99%
- [ ] Beacon v2 `genomicVariations` output