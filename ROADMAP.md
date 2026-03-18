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

## Stage 4 — Whole genome + Beacon v2
- Scale to full GRCh38 (all chromosomes)
- bcftools norm on both VCFs (left-align indels, split multiallelics)
- VEP annotation (gene + consequence)
- Beacon v2 `genomicVariations` output