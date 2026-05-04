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
- Benchmarking: `hap.py v0.3.12` with vcfeval engine against GIAB v4.2.1 truth sets
- Per-sample `truth_vcf` and `truth_bed` carried through samplesheet
- Filter strategy: `bcftools view -f PASS` (replaced earlier `QUAL>20`)

Validated across three technology generations. Headline benchmark table:

| Sample | Technology | Caller + model | SNP F1 | Indel F1 |
|---|---|---|---|---|
| HG001 chr20 (NYGC 30×) | Illumina 2×150bp | DeepVariant 1.5.0 WGS | 0.998 | 0.996 |
| HG002 chr20 (UCSC ultralong) | ONT R9.4.1, Guppy unknown | Clair3 `r941_prom_hac_g360+g422` | 0.987 | 0.215 |
| HG002 chr20 (GIAB v3.4.5) | ONT R9.4.1, Guppy 3.4.5 HAC | Clair3 `r941_prom_hac_g360+g422` | 0.995 | 0.541 |
| HG002 chr20 (GIAB 2025.01) | ONT R10.4.1, Dorado v5.0.0 SUP | Clair3 `r1041_e82_400bps_sup_v500` | 0.999 | 0.903 |

Full results in `BENCHMARKS.md`. Discussion of dataset effects, chemistry effects, and filter calibration in `METHODOLOGY.md` §9 and §16.

### Stage 4a — Short-read validation (HG001) ✓
Real data: NYGC 30× GRCh38 CRAM (ENA ERR3239334), chr20 streamed via HTTP range request. SNP and indel F1 both ~0.99 against GIAB v4.2.1.

### Stage 4b — Long-read R9 validation (HG002 ONT R9.4.1) ✓
Two ultralong R9 datasets tested to disambiguate dataset effects from tool effects. Indel F1 jumped from 0.22 (UCSC release, undocumented basecaller) to 0.54 (GIAB community release, Guppy 3.4.5 HAC) under identical pipeline configuration. Latter result is consistent with published Clair3-on-R9 expectations. Remaining gap to short-read attributed to R9 chemistry limits and mixed-flowcell composition of the GIAB release.

### Stage 4c — Long-read R10 validation (HG002 ONT R10.4.1) ✓
GIAB 2025.01 release, R10.4.1 PromethION, Dorado v5.0.0 SUP basecaller. Single PAW70337 flowcell, chr20 streamed from ONT open data S3. SNP F1 = 0.999 (matches short-read DeepVariant ceiling); indel F1 = 0.903. R10 chemistry + matched Dorado v5 SUP model + Clair3 produces production-grade calling on both variant types — the indel ceiling that limits R9 calling is largely lifted.

## Stage 5 — Beacon v2 output
- [ ] Convert annotated VCF to Beacon Friendly Format (BFF)
- [ ] Produce `genomicVariationsVcf.json.gz` per GA4GH Beacon v2 spec
- [ ] Validate against Beacon v2 reference implementation

## Stage 6 — Framework redesign
- [ ] Per-sample technology metadata in samplesheet (or auto-detected from BAM headers)
- [ ] Externally configurable mapping from data attributes → caller + model
- [ ] Documented mechanism for routing decision updates as tools evolve
- [ ] Caller selection logic moved out of `main.nf` into `conf/routing.config` or equivalent