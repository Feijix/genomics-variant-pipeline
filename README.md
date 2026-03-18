# genomics-variant-pipeline

Nextflow DSL2 pipeline for human variant calling. Processes paired-end short-read FASTQ files and produces filtered VCF files compatible with Beacon v2. Supports multiple samples in parallel — each sample produces an independent VCF.

## Pipeline

```
FASTQ → align (bwa-mem2) → sort & index (samtools) → call variants (DeepVariant) → filter (bcftools) → VCF
```

All bioinformatics tools run inside Docker containers. Reference indexes (`bwt.2bit.64`, `fai`) and BAM indexes (`bai`) are built automatically on first run and cached — subsequent runs skip indexing entirely.

---

## Requirements

Only three things need to be installed on the host machine:

| Dependency | Version tested | Notes |
|---|---|---|
| Java | ≥ 11 | Required by Nextflow |
| Nextflow | 25.10.4 | See install below |
| Docker | any recent | All bioinformatics tools run inside containers |

> **Note:** Use `curl` rather than `wget` to download reference files — some `wget` builds lack FTP support.

### Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/          # or any directory on your $PATH
nextflow -version           # verify
```

### Verify Docker

```bash
docker run --rm hello-world
```

If this fails, Docker is not running or your user lacks permission (`sudo usermod -aG docker $USER`, then log out and back in).

---

## Quickstart

### 1. Get the reference genome

```bash
mkdir -p data/reference
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
gunzip chr20.fa.gz
mv chr20.fa data/reference/chr20.fa
```

### 2. Generate synthetic sample reads

`wgsim` ships inside the samtools container — no extra install needed. The `-S` flag sets the random seed, giving each sample a distinct set of variants.

```bash
mkdir -p data/samples

# HG001 — seed 42
docker run --rm -v $(pwd)/data:$(pwd)/data -w $(pwd)/data staphb/samtools:latest \
  wgsim -S 42 -N 1000000 -1 150 -2 150 -e 0.001 -r 0.001 -R 0.0001 \
  reference/chr20.fa \
  samples/HG001_R1.fastq samples/HG001_R2.fastq

# HG002 — seed 99
docker run --rm -v $(pwd)/data:$(pwd)/data -w $(pwd)/data staphb/samtools:latest \
  wgsim -S 99 -N 1000000 -1 150 -2 150 -e 0.001 -r 0.001 -R 0.0001 \
  reference/chr20.fa \
  samples/HG002_R1.fastq samples/HG002_R2.fastq

# Compress
gzip data/samples/HG001_R1.fastq data/samples/HG001_R2.fastq \
     data/samples/HG002_R1.fastq data/samples/HG002_R2.fastq
```

> **Note on synthetic data:** wgsim places variants randomly with no biological model. Ts/Tv ratio will be ~0.50 (random chance) rather than the ~2.0–2.1 expected from real human WGS data. This is expected and correct — use GIAB reads for biological validation.

### 3. Create the samplesheet

```bash
cat > data/samplesheet.csv << EOF
sample,fastq_1,fastq_2
HG001,data/samples/HG001_R1.fastq.gz,data/samples/HG001_R2.fastq.gz
HG002,data/samples/HG002_R1.fastq.gz,data/samples/HG002_R2.fastq.gz
EOF
```

Add or remove rows to process more or fewer samples. Each row runs independently through the full pipeline in parallel.

### 4. Run

```bash
nextflow run main.nf -profile docker
```

First run pulls Docker images and builds reference indexes (~10 min depending on connection). Subsequent runs skip both.

Expected output:

```
[skipped  ] indexReference  | 1 of 1, stored ✔
[skipped  ] faiReference    | 1 of 1, stored ✔
[xx/xxxxxx] align           | 2 of 2 ✔
[xx/xxxxxx] sortBam         | 2 of 2 ✔
[xx/xxxxxx] callVariants    | 2 of 2 ✔
[xx/xxxxxx] vcfQc           | 2 of 2 ✔
```

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--samplesheet` | `data/samplesheet.csv` | Sample manifest |
| `--reference` | `data/reference/chr20.fa` | Reference FASTA |
| `--outdir` | `data/outputs` | Output directory |

Override at runtime:

```bash
nextflow run main.nf -profile docker \
  --reference data/reference/GRCh38.fa \
  --outdir results/
```

---

## Output

```
data/outputs/
├── bam/
│   ├── <sample>.sorted.bam
│   └── <sample>.sorted.bam.bai
├── vcf/
│   └── <sample>.vcf.gz
└── qc/
    ├── <sample>.final.vcf.gz   ← QUAL > 20 filtered
    └── <sample>.stats.txt
```

One set of outputs per sample, named by the `sample` column in the samplesheet.

---

## Inspecting Results

`bcftools` and `samtools` are not required on the host. Use the same containers the pipeline uses:

```bash
# Add to ~/.bashrc for convenience
alias bcftools='docker run --rm -v $(pwd):$(pwd) -w $(pwd) staphb/bcftools:latest bcftools'
alias samtools='docker run --rm -v $(pwd):$(pwd) -w $(pwd) staphb/samtools:latest samtools'

source ~/.bashrc
```

```bash
# Variant counts per sample
bcftools stats data/outputs/qc/HG001.final.vcf.gz | grep "^SN"
bcftools stats data/outputs/qc/HG002.final.vcf.gz | grep "^SN"

# Ts/Tv ratio
# ~0.50 expected for synthetic data, ~2.0-2.1 for real WGS data
bcftools stats data/outputs/qc/HG001.final.vcf.gz | grep "^TSTV"
bcftools stats data/outputs/qc/HG002.final.vcf.gz | grep "^TSTV"

# BAM alignment stats
samtools flagstat data/outputs/bam/HG001.sorted.bam
samtools flagstat data/outputs/bam/HG002.sorted.bam
```

---

## Containers

| Process | Image |
|---|---|
| indexReference | `dceoy/bwa-mem2:latest` |
| faiReference | `staphb/samtools:latest` |
| align | `dceoy/bwa-mem2:latest` |
| sortBam | `staphb/samtools:latest` |
| callVariants | `google/deepvariant:1.5.0` |
| vcfQc | `staphb/bcftools:latest` |

---

## Troubleshooting

**`ERROR! Unable to open the file: chr20.fa.bwt.2bit.64`**
The reference index is missing. The pipeline builds it automatically — if this error appears, the `indexReference` process failed. Check the `work/` logs.

**`ValueError: NOT_FOUND: could not load fasta and/or fai`**
The `.fai` index was not staged into DeepVariant's work directory. Ensure `faiReference` completed successfully before `callVariants` — handled automatically in the current pipeline version.

**`[E::idx_find_and_load] Could not retrieve index file for '*.bam'`**
The `.bai` was not passed to `callVariants`. Ensure `sortBam` emits both the `.bam` and `.bai` in its output tuple.

**`No non-empty sample name found in the input reads`**
Read group tags are missing from the BAM. The `align` process adds them via `-R "@RG\tID:...\tSM:...\tPL:ILLUMINA"`. If you see this warning, check that the `-R` flag is present in the align script.

**Ts/Tv is ~0.50 instead of ~2.0**
Expected when using synthetic data generated by `wgsim`. This is not a pipeline error — wgsim places variants randomly with no biological model. Run with real GIAB reads to get a biologically valid Ts/Tv ratio.

---

## Status

- [x] Short-read pipeline (bwa-mem2 + DeepVariant WGS)
- [x] Multi-sample support — N samples run in parallel, one VCF each
- [ ] Long-read support (Clair3)
- [ ] VCF normalization and annotation (VEP)
- [ ] Beacon v2 output

---

## Validation

Tested with synthetic chr20 reads generated by `wgsim` from GRCh38 chr20 (UCSC hg38).
Two samples (HG001/HG002) with different random seeds confirmed independent variant sets.

For biological validation, use [GIAB](https://www.nist.gov/programs-projects/genome-bottle) truth sets:
- HG001 (NA12878): `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv3.3.2/GRCh38/`
- HG002 (NA24385): `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/`