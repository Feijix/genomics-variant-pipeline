# genomics-variant-pipeline

Nextflow DSL2 pipeline for human variant calling. Processes paired-end short-read FASTQ files and produces filtered VCF files compatible with Beacon v2.

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
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr20.fna.gz
gunzip chr20.fna.gz
mv chr20.fna data/reference/chr20.fa
```

### 2. Get sample reads

Download real NA12878 reads from ENA (paired-end Illumina WGS, ~2 GB each):

```bash
curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
mv ERR194147_*.fastq.gz data/
```

### 3. Create the samplesheet

```bash
cat > data/samplesheet.csv << EOF
sample,fastq_1,fastq_2
NA12878,data/ERR194147_1.fastq.gz,data/ERR194147_2.fastq.gz
EOF
```

### 4. Run

```bash
nextflow run main.nf -profile docker
```

First run pulls Docker images and builds indexes (~10 min depending on connection). Subsequent runs skip both.

Expected output:

```
[skipped  ] indexReference  | 1 of 1, stored ✔
[skipped  ] faiReference    | 1 of 1, stored ✔
[xx/xxxxxx] align           | 1 of 1 ✔
[xx/xxxxxx] sortBam         | 1 of 1 ✔
[xx/xxxxxx] callVariants    | 1 of 1 ✔
[xx/xxxxxx] vcfQc           | 1 of 1 ✔
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

---

## Inspecting Results

`bcftools` and `samtools` are not required on the host. Use the same containers the pipeline uses:

```bash
# Add to ~/.bashrc for convenience
alias bcftools='docker run --rm -v $(pwd):$(pwd) -w $(pwd) staphb/bcftools:latest bcftools'
alias samtools='docker run --rm -v $(pwd):$(pwd) -w $(pwd) staphb/samtools:latest samtools'

source ~/.bashrc

# Variant summary
bcftools stats data/outputs/qc/NA12878.final.vcf.gz | grep "^SN"

# Ts/Tv ratio (expect 2.0–2.1 for WGS)
bcftools stats data/outputs/qc/NA12878.final.vcf.gz | grep "Ts/Tv"

# BAM alignment stats
samtools flagstat data/outputs/bam/NA12878.sorted.bam
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

**All variants filtered out (0 records in output)**
Expected behaviour with synthetic test data — there are no real variants to call. Use real FASTQ reads from GIAB to produce meaningful output.

---

## Status

- [x] Short-read pipeline (bwa-mem2 + DeepVariant WGS)
- [ ] Multi-sample support
- [ ] Long-read support (Clair3)
- [ ] VCF normalization and annotation (VEP)
- [ ] Beacon v2 output

## Validation

Designed for [GIAB](https://www.nist.gov/programs-projects/genome-bottle) sample NA12878 (HG001) on GRCh38.
Truth VCF: `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv3.3.2/GRCh38/`