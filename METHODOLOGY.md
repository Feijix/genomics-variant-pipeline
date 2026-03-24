# Methodology

## Pipeline Overview

This project implements a reproducible, containerised variant calling pipeline for human genomic data using Nextflow DSL2. The pipeline processes raw sequencing reads — both short-read (Illumina) and long-read (Oxford Nanopore) — and produces filtered, annotated Variant Call Format (VCF) files compatible with the Beacon v2 genomic data discovery protocol.

The core engineering objective is reproducibility: any researcher with Docker and Nextflow installed should be able to reproduce identical results from the same input data, on any machine, without manual environment configuration.

---

## 1. Workflow Orchestration — Nextflow DSL2

**Nextflow** (v25.10.4) was selected as the workflow management system for three reasons relevant to a computer science context:

- **Dataflow programming model** — processes communicate exclusively through typed channels, eliminating shared state and race conditions. Each process is a pure function: given the same inputs, it produces the same outputs.
- **Implicit parallelism** — the runtime automatically parallelises independent processes. When N samples are loaded from a samplesheet, all N execute concurrently without explicit thread management.
- **Resumability** — Nextflow maintains a content-addressable cache of process outputs keyed by input hash. Rerunning a modified pipeline skips unchanged processes, reducing iteration time from hours to minutes.

The pipeline is written in DSL2, which enforces strict process modularity — each process declares its inputs and outputs explicitly, making data dependencies machine-verifiable rather than implicit.

### Channel Branching

From Stage 3 onward, the pipeline routes samples to different variant callers based on sequencing technology. This is implemented using Nextflow's `branch` operator, which partitions a channel into named sub-channels based on a predicate:

```groovy
samples_ch
    .branch {
        shortread: it[1] == 'short'
        longread:  it[1] == 'long'
    }
    .set { branched }
```

The branch labels `shortread` and `longread` were chosen deliberately — `short` and `long` are primitive type names in Groovy (the language Nextflow is built on) and cannot be used as identifiers. This is a non-obvious constraint that manifests as a parser error rather than a semantic one.

The two paths converge at the QC stage using the `.mix()` operator, which merges channels without ordering guarantees — appropriate here since QC is stateless with respect to sequencing technology.

---

## 2. Containerisation — Docker

Every process in the pipeline runs inside a dedicated Docker container. No bioinformatics tool is installed on the host machine. This design decision has two engineering justifications:

- **Dependency isolation** — tool versions are pinned at the container level, not the system level. `google/deepvariant:1.5.0` will always run DeepVariant 1.5.0 regardless of host OS or installed libraries.
- **Portability** — the pipeline runs identically on a laptop, a university HPC node, or a cloud VM, as long as Docker is available.

Containers are declared per-process in Nextflow, so different tools with conflicting dependencies can coexist in the same pipeline without virtual environments or manual path management.

### SELinux Compatibility

On SELinux-enabled systems (Fedora, RHEL, CentOS), Docker volume mounts require an additional relabeling flag. Nextflow's `mountFlags = 'z'` configuration applies the `:z` SELinux label to all container volume mounts, allowing the container process to read host-mounted files without modifying file ownership. This flag has no effect on non-SELinux systems and is safe to include unconditionally.

### Containers Used

| Process | Image | Purpose |
|---|---|---|
| indexReference | `dceoy/bwa-mem2:latest` | BWA-MEM2 FM-index construction |
| faiReference | `staphb/samtools:latest` | FASTA index for random access |
| align | `dceoy/bwa-mem2:latest` | Short-read alignment |
| alignLongRead | `staphb/minimap2:latest` | Long-read alignment |
| sortBam | `staphb/samtools:latest` | BAM coordinate sorting and indexing |
| callVariants | `google/deepvariant:1.5.0` | Short-read variant calling (DNN) |
| callVariantsLong | `hkubal/clair3:latest` | Long-read variant calling (DNN) |
| vcfQc | `staphb/bcftools:latest` | VCF filtering and statistics |
| normaliseVcf | `staphb/bcftools:latest` | VCF normalisation |
| annotateVcf | `ensemblorg/ensembl-vep:latest` | Variant annotation |
| benchmark | `pkrusche/hap.py:latest` | Precision/recall benchmarking |
| beaconConvert | `beacon2ri/beacon2-ri-tools:latest` | Beacon v2 BFF conversion |

---

## 3. Reference Genome

**Reference:** GRCh38 (hg38), chromosome 20 (`chr20.fa`), sourced from UCSC Genome Browser golden path.

GRCh38 is the current human reference assembly maintained by the Genome Reference Consortium. Chromosome 20 (~63 MB) was used for development and validation to reduce compute time while preserving realistic genomic context — repetitive elements, GC variation, and structural complexity are all present in chr20 at representative levels.

### Reference Indexing

Two independent index formats are required by different tools:

**BWA-MEM2 FM-index** (`.bwt.2bit.64`, `.pac`, `.ann`, `.amb`, `.0123`) — a compressed suffix array enabling O(1) seed lookup during alignment. Construction is computationally expensive (~5 min for chr20, ~45 min for full GRCh38) and deterministic, so it is cached using Nextflow's `storeDir` directive and skipped on subsequent runs.

**FASTA index** (`.fai`) — a tabix-style index mapping chromosome names to byte offsets in the FASTA file. Required by DeepVariant and Clair3 for random-access region queries during variant calling.

Both indexes must be staged into each process's isolated working directory as explicit Nextflow inputs, since processes cannot access the host filesystem directly.

---

## 4. Short-Read Alignment — BWA-MEM2

**Tool:** BWA-MEM2 (Burrows-Wheeler Aligner, MEM algorithm, v2)

BWA-MEM2 aligns paired-end Illumina reads (2 × 150 bp) to the reference genome using a seed-and-extend strategy:

1. **Seeding** — exact matches (MEMs, Maximal Exact Matches) are found using the FM-index
2. **Chaining** — seeds are chained into candidate alignment regions
3. **Extension** — Smith-Waterman dynamic programming extends chains into full alignments

BWA-MEM2 is the vectorised successor to BWA-MEM, exploiting AVX-512 SIMD instructions for 2–3× speedup on modern CPUs. The pipeline detected AVX-512 support automatically on the development machine (Intel i7-11800H).

### Read Group Tags

Every alignment includes a read group (`@RG`) header tag:

```
@RG  ID:<sample>  SM:<sample>  PL:ILLUMINA
```

The `SM` (sample name) field propagates through to the VCF genotype column header. Without it, variant callers produce anonymous genotype columns that cannot be disambiguated in multi-sample analyses. This is a mandatory requirement for any downstream tool that parses sample identity, including Beacon v2.

---

## 5. Long-Read Alignment — Minimap2

**Tool:** minimap2 (`staphb/minimap2:latest`)

Long reads (Oxford Nanopore, mean length ~10,000 bp) require a different alignment strategy. BWA-MEM2 is optimised for short, low-error reads and fails on long reads with high indel error rates. Minimap2 uses a minimiser-based seeding strategy that scales to reads of arbitrary length.

The `map-ont` preset was used, which configures:
- Minimiser window size and k-mer length appropriate for ONT error rates
- Gap penalty parameters tuned for the systematic insertion/deletion errors characteristic of nanopore sequencing
- Spliced alignment disabled (genomic DNA, not RNA)

Output is in SAM format (`-a` flag) for compatibility with the downstream `samtools sort` step.

---

## 6. BAM Processing — SAMtools

**Tool:** SAMtools (`staphb/samtools:latest`)

Post-alignment processing involves two operations:

**Coordinate sorting** (`samtools sort`) — reorders alignment records by genomic position. Variant callers require sorted BAMs to perform efficient pileup operations — iterating over all reads covering a given position requires sequential access to a coordinate-sorted file.

**BAM indexing** (`samtools index`) — produces a `.bai` index enabling O(log n) random access to any genomic region. This is required by both DeepVariant and Clair3, which parallelise variant calling by splitting the genome into shards and processing each independently.

The `.bai` file is declared as an explicit Nextflow output and staged into the variant calling process's working directory alongside the BAM. Without explicit staging, Nextflow's process isolation would prevent the caller from finding the index.

---

## 7. Short-Read Variant Calling — DeepVariant

**Tool:** DeepVariant 1.5.0 (`google/deepvariant:1.5.0`, `--model_type=WGS`)

DeepVariant reformulates variant calling as an image classification problem. For each candidate variant site, it:

1. Constructs a pileup image — a matrix encoding read base, base quality, mapping quality, strand, and other features for all reads covering the site
2. Passes the image through a convolutional neural network (Inception v3 architecture) trained on validated variant calls
3. Outputs a genotype likelihood: homozygous reference, heterozygous, or homozygous alternate

This approach captures complex read patterns that rule-based callers (e.g. GATK HaplotypeCaller) cannot model, particularly in low-complexity regions and near structural variants.

Parallelism is controlled via `--num_shards`, set to `task.cpus` (4) to match the Nextflow resource declaration.

---

## 8. Long-Read Variant Calling — Clair3

**Tool:** Clair3 v2 (`hkubal/clair3:latest`, `--platform=ont`)

Clair3 is a deep-learning variant caller designed specifically for long reads. Like DeepVariant, it uses a neural network, but with architecture and training data matched to long-read error profiles.

Key differences from DeepVariant relevant to computer science:

- **Technology-specific models** — separate models are trained for PacBio HiFi, Oxford Nanopore R9, and ONT R10. The `ont` model was used here, trained on nanopore data with characteristic homopolymer errors.
- **Pileup + full-alignment hybrid** — Clair3 uses a two-stage approach: a fast pileup model screens candidate sites, then a slower full-alignment model evaluates candidates using the complete read sequence rather than a fixed window.
- **PyTorch backend** (v2.x) — Clair3 v2 migrated from TensorFlow to PyTorch, making v1 models incompatible with v2. Model and binary versions must be matched.

Clair3 requires absolute paths for input and output. Inside Nextflow's isolated work directories, `$PWD` is used to resolve the absolute path at runtime.

---

## 9. VCF Quality Control and Filtering — BCFtools

**Tool:** BCFtools (`staphb/bcftools:latest`)

Post-calling filtering applies a minimum quality threshold:

```bash
bcftools view -i 'QUAL>20' input.vcf.gz -Oz -o output.final.vcf.gz
```

`QUAL > 20` corresponds to a Phred-scaled probability of error < 1% (10^(−20/10) = 0.01). This removes low-confidence calls while retaining the majority of true variants in well-covered regions.

Quality statistics are computed with `bcftools stats`, which reports:
- Total variant counts (SNPs, indels, MNPs)
- Transition/transversion ratio (Ts/Tv)
- Per-sample genotype counts

### Ts/Tv Ratio as a Quality Metric

The transition/transversion ratio is the primary QC metric for variant calling pipelines. Transitions (purine↔purine, pyrimidine↔pyrimidine: A↔G, C↔T) are biochemically more probable than transversions (purine↔pyrimidine). In real human genomes:

- **Whole genome:** Ts/Tv ≈ 2.0–2.1
- **Whole exome:** Ts/Tv ≈ 2.6–3.0 (coding regions are more constrained)
- **Random noise:** Ts/Tv ≈ 0.5 (4 transitions vs 8 transversions by chance)

Synthetic data generated by `wgsim` and `badread` produces Ts/Tv ≈ 0.5, confirming that the simulators apply no biological mutation model. This is expected and used as a negative control to distinguish pipeline mechanics from biological signal.

---

## 10. VCF Normalisation — BCFtools norm

**Tool:** BCFtools (`staphb/bcftools:latest`, `bcftools norm`)

VCFs produced by different variant callers may represent the same biological variant using different coordinate conventions. Two callers processing the same alignment can produce records that are semantically identical but syntactically different, causing downstream tools that rely on exact coordinate matching to report false positives and false negatives where none exist.

Normalisation resolves this with two operations:

**Left-alignment of indels** — an insertion in a homopolymer run (e.g. `AAAAAAA`) can be placed at any position within the run and remain biologically equivalent. By convention, indels are shifted as far left as possible. Without left-alignment, the same indel at two slightly different coordinates appears as a mismatch during benchmarking.

**Splitting multiallelic sites** — a position with multiple alternate alleles (e.g. `REF=A, ALT=G,T`) may be represented as a single VCF row or as two separate rows depending on the caller. `bcftools norm -m -any` splits all multiallelic records into one row per alternate allele, ensuring consistent representation across callers and downstream tools.

```bash
bcftools norm -f reference.fa -m -any input.vcf.gz -Oz -o normalised.vcf.gz
```

Normalisation is a prerequisite for both `hap.py` benchmarking and Beacon v2 ingestion. Both tools perform coordinate-based matching and require canonical variant representation to function correctly.

---

## 11. Variant Annotation — Ensembl VEP

**Tool:** Ensembl Variant Effect Predictor (`ensemblorg/ensembl-vep:latest`)

Raw VCF files produced by variant callers contain only genomic coordinates and allele information. They carry no biological context — no gene name, no functional consequence, no population frequency. Variant annotation adds this information to the VCF INFO field by cross-referencing each variant against curated databases.

VEP annotates each variant with:

- **Gene** — Ensembl gene ID and HGNC symbol for the gene the variant falls in
- **Consequence** — the predicted functional effect: synonymous, missense, frameshift, splice site, intergenic, etc.
- **Impact** — a severity classification: HIGH (frameshift, stop gained), MODERATE (missense), LOW (synonymous), MODIFIER (intergenic)
- **Existing variant** — cross-reference against dbSNP and ClinVar if the variant has been previously reported
- **Population frequency** — allele frequency in gnomAD, 1000 Genomes, and other population databases

Annotation is mandatory for Beacon v2 ingestion. Without it, the `molecularEffects` fields in the `genomicVariations` entity are empty, producing a technically valid but scientifically limited Beacon record. A researcher querying the Beacon for variants with a specific consequence — for example, all frameshift mutations in BRCA1 — would not find unannotated variants regardless of their coordinates.

---

## 12. Benchmarking — hap.py

**Tool:** hap.py (`pkrusche/hap.py:latest`)

Variant calling accuracy is quantified by comparing pipeline output against a validated truth set. The Genome in a Bottle (GIAB) consortium maintains truth VCFs for a set of well-characterised human samples, including HG001 (NA12878), for which variant calls have been validated through orthogonal sequencing technologies and manual curation.

`hap.py` compares a query VCF (pipeline output) against a truth VCF (GIAB benchmark) and computes:

- **Precision** = TP / (TP + FP) — the fraction of called variants that are real
- **Recall** = TP / (TP + FN) — the fraction of real variants that were called
- **F1 score** = harmonic mean of precision and recall

The comparison uses a graph-based matching algorithm rather than exact string matching. This correctly handles equivalent representations of the same variant — a key reason why normalisation must precede benchmarking. A variant that is left-aligned in the truth set but not in the query would be counted as a false positive and false negative simultaneously without normalisation.

Benchmarking is restricted to the GIAB high-confidence regions (`.bed` file) — genomic intervals where the truth set is considered complete and reliable. Variants outside these regions are excluded from the precision/recall calculation.

Expected performance for DeepVariant on real Illumina WGS data against the GIAB HG001 truth set:
- SNP precision > 99%, recall > 99%
- Indel precision > 95%, recall > 95%

---

## 13. Beacon v2 Output — beacon2-ri-tools

**Tool:** beacon2-ri-tools (`beacon2ri/beacon2-ri-tools:latest`)

Beacon v2 is a GA4GH-approved standard for federated genomic data discovery. It defines an API protocol that allows researchers to query variant databases across institutions without raw data leaving each institution's infrastructure.

The `genomicVariations` entity in the Beacon v2 data model requires more than variant coordinates. Each record links a variant to the individual it came from, the biosample it was extracted from, the sequencing run, and the analysis pipeline. A raw VCF row maps to only a fraction of the required fields.

`beacon2-ri-tools` converts an annotated VCF into Beacon Friendly Format (BFF) — a set of JSON files structured according to the Beacon v2 data model. The primary output for variant data is `genomicVariationsVcf.json.gz`, which contains one JSON record per variant with fields populated from the VCF INFO annotations produced by VEP.

The conversion pipeline is:

```
annotated VCF → beacon2-ri-tools → genomicVariationsVcf.json.gz
```

This output can be loaded directly into a Beacon v2 MongoDB instance using the standard reference implementation ingestion scripts. For this project, producing a correctly formatted and validated `genomicVariationsVcf.json.gz` constitutes the Beacon v2 compatibility deliverable — demonstrating that the pipeline output meets the protocol requirements for federated genomic data sharing.

---

## 14. Synthetic Data Generation

Real sequencing data was not available for pipeline development. Two simulation tools were used:

**wgsim** (bundled with SAMtools) — simulates paired-end short reads from a reference FASTA. Variants are placed randomly at a configurable rate (`-r` for SNPs, `-R` for indels). The random seed (`-S`) controls variant positions, allowing two distinct synthetic individuals to be generated from the same reference.

**badread** (v0.4.1) — simulates long reads with technology-specific error models. The `nanopore2023` error model applies realistic homopolymer errors and quality score distributions derived from real ONT data.

Seed-based reproducibility is system-dependent for simulation tools. The same seed may produce marginally different outputs across operating systems or hardware architectures due to differences in floating-point implementation and random number generation. This is a known limitation of synthetic data validation: results are reproducible within a given environment but may differ slightly across environments. In practice, the differences observed were negligible (0 vs 7 variant calls out of ~64 million base pairs), confirming that the pipeline behaviour is consistent even if the exact simulation output is not.

Parameters used:

| Parameter | Short-read (wgsim) | Long-read (badread) |
|---|---|---|
| Coverage | 1M read pairs (~10x) | 10x |
| Read length | 2 × 150 bp | Variable, mean ~10,000 bp |
| SNP rate | 0.1% | N/A (error model) |
| Indel rate | 0.01% | N/A (error model) |
| Error model | Uniform substitution | nanopore2023 |
| Seed | 42 (HG001), 99 (HG002) | 42 (HG003) |

Synthetic data is appropriate for validating pipeline correctness (processes complete, outputs are well-formed, indexes are staged correctly) but not for biological conclusions. Validation against GIAB truth sets is conducted in Stage 4.

---

## 15. Reproducibility Design Decisions

| Decision | Rationale |
|---|---|
| All tools containerised | Eliminates host environment as a variable |
| Pinned container versions | `deepvariant:1.5.0` not `latest` for production processes |
| `storeDir` for reference indexes | Expensive, deterministic operations cached across runs |
| `mode: 'copy'` in publishDir | Symlinks break after `nextflow clean`; real files persist |
| `task.cpus` not hardcoded | Resource declarations in config, not in process logic |
| Explicit `.bai` and `.fai` staging | Process isolation requires all dependencies declared |
| Read group tags on every BAM | Required for sample identity in all downstream tools |
| `mountFlags = 'z'` | SELinux compatibility without modifying file ownership |
| Normalisation before benchmarking | Canonical variant representation required for coordinate matching |
| Annotation before Beacon ingestion | `molecularEffects` fields require consequence data from VEP |

---

## Validation Status

| Stage | Data | Ts/Tv | Variants called | Status |
|---|---|---|---|---|
| 1 | Synthetic chr20, 1 sample | N/A (0 variants) | 0 | ✓ Pipeline mechanics validated |
| 2 | Synthetic chr20, 2 samples | 0.50 (expected) | ~47,000 per sample | ✓ Multi-sample validated |
| 3 | Synthetic chr20, short + long | 0.50 (short), 0.00 (long) | ~47,000 (short), 0 (long) | ✓ Multi-technology validated |
| 4 | GIAB HG001, GRCh38 chr20 | ~2.0–2.1 (target) | ~68,000 (chr20) | In progress |
