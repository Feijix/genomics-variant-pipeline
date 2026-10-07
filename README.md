# genomics-variant-pipeline

A Nextflow DSL2 variant calling pipeline for human sequencing data. Each sample declares its provenance (platform, chemistry, basecaller, germline or somatic) in the sample sheet. The pipeline uses that declaration to choose the aligner, the variant callers and the caller models. Small variants and structural variants (SVs) are called, filtered, normalised and annotated, then converted to Beacon v2 BFF (`genomicVariations`) records.

The pipeline was developed for the MSc thesis *Bridging EGA and Beacon v2: A Pipeline for Variant-Level Genomic Discovery* (Instituto Superior Técnico, 2026). Every figure in the thesis was produced by the tag `v1.0-thesis`.

## Requirements

| Dependency | Tested version | Notes |
|---|---|---|
| Nextflow | 26.04.6 | needs Java 17 or later |
| Docker | 28.3.3 | every tool runs in a container pinned by version tag |
| Reference FASTA | GRCh38 no-alt analysis set | development used chromosome 20 |
| VEP cache | matching the VEP image | `--vep_cache` |
| Beacon reference data | for bff-tools | `--beacon_refdata`, about 65 GB |

## Running

```bash
nextflow run main.nf -profile docker \
  --samplesheet samplesheet.csv \
  --reference data/reference/chr20.fa \
  --vep_cache /path/to/vep_cache \
  --beacon_refdata /path/to/beacon_refdata \
  --outdir results
```

For whole genomes, add the `wgs` profile (`-profile docker,wgs`). It sets the Mutect2 resources to null, because the default resources cover chromosome 20 only. Supply whole-genome files with `--germline_resource` and `--pon`.

Without `--callers`, each sample runs the default tools of its profile. `--callers deepvariant,clair3` names tools explicitly, and `--callers all` runs every tool the profile lists.

## Sample sheet

| Column | Content |
|---|---|
| `sample` | sample identifier, used in output names |
| `subject` | individual; tumour and normal samples of one subject are paired |
| `status` | `normal` or `tumour` (germline samples are `normal`) |
| `sequencing_type` | `shortread` or `longread`; when empty, inferred from read length (500 bp threshold) |
| `platform` | e.g. `illumina`, `ont`, `pacbio` |
| `chemistry` | e.g. `wgs`, `r10.4.1`, `r9.4.1`, `hifi` |
| `basecaller` | e.g. `dorado_sup`; optional, used when the caller model depends on it |
| `fastq_1`, `fastq_2` | reads; `fastq_2` is empty for single-end and long reads |
| `aligned_bam` | an aligned BAM, used in place of FASTQ to skip alignment |
| `truth_vcf`, `truth_bed` | optional truth set and confident regions; when present, the sample is benchmarked |
| `callers` | optional per-sample tool list, overriding the profile defaults |

## Routing

`conf/routing.config` holds two blocks.

- `profiles`, keyed `platform/chemistry/basecaller`. Each profile lists germline and somatic tools, split into defaults and alternatives. When no profile matches the full key, trailing segments are dropped one at a time and a warning is printed. Alias profiles point to another profile with `falls_back_to`. A sample with no platform is routed to the `shortread` or `longread` family profile from its inferred read length.
- `tools`, a registry of tool properties: `variant_class` (`small` or `sv`), `requires_model` (the profile must name a model), and `requires` (tools that must also be selected, e.g. Strelka2 requires Manta in somatic mode).

A tool whose precondition fails (e.g. `requires_normal` for a tumour-only sample) is withdrawn with a warning. A selection that resolves to zero tools, an unknown tool name, or a missing dependency stops the run with an error. Every sample prints one line, `Routing: <sample> -> profile ..., callers [...]`.

Tools available per data type (`conf/routing.config` is the authoritative list of defaults and alternatives):

| Data | Germline small variants | Germline SVs | Somatic |
|---|---|---|---|
| Illumina short reads | DeepVariant | Manta, DELLY, smoove, GRIDSS | Mutect2, Strelka2 (tumour-normal only) |
| ONT R10.4.1 | DeepVariant (`ONT_R104`), Clair3 | Sniffles2, cuteSV | ClairS |
| ONT R9.4.1 | Clair3 | Sniffles2, cuteSV | ClairS |
| PacBio HiFi | DeepVariant, Clair3 | Sniffles2, cuteSV | ClairS |

## Steps

1. Alignment with bwa-mem2 (short reads) or minimap2 (long reads), then sorting and indexing with samtools. Skipped when `aligned_bam` is given.
2. Variant calling with the selected tools.
3. Small variants: bcftools keeps PASS records, normalises them (`bcftools norm`) and sorts them. VEP annotates the published VCF. bff-tools converts the normalised VCF to BFF and re-annotates it with snpEff.
4. SVs: per-tool VCFs are merged with SURVIVOR, and `scripts/sv_to_bff.py` converts them to BFF.
5. Benchmarking, for samples with a truth set: hap.py (germline small variants), som.py (somatic), Truvari (SVs).

## Layout

| Path | Content |
|---|---|
| `main.nf` | entry point and sample sheet parsing |
| `workflows/` | sub-workflows (variant calling, post-processing) |
| `modules/` | one process per tool |
| `conf/routing.config` | routing profiles and the tool registry |
| `nextflow.config` | parameters, resources, containers, profiles |
| `scripts/` | SV-to-BFF converter and EGA metadata scripts |
| `test/` | regression suite and baselines |

## Mutect2 resources

The defaults in `nextflow.config` are the gnomAD allele-frequency file (`af-only-gnomad`) and the 1000 Genomes panel of normals distributed with GATK, both restricted to chromosome 20. Both are optional. Without them, Mutect2 runs with no germline resource and no panel of normals.

## Regression suite

`test/regression.sh` runs ten configurations (germline short and long reads, somatic tumour-normal and tumour-only, SVs, and routing by inferred read length). It compares each metric with the committed baselines in `test/baselines/` and exits with a non-zero status on any difference. The sequencing data and truth sets are not in the repository. The sample sheets expect them under `data/`, and the thesis lists their public sources.

## EGA metadata

`ega_to_samplesheet.py` builds a sample sheet from an EGA dataset's public metadata. `ega_to_beacon_xlsx.py` maps the same metadata to the Beacon v2 metadata template. From archive metadata alone, three of the six Beacon entities validate (individuals, datasets, cohorts). Biosamples, runs and analyses need fields the archive does not provide.

## Known limitations

- Images are pinned by version tag, not by digest.
- Read-length inference selects only the short- or long-read family, not the chemistry or the caller model.
- Tested on one shared Ubuntu server with the local executor. No cluster executor has been tested.
