#!/usr/bin/env bash
#
# test/regression.sh
#
# Re-run every validated pipeline configuration against the current code, so
# that all figures reported in the dissertation can be attributed to one tagged
# version rather than to the several versions that produced them originally.
#
# Each configuration:
#   - checks its own inputs and SKIPS (does not fail) if any are absent
#   - runs WITHOUT -resume: inheriting cached tasks from an older code version
#     would defeat the purpose of the exercise
#   - writes to its own --outdir, so results cannot overwrite one another
#   - suppresses BAM publication; only benchmark summaries are compared
#
# Usage:
#     bash test/regression.sh                     # all configurations
#     bash test/regression.sh somatic_rep1        # one configuration
#     bash test/regression.sh --list              # print names and exit
#
# Afterwards:
#     python3 test/compare_baselines.py --new test/results --baseline test/baselines
#
# On a disk-constrained machine prefer test/run_sequential.sh, which runs one
# configuration at a time and reclaims intermediates between them.
#
# NOTE ON REFERENCES: every configuration states its reference explicitly.
# The pipeline carries three GRCh38 variants because three data sources were
# aligned against three different flavours of the assembly. The correct
# reference is a property of the input data, not of the pipeline, and choosing
# wrongly does not always fail loudly.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

RESULTS="test/results"
LOGS="test/logs"
mkdir -p "$RESULTS" "$LOGS"

ARG="${1:-}"
PASSED=(); FAILED=(); SKIPPED=()

CONFIG_NAMES=(
    germline_hg001_short
    germline_hg002_long
    somatic_rep1
    somatic_rep2
    somatic_tumouronly
    somatic_longread
    sv_strict
    sv_permissive
    sv_longread
)

if [ "$ARG" = "--list" ]; then
    printf '%s\n' "${CONFIG_NAMES[@]}"
    exit 0
fi

have() {
    for f in "$@"; do [ -e "$f" ] || return 1; done
    return 0
}

run_config() {
    local name="$1"; shift
    local -a inputs=(); local -a nfargs=(); local collecting=1
    for arg in "$@"; do
        if [ "$arg" = "--" ]; then collecting=0; continue; fi
        if [ "$collecting" -eq 1 ]; then inputs+=("$arg"); else nfargs+=("$arg"); fi
    done

    if [ -n "$ARG" ] && [ "$ARG" != "$name" ]; then return 0; fi

    if ! have "${inputs[@]}"; then
        echo "SKIP  $name — missing input(s):"
        for f in "${inputs[@]}"; do [ -e "$f" ] || echo "        $f"; done
        SKIPPED+=("$name"); return 0
    fi

    echo "RUN   $name  ($(date '+%H:%M:%S'))"
    if nextflow run main.nf -profile docker \
         --outdir "$RESULTS/$name" \
         --publish_bams false \
         -with-trace "$LOGS/trace_$name.txt" \
         "${nfargs[@]}" > "$LOGS/$name.log" 2>&1; then
        echo "  ok  $name  ($(date '+%H:%M:%S'))"
        PASSED+=("$name")
    else
        echo "  FAIL $name — see $LOGS/$name.log"
        tail -25 "$LOGS/$name.log" | sed 's/^/        /'
        FAILED+=("$name")
    fi
}

# ── Configurations ───────────────────────────────────────────────────────────
# run_config NAME <required input paths...> -- <nextflow args...>

# 1. Germline, short read, HG001 chr20 (FASTQ input — exercises alignment).
#    No --skip_bwa_index: the index is required here, and chr20 builds in minutes.
run_config germline_hg001_short \
    data/samplesheet.regression.germline_short.csv \
    data/giab/HG001_real_R1.fastq.gz data/giab/HG001_real_R2.fastq.gz \
    data/reference/GRCh38.chr20.fa \
    -- \
    --samplesheet data/samplesheet.regression.germline_short.csv \
    --callers deepvariant \
    --reference data/reference/GRCh38.chr20.fa \
    --reference_germline data/reference/GRCh38.chr20.fa

# 2. Germline, long read, HG002 — Clair3 across ONT R10.4.1, ONT R9.4.1, HiFi.
run_config germline_hg002_long \
    data/samplesheet.regression.germline_long.csv \
    data/samples/HG002_ONT_R10_chr20.fastq.gz \
    data/samples/HG002_HiFi_chr20.bam \
    data/reference/GRCh38.chr20.fa \
    -- \
    --samplesheet data/samplesheet.regression.germline_long.csv \
    --callers clair3 \
    --reference data/reference/GRCh38.chr20.fa \
    --reference_germline data/reference/GRCh38.chr20.fa \
    --skip_bwa_index true

# 3. Somatic, short read, paired, replicate 1.
run_config somatic_rep1 \
    data/samplesheet.seqc2.csv \
    data/seqc2/HCC1395_chr20.bam data/seqc2/HCC1395BL_chr20.bam \
    data/reference/d1vd1/GRCh38.d1.vd1.fa \
    -- \
    --samplesheet data/samplesheet.seqc2.csv \
    --callers mutect2,strelka2,manta_somatic \
    --reference data/reference/d1vd1/GRCh38.d1.vd1.fa \
    --skip_bwa_index true

# 4. Somatic, short read, paired, replicate 2 (~45% shallower than replicate 1).
run_config somatic_rep2 \
    data/samplesheet.seqc2.rep2.csv \
    data/seqc2/HCC1395_rep2_chr20.bam data/seqc2/HCC1395BL_rep2_chr20.bam \
    -- \
    --samplesheet data/samplesheet.seqc2.rep2.csv \
    --callers mutect2,strelka2,manta_somatic \
    --reference data/reference/d1vd1/GRCh38.d1.vd1.fa \
    --skip_bwa_index true

# 5. Somatic, tumour-only. Strelka2 and Manta-somatic are listed deliberately:
#    requires_normal should exclude both, and this run confirms it still does.
run_config somatic_tumouronly \
    data/samplesheet.seqc2.tumoronly.csv \
    data/seqc2/HCC1395_chr20.bam \
    -- \
    --samplesheet data/samplesheet.seqc2.tumoronly.csv \
    --callers mutect2,strelka2,manta_somatic \
    --reference data/reference/d1vd1/GRCh38.d1.vd1.fa \
    --skip_bwa_index true

# 6. Somatic, long read — ClairS on ONT HCC1395/BL.
run_config somatic_longread \
    data/samplesheet.hcc1395_ont.csv \
    data/ont/HCC1395_ONT_chr20.bam data/ont/HCC1395BL_ONT_chr20.bam \
    data/reference/GRCh38.fa \
    -- \
    --samplesheet data/samplesheet.hcc1395_ont.csv \
    --callers clairs,deepvariant \
    --reference data/reference/GRCh38.fa \
    --reference_germline data/reference/GRCh38.fa \
    --skip_bwa_index true

# 7. Structural variants — strict matching (the primary reported figures).
run_config sv_strict \
    data/samplesheet.hg002_sv.csv \
    data/giab/HG002_SV/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz \
    data/giab/HG002_SV/HG002_SV_chr20.benchmark.bed \
    -- \
    --samplesheet data/samplesheet.hg002_sv.csv \
    --callers manta,delly,smoove \
    --reference data/reference/GRCh38.chr20.fa \
    --truvari_pctseq 0.7 --truvari_dup_to_ins false

# 8. Structural variants — permissive matching (sensitivity analysis).
run_config sv_permissive \
    data/samplesheet.hg002_sv.csv \
    -- \
    -resume \
    --samplesheet data/samplesheet.hg002_sv.csv \
    --callers manta,delly,smoove \
    --reference data/reference/GRCh38.chr20.fa \
    --truvari_pctseq 0 --truvari_dup_to_ins true

# 9. Structural variants, long read — Sniffles2 and cuteSV on HG002.
run_config sv_longread \
    data/samplesheet.regression.sv_longread.csv \
    data/samples/HG002_HiFi_chr20.bam \
    data/samples/HG002_ONT_R10_chr20.fastq.gz \
    data/giab/HG002_SV/GRCh38_HG2-T2TQ100-V1.1_stvar.vcf.gz \
    data/giab/HG002_SV/HG002_SV_chr20.benchmark.bed \
    -- \
    --samplesheet data/samplesheet.regression.sv_longread.csv \
    --callers sniffles2,cutesv \
    --reference data/reference/GRCh38.chr20.fa \
    --skip_bwa_index true \
    --truvari_pctseq 0.7 --truvari_dup_to_ins false
# ── Summary ──────────────────────────────────────────────────────────────────
echo
echo "======================================================================"
echo "passed:  ${#PASSED[@]}   ${PASSED[*]:-}"
echo "failed:  ${#FAILED[@]}   ${FAILED[*]:-}"
echo "skipped: ${#SKIPPED[@]}  ${SKIPPED[*]:-}"
echo "======================================================================"
echo
echo "Next: python3 test/compare_baselines.py --new $RESULTS --baseline test/baselines"

[ "${#FAILED[@]}" -eq 0 ] || exit 1