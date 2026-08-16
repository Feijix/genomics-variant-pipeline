#!/usr/bin/env bash
#
# test/run_sequential.sh
#
# Runs the regression configurations one at a time, reclaiming the Nextflow
# work directory between each, so the peak disk requirement is one
# configuration rather than eight. Intended for a machine that cannot hold the
# whole set at once — and, on a shared node, one that must not fill the
# filesystem out from under co-hosted services.
#
# Progress is recorded in test/.regression_state, so an interrupted run can be
# resumed without repeating completed configurations.
#
# Usage:
#     bash test/run_sequential.sh              # run all remaining
#     bash test/run_sequential.sh --reset      # forget progress, start over
#     bash test/run_sequential.sh --status     # show what is done
#
# Safety: aborts before starting a configuration if free disk is below
# MIN_FREE_GB, rather than discovering it mid-run.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

STATE="test/.regression_state"
LOGS="test/logs"
RESULTS="test/results"
MIN_FREE_GB=50          # refuse to start a configuration below this

mkdir -p "$LOGS" "$RESULTS"
touch "$STATE"

CONFIGS=(
    germline_hg001_short
    germline_hg002_long
    somatic_rep1
    somatic_rep2
    somatic_tumouronly
    somatic_longread
    sv_strict
    sv_permissive,
    sv_longread
)

case "${1:-}" in
    --reset)  : > "$STATE"; echo "progress cleared"; exit 0 ;;
    --status)
        echo "completed:"
        while read -r line; do [ -n "$line" ] && echo "  $line"; done < "$STATE"
        echo "remaining:"
        for c in "${CONFIGS[@]}"; do grep -qx "$c" "$STATE" || echo "  $c"; done
        exit 0 ;;
esac

free_gb() { df -BG --output=avail . | tail -1 | tr -dc '0-9'; }

reclaim() {
    # Nothing in work/ is reusable: the regression run deliberately avoids
    # -resume, since inheriting cached tasks from an older code version would
    # invalidate the comparison.
    if pgrep -f 'nextflow' > /dev/null; then
        echo "  a nextflow process is running; not touching work/"
        return
    fi
    local before after
    before=$(free_gb)
    rm -rf work .nextflow.log.* 2>/dev/null
    after=$(free_gb)
    echo "  reclaimed $(( after - before )) GB (now ${after} GB free)"
}

echo "starting: $(free_gb) GB free"
echo

for cfg in "${CONFIGS[@]}"; do
    if grep -qx "$cfg" "$STATE"; then
        echo "done   $cfg (skipping)"
        continue
    fi

    avail=$(free_gb)
    if [ "$avail" -lt "$MIN_FREE_GB" ]; then
        echo
        echo "STOP: only ${avail} GB free, need ${MIN_FREE_GB} GB to start '$cfg'."
        echo "Free space and re-run; completed configurations will be skipped."
        exit 2
    fi

    echo "---- $cfg  (${avail} GB free, $(date '+%H:%M:%S'))"
    if bash test/regression.sh "$cfg"; then
        echo "$cfg" >> "$STATE"
        # Benchmark summaries are already published to test/results; the
        # intermediates are not needed by the comparison.
        #
        # Exception: sv_permissive resumes from sv_strict's work directory so
        # that both matching regimes score the SAME callsets. Re-calling would
        # introduce variation that would then be misattributed to the Truvari
        # parameters.
        if [ "$cfg" != "sv_strict" ]; then reclaim; fi
    else
        echo
        echo "FAILED at '$cfg'. Progress is preserved; fix and re-run to continue."
        echo "Log: $LOGS/$cfg.log"
        exit 1
    fi
    echo
done

echo "======================================================================"
echo "All configurations complete. $(free_gb) GB free."
echo
echo "Now compare against the committed baselines:"
echo "  python3 test/compare_baselines.py --new $RESULTS --baseline test/baselines"
echo "======================================================================"
