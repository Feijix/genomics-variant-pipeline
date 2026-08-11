#!/usr/bin/env python3
"""
test/compare_baselines.py

Compare freshly generated benchmark output against the committed baselines in
test/baselines/. Exits non-zero if any metric has drifted beyond tolerance.

Handles the three comparator output formats used by the pipeline:
  hap.py   <prefix>.summary.csv   (germline small variants)
  som.py   <prefix>.stats.csv     (somatic small variants)
  Truvari  summary.json           (structural variants)

Usage:
    python3 test/compare_baselines.py --new data/outputs/benchmark --baseline test/baselines
    python3 test/compare_baselines.py --new ... --baseline ... --tolerance 0.001

Rationale for the default tolerance of 0.0: hap.py and som.py are
deterministic given identical inputs, so any change at all is a regression
worth seeing. Truvari's --pick single can tie-break differently between
releases, so raise the tolerance if the container version has changed.
"""

import argparse
import csv
import json
import sys
from pathlib import Path


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_happy(path):
    """hap.py summary.csv -> {'SNP/PASS': {'recall': x, 'precision': y}, ...}"""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            key = f"{row['Type']}/{row['Filter']}"
            try:
                out[key] = {
                    'recall':    float(row['METRIC.Recall']),
                    'precision': float(row['METRIC.Precision']),
                }
            except (ValueError, KeyError):
                continue          # rows with empty metrics (no calls of that type)
    return out


def parse_sompy(path):
    """som.py stats.csv -> {'SNVs': {...}, 'indels': {...}, 'records': {...}}"""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t = row.get('type')
            if not t:
                continue
            try:
                out[t] = {
                    'recall':    float(row['recall']),
                    'precision': float(row['precision']),
                    'tp':        float(row['tp']),
                    'fp':        float(row['fp']),
                }
            except (ValueError, KeyError):
                continue
    return out


def parse_truvari(path):
    """Truvari summary.json -> {'SV': {...}}"""
    d = json.load(open(path))
    return {'SV': {
        'recall':    float(d['recall']),
        'precision': float(d['precision']),
        'f1':        float(d['f1']),
        'tp':        float(d['TP-base']),
        'fp':        float(d['FP']),
    }}


def load(path):
    p = Path(path)
    if p.name == 'summary.json':
        return parse_truvari(p)
    if p.name.endswith('.stats.csv'):
        return parse_sompy(p)
    if p.name.endswith('.summary.csv'):
        return parse_happy(p)
    raise ValueError(f"unrecognised comparator output: {p}")


# ── Comparison ───────────────────────────────────────────────────────────────

def relpath_key(root, path):
    """Identify a result by its path relative to the benchmark root, so that a
    baseline and a fresh result for the same configuration line up."""
    return str(Path(path).relative_to(root))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--new', required=True, help='fresh benchmark output root')
    ap.add_argument('--baseline', required=True, help='test/baselines root')
    ap.add_argument('--tolerance', type=float, default=0.0,
                    help='absolute tolerance on each metric (default: exact)')
    args = ap.parse_args()

    new_root, base_root = Path(args.new), Path(args.baseline)

    patterns = ['**/*.summary.csv', '**/*.stats.csv', '**/summary.json']
    baselines = {}
    for pat in patterns:
        for f in base_root.glob(pat):
            baselines[relpath_key(base_root, f)] = f

    if not baselines:
        print(f"ERROR: no baseline files found under {base_root}", file=sys.stderr)
        return 2

    failures, missing, checked = [], [], 0

    for key, base_file in sorted(baselines.items()):
        new_file = new_root / key
        if not new_file.exists():
            missing.append(key)
            continue

        try:
            base, new = load(base_file), load(new_file)
        except Exception as e:                       # noqa: BLE001 - report and continue
            failures.append(f"{key}: could not parse ({e})")
            continue

        for group, metrics in base.items():
            if group not in new:
                failures.append(f"{key} [{group}]: absent from new results")
                continue
            for metric, expected in metrics.items():
                actual = new[group].get(metric)
                if actual is None:
                    failures.append(f"{key} [{group}] {metric}: absent")
                    continue
                checked += 1
                if abs(actual - expected) > args.tolerance:
                    failures.append(
                        f"{key} [{group}] {metric}: "
                        f"baseline {expected:.6f} -> new {actual:.6f} "
                        f"(delta {actual - expected:+.6f})"
                    )

    # ── Report ───────────────────────────────────────────────────────────────
    print(f"compared {checked} metrics across {len(baselines) - len(missing)} results")

    if missing:
        print(f"\nNOT RE-RUN ({len(missing)}) — inputs unavailable, or configuration skipped:")
        for k in missing:
            print(f"  - {k}")

    if failures:
        print(f"\nDRIFT DETECTED ({len(failures)}):")
        for f in failures:
            print(f"  ! {f}")
        return 1

    print("\nAll re-run configurations match their baselines.")
    
    if checked == 0:
        print("\nERROR: no metrics were compared. Baseline and result layouts "
              "do not correspond — see NOT RE-RUN above.", file=sys.stderr)
        return 2
    
    return 0


if __name__ == '__main__':
    sys.exit(main())