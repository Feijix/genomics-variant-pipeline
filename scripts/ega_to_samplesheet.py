#!/usr/bin/env python3
"""
scripts/ega_to_samplesheet.py

Build a pipeline sample sheet from the metadata an EGA dataset exposes through its
public API (the JSON files fetched for ega_to_beacon_xlsx.py), so that routing can be
resolved from archive metadata alone. One row per sequencing run.

  platform, chemistry  <- experiment.platform, experiment.library_strategy
  status               <- sample.phenotype ('Matched Normal' -> normal, else tumor)
  subject              <- stem of sample.subject_id (PD3945a / PD3945b -> PD3945)
  aligned_bam          <- placeholder path; the reads are access-controlled

Usage:
  python3 scripts/ega_to_samplesheet.py --ega-dir data/beacon_metadata/ega \
      --dataset EGAD00001000082 --out data/samplesheet.ega.EGAD00001000082.csv
"""
import argparse, csv, json, re
from pathlib import Path

def load(d, ds, name):
    p = Path(d) / f"{ds}_{name}.json"
    return json.load(open(p)) if p.exists() else []

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ega-dir', required=True)
    ap.add_argument('--dataset', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    samples = {s['accession_id']: s for s in load(a.ega_dir, a.dataset, 'samples')}
    runs    = load(a.ega_dir, a.dataset, 'runs')
    exps    = load(a.ega_dir, a.dataset, 'experiments')

    # experiments reference their sample only in free text
    exp_by_sample = {}
    for e in exps:
        m = re.search(r'(EGAN\d+)', e.get('design_description') or '')
        if m:
            exp_by_sample[m.group(1)] = e

    cols = ['sample','subject','status','sequencing_type','platform','chemistry','basecaller',
            'fastq_1','fastq_2','aligned_bam','truth_vcf','truth_bed','callers']
    n_t = n_n = 0
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh); w.writerow(cols)
        for r in runs:
            s = samples.get(r.get('sample_accession_id'), {})
            e = exp_by_sample.get(r.get('sample_accession_id'), {})
            subject_id = s.get('subject_id') or ''
            m = re.match(r'(PD\d+)', subject_id)           # Sanger naming: PDnnnn + a (tumour) / b (normal)
            subject = m.group(1) if m else subject_id
            status = 'normal' if 'normal' in (s.get('phenotype') or '').lower() else 'tumor'
            n_n += status == 'normal'; n_t += status == 'tumor'
            platform = (e.get('platform') or '').lower()
            chemistry = (e.get('library_strategy') or '').lower()   # 'wgs' for Illumina
            seqtype = 'shortread' if platform == 'illumina' else ''
            w.writerow([r['accession_id'], subject, status, seqtype, platform, chemistry, '',
                        '', '', f"data/ega/{r['accession_id']}.bam", '', '', ''])
    print(f"wrote {a.out}: {len(runs)} runs, {n_t} tumour, {n_n} normal")

if __name__ == '__main__':
    main()
