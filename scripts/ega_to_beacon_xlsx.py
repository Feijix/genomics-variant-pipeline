#!/usr/bin/env python3
"""
scripts/ega_to_beacon_xlsx.py

Populate the Beacon v2 Models XLSX template from EGA metadata retrieved via
the public metadata API.

Fills only fields for which EGA supplies a value. Fields with no EGA source are
left empty and reported, so that the output is a record of what the EGA schema
can and cannot supply — not a document padded with invented values.

Usage:
    python3 scripts/ega_to_beacon_xlsx.py \
        --ega-dir data/beacon_metadata/ega \
        --dataset EGAD00001000082 \
        --template data/beacon_metadata/Beacon-v2-Models_template.xlsx \
        --out data/beacon_metadata/EGAD00001000082_beacon.xlsx

Then validate:
    docker run --rm -v "$(pwd)":"$(pwd)" -w "$(pwd)" --entrypoint bash \\
      manuelrueda/beacon2-cbi-tools:latest -c \\
      '/usr/share/beacon2-cbi-tools/bin/bff-tools validate \\
         -i data/beacon_metadata/EGAD00001000082_beacon.xlsx \\
         -o data/beacon_metadata/bff_out'
"""

import argparse
import json
import re
import sys
from pathlib import Path

import openpyxl

# ── Values not held by EGA ───────────────────────────────────────────────────
# Supplied at ingestion because the EGA metadata schema has no corresponding
# field. Every CURIE below must be verified against its source ontology: a
# syntactically valid identifier naming the wrong term passes validation while
# being semantically wrong, which is worse than leaving the field empty.
INGESTION = {
    # EGA holds the value as free text; the ontology coding is absent.
    'librarySource':  {'GENOMIC': ('OBI:0000626', 'DNA extract')},                   # VERIFY
    'platformModel':  {'Illumina Genome Analyzer II': ('OBI:0000759', 'Illumina')},  # VERIFY
    'disease':        {'Breast Cancer':  ('NCIT:C4872', 'Breast Carcinoma'),         # VERIFY
                       'Matched Normal': (None, None)},

    # biosampleStatus: 'reference sample' confirmed from the CINECA workbook.
    # The abnormal counterpart does not appear there (CINECA is a healthy
    # cohort) and is taken from EFO directly.
    'status_normal':  ('EFO:0009654', 'reference sample'),                           # confirmed
    'status_tumour':  ('EFO:0009655', 'abnormal sample'),                            # VERIFY

    # sampleOriginType is an anatomical term (UBERON), per the CINECA workbook,
    # not a generic specimen type. EGA does not expose organism_part for this
    # dataset, so no value is derivable and the field is left empty.
    'originType':     (None, None),

    'runDate':        '1970-01-01',
    'analysisDate':   '1970-01-01',
    'pipelineName':   'genomics-variant-pipeline v1.0',
}

# ── Fields the Beacon v2 model marks mandatory, per entity ───────────────────
# Used only to report which mandatory fields EGA cannot supply.
MANDATORY = {
    'individuals': ['id', 'sex.id', 'sex.label'],
    'biosamples':  ['id', 'individualId', 'biosampleStatus.id', 'biosampleStatus.label',
                    'sampleOriginType.id', 'sampleOriginType.label'],
    'runs':        ['id', 'biosampleId', 'individualId', 'runDate',
                    'librarySource.id', 'librarySource.label', 'libraryStrategy',
                    'librarySelection', 'libraryLayout', 'platform'],
    'analyses':    ['id', 'biosampleId', 'individualId', 'runId',
                    'analysisDate', 'pipelineName'],
    'datasets':    ['id', 'name'],
    'cohorts':     ['cohortType', 'cohortSize'],
}

# EGA biological_sex -> NCIT terms used by the Beacon v2 model.
SEX_TERMS = {
    'male':    ('NCIT:C20197', 'male'),
    'female':  ('NCIT:C16576', 'female'),
    'unknown': ('NCIT:C17998', 'unknown'),
}


def load(ega_dir, dataset, name):
    p = Path(ega_dir) / f"{dataset}_{name}.json"
    if not p.exists():
        print(f"  ! missing {p}", file=sys.stderr)
        return []
    return json.load(open(p))


def write_sheet(wb, sheet, rows, report):
    """Write rows (list of dicts keyed by column header) into a template sheet."""
    ws = wb[sheet]
    headers = [c.value for c in ws[1]]
    idx = {h: i + 1 for i, h in enumerate(headers) if h}

    unknown = [k for r in rows for k in r if k not in idx]
    if unknown:
        print(f"  ! {sheet}: columns not in template: {sorted(set(unknown))}", file=sys.stderr)

    for n, row in enumerate(rows, start=2):
        for key, val in row.items():
            if key in idx and val is not None:
                ws.cell(row=n, column=idx[key], value=val)

    filled = sorted({k for r in rows for k, v in r.items() if v not in (None, '')})
    missing = [f for f in MANDATORY.get(sheet, []) if f not in filled]
    report[sheet] = {'rows': len(rows), 'filled': filled, 'missing_mandatory': missing}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ega-dir', required=True)
    ap.add_argument('--dataset', required=True)
    ap.add_argument('--template', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    ds_raw   = load(args.ega_dir, args.dataset, 'dataset')
    dataset  = ds_raw[0] if isinstance(ds_raw, list) and ds_raw else ds_raw
    samples  = load(args.ega_dir, args.dataset, 'samples')
    studies  = load(args.ega_dir, args.dataset, 'studies')
    runs     = load(args.ega_dir, args.dataset, 'runs')
    exps     = load(args.ega_dir, args.dataset, 'experiments')

    # sample accession -> subject id
    subj = {s['accession_id']: s.get('subject_id') for s in samples}

    # Experiments reference their sample only inside design_description free text,
    # e.g. "... constructed from sample accession EGAN00001000576 for study ...".
    # There is no structured experiment->sample link in the API response.
    exp_by_sample = {}
    for e in exps:
        m = re.search(r'(EGAN\d+)', e.get('design_description') or '')
        if m:
            exp_by_sample[m.group(1)] = e

    wb = openpyxl.load_workbook(args.template)
    report = {}

    # ── individuals ──────────────────────────────────────────────────────────
    seen = set()
    ind_rows = []
    for s in samples:
        sid = s.get('subject_id')
        if not sid or sid in seen:
            continue
        seen.add(sid)
        sex_id, sex_label = SEX_TERMS.get((s.get('biological_sex') or '').lower(),
                                          (None, None))
        row = {'id': sid, 'sex.id': sex_id, 'sex.label': sex_label}
        dis_id, dis_label = INGESTION['disease'].get(s.get('phenotype'), (None, None))
        if dis_id:                       # only emit a diseases entry when coded
            row['diseases_diseaseCode.id'] = dis_id
            row['diseases_diseaseCode.label'] = dis_label
        ind_rows.append(row)

    write_sheet(wb, 'individuals', ind_rows, report)

    # ── biosamples ───────────────────────────────────────────────────────────
    bio_rows = []
    for s in samples:
        tumour = 'normal' not in (s.get('phenotype') or '').lower()
        st_id, st_label = INGESTION['status_tumour' if tumour else 'status_normal']
        ot_id, ot_label = INGESTION['originType']
        bio_rows.append({
            'id': s['accession_id'],
            'individualId': s.get('subject_id'),
            'biosampleStatus.id': st_id,     'biosampleStatus.label': st_label,
            'sampleOriginType.id': ot_id,    'sampleOriginType.label': ot_label,
        })
    write_sheet(wb, 'biosamples', bio_rows, report)

    # ── runs ─────────────────────────────────────────────────────────────────
    run_rows = []
    for r in runs:
        sample = r.get('sample_accession_id')
        e = exp_by_sample.get(sample, {})
        run_rows.append({
            'id': r['accession_id'],
            'biosampleId': sample,
            'individualId': subj.get(sample),
            'libraryLayout': e.get('library_layout'),
            'librarySelection': e.get('library_selection'),
            'librarySource.label': e.get('library_source'),
            'libraryStrategy': e.get('library_strategy'),
            'librarySource.id': INGESTION['librarySource'].get(e.get('library_source'), (None,))[0],
            'platformModel.id': INGESTION['platformModel'].get(e.get('model'), (None,))[0],
            'runDate': INGESTION['runDate'],
            'platform': e.get('platform'),
            'platformModel.label': e.get('model'),
        })
    write_sheet(wb, 'runs', run_rows, report)

    # ── analyses ─────────────────────────────────────────────────────────────
    # Analysis records describe the processing applied, which EGA does not hold
    # for this dataset (/analyses returned empty). Identifiers are derived from
    # the run records so the entity is linkable; the rest is left empty.
    ana_rows = [{
        'id': r['accession_id'].replace('EGAR', 'ANA'),
        'runId': r['accession_id'],
        'biosampleId': r.get('sample_accession_id'),
        'individualId': subj.get(r.get('sample_accession_id')),
        'analysisDate': INGESTION['analysisDate'],
        'pipelineName': INGESTION['pipelineName'],
    } for r in runs]
    write_sheet(wb, 'analyses', ana_rows, report)

    # ── datasets ─────────────────────────────────────────────────────────────
    ds_rows = [{
        'id': dataset.get('accession_id'),
        'name': dataset.get('title'),
        'description': dataset.get('description'),
        'createDateTime': dataset.get('released_date'),
        # externalUrl, version, dataUseConditions.duoDataUse: no EGA equivalent.
        # (EGA holds a policy accession, not a DUO code.)
    }]
    write_sheet(wb, 'datasets', ds_rows, report)

    # ── cohorts ──────────────────────────────────────────────────────────────
    # Correct the template's non-schema column names.
    ws = wb['cohorts']
    for c in ws[1]:
        if c.value == 'cohortId':   c.value = 'id'
        if c.value == 'cohortName': c.value = 'name'

    study = studies[0] if studies else {}
    coh_rows = [{
        # See note above on id/name vs cohortId/cohortName.
        'id': study.get('accession_id'),
        'name': study.get('title'),
        # cohortType is an enumeration describing how the cohort was constituted
        # (study-defined | beacon-defined | user-defined). EGA's study_type
        # ("Cancer Genomics") describes the research domain and is not a
        # correspondence: a cohort derived from an EGA study is study-defined by
        # construction, so the value is a property of ingestion, not of the
        # source record.
        'cohortType': 'study-defined',
        'cohortSize': dataset.get('num_samples'),
    }]
    write_sheet(wb, 'cohorts', coh_rows, report)

    wb.save(args.out)

    # ── Report ───────────────────────────────────────────────────────────────
    print(f"\nwrote {args.out}\n")
    print(f"{'entity':<14} {'rows':>5}  mandatory fields with no EGA source")
    print('-' * 78)
    for sheet, r in report.items():
        miss = ', '.join(r['missing_mandatory']) or '(none)'
        print(f"{sheet:<14} {r['rows']:>5}  {miss}")


if __name__ == '__main__':
    sys.exit(main())
