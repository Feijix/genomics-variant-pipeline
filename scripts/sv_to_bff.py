#!/usr/bin/env python3
"""
sv_to_bff.py

Convert a typed structural-variant VCF (symbolic <DEL>/<INS>/<DUP>/<INV> with
SVTYPE/SVLEN/END) into Beacon v2 genomicVariations BFF as gzipped NDJSON
(one JSON document per line), schema-valid against the Beacon v2 Models
genomicVariations defaultSchema (validated with bff-validator -gv -gv-vcf).

Schema shape (branch 2 of variation.oneOf, the Beacon-legacy branch):
  variation requires variantType + alternateBases + location.
  - For SVs the change is not a base substitution, so alternateBases is the
    EMPTY string "" (the schema explicitly uses empty alternateBases for
    deletions; symbolic <DEL> would fail alternateBases' IUPAC pattern).
  - variantType carries the SV type (DEL/INS/DUP/INV).
  - location is a VRS SequenceLocation: type=SequenceLocation, a CURIE
    sequence_id (RefSeq accession), and a SequenceInterval with Number
    start/end (0-based, interbase).

Coordinates: VCF POS is 1-based; VRS SequenceInterval is 0-based interbase.
POS -> start = POS-1 ; END (1-based inclusive) -> end = END.
Insertions (no reference span) -> start == end.

Usage:
  sv_to_bff.py <input.vcf[.gz]> <output.json.gz> \
      --biosample <id> --individual <id> --run <id> --analysis <id> \
      --dataset <id> --assembly GRCh38 --caller <name>
"""

import sys
import gzip
import json
import argparse
import re

# GRCh38 chromosome -> RefSeq accession (CURIE value for sequence_id).
REFSEQ_GRCH38 = {
    "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12",
    "chr4": "NC_000004.12", "chr5": "NC_000005.10", "chr6": "NC_000006.12",
    "chr7": "NC_000007.14", "chr8": "NC_000008.11", "chr9": "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
    "chrM": "NC_012920.1",
}

SV_LABELS = {"DEL": "deletion", "INS": "insertion", "DUP": "duplication", "INV": "inversion"}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_info(info):
    d = {}
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def refseq_curie(chrom):
    acc = REFSEQ_GRCH38.get(chrom)
    if acc is None:
        # normalise a bare '20' style name too
        acc = REFSEQ_GRCH38.get("chr" + chrom)
    return "refseq:" + acc if acc else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("out")
    ap.add_argument("--biosample", required=True)
    ap.add_argument("--individual", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--assembly", default="GRCh38")
    ap.add_argument("--caller", default="unknown")
    args = ap.parse_args()

    n = 0
    skipped = 0
    with opener(args.vcf) as fh, gzip.open(args.out, "wt") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = c[:8]
            ii = parse_info(info)

            svtype = ii.get("SVTYPE")
            if svtype is True:
                svtype = None
            if svtype is None:
                m = re.match(r"<(\w+)>", alt)
                svtype = m.group(1) if m else None
            if svtype not in ("DEL", "INS", "DUP", "INV"):
                skipped += 1
                continue

            seqid = refseq_curie(chrom)
            if seqid is None:
                skipped += 1
                continue

            pos = int(pos)
            start0 = pos - 1
            if "END" in ii and ii["END"] is not True:
                end0 = int(ii["END"])
            elif "SVLEN" in ii and ii["SVLEN"] is not True:
                end0 = start0 + abs(int(ii["SVLEN"]))
            else:
                end0 = start0
            if end0 < start0:
                end0 = start0

            internal_id = "{}-{}-{}-{}-{}".format(args.assembly, chrom, pos, end0, svtype)
            svlen = ii.get("SVLEN")
            svlen = None if svlen is True else svlen

            gt = c[9].split(":")[0] if len(c) > 9 else "./."
            zyg = ({"id": "GENO:0000136", "label": "homozygous"}
                   if gt in ("1/1", "1|1")
                   else {"id": "GENO:0000135", "label": "heterozygous"})

            record = {
                "variantInternalId": internal_id,
                "variation": {
                    "variantType": svtype,
                    "alternateBases": "",       # required; empty for symbolic SVs
                    "referenceBases": "",
                    "location": {
                        "type": "SequenceLocation",
                        "sequence_id": seqid,   # CURIE, e.g. refseq:NC_000020.11
                        "interval": {
                            "type": "SequenceInterval",
                            "start": {"type": "Number", "value": start0},
                            "end":   {"type": "Number", "value": end0},
                        },
                    },
                },
                "identifiers": {
                    "genomicHGVSId": "{}:g.{}_{}{}".format(chrom, pos, end0, svtype.lower()),
                },
                "molecularAttributes": {"geneIds": [], "molecularEffects": []},
                "caseLevelData": [
                    {
                        "biosampleId": args.biosample,
                        "individualId": args.individual,
                        "runId": args.run,
                        "analysisId": args.analysis,
                        "zygosity": zyg,
                    }
                ],
                "info": {
                    "caller": args.caller,
                    "svType": svtype,
                    "svTypeLabel": SV_LABELS.get(svtype),
                    "svLen": svlen,
                    "filter": filt,
                    "datasetId": args.dataset,
                    "assemblyId": args.assembly,
                },
            }
            out.write(json.dumps(record) + "\n")
            n += 1

    sys.stderr.write("Wrote {} SV genomicVariations records to {} (skipped {})\n".format(n, args.out, skipped))


if __name__ == "__main__":
    main()
