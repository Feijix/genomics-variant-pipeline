#!/usr/bin/env python3
# parse_trace.py — Nextflow trace.txt -> normalized TSV for §5.3
# Usage: python3 parse_trace.py trace.txt > metrics.tsv
import sys, csv, re

def to_seconds(s):
    # "1h 2m 3s" / "2m 3s" / "3.4s" / "3.4ms" -> float seconds
    if not s or s in ('-','0'): return 0.0
    if s.endswith('ms') and ' ' not in s: return float(s[:-2]) / 1000
    tot = 0.0
    for val, unit in re.findall(r'([\d.]+)(h|m|s|ms)', s):
        tot += float(val) * {'h':3600,'m':60,'s':1,'ms':0.001}[unit]
    return tot

def to_bytes(s):
    # "3.2 GB" / "512 MB" / "0" -> float bytes (Nextflow uses power-of-2)
    if not s or s in ('-','0'): return 0.0
    m = re.match(r'([\d.]+)\s*([KMGT]?B)', s.strip())
    if not m: return 0.0
    v, u = float(m.group(1)), m.group(2)
    return v * {'B':1,'KB':1024,'MB':1024**2,'GB':1024**3,'TB':1024**4}[u]

def to_pct(s):
    if not s or s in ('-','0'): return 0.0
    return float(s.replace('%',''))

rows = list(csv.DictReader(open(sys.argv[1]), delimiter='\t'))
w = csv.writer(sys.stdout, delimiter='\t')
w.writerow(['process','tag','realtime_s','cpu_pct','cpus',
            'cpu_efficiency','peak_rss_gb','peak_vmem_gb',
            'rchar_gb','wchar_gb','status'])

for r in rows:
    cpus = float(r.get('cpus','1') or 1)
    cpu_pct = to_pct(r.get('%cpu',''))
    eff = round(cpu_pct / (cpus * 100), 3) if cpus else 0.0  # 1.0 = full utilisation
    w.writerow([
        r.get('process',''), r.get('tag',''),
        round(to_seconds(r.get('realtime','')), 1),
        cpu_pct, int(cpus), eff,
        round(to_bytes(r.get('peak_rss','')) / 1024**3, 3),
        round(to_bytes(r.get('peak_vmem','')) / 1024**3, 3),
        round(to_bytes(r.get('rchar','')) / 1024**3, 3),
        round(to_bytes(r.get('wchar','')) / 1024**3, 3),
        r.get('status',''),
    ])