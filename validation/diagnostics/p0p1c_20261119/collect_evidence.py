#!/usr/bin/env python3
"""Index of executed runs: inputs, binary, status, steps, cost, hashes, balances.

Usage: collect_evidence.py RESULTS_ROOT OUTPUT_JSON
Every run folder containing run_status.txt is indexed. Balances come from the
independent ledger auditor (scripts/analyze_accepted_conservation.py) using the
accepted fluxes and backward-Euler weights; an initial inventory of zero keeps
its normalization undefined (None) instead of an invented denominator.
External interruptions are not encoded here: a missing termination_reason is
reported as INCOMPLETE, never as a physical failure.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
from analyze_accepted_conservation import audit_ledger  # noqa: E402


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def status(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        if '=' in line:
            k, v = line.split('=', 1)
            out[k.strip()] = v.strip()
    return out


def steps(path: Path) -> dict:
    if not path.is_file():
        return {}
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return dict(accepted=0)
    dts = [float(r['dt_s']) for r in rows]
    return dict(accepted=len(rows), dt_min_s=min(dts), dt_max_s=max(dts),
                dt_mean_s=sum(dts) / len(dts), last_time_s=float(rows[-1]['time_s']),
                rejected_total=int(rows[-1]['rejected_steps_total']))


def main():
    root, output = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    runs = []
    for st in sorted(root.rglob('run_status.txt')):
        folder = st.parent
        s = status(st)
        entry = dict(run=str(folder.relative_to(root)), status=s,
                     termination=s.get('termination_reason', 'INCOMPLETE'), steps=steps(folder / 'steps.csv'))
        deck = folder.parent / (folder.name + '.nml')
        entry['input'] = str(deck.relative_to(root)) if deck.is_file() else None
        entry['input_sha256'] = sha(deck)
        entry['hashes'] = {n: sha(folder / n) for n in ['monitor.csv', 'conservation.csv', 'steps.csv',
                                                     'boundary_conditions.csv', 'fields.csv', 'checkpoint.dat',
                                                     'last_rejected_state.dat', 'run_metadata.txt']
                           if (folder / n).is_file()}
        if (folder / 'conservation.csv').is_file():
            try:
                ledger = audit_ledger(folder / 'conservation.csv')
                entry['ledger_arithmetic_passes'] = ledger['arithmetic_contract_passes']
                entry['junction_integral_mismatch_kg'] = ledger['max_junction_integral_mismatch_kg']
                entry['balances'] = ledger['balances']
                rel = [v['relative_inventory'] for v in ledger['balances'].values() if v['relative_inventory'] is not None]
                entry['max_balance_over_M0'] = max(rel) if rel else None
            except Exception as exc:  # recorded, never silently dropped
                entry['ledger_error'] = f'{type(exc).__name__}: {exc}'
        runs.append(entry)
    binaries = {p.name: sha(p) for p in sorted((root / 'bin').glob('*'))} if (root / 'bin').is_dir() else {}
    output.write_text(json.dumps(dict(results_root=str(root), binaries=binaries, runs=runs), indent=1) + '\n')
    for r in runs:
        print(r['run'], r['termination'], r['status'].get('reached_time_s'), r['steps'].get('accepted'),
              r['steps'].get('rejected_total'), r.get('max_balance_over_M0'))


if __name__ == '__main__':
    main()
