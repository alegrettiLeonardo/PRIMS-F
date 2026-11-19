#!/usr/bin/env python3
"""Short solver checks for accepted ledgers, sparse output and early termination."""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from test_pipe_regime_cli import make_deck, read_run_status
spec = importlib.util.spec_from_file_location('accepted_mass', ROOT / 'scripts/analyze_accepted_conservation.py')
audit = importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)


def terminal_fields(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    return [r for r in rows if r['time_s'] == rows[-1]['time_s']]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): parser.error('Choose a fresh output directory')
    out = args.output.resolve(); out.mkdir(parents=True)
    base = (ROOT / 'input/benchmarks/D_Stable.nml').read_text()
    common = dict(nodes=4, n_riser=4, time_dim=.02, dt_start_s=.001, dt_max_s=.004,
                  max_iter=100, output_interval_s=1., field_interval_s=1., perturbation_amplitude=.01)
    cases = [('sparse', {}), ('frequent', dict(output_interval_s=.001, field_interval_s=.001)),
             ('limited', dict(max_iter=2)),
             ('rejected', dict(newton_max=1, tol_res=1e-30, dt_min_s=.0004))]
    results = []
    for name, extra in cases:
        folder = out / name
        deck = out / (name + '.nml')
        deck.write_text(make_deck(base, dict(common, **extra) | dict(out_dir="'" + str(folder) + "'")))
        process = subprocess.run([str(args.binary.resolve()), str(deck)], cwd=ROOT,
            env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'),
            capture_output=True, text=True, timeout=60)
        (out / (name + '.log')).write_text(process.stdout + process.stderr)
        status = {k: float(v) for k, v in read_run_status(folder / 'run_status.txt').items()}
        ledger = audit.audit_ledger(folder / 'conservation.csv')
        steps = audit.audit_steps(folder / 'steps.csv', ledger)
        expected = 0 if name in ('sparse', 'frequent') else 7000 if name == 'limited' else None
        passed = ledger['arithmetic_contract_passes'] and ledger['accepted_steps'] == status['accepted_steps']
        passed = passed and ledger['final_time_s'] == status['reached_time_s']
        if expected is not None:
            passed = passed and status['status'] == expected and ((process.returncode == 0) == (expected == 0))
        else:
            passed = passed and status['rejected_steps'] > 0 and process.returncode != 0
        if name == 'sparse':
            with (folder / 'monitor.csv').open() as stream: monitor_rows = sum(1 for _ in csv.DictReader(stream))
            passed = passed and ledger['accepted_steps'] > monitor_rows
        results.append(dict(case=name, passed=passed, returncode=process.returncode, status=status, ledger=ledger, steps=steps))
    results.append(dict(case='sampling_does_not_change_terminal_state',
                        passed=terminal_fields(out / 'sparse/fields.csv') == terminal_fields(out / 'frequent/fields.csv')))
    (out / 'summary.json').write_text(json.dumps(dict(binary=str(args.binary.resolve()), cases=results), indent=2) + '\n')
    for result in results: print('PASS' if result['passed'] else 'FAIL', result['case'])
    return 0 if all(r['passed'] for r in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
