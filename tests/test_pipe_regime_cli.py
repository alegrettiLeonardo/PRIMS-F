#!/usr/bin/env python3
"""Standalone CLI smoke: explicit pipe mode, horizontal pipe and domain refusal.

Build the executable first, then run:
  python3 tests/test_pipe_regime_cli.py --binary bin/prims_f
This test launches only tiny four-cell-per-segment, 0.002 s cases.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess


PHYSICAL_MONITOR_COLUMNS = [
    'time_s', 'dt_s', 'p_in_Pa', 'p_riser_base_Pa', 'alpha_riser_base',
    'mdot_l_out_kg_s', 'mdot_g_out_kg_s', 'Ml_kg', 'Mg_kg',
    'Ml_in_integral_kg', 'Mg_in_integral_kg',
    'Ml_out_integral_kg', 'Mg_out_integral_kg',
]


def make_deck(base, values):
    result = base
    for key, value in values.items():
        pattern = rf'(?i)\b{re.escape(key)}\s*=\s*(?:\x27[^\x27]*\x27|[^,\n/]+)'
        replacement = f'{key}={value}'
        if re.search(pattern, result):
            result = re.sub(pattern, lambda _: replacement, result)
        else:
            result = result.replace('\n/', '\n  ' + replacement + ',\n/')
    return result


def read_run_status(path):
    if not path.exists():
        return None
    return dict(
        line.split('=', 1)
        for line in path.read_text().splitlines()
        if '=' in line
    )


def case_passed(returncode, status, output, expected, message):
    if expected == 'PASS':
        return (
            returncode == 0
            and status is not None
            and status.get('status', '').strip() == '0'
        )
    passed = returncode != 0 and message in output
    if expected == 'REJECT_DOMAIN':
        domain_status = re.search(
            r'Steady initialization failed, status=\s*61\b', output
        )
        passed = passed and bool(domain_status)
    return passed


def read_physical_monitor(path):
    with path.open() as stream:
        return [
            {key: row[key] for key in PHYSICAL_MONITOR_COLUMNS}
            for row in csv.DictReader(stream)
        ]


def legacy_output_matches_fields(directory):
    with (directory / 'fields.csv').open() as stream:
        fields = list(csv.DictReader(stream))
    fields = [row for row in fields if row['time_s'] == fields[-1]['time_s']]
    # The ee5810e writer can repeat its terminal frame at the same time.
    # Retain the final record at each physical cell position.
    fields = list({float(row['s_m']): row for row in fields}.values())
    legacy = [
        list(map(float, line.split()))
        for line in (directory / 'results_D_Stable.dat').read_text().splitlines()
        if not line.startswith('#')
    ]
    keys = ['s_m', 'alpha', 'p_Pa', 'ul_m_s', 'ug_m_s']
    if len(fields) != len(legacy):
        return False
    return all(
        math.isclose(row[column + 1], float(field[key]), rel_tol=1e-8, abs_tol=1e-10)
        for row, field in zip(legacy, fields)
        for column, key in enumerate(keys)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--binary', type=Path, default=root / 'bin/prims_f')
    parser.add_argument('--output', type=Path, default=root / 'tmp/pipe_regime_cli')
    args = parser.parse_args()
    binary = args.binary.resolve()
    out = args.output.resolve()
    if not binary.is_file():
        raise SystemExit(f'Build the solver first: {binary}')
    out.mkdir(parents=True, exist_ok=True)
    base = (root / 'input/benchmarks/D_Stable.nml').read_text()
    common = dict(
        nodes=4, n_riser=4, time_dim=.002, dt_start_s=.001, dt_max_s=.001,
        max_iter=100, output_interval_s=.001, field_interval_s=.001,
        perturbation_amplitude=0,
    )
    cases = [
        ('default', {}, 'PASS', None),
        ('explicit_stratified',
         dict(pipe_regime_mode="'FIXED_STRATIFIED'", closure_pipe=1), 'PASS', None),
        ('auto_diagnostic', dict(pipe_regime_mode="'AUTO_BURR15'"), 'PASS', None),
        ('drift_horizontal',
         dict(pipe_regime_mode="'FIXED_DRIFT'", closure_pipe=2, beta_deg=0), 'PASS', None),
        ('drift_downhill',
         dict(pipe_regime_mode="'FIXED_DRIFT'", closure_pipe=2), 'PASS', None),
        ('drift_missing_id', dict(pipe_regime_mode="'FIXED_DRIFT'"),
         'REJECT_CONFIG', 'requires explicit closure_pipe'),
        ('drift_id1', dict(pipe_regime_mode="'FIXED_DRIFT'", closure_pipe=1),
         'REJECT_CONFIG', 'requires explicit closure_pipe'),
        ('drift_id8', dict(pipe_regime_mode="'FIXED_DRIFT'", closure_pipe=8),
         'REJECT_CONFIG', 'requires explicit closure_pipe'),
        ('stratified_drift_id', dict(closure_pipe=2),
         'REJECT_CONFIG', 'require closure_pipe=1'),
        ('auto_drift_id', dict(pipe_regime_mode="'AUTO_BURR15'", closure_pipe=5),
         'REJECT_CONFIG', 'require closure_pipe=1'),
        ('unknown_mode', dict(pipe_regime_mode="'BAD_MODE'"),
         'REJECT_CONFIG', 'Unknown pipe_regime_mode'),
        ('domain_horizontal_countercurrent',
         dict(pipe_regime_mode="'FIXED_DRIFT'", closure_pipe=5,
              closure_riser=2, beta_deg=0, jl=-.1, mul=-1),
         'REJECT_DOMAIN', 'Steady initialization failed'),
    ]
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    results = []
    for name, extra, expected, message in cases:
        dest = out / name
        values = dict(common, out_dir="'" + str(dest) + "'", **extra)
        path = out / (name + '.nml')
        path.write_text(make_deck(base, values))
        process = subprocess.run(
            [str(binary), str(path)], cwd=root, env=env,
            capture_output=True, text=True, timeout=60,
        )
        output = process.stdout + process.stderr
        (out / (name + '.log')).write_text(output)
        status = read_run_status(dest / 'run_status.txt')
        passed = case_passed(process.returncode, status, output, expected, message)
        results.append(dict(
            case=name, expected=expected, returncode=process.returncode,
            passed=passed, status=status,
        ))
    if all(result['passed'] for result in results[:3]):
        default = read_physical_monitor(out / 'default/monitor.csv')
        stratified = read_physical_monitor(out / 'explicit_stratified/monitor.csv')
        auto = read_physical_monitor(out / 'auto_diagnostic/monitor.csv')
        results.append(dict(
            case='default_stratified_auto_physical_identity', expected='PASS',
            passed=default == stratified == auto,
        ))
    for result in results[:5]:
        if not result['passed']:
            continue
        name = result['case']
        results.append(dict(
            case=name + '_legacy_output_matches_operator_fields', expected='PASS',
            passed=legacy_output_matches_fields(out / name),
        ))
    (out / 'summary.json').write_text(
        json.dumps(dict(binary=str(binary), cases=results), indent=2) + '\n'
    )
    for result in results:
        label = 'PASS' if result['passed'] else 'FAIL'
        print(label, result['case'], result['expected'])
    if all(result['passed'] for result in results):
        return 0
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
