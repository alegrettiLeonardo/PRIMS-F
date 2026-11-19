"""Exact-bit terminal probe plus independent 80-digit cocurrent bound."""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import struct
import subprocess
import sys

import mpmath as mp


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--buildmanifest', type=Path, required=True)
    parser.add_argument('--smoke-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    source = args.source_root.resolve()
    build = json.loads(args.buildmanifest.read_text())
    smoke = json.loads(args.smoke_manifest.read_text())
    if smoke['case'] != 'D_SS3' or smoke['state'] != 'FINISHED' or smoke.get('timed_out'):
        raise SystemExit('This replay requires a true completed SS3 terminal, not a timeout snapshot.')
    if smoke['source_commit'] != build['git_head']:
        raise SystemExit('Build and terminal source commits differ.')
    rejection = smoke.get('last_rejected_iteration') or {}
    if rejection.get('proposal_stage') != 'CELL':
        raise SystemExit('Terminal proposal is not a cell rejection; use a junction/boundary diagnostic.')
    focus = int(rejection['proposal_index'])
    objects = Path(build['output']) / 'obj'
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fixture = out / 'fixture'
    fixture.mkdir()
    output = Path(smoke['output_directory'])
    snapshot = output / 'last_rejected_state.dat'
    shutil.copy2(snapshot, fixture / snapshot.name)
    original_parameters = source / 'validation/diagnostics/qualification2_20261119/globalization_flux/fixtures/D_SS3_active_terminal/probe_parameters.dat'
    shutil.copy2(original_parameters, fixture / 'probe_parameters.dat')
    api = source / 'validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90'
    probe = here / 'replay_frozen_terminal.f90'
    chosen = sorted(path for path in objects.glob('*.o') if path.name != 'main.o')
    if not chosen:
        raise SystemExit('No frozen objects available.')
    command = ['gfortran', '-std=f2018', '-O2', '-fopenmp', '-ffree-line-length-none',
               '-I' + str(objects), '-J' + str(out), str(api), str(probe),
               *map(str, chosen), '-llapack', '-lblas', '-o', str(out / 'probe')]
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    compiled = subprocess.run(command, cwd=out, env=environment, text=True, capture_output=True)
    (out / 'compile.log').write_text(compiled.stdout + compiled.stderr)
    compiled.check_returncode()
    replay = subprocess.run([str(out / 'probe'), str(fixture), str(focus)], cwd=out,
                            env=environment, text=True, capture_output=True)
    (out / 'probe.log').write_text(replay.stdout + replay.stderr)
    replay.check_returncode()

    sys.path.insert(0, str(source / 'scripts'))
    reference = load_module('frozen_bound', source / 'scripts/replay_root_bounds.py')
    mp.mp.dps = 80
    controls = {'cl': 1498, 'cg_fixed': 343, 'p_l0': 101325, 'rho_l0': '997.05',
                'd': '.05254', 'eps': '.000046', 'mu_l': '.00089', 'sigma': '.0728'}
    records = []
    with (out / 'proposal_cells_exact.csv').open() as stream:
        candidates = list(csv.DictReader(stream))
    for row in candidates:
        # Whole ray is retained in CSV; high precision targets endpoints and half-step.
        if int(row['k']) not in (-2, -1, 1, 2, 18):
            continue
        for representation in ('serialized_decimal', 'exact_binary64'):
            values = dict(row)
            if representation == 'exact_binary64':
                for key in ('theta', 'U1', 'U2', 'U3'):
                    value = struct.unpack('>d', bytes.fromhex(row[key + '_bits']))[0]
                    values[key] = mp.mpf(value)
            result = reference.bound(values, controls, '.48')
            # Same conservative scale used by the existing sufficient-bound check,
            # converted to SI. This is a numerical separation criterion, not an
            # interval-arithmetic certification of every floating-point operation.
            _, drift, _ = reference.coefficients((
                result['alpha'], result['rho_l_kg_m3'], result['rho_g_kg_m3'],
                0, result['u_g_max_m_s'], values['theta'], controls['d'],
                controls['eps'], controls['mu_l'], controls['sigma']))
            scale = max(mp.mpf('.48'), abs(mp.mpf(result['u_g_max_m_s'])), abs(drift))
            guard = 128 * mp.mpf(2)**(-52) * scale
            lower = mp.mpf(result['minimum_lower_bound_m_s'])
            robust = result['excluded'] and lower > guard
            records.append({
                'method': row['method'], 'representation': representation,
                'physical_status_double': int(row['physical_status']),
                'margin_double_m_s': row['margin_si'], 'norm_double': row['residual_norm'],
                'failure_index': row['failure_index'], 'failure_stage': row['failure_stage'],
                'input_bit_patterns': {key: row[key + '_bits'] for key in ('theta', 'U1', 'U2', 'U3')},
                'roundoff_separation_guard_m_s': mp.nstr(guard, 35),
                'robust_sufficient_exclusion': bool(robust),
                'classification': ('ROBUST_POSITIVE_BOUND_FOR_THIS_PROPOSAL' if robust
                                   else 'MARGINAL_OR_ADMISSIBLE_NOT_EXCLUDED'),
                'high_precision': result,
            })
    report = {
        'schema': 'primsf_active_terminal_80digit_v1', 'source_commit': build['git_head'],
        'focus_cell_from_terminal_sidecar': focus,
        'snapshot_sha256': digest(snapshot), 'build_manifest_sha256': digest(args.buildmanifest),
        'smoke_manifest_sha256': digest(args.smoke_manifest),
        'probe_source_sha256': digest(probe), 'api_sha256': digest(api),
        'probe_binary_sha256': digest(out / 'probe'),
        'object_sha256': {path.name: digest(path) for path in chosen},
        'compile_command': command, 'reference_sha256': digest(source / 'scripts/replay_root_bounds.py'),
        'coefficient_reference_sha256': digest(source / 'scripts/generate_bhagwat_fixtures.py'),
        'decimal_digits': 80, 'controls_nominal_SI': controls,
        'input_interpretations': 'Stored decimal literals and exact binary64 U/theta; physical controls are nominal decimals in both.',
        'limits': 'Only the named local proposals are assessed. Sub-ULP signs and double status1017 alone are not a mathematical exclusion proof; no global nonexistence claim.',
        'records': records,
    }
    (out / 'high_precision.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(out), 'snapshot_sha256': report['snapshot_sha256'],
                      'records': len(records), 'robust_exclusions': sum(row['robust_sufficient_exclusion'] for row in records)}))


if __name__ == '__main__':
    main()
