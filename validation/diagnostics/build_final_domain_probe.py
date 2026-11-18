#!/usr/bin/env python3
"""Diagnose saved SS3 states using each run's mesh and dimensional metadata.

No transient is executed and no production source is modified. Snapshots are
copied separately per case. Example:
  FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/build_final_domain_probe.py \
    validation/results/final_D_SS3_C3_BAROTROPIC_MATCHED_LAGGED_ROE_N51_51_dt0.5
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess


SOURCES = ('mod_kinds mod_rootfinding mod_void_helpers module_flow_regime module_config '
           'module_lapack_small module_closure_law module_find module_catenary module_steady_state '
           'module_jacobian module_matrix_roe module_roe_solver module_boundary_conditions '
           'module_solver_utils module_time_step module_simulate_pipe_riser').split()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metadata(path: Path) -> tuple[dict, dict]:
    data = {}
    text = path.read_text()
    for line in text.splitlines():
        match = re.match(r'\s*(\w+)\s*=\s*(.*?)\s*,?\s*$', line)
        if not match:
            continue
        key, raw = match.groups()
        raw = raw.rstrip(',').strip()
        if raw.startswith(('"', "'")):
            value = raw.strip('"\'').strip()
        elif raw in ('T', 'F', '.TRUE.', '.FALSE.'):
            value = raw in ('T', '.TRUE.')
        else:
            value = float(raw.replace('D', 'E').replace('d', 'e'))
        data[key.lower()] = value
    line = next(line for line in text.splitlines() if line.startswith('Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u:'))
    scale_names = 'Lr CA Cg Gl Gg omega_P omega_rho omega_u'.split()
    values = [float(value) for value in line.split(':', 1)[1].split()]
    if len(values) != len(scale_names):
        raise ValueError('The scaling metadata must have eight values.')
    return data, dict(zip(scale_names, values))


def read_metrics(path: Path) -> dict:
    result = {}
    for line in path.read_text().splitlines():
        key, raw = line.split('=', 1)
        raw = raw.strip()
        if re.fullmatch(r'[-+]?\d+', raw):
            result[key] = int(raw)
        else:
            value = float(raw.replace('D', 'E'))
            result[key] = value if math.isfinite(value) else None
    return result


def main() -> int:
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case', nargs='+', type=Path)
    parser.add_argument('--output-dir', type=Path, default=here / 'final_ss3_cases')
    parser.add_argument('--snapshot-source', choices=('rejected', 'final-fields'), default='rejected',
                        help='Use original rejected state or reconstruct the final physical state from fields.csv.')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    build = output / 'build'
    build.mkdir(parents=True, exist_ok=True)
    sources = [root / 'scr' / f'{name}.f90' for name in SOURCES]
    sources.append(here / 'probe_bg_domain.f90')
    source_hashes = {str(path.relative_to(root)): sha(path) for path in sources}
    exe = build / 'probe_bg_domain'
    command = [*shlex.split(os.environ.get('FC', 'gfortran')), '-std=f2018', '-O0', '-g',
               '-ffree-line-length-none', '-fopenmp', '-fcheck=all,no-array-temps',
               '-ffpe-trap=invalid,zero,overflow', '-J' + str(build), '-I' + str(build),
               *map(str, sources), '-Wl,-z,noexecstack',
               *shlex.split(os.environ.get('LIBS', '-llapack -lblas')), '-o', str(exe)]
    result = subprocess.run(command, cwd=build, text=True, capture_output=True, timeout=60)
    (build / 'build.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        print(result.stdout + result.stderr)
        return result.returncode
    build_provenance = {'utc': datetime.now(timezone.utc).isoformat(), 'command': command,
                        'source_sha256': source_hashes, 'binary_sha256': sha(exe)}
    (build / 'provenance.json').write_text(json.dumps(build_provenance, indent=2) + '\n')
    any_failure = False
    for case_arg in args.case:
        case = case_arg.resolve()
        dest = output / (case.name + ('__final_fields' if args.snapshot_source == 'final-fields' else ''))
        dest.mkdir(exist_ok=True)
        (dest / 'build_provenance.json').write_text(json.dumps(build_provenance, indent=2) + '\n')
        copied = {}
        for name in ('last_rejected_state.dat', 'run_metadata.txt', 'run_status.txt', 'manifest.json'):
            data = (case / name).read_bytes()
            (dest / name).write_bytes(data)
            copied[name] = hashlib.sha256(data).hexdigest()
        pars, scale = metadata(dest / 'run_metadata.txt')
        rows = [line.split() for line in (dest / 'last_rejected_state.dat').read_text().splitlines()[2:] if line.strip()]
        nphysical = len(rows)
        nghost, npipe = int(pars['num_dummy_cells']), int(pars['nodes'])
        nriser = int(pars['n_riser'])
        if nphysical != npipe + nriser:
            raise ValueError(f'{case.name}: snapshot/metadata cell counts disagree.')
        if [int(row[0]) for row in rows] != list(range(nghost + 1, nghost + nphysical + 1)):
            raise ValueError(f'{case.name}: unexpected physical-cell indices.')
        if int(pars['closure_riser']) != 3:
            raise ValueError('This probe diagnoses the ID3 closure only.')
        run_status = read_metrics(dest / 'run_status.txt')
        probe_input = dest / 'last_rejected_state.dat'
        reconstruction = None
        if args.snapshot_source == 'final-fields':
            field_path = case / 'fields.csv'
            all_rows = list(csv.DictReader(field_path.open()))
            final_time = max(float(row['time_s']) for row in all_rows)
            field_rows = [row for row in all_rows if float(row['time_s']) == final_time]
            if len(field_rows) != nphysical:
                raise ValueError('Final field row count differs from the snapshot mesh.')
            omega_t = scale['Lr'] / scale['omega_u']
            header = ['# Reconstructed final fields; dt/norm are zero diagnostic placeholders, not a time step.',
                      f'{final_time / omega_t:.17e} 0.0 0.0 {run_status["status"]}']
            for shape, field in zip(rows, field_rows):
                alpha, rhol, rhog = (float(field[key]) for key in ('alpha', 'rhol_kg_m3', 'rhog_kg_m3'))
                ul, ug = (float(field[key]) for key in ('ul_m_s', 'ug_m_s'))
                ml, mg = (1 - alpha) * rhol / scale['omega_rho'], alpha * rhog / scale['omega_rho']
                state = [ml, mg, (ml * ul + mg * ug) / scale['omega_u']]
                values = [float(field['theta_rad']), float(shape[2]), *state, *state, 0., 0., 0.]
                header.append(str(int(shape[0])) + ' ' + ' '.join(f'{value:.17e}' for value in values))
            probe_input = dest / 'final_state_from_fields.dat'
            probe_input.write_text('\n'.join(header) + '\n')
            reconstruction = {'fields_sha256': sha(field_path), 'final_time_s': final_time,
                              'reconstructed_state_sha256': sha(probe_input),
                              'qualification': 'Conservative state reconstructed from full-precision printed primitives; dt and norm are diagnostic placeholders.'}
        config = [
            [nphysical, nghost, npipe, int(pars['closure_riser']), int(pars['quad_points']), 1000],
            [scale['omega_u'], scale['omega_rho'], scale['omega_P'], scale['Lr'], pars['lp'],
             pars['d'], pars['eps'], pars['g']],
            [pars['mu_l'], pars['mu_g'], pars['sigma'], pars['cl'], scale['Cg'], pars['rho_l0'], pars['p_l0'], pars['ps']],
            [scale['Gl'], scale['Gg'], pars['jl'], pars['jg'], pars['tol'], pars['tola'],
             pars['tol_res'], pars['tol_inc'], pars['jac_dissipation_fd_eps']],
            [f'"{pars[key]}"' for key in ('flux_scheme', 'junction_model', 'dissipation_mode', 'gas_eos')],
        ]
        (dest / 'probe_config.dat').write_text('\n'.join(' '.join(map(str, line)) for line in config) + '\n')
        manifest = json.loads((dest / 'manifest.json').read_text())
        mismatches = [key for key, value in source_hashes.items()
                      if key.startswith('scr/') and manifest['source_sha256'].get(key) != value]
        run = subprocess.run([str(exe), str(probe_input)], cwd=dest, text=True, capture_output=True, timeout=60)
        (dest / 'probe.log').write_text(run.stdout + run.stderr)
        summary = {'case': case.name, 'original_case_path': str(case), 'snapshot_sha256': copied,
                   'metadata': pars, 'scales': scale, 'source_mismatches_against_run': mismatches,
                   'build_provenance': 'build_provenance.json', 'returncode': run.returncode}
        summary['run_status'] = run_status
        summary['input_state_source'] = args.snapshot_source
        summary['reconstruction'] = reconstruction
        metrics_path = dest / 'metrics.txt'
        if run.returncode or not metrics_path.exists():
            summary['classification'] = 'PROBE_FAILED_OR_BASE_OPERATOR_REJECTED'
            any_failure = True
            print(run.stdout + run.stderr)
        else:
            metrics = read_metrics(metrics_path)
            summary['metrics'] = metrics
            summary['snapshot_age_relative_to_run_end_s'] = run_status['reached_time_s'] - metrics['time_s']
            is_liquid_limit = (
                metrics['point_failures'] == 0 and metrics.get('prepare_dissipation_status') == 0 and
                metrics['max_spectral_imag'] == 0 and metrics['critical_ul_m_s'] < 1e-9 and
                metrics.get('candidate_limit_margin_nd', math.inf) < 0 and
                metrics.get('line_search_point_failures', 0) > 0 and
                metrics.get('line_search_common_failure_cell') == metrics['minimum_ul_cell_extended'])
            is_path_limit = (metrics['point_failures'] == 0 and
                             metrics.get('prepare_dissipation_status') == 1214 and
                             metrics.get('path_min_ulzero_margin_nd', 0) < 0)
            if is_path_limit:
                summary['classification'] = 'ROE_PATH_OUTSIDE_COCURRENT_DOMAIN'
            elif is_liquid_limit:
                summary['classification'] = 'CELL_COCURRENT_LIMIT_LIQUID_REVERSAL'
            elif metrics.get('line_search_junction_domain_failures', 0) > 0 and metrics['point_failures'] == 0:
                summary['classification'] = 'MATCHED_JUNCTION_COCURRENT_TRACE_LIMIT'
            else:
                summary['classification'] = 'REQUIRES_LOG_REVIEW'
        (dest / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
        print(json.dumps({'case': case.name, 'classification': summary['classification'],
                          'source_mismatches': mismatches, 'metrics': summary.get('metrics')}, indent=2))
    return 1 if any_failure else 0


if __name__ == '__main__':
    raise SystemExit(main())
