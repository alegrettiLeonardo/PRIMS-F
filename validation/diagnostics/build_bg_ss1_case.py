#!/usr/bin/env python3
"""Freeze and diagnose the SS1 junction for an arbitrary saved case mesh.

Usage: FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/build_bg_ss1_case.py RUN_DIRECTORY
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
from datetime import datetime, timezone

here = Path(__file__).resolve().parent
root = here.parents[1]

def metadata(path):
    text = path.read_text()
    params = {}
    for match in re.finditer(r'^\s*(\w+)\s*=\s*(.*?),?\s*$', text, re.M):
        key, value = match[1].lower(), match[2].strip().rstrip(',').strip()
        if value.startswith('"'):
            params[key] = value.strip('" ').strip()
        elif value in ('T', 'F', '.true.', '.false.'):
            params[key] = value in ('T', '.true.')
        else:
            try:
                params[key] = float(value.replace('D', 'E').replace('d', 'e'))
            except ValueError:
                params[key] = value
    line = next(line for line in text.splitlines() if line.startswith('Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u:'))
    lr, ca, cg, gl, gg, wp, wr, wu = map(float, line.split(':', 1)[1].split())
    return params, (lr, ca, cg, gl, gg, wp, wr, wu)

def summarize(log, params, source, frozen):
    profiles = {}
    summary = {'case': source.name, 'physical_cells': int(params['nodes']+params['n_riser']),
        'junction_extended_cells': [int(params['nodes']+params['num_dummy_cells']),
                                   int(params['nodes']+params['num_dummy_cells']+1)],
        'snapshot_sha256': frozen['last_rejected_state.dat'],
        'evidence_class': 'local_diagnosis_not_experimental_validation',
        'profiles': profiles, 'newton_candidates': [], 'point_failures': []}
    def profile(token):
        return profiles.setdefault(str(int(token)), {})
    def finite(token):
        value = float(token)
        return value if math.isfinite(value) else None
    for line in log.splitlines():
        parts = line.split()
        if line.startswith('snapshot time_s'):
            summary.update(time_s=float(parts[-4]), dt_s=float(parts[-3]),
                           saved_residual_norm=float(parts[-2]), saved_status=int(parts[-1]))
        elif line.startswith('minimum physical alpha'):
            profile(parts[-2])['minimum_cell_alpha'] = float(parts[-1])
        elif line.startswith('point failure'):
            summary['point_failures'].append({'profile':int(parts[-3]),'cell':int(parts[-2]),'status':int(parts[-1])})
        elif line.startswith('inlet boundary') or line.startswith('outlet boundary'):
            profile(parts[-2])[parts[0]+'_boundary_status'] = int(parts[-1])
        elif line.startswith('junction pstar'):
            profile(parts[-4]).update(pressure_Pa=float(parts[-3]), jl_m_s=float(parts[-2]), jg_m_s=float(parts[-1]))
        elif line.startswith('junction gas convective'):
            profile(parts[-4]).update(central_gas_mass_flow_kg_s=float(parts[-3]),
                dissipative_gas_mass_flow_kg_s=float(parts[-2]),gas_mass_flow_kg_s=float(parts[-1]))
        elif line.startswith('riser trace root'):
            profile(parts[-4]).update(trace_status=int(parts[-3]),trace_nroots=int(parts[-2]),alpha_riser_trace=finite(parts[-1]))
        elif line.startswith('neighbor velocities'):
            profile(parts[-7]).update(dict(zip(['ul_pipe_m_s','ug_pipe_m_s','ul_riser_m_s','ug_riser_m_s','alpha_pipe_cell','alpha_riser_cell'],
                map(float,parts[-6:]))))
        elif line.startswith('assembly residual'):
            profile(parts[-3]).update(assembly_status=int(parts[-2]),residual_norm=float(parts[-1]))
        elif line.startswith('Newton band solve'):
            summary.update(linear_solve_status=int(parts[-2]),max_newton_increment=float(parts[-1]))
        elif line.startswith('Newton candidate'):
            summary['newton_candidates'].append({'lambda':float(parts[-3]),'cell_point_failures':int(parts[-2]),'junction_status':int(parts[-1])})
    summary['interpretation'] = ('The snapshot and candidates distinguish cell admissibility from the MATCHED Roe trace domain. '
        'A negative proposed numerical gas flux alone does not prove physical reverse-gas flow.')
    return summary

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_directory', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    source = args.run_directory.resolve()
    out = (args.output or here/'final_ss1_cases'/source.name).resolve()
    out.mkdir(parents=True, exist_ok=True)
    frozen = {}
    for name in ['last_rejected_state.dat', 'run_metadata.txt', 'run_status.txt', 'manifest.json', 'monitor.csv', 'stdout.log']:
        p = source/name
        if p.exists():
            destination = out/name
            if destination.exists():
                if destination.read_bytes() != p.read_bytes():
                    raise SystemExit(f'Refusing to overwrite a changed frozen input: {destination}; use --output with a new directory.')
            else:
                shutil.copyfile(p, destination)
            frozen[name] = hashlib.sha256(destination.read_bytes()).hexdigest()
    params, scales = metadata(out/'run_metadata.txt')
    lr, ca, cg, gl, gg, wp, wr, wu = scales
    nodes = int(params['nodes']);nriser = int(params['n_riser']);nd = int(params['num_dummy_cells'])
    if int(params['closure_riser']) != 3 or params['junction_model'] != 'MATCHED' or params['dissipation_mode'] != 'LAGGED':
        raise SystemExit('This probe qualifies ID3 / MATCHED / LAGGED runs only.')
    values = [
        f'{nodes} {nriser} {nd}',
        ' '.join(map(repr, [wu, wr, wp, lr, params['lp']])),
        ' '.join(map(repr, [params[k] for k in ['d', 'eps', 'g', 'mu_l', 'mu_g', 'sigma']])),
        ' '.join(map(repr, [params['cl']*math.sqrt(wr/wp), cg*math.sqrt(wr/wp), params['rho_l0']/wr, params['p_l0']/wp, params['ps']/wp])),
        ' '.join(map(repr, [gl/(wr*wu), gg/(wr*wu), params['jl']/wu, params['jg']/wu])),
        ' '.join(map(repr, [params[k] for k in ['tol', 'tola', 'tol_res', 'tol_inc']])),
        f'1 3 {int(params["quad_points"])} 1000',
        'T' if params['jac_dissipation_fd'] else 'F',
        str(params['junction_model']), str(params['flux_scheme']), str(params['dissipation_mode']),
    ]
    (out/'case_config.dat').write_text('\n'.join(values)+'\n')
    (out/'snapshot_provenance.json').write_text(json.dumps({'frozen_utc':datetime.now(timezone.utc).isoformat(),
        'source':str(source),'files':frozen,'parameters':params,'dimensional_scales':scales},indent=2)+'\n')
    build = out/'build';build.mkdir(exist_ok=True)
    names = ('mod_kinds mod_rootfinding mod_void_helpers module_flow_regime module_config module_lapack_small '
             'module_closure_law module_find module_catenary module_steady_state module_jacobian module_matrix_roe '
             'module_roe_solver module_boundary_conditions module_solver_utils module_time_step module_simulate_pipe_riser').split()
    files = [root/'scr'/(name+'.f90') for name in names]+[here/'probe_bg_ss1_case.f90']
    command = shlex.split(os.environ.get('FC','gfortran'))+['-std=f2018','-O0','-g','-ffree-line-length-none',
        '-fcheck=all,no-array-temps','-ffpe-trap=invalid,zero,overflow','-finit-real=snan','-fopenmp',
        '-J'+str(build),'-I'+str(build),*map(str,files),'-Wl,-z,noexecstack',
        *shlex.split(os.environ.get('LIBS','-llapack -lblas')),'-o',str(build/'probe_bg_ss1_case')]
    result = subprocess.run(command,capture_output=True,text=True,timeout=180)
    (build/'build.log').write_text(result.stdout+result.stderr)
    provenance = {'command':command,'build_returncode':result.returncode,
        'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    if result.returncode == 0:
        result = subprocess.run([str(build/'probe_bg_ss1_case')],cwd=out,capture_output=True,text=True,timeout=120)
        (out/'probe.log').write_text(result.stdout);(build/'runtime_stderr.log').write_text(result.stderr)
        provenance['probe_returncode'] = result.returncode
        (out/'summary.json').write_text(json.dumps(summarize(result.stdout,params,source,frozen),indent=2,allow_nan=False)+'\n')
        print(result.stdout)
    else:
        print(result.stderr)
    (build/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Evidence directory:',out)
    return result.returncode

if __name__ == '__main__':
    raise SystemExit(main())
