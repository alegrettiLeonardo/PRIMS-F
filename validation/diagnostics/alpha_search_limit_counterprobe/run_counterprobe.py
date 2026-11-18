#!/usr/bin/env python3
"""Synthetic local counterexample to treating finite alpha scans as root proofs.

Production is unchanged.  Pressure and jl are copied from the saved SS1 trace;
jg is deliberately reduced to half its tiny positive value.  This does not
reconstruct the exact rejected Newton candidate of the campaign.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shlex
import subprocess

import mpmath as mp

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
mp.mp.dps = 80
fixture_path = REPO / 'scripts/generate_bhagwat_fixtures.py'
spec = importlib.util.spec_from_file_location('independent_bhagwat', fixture_path)
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = mp.mpf('284053.9719206488')
    jl = mp.mpf('0.5295514577421165')
    jg_saved = mp.mpf('6.445998075334944e-13')
    rl = mp.mpf('997.05') + (p-mp.mpf('101325'))/mp.mpf('1498')**2
    rg = p/mp.mpf('343')**2
    amin = 1/(1+mp.exp(28))
    theta_cell = mp.mpf('4.0731295598015818e-2')
    records = []
    for location, theta in [('junction_trace', mp.mpf(0)), ('first_riser_cell_angle_comparison', theta_cell)]:
        for factor in (mp.mpf(1), mp.mpf('.5')):
            jg = factor*jg_saved

            def evaluate(a):
                c0, ud, groups = reference.coefficients([
                    a, rl, rg, jl/(1-a), jg/a, theta, '.05254', '.000046', '.00089', '.0728'])
                return a*(c0*(jl+jg)+ud)-jg, c0, ud

            lo, hi = mp.mpf('1e-16'), mp.mpf('1e-10')
            flo, fhi = evaluate(lo)[0], evaluate(hi)[0]
            assert flo < 0 < fhi
            for _ in range(210):
                mid = (lo+hi)/2
                if evaluate(mid)[0] > 0:
                    hi = mid
                else:
                    lo = mid
            alpha = (lo+hi)/2
            residual, c0, ud = evaluate(alpha)
            records.append({
                'id': len(records)+1, 'location': location, 'gas_factor_from_saved_trace': str(factor),
                'theta_rad': str(theta), 'p_Pa': str(p), 'jl_m_s': str(jl), 'jg_m_s': str(jg),
                'rho_l_kg_m3': mp.nstr(rl, 45), 'rho_g_kg_m3': mp.nstr(rg, 45),
                'independent_alpha': mp.nstr(alpha, 50), 'alpha_scan_min': mp.nstr(amin, 50),
                'independent_alpha_below_search_min': bool(alpha < amin),
                'C0_at_root': mp.nstr(c0, 45), 'Ud_at_root_m_s': mp.nstr(ud, 45),
                'residual_at_root_m_s': mp.nstr(residual, 20),
                'residual_at_scan_min_m_s': mp.nstr(evaluate(amin)[0], 40),
                'bracket_lower_alpha': '1e-16', 'bracket_upper_alpha': '1e-10',
                'bracket_lower_residual_m_s': mp.nstr(flo, 40),
                'bracket_upper_residual_m_s': mp.nstr(fhi, 40),
            })
    (HERE/'inputs.dat').write_text(''.join(
        ' '.join(str(r[k]) for k in ('theta_rad','p_Pa','jl_m_s','jg_m_s','independent_alpha'))+'\n'
        for r in records))
    build = HERE/'build'
    build.mkdir(exist_ok=True)
    source_files = [REPO/'scr'/(name+'.f90') for name in
                    ('mod_kinds','mod_void_helpers','module_closure_law','module_steady_state')]
    source_files.append(HERE/'probe.f90')
    command = shlex.split(os.environ.get('FC','gfortran'))+['-std=f2018','-O0','-g','-ffree-line-length-none',
               '-fcheck=all,no-array-temps','-ffpe-trap=invalid,zero,overflow',
               '-J'+str(build),'-I'+str(build),*map(str,source_files),'-o',str(build/'probe')]
    compiled = subprocess.run(command,capture_output=True,text=True,check=False)
    (build/'build.log').write_text(compiled.stdout+compiled.stderr)
    if compiled.returncode:
        raise RuntimeError(compiled.stderr)
    ran = subprocess.run([str(build/'probe')],cwd=HERE,capture_output=True,text=True,check=True)
    (HERE/'probe.log').write_text(ran.stdout)
    (build/'runtime_stderr.log').write_text(ran.stderr)
    rows = [line.split() for line in ran.stdout.splitlines() if line.strip()]
    assert len(rows) == len(records)
    for record, row in zip(records, rows):
        record['fortran_solver_status'] = int(row[1])
        record['fortran_nroots'] = int(row[2])
        record['fortran_alpha'] = None if row[3].lower() == 'nan' else float(row[3])
        record['fortran_closure_status_at_independent_root'] = int(row[4])
        record['fortran_residual_at_independent_root_m_s'] = float(row[5])
        record['fortran_status_with_independent_alpha_hint'] = int(row[6])
        assert record['fortran_closure_status_at_independent_root'] == 0
        if record['gas_factor_from_saved_trace'] == '0.5':
            assert record['fortran_solver_status'] == 21
            assert record['independent_alpha_below_search_min']
        else:
            assert record['fortran_solver_status'] == 0
    provenance = {'command':command, 'source_sha256':{str(f.relative_to(REPO)):digest(f) for f in source_files},
                  'independent_reference_sha256':digest(fixture_path),
                  'probe_script_sha256':digest(Path(__file__)), 'production_sources_modified':False}
    (HERE/'build_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    (HERE/'summary.json').write_text(json.dumps({
        'evidence_class':'synthetic_local_counterprobe_not_exact_rejected_campaign_candidate',
        'purpose':'Distinguish a root outside a finite search interval from a proven absence of physical closure roots.',
        'reference_source':reference.SOURCE,
        'reference_arithmetic_digits':mp.mp.dps,
        'production_search':'160 intervals in logit(alpha) from -28 to +28; no production edits.',
        'qualification':'Positive co-current roots below alpha_min exist in these constructed states; this does not establish the cause of every status21/3021 failure and does not permit negative gas flux in the qualified correlation.',
        'cases':records},indent=2)+'\n')
    for r in records:
        print(r['id'],r['location'],'factor',r['gas_factor_from_saved_trace'],
              'alpha',r['independent_alpha'],'solver status',r['fortran_solver_status'])


if __name__ == '__main__':
    main()
