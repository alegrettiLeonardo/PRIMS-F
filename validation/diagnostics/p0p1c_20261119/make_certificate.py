#!/usr/bin/env python3
"""Assemble a blockage certificate from certificate.f90 output and run logs.

Usage: make_certificate.py NAME RAW RUN_DIR OUT_JSON [OUT_MD]

Classification rules (each recorded with the evidence that triggered it):
* where: the accepted state and the last trial are evaluated; if both are
  admissible the blockage is in a Newton PROPOSAL, not an accepted state. A
  failure of the Jacobian evaluation at the trial is an AUXILIARY evaluation.
* code decoding: face 3000+1000+ic -> closure status ic at the junction trace;
  cell 1000+10+pc -> primitive recovery status pc (module_find).
* CONSTITUTIVE_DOMAIN_LIMIT: closure status 12 (CLOSURE_DOMAIN) at the first
  refused point of the Newton ray, with an admissible accepted state.
* NO_ROOT_IN_DECLARED_BRANCH: primitive status 7, proven by the analytic bound
  of the declared co-current branch (not by sampling).
* ROOT_SEARCH_FAILURE: primitive status 4 (bracketing failed without proof).
* INADMISSIBLE_STATE: refused state with domain/EOS status (1, 2, 8, 89...).
* JACOBIAN_FAILURE: Jacobian assembly fails or FD/assembled mismatch > 1e-4.
* GLOBALIZATION_FAILURE: an admissible ray segment reduces the residual by
  more than 1% (descent existed) but the run still stopped.
* UNSUPPORTED_PHYSICAL_BOUNDARY: the trace required by the Newton root lies in
  a physical regime outside every closure of the operator (here: zero/negative
  junction gas flux under ascending liquid, i.e. face gas holdup -> 0).
Absence of a root along local proposals is NEVER reported as proof of the
nonexistence of a global solution.
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

PRIMITIVE = {0: 'OK', 1: 'DOMAIN', 2: 'EOS', 3: 'CLOSURE', 4: 'ROOT_SEARCH_FAILURE', 5: 'SINGULAR',
             6: 'UNSUPPORTED_SWITCH', 7: 'NO_PHYSICAL_ROOT', 8: 'INADMISSIBLE_ROOT'}


def parse(raw):
    rec = {}
    for line in Path(raw).read_text().splitlines():
        if '|' not in line:
            continue
        key, *vals = [x.strip() for x in line.split('|')]
        rec.setdefault(key, []).append(vals)
    return rec


def nums(s):
    return [float(x) for x in s.split()]


def decode(code, stage):
    """assemble(): CELL 1000+point_state, FACE 2000+face_pair, BOUNDARY 3000+boundary."""
    if stage == 'FACE':
        inner = code - 2000
        for lo, site, key in [(1400, 'junction steady trace', 'status'), (1200, 'Roe matrix on the path', 'status'),
                              (1100, 'junction momentum trace', 'status'), (1000, 'junction riser closure', 'closure_status'),
                              (900, 'junction pipe closure', 'closure_status'), (800, 'right trace mapping', 'closure_status')]:
            if inner >= lo:
                return {'site': site, key: inner - lo}
        return dict(site='face', status=inner)
    if stage == 'CELL':
        inner = code - 1000
        if 10 <= inner < 20:
            return dict(site='cell primitive recovery', primitive_status=inner - 10, primitive_name=PRIMITIVE.get(inner - 10))
        return dict(site='cell point state', status=inner)
    return dict(site=str(stage), status=code)


def main():
    name, raw, run, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    r = parse(raw)
    cert = dict(case=name, raw=str(raw), run=str(run))
    states = {int(v[0]): dict(time_s=nums(v[1])[0], dt_s=nums(v[1])[1], status=nums(v[2])[0],
                              index=nums(v[2])[1], stage=nums(v[2])[2], residual_norm=nums(v[3])[0])
              for v in r['STATE']}
    cert['accepted_state'], cert['last_trial'] = states[0], states[1]
    jac = r['JACOBIAN'][0]
    cert['jacobian'] = dict(lapack_info=int(jac[0]), norm1=float(jac[1]), condition_1norm_estimate=float(jac[2]))
    fd = [float(v[3]) for v in r['JACFD'] if float(v[3]) >= 0]
    cert['jacobian']['fd_max_relative_error_local_columns'] = max(fd) if fd else None
    cert['jacobian']['fd_columns_without_admissible_stencil'] = sum(float(v[3]) < 0 for v in r['JACFD'])
    cert['residual_dominant'] = [dict(cell=int(v[1]), component=int(v[2]), scaled=float(v[3])) for v in r['RESID']]
    nw = r['NEWTON'][0]
    cert['newton'] = dict(linear_info=int(nw[0]), max_scaled_step=float(nw[1]), at_cell=int(nw[2]), component=int(nw[3]))
    ray = r['RAY'][0]
    cert['ray'] = dict(admissible_samples=int(ray[0]), best_norm=float(ray[1]), best_lambda=float(ray[2]),
                       best_norm_ratio=float(ray[3]))
    b = r['BOUNDARY'][0]
    if b[0] != 'NONE':
        lam = nums(b[0])
        code, idx, stage = (int(x) for x in nums(b[1]))
        st = {1: 'CELL', 2: 'FACE', 3: 'BOUNDARY'}.get(stage, stage)
        cert['boundary'] = dict(lambda_admissible=lam[0], lambda_refused=lam[1], status=code, index=idx,
                                stage=st, norm_ratio_at_edge=float(b[2]), decoded=decode(code, st))
    if 'LOCATION' in r:
        loc = r['LOCATION'][0]
        cert['location'] = dict(index=int(loc[0]), s_m=float(loc[1]), inclination_deg=float(loc[2]), width_m=float(loc[3]))
    tags = {0: 'accepted', 1: 'trial', 2: 'ray_edge_admissible', 3: 'ray_first_refused'}
    cells = []
    for v in r.get('CELL', []):
        u, p, ev = nums(v[4]), nums(v[5]), nums(v[6])
        cells.append(dict(state=tags[int(v[0])], index=int(v[1]), closure=int(v[2]), status=int(v[3]),
                          U_SI=u, alpha=p[0], rho_l=p[1], rho_g=p[2], p_Pa=p[3], u_l=p[4], u_g=p[5],
                          eigen_real_m_s=ev[:3], eigen_imag_m_s=ev[3:],
                          hyperbolic=all(abs(x) == 0 for x in ev[3:]) if int(v[3]) == 0 else None))
    cert['cells'] = cells
    faces = []
    for v in r.get('FACE', []):
        if len(v) < 3:
            faces.append(dict(state=tags[int(v[0])], status=int(v[1])))
            continue
        f = nums(v[2])
        faces.append(dict(state=tags[int(v[0])], riser_closure_status=int(v[1]), conv_Gl=f[0], conv_Gg=f[1],
                          diss_Gl=f[2], diss_Gg=f[3], Gl=f[4], Gg=f[5], Gg_pipe_cell=f[6], Gg_riser_trace=f[7],
                          p_star_Pa=f[8], alpha_riser_face=float(v[3]) if v[3] != 'NaN' else None,
                          units='kg/m2/s'))
    cert['faces'] = faces
    cert['root_maps'] = [dict(state=tags[int(v[0])], cell=int(v[1]), alpha=float(v[2]), interior_sign_changes=int(v[3]),
                              closure_domain_samples=int(v[4]), phi_ug0_ugmax=nums(v[5]),
                              min_abs_phi_over_alpha=nums(v[6])[0], at_ug_m_s=nums(v[6])[1]) for v in r.get('ROOTMAP', [])]
    m = r['MARGIN'][0]
    cert['margins'] = dict(status=int(m[0]), junction_raw_gas_superficial_nd=float(m[1]),
                           min_riser_drift_margin_m_s=float(m[2]), at_cell=int(m[3]))
    # run history
    steps = list(csv.DictReader((run / 'steps.csv').open()))
    trials = list(csv.DictReader((run / 'rejected_proposals.csv').open()))
    status = dict(l.split('=', 1) for l in (run / 'run_status.txt').read_text().splitlines() if '=' in l)
    dts = [float(s['dt_s']) for s in steps]
    tried = [float(t['attempt_dt_s']) for t in trials]
    cert['run'] = dict(status={k.strip(): v.strip() for k, v in status.items()}, accepted_steps=len(steps),
                       min_dt_accepted_s=min(dts), last_dt_accepted_s=dts[-1],
                       min_dt_tried_s=min(tried) if tried else None,
                       rejected_proposal_codes=dict(Counter(f"{t['assembly_status'].strip()}@{t['stage'].strip()}{t['extended_index'].strip()}"
                                                            for t in trials)),
                       last_trial_codes=sorted({t['assembly_status'].strip() for t in trials[-40:]}),
                       reduction_reason='dt halved after each non-converged Newton sequence (dt_shrink); see rejected_proposals.csv')
    # classification
    cls, why = [], []
    acc_ok = cert['accepted_state']['status'] == 0
    trial_ok = cert['last_trial']['status'] == 0
    where = 'NEWTON_PROPOSAL' if acc_ok and trial_ok else ('ACCEPTED_STATE' if not acc_ok else 'TRIAL_STATE')
    if cert['jacobian']['lapack_info'] != 0 or (cert['jacobian']['fd_max_relative_error_local_columns'] or 0) > 1e-4:
        cls.append('JACOBIAN_FAILURE')
    else:
        why.append('Jacobian assembled at the trial, LU nonsingular, local FD agreement '
                   f"{cert['jacobian']['fd_max_relative_error_local_columns']:.1e}: not a Jacobian failure")
    if cert['ray']['best_norm_ratio'] < .99:
        cls.append('GLOBALIZATION_FAILURE')
        why.append('admissible descent exists on the Newton ray')
    else:
        why.append(f"admissible part of the Newton ray reduces the residual only to {cert['ray']['best_norm_ratio']:.6f}: "
                   'no admissible descent, not a line-search weakness')
    d = cert.get('boundary', {}).get('decoded', {})
    if d.get('closure_status') == 12:
        cls.append('CONSTITUTIVE_DOMAIN_LIMIT')
        why.append('first refused point: CLOSURE_DOMAIN of the junction riser closure')
        ff = [f for f in faces if f['state'] == 'ray_first_refused' and 'Gg' in f]
        fe = [f for f in faces if f['state'] == 'trial' and 'Gg' in f]
        if ff and ff[0]['Gg'] <= 0 and fe and (fe[0]['alpha_riser_face'] or 1) < 1e-6:
            cls.append('UNSUPPORTED_PHYSICAL_BOUNDARY')
            why.append(f"junction gas flux {fe[0]['Gg']:.2e} -> {ff[0]['Gg']:.2e} kg/m2/s with face riser holdup "
                       f"{fe[0]['alpha_riser_face']:.2e}: required mixed->pure liquid transition of the face trace "
                       'while both cells stay two-phase')
            if fe[0]['Gg_pipe_cell'] > 0 and fe[0]['Gg_riser_trace'] > 0:
                why.append(f"both traces still carry gas upward ({fe[0]['Gg_pipe_cell']:.3e}, {fe[0]['Gg_riser_trace']:.3e} "
                           f"kg/m2/s): the vanishing face flux comes from the Roe dissipation (diss_Gg={fe[0]['diss_Gg']:.3e}) "
                           'at the junction holdup jump; a PHYSICAL blockage is NOT established by these states')
                cert.setdefault('flags', []).append('NUMERICAL_FLUX_INDUCED_BOUNDARY')
    if d.get('primitive_status') == 7:
        cls.append('NO_ROOT_IN_DECLARED_BRANCH')
        why.append('primitive recovery status 7: analytic exclusion over the whole co-current interval')
        rm = [x for x in cert['root_maps'] if x['state'] == 'ray_first_refused']
        if rm and rm[0]['interior_sign_changes'] == 0:
            cls.append('CONSTITUTIVE_DOMAIN_LIMIT')
            why.append('root of the trial lies at the co-current endpoint (u_l=0); the proposal requires u_l<0')
    if d.get('primitive_status') == 4:
        cls.append('ROOT_SEARCH_FAILURE')
    if d.get('primitive_status') in (1, 2, 8):
        cls.append('INADMISSIBLE_STATE')
    cert['classification'] = dict(blockage_in=where, classes=cls, flags=cert.get('flags', []), evidence=why,
                                  not_implied='absence of roots along local Newton proposals is not a proof that no '
                                              'global solution exists; the model, not the physics, is certified here')
    out.write_text(json.dumps(cert, indent=1) + '\n')
    if len(sys.argv) > 5:
        md = [f"# Certificado de bloqueio — {name}", '', f"Execução: `{run}`", '',
              f"* Bloqueio em: **{where}**; classes: **{', '.join(cls)}**",
              f"* t = {states[0]['time_s']:.10f} s; último dt aceito {cert['run']['last_dt_accepted_s']:.3e} s; "
              f"menor tentado {cert['run']['min_dt_tried_s']:.3e} s",
              f"* Local: {cert.get('location')}", f"* Fronteira do raio: {cert.get('boundary')}",
              f"* Jacobiano: {cert['jacobian']}", f"* Raio: {cert['ray']}", f"* Margens: {cert['margins']}", '',
              '## Evidências', *[f"* {w}" for w in why], '', f"_{cert['classification']['not_implied']}_", '']
        Path(sys.argv[5]).write_text('\n'.join(md))
    print(name, where, cls)


if __name__ == '__main__':
    main()
