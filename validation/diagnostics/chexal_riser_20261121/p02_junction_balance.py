#!/usr/bin/env python3
"""P0.2 summary of the junction transmission study with the CHEXAL riser.

Inputs: evidence/p02/junction_<fixture>.log (junction_riemann_ref JUNCTION mode,
accepted state `which=0`, CENTER traces), evidence/p02/multistart_<fixture>.log,
and the fixtures (offline_api format) for the distributed-force estimate.

Per fixture:
  T1  (Gl, Gg, F3 continuous): Gg*, p_L*, p_R*, Delta p = p_R* - p_L*, status
  T2  (Gl, Gg, p continuous):  Gg*, J_M = F3_R* - F3_L* (the singular momentum
      source T2 needs; absent from Burr's eqs. 72-75), status
  ROE operator (current PRIMS MATCHED): Gg*
  distributed force of the SAME model: mean |dp/ds| over the first four riser
  cells of the accepted state (BAROTROPIC EOS from U) and the length
  l_eq = |J_M| / |dp/ds| over which the model's own distributed sources would
  produce the T2 point force (also in pipe diameters).
Multistart: number of starts, converged solutions (distinct Gg*) per condition.
Usage: p02_junction_balance.py FIXTURE_ROOT EVIDENCE_P02_DIR OUT_JSON OUT_MD
"""
import json, math, sys
from pathlib import Path
import numpy as np

OMEGA_RHO, OMEGA_P, LR, N = 997.05, 101325.0, 12.546253993197080, 51
NCG, NCL = 343.0 * math.sqrt(OMEGA_RHO / OMEGA_P), 1498.0 * math.sqrt(OMEGA_RHO / OMEGA_P)


def ref_lines(path):
    out = {}
    for l in Path(path).read_text().splitlines():
        f = l.split('|')
        if f[0] == 'REF' and f[1] == '0' and f[2] == 'CENTER':
            v = [float(x) for x in f[5].split()]
            out[f[3]] = dict(status=int(f[4]), Gl=v[0], Gg=v[1], F3L=v[2], F3R=v[3], pL=v[4], pR=v[5], aL=v[6], aR=v[7],
                             lam=v[8:11])
        if f[0] == 'OP' and f[1] == '0' and f[2] == 'CENTER' and f[3] == 'ROE':
            v = [float(x) for x in f[5].split()]
            out['ROE'] = dict(status=int(f[4]), p=v[0], Gl=v[1], Gg=v[2])
    return out


def pressure(u1, u2):
    lo, hi = 1e-4, 100.0
    for _ in range(200):
        p = .5 * (lo + hi); f = u1 / (1 + (p - 1) / NCL ** 2) + u2 / (p / NCG ** 2) - 1
        lo, hi = (p, hi) if f > 0 else (lo, p)
    return p * OMEGA_P


def riser_gradient(fixture):
    d = np.loadtxt(Path(fixture) / 'last_rejected_state.dat', skiprows=2)
    ps = [pressure(r[3], r[4]) for r in d[N:N + 5]]
    ds = LR / N
    return float(np.mean(np.abs(np.diff(ps))) / ds)


def multistart(path):
    if not Path(path).exists():
        return None
    res = {}
    for l in Path(path).read_text().splitlines():
        f = l.split('|')
        if f[0] != 'MS':
            continue
        t, st = f[1], int(f[3]); v = [float(x) for x in f[4].split()]
        r = res.setdefault(t, dict(starts=0, converged=0, gg=set(), alpha0_limit=0))
        r['starts'] += 1
        if st == 0:
            r['converged'] += 1; r['gg'].add(round(v[1], 9))
        if st == 50:
            r['alpha0_limit'] += 1
    for r in res.values():
        r['distinct_Gg'] = sorted(r.pop('gg'))
    return res


def main():
    fx_root, ev, oj, om = map(Path, sys.argv[1:5])
    rows = []
    for log in sorted(ev.glob('junction_*.log')):
        name = log.stem.replace('junction_', '')
        r = ref_lines(log); t1, t2, roe = r.get('T1'), r.get('T2'), r.get('ROE')
        g = riser_gradient(fx_root / name)
        row = dict(fixture=name, T1=t1, T2=t2, ROE=roe, riser_dpds_Pa_m=g, multistart=multistart(ev / f'multistart_{name}.log'))
        if t2 and t2['status'] == 0:
            jm = t2['F3R'] - t2['F3L']
            row.update(J_M_Pa=jm, l_eq_m=abs(jm) / g, l_eq_D=abs(jm) / g / 0.05254)
        if t1 and t1['status'] == 0:
            row['T1_dp_Pa'] = t1['pR'] - t1['pL']
        rows.append(row)
    json.dump(rows, open(oj, 'w'), indent=1)
    f = lambda v, fm='%.4g': '-' if v is None else fm % v
    with open(om, 'w') as m:
        m.write('| state | Gg* T1 | Gg* T2 | Gg* Roe | T1 status | T2 status | T1 p_R*-p_L* (Pa) | T2 J_M = F3_R-F3_L (Pa) | riser |dp/ds| (Pa/m) | l_eq (m) | l_eq/D | multistart T1 / T2 (converged, distinct Gg*) |\n')
        m.write('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n')
        for r in rows:
            ms = r['multistart']
            mst = '-' if not ms else '; '.join('%s %d/%d %s' % (k, v['converged'], v['starts'], ['%.4g' % x for x in v['distinct_Gg']])
                                              for k, v in sorted(ms.items()))
            m.write('| %s | %s | %s | %s | %s | %s | %s | %s | %.0f | %s | %s | %s |\n' % (
                r['fixture'], f(r['T1'] and r['T1']['Gg']), f(r['T2'] and r['T2']['Gg']), f(r['ROE'] and r['ROE']['Gg']),
                r['T1'] and r['T1']['status'], r['T2'] and r['T2']['status'], f(r.get('T1_dp_Pa'), '%.0f'), f(r.get('J_M_Pa'), '%.0f'),
                r['riser_dpds_Pa_m'], f(r.get('l_eq_m'), '%.2f'), f(r.get('l_eq_D'), '%.1f'), mst))
    print(open(om).read())


if __name__ == '__main__':
    main()
