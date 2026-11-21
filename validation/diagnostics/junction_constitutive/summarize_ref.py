#!/usr/bin/env python3
"""Summarize junction_riemann_ref.f90 JUNCTION logs: reference (T1, T2) vs Roe/E1.
Usage: summarize_ref.py OUT_JSON OUT_MD LABEL=LOG ..."""
import json, sys
ST = {0: 'OK', 4: 'NO_CONVERGENCE', 50: 'LIMIT_ALPHA_R0', 60: 'CONFIG_VIOLATED'}

def parse(path):
    d = {}
    for line in open(path):
        p = [x.strip() for x in line.split('|')]
        if p[0] == 'REF' and len(p) > 5:
            v = [float(x) for x in p[5].split()]
            d.setdefault((int(p[1]), p[2]), {})[p[3]] = dict(status=int(p[4]), Gl=v[0], Gg=v[1], F3L=v[2], F3R=v[3], pL=v[4], pR=v[5],
                                                            aL=v[6], aR=v[7], lam1L=v[8], lam2R=v[9], lam3M=v[10], its=v[11])
        elif p[0] == 'OP':
            v = [float(x) for x in p[5].split()]
            d.setdefault((int(p[1]), p[2]), {})[p[3]] = dict(status=int(p[4]), p=v[0], Gl=v[1], Gg=v[2], Gg_raw=v[3], momL=v[4], momR=v[5])
    return d

out = {}
for a in sys.argv[3:]:
    lab, path = a.split('=', 1)
    out[lab] = {'%d/%s' % k: v for k, v in parse(path).items()}
json.dump(out, open(sys.argv[1], 'w'), indent=1)
g = lambda r: '-' if r['status'] not in (0, 50) else '%.3g' % r['Gg']
with open(sys.argv[2], 'w') as m:
    m.write('# P0.2 junction gas flux: independent reference vs Roe and E1 (same cell states)\n\n')
    m.write('Gg in kg/m2s. Reference T1 = conservative layer (Gl, Gg, F3 continuous); T2 = PRIMS convention (Gl, Gg, p continuous).\n')
    m.write('LIMIT_ALPHA_R0 = Newton driven to alpha_R*<1e-14: no weak-wave solution with alpha_R*>0.\n\n')
    m.write('| state | input | REF T1 Gg (status, alpha_R*) | REF T1 pR*-pL* Pa | REF T2 Gg (status) | Roe Gg (status) | E1 Gg (status) | Roe point force momR-momL Pa |\n|---|---|---|---:|---|---|---|---:|\n')
    for lab, d in out.items():
        for key, r in d.items():
            t1, t2, roe, e1 = r['T1'], r['T2'], r['ROE'], r['E1']
            m.write('| %s | %s | %s (%s, %.2g) | %s | %s (%s) | %s (%d) | %s (%d) | %s |\n' % (
                lab, key, g(t1), ST.get(t1['status'], t1['status']), t1['aR'], '%.0f' % (t1['pR'] - t1['pL']) if t1['status'] == 0 else '-',
                g(t2), ST.get(t2['status'], t2['status']), '%.3g' % roe['Gg'] if roe['status'] == 0 else '-', roe['status'],
                '%.3g' % e1['Gg'] if e1['status'] == 0 else ('0 (exact)' if e1['Gg_raw'] != e1['Gg_raw'] and e1['Gl'] == e1['Gl'] else '-'), e1['status'],
                '%.0f' % (roe['momR'] - roe['momL']) if roe['status'] == 0 else '-'))
print(open(sys.argv[2]).read())
