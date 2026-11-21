#!/usr/bin/env python3
"""P0-B/C table: global T1 face (t1_face_check) versus the independent local T1
reference (junction_riemann_ref, which=0 CENTER) and the global Roe/T2 face, on
the SAME accepted conservative states; momentum identity of the assembled residual.
Usage: t1_face_summary.py EVIDENCE_DIR OUT_JSON OUT_MD"""
import json, sys
from pathlib import Path
ev = Path(sys.argv[1]); rows = []
for f in sorted((ev / 't1').glob('face_*.log')):
    name = f.stem.replace('face_', ''); d = {'state': name}
    for l in f.read_text().splitlines():
        p = l.split('|')
        if p[0] == 'FACE' and len(p) >= 4 and p[2].strip().lstrip('-').isdigit():
            v = [float(x) for x in p[3].split()]
            d[p[1]] = dict(status=int(p[2]), Gl=v[0], Gg=v[1], F3L=v[2], F3R=v[3], pL=v[4], pR=v[5], aL=v[6], aR=v[7])
        if p[0] == 'JUMP':
            d['jump_' + p[1]] = [float(x) for x in p[2].split()]
        if p[0].startswith('BAL') and 'identity' in l:
            d['identity_' + p[1]] = float(l.split('identity residual')[1].split('|')[0].split()[2])
    r = (ev / 'p02' / f'junction_{name}.log').read_text().splitlines()
    for l in r:
        p = l.split('|')
        if p[:3] == ['REF', '0', 'CENTER'] and p[3] == 'T1':
            v = [float(x) for x in p[5].split()]
            d['REF_T1'] = dict(status=int(p[4]), Gl=v[0], Gg=v[1], F3=v[2], pL=v[4], pR=v[5], aR=v[7])
    rows.append(d)
json.dump(rows, open(sys.argv[2], 'w'), indent=1)
f = lambda x, fm='%.5g': '-' if x is None else fm % x
with open(sys.argv[3], 'w') as m:
    m.write('| state | Gg* T1 global | Gg* T1 reference | rel. diff | Gl* rel. diff | p_R*-p_L* global (Pa) | reference (Pa) | Gg* Roe/T2 (status) | T1 jump (3 comp.) | T1 momentum identity | T2 momentum identity |\n')
    m.write('|---|---:|---:|---:|---:|---:|---:|---|---|---:|---:|\n')
    for d in rows:
        g, rf, t2 = d.get('T1'), d.get('REF_T1'), d.get('T2')
        rel = abs(g['Gg'] - rf['Gg']) / abs(rf['Gg']) if g and rf else None
        rgl = abs(g['Gl'] - rf['Gl']) / abs(rf['Gl']) if g and rf else None
        m.write('| %s | %s | %s | %s | %s | %s | %s | %s (%s) | %s | %s | %s |\n' % (
            d['state'], f(g and g['Gg']), f(rf and rf['Gg']), f(rel, '%.1e'), f(rgl, '%.1e'),
            f(g and g['pR'] - g['pL'], '%.0f'), f(rf and rf['pR'] - rf['pL'], '%.0f'), f(t2 and t2['Gg']), t2 and t2['status'],
            d.get('jump_T1'), f(d.get('identity_T1'), '%.1e'), f(d.get('identity_T2'), '%.1e')))
print(open(sys.argv[3]).read())
