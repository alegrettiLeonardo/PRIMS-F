#!/usr/bin/env python3
"""Fortran core (chexal_coefficients_si) versus the independent transcription
chexal92_ref.py on a grid of PRIMS-F air-water states (rho_l=998, rho_g 1.2-6,
D=0.05254 m), all riser orientations 0..90 deg, co-current up/down flows, and
inclined/horizontal countercurrent (both must refuse).
Usage: core_crosscheck.py FORTRAN_BINARY OUT_JSON"""
import itertools, json, math, subprocess, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import chexal92_ref as cx
rows = []
for a, rg, thd, (jl, jg) in itertools.product((0.02, 0.2, 0.5, 0.8, 0.97), (1.2, 2.4, 4.2),
        (0.0, 5.0, 30.0, 51.637, 76.47, 89.0, 90.0),
        ((0.547, 1.16), (0.53, 0.11), (1e-4, 0.6), (0.0, 0.3), (2.0, 3.0), (-0.3, -0.5), (-0.05, 0.4), (0.0, 0.0), (-0.2, 0.0))):
    rows.append((a, 998.0, rg, math.radians(thd), jl, jg))
# vertical countercurrent needs an explicit branch in the core (AUTO -> 62);
# the transcription's CC C3 is Eq.3-39/3-40, i.e. the HIGH branch: pass HIGH there.
vcc = lambda r: abs(r[3] - math.pi / 2) < 1e-12 and r[4] < 0 <= r[5]
inp = ''.join('%.17e %.17e %.17e %.17e %.17e %.17e %d\n' % (*r, 1 if vcc(r) else 0) for r in rows)
out = subprocess.run([sys.argv[1]], input=inp, capture_output=True, text=True, check=True).stdout.split('\n')
res = dict(n=len(rows), vertical_cc_high=sum(map(vcc, rows)), both_ok=0, both_refuse=0, mismatch_status=[], max_rel_C0=0.0, max_rel_Vgj=0.0, worst=None)
for r, line in zip(rows, out):
    st, c0f, vf = line.split(); st = int(st)
    try:
        c0p, vp = cx.coefficients(r[0], r[1], r[2], r[3], r[4], r[5])
        pst = 0
    except cx.DomainError:
        pst = 61
    if st != pst:
        res['mismatch_status'].append(dict(state=r, fortran=st, python=pst)); continue
    if st != 0:
        res['both_refuse'] += 1; continue
    res['both_ok'] += 1
    e0 = abs(float(c0f) - c0p) / abs(c0p); ev = abs(float(vf) - vp) / max(abs(vp), 1e-300)
    if max(e0, ev) > max(res['max_rel_C0'], res['max_rel_Vgj']):
        res['worst'] = dict(state=r, C0=(float(c0f), c0p), Vgj=(float(vf), vp))
    res['max_rel_C0'] = max(res['max_rel_C0'], e0); res['max_rel_Vgj'] = max(res['max_rel_Vgj'], ev)
json.dump(res, open(sys.argv[2], 'w'), indent=1)
print({k: v for k, v in res.items() if k != 'mismatch_status'}, 'status mismatches:', len(res['mismatch_status']))
for m in res['mismatch_status'][:10]: print(m)
