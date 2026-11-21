#!/usr/bin/env python3
"""Diagnostic (NOT a correction of the correlation): which single saturated
property set (rho_f, rho_g, sigma) at 1000 psia reproduces the 1986 Table A-1
rows 2, 6, 7 (9 printed values) with the 1986 equations? Row 3 (internally
inconsistent) is predicted out of sample. Viscosities stay those implied by the
printed Re. Least squares on relative errors; result compared with the classical
steam tables (vf 0.02159, vg 0.4456 ft3/lb) and IAPWS sigma.
Usage: np3989_property_fit.py OUT_JSON"""
import json, math, sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, __import__('os').path.dirname(__file__))
import np3989_audit as au
rows = [c for c in au.TABLE if c[1] == 1000]
def model(p, c, ver=1986):
    rf, rg, sig = p
    n, P, Dft, jf, jg, Ref, Reg = c[:7]
    D, jfs, jgs = Dft * au.FT, jf * au.FT, jg * au.FT
    f = lambda a: (lambda q: a * (q[0] * (jfs + jgs) + q[1]) - jgs)(au.coeffs(ver, a, P, D, jfs, jgs, rf, rg, sig, Ref, Reg))
    lo, hi = 0.3, 0.99   # each row has a single root in this interval (checked by the audit scan)
    flo = f(lo)
    for _ in range(80):
        m = .5 * (lo + hi); fm = f(m)
        lo, hi, flo = (m, hi, fm) if flo * fm > 0 else (lo, m, flo)
    a = .5 * (lo + hi)
    c0, v = au.coeffs(ver, a, P, D, jfs, jgs, rf, rg, sig, Ref, Reg)
    return np.array([c0, v / au.FT, a])
fit_rows = [c for c in rows if c[0] in (2, 6, 7)]
res = lambda p: np.concatenate([(model(p, c) - np.array(c[7:10])) / np.array(c[7:10]) for c in fit_rows])
p0 = np.array([au.PROPS[1000]['rf'], au.PROPS[1000]['rg'], au.sigma(au.PROPS[1000]['T'])])
sol = least_squares(lambda x: res(x * p0), np.ones(3), x_scale=np.ones(3), xtol=1e-14, ftol=1e-14)
p = sol.x * p0
out = dict(nominal=dict(rho_f=p0[0], rho_g=p0[1], sigma=p0[2], ratio=p0[1] / p0[0]),
           fitted=dict(rho_f=p[0], rho_g=p[1], sigma=p[2], ratio=p[1] / p[0]),
           relative_change=dict(rho_f=p[0] / p0[0] - 1, rho_g=p[1] / p0[1] - 1, sigma=p[2] / p0[2] - 1),
           rows={})
for c in rows:
    m = model(p, c); m0 = model(p0, c)
    out['rows'][c[0]] = dict(table=list(c[7:10]), fitted_props=m.tolist(), nominal_props=m0.tolist(),
                             max_rel_err_fitted=float(np.max(np.abs(m - np.array(c[7:10])) / np.array(c[7:10]))),
                             max_rel_err_nominal=float(np.max(np.abs(m0 - np.array(c[7:10])) / np.array(c[7:10]))),
                             in_fit=c[0] in (2, 6, 7))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
print(json.dumps(out, indent=1))
