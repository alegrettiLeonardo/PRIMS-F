#!/usr/bin/env python3
"""P1.2 local constitutive study at the ACTUAL blocking states (no solver coupling).

States: the boundary cell of every Chexal run classified BLOCKED_CONSTITUTIVE_DOMAIN
(evidence/DOMAIN_PROBE.json: alpha, p, theta of the accepted state with j_l -> 0+).
For a prescribed signed liquid flux j_l in {0, -1e-3, -1e-2, -3e-2} m/s (liquid falling,
gas rising = inclined countercurrent) the gas velocity u_g = j_g/alpha that satisfies
each candidate drift relation u_g = C0 j + Vgj (j = j_l + j_g) is solved locally:
  CHEXAL_AT_THETA   Chexal 1992 at the real angle: admissible ONLY at j_l = 0 (p.262);
                    for j_l < 0 the core refuses (reported as n/a);
  CHEXAL_VERTICAL_HIGH  Chexal 1992 vertical countercurrent, HIGH branch (Eq.3-40),
                    evaluated at theta_C = 0: SENSITIVITY ONLY, not admissible at theta;
  BENDIKSEN_1984    C0 = 1.05 + 0.15 sin(theta), Vgj = sqrt(gD)(0.35 sin + 0.54 cos)
                    (repository closure 2, Fr<3.5 branch): long-bubble data include
                    stagnant liquid; negative liquid flux is a continuation of the formula;
  BHAGWAT_EDGE      repository Bhagwat transcription is real only for j_l >= 0 (edge only).
Output: u_g, j_g per candidate and the spread max/min of u_g at each j_l: the
epistemic gap that a qualified inclined-countercurrent closure would have to close.
Usage: inclined_cc_local.py DOMAIN_PROBE_JSON OUT_JSON OUT_MD
"""
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'junction_constitutive'))
import chexal92_ref as cx
import ss3_countercurrent_local as bh

D, G, ML, MG, SIG = 0.05254, 9.81, 8.9e-4, 1.85e-5, 0.0728


def rho(p, eos):
    rl = 997.05 + (p - 101325.0) / 1498.0 ** 2
    rg = p / 343.0 ** 2 if eos == 'BAROTROPIC' else p / (287.0 * 293.15)
    return rl, rg


def solve_jg(f, lo=0.0, hi=50.0):
    """bisection on R(jg) = alpha u_g(C0,Vgj) - jg; None when no sign change / refused"""
    try:
        flo, fhi = f(lo), f(hi)
    except cx.DomainError:
        return None
    if flo * fhi > 0:
        return None
    for _ in range(200):
        m = .5 * (lo + hi)
        fm = f(m)
        lo, flo, hi = (m, fm, hi) if flo * fm > 0 else (lo, flo, m)
    return .5 * (lo + hi)


def main():
    probes = [p for p in json.load(open(sys.argv[1])) if 'min_liquid_flux_cell' in p and 'DOMAIN' in (p.get('decoded') or '')]
    out = []
    for pr in probes:
        c = pr['min_liquid_flux_cell']; a, th, p = c['alpha'], math.radians(c['theta_deg']), c['p_Pa']
        eos = 'ISOTHERMAL' if 'ISOTHERMAL' in pr['run'] else 'BAROTROPIC'
        rl, rg = rho(p, eos)
        rows = []
        for jl in (0.0, -1e-3, -1e-2, -3e-2):
            r = dict(jl=jl)
            # Chexal at the real angle
            try:
                jg = solve_jg(lambda q: cx.residual(a, rl, rg, th, jl, q, d=D, g=G, ml=ML, mg=MG, sig=SIG))
            except cx.DomainError:
                jg = None
            r['CHEXAL_AT_THETA'] = jg / a if jg else None
            # Chexal vertical HIGH (sensitivity only)
            jg = solve_jg(lambda q: cx.residual(a, rl, rg, math.pi / 2, jl, q, d=D, g=G, ml=ML, mg=MG, sig=SIG))
            r['CHEXAL_VERTICAL_HIGH'] = jg / a if jg else None
            # Bendiksen (linear): u_g = C0 (jl + a u_g) + Vgj
            c0 = 1.05 + 0.15 * math.sin(th); vgj = math.sqrt(G * D) * (0.35 * math.sin(th) + 0.54 * math.cos(th))
            r['BENDIKSEN_1984'] = (c0 * jl + vgj) / (1 - c0 * a)
            # Bhagwat: edge only
            if jl == 0.0:
                s = dict(bh.S, alpha=a, p=p, rho_l=rl, rho_g=rg, theta_deg=math.degrees(th))
                j = 1.0
                for _ in range(50):
                    c0b, vb, re, w = bh.bhagwat_j0(s, j)
                    ug = vb / (1 - c0b * a); j = a * ug
                r['BHAGWAT_EDGE'] = ug
            else:
                r['BHAGWAT_EDGE'] = None
            vals = [v for k, v in r.items() if k not in ('jl',) and v is not None]
            r['spread_max_over_min'] = max(vals) / min(vals) if len(vals) > 1 and min(vals) > 0 else None
            rows.append(r)
        out.append(dict(run=pr['run'], alpha=a, theta_deg=c['theta_deg'], p_Pa=p, eos=eos, rows=rows))
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
    f = lambda v: '-' if v is None else '%.4f' % v
    with open(sys.argv[3], 'w') as m:
        m.write('| run | alpha | theta deg | j_l m/s | u_g Chexal(theta) | u_g Chexal vertical HIGH (sensitivity) | u_g Bendiksen | u_g Bhagwat (edge) | spread |\n'
                '|---|---:|---:|---:|---:|---:|---:|---:|---:|\n')
        for o in out:
            for r in o['rows']:
                m.write('| %s | %.4f | %.2f | %g | %s | %s | %s | %s | %s |\n' % (o['run'], o['alpha'], o['theta_deg'], r['jl'],
                        f(r['CHEXAL_AT_THETA']), f(r['CHEXAL_VERTICAL_HIGH']), f(r['BENDIKSEN_1984']), f(r['BHAGWAT_EDGE']),
                        f(r['spread_max_over_min'])))
    print(open(sys.argv[3]).read())


if __name__ == '__main__':
    main()
