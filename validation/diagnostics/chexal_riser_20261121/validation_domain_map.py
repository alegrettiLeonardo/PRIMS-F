#!/usr/bin/env python3
"""P1.1: where do the ACCEPTED Chexal riser states of PRIMS-F sit relative to the
air-water data base of Chexal et al. (1992), Table 5-1 (p.280) and Table 5-2 (p.282)?
No void-fraction data points are available locally (only the authors' summary
statistics), so this maps DOMAIN COVERAGE; it is not a validation.

Data families (ranges as printed):
  INCLINED  Beggs 1972, pipe D 25.4-38.1 mm, p 2.56-6.81 bar, G 60-5500 kg/m2s,
            air mass fraction x 0.0008-0.833, 621 points, mean error .003, std .038 (alpha)
  HORIZONTAL (theta=0; Fukano also 10 deg) D 12.7-127 mm, p 1.03-2.25 bar, G 54-4700, x 0.0001-0.447
  UPFLOW (vertical)  D 11-305 mm, p 1.01-4 bar, G 0.18-4060, x 0.00004-1.0
  CCFL (vertical countercurrent)  Bharathan 1979, tubes 25 and 51 mm, 1.01 bar
Classes per state:
  VALIDATED     inside every range of the family of its orientation;
  INTERPOLATED  outside only in diameter, with D inside the span of the families that
                bracket its angle (horizontal 12.7-127 mm and vertical 11-305 mm),
                i.e. the Fr interpolation of Eq.3-2 between covered end points;
  EXTRAPOLATED  p, G or x outside every family bracketing its angle;
  UNSUPPORTED   countercurrent at an inclined angle, horizontal negative flux, gas down/liquid up.
Usage: validation_domain_map.py OUT_JSON OUT_MD RUN_DIR [RUN_DIR ...]
"""
import json, math, re, sys
from pathlib import Path
import numpy as np

FAM = dict(INCLINED=dict(D=(25.4e-3, 38.1e-3), p=(2.56e5, 6.81e5), G=(60, 5500), x=(0.0008, 0.833)),
           HORIZONTAL=dict(D=(12.7e-3, 127e-3), p=(1.03e5, 2.25e5), G=(54, 4700), x=(0.0001, 0.447)),
           UPFLOW=dict(D=(11e-3, 305e-3), p=(1.01e5, 4.0e5), G=(0.18, 4060), x=(0.00004, 1.0)))


def inside(v, r):
    return r[0] <= v <= r[1]


def classify(D, p, G, x, theta, jl, jg):
    if (theta > 0 and jl < 0 <= jg) or (theta == 0 and (jl < 0 or jg < 0)) or (jg < 0 < jl):
        return 'UNSUPPORTED', []
    fam = 'HORIZONTAL' if theta == 0 else 'INCLINED'
    f = FAM[fam]
    out = [k for k, v in (('D', D), ('p', p), ('G', G), ('x', x)) if not inside(v, f[k])]
    if not out:
        return 'VALIDATED', []
    if out == ['D'] and fam == 'INCLINED' and inside(D, FAM['HORIZONTAL']['D']) and inside(D, FAM['UPFLOW']['D']):
        return 'INTERPOLATED', out
    if fam == 'INCLINED':
        # p, G, x covered by at least one bracketing family (horizontal or vertical)?
        cov = all(any(inside(v, FAM[g][k]) for g in ('INCLINED', 'HORIZONTAL', 'UPFLOW')) for k, v in (('p', p), ('G', G), ('x', x)))
        if cov and inside(D, FAM['HORIZONTAL']['D']) and inside(D, FAM['UPFLOW']['D']):
            return 'INTERPOLATED', out
    return 'EXTRAPOLATED', out


def main():
    res = []
    for r in sys.argv[3:]:
        run = Path(r)
        deck = next(run.parent.parent.glob('decks/' + run.name + '.nml')).read_text()
        body = '\n'.join(l for l in deck.splitlines() if not l.lstrip().startswith('!'))
        D = float(re.search(r'\bD\s*=\s*([0-9.eE+-]+)', body).group(1))
        lp = float(re.search(r'\bLp\s*=\s*([0-9.eE+-]+)', body).group(1))
        f = np.genfromtxt(run / 'fields.csv', delimiter=',', names=True)
        rows = f[(f['s_m'] > lp) & (f['theta_rad'] >= 0)]
        counts, reasons = {}, {}
        ranges = dict(p=[np.inf, -np.inf], G=[np.inf, -np.inf], x=[np.inf, -np.inf], theta_deg=[np.inf, -np.inf])
        for q in rows:
            a = float(q['alpha']); jl = (1 - a) * float(q['ul_m_s']); jg = a * float(q['ug_m_s'])
            G = float(q['rhol_kg_m3']) * jl + float(q['rhog_kg_m3']) * jg
            x = float(q['rhog_kg_m3']) * jg / G if G > 0 else math.nan
            th = float(q['theta_rad']); p = float(q['p_Pa'])
            c, why = classify(D, p, abs(G), x, th, jl, jg)
            counts[c] = counts.get(c, 0) + 1
            for w in why:
                reasons[w] = reasons.get(w, 0) + 1
            for k, v in (('p', p), ('G', G), ('x', x), ('theta_deg', math.degrees(th))):
                ranges[k] = [min(ranges[k][0], v), max(ranges[k][1], v)]
        n = int(rows.size)
        res.append(dict(run=run.name, states=n, fractions={k: v / n for k, v in counts.items()}, counts=counts,
                        out_of_range_dimension_counts=reasons, ranges=ranges, D_m=D))
    json.dump(dict(families=FAM, runs=res), open(sys.argv[1], 'w'), indent=1)
    with open(sys.argv[2], 'w') as m:
        m.write('| run | states | VALIDATED | INTERPOLATED | EXTRAPOLATED | UNSUPPORTED | out-of-range dims | p bar | G kg/m2s | x | theta deg |\n'
                '|---|---:|---:|---:|---:|---:|---|---|---|---|---|\n')
        for r in res:
            fr = r['fractions']; rg = r['ranges']
            m.write('| %s | %d | %.3f | %.3f | %.3f | %.3f | %s | %.2f-%.2f | %.0f-%.0f | %.2g-%.2g | %.1f-%.1f |\n' % (
                r['run'], r['states'], fr.get('VALIDATED', 0), fr.get('INTERPOLATED', 0), fr.get('EXTRAPOLATED', 0), fr.get('UNSUPPORTED', 0),
                r['out_of_range_dimension_counts'], rg['p'][0] / 1e5, rg['p'][1] / 1e5, rg['G'][0], rg['G'][1], rg['x'][0], rg['x'][1],
                rg['theta_deg'][0], rg['theta_deg'][1]))
    print(open(sys.argv[2]).read())


if __name__ == '__main__':
    main()
