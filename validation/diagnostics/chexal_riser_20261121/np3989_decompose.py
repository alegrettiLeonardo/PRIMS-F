#!/usr/bin/env python3
"""P0.1: decompose NP-3989 Table A-1 rows into the intermediate coefficients and
locate the FIRST divergence (diagnostic; nothing is fitted or changed).

The table prints Re_f, Re_g (inputs, matched exactly), C0, Vgj and alpha. All
intermediates are evaluated at the PRINTED alpha with the 1986 equations (the
table's own version; upflow identical to 1992):
  C0 chain:  B1 -> K0 = B1 + (1-B1)(rg/rf)^.25 -> r -> L(alpha,p) -> C0
  Vgj chain: prop = [(rf-rg) sigma g / rf^2]^.25 -> C2 -> C3 -> C4 -> (1-alpha)^K1 -> Vgj
Implied quantities from the printed values (the only unknown of each chain
isolated, everything else as computed):
  (rg/rf)_implied   from printed C0 (solving K0 for the density ratio),
  prop_implied      = printed Vgj / (1.41 C2 C3 C4 (1-alpha)^K1).
The first printed quantity that diverges, and which intermediate explains it,
is reported per row; rows at the same pressure must imply the same properties
if properties are the cause.
Usage: np3989_decompose.py OUT_JSON OUT_MD
"""
import json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import np3989_audit as au


def chain(c, ratio=None):
    n, P, Dft, jf, jg, Ref, Reg, C0t, Vt, at = c
    pr = au.PROPS[P]; rf, rg = pr['rf'], pr['rg']
    if ratio is not None:
        rg = rf * ratio
    sig = au.sigma(pr['T']); D = Dft * au.FT
    re = Reg if (Reg > Ref or Reg < 0) else Ref
    b1 = min(0.8, 1 / (1 + math.exp(-re / 60000)))
    k0 = b1 + (1 - b1) * (rg / rf) ** 0.25
    r = (1 + 1.57 * rg / rf) / (1 - b1)
    c1 = 4 * au.PCRIT ** 2 / (P * (au.PCRIT - P)); L = (1 - math.exp(-c1 * at)) / (1 - math.exp(-c1))
    c0 = L / (k0 + (1 - k0) * at ** r)
    c5 = math.sqrt(150 * rg / rf); c2 = 1.0 if c5 >= 1 else 1 / (1 - math.exp(-c5 / (1 - c5)))
    c7 = (0.09144 / D) ** 0.6; c4 = 1.0 if c7 >= 1 else 1 / (1 - math.exp(-c7 / (1 - c7)))
    q = abs(Ref)
    if jf >= 0 and jg >= 0:
        c3 = max(0.5, 2 * math.exp(-q / 60000))
    else:
        c3 = 2 * math.exp((q / 350000) ** 0.4) - 1.75 * q ** 0.03 * math.exp(-q / 50000 * (0.0381 / D) ** 2) + (0.0381 / D) ** 0.25 * q ** 0.001
    k1 = b1 if Reg >= 0 else min(0.65, 0.5 * math.exp(abs(Reg) / 4000))
    prop = ((rf - rg) * sig * au.G / rf ** 2) ** 0.25
    c9 = (1 - at) ** k1
    vgj = 1.41 * prop * c2 * c3 * c4 * c9 / au.FT
    return dict(B1=b1, K0=k0, r=r, L=L, C0=c0, C2=c2, C3=c3, C4=c4, K1=k1, C9=c9, prop=prop, Vgj_fts=vgj,
                identity=at * (c0 * (jf + jg) + vgj) - jg)


def implied_ratio(c):
    lo, hi = 1e-4, 0.2
    for _ in range(200):
        m = math.sqrt(lo * hi)
        lo, hi = (m, hi) if chain(c, m)['C0'] > c[7] else (lo, m)   # C0 decreases with the ratio
    return math.sqrt(lo * hi)


def main():
    rows = []
    for c in au.TABLE:
        ch = chain(c)
        e0 = (ch['C0'] - c[7]) / c[7]; ev = (ch['Vgj_fts'] - c[8]) / c[8]
        prop_imp = c[8] * au.FT / (1.41 * ch['C2'] * ch['C3'] * ch['C4'] * ch['C9'])
        pr = au.PROPS[c[1]]
        row = dict(case=c[0], P_psia=c[1], D_ft=c[2], jf=c[3], chain_at_printed_alpha=ch, C0_rel_err=e0, Vgj_rel_err=ev,
                   ratio_tables=pr['rg'] / pr['rf'], ratio_implied_by_C0=implied_ratio(c), prop_computed=ch['prop'],
                   prop_implied_by_Vgj=prop_imp)
        row['first_divergence'] = ('C0 (density ratio in K0, r)' if abs(e0) > 2e-4 else
                                   'Vgj (property group / C2-C4 factor)' if abs(ev) > 2e-3 else 'none within print/property band')
        rows.append(row)
    json.dump(rows, open(sys.argv[1], 'w'), indent=1)
    with open(sys.argv[2], 'w') as m:
        m.write('| case | P psia | D ft | dir | C0 err | Vgj err | rg/rf tables | rg/rf implied by C0 | prop computed | prop implied by Vgj | C3 | C4 | C9 | first divergence |\n')
        m.write('|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n')
        for r in rows:
            ch = r['chain_at_printed_alpha']
            m.write('| %d | %g | %g | %s | %+.2e | %+.2e | %.5f | %.5f | %.5f | %.5f | %.4f | %.4f | %.4f | %s |\n' % (
                r['case'], r['P_psia'], r['D_ft'], 'up' if r['jf'] > 0 else 'down', r['C0_rel_err'], r['Vgj_rel_err'],
                r['ratio_tables'], r['ratio_implied_by_C0'], r['prop_computed'], r['prop_implied_by_Vgj'], ch['C3'], ch['C4'], ch['C9'],
                r['first_divergence']))
    print(open(sys.argv[2]).read())


if __name__ == '__main__':
    main()
