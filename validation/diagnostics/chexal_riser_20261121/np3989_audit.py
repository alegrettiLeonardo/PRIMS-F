#!/usr/bin/env python3
"""P0.1 audit of EPRI NP-3989-SR Rev.1 (1986) Table A-1 against BOTH versions.

Versions (independent transcriptions, vertical, steam-water fluid parameter):
  V1986  NP-3989 pp.3-1..3-4 (PDF 25-28): C0 = L/[K0+(1-K0)a^r] for all
         directions; Vgj = 1.41[(rf-rg) s g/rf^2]^.25 (1-a)^K1 C2 C3 C4,
         K1 = B1 (Reg>=0) or min(0.65, 0.5 exp(|Reg|/4000)) (Reg<0);
         C2 = 1 / 1/(1-exp(-C5/(1-C5))) with C5 = sqrt(150 rg/rf) (no 0.4757 branch);
         C3 up = max(0.5, 2exp(-|Ref|/60000)); C3 down = C3' =
         2exp[(|Ref|/350000)^0.4] - 1.75|Ref|^0.03 exp(-|Ref|/50000 (D1/D)^2) + (D1/D)^.25 |Ref|^.001.
  V1992  Chexal et al. 1992 eqs.3-4a/b, 3-5/3-6, 3-11..3-17, 3-23..3-38:
         C0 for Reg<0 = max(3-4a, Vgj0 (1-a)^0.2/(|jf|+|jg|)); C9 = (1-a)^B1 or
         min(0.7,(1-a)^0.65); C2 has the rf/rg<=18 branch; C3 down = 2(C10/2)^B2.
  Steam-water L (both): (1-exp(-C1 a))/(1-exp(-C1)), C1 = 4 Pc^2/(p(Pc-p)), p, Pc in psia (absolute).
Properties (declared; the report prints none): saturated densities from the classical
steam tables (Keenan et al.) 14.7 psia: vf 0.01672, vg 26.80 ft3/lb; 1000 psia: vf 0.02159,
vg 0.4456 ft3/lb; sigma IAPWS (1994) at Tsat (212.0 F; 544.6 F); viscosities from the PRINTED
Re of each row (nu = |j| D/|Re|), so Re_f and Re_g are exactly the printed ones.
Each printed row is also checked for internal consistency: alpha (C0 j + Vgj) = jg with the
printed C0, Vgj, alpha (rounding bound from 4 printed digits).
Diagnostic fits (NOT used to change the correlation): the effective Reynolds number Re*
in B1/C3 and the property ratio rg/rf that would reproduce a row.
Usage: np3989_audit.py OUT_JSON OUT_MD
"""
import json, math, sys

FT, LBFT3 = 0.3048, 16.01846337
G = 9.80665
TABLE = [(1, 14.7, 0.05, 5, 10, 78684, 2297, 1.2037, 0.5979, 0.5361),
         (2, 1000, 0.05, 5, 10, 184219, 87257, 1.1116, 0.1410, 0.5947),
         (3, 1000, 1.0, 5, 10, 3684867, 1745362, 1.1119, 0.2234, 0.5914),
         (4, 14.7, 1.0, 5, 10, 1573883, 45947, 1.1922, 0.9054, 0.5323),
         (5, 14.7, 1.0, -5, -10, -1573883, -45947, 1.0400, 5.1172, 0.9538),
         (6, 1000, 1.0, -5, -10, -3684867, -1745362, 1.0498, 4.8444, 0.9171),
         (7, 1000, 0.05, -5, -10, -184219, -87257, 1.2948, 1.8775, 0.5700),
         (8, 14.7, 0.05, -5, -10, -78684, -2297, 1.3036, 4.6590, 0.6714)]
PROPS = {14.7: dict(rf=LBFT3 / 0.01672, rg=LBFT3 / 26.80, T=373.12),
         1000: dict(rf=LBFT3 / 0.02159, rg=LBFT3 / 0.4456, T=557.93)}
PCRIT = 3208.2   # psia (NP-3989 era); 22.064 MPa = 3200.1 psia tested in the band


def sigma(T):
    tau = 1 - T / 647.096
    return 0.2358 * tau ** 1.256 * (1 - 0.625 * tau)


def coeffs(ver, a, P, D, jf, jg, rf, rg, sig, ref, reg, pc=PCRIT):
    re = reg if (reg > ref or reg < 0) else ref
    b1 = min(0.8, 1 / (1 + math.exp(-re / 60000)))
    k0 = b1 + (1 - b1) * (rg / rf) ** 0.25
    r = (1 + 1.57 * rg / rf) / (1 - b1)
    c1 = 4 * pc ** 2 / (P * (pc - P))
    L = (1 - math.exp(-c1 * a)) / (1 - math.exp(-c1))
    c0 = L / (k0 + (1 - k0) * a ** r)
    c5 = math.sqrt(150 * rg / rf)
    if ver == 1992 and rf / rg <= 18:
        c2 = 0.4757 * math.log(rf / rg) ** 0.7
    else:
        c2 = 1.0 if c5 >= 1 else 1 / (1 - math.exp(-c5 / (1 - c5)))
    c7 = (0.09144 / D) ** 0.6
    c4 = 1.0 if c7 >= 1 else 1 / (1 - math.exp(-c7 / (1 - c7)))
    q = abs(ref)
    if jf >= 0 and jg >= 0:
        c3 = max(0.5, 2 * math.exp(-q / 60000))
    else:
        c10 = 2 * math.exp((q / 350000) ** 0.4) - 1.75 * q ** 0.03 * math.exp(-q / 50000 * (0.0381 / D) ** 2) \
            + (0.0381 / D) ** 0.25 * q ** 0.001
        c3 = c10 if ver == 1986 else 2 * (c10 / 2) ** ((1 / (1 + 0.05 * q / 350000)) ** 0.4)
    v0 = 1.41 * ((rf - rg) * sig * G / rf ** 2) ** 0.25 * c2 * c3 * c4
    if ver == 1986:
        k1 = b1 if reg >= 0 else min(0.65, 0.5 * math.exp(abs(reg) / 4000))
        return c0, v0 * (1 - a) ** k1
    if reg >= 0:
        return c0, v0 * (1 - a) ** b1
    c0 = max(c0, v0 * (1 - a) ** 0.2 / (abs(jf) + abs(jg)))
    return c0, v0 * min(0.7, (1 - a) ** 0.65)


def solve(ver, P, D, jf, jg, rf, rg, sig, ref, reg, pc=PCRIT):
    f = lambda a: (lambda c: a * (c[0] * (jf + jg) + c[1]) - jg)(coeffs(ver, a, P, D, jf, jg, rf, rg, sig, ref, reg, pc))
    xs = [i / 20000 for i in range(1, 20000)]
    roots = []
    for x0, x1 in zip(xs, xs[1:]):
        f0, f1 = f(x0), f(x1)
        if f0 == 0 or f0 * f1 < 0:
            lo, hi = x0, x1
            for _ in range(100):
                m = .5 * (lo + hi)
                lo, hi = (m, hi) if f(lo) * f(m) > 0 else (lo, m)
            roots.append(.5 * (lo + hi))
    return roots


def evaluate(ver, c, rscale=1.0, sscale=1.0, pc=PCRIT, re_override=None):
    n, P, Dft, jf, jg, Ref, Reg = c[:7]
    pr = PROPS[P]
    rf, rg = pr['rf'] * rscale, pr['rg'] / rscale
    sig = sigma(pr['T']) * sscale
    if re_override:
        Ref = re_override
    D, jfs, jgs = Dft * FT, jf * FT, jg * FT
    roots = solve(ver, P, D, jfs, jgs, rf, rg, sig, Ref, Reg, pc)
    if len(roots) != 1:
        return dict(nroots=len(roots))
    a = roots[0]
    c0, v = coeffs(ver, a, P, D, jfs, jgs, rf, rg, sig, Ref, Reg, pc)
    return dict(nroots=1, C0=c0, Vgj=v / FT, alpha=a)


def consistency(c):
    """printed alpha (C0 j + Vgj) - jg, relative to jg, and its 4-digit rounding bound"""
    jf, jg, C0, V, a = c[3], c[4], c[7], c[8], c[9]
    j = jf + jg
    val = a * (C0 * j + V) - jg
    bound = 5e-5 * abs(C0 * j + V) + a * 5e-5 * abs(j) + a * 5e-5
    return val / abs(jg), bound / abs(jg)


def main():
    rows = []
    for c in TABLE:
        ref = dict(C0=c[7], Vgj=c[8], alpha=c[9])
        cons, bound = consistency(c)
        row = dict(case=c[0], P_psia=c[1], D_ft=c[2], jf=c[3], jg=c[4], table=ref,
                   table_identity_rel=cons, table_rounding_bound_rel=bound, table_row_consistent=abs(cons) <= bound)
        for ver in (1986, 1992):
            nom = evaluate(ver, c)
            band = [evaluate(ver, c, rs, ss, pc) for rs in (0.995, 1.005) for ss in (0.98, 1.02) for pc in (3200.1, 3208.2)]
            row[str(ver)] = dict(nominal=nom, rel_err={k: (nom[k] - ref[k]) / ref[k] for k in ref} if nom.get('nroots') == 1 else None,
                                 band={k: [min(b[k] for b in band if b.get('nroots') == 1), max(b[k] for b in band if b.get('nroots') == 1)]
                                       for k in ref})
            row[str(ver)]['inside_band'] = {k: row[str(ver)]['band'][k][0] - 5e-5 * abs(ref[k]) <= ref[k] <= row[str(ver)]['band'][k][1] + 5e-5 * abs(ref[k]) for k in ref}
        rows.append(row)
    # diagnostic fit for the 1000 psia upflow rows: effective Reynolds number in B1/C3
    for row, c in zip(rows, TABLE):
        if c[1] == 1000 and c[3] > 0:
            best = None
            for k in range(400):
                re = 30000 * 1.01 ** k
                e = evaluate(1986, c, re_override=re)
                if e.get('nroots') != 1:
                    continue
                err = max(abs(e[q] - row['table'][q]) / row['table'][q] for q in ('C0', 'Vgj', 'alpha'))
                if best is None or err < best[0]:
                    best = (err, re, e)
            row['diagnostic_effective_Ref'] = dict(Re=best[1], max_rel_err=best[0], values=best[2], printed_Ref=c[5])
    # classification (cause/status), decided from the computed numbers
    for row in rows:
        up = row['jf'] > 0
        v86, v92 = row['1986'], row['1992']
        in86, in92 = all(v86['inside_band'].values()), all(v92['inside_band'].values())
        if up and in86 and in92:
            row.update(version='both (upflow form identical)', cause='reproduced', status='PASS')
        elif not up and in86 and not in92:
            row.update(version='1986', cause='downflow revised in 1992 (C3=2(C10/2)^B2, C9, Eq.3-4b); 1986 equations reproduce the row',
                       status='PASS (version difference demonstrated)')
        elif not up and in86 and in92:
            row.update(version='both', cause='reproduced by both versions', status='PASS')
        elif not row['table_row_consistent']:
            row.update(version='both', cause='printed row violates alpha(C0 j+Vgj)=jg beyond its rounding bound (defective reference row)',
                       status='FAIL (reference inconsistent)')
        else:
            row.update(version='both', cause='unexplained: not a version difference (upflow form identical); see diagnostic fit',
                       status='FAIL')
    json.dump(rows, open(sys.argv[1], 'w'), indent=1)
    f = lambda d, k: '%.4f' % d['nominal'][k] if d['nominal'].get('nroots') == 1 else 'n/a'
    with open(sys.argv[2], 'w') as m:
        m.write('# NP-3989 Table A-1: 1986 and 1992 transcriptions (steam-water L, declared properties)\n\n')
        m.write('| case | P psia | D ft | jf, jg | table C0 / Vgj / alpha | identity residual (bound) | V1986 C0 / Vgj / alpha | V1992 C0 / Vgj / alpha | version | cause | status |\n')
        m.write('|---:|---:|---:|---|---|---|---|---|---|---|---|\n')
        for r in rows:
            t = r['table']
            m.write('| %d | %g | %g | %g, %g | %.4f / %.4f / %.4f | %+.1e (%.1e) | %s / %s / %s | %s / %s / %s | %s | %s | %s |\n' % (
                r['case'], r['P_psia'], r['D_ft'], r['jf'], r['jg'], t['C0'], t['Vgj'], t['alpha'], r['table_identity_rel'],
                r['table_rounding_bound_rel'], f(r['1986'], 'C0'), f(r['1986'], 'Vgj'), f(r['1986'], 'alpha'),
                f(r['1992'], 'C0'), f(r['1992'], 'Vgj'), f(r['1992'], 'alpha'), r['version'], r['cause'], r['status']))
        for r in rows:
            if 'diagnostic_effective_Ref' in r:
                d = r['diagnostic_effective_Ref']
                m.write('\nCase %d diagnostic (not a correction): an effective Re_f = %.0f in B1/C3 (printed %d) gives max rel. error %.1e '
                        '(C0 %.4f, Vgj %.4f, alpha %.4f).\n' % (r['case'], d['Re'], d['printed_Ref'], d['max_rel_err'],
                                                               d['values']['C0'], d['values']['Vgj'], d['values']['alpha']))
    print(open(sys.argv[2]).read())


if __name__ == '__main__':
    main()
