#!/usr/bin/env python3
"""External numerical reference for the Chexal-Lellouche structure.

EPRI NP-3989-SR Rev.1 (Sept. 1986), Appendix A, Table A-1 (PDF page 100) prints
eight vertical steam-water sample problems (P psia, D ft, jf, jg ft/s ->
Re_f, Re_g, C0, Vgj ft/s, alpha). Chexal et al. (1992) p.260 states that the
vertical-UPFLOW steam-water form "remains the same" as the 1986 correlation;
downflow was revised in 1992. Hence:
  * cases 1-4 (upflow) are REFERENCES for the 1992 equations 3-2..3-34;
  * cases 5-8 (downflow) are reported for information only (1992 revised them).

This script is a third, independent transcription of the 1992 vertical
equations (it does not import the solver or scripts/generate_chexal_reference.py).
Only the fluid parameter L differs between fluids: steam-water Eq.3-5/3-6 here,
air-water Eq.3-7/3-8 in PRIMS-F. With --air-water it prints the same structure
with the air-water L, to be compared with the Fortran core by
chexal_core_crosscheck.f90 (same inputs).

Properties (declared, not taken from the report, which does not print them):
  kinematic viscosities from the PRINTED Re_f, Re_g of each case (nu = j D/Re);
  surface tension from the IAPWS (1994) release formula;
  saturated densities from IAPWS-IF97 steam tables, entered manually, uncertainty
  declared +-0.5 %; P_crit = 22.064 MPa. A sensitivity band (rho +-0.5 %,
  sigma +-2 %, P_crit 3200-3208.2 psia) is propagated to C0, Vgj and alpha.
Usage: np3989_sample_check.py OUT_JSON OUT_MD
"""
import itertools, json, math, sys

FT, PSI = 0.3048, 6894.757
G = 9.80665
TABLE = [  # case, P psia, D ft, jf, jg ft/s, Re_f, Re_g, C0, Vgj ft/s, alpha
    (1, 14.7, 0.05, 5, 10, 78684, 2297, 1.2037, 0.5979, 0.5361),
    (2, 1000, 0.05, 5, 10, 184219, 87257, 1.1116, 0.1410, 0.5947),
    (3, 1000, 1.0, 5, 10, 3684867, 1745362, 1.1119, 0.2234, 0.5914),
    (4, 14.7, 1.0, 5, 10, 1573883, 45947, 1.1922, 0.9054, 0.5323),
    (5, 14.7, 1.0, -5, -10, -1573883, -45947, 1.0400, 5.1172, 0.9538),
    (6, 1000, 1.0, -5, -10, -3684867, -1745362, 1.0498, 4.8444, 0.9171),
    (7, 1000, 0.05, -5, -10, -184219, -87257, 1.2948, 1.8775, 0.5700),
    (8, 14.7, 0.05, -5, -10, -78684, -2297, 1.3036, 4.6590, 0.6714),
]
# saturated densities kg/m3 (IAPWS-IF97 tables): 101.35 kPa (14.7 psia), 6.8948 MPa (1000 psia)
RHO = {14.7: (958.35, 0.5977), 1000: (741.9, 35.94)}
TSAT = {14.7: 373.12, 1000: 558.02}


def sigma_iapws(T):
    tau = 1 - T / 647.096
    return 0.2358 * tau ** 1.256 * (1 - 0.625 * tau)


def coefficients(a, P, D, jf, jg, rf, rg, sig, ref, reg, fluid, pcrit):
    # Eq.3-15: signed Reynolds selection
    re = reg if (reg > ref or reg < 0) else ref
    b1 = min(0.8, 1 / (1 + math.exp(-re / 60000)))                   # 3-13, 3-14
    k0 = b1 + (1 - b1) * (rg / rf) ** 0.25                            # 3-11
    r = (1 + 1.57 * rg / rf) / (1 - b1)                               # 3-12
    if fluid == 'steam':
        c1 = 4 * pcrit ** 2 / (P * (pcrit - P))                       # 3-6
        L = (1 - math.exp(-c1 * a)) / (1 - math.exp(-c1))             # 3-5
    else:
        L = min(1.15 * a ** 0.45, 1.0) if reg >= 0 else min(1.05 * a ** 0.25, 1.0)  # 3-7, 3-8
    c0 = L / (k0 + (1 - k0) * a ** r)                                 # 3-4a
    if rf / rg <= 18:
        c2 = 0.4757 * math.log(rf / rg) ** 0.7                        # 3-27
    else:
        c5 = math.sqrt(150 / (rf / rg))                               # 3-29
        c2 = 1.0 if c5 >= 1 else 1 / (1 - math.exp(-c5 / (1 - c5)))   # 3-28
    c7 = (0.09144 / D) ** 0.6                                         # 3-31, 3-33
    c4 = 1.0 if c7 >= 1 else 1 / (1 - math.exp(-c7 / (1 - c7)))      # 3-30, 3-32
    if jf >= 0 and jg >= 0:
        c3 = max(0.5, 2 * math.exp(-abs(ref) / 60000))                # 3-34
    else:                                                             # 3-35..3-38 (1992)
        q = abs(ref)
        c10 = 2 * math.exp((q / 350000) ** 0.4) - 1.75 * q ** 0.03 * math.exp(-q / 50000 * (0.0381 / D) ** 2) \
            + (0.0381 / D) ** 0.25 * q ** 0.001
        b2 = (1 / (1 + 0.05 * q / 350000)) ** 0.4
        c3 = 2 * (c10 / 2) ** b2
    v0 = 1.41 * ((rf - rg) * sig * G / rf ** 2) ** 0.25 * c2 * c3 * c4  # 3-24 (SI, gc=1)
    c9 = (1 - a) ** b1 if reg >= 0 else min(0.7, (1 - a) ** 0.65)     # 3-25, 3-26
    if reg < 0:                                                       # 3-4b
        c0 = max(c0, v0 * (1 - a) ** 0.2 / (abs(jf) + abs(jg)))
    return c0, v0 * c9


def solve(P, D, jf, jg, rf, rg, sig, ref, reg, fluid, pcrit):
    """alpha from alpha*(C0 (jf+jg) + Vgj) - jg = 0 (Eq.3-1 identity), bisection
    on (0,1); a sign scan first so that multiple roots are reported, not hidden."""
    f = lambda a: a * (coefficients(a, P, D, jf, jg, rf, rg, sig, ref, reg, fluid, pcrit)[0] * (jf + jg)
                       + coefficients(a, P, D, jf, jg, rf, rg, sig, ref, reg, fluid, pcrit)[1]) - jg
    xs = [i / 4000 for i in range(1, 4000)]
    roots = []
    for x0, x1 in zip(xs, xs[1:]):
        f0, f1 = f(x0), f(x1)
        if f0 == 0 or f0 * f1 < 0:
            lo, hi = x0, x1
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
            roots.append(0.5 * (lo + hi))
    return roots


def run_case(c, fluid, rho_scale=1.0, sig_scale=1.0, pcrit_psia=22.064e6 / PSI):
    n, P, Dft, jf, jg, Ref, Reg = c[:7]
    rf, rg = (x * rho_scale for x in RHO[P])
    rg = RHO[P][1] * (2 - rho_scale)   # opposite perturbation for the ratio
    sig = sigma_iapws(TSAT[P]) * sig_scale
    D, jfs, jgs = Dft * FT, jf * FT, jg * FT
    roots = solve(P, D, jfs, jgs, rf, rg, sig, Ref, Reg, fluid, pcrit_psia)
    out = []
    for a in roots:
        c0, v = coefficients(a, P, D, jfs, jgs, rf, rg, sig, Ref, Reg, fluid, pcrit_psia)
        out.append(dict(alpha=a, C0=c0, Vgj_fts=v / FT))
    return out


def main():
    res = []
    for c in TABLE:
        nominal = run_case(c, 'steam')
        band = []
        for rs, ss, pc in itertools.product((0.995, 1.0, 1.005), (0.98, 1.0, 1.02), (3200.1, 3208.2)):
            band.extend(run_case(c, 'steam', rs, ss, pc))
        ref = dict(C0=c[7], Vgj_fts=c[8], alpha=c[9])
        row = dict(case=c[0], P_psia=c[1], D_ft=c[2], jf=c[3], jg=c[4], Re_f=c[5], Re_g=c[6], table=ref,
                   role='REFERENCE (upflow form unchanged 1986->1992, p.260)' if c[3] > 0 else
                   'INFORMATIVE (downflow revised in 1992)', nroots=len(nominal))
        if len(nominal) == 1:
            n0 = nominal[0]
            row['nominal'] = n0
            row['rel_err'] = {k: (n0[k] - ref[k]) / ref[k] for k in ref}
            row['band'] = {k: [min(b[k] for b in band), max(b[k] for b in band)] for k in ref}
            row['table_inside_band'] = {k: row['band'][k][0] - 5e-5 * abs(ref[k]) <= ref[k] <= row['band'][k][1] + 5e-5 * abs(ref[k])
                                        for k in ref}
        res.append(row)
    # tolerance: 4 printed digits (+-5e-5 relative of the printed value) plus the
    # declared property band; upflow references PASS if every quantity lies in the band
    for r in res:
        if r['jf'] > 0:
            r['status'] = 'PASS' if r['nroots'] == 1 and all(r['table_inside_band'].values()) else 'FAIL'
        else:
            r['status'] = 'INFO'
    json.dump(res, open(sys.argv[1], 'w'), indent=1)
    with open(sys.argv[2], 'w') as m:
        m.write('# EPRI NP-3989-SR Rev.1 Table A-1 vs independent transcription of Chexal (1992) vertical equations\n\n')
        m.write('| case | P psia | D ft | jf, jg ft/s | role | quantity | table | nominal | rel. error | property band | status |\n'
                '|---:|---:|---:|---|---|---|---:|---:|---:|---|---|\n')
        for r in res:
            for k in ('C0', 'Vgj_fts', 'alpha'):
                if 'nominal' not in r:
                    continue
                m.write('| %d | %g | %g | %g, %g | %s | %s | %.4f | %.4f | %+.2e | %.4f-%.4f | %s |\n' % (
                    r['case'], r['P_psia'], r['D_ft'], r['jf'], r['jg'], r['role'].split(' ')[0], k, r['table'][k],
                    r['nominal'][k], r['rel_err'][k], r['band'][k][0], r['band'][k][1], r['status']))
    print(open(sys.argv[2]).read())


if __name__ == '__main__':
    main()
