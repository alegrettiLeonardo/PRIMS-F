#!/usr/bin/env python3
"""P1.2 local constitutive solver for the SS3 critical cell (diagnostic only;
nothing here enters the global solver).

State (accepted SS3 state, cell 69, from evidence/certificates/SS3_closure_map.log):
  alpha, p, rho_l, rho_g, mixture mass flux G = rho_l j_l + rho_g j_g, theta.
Signed drift relation  u_g = C0 j + Vgj,  j = j_l + j_g,  j_g = alpha u_g.
With j_l = (G - rho_g alpha u_g)/rho_l the relation is LINEAR in u_g:
  u_g = (C0 G/rho_l + Vgj) / (1 - C0 alpha (1 - rho_g/rho_l)).
Candidate coefficients, each with its documented basis:
  BENDIKSEN  C0 = 1.05 + 0.15 sin(theta), Vgj = sqrt(gD)(0.35 sin + 0.54 cos)
             (repository closure, Fr<3.5 branch; long bubbles in inclined tubes,
             measurements include stagnant liquid; co-current correlation);
  BHAGWAT_J0 C0 and Vgj of the repository Bhagwat transcription evaluated at the
             STAGNANT-LIQUID edge j_l = 0 (x = 1, (1-x)^1.5 = 0), where it is
             still defined; for j_l < 0 the mass quality x > 1 and (1-x)^1.5 is
             not real: the correlation has NO counter-current continuation.
  Sensitivity band C0 in [0.95, 1.25] with both Vgj (no source claimed: it only
  bounds how much the counter-current root depends on C0).
Output: u_g, u_l, j_l, j_g, sign of the liquid flow, |j_l|/Vgj (distance from the
stagnant-liquid measurements), for each candidate.
"""
import json, math, sys

S = dict(alpha=0.638272, p=228440.2, rho_l=997.1066, rho_g=1.9417, G=0.77404, theta_deg=51.637,
         D=0.05254, g=9.81, sigma=0.0728, mu_l=0.89e-3, eps=4.6e-5)


def colebrook_y(re, rr):
    # Eq.118 of the transcription: y + 4 log10(rr/3.7 + 1.256 y/Re) = 0  (y = 1/sqrt(f))
    lo, hi = 1e-6, 100.0
    for _ in range(200):
        y = 0.5 * (lo + hi)
        h = y + 4 * math.log10(rr / 3.7 + 1.256 * y / re)
        lo, hi = (y, hi) if h < 0 else (lo, y)
    return 0.5 * (lo + hi)


def bhagwat_j0(s, j):
    """C0, Vgj of the transcription at j_l = 0 (beta = 1, x = 1), for mixture flux j."""
    th = math.radians(s['theta_deg']); r = s['rho_g'] / s['rho_l']; a = s['alpha']
    re = s['rho_l'] * j * s['D'] / s['mu_l']
    w = 1 / (1 + (1000 / re) ** 2) if re >= 1000 else (re / 1000) ** 2 / (1 + (re / 1000) ** 2)
    co = math.cos(th)
    bb = ((1 + r * r * co) / (1 + co)) ** (0.2 * (1 - a))
    cc = 0.0  # 0.2(1-sqrt r)((2.6-beta)^0.15 - 1/y)(1-x)^1.5 with x = 1
    C0 = (2 - r * r) * (1 - w) + (bb + cc) * w
    la = math.sqrt(s['sigma'] / (s['g'] * (s['rho_l'] - s['rho_g']))) / s['D']
    c3 = 1.0 if la >= 0.025 else (la / 0.025) ** 0.9
    Vgj = (0.35 * math.sin(th) + 0.45 * co) * math.sqrt(s['g'] * s['D'] * (1 - r)) * math.sqrt(1 - a) * c3
    return C0, Vgj, re, w


def root(s, C0, Vgj):
    a, rl, rg, G = s['alpha'], s['rho_l'], s['rho_g'], s['G']
    ug = (C0 * G / rl + Vgj) / (1 - C0 * a * (1 - rg / rl))
    jg = a * ug; jl = (G - rg * jg) / rl
    return dict(C0=C0, Vgj=Vgj, u_g=ug, j_g=jg, j_l=jl, u_l=jl / (1 - a), counter_current=jl < 0,
                jl_over_Vgj=abs(jl) / Vgj, mass_quality_x=rg * jg / G)


def main():
    s = S; th = math.radians(s['theta_deg'])
    out = dict(state=s, candidates={}, band=[])
    Vb = math.sqrt(s['g'] * s['D']) * (0.35 * math.sin(th) + 0.54 * math.cos(th))
    out['candidates']['BENDIKSEN'] = root(s, 1.05 + 0.15 * math.sin(th), Vb)
    j = 1.0  # fixed point in j for the Re weight of Bhagwat (weight ~1 for Re >> 1000)
    for _ in range(50):
        C0, Vg, re, w = bhagwat_j0(s, j)
        r = root(s, C0, Vg); j = r['j_g'] + r['j_l']
    r.update(Re=re, weight=w); out['candidates']['BHAGWAT_J0'] = r
    for C0 in (0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25):
        for name, V in (('Vgj_Bendiksen', Vb), ('Vgj_Bhagwat_j0', out['candidates']['BHAGWAT_J0']['Vgj'])):
            out['band'].append(dict(Vgj_source=name, **root(s, C0, V)))
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
    with open(sys.argv[2], 'w') as m:
        m.write('# P1.2 SS3 cell 69: signed drift relation near stagnant liquid (local solver)\n\n')
        m.write('State: alpha=%.6f, p=%.1f Pa, rho_l=%.4f, rho_g=%.4f, G=%.5f kg/m2s, theta=%.3f deg.\n\n' % (
            s['alpha'], s['p'], s['rho_l'], s['rho_g'], s['G'], s['theta_deg']))
        m.write('| candidate | C0 | Vgj m/s | u_g m/s | u_l mm/s | j_l mm/s | |j_l|/Vgj | x | counter-current |\n|---|---:|---:|---:|---:|---:|---:|---:|---|\n')
        for n, r in out['candidates'].items():
            m.write('| %s | %.4f | %.4f | %.4f | %.3f | %.3f | %.4f | %.3f | %s |\n' % (n, r['C0'], r['Vgj'], r['u_g'], 1e3 * r['u_l'],
                    1e3 * r['j_l'], r['jl_over_Vgj'], r['mass_quality_x'], r['counter_current']))
        m.write('\nSensitivity band (no source claimed for the C0 values):\n\n| Vgj source | C0 | u_g m/s | u_l mm/s | |j_l|/Vgj |\n|---|---:|---:|---:|---:|\n')
        for r in out['band']:
            m.write('| %s | %.2f | %.4f | %.3f | %.4f |\n' % (r['Vgj_source'], r['C0'], r['u_g'], 1e3 * r['u_l'], r['jl_over_Vgj']))
    print(open(sys.argv[2]).read())


if __name__ == '__main__':
    main()
