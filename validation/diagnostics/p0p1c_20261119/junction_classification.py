#!/usr/bin/env python3
"""P0.1: gas stagnation (A) vs phase disappearance (B) vs numerical null flux (C).

Input: junction_state_compare.f90 output (CELL/FACE records, SI) and the
run's fields.csv (initial snapshot gives the physical holdup scale).
Dimensionless indicators, each against a PHYSICAL scale of the same run:
  r_alpha = alpha_g / alpha_g(t=0, same cell)          (phase presence)
  r_u     = u_g / |U_d(state)|  (Bhagwat drift velocity of that state)
  r_G     = G_g,face / min(G_g of the two cells)  (face flux vs transport)
Classes for a threshold eps:
  B  if r_alpha < eps in either adjacent cell;
  A  if not B and |r_u| < eps in the upstream (pipe) cell;
  C  if not A/B, both cells carry gas upward (r_u >= eps) and r_G < eps;
  N  otherwise (the face transmits the gas transport).
No eps is chosen: the class is listed over eps in [1e-15, 1], so the reader
sees for which range of thresholds each verdict holds.
Usage: junction_classification.py COMPARE_OUT FIELDS_CSV LABEL [state 0|1] [flux] [trace]
"""
import sys
import numpy as np


def read(path, which, flux, trace):
    cells, face = {}, None
    for line in open(path):
        p = [x.strip() for x in line.split('|')]
        if p[0] == 'CELL' and int(p[1]) == which:
            v = [float(x) for x in p[5].split()]
            cells[int(p[2])] = dict(alpha=v[0], rho_g=v[1], u_g=v[2], j_g=v[3], U2=v[4], p=v[5], U3=v[6], Gg=v[7],
                                    theta=v[8], Ud=v[9], closure=int(p[3]))
        if p[0] == 'FACE' and int(p[1]) == which and p[2] == flux and p[3] == trace:
            v = [float(x) for x in p[5].split()]
            face = dict(status=int(p[4]), p=v[0], Gl=v[1], Gg=v[2], a_pipe=v[3], a_riser=v[4])
    return cells, face


def main():
    out, fields, label = sys.argv[1:4]
    which = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    flux = sys.argv[5] if len(sys.argv) > 5 else 'ROE'
    trace = sys.argv[6] if len(sys.argv) > 6 else 'CENTER'
    cells, face = read(out, which, flux, trace)
    left, right = sorted(cells)
    f = np.genfromtxt(fields, delimiter=',', names=True)
    t0 = f['time_s'] == f['time_s'].min()
    a0 = f['alpha'][t0]
    ra = {k: cells[k]['alpha'] / a0[k - 3] for k in cells}      # physical cell k -> row k-3
    ru = {k: cells[k]['u_g'] / abs(cells[k]['Ud']) for k in cells}
    rg = face['Gg'] / min(cells[left]['Gg'], cells[right]['Gg']) if face['status'] == 0 else float('nan')
    print('%s state=%d face=%s/%s status=%d' % (label, which, flux, trace, face['status']))
    for k in (left, right):
        c = cells[k]
        print('  cell %d: alpha=%.4g (r_alpha=%.3g) rho_g=%.4g u_g=%.4g m/s (|Ud|=%.3g, r_u=%.3g) j_g=%.4g G_g=%.4g p=%.1f' % (
            k, c['alpha'], ra[k], c['rho_g'], c['u_g'], abs(c['Ud']), ru[k], c['j_g'], c['Gg'], c['p']))
    print('  face: G_g=%.4g kg/m2s  r_G=%.3g  alpha*_pipe=%.4g alpha*_riser=%.4g' % (face['Gg'], rg, face['a_pipe'], face['a_riser']))
    lo_b, lo_a = min(ra.values()), abs(ru[left])
    up = min(ru.values())

    def klass(eps):
        if lo_b < eps:
            return 'B'
        if lo_a < eps:
            return 'A'
        if np.isfinite(rg) and rg < eps and up >= eps:
            return 'C'
        return 'N'   # face transmits the transport (no blockage at this eps)
    grid = 10.0 ** np.linspace(-15, 0, 1501)
    ks = [klass(e) for e in grid]
    spans, start = [], 0
    for i in range(1, len(grid) + 1):
        if i == len(grid) or ks[i] != ks[start]:
            spans.append((ks[start], grid[start], grid[i - 1]))
            start = i
    print('  classes vs eps: ' + '; '.join('%s for eps in [%.2g, %.2g]' % sp for sp in spans))
    print('  distance to A: r_u(pipe)=%.3g ; distance to B: min r_alpha=%.3g ; null-face ratio r_G=%.3g' % (lo_a, lo_b, rg))

if __name__ == '__main__':
    main()
