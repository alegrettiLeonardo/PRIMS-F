#!/usr/bin/env python3
"""Mesh-independent physical perturbation of a discrete stationary state (P0.3).

  p(s,0) = p_ss,h(s) + DeltaP * f(s),   f(s) = (1 - tanh((s - Lp)/w)) / 2

The stationary state p_ss,h is the DISCRETE steady state of the same operator
on the same mesh (checkpoint of a large-dt backward-Euler run), so the only
seed of the unstable mode is the declared perturbation (plus round-off).
The perturbation is thermodynamically admissible and EOS-consistent: at fixed
alpha and phase velocities, rho_l = rho_l0 + (p - p_l0)/c_l^2 and
rho_g = p/c_g^2 are affine in p, hence
  dU1 = (1-alpha) dp/c_l^2,  dU2 = alpha dp/c_g^2,  dU3 = dU1 u_l + dU2 u_g
is EXACTLY linear in dp; with the analytic cell mean of f the conservative
projection preserves the integral of dp exactly (no quadrature error).
Phase velocities are those of the stationary fields (fields.csv, last
snapshot). For DeltaP = 0 the written U equals the checkpoint U bitwise.

Usage: make_perturbed_ic.py STEADY_RUN_DIR DeltaP_Pa w_m c_gas_m_s OUT_U OUT_JSON
"""
import json
import math
import sys
import numpy as np

CL = 1498.0
OMEGA_RHO, OMEGA_P = 997.05, 101325.0


def logcosh(x):
    ax = abs(x)
    return ax + math.log1p(math.exp(-2 * ax)) - math.log(2.0)


def F(s, lp, w):  # primitive of f
    return 0.5 * (s - w * logcosh((s - lp) / w))


def main():
    run, dp, w, cg, out_u, out_json = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5], sys.argv[6]
    lines = open(run + '/checkpoint.dat').read().split('\n')
    head = lines[1].split()
    ncell, nd, npipe = int(head[0]), int(head[1]), int(head[2])
    omega_t, omega_u = float(head[7]), float(head[8])  # 5 integers, then t, dt, omega_t, omega_u
    lr = omega_t * omega_u
    rows = np.array([[float(v) for v in l.split()] for l in lines[2:2 + ncell]])
    u = rows[nd:ncell - nd, 4:7]
    width = rows[nd:ncell - nd, 3] * lr
    faces = np.concatenate([[0.0], np.cumsum(width)])
    lp = faces[npipe]
    f = np.genfromtxt(run + '/fields.csv', delimiter=',', names=True)
    # The terminal snapshot may be written twice by the finalization; take the
    # last complete mesh of the last time.
    idx = np.nonzero(f['time_s'] == f['time_s'].max())[0][-len(u):]
    alpha, ul, ug = f['alpha'][idx], f['ul_m_s'][idx] / omega_u, f['ug_m_s'][idx] / omega_u
    assert np.allclose(f['s_m'][idx], (faces[:-1] + faces[1:]) / 2, rtol=0, atol=1e-9), 'snapshot/mesh mismatch'
    assert len(alpha) == len(u), 'fields snapshot does not match the checkpoint mesh'
    ncl2 = (CL * math.sqrt(OMEGA_RHO / OMEGA_P)) ** 2
    ncg2 = (cg * math.sqrt(OMEGA_RHO / OMEGA_P)) ** 2
    fbar = np.array([(F(b, lp, w) - F(a, lp, w)) / (b - a) for a, b in zip(faces[:-1], faces[1:])])
    dpn = dp * fbar / OMEGA_P
    du = np.zeros_like(u)
    du[:, 0] = (1 - alpha) * dpn / ncl2
    du[:, 1] = alpha * dpn / ncg2
    du[:, 2] = du[:, 0] * ul + du[:, 1] * ug
    unew = u + du if dp != 0 else u.copy()
    with open(out_u, 'w') as h:
        for r in unew:
            h.write('%.17e %.17e %.17e\n' % tuple(r))
    # projection diagnostics: realized integral, L2 distance between the
    # piecewise-constant projection and the continuous shape (fine sampling)
    integral_exact = F(faces[-1], lp, w) - F(faces[0], lp, w)
    integral_proj = float((fbar * width).sum())
    xs = np.linspace(0, faces[-1], 200001)
    fc = 0.5 * (1 - np.tanh((xs - lp) / w))
    fp = fbar[np.clip(np.searchsorted(faces, xs, side='right') - 1, 0, len(fbar) - 1)]
    l2 = float(np.sqrt(np.trapz((fc - fp) ** 2, xs) / np.trapz(fc ** 2, xs)))
    gas_mass_rel = float((du[:, 1] * width).sum() / (u[:, 1] * width).sum())
    json.dump(dict(steady_run=run, DeltaP_Pa=dp, w_m=w, Lp_m=float(lp), c_gas_m_s=cg, cells=len(u),
                   integral_f_exact_m=integral_exact, integral_f_projected_m=integral_proj,
                   projection_L2_relative=l2, gas_inventory_relative_change=gas_mass_rel,
                   max_dp_cell_Pa=float(dp * fbar.max()),
                   definition=__doc__.split('Usage')[0]), open(out_json, 'w'), indent=1)
    print(json.dumps(dict(cells=len(u), Lp=lp, integral_rel_err=abs(integral_proj - integral_exact) / integral_exact,
                          projection_L2=l2, gas_rel=gas_mass_rel)))


if __name__ == '__main__':
    main()
