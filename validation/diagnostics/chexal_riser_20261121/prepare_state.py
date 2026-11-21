#!/usr/bin/env python3
"""Build the inputs of global_jacobian_check from an ACCEPTED saved state of a
run: fields.csv row block nearest to TIME (primitives alpha, rho_l, rho_g, u_l,
u_g, theta) -> nondimensional conservative U per physical cell, exactly as
module_solver_utils.compute_conservative_variables, with the deck's scales
(omega_rho=rho_l0, omega_u=max(jl,jg), omega_P=P_l0). Mesh: mesh_ratio=1 gives
uniform sections (checked against the saved cell centres).
Usage: prepare_state.py RUN_DIR TIME OUT_PREFIX [DT_S]
       TIME = DUMP_OLD | DUMP_TRIAL: take Uold / refused trial from last_rejected_state.dat
"""
import json, math, re, sys
from pathlib import Path
import numpy as np

LR = 12.546253993197080   # catenary arc length of the D geometry (X=6.435, Z=9.886), as offline_api


def deck_value(text, key, cast=float):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))
    m = re.search(r'\b' + key + r'\s*=\s*([^,\n]+)', text, re.I)
    return cast(m.group(1).strip().strip("'"))


def write_params(out, ou, deck, jl, jg, cg, dt_s, ncell, eos):
    with open(out + '.params', 'w') as p:
        p.write('%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %d\n' % (
            ou, deck_value(deck, 'Ps'), jl, jg, deck_value(deck, 'mul'), deck_value(deck, 'mug'), cg, LR, dt_s, ncell))
        p.write(eos + '\n')


def main():
    run, out = Path(sys.argv[1]), sys.argv[3]
    dt_s = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    deck = next(run.parent.parent.glob('decks/' + run.name + '.nml')).read_text()
    jl, jg = deck_value(deck, 'jl'), deck_value(deck, 'jg')
    eos = deck_value(deck, 'gas_eos', str)
    cg = 343.0 if eos == 'BAROTROPIC' else math.sqrt(deck_value(deck, 'Rg') * deck_value(deck, 'Tgas'))
    n = int(deck_value(deck, 'nodes')); lp = deck_value(deck, 'Lp')
    ou, orho = max(jl, jg), 997.05
    if sys.argv[2].startswith('DUMP_'):
        d = np.loadtxt(run / 'last_rejected_state.dat', skiprows=2)
        cols = slice(3, 6) if sys.argv[2] == 'DUMP_OLD' else slice(6, 9)
        with open(out + '.state', 'w') as s:
            for row in d:
                s.write('%.17e %.17e %.17e %.17e %.17e\n' % (row[1], row[2], *row[cols]))
        write_params(out, ou, deck, jl, jg, cg, dt_s, 2 * n, eos)
        json.dump(dict(run=str(run), source=sys.argv[2], eos=eos, cg=cg, omega_u=ou), open(out + '.json', 'w'))
        return
    f = np.genfromtxt(run / 'fields.csv', delimiter=',', names=True)
    times = np.unique(f['time_s'])
    t = float(sys.argv[2])
    tsel = times[np.argmin(np.abs(times - t))]
    rows = f[f['time_s'] == tsel]
    assert rows.size == 2 * n, rows.size
    wp, wr = lp / n / LR, 1.0 / n
    centers = np.concatenate([(np.arange(n) + .5) * lp / n, lp + (np.arange(n) + .5) * LR / n])
    assert np.max(np.abs(centers - rows['s_m'])) < 1e-9, 'non-uniform mesh'
    with open(out + '.state', 'w') as s:
        for k, r in enumerate(rows):
            a, rl, rg, ul, ug = (float(r[c]) for c in ('alpha', 'rhol_kg_m3', 'rhog_kg_m3', 'ul_m_s', 'ug_m_s'))
            U = [(1 - a) * rl / orho, a * rg / orho, ((1 - a) * rl * ul + a * rg * ug) / (orho * ou)]
            s.write('%.17e %.17e %.17e %.17e %.17e\n' % (float(r['theta_rad']), wp if k < n else wr, *U))
    write_params(out, ou, deck, jl, jg, cg, dt_s, 2 * n, eos)
    json.dump(dict(run=str(run), requested_time=t, state_time=float(tsel), eos=eos, cg=cg, omega_u=ou), open(out + '.json', 'w'))
    print('state time', tsel)


if __name__ == '__main__':
    main()
