#!/usr/bin/env python3
"""Classify a terminated run from its own dump (last_rejected_state.dat: the
LAST ACCEPTED state Uold and the refused Newton trial, per physical cell) and
run_status.txt / last_rejected_iteration.txt.

Status decoding (point_state / assemble composition): CELL failure
1000+10+c where c is the primitive-recovery code; for CHEXAL_1992 c is the core
status (60 invalid, 61 DOMAIN, 62 branch required, 63 no root, 64 CCFL exceeded,
65 singular, 66 numeric). FACE failures are 2000+face code.

For every riser cell of the accepted state the independent transcription
chexal92_ref.py recovers the co-current root on the liquid-flux interval
jl in [0, M/rho_l] (bisection on R(jl) with jg=(M-rho_l jl)/rho_g); the
liquid flux and its distance to zero locate the stagnant-liquid boundary.
Inclined countercurrent (jl<0<jg at theta_C>0) is OUTSIDE the 1992 source
(p.262): it is reported, never evaluated.
Usage: domain_boundary_probe.py OUT_JSON RUN_DIR [RUN_DIR ...]
"""
import json, math, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chexal92_ref as cx

CHEXAL = {60: 'INVALID', 61: 'DOMAIN', 62: 'BRANCH_REQUIRED', 63: 'NO_ROOT', 64: 'CCFL_EXCEEDED', 65: 'SINGULAR', 66: 'NUMERIC'}


def kv(path):
    d = {}
    if Path(path).exists():
        for line in Path(path).read_text().splitlines():
            if '=' in line:
                k, v = line.split('=', 1); d[k.strip()] = v.strip()
    return d


def deck_value(text, key, cast=float):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))
    m = re.search(r'\b' + key + r'\s*=\s*([^,\n]+)', text, re.I)
    return cast(m.group(1).strip().strip("'"))


def probe(run):
    run = Path(run)
    deck = next(run.parent.parent.glob('decks/' + run.name + '.nml')).read_text()
    eos = deck_value(deck, 'gas_eos', str)
    cg = 343.0 if eos == 'BAROTROPIC' else math.sqrt(deck_value(deck, 'Rg') * deck_value(deck, 'Tgas'))
    rho0, P0 = deck_value(deck, 'rho_l0'), deck_value(deck, 'P_l0')
    ou = max(deck_value(deck, 'jl'), deck_value(deck, 'jg'))
    nCg, nCl = cg * math.sqrt(rho0 / P0), 1498.0 * math.sqrt(rho0 / P0)
    st, it = kv(run / 'run_status.txt'), kv(run / 'last_rejected_iteration.txt')
    out = dict(run=run.name, closure_riser=deck_value(deck, 'closure_riser', int), termination=st.get('termination_reason'),
               reached_time_s=float(st['reached_time_s']) if 'reached_time_s' in st else None,
               accepted_steps=st.get('accepted_steps'), rejected_steps=st.get('rejected_steps'),
               proposal_status=it.get('proposal_status'), proposal_index=it.get('proposal_index'), proposal_stage=it.get('proposal_stage'))
    ps = int(it['proposal_status']) if it.get('proposal_status', '').lstrip('-').isdigit() else None
    if ps is not None and it.get('proposal_stage') == 'CELL' and 1010 <= ps < 1100:
        c = ps - 1010
        out['decoded'] = 'CELL primitive recovery code %d (%s)' % (c, CHEXAL.get(c, 'non-Chexal primitive code') if out['closure_riser'] >= 5 else 'primitive code')
    elif ps is not None and ps >= 2000:
        out['decoded'] = 'FACE code %d' % (ps - 2000)
    dump = run / 'last_rejected_state.dat'
    if not dump.exists() or out['closure_riser'] < 5:
        return out   # the Chexal recovery below is meaningless for other riser closures
    d = np.loadtxt(dump, skiprows=2)
    n = d.shape[0] // 2
    lp = deck_value(deck, 'Lp')
    cells = []
    for row in d[n:]:
        i, th = int(row[0]), float(row[1]); u1, u2, u3 = row[3:6]
        lo, hi = 1e-6, 100.0
        for _ in range(200):
            p = .5 * (lo + hi); f = u1 / (1 + (p - 1) / nCl ** 2) + u2 / (p / nCg ** 2) - 1
            lo, hi = (p, hi) if f > 0 else (lo, p)
        rl, rg = (1 + (p - 1) / nCl ** 2) * rho0, p / nCg ** 2 * rho0
        a, M = u2 / (p / nCg ** 2), u3 * rho0 * ou
        rec = dict(index=i, theta_deg=math.degrees(th), alpha=a, p_Pa=p * P0, M_kg_m2s=M)
        if M > 0:
            F = lambda jl: cx.residual(a, rl, rg, th, jl, (M - rl * jl) / rg)
            x0, x1 = 0.0, M / rl
            try:
                # end point jl=M/rho_l means jg=0 EXACTLY (avoid a -1e-17 round-off sign)
                f0, f1 = F(x0), cx.residual(a, rl, rg, th, x1, 0.0)
            except cx.DomainError as e:
                f0 = f1 = None; rec['note'] = str(e)
            if f0 is not None and f0 * f1 <= 0:
                for _ in range(200):
                    xm = .5 * (x0 + x1); fm = F(xm)
                    x0, f0, x1 = (xm, fm, x1) if f0 * fm > 0 else (x0, f0, xm)
                jl = .5 * (x0 + x1); rec.update(jl_mm_s=1e3 * jl, jg_m_s=(M - rl * jl) / rg)
            elif f0 is not None:
                rec['cocurrent_root_on_jl>=0'] = False
                rec['R_at_jl0'] = f0
        cells.append(rec)
    out['riser_cells_accepted_state'] = cells
    lows = [c for c in cells if 'jl_mm_s' in c]
    if lows:
        k = min(lows, key=lambda c: c['jl_mm_s'])
        out['min_liquid_flux_cell'] = k
    out['cells_without_cocurrent_root'] = [c['index'] for c in cells if 'jl_mm_s' not in c]
    return out


def main():
    res = [probe(r) for r in sys.argv[2:]]
    json.dump(res, open(sys.argv[1], 'w'), indent=1)
    for r in res:
        print(r['run'], r['termination'], r.get('reached_time_s'), r.get('decoded'), 'index', r.get('proposal_index'))
        if 'min_liquid_flux_cell' in r:
            c = r['min_liquid_flux_cell']
            print('   min jl cell %d theta=%.2f deg alpha=%.4f jl=%.4g mm/s; cells without co-current root: %s' % (
                c['index'], c['theta_deg'], c['alpha'], c['jl_mm_s'], r['cells_without_cocurrent_root']))


if __name__ == '__main__':
    main()
