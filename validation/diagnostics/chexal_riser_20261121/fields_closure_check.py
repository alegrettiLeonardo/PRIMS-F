#!/usr/bin/env python3
"""End-to-end check that the ACCEPTED riser states of a run satisfy the
Chexal (1992) drift relation, recomputed by the independent transcription
chexal92_ref.py from fields.csv (alpha, rho_l, rho_g, u_l, u_g, theta).

For every saved time and every riser cell (s > Lp, theta >= 0):
  jl = (1-alpha) ul, jg = alpha ug,
  R = alpha (C0 (jl+jg) + Vgj) - jg,
  R_rel = |R| / max(|jg|, |alpha C0 (jl+jg)|, |alpha Vgj|).
A Chexal run (closure_riser=5) must give R_rel at round-off/solver tolerance in
every cell; a Bhagwat run (closure_riser=3) is the discriminating control
(its states satisfy a different relation, so R_rel is O(1e-2..1e-1)).
States outside the source domain (inclined countercurrent, negative horizontal
flux) are COUNTED, never evaluated by extrapolation.
Usage: fields_closure_check.py OUT_JSON RUN_DIR [RUN_DIR ...]  (deck values: D, g, mu, sigma, Lp)
"""
import json, math, re, sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chexal92_ref as cx


def deck_value(text, key):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))
    m = re.search(r'\b' + key + r'\s*=\s*([^,\n]+)', text, re.I)
    return float(m.group(1).strip().strip("'"))


def check(run):
    run = Path(run)
    deck = Path(json.loads((run / 'GUARD_STATUS.json').read_text())['deck']).read_text() if (run / 'GUARD_STATUS.json').exists() \
        else None
    if deck is None:
        deck = next(run.parent.parent.glob('decks/' + run.name + '.nml')).read_text()
    kw = dict(d=deck_value(deck, 'D'), g=deck_value(deck, 'G'), ml=deck_value(deck, 'mu_l'), mg=deck_value(deck, 'mu_g'),
              sig=deck_value(deck, 'sigma'))
    Lp = deck_value(deck, 'Lp')
    cid = int(deck_value(deck, 'closure_riser'))
    f = np.genfromtxt(run / 'fields.csv', delimiter=',', names=True)
    sel = (f['s_m'] > Lp) & (f['theta_rad'] >= 0)
    rows = f[sel]
    rel, dom, worst = [], 0, None
    for r in rows:
        a = float(r['alpha']); jl = (1 - a) * float(r['ul_m_s']); jg = a * float(r['ug_m_s'])
        try:
            c0, v = cx.coefficients(a, float(r['rhol_kg_m3']), float(r['rhog_kg_m3']), float(r['theta_rad']), jl, jg, **kw)
        except cx.DomainError:
            dom += 1
            continue
        R = a * (c0 * (jl + jg) + v) - jg
        sc = max(abs(jg), abs(a * c0 * (jl + jg)), abs(a * v))
        q = abs(R) / sc
        rel.append(q)
        if worst is None or q > worst[0]:
            worst = (q, float(r['time_s']), float(r['s_m']), a, jl, jg)
    rel = np.array(rel)
    return dict(run=run.name, closure_riser=cid, riser_cell_states=int(sel.sum()), evaluated=int(rel.size), outside_source_domain=dom,
                times=int(np.unique(rows['time_s']).size), t_last=float(rows['time_s'].max()) if rows.size else None,
                R_rel_max=float(rel.max()) if rel.size else None, R_rel_median=float(np.median(rel)) if rel.size else None,
                R_rel_p99=float(np.quantile(rel, 0.99)) if rel.size else None,
                worst=dict(R_rel=worst[0], t=worst[1], s=worst[2], alpha=worst[3], jl=worst[4], jg=worst[5]) if worst else None)


def main():
    out = [check(r) for r in sys.argv[2:]]
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
    for o in out:
        print('%-36s C%d states=%7d eval=%7d outside=%d times=%4d R_rel max=%.2e p99=%.2e median=%.2e' % (
            o['run'], o['closure_riser'], o['riser_cell_states'], o['evaluated'], o['outside_source_domain'], o['times'],
            o['R_rel_max'] or float('nan'), o['R_rel_p99'] or float('nan'), o['R_rel_median'] or float('nan')))


if __name__ == '__main__':
    main()
