#!/usr/bin/env python3
"""P1 summary: CHEXAL_1992 (C5) versus Bhagwat (C3), one row per run.

Per run: termination (run_status.txt), guard verdict, reached time, accepted /
rejected steps, Newton iterations, CPU; classification of the terminal event
from the solver's own composite code (domain_boundary_probe decoding):
  COMPLETED
  BLOCKED_CONSTITUTIVE_DOMAIN   CELL primitive recovery -> CHEXAL_DOMAIN (61)
                                (inclined countercurrent proposal, p.262)
  BLOCKED_JUNCTION_TRACE_DOMAIN FACE 1061: riser trace at theta=0 refused
                                (horizontal reverse gas flux, pp.264-265)
  FAIL:<reason>/<code>          anything else (e.g. Bhagwat 1012, primitive 7)
Six phase balances (2 phases x system/pipe/riser) from the ACCEPTED fluxes by
scripts/analyze_accepted_conservation.audit_ledger; event E_g(0.01) = first time
the base gas mass flow drops below 1 % of its initial value (a FLUX threshold,
not phase disappearance); inventories; for COMPLETED runs the period and
semi-amplitude of p_in over the last ten complete cycles (up-crossings of the
mean of the last third), and the relative oscillation for Stable.
Usage: summarize_campaign.py OUT_JSON OUT_MD RUN_DIR [RUN_DIR ...]
"""
import json, re, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(HERE))
from analyze_accepted_conservation import audit_ledger  # noqa: E402
import domain_boundary_probe as dbp  # noqa: E402


def cycles(t, p):
    sel = t >= t[-1] - (t[-1] - t[0]) / 3
    m = p[sel].mean()
    rel = float((p[sel].max() - p[sel].min()) / m)
    if rel < 1e-6:   # decayed to a steady state: crossings would be round-off noise
        return dict(decayed_to_steady=True, rel_oscillation_last_third=rel, p_in_mean_Pa=float(m))
    x = p - m
    up = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = [t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i]) for i in up]
    cyc = [(a, b, 0.5 * (p[i:j + 1].max() - p[i:j + 1].min())) for (i, j, a, b) in zip(up[:-1], up[1:], tc[:-1], tc[1:])]
    if len(cyc) < 10:
        return dict(complete_cycles=len(cyc))
    last = cyc[-10:]
    per = [b - a for a, b, _ in last]; amp = [c for _, _, c in last]
    return dict(complete_cycles=len(cyc), window_s=[last[0][0], last[-1][1]], period_s=float(np.mean(per)),
                period_spread_rel=float((max(per) - min(per)) / np.mean(per)), semi_amplitude_Pa=float(np.mean(amp)),
                amplitude_spread_rel=float((max(amp) - min(amp)) / np.mean(amp)), p_in_mean_Pa=float(m))


def summarize(run):
    run = Path(run)
    pr = dbp.probe(run)
    g = json.loads((run / 'GUARD_STATUS.json').read_text()) if (run / 'GUARD_STATUS.json').exists() else {}
    st = dbp.kv(run / 'run_status.txt')
    term = st.get('termination_reason', 'RUNNING' if not g else 'NO_STATUS')
    dec = pr.get('decoded', '')
    if term == 'COMPLETED':
        cls = 'COMPLETED'
    elif 'code 61' in dec and 'CELL' in dec:
        cls = 'BLOCKED_CONSTITUTIVE_DOMAIN'
    elif dec == 'FACE code 1061' or (term == 'MAX_NEWTON' and last_face_code(run) == 3061):
        cls = 'BLOCKED_JUNCTION_TRACE_DOMAIN'
    elif term == 'RUNNING':
        cls = 'RUNNING'
    else:
        code = pr.get('proposal_status') or last_face_code(run) or st.get('status')
        cls = 'FAIL:%s/%s' % (term, code)
        if term == 'PREPARE_DISSIPATION' and str(st.get('status')) == '9323':
            # 8000 + face 1200 + Roe path 100 + DFduj 10 + CLOSURE_DERIVATIVE 13:
            # a Gauss node of the lagged Roe path sits on a Chexal piece boundary
            cls = 'FAIL:ROE_PATH_DERIVATIVE_AT_CONSTITUTIVE_KINK/9323'
    m = np.genfromtxt(run / 'monitor.csv', delimiter=',', names=True)
    t = m['time_s']; gb = m['mdot_g_base_kg_s']
    below = np.nonzero(gb < 0.01 * gb[0])[0]
    out = dict(run=run.name, closure_riser=pr['closure_riser'], termination=term, guard=g.get('verdict'), classification=cls,
               decoded=dec or None, proposal_index=pr.get('proposal_index'), reached_time_s=float(t[-1]),
               accepted_steps=int(st.get('accepted_steps', -1)), rejected_steps=int(st.get('rejected_steps', -1)),
               newton_iterations=int(st.get('newton_iterations_total', -1)), cpu_s=float(st['cpu_seconds']) if 'cpu_seconds' in st else None,
               E_g_001_s=float(t[below[0]]) if below.size else None, mdot_g_base_min_rel=float(gb.min() / gb[0]),
               p_in_range_Pa=[float(m['p_in_Pa'].min()), float(m['p_in_Pa'].max())],
               alpha_riser_base_range=[float(m['alpha_riser_base'].min()), float(m['alpha_riser_base'].max())],
               Ml_kg=[float(m['Ml_kg'][0]), float(m['Ml_kg'][-1])], Mg_kg=[float(m['Mg_kg'][0]), float(m['Mg_kg'][-1])],
               Mg_riser_kg=[float(m['Mg_riser_kg'][0]), float(m['Mg_riser_kg'][-1])])
    if 'min_liquid_flux_cell' in pr:
        c = pr['min_liquid_flux_cell']; out['boundary_cell'] = dict(index=c['index'], theta_deg=c['theta_deg'], alpha=c['alpha'], jl_m_s=c['jl_mm_s'] * 1e-3)
    if term == 'COMPLETED':
        out['cycles'] = cycles(t, m['p_in_Pa'])
    try:
        led = audit_ledger(run / 'conservation.csv')
        b = led['balances']
        out['balances'] = {k: dict(max_abs_kg=v['max_abs_balance_kg'], rel_inventory=v['relative_inventory'],
                                   rel_processed=v['relative_processed']) for k, v in b.items()}
        out['ledger_arithmetic_pass'] = led['arithmetic_contract_passes']
        out['junction_integral_mismatch_kg'] = led['max_junction_integral_mismatch_kg']
        out['max_rel_balance'] = max((v['relative_inventory'] or 0.0) for v in b.values())
    except Exception as e:   # recorded, never hidden
        out['balances_error'] = repr(e)
    return out


def last_face_code(run):
    p = run / 'rejected_proposals.csv'
    if not p.exists():
        return None
    lines = p.read_text().strip().splitlines()
    if len(lines) < 2:
        return None
    try:
        return int(lines[-1].split(',')[7])
    except (ValueError, IndexError):
        return None


def main():
    rs = [summarize(r) for r in sys.argv[3:]]
    json.dump(rs, open(sys.argv[1], 'w'), indent=1)
    f = lambda v, fmt='%.4g': '-' if v is None else fmt % v
    with open(sys.argv[2], 'w') as md:
        md.write('| run | C | classification | t reached s | acc/rej | Newton | CPU s | E_g(0.01) s | min mdot_g,base/0 | boundary cell (θ°, j_l m/s) | period s | semi-amp. Pa | max |balance|/M0 | ledger |\n')
        md.write('|---|---:|---|---:|---|---:|---:|---:|---:|---|---:|---:|---:|---|\n')
        for r in rs:
            c = r.get('cycles', {}); bc = r.get('boundary_cell')
            md.write('| %s | %d | %s | %.2f | %d/%d | %d | %s | %s | %.3g | %s | %s | %s | %s | %s |\n' % (
                r['run'], r['closure_riser'], r['classification'], r['reached_time_s'], r['accepted_steps'], r['rejected_steps'],
                r['newton_iterations'], f(r['cpu_s'], '%.0f'), f(r['E_g_001_s'], '%.2f'), r['mdot_g_base_min_rel'],
                '%d (%.2f, %.1e)' % (bc['index'], bc['theta_deg'], bc['jl_m_s']) if bc else '-',
                f(c.get('period_s'), '%.3f'), f(c.get('semi_amplitude_Pa'), '%.0f'), f(r.get('max_rel_balance'), '%.1e'),
                'PASS' if r.get('ledger_arithmetic_pass') else r.get('balances_error', 'FAIL')))
    print(open(sys.argv[2]).read())


if __name__ == '__main__':
    main()
