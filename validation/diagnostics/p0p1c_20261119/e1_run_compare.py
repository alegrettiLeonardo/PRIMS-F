#!/usr/bin/env python3
"""P0.2: dynamic comparison Roe (qualified) vs E1 (experimental) junction vector.

Same deck, same frozen operator except junction_flux/junction_trace. For each
run: termination, reached time, accepted/rejected steps, minimum accepted dt,
CPU, last-state junction quantities. Each E1 run is compared with its Roe
reference on the common window: max/RMS |dp_in|/mean, max |dMl|, |dMg| (kg),
and the per-cycle amplitude vs gas-trough map (iso_event_analysis.trough_map).
Usage: e1_run_compare.py OUT_JSON OUT_MD REF_DIR:RUN_DIR [REF_DIR:RUN_DIR ...]
"""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from iso_event_analysis import status, trough_map  # noqa: E402


def summary(d):
    st = status(d)
    s = np.genfromtxt(os.path.join(d, 'steps.csv'), delimiter=',', names=True)
    m = np.genfromtxt(os.path.join(d, 'monitor.csv'), delimiter=',', names=True)
    t = m['time_s']
    gin = m['mdot_g_in_kg_s'][(t >= 20) & (t <= 60)].mean() if t[-1] > 60 else m['mdot_g_in_kg_s'].mean()
    rows, astar = trough_map(t, m['p_in_Pa'], m['mdot_g_base_kg_s'], m['holdup_liquid_base'], gin)
    reason = st.get('termination_reason', 'RUNNING')
    if reason == 'RUNNING' and os.path.exists(os.path.join(d, 'OPERATOR_TERMINATED.txt')):
        reason = 'OPERATOR_TERMINATED_COST_STAGNATION'
    return dict(run=d, reason=reason, t_reached_s=float(t[-1]),
                accepted=int(st.get('accepted_steps', len(s))), rejected=int(st.get('rejected_steps', -1)),
                dt_min_accepted_s=float(s['dt_s'].min()), cpu_s=float(st['cpu_seconds']) if 'cpu_seconds' in st else None,
                last=dict(p_in_Pa=float(m['p_in_Pa'][-1]), mdot_g_base_kg_s=float(m['mdot_g_base_kg_s'][-1]),
                          alpha_pipe_end=float(m['alpha_pipe_end'][-1]), alpha_riser_base=float(m['alpha_riser_base'][-1]),
                          holdup_base=float(m['holdup_liquid_base'][-1])),
                cycles=rows, amplitude_events_Pa=astar), m


def main():
    out_json, out_md = sys.argv[1], sys.argv[2]
    res = []
    for pair in sys.argv[3:]:
        ref, run = pair.split(':')
        a, ma = summary(ref)
        b, mb = summary(run)
        tc = min(a['t_reached_s'], b['t_reached_s'])
        grid = np.arange(0, tc, 0.1)
        ia = lambda m, k: np.interp(grid, m['time_s'], m[k])
        dp = np.abs(ia(ma, 'p_in_Pa') - ia(mb, 'p_in_Pa')) / ia(ma, 'p_in_Pa').mean()
        res.append(dict(reference=a, e1=b, common_window_s=[0, float(tc)], p_in_max_rel=float(dp.max()),
                        p_in_rms_rel=float(np.sqrt((dp ** 2).mean())),
                        Ml_max_abs_kg=float(np.abs(ia(ma, 'Ml_kg') - ia(mb, 'Ml_kg')).max()),
                        Mg_max_abs_kg=float(np.abs(ia(ma, 'Mg_kg') - ia(mb, 'Mg_kg')).max())))
    json.dump(res, open(out_json, 'w'), indent=1)
    f = lambda v: '-' if v is None else ('%.4g' % v)
    with open(out_md, 'w') as m:
        m.write('# Roe vs E1: dynamic replays (frozen operator, same deck)\n\n')
        m.write('| reference | E1 run | t ref (s) | t E1 (s) | reason E1 | acc/rej ref | acc/rej E1 | dt_min E1 | cpu ref | cpu E1 | max/RMS dp_in | max dMl, dMg (kg) | A*(0.25) ref | A*(0.25) E1 |\n')
        m.write('|---|---|---:|---:|---|---|---|---:|---:|---:|---|---|---:|---:|\n')
        for r in res:
            a, b = r['reference'], r['e1']
            m.write('| %s | %s | %.3f | %.3f | %s | %d/%d | %d/%d | %.2e | %s | %s | %.4f/%.4f | %.3f, %.4f | %s | %s |\n' % (
                os.path.basename(a['run']), os.path.basename(b['run']), a['t_reached_s'], b['t_reached_s'], b['reason'],
                a['accepted'], a['rejected'], b['accepted'], b['rejected'], b['dt_min_accepted_s'], f(a['cpu_s']), f(b['cpu_s']),
                r['p_in_max_rel'], r['p_in_rms_rel'], r['Ml_max_abs_kg'], r['Mg_max_abs_kg'],
                f(a['amplitude_events_Pa']['A_star_0.25']), f(b['amplitude_events_Pa']['A_star_0.25'])))
        m.write('\nLast accepted junction state:\n\n| run | p_in Pa | mdot_g junction kg/s | alpha pipe end | alpha riser base | holdup base |\n|---|---:|---:|---:|---:|---:|\n')
        seen = set()
        for r in res:
            for x in (r['reference'], r['e1']):
                if x['run'] in seen:
                    continue
                seen.add(x['run'])
                l = x['last']
                m.write('| %s | %.0f | %.3e | %.4f | %.4f | %.4f |\n' % (os.path.basename(x['run']), l['p_in_Pa'], l['mdot_g_base_kg_s'],
                        l['alpha_pipe_end'], l['alpha_riser_base'], l['holdup_base']))
    print(open(out_md).read())


if __name__ == '__main__':
    main()
