#!/usr/bin/env python3
"""P0.3: OSC ISOTHERMAL from the discrete stationary state + fixed physical perturbation.

Reuses the pre-registered definitions of iso_event_analysis.py (events E_g/E_h
after 20 s, half-cycle amplitudes between up-crossings of p_in about its
initial mean, A*(f) per-cycle trough map). Specific to P0.3:
  * reference mean of p_in = the stationary value at t=0 (exact rest state);
  * exponential law A_k = A0 exp(sigma t_k) fitted on the complete cycles
    whose trough ratio is still >= 0.5 (pre-registered linear regime), with
    the maximum log residual reported (the law is used only if <= 0.05);
  * predicted event time t_E = t_mid,1 + ln(A*(0.25)_ref / A_1)/sigma, compared
    with the observed E_g(0.25), where A*(0.25)_ref is the amplitude-space
    event of the legacy campaign (11.1-11.6 kPa).
Usage: p03_analysis.py OUT_JSON OUT_MD RUN_DIR [RUN_DIR ...]
"""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from iso_event_analysis import first_cross, status  # noqa: E402

A_STAR_025_LEGACY = (11056.0, 11594.0)


def cycles(t, p, g, h, gin):
    x = p - p[0]
    up = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    out = []
    for ia, ib in zip(up[:-1], up[1:]):
        s = slice(ia, ib + 1)
        out.append(dict(t_start=float(t[ia]), t_end=float(t[ib]), period_s=float(t[ib] - t[ia]),
                        amplitude_Pa=float(0.5 * (p[s].max() - p[s].min())),
                        gas_trough_ratio=float(g[s].min() / gin), holdup_base_max=float(h[s].max())))
    return out


def main():
    out_json, out_md = sys.argv[1], sys.argv[2]
    res = {}
    for d in sys.argv[3:]:
        name = os.path.basename(d.rstrip('/'))
        a = np.genfromtxt(os.path.join(d, 'monitor.csv'), delimiter=',', names=True)
        t, p, g, h = a['time_s'], a['p_in_Pa'], a['mdot_g_base_kg_s'], a['holdup_liquid_base']
        st = status(d)
        gin = float(g[0])   # stationary junction gas flux = feed at the exact rest state
        cyc = cycles(t, p, g, h, gin)
        lin = [c for c in cyc if c['gas_trough_ratio'] >= 0.5]
        fit = dict(sigma=None)
        if len(lin) >= 3:
            tm = np.array([0.5 * (c['t_start'] + c['t_end']) for c in lin]);
            la = np.log([c['amplitude_Pa'] for c in lin])
            sigma, lna0 = np.polyfit(tm, la, 1)
            resid = float(np.max(np.abs(la - (lna0 + sigma * tm))))
            fit = dict(sigma=float(sigma), A0_Pa=float(np.exp(lna0)), cycles=len(lin), max_log_residual=resid,
                       exponential_law_supported=resid <= 0.05)
        ev = {('E_g_%g' % f): first_cross(t, g, f * gin, above=False) for f in (0.5, 0.25, 0.1, 0.01)}
        ev.update({('E_h_%g' % hh): first_cross(t, h, hh, above=True) for hh in (0.8, 0.9, 0.98)})
        astar = None
        for c0, c1 in zip(cyc[:-1], cyc[1:]):
            if c0['gas_trough_ratio'] > 0.25 >= c1['gas_trough_ratio']:
                w = (c0['gas_trough_ratio'] - .25) / (c0['gas_trough_ratio'] - c1['gas_trough_ratio'])
                astar = float(np.exp((1 - w) * np.log(c0['amplitude_Pa']) + w * np.log(c1['amplitude_Pa'])))
        pred = None
        if fit['sigma'] and cyc:
            t1 = 0.5 * (cyc[0]['t_start'] + cyc[0]['t_end'])
            pred = [t1 + np.log(a_ / cyc[0]['amplitude_Pa']) / fit['sigma'] for a_ in A_STAR_025_LEGACY]
        res[name] = dict(termination_reason=st.get('termination_reason'), t_last_s=float(t[-1]),
                         t_fail_s=float(st['reached_time_s']) if 'reached_time_s' in st else None,
                         p_in_rest_Pa=float(p[0]), feed_gas_kg_s=gin, cycles=cyc, fit=fit, events_s=ev,
                         A_star_025_Pa=astar, predicted_E_g025_from_legacy_A_star_s=pred,
                         period_mean_s=float(np.mean([c['period_s'] for c in cyc[1:]])) if len(cyc) > 1 else None)
    json.dump(res, open(out_json, 'w'), indent=1)
    fmt = lambda v: '-' if v is None else '%.3f' % v
    with open(out_md, 'w') as m:
        m.write('# P0.3 OSC ISOTHERMAL from the discrete stationary state + fixed perturbation\n\n')
        m.write('| run | reason | t_fail | A1 Pa | sigma 1/s | log resid | period s | E_g 0.25 | predicted E_g 0.25 | E_g 0.01 | E_h 0.9 | A*(0.25) |\n')
        m.write('|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|\n')
        for n, r in res.items():
            f = r['fit']
            m.write('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' % (
                n, r['termination_reason'], fmt(r['t_fail_s'] or r['t_last_s']),
                '%.0f' % r['cycles'][0]['amplitude_Pa'] if r['cycles'] else '-',
                '-' if f['sigma'] is None else '%.5f' % f['sigma'],
                '-' if f['sigma'] is None else '%.3f' % f['max_log_residual'], fmt(r['period_mean_s']),
                fmt(r['events_s']['E_g_0.25']),
                '-' if r['predicted_E_g025_from_legacy_A_star_s'] is None else '%.1f-%.1f' % tuple(r['predicted_E_g025_from_legacy_A_star_s']),
                fmt(r['events_s']['E_g_0.01']), fmt(r['events_s']['E_h_0.9']),
                '-' if r['A_star_025_Pa'] is None else '%.0f' % r['A_star_025_Pa']))
        m.write('\nCycle amplitudes (Pa @ start s, trough ratio):\n\n')
        for n, r in res.items():
            m.write('- %s: %s\n' % (n, ', '.join('%.0f@%.1f(%.2f)' % (c['amplitude_Pa'], c['t_start'], c['gas_trough_ratio']) for c in r['cycles'])))
    print(open(out_md).read())


if __name__ == '__main__':
    main()
