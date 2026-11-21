#!/usr/bin/env python3
"""P1.1: OSC ISOTHERMAL, fixed physical perturbation (DeltaP=250 Pa, P0.3 IC),
time and space refinement with PRE-REGISTERED windows (declared before the new
runs were analysed):

  * reference pressure p_rest = p_in(0) - DeltaP  (inlet cell mean of the
    perturbation equals DeltaP exactly: f(s_in)=1 to round-off);
  * cycles = consecutive up-crossings of p_in - p_rest;
  * exponential window W = complete cycles with start time in [15, 100] s
    (the same four cycles for every discretisation; linear regime checked by
    the minimum gas-trough ratio inside W);
  * sigma, A0 from log(half-amplitude) vs cycle mid-time on W; max log residual;
  * modal fit on t in [10, 100] s: p_in - p_rest = e^{sigma t}(a cos wt + b sin wt)
    (grid in sigma, omega; linear least squares in a, b): A_modal=hypot(a,b),
    phase phi, and A_modal*cos(phi) to compare with the prescribed DeltaP;
  * phase: time of the up-crossing closest to 95 s; period = mean over W;
  * events, reported separately: E_g(0.01) = first time the junction gas mass
    flow drops below 1% of its rest value (a FLUX threshold, not phase
    disappearance); t_fail = solver termination time (Newton), if any;
  * inventories: total gas mass at t=100 s relative to t=0.
Convergence: observed order from three levels with ratio 2,
p = log2(|q_4h - q_2h| / |q_2h - q_h|), and Richardson extrapolation only when
the order is within [0.5, 2] and the differences keep their sign.
Usage: p11_time_convergence.py OUT_JSON OUT_MD RUN_DIR [RUN_DIR ...]
"""
import json, math, os, re, sys
import numpy as np

DP = 250.0


def status(d):
    s = {}
    p = os.path.join(d, 'run_status.txt')
    if os.path.exists(p):
        for line in open(p):
            if '=' in line:
                k, v = line.split('=', 1); s[k.strip()] = v.strip()
    return s


def analyse(d):
    m = np.genfromtxt(os.path.join(d, 'monitor.csv'), delimiter=',', names=True)
    t, p, g = m['time_s'], m['p_in_Pa'], m['mdot_g_base_kg_s']
    prest = p[0] - DP
    x = p - prest
    up = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = []
    for i in up:   # linear interpolation of the crossing time
        tc.append(t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i]))
    cyc = []
    for a, b, ta, tb in zip(up[:-1], up[1:], tc[:-1], tc[1:]):
        s = slice(a, b + 1)
        cyc.append(dict(t_start=ta, t_end=tb, A=0.5 * (p[s].max() - p[s].min()), trough=g[s].min() / g[0]))
    W = [c for c in cyc if 15 <= c['t_start'] <= 100 and c['t_end'] <= t[-1]]
    res = dict(run=os.path.basename(d.rstrip('/')), n_cycles=len(cyc), window_cycles=len(W))
    if len(W) >= 3:
        tm = np.array([0.5 * (c['t_start'] + c['t_end']) for c in W]); la = np.log([c['A'] for c in W])
        sg, l0 = np.polyfit(tm, la, 1)
        res.update(sigma=float(sg), A0=float(math.exp(l0)), max_log_resid=float(np.max(np.abs(la - (l0 + sg * tm)))),
                   period=float(np.mean([c['t_end'] - c['t_start'] for c in W])), min_trough_in_W=float(min(c['trough'] for c in W)),
                   A_cycle4=float(W[-1]['A']) if len(W) >= 4 else None)
    # modal fit
    sel = (t >= 10) & (t <= 100)
    ts, xs = t[sel], x[sel]
    best = None
    if sel.sum() > 50 and 'sigma' in res:
        w0 = 2 * math.pi / res['period']
        for sgm in np.linspace(res['sigma'] - 0.003, res['sigma'] + 0.003, 61):
            for om in np.linspace(0.97 * w0, 1.03 * w0, 61):
                e = np.exp(sgm * ts)
                M = np.column_stack([e * np.cos(om * ts), e * np.sin(om * ts)])
                coef, *_ = np.linalg.lstsq(M, xs, rcond=None)
                r = float(np.sqrt(np.mean((M @ coef - xs) ** 2)))
                if best is None or r < best[0]:
                    best = (r, sgm, om, coef)
        r, sgm, om, (a, b) = best
        amp = math.hypot(a, b); phi = math.atan2(-b, a)
        res.update(modal=dict(sigma=float(sgm), omega=float(om), A_modal=amp, phase_rad=phi, A_modal_cos_phi=amp * math.cos(phi),
                              rms_resid_Pa=r, rms_resid_rel=r / float(np.sqrt(np.mean(xs ** 2)))))
    if tc:
        k = int(np.argmin(np.abs(np.array(tc) - 95)))
        res['upcross_near_95s'] = float(tc[k])
    below = np.nonzero(g < 0.01 * g[0])[0]
    res['E_g_0.01'] = float(t[below[0]]) if below.size else None
    st = status(d)
    res['termination'] = st.get('termination_reason')
    res['t_fail'] = float(st['reached_time_s']) if st.get('termination_reason') not in (None, 'COMPLETED') and 'reached_time_s' in st else None
    i100 = int(np.argmin(np.abs(t - 100)))
    res['Mg_rel_change_100s'] = float(m['Mg_kg'][i100] / m['Mg_kg'][0] - 1)
    res['accepted_steps'] = int(st.get('accepted_steps', -1)); res['rejected_steps'] = int(st.get('rejected_steps', -1))
    res['cpu_s'] = float(st['cpu_seconds']) if 'cpu_seconds' in st else None
    mm = re.search(r'N(\d+)_dt([0-9.]+)_', res['run'])
    res['N'] = int(mm.group(1)); res['dt'] = float(mm.group(2))
    return res


def order(q):
    """q: values at h, 2h, 4h -> (p, extrapolated, error of finest)"""
    qh, q2, q4 = q
    d1, d2 = q4 - q2, q2 - qh
    if d1 == 0 or d2 == 0 or d1 * d2 < 0:
        return None, None, None
    p = math.log2(abs(d1) / abs(d2))
    if not 0.5 <= p <= 2.0:
        return p, None, None
    ext = qh + (qh - q2) / (2 ** p - 1)
    return p, ext, abs(qh - ext) / abs(ext)


def main():
    out_json, out_md = sys.argv[1:3]
    rs = [analyse(d) for d in sys.argv[3:]]
    rs.sort(key=lambda r: (r['N'], -r['dt']))
    seqs = {}
    for N in sorted(set(r['N'] for r in rs)):
        s = sorted([r for r in rs if r['N'] == N and 'sigma' in r], key=lambda r: -r['dt'])
        seqs['N%d_time' % N] = s
    for dt in sorted(set(r['dt'] for r in rs)):
        s = sorted([r for r in rs if r['dt'] == dt and 'sigma' in r], key=lambda r: r['N'])
        seqs['dt%g_space' % dt] = s
    conv = {}
    for name, s in seqs.items():
        for qn in ('sigma', 'period', 'E_g_0.01', 't_fail', 'upcross_near_95s', 'A_cycle4'):
            vals = [(r['dt'] if 'time' in name else 1.0 / r['N'], r.get(qn)) for r in s if r.get(qn) is not None]
            out = []
            for i in range(len(vals) - 2):
                (h4, q4), (h2, q2), (h1, q1) = vals[i], vals[i + 1], vals[i + 2]
                if abs(h4 / h2 - 2) < 1e-9 and abs(h2 / h1 - 2) < 1e-9:
                    p, ext, err = order((q1, q2, q4))
                    out.append(dict(levels=[h4, h2, h1], values=[q4, q2, q1], order=p, extrapolated=ext, finest_rel_error=err))
            if out:
                conv.setdefault(name, {})[qn] = out
    json.dump(dict(runs=rs, convergence=conv), open(out_json, 'w'), indent=1)
    f = lambda v, fmt='%.4g': '-' if v is None else fmt % v
    with open(out_md, 'w') as m:
        m.write('# P1.1 OSC ISOTHERMAL fixed perturbation: time/space refinement (pre-registered window W = cycles starting in [15,100] s)\n\n')
        m.write('| run | sigma 1/s | log resid | A0 Pa (cycle fit) | A_modal Pa | phase rad | A_modal cos(phi) Pa | period s | up-crossing ~95 s | min trough in W | E_g(0.01) s | t_fail s | termination | steps acc/rej | cpu s |\n')
        m.write('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|\n')
        for r in rs:
            md = r.get('modal', {})
            m.write('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d/%d | %s |\n' % (
                r['run'], f(r.get('sigma'), '%.5f'), f(r.get('max_log_resid'), '%.3f'), f(r.get('A0'), '%.0f'), f(md.get('A_modal'), '%.0f'),
                f(md.get('phase_rad'), '%.3f'), f(md.get('A_modal_cos_phi'), '%.0f'), f(r.get('period'), '%.3f'), f(r.get('upcross_near_95s'), '%.3f'),
                f(r.get('min_trough_in_W'), '%.2f'), f(r.get('E_g_0.01'), '%.2f'), f(r.get('t_fail'), '%.2f'), r['termination'],
                r['accepted_steps'], r['rejected_steps'], f(r.get('cpu_s'), '%.0f')))
        m.write('\n## Observed orders (three levels, ratio 2)\n\n| sequence | quantity | levels | values | order | extrapolated | finest rel. error |\n|---|---|---|---|---:|---:|---:|\n')
        for name, qs in conv.items():
            for qn, lst in qs.items():
                for c in lst:
                    m.write('| %s | %s | %s | %s | %s | %s | %s |\n' % (name, qn, ', '.join('%g' % v for v in c['levels']),
                            ', '.join('%.5g' % v for v in c['values']), f(c['order'], '%.2f'), f(c['extrapolated'], '%.5g'), f(c['finest_rel_error'], '%.2e')))
    print(open(out_md).read())


if __name__ == '__main__':
    main()
