#!/usr/bin/env python3
"""P0.4: metrics of the isolated benchmark runs (one process at a time, pinned core).
Per run: termination, reached time, accepted/rejected steps, Newton iterations,
Jacobian and residual assemblies, backtracks, junction face solves and their CPU,
implicit face Jacobians, total CPU, CPU per simulated second, max RSS; solution
difference of p_in against a reference run (same horizon, interpolated) and the six
mass balances (max |error| / initial mass, conservation.csv).
Usage: bench_summary.py RUNS_DIR OUT_JSON OUT_MD [REF_PAIRS e.g. ss1_C16fd:ss1_C16ij ...]
"""
import csv, json, re, sys
from pathlib import Path


def kv(path):
    d = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if '=' in line:
                k, v = line.split('=', 1); d[k.strip()] = v.strip()
    return d


def series(run, col):
    r = list(csv.DictReader(open(run / 'monitor.csv')))
    return [float(x['time_s']) for x in r], [float(x[col]) for x in r]


def interp(t, x, y):
    if t <= x[0]: return y[0]
    for i in range(1, len(x)):
        if x[i] >= t:
            w = (t - x[i - 1]) / (x[i] - x[i - 1]); return y[i - 1] + w * (y[i] - y[i - 1])
    return y[-1]


def summarize(run):
    st = kv(run / 'run_status.txt')
    steps = list(csv.DictReader(open(run / 'steps.csv')))
    n = len(steps); last = steps[-1] if steps else {}
    tm = (run.parent / (run.name + '.time'))
    rss = cpu_u = None
    if tm.exists():
        txt = tm.read_text()
        m = re.search(r'Maximum resident set size \(kbytes\): (\d+)', txt); rss = int(m.group(1)) if m else None
        m = re.search(r'User time \(seconds\): ([\d.]+)', txt); cpu_u = float(m.group(1)) if m else None
    t_end = float(last['time_s']) if last else 0.0
    cons = list(csv.DictReader(open(run / 'conservation.csv')))
    m0 = [float(cons[0][k]) for k in ('Ml_kg', 'Mg_kg', 'Ml_pipe_kg', 'Mg_pipe_kg', 'Ml_riser_kg', 'Mg_riser_kg')]
    bal = ['balance_l_kg', 'balance_g_kg', 'balance_l_pipe_kg', 'balance_g_pipe_kg', 'balance_l_riser_kg', 'balance_g_riser_kg']
    mb = max(abs(float(r[c])) / max(abs(m0[i]), 1e-300) for r in cons for i, c in enumerate(bal))
    g = lambda k: int(float(last.get(k, 0) or 0)) if last else 0
    out = dict(run=run.name, termination=st.get('termination_reason'), reached_s=t_end, accepted_steps=n,
               rejected_steps=g('rejected_steps_total'), newton_iterations=g('newton_iterations_total'),
               jacobian_assemblies=g('jacobian_evaluations_total'), residual_assemblies=g('residual_evaluations_total'),
               backtracks=g('backtrack_trials_total'),
               junction_face_solves=int(st.get('junction_face_solves_total', 0) or 0),
               junction_face_cpu_s=float(st.get('junction_face_cpu_s', 0) or 0),
               implicit_face_jacobians=int(st.get('junction_face_implicit_jacobians_total', 0) or 0),
               cpu_total_s=float(last.get('cpu_total_s', 0) or 0) if last else 0.0, cpu_user_time_s=cpu_u, max_rss_kb=rss,
               max_mass_balance_over_M0=mb)
    acc = max(n, 1)
    out.update(newton_per_step=out['newton_iterations'] / acc, backtracks_per_newton=out['backtracks'] / max(out['newton_iterations'], 1),
               face_solves_per_step=out['junction_face_solves'] / acc,
               cpu_per_face_s=out['junction_face_cpu_s'] / max(out['junction_face_solves'], 1),
               cpu_per_sim_s=out['cpu_total_s'] / max(t_end, 1e-300), mean_dt_s=t_end / acc)
    return out


def main():
    runs = Path(sys.argv[1]); res = {}
    for d in sorted(p for p in runs.iterdir() if p.is_dir() and (p / 'steps.csv').exists()):
        res[d.name] = summarize(d)
    pairs = {}
    for spec in sys.argv[4:]:
        a, b = spec.split(':')
        ta, pa = series(runs / a, 'p_in_Pa'); tb, pb = series(runs / b, 'p_in_Pa')
        T = min(ta[-1], tb[-1]); grid = [T * k / 200 for k in range(201)]
        pairs[spec] = dict(horizon_s=T, max_abs_dp_in_Pa=max(abs(interp(t, ta, pa) - interp(t, tb, pb)) for t in grid))
    Path(sys.argv[2]).write_text(json.dumps(dict(runs=res, solution_differences=pairs), indent=1) + '\n')
    cols = ['termination', 'reached_s', 'accepted_steps', 'rejected_steps', 'newton_per_step', 'backtracks_per_newton',
            'jacobian_assemblies', 'residual_assemblies', 'face_solves_per_step', 'cpu_per_face_s', 'cpu_total_s',
            'cpu_per_sim_s', 'mean_dt_s', 'max_rss_kb', 'max_mass_balance_over_M0']
    with open(sys.argv[3], 'w') as md:
        md.write('| run | ' + ' | '.join(cols) + ' |\n|' + '---|' * (len(cols) + 1) + '\n')
        for k, v in res.items():
            md.write('| %s | ' % k + ' | '.join(('%.4g' % v[c]) if isinstance(v[c], float) else str(v[c]) for c in cols) + ' |\n')
        md.write('\n| pair | horizon s | max abs diff p_in (Pa) |\n|---|---|---|\n')
        for k, v in pairs.items():
            md.write('| %s | %.3f | %.3e |\n' % (k, v['horizon_s'], v['max_abs_dp_in_Pa']))
    print(open(sys.argv[3]).read())


if __name__ == '__main__':
    main()
