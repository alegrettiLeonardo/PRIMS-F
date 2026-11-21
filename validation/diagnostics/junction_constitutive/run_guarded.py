#!/usr/bin/env python3
"""Run one PRIMS-F deck under pre-declared cost limits (no physics is touched).

Classification written to GUARD_STATUS.json in the run folder:
  COMPLETED / FAIL       the solver's own termination_reason (run_status.txt);
  TIMEOUT                wall-time cap reached while still advancing;
  COST_STAGNATION        simulated time advanced less than --min-advance-s over
                         the last --window-s of wall time (the run is stopped with
                         SIGTERM; its last accepted state, checkpoint and logs stay).
Limits are fixed BEFORE the run and recorded with it.
Usage: run_guarded.py BINARY DECK.nml RUN_DIR [--max-wall-s S] [--window-s S] [--min-advance-s S]
"""
import argparse, json, os, signal, subprocess, time
from pathlib import Path


def last_time(run):
    m = Path(run) / 'monitor.csv'
    if not m.exists():
        return 0.0
    try:
        with open(m, 'rb') as f:
            f.seek(0, 2)
            n = f.tell()
            f.seek(max(0, n - 4096))
            lines = f.read().decode(errors='replace').strip().splitlines()
        return float(lines[-1].split(',')[0])
    except (ValueError, IndexError):
        return 0.0


def status(run):
    p = Path(run) / 'run_status.txt'
    d = {}
    if p.exists():
        for line in p.read_text().splitlines():
            if '=' in line:
                k, v = line.split('=', 1)
                d[k.strip()] = v.strip()
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('binary'); ap.add_argument('deck'); ap.add_argument('run')
    ap.add_argument('--max-wall-s', type=float, default=14400.0)
    ap.add_argument('--window-s', type=float, default=1200.0)
    ap.add_argument('--min-advance-s', type=float, default=0.5)
    a = ap.parse_args()
    run = Path(a.run); run.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    # Monotonic clock: a machine suspension must not count as stagnation or wall time.
    t0 = time.monotonic(); t0_wall = time.time(); hist = []
    with open(str(run) + '.stdout', 'w') as out:
        p = subprocess.Popen([a.binary, a.deck], stdout=out, stderr=subprocess.STDOUT, env=env, cwd=run.parent)
        verdict = None
        while p.poll() is None:
            time.sleep(30)
            now = time.monotonic(); ts = last_time(run); hist.append((now, ts))
            old = [h for h in hist if now - h[0] >= a.window_s]
            if now - t0 > a.max_wall_s:
                verdict = 'TIMEOUT'
            elif old and ts - old[-1][1] < a.min_advance_s:
                verdict = 'COST_STAGNATION'
            if verdict:
                p.send_signal(signal.SIGTERM)
                try:
                    p.wait(60)
                except subprocess.TimeoutExpired:
                    p.kill(); p.wait()
                break
    st = status(run)
    res = dict(deck=a.deck, binary=a.binary, limits=dict(max_wall_s=a.max_wall_s, window_s=a.window_s, min_advance_s=a.min_advance_s),
               wall_s=time.monotonic() - t0, wall_clock_incl_suspension_s=time.time() - t0_wall,
               last_monitor_time_s=last_time(run), returncode=p.returncode,
               verdict=verdict or ('COMPLETED' if st.get('termination_reason') == 'COMPLETED' else 'FAIL:' + st.get('termination_reason', 'NO_STATUS')),
               solver_status=st)
    if verdict == 'COST_STAGNATION':
        w = [h for h in hist if time.monotonic() - h[0] <= a.window_s + 60]
        res['stagnation'] = dict(advance_over_window_s=(w[-1][1] - w[0][1]) if len(w) > 1 else None)
        (run / 'OPERATOR_TERMINATED.txt').write_text(
            'COST_STAGNATION (pre-declared guard): simulated time advanced %.4g s in the last %.0f s of wall time; '
            'neither TIMEOUT nor solver FAIL. Last accepted state and checkpoint preserved.\n'
            % (res['stagnation']['advance_over_window_s'] or 0.0, a.window_s))
    (run / 'GUARD_STATUS.json').write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ('verdict', 'wall_s', 'last_monitor_time_s')}))


if __name__ == '__main__':
    main()
