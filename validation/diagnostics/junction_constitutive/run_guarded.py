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


def _tail_time(path, column):
    with open(path, 'rb') as f:
        f.seek(0, 2)
        n = f.tell()
        f.seek(max(0, n - 4096))
        lines = f.read().decode(errors='replace').strip().splitlines()
    return float(lines[-1].split(',')[column])


def last_time(run):
    """Simulated time of the last ACCEPTED step: steps.csv (one row per accepted
    step, column time_s); monitor.csv (written at output_interval_s) is the fallback.
    Round 20261122: before, only monitor.csv was read, which lags by up to one
    output interval."""
    s = Path(run) / 'steps.csv'
    try:
        if s.exists():
            return _tail_time(s, 1)
    except (ValueError, IndexError):
        pass
    m = Path(run) / 'monitor.csv'
    if not m.exists():
        return 0.0
    try:
        return _tail_time(m, 0)
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
        verdict = None; trigger = None
        while p.poll() is None:
            time.sleep(30)
            now = time.monotonic(); ts = last_time(run); hist.append((now, ts))
            old = [h for h in hist if now - h[0] >= a.window_s]
            if now - t0 > a.max_wall_s:
                verdict = 'TIMEOUT'
            elif old and ts - old[-1][1] < a.min_advance_s:
                verdict = 'COST_STAGNATION'
                # the exact comparison that fired (audit, round 20261122)
                trigger = dict(window_start_monotonic_s=old[-1][0] - t0, window_end_monotonic_s=now - t0,
                               window_length_s=now - old[-1][0], sim_time_start_s=old[-1][1], sim_time_end_s=ts,
                               advance_s=ts - old[-1][1], min_advance_s=a.min_advance_s,
                               comparison='advance_s < min_advance_s')
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
    res['history_monotonic_s_sim_time_s'] = [[round(h[0] - t0, 3), h[1]] for h in hist]
    if verdict == 'COST_STAGNATION':
        res['stagnation'] = trigger
        (run / 'OPERATOR_TERMINATED.txt').write_text(
            'COST_STAGNATION (pre-declared guard): simulated time advanced %.6g s between monotonic %.1f s and %.1f s '
            '(window %.1f s >= %.0f s); %.6g < min_advance_s=%.4g. Neither TIMEOUT nor solver FAIL. '
            'Last accepted state and checkpoint preserved.\n'
            % (trigger['advance_s'], trigger['window_start_monotonic_s'], trigger['window_end_monotonic_s'],
               trigger['window_length_s'], a.window_s, trigger['advance_s'], a.min_advance_s))
    (run / 'GUARD_STATUS.json').write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ('verdict', 'wall_s', 'last_monitor_time_s')}))


if __name__ == '__main__':
    main()
