#!/usr/bin/env python3
"""P1 campaign: Stable/OSC/SS1/SS3 with CHEXAL_1992 (closure_riser=5) versus the
legacy Bhagwat baseline (closure_riser=3), everything else identical.

Each deck is input/benchmarks/<CASE>[_ISO].nml (primary decks, now ID5) with
ONLY these overrides, written to the deck and to PROTOCOL.json before running:
  closure_riser (5 or 3), tol_res=1e-10 (reference tolerance), out_dir, regime.
Geometry, properties, EOS, feed, mesh N51+51, dt_max 0.5 s, 2400 s, Roe,
LAGGED, MATCHED are the benchmark values. The steady initial state is computed
by the solver at start-up with the selected riser closure (no reuse of a
Bhagwat initial state for Chexal).
Each run goes through junction_constitutive/run_guarded.py with pre-declared
cost limits (COMPLETED / FAIL / TIMEOUT / COST_STAGNATION).
Usage: run_campaign.py BINARY RESULTS_ROOT [--workers 6] [--max-wall-s 14400] [--grid N:dt ...]
"""
import argparse, concurrent.futures as cf, hashlib, json, os, re, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
GUARD = ROOT / 'validation/diagnostics/junction_constitutive/run_guarded.py'
CASES = ['D_Stable', 'D_OSC', 'D_SS1', 'D_SS3']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def assign(text, key, value):
    """Replace key=value on NON-COMMENT lines only (a provenance comment that
    mentions a key must never be the one edited); exactly one assignment."""
    pat = re.compile(r'(\b' + re.escape(key) + r'\s*=)\s*(\'[^\']*\'|[^,\n]+)', re.I)
    lines, hits = text.splitlines(keepends=True), 0
    for i, ln in enumerate(lines):
        if ln.lstrip().startswith('!'):
            continue
        new, n = pat.subn(lambda m: m[1] + value, ln)
        if n:
            lines[i], hits = new, hits + n
    if hits > 1:
        raise RuntimeError(f'{key} assigned {hits} times')
    if hits == 1:
        return ''.join(lines)
    return re.sub(r'(?m)^\s*/\s*$', f'  {key}={value}\n/', text)


def main():
    a = argparse.ArgumentParser()
    a.add_argument('binary', type=Path)
    a.add_argument('root', type=Path)
    a.add_argument('--workers', type=int, default=6)
    a.add_argument('--max-wall-s', type=float, default=14400.0)
    a.add_argument('--window-s', type=float, default=1800.0)
    a.add_argument('--min-advance-s', type=float, default=1.0)
    a.add_argument('--eos', nargs='+', default=['BAROTROPIC', 'ISOTHERMAL'])
    a.add_argument('--closures', nargs='+', type=int, default=[5, 3])
    a.add_argument('--cases', nargs='+', default=CASES)
    a.add_argument('--grid', nargs='+', default=['51:0.5'], help='N:dt_max pairs (cells per section : seconds)')
    args = a.parse_args()
    root = args.root.resolve()
    (root / 'decks').mkdir(parents=True, exist_ok=True)
    (root / 'runs').mkdir(exist_ok=True)
    jobs = []
    for eos in args.eos:
        for case in args.cases:
            src = ROOT / 'input/benchmarks' / (case + ('_ISO' if eos == 'ISOTHERMAL' else '') + '.nml')
            for cid, g in ((c, g) for c in args.closures for g in args.grid):
                n, dt = g.split(':')
                name = f'{case}_C{cid}_{eos}_N{n}_dt{dt}'
                run = root / 'runs' / name
                text = src.read_text()
                over = [('closure_riser', str(cid)), ('tol_res', '1e-10'), ('out_dir', repr(str(run))), ('regime', repr(name))]
                if g != '51:0.5':   # benchmark grid otherwise; refinement changes only N and dt_max
                    over += [('nodes', n), ('n_riser', n), ('dt_max_s', dt), ('output_interval_s', str(min(0.1, float(dt))))]
                for k, v in over:
                    text = assign(text, k, v)
                deck = root / 'decks' / (name + '.nml')
                deck.write_text(text)
                body = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))
                got = re.search(r'\bclosure_riser\s*=\s*(\d+)', body).group(1)
                if int(got) != cid:
                    raise RuntimeError(f'{deck}: closure_riser={got}, expected {cid}')
                jobs.append(dict(name=name, case=case, eos=eos, closure_riser=cid, grid=g, source_deck=str(src.relative_to(ROOT)),
                                 source_sha256=sha(src), deck=str(deck), deck_sha256=sha(deck), run=str(run)))
    protocol = dict(created=time.strftime('%Y-%m-%dT%H:%M:%S%z'), binary=str(args.binary.resolve()), binary_sha256=sha(args.binary),
                    git_head=subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip(),
                    git_dirty=bool(subprocess.check_output(['git', '-C', str(ROOT), 'status', '--porcelain', 'scr'], text=True).strip()),
                    limits=dict(max_wall_s=args.max_wall_s, window_s=args.window_s, min_advance_s=args.min_advance_s),
                    env=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'), jobs=jobs)
    (root / f'PROTOCOL_{int(time.time())}.json').write_text(json.dumps(protocol, indent=1) + '\n')
    env = dict(os.environ, **protocol['env'])

    def go(j):
        if (Path(j['run']) / 'GUARD_STATUS.json').exists():
            return j['name'], 'skipped'
        Path(j['run']).mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(GUARD), str(args.binary.resolve()), j['deck'], j['run'], '--max-wall-s', str(args.max_wall_s),
               '--window-s', str(args.window_s), '--min-advance-s', str(args.min_advance_s)]
        with open(Path(j['run']) / 'guard.log', 'w') as f:
            p = subprocess.run(cmd, env=env, stdout=f, stderr=subprocess.STDOUT)
        return j['name'], p.returncode

    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        for name, rc in ex.map(go, jobs):
            print(name, rc, flush=True)


if __name__ == '__main__':
    main()
