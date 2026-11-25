#!/usr/bin/env python3
"""Frozen BAROTROPIC OSC refinement campaign (common initial value, one binary).

Usage: run_baro_campaign.py --binary FROZEN_RESTART --output DIR --jobs N:DT[,N:DT...]
       [--duration 900] [--wall-cap 21600] [--workers 4]

The binary must be the restart_refinement driver of ONE operator; it is hashed
into the manifest and never rebuilt here. The common initial U (51+51 cells,
validation/diagnostics/p0p1_20261119/common_initial_U.txt) is projected by exact
cell-overlap averaging, so integral phase masses and momentum are preserved to
round-off. A wall-clock cap ends a run as TIMEOUT (external interruption), which
is recorded separately from any solver status.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
IC = ROOT / 'validation/diagnostics/p0p1_20261119/common_initial_U.txt'
LP, LR = 57.4, 12.54625399319708


def edges(n):
    return np.r_[np.linspace(0, LP, n + 1), LP + np.linspace(0, LR, n + 1)[1:]]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    a = argparse.ArgumentParser()
    a.add_argument('--binary', type=Path, required=True)
    a.add_argument('--output', type=Path, required=True)
    a.add_argument('--jobs', required=True)
    a.add_argument('--duration', type=float, default=900.)
    a.add_argument('--wall-cap', type=float, default=21600.)
    a.add_argument('--workers', type=int, default=4)
    args = a.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    u0 = np.loadtxt(IC)
    assert u0.shape == (102, 3)
    e0 = edges(51)
    jobs = [(int(x.split(':')[0]), float(x.split(':')[1])) for x in args.jobs.split(',')]
    protocol = dict(binary=str(args.binary.resolve()), binary_sha256=sha(args.binary), initial_U=str(IC),
                    initial_U_sha256=sha(IC), duration_s=args.duration, wall_cap_s=args.wall_cap, jobs=jobs,
                    env=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1'),
                    host=platform.node(), python=platform.python_version(),
                    git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip())
    (out / f'PROTOCOL_{int(time.time())}.json').write_text(json.dumps(protocol, indent=1) + '\n')
    env = dict(os.environ, **protocol['env'])

    def run(job):
        n, dt = job
        folder = out / f'N{n}_dt{dt}'
        if (folder / 'run_status.txt').exists():
            return dict(run=folder.name, skipped='already completed')
        folder.mkdir(exist_ok=True)
        e = edges(n)
        new = np.array([np.maximum(0, np.minimum(e0[1:], r) - np.maximum(e0[:-1], l)) @ u0 / (r - l)
                        for l, r in zip(e[:-1], e[1:])])
        defect = np.diff(e) @ new - np.diff(e0) @ u0
        assert np.max(np.abs(defect)) < 1e-12
        np.savetxt(folder / 'initial_U.txt', new, fmt='%.17e')
        start = time.monotonic()
        timeout = False
        with (folder / 'stdout.log').open('w') as f:
            try:
                p = subprocess.run([str(args.binary.resolve()), str(folder), str(folder / 'initial_U.txt'), str(n), str(dt),
                                    repr(args.duration)], cwd=folder, env=env, stdout=f, stderr=subprocess.STDOUT,
                                   timeout=args.wall_cap)
                code = p.returncode
            except subprocess.TimeoutExpired:
                timeout, code = True, None
        m = dict(run=folder.name, n_pipe=n, n_riser=n, dt_max_s=dt, returncode=code, timeout_external=timeout,
                 wall_s=time.monotonic() - start, projection_defect_nd=defect.tolist(),
                 status=(folder / 'run_status.txt').read_text() if (folder / 'run_status.txt').exists() else None,
                 files={p.name: sha(p) for p in folder.iterdir() if p.is_file() and p.name != 'manifest.json'})
        (folder / 'manifest.json').write_text(json.dumps(m, indent=1) + '\n')
        print(folder.name, code, 'TIMEOUT' if timeout else '', flush=True)
        return m

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run, jobs))
    (out / f'manifest_{int(time.time())}.json').write_text(json.dumps(dict(protocol=protocol, runs=results), indent=1) + '\n')


if __name__ == '__main__':
    main()
