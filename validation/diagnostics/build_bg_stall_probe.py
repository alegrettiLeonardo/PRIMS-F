#!/usr/bin/env python3
"""Build and reproduce the frozen Bhagwat SS3 domain diagnosis (no transient).

Run from any directory, for example:
  FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/build_bg_stall_probe.py

The default input is the adjacent ss3_bg_stall_snapshot.dat. Compiler and
library selection use FC/LIBS; paths are resolved relative to this script.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess


def main() -> int:
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--snapshot', type=Path, default=here / 'ss3_bg_stall_snapshot.dat')
    parser.add_argument('--build-dir', type=Path, default=here / 'build_bg_stall')
    args = parser.parse_args()
    build = args.build_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    snapshot = args.snapshot.resolve()
    names = ('mod_kinds mod_rootfinding mod_void_helpers module_flow_regime module_config '
             'module_lapack_small module_closure_law module_find module_catenary module_steady_state '
             'module_jacobian module_matrix_roe module_roe_solver module_boundary_conditions '
             'module_solver_utils module_time_step module_simulate_pipe_riser').split()
    sources = [root / 'scr' / f'{name}.f90' for name in names]
    sources.append(here / 'probe_ss3_bg_stall.f90')
    command = [*shlex.split(os.environ.get('FC', 'gfortran')),
               '-std=f2018', '-O0', '-g', '-fcheck=all,no-array-temps',
               '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none', '-fopenmp',
               '-J' + str(build), '-I' + str(build), *map(str, sources), '-Wl,-z,noexecstack',
               *shlex.split(os.environ.get('LIBS', '-llapack -lblas')),
               '-o', str(build / 'probe_ss3_bg_stall')]
    hashes = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
              for path in sources}
    result = subprocess.run(command, cwd=build, capture_output=True, text=True, timeout=60)
    (build / 'build.log').write_text(result.stdout + result.stderr)
    if result.returncode:
        print(result.stdout + result.stderr)
        return result.returncode
    provenance = {'compiler_command': command, 'source_sha256': hashes,
                  'snapshot': str(snapshot),
                  'snapshot_sha256': hashlib.sha256(snapshot.read_bytes()).hexdigest()}
    (build / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    result = subprocess.run([str(build / 'probe_ss3_bg_stall'), str(snapshot)],
                            cwd=here, capture_output=True, text=True, timeout=60)
    (here / 'ss3_bg_stall_probe.log').write_text(result.stdout + result.stderr)
    print(result.stdout + result.stderr)
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
