#!/usr/bin/env python3
"""Reproduce the frozen SS1 junction-domain diagnosis with current sources.

FC=gfortran LIBS='-llapack -lblas' python3 validation/diagnostics/bg_ss1_stall/build_probe.py
No production source, input state or requested physical model is modified.
"""
from pathlib import Path
import hashlib
import json
import os
import shlex
import subprocess

here = Path(__file__).resolve().parent
root = here.parents[2]
build = here / 'build'
build.mkdir(exist_ok=True)
names = ('mod_kinds mod_rootfinding mod_void_helpers module_flow_regime module_config '
         'module_lapack_small module_closure_law module_find module_catenary module_steady_state '
         'module_jacobian module_matrix_roe module_roe_solver module_boundary_conditions '
         'module_solver_utils module_time_step module_simulate_pipe_riser').split()
files = [root / 'scr' / (name + '.f90') for name in names] + [here / 'probe_bg_ss1_stall.f90']
command = shlex.split(os.environ.get('FC', 'gfortran')) + [
    '-std=f2018', '-O0', '-g', '-ffree-line-length-none', '-fcheck=all,no-array-temps',
    '-ffpe-trap=invalid,zero,overflow', '-finit-real=snan', '-fopenmp',
    '-J' + str(build), '-I' + str(build), *map(str, files), '-Wl,-z,noexecstack',
    *shlex.split(os.environ.get('LIBS', '-llapack -lblas')), '-o', str(build / 'probe_bg_ss1_stall')]
result = subprocess.run(command, capture_output=True, text=True, timeout=180)
(build / 'build.log').write_text(result.stdout + result.stderr)
provenance = {
    'command': command, 'build_returncode': result.returncode,
    'source_sha256': {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
    'snapshot_sha256': hashlib.sha256((here / 'last_rejected_state.dat').read_bytes()).hexdigest(),
}
if result.returncode == 0:
    result = subprocess.run([str(build / 'probe_bg_ss1_stall')], cwd=here,
                            capture_output=True, text=True, timeout=60)
    (here / 'probe.log').write_text(result.stdout)
    (build / 'runtime_stderr.log').write_text(result.stderr)
    provenance['probe_returncode'] = result.returncode
    print(result.stdout)
else:
    print(result.stderr)
(build / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
raise SystemExit(result.returncode)
