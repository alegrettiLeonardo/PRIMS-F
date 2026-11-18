"""Compile and run the frozen SS3 algebraic-gap probe with source provenance."""
from pathlib import Path
import hashlib
import json
import os
import subprocess

here = Path(__file__).resolve().parent
root = here.parents[1]
out = here / 'build_gap'
out.mkdir(exist_ok=True)
names = 'mod_kinds mod_rootfinding mod_void_helpers module_closure_law module_find'.split()
sources = [root / 'scr' / (name + '.f90') for name in names]
sources.append(here / 'probe_ss3_bendiksen_gap.f90')
cmd = [os.environ.get('FC', 'gfortran'), '-std=f2018', '-O2', '-g',
       '-fcheck=all,no-array-temps', '-ffpe-trap=invalid,zero,overflow',
       '-ffree-line-length-none', '-J', str(out), '-I', str(out)]
cmd += list(map(str, sources)) + ['-o', str(out / 'probe_gap')]
p = subprocess.run(cmd, capture_output=True, text=True)
(out / 'build.log').write_text(p.stdout + p.stderr)
if p.returncode:
    print(p.stderr)
    raise SystemExit(p.returncode)
(out / 'provenance.json').write_text(json.dumps({
    'command': cmd,
    'source_sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
    'snapshot_sha256': hashlib.sha256((here / 'ss3_late_snapshot.dat').read_bytes()).hexdigest(),
}, indent=2) + '\n')
p = subprocess.run([str(out / 'probe_gap')], cwd=here, capture_output=True, text=True)
(here / 'ss3_bendiksen_gap_probe.log').write_text(p.stdout + p.stderr)
print(p.stdout + p.stderr)
raise SystemExit(p.returncode)
