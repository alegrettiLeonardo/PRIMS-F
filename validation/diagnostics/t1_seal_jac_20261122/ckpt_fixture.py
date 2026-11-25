#!/usr/bin/env python3
"""Build an offline_api fixture (probe_parameters.dat + last_rejected_state.dat) from a
PRIMS-F checkpoint.dat: Uold = trial = accepted U, header t_nd dt_nd (dt = next attempted
step), residual 0. probe_parameters.dat is copied from a template fixture of the same case.
Usage: ckpt_fixture.py CHECKPOINT TEMPLATE_FIXTURE OUT_DIR"""
import shutil, sys
from pathlib import Path
ck, tpl, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)
shutil.copy(tpl / 'probe_parameters.dat', out / 'probe_parameters.dat')
L = ck.read_text().splitlines()
h = L[1].split(); n, nd = int(h[0]), int(h[1]); t, dt = float(h[5]), float(h[6])
rows = [l.split() for l in L[2:2 + n]]
with open(out / 'last_rejected_state.dat', 'w') as f:
    f.write('# from %s: t_nd dt_nd residual code; then i,theta,width,Uold(3),Utrial(3),residual(3)\n' % ck)
    f.write('%25.16e%25.16e%25.16e%10d\n' % (t, dt, 0.0, 0))
    for r in rows[nd:n - nd]:
        i = int(r[0]); th, w = float(r[2]), float(r[3]); u = [float(v) for v in r[4:7]]
        f.write('%8d' % i + ''.join('%25.16e' % v for v in [th, w] + u + u + [0.0, 0.0, 0.0]) + '\n')
print('fixture', out, 't_nd', t, 'dt_nd', dt)
