from pathlib import Path
import json
import subprocess
import hashlib
import os
here=Path(__file__).resolve().parent
dev=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-20261119')
root=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')
out=root/'minimax_replays_v3'
out.mkdir(exist_ok=True)
objects=dev/'tmp/release'
names=['mod_kinds','mod_rootfinding','mod_void_helpers','module_flow_regime','module_config',
       'module_lapack_small','module_chexal_1992','module_closure_law','module_find','module_catenary',
       'module_steady_state','module_jacobian','module_matrix_roe','module_roe_solver',
       'module_boundary_conditions','module_solver_utils','module_time_step','module_active_set','module_simulate_pipe_riser']
env=os.environ.copy()
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    env[name]='1'
cmd=['gfortran','-std=f2018','-O2','-fopenmp','-ffree-line-length-none','-I'+str(objects),'-J'+str(out),
     str(here/'offline_api.f90'),str(here/'replay_minimax.f90'),*[str(objects/(name+'.o')) for name in names],
     '-llapack','-lblas','-o',str(out/'replay')]
r=subprocess.run(cmd,cwd=out,env=env,capture_output=True,text=True)
(out/'compile.log').write_text(r.stdout+r.stderr)
if r.returncode:
    print(r.stderr)
    raise SystemExit(r.returncode)
manifest={'sources':{name:hashlib.sha256((dev/'scr'/(name+'.f90')).read_bytes()).hexdigest() for name in names},'runs':[]}
for case in ('D_SS3_active_terminal','D_SS3_active_terminal_v2'):
    fixture=root/'offline_replay'/case
    r=subprocess.run([str(out/'replay'),str(fixture)],env=env,text=True,capture_output=True)
    (out/(case+'.log')).write_text(r.stdout+r.stderr)
    manifest['runs'].append({'case':case,'returncode':r.returncode,'snapshot_sha256':hashlib.sha256((fixture/'last_rejected_state.dat').read_bytes()).hexdigest(),'log':r.stdout+r.stderr})
    print(case,r.returncode)
    print('\n'.join(row for row in r.stdout.splitlines() if not row.startswith('proposal') or ' 1 ' in row[:18]))
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
