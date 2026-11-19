from pathlib import Path
import hashlib
import json
import os
import subprocess
here=Path(__file__).resolve().parent
dev=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-20261119')
base=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')
out=base/'minimax_steps_v3'
out.mkdir(exist_ok=True)
obj=dev/'tmp/release'
names=['mod_kinds','mod_rootfinding','mod_void_helpers','module_flow_regime','module_config',
       'module_lapack_small','module_chexal_1992','module_closure_law','module_find','module_catenary',
       'module_steady_state','module_jacobian','module_matrix_roe','module_roe_solver',
       'module_boundary_conditions','module_solver_utils','module_time_step','module_active_set','module_simulate_pipe_riser']
env=os.environ.copy()
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    env[name]='1'
cmd=['gfortran','-std=f2018','-O2','-fopenmp','-ffree-line-length-none','-I'+str(obj),'-J'+str(out),
     str(here/'offline_api.f90'),str(here/'replay_one_step.f90'),*[str(obj/(name+'.o')) for name in names],
     '-llapack','-lblas','-o',str(out/'replay')]
r=subprocess.run(cmd,cwd=out,env=env,capture_output=True,text=True)
(out/'compile.log').write_text(r.stdout+r.stderr)
if r.returncode:
    print(r.stderr)
    raise SystemExit(r.returncode)
report={'source_hashes':{name:hashlib.sha256((dev/'scr'/(name+'.f90')).read_bytes()).hexdigest() for name in names},'runs':[]}
for case in ('D_SS1_active_terminal','D_SS3_active_terminal','D_SS3_active_terminal_v2'):
    fixture=base/'offline_replay'/case
    for mode in ('LINE_SEARCH','ACTIVE_SET'):
        folder=out/(case+'_'+mode)
        folder.mkdir(exist_ok=True)
        r=subprocess.run([str(out/'replay'),str(fixture),str(folder),mode],env=env,text=True,capture_output=True)
        (folder/'replay.log').write_text(r.stdout+r.stderr)
        item={'case':case,'mode':mode,'returncode':r.returncode,'snapshot_sha256':hashlib.sha256((fixture/'last_rejected_state.dat').read_bytes()).hexdigest(),
              'run_status':(folder/'run_status.txt').read_text()}
        report['runs'].append(item)
        (out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(item),flush=True)
