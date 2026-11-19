from pathlib import Path
import hashlib
import json
import os
import subprocess

HERE=Path(__file__).resolve().parent
DEV=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-20261119')
OUT=Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119/offline_replay/fortran')
OUT.mkdir(exist_ok=True)
objects=DEV/'tmp/release'
order=['mod_kinds','mod_rootfinding','mod_void_helpers','module_flow_regime','module_config',
       'module_lapack_small','module_chexal_1992','module_closure_law','module_find','module_catenary',
       'module_steady_state','module_jacobian','module_matrix_roe','module_roe_solver',
       'module_boundary_conditions','module_solver_utils','module_time_step','module_active_set','module_simulate_pipe_riser']
env=os.environ.copy()
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    env[name]='1'
cmd=['gfortran','-std=f2018','-O2','-fopenmp','-ffree-line-length-none','-I'+str(objects),
     '-J'+str(OUT),str(HERE/'offline_api.f90'),str(HERE/'replay_active_fortran.f90'),
     *[str(objects/(name+'.o')) for name in order],'-llapack','-lblas','-o',str(OUT/'replay')]
compile_result=subprocess.run(cmd,cwd=OUT,env=env,text=True,capture_output=True)
(OUT/'compile.log').write_text(compile_result.stdout+compile_result.stderr)
if compile_result.returncode:
    print(compile_result.stderr)
    raise SystemExit(compile_result.returncode)
manifest={'sources':{},'compile_command':cmd,'runs':{}}
for name in order:
    manifest['sources'][name]=hashlib.sha256((DEV/'scr'/(name+'.f90')).read_bytes()).hexdigest()
for case in ('D_SS1','D_SS3'):
    result=subprocess.run([str(OUT/'replay'),str(OUT.parent/case)],env=env,text=True,capture_output=True)
    (OUT/(case+'.log')).write_text(result.stdout+result.stderr)
    manifest['runs'][case]={'returncode':result.returncode,'log':result.stdout+result.stderr}
    print(case,result.returncode,result.stdout,result.stderr)
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
