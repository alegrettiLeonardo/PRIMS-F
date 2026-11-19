from pathlib import Path
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/home/leonardoab/Documentos/PRIMS-F-numerics-run-20261119')
OUT = Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119/offline_replay')
HERE = Path(__file__).resolve().parent
NAMES = ['mod_kinds','mod_rootfinding','mod_void_helpers','module_flow_regime','module_config',
         'module_lapack_small','module_chexal_1992','module_closure_law','module_find',
         'module_catenary','module_steady_state','module_jacobian','module_matrix_roe',
         'module_roe_solver','module_boundary_conditions','module_solver_utils',
         'module_time_step','module_simulate_pipe_riser']
OUT.mkdir(parents=True, exist_ok=True)
build = OUT / 'build'
build.mkdir(exist_ok=True)
manifest = {'source_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(), 'sources':{}}
for name in NAMES:
    source = ROOT / 'scr' / (name+'.f90')
    manifest['sources'][name] = hashlib.sha256(source.read_bytes()).hexdigest()
    shutil.copy2(source,build/source.name)
api = HERE / 'offline_api.f90'
shutil.copy2(api,build/api.name)
command = ['gfortran','-std=f2018','-O2','-fPIC','-shared','-fopenmp','-ffree-line-length-none',
           *[name+'.f90' for name in NAMES],api.name,'-llapack','-lblas','-o','offline_api.so']
completed = subprocess.run(command,cwd=build,text=True,capture_output=True)
(OUT/'build.log').write_text(completed.stdout+completed.stderr)
print(completed.stdout+completed.stderr)
manifest['command']=command
manifest['returncode']=completed.returncode
if completed.returncode:
    raise SystemExit(completed.returncode)
manifest['api_sha256']=hashlib.sha256(api.read_bytes()).hexdigest()
manifest['library_sha256']=hashlib.sha256((build/'offline_api.so').read_bytes()).hexdigest()
for case,params in [('D_SS1','.53 208000 .53 .11 528.4365 .1324434043053912 0'),
                    ('D_SS3','.48 201000 .28 .48 279.174 .577934855150798 0')]:
    folder=OUT/case
    folder.mkdir(exist_ok=True)
    original=Path('/home/leonardoab/Documentos/PRIMS-F-numerics-results-20261119/v1_2400/runs')/case
    for name in ['last_rejected_state.dat','last_rejected_iteration.txt','run_metadata.txt']:
        shutil.copy2(original/name,folder/name)
    (folder/'probe_parameters.dat').write_text(params+'\n')
    manifest[case+'_snapshot_sha256']=hashlib.sha256((folder/'last_rejected_state.dat').read_bytes()).hexdigest()
(OUT/'build_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
