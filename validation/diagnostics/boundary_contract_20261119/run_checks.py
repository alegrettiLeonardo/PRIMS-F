from pathlib import Path
import subprocess, time, os, json
from audit_paths import ROOT as root
from audit_paths import OUT
out=OUT/'checks';out.mkdir(parents=True,exist_ok=True)
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONPATH=str(root/'scripts'))
entries=[]
commands=[['make','-j2','FC=gfortran','LIBS=-llapack -lblas'],['make','check_all','FC=gfortran','LIBS=-llapack -lblas'],['make','test','FC=gfortran','LIBS=-llapack -lblas'],['python3','-m','unittest','discover','-s','tests','-p','test_*.py','-v'],['python3','tests/test_pipe_regime_cli.py','--binary','bin/prims_f_check_all','--output','tmp/bc_cli'],['python3','tests/test_accepted_telemetry_cli.py','--binary','bin/prims_f','--output','tmp/bc_telemetry']]
for i,cmd in enumerate(commands):
 start=time.perf_counter()
 with (out/f'{i}.log').open('w') as f:
  try:r=subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=240);code=r.returncode
  except subprocess.TimeoutExpired:code='TIMEOUT'
 entries.append({'command':cmd,'returncode':code,'wall_s':time.perf_counter()-start,'log':str(out/f'{i}.log')})
 (out/'manifest.json').write_text(json.dumps(entries,indent=2)+'\n')
 print(i,code,flush=True)
