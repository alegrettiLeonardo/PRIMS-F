from pathlib import Path
import hashlib,json,subprocess,os,time,shutil
from audit_paths import ROOT
from audit_paths import OUT
HERE=ROOT/'validation/diagnostics/boundary_contract_20261119'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
 commands=[['python3',str(HERE/'test_audit_prototypes.py'),'-v'],['python3',str(HERE/'qualify_characteristics.py')]]
 result={'commands':[],'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0]}
 for i,cmd in enumerate(commands):
  start=time.perf_counter();r=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True,env=env,timeout=60)
  (OUT/f'final_check_{i}.log').write_text(r.stdout+r.stderr);result['commands'].append({'command':cmd,'returncode':r.returncode,'wall_s':time.perf_counter()-start})
  if r.returncode:raise RuntimeError(r.stdout+r.stderr)
 result['git_status_at_freeze']=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True)
 result['head']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 result['base']='2c3e3af8c99d490ffc1dcbf9357da6c0d59d633d'
 result['production_diff_from_base']=subprocess.check_output(['git','diff',result['base'],'--','scr','Makefile','.github'],cwd=ROOT,text=True)
 result['tracked_source_hashes']={str(p.relative_to(ROOT)):sha(p) for p in list((ROOT/'scr').glob('*.f90'))+list(HERE.glob('*.py'))+list(HERE.glob('*.f90'))}
 result['binaries']={str(p.relative_to(OUT)):sha(p) for glob in ['variants/*/window','variants/*/audit.so','variants_v2*/*/window','variants_v2*/*/audit.so','release/audit.so','debug/audit.so'] for p in OUT.glob(glob) if p.name!='evidence_manifest.json'}
 result['main_status']=subprocess.check_output(['git','status','--short','--branch'],cwd=ROOT.parent/'PRIMS-F',text=True)
 result['worktrees']=subprocess.check_output(['git','worktree','list','--porcelain'],cwd=ROOT,text=True)
 result['dependencies']=subprocess.check_output(['ldd',str(ROOT/'bin/prims_f')],text=True)
 result['dependency_package_versions']=subprocess.check_output(['dpkg-query','-W','gfortran-13','liblapack3','libblas3'],text=True)
 result['production_binary_sha256']=sha(ROOT/'bin/prims_f');result['debug_binary_sha256']=sha(ROOT/'bin/prims_f_check_all')
 result['script_notes']='Original campaigns v1/v2 identify generated source bytes and binaries in each manifest. SHA here freezes the final audit harness; production unchanged.'
 result['data_hashes']={str(p.relative_to(OUT)):sha(p) for glob in ['*.json','windows/*.json','windows_v2/*.json','controls/*.json','checks/*.json','variants*/*/manifest.json'] for p in OUT.glob(glob) if p.name!='evidence_manifest.json'}
 (OUT/'evidence_manifest.json').write_text(json.dumps(result,indent=2)+'\n');print('source differences:',repr(result['production_diff_from_base']),'head:',result['head'],'main:',result['main_status'].strip())
if __name__=='__main__':main()
