from pathlib import Path
import json, subprocess, time, os, math
import numpy as np
from audit_local import Audit,scalar
from audit_paths import OUT
HERE=Path(__file__).parent
def main():
 records=[]
 for v in ['standard','volume','partition']:
  b=OUT/'variants'/v
  cmd=['gfortran','-std=f2018','-O2','-I'+str(b),str(HERE/'test_boundary_contract.f90'),str(b/'audit.so'),'-Wl,-rpath,'+str(b),'-o',str(b/'regression')]
  subprocess.run(cmd,check=True,cwd=b)
  for mode in ['contract','partition']:
   r=subprocess.run([str(b/'regression'),mode],capture_output=True,text=True)
   (b/f'regression_{mode}.log').write_text(r.stdout+r.stderr)
   records.append({'variant':v,'mode':mode,'returncode':r.returncode,'output':r.stdout,'compile':cmd})
 (OUT/'regressions.json').write_text(json.dumps(records,indent=2)+'\n')
 a=Audit(OUT/'release/audit.so');saved=a.par.copy();a.par[17:20]=0;a.setup()
 # Pure-liquid analytic hydrostatic EOS reference, zero physical feed; no internal face.
 # Threshold BEFORE evaluation: |R| <= dt*|midpoint_source-exact_source| + roundoff.
 hydro=[]
 for ds in [1.,.5,.25]:
  length=a.par[4]*ds;theta=.7;cl2=1498.**2;pt=a.par[16]*a.par[2]
  pc=pt+(pt-101325+997.05*cl2)*math.expm1(9.81*math.sin(theta)*length/(2*cl2))
  pb=pt+(pt-101325+997.05*cl2)*math.expm1(9.81*math.sin(theta)*length/cl2)
  u=a.state(p=pc/a.par[2],alpha=0,mom=0)[None,:];dt=.001
  q=a.assemble(u,theta=[theta],width=[ds],npip=1,jac=0,dt=dt)
  # Midpoint-source residual if an exact hydrostatic lower-face pressure were used.
  correct=q['r'][2,0]-dt*(pb-pc)/(a.par[1]*a.par[0]**2*ds)
  expected_bound=abs(correct)+1e-11
  hydro.append({'ds_nd':ds,'length_m':length,'status':q['status'],'momentum_residual_nd':q['r'][2,0],
   'hydrostatic_lower_face_reference_Pa':pb,'executed_lower_face_pressure_Pa':pc,
   'exact_lower_face_midpoint_source_residual_nd':correct,'predefined_truncation_bound':expected_bound,
   'stationary_test_pass':bool(abs(q['r'][2,0])<=expected_bound),'phase_fluxes_in':q['fin'][:2],
   'phase_fluxes_out':q['fout'][:2]})
 (OUT/'hydrostatic.json').write_text(json.dumps(scalar(hydro),indent=2)+'\n')
 # Declared short controls, distinct from the pipeline/riser benchmark.
 # Same physical setup in each pair: closure2, constant angle, 20 physical FV cells,
 # same feed and fixed pressure. Only initial momentum distinguishes fill/drain.
 a.par=saved;a.setup();p=a.par
 folder=OUT/'controls';folder.mkdir(exist_ok=True)
 results=[]
 for name,mom,alpha,mg in [('filling',.2,.25,.001),('draining',.8,.25,.001),('near_absent_gas',.6,1e-8,.001)]:
  for variant in ['standard','volume']:
   dest=folder/(name+'_'+variant);dest.mkdir(exist_ok=False)
   n=24;params=p.copy();params[18]=mg
   a.par=params;a.setup();u=a.state(alpha=alpha,mom=mom)
   lines=[f'{n} 2 10 2 2 2.0 1e-10 1e-10 25 F 2000',
       ' '.join(format(x,'.17e') for x in params),
       ' '.join(format(x/params[3],'.17e') for x in [.001,.02,1e-8])]
   lines += [str(i+1)+' '+' '.join(format(x,'.17e') for x in [-.035,.1,*u]) for i in range(n)]
   inp=dest/'input.dat';inp.write_text('\n'.join(lines)+'\n')
   cmd=[str(OUT/'variants'/variant/'window'),str(inp),str(dest)];start=time.perf_counter()
   with (dest/'stdout.log').open('w') as f:
    try:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=60,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));code=r.returncode
    except subprocess.TimeoutExpired:code='TIMEOUT'
   record={'control':name,'variant':variant,'returncode':code,'command':cmd,'wall_s':time.perf_counter()-start}
   if (dest/'run_status.txt').exists():record['status']=(dest/'run_status.txt').read_text()
   results.append(record);print(name,variant,code,record.get('status','')[:100],flush=True)
 (OUT/'controls/manifest.json').write_text(json.dumps(results,indent=2)+'\n')
 print('regressions',[(r['variant'],r['mode'],r['returncode']) for r in records])
 print('hydrostatic',scalar(hydro))
if __name__=='__main__':main()
