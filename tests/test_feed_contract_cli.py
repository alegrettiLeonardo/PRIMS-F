"""EOS/feed integration checks. Outputs are new; old evidence is never overwritten."""
from pathlib import Path
import argparse,csv,hashlib,json,math,os,subprocess,time
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--baseline',type=Path);a=ap.parse_args()
 binary=a.binary.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 base=dict(nodes=7,n_riser=7,X=6.435,Z=9.886,Lp=57.4,beta_deg=2.,D=.05254,eps=.000046,rho_l0=997.05,P_l0=101325.,Cl=1498.,Rg=287.,mu_l=.00089,mu_g=.0000185,sigma=.0728,G=9.81,Ps=206000.,jl=.3,jg=.6,Tgas=293.15,P_standard=101300.,T_standard=293.15,Cg_fixed=343.,gas_flow_basis='STANDARD',time_dim=.2,tol=1e-12,tola=1e-13,tol_res=1e-10,tol_inc=1e-10,newton_max=25,dt_max_s=.01,dt_start_s=.001,dt_min_s=1e-8,perturbation_amplitude=0.,output_interval_s=.01,field_interval_s=.01,regime='FEED_TEST',closure_riser=2)
 records=[]
 def run(name,patch,binary=binary):
  dest=out/name;dest.mkdir();params=dict(base,**patch);params['out_dir']=str(dest)
  def value(v):return "'"+v+"'" if isinstance(v,str) else format(v,'.17g') if isinstance(v,float) else str(v)
  inp=dest/'input.nml';inp.write_text('&PARAMS\n'+''.join(' '+k+'='+value(v)+',\n' for k,v in params.items())+'/\n')
  start=time.perf_counter()
  with (dest/'stdout.log').open('w') as f:r=subprocess.run([str(binary),str(inp)],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),timeout=60)
  item=dict(name=name,returncode=r.returncode,wall_s=time.perf_counter()-start,input_sha256=sha(inp),binary_sha256=sha(binary));records.append(item)
  if (dest/'run_status.txt').exists():item['status']=(dest/'run_status.txt').read_text()
  return dest,item
 for eos in ['BAROTROPIC','ISOTHERMAL']:
  legacy,item=run(eos+'_legacy',dict(gas_eos=eos));assert item['returncode']==0
  explicit,item=run(eos+'_mass',dict(gas_eos=eos,liquid_flow_basis='MASS',gas_flow_basis='MASS',mul=997.05*.3,mug=101300/(287*293.15)*.6));assert item['returncode']==0
  for n in ['monitor.csv','fields.csv','conservation.csv']:assert sha(legacy/n)==sha(explicit/n),(eos,n)
  if a.baseline:
   old,item=run(eos+'_old_binary',dict(gas_eos=eos),a.baseline.resolve());assert item['returncode']==0
   for n in ['monitor.csv','fields.csv','conservation.csv']:assert sha(legacy/n)==sha(old/n),('legacy changed',eos,n)
  for gasbasis in ['STANDARD','LOCAL_VOLUME']:
   dest,item=run(eos+'_volume_'+gasbasis,dict(gas_eos=eos,liquid_flow_basis='LOCAL_VOLUME',gas_flow_basis=gasbasis));assert item['returncode']==0
   c=dict(s.split('=',1) for s in (dest/'inlet_contract.txt').read_text().splitlines());assert float(c['steady_inlet_contract_relative_error'])<=1e-11
   p=float(c['initial_inlet_face_pressure_Pa']);expected=(997.05+(p-101325)/1498**2)*.3
   assert abs(float(c['initial_liquid_mass_flux_kg_m2_s'])/expected-1)<1e-10
   if gasbasis=='LOCAL_VOLUME':
    cg2=343**2 if eos=='BAROTROPIC' else 287*293.15
    assert abs(float(c['initial_gas_mass_flux_kg_m2_s'])/(p/cg2*.6)-1)<1e-10
 for name,patch in [('liquid_conflict',dict(liquid_flow_basis='LOCAL_VOLUME',mul=1.)),('gas_conflict',dict(gas_flow_basis='LOCAL_VOLUME',mug=1.)),('missing_mass',dict(liquid_flow_basis='MASS'))]:
  dest,item=run(name,patch);assert item['returncode']!=0;assert 'inlet contract' in (dest/'stdout.log').read_text()
 (out/'manifest.json').write_text(json.dumps(records,indent=2)+'\n');print('PASS',len(records),'EOS/feed CLI executions')
if __name__=='__main__':main()
