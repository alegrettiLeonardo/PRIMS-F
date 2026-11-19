from pathlib import Path
import subprocess,sys,hashlib,json,time,os,concurrent.futures,shutil
import argparse
p=argparse.ArgumentParser();p.add_argument('--frozen',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
b=Path(__file__).resolve().parents[3];frozen=args.frozen.resolve();r=args.output.resolve().parent
sys.path.insert(0,str(b/'tests'));from test_pipe_regime_cli import make_deck
out=args.output.resolve();out.mkdir();shutil.copy2(b/'docs/QUALIFICATION_P0P1.md',out/'PROTOCOL.md')
binary=frozen/'prims_f';binhash=hashlib.sha256(binary.read_bytes()).hexdigest()
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
jobs=[]
for case in ['D_SS1','D_SS3']:
 for mode in ['LEGACY','STEADY_TRACE']:
  jobs.append((case+'_'+mode,case,dict(boundary_method="'"+mode+"'",time_dim=110),180))
jobs.extend([('D_Stable','D_Stable',dict(time_dim=2400),900),('D_OSC','D_OSC',dict(time_dim=2400),900),('OSC_SOURCE_120','D_OSC',dict(time_dim=120,boundary_method="'STEADY_TRACE'"),180)])
for n in [8,16,32,64,201]:
 for g in [0,4.905,9.81,19.62]:
  jobs.append((f'HYDRO_N{n}_G{g}','D_OSC',dict(nodes=n,n_riser=n,G=g,time_dim=1,jl=0,jg=0,mul=0,mug=0,gas_eos="'ISOTHERMAL'",single_phase_mode="'REDUCED'",source_scheme="'HYDROSTATIC'",perturbation_amplitude=0,dt_max_s=.01,field_interval_s=.25),90))
build=json.loads((frozen/'build_manifest.json').read_text());operator=build.get('operator',build.get('sha'))
protocol=dict(operator=operator,binary_sha256=binhash,environment={k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},cases=[dict(name=n,benchmark=c,changes=d,wall_cap_s=t)for n,c,d,t in jobs])
(out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
def run(job):
 name,case,extra,limit=job;folder=out/name;folder.mkdir()
 common=dict(nodes=51,n_riser=51,time_dim=2400,dt_min=-1,dt_max=-1,dt_start=-1,dt_start_s=.001,dt_min_s=1e-7,dt_max_s=.125,tol_res=1e-10,tol_inc=1e-10,newton_max=25,field_interval_s=2,output_interval_s=.125,checkpoint_interval_s=10,max_iter=10000000,out_dir="'"+str(folder)+"'")
 common.update(extra)
 if name.startswith('HYDRO') and not operator.startswith('18aa7f7'):common['outlet_reservoir_phase']=1
 deck=make_deck((b/f'input/benchmarks/{case}.nml').read_text(),common);(folder/'input.nml').write_text(deck)
 start=time.monotonic();timed=False
 with (folder/'stdout.log').open('w')as log:
  try:p=subprocess.run([str(binary),str(folder/'input.nml')],cwd=b,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=limit);code=p.returncode
  except subprocess.TimeoutExpired:timed=True;code=None
 status=(folder/'run_status.txt').read_text()if(folder/'run_status.txt').exists()else None
 manifest=dict(name=name,operator=protocol['operator'],binary_sha256=binhash,command=[str(binary),str(folder/'input.nml')],returncode=code,timeout=timed,wall_s=time.monotonic()-start,status=status,files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in folder.iterdir()if p.is_file()})
 (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(name,code,'TIMEOUT'if timed else '',(status or '')[:100].replace('\n',' '),flush=True);return manifest
with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:
 results=list(pool.map(run,jobs))
(out/'manifest.json').write_text(json.dumps(dict(protocol=protocol,jobs=results),indent=2)+'\n')
