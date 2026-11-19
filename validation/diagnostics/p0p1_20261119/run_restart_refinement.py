from pathlib import Path
import csv,json,hashlib,subprocess,os,time,concurrent.futures,shutil,sys
import numpy as np
import argparse
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--initial-fields',type=Path,required=True);args=p.parse_args()
b=Path(__file__).resolve().parents[3];out=args.output.resolve();out.mkdir();r=out.parent;v=b/'validation/diagnostics/p0p1_20261119'
source=r'''program restart_refinement
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::u(:,:)
 real(wp)::s,dtmax,duration,dummy(6)
 integer::n,i,k,iu,st
 character(512)::folder,input,arg
 call get_command_argument(1,folder);call get_command_argument(2,input)
 call get_command_argument(3,arg);read(arg,*)n
 call get_command_argument(4,arg);read(arg,*)dtmax
 call get_command_argument(5,arg);read(arg,*)duration
 cfg%omega_u=1.16_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325
 cfg%Lr=12.54625399319708_wp;cfg%Lp=57.4_wp;cfg%CA=3.01818193717446_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
 cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%beta=2*acos(-1._wp)/180
 cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
 cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=206000/cfg%omega_P
 cfg%mul=997.05_wp*.547_wp/(cfg%omega_rho*cfg%omega_u)
 cfg%mug=1.396675899947762_wp/(cfg%omega_rho*cfg%omega_u);cfg%njl=.547_wp/cfg%omega_u;cfg%njg=1.16_wp/cfg%omega_u
 cfg%gas_eos='BAROTROPIC';cfg%closure_pipe=1;cfg%closure_riser=3
 cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp;cfg%newton_max=25
 cfg%dissipation_mode='LAGGED';cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%jac_dissipation_fd=.true.
 cfg%nCells=2*n+4;cfg%N_pipe=n;cfg%num_dummy_cells=2;cfg%quad_points=3
 cfg%dt_start=.001_wp/cfg%omega_t;cfg%dt_min=1e-7_wp/cfg%omega_t;cfg%dt_max=dtmax/cfg%omega_t
 cfg%out_dir=trim(folder);cfg%output_interval_s=.125_wp;cfg%field_interval_s=1;cfg%checkpoint_interval_s=30
 allocate(u(cfg%nCells,3),cfg%vns(cfg%nCells),cfg%thetav(cfg%nCells),cfg%cell_width(cfg%nCells),cfg%face_positions(0:2*n))
 do i=0,n
  cfg%face_positions(i)=cfg%Lp/cfg%Lr*i/n
  cfg%face_positions(n+i)=cfg%Lp/cfg%Lr+real(i,wp)/n
 enddo
 do i=1,2*n
  cfg%vns(i+2)=.5_wp*(cfg%face_positions(i-1)+cfg%face_positions(i))
  cfg%cell_width(i+2)=cfg%face_positions(i)-cfg%face_positions(i-1)
  cfg%thetav(i+2)=-cfg%beta
  if(i>n)cfg%thetav(i+2)=atan((cfg%vns(i+2)*cfg%Lr-cfg%Lp)/cfg%CA)
 enddo
 open(newunit=iu,file=trim(input),status='old')
 do i=3,2*n+2
  read(iu,*)u(i,:)
 enddo
 close(iu)
 do i=1,2
  u(i,:)=u(3,:);u(2*n+2+i,:)=u(2*n+2,:)
  cfg%vns(i)=cfg%vns(3)-(3-i)*cfg%cell_width(3);cfg%vns(2*n+2+i)=cfg%vns(2*n+2)+i*cfg%cell_width(2*n+2)
  cfg%cell_width(i)=cfg%cell_width(3);cfg%cell_width(2*n+2+i)=cfg%cell_width(2*n+2)
  cfg%thetav(i)=cfg%thetav(3);cfg%thetav(2*n+2+i)=cfg%thetav(2*n+2)
 enddo
 call sim%init(cfg)
 call sim%simulate(u,cfg%vns,cfg%thetav,duration/cfg%omega_t,max_iter=10000000,status=st)
 if(st/=0)stop 3
end program
'''
(v/'restart_refinement.f90').write_text(source)
objs=list((b/'tmp/release').glob('*.o'));objs=[p for p in objs if p.name!='main.o']
cmd=['gfortran','-std=f2018','-O2','-ffree-line-length-none','-fopenmp','-I'+str(b/'tmp/release'),str(v/'restart_refinement.f90'),*map(str,objs),'-llapack','-lblas','-o',str(out/'restart')]
with (out/'build.log').open('w')as f:p=subprocess.run(cmd,cwd=out,stdout=f,stderr=subprocess.STDOUT)
assert p.returncode==0,(out/'build.log').read_text()[-2000:]
src=args.initial_fields.resolve()
with src.open()as f:fields=list(csv.DictReader(f))
fields=[x for x in fields if x['time_s']==fields[-1]['time_s']][-102:]
edges=np.r_[np.linspace(0,57.4,52),57.4+np.linspace(0,12.54625399319708,52)[1:]]
U=[]
for x in fields:
 a=float(x['alpha']);rl=float(x['rhol_kg_m3']);rg=float(x['rhog_kg_m3']);ul=float(x['ul_m_s']);ug=float(x['ug_m_s'])
 U.append([(1-a)*rl/997.05,a*rg/997.05,((1-a)*rl*ul+a*rg*ug)/(997.05*1.16)])
U=np.array(U);np.savetxt(out/'common_initial_U.txt',U,fmt='%.17e')
jobs=[(51,.125),(101,.125),(201,.125),(401,.125),(51,.0625),(51,.03125),(101,.0625)]
protocol=dict(operator=subprocess.check_output(['git','rev-parse','HEAD'],cwd=b,text=True).strip(),driver_sha256=hashlib.sha256((v/'restart_refinement.f90').read_bytes()).hexdigest(),binary_sha256=hashlib.sha256((out/'restart').read_bytes()).hexdigest(),module_objects={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in objs},initial_field_source=str(src),initial_field_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),initial_U_sha256=hashlib.sha256((out/'common_initial_U.txt').read_bytes()).hexdigest(),initial_definition='New initial-value campaign: conservative overlap projection of the same piecewise constant FV U reconstructed from 1365094 field at 2400s. Not a bitwise continuation; time starts at zero. Geometry, feed, EOS, solver and tolerances held fixed.',duration_s=600,settling_limit=.005,convergence_target=.01,last_complete_cycles=10,wall_cap_s=2400,jobs=jobs)
(out/'PROTOCOL.json').write_text(json.dumps(protocol,indent=2)+'\n')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
def run(job):
 n,dt=job;folder=out/f'N{n}_dt{dt}';folder.mkdir();newedges=np.r_[np.linspace(0,57.4,n+1),57.4+np.linspace(0,12.54625399319708,n+1)[1:]]
 new=[]
 for left,right in zip(newedges[:-1],newedges[1:]):
  overlap=np.maximum(0,np.minimum(edges[1:],right)-np.maximum(edges[:-1],left));new.append(overlap@U/(right-left))
 new=np.array(new);np.savetxt(folder/'initial_U.txt',new,fmt='%.17e');defect=np.diff(newedges)@new-np.diff(edges)@U
 assert np.max(np.abs(defect))<1e-12
 start=time.monotonic();timed=False
 with(folder/'stdout.log').open('w')as f:
  try:p=subprocess.run([str(out/'restart'),str(folder),str(folder/'initial_U.txt'),str(n),str(dt),'600'],cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=2400);code=p.returncode
  except subprocess.TimeoutExpired:timed=True;code=None
 status=(folder/'run_status.txt').read_text()if(folder/'run_status.txt').exists()else None
 manifest=dict(operator=protocol['operator'],n_pipe=n,n_riser=n,dt_max_s=dt,time_s=600,projection_defect_nd_integral=defect.tolist(),returncode=code,timeout=timed,wall_s=time.monotonic()-start,status=status,files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in folder.iterdir()if p.is_file()})
 (folder/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(folder.name,code,'TIMEOUT'if timed else'',flush=True);return manifest
with concurrent.futures.ThreadPoolExecutor(max_workers=2)as pool:results=list(pool.map(run,jobs))
(out/'manifest.json').write_text(json.dumps(dict(protocol=protocol,runs=results),indent=2)+'\n')
