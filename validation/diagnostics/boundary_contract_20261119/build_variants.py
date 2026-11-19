"""Generate narrowly modified copies OUTSIDE the source tree; no production edits."""
from pathlib import Path
import argparse, hashlib, json, subprocess, sys, time
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def once(text,old,new):
 assert text.count(old)==1,(old,text.count(old))
 return text.replace(old,new)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
 ap.add_argument('--variants',nargs='+',default=['plain','standard','volume','partition'],choices=['plain','standard','volume','partition']);ap.add_argument('--debug',action='store_true');args=ap.parse_args()
 root=args.root.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True);here=Path(__file__).parent
 sys.path.insert(0,str(root/'scripts'));from run_unit_tests import SOURCES
 for variant in args.variants:
  build=out/variant;build.mkdir(exist_ok=True)
  mf={'variant':variant,'sha':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'original_sources':{},'generated_sources':{}}
  for name in SOURCES:
   source=root/'scr'/(name+'.f90');text=source.read_text();mf['original_sources'][source.name]=sha(source)
   if name=='module_boundary_conditions' and variant=='volume':
    # One model change only: preserve Q_l rather than the reference-density mass flow.
    start=text.index('  subroutine boundary_flux(');end=text.index('  subroutine apply_boundary_conditions(')
    before,body,after=text[:start],text[start:end],text[end:]
    body=body.replace('cfg%mul','ml_prescribed')
    body=once(body,'real(wp)::p,rhol','real(wp)::ml_prescribed\n    real(wp)::p,rhol')
    body=once(body,'    if(is_inlet)then','    ml_prescribed=cfg%mul\n    if(is_inlet)then')
    body=once(body,'      alpha_interior=Uinterior(2)/rhog','      ml_prescribed=rhol*cfg%njl\n      alpha_interior=Uinterior(2)/rhog')
    text=before+body+after
   if name=='module_simulate_pipe_riser':
    if variant=='partition':
     text=once(text,'    match=match.and.(trim(this%cfg%junction_model)', '    match=match.and.(id_left/=id_right.or.thetaL/=thetaR)\n    match=match.and.(trim(this%cfg%junction_model)')
    if variant!='plain':
     text=once(text,'  use mod_kinds, only: wp,i4','  use audit_observer, only:audit_record\n  use mod_kinds, only: wp,i4')
     text=once(text,'    call write_mass_record(0.0_wp,0.0_wp,initial_flux_status)','    call write_mass_record(0.0_wp,0.0_wp,initial_flux_status)\n    if(initial_flux_status==0)call audit_record(this%cfg,U,0,0.0_wp,0.0_wp,fin,fout)')
     text=once(text,'      call write_mass_record(dt*this%cfg%omega_t,norm,0)','      call write_mass_record(dt*this%cfg%omega_t,norm,0)\n      call audit_record(this%cfg,U,step,tcur*this%cfg%omega_t,dt*this%cfg%omega_t,fin,fout)')
   target=build/source.name;target.write_text(text);mf['generated_sources'][target.name]=sha(target)
  (build/'audit_observer.f90').write_bytes((here/'audit_observer.f90').read_bytes())
  (build/'window_main.f90').write_bytes((here/'window_main.f90').read_bytes())
  names=[x+'.f90' for x in SOURCES]
  if variant!='plain':names.insert(names.index('module_simulate_pipe_riser.f90'),'audit_observer.f90')
  flags=['-std=f2018','-O0' if args.debug else '-O2','-fPIC','-fopenmp','-ffree-line-length-none']
  if args.debug:flags+=['-g','-fcheck=all','-fbacktrace','-ffpe-trap=invalid,zero,overflow']
  cmd=['gfortran',*flags,*names,'window_main.f90','-llapack','-lblas','-o','window']
  st=time.perf_counter();r=subprocess.run(cmd,cwd=build,text=True,capture_output=True);(build/'build.log').write_text(r.stdout+r.stderr)
  if r.returncode:print(r.stderr[-4000:]);raise SystemExit(r.returncode)
  mf.update(command=cmd,returncode=r.returncode,wall_s=time.perf_counter()-st,binary_sha256=sha(build/'window'),observer_sha256=sha(build/'audit_observer.f90'),driver_sha256=sha(build/'window_main.f90'),generator_sha256=sha(__file__))
  # Same generated variant also provides a C ABI for exact snapshot/operator comparisons.
  api=root/'validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90'
  cmd=['gfortran',*flags,'-shared',*names,str(api),str(here/'boundary_audit_api.f90'),'-llapack','-lblas','-o','audit.so']
  r=subprocess.run(cmd,cwd=build,text=True,capture_output=True);(build/'library.log').write_text(r.stdout+r.stderr)
  if r.returncode:print(r.stderr[-4000:]);raise SystemExit(r.returncode)
  mf['library_command']=cmd;mf['library_sha256']=sha(build/'audit.so')
  (build/'manifest.json').write_text(json.dumps(mf,indent=2)+'\n');print(variant,mf['binary_sha256'],flush=True)
if __name__=='__main__':main()
