"""Build isolated audit libraries; never edits production sources."""
from pathlib import Path
import argparse, hashlib, json, os, subprocess, sys, time

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
 ap.add_argument('--debug',action='store_true');args=ap.parse_args()
 root=args.root.resolve();out=args.out.resolve();out.mkdir(parents=True,exist_ok=True)
 sys.path.insert(0,str(root/'scripts'));from run_unit_tests import SOURCES
 mf={'base_sha':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
     'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
     'production_sources':{str(p.relative_to(root)):sha(p) for p in (root/'scr').glob('*.f90')},'commands':[]}
 flags=['-std=f2018','-O0' if args.debug else '-O2','-fPIC','-fopenmp','-ffree-line-length-none','-Wall','-Wextra']
 if args.debug:flags+=['-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow','-fbacktrace']
 here=Path(__file__).parent
 # All copies, objects and executables live outside version control.
 sources=[root/'scr'/(x+'.f90') for x in SOURCES]+[root/'validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90',here/'boundary_audit_api.f90']
 cmd=['gfortran',*flags,'-shared',*[str(x) for x in sources],'-llapack','-lblas','-o','audit.so']
 start=time.perf_counter();r=subprocess.run(cmd,cwd=out,text=True,capture_output=True)
 (out/'build.log').write_text(r.stdout+r.stderr)
 mf.update(command=cmd,returncode=r.returncode,elapsed_s=time.perf_counter()-start)
 if r.returncode:print((out/'build.log').read_text()[-5000:]);raise SystemExit(r.returncode)
 mf['harness_sources']={str(x):sha(x) for x in sources[-2:]};mf['library_sha256']=sha(out/'audit.so')
 (out/'manifest.json').write_text(json.dumps(mf,indent=2)+'\n')
 print(out/'manifest.json',mf['library_sha256'])
if __name__=='__main__':main()
