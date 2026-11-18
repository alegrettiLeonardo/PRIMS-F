from pathlib import Path
import subprocess,os,sys,json,hashlib
here=Path(__file__).resolve().parent;root=here.parents[1]
label=sys.argv[1] if len(sys.argv)>1 else 'baseline'
src=here/'baseline_sources' if label=='baseline' else root/'scr'
dst=here/f'build_branch_{label}';dst.mkdir(exist_ok=True)
names='mod_kinds mod_rootfinding mod_void_helpers module_flow_regime module_config module_lapack_small module_closure_law module_find module_jacobian module_matrix_roe module_roe_solver'.split()
files=[src/f'{n}.f90' for n in names]+[here/'probe_ss3_branch.f90']
cmd=[os.environ.get('FC','gfortran'),'-std=f2018','-O2','-g','-fcheck=all,no-array-temps','-ffpe-trap=invalid,zero,overflow','-ffree-line-length-none','-ftrampoline-impl=heap','-J',str(dst),'-I',str(dst)]
cmd += list(map(str,files))+['-l:liblapack.so.3','-l:libblas.so.3','-o',str(dst/'probe_ss3_branch')]
p=subprocess.run(cmd,capture_output=True,text=True);(dst/'build.log').write_text(p.stdout+p.stderr);print('build',label,p.returncode)
if p.returncode:print(p.stderr);raise SystemExit(p.returncode)
(dst/'provenance.json').write_text(json.dumps({'source_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},'command':cmd},indent=2)+'\n')
p=subprocess.run([str(dst/'probe_ss3_branch'),label],cwd=here,capture_output=True,text=True)
(here/f'ss3_branch_{label}.log').write_text(p.stdout+p.stderr);print(p.stdout+p.stderr);raise SystemExit(p.returncode)
