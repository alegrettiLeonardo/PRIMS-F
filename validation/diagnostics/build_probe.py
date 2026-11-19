from pathlib import Path
import os,subprocess,hashlib,json
here=Path(__file__).resolve().parent
root=here.parents[1]
dst=here/'build';dst.mkdir(exist_ok=True)
names='mod_kinds mod_rootfinding mod_void_helpers module_flow_regime module_config module_lapack_small module_closure_law module_find module_catenary module_steady_state module_jacobian module_matrix_roe module_roe_solver module_boundary_conditions module_solver_utils module_time_step module_simulate_pipe_riser'.split()
files=[root/'scr'/f'{n}.f90' for n in names]+[here/'probe_ss3_stall.f90']
cmd=[os.environ.get('FC','gfortran'),'-std=f2018','-O2','-g','-fcheck=all,no-array-temps','-ffpe-trap=invalid,zero,overflow','-ffree-line-length-none','-ftrampoline-impl=heap','-J',str(dst),'-I',str(dst)]
cmd += list(map(str,files))+['-l:liblapack.so.3','-l:libblas.so.3','-o',str(dst/'probe_ss3_stall')]
p=subprocess.run(cmd,capture_output=True,text=True);(dst/'build.log').write_text(p.stdout+p.stderr)
print('build',p.returncode)
if p.returncode: print(p.stderr);raise SystemExit(p.returncode)
meta={'source_sha256':{str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},'compiler_command':cmd}
(dst/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n')
p=subprocess.run([str(dst/'probe_ss3_stall')],cwd=here,capture_output=True,text=True)
(here/'ss3_stall_probe.log').write_text(p.stdout);(here/'build/runtime_stderr.log').write_text(p.stderr);print(p.stdout);print('runtime stderr lines',len(p.stderr.splitlines()));raise SystemExit(p.returncode)
