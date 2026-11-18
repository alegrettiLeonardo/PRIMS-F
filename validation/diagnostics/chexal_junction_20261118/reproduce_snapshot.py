#!/usr/bin/env python3
"""Replay frozen SS1 Newton snapshot and test resulting traces with Chexal IDs5-7.

The full state/linearization retains Bhagwat. Chexal is evaluated only on its
horizontal proposed traces. This is not a Chexal transient or experiment.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,os,re,shlex,shutil,subprocess

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo',type=Path,default=Path(__file__).resolve().parents[3])
parser.add_argument('--output',type=Path,default=None)
args=parser.parse_args()
root=args.repo.resolve()
output=(args.output or root/'tmp/chexal_matched_snapshot_20261118').resolve()
output.mkdir(parents=True,exist_ok=True)
build=output/'build';build.mkdir(exist_ok=True)
original=root/'validation/diagnostics/bg_ss1_stall/probe_bg_ss1_stall.f90'
snapshot=root/'validation/diagnostics/bg_ss1_stall/last_rejected_state.dat'
text=original.read_text()
text=text.replace('integer::st,nroots','integer::st,nroots,chexal_model')
anchor="    write(*,'(a,i3,6es25.16)')'neighbor velocities ulL,ugL,ulR,ugR,alphaL,alphaR profile ',iprofile, &"
addition="""    do chexal_model=5,7
      astr=solve_alpha_from_superficials(jli,jga,rl,rg,0.0_wp,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
        cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=st,closure_id=chexal_model,nroots=nroots)
      write(*,'(a,4i5,es25.16)')'CHEXAL_REPLAY profile,model,status,nroots,jgSI ',iprofile,chexal_model,st,nroots,jga*cfg%omega_u
    enddo
"""
assert text.count(anchor)==1
text=text.replace(anchor,addition+anchor)
probe=output/'probe_matched_chexal.f90';probe.write_text(text)
shutil.copyfile(snapshot,output/'last_rejected_state.dat')
spec=importlib.util.spec_from_file_location('unit_runner',root/'scripts/run_unit_tests.py')
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
sources=[root/'scr'/f'{name}.f90' for name in runner.SOURCES]
compiler=shlex.split(os.environ.get('FC','gfortran'))
command=compiler+['-std=f2018','-O0','-g','-ffree-line-length-none','-fcheck=all,no-array-temps',
                 '-ffpe-trap=invalid,zero,overflow','-fopenmp',f'-J{build}',f'-I{build}',
                 *map(str,sources),str(probe),'-llapack','-lblas','-o',str(build/'probe')]
compile=subprocess.run(command,cwd=build,text=True,capture_output=True)
(output/'compile.log').write_text(compile.stdout+compile.stderr)
assert compile.returncode==0,compile.stderr[-1500:]
run=subprocess.run([str(build/'probe')],cwd=output,text=True,capture_output=True)
log=run.stdout+run.stderr
(output/'probe.log').write_text(log)
assert run.returncode==0,log[-1500:]
rows=[]
for line in log.splitlines():
    if line.startswith('CHEXAL_REPLAY'):
        values=line.split()[2:]
        rows.append(dict(zip(['profile','model','status','nroots'],map(int,values[:4])))|{'jg_m_s':float(values[4])})
negative=[r for r in rows if r['jg_m_s']<0]
assert negative and all(r['status']==61 and r['nroots']==0 for r in negative)
newton=[]
for line in log.splitlines():
    if line.startswith('Newton candidate lambda,point failures,junction status'):
        values=line.split()[-3:]
        newton.append({'lambda':float(values[0]),'point_failures':int(values[1]),'junction_status':int(values[2])})
assert len(newton)==5 and all(r['point_failures']==0 and r['junction_status']==1012 for r in newton)
summary={'status':'PASS','git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
         'scope':'Frozen Bhagwat SS1 Newton relinearization; horizontal trace replay through Chexal only',
         'compiler':subprocess.check_output(compiler+['--version'],text=True).splitlines()[0],
         'source_sha256':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*sources,original,snapshot]},
         'newton_candidates':newton,'chexal_trace_replays':rows,
         'negative_replays_count':len(negative),'negative_replays_all_domain61':True,
         'physical_reverse_gas_claim':False,'chexal_transient_claim':False}
(output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','chexal_trace_replays']},indent=2))
