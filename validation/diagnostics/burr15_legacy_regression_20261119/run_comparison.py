#!/usr/bin/env python3
"""Real short default-policy regression, isolated frozen worktrees, one worker."""
from pathlib import Path
import csv,hashlib,json,os,re,subprocess,time
repo=Path('/home/leonardoab/Documentos/PRIMS-F-burr15-20261119')
report=repo/'validation/diagnostics/burr15_legacy_regression_20261119';report.mkdir(parents=True,exist_ok=True)
versions=[('legacy','ee5810e'),('fixed_policy','8a85e37')]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_DYNAMIC='FALSE')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def run(cmd,cwd,log,timeout=180):
    start=time.monotonic();p=subprocess.run(cmd,cwd=cwd,env=env,capture_output=True,text=True,timeout=timeout)
    log.write_text(p.stdout+p.stderr)
    return dict(command=cmd,returncode=p.returncode,wall_seconds=time.monotonic()-start,log=str(log))
def kv(p):return {k:v.strip() for k,v in (l.split('=',1) for l in p.read_text().splitlines() if '=' in l)}
def set_value(text,key,value):
    pat=rf'(?i)\b{re.escape(key)}\s*=\s*(?:\x27[^\x27]*\x27|[^,\n/]+)'
    assert len(re.findall(pat,text))==1,(key,len(re.findall(pat,text)))
    return re.sub(pat,lambda _:f'{key}={value}',text)

worktrees={};builds={};runs={}
for label,commit in versions:
    work=repo.parent/f'PRIMS-F-burr15-regression-{commit}-20261119'
    if not work.exists():
        subprocess.run(['git','worktree','add','--detach',str(work),commit],cwd=repo,check=True,capture_output=True,text=True)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=work,text=True).strip();assert head.startswith(commit)
    status=subprocess.check_output(['git','status','--porcelain','--untracked-files=no'],cwd=work,text=True);assert not status
    build=run(['make','-j1','FC=gfortran','LIBS=-llapack -lblas'],work,report/(label+'_build_v2.log'))
    assert build['returncode']==0,build
    build.update(commit=head,binary_sha256=sha(work/'bin/prims_f'),
                 sources_sha256={str(p.relative_to(work)):sha(p) for p in (work/'scr').glob('*.f90')},
                 compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0])
    builds[label]=build;worktrees[label]=work

for case in ['D_Stable','D_OSC']:
    base=(worktrees['legacy']/'input/benchmarks'/(case+'.nml')).read_text()
    common=set_value(base,'time_dim','1.0')
    canonical=set_value(common,'out_dir',"'<OUTPUT>'")
    record={}
    for label,_ in versions:
        work=worktrees[label];output=work/'tmp/reg2'/case
        assert not output.exists(),str(output)
        deck=work/'tmp'/(case+'_default_policy_1s_v2.nml');deck.parent.mkdir(exist_ok=True)
        text=set_value(common,'out_dir',"'"+str(output)+"'");deck.write_text(text)
        assert set_value(text,'out_dir',"'<OUTPUT>'")==canonical
        result=run([str(work/'bin/prims_f'),str(deck)],work,report/(label+'_'+case+'_v2.log'),timeout=180)
        result.update(input_sha256=sha(deck),canonical_input_sha256=hashlib.sha256(canonical.encode()).hexdigest(),
                      output=str(output),run_status=kv(output/'run_status.txt') if (output/'run_status.txt').exists() else None)
        assert result['returncode']==0 and result['run_status']['status']=='0',result
        result['artifacts_sha256']={name:sha(output/name) for name in ['monitor.csv','fields.csv','results_'+case+'.dat','run_status.txt']}
        record[label]=result
    old=Path(record['legacy']['output']);new=Path(record['fixed_policy']['output'])
    equal={name:(old/name).read_bytes()==(new/name).read_bytes() for name in ['monitor.csv','fields.csv']}
    a=[list(map(float,l.split())) for l in (old/('results_'+case+'.dat')).read_text().splitlines() if not l.startswith('#')]
    b=[list(map(float,l.split())) for l in (new/('results_'+case+'.dat')).read_text().splitlines() if not l.startswith('#')]
    assert len(a)==len(b)==102
    sections={}
    for sec,rr in [('pipe',range(51)),('riser',range(51,102))]:
        sections[sec]={key:max(abs(a[i][j]-b[i][j]) for i in rr) for j,key in enumerate(['s_nd','s_m','alpha','P_Pa','ul_m_s','ug_m_s'])}
    record.update(canonical_inputs_equal=record['legacy']['canonical_input_sha256']==record['fixed_policy']['canonical_input_sha256'],
                  byte_equal=equal,legacy_output_max_absolute_difference=sections)
    runs[case]=record
    print(case,json.dumps(dict(byte_equal=equal,results_difference=sections)),flush=True)
    (report/'summary.json').write_text(json.dumps(dict(schema='PRIMS-F.default_pipe_policy_regression.v1',
        physical_time_s=1,workers=1,environment={k:env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_DYNAMIC']},
        process_nice=os.getpriority(os.PRIO_PROCESS,0),builds=builds,cases=runs),indent=2)+'\n')
assert all(all(r['byte_equal'].values()) for r in runs.values())
print('PASS: both complete monitor/fields are byte-identical; result differences quantified by section.',flush=True)
