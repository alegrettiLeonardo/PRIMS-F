#!/usr/bin/env python3
"""Diagnose saved MATCHED/LAGGED ID3 or ID5 states without changing the solver.

Uses the frozen release objects from the cohort, and the existing generic
junction probe with its trace closure selector replaced by the case's selector.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,os,subprocess,shutil,csv

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--repo',type=Path,required=True)
parser.add_argument('--objects-repo',type=Path,help='Separate source of release objects; snapshot inputs remain frozen in --repo.')
parser.add_argument('--closure',type=int,choices=[3,5])
parser.add_argument('--label',default='comparison_20261118')
parser.add_argument('--tag',default='diagnostics_cell_detail')
parser.add_argument('--case',choices=['D_SS1','D_SS3','D_OSC','D_Stable'])
args=parser.parse_args();root=args.repo.resolve()
objects_root=(args.objects_repo or root).resolve()
evidence=root/'validation/comparisons'/args.label/args.tag
evidence.mkdir(exist_ok=True)
spec=importlib.util.spec_from_file_location('historical_case_probe',root/'validation/diagnostics/build_bg_ss1_case.py')
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
original=root/'validation/diagnostics/probe_bg_ss1_case.f90'
code=original.read_text()
assert code.count('status=st,closure_id=3,nroots=nroots')==1
code=code.replace('status=st,closure_id=3,nroots=nroots','status=st,closure_id=cfg%closure_riser,nroots=nroots')
old='          if(code/=0)pointfails=pointfails+1'
new='''          if(code/=0)then
            pointfails=pointfails+1
            write(*,'(a,3i8,es25.16)')'Newton cell failure step,cell,status,theta ',k,i,code,cfg%thetav(i)
          endif'''
assert code.count(old)==1
code=code.replace(old,new)
probe=evidence/'probe_selected_closure.f90';probe.write_text(code)
objects=sorted(p for p in (objects_root/'tmp/release').glob('*.o') if p.name!='main.o')
assert objects
command=['gfortran','-std=f2018','-O0','-g','-ffree-line-length-none','-fcheck=all,no-array-temps',
         '-ffpe-trap=invalid,zero,overflow','-fopenmp','-I'+str(objects_root/'tmp/release'),str(probe),
         *map(str,objects),'-llapack','-lblas','-o',str(evidence/'probe')]
compile=subprocess.run(command,cwd=root,text=True,capture_output=True)
(evidence/'compile.log').write_text(compile.stdout+compile.stderr)
assert compile.returncode==0,compile.stderr
common={'probe_sha256':hashlib.sha256(probe.read_bytes()).hexdigest(),
        'snapshot_repo':str(root),'objects_repo':str(objects_root),
        'objects_git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=objects_root,text=True).strip(),
        'objects_git_status':subprocess.check_output(['git','status','--short'],cwd=objects_root,text=True).strip(),
        'objects_source_sha256':{str(p.relative_to(objects_root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((objects_root/'scr').glob('*.f90'))},
        'release_object_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in objects},
        'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0]}
(evidence/'build_provenance.json').write_text(json.dumps(common,indent=2)+'\n')
for source in sorted((root/'validation/results').glob(args.label+'_*')):
    manifest=json.loads((source/'manifest.json').read_text())
    if args.case and manifest['controls']['case']!=args.case:continue
    if args.closure is not None and manifest['controls']['closure']!=args.closure:continue
    if manifest['status']!='FAILED' or not (source/'last_rejected_state.dat').exists():continue
    assert all(hashlib.sha256((root/k).read_bytes()).hexdigest()==v for k,v in manifest['source_sha256'].items())
    out=evidence/source.name;out.mkdir(exist_ok=False)
    frozen={}
    for name in ['last_rejected_state.dat','run_metadata.txt','run_status.txt','manifest.json']:
        shutil.copyfile(source/name,out/name)
        frozen[name]=hashlib.sha256((out/name).read_bytes()).hexdigest()
    p,scales=helper.metadata(out/'run_metadata.txt')
    lr,ca,cg,gl,gg,wp,wr,wu=scales
    assert int(p['closure_riser']) in (3,5) and p['junction_model']=='MATCHED' and p['dissipation_mode']=='LAGGED'
    values=[f"{int(p['nodes'])} {int(p['n_riser'])} {int(p['num_dummy_cells'])}",
            ' '.join(map(repr,[wu,wr,wp,lr,p['lp']])),
            ' '.join(map(repr,[p[k] for k in ['d','eps','g','mu_l','mu_g','sigma']])),
            ' '.join(map(repr,[p['cl']*math.sqrt(wr/wp),cg*math.sqrt(wr/wp),p['rho_l0']/wr,p['p_l0']/wp,p['ps']/wp])),
            ' '.join(map(repr,[gl/(wr*wu),gg/(wr*wu),p['jl']/wu,p['jg']/wu])),
            ' '.join(map(repr,[p[k] for k in ['tol','tola','tol_res','tol_inc']])),
            f"1 {int(p['closure_riser'])} {int(p['quad_points'])} 1000",
            'T' if p['jac_dissipation_fd'] else 'F',p['junction_model'],p['flux_scheme'],p['dissipation_mode']]
    (out/'case_config.dat').write_text('\n'.join(values)+'\n')
    run=subprocess.run([str(evidence/'probe')],cwd=out,text=True,capture_output=True,timeout=300)
    (out/'probe.log').write_text(run.stdout+run.stderr)
    summary=helper.summarize(run.stdout,p,source,frozen)
    summary['probe_returncode']=run.returncode
    summary['closure_riser']=int(p['closure_riser'])
    summary['cohort_source_hashes_verified']=True
    summary['newton_cell_failures']=[]
    for line in run.stdout.splitlines():
      if line.startswith('Newton cell failure'):
        f=line.split()[-4:]
        summary['newton_cell_failures'].append({'step_halving_index':int(f[0]),'cell':int(f[1]),'point_state_status':int(f[2]),'theta_rad':float(f[3])})
    summary['profiles_minima']={}
    point_states=out/'point_states.csv'
    if point_states.exists():
      rows=list(csv.DictReader(point_states.open()))
      for profile in (1,2):
        valid=[r for r in rows if int(r['profile'])==profile and int(r['status'])==0]
        if valid:
          summary['profiles_minima'][str(profile)]={'valid_cells':len(valid),
            'min_ul_m_s':min(float(r['ul_m_s']) for r in valid),
            'min_ug_m_s':min(float(r['ug_m_s']) for r in valid)}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k in ['case','time_s','saved_status','closure_riser','probe_returncode','newton_cell_failures','profiles_minima']}),flush=True)
