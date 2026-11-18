#!/usr/bin/env python3
"""Reproducible independent mesh/time studies; no physical parameter calibration.

Examples:
  python3 scripts/run_validation.py --case D_SS3 --cells 51 --dt 0.5 --time 2400
  python3 scripts/run_validation.py --case D_SS3 --cells 101 --dt 0.5 --time 2400
  python3 scripts/run_validation.py --case D_SS3 --cells 51 --dt 0.25 --time 2400
Numerical and model variants are explicit in every case name and manifest.
"""
from __future__ import annotations
import argparse,hashlib,json,re,subprocess,time
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def assign(text,key,value):
    pattern=r'(\b'+re.escape(key)+r'\s*=)\s*(\'[^\']*\'|[^,\n]+)'
    if re.search(pattern,text,re.I):
        return re.sub(pattern,lambda m:m[1]+value,text,flags=re.I)
    return re.sub(r'(?m)^\s*/\s*$',f'  {key}={value}\n/',text)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',choices=['D_Stable','D_OSC','D_SS3','D_SS1'],required=True)
    p.add_argument('--cells',type=int,default=51,help='Cells in EACH section; nphysical=2*cells')
    p.add_argument('--riser-cells',type=int)
    p.add_argument('--dt',type=float,default=.5,help='Maximum physical time step, seconds')
    p.add_argument('--time',type=float,default=2400)
    p.add_argument('--closure',type=int,default=3)
    p.add_argument('--min-dt',type=float,default=1e-7,help='Failure floor for physical time step, seconds')
    p.add_argument('--eos',choices=['BAROTROPIC','ISOTHERMAL'],default='BAROTROPIC')
    p.add_argument('--junction',choices=['MATCHED','BLOCKING','CONSERVATIVE'],default='MATCHED')
    p.add_argument('--dissipation',choices=['LAGGED','FULLY_IMPLICIT'],default='LAGGED')
    p.add_argument('--flux',choices=['ROE','RUSANOV'],default='ROE')
    p.add_argument('--perturbation',type=float,default=.01)
    p.add_argument('--timeout',type=float,default=1800,help='Wall-clock safety limit per run, seconds')
    p.add_argument('--label',default='final')
    p.add_argument('--exe',type=Path,default=ROOT/'bin/prims_f')
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--prepare-only',action='store_true')
    args=p.parse_args()
    nr=args.riser_cells or args.cells
    if min(args.cells,nr)<2 or min(args.dt,args.time,args.timeout,args.min_dt)<=0 or args.min_dt>=args.dt:
        p.error('positive duration/time step and at least two cells in each section are required')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',args.label): p.error('label must contain letters, digits, _ or -')
    name=f'{args.label}_{args.case}_C{args.closure}_{args.eos}_{args.junction}_{args.dissipation}_{args.flux}_N{args.cells}_{nr}_dt{args.dt:g}'
    out=ROOT/'validation/results'/name
    if (out/'manifest.json').exists() and not args.overwrite:
        p.error(f'{out} already exists; choose a new label or explicitly --overwrite')
    out.mkdir(parents=True,exist_ok=True)
    decks=ROOT/'validation/cases';decks.mkdir(parents=True,exist_ok=True)
    deck=decks/(name+'.nml')
    source=ROOT/'input/benchmarks'/(args.case+'.nml')
    text=source.read_text()
    controls={'nodes':args.cells,'n_riser':nr,'time_dim':args.time,'dt_start_s':min(.001,args.dt),
              'dt_max_s':args.dt,'dt_min_s':args.min_dt,'dt_min':-1,'tol':1e-12,'tola':1e-13,'tol_res':1e-8,'tol_inc':1e-8,
              'closure_riser':args.closure,'gas_eos':repr(args.eos),'junction_model':repr(args.junction),
              'dissipation_mode':repr(args.dissipation),'flux_scheme':repr(args.flux),
              'perturbation_amplitude':args.perturbation,'output_interval_s':min(.1,args.dt),
              'field_interval_s':15.,'use_muscl':'.false.','save_every':0,
              'out_dir':repr(str(out.relative_to(ROOT)))}
    for k,v in controls.items(): text=assign(text,k,str(v))
    text=text.rstrip()+"\n"
    deck.write_text(text)
    # Reject the malformed quoted-path-plus-slash construction seen in an early pilot generator.
    if re.search(r"'[^']*'/[^\n]*'",text): raise RuntimeError('Malformed namelist path')
    if args.prepare_only:
        print(deck);return 0
    exe=args.exe.resolve()
    if not exe.is_file():p.error(f'Build the executable first: {exe}')
    command=[str(exe),str(deck)]
    sources={str(q.relative_to(ROOT)):sha(q) for q in sorted((ROOT/'scr').glob('*.f90'))}
    manifest={'name':name,'started_utc':datetime.now(timezone.utc).isoformat(),'command':command,
              'binary_sha256':sha(exe),'source_sha256':sources,'namelist_sha256':sha(deck),
              'reference_sha256':sha(ROOT/'docs/benchmark_sources.json'),
              'controls':vars(args)|{'exe':str(exe),'riser_cells':nr},'status':'RUNNING'}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    start=time.monotonic()
    with (out/'stdout.log').open('w') as log:
        try:
            run=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
            manifest['returncode']=run.returncode
            manifest['status']='COMPLETED' if run.returncode==0 else 'FAILED'
        except subprocess.TimeoutExpired:
            manifest['returncode']=None;manifest['status']='TIMEOUT'
            log.write('\nWall-clock limit reached. Simulation NOT completed.\n')
    manifest['wall_seconds']=time.monotonic()-start
    manifest['ended_utc']=datetime.now(timezone.utc).isoformat()
    manifest['sources_unchanged']=all((ROOT/k).is_file() and sha(ROOT/k)==v for k,v in sources.items())
    # The executable already loaded by the OS is immutable; a later rebuild
    # cannot change the completed run, but source drift must remain visible.
    manifest['binary_unchanged']=sha(exe)==manifest['binary_sha256']
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps({k:manifest[k] for k in ['name','status','returncode','wall_seconds','sources_unchanged']}))
    return 0 if manifest['status']=='COMPLETED' else 1

if __name__=='__main__':raise SystemExit(main())
