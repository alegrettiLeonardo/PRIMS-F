#!/usr/bin/env python3
"""Run the declared catenary verification matrix, one numerical factor at a time.

24 attempts: four reference cases, three meshes at dt=0.5 s, three time-step
caps on the 51+51 mesh (the base run is shared), and four isothermal sensitivity
runs. A failed physical/numerical domain is recorded, never counted as a pass.
"""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--time',type=float,default=2400.)
    p.add_argument('--timeout',type=float,default=7200.)
    p.add_argument('--label',default='final')
    p.add_argument('--study',choices=['baseline','osc-refinement'],default='baseline')
    p.add_argument('--resume',action='store_true',help='Reuse finished attempts only when source and binary hashes agree')
    args=p.parse_args()
    if not 1<=args.workers<=8:p.error('workers must be between 1 and 8')
    sources={str(q.relative_to(ROOT)):sha(q) for q in sorted((ROOT/'scr').glob('*.f90'))}
    binary_hash=sha(ROOT/'bin/prims_f')
    plan=[]
    cases=('D_Stable','D_OSC','D_SS3','D_SS1') if args.study=='baseline' else ('D_OSC',)
    variants=((51,.5,'BAROTROPIC'),(101,.5,'BAROTROPIC'),(201,.5,'BAROTROPIC'),
              (51,.25,'BAROTROPIC'),(51,.125,'BAROTROPIC'),(51,.5,'ISOTHERMAL'))
    if args.study=='osc-refinement':
        variants=((51,.0625,'BAROTROPIC'),(51,.03125,'BAROTROPIC'),(101,.125,'BAROTROPIC'),(201,.125,'BAROTROPIC'))
    for case in cases:
        for n,dt,eos in variants:
            name=f'{args.label}_{case}_C3_{eos}_MATCHED_LAGGED_ROE_N{n}_{n}_dt{dt:g}'
            command=[sys.executable,str(ROOT/'scripts/run_validation.py'),'--case',case,'--cells',str(n),
                     '--dt',str(dt),'--time',str(args.time),'--timeout',str(args.timeout),'--label',args.label,
                     '--closure','3','--eos',eos,'--min-dt','1e-7']
            plan.append({'name':name,'command':command,'case':case,'cells_per_section':n,
                         'dt_max_s':dt,'gas_eos':eos,'duration_s':args.time})
    dest=ROOT/'validation';dest.mkdir(exist_ok=True)
    (dest/f'{args.label}_study_plan.json').write_text(json.dumps({'schema':'PRIMS-F.study.v1',
        'study':args.study,'binary_sha256':binary_hash,'source_sha256':sources,'attempts':plan},indent=2))

    def run(spec):
        path=ROOT/'validation/results'/spec['name']/'manifest.json'
        if path.exists():
            old=json.loads(path.read_text())
            identical=old.get('binary_sha256')==binary_hash and old.get('source_sha256')==sources
            finished=old.get('status') in ('COMPLETED','FAILED','TIMEOUT')
            same_duration=old.get('controls',{}).get('time')==args.time
            if args.resume and identical and finished and same_duration:
                return {'name':spec['name'],'action':'REUSED','status':old['status']}
            return {'name':spec['name'],'action':'REFUSED_EXISTING','status':old.get('status')}
        log=dest/(spec['name']+'.log')
        with log.open('w') as stream:
            result=subprocess.run(spec['command'],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
        manifest=json.loads(path.read_text()) if path.exists() else {}
        return {'name':spec['name'],'action':'EXECUTED','status':manifest.get('status','RUNNER_ERROR'),
                'returncode':result.returncode}

    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for result in pool.map(run,plan):
            results.append(result);print(json.dumps(result),flush=True)
    (dest/f'{args.label}_study_execution.json').write_text(json.dumps(results,indent=2))
    return 1 if any(r['status']!='COMPLETED' for r in results) else 0


if __name__=='__main__':raise SystemExit(main())
