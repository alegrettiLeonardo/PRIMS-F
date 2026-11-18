#!/usr/bin/env python3
"""Matched physical inputs; independent closure/mesh/time cohorts, no calibration."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--label',required=True)
    p.add_argument('--cases',nargs='+',default=['D_Stable','D_OSC','D_SS3','D_SS1'],
                   choices=['D_Stable','D_OSC','D_SS3','D_SS1'])
    p.add_argument('--closures',nargs='+',type=int,default=[3,5],choices=[3,5,6,7])
    p.add_argument('--cells',type=int,default=51)
    p.add_argument('--dt',type=float,default=.125)
    p.add_argument('--time',type=float,default=2400)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--timeout',type=float,default=86400)
    a=p.parse_args()
    if not 1<=a.workers<=8:p.error('workers must be 1..8')
    dest=ROOT/'validation'/'comparisons'/a.label
    dest.mkdir(parents=True,exist_ok=False)
    sources={str(q.relative_to(ROOT)):sha(q) for q in sorted((ROOT/'scr').glob('*.f90'))}
    plan={'started_utc':datetime.now(timezone.utc).isoformat(),
          'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
          'binary_sha256':sha(ROOT/'bin/prims_f'),'source_sha256':sources,
          'controls':vars(a),'attempts':[]}
    for case in a.cases:
        for closure in a.closures:
            command=[sys.executable,str(ROOT/'scripts/run_validation.py'),'--label',a.label,
                     '--case',case,'--closure',str(closure),'--cells',str(a.cells),
                     '--dt',str(a.dt),'--time',str(a.time),'--timeout',str(a.timeout),
                     '--eos','BAROTROPIC','--junction','MATCHED','--flux','ROE',
                     '--dissipation','LAGGED','--min-dt','1e-7']
            plan['attempts'].append({'case':case,'closure':closure,'command':command})
    (dest/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    def run(spec):
        log=dest/f"C{spec['closure']}_{spec['case']}.log"
        with log.open('w') as f:
            r=subprocess.run(spec['command'],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        return dict(case=spec['case'],closure=spec['closure'],returncode=r.returncode,
                    log=str(log.relative_to(ROOT)))
    results=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        pending=[pool.submit(run,spec) for spec in plan['attempts']]
        for task in as_completed(pending):
            result=task.result();results.append(result)
            (dest/'execution.json').write_text(json.dumps(results,indent=2)+'\n')
            print(json.dumps(result),flush=True)
    return 1 if any(r['returncode'] for r in results) else 0

if __name__=='__main__':raise SystemExit(main())
