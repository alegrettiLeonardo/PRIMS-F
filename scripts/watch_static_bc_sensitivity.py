#!/usr/bin/env python3
"""Analyze terminal static-BC runs with sealed common analysis tools; no solver."""
import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--analysis-tools',type=Path,required=True)
    parser.add_argument('--once',action='store_true')
    args=parser.parse_args()
    plan=json.loads(args.plan.read_text());plan_hash=sha(args.plan)
    seal_path=args.plan.parent/'seal.json';seal=json.loads(seal_path.read_text());seal_hash=sha(seal_path)
    tools_manifest=args.analysis_tools/'manifest.json';tool_hash=sha(tools_manifest)
    tools=json.loads(tools_manifest.read_text())
    output=args.plan.parent/'analysis'
    output.mkdir(exist_ok=True)
    while True:
        if sha(args.plan)!=plan_hash or sha(seal_path)!=seal_hash or sha(tools_manifest)!=tool_hash:
            raise ValueError('Plan, seal or analysis-tool manifest changed')
        for name,info in tools.items():
            if sha(args.analysis_tools/name)!=info['sha256']:raise ValueError('Analysis tool changed: '+name)
        terminal=[]
        for job in plan['jobs']:
            root=Path(job['output']);mp=root/'manifest.json'
            if not mp.exists():continue
            manifest=json.loads(mp.read_text())
            if manifest['status'] not in ('COMPLETED','FAILED'):continue
            if manifest['plan_sha256']!=plan_hash or manifest['seal_sha256']!=seal_hash:
                raise ValueError('Run belongs to a different plan/seal')
            if manifest['binary_sha256']!=seal['binary_sha256'] or manifest['source_sha256']!=seal['source_sha256']:
                raise ValueError('Run freeze differs')
            if not manifest['binary_unchanged'] or not manifest['sources_unchanged']:
                raise ValueError('Run freeze was modified during execution')
            dest=output/job['name'];dest.mkdir(exist_ok=True)
            results={}
            for script,label in [('analyze_accepted_conservation.py','accepted_conservation'),
                                 ('analyze_qualification_cross.py','cross_summary')]:
                target=dest/(label+'.json')
                if not target.exists():
                    command=[sys.executable,str(args.analysis_tools/script),'--run-dir',str(root),'--output',str(target)]
                    proc=subprocess.run(command,capture_output=True,text=True)
                    (dest/(label+'.log')).write_text(proc.stdout+'\n'+proc.stderr)
                    if not target.exists():raise RuntimeError(f'{script} failed for {job["name"]}: {proc.returncode}')
                results[label]=dict(path=str(target),sha256=sha(target))
            cons=json.loads(Path(results['accepted_conservation']['path']).read_text())
            cross=json.loads(Path(results['cross_summary']['path']).read_text())['runs'][0]
            # Strict physical balances come from accepted-step P2, not the
            # cycle analyzer's historical absolute floor in its own ledger.
            summary=dict(name=job['name'],axis=job['axis'],static_offset=job['static_offset'],
                         manifest_sha256=sha(mp),status=manifest['status'],
                         completed=cross['completed'],analysis_available=cross.get('analysis_available',False),
                         analysis_reason=cross.get('reason'),cycles_available=cross.get('cycles_available',False),
                         accommodated=cross.get('accommodated',False),
                         needs_extension=bool(cross['completed'] and not cross.get('accommodated',False)),
                         arithmetic_contract_passes=cons['ledger']['arithmetic_contract_passes'],
                         balances=cons['ledger']['balances'],analyses=results)
            record=dest/'summary.json'
            if not record.exists():
                with record.open('x') as stream:json.dump(summary,stream,indent=2,allow_nan=False);stream.write('\n')
                print(json.dumps({k:summary[k] for k in ['name','status','cycles_available','accommodated','needs_extension']}),flush=True)
            elif json.loads(record.read_text())!=summary:
                raise ValueError('Existing terminal summary identity differs')
            terminal.append(summary)
        if len(terminal)==len(plan['jobs']):
            result=dict(schema='PRIMS-F.static_bc_analysis.v1',plan_sha256=plan_hash,
                        seal_sha256=seal_hash,analysis_tools_manifest_sha256=tool_hash,
                        watcher_sha256=sha(__file__),runs=terminal,
                        limits=['Static sensitivities are not experimental time signals or a probability sample.',
                                'Absence of resolved cycles does not prove a physical stable regime.',
                                'Extension flags preserve existing results; this watcher never launches solvers.'])
            final=output/'aggregate.json'
            if not final.exists():
                with final.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
            return
        if args.once:
            print(json.dumps(dict(terminal=len(terminal),expected=len(plan['jobs']))));return
        time.sleep(60)


if __name__=='__main__':main()
