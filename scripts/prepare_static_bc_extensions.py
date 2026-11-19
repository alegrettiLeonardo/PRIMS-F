#!/usr/bin/env python3
"""Prepare duration-only reruns after every initial BC case has been audited."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
from prepare_static_bc_sensitivity import FREEZE,HARNESS_SHA256
from summarize_static_bc_sensitivity import summarize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harness',type=Path,required=True)
    parser.add_argument('--reference-plan',type=Path,required=True)
    parser.add_argument('--duration-s',type=float,default=7200.)
    parser.add_argument('--workspace',type=Path,required=True)
    args=parser.parse_args()
    import hashlib
    if hashlib.sha256(args.harness.read_bytes()).hexdigest()!=HARNESS_SHA256:
        raise ValueError('Executor identity differs')
    spec=importlib.util.spec_from_file_location('sealed_executor',args.harness)
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    old,seal=h.verify_seal(args.reference_plan.resolve())
    if seal['git_head']!=FREEZE:raise ValueError('Operator freeze differs')
    if not math.isfinite(args.duration_s) or args.duration_s<=old['duration_s']:
        raise ValueError('Extension duration must exceed the previous duration')
    root=Path(old['root']);workspace=args.workspace.resolve()
    label='bc_static_ext'+format(args.duration_s,'.0f');target=workspace/label
    if target.exists():raise ValueError('Preserve prior evidence: extension target exists')
    if root==workspace or root in workspace.parents:raise ValueError('Results must stay external')
    selected=[];inputs=dict(old['input_sha256']);evidence=[]
    for job in old['jobs']:
        directory=Path(job['output']);mp=directory/'manifest.json'
        analysis=args.reference_plan.parent/'analysis'/job['name']/'summary.json'
        boundary=args.reference_plan.parent/'boundary_audit'/job['name']/'summary.json'
        cycle=args.reference_plan.parent/'analysis'/job['name']/'cycle_conservation.json'
        if not all(p.exists() for p in (mp,analysis,boundary,cycle)):
            raise ValueError('All initial terminal analyses must exist before extensions: '+job['name'])
        m=json.loads(mp.read_text());a=json.loads(analysis.read_text())
        if m['status'] not in ('COMPLETED','FAILED') or a['manifest_sha256']!=h.sha(mp):
            raise ValueError('Terminal analysis identity differs')
        if m.get('plan_sha256')!=h.sha(args.reference_plan) or m.get('seal_sha256')!=h.sha(args.reference_plan.parent/'seal.json'):
            raise ValueError('Reference run belongs to a different plan/seal')
        if m.get('git_head')!=FREEZE:raise ValueError('Reference run git freeze differs')
        for label in ('cross_summary','accepted_conservation'):
            item=a['analyses'][label];source=analysis.parent/(label+'.json')
            if Path(item['path'])!=source or item['sha256']!=h.sha(source):
                raise ValueError('Terminal analysis content identity differs')
        verified=summarize(directory,analysis.parent,boundary)
        if verified['completed']!=a['completed'] or verified['accommodated']!=a['accommodated']:
            raise ValueError('Terminal summary flags disagree with audited analyses')
        if m['binary_sha256']!=seal['binary_sha256'] or m['source_sha256']!=seal['source_sha256']:
            raise ValueError('Reference run freeze differs')
        if not m['binary_unchanged'] or not m['sources_unchanged']:
            raise ValueError('Reference freeze changed during a run')
        evidence.append(dict(name=job['name'],manifest_sha256=h.sha(mp),analysis_sha256=h.sha(analysis),
                             completed=a['completed'],accommodated=a['accommodated']))
        if a['completed'] and not a['accommodated']:
            selected.append(job)
            for p in (mp,analysis,boundary,cycle):inputs[str(p)]=h.sha(p)
    if not selected:
        print(json.dumps(dict(jobs=0,reason='No completed unaccommodated cases')));return
    jobs=[];decks={}
    for reference in selected:
        original=Path(reference['deck'])
        if h.sha(original)!=reference['deck_sha256']:raise ValueError('Original deck changed')
        text0=original.read_text();before=h.parameters(text0)
        name=reference['name']+'_T'+format(args.duration_s,'.0f')
        output=target/'runs'/name;deck=target/'decks'/(name+'.nml');relative=os.path.relpath(output,root)
        text=h.assign(h.assign(text0,'time_dim',format(args.duration_s,'.17g')),'out_dir',repr(relative))
        after=h.parameters(text);changes=h.parameter_changes(before,after)
        if changes!=['out_dir','time_dim']:raise ValueError('Extension changed a different control')
        if any(after[k]!=before[k] for k in h.PHYSICS):raise ValueError('Physical parameter changed')
        decks[deck]=text
        jobs.append(dict(name=name,case='D_OSC',origin_name=reference['name'],axis=reference['axis'],
                         static_offset=reference['static_offset'],pipe_regime_mode=reference['pipe_regime_mode'],
                         deck=str(deck),output=str(output),output_relative_to_root=relative,
                         parameters=after,physics=reference['physics'],reference_parameter_changes=changes,
                         previous_run=reference['output'],previous_duration_s=old['duration_s']))
        inputs[str(original)]=h.sha(original)
    for path,text in decks.items():
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x') as stream:stream.write(text)
    for job in jobs:job['deck_sha256']=h.sha(job['deck'])
    plan=dict(schema='PRIMS-F.static_bc_extension_plan.v1',state='PREPARED_UNSEALED',created_utc=h.now(),
              root=str(root),workspace=str(workspace),label=label,prepared_git_head=FREEZE,
              input_sha256=inputs,jobs=jobs,environment=h.ENV,duration_s=args.duration_s,
              maximum_total_solvers=4,maximum_owned_workers=2,runner_sha256=h.sha(args.harness),preparer_sha256=h.sha(__file__),
              reference_plan=dict(path=str(args.reference_plan.resolve()),sha256=h.sha(args.reference_plan)),
              reference_manifest=old['reference_manifest'],reference_deck=old['reference_deck'],
              source=old['source'],criteria=old['criteria'],previous_cases=evidence,
              execution_gate='Parent authorized duration-only extensions after six initial cases; at most two owned workers, one thread each, four solvers globally.',
              interpretation='Fresh integration from the same initial prescription, not a checkpoint restart; no forcing or parameter changes.')
    h.dump_new(target/'plan.json',plan)
    print(json.dumps(dict(plan=str(target/'plan.json'),jobs=len(jobs),duration_s=args.duration_s,solvers_started=0)))


if __name__=='__main__':main()
