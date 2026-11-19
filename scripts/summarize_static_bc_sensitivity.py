#!/usr/bin/env python3
"""Index terminal BC sensitivities, hypotheses and common-window metrics."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read(path):return json.loads(Path(path).read_text())


def channel_metrics(rows):
    if not rows:return None
    if any(r.get('coverage_fraction')!=1.0 for r in rows):
        raise ValueError('Incomplete channel coverage cannot be summarized as complete cycles')
    duration=sum(r['covered_duration_s'] for r in rows)
    return dict(time_weighted_mean=sum(r['mean']*r['covered_duration_s'] for r in rows)/duration,
                mean_cycle_peak_to_peak=sum(r['peak_to_peak'] for r in rows)/len(rows),
                minimum=min(r['minimum'] for r in rows),maximum=max(r['maximum'] for r in rows))


def summarize(run, analysis, boundary=None):
    manifest_path=run/'manifest.json';manifest=read(manifest_path)
    cross_path=analysis/'cross_summary.json';cross=read(cross_path)['runs'][0]
    cons_path=analysis/'accepted_conservation.json';cons=read(cons_path)
    if cross['manifest']!=manifest:raise ValueError('Cycle summary belongs to a different run manifest')
    if cons['artifacts']['manifest.json']['sha256']!=sha(manifest_path):
        raise ValueError('Accepted-step summary belongs to a different manifest')
    channels={k:channel_metrics(v) for k,v in cross.get('cycle_channels',{}).items()}
    result=dict(name=manifest['name'],completed=cross['completed'],status=manifest['status'],
                reached_time_s=manifest['run_status']['reached_time_s'],
                accepted_steps=manifest['run_status']['accepted_steps'],
                rejected_steps=manifest['run_status']['rejected_steps'],
                cpu_seconds=manifest['run_status']['cpu_seconds'],wall_seconds=manifest['wall_seconds'],
                analysis_available=cross.get('analysis_available',False),
                cycles_available=cross.get('cycles_available',False),
                accommodated=cross.get('accommodated',False),
                cycle_edges_s=cross.get('cycle_edges_s',[]),pressure_column=cross.get('pressure_column'),
                period_mean_s=cross.get('period_mean_s'),period_cv=cross.get('period_cv'),
                cycle_channels=channels,
                arithmetic_contract_passes=cons['ledger']['arithmetic_contract_passes'],
                phase_region_balances=cons['ledger']['balances'],
                provenance={str(p):sha(p) for p in [manifest_path,cross_path,cons_path]})
    cycle_path=analysis/'cycle_conservation.json'
    if cycle_path.exists():
        cycle=read(cycle_path)
        if cycle['mass_summary']['sha256']!=sha(cons_path) or cycle['cycle_summary']['sha256']!=sha(cross_path):
            raise ValueError('Cycle ledger is linked to different analysis versions')
        intervals=cycle['cycles']
        cycle_maxima={}
        for key in cons['ledger']['balances']:
            values=[c['balances'][key] for c in intervals]
            cycle_maxima[key]=dict(
                max_relative_inventory=max((v['relative_inventory'] for v in values if v['relative_inventory'] is not None),default=None),
                max_relative_processed=max((v['relative_processed'] for v in values if v['relative_processed'] is not None),default=None),
                passes_inventory_1e6=all(v['passes_inventory_1e6'] is True for v in values) if values else None,
                passes_processed_1e6=all(v['passes_processed_1e6'] is True for v in values) if values else None)
        result['cycle_balance_audit']=dict(intervals=len(intervals),phase_regions=cycle_maxima)
        result['provenance'][str(cycle_path)]=sha(cycle_path)
    if boundary is not None:
        b=read(boundary)
        if b['provenance'].get(str(manifest_path))!=sha(manifest_path):
            raise ValueError('Boundary audit belongs to a different run manifest')
        buckets=b['boundary_summary'];frames=b['input_frames']
        expected={'1_1':2,'1_2':2,'2_1':1,'2_2':1}
        result['boundary_audit']=dict(input_frames=frames,
             four_buckets_complete=set(buckets)==set(expected) and all(v['samples']==frames for v in buckets.values()),
             all_status_zero=all(v['status_counts']=={'0':frames} for v in buckets.values()),
             expected_incoming_count=all(buckets.get(k,{}).get('incoming_counts')=={str(n):frames} for k,n in expected.items()),
             negative_ul_samples=sum(v['negative_ul_samples'] for v in buckets.values()),
             negative_ug_samples=sum(v['negative_ug_samples'] for v in buckets.values()))
        result['provenance'][str(boundary)]=sha(boundary)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--allow-partial',action='store_true')
    a=p.parse_args();plan=read(a.plan);plan_hash=sha(a.plan);base=a.plan.parent
    nominal=summarize(Path(plan['reference_manifest']['path']).parent,base/'analysis/nominal')
    rows=[];pending=[]
    for job in plan['jobs']:
        run=Path(job['output']);analysis=base/'analysis'/job['name']
        boundary=base/'boundary_audit'/job['name']/'summary.json'
        paths=[run/'manifest.json',analysis/'cross_summary.json',analysis/'accepted_conservation.json',analysis/'cycle_conservation.json',boundary]
        if not all(path.exists() for path in paths):pending.append(job['name']);continue
        m=read(run/'manifest.json')
        if m['status']=='RUNNING':pending.append(job['name']);continue
        if m['plan_sha256']!=plan_hash:raise ValueError('Mixed campaigns')
        r=summarize(run,analysis,boundary);r.update(axis=job['axis'],static_offset=job['static_offset'],controls=job['parameters'])
        r['relative_difference_from_nominal']={}
        for key,left,right in [('period_mean_s',nominal['period_mean_s'],r['period_mean_s'])]:
            r['relative_difference_from_nominal'][key]=(right-left)/abs(left) if left and right is not None else None
        for metric in ['time_weighted_mean','mean_cycle_peak_to_peak']:
            left=nominal['cycle_channels'].get('pressure',{}).get(metric)
            right=r['cycle_channels'].get('pressure',{}).get(metric)
            r['relative_difference_from_nominal']['pressure_'+metric]=(right-left)/abs(left) if left and right is not None else None
        phase=base/'phase_vs_nominal'/(job['name'].removeprefix('D_OSC_bc_')+'.json')
        if phase.exists():r['phase_comparison']=dict(path=str(phase),sha256=sha(phase))
        rows.append(r)
    if pending and not a.allow_partial:raise ValueError('Terminal analyses pending: '+', '.join(pending))
    result=dict(schema='PRIMS-F.static_bc_metrics.v1',complete_campaign=not pending,
                pending_runs=pending,plan_sha256=plan_hash,seal_sha256=sha(base/'seal.json'),
                summarizer_sha256=sha(__file__),nominal=nominal,runs=rows,
                hypotheses=['Offsets are static, one axis at a time; published signal deviations have no stated statistical definition.',
                            'Gas superficial flow remains STANDARD at 101300 Pa and 293.15 K; its mass override changes only with jg.',
                            'The same 69813ee operator, mesh, timestep controls, EOS, geometry and closure are retained.',
                            'Cycle channels use the same ten complete intervals and recorded coverage; pressure location is explicit.',
                            'Physical mass limits are tested separately from arithmetic identities and boundary counts.',
                            'No forcing frequency, probability distribution, physical regime validation or automatic switching is inferred.'])
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(terminal_analyses=len(rows),pending=len(pending),output=str(a.output))))


if __name__=='__main__':main()
