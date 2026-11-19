#!/usr/bin/env python3
"""Supplementary fixed-junction holdup; preserve the thirteen primary gates."""
import argparse,bisect,hashlib,json,math,statistics
from pathlib import Path
TOOLS=Path(__file__).resolve().parent
from cycle_waveform import difference,waveform
from validation_evidence import load_monitor,time_statistics

def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def calculate(data,edges):
    if len(edges)!=11 or any(not math.isfinite(x) for x in edges) or any(b<=a for a,b in zip(edges,edges[1:])):
        raise ValueError('Exactly ten increasing complete intervals required')
    times=data['time_s'];alpha=data['alpha_riser_trace']
    if len(times)!=len(alpha) or edges[0]<times[0] or edges[-1]>times[-1]:
        raise ValueError('Trace coverage incomplete')
    if any(b<=a for a,b in zip(times,times[1:])):raise ValueError('Trace time must advance')
    left=max(0,bisect.bisect_right(times,edges[0])-1)
    right=min(len(times),bisect.bisect_left(times,edges[-1])+1)
    if any(x is None or not math.isfinite(x) or not 0<=x<=1 for x in alpha[left:right]):
        raise ValueError('Missing or nonphysical trace fraction in declared window')
    derived=dict(data);derived['holdup_fixed_trace']=[None if x is None else 1-x for x in alpha]
    cycles=[time_statistics(derived,'holdup_fixed_trace',a,b) for a,b in zip(edges,edges[1:])]
    if any(x['coverage_fraction']<1-1e-12 for x in cycles):raise ValueError('A cycle lacks complete trace coverage')
    checks={key:difference(statistics.fmean(x[key] for x in cycles[:5]),
                           statistics.fmean(x[key] for x in cycles[5:]),1e-6)
            for key in ('mean','peak_to_peak')}
    return dict(cycles=cycles,cycle_mean_statistics={key:statistics.fmean(x[key] for x in cycles)
                for key in ('mean','minimum','maximum','peak_to_peak')},
                auxiliary_accommodation=checks,
                auxiliary_accommodated=all(x['passes_1percent_plus_floor'] for x in checks.values()),
                waveform=waveform(derived,'holdup_fixed_trace',edges,256))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--cycle-summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('Preserve previous evidence; choose a new output')
    folder=a.run_dir.resolve();summary=json.loads(a.cycle_summary.read_text())
    matched=[r for r in summary['runs'] if r['run']==folder.name]
    if len(matched)!=1:raise ValueError('Cycle summary must identify exactly this run')
    row=matched[0]
    for name in ('monitor.csv','manifest.json','run_metadata.txt','run_status.txt'):
        if row['conservation']['provenance'][name]!=sha(folder/name):raise ValueError('Stale cycle summary')
    result=dict(schema='PRIMS-F.fixed_trace_holdup.v1',run=folder.name,
        cycle_summary=dict(path=str(a.cycle_summary.resolve()),sha256=sha(a.cycle_summary)),
        source_column='alpha_riser_trace',transformation='1-alpha_riser_trace, without clipping',
        location='Matched junction trace at fixed s=Lp; recomputed from the accepted state, distinct from the accepted mass flux.',
        completed=row['completed'],primary_thirteen_channels_accommodated=row.get('accommodated'),
        supplemental_only=True,phase_grid_points=256,phase_origin='Same declared 20/80 pressure markers as the primary cycle summary.',
        script_sha256=sha(Path(__file__)),dependencies_sha256={name:sha(TOOLS/name) for name in
            ('cycle_waveform.py','validation_evidence.py')})
    if row['completed'] and row.get('cycles_available'):
        data,warnings=load_monitor(folder/'monitor.csv')
        if data is None:raise ValueError('Monitor has no usable observations')
        result.update(cycle_edges_s=row['cycle_edges_s'],monitor_warnings=warnings,
                      **calculate(data,row['cycle_edges_s']))
    else:result.update(cycle_edges_s=[],cycles=[],reason='No ten complete intervals in a completed run')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(output=str(a.output),cycles=len(result['cycles']))))

if __name__=='__main__':main()
