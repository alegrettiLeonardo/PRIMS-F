#!/usr/bin/env python3
"""Six phase/segment balances on declared cycles, using every accepted step."""
import argparse,csv,hashlib,json,math
from pathlib import Path

PHASES=('l','g')
REGIONS=(('', 'in','out'),('_pipe','in','base_left'),('_riser','base_right','out'))


def interpolate(left,right,time):
    if right['time_s']<=left['time_s']:raise ValueError('Accepted ledger time must advance')
    weight=(time-left['time_s'])/(right['time_s']-left['time_s'])
    if not -1e-12<=weight<=1+1e-12:raise ValueError('Cycle boundary lies outside the accepted-step bracket')
    return {key:(time if key=='time_s' else value+weight*(right[key]-value)) for key,value in left.items()}


def balance(start,row,phase,region,incoming,outgoing):
    return row[f'M{phase}{region}_kg']-start[f'M{phase}{region}_kg']- (
        row[f'M{phase}_{incoming}_integral_kg']-start[f'M{phase}_{incoming}_integral_kg'])+ (
        row[f'M{phase}_{outgoing}_integral_kg']-start[f'M{phase}_{outgoing}_integral_kg'])


def cycle_balances(rows,edges):
    if len(edges)<2 or not all(math.isfinite(t) for t in edges) or any(b<=a for a,b in zip(edges,edges[1:])):
        raise ValueError('Increasing finite cycle edges required')
    previous=None;boundary_index=0;start=None;maximum={};inventory_min={};inventory_max={};cycles=[];interior=0

    def begin(row):
        nonlocal start,maximum,inventory_min,inventory_max,interior
        start=row;interior=0
        maximum={p+r:0. for p in PHASES for r,_,_ in REGIONS}
        inventory_min={p+r:row[f'M{p}{r}_kg'] for p in PHASES for r,_,_ in REGIONS};inventory_max=inventory_min.copy()

    def observe(row):
        for p in PHASES:
            for r,inlet,outlet in REGIONS:
                key=p+r;value=row[f'M{p}{r}_kg']
                maximum[key]=max(maximum[key],abs(balance(start,row,p,r,inlet,outlet)))
                inventory_min[key]=min(inventory_min[key],value);inventory_max[key]=max(inventory_max[key],value)

    def end(row):
        observe(row);results={}
        for p in PHASES:
            for r,inlet,outlet in REGIONS:
                key=p+r;mass0=start[f'M{p}{r}_kg']
                processed=sum(row[f'M{p}_{face}_abs_integral_kg']-start[f'M{p}_{face}_abs_integral_kg'] for face in (inlet,outlet))
                error=balance(start,row,p,r,inlet,outlet)
                residual=row[f'mass_residual_sum_{p}{r}_kg']-start[f'mass_residual_sum_{p}{r}_kg']
                junction=0.
                if not r:
                    junction=(row[f'M{p}_base_left_integral_kg']-start[f'M{p}_base_left_integral_kg'])-(
                        row[f'M{p}_base_right_integral_kg']-start[f'M{p}_base_right_integral_kg'])
                results[p+(r or '_system')]=dict(initial_inventory_kg=mass0,minimum_inventory_kg=inventory_min[key],
                    maximum_inventory_kg=inventory_max[key],absolute_processed_mass_kg=processed,
                    final_balance_kg=error,max_abs_balance_kg=maximum[key],newton_residual_sum_kg=residual,
                    junction_defect_integral_kg=junction,raw_balance_vs_residual_difference_kg=error+junction-residual,
                    relative_inventory=maximum[key]/mass0 if mass0>0 else None,
                    relative_processed=maximum[key]/processed if processed>0 else None,
                    passes_inventory_1e6=maximum[key]<=1e-6*mass0 if mass0>0 else None,
                    passes_processed_1e6=maximum[key]<=1e-6*processed if processed>0 else None)
        cycles.append(dict(start_s=start['time_s'],end_s=row['time_s'],period_s=row['time_s']-start['time_s'],
                           accepted_rows_strictly_inside=interior,balances=results))

    for row in rows:
        if previous is None:
            if edges[0]<row['time_s']:raise ValueError('Cycle starts before ledger coverage')
            previous=row
            if edges[0]==row['time_s']:
                begin(row);boundary_index=1
            continue
        if row['time_s']<=previous['time_s']:raise ValueError('Nonconsecutive or repeated ledger time')
        while boundary_index<len(edges) and edges[boundary_index]<=row['time_s']:
            point=interpolate(previous,row,edges[boundary_index])
            if start is None:begin(point)
            else:
                end(point)
                if boundary_index<len(edges)-1:begin(point)
                else:start=None
            boundary_index+=1
        if start is not None and row['time_s']>start['time_s']:
            observe(row);interior+=1
        previous=row
        if boundary_index==len(edges):break
    if len(cycles)!=len(edges)-1:raise ValueError('Ledger does not cover every declared cycle')
    return cycles


def ledger_rows(path):
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            yield {key:float(value) for key,value in row.items()}


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-dir',type=Path,required=True)
    p.add_argument('--cycle-summary',type=Path,required=True);p.add_argument('--mass-summary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():p.error('Preserve existing evidence; choose another output')
    folder=a.run_dir.resolve();payload=json.loads(a.cycle_summary.read_text())
    matches=[r for r in payload['runs'] if r['run']==folder.name]
    if len(matches)!=1:raise ValueError('Cycle summary must identify this unique run')
    run=matches[0]
    for name in ('manifest.json','monitor.csv','run_status.txt','run_metadata.txt'):
        if run['conservation']['provenance'][name]!=sha(folder/name):raise ValueError('Cycle summary belongs to another run/version')
    mass=json.loads(a.mass_summary.read_text())
    if mass.get('schema')!='PRIMS-F.accepted_conservation.v2' or Path(mass['run_dir']).resolve()!=folder:
        raise ValueError('Accepted-conservation audit must identify this run')
    for name in ('conservation.csv','steps.csv','manifest.json','run_status.txt','run_metadata.txt'):
        if mass['artifacts'][name]['sha256']!=sha(folder/name):raise ValueError('Accepted-conservation audit is stale')
    result=dict(schema='PRIMS-F.accepted_cycle_balances.v2',run=folder.name,completed=run['completed'],
        accommodated=run.get('accommodated'),cycles_available=run.get('cycles_available'),pressure_column=run.get('pressure_column'),
        cycle_summary=dict(path=str(a.cycle_summary.resolve()),sha256=sha(a.cycle_summary)),
        mass_summary=dict(path=str(a.mass_summary.resolve()),sha256=sha(a.mass_summary),
                          arithmetic_contract_passes=mass['ledger']['arithmetic_contract_passes']),
        ledger=dict(path=str(folder/'conservation.csv'),sha256=sha(folder/'conservation.csv')),
        script_sha256=sha(Path(__file__)),
        method='Use all accepted ledger rows; linearly interpolate inventories and cumulative integrals at declared cycle edges. Partial-step boundaries are an explicit observation convention.',
        balance_policy='Errors are never corrected by residuals; each phase/segment uses strict 1e-6 initial-cycle inventory and absolute processed-mass criteria.')
    if run['completed'] and run.get('cycles_available'):
        edges=run['cycle_edges_s']
        if len(edges)!=11:raise ValueError('Expected ten complete intervals')
        result['cycles']=cycle_balances(ledger_rows(folder/'conservation.csv'),edges)
    else:result['cycles']=[]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(dict(output=str(a.output),cycles=len(result['cycles']),accommodated=result['accommodated'])))


if __name__=='__main__':main()
