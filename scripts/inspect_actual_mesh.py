#!/usr/bin/env python3
"""Recover FV faces from the saved initial centers without reading full fields."""
from pathlib import Path
import csv,json,hashlib,argparse


def inspect(folder,parameters):
    path=folder/'fields.csv';rows=[]
    with path.open(newline='') as stream:
        reader=csv.DictReader(stream)
        for row in reader:
            t=float(row['time_s'])
            if t!=0.:break
            rows.append(row)
    np=int(parameters['nodes']);nr=int(parameters['n_riser'])
    if len(rows)!=np+nr:raise ValueError('Initial field does not cover both complete physical segments')
    if parameters['mesh_ratio']!=1.:raise ValueError('This audit assumes the declared uniform mesh within each segment')
    centers=[float(row['s_m']) for row in rows];faces=[0.]
    for center in centers:faces.append(2*center-faces[-1])
    widths=[b-a for a,b in zip(faces,faces[1:])]
    lp=float(parameters['lp']);lr=faces[-1]-lp
    if abs(faces[np]-lp)>1e-9:raise ValueError('Reconstructed interface disagrees with prescribed pipe length')
    for group,length,count in ((widths[:np],lp,np),(widths[np:],lr,nr)):
        if min(group)<=0 or max(abs(value-length/count) for value in group)>1e-9:
            raise ValueError('Observed FV widths are inconsistent with length/N')
    raw='\n'.join(','.join(row.values()) for row in rows).encode()
    return dict(folder=str(folder),initial_field_rows=len(rows),initial_frame_canonical_sha256=hashlib.sha256(raw).hexdigest(),
        cells_pipe=np,cells_riser=nr,length_pipe_m=lp,length_riser_m=lr,
        h_pipe_m=lp/np,h_riser_m=lr/nr,max_width_m=max(widths),min_width_m=min(widths),
        pipe_first_center_m=centers[0],riser_first_center_m=centers[np],
        junction_position_m=faces[np],right_boundary_m=faces[-1],
        interpretation='Frozen malha_progressiva normalizes N cell widths; h=L/N, not L/(N-1). Faces independently reconstructed as f[i+1]=2*c[i]-f[i].')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();plan=json.loads(a.plan.read_text());result=[]
    for job in plan['jobs']:
        folder=Path(job['output'])
        if not (folder/'fields.csv').exists():continue
        result.append(inspect(folder,job['parameters']))
    with a.output.open('x') as stream:json.dump(dict(meshes=result,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),stream,indent=2);stream.write('\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
