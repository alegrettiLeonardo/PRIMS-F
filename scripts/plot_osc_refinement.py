#!/usr/bin/env python3
"""Plot the audited independent h/dt triples without a convergence claim."""
import argparse
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--spatial-coarse',type=Path,required=True)
    p.add_argument('--spatial-fine',type=Path,required=True)
    p.add_argument('--temporal',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    files=[a.spatial_coarse,a.spatial_fine,a.temporal]
    coarse,fine,temporal=[json.loads(path.read_text()) for path in files]
    if any(not item['eligible'] for item in (coarse,fine,temporal)):
        raise ValueError('The audited triples must be eligible')
    if coarse['scales'][1:]!=fine['scales'][:2]:
        raise ValueError('Spatial triples do not share the intermediate levels')
    if a.output.exists():
        raise ValueError('Preserve previous evidence')
    a.output.mkdir(parents=True)
    fig,axes=plt.subplots(3,2,figsize=(11,9),layout='constrained')
    observables=[('period_s',1,'Período (s)'),
        ('pressure.mean_semi_amplitude_Pa',1e-5,'Semi-amplitude na junção (bar)'),
        ('pressure.mean',1e-5,'Pressão média na junção (bar)')]
    table=[]
    for row,(key,factor,label) in enumerate(observables):
        hc=coarse['orders'][key]['h']
        hf=fine['orders'][key]['h']
        yc=coarse['orders'][key]['observable']
        yf=fine['orders'][key]['observable']
        if yc[1:]!=yf[:2]:
            raise ValueError('Shared spatial observations differ')
        xh=hc[:1]+hf
        yh=yc[:1]+yf
        xt=temporal['orders'][key]['h']
        yt=temporal['orders'][key]['observable']
        for ax,x,y,color in ((axes[row,0],xh,yh,'#235789'),(axes[row,1],xt,yt,'#b45128')):
            ax.plot(x,[v*factor for v in y],marker='o',color=color,lw=1.3)
            ax.set_xscale('log',base=2)
            ax.invert_xaxis()
            ax.set_ylabel(label)
            ax.grid(alpha=.25)
        table.append(dict(observable=key,h=xh,spatial_values=yh,dt_s=xt,temporal_values=yt))
    axes[0,0].set_title('Malha · dt máximo 0,125 s')
    axes[0,1].set_title('Passo · 51+51 volumes')
    axes[-1,0].set_xlabel('Δs do oleoduto (m) · mais fino à direita')
    axes[-1,1].set_xlabel('dt (s) · mais fino à direita')
    fig.suptitle('OSC · dez intervalos acomodados; ordem observada não demonstra resolução do outro eixo',fontsize=10)
    target=a.output/'osc_refinamento.png'
    fig.savefig(target,dpi=180)
    plt.close(fig)
    (a.output/'index.json').write_text(json.dumps(dict(
        inputs={str(path):sha(path) for path in files},script_sha256=sha(Path(__file__)),
        figure_sha256=sha(target),data=table,
        fixed_other_axis_proven_resolved=all(item['fixed_other_axis_proven_resolved'] for item in (coarse,fine,temporal)),
        interpretation='Observed discrete trends at fixed other axes. No GCI, asymptotic or experimental validation claimed.'
    ),indent=2)+'\n')


if __name__=='__main__':
    main()
