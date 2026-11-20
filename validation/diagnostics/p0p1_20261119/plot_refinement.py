#!/usr/bin/env python3
"""Scientific plots of the measured cycle waveforms and N401 phase fields."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--summary',type=Path,required=True);p.add_argument('--spatial',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
 runs=json.loads(a.summary.read_text())['runs'];phase=np.arange(256)/256
 selected=[x for x in runs if (x['dt_max_s']==.125 and x['n']in[201,401])or(x['n']==51 and x['dt_max_s']in[.0625,.03125])]
 fig,axes=plt.subplots(3,2,figsize=(11,10),sharex=True)
 channels=[('p_riser_base_Pa','Base pressure (bar absolute)',1e-5),('holdup_liquid_base','Liquid fraction at base',1),('Ml_kg','Liquid inventory (kg)',1),('Mg_kg','Gas inventory (kg)',1),('mdot_l_base_kg_s','Liquid mass flow at base (kg/s)',1),('mdot_g_base_kg_s','Gas mass flow at base (kg/s)',1)]
 for ax,(key,label,scale)in zip(axes.flat,channels):
  for run in selected:
   ax.plot(phase,np.array(run['phase_waveforms'][key]['mean'])*scale,label=f"N={run['n']}+{run['n']}, dt={run['dt_max_s']}")
  ax.set_ylabel(label);ax.set_xlabel('Cycle phase');ax.grid(alpha=.2)
 axes[0,0].legend(fontsize=7);fig.suptitle('Last ten settled cycles; same physical model, separate h/dt cohorts');fig.tight_layout();fig.savefig(a.output/'osc_phase_waveforms.png',dpi=150);plt.close(fig)
 x=np.load(a.spatial/'N401_dt0.125_phase_fields.npz');n=x['fields'].shape[-1]//2
 edges=np.r_[np.linspace(0,57.4,n+1),57.4+np.linspace(0,12.54625399319708,n+1)[1:]]
 fig,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
 for ax,(k,title,scale)in zip(axes.flat,[(0,'Pressure (bar absolute)',1e-5),(1,'Void fraction',1),(4,'Signed liquid superficial velocity (m/s)',1),(5,'Signed gas superficial velocity (m/s)',1)]):
  im=ax.pcolormesh(np.linspace(0,1,257),edges,x['fields'][k].T*scale,shading='flat');ax.axhline(57.4,color='white',lw=.8);fig.colorbar(im,ax=ax);ax.set_title(title);ax.set_xlabel('Cycle phase');ax.set_ylabel('Axial position (m)')
 fig.suptitle('N401+401, dt_max=0.125 s: phase-averaged fields, not experimental regime labels');fig.tight_layout();fig.savefig(a.output/'osc_N401_phase_fields.png',dpi=150);plt.close(fig)
if __name__=='__main__':main()
