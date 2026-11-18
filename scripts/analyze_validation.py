#!/usr/bin/env python3
"""Measure numerical evidence without equating oscillation to validation.

An experimental period is not a forcing or calibration target. Amplitude below
100 Pa is below this protocol's cycle-detection threshold, not proof of physical
stability. CSV times are actual accepted solver times; no uniformity is assumed.
"""
from __future__ import annotations
import argparse,csv,json,re,math
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks,periodogram,detrend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
# NumPy 1.x names this same integration rule trapz; NumPy 2.x uses trapezoid.
trapezoidal_integral=np.trapezoid if hasattr(np,'trapezoid') else np.trapz


def number(x):
    x=float(x)
    return x if math.isfinite(x) else None

def time_statistics(t,values):
    """Trapezoidal time weights avoid bias from adaptive accepted steps."""
    if len(t)<2 or t[-1]<=t[0]:return float(values[-1]),0.0
    mean=float(trapezoidal_integral(values,t)/(t[-1]-t[0]))
    variance=float(trapezoidal_integral((values-mean)**2,t)/(t[-1]-t[0]))
    return mean,float(np.sqrt(max(variance,0.0)))

def read_metadata(path):
    text=path.read_text(errors='replace') if path.exists() else ''
    pars={}
    for m in re.finditer(r'(\w+)\s*=\s*(\"[^\"]*\"|\'[^\']*\'|[^,\n]+)',text):
        key,val=m[1].lower(),m[2].strip()
        if val[0:1] in ['\"',"'"]: pars[key]=val.strip('\"\' ').strip();continue
        try: pars[key]=float(val.replace('D','E').replace('d','e'))
        except ValueError:pars[key]=val
    return pars

def load_monitor(path):
    with path.open() as f:
        rows=list(csv.reader(f))
    if len(rows)<2:return None
    names=rows[0];good=[]
    for row in rows[1:]:
        if len(row)!=len(names):continue
        try:good.append([float(v) for v in row])
        except ValueError:continue
    if not good:return None
    raw=np.array(good);raw=raw[np.isfinite(raw[:,0])]
    _,unique=np.unique(raw[:,0],return_index=True);raw=raw[np.sort(unique)]
    return {k:raw[:,i] for i,k in enumerate(names)}

def measure(folder,start,amplitude_floor):
    data=load_monitor(folder/'monitor.csv')
    if data is None:return None,None
    pars=read_metadata(folder/'run_metadata.txt')
    status={}
    if (folder/'run_status.txt').exists():
        for line in (folder/'run_status.txt').read_text().splitlines():
            if '=' in line:
                k,v=line.split('=',1);status[k]=v.strip()
    manifest=json.loads((folder/'manifest.json').read_text()) if (folder/'manifest.json').exists() else {}
    t=data['time_s'];end=float(t[-1]);request=float(status.get('requested_time_s',pars.get('time_dim',end)))
    complete=status.get('status')=='0' and abs(end-request)<=max(1e-6,request*1e-10)
    actual_start=start if end>start+10 else end/3
    selection=(t>=actual_start);tc=t[selection]
    pkey='p_junction_Pa' if 'p_junction_Pa' in data and np.all(np.isfinite(data['p_junction_Pa'])) else 'p_riser_base_Pa'
    press=data[pkey][selection]
    out={'run':folder.name,'case':str(pars.get('regime','unknown')).strip(),'completed':complete,
         'run_status':status.get('status',manifest.get('status','INCOMPLETE')),'reached_time_s':end,'requested_time_s':request,
         'window_start_s':actual_start,'window_end_s':end,'window_policy':'900s discard' if actual_start==start else 'pilot: final two thirds',
         'pressure_location':pkey,'parameters':pars,'manifest':manifest,'mass':{},'cycles':{},'phase_indicators':{}}
    for key,phase in [('l','liquid'),('g','gas')]:
        mass=data['M'+key+'_kg'];inflow=data['M'+key+'_in_integral_kg'];error=data['balance_'+key+'_kg']
        denom=max(abs(float(mass[0])),abs(float(inflow[-1])),1e-12)
        out['mass'][phase]={'max_absolute_error_kg':number(np.max(abs(error))),
            'final_error_kg':number(error[-1]),'normalization_kg':denom,
            'max_relative_error':number(np.max(abs(error))/denom),
            'passes_1e_minus_6_relative':bool(np.max(abs(error))<=1e-8+1e-6*denom)}
    if len(tc)>3:
        dt=float(np.median(np.diff(tc)))
        ptp=float(np.ptp(press));pmean,rms=time_statistics(tc,press);prom=max(amplitude_floor,0.10*ptp)
        # Interpolation is only for uniform-sample spectral processing.
        tu=np.arange(tc[0],tc[-1]+dt*.1,dt);pu=np.interp(tu,tc,press)
        peaks,props=find_peaks(pu,prominence=prom,distance=2)
        times=tu[peaks];periods=np.diff(times)
        cycles={'pressure_mean_Pa':pmean,'statistics':'trapezoidal time weighting','pressure_peak_to_peak_Pa':ptp,
                'pressure_rms_Pa':rms,'detection_floor_Pa':amplitude_floor,'prominence_threshold_Pa':prom,
                'sample_interval_s':dt,'peaks_found':len(peaks),'complete_intervals_found':len(periods),
                'peak_times_s':times.tolist(),'periods_s':periods.tolist(),'median_period_s':None,'period_cv':None,
                'spectral_period_s':None,'ten_intervals_observed':len(periods)>=10,
                'classification':'no_resolved_oscillation' if ptp<amplitude_floor else 'transient_or_oscillation'}
        if len(periods)>=2:
            cycles['median_period_s']=float(np.median(periods));cycles['period_cv']=float(np.std(periods)/np.mean(periods))
            cycles['samples_per_period']=cycles['median_period_s']/dt
            cycles['classification']='sustained_cycle_candidate' if len(periods)>=10 and cycles['period_cv']<.2 else 'nonstationary_or_modulated'
        if rms>=1 and ptp>=amplitude_floor and len(pu)>=16:
            frequencies,power=periodogram(detrend(pu),fs=1/dt,window='hann')
            keep=(frequencies>0)&(frequencies>=2/(tu[-1]-tu[0]))
            if keep.any(): cycles['spectral_period_s']=float(1/frequencies[keep][np.argmax(power[keep])])
        blocks=[]
        for ids in np.array_split(np.arange(len(tc)),3):
            bm,bs=time_statistics(tc[ids],press[ids])
            blocks.append({'from_s':float(tc[ids[0]]),'to_s':float(tc[ids[-1]]),'mean_pressure_Pa':bm,
                           'peak_to_peak_Pa':float(np.ptp(press[ids])),'rms_Pa':bs})
        cycles['stationarity_blocks']=blocks
        out['cycles']=cycles
        alpha=data['alpha_riser_base'][selection]
        indicators={'min_alpha_riser_cell':number(np.min(alpha)),'max_alpha_riser_cell':number(np.max(alpha)),
                    'liquid_inventory_range_kg':number(np.ptp(data['Ml_kg'][selection])),
                    'gas_inventory_range_kg':number(np.ptp(data['Mg_kg'][selection])),
                    'min_liquid_outlet_kg_s':number(np.min(data['mdot_l_out_kg_s'][selection])),
                    'max_liquid_outlet_kg_s':number(np.max(data['mdot_l_out_kg_s'][selection])),
                    'min_gas_outlet_kg_s':number(np.min(data['mdot_g_out_kg_s'][selection])),
                    'max_gas_outlet_kg_s':number(np.max(data['mdot_g_out_kg_s'][selection]))}
        if 'blocking_active' in data:
            indicators['experimental_barrier_active_fraction']=number(np.mean(data['blocking_active'][selection]))
            indicators['min_gas_base_kg_s']=number(np.min(data['mdot_g_base_kg_s'][selection]))
            indicators['max_gas_base_kg_s']=number(np.max(data['mdot_g_base_kg_s'][selection]))
        # Indicators record evidence, not a predetermined ordered cycle label.
        out['phase_indicators']=indicators
    out['interpretation']='Numerical evidence only; period/oscillation alone does not establish physical slugging validation.'
    if pars.get('junction_model')=='BLOCKING':out['interpretation']+=' EXPERIMENTAL gas-barrier projection: not the reference moving-front model.'
    return out,data


def plots(metric,data,destination):
    t=data['time_s'];pkey=metric['pressure_location'];start=metric['window_start_s']
    fig,axes=plt.subplots(3,2,figsize=(11.7,8.3),sharex=True,layout='constrained')
    pairs=[(axes[0,0],pkey,'Pressão na junção / célula de base','bar',1e-5),
           (axes[0,1],'alpha_riser_base','Fração de gás na primeira célula do riser','–',1),
           (axes[1,0],'mdot_l_out_kg_s','Vazão mássica de líquido na saída','kg/s',1),
           (axes[1,1],'mdot_g_out_kg_s','Vazão mássica de gás na saída','kg/s',1),
           (axes[2,0],'Ml_kg','Inventário de líquido','kg',1),
           (axes[2,1],'Mg_kg','Inventário de gás','kg',1)]
    for ax,key,title,unit,factor in pairs:
        ax.plot(t,data[key]*factor,color='#175c8c',lw=.8)
        ax.axvline(start,color='#be6b25',ls='--',lw=.7)
        ax.set_title(title,fontsize=10);ax.set_ylabel(unit);ax.grid(alpha=.22)
    for ax in axes[-1]:ax.set_xlabel('Tempo físico (s)')
    pars=metric['parameters']
    title=f"{metric['case']} | fechamento {pars.get('closure_riser',0):g} | {pars.get('gas_eos','')} | {pars.get('junction_model','MATCHED')}"
    details=f"{pars.get('nodes',0):g} + {pars.get('n_riser',0):g} células, Δt máximo {pars.get('dt_max_s',0):g} s — "+('Execução completa' if metric['completed'] else 'Execução incompleta')
    fig.suptitle(title+'\n'+details,fontsize=11)
    fig.savefig(destination/(metric['run']+'.png'),dpi=150)
    fig.savefig(destination/(metric['run']+'.pdf'))
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=ROOT/'validation/results')
    p.add_argument('--output',type=Path,default=ROOT/'validation/analysis')
    p.add_argument('--pattern',default='final_*')
    p.add_argument('--window-start',type=float,default=900)
    p.add_argument('--amplitude-floor',type=float,default=100)
    p.add_argument('--plots',action='store_true')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    metrics=[]
    for folder in sorted(args.results.glob(args.pattern)):
        if not (folder/'monitor.csv').exists():continue
        m,d=measure(folder,args.window_start,args.amplitude_floor)
        if m is None:continue
        metrics.append(m)
        if args.plots:plots(m,d,args.output)
    (args.output/'metrics.json').write_text(json.dumps({'schema':'PRIMS-F.validation.v1','runs':metrics},indent=2,allow_nan=False))
    cols=['run','completed','reached_time_s','pressure_peak_to_peak_Pa','median_period_s','period_cv','complete_intervals_found',
          'mass_error_liquid_kg','mass_error_gas_kg','classification']
    with (args.output/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        for m in metrics:
            row={k:m.get(k,m['cycles'].get(k,'')) for k in cols}
            row['mass_error_liquid_kg']=m['mass']['liquid']['max_absolute_error_kg']
            row['mass_error_gas_kg']=m['mass']['gas']['max_absolute_error_kg'];w.writerow(row)
    print(json.dumps({'runs':len(metrics),'completed':sum(m['completed'] for m in metrics),'output':str(args.output)}))
if __name__=='__main__':main()
