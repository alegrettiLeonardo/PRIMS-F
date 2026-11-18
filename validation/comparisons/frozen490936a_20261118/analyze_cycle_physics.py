#!/usr/bin/env python3
"""Cycle-conditioned inventory/phase-flow observations, without regime validation."""
import argparse,csv,hashlib,json,statistics,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--run-dir',type=Path,required=True)
p.add_argument('--envelopes',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();a.output.mkdir(exist_ok=True,parents=True)
raw=(a.run_dir/'monitor.csv').read_bytes()
data=np.genfromtxt(raw.decode().splitlines(),delimiter=',',names=True)
case=next(r for r in json.loads(a.envelopes.read_text())['cases'] if Path(r['source_directory']).name==a.run_dir.name)
cycles=case['last_ten_cycles']
assert cycles['available'] and cycles['n_complete_cycles']==10
peak_times=np.array([r['peak_time_s'] for r in cycles['peak_records']])
assert len(peak_times)==11 and np.all(np.diff(peak_times)>0)
phase=np.arange(256)/256
columns=['p_junction_Pa','holdup_liquid_base','holdup_liquid_out','Ml_riser_kg','Mg_pipe_kg',
         'mdot_l_base_kg_s','mdot_l_out_kg_s','mdot_g_base_kg_s','mdot_g_out_kg_s']
waves={key:np.array([np.interp(t0+phase*(t1-t0),data['time_s'],data[key])
                    for t0,t1 in zip(peak_times[:-1],peak_times[1:])]) for key in columns}
means={key:values.mean(axis=0) for key,values in waves.items()}
period=np.diff(peak_times).mean()
sys.path.insert(0,str(a.run_dir.resolve().parents[2]/'scripts'))
from validation_evidence import load_monitor,time_statistics
monitor,warnings=load_monitor(a.run_dir/'monitor.csv')
channels={'p_junction_Pa':'Pa','holdup_liquid_base':'1','holdup_liquid_out':'1',
 'mdot_l_base_kg_s':'kg/s','mdot_g_base_kg_s':'kg/s','mdot_l_out_kg_s':'kg/s','mdot_g_out_kg_s':'kg/s',
 'Ml_kg':'kg','Mg_kg':'kg','Ml_riser_kg':'kg','Mg_pipe_kg':'kg','Ml_pipe_kg':'kg','Mg_riser_kg':'kg'}
accommodation={}
for key,unit in channels.items():
  blocks=[time_statistics(monitor,key,t0,t1) for t0,t1 in zip(peak_times[:-1],peak_times[1:])]
  assert all(x is not None and x['coverage_fraction']>=1-1e-8 for x in blocks),key
  means_each=[x['mean'] for x in blocks];ranges_each=[x['peak_to_peak'] for x in blocks]
  first_mean,last_mean=statistics.mean(means_each[:5]),statistics.mean(means_each[-5:])
  first_range,last_range=statistics.mean(ranges_each[:5]),statistics.mean(ranges_each[-5:])
  mean_scale=max(abs(first_mean),abs(last_mean),1e-30)
  range_scale=max(abs(first_range),abs(last_range),1e-30)
  accommodation[key]={'unit':unit,'first_five_mean':first_mean,'last_five_mean':last_mean,
    'mean_absolute_change':last_mean-first_mean,'mean_normalization':mean_scale,
    'mean_relative_change':(last_mean-first_mean)/mean_scale,
    'first_five_mean_excursion':first_range,'last_five_mean_excursion':last_range,
    'excursion_absolute_change':last_range-first_range,'excursion_normalization':range_scale,
    'excursion_relative_change':(last_range-first_range)/range_scale,
    'cycle_means':means_each,'cycle_excursions':ranges_each}
def derivative(y):return (np.roll(y,-1)-np.roll(y,1))/(2*period/len(y))
dp=derivative(means['p_junction_Pa']);dmg=derivative(means['Mg_pipe_kg']);dml=derivative(means['Ml_riser_kg'])
def features(y):
  return {'mean':float(y.mean()),'minimum':float(y.min()),'maximum':float(y.max()),
          'peak_to_peak':float(np.ptp(y)),
          'maximum_phase_after_pressure_peak':float(phase[np.argmax(y)]),
          'minimum_phase_after_pressure_peak':float(phase[np.argmin(y)])}
result={'run':a.run_dir.name,'evidence_class':'numerical_cycle_conditioned_observations_not_regime_validation',
 'analysis_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'monitor_sha256':hashlib.sha256(raw).hexdigest(),'envelopes_sha256':hashlib.sha256(a.envelopes.read_bytes()).hexdigest(),
 'cycles':10,'window_start_s':float(peak_times[0]),'window_end_s':float(peak_times[-1]),'period_mean_s':float(period),
 'phase_reference':'pressure maximum = 0; phases increase to next pressure maximum',
 'last_ten_channel_accommodation':accommodation,
 'accommodation_normalization':'Each mean/excursion change divided by max(abs(first-five),abs(last-five),1e-30 in that channel unit); exact zero/zero gives zero. Absolute differences and scales are retained.',
 'accommodation_statistics':'Original monitor, piecewise-linear time weighting on each of ten pressure-peak intervals; same definitions as frozen baseline13-channel summary.',
 'one_percent_channel_checks':{'all_means':all(abs(v['mean_relative_change'])<=.01 for v in accommodation.values()),
                              'all_excursions':all(abs(v['excursion_relative_change'])<=.01 for v in accommodation.values())},
 'data_quality_warnings':warnings,
 'method':'256 equally spaced phase samples per complete pressure cycle; interpolate each named signal then average ten cycles',
 'signals':{key:features(y) for key,y in means.items()},
 'pressure_gas_pipe_inventory_correlation':float(np.corrcoef(means['p_junction_Pa'],means['Mg_pipe_kg'])[0,1]),
 'pressure_liquid_riser_inventory_correlation':float(np.corrcoef(means['p_junction_Pa'],means['Ml_riser_kg'])[0,1]),
 'fractions_of_average_cycle':{
   'rising_pressure_and_gas_pipeline_accumulation':float(np.mean((dp>0)&(dmg>0))),
   'rising_pressure_and_liquid_riser_accumulation':float(np.mean((dp>0)&(dml>0))),
   'liquid_outflow_exceeds_base_inflow_and_liquid_riser_depletes':float(np.mean((means['mdot_l_out_kg_s']>means['mdot_l_base_kg_s'])&(dml<0))),
   'falling_pressure_and_gas_pipeline_depletion':float(np.mean((dp<0)&(dmg<0)))},
 'gas_base_remains_positive_in_all_resampled_cycles':bool(np.all(waves['mdot_g_base_kg_s']>0)),
 'limits':['Phase averages describe observed numerical signals; they do not validate a flow regime.',
           'Pressure/inventory synchrony supports accumulation/compression and depletion proxies; no local front penetration was measured.',
           'No complete gas blockage is required or imposed for OSC/SS3; positive gas at the base is retained.',
           'Integrated conservation uses original monitor integrals, not derivatives of phase averages.']}
(a.output/'cycle_physics.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
with (a.output/'mean_cycle.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['phase',*columns]);w.writerows(zip(phase,*[means[k] for k in columns]))
fig,ax=plt.subplots(3,2,figsize=(11,9),sharex=True,constrained_layout=True)
groups=[(['p_junction_Pa'],'Pressure [Pa]'),(['holdup_liquid_base','holdup_liquid_out'],'Liquid fraction [-]'),
        (['Ml_riser_kg'],'Riser liquid inventory [kg]'),(['Mg_pipe_kg'],'Pipeline gas inventory [kg]'),
        (['mdot_l_base_kg_s','mdot_l_out_kg_s'],'Liquid mass flow [kg/s]'),
        (['mdot_g_base_kg_s','mdot_g_out_kg_s'],'Gas mass flow [kg/s]')]
for axis,(keys,ylabel) in zip(ax.flat,groups):
  for key in keys:
    y=means[key];axis.plot(np.r_[phase,1],np.r_[y,y[0]],label=key)
    axis.fill_between(phase,waves[key].min(axis=0),waves[key].max(axis=0),alpha=.15)
  axis.set_ylabel(ylabel);axis.grid(alpha=.25)
  if len(keys)>1:axis.legend(fontsize=7)
for axis in ax[-1]:axis.set_xlabel('Phase after pressure maximum [cycles]')
fig.suptitle(a.run_dir.name+'\nTen final complete cycles; shaded range across cycles',fontsize=10)
fig.savefig(a.output/'mean_cycle.png',dpi=180);fig.savefig(a.output/'mean_cycle.pdf');plt.close(fig)
print(json.dumps({'run':a.run_dir.name,'period_mean_s':period,'output':str(a.output),'phase_fractions':result['fractions_of_average_cycle']}))
