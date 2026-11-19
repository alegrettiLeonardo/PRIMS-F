"""Piecewise-linear cycle observations and diagnostic phase comparisons."""
import math
import statistics
from validation_evidence import value_at


def difference(a,b,atol):
    if a is None or b is None:return {'available':False}
    d=abs(b-a);scale=max(abs(a),abs(b))
    return {'available':True,'coarse':a,'fine':b,'absolute_difference':d,
            'relative_difference':d/scale if scale else None,'absolute_tolerance':atol,
            'passes_1percent_plus_floor':d<=atol+.01*scale}


def waveform(data,key,edges,points=128):
    rows=[]
    for left,right in zip(edges,edges[1:]):
        row=[value_at(data['time_s'],data[key],left+(right-left)*i/points) for i in range(points)]
        if any(v is None for v in row):return None
        rows.append(row)
    if not rows:return None
    return {'phase':[i/points for i in range(points)],'n_cycles':len(rows),
            'mean':[statistics.mean(x) for x in zip(*rows)],
            'std_across_cycles':[statistics.pstdev(x) for x in zip(*rows)],
            'alignment':'each pressure marker interval mapped linearly to [0,1); period reported separately'}


def wave_difference(a,b,period):
    if not a or not b:return {'available':False}
    a,b=a['mean'],b['mean'];n=len(a)
    if len(b)!=n:raise ValueError('Different phase grids')
    ma,mb=statistics.mean(a),statistics.mean(b)
    aa=[x-ma for x in a];bb=[x-mb for x in b]
    scale=math.sqrt(math.fsum(x*x for x in aa)*math.fsum(x*x for x in bb))
    out={'available':True,'phase_fixed_rmse':math.sqrt(math.fsum((x-y)**2 for x,y in zip(a,b))/n),
         'mean_difference':mb-ma,'phase_fixed_max_absolute_error':max(abs(x-y) for x,y in zip(a,b))}
    if scale<=0:return out|{'best_lag_fraction':None,'reason':'constant_waveform_lag_undefined'}
    shifts=range(-n//2,n//2)
    lag=max(shifts,key=lambda k:math.fsum(aa[i]*bb[(i+k)%n] for i in range(n)))
    out.update(best_lag_fraction=lag/n,best_lag_s=lag/n*period if period else None,
               best_lag_correlation=math.fsum(aa[i]*bb[(i+lag)%n] for i in range(n))/scale,
               lag_adjusted_rmse=math.sqrt(math.fsum((a[i]-b[(i+lag)%n])**2 for i in range(n))/n),
               lag_convention='fine(phi+lag) compared with coarse(phi); diagnostic only, fixed-phase error retained')
    return out
