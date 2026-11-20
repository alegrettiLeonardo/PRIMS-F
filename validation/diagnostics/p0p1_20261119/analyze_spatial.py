#!/usr/bin/env python3
"""Phase-aligned field diagnostics and geometric isovoid crossings, no regime labels."""
import argparse,csv,json,sys,hashlib
from itertools import groupby
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'scripts'))
from analyze_burr15_campaign import front_candidates,unique_centres

def read(folder,left,right,frontfile):
 times=[];values=[];positions=None;counts={str(x):0 for x in [.4,.5,.6]}
 with(folder/'fields.csv').open()as f,frontfile.open('w')as output:
  writer=csv.writer(output);writer.writerow(['time_s','alpha_level','s_m','orientation'])
  for t,group in groupby(csv.DictReader(f),key=lambda x:float(x['time_s'])):
   if t<left or t>right:continue
   rows=unique_centres(group);s=np.array([float(x['s_m'])for x in rows])
   if positions is None:positions=s
   else:assert np.array_equal(s,positions)
   u=np.array([[float(x[k])for x in rows]for k in ['p_Pa','alpha','ul_m_s','ug_m_s']])
   times.append(t);values.append(np.r_[u,((1-u[1])*u[2])[None,:],(u[1]*u[3])[None,:]])
   for level,x,orientation in front_candidates(rows,57.4,[.4,.5,.6]):
    writer.writerow([t,level,x,orientation]);counts[str(level)]+=1
 return np.array(times),positions,np.array(values),counts

def main():
 p=argparse.ArgumentParser();p.add_argument('--summary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
 summaries=json.loads(a.summary.read_text())['runs'];records=[]
 # Common observation grid on the intersection of all cell-centre domains.
 positions={}
 for s in summaries:
  folder=Path(s['run']);lines=(folder/'checkpoint.dat').read_text().splitlines();head=lines[1].split();nd=int(head[1]);n=int(head[0]);coords=np.array([float(x.split()[1])*12.54625399319708 for x in lines[2:]])[nd:n-nd];positions[s['base_run']]=coords
 lo=max(x[0]for x in positions.values());hi=min(x[-1]for x in positions.values());pmax=min(x[x<57.4][-1]for x in positions.values());rmin=max(x[x>57.4][0]for x in positions.values())
 grid=np.r_[np.linspace(lo,pmax,256),np.linspace(rmin,hi,256)]
 profiles={};channels=['p_Pa','alpha','ul_m_s','ug_m_s','jl_m_s','jg_m_s']
 for s in summaries:
  if 'window_s'not in s:continue
  folder=Path(s['run']);name=Path(s['base_run']).name;left,right=s['window_s'];t,x,u,counts=read(folder,left-2,right+2,a.output/(name+'_fronts.csv'))
  edges=sorted(set([v['peak_time_s']for v in s['cycles']['peak_records']]+[right]));phase=np.arange(256)/256
  # Field interpolation is diagnostic; phase mass balances use accepted ledgers elsewhere.
  target=np.array([l+(r-l)*phase for l,r in zip(edges[:-1],edges[1:])]);assert t[0]<=target.min()and t[-1]>=target.max()
  mean=np.empty((6,256,len(x)))
  for c in range(6):
   mean[c]=np.array([np.interp(target.ravel(),t,u[:,c,j]).reshape(10,256).mean(0)for j in range(len(x))]).T
  mapped=np.empty((6,256,len(grid)))
  for c in range(6):
   for ph in range(256):
    # Interpolate within each physical segment; never smear the closure jump.
    mapped[c,ph,:256]=np.interp(grid[:256],x[x<57.4],mean[c,ph,x<57.4])
    mapped[c,ph,256:]=np.interp(grid[256:],x[x>57.4],mean[c,ph,x>57.4])
  profiles[name]=mapped;np.savez_compressed(a.output/(name+'_phase_fields.npz'),phase=phase,s_m=x,fields=mean,common_s_m=grid,common_fields=mapped)
  records.append(dict(name=name,run=str(folder),window_s=[left,right],frames=len(t),field_support_time_s=[float(t[0]),float(t[-1])],crossing_count_scope='raw field support includes up to2s interpolation padding outside cycle markers',field_dt_max_s=float(np.max(np.diff(t))),crossing_counts=counts,fields_sha256=hashlib.sha256((folder/'fields.csv').read_bytes()).hexdigest()))
 comparisons=[]
 for c in json.loads(a.summary.read_text())['comparisons']:
  one=profiles[c['coarse']];two=profiles[c['fine']];entry=dict(coarse=c['coarse'],fine=c['fine'],axis=c['axis'],channels={})
  for j,key in enumerate(channels):
   err=two[j]-one[j];entry['channels'][key]=dict(rms_absolute=float(np.sqrt(np.mean(err**2))),max_absolute=float(np.max(abs(err))),fine_rms=float(np.sqrt(np.mean(two[j]**2))))
  comparisons.append(entry)
 (a.output/'index.json').write_text(json.dumps(dict(runs=records,comparisons=comparisons,common_grid_domains_m=[[float(lo),float(pmax)],[float(rmin),float(hi)]],scope='Geometric alpha-level candidates and phase averages on common interior observation domains. No material-front speed, complete blocking criterion or experimental regime classification. 1s field sampling does not resolve acoustic fronts.'),indent=2)+'\n')
 print('PASS diagnostic fields',len(records),'comparisons',len(comparisons))
if __name__=='__main__':main()
