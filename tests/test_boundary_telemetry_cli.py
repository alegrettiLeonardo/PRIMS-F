#!/usr/bin/env python3
"""Independent accepted-boundary/atomic-checkpoint contract checks (short CLI)."""
import argparse,csv,json,math,os,re,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
from test_pipe_regime_cli import make_deck

def rows(p):
    with p.open() as f:return list(csv.DictReader(f))

def check_case(binary,out,name,extra):
    folder=out/name;deck=out/(name+'.nml')
    params=dict(nodes=4,n_riser=4,time_dim=.02,dt_start_s=.001,dt_max_s=.004,
        tol_res=1e-10,tol_inc=1e-10,newton_max=25,output_interval_s=1,field_interval_s=1,
        out_dir="'"+str(folder)+"'",**extra)
    base=(ROOT/'input/benchmarks/D_Stable.nml').read_text()
    for key in ('mul','mug'):
        basis='liquid_flow_basis' if key=='mul' else 'gas_flow_basis'
        if 'LOCAL_VOLUME' in str(extra.get(basis,'')):
            base=re.sub(r'(?i)\b'+key+r'\s*=\s*[^,\n/]+,?', '',base)
    deck.write_text(make_deck(base,params))
    p=subprocess.run([str(binary),str(deck)],cwd=ROOT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),capture_output=True,text=True,timeout=90)
    (out/(name+'.log')).write_text(p.stdout+p.stderr)
    bc=rows(folder/'boundary_conditions.csv');led=rows(folder/'conservation.csv')
    assert len(bc)==len(led)>1
    for b,l in zip(bc,led):
        assert int(b['step'])==int(l['step'])
        assert float(b['time_s'])==float(l['time_s']) and float(b['dt_s'])==float(l['dt_s'])
        assert b['status_in']==b['status_out']=='0'
        assert float(b['recompute_max_nd'])<1e-12
        for phase in ('l','g'):
            eff=float(b[f'mdot_{phase}_eff_kg_s']);target=float(b[f'mdot_{phase}_target_kg_s'])
            assert eff==float(l[f'mdot_{phase}_in_kg_s'])
            assert abs(float(b[f'rm_{phase}_kg_s'])-(eff-target))<1e-14
            rel=float(b[f'relative_{phase}'])
            if target==0:assert math.isnan(rel) and eff==0
            else:assert abs(rel)<1e-12
        assert abs(float(b['p_out_flux_Pa'])-float(b['Pt_Pa']))<2e-5
        assert abs(float(b['rp_Pa'])-(float(b['p_out_flux_Pa'])-float(b['Pt_Pa'])))<1e-8
        assert b['spectrum_status_in']==b['spectrum_status_out']=='0'
        if int(b['reduced_phase']):
            assert b['nwaves_in']==b['nwaves_out']=='2'
            assert b['nin_in']==b['nin_out']=='1'
        else:
            assert b['nwaves_in']==b['nwaves_out']=='3'
            assert b['nin_in']=='2' and b['nin_out']=='1'
    lines=(folder/'checkpoint.dat').read_text().splitlines();head=lines[1].split()
    assert len(head)==11 and len(lines)==int(head[0])+2
    assert int(head[3])==int(bc[-1]['step'])
    assert abs(float(head[5])*float(head[7])-float(bc[-1]['time_s']))<1e-12
    assert not (folder/'checkpoint.next').exists()
    fields=rows(folder/'fields.csv');fields=[f for f in fields if f['time_s']==fields[-1]['time_s']]
    states=[list(map(float,line.split())) for line in lines[2:]]
    nd=int(head[1]);rho_scale=float(head[9]);u_scale=float(head[8])
    n=int(head[0])-2*nd
    # Legacy finalization can write the terminal frame twice; compare the last full frame.
    assert len(fields)%n==0
    fields=fields[-n:]
    for f,u in zip(fields,states[nd:-nd]):
        alpha=float(f['alpha']);rl=float(f['rhol_kg_m3']);rg=float(f['rhog_kg_m3'])
        assert abs(u[4]*rho_scale-(1-alpha)*rl)<1e-10
        assert abs(u[5]*rho_scale-alpha*rg)<1e-10
        ul=float(f['ul_m_s']) if alpha<1 else 0
        ug=float(f['ug_m_s']) if alpha>0 else 0
        assert abs(u[6]*rho_scale*u_scale-((1-alpha)*rl*ul+alpha*rg*ug))<1e-9
    return dict(case=name,returncode=p.returncode,accepted=len(bc)-1,time_s=float(bc[-1]['time_s']),passed=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    cases=[('mass',{}),('local',dict(liquid_flow_basis="'LOCAL_VOLUME'",gas_flow_basis="'LOCAL_VOLUME'")),
       ('source',dict(boundary_method="'STEADY_TRACE'")),
       ('hydro',dict(jl=0,jg=0,mul=0,mug=0,single_phase_mode="'REDUCED'",source_scheme="'HYDROSTATIC'",perturbation_amplitude=0)),
       ('limited',dict(max_iter=2))]
    results=[check_case(a.binary.resolve(),a.output.resolve(),n,e) for n,e in cases]
    (a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
    for r in results:print('PASS',r['case'],r['accepted'])
if __name__=='__main__':main()
