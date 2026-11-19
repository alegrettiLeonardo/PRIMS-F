"""Independent 60-digit references for the observed near-zero velocity case.

Uses the primary-equation reference evaluator, never the production Fortran.
Checks arbitrary-precision differentiation against a local analytic derivative.
Sources: Chexal1992 pp.261-269, especially Eqs.3-2,3-3,3-15,3-21,3-24,3-34.
"""
import json
import mpmath as mp
from generate_chexal_reference import coeff

mp.mp.dps = 60
m = mp.mpf
p = dict(rl=m('998'),rg=m('2.3'),d=m('.05254'),g=m('9.81'),
         ml=m('.001'),mg=m('.000018'),sigma=m('.072'))
a=m('.4');ug=m('.5');theta=m('1.151159')

def values(ul):
    jl=(1-a)*ul;jg=a*ug
    cv,vv=coeff(a,jl,jg,p=p)
    ch,vh=coeff(a,jl,jg,p=p,horizontal=True)
    fr=2*theta/mp.pi
    c=fr*cv+(1-fr)*ch;v=fr*vv+(1-fr)*vh
    return c,v,a*(c*(jl+jg)+v)-jg

def text(x):return mp.nstr(x,55)

def main():
    out={'precision_digits':60,'mpmath_version':mp.__version__,
         'meaning':'partial R_SI / partial ul_SI at fixed alpha,rhol,rhog,ug,theta and properties',
         'residual':'alpha*(C0*((1-alpha)*ul+alpha*ug)+Vgj)-alpha*ug',
         'parameters_SI':{k:str(v) for k,v in p.items()},'alpha':str(a),'ug':str(ug),
         'theta_horizontal_rad':str(theta),'cases':[]}
    gas_re=p['rg']*(a*ug)*p['d']/p['mg']
    cross=gas_re*p['ml']/(p['rl']*(1-a)*p['d'])
    floor=60000*mp.log(4)*p['ml']/(p['rl']*(1-a)*p['d'])
    out['positive_velocity_piece_boundaries']={'Rel_equals_Reg_ul':text(cross),'C3_upflow_floor_ul':text(floor)}
    for s in ['2.2586e-8','1e-14']:
        ul=m(s);c,v,r=values(ul)
        numeric=mp.diff(lambda u:values(u)[2],ul)
        analytic=a*(1-a)*(c-v*p['rl']*p['d']/(p['ml']*60000))
        assert abs(numeric-analytic)<m('1e-50')
        assert 0<ul<cross and ul<floor
        out['cases'].append({'ul':str(ul),'jl':str((1-a)*ul),'Re_l':text(p['rl']*(1-a)*ul*p['d']/p['ml']),
            'Re_g':text(gas_re),'C0':text(c),'Vgj':text(v),'R':text(r),
            'dR_dul_mpmath_diff':text(numeric),'dR_dul_analytic':text(analytic),
            'independent_derivative_difference':text(numeric-analytic)})
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
