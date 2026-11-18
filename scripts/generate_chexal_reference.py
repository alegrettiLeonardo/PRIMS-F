"""Independent high-precision evaluation of the 1992 equations, not solver output.

Sources: Chexal et al. (1992), equations 3-2--41, printed pp261--269;
RELAP5/MOD3 (1995) vol4 p6-10 clarifies only C10's printed parenthesization.
No copyrighted PDFs/text accompany generated numerical fixtures.
"""
import json
import mpmath as mp

mp.mp.dps = 60
m = mp.mpf
P = dict(rl=m('998'), rg=m('1.2'), d=m('0.05254'), g=m('9.81'),
         ml=m('0.001'), mg=m('0.000018'), sigma=m('0.072'))

def coeff(a, jl, jg, p=P, horizontal=False, branch='high', flood=None):
    rl,rg,d,g,ml,mg,sigma=[p[k] for k in ['rl','rg','d','g','ml','mg','sigma']]
    rel,reg=rl*jl*d/ml,rg*jg*d/mg
    re=reg if (reg>rel or reg<0) else rel
    b1=min(m('.8'),1/(1+mp.exp(-re/60000)))
    k0=b1+(1-b1)*(rg/rl)**m('.25')
    r=(1+m('1.57')*rg/rl)/(1-b1)
    L=(m('1.125')*a**m('.6') if horizontal else
       m('1.15')*a**m('.45') if reg>=0 else m('1.05')*a**m('.25'))
    c0=min(L,m(1))/(k0+(1-k0)*a**r)
    if rl/rg<=18:
        c2=m('.4757')*mp.log(rl/rg)**m('.7')
    else:
        c5=mp.sqrt(150*rg/rl)
        c2=1 if c5>=1 else -1/mp.expm1(-c5/(1-c5))
    c7=(m('.09144')/d)**m('.6')
    c4=1 if c7>=1 else -1/mp.expm1(-c7/(1-c7))
    if jl>=0 and jg>=0:
        c3=max(m('.5'),2*mp.exp(-abs(rel)/60000))
    else:
        q=abs(rel)
        c10=2*mp.exp((q/350000)**m('.4'))-m('1.75')*q**m('.03')*mp.exp(-q/50000*(m('.0381')/d)**2)
        c10+=(m('.0381')/d)**m('.25')*q**m('.001')
        b2=(1/(1+m('.05')*q/350000))**m('.4')
        c3=2*(c10/2)**b2
        if branch=='low' and jl<0 and jg>=0:
            q=jl/flood
            c3=min(c3*q+2*(1+abs(rel)/60000)*(1-q),c3)
    v0=m('1.41')*((rl-rg)*sigma*g/rl**2)**m('.25')*c2*c3*c4
    c9=(1-a)**b1 if reg>=0 else min(m('.7'),(1-a)**m('.65'))
    if reg<0:c0=max(c0,v0*(1-a)**m('.2')/(abs(jl)+abs(jg)))
    if horizontal:c0*=1+a**m('.05')*(1-a)**2
    return c0,v0*c9

def R(a,jl,jg,branch='high',flood=None):
    c,v=coeff(a,jl,jg,branch=branch,flood=flood)
    return a*(c*(jl+jg)+v)-jg

def alpha_roots(jl,jg,branch='high',flood=None):
    fun=lambda a:R(a,jl,jg,branch,flood)
    roots=[]
    xs=[m(i)/1000 for i in range(1001)]
    for x,y in zip(xs,xs[1:]):
        fx,fy=fun(x),fun(y)
        if fx==0:roots.append(x)
        if fx*fy<0:roots.append(mp.findroot(fun,(x,y),solver='anderson'))
    if fun(1)==0:roots.append(m(1))
    return roots

def main():
    data={'precision_digits':mp.mp.dps,'parameters_SI':{k:str(v) for k,v in P.items()},'coefficients':[]}
    for name,a,jl,jg,h in [('vertical_up','.2','.4','.7',False),('horizontal','.2','.4','.7',True),
                         ('vertical_down','.6','-.5','-.2',False),('cc_high','.4','-.1','.1',False),
                         ('liquid_limit','0','.4','0',False),('gas_limit','1','0','.7',False)]:
        c,v=coeff(m(a),m(jl),m(jg),horizontal=h)
        data['coefficients'].append(dict(name=name,alpha=a,jl=jl,jg=jg,horizontal=h,c0=mp.nstr(c,55),vgj=mp.nstr(v,55)))
    data['countercurrent']=[]
    for gas in ['0.05','0.1','0.2']:
        jg=m(gas)
        fn=lambda a,l:(R(a,l,jg),mp.diff(lambda aa:R(aa,l,jg),a))
        a_star,l_star=mp.findroot(fn,(m('.5'),m('-.2')),tol=m('1e-50'),maxsteps=100)
        jl=l_star/2
        high=alpha_roots(jl,jg)[-1]
        low=alpha_roots(jl,jg,'low',l_star)[0]
        data['countercurrent'].append(dict(jg=gas,jl_ccfl=mp.nstr(l_star,55),alpha_ccfl=mp.nstr(a_star,55),
            jl=mp.nstr(jl,55),alpha_high=mp.nstr(high,55),alpha_low=mp.nstr(low,55),
            ccfl_R=mp.nstr(R(a_star,l_star,jg),4),ccfl_dR=mp.nstr(mp.diff(lambda a:R(a,l_star,jg),a_star),4)))
    data['cocurrent_root']=mp.nstr(alpha_roots(m('.4'),m('.7'))[0],55)
    data['zero_total_flux_roots_high_coefficient']=[mp.nstr(v,55) for v in alpha_roots(m('-.1'),m('.1'))]
    c,v=coeff(m('.4'),m('-1e-25'),m('.1'),branch='low',flood=m(data['countercurrent'][1]['jl_ccfl']))
    data['coefficients'].append(dict(name='cc_low_active_eq3_41',alpha='.4',jl='-1e-25',jg='.1',
        branch='low',c0=mp.nstr(c,55),vgj=mp.nstr(v,55)))
    print(json.dumps(data,indent=2))

if __name__=='__main__':main()
