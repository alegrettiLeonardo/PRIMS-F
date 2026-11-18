#!/usr/bin/env python3
"""Independent high-precision evaluation of Azevedo et al. (IMAACA 2015), p.121.

This reference calculation uses the printed equations directly and 80-digit
arithmetic. Colebrook is solved by bisection, independently of the safeguarded
Newton implementation in Fortran. Requires mpmath; not a runtime dependency.
"""
from pathlib import Path
import json
import mpmath as mp

mp.mp.dps = 80
SOURCE = {
    'authors': 'G. R. Azevedo; J. L. Baliño; K. P. Burr',
    'title': 'Linear stability analysis for severe slugging: sensitivity to void fraction correlations',
    'year': 2015, 'printed_page': 121, 'pdf_page': 135, 'equations': '116–128',
    'url': 'https://www.msc-les.org/proceedings/imaaca/2015/IMAACA2015.pdf',
    'pdf_sha256': '51891ab67e550804e732603faf37bfac53504ee2d1e36b4f589a4921aefbd171',
    'qualification': 'Transcription in this primary application paper; not a full audit of the 2014 original.',
}

def coefficients(values):
    a, rl, rg, ul, ug, th, d, eps, mu, sigma = map(mp.mpf, values)
    g = mp.mpf('9.81')
    jl, jg = (1-a)*ul, a*ug
    j = jl+jg
    re = rl*j*d/mu
    beta, quality = jg/j, rg*jg/(rg*jg+rl*jl)
    lo, hi = mp.mpf(0), mp.mpf(128)
    for _ in range(320):
        mid = (lo+hi)/2
        h = mid+4*mp.log10(eps/(mp.mpf('3.7')*d)+mp.mpf('1.256')*mid/re)
        if h > 0:
            hi = mid
        else:
            lo = mid
    y = (lo+hi)/2
    ftp = 1/y**2
    ratio = rg/rl
    c1 = mp.mpf('.2')*(1-mp.sqrt(ratio))*((mp.mpf('2.6')-beta)**mp.mpf('.15')-mp.sqrt(ftp))*(1-quality)**mp.mpf('1.5')
    c0 = (2-ratio**2)/(1+(re/1000)**2)
    c0 += (((1+ratio**2*mp.cos(th))/(1+mp.cos(th)))**((1-a)/5)+c1)/(1+(1000/re)**2)
    mustar = mu/mp.mpf('.001')
    c2 = (mp.mpf('.434')/mp.log10(mustar))**mp.mpf('.15') if mustar > 10 else mp.mpf(1)
    la = mp.sqrt(sigma/(g*(rl-rg)))/d
    c3 = (la/mp.mpf('.025'))**mp.mpf('.9') if la < mp.mpf('.025') else mp.mpf(1)
    frsg = mp.sqrt(rg/(rl-rg))*jg/mp.sqrt(g*d*mp.cos(th)) if th < 0 else None
    c4 = -1 if th < 0 and frsg <= mp.mpf('.1') else 1
    ud = (mp.mpf('.35')*mp.sin(th)+mp.mpf('.45')*mp.cos(th))*mp.sqrt(g*d*(rl-rg)/rl)*mp.sqrt(1-a)*c2*c3*c4
    return c0, ud, {'Re_tp': re, 'beta': beta, 'x': quality, 'f_tp': ftp,
                    'La': la, 'C2': c2, 'C3': c3, 'C4': c4, 'Fr_sg_downhill': frsg}

def main():
    cases = [
        ('.3','1000','3','.5','2','.6','.0508','.000046','.00089','.0728'),
        ('.7','998','20','.1','.8',mp.pi/2,'.0508','.000046','.00089','.0728'),
        ('.1','850','5','.02','.2','0','.5','.000046','.02','.03'),
        ('.8','1000','3','.2','3','-.3','.0508','.000046','.00089','.0728'),
        ('.2','1000','3','.1','.1','-.3','.0508','.000046','.00089','.0728'),
        ('.95','1000','3','0','.1',mp.pi/2,'.0508','.000046','.00089','.0728'),
        ('.05','1000','3','.3','0','.2','.0508','.000046','.00089','.0728'),
    ]
    records = []
    for number, values in enumerate(cases, 1):
        c0, ud, groups = coefficients(values)
        record = {'case': number, 'inputs_SI': dict(zip(
            ['alpha','rho_l','rho_g','u_l','u_g','theta_rad','D','epsilon','mu_l','sigma'],
            (mp.nstr(mp.mpf(v), 50) for v in values))),
            'C0': mp.nstr(c0, 50), 'Ud_m_s': mp.nstr(ud, 50),
            'groups': {k: None if v is None else mp.nstr(v, 50) for k, v in groups.items()}}
        records.append(record)
        print(number, mp.nstr(c0, 19), mp.nstr(ud, 19))
    path = Path(__file__).resolve().parents[1]/'tests'/'fixtures'/'bhagwat_azevedo2015.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({'source': SOURCE, 'precision_digits': mp.mp.dps, 'cases': records}, indent=2, ensure_ascii=False)+'\n')
    print(path)

if __name__ == '__main__':
    main()
