#!/usr/bin/env python3
"""Summarize junction_limits.f90 logs (P0.1) into JSON + Markdown tables.
Usage: summarize_limits.py OUT_JSON OUT_MD LABEL=LOG [LABEL=LOG ...]"""
import json, sys

def parse(path):
    d = dict(cells={}, faces={}, lim={})
    side = None
    for line in open(path):
        p = [x.strip() for x in line.rstrip('\n').split('|')]
        if p[0] == 'CELL':
            w, i, cid, st = int(p[1]), int(p[2]), int(p[3]), int(p[4])
            v = [float(x) for x in p[6].split()]
            keys = ['alpha', 'rho_g', 'u_g', 'j_g', 'U_g', 'p', 'u_l', 'j_l', 'C0', 'Vgj', 'lam1', 'lam2', 'lam3', 'lam_im']
            d['cells'].setdefault(w, {})[i] = dict(closure=cid, status=st, theta_deg=float(p[5]), **dict(zip(keys, v)))
        elif p[0] == 'FACE':
            w, fl, tr, st = int(p[1]), p[2], p[3], int(p[4])
            v = [float(x) for x in p[5].split()]
            rec = dict(status=st, p=v[0], Gl=v[1], Gg=v[2], a_pipe=v[3], a_riser=v[4])
            for blk in p[6:]:
                name, *vals = blk.split()
                rec[name] = [float(x) for x in vals]
            d['faces'].setdefault(w, {})['%s/%s' % (fl, tr)] = rec
        elif p[0] == 'LIMSIDE':
            side = p[1]
            d['lim'][side] = dict(p_Pa=float(p[4].split()[-1]), j_l=float(p[5].split()[-1]), lim1=[], lim2=[], lim3=[])
        elif p[0] == 'LIM1' and 'root alpha' in p[2]:
            d['lim'][side]['lim1'].append([float(x) for x in p[2].split()[-2:]])
        elif p[0] == 'LIM2':
            d['lim'][side]['lim2'].append(dict(alpha=float(p[2].split()[-1]), status=int(p[3].split()[-1]),
                                               u_g=float(p[4].split()[-1]), C0=float(p[6].split()[-2]), Vgj=float(p[6].split()[-1])))
        elif p[0] == 'LIM3':
            d['lim'][side]['lim3'].append(dict(j_g=float(p[2].split()[-1]), status=int(p[3].split()[-1]),
                                               nroots=int(p[4].split()[-1]), alpha=float(p[5].split()[-1])))
    return d

def main():
    out_json, out_md = sys.argv[1:3]
    res = {}
    for arg in sys.argv[3:]:
        lab, path = arg.split('=', 1)
        res[lab] = parse(path)
    json.dump(res, open(out_json, 'w'), indent=1)
    f = lambda x: '%.3g' % x
    with open(out_md, 'w') as m:
        m.write('# P0.1 junction constitutive limits at the critical states\n\n')
        m.write('## Adjacent cells (accepted state 0)\n\n| state | cell | closure | alpha | u_g m/s | j_g m/s | j_l m/s | C0 | Vgj m/s | eigenvalues m/s |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n')
        for lab, d in res.items():
            for i, c in sorted(d['cells'][0].items()):
                m.write('| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s, %s, %s |\n' % (lab, i, c['closure'], f(c['alpha']), f(c['u_g']), f(c['j_g']),
                        f(c['j_l']), f(c['C0']), f(c['Vgj']), f(c['lam1']), f(c['lam2']), f(c['lam3'])))
        m.write('\n(For the stratified pipe, C0 and Vgj are the apparent u_g/j and u_g-j, not closure inputs.)\n\n')
        for w in (0, 1):
            m.write('## Junction face, state %d (%s)\n\n' % (w, 'accepted' if w == 0 else 'last refused Newton trial'))
            m.write('| state | operator | status | Gg* kg/m2s | alpha* pipe | alpha* riser | Fconv_g | Fdiss_g | Roe raw Gg | transmitted Gg |\n|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n')
            for lab, d in res.items():
                for op, r in d['faces'].get(w, {}).items():
                    m.write('| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |\n' % (lab, op, r['status'], f(r['Gg']), f(r['a_pipe']), f(r['a_riser']),
                            f(r['Fconv'][1]) if r['Fconv'][1] > -1e5 else '-', f(r['Fdiss'][1]) if r['Fdiss'][1] > -1e5 else '-',
                            f(r['Froe'][1]) if r['Froe'][1] > -1e5 else '-', f(r['Ftrans'][1])))
            m.write('\n')
        m.write('## Limits (cell p and j_l of each side)\n\n| state | side | stagnant-gas roots alpha (u_g=0) | alpha as j_g->0+ (1e-8 m/s) | inversion status j_g=0 | status j_g<0 | lim u_g(alpha->0) m/s | C0, Vgj at alpha->0 |\n|---|---|---|---:|---:|---:|---:|---|\n')
        for lab, d in res.items():
            for s, L in d['lim'].items():
                r0 = [x for x in L['lim3'] if x['j_g'] == 0][0]
                rn = [x for x in L['lim3'] if x['j_g'] < 0]
                rp = [x for x in L['lim3'] if abs(x['j_g'] - 1e-8) < 1e-12][0]
                l2 = L['lim2'][-1]
                m.write('| %s | %s | %s | %s (n=%d) | %d (alpha=%s) | %s | %s | %s, %s |\n' % (lab, s, ', '.join('%.4g' % (0.5 * (a + b)) for a, b in L['lim1']) or 'none',
                        f(rp['alpha']), rp['nroots'], r0['status'], f(r0['alpha']), '/'.join(sorted(set(str(x['status']) for x in rn))),
                        f(l2['u_g']) if l2['status'] == 0 else '-', f(l2['C0']), f(l2['Vgj'])))
    print(open(out_md).read())

if __name__ == '__main__':
    main()
