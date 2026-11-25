#!/usr/bin/env python3
"""P0.5: integral momentum balance of the accepted steps (momentum_ledger.csv).

Backward Euler, per accepted step n (dt_n, every term at t_{n+1}, SI per unit area):
  pipe : I_p^{n+1} - I_p^n = dt (F3_in - F3_jL + S_grav,p + S_fric,p) + dt R_p
  riser: I_r^{n+1} - I_r^n = dt (F3_jR - F3_out + S_grav,r + S_fric,r) + dt R_r
  system = pipe + riser (+ dt (F3_jL - F3_jR): zero for T1, one flux at the face)
where R = sum ds R3/dt is the residual of the accepted Newton iterate (|R| <= tol_res).
E_P (cumulative) = [I]_{t0}^{tn} - sum dt*(fluxes + sources); the ledger row of the first
step needs I^0, which is reconstructed from the first row's own balance (excluded).
Normalization: N = sum dt*(|F3_in| + |F3_out| + |S_grav| + |S_fric|) over the same steps
(the gross momentum transported/produced), and also the hydrostatic scale rho_l g Z.
Usage: momentum_balance.py RUN_DIR [OUT_JSON]
"""
import csv, json, sys
from pathlib import Path


def main():
    run = Path(sys.argv[1])
    rows = list(csv.DictReader(open(run / 'momentum_ledger.csv')))
    f = lambda r, k: float(r[k])
    out = dict(run=str(run), steps=len(rows))
    if len(rows) < 3:
        print(json.dumps(out)); return
    acc = dict(pipe=0.0, riser=0.0, system=0.0, junction_jump=0.0)
    res = dict(pipe=0.0, riser=0.0)
    gross = 0.0; maxstep = dict(pipe=0.0, riser=0.0, system=0.0)
    terms = dict(F3_in=0.0, F3_out=0.0, grav_pipe=0.0, fric_pipe=0.0, grav_riser=0.0, fric_riser=0.0)
    bad = 0
    for a, b in zip(rows[:-1], rows[1:]):
        dt = f(b, 'dt_s')
        if int(b['status']) != 0:
            bad += 1
        dIp = f(b, 'I3_pipe_kg_m_s') - f(a, 'I3_pipe_kg_m_s')
        dIr = f(b, 'I3_riser_kg_m_s') - f(a, 'I3_riser_kg_m_s')
        rp = dt * (f(b, 'F3_in_Pa') - f(b, 'F3_junction_left_Pa') + f(b, 'S_grav_pipe_Pa') + f(b, 'S_fric_pipe_Pa'))
        rr = dt * (f(b, 'F3_junction_right_Pa') - f(b, 'F3_out_Pa') + f(b, 'S_grav_riser_Pa') + f(b, 'S_fric_riser_Pa'))
        ep = dIp - rp; er = dIr - rr
        acc['pipe'] += ep; acc['riser'] += er; acc['system'] += ep + er
        acc['junction_jump'] += dt * (f(b, 'F3_junction_left_Pa') - f(b, 'F3_junction_right_Pa'))
        res['pipe'] += dt * f(b, 'R3_pipe_Pa'); res['riser'] += dt * f(b, 'R3_riser_Pa')
        g = dt * (abs(f(b, 'F3_in_Pa')) + abs(f(b, 'F3_out_Pa')) + abs(f(b, 'S_grav_pipe_Pa')) + abs(f(b, 'S_fric_pipe_Pa'))
                  + abs(f(b, 'S_grav_riser_Pa')) + abs(f(b, 'S_fric_riser_Pa')))
        gross += g
        for k, v in (('pipe', ep), ('riser', er), ('system', ep + er)):
            maxstep[k] = max(maxstep[k], abs(v) / max(g, 1e-300))
        terms['F3_in'] += dt * f(b, 'F3_in_Pa'); terms['F3_out'] += dt * f(b, 'F3_out_Pa')
        terms['grav_pipe'] += dt * f(b, 'S_grav_pipe_Pa'); terms['fric_pipe'] += dt * f(b, 'S_fric_pipe_Pa')
        terms['grav_riser'] += dt * f(b, 'S_grav_riser_Pa'); terms['fric_riser'] += dt * f(b, 'S_fric_riser_Pa')
    out.update(t_start_s=f(rows[0], 'time_s'), t_end_s=f(rows[-1], 'time_s'), rows_with_status=bad,
               E_P_kg_m_s=acc, residual_integral_kg_m_s=res, gross_normalization_kg_m_s=gross,
               E_P_over_gross={k: acc[k] / gross for k in ('pipe', 'riser', 'system')},
               E_P_minus_residual_over_gross={k: (acc[k] - res[k]) / gross for k in ('pipe', 'riser')},
               max_step_E_over_step_gross=maxstep, integrated_terms_kg_m_s=terms,
               storage_change_kg_m_s=dict(pipe=f(rows[-1], 'I3_pipe_kg_m_s') - f(rows[0], 'I3_pipe_kg_m_s'),
                                          riser=f(rows[-1], 'I3_riser_kg_m_s') - f(rows[0], 'I3_riser_kg_m_s')))
    txt = json.dumps(out, indent=1)
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
