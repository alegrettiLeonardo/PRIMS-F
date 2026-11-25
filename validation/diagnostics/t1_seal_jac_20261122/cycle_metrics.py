#!/usr/bin/env python3
"""P1.2: severe-slugging cycles of a T1 + LIQUID_BRIDGE run (monitor.csv, SI).

Events (pre-declared, no tuning):
  SEAL   onset of the liquid bridge at the junction: base gas mass flow becomes exactly
         0 (configuration H7; the H4 flux is > 0) after a sample with gas flowing;
  PEN    end of the bridge (gas penetration): first sample after SEAL with base gas > 0;
  BLOW   blowdown: outlet void fraction crosses ALPHA_OUT_BLOW upwards after PEN;
  PMAX   maximum of p_in between consecutive SEAL events.
Cycle k = [SEAL_k, SEAL_{k+1}). Period T_k = t(SEAL_{k+1}) - t(SEAL_k); also from PMAX.
Per cycle: time-mean p_in (trapezoid), p_max, p_min, excursion dp = p_max - p_min,
semi-amplitude A_p = dp/2, bridge duration, SEAL->BLOW, liquid/gas inventory at SEAL.
Usage: cycle_metrics.py RUN_DIR [OUT_JSON]
"""
import csv, json, sys
from pathlib import Path

ALPHA_OUT_BLOW = 0.05


def main():
    run = Path(sys.argv[1])
    r = list(csv.DictReader(open(run / 'monitor.csv')))
    t = [float(x['time_s']) for x in r]; p = [float(x['p_in_Pa']) for x in r]
    gb = [float(x['mdot_g_base_kg_s']) for x in r]; ao = [float(x['alpha_out']) for x in r]
    ml = [float(x['Ml_kg']) for x in r]; mg = [float(x['Mg_kg']) for x in r]
    ev = []
    state = 'H4' if gb[0] > 0 else 'H7'
    blown = state == 'H4'   # a run starting sealed (replay) still records its blowdown
    for i in range(1, len(t)):
        if gb[i] == 0.0 and gb[i - 1] > 0.0:
            ev.append(('SEAL', t[i], i)); state = 'H7'; blown = False
        elif state == 'H7' and gb[i] > 0.0:
            ev.append(('PEN', t[i], i)); state = 'H4'
        elif state == 'H4' and not blown and ao[i - 1] < ALPHA_OUT_BLOW <= ao[i]:
            ev.append(('BLOW', t[i], i)); blown = True
    seals = [e for e in ev if e[0] == 'SEAL']
    cycles = []
    for a, b in zip(seals[:-1], seals[1:]):
        i0, i1 = a[2], b[2]
        integ = sum(0.5 * (p[k] + p[k + 1]) * (t[k + 1] - t[k]) for k in range(i0, i1))
        seg = p[i0:i1 + 1]; kmax = max(range(i0, i1 + 1), key=lambda k: p[k])
        pen = [e for e in ev if e[0] == 'PEN' and a[1] <= e[1] < b[1]]
        blw = [e for e in ev if e[0] == 'BLOW' and a[1] <= e[1] < b[1]]
        cycles.append(dict(t_seal_s=a[1], period_s=b[1] - a[1], p_mean_Pa=integ / (b[1] - a[1]),
                           p_max_Pa=max(seg), p_min_Pa=min(seg), t_pmax_s=t[kmax], dp_Pa=max(seg) - min(seg),
                           semi_amplitude_Pa=0.5 * (max(seg) - min(seg)),
                           bridge_duration_s=(pen[0][1] - a[1]) if pen else None,
                           seal_to_blowdown_s=(blw[0][1] - a[1]) if blw else None,
                           Ml_at_seal_kg=ml[i0], Mg_at_seal_kg=mg[i0]))
    out = dict(run=str(run), t_end_s=t[-1], events=[(e[0], e[1]) for e in ev], n_complete_cycles=len(cycles),
               cycles=cycles)
    if len(cycles) >= 2:
        T = [c['period_s'] for c in cycles]; A = [c['semi_amplitude_Pa'] for c in cycles]
        m = lambda v: sum(v) / len(v)
        sd = lambda v: (sum((x - m(v)) ** 2 for x in v) / (len(v) - 1)) ** 0.5
        out['summary'] = dict(T_mean_s=m(T), T_sd_s=sd(T), A_mean_Pa=m(A), A_sd_Pa=sd(A),
                              dp_mean_Pa=m([c['dp_Pa'] for c in cycles]),
                              p_mean_Pa=m([c['p_mean_Pa'] for c in cycles]),
                              Ml_drift_kg=cycles[-1]['Ml_at_seal_kg'] - cycles[0]['Ml_at_seal_kg'],
                              Mg_drift_kg=cycles[-1]['Mg_at_seal_kg'] - cycles[0]['Mg_at_seal_kg'])
    txt = json.dumps(out, indent=1)
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
