#!/usr/bin/env python3
"""BAROTROPIC OSC refinement under the pre-registered v2 criteria (813da75).

Separate axes (space at dt=0.125, time at N51) and the combined pair
(N51/dt0.125 -> N101/dt0.0625), each run judged on its last ten complete
cycles; the per-run summaries and the pairwise v2 comparisons reuse
analyze_refinement.summarize and settling_criteria.compare unchanged.
Usage: baro_refinement_v2.py OUT_JSON OUT_MD RUN_DIR [RUN_DIR ...]
  Each RUN_DIR is named N<cells>_dt<dt>; runs without run_status.txt with
  termination_reason=COMPLETED are listed as not evaluable.
"""
import json, sys, hashlib
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'p0p1_20261119'))
sys.path.insert(0, str(HERE))
import analyze_refinement as ar  # noqa: E402
import settling_criteria as v2  # noqa: E402


def key(folder):
    n, dt = folder.name.split('_dt')
    return int(n[1:]), float(dt)


def main():
    out_json, out_md = Path(sys.argv[1]), Path(sys.argv[2])
    runs = {}
    for f in map(Path, sys.argv[3:]):
        st = (f / 'run_status.txt').read_text() if (f / 'run_status.txt').exists() else ''
        k = key(f)
        if 'termination_reason=COMPLETED' not in st:
            runs[k] = dict(run=str(f), evaluable=False, reason='not completed')
            continue
        s = ar.summarize(f)
        s.update(evaluable=True, n=k[0], dt_max_s=k[1])
        runs[k] = s
    pairs = []
    space = sorted(k for k in runs if k[1] == .125)
    time = sorted((k for k in runs if k[0] == 51), key=lambda k: -k[1])
    for axis, coh in [('space', space), ('time', time)]:
        pairs += [(axis, a, b) for a, b in zip(coh, coh[1:])]
    if (51, .125) in runs and (101, .0625) in runs:
        pairs.append(('combined', (51, .125), (101, .0625)))
    comps = []
    for axis, a, b in pairs:
        ra, rb = runs[a], runs[b]
        d = dict(axis=axis, coarse='N%d_dt%g' % a, fine='N%d_dt%g' % b)
        if not (ra.get('evaluable') and rb.get('evaluable') and 'period_s' in ra and 'period_s' in rb):
            d['status'] = 'NOT_EVALUABLE'
        else:
            for q in ['period_s', 'semi_amplitude_Pa', 'mean_Pa']:
                d[q + '_rel'] = abs(rb[q] - ra[q]) / max(abs(ra[q]), abs(rb[q]))
            d['both_settled'] = ra['settled'] and rb['settled']
            d['convergence_v2'] = v2.compare(ra['settling_v2'], rb['settling_v2'])
            d['status'] = 'PASS' if (d['period_s_rel'] <= .01 and d['semi_amplitude_Pa_rel'] <= .01
                                     and d['convergence_v2'].get('converged_v2_channels')) else 'FAIL'
        comps.append(d)
    slim = {('N%d_dt%g' % k): {x: v for x, v in r.items() if x not in ('phase_waveforms',)} for k, r in runs.items()}
    out_json.write_text(json.dumps(dict(criteria_v2_protocol_sha256=v2.protocol_sha256(),
                                        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                        runs=slim, comparisons=comps), indent=1, default=float))
    with out_md.open('w') as m:
        m.write('# BAROTROPIC OSC refinement, criteria v2\n\n| run | settled v2 | period s | semi-amplitude Pa | mean Pa | window s |\n|---|---|---:|---:|---:|---|\n')
        for k in sorted(runs):
            r = runs[k]
            if not r.get('evaluable') or 'period_s' not in r:
                m.write('| N%d_dt%g | NOT_EVALUABLE | - | - | - | - |\n' % k)
                continue
            m.write('| N%d_dt%g | %s | %.6f | %.2f | %.1f | %.1f-%.1f |\n' % (k + (r['settled'], r['period_s'], r['semi_amplitude_Pa'], r['mean_Pa'], *r['window_s'])))
        m.write('\n| axis | coarse | fine | period | amplitude | mean | v2 converged | status |\n|---|---|---|---:|---:|---:|---|---|\n')
        for d in comps:
            if d['status'] == 'NOT_EVALUABLE':
                m.write('| %s | %s | %s | - | - | - | - | NOT_EVALUABLE |\n' % (d['axis'], d['coarse'], d['fine']))
                continue
            m.write('| %s | %s | %s | %.4f | %.4f | %.4f | %s | %s |\n' % (d['axis'], d['coarse'], d['fine'], d['period_s_rel'],
                    d['semi_amplitude_Pa_rel'], d['mean_Pa_rel'], d['convergence_v2'].get('converged_v2_channels'), d['status']))
    print(out_md.read_text())


if __name__ == '__main__':
    main()
