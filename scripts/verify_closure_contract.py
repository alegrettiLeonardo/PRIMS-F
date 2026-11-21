#!/usr/bin/env python3
"""Audit that the closure REQUESTED for a run is the closure EXECUTED.

Cross-checks, for one run directory:
  deck       closure_riser / closure_pipe assignments on NON-COMMENT lines of the
             input deck (absent -> DEFAULT);
  label      the `_C<k>_` token of the run name, if any (or --expect-riser k);
  contract   closure_contract.txt written by the solver before the transient:
             requested (DEFAULT or value), effective, and the number of cells to
             which that closure was actually dispatched (must be all of them);
  metadata   CLOSURE_RISER / CLOSURE_PIPE in run_metadata.txt;
  checkpoint closure ids in the checkpoint.dat header, if present;
  data       accepted riser states in fields.csv: the Chexal (1992) relation is
             recomputed with the independent transcription
             validation/diagnostics/chexal_riser_20261121/chexal92_ref.py. A
             Chexal run (5-7) must satisfy it in every state (R_rel <= 1e-9); any
             other riser closure must NOT satisfy it (<= 1 % of the states), so a
             mislabelled run is detected from the data themselves.
Writes CLOSURE_CONTRACT_STATUS.json in the run directory (or --out) and exits
1 on any mismatch. Defaults resolved by the solver: closure_riser 5, pipe 1.
Usage: verify_closure_contract.py RUN_DIR [--deck DECK] [--expect-riser K] [--out JSON] [--no-data]
"""
import argparse, hashlib, json, math, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'validation/diagnostics/chexal_riser_20261121'))
DEFAULT_RISER, DEFAULT_PIPE = 5, 1


def body(text):
    return '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))


def deck_ints(text, key):
    return [int(x) for x in re.findall(r'\b' + key + r'\s*=\s*(-?\d+)', body(text), re.I)]


def deck_float(text, key, default=None):
    m = re.search(r'\b' + key + r'\s*=\s*([-+0-9.eEdD]+)', body(text), re.I)
    return float(m.group(1).replace('d', 'e').replace('D', 'e')) if m else default


def kv(path):
    d = {}
    for line in Path(path).read_text().splitlines():
        if '=' in line:
            k, v = line.split('=', 1); d[k.strip()] = v.strip()
    return d


def data_check(run, deck_text, effective):
    import numpy as np
    import chexal92_ref as cx
    f = np.genfromtxt(run / 'fields.csv', delimiter=',', names=True)
    lp = deck_float(deck_text, 'Lp', 57.4)
    kw = dict(d=deck_float(deck_text, 'D', .05254), g=deck_float(deck_text, 'G', 9.81), ml=deck_float(deck_text, 'mu_l', 8.9e-4),
              mg=deck_float(deck_text, 'mu_g', 1.85e-5), sig=deck_float(deck_text, 'sigma', .0728))
    rows = f[(f['s_m'] > lp) & (f['theta_rad'] >= 0)]
    rel = []
    for r in rows:
        a = float(r['alpha']); jl = (1 - a) * float(r['ul_m_s']); jg = a * float(r['ug_m_s'])
        try:
            c0, v = cx.coefficients(a, float(r['rhol_kg_m3']), float(r['rhog_kg_m3']), float(r['theta_rad']), jl, jg, **kw)
        except cx.DomainError:
            continue
        sc = max(abs(jg), abs(a * c0 * (jl + jg)), abs(a * v))
        if sc > 0:
            rel.append(abs(a * (c0 * (jl + jg) + v) - jg) / sc)
    if not rel:
        return dict(ok=False, reason='no evaluable riser states')
    frac = sum(q <= 1e-9 for q in rel) / len(rel)
    chexal = effective in (5, 6, 7)
    ok = frac == 1.0 if chexal else frac <= 0.01
    return dict(ok=ok, states=len(rel), fraction_satisfying_chexal=frac, max_rel=max(rel),
                rule='all states satisfy Chexal' if chexal else 'at most 1% of states satisfy Chexal')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run', type=Path)
    ap.add_argument('--deck', type=Path)
    ap.add_argument('--expect-riser', type=int)
    ap.add_argument('--out', type=Path)
    ap.add_argument('--no-data', action='store_true')
    a = ap.parse_args()
    run = a.run.resolve()
    checks, fail = [], []

    def check(name, ok, detail):
        checks.append(dict(check=name, ok=bool(ok), detail=detail))
        if not ok:
            fail.append(name)

    cpath = run / 'closure_contract.txt'
    check('contract_present', cpath.exists(), str(cpath))
    if not cpath.exists():
        return finish(a, run, checks, fail, {})
    c = kv(cpath)
    deck_path = a.deck or Path(c.get('input_file', ''))
    deck_text = deck_path.read_text() if deck_path and deck_path.exists() else None
    check('deck_readable', deck_text is not None, str(deck_path))
    eff_r, eff_p = int(c['closure_riser_effective']), int(c['closure_pipe_effective'])
    if deck_text is not None:
        dr, dp = deck_ints(deck_text, 'closure_riser'), deck_ints(deck_text, 'closure_pipe')
        check('deck_single_riser_assignment', len(dr) <= 1, dr)
        req_r = str(dr[0]) if dr else 'DEFAULT'
        check('deck_requested_riser == contract_requested', req_r == c['closure_riser_requested'], [req_r, c['closure_riser_requested']])
        req_p = str(dp[0]) if dp else 'DEFAULT'
        check('deck_requested_pipe == contract_requested', req_p == c['closure_pipe_requested'], [req_p, c['closure_pipe_requested']])
        expected_r = int(dr[0]) if dr else DEFAULT_RISER
        check('resolved_riser == contract_effective', expected_r == eff_r, [expected_r, eff_r])
    check('all_riser_cells_dispatched_effective', c['riser_cells'] == c['riser_cells_dispatched_effective'],
          [c['riser_cells'], c['riser_cells_dispatched_effective']])
    check('all_pipe_cells_dispatched_effective', c['pipe_cells'] == c['pipe_cells_dispatched_effective'],
          [c['pipe_cells'], c['pipe_cells_dispatched_effective']])
    m = re.search(r'_C(\d+)(?:_|$)', run.name)
    label = a.expect_riser if a.expect_riser is not None else (int(m.group(1)) if m else None)
    if label is not None:
        check('label_riser == contract_effective', label == eff_r, [label, eff_r])
    meta = run / 'run_metadata.txt'
    if meta.exists():
        mm = re.search(r'CLOSURE_RISER\s*=\s*(-?\d+)', meta.read_text(), re.I)
        check('metadata_riser == contract_effective', mm is not None and int(mm.group(1)) == eff_r, mm.group(1) if mm else None)
    ck = run / 'checkpoint.dat'
    if ck.exists():
        head = ck.read_text().splitlines()[0]
        mk = re.search(r'closure_riser=(\d+)', head)
        if mk:
            check('checkpoint_riser == contract_effective', int(mk.group(1)) == eff_r, mk.group(1))
            mp = re.search(r'closure_pipe=(\d+)', head)
            check('checkpoint_pipe == contract_effective', mp is not None and int(mp.group(1)) == eff_p, mp.group(1) if mp else None)
        else:
            check('checkpoint_header_has_contract', False, 'legacy checkpoint header without closure ids')
    data = None
    if not a.no_data and (run / 'fields.csv').exists() and deck_text is not None:
        data = data_check(run, deck_text, eff_r)
        check('data_identifies_effective_closure', data['ok'], data)
    return finish(a, run, checks, fail, c, data)


def finish(a, run, checks, fail, c, data=None):
    out = dict(run=str(run), status='PASS' if not fail else 'FAIL', failed=fail, checks=checks,
               contract_sha256=hashlib.sha256((run / 'closure_contract.txt').read_bytes()).hexdigest()
               if (run / 'closure_contract.txt').exists() else None,
               closure_riser_effective=c.get('closure_riser_effective'), closure_riser_name=c.get('closure_riser_name'))
    path = a.out or (run / 'CLOSURE_CONTRACT_STATUS.json')
    path.write_text(json.dumps(out, indent=1) + '\n')
    print(out['status'], run.name, ('failed: ' + ', '.join(fail)) if fail else '')
    return 0 if not fail else 1


if __name__ == '__main__':
    sys.exit(main())
