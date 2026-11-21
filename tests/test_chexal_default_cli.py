"""CHEXAL_1992 as development default of the WHOLE riser: CLI contract.

  python3 tests/test_chexal_default_cli.py --binary bin/prims_f_check_all --output tmp/ci_chexal_default

Checks (outputs are new; nothing is overwritten):
  * a deck WITHOUT closure_riser runs with CHEXAL_1992 (model banner and
    run_metadata) and its saved riser states satisfy the Chexal relation of the
    independent transcription (validation/diagnostics/chexal_riser_20261121);
  * closure_riser=3 (Bhagwat) remains selectable and runs;
  * options derived/qualified only with Bhagwat are REFUSED with Chexal
    (E1 junction_flux, E2 junction_riser_limit, ACTIVE_SET), never silently inactive;
  * every primary benchmark deck assigns closure_riser=5 exactly once outside
    comments, and the archived Bhagwat decks assign 3.
"""
from pathlib import Path
import argparse, json, os, re, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'validation/diagnostics/chexal_riser_20261121'))
import chexal92_ref as cx  # noqa: E402

BASE = dict(nodes=7, n_riser=7, X=6.435, Z=9.886, Lp=57.4, beta_deg=2., D=.05254, eps=.000046, rho_l0=997.05,
            P_l0=101325., Cl=1498., Rg=287., mu_l=.00089, mu_g=.0000185, sigma=.0728, G=9.81, Ps=206000., jl=.547,
            jg=1.16, Tgas=293.15, P_standard=101300., T_standard=293.15, Cg_fixed=343., gas_eos='BAROTROPIC',
            gas_flow_basis='STANDARD', time_dim=.2, tol=1e-12, tola=1e-13, tol_res=1e-10, tol_inc=1e-10, newton_max=25,
            dt_max_s=.01, dt_start_s=.001, dt_min_s=1e-8, perturbation_amplitude=0., output_interval_s=.01,
            field_interval_s=.05, regime='CHEXAL_DEFAULT_TEST')


def body(text):
    return '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    binary, out = a.binary.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    records = []

    def run(name, patch):
        dest = out / name; dest.mkdir()
        params = dict(BASE, **patch); params['out_dir'] = str(dest)
        val = lambda v: "'" + v + "'" if isinstance(v, str) else format(v, '.17g') if isinstance(v, float) else str(v)
        deck = dest / 'input.nml'
        deck.write_text('&PARAMS\n' + ''.join(' ' + k + '=' + val(v) + ',\n' for k, v in params.items()) + '/\n')
        with (dest / 'stdout.log').open('w') as f:
            r = subprocess.run([str(binary), str(deck)], stdout=f, stderr=subprocess.STDOUT, timeout=300,
                               env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1'))
        records.append(dict(name=name, returncode=r.returncode, patch=patch))
        return dest, r.returncode, (dest / 'stdout.log').read_text()

    # 1. default (no closure_riser in the deck) is CHEXAL_1992 everywhere in the riser
    dest, rc, log = run('default_no_closure_key', {})
    assert rc == 0, log[-2000:]
    assert '[MODEL] CHEXAL_1992 air-water' in log, 'default riser closure banner'
    meta = (dest / 'run_metadata.txt').read_text()
    assert re.search(r'CLOSURE_RISER\s*=\s*5\b', meta, re.I), 'run_metadata closure_riser'
    import numpy as np
    f = np.genfromtxt(dest / 'fields.csv', delimiter=',', names=True)
    rows = f[(f['s_m'] > BASE['Lp']) & (f['theta_rad'] >= 0)]
    assert rows.size > 0
    worst = 0.0
    for r in rows:
        al = float(r['alpha']); jl = (1 - al) * float(r['ul_m_s']); jg = al * float(r['ug_m_s'])
        c0, v = cx.coefficients(al, float(r['rhol_kg_m3']), float(r['rhog_kg_m3']), float(r['theta_rad']), jl, jg,
                                d=BASE['D'], g=BASE['G'], ml=BASE['mu_l'], mg=BASE['mu_g'], sig=BASE['sigma'])
        worst = max(worst, abs(al * (c0 * (jl + jg) + v) - jg) / max(abs(jg), abs(al * c0 * (jl + jg)), abs(al * v)))
    assert worst < 1e-9, ('saved riser states do not satisfy Chexal', worst)
    records[-1]['riser_states'] = int(rows.size); records[-1]['chexal_relation_max_rel'] = worst

    # 2. legacy Bhagwat remains explicitly selectable
    dest, rc, log = run('explicit_bhagwat_3', dict(closure_riser=3))
    assert rc == 0, log[-2000:]
    assert '[MODEL] CHEXAL' not in log

    # 3. Bhagwat-only experimental options are refused with Chexal (5, 6, 7)
    refusals = [('E1_phase_upwind', dict(junction_flux='PHASE_UPWIND'), 'junction_flux (E1) is not qualified with CHEXAL_1992'),
                ('E2_alpha0', dict(junction_riser_limit='ALPHA0'), 'junction_riser_limit (E2) is not qualified with CHEXAL_1992'),
                ('active_set', dict(newton_globalization='ACTIVE_SET'), 'ACTIVE_SET admissibility margin is qualified only for closure_riser=3')]
    for name, patch, msg in refusals:
        for cid in (None, 6, 7):
            p = dict(patch) if cid is None else dict(patch, closure_riser=cid)
            dest, rc, log = run(f'{name}_C{cid or "default"}', p)
            assert rc != 0 and msg in log, (name, cid, log[-800:])
            assert not (dest / 'fields.csv').exists(), 'refused before any transient output'

    # 4. deck contract: primary decks select 5 exactly once outside comments; archive keeps 3
    for d in sorted((ROOT / 'input/benchmarks').glob('D_*.nml')):
        hits = re.findall(r'\bclosure_riser\s*=\s*(\d+)', body(d.read_text()), re.I)
        assert hits == ['5'], (d.name, hits)
        old = re.findall(r'\bclosure_riser\s*=\s*(\d+)', body((ROOT / 'input/benchmarks/bhagwat' / d.name).read_text()), re.I)
        assert old == ['3'], (d.name, old)
    (out / 'manifest.json').write_text(json.dumps(records, indent=1) + '\n')
    print('PASS', len(records), 'CHEXAL default CLI executions')


if __name__ == '__main__':
    main()
