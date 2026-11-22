"""Closure contract: the closure requested is the closure executed (P0.3).

  python3 tests/test_closure_contract_cli.py --binary bin/prims_f_check_all --output tmp/ci_closure_contract

Positive and NEGATIVE cases, verified by scripts/verify_closure_contract.py:
  default deck (no key)       requested DEFAULT -> effective 5, data are Chexal        PASS
  closure_riser=3             effective 3, data are NOT Chexal                          PASS
  comment 'closure_riser=3' + body 5, run labelled _C3_   (the 21/11 campaign bug)    detected FAIL
  same deck labelled _C5_                                                              PASS
  closure_riser=9             unknown model refused by the solver before any output
  E1 + Chexal                 incompatible configuration refused
  contract tampered after the run (effective edited to 3)                            detected FAIL
"""
from pathlib import Path
import argparse, json, os, shutil, subprocess, sys

ROOT = Path(__file__).resolve().parents[1]
VERIFY = ROOT / 'scripts/verify_closure_contract.py'
BASE = dict(nodes=7, n_riser=7, X=6.435, Z=9.886, Lp=57.4, beta_deg=2., D=.05254, eps=.000046, rho_l0=997.05,
            P_l0=101325., Cl=1498., Rg=287., mu_l=.00089, mu_g=.0000185, sigma=.0728, G=9.81, Ps=206000., jl=.547,
            jg=1.16, Tgas=293.15, P_standard=101300., T_standard=293.15, Cg_fixed=343., gas_eos='BAROTROPIC',
            gas_flow_basis='STANDARD', time_dim=.1, tol=1e-12, tola=1e-13, tol_res=1e-10, tol_inc=1e-10, newton_max=25,
            dt_max_s=.01, dt_start_s=.001, dt_min_s=1e-8, perturbation_amplitude=0., output_interval_s=.01,
            field_interval_s=.05, regime='CONTRACT_TEST')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    binary, out = a.binary.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    records = []

    def run(name, patch, header=''):
        dest = out / name; dest.mkdir()
        params = dict(BASE, **patch); params['out_dir'] = str(dest)
        val = lambda v: "'" + v + "'" if isinstance(v, str) else format(v, '.17g') if isinstance(v, float) else str(v)
        deck = out / (name + '.nml')
        deck.write_text(header + '&PARAMS\n' + ''.join(' ' + k + '=' + val(v) + ',\n' for k, v in params.items()) + '/\n')
        r = subprocess.run([str(binary), str(deck)], capture_output=True, text=True, timeout=300, env=env)
        (out / (name + '.log')).write_text(r.stdout + r.stderr)
        return dest, deck, r.returncode, r.stdout + r.stderr

    def verify(dest, deck, *extra):
        r = subprocess.run([sys.executable, str(VERIFY), str(dest), '--deck', str(deck), *extra], capture_output=True, text=True)
        status = json.loads((dest / 'CLOSURE_CONTRACT_STATUS.json').read_text()) if (dest / 'CLOSURE_CONTRACT_STATUS.json').exists() else {}
        return r.returncode, status

    # 1. default -> requested DEFAULT, effective 5
    d, deck, rc, log = run('default_C5', {})
    assert rc == 0, log[-1500:]
    c = (d / 'closure_contract.txt').read_text()
    assert 'closure_riser_requested=DEFAULT' in c and 'closure_riser_effective=5' in c and 'CHEXAL_1992_AUTO' in c
    vr, st = verify(d, deck); assert vr == 0 and st['status'] == 'PASS', st
    records.append(dict(case='default', verify=st['status']))
    # 2. explicit Bhagwat
    d, deck, rc, log = run('bhagwat_C3', dict(closure_riser=3))
    assert rc == 0, log[-1500:]
    vr, st = verify(d, deck); assert vr == 0 and st['status'] == 'PASS', st
    records.append(dict(case='bhagwat', verify=st['status']))
    # 3. the historical bug: a comment mentions closure_riser=3, the body selects 5, label says C3
    hdr = '! provenance: legacy Bhagwat = closure_riser=3.\n'
    d, deck, rc, log = run('commented_C3', dict(closure_riser=5), header=hdr)
    assert rc == 0, log[-1500:]
    vr, st = verify(d, deck)
    # executed and data are Chexal (consistent); the LABEL is what lies, as on 21/11
    assert vr == 1 and st['failed'] == ['label_riser == contract_effective'], st
    records.append(dict(case='comment_trap_mislabelled', verify=st['status'], failed=st['failed']))
    # 4. the same deck correctly labelled
    d, deck, rc, log = run('commented_C5', dict(closure_riser=5), header=hdr)
    vr, st = verify(d, deck); assert vr == 0 and st['status'] == 'PASS', st
    records.append(dict(case='comment_trap_labelled', verify=st['status']))
    # 5. unknown model
    d, deck, rc, log = run('unknown_C9', dict(closure_riser=9))
    assert rc != 0 and 'Unsupported riser closure ID' in log and not (d / 'closure_contract.txt').exists()
    records.append(dict(case='unknown_model', refused=True))
    # 6. incompatible configuration
    d, deck, rc, log = run('incompatible_E1_C5', dict(junction_flux='PHASE_UPWIND'))
    assert rc != 0 and 'not qualified with CHEXAL_1992' in log and not (d / 'closure_contract.txt').exists()
    records.append(dict(case='incompatible_E1', refused=True))
    # 6b. junction transmission T1 is recorded and verified; unknown value refused
    d, deck, rc, log = run('t1_C5', dict(junction_transmission='T1'))
    assert rc == 0, log[-1500:]
    assert 'junction_transmission=T1' in (d / 'closure_contract.txt').read_text()
    vr, st = verify(d, deck); assert vr == 0 and st['status'] == 'PASS', st
    assert any(c['check'] == 'deck_transmission == contract' and c['ok'] for c in st['checks'])
    records.append(dict(case='transmission_T1', verify=st['status']))
    d, deck, rc, log = run('t3_C5', dict(junction_transmission='T3'))
    assert rc != 0 and 'Unknown junction_transmission' in log and not (d / 'closure_contract.txt').exists()
    records.append(dict(case='unknown_transmission', refused=True))
    # 6c. seal extension (experimental): recorded and verified with T1; refused with T2
    d, deck, rc, log = run('t1seal_C5', dict(junction_transmission='T1', junction_seal='LIQUID_BRIDGE', time_dim=0.01))
    assert rc == 0, log[-1500:]
    assert 'junction_seal=LIQUID_BRIDGE' in (d / 'closure_contract.txt').read_text()
    vr, st = verify(d, deck); assert vr == 0 and st['status'] == 'PASS', st
    assert any(c['check'] == 'deck_seal == contract' and c['ok'] for c in st['checks'])
    records.append(dict(case='seal_T1', verify=st['status']))
    d, deck, rc, log = run('t2seal_C5', dict(junction_seal='LIQUID_BRIDGE'))
    assert rc != 0 and 'requires junction_transmission=T1' in log and not (d / 'closure_contract.txt').exists()
    records.append(dict(case='seal_with_T2', refused=True))
    # 7. contract tampered after the run
    src = out / 'default_C5'; tam = out / 'tampered_C5'
    shutil.copytree(src, tam)
    t = (tam / 'closure_contract.txt').read_text().replace('closure_riser_effective=5', 'closure_riser_effective=3')
    (tam / 'closure_contract.txt').write_text(t)
    vr, st = verify(tam, out / 'default_C5.nml')
    assert vr == 1 and 'resolved_riser == contract_effective' in st['failed'] and 'metadata_riser == contract_effective' in st['failed'], st
    records.append(dict(case='tampered_contract', verify=st['status'], failed=st['failed']))
    (out / 'manifest.json').write_text(json.dumps(records, indent=1) + '\n')
    print('PASS', len(records), 'closure-contract cases (3 negative detections, 4 refusals)')


if __name__ == '__main__':
    main()
