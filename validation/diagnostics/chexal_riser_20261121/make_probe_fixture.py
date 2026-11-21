#!/usr/bin/env python3
"""Build an offline_api fixture (probe_parameters.dat + last_rejected_state.dat)
from a CHEXAL run, for the junction Riemann reference (P0.2).

  source = DUMP     the solver's own last_rejected_state.dat (accepted Uold and
                    the refused trial), copied byte for byte;
  source = <time>   an ACCEPTED saved state of fields.csv (prepare_state.py
                    reconstruction), written with old = trial = U.
probe_parameters.dat: omega_u, Ps, jl, jg, mul, mug, 0 (deck values; offline_api
fixes the D geometry, 51+51 cells and BAROTROPIC c_g=343 m/s, so only BAROTROPIC
51+51 runs are accepted here).
Usage: make_probe_fixture.py RUN_DIR SOURCE OUT_DIR
"""
import re, shutil, subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def deck_value(text, key, cast=float):
    text = '\n'.join(l for l in text.splitlines() if not l.lstrip().startswith('!'))
    return cast(re.search(r'\b' + key + r'\s*=\s*([^,\n]+)', text, re.I).group(1).strip().strip("'"))


def main():
    run, src, out = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
    deck = next(run.parent.parent.glob('decks/' + run.name + '.nml')).read_text()
    assert deck_value(deck, 'gas_eos', str) == 'BAROTROPIC' and int(deck_value(deck, 'nodes')) == 51, 'offline_api fixes BARO 51+51'
    out.mkdir(parents=True, exist_ok=True)
    jl, jg = deck_value(deck, 'jl'), deck_value(deck, 'jg')
    (out / 'probe_parameters.dat').write_text('%.17e %.17e %.17e %.17e %.17e %.17e 0\n' % (
        max(jl, jg), deck_value(deck, 'Ps'), jl, jg, deck_value(deck, 'mul'), deck_value(deck, 'mug')))
    if src == 'DUMP':
        shutil.copyfile(run / 'last_rejected_state.dat', out / 'last_rejected_state.dat')
    else:
        tmp = out / 'state'
        subprocess.run([sys.executable, str(HERE / 'prepare_state.py'), str(run), src, str(tmp)], check=True, capture_output=True)
        rows = [l.split() for l in (out / 'state.state').read_text().splitlines()]
        with open(out / 'last_rejected_state.dat', 'w') as f:
            f.write('# t_nd dt_nd residual code; then i,theta,width,Uold(3),Utrial(3),residual(3) [accepted saved state, old=trial]\n')
            f.write('0.0 1.0e-3 0.0 0\n')
            for k, r in enumerate(rows):
                f.write('%d %s %s %s %s %s %s %s %s 0 0 0\n' % (k + 3, r[0], r[1], r[2], r[3], r[4], r[2], r[3], r[4]))
    print(out)


if __name__ == '__main__':
    main()
