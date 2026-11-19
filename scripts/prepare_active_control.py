#!/usr/bin/env python3
"""Prepare the qualified stable control with only opt-in Newton globalization changed."""
import argparse
import importlib.util
import json
import os
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--harness', type=Path, required=True)
    p.add_argument('--reference-plan', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--expected-commit', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    import hashlib
    runner_hash = hashlib.sha256(a.harness.read_bytes()).hexdigest()
    if runner_hash != 'ea973346fe0c9bec8a96df88abc98bd205d54b708668b717a2e47de4857e7a24':
        raise ValueError('Immutable executor changed')
    spec = importlib.util.spec_from_file_location('sealed_runner', a.harness)
    h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    root, out = a.root.resolve(), a.output.resolve()
    if out.exists() or root == out or root in out.parents:
        raise ValueError('Use a new destination outside the source tree')
    if h.command(['git', 'rev-parse', 'HEAD'], root) != a.expected_commit:
        raise ValueError('Source freeze differs')
    if h.command(['git', 'status', '--porcelain', '--untracked-files=no'], root):
        raise ValueError('Source tree is not clean')
    old, seal = h.verify_seal(a.reference_plan.resolve())
    job = next(job for job in old['jobs'] if job['name'] == 'D_Stable_tol1e10')
    reference = Path(job['output'])/'manifest.json'
    manifest = json.loads(reference.read_text())
    if not h.normal_completion(manifest['returncode'], manifest['run_status'], 2400):
        raise ValueError('Stable control incomplete')
    template = Path(job['deck'])
    if h.sha(template) != job['deck_sha256']:
        raise ValueError('Reference deck changed')
    text = template.read_text()
    before = h.parameters(text)
    if before['tol_res'] != 1e-10 or before['time_dim'] != 2400:
        raise ValueError('Unqualified stable reference')
    folder = out/'runs'/'D_Stable_ACTIVE_tol1e10'
    deck = out/'decks'/'D_Stable_ACTIVE_tol1e10.nml'
    relative = os.path.relpath(folder, root)
    text = h.assign(text, 'out_dir', repr(relative))
    text = h.assign(text, 'newton_globalization', repr('ACTIVE_SET'))
    # Adding a missing key through the historical assign() consumes the final
    # newline. A terminated Fortran namelist record must be emitted explicitly.
    text = text.rstrip()+'\n'
    values = h.parameters(text)
    if h.parameter_changes(before, values) != ['newton_globalization', 'out_dir']:
        raise ValueError('Control or physics changed')
    for key in h.PHYSICS:
        if values[key] != job['physics'][key]:
            raise ValueError('Physical input changed: '+key)
    deck.parent.mkdir(parents=True)
    with deck.open('x') as stream:
        stream.write(text)
    newjob = dict(name=folder.name, case='D_Stable', pipe_regime_mode='FIXED_STRATIFIED',
        deck=str(deck), deck_sha256=h.sha(deck), output=str(folder), output_relative_to_root=relative,
        parameters=values, physics=job['physics'], reference_parameter_changes=['newton_globalization', 'out_dir'])
    plan = dict(schema='PRIMS-F.active_control_plan.v1', state='PREPARED_UNSEALED', created_utc=h.now(),
        root=str(root), workspace=str(out.parent), label=out.name, prepared_git_head=a.expected_commit,
        jobs=[newjob], environment=h.ENV, duration_s=2400, maximum_total_solvers=4,
        input_sha256={str(template): h.sha(template)}, runner_sha256=runner_hash,
        preparer_sha256=h.sha(__file__),
        reference=dict(plan=str(a.reference_plan), plan_sha256=h.sha(a.reference_plan),
                       manifest=str(reference), manifest_sha256=h.sha(reference),
                       binary_sha256=seal['binary_sha256'], git_head=seal['git_head']),
        criteria=dict(inventory_relative_limit=1e-6, regions=['system', 'pipe', 'riser'], phases=['l', 'g']),
        interpretation='Same qualified tol_res=1e-10 stable physics and controls; new opt-in globalization; no calibration.')
    h.dump_new(out/'plan.json', plan)
    print(json.dumps(dict(plan=str(out/'plan.json'), solvers_started=0, changes=newjob['reference_parameter_changes'])))


if __name__ == '__main__':
    main()
