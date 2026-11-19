#!/usr/bin/env python3
"""Prepare one tolerance-only study of a terminal, frozen jg+ extension."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path

from prepare_static_bc_sensitivity import FREEZE, HARNESS_SHA256


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--harness', type=Path, required=True)
    p.add_argument('--reference-plan', type=Path, required=True)
    p.add_argument('--workspace', type=Path, required=True)
    args = p.parse_args()
    if hashlib.sha256(args.harness.read_bytes()).hexdigest() != HARNESS_SHA256:
        raise ValueError('Executor identity differs')
    spec = importlib.util.spec_from_file_location('sealed_executor', args.harness)
    h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    old, seal = h.verify_seal(args.reference_plan.resolve())
    if seal['git_head'] != FREEZE or old['duration_s'] != 7200.:
        raise ValueError('Reference operator or duration differs')
    reference = next(j for j in old['jobs'] if j['name'] == 'D_OSC_bc_jg_plus_T7200')
    prior = Path(reference['output'])
    manifest_path = prior / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    mass_path = args.reference_plan.parent / 'analysis' / reference['name'] / 'accepted_conservation.json'
    mass = json.loads(mass_path.read_text())
    if manifest['status'] != 'COMPLETED' or not h.normal_completion(manifest['returncode'], manifest['run_status'], 7200.):
        raise ValueError('Prior run is not terminal at the requested duration')
    if manifest['plan_sha256'] != h.sha(args.reference_plan) or manifest['git_head'] != FREEZE:
        raise ValueError('Prior run belongs to another plan or operator')
    if (manifest['binary_sha256'] != seal['binary_sha256'] or
            manifest['source_sha256'] != seal['source_sha256'] or
            not manifest['binary_unchanged'] or not manifest['sources_unchanged']):
        raise ValueError('Prior operator identity differs')
    if mass['artifacts']['manifest.json']['sha256'] != h.sha(manifest_path):
        raise ValueError('Mass audit belongs to another run')
    if not mass['ledger']['arithmetic_contract_passes']:
        raise ValueError('Resolve the arithmetic identity before a tolerance study')
    if all(b['passes_inventory_1e6'] for b in mass['ledger']['balances'].values()):
        raise ValueError('This study requires the documented inventory-balance failure')
    original = Path(reference['deck'])
    if h.sha(original) != reference['deck_sha256']:
        raise ValueError('Reference deck changed')
    root = Path(old['root'])
    workspace = args.workspace.resolve()
    if root == workspace or root in workspace.parents:
        raise ValueError('Results must remain external')
    target = workspace / 'bc_static_tolres1e10_7200'
    if target.exists():
        raise ValueError('Preserve prior evidence: output cohort already exists')
    name = 'D_OSC_bc_jg_plus_T7200_tolres1e10'
    output = target / 'runs' / name
    deck = target / 'decks' / (name + '.nml')
    relative = os.path.relpath(output, root)
    text0 = original.read_text()
    before = h.parameters(text0)
    if before['tol_res'] != 1e-8 or before['tol_inc'] != 1e-8:
        raise ValueError('Reference tolerances differ from the authorized study')
    text = h.assign(h.assign(text0, 'tol_res', '1e-10'), 'out_dir', repr(relative))
    after = h.parameters(text)
    changes = h.parameter_changes(before, after)
    if changes != ['out_dir', 'tol_res'] or after['tol_res'] != 1e-10 or after['tol_inc'] != 1e-8:
        raise ValueError('Only residual tolerance and output directory may change')
    if any(before[k] != after[k] for k in h.PHYSICS):
        raise ValueError('A physical parameter changed')
    deck.parent.mkdir(parents=True)
    with deck.open('x') as stream:
        stream.write(text)
    job = dict(name=name, case='D_OSC', axis=reference['axis'], static_offset=reference['static_offset'],
               pipe_regime_mode=reference['pipe_regime_mode'], deck=str(deck), output=str(output),
               output_relative_to_root=relative, parameters=after, physics=reference['physics'],
               deck_sha256=h.sha(deck), reference_parameter_changes=changes,
               prior_run=str(prior), study_axis='tol_res', study_reference_value=1e-8, study_value=1e-10)
    inputs = dict(old['input_sha256'])
    for path in (args.reference_plan, args.reference_plan.parent / 'seal.json', original, manifest_path, mass_path):
        inputs[str(path.resolve())] = h.sha(path)
    plan = dict(schema='PRIMS-F.static_bc_residual_tolerance_plan.v1', state='PREPARED_UNSEALED',
                created_utc=h.now(), root=str(root), workspace=str(workspace), label=target.name,
                prepared_git_head=FREEZE, input_sha256=inputs, jobs=[job], environment=h.ENV,
                duration_s=7200., maximum_total_solvers=4, maximum_owned_workers=1,
                runner_sha256=h.sha(args.harness), preparer_sha256=h.sha(__file__),
                reference_plan=dict(path=str(args.reference_plan.resolve()), sha256=h.sha(args.reference_plan)),
                reference_manifest=old['reference_manifest'], reference_deck=old['reference_deck'],
                source=old['source'], criteria=old['criteria'],
                study_reference_manifest=dict(path=str(manifest_path), sha256=h.sha(manifest_path)),
                study_reference_mass_audit=dict(path=str(mass_path), sha256=h.sha(mass_path)),
                execution_gate='Parent authorized one tolerance-only study; wait for coordinated SS3 slot release.',
                interpretation='Independent numerical tolerance study; preserve tol1e-8 FAIL and do not replace the common-control sensitivity cohort.')
    h.dump_new(target / 'plan.json', plan)
    print(json.dumps(dict(plan=str(target / 'plan.json'), jobs=1, changed=changes, solvers_started=0)))


if __name__ == '__main__':
    main()
