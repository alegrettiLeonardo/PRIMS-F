#!/usr/bin/env python3
"""Prepare six static one-axis OSC sensitivities; never launch a solver.

Baliño2010 Table3/case13 signal deviations set offsets, not probabilities or
time-dependent excitations. The sealed nominal operator and controls persist.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import math
import os

FREEZE = '69813ee70551fa4a3f075919aeafc45de9693f11'
HARNESS_SHA256 = 'ea973346fe0c9bec8a96df88abc98bd205d54b708668b717a2e47de4857e7a24'
VARIANTS = (('Ps_minus', 'ps', 201000.), ('Ps_plus', 'ps', 211000.),
            ('jg_minus', 'jg', 1.10), ('jg_plus', 'jg', 1.22),
            ('jl_minus', 'jl', .544), ('jl_plus', 'jl', .550))


def validate_nominal(values):
    required = dict(nodes=51., n_riser=51., dt_max_s=.125, time_dim=3600.,
                    tol_res=1e-8, tol_inc=1e-8, ps=206000., jg=1.16, jl=.547,
                    tgas=293.15, gas_flow_basis='STANDARD', p_standard=101300.,
                    t_standard=293.15, rg=287., rho_l0=997.05,
                    pipe_regime_mode='FIXED_STRATIFIED', closure_riser=3.,
                    gas_eos='BAROTROPIC', dissipation_mode='LAGGED',
                    junction_model='MATCHED', flux_scheme='ROE')
    for key, value in required.items():
        if values.get(key) != value:
            raise ValueError(f'Nominal input differs for {key}')
    for key, value in (('mul', values['rho_l0']*values['jl']),
                       ('mug', values['p_standard']/(values['rg']*values['t_standard'])*values['jg'])):
        if not math.isclose(values[key], value, rel_tol=2e-15, abs_tol=0):
            raise ValueError(f'Nominal SI mass flux inconsistent: {key}')


def offsets(values):
    validate_nominal(values)
    variants=[]
    for suffix, axis, value in VARIANTS:
        change={axis:value}
        if axis == 'jl':
            change['mul']=values['rho_l0']*value
        elif axis == 'jg':
            change['mug']=values['p_standard']/(values['rg']*values['t_standard'])*value
        variants.append((suffix,axis,change))
    return variants


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--harness',type=Path,required=True)
    parser.add_argument('--reference-plan',type=Path,required=True)
    parser.add_argument('--workspace',type=Path,required=True)
    args=parser.parse_args()
    import hashlib
    if hashlib.sha256(args.harness.read_bytes()).hexdigest()!=HARNESS_SHA256:
        raise ValueError('Immutable executor identity differs')
    spec=importlib.util.spec_from_file_location('sealed_executor',args.harness)
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    old,seal=h.verify_seal(args.reference_plan.resolve())
    if seal['git_head']!=FREEZE:raise ValueError('Reference operator freeze differs')
    root=Path(old['root']);workspace=args.workspace.resolve();target=workspace/'bc_static'
    if target.exists():raise ValueError('Preserve existing campaign; destination exists')
    if root==workspace or root in workspace.parents:raise ValueError('Results must be external')
    reference=next(job for job in old['jobs'] if job['case']=='D_OSC')
    manifest_path=Path(reference['output'])/'manifest.json'
    manifest=json.loads(manifest_path.read_text())
    if manifest['status']!='COMPLETED' or manifest['returncode']!=0:raise ValueError('Nominal incomplete')
    if manifest['binary_sha256']!=seal['binary_sha256'] or manifest['source_sha256']!=seal['source_sha256']:
        raise ValueError('Nominal source/build provenance differs')
    if not manifest['binary_unchanged'] or not manifest['sources_unchanged']:raise ValueError('Nominal freeze changed')
    template=Path(reference['deck']);text0=template.read_text()
    if h.sha(template)!=reference['deck_sha256']:raise ValueError('Nominal deck changed')
    before=h.parameters(text0)
    metadata_path=Path(reference['output'])/'run_metadata.txt'
    observed=h.parameters(metadata_path.read_text())
    def canonical(v):
        if isinstance(v,str) and v.lower() in ('.true.','t'):return True
        if isinstance(v,str) and v.lower() in ('.false.','f'):return False
        return v
    if any(canonical(value)!=canonical(observed.get(key)) for key,value in before.items()):
        raise ValueError('Nominal effective deck differs from run metadata')
    jobs=[];decks={}
    for suffix,axis,changes in offsets(before):
        name='D_OSC_bc_'+suffix;folder=target/'runs'/name;deck=target/'decks'/(name+'.nml')
        relative=os.path.relpath(folder,root);text=text0
        for key,value in changes.items():text=h.assign(text,key,format(value,'.17g'))
        text=h.assign(text,'out_dir',repr(relative))
        values=h.parameters(text);actual=h.parameter_changes(before,values)
        if actual!=sorted([*changes,'out_dir']):raise ValueError('Mixed axes: '+str(actual))
        for key,value in before.items():
            if key not in changes and key!='out_dir' and values[key]!=value:
                raise ValueError('Unrequested control changed: '+key)
        decks[deck]=text
        jobs.append(dict(name=name,case='D_OSC',pipe_regime_mode='FIXED_STRATIFIED',axis=axis,
                         static_offset=changes[axis]-before[axis],dependent_si_overrides=changes,
                         deck=str(deck),output=str(folder),output_relative_to_root=relative,
                         parameters=values,physics={k:values[k] for k in h.PHYSICS},
                         reference_parameter_changes=actual))
    for path,text in decks.items():
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('x') as stream:stream.write(text)
    for job in jobs:job['deck_sha256']=h.sha(job['deck'])
    inputs=dict(old['input_sha256']);inputs[str(template)]=h.sha(template)
    plan=dict(schema='PRIMS-F.static_bc_plan.v1',state='PREPARED_UNSEALED',created_utc=h.now(),
              root=str(root),workspace=str(workspace),label='bc_static',prepared_git_head=FREEZE,
              input_sha256=inputs,jobs=jobs,environment=h.ENV,duration_s=3600.,
              maximum_total_solvers=4,maximum_owned_workers=1,runner_sha256=h.sha(args.harness),
              preparer_sha256=h.sha(__file__),
              reference_plan=dict(path=str(args.reference_plan.resolve()),sha256=h.sha(args.reference_plan)),
              reference_manifest=dict(path=str(manifest_path),sha256=h.sha(manifest_path)),
              reference_deck=dict(path=str(template),sha256=h.sha(template)),
              reference_metadata=dict(path=str(metadata_path),sha256=h.sha(metadata_path)),
              source=dict(doi='10.1016/j.ijmultiphaseflow.2010.04.003',table=3,case=13,page=653,
                          deviation_definition='Calculated from experimental signals (p.652); statistic not specified'),
              criteria=dict(complete_intervals=10,accommodation_relative_limit=.01,
                            inventory_relative_limit=1e-6,processed_mass_relative_limit=1e-6),
              execution_gate='Parent authorized six sequential runs after verification/seal; workers=1 only.',
              extension_policy='If complete but unaccommodated, preserve and create a duration-only extension in a new cohort.',
              interpretation='Static one-axis sensitivity, not inferred experimental forcing, frequency, probability or calibration.')
    h.dump_new(target/'plan.json',plan)
    print(json.dumps(dict(plan=str(target/'plan.json'),jobs=len(jobs),solvers_started=0,
                          binary_sha256=seal['binary_sha256'])))


if __name__=='__main__':main()
