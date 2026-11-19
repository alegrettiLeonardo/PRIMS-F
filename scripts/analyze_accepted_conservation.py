#!/usr/bin/env python3
"""Audit accepted-step mass ledgers; do not infer transport from sampled flux.

Only CSV arithmetic identities are tested here. The independent inventory and
processed-mass criteria are reported separately and never repaired.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import struct
import sys

PHASES = ('l', 'g')
REGIONS = ('', '_pipe', '_riser')
FACES = ('in', 'out', 'base_left', 'base_right')
EPS = sys.float_info.epsilon


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class Accumulator:
    def __init__(self):
        self.value = self.correction = 0.0

    def add(self, value):
        delta = value - self.correction
        result = self.value + delta
        self.correction = (result - self.value) - delta
        self.value = result
        return result


def audit_ledger(path):
    """Stream one initial row and consecutive accepted rows, retaining no history."""
    signed = {(p, f): Accumulator() for p in PHASES for f in FACES}
    absolute = {(p, f): Accumulator() for p in PHASES for f in FACES}
    residual = {(p, r): Accumulator() for p in PHASES for r in REGIONS}
    maximum = {(p, r): 0.0 for p in PHASES for r in REGIONS}
    mismatches = []
    identity_errors = {}
    junction_difference = dict.fromkeys(PHASES, 0.0)
    first = previous = None
    count = 0
    accepted_time_hash = hashlib.sha256()

    def identity(label, value, expected, scale):
        error = abs(value - expected)
        identity_errors[label] = max(identity_errors.get(label, 0.0), error)
        # This roundoff envelope verifies an identity; it is not the physical
        # balance criterion. The latter remains 1e-6 on named mass scales.
        tolerance = 256 * EPS * max(scale, abs(value), abs(expected), 1e-300)
        if error > tolerance and len(mismatches) < 20:
            mismatches.append(dict(row=count, identity=label, error=error, tolerance=tolerance))

    with Path(path).open(newline='') as stream:
        for raw in csv.DictReader(stream):
            if None in raw:
                raise ValueError('CSV row length differs from header')
            row = {key: float(value) for key, value in raw.items()}
            step = row['step']
            if step != count or not step.is_integer():
                raise ValueError('Expected initial step 0 and consecutive accepted steps')
            if row['time_s'] < 0 or row['dt_s'] < 0:
                raise ValueError('Negative time or timestep')
            for key, value in row.items():
                unavailable_initial_flux = count == 0 and row['flux_status'] != 0 and key.startswith('mdot_')
                if not math.isfinite(value) and not unavailable_initial_flux:
                    raise ValueError(f'Nonfinite {key} in accepted ledger row {count}')
            if count == 0:
                if row['time_s'] != 0 or row['dt_s'] != 0:
                    raise ValueError('Initial row must have time=dt=0')
                first = row
            else:
                if row['flux_status'] != 0 or row['dt_s'] <= 0 or row['time_s'] <= previous['time_s']:
                    raise ValueError('Accepted row must have valid fluxes and advancing time')
                identity('accepted_time', row['time_s'] - previous['time_s'], row['dt_s'], row['time_s'])
                accepted_time_hash.update(struct.pack('!qddd', count, row['time_s'], row['dt_s'], row['residual_norm']))
            for phase in PHASES:
                for face in FACES:
                    flow = row[f'mdot_{phase}_{face}_kg_s']
                    if count:
                        signed[phase, face].add(row['dt_s'] * flow)
                        absolute[phase, face].add(row['dt_s'] * abs(flow))
                    scale = absolute[phase, face].value
                    identity(f'{phase}_{face}_signed_integral', row[f'M{phase}_{face}_integral_kg'], signed[phase, face].value, scale)
                    identity(f'{phase}_{face}_absolute_integral', row[f'M{phase}_{face}_abs_integral_kg'], absolute[phase, face].value, scale)
                inlet, outlet, left, right = (signed[phase, f].value for f in FACES)
                junction_difference[phase] = max(junction_difference[phase], abs(left - right))
                for region, incoming, outgoing in (('', 'in', 'out'), ('_pipe', 'in', 'base_left'), ('_riser', 'base_right', 'out')):
                    mass = row[f'M{phase}{region}_kg']
                    if mass < 0:
                        raise ValueError('Negative phase inventory')
                    mass0 = first[f'M{phase}{region}_kg']
                    transport = signed[phase, outgoing].value - signed[phase, incoming].value
                    error = mass - mass0 + transport
                    processed = absolute[phase, incoming].value + absolute[phase, outgoing].value
                    scale = abs(mass) + abs(mass0) + processed
                    stored_error = row[f'balance_{phase}{region}_kg']
                    identity(f'{phase}{region}_balance', stored_error, error, scale)
                    maximum[phase, region] = max(maximum[phase, region], abs(stored_error))
                    defect = row[f'mass_residual_{phase}{region}_kg']
                    if count:
                        residual[phase, region].add(defect)
                        transport_step = row['dt_s'] * (row[f'mdot_{phase}_{outgoing}_kg_s'] - row[f'mdot_{phase}_{incoming}_kg_s'])
                        step_error = mass - previous[f'M{phase}{region}_kg'] + transport_step
                        identity(f'{phase}{region}_step_balance', row[f'step_balance_{phase}{region}_kg'], step_error, scale)
                        # The global residual also contains any noncanceling
                        # interface mass flux; retaining both sides exposes it.
                        if not region:
                            step_error += row['dt_s'] * (row[f'mdot_{phase}_base_left_kg_s'] - row[f'mdot_{phase}_base_right_kg_s'])
                        identity(f'{phase}{region}_accepted_residual', defect, step_error, scale)
                    else:
                        identity(f'{phase}{region}_initial_residual', defect, 0.0, scale)
                    cumulative = residual[phase, region].value
                    identity(f'{phase}{region}_residual_accumulator', row[f'mass_residual_sum_{phase}{region}_kg'], cumulative, scale)
                    identity(f'{phase}{region}_cumulative_defect', cumulative, error + (left - right if not region else 0.0), scale)
                identity(f'{phase}_inventory_partition', row[f'M{phase}_kg'], row[f'M{phase}_pipe_kg'] + row[f'M{phase}_riser_kg'], abs(row[f'M{phase}_kg']))
            previous = row
            count += 1
    if first is None:
        raise ValueError('Empty conservation ledger')
    balances = {}
    for phase in PHASES:
        for region, incoming, outgoing in (('', 'in', 'out'), ('_pipe', 'in', 'base_left'), ('_riser', 'base_right', 'out')):
            initial = abs(first[f'M{phase}{region}_kg'])
            processed = absolute[phase, incoming].value + absolute[phase, outgoing].value
            maximum_error = maximum[phase, region]
            balances[phase + (region or '_system')] = dict(
                max_abs_balance_kg=maximum_error,
                initial_inventory_kg=initial,
                absolute_processed_mass_kg=processed,
                relative_inventory=maximum_error / initial if initial else None,
                relative_processed=maximum_error / processed if processed else None,
                passes_inventory_1e6=maximum_error <= 1e-6 * initial if initial else None,
                passes_processed_1e6=maximum_error <= 1e-6 * processed if processed else None,
                final_residual_sum_kg=residual[phase, region].value,
            )
    return dict(accepted_steps=count - 1, final_time_s=previous['time_s'],
                accepted_time_dt_norm_sha256=accepted_time_hash.hexdigest(),
                arithmetic_contract_passes=not mismatches, mismatches=mismatches,
                max_identity_error_kg_or_s=identity_errors,
                max_junction_integral_mismatch_kg=junction_difference, balances=balances)


def audit_steps(path, ledger):
    counters = ('newton_iterations_total', 'jacobian_evaluations_total', 'residual_evaluations_total',
                'linear_solves_total', 'backtrack_trials_total', 'prepare_calls_total', 'rejected_steps_total')
    previous = None
    dt_min = math.inf
    dt_max = 0.0
    count = 0
    accepted_time_hash = hashlib.sha256()
    with Path(path).open(newline='') as stream:
        for raw in csv.DictReader(stream):
            row = {key: float(value) for key, value in raw.items()}
            count += 1
            if row['step'] != count or row['dt_s'] <= 0 or not all(math.isfinite(v) for v in row.values()):
                raise ValueError('Invalid accepted-step record')
            for key in counters:
                if row[key] < 0 or not row[key].is_integer() or (previous and row[key] < previous[key]):
                    raise ValueError(f'Invalid cumulative counter {key}')
            if previous and row['time_s'] <= previous['time_s']:
                raise ValueError('Step time did not advance')
            accepted_time_hash.update(struct.pack('!qddd', count, row['time_s'], row['dt_s'], row['residual_norm']))
            if row['cpu_since_previous_accept_s'] < 0 or row['cpu_total_s'] < 0 or row['wall_total_s'] < 0:
                raise ValueError('Negative cost')
            dt_min = min(dt_min, row['dt_s'])
            dt_max = max(dt_max, row['dt_s'])
            previous = row
    if count != ledger['accepted_steps'] or (previous and previous['time_s'] != ledger['final_time_s']):
        raise ValueError('Steps and conservation ledgers disagree')
    if accepted_time_hash.hexdigest() != ledger['accepted_time_dt_norm_sha256']:
        raise ValueError('Step/time/dt/residual histories disagree between ledgers')
    return dict(accepted_steps=count, min_dt_s=dt_min if count else None, max_dt_s=dt_max if count else None,
                last_accepted_counters={key: previous[key] for key in counters} if previous else {},
                cpu_to_last_accept_s=previous['cpu_total_s'] if previous else None,
                wall_to_last_accept_s=previous['wall_total_s'] if previous else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output exists; preserve earlier evidence and choose another path')
    root = args.run_dir.resolve()
    ledger = audit_ledger(root / 'conservation.csv')
    steps = audit_steps(root / 'steps.csv', ledger)
    status = dict(line.split('=', 1) for line in (root / 'run_status.txt').read_text().splitlines() if '=' in line)
    status = {k: float(v) for k, v in status.items()}
    if status['accepted_steps'] != ledger['accepted_steps'] or status['reached_time_s'] != ledger['final_time_s']:
        raise ValueError('Final solver status disagrees with accepted ledgers')
    for key, value in steps['last_accepted_counters'].items():
        status_key = 'rejected_steps' if key == 'rejected_steps_total' else key
        if status[status_key] < value or (status['status'] == 0 and status[status_key] != value):
            raise ValueError('Final solver counters disagree with accepted-step counters')
    artifacts = {}
    for name in ('conservation.csv', 'steps.csv', 'run_status.txt', 'run_metadata.txt', 'manifest.json'):
        path = root / name
        if path.is_file():
            artifacts[name] = dict(path=str(path), sha256=sha256(path), size_bytes=path.stat().st_size)
    result = dict(schema='PRIMS-F.accepted_conservation.v1', run_dir=str(root),
                  analyzer_sha256=sha256(__file__), artifacts=artifacts, status=status, ledger=ledger, steps=steps,
                  limits=['The balance criteria are separate from arithmetic-contract verification.',
                          'Mass residuals are measured nonlinear-solve defects, never balance corrections.',
                          'No regime classification or experimental validation is inferred.'])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(accepted_steps=ledger['accepted_steps'], final_time_s=ledger['final_time_s'],
                          arithmetic_contract_passes=ledger['arithmetic_contract_passes'], output=str(args.output))))
    return 0 if ledger['arithmetic_contract_passes'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
