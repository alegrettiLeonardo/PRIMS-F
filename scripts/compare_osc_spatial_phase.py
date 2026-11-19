#!/usr/bin/env python3
"""Compare OSC primitive fields on fixed physical positions and pressure phase.

Isovoid crossings are geometric level sets, not identified material fronts.
Cell-centred transport is kept separate from the accepted boundary ledger.
"""
import argparse
import csv
import hashlib
import json
from itertools import groupby
from pathlib import Path

import numpy as np

FIELDS = ('p_Pa', 'alpha', 'ul_m_s', 'ug_m_s', 'rhol_kg_m3', 'rhog_kg_m3')
CHANNELS = ('p_Pa', 'liquid_fraction', 'ul_m_s', 'ug_m_s', 'jl_m_s', 'jg_m_s',
            'mdot_l_local_kg_s', 'mdot_g_local_kg_s')
LEVELS = (.4, .5, .6)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def sample_positions(position, values, targets):
    """Linear interpolation inside a segment; never extrapolate or cross junction."""
    position = np.asarray(position)
    targets = np.asarray(targets)
    if (len(position) < 2 or np.any(np.diff(position) <= 0) or
            targets[0] < position[0] or targets[-1] > position[-1]):
        raise ValueError('Unordered coordinates or target outside cell-centre coverage')
    return np.stack([np.interp(targets, position, row) for row in np.asarray(values).T], axis=1)


def crossings(position, alpha, level):
    """Half-open crossing convention counts a node exactly on the level once."""
    result = []
    for i in range(len(position)-1):
        a, b = alpha[i:i+2]
        if a < level <= b or b < level <= a:
            x = position[i] + (position[i+1]-position[i])*(level-a)/(b-a)
            result.append((float(x), 'increasing' if b > a else 'decreasing'))
    return result


def read_window(path, bounds, pipe_length, area, fixed=None):
    """Keep one bracketing frame on either side of the complete-cycle window."""
    left, right = bounds[0], bounds[-1]
    previous = None
    retained = []
    seen_right = False
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        for time, rows in groupby(reader, key=lambda row: float(row['time_s'])):
            if time < left:
                previous = (time, list(rows))
                continue
            if not retained and previous is not None:
                retained.append(previous)
            retained.append((time, list(rows)))
            if time >= right:
                seen_right = True
                break
    if not retained or retained[0][0] > left or not seen_right:
        raise ValueError('Fields do not bracket the declared complete-cycle window')
    times = []
    frames = {'pipe': [], 'riser': []}
    crossing_rows = []
    native = None
    for time, rows in retained:
        by_position = {}
        for row in rows:
            s = float(row['s_m'])
            if s in by_position and row != by_position[s]:
                raise ValueError('Conflicting duplicate field record')
            by_position[s] = row
        position = np.array(sorted(by_position))
        values = np.array([[float(by_position[s][key]) for key in FIELDS] for s in position])
        if not np.isfinite(values).all():
            raise ValueError('Nonfinite field')
        if (np.any(values[:, 0] <= 0) or np.any(values[:, 1] < 0) or
                np.any(values[:, 1] > 1) or np.any(values[:, 4:] <= 0)):
            raise ValueError('Inadmissible pressure, fraction or density in field')
        if native is None:
            native = position
            if fixed is None:
                fixed = {name: position[mask] for name, mask in
                         (('pipe', position < pipe_length), ('riser', position > pipe_length))}
        elif not np.array_equal(native, position):
            raise ValueError('Mesh changed within the declared window')
        times.append(time)
        for name, mask in (('pipe', position < pipe_length), ('riser', position > pipe_length)):
            x, states = position[mask], values[mask]
            p, a, ul, ug, rl, rg = states.T
            jl, jg = (1-a)*ul, a*ug
            derived = np.stack((p, 1-a, ul, ug, jl, jg, area*rl*jl, area*rg*jg), axis=1)
            frames[name].append(sample_positions(x, derived, fixed[name]))
            if left <= time <= right:
                for level in LEVELS:
                    crossing_rows.extend((time, name, level, *item) for item in crossings(x, a, level))
    times = np.array(times)
    if np.any(np.diff(times) <= 0):
        raise ValueError('Fields time ordering is inconsistent')
    return times, fixed, {name: np.array(value) for name, value in frames.items()}, crossing_rows


def phase_average(times, frames, bounds, points=256):
    phase = np.arange(points)/points
    sampled = []
    for left, right in zip(bounds, bounds[1:]):
        targets = left + (right-left)*phase
        sampled.append(np.stack([np.interp(targets, times, signal)
                                 for signal in frames.reshape(len(times), -1).T], axis=1)
                       .reshape(points, *frames.shape[1:]))
    sampled = np.array(sampled)
    mean = sampled.mean(axis=0)
    scatter = np.sqrt(np.mean((sampled-mean)**2, axis=(0, 1)))
    return mean, scatter


def compare_fields(first, second):
    result = {}
    for segment in ('pipe', 'riser'):
        result[segment] = {}
        for i, key in enumerate(CHANNELS):
            a, b = first[segment][:, :, i], second[segment][:, :, i]
            rms = float(np.sqrt(np.mean((b-a)**2)))
            excursion = float(np.ptp(a))
            result[segment][key] = dict(rms_absolute=rms,
                maximum_absolute=float(np.max(abs(b-a))),
                rms_relative_to_reference_space_phase_range=rms/excursion if excursion else None,
                mean_difference=float((b-a).mean()))
    return result


def compare_isovoid(first, second, positions):
    """Ordered geometric crossings at each phase, with explicit topology mismatch."""
    result = {}
    for segment in ('pipe', 'riser'):
        result[segment] = {}
        for level in LEVELS:
            errors = []
            mismatch = 0
            empty = 0
            for a, b in zip(first[segment], second[segment]):
                left = crossings(positions[segment], 1-a[:, 1], level)
                right = crossings(positions[segment], 1-b[:, 1], level)
                if [item[1] for item in left] != [item[1] for item in right]:
                    mismatch += 1
                    continue
                if not left:
                    empty += 1
                errors.extend(last[0]-initial[0] for initial, last in zip(left, right))
            result[segment][str(level)] = dict(matched_crossings=len(errors),
                topology_mismatch_phases=mismatch, no_crossing_in_either_phases=empty,
                rms_position_difference_m=float(np.sqrt(np.mean(np.array(errors)**2))) if errors else None,
                maximum_position_difference_m=max(map(abs, errors)) if errors else None,
                interpretation='Ordered isovoid crossings of phase-mean fields; no material front identity')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase-report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('runs', type=Path, nargs='+')
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Preserve existing output; choose a new directory')
    report = json.loads(args.phase_report.read_text())
    summaries = {run['run']: run for run in report['runs']}
    args.output.mkdir(parents=True)
    index = dict(schema='PRIMS-F.spatial_phase_comparison.v1', phase_points=256,
        phase_marker='p_junction_Pa peak; same ten enclosed intervals as phase report',
        spatial_rule='Fixed physical reference cell centres; linear interpolation inside each segment; no extrapolation',
        front_rule='Piecewise-linear isovoid crossings, half-open levels; no material identity or inferred velocity',
        local_transport_rule='Cell-centred A*rho*j; accepted face ledger remains conservation authority',
        phase_report=dict(path=str(args.phase_report), sha256=sha(args.phase_report)),
        script_sha256=sha(__file__), runs=[], comparisons=[])
    reference = None
    fixed = None
    freeze = None
    for path in args.runs:
        summary = summaries[path.name]
        if summary['status'] not in ('SETTLED', 'UNSETTLED'):
            raise ValueError('Resolved cycles required')
        identity = (summary['freeze'], summary['binary_sha256'])
        if freeze is None:
            freeze = identity
        elif identity != freeze:
            raise ValueError('Operator changed between runs')
        if sha(path/'fields.csv') != summary['hashes']['fields.csv']:
            raise ValueError('Field hash differs from phase report')
        bounds = [row['peak_time_s'] for row in summary['cycles']['peak_records']]
        controls = summary['controls']
        area = np.pi*controls['d']**2/4
        times, fixed, frames, levels = read_window(path/'fields.csv', bounds, controls['lp'], area, fixed)
        averaged, scatter = {}, {}
        for name in frames:
            averaged[name], scatter[name] = phase_average(times, frames[name], bounds)
        npz = args.output/(path.name+'_phase_fields.npz')
        np.savez_compressed(npz, phase=np.arange(256)/256, **averaged,
            **{name+'_s_m': values for name, values in fixed.items()},
            **{name+'_scatter_rms': values for name, values in scatter.items()})
        csv_path = args.output/(path.name+'_isovoid_crossings.csv')
        with csv_path.open('x', newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(('time_s', 'segment', 'alpha_level', 's_m', 'orientation'))
            writer.writerows(levels)
        entry = dict(run=path.name, status=summary['status'], window_s=[bounds[0], bounds[-1]],
            frames=len(times), positions_per_segment={name: len(value) for name, value in fixed.items()},
            isovoid_counts={str(level): sum(row[2] == level for row in levels) for level in LEVELS},
            output_hashes={item.name: sha(item) for item in (npz, csv_path)},
            field_source_sha256=summary['hashes']['fields.csv'])
        index['runs'].append(entry)
        if reference is None:
            reference = (path.name, averaged)
        else:
            index['comparisons'].append(dict(reference=reference[0], run=path.name,
                channels=compare_fields(reference[1], averaged),
                geometric_isovoid=compare_isovoid(reference[1], averaged, fixed)))
    with (args.output/'index.json').open('x') as stream:
        json.dump(index, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
