#!/usr/bin/env python3
"""Pre-registered settling/convergence criteria v2 (SETTLING_CRITERIA_v2.json).

Pure functions on piecewise-linear monitor series. The window is fixed by
eleven pressure markers; every channel is judged on cycle-mean drift, cycle
range variation and cycle-to-cycle waveform shape, normalized by
max(mean cycle range, dimensioned absolute floor). Near-zero channels are
never divided by a near-zero mean. Optimal phase alignment is computed only
as a diagnostic and never replaces the phase-fixed error.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import numpy as np

PROTOCOL = Path(__file__).with_name('SETTLING_CRITERIA_v2.json')
POINTS = 256
DRIFT, RANGE, SHAPE = .05, .01, .02
CONV_LEVEL, CONV_RANGE, CONV_SHAPE = .01, .01, .02
PRESSURE_FLOOR_PA, HOLDUP_FLOOR, FEED_FRACTION = 1.0, 1e-4, 1e-4
CHANNELS = {
    'p_riser_base_Pa': 'pressure', 'p_junction_Pa': 'pressure',
    'holdup_liquid_base': 'holdup', 'holdup_liquid_out': 'holdup',
    'mdot_l_base_kg_s': 'flux_l', 'mdot_l_out_kg_s': 'flux_l',
    'mdot_g_base_kg_s': 'flux_g', 'mdot_g_out_kg_s': 'flux_g',
    'Ml_kg': 'inventory_l', 'Ml_pipe_kg': 'inventory_l', 'Ml_riser_kg': 'inventory_l',
    'Mg_kg': 'inventory_g', 'Mg_pipe_kg': 'inventory_g', 'Mg_riser_kg': 'inventory_g',
}


def protocol_sha256() -> str:
    return hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()


def validate_markers(edges) -> np.ndarray:
    e = np.asarray(edges, dtype=float)
    if e.ndim != 1 or len(e) != 11 or not np.all(np.isfinite(e)) or np.any(np.diff(e) <= 0):
        raise ValueError('ten complete cycles require eleven distinct strictly increasing markers')
    return e


def window_mean(t, y, left, right) -> float:
    inside = (t > left) & (t < right)
    x = np.r_[left, t[inside], right]
    v = np.r_[np.interp(left, t, y), y[inside], np.interp(right, t, y)]
    return float(np.sum(np.diff(x) * (v[1:] + v[:-1]) / 2) / (right - left))


def cycle_waveforms(t, y, edges, points=POINTS) -> np.ndarray:
    return np.array([np.interp(l + (r - l) * np.arange(points) / points, t, y)
                     for l, r in zip(edges[:-1], edges[1:])])


def feed_rates(conservation: dict, edges) -> dict:
    """Mean accepted inlet mass flow of each phase over the window (kg/s)."""
    out = {}
    for phase in 'lg':
        key = f'M{phase}_in_integral_kg'
        if key not in conservation:
            out[phase] = None
            continue
        a = np.interp(edges[0], conservation['time_s'], conservation[key])
        b = np.interp(edges[-1], conservation['time_s'], conservation[key])
        out[phase] = float(abs(b - a) / (edges[-1] - edges[0]))
    return out


def floor_for(kind: str, feeds: dict, period: float):
    if kind == 'pressure':
        return PRESSURE_FLOOR_PA
    if kind == 'holdup':
        return HOLDUP_FLOOR
    phase = kind[-1]
    feed = feeds.get(phase)
    if feed is None:
        return None
    if kind.startswith('flux'):
        return FEED_FRACTION * feed
    return FEED_FRACTION * feed * period


def channel_metrics(t, y, edges, floor) -> dict:
    edges = validate_markers(edges)
    waves = cycle_waveforms(t, y, edges)
    means = np.array([window_mean(t, y, l, r) for l, r in zip(edges[:-1], edges[1:])])
    ranges = waves.max(axis=1) - waves.min(axis=1)
    mean_wave = waves.mean(axis=0)
    rbar, mbar = float(ranges.mean()), float(means.mean())
    scale = max(rbar, floor)
    k = np.arange(len(means), dtype=float)
    slope = float(np.polyfit(k, means, 1)[0])
    drift = max(float(np.ptp(means)), 9 * abs(slope))
    range_var = float(np.ptp(ranges))
    shape = float(np.max(np.sqrt(np.mean((waves - mean_wave) ** 2, axis=1))))
    if not (math.isfinite(scale) and scale > 0):
        # No oscillation and no dimensioned floor (e.g. zero feed): not evaluable,
        # never normalized by an invented denominator.
        return dict(evaluable=False, mean=mbar, mean_range=rbar, floor=floor, passed=False,
                    reason='zero range and zero absolute floor')
    ratios = dict(drift=drift / scale, range=range_var / scale, shape=shape / scale)
    passes = dict(drift=ratios['drift'] <= DRIFT, range=ratios['range'] <= RANGE,
                  shape=ratios['shape'] <= SHAPE)
    return dict(mean=mbar, mean_range=rbar, floor=floor, scale=scale,
                normalization='max(mean cycle range, absolute floor)',
                cycle_means=means.tolist(), cycle_ranges=ranges.tolist(),
                drift=drift, drift_slope_per_cycle=slope, range_variation=range_var,
                max_cycle_rms_shape_deviation=shape, ratios=ratios, passes=passes,
                passed=all(passes.values()), mean_waveform=mean_wave.tolist())


def pressure_c0(periods, amplitudes, means) -> dict:
    variation = {name: float(np.ptp(x) / abs(np.mean(x)))
                 for name, x in [('period', periods), ('amplitude', amplitudes), ('mean', means)]}
    return dict(variation=variation, passed=all(v <= .005 for v in variation.values()))


def settle(t, monitor: dict, edges, feeds: dict, c0: dict) -> dict:
    edges = validate_markers(edges)
    period = float(np.mean(np.diff(edges)))
    channels, missing = {}, []
    for name, kind in CHANNELS.items():
        floor = floor_for(kind, feeds, period)
        if name not in monitor or floor is None:
            missing.append(name)
            continue
        metrics = channel_metrics(t, np.asarray(monitor[name], float), edges, floor)
        if not metrics.get('evaluable', True):
            missing.append(name)
            continue
        channels[name] = metrics
    failed = sorted(f'{n}:{c}' for n, m in channels.items() for c, ok in m['passes'].items() if not ok)
    settled = bool(c0['passed'] and not missing and not failed)
    return dict(protocol_sha256=protocol_sha256(), pressure_C0=c0, channels=channels,
                missing_channels=missing, failed_criteria=failed, settled=settled,
                feeds_kg_s=feeds, mean_period_s=period)


def best_lag(a, b):
    a = np.asarray(a) - np.mean(a)
    b = np.asarray(b) - np.mean(b)
    n = len(a)
    scale = math.sqrt(float(np.dot(a, a) * np.dot(b, b)))
    if scale <= 0:
        return None, None
    corr = [float(np.dot(a, np.roll(b, -k))) for k in range(-n // 2, n // 2)]
    k = int(np.argmax(corr)) - n // 2
    return k / n, corr[k + n // 2] / scale


def compare(coarse: dict, fine: dict) -> dict:
    """Convergence of two settled analyses of the same operator/EOS/build."""
    out = dict(both_settled=bool(coarse['settled'] and fine['settled']), channels={})
    for name in CHANNELS:
        if name not in coarse['channels'] or name not in fine['channels']:
            out['channels'][name] = dict(available=False)
            continue
        c, f = coarse['channels'][name], fine['channels'][name]
        floor = max(c['floor'], f['floor'])
        rbar = max(c['mean_range'], f['mean_range'])
        level_scale = max(abs(c['mean']), abs(f['mean']), rbar, floor)
        shape_scale = max(rbar, floor)
        rms = float(np.sqrt(np.mean((np.asarray(f['mean_waveform']) - np.asarray(c['mean_waveform'])) ** 2)))
        lag, corr = best_lag(c['mean_waveform'], f['mean_waveform'])
        ratios = dict(level=abs(f['mean'] - c['mean']) / level_scale,
                      range=abs(f['mean_range'] - c['mean_range']) / shape_scale,
                      shape_phase_fixed=rms / shape_scale)
        passes = dict(level=ratios['level'] <= CONV_LEVEL, range=ratios['range'] <= CONV_RANGE,
                      shape_phase_fixed=ratios['shape_phase_fixed'] <= CONV_SHAPE)
        out['channels'][name] = dict(available=True, ratios=ratios, passes=passes, passed=all(passes.values()),
                                     phase_error_fraction_diagnostic=lag, lag_correlation_diagnostic=corr)
    failed = sorted(f'{n}:{c}' for n, m in out['channels'].items() if m.get('available')
                    for c, ok in m['passes'].items() if not ok)
    missing = sorted(n for n, m in out['channels'].items() if not m.get('available'))
    out.update(failed_criteria=failed, missing_channels=missing,
               converged_v2_channels=bool(out['both_settled'] and not failed and not missing))
    return out


if __name__ == '__main__':
    print(json.dumps(dict(protocol=str(PROTOCOL), sha256=protocol_sha256())))
