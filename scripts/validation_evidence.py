"""Named-column, standard-library helpers for numerical evidence reports.

This module does not fit a physical model or classify a flow regime. Time
statistics integrate piecewise-linear observations over the stated window;
missing observations reduce coverage instead of becoming zeros.
"""
from __future__ import annotations

import bisect
import csv
import fnmatch
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRESSURE_KEYS = ('p_junction_Pa', 'p_riser_base_Pa', 'p_pipe_end_Pa')


def finite(value):
    try:
        x = float(str(value).replace('D', 'E').replace('d', 'e'))
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def metadata(path):
    if not path.is_file():
        return {}
    text = path.read_text(errors='replace')
    out = {}
    for match in re.finditer(r'(\w+)\s*=\s*("[^"]*"|\x27[^\x27]*\x27|[^,\n]+)', text):
        key, value = match[1].lower(), match[2].strip()
        if value[:1] in ('"', "'"):
            out[key] = value.strip('"\' ').strip()
        elif value.upper() in ('T', '.TRUE.', 'F', '.FALSE.'):
            out[key] = value.upper() in ('T', '.TRUE.')
        else:
            out[key] = finite(value) if finite(value) is not None else value
    return out


def load_monitor(path):
    """Return sorted unique accepted times, named channels and data-quality flags."""
    if not path.is_file():
        return None, ['monitor_missing']
    warnings = []
    by_time = {}
    skipped = duplicates = 0
    previous = -math.inf
    reordered = False
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            return None, ['monitor_header_missing']
        if 'time_s' not in [key.strip() for key in reader.fieldnames]:
            return None, ['named_time_column_missing']
        for row in reader:
            if None in row:
                skipped += 1
                continue
            row = {key.strip(): finite(value) for key, value in row.items()}
            time = row.get('time_s')
            if time is None:
                skipped += 1
                continue
            reordered |= time < previous
            previous = time
            duplicates += time in by_time
            by_time[time] = row  # A final duplicate replaces the earlier observation.
    if skipped:
        warnings.append(f'malformed_or_nonfinite_time_rows_skipped:{skipped}')
    if duplicates:
        warnings.append(f'duplicate_times_last_row_retained:{duplicates}')
    if reordered:
        warnings.append('nonmonotone_input_times_sorted')
    if not by_time:
        return None, warnings + ['no_valid_monitor_rows']
    times = sorted(by_time)
    keys = sorted(set().union(*(row.keys() for row in by_time.values())))
    data = {key: [by_time[time].get(key) for time in times] for key in keys}
    data['time_s'] = times
    return data, warnings


def read_json(path, default=None):
    if not path.is_file():
        return default
    return json.loads(path.read_text())


def load_runs(results, pattern='final_*', metrics_path=None, run_dirs=()):
    """Metrics supply period evidence; raw monitors supply aligned time statistics."""
    metric_map = {}
    metrics_digest = None
    if metrics_path and Path(metrics_path).is_file():
        raw = Path(metrics_path).read_bytes()
        metrics_digest = hashlib.sha256(raw).hexdigest()
        payload = json.loads(raw)
        rows = payload if isinstance(payload, list) else payload.get('runs', [])
        metric_map = {str(row['run']): row for row in rows if isinstance(row, dict) and 'run' in row}
    folders = {p.name: p for p in Path(results).glob(pattern) if p.is_dir()}
    for name in metric_map:
        if fnmatch.fnmatch(name, pattern):
            folders.setdefault(name, Path(results) / name)
    for folder in run_dirs:
        folders[Path(folder).name] = Path(folder)
    records = []
    for name, folder in sorted(folders.items()):
        metric = metric_map.get(name, {})
        data, warnings = load_monitor(folder / 'monitor.csv')
        if data is None and not metric:
            continue
        pars = dict(metric.get('parameters', {}))
        pars.update(metadata(folder / 'run_metadata.txt'))
        manifest = read_json(folder / 'manifest.json', metric.get('manifest', {})) or {}
        controls = manifest.get('controls', {})
        for source, target in [('cells', 'nodes'), ('riser_cells', 'n_riser'), ('dt', 'dt_max_s'),
                               ('case', 'regime'), ('closure', 'closure_riser'), ('eos', 'gas_eos'),
                               ('junction', 'junction_model'), ('dissipation', 'dissipation_mode'),
                               ('flux', 'flux_scheme'), ('time', 'time_dim')]:
            if target not in pars and source in controls:
                pars[target] = controls[source]
        status = metadata(folder / 'run_status.txt')
        end = data['time_s'][-1] if data else finite(metric.get('reached_time_s'))
        requested = finite(status.get('requested_time_s', pars.get('time_dim', metric.get('requested_time_s'))))
        termination = status.get('status', manifest.get('status', metric.get('run_status', 'UNKNOWN')))
        reports_success = termination in (0, 0.0, '0', 'COMPLETED')
        completed = bool(reports_success and end is not None and requested is not None and
                         abs(end - requested) <= max(1e-6, abs(requested) * 1e-10))
        if not status and not manifest and metric.get('completed') is True:
            completed = bool(end is not None and requested is not None and
                             abs(end - requested) <= max(1e-6, abs(requested) * 1e-10))
        if not completed:
            warnings.append('run_not_completed')
        records.append({'run': name, 'case': str(pars.get('regime', metric.get('case', 'unknown'))),
                        'folder': str(folder.resolve()), 'parameters': pars, 'manifest': manifest,
                        'termination_status': termination, 'completed': completed,
                        'reached_time_s': end, 'requested_time_s': requested,
                        'data': data, 'upstream_metrics': metric, 'warnings': warnings,
                        'provenance': {'monitor_sha256': file_hash(folder / 'monitor.csv'),
                                       'metadata_sha256': file_hash(folder / 'run_metadata.txt'),
                                       'metrics_sha256': metrics_digest}})
    return records


def value_at(times, values, time):
    if not times or time < times[0] or time > times[-1]:
        return None
    i = bisect.bisect_left(times, time)
    if i < len(times) and times[i] == time:
        return values[i]
    if i == 0 or i == len(times) or values[i - 1] is None or values[i] is None:
        return None
    weight = (time - times[i - 1]) / (times[i] - times[i - 1])
    return (1 - weight) * values[i - 1] + weight * values[i]


def clipped_series(data, key, start, end):
    if not data or key not in data or end <= start:
        return [], []
    times, values = data['time_s'], data[key]
    left, right = max(start, times[0]), min(end, times[-1])
    if right <= left:
        return [], []
    inside = [i for i, time in enumerate(times) if left < time < right]
    return ([left] + [times[i] for i in inside] + [right],
            [value_at(times, values, left)] + [values[i] for i in inside] + [value_at(times, values, right)])


def time_statistics(data, key, start, end):
    """Exact integrals of the linear interpolant, with a centered second moment."""
    times, values = clipped_series(data, key, start, end)
    valid = [value for value in values if value is not None]
    if len(times) < 2 or not valid:
        return None
    origin = valid[0]
    pieces, first, second = [], [], []
    for i in range(len(times) - 1):
        if values[i] is None or values[i + 1] is None:
            continue
        dt = times[i + 1] - times[i]
        a, b = values[i] - origin, values[i + 1] - origin
        pieces.append(dt)
        first.append(dt * (a + b) / 2)
        second.append(dt * (a * a + a * b + b * b) / 3)
    covered = math.fsum(pieces)
    if covered <= 0:
        return None
    centered_mean = math.fsum(first) / covered
    variance = max(0.0, math.fsum(second) / covered - centered_mean ** 2)
    return {'mean': origin + centered_mean, 'std': math.sqrt(variance),
            'minimum': min(valid), 'maximum': max(valid), 'peak_to_peak': max(valid) - min(valid),
            'initial': values[0], 'final': values[-1],
            'change': values[-1] - values[0] if values[0] is not None and values[-1] is not None else None,
            'sample_count': len(valid), 'covered_duration_s': covered,
            'coverage_fraction': covered / (end - start),
            'method': 'time-weighted piecewise-linear observations; missing segments excluded'}


def pressure_key(data, start, end):
    for key in PRESSURE_KEYS:
        stats = time_statistics(data, key, start, end)
        if stats and stats['coverage_fraction'] >= 1 - 1e-8:
            return key
    return None


def median_interval(times):
    intervals = [b - a for a, b in zip(times, times[1:]) if b > a]
    return statistics.median(intervals) if intervals else None


def period_summary(markers, interval, method, amplitude, floor):
    periods = [b - a for a, b in zip(markers, markers[1:]) if b > a]
    median = statistics.median(periods) if len(periods) >= 2 else None
    mean = statistics.mean(periods) if periods else None
    cv = statistics.pstdev(periods) / mean if len(periods) >= 2 and mean else None
    return {'method': method, 'marker_times_s': markers, 'periods_s': periods,
            'complete_intervals_found': len(periods), 'median_period_s': median, 'period_cv': cv,
            'samples_per_period': median / interval if median and interval else None,
            'ten_intervals_observed': len(periods) >= 10,
            'pressure_amplitude_above_detection_floor': amplitude is not None and amplitude >= floor,
            'detection_floor_Pa': floor,
            'interpretation': 'Period candidates in a numerical pressure signal; no flow-regime validation.'}


def pressure_cycles(record, key, start, end, floor=100.0):
    """Prefer upstream peak markers, with an explicit threshold-crossing fallback."""
    data = record.get('data')
    stats = time_statistics(data, key, start, end) if key else None
    amplitude = stats['peak_to_peak'] if stats else None
    times, values = clipped_series(data, key, start, end) if key else ([], [])
    interval = median_interval(times)
    upstream = record.get('upstream_metrics', {})
    cycles = upstream.get('cycles', {})
    peak_times = cycles.get('peak_times_s')
    upstream_start = finite(upstream.get('window_start_s'))
    upstream_end = finite(upstream.get('window_end_s'))
    covers_window = (upstream_start is not None and upstream_end is not None and
                     upstream_start <= start + 1e-7 and upstream_end >= end - 1e-7)
    rms = stats['std'] if stats else finite(cycles.get('pressure_rms_Pa'))
    if rms is not None and rms < 1.0:
        result = period_summary([], interval, 'no resolved pressure variation (time std<1Pa)', amplitude, floor)
        result['pressure_std_floor_Pa'] = 1.0
        return result
    if peak_times is not None and covers_window and (not key or upstream.get('pressure_location') == key):
        markers = [float(t) for t in peak_times if finite(t) is not None and start <= float(t) <= end]
        if amplitude is None:
            amplitude = finite(cycles.get('pressure_peak_to_peak_Pa'))
        if amplitude is not None and amplitude < floor:
            markers = []
        return period_summary(markers, interval or finite(cycles.get('sample_interval_s')),
                              'upstream analyze_validation.py peak markers, restricted to comparison window',
                              amplitude, floor)
    if not stats or stats['coverage_fraction'] < 1 - 1e-8 or amplitude < floor:
        return period_summary([], interval, 'no resolved pressure excursion', amplitude, floor)
    low = stats['minimum'] + .2 * amplitude
    high = stats['minimum'] + .8 * amplitude
    armed, markers = False, []
    for i, value in enumerate(values):
        if value is None:
            armed = False
            continue
        if value <= low:
            armed = True
        if armed and value >= high:
            time = times[i]
            if i and values[i - 1] is not None and value > values[i - 1] and values[i - 1] < high:
                time = times[i - 1] + (times[i] - times[i - 1]) * (high - values[i - 1]) / (value - values[i - 1])
            markers.append(time)
            armed = False
    out = period_summary(markers, interval, 'upcrossing80% after observation below20% of window pressure range',
                         amplitude, floor)
    out['thresholds_Pa'] = {'low': low, 'high': high}
    return out


def analysis_window(record, start, end=None):
    reached = finite(record.get('reached_time_s'))
    if reached is None:
        return None
    finish = min(reached, end) if end is not None else reached
    data = record.get('data')
    begin = max(start, data['time_s'][0]) if data else start
    return (begin, finish) if finish > begin else None


def json_ready(record):
    return {key: value for key, value in record.items() if key not in ('data', 'upstream_metrics')}


def fmt(value, digits=6):
    if value is None:
        return '—'
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    if isinstance(value, (float, int)):
        return f'{value:.{digits}g}'
    return str(value).replace('|', '\\|').replace('\n', ' ')


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
