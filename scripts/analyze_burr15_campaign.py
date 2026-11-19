#!/usr/bin/env python3
"""Analyze a frozen Burr 15 diagnostic campaign. Pressure cycles are OSC metrics."""
import argparse
import csv
from itertools import groupby
import json
import math
from pathlib import Path

import numpy as np

from analyze_osc_envelopes import last_ten_cycles, pressure_envelope
from analyze_phase_sequence import analyze as phase_sequence
from validation_evidence import load_runs, time_statistics, value_at

PRESSURE = "p_riser_base_Pa"
CHANNELS = (
    PRESSURE, "holdup_liquid_base", "holdup_liquid_out", "mdot_l_base_kg_s",
    "mdot_g_base_kg_s", "mdot_l_out_kg_s", "mdot_g_out_kg_s",
    "Ml_kg", "Mg_kg", "Ml_pipe_kg", "Mg_pipe_kg", "Ml_riser_kg", "Mg_riser_kg",
)


def mass_balance(data, cycle_bounds):
    times = data["time_s"]
    results = {}
    for phase in ("l", "g"):
        inventory = f"M{phase}_kg"
        incoming = f"M{phase}_in_integral_kg"
        outgoing = f"M{phase}_out_integral_kg"
        error = f"balance_{phase}_kg"
        peak = max(abs(v) for v in data[error] if v is not None)
        processed = max(abs(data[incoming][-1]), abs(data[outgoing][-1]))
        initial = data[inventory][0]
        cycles = []
        for left, right in zip(cycle_bounds, cycle_bounds[1:]):
            def increment(channel):
                return value_at(times, data[channel], right)-value_at(times, data[channel], left)
            dm, min_, mout = (increment(channel) for channel in (inventory, incoming, outgoing))
            closure = dm-min_+mout
            scale = max(abs(min_), abs(mout))
            cycles.append({"start_s": left, "end_s": right, "error_kg": closure,
                           "relative_to_net_throughput": abs(closure)/scale if scale else None})
        results[phase] = {
            "initial_kg": initial, "final_kg": data[inventory][-1],
            "maximum_abs_error_kg": peak, "final_error_kg": data[error][-1],
            "relative_to_initial_inventory": peak/initial if initial else None,
            "relative_to_net_throughput": peak/processed if processed else None,
            "cycles": cycles,
        }
    return results


def front_candidates(snapshot, pipe_length, levels):
    """Piecewise-linear void-fraction level crossings within each segment."""
    for left, right in zip(snapshot, snapshot[1:]):
        s0, a0 = float(left["s_m"]), float(left["alpha"])
        s1, a1 = float(right["s_m"]), float(right["alpha"])
        if s0 < pipe_length < s1:
            continue
        for level in levels:
            if a0 < level <= a1 or a1 < level <= a0:
                position = s0 + (s1-s0)*(level-a0)/(a1-a0)
                yield (level, position, "rising_alpha" if a1 > a0 else "falling_alpha")


def local_phase_transport(alpha, ul, ug, rhol, rhog, area):
    jl, jg = (1-alpha)*ul, alpha*ug
    return jl, jg, area*rhol*jl, area*rhog*jg


def unique_centres(rows):
    centres = {}
    for row in rows:
        position = float(row["s_m"])
        if position in centres and centres[position] != row:
            raise ValueError(f"Conflicting field records at t={row['time_s']}, s={position}")
        centres[position] = row
    return list(centres.values())


def scan_fields(folder, parameters, front_path, flux_path):
    pipe_length = parameters["lp"]
    area = math.pi * parameters["d"]**2 / 4
    pipe_id, riser_id = parameters.get("closure_pipe", 1), parameters["closure_riser"]
    alpha_min, alpha_max = 1.0, 0.0
    velocities = {"ul_m_s": 0, "ug_m_s": 0}
    levels = (.4, .5, .6)
    counts = {f"{value:.1f}": 0 for value in levels}
    samples = 0
    with ((folder/"fields.csv").open() as stream,
          front_path.open("w", newline="") as out,
          flux_path.open("w", newline="") as flux_stream):
        writer = csv.writer(out)
        writer.writerow(("time_s", "alpha_level", "s_m", "orientation"))
        flux = csv.writer(flux_stream)
        flux.writerow(("time_s", "s_m", "closure_id", "jl_m_s", "jg_m_s", "mdot_l_kg_s", "mdot_g_kg_s"))
        for time, group in groupby(csv.DictReader(stream), key=lambda row: float(row["time_s"])):
            # The terminal write may duplicate the last snapshot.
            snapshot = unique_centres(group)
            samples += len(snapshot)
            for row in snapshot:
                alpha = float(row["alpha"])
                if np.isfinite(alpha):
                    alpha_min, alpha_max = min(alpha_min, alpha), max(alpha_max, alpha)
                for key in velocities:
                    velocities[key] += float(row[key]) < 0
                position = float(row["s_m"])
                transport = local_phase_transport(alpha, float(row["ul_m_s"]), float(row["ug_m_s"]),
                                                  float(row["rhol_kg_m3"]), float(row["rhog_kg_m3"]), area)
                flux.writerow((time, position, int(pipe_id if position < pipe_length else riser_id), *transport))
            for level, position, orientation in front_candidates(snapshot, pipe_length, levels):
                counts[f"{level:.1f}"] += 1
                writer.writerow((time, level, position, orientation))
    return {"sample_count": samples, "alpha_min": alpha_min, "alpha_max": alpha_max,
            "negative_sample_counts": velocities, "crossings_by_level": counts,
            "interpretation": "spatial level crossings, without assigning front identity or regime"}


def summarize(record, output):
    folder = Path(record["folder"])
    out = {key: record[key] for key in ("run", "case", "completed", "reached_time_s", "termination_status", "provenance")}
    data = record["data"]
    if not data:
        return out
    stop = data["time_s"][-1]
    bounds = []
    if record["completed"]:
        out["late_window"] = {key: time_statistics(data, key, 900., stop) for key in CHANNELS}
        if record["case"].upper().endswith("OSC"):
            time, pressure = np.asarray(data["time_s"]), np.asarray(data[PRESSURE])
            cycles = last_ten_cycles(time, pressure, pressure_envelope(time, pressure))
            out["last_ten_pressure_cycles"] = cycles
            if cycles["available"] and cycles["useful_above_declared_detection_thresholds"]:
                bounds = [peak["peak_time_s"] for peak in cycles["peak_records"]]
                per_channel = {}
                for key in CHANNELS:
                    stats = [time_statistics(data, key, a, b) for a, b in zip(bounds, bounds[1:])]
                    values = {}
                    for metric in ("mean", "peak_to_peak"):
                        first = np.mean([cycle[metric] for cycle in stats[:5]])
                        last = np.mean([cycle[metric] for cycle in stats[5:]])
                        values[metric] = {"first_five": float(first), "last_five": float(last),
                                          "relative_change": float((last-first)/abs(first)) if first else None}
                    per_channel[key] = values
                out["last_ten_channel_comparison"] = per_channel
        sequence = phase_sequence(record, 900., stop)
        output.joinpath(f"{folder.name}_phase_sequence.json").write_text(json.dumps(sequence, allow_nan=False, indent=2)+"\n")
    steps = [dt for dt in data["dt_s"] if dt is not None and dt > 0]
    if steps:
        out["sampled_accepted_dt_s"] = {"count": len(steps), "min": min(steps),
                                       "median": float(np.median(steps)), "max": max(steps)}
    out["phase_balance"] = mass_balance(data, bounds)
    out["base_raw_gas_flux_min_kg_s"] = min((v for v in data["mdot_g_unconstrained_kg_s"] if v is not None), default=None)
    fields = folder/"fields.csv"
    if fields.is_file():
        parameters = dict(record["manifest"].get("parameters", {}))
        parameters.update(record["parameters"])
        out["fields"] = scan_fields(folder, parameters, output/f"{folder.name}_front_candidates.csv",
                                    output/f"{folder.name}_phase_transport.csv")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    records = load_runs(args.output, run_dirs=args.runs)
    result = [summarize(record, args.output) for record in records]
    (args.output/"campaign_metrics.json").write_text(json.dumps(result, allow_nan=False, indent=2)+"\n")


if __name__ == "__main__":
    main()
