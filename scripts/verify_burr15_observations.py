#!/usr/bin/env python3
"""Recompute Burr Eq.15 from SI fields, without importing the solver.

The height inversion uses circular-cap area as a function of height. The
observed criterion is a loss of stratified-candidate stability, not a flow label.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

STATUS = {0: "OK", 1: "INVALID_STATE", 2: "PHASE_ABSENT", 3: "DIRECTION",
          4: "ANGLE", 5: "NOT_EVALUATED", 6: "GEOMETRY", 7: "NUMERIC"}


def cap_fraction(height):
    """Area fraction of a circular cap, 0 <= height/D <= 1/2.

    Integrate sqrt(t*(1-t)) by its binomial series in the small cap to
    avoid cancellation in acos(x)-x*sqrt(1-x*x). No alpha/gamma inverse.
    """
    if height == 0:
        return 0.0
    if height >= 0.05:
        x = 1 - 2*height
        return (math.acos(x) - x*2*math.sqrt(height*(1-height)))/math.pi
    total, coefficient, power = 2/3, 1.0, 1.0
    for n in range(1, 40):
        coefficient *= (n - 1.5)/n
        power *= height
        term = coefficient*power/(n + 1.5)
        previous = total
        total += term
        if total == previous:
            break
    return (8/math.pi)*height**1.5*total


def height_geometry(alpha):
    """Return liquid height/D, gas headroom/D, liquid wall angle in radians."""
    if not 0 < alpha < 1:
        raise ValueError("Biphasic alpha required")
    fraction = min(alpha, 1-alpha)
    low, high = 0.0, 0.5
    for _ in range(1100):
        middle = (low+high)/2
        if middle == low or middle == high:
            break
        if cap_fraction(middle) < fraction:
            low = middle
        else:
            high = middle
    cap = (low+high)/2
    h, room = (cap, 1-cap) if alpha >= 0.5 else (1-cap, cap)
    angle = math.atan2(2*math.sqrt(h*room), room-h)
    return h, room, angle


def reference(field, diameter, gravity, area=None, closure_id=1):
    """Declared implementation domain; this does not qualify the source domain."""
    if closure_id != 1 or int(float(field.get("state_status", 0))) != 0:
        return {"status": 5}
    a, ul, ug, rl, rg, theta = (float(field[k]) for k in
        ("alpha", "ul_m_s", "ug_m_s", "rhol_kg_m3", "rhog_kg_m3", "theta_rad"))
    area = (math.pi/4)*diameter*diameter if area is None else area
    if not all(math.isfinite(x) for x in (a, ul, ug, rl, rg, theta, diameter, gravity, area)):
        return {"status": 1}
    if not 0 <= a <= 1 or not rl > rg > 0 or gravity <= 0:
        return {"status": 1}
    if diameter <= 0 or area <= 0 or not math.isclose(area, (math.pi/4)*diameter*diameter, rel_tol=32*math.ulp(1.0), abs_tol=0):
        return {"status": 6}
    if a in (0, 1):
        return {"status": 2}
    if not -math.pi/2 < theta <= 0:
        return {"status": 4}
    if ug < 0 or ul < 0:
        return {"status": 3}
    try:
        h, room, angle = height_geometry(a)
        sl = diameter*angle
        critical = room**2*((rl-rg)/rl)*gravity*math.cos(theta)*area*a/sl
        ug2 = ug*ug
        margin = ug2-critical
        if not all(math.isfinite(x) for x in (h, critical, ug2, margin)) or critical <= 0:
            return {"status": 7}
    except (ArithmeticError, ValueError):
        return {"status": 7}
    return dict(status=0, hL_over_D=h, ug2_m2_s2=ug2,
                ug2_critical_m2_s2=critical, margin_m2_s2=margin,
                signed_margin_m_s=ug-math.sqrt(critical), ug_m_s=ug,
                critical_ug_m_s=math.sqrt(critical), unstable=margin >= 0)


def frames(path, counts):
    """Stream one timestamp, collapsing only exactly identical final duplicates."""
    current, rows = None, {}
    with Path(path).open(newline="", encoding="utf-8") as source:
        for row in csv.DictReader(source):
            t, s = float(row["time_s"]), float(row["s_m"])
            if not math.isfinite(t) or not math.isfinite(s):
                raise ValueError(f"Nonfinite coordinates in {path}")
            if current is not None and t < current:
                raise ValueError(f"Unordered timestamps in {path}")
            if current is not None and t != current:
                yield current, rows
                rows = {}
            current = t
            if s in rows:
                if rows[s] != row:
                    raise ValueError(f"Conflicting duplicate at {(t,s)} in {path}")
                counts["identical_duplicates"] += 1
            else:
                rows[s] = row
        if current is not None:
            yield current, rows


class Intervals:
    """Observed intervals; event times remain bracketed by the sampling times."""
    def __init__(self):
        self.previous = None
        self.active = None
        self.intervals = []

    def update(self, t, value):
        if value is True:
            if self.active is None:
                self.active = dict(first_observed_s=t, last_observed_s=t,
                    onset_bracket_s=[self.previous, t], left_censored=self.previous is None,
                    samples=0)
            self.active["last_observed_s"] = t
            self.active["samples"] += 1
        elif self.active is not None:
            self.active.update(end_bracket_s=[self.active["last_observed_s"], t if value is False else None],
                               right_censored=value is None)
            self.intervals.append(self.active)
            self.active = None
        self.previous = t if value is not None else None

    def finish(self):
        self.update(None, None)
        return self.intervals


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_cycle_edges(path, run_root=None):
    if path is None:
        return []
    data = json.loads(Path(path).read_text())
    if "runs" in data:
        if len(data["runs"]) != 1:
            raise ValueError("Cycle summary must contain exactly one run")
        run = data["runs"][0]
        if not run.get("cycles_available", False):
            raise ValueError("Cycle summary does not provide detected cycles")
        if run_root is not None:
            run_root = Path(run_root)
            if run.get("run") != run_root.name:
                raise ValueError("Cycle summary belongs to another run")
            provenance = run.get("conservation", {}).get("provenance", {})
            for name in ("monitor.csv", "manifest.json", "run_metadata.txt", "run_status.txt"):
                if name not in provenance or not (run_root/name).is_file() or file_sha256(run_root/name) != provenance[name]:
                    raise ValueError(f"Cycle provenance mismatch or missing: {name}")
        edges = run["cycle_edges_s"]
    else:
        if run_root is not None:
            raise ValueError("CLI cycle summary requires run provenance")
        edges = data["cycle_edges_s"]
    edges = [float(x) for x in edges]
    if len(edges) < 2 or not all(math.isfinite(x) for x in edges) or any(b <= a for a, b in zip(edges, edges[1:])):
        raise ValueError("Cycle edges must be finite and strictly increasing")
    return edges


def verify(fields, burr, diameter, gravity=9.81, area=None, expected_cells=None,
           start=None, end=None, max_gap=None, cycle_edges=(), rtol=2e-10, atol=2e-12):
    paths = [Path(fields), Path(burr)]
    snapshots = [(p.stat().st_size, p.stat().st_mtime_ns) for p in paths]
    counts, mismatches, domain, errors = Counter(), Counter(), Counter(), Counter()
    examples, per_cell, grid, previous_t = [], {}, None, None
    reference_bounds = {}
    global_loss = Intervals()
    field_iter = iter(frames(fields, counts))
    frame = next(field_iter, None)
    first_loss = None
    cycle_counts = [Counter() for _ in range(max(0, len(cycle_edges)-1))]
    def missing_burr_frame(t):
        if (start is not None and t < start) or (end is not None and t > end):
            return
        counts["frames"] += 1
        counts["incomplete_frames"] += 1
        counts["field_only_frames"] += 1
        mismatches["coverage"] += 1
        global_loss.update(t, None)
        for track in per_cell.values():
            track.update(t, None)
    for t, observations in frames(burr, counts):
        while frame is not None and frame[0] < t:
            missing_burr_frame(frame[0])
            frame = next(field_iter, None)
        if start is not None and t < start:
            continue
        if end is not None and t > end:
            break
        counts["frames"] += 1
        matched = frame[1] if frame is not None and frame[0] == t else {}
        if grid is None and (expected_cells is None or len(observations) == expected_cells):
            grid = set(observations)
        complete = grid is not None and set(observations) == grid and set(observations) <= set(matched)
        if expected_cells is not None:
            complete = complete and len(observations) == expected_cells
        if max_gap is not None and previous_t is not None and t-previous_t > max_gap*(1+1e-12):
            counts["time_gaps"] += 1
            global_loss.update(t, None)
            for track in per_cell.values():
                track.update(t, None)
        previous_t = t
        counts["complete_frames" if complete else "incomplete_frames"] += 1
        if not complete:
            mismatches["coverage"] += 1
        flags = []
        for s in sorted((grid or set(observations)) | set(observations)):
            track = per_cell.setdefault(s, Intervals())
            if s not in observations or s not in matched:
                counts["unmatched_samples"] += 1
                track.update(t, None)
                continue
            row, field = observations[s], matched[s]
            counts["samples"] += 1
            ref = reference(field, diameter, gravity, area, int(float(row["closure_id"])))
            domain[STATUS[ref["status"]]] += 1
            observed_status = int(float(row["status"]))
            local_errors = []
            if observed_status != ref["status"]:
                local_errors.append("status")
            if ref["status"] == 0:
                counts["evaluable_samples"] += 1
                for quantity in ("hL_over_D", "ug_m_s", "critical_ug_m_s", "ug2_m2_s2", "ug2_critical_m2_s2", "margin_m2_s2", "signed_margin_m_s"):
                    value = ref[quantity]
                    bounds = reference_bounds.setdefault(quantity, {"min": value, "max": value})
                    bounds["min"] = min(bounds["min"], value)
                    bounds["max"] = max(bounds["max"], value)
                for key in ("hL_over_D", "ug2_m2_s2", "ug2_critical_m2_s2", "margin_m2_s2"):
                    actual = float(row[key])
                    if math.isfinite(actual):
                        errors[key] = max(errors[key], abs(actual-ref[key]))
                    else:
                        counts["nonfinite_recorded_values"] += 1
                    if not math.isclose(actual, ref[key], rel_tol=rtol, abs_tol=atol):
                        local_errors.append(key)
                if bool(int(float(row["unstable"]))) != ref["unstable"]:
                    local_errors.append("unstable")
                if abs(ref["margin_m2_s2"]) <= atol+rtol*max(ref["ug2_m2_s2"], ref["ug2_critical_m2_s2"]):
                    counts["near_threshold_samples"] += 1
                flag = ref["unstable"]
                flags.append(flag)
                track.update(t, flag if complete else None)
                if flag:
                    counts["loss_samples"] += 1
                    if first_loss is None:
                        first_loss = dict(time_s=t, s_m=s, **ref)
            else:
                if bool(int(float(row["unstable"]))):
                    local_errors.append("unstable_outside_domain")
                track.update(t, None)
                flags.append(None)
            mismatches.update(local_errors)
            if local_errors and len(examples) < 20:
                examples.append(dict(time_s=t, s_m=s, mismatches=local_errors,
                                     recorded_status=observed_status, reference=ref))
        any_loss = any(value is True for value in flags)
        all_known = complete and bool(flags) and all(value is not None for value in flags)
        global_loss.update(t, any_loss if all_known else None)
        counts["loss_frames"] += int(any_loss)
        for i, (left, right) in enumerate(zip(cycle_edges, cycle_edges[1:])):
            if left <= t < right:
                cycle_counts[i]["sampled_frames"] += 1
                cycle_counts[i]["complete_evaluable_frames"] += int(all_known)
                cycle_counts[i]["frames_with_observed_loss"] += int(any_loss)
                break
        if frame is not None and frame[0] == t:
            frame = next(field_iter, None)
    if counts["samples"] == 0:
        raise ValueError("No matched observations in the requested interval")
    # Extra field frames cannot silently disappear from coverage or tracking.
    if frame is not None:
        missing_burr_frame(frame[0])
    for t, _ in field_iter:
        if end is not None and t > end:
            break
        missing_burr_frame(t)
    hashes = [file_sha256(p) for p in paths]
    if snapshots != [(p.stat().st_size, p.stat().st_mtime_ns) for p in paths]:
        raise ValueError("Input changed during analysis; use a completed run or frozen snapshot")
    return dict(schema="PRIMS-F.burr15_observations.v1", consistent=not any(mismatches.values()),
        inputs=[dict(path=str(p.resolve()), size_bytes=size, mtime_ns=mtime, sha256=digest)
                for p, (size,mtime), digest in zip(paths,snapshots,hashes)],
        verifier_sha256=file_sha256(__file__),
        run_provenance={name:file_sha256(paths[0].parent/name) for name in
                        ("manifest.json", "run_metadata.txt", "run_status.txt") if (paths[0].parent/name).is_file()},
        formula="Burr Eq15/p3: rhol*Sl denominator; Sl is liquid-wetted wall",
        method="Independent circular-cap height inversion; SI fields; signed ug and ul domain guards",
        domain_contract="Implementation contract theta in (-pi/2,0], ul>=0, ug>=0; not physical-domain qualification",
        coverage_limit="Expected cell count supplied" if expected_cells else "First frame reference only; initially truncated grid cannot be excluded",
        sampling_limit="Intervals contain sampled losses; onset/end bracketed, unobserved crossings not excluded; no slugging classification",
        parameters=dict(diameter_m=diameter, gravity_m_s2=gravity, area_m2=area, rtol=rtol, atol=atol,
                        start_s=start, end_s=end, max_gap_s=max_gap, expected_pipe_cells=expected_cells),
        counts=dict(counts), domain_counts=dict(domain), mismatches=dict(mismatches),
        maximum_absolute_errors=dict(errors), reference_bounds=reference_bounds,
        bounds_universe="Matched samples with independent status OK only; invalid/absent observations excluded",
        mismatch_examples=examples, first_observed_loss=first_loss,
        observed_loss_intervals=global_loss.finish(),
        cell_loss_intervals=[dict(s_m=s, intervals=track.finish()) for s,track in sorted(per_cell.items())],
        supplied_cycles=[dict(index=i+1, start_s=cycle_edges[i], end_s=cycle_edges[i+1], **dict(c)) for i,c in enumerate(cycle_counts)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fields", type=Path, required=True)
    parser.add_argument("--burr15", type=Path, required=True)
    parser.add_argument("--diameter", type=float, required=True)
    parser.add_argument("--gravity", type=float, default=9.81)
    parser.add_argument("--area", type=float)
    parser.add_argument("--expected-pipe-cells", type=int)
    parser.add_argument("--start", type=float)
    parser.add_argument("--end", type=float)
    parser.add_argument("--max-gap", type=float)
    parser.add_argument("--cycle-summary", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists; choose a new versioned output")
    if args.fields.parent.resolve() != args.burr15.parent.resolve():
        parser.error("fields and burr15 must belong to the same run directory")
    result = verify(args.fields, args.burr15, args.diameter, args.gravity, args.area,
                    args.expected_pipe_cells, args.start, args.end, args.max_gap,
                    read_cycle_edges(args.cycle_summary, args.fields.parent))
    if args.cycle_summary:
        result["cycle_source"] = dict(path=str(args.cycle_summary.resolve()),
            sha256=file_sha256(args.cycle_summary),
            limitation="Supplied cycle edges label samples; waveform detection is not repeated here")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"consistent":result["consistent"], "counts":result["counts"], "mismatches":result["mismatches"]}))
    return 0 if result["consistent"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
