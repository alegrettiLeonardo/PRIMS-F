#!/usr/bin/env python3
"""ISOTHERMAL OSC: common-window comparison and continuous-variable events.

Pre-registered definitions (fixed before reading the refined runs):
  * Common window  W = [0, T_c], T_c = min over runs of the last accepted time
    (a run still in progress contributes its current time and is marked RUNNING).
  * Oscillation    p_in about its mean on [20, 60] s; half-cycle amplitudes
    A_k = (max - min)/2 between successive up-crossings; period from up-crossings;
    growth rate sigma = slope of ln A_k vs cycle mid-time over the complete
    cycles inside W that start in (60, 150) s. (Amended after N401/dt0.015625
    shortened W: the first rule, "start after 60 s", fitted 7 cycles in one run and
    6 in the others; the amended rule fits the same four cycle indices everywhere.)
  * Events from continuous monitor quantities (never from the solver status):
      E_g(f)  first t > 20 s with mdot_g_base <= f * mean(mdot_g_in on [20,60]),
              f in {0.5, 0.25, 0.10, 0.01}, by linear interpolation;
      E_h(h)  first t > 20 s with holdup_liquid_base >= h, h in {0.80, 0.90, 0.98}.
  * Newton failure time  t_fail = reached_time_s of run_status (termination_reason).
  * Amplitude-aligned clock: t_A(L) = time at which the fitted envelope
    ln A = ln A_0 + sigma t (cycles after 60 s inside W) reaches L = 3000 Pa;
    an event is discretization-consistent if E - t_A(L) converges even when
    E itself does not (seed amplitude A_0 depends on the initial-state error).
  * Amplitude-space event (seed-independent): for every complete cycle of the
    whole record, (A_k, min mdot_g_base / mean mdot_g_in, max holdup_base);
    A*(f) = amplitude at which the per-cycle gas trough reaches f, by linear
    interpolation of the trough ratio in ln A between successive complete cycles.
  * Waveform distance on W: max and RMS of |p_in^a - p_in^b| / mean(p_in) on a
    0.1 s grid (linear interpolation of the accepted monitor records).
Usage: iso_event_analysis.py RESULTS_ISO_DIR OUT_JSON OUT_MD
"""
import json, os, sys
import numpy as np

RUNS = {  # name: (N_pipe_cells, dt_max_s, mass/volume flux convention)
    "B_ISO_MASS": (51, 0.125, "MASS"),
    "C_ISO_VOLUME": (51, 0.125, "VOLUME"),
    "B_ISO_MASS_dt0.0625": (51, 0.0625, "MASS"),
    "B_ISO_MASS_dt0.03125": (51, 0.03125, "MASS"),
    "B_ISO_MASS_N101": (101, 0.125, "MASS"),
    "B_ISO_MASS_N201": (201, 0.125, "MASS"),
    "B_ISO_MASS_N101_dt0.0625": (101, 0.0625, "MASS"),
    "B_ISO_MASS_N401": (401, 0.125, "MASS"),
    "B_ISO_MASS_dt0.015625": (51, 0.015625, "MASS"),
}
FG = [0.5, 0.25, 0.10, 0.01]
HH = [0.80, 0.90, 0.98]


def status(d):
    f = os.path.join(d, "run_status.txt")
    out = {"termination_reason": "RUNNING"}
    if os.path.exists(f):
        for line in open(f):
            if "=" in line:
                k, v = line.strip().split("=", 1)
                out[k] = v.strip()
    return out


def first_cross(t, y, level, above, t0=20.0):
    m = t > t0
    t, y = t[m], y[m]
    ok = (y >= level) if above else (y <= level)
    idx = np.nonzero(ok)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return float(t[0])
    y0, y1 = y[i - 1], y[i]
    return float(t[i - 1] + (level - y0) * (t[i] - t[i - 1]) / (y1 - y0))


def cycles(t, p, tc):
    m0 = (t >= 20) & (t <= 60)
    x = p - p[m0].mean()
    up = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    tu = t[up] - x[up] * (t[up + 1] - t[up]) / (x[up + 1] - x[up])
    cyc = []
    for a, b, ia, ib in zip(tu[:-1], tu[1:], up[:-1], up[1:]):
        if b > tc:
            break
        seg = p[ia:ib + 2]
        cyc.append({"t_start": float(a), "t_end": float(b), "period_s": float(b - a),
                    "amplitude_Pa": float(0.5 * (seg.max() - seg.min()))})
    fit = [c for c in cyc if 60 < c["t_start"] < 150]  # same four cycle indices in every run
    sigma = None
    if len(fit) >= 3:
        tm = np.array([0.5 * (c["t_start"] + c["t_end"]) for c in fit])
        la = np.log([c["amplitude_Pa"] for c in fit])
        sigma, lna0 = np.polyfit(tm, la, 1)
        return cyc, float(sigma), float(lna0)
    return cyc, None, None


def trough_map(t, p, g, h, gin):
    m0 = (t >= 20) & (t <= 60)
    x = p - p[m0].mean()
    up = list(np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0])
    rows = []
    for ia, ib in zip(up[:-1], up[1:]):
        s = slice(ia, ib + 1)
        rows.append({"t_start": float(t[ia]), "amplitude_Pa": float(0.5 * (p[s].max() - p[s].min())),
                     "gas_trough_ratio": float(g[s].min() / gin), "holdup_base_max": float(h[s].max())})
    astar = {}
    for f in FG:
        v = None
        for r0, r1 in zip(rows[:-1], rows[1:]):
            if r0["gas_trough_ratio"] > f >= r1["gas_trough_ratio"]:
                w = (r0["gas_trough_ratio"] - f) / (r0["gas_trough_ratio"] - r1["gas_trough_ratio"])
                v = float(np.exp((1 - w) * np.log(r0["amplitude_Pa"]) + w * np.log(r1["amplitude_Pa"])))
                break
        astar["A_star_%g" % f] = v
    return rows, astar


def main():
    root, out_json, out_md = sys.argv[1:4]
    data, res = {}, {}
    for name, (n, dt, conv) in RUNS.items():
        d = os.path.join(root, name)
        f = os.path.join(d, "monitor.csv")
        if not os.path.exists(f):
            continue
        a = np.genfromtxt(f, delimiter=",", names=True)
        st = status(d)
        data[name] = a
        res[name] = {"N_pipe": n, "dt_max_s": dt, "flux": conv,
                     "termination_reason": st.get("termination_reason"),
                     "t_last_s": float(a["time_s"][-1]),
                     "t_fail_s": float(st["reached_time_s"]) if "reached_time_s" in st else None}
    tc = min(r["t_last_s"] for r in res.values())
    grid = np.arange(0.0, tc + 1e-12, 0.1)
    for name, a in data.items():
        t = a["time_s"]
        gin = a["mdot_g_in_kg_s"][(t >= 20) & (t <= 60)].mean()
        ev = {}
        for f in FG:
            ev["E_g_%g" % f] = first_cross(t, a["mdot_g_base_kg_s"], f * gin, above=False)
        for h in HH:
            ev["E_h_%g" % h] = first_cross(t, a["holdup_liquid_base"], h, above=True)
        cyc, sigma, lna0 = cycles(t, a["p_in_Pa"], tc)
        ta = (np.log(3000.0) - lna0) / sigma if sigma else None
        rows, astar = trough_map(t, a["p_in_Pa"], a["mdot_g_base_kg_s"], a["holdup_liquid_base"], gin)
        r = res[name]
        r.update({"mean_mdot_g_in_kg_s": float(gin), "events_s": ev, "cycles_in_window": cyc,
                  "growth_rate_1_s": sigma,
                  "growth_per_cycle": float(np.exp(sigma * np.mean([c["period_s"] for c in cyc[-3:]])))
                  if sigma else None,
                  "period_mean_s": float(np.mean([c["period_s"] for c in cyc[2:]])),
                  "seed_amplitude_first_cycle_Pa": cyc[0]["amplitude_Pa"],
                  "t_A3000_s": ta, "complete_cycles_full_record": rows, "amplitude_events_Pa": astar,
                  "aligned_events_s": {k: (v - ta if (v is not None and ta) else None) for k, v in ev.items()},
                  "min_mdot_g_base_in_window_kg_s": float(a["mdot_g_base_kg_s"][t <= tc].min())})
        if r["t_fail_s"] is not None:
            r["fail_minus_E_g_0.01_s"] = (r["t_fail_s"] - ev["E_g_0.01"]) if ev["E_g_0.01"] else None
    ref = {}
    for name, a in data.items():
        ref[name] = np.interp(grid, a["time_s"], a["p_in_Pa"])
    pairs = [("time", "B_ISO_MASS", "B_ISO_MASS_dt0.0625"), ("time", "B_ISO_MASS_dt0.0625", "B_ISO_MASS_dt0.03125"),
             ("time", "B_ISO_MASS_dt0.03125", "B_ISO_MASS_dt0.015625"),
             ("space", "B_ISO_MASS", "B_ISO_MASS_N101"), ("space", "B_ISO_MASS_N101", "B_ISO_MASS_N201"),
             ("space", "B_ISO_MASS_N201", "B_ISO_MASS_N401"),
             ("combined", "B_ISO_MASS", "B_ISO_MASS_N101_dt0.0625"),
             ("convention", "B_ISO_MASS", "C_ISO_VOLUME")]
    dist = []
    for axis, a, b in pairs:
        if a in ref and b in ref:
            e = np.abs(ref[a] - ref[b]) / ref[a].mean()
            dist.append({"axis": axis, "a": a, "b": b, "max_rel": float(e.max()), "rms_rel": float(np.sqrt((e ** 2).mean()))})
    out = {"common_window_s": [0.0, float(tc)], "definitions": __doc__, "runs": res, "waveform_distance_p_in": dist}
    json.dump(out, open(out_json, "w"), indent=1)
    with open(out_md, "w") as m:
        m.write("# ISOTHERMAL OSC: common window and continuous events\n\n")
        m.write("Common window [0, %.3f] s (shortest run). Events are interpolated first crossings after 20 s.\n\n" % tc)
        m.write("| run | N | dt | reason | t_fail s | E_g 0.25 | E_g 0.10 | E_g 0.01 | E_h 0.90 | E_h 0.98 | sigma 1/s | growth/cycle | cycles |\n")
        m.write("|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        fmt = lambda v: "-" if v is None else "%.3f" % v
        for name, r in res.items():
            e = r["events_s"]
            m.write("| %s | %d | %g | %s | %s | %s | %s | %s | %s | %s | %s | %s | %d |\n" % (
                name, r["N_pipe"], r["dt_max_s"], r["termination_reason"], fmt(r["t_fail_s"] or r["t_last_s"]),
                fmt(e["E_g_0.25"]), fmt(e["E_g_0.1"]), fmt(e["E_g_0.01"]), fmt(e["E_h_0.9"]), fmt(e["E_h_0.98"]),
                "-" if r["growth_rate_1_s"] is None else "%.5f" % r["growth_rate_1_s"],
                fmt(r["growth_per_cycle"]), len(r["cycles_in_window"])))
        m.write("\n## Converging observables and amplitude-aligned events\n\n| run | period s | seed A1 Pa | t(A=3000 Pa) | E_g0.25 - t_A | E_g0.10 - t_A | E_h0.90 - t_A |\n|---|---:|---:|---:|---:|---:|---:|\n")
        for name, r in res.items():
            e = r["aligned_events_s"]
            m.write("| %s | %.3f | %.0f | %s | %s | %s | %s |\n" % (name, r["period_mean_s"], r["seed_amplitude_first_cycle_Pa"],
                    fmt(r["t_A3000_s"]), fmt(e["E_g_0.25"]), fmt(e["E_g_0.1"]), fmt(e["E_h_0.9"])))
        m.write("\n## Amplitude-space events (complete cycles of the whole record)\n\n| run | A*(0.5) Pa | A*(0.25) Pa | A*(0.10) Pa | A*(0.01) Pa |\n|---|---:|---:|---:|---:|\n")
        fa = lambda v: "-" if v is None else "%.0f" % v
        for name, r in res.items():
            e = r["amplitude_events_Pa"]
            m.write("| %s | %s | %s | %s | %s |\n" % (name, fa(e["A_star_0.5"]), fa(e["A_star_0.25"]), fa(e["A_star_0.1"]), fa(e["A_star_0.01"])))
        m.write("\n## p_in waveform distance on the common window\n\n| axis | a | b | max rel | rms rel |\n|---|---|---|---:|---:|\n")
        for x in dist:
            m.write("| %s | %s | %s | %.4f | %.4f |\n" % (x["axis"], x["a"], x["b"], x["max_rel"], x["rms_rel"]))
        m.write("\n## Cycle amplitudes (Pa) inside the window\n\n")
        for name, r in res.items():
            m.write("- %s: %s\n" % (name, ", ".join("%.0f@%.1f" % (c["amplitude_Pa"], c["t_start"]) for c in r["cycles_in_window"])))
    print(open(out_md).read())


if __name__ == "__main__":
    main()
