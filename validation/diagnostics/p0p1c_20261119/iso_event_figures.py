#!/usr/bin/env python3
"""Figures for the ISOTHERMAL OSC refinement (reads iso_event_analysis JSON + raw runs).

Usage: iso_event_figures.py RESULTS_ISO_DIR ANALYSIS_JSON OUT_DIR
  iso_temporal.png  p_in, junction gas flux / feed and base holdup vs time;
  iso_growth.png    per-cycle amplitude (log) vs time and gas trough vs amplitude;
  iso_spatial.png   alpha(s) and u_g(s) snapshots of the coarsest and finest meshes.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    root, js, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    d = json.load(open(js))
    names = list(d["runs"])
    fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for n in names:
        a = np.genfromtxt(os.path.join(root, n, "monitor.csv"), delimiter=",", names=True)
        gin = d["runs"][n]["mean_mdot_g_in_kg_s"]
        ax[0].plot(a["time_s"], a["p_in_Pa"] / 1e3, lw=.7, label=n)
        ax[1].plot(a["time_s"], a["mdot_g_base_kg_s"] / gin, lw=.7)
        ax[2].plot(a["time_s"], a["holdup_liquid_base"], lw=.7)
    tc = d["common_window_s"][1]
    for x in ax:
        x.axvline(tc, color="k", ls=":", lw=.8)
    ax[0].set_ylabel("p_in (kPa)"); ax[1].set_ylabel("mdot_g junction / feed"); ax[2].set_ylabel("liquid holdup, riser base")
    ax[2].set_xlabel("t (s); dotted: end of common window"); ax[0].legend(fontsize=6, ncol=3)
    fig.tight_layout(); fig.savefig(os.path.join(out, "iso_temporal.png"), dpi=150); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    for n in names:
        r = d["runs"][n]["complete_cycles_full_record"]
        ax[0].semilogy([c["t_start"] for c in r], [c["amplitude_Pa"] for c in r], "o-", ms=3, lw=.8, label=n)
        ax[1].semilogx([c["amplitude_Pa"] for c in r], [c["gas_trough_ratio"] for c in r], "o-", ms=3, lw=.8)
    ax[0].set_xlabel("cycle start (s)"); ax[0].set_ylabel("p_in half-cycle amplitude (Pa)"); ax[0].legend(fontsize=6)
    ax[1].set_xlabel("p_in amplitude (Pa)"); ax[1].set_ylabel("min junction gas flux / feed in cycle")
    for f in (0.5, 0.25, 0.1):
        ax[1].axhline(f, color="k", ls=":", lw=.6)
    fig.tight_layout(); fig.savefig(os.path.join(out, "iso_growth.png"), dpi=150); plt.close(fig)
    fine = max((n for n in names if d["runs"][n]["dt_max_s"] == 0.125), key=lambda n: d["runs"][n]["N_pipe"])
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for j, n in enumerate(["B_ISO_MASS", fine]):
        f = np.genfromtxt(os.path.join(root, n, "fields.csv"), delimiter=",", names=True)
        times = np.unique(f["time_s"])
        for t in times[::max(1, len(times) // 6)]:
            m = f["time_s"] == t
            ax[0, j].plot(f["s_m"][m], f["alpha"][m], lw=.8, label="%.0f s" % t)
            ax[1, j].plot(f["s_m"][m], f["ug_m_s"][m], lw=.8)
        ax[0, j].set_title(n); ax[0, j].set_ylabel("alpha"); ax[1, j].set_ylabel("u_g (m/s)"); ax[1, j].set_xlabel("s (m)")
        ax[0, j].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(os.path.join(out, "iso_spatial.png"), dpi=150); plt.close(fig)


if __name__ == "__main__":
    main()
