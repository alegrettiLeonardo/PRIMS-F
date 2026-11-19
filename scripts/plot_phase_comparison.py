#!/usr/bin/env python3
"""Plot stored phase-aligned cycle comparisons without shifting channels."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(report,output):
    runs=[run for run in report["runs"] if "waveforms" in run]
    phase=np.arange(report["phase_grid_points"])/report["phase_grid_points"]
    fig,axes=plt.subplots(3,2,figsize=(12,9),sharex=True,layout="constrained")
    panels=[("p_junction_Pa","Pressão na junção (bar)",1e-5),
            ("holdup_liquid_base","Fração líquida · primeiro centro do riser",1),
            ("mdot_l_base_kg_s","Vazão líquida na base (kg/s)",1),
            ("mdot_g_base_kg_s","Vazão gasosa na base (kg/s)",1),
            ("Ml_riser_kg","Inventário líquido no riser (kg)",1),
            ("Mg_pipe_kg","Inventário gasoso no oleoduto (kg)",1)]
    for ax,(key,label,scale) in zip(axes.flat,panels):
        for run in runs:
            ax.plot(phase,np.array(run["waveforms"][key]["mean"])*scale,label=run["run"],lw=1.2)
        ax.set_ylabel(label)
        ax.grid(alpha=.25)
    for ax in axes[-1,:]:
        ax.set_xlabel("Fase do ciclo · pico de pressão em 0")
    axes[0,0].legend(fontsize=8,frameon=False)
    labels=[run["run"]+": "+run.get("status","NOT_VERIFIED") for run in report["runs"]]
    fig.suptitle(" | ".join(labels)+"\nSéries sem ciclos úteis não recebem curva de fase",fontsize=9)
    fig.savefig(output,dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    with args.output.open("xb"):
        pass
    plot(json.loads(args.report.read_text()),args.output)


if __name__=="__main__":
    main()
