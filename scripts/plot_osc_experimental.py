#!/usr/bin/env python3
"""Compare the accepted OSC run with the published case-13 pressure trace."""
import argparse
import csv
import hashlib
from itertools import groupby
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

from analyze_burr15_campaign import unique_centres
from analyze_osc_envelopes import last_ten_cycles, pressure_envelope
from validation_evidence import load_runs, time_statistics


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def digitized_pressure(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    columns = ("time_s", "p_base_Pa", "p_lower_Pa", "p_upper_Pa")
    array = np.array([[float(row[key]) for key in columns] for row in rows])
    if array.shape[0] < 2 or not np.isfinite(array).all():
        raise ValueError("Published pressure trace must contain finite coordinates and envelopes")
    time, pressure, lower, upper = array.T
    if np.any(np.diff(time) <= 0) or np.any(lower > pressure) or np.any(pressure > upper):
        raise ValueError("Published trace time order or digitization envelope is inconsistent")
    return time, pressure, lower, upper


def field_window(path, left, right):
    time, data = [], []
    positions = None
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream)
        for t, group in groupby(rows, key=lambda row: float(row["time_s"])):
            if not left <= t <= right:
                continue
            cells = unique_centres(group)
            position = np.array([float(cell["s_m"]) for cell in cells])
            if positions is None:
                positions = position
            if not np.array_equal(position, positions):
                raise ValueError("The field mesh changed within the requested window")
            time.append(t)
            data.append([[float(cell[key]) for key in ("p_Pa", "alpha", "ug_m_s")] for cell in cells])
    fields = np.array(data)
    if len(time) < 2 or not np.isfinite(fields).all():
        raise ValueError("Resolved space–time fields are required")
    return np.array(time), positions, fields


def overlay(experiment, model, target):
    te, pe, lower, upper = experiment
    tm, first_cell, junction = model
    fig, ax = plt.subplots(figsize=(11, 4.6), layout="constrained")
    ax.fill_between(te-te[0],lower/1e5,upper/1e5,color="#53789d",alpha=.25,label="Envelope da digitalização")
    ax.plot(te-te[0],pe/1e5,color="#235789",lw=1.35,label="Experimento, Fig. 18")
    ax.plot(tm-tm[0],junction/1e5,color="#bb4a18",lw=1.45,label="PRIMS-F, traço na junção")
    ax.plot(tm-tm[0],first_cell/1e5,color="#87523e",ls="--",lw=.9,label="PRIMS-F, primeira célula do riser")
    ax.set(xlabel="Tempo decorrido na janela (s)",ylabel="Pressão absoluta (bar)",
           title="Caso 13 · pressão de base com origens temporais independentes")
    ax.grid(alpha=.25)
    ax.legend(frameon=False,ncol=2)
    fig.savefig(target,dpi=180)
    plt.close(fig)


def spectrum(time, pressure):
    """Hann-windowed relative power, on the observed approximately uniform cadence."""
    samples = np.linspace(time[0],time[-1],len(time))
    uniform = np.interp(samples,time,pressure)
    frequencies = np.fft.rfftfreq(len(samples),samples[1]-samples[0])
    power = np.abs(np.fft.rfft((uniform-uniform.mean())*np.hanning(len(samples))))**2
    if power[1:].max() <= 0:
        raise ValueError("Pressure spectrum is unresolved")
    return frequencies,power/power[1:].max()


def spectra(experiment, model, reference_period, target):
    fig,ax = plt.subplots(figsize=(8.5,4.2),layout="constrained")
    summaries = {}
    for (time,pressure),name,color in zip(
            (experiment,model),("Experimento digitalizado","PRIMS-F na junção"),("#235789","#bb4a18")):
        frequencies,power = spectrum(time,pressure)
        band = frequencies <= .10
        ax.plot(frequencies[band],power[band],color=color,marker=".",lw=1.15,label=name)
        peaks = [idx for idx in find_peaks(power)[0] if 0 < frequencies[idx] <= .10]
        summaries[name] = [{"f_Hz":float(frequencies[idx]),"relative_power":float(power[idx])}
                           for idx in sorted(peaks,key=lambda idx:power[idx],reverse=True)[:3]]
    for period in reference_period:
        ax.axvline(1/period,color="#666666",ls=":",lw=.9)
    ax.set(xlabel="Frequência (Hz)",ylabel="Potência relativa ao máximo de cada série",
           title="Caso 13 · conteúdo espectral em janelas de 300 s")
    ax.grid(alpha=.25)
    ax.legend(frameon=False)
    fig.savefig(target,dpi=180)
    plt.close(fig)
    return summaries


def space_time(fields, pipe_length, target):
    tm, positions, states = fields
    pressure, alpha, ug = states.transpose(2, 0, 1)
    fig, axes = plt.subplots(3,1,figsize=(11,8.4),sharex=True,layout="constrained")
    for ax, values, label, cmap in zip(
            axes, (pressure/1e5,1-alpha,alpha*ug),
            ("Pressão absoluta (bar)","Fração líquida","Velocidade superficial do gás (m/s)"),
            ("viridis","cividis","magma")):
        mesh = ax.pcolormesh(positions,tm,values,shading="nearest",cmap=cmap,rasterized=True)
        ax.axvline(pipe_length,color="#e8e8e8",ls="--",lw=.9)
        ax.set(ylabel="Tempo (s)",title=label)
        fig.colorbar(mesh,ax=ax)
    for where in (positions < pipe_length,positions > pipe_length):
        holdup = 1-alpha[:, where]
        levels = [value for value in (.4,.5,.6) if np.min(holdup)<value<np.max(holdup)]
        if levels and np.count_nonzero(where)>1:
            axes[1].contour(positions[where],tm,holdup,levels=levels,colors="white",linewidths=.45)
    axes[-1].set_xlabel("Posição ao longo do tubo (m)")
    fig.savefig(target,dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run",type=Path,required=True,help="completed D_OSC reference calculation")
    parser.add_argument("--digitization",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    root = args.run.resolve()
    record, = load_runs(root.parent,pattern=root.name)
    if not record["completed"] or record["case"]!="D_OSC":
        raise SystemExit("Require a completed D_OSC run")
    data = record["data"]
    te, pe, lower, upper = digitized_pressure(args.digitization)
    time = np.asarray(data["time_s"])
    pcell = np.asarray(data["p_riser_base_Pa"],dtype=float)
    pjunction = np.asarray(data["p_junction_Pa"],dtype=float)
    left = float(time[-1]-(te[-1]-te[0]))
    use = time>=left
    if not np.isfinite(pjunction[use]).all():
        raise ValueError("Junction pressure unavailable in the comparison window")
    cycles = last_ten_cycles(time,pcell,pressure_envelope(time,pcell))
    if not cycles["available"] or not cycles["useful_above_declared_detection_thresholds"]:
        raise ValueError("Ten fully enclosed pressure intervals are required")

    out = args.output.resolve()
    out.mkdir(parents=True)
    overlay((te,pe,lower,upper),(time[use],pcell[use],pjunction[use]),out/"osc_pressao_fig18.png")
    spatial = field_window(root/"fields.csv",cycles["window_start_s"],cycles["window_end_s"])
    space_time(spatial,record["parameters"]["lp"],out/"osc_campos_espaco_tempo.png")
    reference_period = (19.3,33.8)
    frequency_peaks = spectra((te,pe),(time[use],pjunction[use]),reference_period,out/"osc_espectro_fig18.png")
    summary = {
        "experiment": "Baliño, Burr e Nemoto 2010, Fig. 18, caso 13, curva publicada de Pb",
        "table3_reported_period_s": reference_period,
        "table3_DeltaP_reported_max_Pa": 69000,
        "experiment_time_s": [float(te[0]),float(te[-1])],
        "digitized_pressure_sample_range_Pa": [float(pe.min()),float(pe.max())],
        "digitized_pressure_envelope_range_Pa": [float(lower.min()),float(upper.max())],
        "digitized_300s_excursion_interval_Pa": [
            max(0.0,float(lower.max()-upper.min())),float(upper.max()-lower.min())],
        "spectral_bins": frequency_peaks,
        "model_time_s": [float(time[use][0]),float(time[-1])],
        "model_junction_pressure": time_statistics(data,"p_junction_Pa",left,float(time[-1])),
        "model_first_riser_cell_pressure": time_statistics(data,"p_riser_base_Pa",left,float(time[-1])),
        "model_last_ten_first_riser_cell_cycles": cycles,
        "provenance": {str(path):sha(path) for path in (
            root/"manifest.json",root/"monitor.csv",root/"fields.csv",args.digitization,Path(__file__))},
        "qualification": "A raster-derived comparison with independent time origins and a fixed model separator. "
                         "No raw synchronized sensor records, experimental error estimate or regime validation.",
    }
    (out/"comparison.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
