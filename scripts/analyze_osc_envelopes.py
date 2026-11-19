#!/usr/bin/env python3
"""Analyze OSC pressure envelopes without declaring physical validation.

Reads immutable snapshots of the bytes currently present in monitor.csv.  A run
without a successful run_status.txt remains explicitly incomplete.  The main
observable is the slope of log(semi-amplitude) in the fixed 40--400 s window.
This is a signal fit, not a proven eigenvalue of the spatial discrete operator.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import linregress
from osc_run_selection import DEFAULT_PATTERNS, cohort_metadata, select_runs


WINDOWS = ((40., 200.), (200., 400.), (40., 400.), (900., 1400.),
           (1400., 1900.), (1900., 2400.), (2100., 2400.))
BLOCKS = ((900., 1400.), (1400., 1900.), (1900., 2400.))
PATTERNS = DEFAULT_PATTERNS
THRESHOLD_AMPLITUDE_PA = 1.0  # Fit floor; not the physical-cycle detector.


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def finite_or_none(value: float) -> float | None:
    return float(value) if np.isfinite(value) else None


def metadata_number(text: str, key: str) -> float:
    match = re.search(r"^\s*" + re.escape(key) + r"\s*=\s*([^,\n]+)",
                      text, flags=re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise ValueError(f"Missing metadata key {key}")
    return float(match.group(1).strip().replace("D", "E").replace("d", "e"))


def parse_status(text: str) -> dict:
    result = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = (x.strip() for x in line.split("=", 1))
        try:
            result[key] = float(value.replace("D", "E"))
        except ValueError:
            result[key] = value
    return result


def refined_extremum(t: np.ndarray, y: np.ndarray, index: int,
                     maximum: bool) -> tuple[float, float]:
    """Three-point local quadratic vertex; never extrapolate outside neighbors."""
    if index == 0 or index == len(t) - 1:
        return float(t[index]), float(y[index])
    x = t[index - 1:index + 2] - t[index]
    z = y[index - 1:index + 2] - y[index]
    a, b, c = np.polyfit(x, z, 2)
    if a == 0 or (maximum and a >= 0) or (not maximum and a <= 0):
        return float(t[index]), float(y[index])
    vertex = -b / (2 * a)
    if not x[0] <= vertex <= x[-1]:
        return float(t[index]), float(y[index])
    return float(t[index] + vertex), float(y[index] + (a * vertex + b) * vertex + c)


def pressure_envelope(t: np.ndarray, y: np.ndarray) -> list[dict]:
    sample_step = float(np.median(np.diff(t)))
    distance = max(1, int(10. / sample_step))
    peaks = find_peaks(y, distance=distance, prominence=1.e-5)[0]
    troughs = find_peaks(-y, distance=distance, prominence=1.e-5)[0]
    minima = [refined_extremum(t, y, int(i), False) for i in troughs]
    result = []
    for index in peaks:
        k = int(np.searchsorted(troughs, index))
        if k == 0 or k == len(troughs):
            continue
        left, right = minima[k - 1], minima[k]
        if not 8. < right[0] - left[0] < 40.:
            continue
        tp, yp = refined_extremum(t, y, int(index), True)
        baseline = left[1] + (tp - left[0]) / (right[0] - left[0]) * (right[1] - left[1])
        amplitude = .5 * (yp - baseline)
        if not np.isfinite(amplitude) or amplitude <= 0:
            continue
        result.append({"peak_time_s": tp, "peak_pressure_Pa": yp,
                       "left_trough_time_s": left[0], "right_trough_time_s": right[0],
                       "trough_baseline_Pa": baseline,
                       "semi_amplitude_Pa": amplitude,
                       "peak_to_interpolated_trough_Pa": 2 * amplitude})
    return result


def fit_envelope(envelope: list[dict], lo: float, hi: float,
                 observed_end: float) -> dict:
    rows = [r for r in envelope if lo <= r["peak_time_s"] <= hi
            and r["semi_amplitude_Pa"] >= THRESHOLD_AMPLITUDE_PA]
    result = {"window_start_s": lo, "window_end_s": hi,
              "full_window_observed": observed_end >= hi, "n_peaks": len(rows),
              "amplitude_floor_Pa": THRESHOLD_AMPLITUDE_PA}
    if len(rows) < 4:
        result["available"] = False
        return result
    t = np.array([r["peak_time_s"] for r in rows])
    amplitude = np.array([r["semi_amplitude_Pa"] for r in rows])
    regression = linregress(t, np.log(amplitude))
    result.update(available=True,
                  growth_rate_s_inv=float(regression.slope),
                  log_amplitude_intercept=float(regression.intercept),
                  r_squared=finite_or_none(regression.rvalue ** 2),
                  slope_standard_error_s_inv=float(regression.stderr),
                  first_peak_time_s=float(t[0]), last_peak_time_s=float(t[-1]),
                  first_amplitude_Pa=float(amplitude[0]), last_amplitude_Pa=float(amplitude[-1]),
                  mean_peak_interval_s=float(np.mean(np.diff(t))),
                  peak_interval_cv=float(np.std(np.diff(t), ddof=1) / np.mean(np.diff(t))))
    return result


def window_statistics(t: np.ndarray, y: np.ndarray, lo: float, hi: float) -> dict:
    """Exact moments of the piecewise-linear sampled signal, avoiding cancellation."""
    result = {"window_start_s": lo, "window_end_s": hi,
              "full_window_observed": bool(t[0] <= lo and t[-1] >= hi)}
    if not result["full_window_observed"]:
        return result
    use = (t > lo) & (t < hi)
    tx = np.r_[lo, t[use], hi]
    yy = np.r_[np.interp(lo, t, y), y[use], np.interp(hi, t, y)]
    mean = float(np.sum(.5 * (yy[:-1] + yy[1:]) * np.diff(tx)) / (hi - lo))
    z = yy - mean
    variance = float(np.sum(np.diff(tx) * (z[:-1] ** 2 + z[:-1] * z[1:] + z[1:] ** 2) / 3)
                     / (hi - lo))
    result.update(mean_Pa=mean, rms_about_mean_Pa=math.sqrt(max(0., variance)),
                  peak_to_peak_Pa=float(np.ptp(yy)), minimum_Pa=float(np.min(yy)),
                  maximum_Pa=float(np.max(yy)))
    return result


def last_ten_cycles(t: np.ndarray, y: np.ndarray, envelope: list[dict]) -> dict:
    """Eleven consecutive enclosed peaks define ten complete peak-to-peak cycles."""
    if len(envelope) < 11:
        return {"available": False, "reason": "fewer than eleven enclosed peaks"}
    last = envelope[-11:]
    times = np.array([r["peak_time_s"] for r in last])
    amplitudes = np.array([r["semi_amplitude_Pa"] for r in last])
    periods = np.diff(times)
    if np.any((periods < 8) | (periods > 40)):
        return {"available": False, "reason": "nonconsecutive oscillations"}
    statistics = window_statistics(t, y, float(times[0]), float(times[-1]))
    regression = linregress(times, np.log(amplitudes))
    relative_change = float(amplitudes[-1] / amplitudes[0] - 1)
    useful = bool(np.min(2 * amplitudes) >= 100 and statistics["rms_about_mean_Pa"] >= 1)
    return {"available": True, "n_complete_cycles": 10,
            "window_start_s": float(times[0]), "window_end_s": float(times[-1]),
            "window_is_end_of_run": bool(t[-1] - times[-1] < np.max(periods)),
            "mean_period_s": float(np.mean(periods)),
            "period_cv": float(np.std(periods, ddof=1) / np.mean(periods)),
            "first_semi_amplitude_Pa": float(amplitudes[0]),
            "last_semi_amplitude_Pa": float(amplitudes[-1]),
            "mean_semi_amplitude_Pa": float(np.mean(amplitudes)),
            "amplitude_cv": float(np.std(amplitudes, ddof=1) / np.mean(amplitudes)),
            "relative_amplitude_change": relative_change,
            "relative_difference_last_five_vs_first_five_amplitudes":
                float(np.mean(amplitudes[-5:]) / np.mean(amplitudes[:5]) - 1),
            "growth_rate_s_inv": float(regression.slope),
            "log_amplitude_r_squared": finite_or_none(regression.rvalue ** 2),
            "useful_above_declared_detection_thresholds": useful,
            "amplitude_change_within_one_percent": bool(abs(relative_change) <= .01),
            "amplitude_change_within_five_percent": bool(abs(relative_change) <= .05),
            "settled_or_physically_validated": False,
            "pressure_statistics": statistics,
            "peak_records": last}


def read_run(path: Path) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    raw = (path / "monitor.csv").read_bytes()
    if not raw.endswith(b"\n"):
        raw = raw[:raw.rfind(b"\n") + 1]
    meta_raw = (path / "run_metadata.txt").read_bytes()
    status_raw = (path / "run_status.txt").read_bytes() if (path / "run_status.txt").exists() else b""
    data = np.genfromtxt(io.BytesIO(raw), delimiter=",", names=True)
    data = np.atleast_1d(data)
    t = data["time_s"]
    pressure = data["p_junction_Pa"]
    if len(t) < 4 or np.any(np.diff(t) <= 0) or np.any(~np.isfinite(pressure)):
        raise ValueError(f"Invalid monitor signal: {path}")
    meta = meta_raw.decode()
    status = parse_status(status_raw.decode())
    requested = metadata_number(meta, "TIME_DIM")
    completed = bool(status.get("status") == 0 and t[-1] >= requested - 1.e-7)
    envelope = pressure_envelope(t, pressure)
    result = {"case": path.name, "source_directory": str(path.resolve()),
              "monitor_sha256_bytes_read": sha256(raw),
              "metadata_sha256": sha256(meta_raw),
              "run_status_sha256": sha256(status_raw) if status_raw else None,
              "run_status": status, "completed": completed,
              "run_status_present": bool(status_raw),
              "label": "completed" if completed else "incomplete_at_analysis",
              "observed_end_s": float(t[-1]), "requested_end_s": requested,
              "nodes_pipe": int(metadata_number(meta, "NODES")),
              "nodes_riser": int(metadata_number(meta, "N_RISER")),
              "dt_max_s": metadata_number(meta, "DT_MAX_S"),
              "requested_output_interval_s": metadata_number(meta, "OUTPUT_INTERVAL_S"),
              "sample_interval_median_s": float(np.median(np.diff(t))),
              "actual_step_range_40_400_s":
                  [float(np.min(data["dt_s"][(t >= 40) & (t <= 400)])),
                   float(np.max(data["dt_s"][(t >= 40) & (t <= 400)]))]
                  if t[-1] >= 40 else None,
              "pressure_column": "p_junction_Pa", "sample_count": len(t),
              "envelope": envelope,
              "fits": [fit_envelope(envelope, lo, hi, float(t[-1])) for lo, hi in WINDOWS],
              "pressure_blocks": [window_statistics(t, pressure, lo, hi) for lo, hi in BLOCKS],
              "pressure_900_2400": window_statistics(t, pressure, 900., 2400.),
              "last_ten_cycles": last_ten_cycles(t, pressure, envelope)}
    # Refined runs save more frequently.  Repeat the initial-window fit on a
    # common 0.125 s grid to expose the effect of sampling separately from h.
    if result["dt_max_s"] <= .125 and t[-1] >= 450:
        common_t = np.arange(0., 450.00001, .125)
        common_pressure = np.interp(common_t, t, pressure)
        common_fit = fit_envelope(pressure_envelope(common_t, common_pressure), 40., 400., 450.)
        result["common_sample_step_s"] = .125
        result["initial_fit_common_sampling"] = common_fit
    return result, t, pressure, data["dt_s"]


def initial_fit(case: dict) -> dict | None:
    return next((x for x in case["fits"] if x["window_start_s"] == 40 and x["window_end_s"] == 400
                 and x.get("available") and x["full_window_observed"]), None)


def backward_euler_model(cases: list[dict]) -> dict:
    baseline = []
    for h in (.5, .25, .125):
        case = next((c for c in cases if c["nodes_pipe"] == 51 and c["nodes_riser"] == 51
                     and math.isclose(c["dt_max_s"], h) and initial_fit(c)), None)
        if case is None:
            return {"available": False}
        baseline.append(case)
    h = np.array([c["dt_max_s"] for c in baseline])
    growth = np.array([initial_fit(c)["growth_rate_s_inv"] for c in baseline])
    # g(h)=-log((1-h*sigma)^2+(h*omega)^2)/(2h)
    # expm1(-2*h*g)/h = -2*sigma+h*(sigma^2+omega^2).
    k, intercept = np.polyfit(h, np.expm1(-2 * h * growth) / h, 1)
    sigma = -intercept / 2
    omega2 = k - sigma ** 2
    if sigma <= 0 or omega2 <= 0:
        return {"available": False, "reason": "No positive oscillatory descriptive fit"}
    omega = math.sqrt(omega2)
    predictions = []
    for case in sorted(cases, key=lambda c: -c["dt_max_s"]):
        fit = initial_fit(case)
        if case["nodes_pipe"] != 51 or case["nodes_riser"] != 51 or fit is None:
            continue
        step = case["dt_max_s"]
        predicted_g = -math.log((1 - step * sigma) ** 2 + (step * omega) ** 2) / (2 * step)
        predicted_frequency = math.atan2(step * omega, 1 - step * sigma) / step
        predicted_period = 2 * math.pi / predicted_frequency
        predictions.append({"case": case["case"], "dt_s": step,
                            "used_to_fit": step in (.5, .25, .125),
                            "observed_growth_s_inv": fit["growth_rate_s_inv"],
                            "predicted_growth_s_inv": predicted_g,
                            "growth_error_s_inv": predicted_g - fit["growth_rate_s_inv"],
                            "observed_period_s": fit["mean_peak_interval_s"],
                            "predicted_period_s": predicted_period,
                            "period_relative_error": predicted_period / fit["mean_peak_interval_s"] - 1})
    order = None
    fine = []
    for step in (.125, .0625, .03125):
        item = next((x for x in predictions if math.isclose(x["dt_s"], step)), None)
        if item is None:
            break
        fine.append(item["observed_growth_s_inv"])
    if len(fine) == 3:
        delta1, delta2 = fine[1] - fine[0], fine[2] - fine[1]
        if delta1 * delta2 > 0:
            order = {"dt_s": [.125, .0625, .03125], "growth_rates_s_inv": fine,
                     "increment_ratio": delta1 / delta2,
                     "observed_order_for_growth_rate_only": math.log2(delta1 / delta2),
                     "linear_extrapolated_growth_s_inv": 2 * fine[2] - fine[1],
                     "global_convergence_or_gci_claim": False}
    return {"available": True, "interpretation": "descriptive signal fit; not a proven operator eigenvalue",
            "fit_window_s": [40, 400], "fit_cases": [c["case"] for c in baseline],
            "sigma_effective_s_inv": float(sigma), "omega_effective_rad_s": omega,
            "continuous_period_of_descriptive_model_s": 2 * math.pi / omega,
            "critical_step_of_descriptive_model_s": float(2 * sigma / k),
            "predictions": predictions, "temporal_order_indicator": order,
            "limitation": "LAGGED freezes the previous-step Roe dissipation: its full linearization is a split operator. "
                          "A scalar backward-Euler fit to amplitudes is therefore illustrative; its period mismatch is reported."}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, keys)
        writer.writeheader()
        writer.writerows(rows)


def figures(output: Path, cases: list[dict], signals: dict, model: dict) -> None:
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .22,
                         "savefig.dpi": 180, "pdf.fonttype": 42})
    temporal = sorted([c for c in cases if c["nodes_pipe"] == 51 and c["nodes_riser"] == 51],
                      key=lambda c: -c["dt_max_s"])
    fig, ax = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True)
    for i, case in enumerate(temporal):
        rows = case["envelope"]
        use = [r for r in rows if r["semi_amplitude_Pa"] >= 1.e-3]
        label = f"Δt máx. = {case['dt_max_s']:g} s" + (" (em curso)" if not case["completed"] else "")
        ax[0].semilogy([r["peak_time_s"] for r in use], [r["semi_amplitude_Pa"] for r in use],
                       color=f"C{i}", label=label, linewidth=1.8)
        fit = initial_fit(case)
        if fit:
            ax[1].plot(case["dt_max_s"], fit["growth_rate_s_inv"], "o", color=f"C{i}", markersize=7)
    ax[0].axvspan(40, 400, color="gray", alpha=.10, label="Janela do ajuste: 40–400 s")
    ax[0].set(xlabel="Tempo (s)", ylabel="Semi-amplitude da pressão na junção (Pa)",
              title="Caso OSC, N = 51: amortecimento ou crescimento dependem do passo")
    ax[0].legend(ncol=2, fontsize=9, loc="lower right")
    if model.get("available"):
        h = np.linspace(.001, .52, 300)
        s, w = model["sigma_effective_s_inv"], model["omega_effective_rad_s"]
        g = -np.log((1 - h * s) ** 2 + (h * w) ** 2) / (2 * h)
        ax[1].plot(h, g, color="black", lw=1.2, label="Modelo escalar Euler implícito: ajuste descritivo")
        ax[1].axvline(model["critical_step_of_descriptive_model_s"], ls=":", color="gray")
    ax[1].axhline(0, color="gray", linewidth=.8)
    ax[1].set(xlabel="Passo máximo de tempo (s)", ylabel="Taxa ajustada de log-amplitude (s⁻¹)",
              title="Ajuste dos picos em 40–400 s; os pontos não são autovalores demonstrados")
    ax[1].legend(fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"osc_envelopes.{suffix}")
    plt.close(fig)

    completed_oscillatory = [c for c in temporal if c["completed"]
                            and c["last_ten_cycles"].get("available")
                            and c["last_ten_cycles"].get("useful_above_declared_detection_thresholds")]
    if completed_oscillatory:
        fig, axes = plt.subplots(len(completed_oscillatory), 2,
                                 figsize=(12, 3.3 * len(completed_oscillatory)),
                                 squeeze=False, constrained_layout=True)
        for row, case in enumerate(completed_oscillatory):
            info = case["last_ten_cycles"]
            t, p, _ = signals[case["case"]]
            use = (t >= info["window_start_s"]) & (t <= info["window_end_s"])
            axes[row, 0].plot(t[use], p[use] / 1000, linewidth=1)
            axes[row, 0].set(xlabel="Tempo (s)", ylabel="Pressão na junção (kPa)",
                             title=f"N = 51, Δt máx. = {case['dt_max_s']:g} s: últimos 10 ciclos completos")
            records = info["peak_records"]
            axes[row, 1].plot(range(11), [r["semi_amplitude_Pa"] / 1000 for r in records], "o-")
            change = 100 * info["relative_amplitude_change"]
            axes[row, 1].set(xlabel="Índice do pico (11 picos definem 10 ciclos)",
                             ylabel="Semi-amplitude (kPa)",
                             title=f"Variação: {change:+.3f}%; período médio: {info['mean_period_s']:.4f} s")
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"osc_envelopes_last_cycles.{suffix}")
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for n in (51, 101, 201):
        subset = sorted([c for c in cases if c["nodes_pipe"] == n and c["nodes_riser"] == n and initial_fit(c)],
                        key=lambda c: c["dt_max_s"])
        if subset:
            ax.plot([c["dt_max_s"] for c in subset], [initial_fit(c)["growth_rate_s_inv"] for c in subset],
                    "o-", label=f"N = {n}")
    ax.axhline(0, color="gray", linewidth=.8)
    ax.set(xlabel="Passo máximo de tempo (s)", ylabel="Taxa ajustada de log-amplitude (s⁻¹)",
           title="A malha espacial também altera o crescimento observado (40–400 s)")
    ax.legend()
    for suffix in ("png", "pdf"):
        fig.savefig(output / f"osc_envelopes_spatial.{suffix}")
    plt.close(fig)


def report(cases: list[dict], model: dict, generated: str) -> str:
    lines = ["# Envelopes do caso OSC e refinamento temporal", "", f"Gerado em {generated}.", "",
             "A classificação do sinal muda com o passo de tempo. Os resultados abaixo são observáveis numéricos "
             "da pressão na junção; não demonstram autovalores do operador, convergência global, um ciclo limite "
             "assintótico ou validação física de severe slugging.", "",
             "## Método e janelas", "",
             "Os máximos e mínimos locais são refinados por uma parábola de três amostras, sem extrapolação. "
             "Cada pico usa os dois mínimos vizinhos para interpolar a linha de base. A semi-amplitude é metade "
             "da diferença entre o pico e essa linha de base. O ajuste principal é a regressão linear de "
             "log(amplitude) em **40–400 s**, usando amplitudes ≥ 1 Pa. A separação mínima entre extremos do "
             "mesmo tipo é 10 s e somente pares de mínimos separados por 8–40 s são associados a um ciclo. "
             "Esses parâmetros medem o sinal observado; não definem um regime físico.", "",
             "O diagnóstico de amortecimento pode usar sinais pequenos. Para rotular dez ciclos finais como "
             "úteis, exige-se excursão pico–mínimo ≥ 100 Pa em todos eles e RMS da pressão em torno da média "
             "≥ 1 Pa. Onze picos consecutivos definem exatamente dez intervalos. O RMS é integrado no tempo "
             "sobre a interpolação linear das amostras. Dados ainda em execução são explicitamente identificados.", "",
             "Os refinados 0,0625 e 0,03125 s também gravam o monitor com intervalo menor. O JSON registra "
             "essa cadência e repete o ajuste inicial em uma grade comum de amostragem de 0,125 s, por "
             "interpolação linear. Isso permite distinguir o efeito da gravação do efeito do passo do solver; "
             "a tabela principal usa as amostras nativas.", "",
             "## Crescimento inicial observado", "",
             "| N | Δt máximo (s) | Fim observado (s) | Execução completa | Taxa (s⁻¹) | R² | Picos | Intervalo médio entre picos (s) |",
             "|---:|---:|---:|:---:|---:|---:|---:|---:|"]
    for case in sorted(cases, key=lambda c: (c["nodes_pipe"], -c["dt_max_s"])):
        fit = initial_fit(case)
        if fit:
            r_squared = f"{fit['r_squared']:.9f}" if fit['r_squared'] is not None else "indefinido"
            lines.append(f"| {case['nodes_pipe']} | {case['dt_max_s']:g} | {case['observed_end_s']:.3f} | "
                         f"{'sim' if case['completed'] else 'não'} | {fit['growth_rate_s_inv']:+.9f} | "
                         f"{r_squared} | {fit['n_peaks']} | {fit['mean_peak_interval_s']:.6f} |")
    lines += ["", "O alto R² descreve a regularidade exponencial desta janela. O erro padrão da regressão não "
              "é uma incerteza de modelo: os ciclos vêm de uma única trajetória numérica, e não de réplicas independentes.", "",
              "## Acomodação nas execuções completas", "",
              "| N | Δt (s) | RMS 900–1400 / 1400–1900 / 1900–2400 (Pa) | Últimos 10 ciclos úteis | Variação de amplitude | Período médio (s) | CV do período |",
              "|---:|---:|---|:---:|---:|---:|---:|"]
    for case in sorted(cases, key=lambda c: (c["nodes_pipe"], -c["dt_max_s"])):
        if not case["completed"]:
            continue
        rms = " / ".join(f"{b['rms_about_mean_Pa']:.6g}" if b["full_window_observed"] else "—"
                         for b in case["pressure_blocks"])
        last = case["last_ten_cycles"]
        useful = last.get("useful_above_declared_detection_thresholds", False)
        change = f"{100 * last['relative_amplitude_change']:+.4f}%" if useful else "—"
        period = f"{last['mean_period_s']:.6f}" if useful else "—"
        cv = f"{last['period_cv']:.6g}" if useful else "—"
        lines.append(f"| {case['nodes_pipe']} | {case['dt_max_s']:g} | {rms} | "
                     f"{'sim' if useful else 'não'} | {change} | {period} | {cv} |")
    lines += ["", "A estabilidade do período não garante a acomodação da amplitude. O JSON registra separadamente "
              "se a amplitude dos últimos dez ciclos varia menos de 1% ou de 5%; essas duas verificações são "
              "descrições, não critérios universais de validação ou demonstração de um ciclo limite.", "",
              "## Por que Euler implícito pode ocultar o crescimento", "",
              "Para o problema escalar z′ = (σ + iω)z, a amplificação de Euler implícito é", "",
              "```text", "q(h) = 1 / [1 − h(σ + iω)]",
              "g_BE(h) = log|q| / h = −log[(1 − hσ)² + (hω)²] / (2h)",
              "ω_BE(h) = atan2(hω, 1 − hσ) / h", "```", "",
              "Mesmo para σ > 0, há crescimento somente se h < 2σ/(σ² + ω²). A estabilidade absoluta no "
              "semiplano esquerdo não garante a reprodução de instabilidades no semiplano direito. Para h "
              "pequeno, g_BE ≈ σ + h(σ² − ω²)/2: em um modo oscilatório, o termo temporal pode acrescentar "
              "amortecimento suficiente para inverter o sinal observado.", ""]
    if model.get("available"):
        lines += [f"O ajuste descritivo dos envelopes N = 51 com h = 0,5; 0,25; 0,125 s produz "
                  f"σ efetivo = {model['sigma_effective_s_inv']:.9g} s⁻¹ e ω efetivo = "
                  f"{model['omega_effective_rad_s']:.9g} rad/s, com passo de mudança de sinal "
                  f"≈ {model['critical_step_of_descriptive_model_s']:.6f} s. **Esses parâmetros não são "
                  "autovalores demonstrados.** O modo LAGGED congela a dissipação de Roe do passo anterior; "
                  "sua linearização completa é uma separação de operadores, e não necessariamente a fórmula "
                  "escalar de um único operador fixo.", "",
                  "| Δt (s) | Usado no ajuste | Taxa observada (s⁻¹) | Previsão BE (s⁻¹) | Erro na previsão do período |",
                  "|---:|:---:|---:|---:|---:|"]
        for prediction in model["predictions"]:
            lines.append(f"| {prediction['dt_s']:g} | {'sim' if prediction['used_to_fit'] else 'não'} | "
                         f"{prediction['observed_growth_s_inv']:+.9f} | {prediction['predicted_growth_s_inv']:+.9f} | "
                         f"{100 * prediction['period_relative_error']:+.3f}% |")
        order = model["temporal_order_indicator"]
        if order:
            lines += ["", f"Nos passos 0,125; 0,0625; 0,03125 s, a razão entre os incrementos da taxa "
                      f"é {order['increment_ratio']:.6f}, equivalente à ordem observada "
                      f"{order['observed_order_for_growth_rate_only']:.6f} para **esta grandeza e janela**. "
                      f"A extrapolação linear dos dois passos menores dá "
                      f"{order['linear_extrapolated_growth_s_inv']:.9f} s⁻¹. Isso é um indício de erro "
                      "temporal de primeira ordem no crescimento inicial, sem constituir GCI global ou "
                      "convergência da amplitude final, do período, do campo transiente ou do regime físico."]
    lines += ["", "## Reprodutibilidade e limites", "",
              "O JSON contém hashes dos bytes lidos de monitor.csv, run_metadata.txt e run_status.txt, "
              "janelas completas/incompletas, R², pisos, valores de cada pico e resultados de todos os ajustes. "
              "O CSV osc_envelopes.csv permite reconstruir as curvas de amplitude; osc_envelopes_fits.csv "
              "lista as regressões. Os PNG e PDF são figuras independentes para exportação.", "",
              "As execuções refinadas em curso podem mudar suas estatísticas tardias. Reexecute este script "
              "após a conclusão. Nenhum alvo experimental é usado para ajustar as taxas, o modelo descritivo "
              "ou os limiares. Captar uma oscilação de pressão não comprova a sequência física de acumulação, "
              "compressão, produção, penetração de gás e blowdown. Bloqueio completo do gás não é um "
              "requisito universal de OSC/SS3; a análise dessas fases exige fluxos e inventários.", ""]
    return "\n".join(lines)


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=repo / "validation/results")
    parser.add_argument("--output", type=Path, default=repo / "validation/analysis_final")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--pattern", action="append", help="Repeatable glob; replaces historical defaults")
    selection.add_argument("--run-dir", action="append", type=Path,
                           help="Explicit run directory; repeat for independent mesh/time refinements")
    args = parser.parse_args()
    try:
        paths = select_runs(args.results_root, args.pattern, args.run_dir)
        cohort = cohort_metadata(paths)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    args.output.mkdir(parents=True, exist_ok=True)
    cases, signals, errors = [], {}, []
    for path in paths:
        try:
            case, t, pressure, steps = read_run(path)
        except (ValueError, OSError) as exc:
            errors.append({"case": path.name, "error": str(exc)})
            continue
        cases.append(case)
        signals[case["case"]] = (t, pressure, steps)
    if not cases:
        raise SystemExit("No usable OSC monitors found")
    generated = datetime.now(timezone.utc).isoformat()
    model = backward_euler_model(cases)
    summary = {"generated_utc": generated, "script_sha256": sha256(Path(__file__).read_bytes()),
               "selection": {"patterns": args.pattern if not args.run_dir else None,
                             "used_historical_defaults": not (args.pattern or args.run_dir),
                             "run_directories": [str(p) for p in paths]},
               "cohort": {"parameters": cohort[0]["parameters"],
                          "source_sha256": cohort[0]["manifest"]["source_sha256"]},
               "interpretation": "numerical signal diagnostics; no physical validation claim",
               "main_fit_window_s": [40, 400], "fit_amplitude_floor_Pa": THRESHOLD_AMPLITUDE_PA,
               "cycle_detection_thresholds": {"minimum_peak_to_trough_Pa": 100, "minimum_rms_Pa": 1},
               "cases": cases, "backward_euler_descriptive_model": model, "read_errors": errors}
    (args.output / "osc_envelopes.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    write_csv(args.output / "osc_envelopes.csv",
              [{"case": c["case"], "N": c["nodes_pipe"], "dt_max_s": c["dt_max_s"],
                "completed": c["completed"], **row} for c in cases for row in c["envelope"]])
    write_csv(args.output / "osc_envelopes_fits.csv",
              [{"case": c["case"], "N": c["nodes_pipe"], "dt_max_s": c["dt_max_s"],
                "completed": c["completed"], **row} for c in cases for row in c["fits"]])
    figures(args.output, cases, signals, model)
    (args.output / "osc_envelopes.md").write_text(report(cases, model, generated))
    print(json.dumps({"cases": len(cases), "completed": sum(c["completed"] for c in cases),
                      "errors": errors, "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    main()
