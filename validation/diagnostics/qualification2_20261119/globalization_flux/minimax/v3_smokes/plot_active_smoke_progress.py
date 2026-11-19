"""Scientific diagnostic of SS1 progress; no experimental-regime assertion."""
from pathlib import Path
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def complete_rows(path):
    with path.open() as stream:
        lines = (line for line in stream if line.endswith('\n'))
        return [row for row in csv.DictReader(lines)
                if None not in row and None not in row.values()]


root = Path('/home/leonardoab/Documentos/PRIMS-F-qualification2-results-20261119')
out = root / 'active_set_smokes_v3_a179e51'
reference_time = 78.57178847673931
cases = [('v1: L2 tardio (3940bb7)', root / 'active_set_smokes_v1/D_SS1'),
         ('v2: L2 antecipado', root / 'active_set_smokes_v2/D_SS1'),
         ('v3: L2 + L∞ (a179e51)', out / 'D_SS1/output')]
fig, axes = plt.subplots(3, 1, figsize=(9, 8.5), sharex=True, constrained_layout=True)
for label, folder in cases:
    steps = [row for row in complete_rows(folder / 'steps.csv')
             if float(row['time_s']) >= reference_time]
    cpu = [float(row['cpu_total_s']) for row in steps]
    delta = [(float(row['time_s']) - reference_time) * 1000 for row in steps]
    dt = [float(row['dt_s']) for row in steps]
    axes[0].plot(cpu, delta, label=label, linewidth=1.2)
    axes[1].plot(cpu, dt, linewidth=1.2)
    ledgers = {int(row['step']): row for row in complete_rows(folder / 'conservation.csv')}
    selected = [(float(row['cpu_total_s']), ledgers[int(row['step'])])
                for row in steps if int(row['step']) in ledgers]
    axes[2].plot([item[0] for item in selected],
                 [float(item[1]['mdot_g_base_left_kg_s']) for item in selected], linewidth=1.2)
axes[0].set_ylabel('Avanço após parada original (ms)')
axes[0].legend(loc='upper left', fontsize=9)
axes[0].set_title('SS1 — progresso aceito próximo ao limite do traço MATCHED')
axes[1].set_ylabel('dt aceito (s)')
axes[1].set_yscale('log')
axes[2].set_ylabel('Fluxo gasoso na junção (kg/s)')
axes[2].set_yscale('symlog', linthresh=1e-13)
axes[2].set_xlabel('CPU acumulada registrada no último passo aceito (s)')
for ax in axes:
    ax.grid(True, alpha=.25)
fig.savefig(out / 'SS1_progress.png', dpi=170)
print(out / 'SS1_progress.png')
