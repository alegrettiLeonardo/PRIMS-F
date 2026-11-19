#!/usr/bin/env python3
"""Compare static one-axis sensitivities with the published finite pressure excerpt."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_osc_experimental import digitized_pressure, spectrum
from summarize_osc_spectra import sha, spectral_peaks
from validation_evidence import load_runs, time_statistics


def load(folder, duration):
    record, = load_runs(folder.parent, pattern=folder.name)
    if not record['completed']:
        raise ValueError('Require completed simulation: '+str(folder))
    t = np.array(record['data']['time_s'])
    p = np.array(record['data']['p_junction_Pa'])
    start = record['reached_time_s']-duration
    use = t >= start
    if not np.isfinite(p[use]).all():
        raise ValueError('Nonfinite pressure')
    manifest = json.loads((folder/'manifest.json').read_text())
    return dict(record=record, time=t[use], pressure=p[use], controls=manifest['controls'],
        statistics=time_statistics(record['data'], 'p_junction_Pa', start, record['reached_time_s']),
        spectrum=spectral_peaks(t[use], p[use]), manifest=manifest,
        input_hashes={name: sha(folder/name) for name in ('manifest.json', 'monitor.csv', 'run_metadata.txt')})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nominal', type=Path, required=True)
    parser.add_argument('--digitization', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('runs', type=Path, nargs='+')
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Choose new output; preserve previous evidence')
    te, pe, lo, hi = digitized_pressure(args.digitization)
    nominal = load(args.nominal, te[-1]-te[0])
    groups = {key: [] for key in ('ps', 'jg', 'jl')}
    for folder in args.runs:
        run = load(folder, te[-1]-te[0])
        if (run['manifest']['git_head'] != nominal['manifest']['git_head'] or
                run['manifest']['binary_sha256'] != nominal['manifest']['binary_sha256']):
            raise ValueError('Operator or binary changed')
        axes = [key for key in groups if run['controls'][key] != nominal['controls'][key]]
        if len(axes) != 1:
            raise ValueError('Require one changed experimental control')
        groups[axes[0]].append(run)
    if any(len(group) != 2 for group in groups.values()):
        raise ValueError('Exactly two offsets per axis are required')
    args.output.mkdir(parents=True)
    fp, ap = plt.subplots(3, 1, figsize=(12, 10), sharex=True, sharey=True, layout='constrained')
    ff, af = plt.subplots(3, 1, figsize=(10, 10), sharex=True, layout='constrained')
    labels = {'ps': 'Ps absoluto (Pa)', 'jg': 'jg padrão (m/s)', 'jl': 'jl (m/s)'}
    expf, exppower = spectrum(te, pe)
    summaries = []
    for axis, axp, axf in zip(groups, ap, af):
        axp.fill_between(te-te[0], lo/1e5, hi/1e5, color='#a8bbd0', alpha=.5, label='Envelope gráfico experimental')
        axp.plot(te-te[0], pe/1e5, color='#235789', lw=1, label='Experimento, Fig.18(a)')
        use = expf <= .1
        axf.plot(expf[use], exppower[use], color='#235789', label='Experimento digitalizado')
        for run, color, label in [(nominal, '#555555', 'Nominal')]+[
                (run, color, labels[axis]+'='+format(run['controls'][axis], '.6g'))
                for run, color in zip(sorted(groups[axis], key=lambda r: r['controls'][axis]), ('#bc4a20', '#16836d'))]:
            t, p = run['time'], run['pressure']
            excursion = run['statistics']['peak_to_peak']/1e5
            axp.plot(t-t[0], p/1e5, color=color, lw=1.1, label=label+f' · excursão {excursion:.4g} bar')
            if run['spectrum']['window_excursion_at_least_100Pa']:
                freq, power = spectrum(t, p)
                use = freq <= .1
                axf.plot(freq[use], power[use], color=color, marker='.', ms=2, lw=1, label=label)
            else:
                axf.plot([], [], color=color, label=label+' · excursão <100 Pa; sem período qualificado')
        for ax in (axp, axf):
            ax.grid(alpha=.25)
            ax.legend(fontsize=8, frameon=False)
        axp.set(title=labels[axis]+' · contornos constantes', ylabel='Pressão absoluta (bar)')
        axf.set(title=labels[axis], ylabel='Potência relativa de cada série')
        for run in groups[axis]:
            summaries.append(dict(run=run['record']['run'], axis=axis, value=run['controls'][axis],
                time_window_s=[float(run['time'][0]), float(run['time'][-1])],
                statistics=run['statistics'], spectrum=run['spectrum'], input_hashes=run['input_hashes']))
    ap[-1].set_xlabel('Tempo na janela final de 300 s · origens independentes')
    af[-1].set_xlabel('Frequência (Hz)')
    fp.suptitle('Sensibilidade estática; os offsets não reconstruem o forçamento experimental', fontsize=11)
    ff.suptitle('Hann e média removida; bins finitos, sem inferência de modos físicos', fontsize=11)
    paths = (args.output/'pressao_contornos_fig18.png', args.output/'espectros_contornos_fig18.png')
    fp.savefig(paths[0], dpi=170)
    ff.savefig(paths[1], dpi=170)
    plt.close(fp)
    plt.close(ff)
    index = dict(schema='PRIMS-F.static_bc_experimental_figures.v1', script_sha256=sha(__file__),
        digitization_sha256=sha(args.digitization), nominal_hashes=nominal['input_hashes'],
        runs=summaries, figure_hashes={path.name: sha(path) for path in paths},
        limitations='Independent time origins; static offsets of measured-signal deviations; no synchronized forcing, confidence interval or experimental validation inferred.')
    with (args.output/'index.json').open('x') as stream:
        json.dump(index, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
