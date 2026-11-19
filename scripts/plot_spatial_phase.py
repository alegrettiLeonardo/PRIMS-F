#!/usr/bin/env python3
"""Plot previously computed fields at fixed positions and a common pressure phase."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    index = json.loads((args.folder/'index.json').read_text())
    if len(index['runs']) != 2:
        raise ValueError('A reference and one comparison are required')
    if args.output.exists():
        raise ValueError('Choose a new output')
    a, b = [np.load(args.folder/(run['run']+'_phase_fields.npz')) for run in index['runs']]
    fig, axes = plt.subplots(3, 3, figsize=(13, 9), sharex=True, sharey=True, layout='constrained')
    for row, (channel, factor, label) in enumerate(((0, 1e-5, 'Pressão absoluta (bar)'),
            (1, 1, 'Fração líquida'), (5, 1, 'Velocidade superficial gasosa (m/s)'))):
        values = [np.concatenate([item[segment][:, :, channel] for segment in ('pipe', 'riser')], axis=1)*factor
                  for item in (a, b)]
        low, high = min(v.min() for v in values), max(v.max() for v in values)
        delta = values[1]-values[0]
        span = max(abs(delta.min()), abs(delta.max()))
        for column in range(3):
            ax = axes[row, column]
            for segment in ('pipe', 'riser'):
                left = a[segment][:, :, channel]*factor
                right = b[segment][:, :, channel]*factor
                field = (left, right, right-left)[column]
                mesh = ax.pcolormesh(a[segment+'_s_m'], a['phase'], field,
                    shading='nearest', cmap='coolwarm' if column == 2 else 'viridis',
                    vmin=-span if column == 2 else low, vmax=span if column == 2 else high,
                    rasterized=True)
            ax.set(title=label if column < 2 else 'Diferença: '+label)
            fig.colorbar(mesh, ax=ax, shrink=.85)
    for ax in axes[:, 0]:
        ax.set_ylabel('Fase (ciclo)')
    for ax in axes[-1]:
        ax.set_xlabel('Posição física ao longo do tubo (m)')
    for column, text in enumerate((index['runs'][0]['run'], index['runs'][1]['run'], 'Comparação em posições comuns')):
        axes[0, column].annotate(text, xy=(.5, 1.22), xycoords='axes fraction', ha='center', fontsize=9)
    fig.suptitle('Fase zero: pico de pressão da junção; média de dez intervalos completos', fontsize=11)
    fig.savefig(args.output, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    main()
