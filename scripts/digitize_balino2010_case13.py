#!/usr/bin/env python3
"""Limited digitization of experimental Fig.18(a), Baliño et al. 2010.

Requires PyMuPDF, NumPy and matplotlib. Outputs belong outside the repository.
The input PDF is private and is never copied to the output directory.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

PDF_SHA256 = 'ce04a2d1b79d3a10d442fbaeb3863e4639c8508bd34ba9e1a414ae03e5479e2b'


def main():
    import pymupdf
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pdf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--historical-period-s', type=float, default=23.10864)
    args = parser.parse_args()
    if not np.isfinite(args.historical_period_s) or args.historical_period_s <= 0:
        parser.error('Historical period must be positive and finite')
    if hashlib.sha256(args.pdf.read_bytes()).hexdigest() != PDF_SHA256:
        parser.error('PDF identity differs from the inspected source')
    if args.output_dir.exists():
        parser.error('Output directory exists; preserve earlier evidence')
    doc = pymupdf.open(args.pdf)
    item = doc.extract_image(263)  # image on PDF15 = printed p.657, Fig.18
    pix = pymupdf.Pixmap(item['image'])
    if (pix.width, pix.height, pix.n) != (1256, 434, 3):
        raise ValueError('Native raster layout differs from inspected source')
    rgb = np.frombuffer(pix.samples, dtype=np.uint8).reshape(434, 1256, 3).astype(int)
    # Native image coordinates, manually inspected outer axes and tick labels.
    x0, x1, y0, y1 = 36, 579, 13, 371
    t0, t1, p0, p1 = 1200.0, 1500.0, 3.5e5, 1.5e5
    dt, dp = (t1-t0)/(x1-x0), (p0-p1)/(y1-y0)
    blue = (rgb[:,:,2] > 150) & (rgb[:,:,0] < 120) & (rgb[:,:,1] < 200)
    blue &= (rgb[:,:,2]-rgb[:,:,0] > 80) & (rgb[:,:,2]-rgb[:,:,1] > 40)
    rows = []
    for x in range(x0, x1+1):
        ys = np.flatnonzero(blue[y0:y1+1, x]) + y0
        time = t0+(x-x0)*dt
        # Pixel cell half-width + one pixel manual axis calibration allowance.
        row = dict(native_pixel_x=x, time_s=time, time_lower_s=time-1.5*dt,
                   time_upper_s=time+1.5*dt, status='missing' if not len(ys) else 'observed')
        if len(ys):
            low = p0-(ys[-1]-y0+1.5)*dp
            high = p0-(ys[0]-y0-1.5)*dp
            row.update(p_base_Pa=(low+high)/2, p_lower_Pa=low, p_upper_Pa=high,
                       blue_pixels=int(len(ys)))
        else:
            row.update(p_base_Pa=None,p_lower_Pa=None,p_upper_Pa=None,blue_pixels=0)
        rows.append(row)
    args.output_dir.mkdir(parents=True)
    csv_path=args.output_dir/'balino2010_case13_fig18_digitized.csv'
    with csv_path.open('x',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    summary=dict(schema='PRIMS-F.balino2010_digitization.v1',
                 source=dict(doi='10.1016/j.ijmultiphaseflow.2010.04.003',case=13,
                             printed_page=657,pdf_page=15,figure='18(a)',pdf_sha256=PDF_SHA256,
                             native_image_xref=263,native_image_sha256=hashlib.sha256(item['image']).hexdigest(),
                             native_image_shape=[434,1256,3]),
                 analyzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 calibration=dict(x_pixel=[x0,x1],time_s=[t0,t1],y_pixel=[y0,y1],pressure_Pa=[p0,p1],
                                  manual_axis_allowance_px=1.0,pixel_half_width_px=0.5,
                                  time_resolution_s=dt,pressure_resolution_Pa=dp,
                                  blue_mask='B>150, R<120, G<200, B-R>80, B-G>40'),
                 rows=len(rows),missing_columns=sum(r['status']=='missing' for r in rows),
                 nominal_time_coverage_s=[t0,t1],
                 csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                 conditional_period_comparison=dict(experimental_interval_s=[19.3,33.8],
                                                    historical_numerical_period_s=args.historical_period_s,
                                                    numerical_freeze='8a85e37319995b37795d2c798d6f865b863b1834',
                                                    inside_published_range=19.3<=args.historical_period_s<=33.8),
                 limits=['Digitized published raster, not the original experimental 10 Hz record.',
                         'Central values are midpoints of the blue vertical span per raster column, not recovered raw samples.',
                         'Bands cover unresolved fluctuations/line thickness, pixel cell and declared axis allowance; not statistical confidence.',
                         'Instrument and plotting uncertainties beyond that allowance remain unknown.',
                         'No interpolation fills absent columns; no absolute-time alignment with a simulation.',
                         'The plotted 300 s excerpt is not the complete approximately 1500 s acquisition.',
                         'Table3 Max 0.69 bar is not replaced by the excursion of this digitized excerpt.',
                         'Period inclusion is conditional consistency, not model validation.'])
    (args.output_dir/'balino2010_case13_fig18_digitization.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    time=np.array([r['time_s'] for r in rows])
    pressure=np.array([r['p_base_Pa'] if r['p_base_Pa'] is not None else np.nan for r in rows])/1e5
    lower=np.array([r['p_lower_Pa'] if r['p_lower_Pa'] is not None else np.nan for r in rows])/1e5
    upper=np.array([r['p_upper_Pa'] if r['p_upper_Pa'] is not None else np.nan for r in rows])/1e5
    fig,ax=plt.subplots(figsize=(10,4.1),layout='constrained')
    ax.fill_between(time,lower,upper,color='#5e87a2',alpha=.4,label='Envelope grafico declarado')
    ax.plot(time,pressure,color='#164e73',lw=.9,label='Centro do traco experimental digitalizado')
    ax.set(xlabel='Tempo experimental publicado (s)',ylabel='Pressao absoluta de base (bar)',
           title='Baliño et al. 2010 — caso 13, Fig.18(a), p.657')
    ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.text(.5,-.025,'Raster publicado: ~0,55 s/px; ~559 Pa/px. Nao e a serie bruta de 10 Hz.',ha='center',fontsize=9)
    fig.savefig(args.output_dir/'balino2010_case13_experimental_digitized.png',dpi=180,bbox_inches='tight')
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,2.8),layout='constrained')
    ax.plot([19.3,33.8],[1,1],lw=8,color='#86aec4',solid_capstyle='butt',label='Experimental: Tabela 3, caso 13')
    ax.plot([19.3,33.8],[1,1],'|',color='#164e73',markersize=18)
    ax.plot(args.historical_period_s,0,'o',color='#8e3b32',label='Historico numerico: freeze 8a85e37')
    ax.annotate(f'{args.historical_period_s:.5f} s',(args.historical_period_s,0),xytext=(0,12),textcoords='offset points',ha='center')
    ax.set(xlabel='Periodo (s)',yticks=[0,1],yticklabels=['Simulacao historica','Faixa experimental'],ylim=[-.45,1.45],xlim=[17,36],
           title='Comparacao condicional — inclusao na faixa nao e validacao')
    ax.grid(axis='x',alpha=.2)
    fig.savefig(args.output_dir/'balino2010_case13_period_conditional.png',dpi=180)
    plt.close(fig)
    print(json.dumps(dict(rows=len(rows),missing_columns=summary['missing_columns'],output_dir=str(args.output_dir))))


if __name__ == '__main__':
    main()
