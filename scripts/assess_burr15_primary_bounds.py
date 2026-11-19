#!/usr/bin/env python3
"""Conservative graph-reading bounds, not raw experimental measurements.

Carraretto et al. (2020), DOI 10.3390/fluids5030101, Fig.9/p.10.
This calculation never imports the PRIMS-F implementation or closure.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

POINTS = (
    dict(id='figure9a_third_marker_from_left', jl_nominal_m_s=0.03,
         xv_reading_interval=[0.960, 0.967], alpha_reading_interval=[0.76, 0.80]),
    dict(id='figure9b_third_marker_from_left', jl_nominal_m_s=0.04,
         xv_reading_interval=[0.948, 0.954], alpha_reading_interval=[0.71, 0.76]),
)


def geometry(alpha, diameter):
    """Invert circular liquid area by height, independently of Fortran."""
    if not 0 < alpha < 1 or not diameter > 0:
        raise ValueError('Interior gas fraction and positive diameter required')
    lo, hi = 0.0, diameter
    for _ in range(100):
        h = (lo + hi) / 2
        phi = math.acos(1 - 2 * h / diameter)
        gas_area_fraction = 1 - (phi - math.sin(phi) * math.cos(phi)) / math.pi
        if gas_area_fraction > alpha:
            lo = h
        else:
            hi = h
    h = (lo + hi) / 2
    return h, diameter * math.acos(1 - 2 * h / diameter)


def assess(point):
    d, g = 0.060, 9.81
    # |d alpha/dh| <= 4/(pi D). Add height accuracy to the deliberately
    # coarse graph rectangle. This is an envelope, not a confidence interval.
    delta_alpha = 4 * 0.0005 / (math.pi * d)
    amin, amax = (point['alpha_reading_interval'][0] - delta_alpha,
                 point['alpha_reading_interval'][1] + delta_alpha)
    xmin, xmax = point['xv_reading_interval']
    jl = point['jl_nominal_m_s']
    # Applying both maximum reported relative flow uncertainties enlarges
    # the inferred range; no independence/statistical interpretation used.
    jgmin = jl * xmin / (1 - xmin) * (1 - 0.03) * (1 - 0.07)
    jgmax = jl * xmax / (1 - xmax) * (1 + 0.03) * (1 + 0.07)
    h, sl = geometry(amax, d)
    # Ucrit increases with alpha: (1-h/D), alpha and 1/Sl all increase.
    # Use (rho_l-rho_g)/rho_l <= 1, so no missing density is invented.
    critical_upper = (1 - h / d) * math.sqrt(g * math.pi * d*d / 4 * amax / sl)
    ug_lower = jgmin / amax
    return dict(**point, diameter_m=d, gravity_m_s2=g,
                alpha_envelope=[amin, amax], jg_envelope_m_s=[jgmin, jgmax],
                ug_lower_m_s=ug_lower, burr15_critical_upper_m_s=critical_upper,
                signed_margin_lower_m_s=ug_lower-critical_upper,
                squared_margin_lower_m2_s2=ug_lower**2-critical_upper**2,
                literal_detector_loss_for_entire_envelope=ug_lower >= critical_upper,
                observed_pattern='stratified, smooth interface',
                diameter_needed_to_remove_this_bound_m=d*(ug_lower/critical_upper)**2)


def result():
    return dict(schema='PRIMS-F.burr15_primary_bounds.v1',
                source=dict(doi='10.3390/fluids5030101', page=10, figure='9(a,b)',
                            pdf_sha256='960446bbdecef9acd5687978d8918390dde00ac58f045c2a00968144abbb7887',
                            uncertainty_source='pp.5-6, Table4; smooth observation p.7'),
                analyst_assumptions=dict(height_accuracy_m=0.0005,
                                         maximum_relative_jg_uncertainty=0.07,
                                         maximum_relative_jl_uncertainty=0.03,
                                         density_ratio_upper=1.0,
                                         uncertainty_kind='conservative envelope; not statistical CI'),
                results=[assess(p) for p in POINTS],
                limits=['Graph rectangles are manual readings, not tabulated raw points.',
                        'No reversal threshold, transition location or regime switching is inferred.',
                        'The geometric diameter is the reported nominal 60 mm.',
                        'This is a contradiction check of literal Burr15 with liquid wall perimeter.',
                        'It neither selects a replacement formula nor validates the full simulator.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    payload = result()
    payload['analyzer_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
