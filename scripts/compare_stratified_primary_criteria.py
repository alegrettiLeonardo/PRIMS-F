#!/usr/bin/env python3
"""Compare literal Burr15 and the primary 1977 stratified-loss criterion.

No regime selection or solver change. Missing experimental density is kept as
a conditional density-ratio threshold rather than replaced with a fitted value.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
from assess_burr15_primary_bounds import POINTS, assess, geometry


def td1977_factor(alpha, diameter, gravity, theta=0.0):
    """Ucrit**2 = (rho_l/rho_g - 1) * factor, horizontal Eq1/p6."""
    if theta != 0.0:
        raise ValueError('Only the source horizontal case is compared here')
    if not (math.isfinite(gravity) and gravity>0):
        raise ValueError('Positive finite gravity required')
    h, sl = geometry(alpha, diameter)
    phi=sl/diameter
    si=diameter*math.sin(phi)  # dA_l/dh, Eq11/p9; not the liquid wall perimeter.
    return (1-h/diameter)**2*gravity*(math.pi*diameter**2/4)*alpha/si


def compare(point):
    burr=assess(point)
    amin,amax=burr['alpha_envelope']
    if amin<=0.5:
        raise ValueError('This simple monotone interval proof requires alpha>0.5')
    # In this interval h and Si decrease with alpha; every factor in K grows.
    low=td1977_factor(amin,burr['diameter_m'],burr['gravity_m_s2'])
    ugmax=burr['jg_envelope_m_s'][1]/amin
    required=1+ugmax**2/low
    return dict(point_id=point['id'],burr15_margin_lower_m_s=burr['signed_margin_lower_m_s'],
                burr15_loss_entire_envelope=burr['literal_detector_loss_for_entire_envelope'],
                td1977_factor_lower_m2_s2=low,ug_upper_m_s=ugmax,
                td1977_no_loss_entire_envelope_if_rhol_over_rhog_strictly_greater_than=required,
                measured_density_ratio=None,td1977_experimental_verdict='CONDITIONAL_MISSING_DENSITY',
                observed_pattern='stratified smooth',
                lower_gas_to_liquid_phase_velocity_ratio=burr['ug_lower_m_s']/(
                    point['jl_nominal_m_s']*1.03/(1-amax)))


def result():
    return dict(schema='PRIMS-F.stratified_primary_comparison.v1',
        source=dict(title='Taitel, Lee and Dukler (1977), DOE/ID/01571--T1',
                    url='https://www.osti.gov/servlets/purl/5902200',equation=1,printed_page=6,pdf_page=10,
                    geometry='Eq11, printed p9 / PDF13: dAl/dh=Si',
                    pdf_sha256='c953722d3abb7e22369b4c89b6fd15c6a28482d896a3bbdfc8f654f41f9854fa',
                    relation_to_1976='Author report summarizes 1976 steady theory; not a claim of reading 1976 full text'),
        controls=[compare(p) for p in POINTS],
        nonstratified_record=dict(source='Balino2010 Table3/case14, p653; Fig14/p655',
                                 observed_pattern='SF',alpha_p_reported_mean=.74,
                                 jg_standard_m_s=1.86,jl_m_s=.82,
                                 status='NOT_APPLICABLE_GEOMETRY_NOT_ESTABLISHED',
                                 critical_velocity=None,
                                 reason='System flow label and mean void fraction do not establish a local flat stratified candidate; standard gas flux is not local gas velocity'),
        limits=['The comparison is conditional on density ratio; no experimental point density was invented.',
                'The author report assumes gas appreciably faster than liquid in the stratified model; the envelope ratios are reported, not a validation of that approximation.',
                'Loss of stratification is not equivalent to slug onset: the source distinguishes intermittent/annular and smooth/wavy transitions.',
                'No return law, hysteresis coefficient, angled extension or automatic switching is implemented.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    value=result();value['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    value['bounds_script_sha256']=hashlib.sha256(Path(__file__).with_name('assess_burr15_primary_bounds.py').read_bytes()).hexdigest()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


if __name__=='__main__':main()
