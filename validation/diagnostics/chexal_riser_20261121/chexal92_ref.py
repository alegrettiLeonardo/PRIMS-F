"""Independent transcription of Chexal et al. (1992), air-water, COCURRENT and
vertical-countercurrent-HIGH coefficients only (no root search policy, no
solver code imported). Used to (a) cross-check the Fortran core at PRIMS-F
states and (b) recompute the closure residual from accepted solver fields.

Angle: PRIMS-F theta in radians from the horizontal (riser theta >= 0);
the source angle is degrees from the vertical, theta_C = 90 - theta*180/pi.
Equations: 3-2/3-3 (Fr), 3-4a/b, 3-7/3-8 (L air-water), 3-11..3-17, 3-18/3-19
(horizontal), 3-21/3-22, 3-23..3-34, 3-35..3-38 (downflow / CCFL C3).
Countercurrent (jl<0<=jg) is admitted ONLY at theta_C == 0 (p.262); inclined
countercurrent raises DomainError, as does horizontal flow with negative fluxes.
"""
import math


class DomainError(ValueError):
    pass


def _component(a, rl, rg, d, g, ml, mg, sig, jl, jg, horizontal):
    rel, reg = rl * jl * d / ml, rg * jg * d / mg
    re = reg if (reg > rel or reg < 0) else rel                        # 3-15
    b1 = min(0.8, 1.0 / (1.0 + math.exp(-re / 60000.0)))               # 3-13/14
    k0 = b1 + (1 - b1) * (rg / rl) ** 0.25                             # 3-11
    r = (1 + 1.57 * rg / rl) / (1 - b1)                                # 3-12
    if horizontal:
        L = min(1.125 * a ** 0.6, 1.0)                                 # 3-19
    elif reg >= 0:
        L = min(1.15 * a ** 0.45, 1.0)                                 # 3-7
    else:
        L = min(1.05 * a ** 0.25, 1.0)                                 # 3-8
    c0 = L / (k0 + (1 - k0) * a ** r)                                  # 3-4a
    if rl / rg <= 18:
        c2 = 0.4757 * math.log(rl / rg) ** 0.7                         # 3-27
    else:
        c5 = math.sqrt(150 * rg / rl)                                  # 3-29
        c2 = 1.0 if c5 >= 1 else 1 / (1 - math.exp(-c5 / (1 - c5)))   # 3-28
    c7 = (0.09144 / d) ** 0.6                                          # 3-31/33
    c4 = 1.0 if c7 >= 1 else 1 / (1 - math.exp(-c7 / (1 - c7)))        # 3-30/32
    if jl >= 0 and jg >= 0:
        c3 = max(0.5, 2 * math.exp(-abs(rel) / 60000.0))               # 3-34
    else:
        q = abs(rel)                                                   # 3-35..38 (C10 per RELAP5/MOD3 1995 p.6-10)
        c10 = 2 * math.exp((q / 350000.0) ** 0.4) - 1.75 * q ** 0.03 * math.exp(-q / 50000.0 * (0.0381 / d) ** 2) \
            + (0.0381 / d) ** 0.25 * q ** 0.001
        b2 = (1 / (1 + 0.05 * q / 350000.0)) ** 0.4
        c3 = 2 * (c10 / 2) ** b2
    v0 = 1.41 * ((rl - rg) * sig * g / rl ** 2) ** 0.25 * c2 * c3 * c4  # 3-24
    if reg >= 0:
        c9 = (1 - a) ** b1                                             # 3-25
    else:
        c9 = min(0.7, (1 - a) ** 0.65)                                 # 3-26
        c0 = max(c0, v0 * (1 - a) ** 0.2 / (abs(jl) + abs(jg)))        # 3-4b
    if horizontal:
        c0 = (1 + a ** 0.05 * (1 - a) ** 2) * c0                       # 3-18
    return c0, v0 * c9


def coefficients(a, rl, rg, theta, jl, jg, d=0.05254, g=9.81, ml=0.89e-3, mg=1.85e-5, sig=0.0728):
    """(C0, Vgj) in the PRIMS-F axis (Vgj sign follows the axis orientation)."""
    direction = -1.0 if theta < 0 else 1.0
    fl, fg = direction * jl, direction * jg
    tc = 90.0 - abs(theta) * 180.0 / math.pi
    if abs(abs(theta) - math.pi / 2) <= 32 * 2.220446049250313e-16:
        tc = 0.0
    if theta == 0.0 and (fl < 0 or fg < 0):
        raise DomainError('horizontal flow with a negative flux (pp.264-265)')
    if fg < 0 < fl:
        raise DomainError('gas down / liquid up')
    cc = fl < 0 <= fg
    if cc and tc != 0.0:
        raise DomainError('inclined countercurrent (p.262: vertical only)')
    cv, vv = _component(a, rl, rg, d, g, ml, mg, sig, fl, fg, False)
    if tc == 0.0:
        return cv, direction * vv
    ch, vh = _component(a, rl, rg, d, g, ml, mg, sig, abs(fl), abs(fg), True)
    if fg >= 0:
        fr = (90 - tc) / 90                                            # 3-3a
        v = fr * vv + (1 - fr) * vh                                    # 3-21
    else:
        fr = 1.0 if tc <= 80 else (90 - tc) / 10                       # 3-3b
        v = fr * vv + (fr - 1) * vh                                    # 3-22
    return fr * cv + (1 - fr) * ch, direction * v                     # 3-2


def residual(a, rl, rg, theta, jl, jg, **kw):
    """R = alpha (C0 (jl+jg) + Vgj) - jg  [m/s]  (Eq.3-1 identity, never divided by j)."""
    c0, v = coefficients(a, rl, rg, theta, jl, jg, **kw)
    return a * (c0 * (jl + jg) + v) - jg
