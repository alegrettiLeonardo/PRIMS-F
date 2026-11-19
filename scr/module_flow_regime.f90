module module_flow_regime
  !! Flow-regime utilities and transition criteria (Montgomery, 2002).
  !!
  !! The goal is to make regime/closure selection "pluggable".
  !!
  !! Modes:
  !!   MODE_STRAT   : stratified (local equilibrium) closure
  !!   MODE_DF      : drift-flux closure using actual inclination theta
  !!   MODE_DF_PIPE : drift-flux closure in pipeline (theta_eff = 0)
  !!
  !! Notes on CriteriaStabFlow:
  !!   This is a translation of the legacy MATLAB/Python implementation used
  !!   in the project. The Python draft had a typo (it compared "auxa >= 1"
  !!   instead of the final aux value). Here we apply the corrected check.
  !!
  use mod_kinds, only: wp, i4
  use mod_void_helpers, only: alpha2gamma
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  public

  integer(i4), parameter :: MODE_STRAT   = 0_i4
  integer(i4), parameter :: MODE_DF      = 1_i4
  integer(i4), parameter :: MODE_DF_PIPE = 2_i4

  integer, parameter :: BURR15_OK = 0, BURR15_INVALID_STATE = 1, BURR15_PHASE_ABSENT = 2
  integer, parameter :: BURR15_DIRECTION = 3, BURR15_ANGLE = 4, BURR15_NOT_EVALUATED = 5
  integer, parameter :: BURR15_GEOMETRY = 6, BURR15_NUMERIC = 7

  type :: burr15_result_t
    integer :: status
    real(wp) :: hL_over_D, ug2, critical_ug2, margin
    logical :: unstable
  end type

contains

  function burr15_stability(alpha, ul, ug, rhol, rhog, theta, diameter, area, g) result(test)
    !! Burr, sistema_oleoduto_riser.pdf, p. 3, Eq. (15). Inputs in SI.
    !! The liquid headroom multiplies the root in the published inequality.
    !! Use the squared comparison only after checking the velocity directions.
    !! This test describes a *stratified candidate* and its loss of stability.
    real(wp), intent(in) :: alpha, ul, ug, rhol, rhog, theta, diameter, area, g
    type(burr15_result_t) :: test
    real(wp), parameter :: pi = acos(-1.0_wp)
    real(wp) :: gamma, half_angle, headroom, area_per_d, length_factor, nan
    integer :: code

    nan = ieee_value(0.0_wp, ieee_quiet_nan)
    test = burr15_result_t(BURR15_INVALID_STATE, nan, nan, nan, nan, .false.)
    if (.not. all(ieee_is_finite([alpha, ul, ug, rhol, rhog, theta, diameter, area, g]))) return
    if (alpha < 0.0_wp .or. alpha > 1.0_wp .or. rhog <= 0.0_wp .or. rhol <= rhog .or. g <= 0.0_wp) return
    test%status = BURR15_GEOMETRY
    if (diameter <= 0.0_wp .or. area <= 0.0_wp) return
    if (diameter < 1.0_wp) then
      if (area > huge(area)*diameter) return
    endif
    area_per_d = area/diameter
    if (abs(area_per_d - (pi/4)*diameter) > 32*epsilon(area)*(pi/4)*diameter) return
    if (alpha == 0.0_wp .or. alpha == 1.0_wp) then
      test%status = BURR15_PHASE_ABSENT
      return
    endif
    test%status = BURR15_ANGLE
    if (theta > 0.0_wp .or. theta <= -pi/2) return
    test%status = BURR15_DIRECTION
    if (ug < 0.0_wp .or. ul < 0.0_wp) return

    test%status = BURR15_NUMERIC
    gamma = alpha2gamma(alpha, 1e-14_wp, code)
    if (code /= 0) return
    if (gamma <= 0.0_wp .or. gamma >= 1.0_wp) return
    half_angle = (pi/2)*gamma
    test%hL_over_D = sin(half_angle)**2
    headroom = sin((pi/2)*(1.0_wp-gamma))**2
    ! S_l = pi*D*gamma (wall wetted by liquid). A/S_l has units of length.
    length_factor = (alpha*headroom**2/(pi*gamma))*area_per_d
    length_factor = length_factor*((rhol-rhog)/rhol)*cos(theta)
    if (length_factor > 1.0_wp) then
      if (g > huge(g)/length_factor) return
    endif
    if (ug >= sqrt(huge(ug))) return
    test%ug2 = ug**2
    test%critical_ug2 = g*length_factor
    if (test%critical_ug2 <= 0.0_wp) return
    test%margin = test%ug2-test%critical_ug2
    test%unstable = test%margin >= 0.0_wp
    test%status = BURR15_OK
  end function burr15_stability

  pure real(wp) function clamp(x, lo, hi) result(y)
    real(wp), intent(in) :: x, lo, hi
    y = min(max(x, lo), hi)
  end function clamp


  pure real(wp) function alphap_from_h(h) result(alphap)
    !! GAS area fraction above a liquid level h=H_l/D (circular pipe).
    !! This function decreases from 1 to 0 as h increases from 0 to 1.
    !! alpha_p = (acos(2h-1) - (2h-1)*sqrt(1-(2h-1)^2)) / pi
    real(wp), intent(in) :: h
    real(wp) :: auxa, auxb, auxc
    auxa = 2.0_wp*h - 1.0_wp
    auxa = clamp(auxa, -1.0_wp, 1.0_wp)
    auxb = acos(auxa)
    auxc = sqrt(max(0.0_wp, 1.0_wp - auxa*auxa))
    alphap = (auxb - auxa*auxc) / acos(-1.0_wp)
  end function alphap_from_h


  real(wp) function h_from_alphap(alphap, tol) result(h)
    !! Invert alphap_from_h by bisection.
    real(wp), intent(in) :: alphap
    real(wp), intent(in) :: tol
    real(wp) :: h_lo, h_hi, h_mid, a_mid
    integer :: it

    h_lo = 0.0_wp
    h_hi = 1.0_wp
    h_mid = 0.5_wp

    h=ieee_value(0.0_wp,ieee_quiet_nan)
    if (.not. ieee_is_finite(alphap)) return
    if (alphap<0.0_wp .or. alphap>1.0_wp) return
    if (alphap == 0.0_wp) then
      h = 1.0_wp
      return
    end if
    if (alphap == 1.0_wp) then
      h = 0.0_wp
      return
    end if

    do it = 1, 80
      h_mid = 0.5_wp*(h_lo + h_hi)
      a_mid = alphap_from_h(h_mid)
      if (abs(a_mid - alphap) <= max(tol, 1.0e-14_wp)) exit
      if (a_mid > alphap) then
        h_lo = h_mid
      else
        h_hi = h_mid
      end if
      if (abs(h_hi - h_lo) <= max(tol, 1.0e-14_wp)) exit
    end do

    h = h_mid
  end function h_from_alphap

  real(wp) function criteria_stab_flow_index(alpha_p, jg, jl0, area, D, g, rhol, rhog, beta, tol) result(SI)
    !! LEGACY diagnostic only; NOT used to select the solver closure.
    !! The original normalization is dimensionally ambiguous and the supplied
    !! Burr equations disagree about the density in the threshold. Retained
    !! for comparison with historical outputs, not a qualified SI criterion.
    !!
    !! Translation of CriteriaStabFlow.m (legacy MATLAB).
    !!
    !! Returns a dimensionless metric SI:
    !!   SI >= 1  => stratified flow is unstable (switch to non-stratified closure)
    !!
    !! Inputs are DIMENSIONAL (but density ratios are scale-invariant, so rhol/rhog
    !! may also be provided nondimensional as long as both use the same scale):
    !!   alpha_p : pipeline void fraction [-]
    !!   jg      : superficial gas velocity [m/s]
    !!   jl0     : reference superficial liquid velocity at pipeline inlet [m/s]
    !!   area    : cross-sectional area [m^2]
    !!   D       : pipe diameter [m]
    !!   g       : gravity [m/s^2]
    !!   rhol,g  : densities [kg/m^3]
    !!   beta    : pipeline inclination magnitude [rad]
    !!
    real(wp), intent(in) :: alpha_p, jg, jl0, area, D, g, rhol, rhog, beta, tol

    real(wp) :: h
    real(wp) :: auxa, auxb, F
    real(wp) :: denom
    real(wp) :: sIgeom
    real(wp) :: crit

    SI = -1.0_wp

    if (.not. ( &
         ieee_is_finite(alpha_p) .and. ieee_is_finite(jg)   .and. ieee_is_finite(jl0)  .and. &
         ieee_is_finite(area)    .and. ieee_is_finite(D)    .and. ieee_is_finite(g)    .and. &
         ieee_is_finite(rhol)    .and. ieee_is_finite(rhog) .and. ieee_is_finite(beta) )) return

    if (alpha_p <= 0.0_wp .or. alpha_p >= 1.0_wp) return
    if (jg <= 0.0_wp) return
    if (jl0 <= 0.0_wp) return

    denom = max(rhol - rhog, 1.0e-12_wp)

    auxa = sqrt(max(0.0_wp, rhog / denom))
    auxb = sqrt(max(0.0_wp, D / (max(area, 1.0e-16_wp) * max(g, 1.0e-16_wp) * max(cos(beta), 1.0e-6_wp))))

    ! MATLAB: F = sqrt(rhog/(rhol-rhog))*sqrt(D/(A*g*cosB))*(QL0/A)
    F = auxa * auxb * jl0

    ! pipe liquid level from void fraction
    h = h_from_alphap(alpha_p, tol)
    h = clamp(h, 0.0_wp, 1.0_wp)

    ! non-dimensional interface perimeter term: sqrt(1-(2h-1)^2)
    sIgeom = sqrt(max(0.0_wp, 1.0_wp - (2.0_wp*h - 1.0_wp)**2))

    if (alpha_p < 1.0e-12_wp) return
    if (1.0_wp - h < 1.0e-12_wp) return

    crit = sIgeom / max(alpha_p**3, 1.0e-24_wp)
    crit = sqrt(max(0.0_wp, crit))
    crit = jg * crit / max(1.0_wp - h, 1.0e-12_wp)
    crit = crit * F

    SI = crit
  end function criteria_stab_flow_index


  pure logical function hysteresis_bool(prev_state, metric, threshold_on, threshold_off) result(state)
    !! Generic boolean hysteresis:
    !!   prev_state = .false.  -> switch ON when metric >= threshold_on
    !!   prev_state = .true.   -> switch OFF when metric <  threshold_off
    !!
    logical, intent(in) :: prev_state
    real(wp), intent(in) :: metric, threshold_on, threshold_off

    if (.not. ieee_is_finite(metric)) then
      state = prev_state
      return
    end if


    if (metric < 0.0_wp) then
      state = prev_state
      return
    end if
    if (.not. prev_state) then
      state = (metric >= threshold_on)
    else
      state = (metric >= threshold_off)
    end if
  end function hysteresis_bool


  logical function criteria_stab_flow(alpha_p, jg, jl0, area, D, g, rhol, rhog, beta, tol) result(is_unstable)
    !! Backward-compatible boolean form (SI >= 1).
    real(wp), intent(in) :: alpha_p, jg, jl0, area, D, g, rhol, rhog, beta, tol
    real(wp) :: SI
    SI = criteria_stab_flow_index(alpha_p, jg, jl0, area, D, g, rhol, rhog, beta, tol)
    is_unstable = ieee_is_finite(SI) .and. (SI >= 1.0_wp)
  end function criteria_stab_flow


  ! ---------------------------------------------------------------------------
  ! Montgomery (2002) severe-slugging correlations
  ! ---------------------------------------------------------------------------

  pure real(wp) function montgomery_eotvos(rhol, rhog, g, D, sigma) result(Eo)
    !! Eotvos number (typical definition)
    real(wp), intent(in) :: rhol, rhog, g, D, sigma
    Eo = max(0.0_wp, (rhol - rhog) * g * D * D / max(sigma, 1.0e-16_wp))
  end function montgomery_eotvos


  pure real(wp) function montgomery_eg_tb(Eo) result(Eg)
    !! Eq. (5.39) as reconstructed from the thesis layout:
    !!   Eg_TB = 0.59 + 0.3031 * (1 - exp(-Eo/90))^0.2308
    !! Ensures Eg_TB -> 0.8931 for large Eo.
    real(wp), intent(in) :: Eo
    real(wp) :: term
    term = 1.0_wp - exp(-max(Eo, 0.0_wp) / 90.0_wp)
    term = max(0.0_wp, min(1.0_wp, term))
    Eg = 0.59_wp + 0.3031_wp * term**0.2308_wp
  end function montgomery_eg_tb


  pure real(wp) function montgomery_c0(theta, jl, D, g) result(C0)
    !! Eq. (5.43) with Eq. (5.45)
    !!   If Fr_L < 3.5:
    !!     C0 = C0H + (C0V - C0H) * sin(theta)^2
    !!   Else:
    !!     C0 = C0V
    !! with C0H=1.05, C0V=1.2
    real(wp), intent(in) :: theta, jl, D, g
    real(wp) :: FrL
    real(wp), parameter :: C0H=1.05_wp, C0V=1.2_wp

    FrL = jl / sqrt(max(g*D, 1.0e-16_wp))
    if (FrL < 3.5_wp) then
      C0 = C0H + (C0V - C0H) * (sin(theta)**2)
    else
      C0 = C0V
    end if
  end function montgomery_c0


  pure real(wp) function montgomery_ud(theta, D, g) result(Ud)
    !! Eq. (5.46)-(5.48)
    !!   Ud = UdH cos(theta) + UdV sin(theta)
    !!   UdH = 0.54 sqrt(g D)
    !!   UdV = 0.35 sqrt(g D)
    real(wp), intent(in) :: theta, D, g
    real(wp) :: UdH, UdV
    UdH = 0.54_wp * sqrt(max(g*D, 1.0e-16_wp))
    UdV = 0.35_wp * sqrt(max(g*D, 1.0e-16_wp))
    Ud  = UdH*cos(theta) + UdV*sin(theta)
  end function montgomery_ud


  pure real(wp) function montgomery_ugsb(eg, c0, ud, ulsb) result(ugsb)
    !! Eq. (5.36) derived from drift-flux definition:
    !!   U_GSB = eg*(c0*U_LSB + Ud) / (1 - c0*eg)
    real(wp), intent(in) :: eg, c0, ud, ulsb
    real(wp) :: denom
    denom = 1.0_wp - c0*eg
    denom = sign(max(abs(denom), 1.0e-10_wp), denom)
    ugsb = eg*(c0*ulsb + ud) / denom
  end function montgomery_ugsb



pure real(wp) function montgomery_ugsb_crit(theta, jl, rhol, rhog, D, g, sigma, EgP) result(ugsb)
  !! Convenience wrapper:
  !!   Eo -> Eg_TB -> (optionally limited by pipeline void fraction, Eq. 5.41) -> C0, Ud -> U_GSB
  !!
  !! Montgomery (2002):
  !!   - Eq. (5.36): U_GSB = Eg*(C0*jl + Ud)/(1 - C0*Eg)
  !!   - Eq. (5.41): Eg behind the bubble nose at the riser base is bounded by the
  !!                pipeline void fraction: Eg = min(Eg_P, Eg_TB)
  !!
  !! Notes:
  !!   - For numerical robustness, we also enforce the "meaningful results" constraint
  !!     Eg < 1/C0 (Eq. 5.42) with a small margin.
  real(wp), intent(in) :: theta, jl, rhol, rhog, D, g, sigma
  real(wp), intent(in), optional :: EgP
  real(wp) :: Eo, EgTB, Eg, C0, Ud, Eg_lim

  Eo   = montgomery_eotvos(rhol, rhog, g, D, sigma)
  EgTB = montgomery_eg_tb(Eo)

  Eg = EgTB
  if (present(EgP)) then
    Eg = min(EgTB, max(0.0_wp, min(1.0_wp, EgP)))
  end if

  C0 = montgomery_c0(theta, jl, D, g)

  ! Meaningful results require Eg < 1/C0 (Montgomery Eq. 5.42)
  Eg_lim = 1.0_wp / max(C0, 1.0e-12_wp) - 1.0e-6_wp
  Eg = max(0.0_wp, min(Eg, max(0.0_wp, Eg_lim)))

  Ud   = montgomery_ud(theta, D, g)
  ugsb = montgomery_ugsb(Eg, C0, Ud, jl)
end function montgomery_ugsb_crit


end module module_flow_regime
