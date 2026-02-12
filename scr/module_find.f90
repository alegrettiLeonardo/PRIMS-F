module module_find
  !! Fortran translation of find.py
  !!
  !! Core responsibility: invert conservative variables (u1,u2,u3)
  !! to primitive variables (alpha, rhol, rhog, P, ul, ug).
  !!
  !! Notes:
  !! - Pressure is obtained from the quadratic relation in pressure_fun_uj.
  !! - Void fraction is obtained from u2 = alpha * rhog.
  !! - Velocities are obtained by solving either:
  !!     (a) drift-flux closure (theta > 0, non-stratified)
  !!     (b) local equilibrium (theta <= 0, stratified)
  !!   using Brent's method from mod_rootfinding.
  !! - Optional warm-start seeds (ul_prev/ug_prev) are used in drift-flux mode.
  !!
  !! IMPORTANT Fortran rule:
  !!   If a procedure has internal procedures, all executable statements
  !!   must appear BEFORE the CONTAINS statement.

  use mod_kinds,          only: wp, i4
  use mod_rootfinding,    only: brent_root
  use module_closure_law, only: drift_flux_swananda_ndim, RelEquilLocalPipe_comp
  use module_flow_regime, only: criteria_stab_flow_index, hysteresis_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  ! ---------------------------------------------------------------------------
  ! Pipeline stability switch (Stratified <-> Drift-Flux) with hysteresis
  !
  ! SI hysteresis (Schmitt trigger):
  !   - switch ON  (stratified -> drift-flux) when SI >= si_on
  !   - switch OFF (drift-flux -> stratified) when SI <  si_off
  ! Defaults follow the MATLAB example: SI_ON=1.0, SI_OFF=0.9.
  ! ---------------------------------------------------------------------------
  real(wp), save :: si_on  = 1.0_wp
  real(wp), save :: si_off = 0.9_wp
  logical,  save :: use_stab_hysteresis = .true.

  public :: set_stability_params
  private

  public :: compute_primitive_variables
  public :: alpha_fun_uj
  public :: pressure_fun_uj

  integer(i4), parameter :: max_roots_default = 16
  integer(i4), parameter :: nscan_default     = 200

contains

  subroutine set_stability_params(si_on_in, si_off_in, enable_hysteresis)
    !! Set pipeline stability switch hysteresis parameters.
    !!
    !! If you do not call this routine, defaults are used:
    !!   si_on  = 1.0
    !!   si_off = 0.9
    !!   use_stab_hysteresis = .true.
    real(wp), intent(in), optional :: si_on_in, si_off_in
    logical,  intent(in), optional :: enable_hysteresis

    if (present(si_on_in)) then
      if (ieee_is_finite(si_on_in) .and. si_on_in > 0.0_wp) si_on = si_on_in
    end if
    if (present(si_off_in)) then
      if (ieee_is_finite(si_off_in) .and. si_off_in > 0.0_wp) si_off = si_off_in
    end if

    ! Keep consistency: off threshold must be <= on threshold.
    if (si_off > si_on) si_off = 0.9_wp * si_on

    if (present(enable_hysteresis)) use_stab_hysteresis = enable_hysteresis
  end subroutine set_stability_params

  pure real(wp) function pressure_fun_uj(u1, u2, Cg, Cl, rho_l0, P_l0) result(P)
    real(wp), intent(in) :: u1, u2, Cg, Cl, rho_l0, P_l0
    real(wp) :: auxa, auxb, auxc, delta, sqrt_delta
    real(wp) :: P1, P2

    auxa = (Cg/Cl)**2
    auxb = rho_l0 - P_l0/(Cl*Cl)
    auxc = auxb - u1 - u2*auxa

    delta = auxc*auxc + 4.0_wp*u2*auxa*auxb
    if (delta < 0.0_wp) delta = 0.0_wp
    sqrt_delta = sqrt(delta)

    P1 = ((-auxc - sqrt_delta)/2.0_wp) * (Cl*Cl)
    P2 = ((-auxc + sqrt_delta)/2.0_wp) * (Cl*Cl)

    if (P1 > 0.0_wp) then
      P = P1
    else
      P = P2
    end if
  end function pressure_fun_uj


  pure real(wp) function alpha_fun_uj(u1, u2, Cg, Cl, rho_l0, P_l0, tol) result(alpha)
    real(wp), intent(in) :: u1, u2, Cg, Cl, rho_l0, P_l0, tol
    real(wp) :: P, rhog

    if (abs(u2) < tol) then
      alpha = 0.0_wp
      return
    end if

    P = pressure_fun_uj(u1, u2, Cg, Cl, rho_l0, P_l0)
    rhog = P/(Cg*Cg)

    if (rhog < tol) then
      alpha = 0.0_wp
      return
    end if

    alpha = u2 / rhog
    if (alpha < 0.0_wp) alpha = 0.0_wp
    if (alpha > 1.0_wp) alpha = 1.0_wp
  end function alpha_fun_uj


  subroutine compute_primitive_variables(u1, u2, u3, theta, Cg, Cl, rho_l0, P_l0, D, AREA, EPS, G, MUL, MUG, sigma, njl_ref, &
                                         w_u, w_rho, tol, tola, maxit, alpha, rhol, rhog, P, ul, ug, index, idx_prev, &
                                         ul_prev, ug_prev)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: theta
    real(wp), intent(in) :: Cg, Cl, rho_l0, P_l0
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma, njl_ref
    real(wp), intent(in) :: w_u, w_rho
    real(wp), intent(in) :: tol, tola
    integer(i4), intent(in) :: maxit

    real(wp), intent(out) :: alpha, rhol, rhog, P, ul, ug
    integer(i4), intent(out) :: index
    integer(i4), intent(in), optional :: idx_prev

    real(wp), intent(in), optional :: ul_prev, ug_prev
    real(wp), parameter :: alpha_eps = 1.0e-12_wp
    real(wp), parameter :: rho_eps   = 1.0e-16_wp
    real(wp), parameter :: P_eps     = 1.0e-12_wp

    alpha = alpha_fun_uj(u1, u2, Cg, Cl, rho_l0, P_l0, tol)
    if (.not. ieee_is_finite(alpha)) alpha = 0.5_wp
    alpha = max(min(alpha, 1.0_wp - alpha_eps), alpha_eps)

    P = pressure_fun_uj(u1, u2, Cg, Cl, rho_l0, P_l0)
    if (.not. ieee_is_finite(P)) P = P_eps
    P = max(P, P_eps)

    rhog = max(P/(Cg*Cg), rho_eps)
    rhol = max(rho_l0 + (P - P_l0)/(Cl*Cl), rho_eps)

    call find_ulug_from_uj_genflw_simp(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, njl_ref, &
                                      w_u, w_rho, tol, tola, maxit, ul, ug, index, idx_prev=idx_prev, ul_prev=ul_prev, ug_prev=ug_prev)
  end subroutine compute_primitive_variables


  subroutine find_ulug_from_uj_genflw_simp(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, njl_ref, &
                                          w_u, w_rho, tol, tola, maxit, ul, ug, index, idx_prev, ul_prev, ug_prev)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, theta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma, njl_ref
    real(wp), intent(in) :: w_u, w_rho
    real(wp), intent(in) :: tol, tola
    integer(i4), intent(in) :: maxit

    real(wp), intent(out) :: ul, ug
    integer(i4), intent(out) :: index
    integer(i4), intent(in), optional :: idx_prev

    real(wp), intent(in), optional :: ul_prev, ug_prev
    ! --- locals for pipeline stability criterion
    real(wp) :: ul_s, ug_s
    real(wp) :: jg_dim, jl0_dim
    real(wp) :: rhol_dim, rhog_dim
    real(wp) :: SI
    real(wp) :: beta_mag
    logical :: use_df
    if (theta > 0.0_wp) then
      ! Riser / non-stratified: drift-flux with actual inclination
      index = 1
      call find_ul_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, &
                                 w_u, w_rho, tol, tola, maxit, ul, ug, ul_prev=ul_prev, ug_prev=ug_prev)
    else
      ! Pipeline branch: start from stratified local-equilibrium and apply
      ! stability criterion (CriteriaStabFlow) to decide transition to
      ! non-stratified drift-flux in the pipeline.
      beta_mag = abs(theta)

      ! 1) Stratified solution (local equilibrium)
      call find_ul_ug_from_uj(u1, u2, u3, rhol, rhog, alpha, beta_mag, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, tola, ul_s, ug_s)

      ! 2) Evaluate stability index using dimensional superficial velocities
      jg_dim   = alpha * ug_s * w_u
      jl0_dim  = njl_ref * w_u
      rhol_dim = rhol * w_rho
      rhog_dim = rhog * w_rho

      SI = criteria_stab_flow_index(alpha, jg_dim, jl0_dim, AREA, D, G, rhol_dim, rhog_dim, beta_mag, tol)

      ! 3) Decision with optional hysteresis (avoid chattering near SI≈1)
      !    Thresholds are configurable via set_stability_params():
      !      si_on  (default 1.0) and si_off (default 0.9)
      if (present(idx_prev) .and. use_stab_hysteresis) then
        use_df = hysteresis_bool(idx_prev /= 0_i4, SI, si_on, si_off)
      else
        use_df = ieee_is_finite(SI) .and. (SI >= si_on)
      end if

      if (use_df) then
        ! Unstable stratified (or latched DF): switch to drift-flux in pipeline
        index = 1
        call find_ul_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, beta_mag, D, AREA, EPS, G, MUL, MUG, sigma, &
                                   w_u, w_rho, tol, tola, maxit, ul, ug, ul_prev=ul_prev, ug_prev=ug_prev)
      else
        ! Stable stratified
        index = 0
        ul = ul_s
        ug = ug_s
      end if
    end if
  end subroutine find_ulug_from_uj_genflw_simp


  subroutine find_ul_ug_from_uj(u1, u2, u3, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, tola, ul, ug)
    !! Stratified branch (local equilibrium): solve RelEquilLocalPipe_comp.
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, beta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG
    real(wp), intent(in) :: w_u, w_rho, tol, tola
    real(wp), intent(out) :: ul, ug

    if (abs(u2) < tola) then
      ug = 0.0_wp
      ul = u3 / u1
      return
    end if
    if (abs(u1) < tola) then
      ul = 0.0_wp
      ug = u3 / u2
      return
    end if

    if (u2 <= u1) then
      ug = find_ug_from_uj(u1, u2, u3, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol)
      ul = u3 / u1 - (u2/u1) * ug
    else
      ul = find_ul_from_uj(u1, u2, u3, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol)
      ug = u3 / u2 - (u1/u2) * ul
    end if
  end subroutine find_ul_ug_from_uj


  function find_ug_from_uj(u1, u2, u3, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol) result(ug_root)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, beta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG
    real(wp), intent(in) :: w_u, w_rho, tol

    real(wp) :: ug_root
    real(wp) :: a, b, fa, fb
    integer(i4) :: info

    a = -(u3/u2)
    b = +(u3/u2)

    fa = f(a)
    fb = f(b)

    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
      ug_root = 0.0_wp
      return
    end if
    if (fa*fb > 0.0_wp) then
      ug_root = 0.0_wp
      return
    end if

    ug_root = brent_root(f, a, b, tol, 1000_i4, info)
    if (info /= 0 .or. .not. ieee_is_finite(ug_root)) ug_root = 0.0_wp

  contains
    function f(ug) result(val)
      real(wp), intent(in) :: ug
      real(wp) :: val
      real(wp) :: ul
      ul = u3/u1 - (u2/u1)*ug
      call RelEquilLocalPipe_comp(ul, ug, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, val)
    end function f
  end function find_ug_from_uj


  function find_ul_from_uj(u1, u2, u3, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol) result(ul_root)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, beta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG
    real(wp), intent(in) :: w_u, w_rho, tol

    real(wp) :: ul_root
    real(wp) :: a, b, fa, fb
    integer(i4) :: info

    a = -(u3/u1)
    b = +(u3/u1)

    fa = f(a)
    fb = f(b)

    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
      ul_root = 0.0_wp
      return
    end if
    if (fa*fb > 0.0_wp) then
      ul_root = 0.0_wp
      return
    end if

    ul_root = brent_root(f, a, b, tol, 1000_i4, info)
    if (info /= 0 .or. .not. ieee_is_finite(ul_root)) ul_root = 0.0_wp

  contains
    function f(ul) result(val)
      real(wp), intent(in) :: ul
      real(wp) :: val
      real(wp) :: ug
      ug = u3/u2 - (u1/u2)*ul
      call RelEquilLocalPipe_comp(ul, ug, rhol, rhog, alpha, beta, D, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, val)
    end function f
  end function find_ul_from_uj


  subroutine find_ul_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, &
                                    w_u, w_rho, tol, tola, maxit, ul, ug, ul_prev, ug_prev)
    !! Non-stratified branch: solve drift_flux_swananda_ndim.
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, theta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in) :: w_u, w_rho
    real(wp), intent(in) :: tol, tola
    integer(i4), intent(in) :: maxit

    real(wp), intent(out) :: ul, ug
    real(wp), intent(in), optional :: ul_prev, ug_prev

    real(wp) :: Fr, denom, cos_t
    real(wp) :: roots(max_roots_default)
    integer(i4) :: nroots

    ! Degenerate cases
    if (abs(u2) < tola) then
      ug = 0.0_wp
      ul = u3/u1
      return
    end if
    if (abs(u1) < tola) then
      ul = 0.0_wp
      ug = u3/u2
      return
    end if

    cos_t = cos(theta)
    denom = (rhol - rhog) * cos_t
    denom = denom / max(rhog, 1.0e-30_wp)
    denom = max(denom, 1.0e-30_wp)

    Fr = w_u / sqrt(max(D*G, 1.0e-30_wp))
    Fr = Fr / sqrt(denom)

    if (u2 <= u1) then
      if (Fr > 0.1_wp) then
        ug = find_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, ug_prev)
        ul = u3/u1 - (u2/u1)*ug
      else
        call find_all_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, roots, nroots)
        if (nroots > 0) then
          ug = roots(1)
          ul = u3/u1 - (u2/u1)*ug
        else
          call find_all_ul_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, roots, nroots)
          if (nroots > 0) then
            ul = roots(1)
            ug = u3/u2 - (u1/u2)*ul
          else
            ug = 0.0_wp
            ul = u3/u1
          end if
        end if
      end if
    else
      if (Fr > 0.1_wp) then
        ul = find_ul_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, ul_prev)
        ug = u3/u2 - (u1/u2)*ul
      else
        call find_all_ul_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, roots, nroots)
        if (nroots > 0) then
          ul = roots(1)
          ug = u3/u2 - (u1/u2)*ul
        else
          call find_all_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, roots, nroots)
          if (nroots > 0) then
            ug = roots(1)
            ul = u3/u1 - (u2/u1)*ug
          else
            ug = 0.0_wp
            ul = u3/u1
          end if
        end if
      end if
    end if
  end subroutine find_ul_ug_from_uj_dfr


  function find_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, ug_prev) result(ug_root)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, theta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in) :: w_u, w_rho, tol
    integer(i4), intent(in) :: maxit
    real(wp), intent(in), optional :: ug_prev

    real(wp) :: ug_root
    real(wp) :: ug_min, ug_max
    real(wp) :: a, b, fa, fb, span, ugp
    integer(i4) :: info

    ug_min = 0.0_wp
    ug_max = u3/u2
    if (ug_max < ug_min) ug_max = ug_min

    ! Warm-start bracket
    if (present(ug_prev)) then
      ugp = ug_prev
      if (ieee_is_finite(ugp) .and. ugp > ug_min .and. ugp < ug_max) then
        span = 0.25_wp * max(ugp, 1.0e-12_wp)
        a = max(ug_min, ugp - span)
        b = min(ug_max, ugp + span)
        fa = f(a)
        fb = f(b)
        if (ieee_is_finite(fa) .and. ieee_is_finite(fb) .and. fa*fb <= 0.0_wp) then
          ug_root = brent_root(f, a, b, tol, maxit, info)
          if (info == 0 .and. ieee_is_finite(ug_root)) return
        end if
      end if
    end if

    fa = f(ug_min)
    fb = f(ug_max)
    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
      ug_root = 0.0_wp
      return
    end if
    if (fa*fb > 0.0_wp) then
      ug_root = 0.0_wp
      return
    end if

    ug_root = brent_root(f, ug_min, ug_max, tol, maxit, info)
    if (info /= 0 .or. .not. ieee_is_finite(ug_root)) ug_root = 0.0_wp

  contains
    function f(ug) result(val)
      real(wp), intent(in) :: ug
      real(wp) :: val
      real(wp) :: ul
      ul = u3/u1 - (u2/u1)*ug
      call drift_flux_swananda_ndim(ul, ug, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, val)
    end function f
  end function find_ug_from_uj_dfr


  function find_ul_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, maxit, ul_prev) result(ul_root)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, theta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in) :: w_u, w_rho, tol
    integer(i4), intent(in) :: maxit
    real(wp), intent(in), optional :: ul_prev

    real(wp) :: ul_root
    real(wp) :: ul_min, ul_max
    real(wp) :: a, b, fa, fb, span, ulp
    integer(i4) :: info
    ul_min = min(0.0_wp, u3/u1)
    ul_max = max(0.0_wp, u3/u1)
    if (abs(ul_max-ul_min) < 1.0e-12_wp) then
      ul_min = ul_min - 1.0_wp
      ul_max = ul_max + 1.0_wp
    end if
    ! Warm-start bracket
    if (present(ul_prev)) then
      ulp = ul_prev
      if (ieee_is_finite(ulp) .and. ulp > ul_min .and. ulp < ul_max) then
        span = 0.25_wp * max(ulp, 1.0e-12_wp)
        a = max(ul_min, ulp - span)
        b = min(ul_max, ulp + span)
        fa = f(a)
        fb = f(b)
        if (ieee_is_finite(fa) .and. ieee_is_finite(fb) .and. fa*fb <= 0.0_wp) then
          ul_root = brent_root(f, a, b, tol, maxit, info)
          if (info == 0 .and. ieee_is_finite(ul_root)) return
        end if
      end if
    end if

    fa = f(ul_min)
    fb = f(ul_max)
    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb)) then
      ul_root = 0.0_wp
      return
    end if
    if (fa*fb > 0.0_wp) then
      ul_root = 0.0_wp
      return
    end if

    ul_root = brent_root(f, ul_min, ul_max, tol, maxit, info)
    if (info /= 0 .or. .not. ieee_is_finite(ul_root)) ul_root = 0.0_wp

  contains
    function f(ul) result(val)
      real(wp), intent(in) :: ul
      real(wp) :: val
      real(wp) :: ug
      ug = u3/u2 - (u1/u2)*ul
      call drift_flux_swananda_ndim(ul, ug, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, val)
    end function f
  end function find_ul_from_uj_dfr


  subroutine find_all_ug_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, &
                                     w_u, w_rho, tol, maxit, roots, nroots)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, theta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in) :: w_u, w_rho, tol
    integer(i4), intent(in) :: maxit

    real(wp), intent(out) :: roots(max_roots_default)
    integer(i4), intent(out) :: nroots

    real(wp) :: ugmin, ugmax, dx, a, b, fa, fb, root
    integer(i4) :: j, info

    roots = 0.0_wp
    nroots = 0

    ugmin = 0.0_wp
    ugmax = u3/u2
    if (ugmax <= ugmin) return
    dx = (ugmax-ugmin)/real(nscan_default, wp)

    a = ugmin
    fa = f(a)

    do j = 1, nscan_default
      b = ugmin + real(j, wp)*dx
      fb = f(b)
      if (ieee_is_finite(fa) .and. ieee_is_finite(fb)) then
        if (fa*fb <= 0.0_wp) then
          root = brent_root(f, a, b, tol, maxit, info)
          if (info == 0 .and. ieee_is_finite(root)) then
            call push_unique_root(root, tol, roots, nroots)
            if (nroots >= max_roots_default) exit
          end if
        end if
      end if
      a = b
      fa = fb
    end do

  contains
    function f(ug) result(val)
      real(wp), intent(in) :: ug
      real(wp) :: val
      real(wp) :: ul
      ul = u3/u1 - (u2/u1)*ug
      call drift_flux_swananda_ndim(ul, ug, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, val)
    end function f
  end subroutine find_all_ug_from_uj_dfr


  subroutine find_all_ul_from_uj_dfr(u1, u2, u3, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, &
                                     w_u, w_rho, tol, maxit, roots, nroots)
    real(wp), intent(in) :: u1, u2, u3
    real(wp), intent(in) :: rhol, rhog, alpha, theta
    real(wp), intent(in) :: D, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in) :: w_u, w_rho, tol
    integer(i4), intent(in) :: maxit

    real(wp), intent(out) :: roots(max_roots_default)
    integer(i4), intent(out) :: nroots

    real(wp) :: ulmin, ulmax, dx, a, b, fa, fb, root
    integer(i4) :: j, info

    roots = 0.0_wp
    nroots = 0

    ulmin = 0.0_wp
    ulmax = u3/u1
    if (ulmax <= ulmin) return
    dx = (ulmax-ulmin)/real(nscan_default, wp)

    a = ulmin
    fa = f(a)

    do j = 1, nscan_default
      b = ulmin + real(j, wp)*dx
      fb = f(b)
      if (ieee_is_finite(fa) .and. ieee_is_finite(fb)) then
        if (fa*fb <= 0.0_wp) then
          root = brent_root(f, a, b, tol, maxit, info)
          if (info == 0 .and. ieee_is_finite(root)) then
            call push_unique_root(root, tol, roots, nroots)
            if (nroots >= max_roots_default) exit
          end if
        end if
      end if
      a = b
      fa = fb
    end do

  contains
    function f(ul) result(val)
      real(wp), intent(in) :: ul
      real(wp) :: val
      real(wp) :: ug
      ug = u3/u2 - (u1/u2)*ul
      call drift_flux_swananda_ndim(ul, ug, rhol, rhog, alpha, theta, D, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, val)
    end function f
  end subroutine find_all_ul_from_uj_dfr


  subroutine push_unique_root(x, tol, roots, nroots)
    real(wp), intent(in) :: x, tol
    real(wp), intent(inout) :: roots(max_roots_default)
    integer(i4), intent(inout) :: nroots
    integer(i4) :: k
    real(wp) :: thr

    thr = 10.0_wp*tol
    do k = 1, nroots
      if (abs(x - roots(k)) <= thr) return
    end do
    if (nroots < max_roots_default) then
      nroots = nroots + 1
      roots(nroots) = x
    end if
  end subroutine push_unique_root

end module module_find
