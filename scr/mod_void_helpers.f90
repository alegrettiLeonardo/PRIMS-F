module mod_void_helpers
  !! Utilities originally from voidFraction.py:
  !! - numerical guards (clamp/safe_div)
  !! - Fanning friction factor ffan and its partial derivatives
  !! - alpha <-> gamma relation via bisection
  use mod_kinds,       only: wp
  use mod_rootfinding, only: bisect_root
  implicit none
  private

  public :: clamp, safe_div
  public :: calcular_fatrito, ffan, Dffan
  public :: fvgamma, alpha2gamma

  real(wp), parameter :: PI = acos(-1.0_wp)

contains

  pure real(wp) function clamp(x, lo, hi) result(y)
    real(wp), intent(in) :: x, lo, hi
    y = min(max(x, lo), hi)
  end function clamp

  pure real(wp) function safe_div(a, b, eps) result(y)
    real(wp), intent(in) :: a, b
    real(wp), intent(in), optional :: eps
    real(wp) :: e, bb

    if (present(eps)) then
      e = eps
    else
      e = 1.0e-30_wp
    end if

    bb = b
    if (abs(bb) <= e) bb = sign(e, bb)
    y = a / bb
  end function safe_div


  ! =========================
  ! Fator de atrito (Fanning)
  ! =========================

  real(wp) function calcular_fatrito(Re, epsD) result(f)
    !! Explicit approximation (same algebra as in Python).
    real(wp), intent(in) :: Re, epsD
    real(wp) :: Re2, termo_interno, f0, f1

    Re2 = max(Re, 1.0e-12_wp)

    termo_interno = 5.8506_wp / (Re2 ** 0.8981_wp) + (epsD ** 1.1098_wp) / 2.8257_wp
    termo_interno = max(termo_interno, 1.0e-30_wp)

    f0 = log10(termo_interno)
    f1 = (epsD / 3.7065_wp) - (5.0452_wp / Re2) * f0
    f1 = max(f1, 1.0e-30_wp)

    f = (-4.0_wp * log10(f1)) ** (-2.0_wp)
  end function calcular_fatrito


  real(wp) function ffan(epsD, Re) result(f)
    !! Fanning friction factor piecewise:
    !!   - Laminar (Re < 2000): f = 16/Re
    !!   - Transitional (2000..2300): linear blend
    !!   - Turbulent (Re > 2300): calcular_fatrito(Re, epsD)
    real(wp), intent(in) :: epsD, Re
    real(wp) :: Re2, fl, ft

    Re2 = max(Re, 1.0e-12_wp)

    if (Re2 < 2000.0_wp) then
      f = 16.0_wp / Re2
    elseif (Re2 > 2300.0_wp) then
      f = calcular_fatrito(Re2, epsD)
    else
      fl = 16.0_wp / 2000.0_wp
      ft = calcular_fatrito(2300.0_wp, epsD)
      f = ((Re2 - 2000.0_wp) * ft + (2300.0_wp - Re2) * fl) / 300.0_wp
    end if
  end function ffan


  subroutine Dffan(Rem, epsD, D1ffan, D2ffan)
    !! Partial derivatives of ffan:
    !!   D1ffan = df/dRem
    !!   D2ffan = df/depsD
    real(wp), intent(in)  :: Rem, epsD
    real(wp), intent(out) :: D1ffan, D2ffan

    real(wp) :: Re2, ED
    real(wp) :: a, b, c, d, e
    real(wp) :: ln10
    real(wp) :: inner, L, Aval, logA10
    real(wp) :: df_dA
    real(wp) :: dL_dRem, dA_dRem
    real(wp) :: d_inner_dED, dL_dED, dA_dED

    Re2 = max(Rem, 1.0e-12_wp)
    ED  = epsD

    if (Re2 < 2000.0_wp) then
      D1ffan = -16.0_wp / (Re2 * Re2)
      D2ffan = 0.0_wp
      return
    end if

    if (Re2 > 2300.0_wp) then
      a = ED / 3.7065_wp
      b = 5.0452_wp
      c = (ED ** 1.1098_wp) / 2.8257_wp
      d = 5.8506_wp
      e = 0.8981_wp
      ln10 = log(10.0_wp)

      inner = c + d * Re2 ** (-e)
      inner = max(inner, 1.0e-30_wp)
      L = log10(inner)

      Aval = a - (b / Re2) * L
      Aval = max(Aval, 1.0e-30_wp)

      logA10 = log10(Aval)

      ! df/dA  where f = 1/[16 (log10 Aval)^2]
      df_dA = -2.0_wp * (1.0_wp / 16.0_wp) * (logA10 ** (-3.0_wp)) * (1.0_wp / (Aval * ln10))

      ! dA/dRe
      dL_dRem = (1.0_wp / (inner * ln10)) * (d * (-e) * Re2 ** (-e - 1.0_wp))
      dA_dRem = - ( -b / (Re2 ** 2) * L + (b / Re2) * dL_dRem )

      ! dA/dED
      d_inner_dED = (1.1098_wp * (ED ** 0.1098_wp)) / 2.8257_wp
      dL_dED = (1.0_wp / (inner * ln10)) * d_inner_dED
      dA_dED = (1.0_wp / 3.7065_wp) - (b / Re2) * dL_dED

      D1ffan = df_dA * dA_dRem
      D2ffan = df_dA * dA_dED
      return
    end if

    ! Transition region: linear blend.
    ! ft evaluated at Re=2300 -> treated constant w.r.t Re in this branch.
    D1ffan = (calcular_fatrito(2300.0_wp, ED) - 16.0_wp/2000.0_wp) / 300.0_wp
    D2ffan = 0.0_wp
  end subroutine Dffan


  ! =========================
  ! Relação alpha <-> gamma
  ! =========================

  real(wp) function fvgamma(alpha, gamma) result(val)
    real(wp), intent(in) :: alpha, gamma
    val = alpha - 1.0_wp + gamma - sin(2.0_wp*PI*gamma)/(2.0_wp*PI)
  end function fvgamma


  real(wp) function alpha2gamma(alpha, tol) result(gamma)
    !! Solve gamma in [0,1] given alpha.
    !! Uses bisection on fvgamma(alpha, g) == 0.
    real(wp), intent(in) :: alpha, tol

    real(wp) :: root
    integer :: ier

    if (abs(alpha) <= tol) then
      gamma = 1.0_wp
      return
    end if
    if (abs(1.0_wp - alpha) <= tol) then
      gamma = 0.0_wp
      return
    end if

    ! NOTE: bisect_root is a FUNCTION (returns the root) with an integer INFO output.
    root = bisect_root(fun_gamma, 0.0_wp, 1.0_wp, tol, 200, ier)
    if (ier /= 0) then
      gamma = 0.5_wp
    else
      gamma = root
    end if

  contains
    real(wp) function fun_gamma(g) result(val)
      real(wp), intent(in) :: g
      val = fvgamma(alpha, g)
    end function fun_gamma
  end function alpha2gamma

end module mod_void_helpers
