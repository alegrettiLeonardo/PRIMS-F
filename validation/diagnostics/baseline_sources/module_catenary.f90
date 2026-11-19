module module_catenary
  !! Fortran translation of catenary.py
  !!
  !! Replaces scipy.optimize.fsolve with a damped Newton method and a robust
  !! bisection fallback (using mod_rootfinding).
  use mod_kinds,       only: wp
  use mod_rootfinding, only: bisect_root
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private

  public :: eq_transcendental, solve_A
  public :: altura_catenaria, comprimento_arco, angulo_inclinacao
  public :: catenaria_geometria
  public :: fun_or_geo, catenary_constant
  public :: fungeo, Dfungeo

contains

  pure real(wp) function eq_transcendental(A, Z, X) result(val)
    !! Z - A*(cosh(X/A) - 1)
    real(wp), intent(in) :: A, Z, X
    val = Z - A * (cosh(X / A) - 1.0_wp)
  end function eq_transcendental


  real(wp) function solve_A(Z, X) result(A)
    !! Solve A>0 such that Z - A*(cosh(X/A)-1) = 0.
    !! Faithful to Python's behavior:
    !!   - Uses a physically informed initial guess
    !!   - Tight tolerance (1e-10)
    !!   - Raises error if no convergence (here: returns A=0 and stops via info)
    real(wp), intent(in) :: Z, X

    integer :: it
    real(wp) :: A_guess, f, df, y, c, s
    real(wp) :: tol
    real(wp) :: damp
    real(wp) :: Anew
    real(wp) :: root
    integer :: ier

    tol = 1.0e-10_wp

    if (Z <= 0.0_wp .or. X <= 0.0_wp) then
      error stop 'solve_A: Z e X devem ser positivos.'
    end if

    ! Initial guess (same as Python)
    A_guess = X / log(Z / X + 1.0_wp)
    if (A_guess <= 0.0_wp .or. ieee_is_nan(A_guess)) then
      A_guess = max(Z, X)
    end if

    A = A_guess
    damp = 1.0_wp

    do it = 1, 100
      f = eq_transcendental(A, Z, X)
      if (abs(f) <= tol) return

      y = X / A
      c = cosh(y) - 1.0_wp
      s = sinh(y)

      ! df/dA = -(cosh(y)-1) + y*sinh(y)
      df = -c + y*s
      if (abs(df) <= 1.0e-30_wp .or. ieee_is_nan(df)) exit

      Anew = A - damp * f / df

      ! Keep positive; if it goes bad, reduce damping
      if (Anew <= 0.0_wp .or. ieee_is_nan(Anew)) then
        damp = 0.5_wp * damp
        if (damp < 1.0e-6_wp) exit
        cycle
      end if

      ! accept step; if residual worsened, damp
      if (abs(eq_transcendental(Anew, Z, X)) > abs(f)) then
        damp = 0.5_wp * damp
        if (damp < 1.0e-6_wp) exit
      else
        A = Anew
        damp = min(1.0_wp, 2.0_wp*damp)
      end if
    end do

    ! Robust fallback: bisection on [A_lo, A_hi]
    root = bisect_root(funA, 1.0e-12_wp, max(Z,X)*1.0e3_wp, tol, 5000, ier)
    if (ier /= 0) then
      error stop 'solve_A: Falha ao resolver A (Newton+bisection).'
    end if
    A = root

  contains
    real(wp) function funA(a0) result(v)
      real(wp), intent(in) :: a0
      v = eq_transcendental(a0, Z, X)
    end function funA
  end function solve_A


  pure real(wp) function altura_catenaria(x, A) result(zx)
    real(wp), intent(in) :: x, A
    zx = A * (cosh(x / A) - 1.0_wp)
  end function altura_catenaria


  pure real(wp) function comprimento_arco(x, A) result(sx)
    real(wp), intent(in) :: x, A
    sx = A * sinh(x / A)
  end function comprimento_arco


  pure real(wp) function angulo_inclinacao(s, A) result(theta)
    !! theta(s) = 0.5*log((1+s/A)/(1-s/A)), returns +Inf if |s/A|>=1
    real(wp), intent(in) :: s, A
    real(wp) :: r
    r = s / A
    if (abs(r) >= 1.0_wp) then
      theta = huge(1.0_wp)
    else
      theta = 0.5_wp * log((1.0_wp + r) / (1.0_wp - r))
    end if
  end function angulo_inclinacao


  subroutine catenaria_geometria(Z, X, x_pos, z_x, s_x, theta_x, A)
    !! Equivalent to catenaria_geometria() in Python.
    real(wp), intent(in)  :: Z, X, x_pos
    real(wp), intent(out) :: z_x, s_x, theta_x
    real(wp), intent(out) :: A

    A = solve_A(Z, X)
    z_x = altura_catenaria(x_pos, A)
    s_x = comprimento_arco(x_pos, A)
    theta_x = angulo_inclinacao(s_x, A)
  end subroutine catenaria_geometria


  pure real(wp) function fun_or_geo(s, Lp, beta, AC) result(theta)
    !! Local inclination angle in a pipeline–riser system (catenary).
    !! Convention:
    !!   - Pipe (s < Lp): theta = -beta
    !!   - Riser (s >= Lp): theta = atan( sinh(asinh((s-Lp)/A)) )
    real(wp), intent(in) :: s, Lp, beta, AC
    real(wp) :: xa

    if (s < Lp) then
      theta = -beta
    else
      xa = asinh((s - Lp) / AC)
      theta = atan(sinh(xa))
    end if
  end function fun_or_geo


  real(wp) function catenary_constant(X, Z) result(A)
    !! Equivalent to catenary_constant() in Python (uses fsolve there).
    real(wp), intent(in) :: X, Z
    A = solve_A(Z, X)
  end function catenary_constant


  pure real(wp) function fungeo(s, VPARAM) result(theta)
    !! Generates riser inclination angle at position s.
    real(wp), intent(in) :: s
    real(wp), intent(in) :: VPARAM(2)
    real(wp) :: A, LR

    A  = VPARAM(1)
    LR = VPARAM(2)
    theta = atan2(LR*s, A)
  end function fungeo


  pure real(wp) function Dfungeo(s, VPARAM) result(dthetads)
    !! Derivative dtheta/ds for the catenary-based geometry used in fungeo.
    real(wp), intent(in) :: s
    real(wp), intent(in) :: VPARAM(2)
    real(wp) :: A, LR, aux

    A  = VPARAM(1)
    LR = VPARAM(2)
    aux = LR*s/A
    dthetads = (1.0_wp/(1.0_wp + aux*aux))*(LR/A)
  end function Dfungeo

end module module_catenary
