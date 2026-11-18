module module_catenary
  !! Fortran translation of catenary.py
  !!
  !! Catenary geometry with a dimensionless, bracketed positive root.
  use mod_kinds,       only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
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


  real(wp) function solve_A(Z,X,status) result(A)
    !! Solve sinh(y/2)^2*2/y=Z/X, y=X/A. No procedure callbacks,
    !! executable-stack trampolines, or overflowing trial X/A are needed.
    real(wp),intent(in)::Z,X
    integer,intent(out),optional::status
    real(wp)::lo,hi,y,ratio,value
    integer::k
    A=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(status))status=1
    if(.not.all(ieee_is_finite([Z,X])))return
    if(min(Z,X)<=0.0_wp)return
    ratio=Z/X
    if(.not.ieee_is_finite(ratio).or.ratio<=0.0_wp)return
    lo=0.0_wp;hi=1.0_wp
    do k=1,10
      value=2.0_wp*sinh(hi/2.0_wp)**2/hi
      if(value>=ratio)exit
      hi=min(2.0_wp*hi,700.0_wp)
    enddo
    if(value<ratio)return
    do k=1,1100
      y=0.5_wp*(lo+hi)
      if(y<=0.0_wp)return
      value=2.0_wp*sinh(y/2.0_wp)**2/y
      if(value<ratio)then
        lo=y
      else
        hi=y
      endif
      if(hi-lo<=8.0_wp*epsilon(1.0_wp)*y)exit
    enddo
    A=X/(0.5_wp*(lo+hi))
    if(.not.ieee_is_finite(A))return
    if(present(status))status=0
  end function solve_A


  pure real(wp) function altura_catenaria(x, A) result(zx)
    real(wp), intent(in) :: x, A
    zx = A * (cosh(x / A) - 1.0_wp)
  end function altura_catenaria


  pure real(wp) function comprimento_arco(x, A) result(sx)
    real(wp), intent(in) :: x, A
    sx = A * sinh(x / A)
  end function comprimento_arco


  pure real(wp) function angulo_inclinacao(s,A) result(theta)
    !! tan(theta)=s/A for catenary arc length s; this is not atanh(s/A).
    real(wp),intent(in)::s,A
    theta=ieee_value(0.0_wp,ieee_quiet_nan)
    if(.not.all(ieee_is_finite([s,A])).or.A<=0.0_wp)return
    theta=atan2(s,A)
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
