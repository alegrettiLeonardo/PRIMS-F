module mod_rootfinding
  !! Simple 1D rootfinding utilities (bisection + Brent).
  !!
  !! Robustness fixes (per gfortran warnings):
  !! - Avoid direct equality comparisons on reals.
  !! - Use IEEE_IS_NAN for NaN checks.

  use mod_kinds, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  implicit none
  private

  public :: func_1d
  public :: bisect_root
  public :: brent_root

  abstract interface
    function func_1d(x) result(f)
      import :: wp
      real(wp), intent(in) :: x
      real(wp) :: f
    end function func_1d
  end interface

contains

  pure logical function nearly_zero(x) result(is0)
    real(wp), intent(in) :: x
    real(wp) :: eps0
    eps0 = 10.0_wp * epsilon(1.0_wp)
    is0 = abs(x) <= eps0
  end function nearly_zero


  function bisect_root(f, a, b, xtol, maxit, info) result(root)
    procedure(func_1d) :: f
    real(wp), intent(in) :: a, b
    real(wp), intent(in) :: xtol
    integer, intent(in) :: maxit
    integer, intent(out) :: info
    real(wp) :: root

    real(wp) :: lo, hi, mid
    real(wp) :: flo, fhi, fmid
    integer :: it

    info = 0
    lo = a
    hi = b
    flo = f(lo)
    fhi = f(hi)

    if (ieee_is_nan(flo) .or. ieee_is_nan(fhi)) then
      info = 3
      root = 0.5_wp*(a+b)
      return
    end if

    if (nearly_zero(flo)) then
      root = lo
      return
    end if
    if (nearly_zero(fhi)) then
      root = hi
      return
    end if

    if (flo*fhi > 0.0_wp) then
      info = 1
      root = 0.5_wp*(a+b)
      return
    end if

    do it = 1, maxit
      mid = 0.5_wp*(lo + hi)
      fmid = f(mid)

      if (ieee_is_nan(fmid)) then
        info = 3
        root = mid
        return
      end if

      if (nearly_zero(fmid) .or. abs(hi-lo) <= xtol) exit

      if (flo*fmid < 0.0_wp) then
        hi = mid
        fhi = fmid
      else
        lo = mid
        flo = fmid
      end if
    end do

    root = 0.5_wp*(lo + hi)

    if (.not. ieee_is_finite(root)) info = 3

  end function bisect_root


  function brent_root(f, a, b, xtol, maxit, info) result(root)
    !! Brent's method for a single root in [a,b].
    !! Returns info=0 on success, info=1 if no bracket, info=2 if maxit,
    !! info=3 if NaN/Inf encountered.
    procedure(func_1d) :: f
    real(wp), intent(in) :: a, b
    real(wp), intent(in) :: xtol
    integer, intent(in) :: maxit
    integer, intent(out) :: info
    real(wp) :: root

    real(wp) :: aa, bb, cc
    real(wp) :: fa, fb, fc
    real(wp) :: d, e
    real(wp) :: p, q, r, s
    real(wp) :: tol1, xm
    real(wp) :: eps0
    real(wp) :: tmp
    integer :: it

    info = 0
    eps0 = epsilon(1.0_wp)

    aa = a
    bb = b
    fa = f(aa)
    fb = f(bb)

    if (ieee_is_nan(fa) .or. ieee_is_nan(fb)) then
      info = 3
      root = 0.5_wp*(a+b)
      return
    end if

    if (nearly_zero(fa)) then
      root = aa
      return
    end if
    if (nearly_zero(fb)) then
      root = bb
      return
    end if

    if (fa*fb > 0.0_wp) then
      info = 1
      root = 0.5_wp*(a+b)
      return
    end if

    cc = aa
    fc = fa
    d = bb - aa
    e = d

    do it = 1, maxit

      ! Ensure |fb| <= |fc|
      if (abs(fc) < abs(fb)) then
        tmp = aa;  aa = bb;  bb = cc;  cc = tmp
        tmp = fa;  fa = fb;  fb = fc;  fc = tmp
      end if

      tol1 = 2.0_wp*eps0*abs(bb) + 0.5_wp*xtol
      xm   = 0.5_wp*(cc - bb)

      if (abs(xm) <= tol1 .or. nearly_zero(fb)) exit

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        ! Attempt inverse quadratic interpolation
        s = fb/fa
        if (abs(aa-cc) <= eps0*max(1.0_wp, abs(aa), abs(cc))) then
          p = 2.0_wp*xm*s
          q = 1.0_wp - s
        else
          q = fa/fc
          r = fb/fc
          p = s*(2.0_wp*xm*q*(q-r) - (bb-aa)*(r-1.0_wp))
          q = (q-1.0_wp)*(r-1.0_wp)*(s-1.0_wp)
        end if

        if (p > 0.0_wp) q = -q
        p = abs(p)

        if (2.0_wp*p < min(3.0_wp*xm*q - abs(tol1*q), abs(e*q))) then
          e = d
          d = p/q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if

      aa = bb
      fa = fb

      if (abs(d) > tol1) then
        bb = bb + d
      else
        bb = bb + sign(tol1, xm)
      end if

      fb = f(bb)

      if (ieee_is_nan(fb)) then
        info = 3
        root = bb
        return
      end if

      if ((fb > 0.0_wp .and. fc > 0.0_wp) .or. (fb < 0.0_wp .and. fc < 0.0_wp)) then
        cc = aa
        fc = fa
        d = bb - aa
        e = d
      end if

    end do

    root = bb
    if (it >= maxit) info = 2

    if (.not. ieee_is_finite(root)) info = 3

  end function brent_root

end module mod_rootfinding
