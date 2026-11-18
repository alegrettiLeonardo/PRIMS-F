module mod_void_helpers
  !! Geometry and Fanning friction. No hidden root fallback or negative-Re clipping.
  use mod_kinds, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: clamp, safe_div, calcular_fatrito, ffan, Dffan, fvgamma, alpha2gamma
  real(wp), parameter :: PI=acos(-1.0_wp)
contains
  pure real(wp) function clamp(x,lo,hi) result(y)
    real(wp),intent(in)::x,lo,hi
    y=min(max(x,lo),hi)
  end function
  pure real(wp) function safe_div(a,b,eps) result(y)
    real(wp),intent(in)::a,b
    real(wp),intent(in),optional::eps
    real(wp)::e,bb
    e=1e-30_wp
    if(present(eps)) e=eps
    bb=b
    if(abs(bb)<=e)bb=sign(e,bb)
    y=a/bb
  end function
  real(wp) function calcular_fatrito(Re,epsD,status) result(f)
    !! Explicit turbulent Fanning approximation from the supplied implementation.
    !! The -4 logarithm coefficient is Fanning, not Darcy.
    real(wp),intent(in)::Re,epsD
    integer,intent(out),optional::status
    real(wp)::q,b
    if(present(status))status=10
    f=ieee_value(0.0_wp,ieee_quiet_nan)
    if(.not.all(ieee_is_finite([Re,epsD])))return
    if(Re<=0.0_wp.or.epsD<0.0_wp)return
    q=5.8506_wp/Re**0.8981_wp+epsD**1.1098_wp/2.8257_wp
    if(q<=0.0_wp)return
    b=epsD/3.7065_wp-(5.0452_wp/Re)*log10(q)
    if(b<=0.0_wp)return
    q=-4.0_wp*log10(b)
    if(abs(q)<=tiny(1.0_wp))return
    f=1.0_wp/(q*q)
    if(present(status))status=0
  end function
  real(wp) function ffan(epsD,Re,status) result(f)
    real(wp),intent(in)::epsD,Re
    integer,intent(out),optional::status
    integer::ierr
    real(wp)::ft
    f=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(status))status=10
    if(.not.all(ieee_is_finite([Re,epsD])))return
    if(Re<=0.0_wp.or.epsD<0.0_wp)return
    if(Re<2000.0_wp)then
      f=16.0_wp/Re
    else if(Re>2300.0_wp)then
      f=calcular_fatrito(Re,epsD,ierr)
      if(ierr/=0)return
    else
      ft=calcular_fatrito(2300.0_wp,epsD,ierr)
      if(ierr/=0)return
      f=((Re-2000.0_wp)*ft+(2300.0_wp-Re)*(16.0_wp/2000.0_wp))/300.0_wp
    endif
    if(present(status))status=0
  end function
  subroutine Dffan(Re,epsD,dRe,dE,status)
    real(wp),intent(in)::Re,epsD
    real(wp),intent(out)::dRe,dE
    integer,intent(out),optional::status
    real(wp)::re_eval,inner,ll,b,dbdr,dbde,fac,ft,der,dee
    integer::ierr
    dRe=ieee_value(0.0_wp,ieee_quiet_nan);dE=dRe
    if(present(status))status=10
    if(.not.all(ieee_is_finite([Re,epsD])))return
    if(Re<=0.0_wp.or.epsD<0.0_wp)return
    if(Re<2000.0_wp)then
      dRe=-16.0_wp/Re**2;dE=0.0_wp
    else
      re_eval=max(Re,2300.0_wp)
      inner=5.8506_wp/re_eval**0.8981_wp+epsD**1.1098_wp/2.8257_wp
      ll=log10(inner)
      b=epsD/3.7065_wp-5.0452_wp*ll/re_eval
      if(b<=0.0_wp)return
      fac=-1.0_wp/(8.0_wp*log10(b)**3*b*log(10.0_wp))
      dbdr=5.0452_wp*ll/re_eval**2+ &
        (5.0452_wp/re_eval)*(0.8981_wp*5.8506_wp/re_eval**1.8981_wp)/(inner*log(10.0_wp))
      dbde=1.0_wp/3.7065_wp-(5.0452_wp/re_eval)* &
        (1.1098_wp*epsD**0.1098_wp/2.8257_wp)/(inner*log(10.0_wp))
      der=fac*dbdr;dee=fac*dbde
      if(Re>2300.0_wp)then
        dRe=der;dE=dee
      else
        ft=calcular_fatrito(2300.0_wp,epsD,ierr)
        if(ierr/=0)return
        dRe=(ft-16.0_wp/2000.0_wp)/300.0_wp
        ! Unlike the previous version, the transition retains roughness sensitivity.
        dE=(Re-2000.0_wp)*dee/300.0_wp
      endif
    endif
    if(present(status))status=0
  end subroutine
  pure real(wp) function fvgamma(alpha,gamma) result(v)
    real(wp),intent(in)::alpha,gamma
    v=alpha-1.0_wp+gamma-sin(2.0_wp*PI*gamma)/(2.0_wp*PI)
  end function
  real(wp) function alpha2gamma(alpha,tol,status) result(gamma)
    real(wp),intent(in)::alpha,tol
    integer,intent(out),optional::status
    real(wp)::lo,hi,mid,f,xtol
    integer::k
    gamma=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(status))status=10
    if(.not.all(ieee_is_finite([alpha,tol])))return
    if(alpha<0.0_wp.or.alpha>1.0_wp.or.tol<=0.0_wp)return
    if(alpha==0.0_wp)then
      gamma=1.0_wp
    else if(alpha==1.0_wp)then
      gamma=0.0_wp
    else
      lo=0.0_wp;hi=1.0_wp
      xtol=max(4.0_wp*epsilon(1.0_wp),min(tol,1e-12_wp))
      do k=1,100
        mid=0.5_wp*(lo+hi)
        f=fvgamma(alpha,mid)
        if(f>0.0_wp)then
          hi=mid
        else
          lo=mid
        endif
        if(hi-lo<=xtol)exit
      enddo
      gamma=0.5_wp*(lo+hi)
    endif
    if(present(status))status=0
  end function
end module mod_void_helpers
