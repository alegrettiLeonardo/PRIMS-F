module module_closure_law
  !! Verified closure contracts. Input rho/u are nondimensional; geometry,
  !! viscosity, gravity and surface tension are SI. Operators are evaluated in SI.
  !! STRATIFIED: Burr, Sistema Oleoduto-Riser, eqs. 8--13; Fanning friction.
  !! BENDIKSEN: Balino, Burr & Nemoto (2010), IJMF, DOI 10.1016/j.ijmultiphaseflow.2010.04.003,
  !! equations 32--36. This is not the different Montgomery (2002), p.151 form.
  !! BHAGWAT is reserved but unsupported until the complete primary correlation
  !! is verified. Legacy names containing 'swananda' never silently select Bendiksen.
  use mod_kinds, only: wp
  use mod_void_helpers, only: ffan, alpha2gamma
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  integer,parameter,public :: STRATIFIED=1, BENDIKSEN=2, BHAGWAT=3, MONTGOMERY=4
  integer,parameter,public :: CLOSURE_OK=0, CLOSURE_INVALID=10, CLOSURE_UNSUPPORTED=11
  integer,parameter,public :: CLOSURE_DOMAIN=12, CLOSURE_DERIVATIVE=13
  real(wp),parameter :: PI=acos(-1.0_wp), FI_FANNING=0.0142_wp
  public :: CdUd_bendiksen, CdUd_drift, drift_flux_ndim, fanning_wall_shear
  public :: CdUd_swananda, DCdUd_swananda, drift_flux_swananda_ndim, Ddrift_flux_swananda_ndim
  public :: RelEquilLocalPipe_comp, DRelEquiLocal_pipe_comp
  public :: lei_fechamento_or_ndim_simp, Dlei_fechamento_or_ndim_simp
contains
  pure real(wp) function invalid_real() result(x)
    x=ieee_value(0.0_wp,ieee_quiet_nan)
  end function
  real(wp) function fanning_wall_shear(rho,vel,diam,mu,rough,status) result(tau)
    real(wp),intent(in)::rho,vel,diam,mu,rough
    integer,intent(out)::status
    real(wp)::re,f
    status=CLOSURE_INVALID;tau=invalid_real()
    if(.not.all(ieee_is_finite([rho,vel,diam,mu,rough])))return
    if(rho<=0.0_wp.or.diam<=0.0_wp.or.mu<=0.0_wp.or.rough<0.0_wp)return
    re=rho*abs(vel)*diam/mu
    if(re<2000.0_wp)then
      ! Exact signed laminar limit, including zero speed; no 0*infinity.
      tau=8.0_wp*mu*vel/diam
    else
      f=ffan(rough/diam,re,status)
      if(status/=0)return
      tau=0.5_wp*f*rho*vel*abs(vel)
    endif
    status=CLOSURE_OK
  end function
  subroutine drift_coefficients_si(alpha,rhol,rhog,ul,ug,theta,D,G,model,Cd,Ud,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,D,G
    integer,intent(in)::model
    real(wp),intent(out)::Cd,Ud
    integer,intent(out)::status
    real(wp)::vref,fr,j
    Cd=invalid_real();Ud=Cd;status=CLOSURE_INVALID
    if(.not.all(ieee_is_finite([alpha,rhol,rhog,ul,ug,theta,D,G])))return
    if(alpha<=0.0_wp.or.alpha>=1.0_wp.or.rhol<=0.0_wp.or.rhog<=0.0_wp.or.D<=0.0_wp.or.G<0.0_wp)return
    if(model/=BENDIKSEN.and.model/=MONTGOMERY)then
      status=CLOSURE_UNSUPPORTED;return
    endif
    ! Reported inclination coverage of the long-bubble formulation. Signed phase
    ! speeds are retained: opposite phase directions are not clipped to co-current.
    ! Countercurrent-film validation is a separate physical qualification.
    if(theta < -PI/6.0_wp-1e-12_wp.or.theta > PI/2.0_wp+1e-12_wp)then
      status=CLOSURE_DOMAIN;return
    endif
    if(G==0.0_wp)then
      Cd=1.2_wp;Ud=0.0_wp;status=0;return
    endif
    vref=sqrt(G*D)
    if(model==BENDIKSEN)then
      j=(1.0_wp-alpha)*ul+alpha*ug
      fr=abs(j)/vref
      if(fr<3.5_wp)then
        Cd=1.05_wp+0.15_wp*sin(theta)
        Ud=vref*(0.35_wp*sin(theta)+0.54_wp*cos(theta))
      else
        Cd=1.2_wp
        Ud=0.35_wp*vref*sin(theta)
      endif
    else
      ! Montgomery (2002), printed p.151, equations 5.43--5.48.
      fr=abs((1.0_wp-alpha)*ul)/vref
      if(fr<3.5_wp)then
        Cd=1.05_wp+0.15_wp*sin(theta)**2
      else
        Cd=1.2_wp
      endif
      Ud=vref*(0.35_wp*sin(theta)+0.54_wp*cos(theta))
    endif
    status=CLOSURE_OK
  end subroutine
  subroutine residual_si(alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug,model,value,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug
    integer,intent(in)::model
    real(wp),intent(out)::value
    integer,intent(out)::status
    real(wp)::gamma,gammai,dl,dg,rel,ui,tl,tg,ti,Cd,Ud
    value=invalid_real();status=CLOSURE_INVALID
    if(.not.all(ieee_is_finite([alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug])))return
    if(alpha<=0.0_wp.or.alpha>=1.0_wp.or.rhol<=0.0_wp.or.rhog<=0.0_wp)return
    if(D<=0.0_wp.or.area<=0.0_wp.or.rough<0.0_wp.or.G<0.0_wp.or.mul<=0.0_wp.or.mug<=0.0_wp)return
    if(abs(area-PI*D**2/4.0_wp)>1e-6_wp*area)return
    if(model==STRATIFIED)then
      if(abs(theta)>=PI/2.0_wp)then
        status=CLOSURE_DOMAIN;return
      endif
      gamma=alpha2gamma(alpha,1e-13_wp,status)
      if(status/=0)return
      gammai=sin(PI*gamma)/PI
      dl=(1.0_wp-alpha)*D/gamma
      dg=alpha*D/(1.0_wp-gamma+gammai)
      rel=rhol*abs(ul)*dl/mul
      if(rel<=2000.0_wp)then
        ui=1.8_wp*ul
      else if(rel>=2200.0_wp)then
        ui=ul
      else
        ui=ul*(1.8_wp*(2200.0_wp-rel)+(rel-2000.0_wp))/200.0_wp
      endif
      tl=fanning_wall_shear(rhol,ul,dl,mul,rough,status)
      if(status/=0)return
      tg=fanning_wall_shear(rhog,ug,dg,mug,rough,status)
      if(status/=0)return
      ti=0.5_wp*FI_FANNING*rhog*(ug-ui)*abs(ug-ui)
      ! Burr eq.8 divided by pi*D; beta=-theta (downhill beta positive).
      value=tg*(1.0_wp-gamma)/alpha-tl*gamma/(1.0_wp-alpha)+ &
        ti*gammai*(1.0_wp/(1.0_wp-alpha)+1.0_wp/alpha)- &
        (rhol-rhog)*G*area*sin(theta)/(PI*D)
    else
      call drift_coefficients_si(alpha,rhol,rhog,ul,ug,theta,D,G,model,Cd,Ud,status)
      if(status/=0)return
      value=-alpha*ug+alpha*(Cd*((1.0_wp-alpha)*ul+alpha*ug)+Ud)
    endif
    if(.not.ieee_is_finite(value))then
      value=invalid_real();status=CLOSURE_DOMAIN;return
    endif
    status=CLOSURE_OK
  end subroutine
  subroutine lei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol, &
      value,closure_id,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::value
    integer,intent(in),optional::closure_id
    integer,intent(out),optional::status
    integer::model,ierr
    value=invalid_real();ierr=CLOSURE_INVALID
    if(present(status))status=ierr
    if(.not.all(ieee_is_finite([w_u,w_rho,tol,sigma])))return
    if(w_u<=0.0_wp.or.w_rho<=0.0_wp.or.tol<=0.0_wp.or.sigma<0.0_wp)return
    model=BENDIKSEN
    if(theta<0.0_wp)model=STRATIFIED
    if(present(closure_id))model=closure_id
    call residual_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,AREA,EPS,G,MUL,MUG,model,value,ierr)
    if(ierr==0)then
      if(model==STRATIFIED)then
        value=value/(w_rho*w_u*w_u)
      else
        value=value/w_u
      endif
    endif
    if(present(status))status=ierr
  end subroutine
  subroutine Dlei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol, &
      dfda,dfdrhog,dfdrhol,dfdug,dfdul,closure_id,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::dfda,dfdrhog,dfdrhol,dfdug,dfdul
    integer,intent(in),optional::closure_id
    integer,intent(out),optional::status
    real(wp)::x(5),xp(5),xm(5),der(5),scale(5),h,f0,fp,fm,fo,cd,ud,fr,den,j
    integer::model,ierr,k
    dfda=invalid_real();dfdrhog=dfda;dfdrhol=dfda;dfdug=dfda;dfdul=dfda
    if(present(status))status=CLOSURE_INVALID
    model=BENDIKSEN
    if(theta<0.0_wp)model=STRATIFIED
    if(present(closure_id))model=closure_id
    call lei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,f0,model,ierr)
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    if(model==BENDIKSEN.or.model==MONTGOMERY)then
      call drift_coefficients_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,G,model,cd,ud,ierr)
      if(G>0.0_wp)then
        den=sqrt(G*DH)
        j=(1.0_wp-alpha)*ul+alpha*ug
        if(model==BENDIKSEN)then
          fr=abs(j*w_u)/den
        else
          fr=abs((1.0_wp-alpha)*ul*w_u)/den
        endif
        if(abs(fr-3.5_wp)<1e-10_wp)then
          if(present(status))status=CLOSURE_DERIVATIVE
          return
        endif
      endif
      ! Coefficients are constant with respect to the five state variables
      ! within each Froude branch. Derivative at the branch jump is undefined.
      j=(1.0_wp-alpha)*ul+alpha*ug
      dfda=-ug+cd*j+ud/w_u+alpha*cd*(ug-ul)
      dfdrhog=0.0_wp;dfdrhol=0.0_wp
      dfdug=-alpha+alpha*alpha*cd
      dfdul=alpha*cd*(1.0_wp-alpha)
      if(present(status))status=CLOSURE_OK
      return
    endif
    ! Finite differences of the same SI residual, then chain-rule conversion.
    ! This avoids dependence of perturbation size on arbitrary rho/u scales.
    x=[alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u]
    scale=[1.0_wp,max(abs(x(2)),1.0_wp),max(abs(x(3)),1.0_wp), &
      max(abs(x(4)),1.0_wp),max(abs(x(5)),1.0_wp)]
    do k=1,5
      h=epsilon(1.0_wp)**(1.0_wp/3.0_wp)*scale(k)
      if(k==2.or.k==3)h=min(h,0.25_wp*x(k))
      xp=x;xm=x
      xp(k)=x(k)+h;xm(k)=x(k)-h
      if(k==1.and.xm(1)<=0.0_wp)then
        call eval(x,fo);if(ierr/=0)exit
        call eval(xp,fp);if(ierr/=0)exit
        xp(k)=x(k)+2.0_wp*h
        call eval(xp,fm);if(ierr/=0)exit
        der(k)=(-3.0_wp*fo+4.0_wp*fp-fm)/(2.0_wp*h)
      else if(k==1.and.xp(1)>=1.0_wp)then
        call eval(x,fo);if(ierr/=0)exit
        call eval(xm,fm);if(ierr/=0)exit
        xm(k)=x(k)-2.0_wp*h
        call eval(xm,fp);if(ierr/=0)exit
        der(k)=(3.0_wp*fo-4.0_wp*fm+fp)/(2.0_wp*h)
      else
        call eval(xp,fp);if(ierr/=0)exit
        call eval(xm,fm);if(ierr/=0)exit
        der(k)=(fp-fm)/(2.0_wp*h)
      endif
    enddo
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    den=w_rho*w_u*w_u
    dfda=der(1)/den;dfdrhol=der(2)*w_rho/den;dfdrhog=der(3)*w_rho/den
    dfdul=der(4)*w_u/den;dfdug=der(5)*w_u/den
    if(present(status))status=CLOSURE_OK
  contains
    subroutine eval(xx,ff)
      real(wp),intent(in)::xx(5)
      real(wp),intent(out)::ff
      call residual_si(xx(1),xx(2),xx(3),xx(4),xx(5),theta,DH,AREA,EPS,G,MUL,MUG,model,ff,ierr)
    end subroutine
  end subroutine
  subroutine CdUd_drift(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,closure_id,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::Cd,Ud
    integer,intent(in),optional::closure_id
    integer,intent(out),optional::status
    integer::model,ierr
    Cd=invalid_real();Ud=Cd
    if(present(status))status=CLOSURE_INVALID
    if(.not.all(ieee_is_finite([w_u,w_rho])).or.w_u<=0.0_wp.or.w_rho<=0.0_wp)return
    model=BENDIKSEN
    if(present(closure_id))model=closure_id
    call drift_coefficients_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,G,model,Cd,Ud,ierr)
    if(present(status))status=ierr
  end subroutine
  subroutine CdUd_bendiksen(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::Cd,Ud
    integer,intent(out),optional::status
    call CdUd_drift(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,BENDIKSEN,status)
  end subroutine
  subroutine drift_flux_ndim(ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,value,closure_id,status)
    real(wp),intent(in)::ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::value
    integer,intent(in),optional::closure_id
    integer,intent(out),optional::status
    integer::model
    model=BENDIKSEN
    if(present(closure_id))model=closure_id
    call lei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,value,model,status)
  end subroutine
  subroutine RelEquilLocalPipe_comp(ul,ug,rhol,rhog,alpha,BETA,DH,AREA,EPS,G,MUL,MUG,w_u,w_rho,tol,value,status)
    real(wp),intent(in)::ul,ug,rhol,rhog,alpha,BETA,DH,AREA,EPS,G,MUL,MUG,w_u,w_rho,tol
    real(wp),intent(out)::value
    integer,intent(out),optional::status
    call lei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,-BETA,DH,AREA,EPS,G,MUL,MUG,0.0_wp,w_u,w_rho,tol,value,STRATIFIED,status)
  end subroutine
  subroutine DRelEquiLocal_pipe_comp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,w_u,w_rho,tol, &
      dfda,dfdrhog,dfdrhol,dfdug,dfdul,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,w_u,w_rho,tol
    real(wp),intent(out)::dfda,dfdrhog,dfdrhol,dfdug,dfdul
    integer,intent(out),optional::status
    call Dlei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,-theta,DH,AREA,EPS,G,MUL,MUG,0.0_wp,w_u,w_rho,tol, &
      dfda,dfdrhog,dfdrhol,dfdug,dfdul,STRATIFIED,status)
  end subroutine
  subroutine CdUd_swananda(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::Cd,Ud
    integer,intent(out),optional::status
    call CdUd_drift(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,BHAGWAT,status)
  end subroutine
  subroutine DCdUd_swananda(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol, &
      dC0da,dC0drhol,dC0drhog,dC0dul,dC0dug,dUdda,dUddrhol,dUddrhog,dUddul,dUddug,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::dC0da,dC0drhol,dC0drhog,dC0dul,dC0dug,dUdda,dUddrhol,dUddrhog,dUddul,dUddug
    integer,intent(out),optional::status
    dC0da=invalid_real();dC0drhol=dC0da;dC0drhog=dC0da;dC0dul=dC0da;dC0dug=dC0da
    dUdda=dC0da;dUddrhol=dC0da;dUddrhog=dC0da;dUddul=dC0da;dUddug=dC0da
    if(present(status))status=CLOSURE_UNSUPPORTED
  end subroutine
  subroutine drift_flux_swananda_ndim(ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,value,status)
    real(wp),intent(in)::ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::value
    integer,intent(out),optional::status
    call drift_flux_ndim(ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,value,BHAGWAT,status)
  end subroutine
  subroutine Ddrift_flux_swananda_ndim(ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol, &
      dfda,dfdrhog,dfdrhol,dfdug,dfdul,status)
    real(wp),intent(in)::ul,ug,rhol,rhog,alpha,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::dfda,dfdrhog,dfdrhol,dfdug,dfdul
    integer,intent(out),optional::status
    call Dlei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol, &
      dfda,dfdrhog,dfdrhol,dfdug,dfdul,BHAGWAT,status)
  end subroutine
end module module_closure_law
