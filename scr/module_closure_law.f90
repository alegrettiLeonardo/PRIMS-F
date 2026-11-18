module module_closure_law
  !! Verified closure contracts. Input rho/u are nondimensional; geometry,
  !! viscosity, gravity and surface tension are SI. Operators are evaluated in SI.
  !! STRATIFIED: Burr, Sistema Oleoduto-Riser, eqs. 8--13; Fanning friction.
  !! BENDIKSEN: Balino, Burr & Nemoto (2010), IJMF, DOI 10.1016/j.ijmultiphaseflow.2010.04.003,
  !! equations 32--36. This is not the different Montgomery (2002), p.151 form.
  !! BHAGWAT: the explicit transcription in Azevedo, Balino & Burr (IMAACA 2015),
  !! printed p.121, eqs.116--128. The 2014 original was not available in full.
  !! Legacy names containing 'swananda' select this identified transcription.
  use module_chexal_1992, only: chexal_coefficients_si, chexal_residual_si
  use mod_kinds, only: wp
  use mod_void_helpers, only: ffan, Dffan, alpha2gamma
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  integer,parameter,public :: STRATIFIED=1, BENDIKSEN=2, BHAGWAT=3, MONTGOMERY=4
  ! IDs 6/7 select the independently specified high/low void branch of 1992.
  ! ID5 never silently chooses between countercurrent branches.
  integer,parameter,public :: CHEXAL_1992=5, CHEXAL_1992_HIGH=6, CHEXAL_1992_LOW=7
  public :: is_chexal, chexal_branch_id
  integer,parameter,public :: CLOSURE_OK=0, CLOSURE_INVALID=10, CLOSURE_UNSUPPORTED=11
  integer,parameter,public :: CLOSURE_DOMAIN=12, CLOSURE_DERIVATIVE=13
  real(wp),parameter :: PI=acos(-1.0_wp), FI_FANNING=0.0142_wp
  public :: CdUd_bendiksen, CdUd_bhagwat, CdUd_drift, drift_flux_ndim, fanning_wall_shear
  public :: CdUd_swananda, DCdUd_swananda, drift_flux_swananda_ndim, Ddrift_flux_swananda_ndim
  public :: RelEquilLocalPipe_comp, DRelEquiLocal_pipe_comp
  public :: lei_fechamento_or_ndim_simp, Dlei_fechamento_or_ndim_simp
contains
  pure logical function is_chexal(model)
    integer,intent(in)::model
    is_chexal=model>=CHEXAL_1992.and.model<=CHEXAL_1992_LOW
  end function
  pure integer function chexal_branch_id(model) result(branch)
    integer,intent(in)::model
    branch=0
    if(model==CHEXAL_1992_HIGH)branch=1
    if(model==CHEXAL_1992_LOW)branch=-1
  end function
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
  subroutine shear_derivatives_si(rho,vel,diam,mu,rough,tau,drho,dvel,ddiam,status)
    !! Exact derivatives within the branch selected by fanning_wall_shear.
    !! The interpolation at Re=2000..2300 remains unchanged and is only C0.
    real(wp),intent(in)::rho,vel,diam,mu,rough
    real(wp),intent(out)::tau,drho,dvel,ddiam
    integer,intent(out)::status
    real(wp)::re,f,fre,fe
    tau=invalid_real();drho=tau;dvel=tau;ddiam=tau
    status=CLOSURE_INVALID
    if(.not.all(ieee_is_finite([rho,vel,diam,mu,rough])))return
    if(rho<=0.0_wp.or.diam<=0.0_wp.or.mu<=0.0_wp.or.rough<0.0_wp)return
    re=rho*abs(vel)*diam/mu
    if(re<2000.0_wp)then
      tau=8.0_wp*mu*vel/diam
      drho=0.0_wp;dvel=8.0_wp*mu/diam;ddiam=-tau/diam
    else
      f=ffan(rough/diam,re,status)
      if(status/=0)return
      call Dffan(re,rough/diam,fre,fe,status)
      if(status/=0)return
      tau=0.5_wp*f*rho*vel*abs(vel)
      drho=0.5_wp*vel*abs(vel)*(f+fre*re)
      dvel=rho*abs(vel)*(f+0.5_wp*fre*re)
      ddiam=0.5_wp*rho*vel*abs(vel)*(fre*re/diam-fe*rough/diam**2)
    endif
    if(.not.all(ieee_is_finite([tau,drho,dvel,ddiam])))then
      status=CLOSURE_DERIVATIVE;return
    endif
    status=CLOSURE_OK
  end subroutine
  subroutine stratified_derivatives_si(alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug,der,status)
    !! d(Phi_SI)/d[alpha,rhol,rhog,ul,ug], with Phi divided by pi*D.
    !! No finite-difference perturbation may cross a Reynolds branch here.
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug
    real(wp),intent(out)::der(5)
    integer,intent(out)::status
    real(wp)::gamma,gammai,gamma_a,gammai_a,sg,den,dl,dg,dl_a,dg_a
    real(wp)::tl,tg,ti,tr,tu,td,re,fac,fac_re,ui,rel
    real(wp)::dtl(5),dtg(5),dti(5),drel(5),dui(5)
    real(wp)::wl,wg,wi,wl_a,wg_a,wi_a,grav
    der=invalid_real();status=CLOSURE_DERIVATIVE
    gamma=alpha2gamma(alpha,1e-13_wp,status)
    if(status/=0)return
    sg=sin(PI*gamma)
    if(sg==0.0_wp)then
      status=CLOSURE_DERIVATIVE;return
    endif
    gammai=sg/PI
    gamma_a=-0.5_wp/(sg*sg)
    gammai_a=cos(PI*gamma)*gamma_a
    den=1.0_wp-gamma+gammai
    dl=(1.0_wp-alpha)*D/gamma
    dg=alpha*D/den
    dl_a=-D/gamma-(1.0_wp-alpha)*D*gamma_a/gamma**2
    dg_a=D/den-alpha*D*(-gamma_a+gammai_a)/den**2
    call shear_derivatives_si(rhol,ul,dl,mul,rough,tl,tr,tu,td,status)
    if(status/=0)return
    dtl=[td*dl_a,tr,0.0_wp,tu,0.0_wp]
    call shear_derivatives_si(rhog,ug,dg,mug,rough,tg,tr,tu,td,status)
    if(status/=0)return
    dtg=[td*dg_a,0.0_wp,tr,0.0_wp,tu]
    re=rhol*abs(ul)*dl/mul
    if(re<=2000.0_wp)then
      fac=1.8_wp;fac_re=0.0_wp
    else if(re>=2200.0_wp)then
      fac=1.0_wp;fac_re=0.0_wp
    else
      fac=(1.8_wp*(2200.0_wp-re)+(re-2000.0_wp))/200.0_wp
      fac_re=-0.004_wp
    endif
    drel=[rhol*abs(ul)*dl_a/mul,abs(ul)*dl/mul,0.0_wp,0.0_wp,0.0_wp]
    if(ul/=0.0_wp)drel(4)=rhol*sign(1.0_wp,ul)*dl/mul
    ui=fac*ul
    dui=ul*fac_re*drel
    dui(4)=dui(4)+fac
    rel=ug-ui
    ti=0.5_wp*FI_FANNING*rhog*rel*abs(rel)
    dti=-FI_FANNING*rhog*abs(rel)*dui
    dti(3)=dti(3)+0.5_wp*FI_FANNING*rel*abs(rel)
    dti(5)=dti(5)+FI_FANNING*rhog*abs(rel)
    wg=(1.0_wp-gamma)/alpha
    wl=gamma/(1.0_wp-alpha)
    wi=gammai*(1.0_wp/alpha+1.0_wp/(1.0_wp-alpha))
    wg_a=-gamma_a/alpha-(1.0_wp-gamma)/alpha**2
    wl_a=gamma_a/(1.0_wp-alpha)+gamma/(1.0_wp-alpha)**2
    wi_a=gammai_a*(1.0_wp/alpha+1.0_wp/(1.0_wp-alpha))+ &
      gammai*(-1.0_wp/alpha**2+1.0_wp/(1.0_wp-alpha)**2)
    der=wg*dtg-wl*dtl+wi*dti
    der(1)=der(1)+wg_a*tg-wl_a*tl+wi_a*ti
    grav=G*area*sin(theta)/(PI*D)
    der(2)=der(2)-grav
    der(3)=der(3)+grav
    if(.not.all(ieee_is_finite(der)))then
      der=invalid_real();status=CLOSURE_DERIVATIVE;return
    endif
    status=CLOSURE_OK
  end subroutine
  subroutine colebrook_bhagwat(re,relative_rough,y,dy_dre,status)
    !! Eq.118: y=1/sqrt(f_tp), H=y+4*log10(eps/(3.7D)+1.256*y/Re)=0.
    !! H is strictly increasing; safeguarded Newton has one positive root.
    !! This coefficient is the literal correlation input, not the wall-friction
    !! routine with its separate laminar/transition model.
    real(wp),intent(in)::re,relative_rough
    real(wp),intent(out)::y,dy_dre
    integer,intent(out)::status
    real(wp)::lo,hi,q,h,hy,trial,t,k
    integer::i
    y=invalid_real();dy_dre=y;status=CLOSURE_DOMAIN
    if(re<=0.0_wp.or.relative_rough<0.0_wp.or.relative_rough>=3.7_wp)return
    k=4.0_wp/log(10.0_wp);t=1.256_wp/re
    lo=0.0_wp;hi=1.0_wp
    do i=1,64
      q=relative_rough/3.7_wp+t*hi
      h=hi+4.0_wp*log10(q)
      if(h>=0.0_wp)exit
      hi=2.0_wp*hi
    enddo
    if(i>64)return
    y=0.5_wp*hi
    do i=1,100
      q=relative_rough/3.7_wp+t*y
      h=y+4.0_wp*log10(q)
      hy=1.0_wp+k*t/q
      if(abs(h)<=64.0_wp*epsilon(y)*max(1.0_wp,y))exit
      if(h>0.0_wp)then
        hi=y
      else
        lo=y
      endif
      trial=y-h/hy
      if(trial<=lo.or.trial>=hi)trial=0.5_wp*(lo+hi)
      y=trial
    enddo
    if(i>100.or.y<=0.0_wp)then
      y=invalid_real();status=CLOSURE_DERIVATIVE;return
    endif
    dy_dre=(k*t*y/(re*q))/hy
    if(.not.all(ieee_is_finite([y,dy_dre])))then
      y=invalid_real();dy_dre=y;return
    endif
    status=CLOSURE_OK
  end subroutine
  subroutine bhagwat_coefficients_si(alpha,rhol,rhog,ul,ug,theta,D,rough,G,mul,sigma,Cd,Ud,status,dCd,dUd)
    !! Azevedo, Balino & Burr (2015), p.121, eqs.116--128, circular duct.
    !! dCd/dUd, when requested, are with respect to [alpha,rhol,rhog,ul,ug] SI.
    !! Signed velocities are checked, never replaced by absolute values. This
    !! implementation accepts co-current nonnegative phase fluxes only; the
    !! transcription does not establish a countercurrent-film extension.
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,D,rough,G,mul,sigma
    real(wp),intent(out)::Cd,Ud
    integer,intent(out)::status
    real(wp),intent(out),optional::dCd(5),dUd(5)
    real(wp)::jl,jg,j,r,delta,mass,x,beta,re,y,yre,weight,la,c2,c3,c4,frsg,mu_star
    real(wp)::aa,bb,cc,base,logbase,power,qq,tt,vref,co,cd_der(5),ud_der(5)
    real(wp)::djl(5),djg(5),dj(5),dr(5),ddelta(5),dmass(5),dx(5),dbeta(5),dre(5)
    real(wp)::dweight(5),daa(5),dbb(5),dcc(5),dqq(5),dtt(5),dy(5)
    logical::want_derivatives
    Cd=invalid_real();Ud=Cd;status=CLOSURE_INVALID
    if(present(dCd))dCd=Cd
    if(present(dUd))dUd=Cd
    if(.not.all(ieee_is_finite([alpha,rhol,rhog,ul,ug,theta,D,rough,G,mul,sigma])))return
    if(alpha<=0.0_wp.or.alpha>=1.0_wp.or.rhog<=0.0_wp.or.rhol<=rhog)return
    if(D<=0.0_wp.or.G<=0.0_wp.or.mul<=0.0_wp.or.sigma<=0.0_wp.or.rough<0.0_wp)return
    status=CLOSURE_DOMAIN
    if(theta< -50.0_wp*PI/180.0_wp.or.theta>PI/2.0_wp)return
    if(ul<0.0_wp.or.ug<0.0_wp)return
    jl=(1.0_wp-alpha)*ul;jg=alpha*ug;j=jl+jg
    if(j<=0.0_wp)return
    r=rhog/rhol;delta=rhol-rhog
    mass=rhol*jl+rhog*jg
    x=rhog*jg/mass;beta=jg/j
    re=rhol*j*D/mul
    call colebrook_bhagwat(re,rough/D,y,yre,status)
    if(status/=0)return
    ! Eq.116 written with a bounded weight to avoid squaring large Re.
    if(re<1000.0_wp)then
      weight=(re/1000.0_wp)**2/(1.0_wp+(re/1000.0_wp)**2)
    else
      weight=1.0_wp/(1.0_wp+(1000.0_wp/re)**2)
    endif
    co=cos(theta)
    base=(1.0_wp+r*r*co)/(1.0_wp+co)
    logbase=log(base);power=0.2_wp*(1.0_wp-alpha)
    aa=2.0_wp-r*r;bb=exp(power*logbase)
    qq=(2.6_wp-beta)**0.15_wp-1.0_wp/y
    tt=(1.0_wp-x)**1.5_wp
    cc=0.2_wp*(1.0_wp-sqrt(r))*qq*tt
    Cd=aa*(1.0_wp-weight)+(bb+cc)*weight
    mu_star=mul/0.001_wp;c2=1.0_wp
    if(mu_star>10.0_wp)c2=(0.434_wp/log10(mu_star))**0.15_wp
    la=sqrt(sigma/(G*delta))/D;c3=1.0_wp
    if(la<0.025_wp)c3=(la/0.025_wp)**0.9_wp
    c4=1.0_wp;frsg=0.0_wp
    if(theta<0.0_wp)then
      frsg=sqrt(rhog/delta)*jg/sqrt(G*D*co)
      if(frsg<=0.1_wp)c4=-1.0_wp
    endif
    vref=sqrt(G*D*(1.0_wp-r))
    Ud=(0.35_wp*sin(theta)+0.45_wp*co)*vref*sqrt(1.0_wp-alpha)*c2*c3*c4
    if(.not.all(ieee_is_finite([Cd,Ud])))then
      Cd=invalid_real();Ud=Cd;status=CLOSURE_DOMAIN;return
    endif
    want_derivatives=present(dCd).or.present(dUd)
    if(.not.want_derivatives)then
      status=CLOSURE_OK;return
    endif
    ! C4 jumps at Fr_sg=.1 for downhill angles; no derivative is invented.
    if(theta<0.0_wp.and.abs(frsg-0.1_wp)<=64.0_wp*epsilon(frsg))then
      status=CLOSURE_DERIVATIVE;return
    endif
    djl=[-ul,0.0_wp,0.0_wp,1.0_wp-alpha,0.0_wp]
    djg=[ug,0.0_wp,0.0_wp,0.0_wp,alpha];dj=djl+djg
    dr=[0.0_wp,-r/rhol,1.0_wp/rhol,0.0_wp,0.0_wp]
    ddelta=[0.0_wp,1.0_wp,-1.0_wp,0.0_wp,0.0_wp]
    dmass=rhol*djl+rhog*djg;dmass(2)=dmass(2)+jl;dmass(3)=dmass(3)+jg
    dx=rhog*djg;dx(3)=dx(3)+jg;dx=(dx-x*dmass)/mass
    dbeta=(djg-beta*dj)/j
    dre=rhol*D/mul*dj;dre(2)=dre(2)+j*D/mul
    dy=yre*dre;dweight=2.0_wp*weight*(1.0_wp-weight)/re*dre
    daa=-2.0_wp*r*dr
    dbb=bb*power*(2.0_wp*r*co/(1.0_wp+r*r*co))*dr
    dbb(1)=dbb(1)-0.2_wp*bb*logbase
    dqq=-0.15_wp*(2.6_wp-beta)**(-0.85_wp)*dbeta+dy/y**2
    dtt=-1.5_wp*sqrt(1.0_wp-x)*dx
    dcc=0.2_wp*(-0.5_wp*dr/sqrt(r)*qq*tt+(1.0_wp-sqrt(r))*(dqq*tt+qq*dtt))
    cd_der=(1.0_wp-weight)*daa+weight*(dbb+dcc)+(bb+cc-aa)*dweight
    ud_der=-0.5_wp*Ud*dr/(1.0_wp-r)
    ud_der(1)=ud_der(1)-0.5_wp*Ud/(1.0_wp-alpha)
    if(la<0.025_wp)ud_der=ud_der-0.45_wp*Ud*ddelta/delta
    if(.not.all(ieee_is_finite(cd_der)).or..not.all(ieee_is_finite(ud_der)))then
      status=CLOSURE_DERIVATIVE;return
    endif
    if(present(dCd))dCd=cd_der
    if(present(dUd))dUd=ud_der
    status=CLOSURE_OK
  end subroutine
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
  subroutine residual_si(alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug,sigma,model,value,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,D,area,rough,G,mul,mug,sigma
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
      if(is_chexal(model))then
        call chexal_residual_si(alpha,rhol,rhog,(1.0_wp-alpha)*ul,alpha*ug,theta,D,G,mul,mug,sigma, &
          value,status,branch=chexal_branch_id(model))
        return
      else if(model==BHAGWAT)then
        call bhagwat_coefficients_si(alpha,rhol,rhog,ul,ug,theta,D,rough,G,mul,sigma,Cd,Ud,status)
      else
        call drift_coefficients_si(alpha,rhol,rhog,ul,ug,theta,D,G,model,Cd,Ud,status)
      endif
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
    call residual_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,AREA,EPS,G,MUL,MUG,sigma,model,value,ierr)
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
    real(wp)::der(5),f0,cd,ud,fr,den,j,dcd(5),dud(5),dj(5),rsi(5)
    integer::model,ierr
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
    if(is_chexal(model))then
      rsi=[alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u]
      call chexal_derivatives_si(rsi,theta,DH,G,MUL,MUG,sigma,chexal_branch_id(model),der,ierr)
      if(present(status))status=ierr
      if(ierr/=0)return
      dfda=der(1)/w_u;dfdrhol=der(2)*w_rho/w_u;dfdrhog=der(3)*w_rho/w_u
      dfdul=der(4);dfdug=der(5)
      return
    endif
    if(model==BHAGWAT)then
      rsi=[alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u]
      call bhagwat_coefficients_si(rsi(1),rsi(2),rsi(3),rsi(4),rsi(5),theta,DH,EPS,G,MUL,sigma, &
        cd,ud,ierr,dcd,dud)
      if(ierr/=0)then
        if(present(status))status=ierr
        return
      endif
      j=(1.0_wp-alpha)*rsi(4)+alpha*rsi(5)
      dj=[rsi(5)-rsi(4),0.0_wp,0.0_wp,1.0_wp-alpha,alpha]
      der=alpha*(j*dcd+cd*dj+dud)
      der(1)=der(1)-rsi(5)+cd*j+ud
      der(5)=der(5)-alpha
      dfda=der(1)/w_u;dfdrhol=der(2)*w_rho/w_u;dfdrhog=der(3)*w_rho/w_u
      dfdul=der(4);dfdug=der(5)
      if(present(status))status=CLOSURE_OK
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
    ! Analytic derivatives of the same SI operator, in its selected Reynolds
    ! branches. Central differences near Re=2000 mixed incompatible slopes.
    call stratified_derivatives_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,AREA,EPS,G,MUL,MUG,der,ierr)
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    den=w_rho*w_u*w_u
    dfda=der(1)/den;dfdrhol=der(2)*w_rho/den;dfdrhog=der(3)*w_rho/den
    dfdul=der(4)*w_u/den;dfdug=der(5)*w_u/den
    if(present(status))status=CLOSURE_OK
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
    if(is_chexal(model))then
      if(.not.all(ieee_is_finite([AREA,MUG,tol])).or.AREA<=0.0_wp.or.MUG<=0.0_wp.or.tol<=0.0_wp)return
      if(DH<=0.0_wp)return
      if(abs(AREA-PI*DH**2/4.0_wp)>1e-6_wp*AREA)return
      call chexal_coefficients_si(alpha,rhol*w_rho,rhog*w_rho,(1.0_wp-alpha)*ul*w_u,alpha*ug*w_u, &
        theta,DH,G,MUL,MUG,sigma,Cd,Ud,ierr,branch=chexal_branch_id(model))
    else if(model==BHAGWAT)then
      if(.not.all(ieee_is_finite([AREA,MUG,tol])).or.AREA<=0.0_wp.or.MUG<=0.0_wp.or.tol<=0.0_wp)return
      if(DH<=0.0_wp)return
      if(abs(AREA-PI*DH**2/4.0_wp)>1e-6_wp*AREA)return
      call bhagwat_coefficients_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,EPS,G,MUL,sigma,Cd,Ud,ierr)
    else
      call drift_coefficients_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,G,model,Cd,Ud,ierr)
    endif
    if(present(status))status=ierr
  end subroutine
  subroutine CdUd_bhagwat(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,status)
    real(wp),intent(in)::alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol
    real(wp),intent(out)::Cd,Ud
    integer,intent(out),optional::status
    call CdUd_drift(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,Cd,Ud,BHAGWAT,status)
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
    real(wp)::cd,ud,dc(5),du(5)
    integer::ierr
    dC0da=invalid_real();dC0drhol=dC0da;dC0drhog=dC0da;dC0dul=dC0da;dC0dug=dC0da
    dUdda=dC0da;dUddrhol=dC0da;dUddrhog=dC0da;dUddul=dC0da;dUddug=dC0da
    call CdUd_drift(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma,w_u,w_rho,tol,cd,ud,BHAGWAT,ierr)
    if(ierr==0)call bhagwat_coefficients_si(alpha,rhol*w_rho,rhog*w_rho,ul*w_u,ug*w_u,theta,DH,EPS,G,MUL,sigma, &
      cd,ud,ierr,dc,du)
    if(ierr==0)then
      dC0da=dc(1);dC0drhol=dc(2)*w_rho;dC0drhog=dc(3)*w_rho;dC0dul=dc(4)*w_u;dC0dug=dc(5)*w_u
      ! Ud is dimensional in this public API; derivatives are with respect to
      ! the nondimensional arguments, so the chain factors are retained.
      dUdda=du(1);dUddrhol=du(2)*w_rho;dUddrhog=du(3)*w_rho;dUddul=du(4)*w_u;dUddug=du(5)*w_u
    endif
    if(present(status))status=ierr
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
  subroutine chexal_derivatives_si(x,theta,D,G,mul,mug,sigma,branch,der,status)
    !! Fourth-order differences of ONE selected residual. A velocity very
    !! close to zero uses a one-sided stencil directed away from the sign
    !! boundary. A local two-sided signature probe first excludes real kinks.
    !! Exact zero is not assigned an invented bilateral derivative.
    real(wp),intent(in)::x(5),theta,D,G,mul,mug,sigma
    integer,intent(in)::branch
    real(wp),intent(out)::der(5)
    integer,intent(out)::status
    real(wp)::xp(5),h,scale,values(6),f0,dc,df,offsets(6),probe,sgn
    integer::sig0,sig,code,k,n,attempt,np,guard_attempt
    logical::valid,one_sided
    der=invalid_real();status=CLOSURE_DERIVATIVE
    call eval(x,f0,code,sig0)
    if(code/=0)then
      status=code;return
    endif
    do k=1,5
      scale=max(abs(x(k)),1.0_wp)
      if(k==1)scale=min(x(1),1.0_wp-x(1))
      h=epsilon(1.0_wp)**0.2_wp*scale
      one_sided=k>=4.and.x(k)/=0.0_wp.and.abs(x(k))<4*h
      sgn=1.0_wp;np=4
      offsets(:4)=[1.0_wp,-1.0_wp,0.5_wp,-0.5_wp]
      if(one_sided)then
        ! A nearby switch need not be a kink AT the state. Reduce the
        ! bilateral probe until both samples use its piece; never accept
        ! an unresolved switch or an unrepresentable probe. The finite
        ! limit keeps exact switches well above rounding-width plateaus.
        probe=0.125_wp*abs(x(k))
        do guard_attempt=1,16
          valid=.true.
          do n=-1,1,2
            xp=x;xp(k)=x(k)+n*probe
            if(xp(k)==x(k))return
            call eval(xp,values(1),code,sig)
            if(code/=0.or.sig/=sig0)valid=.false.
          enddo
          if(valid)exit
          probe=0.5_wp*probe
        enddo
        if(guard_attempt>16)return
        sgn=sign(1.0_wp,x(k));np=6
        offsets=sgn*[0.5_wp,1.0_wp,1.5_wp,2.0_wp,3.0_wp,4.0_wp]
      endif
      do attempt=1,16
        valid=.true.
        do n=1,np
          xp=x;xp(k)=x(k)+offsets(n)*h
          if(xp(k)==x(k))then
            valid=.false.;exit
          endif
          call eval(xp,values(n),code,sig)
          if(code/=0.or.sig/=sig0)then
            valid=.false.;exit
          endif
        enddo
        if(valid)then
          if(one_sided)then
            dc=(-25*f0+48*values(2)-36*values(4)+16*values(5)-3*values(6))/(12*h*sgn)
            df=(-25*f0+48*values(1)-36*values(2)+16*values(3)-3*values(4))/(6*h*sgn)
          else
            dc=(values(1)-values(2))/(2.0_wp*h)
            df=(values(3)-values(4))/h
          endif
          if(abs(dc-df)<=2e-6_wp*max(abs(df),abs(f0)/scale,1e-9_wp))exit
        endif
        h=0.5_wp*h
      enddo
      if(attempt>16)return
      if(one_sided)then
        der(k)=(16.0_wp*df-dc)/15.0_wp
      else
        der(k)=(4.0_wp*df-dc)/3.0_wp
      endif
    enddo
    if(.not.all(ieee_is_finite(der)))return
    status=CLOSURE_OK
  contains
    subroutine eval(q,f,code,sig)
      real(wp),intent(in)::q(5)
      real(wp),intent(out)::f
      integer,intent(out)::code,sig
      call chexal_residual_si(q(1),q(2),q(3),(1.0_wp-q(1))*q(4),q(1)*q(5), &
        theta,D,G,mul,mug,sigma,f,code,branch=branch,signature=sig)
    end subroutine
  end subroutine
end module module_closure_law
