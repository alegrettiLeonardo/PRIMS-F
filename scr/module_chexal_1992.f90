module module_chexal_1992
  !! Chexal et al., Progress in Nuclear Energy 27 (1992), pp.255-295.
  !! Air-water ONLY; SI input. The source's angle is degrees FROM VERTICAL.
  !! Countercurrent use is vertical only (p.262); it needs an explicit branch
  !! (p.268, eqs.3-39--41). No 1997 coefficients or inclined CC extension.
  use mod_kinds, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  integer,parameter,public :: CHEXAL_OK=0,CHEXAL_INVALID=60,CHEXAL_DOMAIN=61
  integer,parameter,public :: CHEXAL_BRANCH_REQUIRED=62,CHEXAL_NO_ROOT=63,CHEXAL_CCFL_EXCEEDED=64
  integer,parameter,public :: CHEXAL_SINGULAR=65,CHEXAL_NUMERIC=66,CHEXAL_SOURCE_UNRESOLVED=67
  integer,parameter,public :: CHEXAL_BRANCH_AUTO=0,CHEXAL_BRANCH_HIGH=1,CHEXAL_BRANCH_LOW=-1
  ! Fluid parameter L (p.262-265). AIR_WATER (default, the only one PRIMS-F
  ! selects): Eqs.3-7/3-8/3-19. STEAM_WATER (verification against the EPRI
  ! sample problems only): Eq.3-5/3-6, all orientations, needs 0<p<p_crit
  ! in any consistent absolute unit. Refrigerant L is not implemented.
  integer,parameter,public :: CHEXAL_FLUID_AIR_WATER=0,CHEXAL_FLUID_STEAM_WATER=1
  real(wp),parameter :: PI=acos(-1.0_wp),ANGLE_TOL=32*epsilon(1.0_wp)
  public :: chexal_coefficients_si,chexal_residual_si,chexal_ccfl
  public :: chexal_solve_alpha,chexal_recover_velocities
  type :: parameters
    real(wp) :: rl,rg,d,g,ml,mg,sigma,theta
    integer :: fluid=CHEXAL_FLUID_AIR_WATER
    real(wp) :: c1=0.0_wp
  end type
contains
  pure real(wp) function nan() result(x)
    x=ieee_value(0.0_wp,ieee_quiet_nan)
  end function

  logical function valid(p) result(ok)
    type(parameters),intent(in)::p
    ok=.false.
    if(.not.all(ieee_is_finite([p%rl,p%rg,p%d,p%g,p%ml,p%mg,p%sigma,p%theta])))return
    ok=p%rl>p%rg.and.p%rg>0.and.p%d>0.and.p%g>=0.and.p%ml>0.and.p%mg>0.and.p%sigma>=0 &
      .and.abs(p%theta)<=PI/2+ANGLE_TOL
  end function

  subroutine set_fluid(p,fluid,pressure,pcrit,status)
    !! Eq.3-6: C1 = 4 p_crit^2/[p (p_crit - p)], absolute pressures, same unit.
    type(parameters),intent(inout)::p
    integer,intent(in),optional::fluid
    real(wp),intent(in),optional::pressure,pcrit
    integer,intent(out)::status
    status=CHEXAL_OK
    if(.not.present(fluid))return
    status=CHEXAL_INVALID
    select case(fluid)
    case(CHEXAL_FLUID_AIR_WATER)
      p%fluid=fluid
    case(CHEXAL_FLUID_STEAM_WATER)
      if(.not.(present(pressure).and.present(pcrit)))return
      if(.not.all(ieee_is_finite([pressure,pcrit])))return
      if(pressure<=0.or.pressure>=pcrit)return
      p%fluid=fluid;p%c1=4*pcrit**2/(pressure*(pcrit-pressure))
    case default
      return
    end select
    status=CHEXAL_OK
  end subroutine

  subroutine source_coordinates(p,jl,jg,fl,fg,theta_deg,direction,status)
    type(parameters),intent(in)::p
    real(wp),intent(in)::jl,jg
    real(wp),intent(out)::fl,fg,theta_deg,direction
    integer,intent(out)::status
    status=CHEXAL_INVALID
    if(.not.valid(p).or..not.all(ieee_is_finite([jl,jg])))return
    direction=1.0_wp
    if(p%theta<0)direction=-1.0_wp
    fl=direction*jl;fg=direction*jg
    theta_deg=90.0_wp-abs(p%theta)*180.0_wp/PI
    ! Exact geometry up to roundoff; never regard the catenary top as vertical.
    if(abs(abs(p%theta)-PI/2)<=ANGLE_TOL)theta_deg=0.0_wp
    status=CHEXAL_DOMAIN
    if(p%theta==0.0_wp.and.(fl<0.or.fg<0))return
    if((fl<0.and.fg>0).or.(fl>0.and.fg<0))then
      if(theta_deg/=0.0_wp.or.fg<0)return
    endif
    if(fg<0.and.fl>0)return
    status=CHEXAL_OK
  end subroutine

  subroutine downward_c3(rel,d,c3,status)
    !! Eq.3-35--38, p.268. Eq.3-36 exp[(|Rel|/350000)^0.4] is clarified
    !! by primary RELAP5/MOD3 (1995), vol4 p.6-10, BEFORE the 1997 variant.
    !! Only that parenthesization is cross-checked; coefficients stay 1992.
    real(wp),intent(in)::rel,d
    real(wp),intent(out)::c3
    integer,intent(out)::status
    real(wp)::rr,arg,b2,c10
    c3=nan();status=CHEXAL_NUMERIC
    rr=abs(rel);arg=(rr/350000.0_wp)**0.4_wp
    if(arg>log(huge(1.0_wp)/4))return
    c10=2*exp(arg)-1.75_wp*rr**0.03_wp*exp(-rr/50000.0_wp*(0.0381_wp/d)**2) &
      +(0.0381_wp/d)**0.25_wp*rr**0.001_wp
    if(c10<=0.or..not.ieee_is_finite(c10))return
    b2=(1/(1+0.05_wp*(rr/350000.0_wp)))**0.4_wp
    c3=2*(c10/2)**b2
    if(.not.ieee_is_finite(c3))return
    status=CHEXAL_OK
  end subroutine

  subroutine component(a,p,fl,fg,horizontal,c3_override,c0,v0,v,status,signature)
    !! One source-coordinate component. Horizontal evaluation uses positive
    !! fluxes, as prescribed in pp.264-265/269, not to extend countercurrent use.
    real(wp),intent(in)::a,fl,fg,c3_override
    type(parameters),intent(in)::p
    logical,intent(in)::horizontal
    real(wp),intent(out)::c0,v0,v
    integer,intent(out)::status,signature
    real(wp)::rel,reg,re,b1,k0,r,l,c2,c3,c4,c5,c7,c9,ratio,alt,x
    c0=nan();v0=c0;v=c0;signature=0;status=CHEXAL_NUMERIC
    rel=p%rl*fl*p%d/p%ml;reg=p%rg*fg*p%d/p%mg
    if(.not.all(ieee_is_finite([rel,reg])))return
    ! Signed Reynolds: Eq.3-15, with the negative gas branch taking precedence.
    if(reg>rel.or.fg<0)then
      re=reg;signature=ibset(signature,0)
    else
      re=rel
    endif
    ! Stable logistic, mathematically identical to Eq.3-14.
    if(re>=0)then
      b1=1/(1+exp(-re/60000.0_wp))
    else
      x=exp(re/60000.0_wp);b1=x/(1+x)
    endif
    if(b1>=0.8_wp)then
      b1=0.8_wp;signature=ibset(signature,1)
    endif
    ratio=p%rg/p%rl;k0=b1+(1-b1)*ratio**0.25_wp
    r=(1+1.57_wp*ratio)/(1-b1)
    if(p%fluid==CHEXAL_FLUID_STEAM_WATER)then
      l=(1-exp(-p%c1*a))/(1-exp(-p%c1))            ! Eq.3-5, C1 of Eq.3-6
    else if(horizontal)then
      l=1.125_wp*a**0.6_wp                         ! Eq.3-19
    else if(fg>=0)then
      l=1.15_wp*a**0.45_wp                         ! Eq.3-7
    else
      l=1.05_wp*a**0.25_wp                         ! Eq.3-8
      signature=ibset(signature,2)
    endif
    if(l>=1)then
      l=1.0_wp;signature=ibset(signature,3)
    endif
    c0=l/(k0+(1-k0)*a**r)                         ! Eq.3-4a
    if(1/ratio<=18)then
      c2=0.4757_wp*log(1/ratio)**0.7_wp           ! Eq.3-27
      signature=ibset(signature,4)
    else
      c5=sqrt(150*ratio)
      c2=1.0_wp
      if(c5<1)then
        c2=1/(1-exp(-c5/(1-c5)))                 ! Eq.3-28--29
        signature=ibset(signature,5)
      endif
    endif
    c7=(0.09144_wp/p%d)**0.6_wp;c4=1.0_wp         ! Eq.3-30--33
    if(c7<1)then
      c4=1/(1-exp(-c7/(1-c7)))
      signature=ibset(signature,6)
    endif
    if(fl>=0.and.fg>=0)then
      c3=2*exp(-abs(rel)/60000.0_wp)              ! Eq.3-34
      if(c3<=0.5_wp)then
        c3=0.5_wp;signature=ibset(signature,7)
      endif
    else
      call downward_c3(rel,p%d,c3,status)
      if(status/=CHEXAL_OK)return
      signature=ibset(signature,8)
    endif
    if(c3_override>=0)then
      c3=c3_override;signature=ibset(signature,9)
    endif
    v0=1.41_wp*((p%rl-p%rg)*p%sigma*p%g/p%rl**2)**0.25_wp*c2*c3*c4 ! Eq.3-24, gc=1 SI
    if(fg>=0)then
      c9=(1-a)**b1                               ! Eq.3-25
    else
      c9=(1-a)**0.65_wp
      if(c9>=0.7_wp)then
        c9=0.7_wp;signature=ibset(signature,10)
      endif
      if(abs(fl)+abs(fg)==0)then
        status=CHEXAL_SINGULAR;return
      endif
      alt=v0*(1-a)**0.2_wp/(abs(fl)+abs(fg))      ! Eq.3-4b
      if(alt>c0)then
        c0=alt;signature=ibset(signature,11)
      endif
    endif
    if(horizontal)c0=(1+a**0.05_wp*(1-a)**2)*c0 ! Eq.3-18
    v=v0*c9                                     ! Eq.3-23
    if(.not.all(ieee_is_finite([c0,v0,v])))then
      status=CHEXAL_NUMERIC;return
    endif
    status=CHEXAL_OK
  end subroutine

  recursive subroutine evaluate(a,p,jl,jg,branch,c0,v,status,signature,flood)
    real(wp),intent(in)::a,jl,jg
    type(parameters),intent(in)::p
    integer,intent(in)::branch
    real(wp),intent(out)::c0,v
    integer,intent(out)::status,signature
    real(wp),intent(in),optional::flood
    real(wp)::fl,fg,theta_deg,direction,fr,cv,ch,vv,vh,v0,c3,c3low,fs,af,rr
    integer::sv,sh
    logical::cc
    c0=nan();v=c0;signature=0;status=CHEXAL_INVALID
    if(.not.ieee_is_finite(a))return
    if(a<0.or.a>1)return
    if(branch/=CHEXAL_BRANCH_AUTO.and.branch/=CHEXAL_BRANCH_HIGH.and.branch/=CHEXAL_BRANCH_LOW)return
    call source_coordinates(p,jl,jg,fl,fg,theta_deg,direction,status)
    if(status/=CHEXAL_OK)return
    cc=fl<0.and.fg>=0
    if(cc.and.theta_deg/=0)then
      status=CHEXAL_DOMAIN;return
    endif
    if(cc.and.branch==CHEXAL_BRANCH_AUTO)then
      status=CHEXAL_BRANCH_REQUIRED;return
    endif
    c3=-1.0_wp
    if(cc.and.branch==CHEXAL_BRANCH_LOW)then
      call downward_c3(p%rl*fl*p%d/p%ml,p%d,c3,status)
      if(status/=CHEXAL_OK)return
      c3low=2*(1+abs(p%rl*fl*p%d/p%ml)/60000.0_wp)
      ! Eq.3-41 is a min of a convex interpolation and C3_HIGH. If its
      ! other endpoint is >= C3_HIGH, it equals C3_HIGH throughout the
      ! sub-CCFL domain, independent of jl*. Do not solve an irrelevant CCFL
      ! problem here. Final roots are still checked against physical CCFL.
      if(c3low<c3)then
        if(fg==0)then
          ! At jg=0 the positive-alpha maximum persists for every finite
          ! downward liquid flux: jl* tends to -infinity, hence jl/jl* -> 0.
          c3=c3low;signature=ibset(signature,26)
        else
        if(present(flood))then
          fs=flood
        else
          call ccfl_internal(p,fg,fs,af,status)
          if(status/=CHEXAL_OK)return
        endif
        if(.not.ieee_is_finite(fs).or.fs>=0)then
          status=CHEXAL_INVALID;return
        endif
        if(fl<fs*(1+1e-10_wp))then
          status=CHEXAL_CCFL_EXCEEDED;return
        endif
        rr=fl/fs
        c3=c3*rr+c3low*(1-rr)
        signature=ibset(signature,26)
        endif
      endif
      signature=ibset(signature,27)
    endif
    call component(a,p,fl,fg,.false.,c3,cv,v0,vv,status,sv)
    if(status/=CHEXAL_OK)return
    signature=ior(signature,sv)
    if(theta_deg==0)then
      c0=cv;v=direction*vv;return
    endif
    ! Vertical CC is already handled; the horizontal component is evaluated
    ! in its own positive-flow convention for cocurrent interpolation only.
    call component(a,p,abs(fl),abs(fg),.true.,-1.0_wp,ch,v0,vh,status,sh)
    if(status/=CHEXAL_OK)return
    signature=ior(signature,ishft(sh,13))
    if(fg>=0)then
      fr=(90-theta_deg)/90
      v=fr*vv+(1-fr)*vh                          ! Eq.3-21
    else
      if(theta_deg<=80)then
        fr=1.0_wp;signature=ibset(signature,28)
      else
        fr=(90-theta_deg)/10
      endif
      v=fr*vv+(fr-1)*vh                          ! Eq.3-22
      signature=ibset(signature,29)
    endif
    c0=fr*cv+(1-fr)*ch;v=direction*v             ! Eq.3-2
  end subroutine

  recursive subroutine chexal_coefficients_si(alpha,rhol,rhog,jl,jg,theta,D,G,mul,mug,sigma,c0,vgj,status, &
      branch,signature,jl_ccfl,fluid,pressure,pcrit)
    real(wp),intent(in)::alpha,rhol,rhog,jl,jg,theta,D,G,mul,mug,sigma
    real(wp),intent(out)::c0,vgj
    integer,intent(out)::status
    integer,intent(in),optional::branch
    integer,intent(out),optional::signature
    real(wp),intent(in),optional::jl_ccfl
    integer,intent(in),optional::fluid
    real(wp),intent(in),optional::pressure,pcrit
    type(parameters)::p
    integer::b,sig
    real(wp)::fl,fg,theta_deg,direction,fs,af
    p=parameters(rhol,rhog,D,G,mul,mug,sigma,theta);b=CHEXAL_BRANCH_AUTO
    if(present(branch))b=branch
    c0=nan();vgj=nan();sig=0;status=CHEXAL_INVALID
    if(present(signature))signature=0
    call set_fluid(p,fluid,pressure,pcrit,status)
    if(status/=CHEXAL_OK)return
    status=CHEXAL_INVALID
    if(.not.ieee_is_finite(alpha))return
    if(alpha<0.or.alpha>1)return
    call source_coordinates(p,jl,jg,fl,fg,theta_deg,direction,status)
    if(status/=CHEXAL_OK)return
    if(fl<0.and.fg>0.and.(b==CHEXAL_BRANCH_HIGH.or.b==CHEXAL_BRANCH_LOW))then
      ! Public coefficients describe an admissible flow pair. Internal
      ! residual trials can use the algebraic continuation to locate roots;
      ! those trials are never themselves accepted as physical solutions.
      if(present(jl_ccfl))then
        fs=jl_ccfl
        if(.not.ieee_is_finite(fs).or.fs>=0)then
          status=CHEXAL_INVALID;return
        endif
      else
        call ccfl_internal(p,fg,fs,af,status)
        if(status/=CHEXAL_OK)return
      endif
      if(fl<fs*(1+1e-10_wp))then
        status=CHEXAL_CCFL_EXCEEDED;return
      endif
      call evaluate(alpha,p,jl,jg,b,c0,vgj,status,sig,fs)
    else
      call evaluate(alpha,p,jl,jg,b,c0,vgj,status,sig,jl_ccfl)
    endif
    if(present(signature))signature=sig
  end subroutine

  subroutine chexal_residual_si(alpha,rhol,rhog,jl,jg,theta,D,G,mul,mug,sigma,residual,status, &
      branch,signature,jl_ccfl,fluid,pressure,pcrit)
    real(wp),intent(in)::alpha,rhol,rhog,jl,jg,theta,D,G,mul,mug,sigma
    real(wp),intent(out)::residual
    integer,intent(out)::status
    integer,intent(in),optional::branch
    integer,intent(out),optional::signature
    real(wp),intent(in),optional::jl_ccfl
    integer,intent(in),optional::fluid
    real(wp),intent(in),optional::pressure,pcrit
    real(wp)::c0,v
    call chexal_coefficients_si(alpha,rhol,rhog,jl,jg,theta,D,G,mul,mug,sigma,c0,v,status, &
      branch,signature,jl_ccfl,fluid,pressure,pcrit)
    residual=nan()
    if(status==CHEXAL_OK)residual=alpha*(c0*(jl+jg)+v)-jg
  end subroutine

  recursive real(wp) function residual(a,p,jl,jg,b,status,flood) result(f)
    real(wp),intent(in)::a,jl,jg
    type(parameters),intent(in)::p
    integer,intent(in)::b
    integer,intent(out)::status
    real(wp),intent(in),optional::flood
    real(wp)::c,v
    integer::sig
    call evaluate(a,p,jl,jg,b,c,v,status,sig,flood)
    f=nan()
    if(status==CHEXAL_OK)f=a*(c*(jl+jg)+v)-jg
  end function

  subroutine max_residual(p,jl,jg,amax,fmax,status)
    !! Locate the first-fold residual maximum, including mesh-local maxima.
    !! Golden section refines every sampled maximum, not only the largest bin.
    type(parameters),intent(in)::p
    real(wp),intent(in)::jl,jg
    real(wp),intent(out)::amax,fmax
    integer,intent(out)::status
    real(wp)::xs(0:128),ys(0:128),lo,hi,x1,x2,y1,y2
    real(wp),parameter::gr=0.6180339887498948482_wp
    integer::i,k
    status=CHEXAL_OK;amax=nan();fmax=-huge(1.0_wp)
    do i=0,128
      xs(i)=0.5_wp*(1-cos(PI*real(i,wp)/128))
      ys(i)=residual(xs(i),p,jl,jg,CHEXAL_BRANCH_HIGH,status)
      if(status/=CHEXAL_OK)return
      if(ys(i)>fmax)then
        fmax=ys(i);amax=xs(i)
      endif
    enddo
    do i=1,127
      if(ys(i)<ys(i-1).or.ys(i)<ys(i+1))cycle
      lo=xs(i-1);hi=xs(i+1);x1=hi-gr*(hi-lo);x2=lo+gr*(hi-lo)
      y1=residual(x1,p,jl,jg,CHEXAL_BRANCH_HIGH,status)
      if(status/=CHEXAL_OK)return
      y2=residual(x2,p,jl,jg,CHEXAL_BRANCH_HIGH,status)
      if(status/=CHEXAL_OK)return
      do k=1,64
        if(y1<y2)then
          lo=x1;x1=x2;y1=y2;x2=lo+gr*(hi-lo)
          y2=residual(x2,p,jl,jg,CHEXAL_BRANCH_HIGH,status)
        else
          hi=x2;x2=x1;y2=y1;x1=hi-gr*(hi-lo)
          y1=residual(x1,p,jl,jg,CHEXAL_BRANCH_HIGH,status)
        endif
        if(status/=CHEXAL_OK)return
      enddo
      if(y1>fmax)then
        fmax=y1;amax=x1
      endif
      if(y2>fmax)then
        fmax=y2;amax=x2
      endif
    enddo
  end subroutine

  subroutine ccfl_internal(p,jg,jlstar,astar,status)
    type(parameters),intent(in)::p
    real(wp),intent(in)::jg
    real(wp),intent(out)::jlstar,astar
    integer,intent(out)::status
    type(parameters)::pv
    real(wp)::lo,hi,mid,fl,fh,fm,a,v,c,scale,rscale,ftol,delta,fminus,fplus
    integer::i,sig
    jlstar=nan();astar=nan();status=CHEXAL_INVALID
    if(.not.valid(p).or..not.ieee_is_finite(jg))return
    if(jg<=0)return
    pv=p;pv%theta=PI/2
    call evaluate(0.5_wp,pv,0.0_wp,jg,CHEXAL_BRANCH_HIGH,c,v,status,sig)
    if(status/=CHEXAL_OK)return
    scale=max(jg,abs(v),tiny(1.0_wp))
    lo=0;hi=-0.01_wp*scale
    call max_residual(pv,lo,jg,a,fl,status)
    if(status/=CHEXAL_OK)return
    if(fl<=0)then
      status=CHEXAL_NO_ROOT;return
    endif
    do i=1,80
      call max_residual(pv,hi,jg,a,fh,status)
      if(status/=CHEXAL_OK)return
      if(fh<=0)exit
      lo=hi;hi=2*hi
    enddo
    if(i>80)then
      status=CHEXAL_NUMERIC;return
    endif
    do i=1,70
      mid=0.5_wp*(lo+hi)
      call max_residual(pv,mid,jg,a,fm,status)
      if(status/=CHEXAL_OK)return
      if(fm>0)then
        lo=mid
      else
        hi=mid
      endif
      if(abs(hi-lo)<=32*epsilon(1.0_wp)*max(abs(mid),scale))exit
    enddo
    jlstar=0.5_wp*(lo+hi)
    call max_residual(pv,jlstar,jg,astar,fm,status)
    if(status/=CHEXAL_OK)return
    ! A disappearing sampled maximum can merely fall BETWEEN the mesh nodes.
    ! Never report that sampling event as a physical CCFL fold. Check the
    ! candidate's interiority, actual residual and local maximality first.
    if(astar<=0.0_wp.or.astar>=1.0_wp)then
      jlstar=nan();astar=nan();status=CHEXAL_NUMERIC;return
    endif
    call evaluate(astar,pv,jlstar,jg,CHEXAL_BRANCH_HIGH,c,v,status,sig)
    if(status/=CHEXAL_OK)return
    rscale=max(abs(astar*c*(jlstar+jg)),abs(astar*v),abs(jg))
    ftol=2048*epsilon(1.0_wp)*rscale
    if(abs(fm)>ftol)then
      jlstar=nan();astar=nan();status=CHEXAL_NUMERIC;return
    endif
    delta=epsilon(1.0_wp)**(1.0_wp/3.0_wp)*min(astar,1-astar)
    if(delta<=0.or.astar-delta==astar.or.astar+delta==astar)then
      jlstar=nan();astar=nan();status=CHEXAL_SINGULAR;return
    endif
    fminus=residual(astar-delta,pv,jlstar,jg,CHEXAL_BRANCH_HIGH,status)
    if(status/=CHEXAL_OK)return
    fplus=residual(astar+delta,pv,jlstar,jg,CHEXAL_BRANCH_HIGH,status)
    if(status/=CHEXAL_OK)return
    ! One-sided inequalities also permit a genuine maximum at a constitutive
    ! kink; demanding a zero bilateral derivative there would be incorrect.
    if(fminus>fm+ftol.or.fplus>fm+ftol)then
      jlstar=nan();astar=nan();status=CHEXAL_NUMERIC;return
    endif
  end subroutine

  subroutine chexal_ccfl(jg,rhol,rhog,D,G,mul,mug,sigma,jl_ccfl,alpha_ccfl,status)
    real(wp),intent(in)::jg,rhol,rhog,D,G,mul,mug,sigma
    real(wp),intent(out)::jl_ccfl,alpha_ccfl
    integer,intent(out)::status
    type(parameters)::p
    p=parameters(rhol,rhog,D,G,mul,mug,sigma,PI/2)
    call ccfl_internal(p,jg,jl_ccfl,alpha_ccfl,status)
  end subroutine

  subroutine bracket_root(p,jl,jg,b,lo_in,hi_in,root,status,flood)
    type(parameters),intent(in)::p
    real(wp),intent(in)::jl,jg,lo_in,hi_in
    integer,intent(in)::b
    real(wp),intent(out)::root
    integer,intent(out)::status
    real(wp),intent(in),optional::flood
    real(wp)::lo,hi,fl,fh,x,fx
    integer::i
    root=nan();lo=lo_in;hi=hi_in
    fl=residual(lo,p,jl,jg,b,status,flood);if(status/=CHEXAL_OK)return
    fh=residual(hi,p,jl,jg,b,status,flood);if(status/=CHEXAL_OK)return
    if(fl==0)then
      root=lo;return
    endif
    if(fh==0)then
      root=hi;return
    endif
    if(sign(1.0_wp,fl)==sign(1.0_wp,fh))then
      status=CHEXAL_NO_ROOT;return
    endif
    do i=1,160
      x=lo+(hi-lo)/2
      if(x==lo.or.x==hi)exit
      fx=residual(x,p,jl,jg,b,status,flood);if(status/=CHEXAL_OK)return
      if(fx==0)then
        lo=x;hi=x;exit
      endif
      if(sign(1.0_wp,fx)==sign(1.0_wp,fl))then
        lo=x;fl=fx
      else
        hi=x;fh=fx
      endif
    enddo
    root=lo+(hi-lo)/2;status=CHEXAL_OK
  end subroutine

  subroutine add_root(x,rr,n,status)
    real(wp),intent(in)::x
    real(wp),intent(inout)::rr(:)
    integer,intent(inout)::n
    integer,intent(out)::status
    integer::i
    status=CHEXAL_OK
    do i=1,n
      if(abs(rr(i)-x)<=64*epsilon(1.0_wp)*max(abs(x),abs(rr(i)),tiny(1.0_wp)))return
    enddo
    if(n==size(rr))then
      status=CHEXAL_NUMERIC;return
    endif
    n=n+1;rr(n)=x
  end subroutine

  subroutine enumerate_alpha(p,jl,jg,b,rr,n,status,flood)
    type(parameters),intent(in)::p
    real(wp),intent(in)::jl,jg
    integer,intent(in)::b
    real(wp),intent(out)::rr(8)
    integer,intent(out)::n,status
    real(wp),intent(in),optional::flood
    real(wp)::xx(0:1010),x,fl,fh,r,amax,fmax,tmp
    integer::i,k,count
    n=0;rr=nan();status=CHEXAL_OK
    ! Linear interior plus representable positive tails: no alpha clipping.
    count=0;xx(0)=0
    do i=1,256
      count=count+1;xx(count)=real(i,wp)/256
    enddo
    do i=1,708
      x=exp(-real(i,wp))
      count=count+1;xx(count)=x
    enddo
    do i=1,36
      x=1-exp(-real(i,wp))
      count=count+1;xx(count)=x
    enddo
    if(jl<0.and.jg>0.and.b==CHEXAL_BRANCH_HIGH)then
      call max_residual(p,jl,jg,amax,fmax,status)
      if(status/=CHEXAL_OK)return
      count=count+1;xx(count)=amax
    endif
    do i=1,count
      tmp=xx(i);k=i-1
      do while(k>=0)
        if(xx(k)<=tmp)exit
        xx(k+1)=xx(k);k=k-1
      enddo
      xx(k+1)=tmp
    enddo
    fl=residual(xx(0),p,jl,jg,b,status,flood);if(status/=CHEXAL_OK)return
    if(fl==0)call add_root(xx(0),rr,n,status)
    do i=1,count
      if(xx(i)==xx(i-1))cycle
      fh=residual(xx(i),p,jl,jg,b,status,flood);if(status/=CHEXAL_OK)return
      if(fh==0)then
        call add_root(xx(i),rr,n,status)
      else if(fl/=0.and.sign(1.0_wp,fl)/=sign(1.0_wp,fh))then
        call bracket_root(p,jl,jg,b,xx(i-1),xx(i),r,status,flood)
        if(status/=CHEXAL_OK)return
        call add_root(r,rr,n,status)
      endif
      if(status/=CHEXAL_OK)return
      fl=fh
    enddo
    if(n==0)status=CHEXAL_NO_ROOT
  end subroutine

  subroutine chexal_solve_alpha(jl,jg,rhol,rhog,theta,D,G,mul,mug,sigma,alpha,status,branch,nroots,alpha_hint,roots, &
      fluid,pressure,pcrit)
    real(wp),intent(in)::jl,jg,rhol,rhog,theta,D,G,mul,mug,sigma
    real(wp),intent(out)::alpha
    integer,intent(out)::status
    integer,intent(in),optional::branch
    integer,intent(out),optional::nroots
    real(wp),intent(in),optional::alpha_hint
    real(wp),intent(out),optional::roots(:)
    integer,intent(in),optional::fluid
    real(wp),intent(in),optional::pressure,pcrit
    type(parameters)::p,pv
    real(wp)::fl,fg,angle,direction,rr(8),rlow(8),fs,af
    integer::b,n,nl,k,idx
    alpha=nan();status=CHEXAL_INVALID
    if(present(nroots))nroots=0
    if(present(roots))roots=nan()
    p=parameters(rhol,rhog,D,G,mul,mug,sigma,theta);b=CHEXAL_BRANCH_AUTO
    call set_fluid(p,fluid,pressure,pcrit,status)
    if(status/=CHEXAL_OK)return
    status=CHEXAL_INVALID
    if(present(branch))b=branch
    if(b/=CHEXAL_BRANCH_AUTO.and.b/=CHEXAL_BRANCH_HIGH.and.b/=CHEXAL_BRANCH_LOW)return
    if(present(alpha_hint))then
      if(.not.ieee_is_finite(alpha_hint))return
      if(alpha_hint<0.or.alpha_hint>1)return
    endif
    call source_coordinates(p,jl,jg,fl,fg,angle,direction,status)
    if(status/=CHEXAL_OK)return
    if(fl<0.and.fg>0)then
      pv=p;pv%theta=PI/2
      call ccfl_internal(pv,fg,fs,af,status)
      if(status/=CHEXAL_OK)return
      if(fl<fs*(1+1e-10_wp))then
        status=CHEXAL_CCFL_EXCEEDED;return
      endif
      if(abs(fl-fs)<=1e-10_wp*abs(fs))then
        alpha=af;status=CHEXAL_SINGULAR
        if(present(nroots))nroots=1
        if(present(roots))then
          if(size(roots)>0)roots(1)=af
        endif
        return
      endif
      call enumerate_alpha(pv,fl,fg,CHEXAL_BRANCH_HIGH,rr,n,status,fs)
      if(status/=CHEXAL_OK)return
      call enumerate_alpha(pv,fl,fg,CHEXAL_BRANCH_LOW,rlow,nl,status,fs)
      if(status/=CHEXAL_OK)return
      rr(2)=maxval(rr(:n));rr(1)=minval(rlow(:nl));n=2
      if(b==CHEXAL_BRANCH_HIGH)then
        alpha=rr(2)
      else if(b==CHEXAL_BRANCH_LOW)then
        alpha=rr(1)
      else
        status=CHEXAL_BRANCH_REQUIRED
      endif
    else
      call enumerate_alpha(p,jl,jg,b,rr,n,status)
      if(status/=CHEXAL_OK)return
      if(n==1)then
        alpha=rr(1)
      else if(present(alpha_hint))then
        if(.not.ieee_is_finite(alpha_hint))then
          status=CHEXAL_INVALID;return
        endif
        idx=minloc(abs(rr(:n)-alpha_hint),dim=1);alpha=rr(idx)
      else
        status=CHEXAL_BRANCH_REQUIRED
      endif
    endif
    if(present(nroots))nroots=n
    if(present(roots))then
      k=min(size(roots),n);roots(:k)=rr(:k)
    endif
  end subroutine

  subroutine chexal_recover_velocities(alpha,rhol,rhog,momentum,theta,D,G,mul,mug,sigma,ul,ug,status, &
      branch,ul_hint,ug_hint)
    real(wp),intent(in)::alpha,rhol,rhog,momentum,theta,D,G,mul,mug,sigma
    real(wp),intent(out)::ul,ug
    integer,intent(out)::status
    integer,intent(in),optional::branch
    real(wp),intent(in),optional::ul_hint,ug_hint
    type(parameters)::p
    real(wp)::lo,hi,x,jl,fl,fh,fx,direction,scale
    integer::b,i,ncc
    logical::have_prior
    ul=nan();ug=nan();status=CHEXAL_INVALID
    p=parameters(rhol,rhog,D,G,mul,mug,sigma,theta);b=CHEXAL_BRANCH_AUTO
    if(present(branch))b=branch
    if(b/=CHEXAL_BRANCH_AUTO.and.b/=CHEXAL_BRANCH_HIGH.and.b/=CHEXAL_BRANCH_LOW)return
    if(.not.valid(p).or..not.all(ieee_is_finite([alpha,momentum])))return
    if(alpha<0.or.alpha>1)return
    ! Zero velocity for an absent phase is storage convention, not a root.
    if(alpha==0)then
      ul=momentum/rhol;ug=0;status=CHEXAL_OK;return
    else if(alpha==1)then
      ul=0;ug=momentum/rhog;status=CHEXAL_OK;return
    endif
    direction=1;if(theta<0)direction=-1
    have_prior=.false.
    ! Continuation information takes precedence: never silently replace a CC
    ! prior by a cocurrent root of the same conservative state.
    if(present(ul_hint).and.present(ug_hint))then
      if(.not.all(ieee_is_finite([ul_hint,ug_hint])))then
        status=CHEXAL_INVALID;return
      endif
      have_prior=.true.
      if(direction*ul_hint<0.and.direction*ug_hint>0)then
        call recover_countercurrent(alpha,p,momentum,b,ug_hint,ul,ug,status)
        return
      endif
    endif
    if(abs(abs(theta)-PI/2)<=ANGLE_TOL.and.b/=CHEXAL_BRANCH_AUTO)then
      if(.not.have_prior)then
        call enumerate_cc_velocity(alpha,p,momentum,b,ul,ug,ncc,status)
        if(status/=CHEXAL_NO_ROOT)return
        ! No admissible CC root was detected in the derived bounded interval;
        ! finite sampling is not a global certificate of root completeness.
      endif
    endif
    if(abs(abs(theta)-PI/2)<=ANGLE_TOL.and.b==CHEXAL_BRANCH_AUTO.and..not.have_prior)then
      ! AUTO tests cocurrent inversion only after both source-branch searches
      ! detect no vertical CC alternative. Their finite sampling is documented.
      call enumerate_cc_velocity(alpha,p,momentum,CHEXAL_BRANCH_HIGH,ul,ug,ncc,status)
      if(status/=CHEXAL_NO_ROOT)then
        if(status==CHEXAL_OK)status=CHEXAL_BRANCH_REQUIRED
        ul=nan();ug=nan();return
      endif
      call enumerate_cc_velocity(alpha,p,momentum,CHEXAL_BRANCH_LOW,ul,ug,ncc,status)
      if(status/=CHEXAL_NO_ROOT)then
        if(status==CHEXAL_OK)status=CHEXAL_BRANCH_REQUIRED
        ul=nan();ug=nan();return
      endif
    endif
    lo=min(0.0_wp,momentum/rhog);hi=max(0.0_wp,momentum/rhog)
    ! Evaluate one-sided cocurrent endpoint limits in the representable
    ! interval. In particular jg=0 with descending liquid is NOT evidence for
    ! an inclined countercurrent state. No velocity floor is imposed.
    if(momentum/=0)then
      lo=nearest(lo,1.0_wp);hi=nearest(hi,-1.0_wp)
    endif
    fl=residual(alpha,p,(momentum-rhog*lo)/rhol,lo,b,status)
    if(status/=CHEXAL_OK)return
    fh=residual(alpha,p,(momentum-rhog*hi)/rhol,hi,b,status)
    if(status/=CHEXAL_OK)return
    if(fl/=0.and.fh/=0.and.sign(1.0_wp,fl)==sign(1.0_wp,fh))then
      ! A missing cocurrent solution must not be continued into an inclined
      ! countercurrent state. The source permits CC only in a vertical tube.
      if(abs(abs(theta)-PI/2)>ANGLE_TOL)then
        status=CHEXAL_DOMAIN;return
      endif
      if(b==CHEXAL_BRANCH_AUTO)then
        status=CHEXAL_BRANCH_REQUIRED;return
      endif
      status=CHEXAL_BRANCH_REQUIRED;return
    endif
    do i=1,160
      if(fl==0)then
        hi=lo;exit
      else if(fh==0)then
        lo=hi;exit
      endif
      x=lo+(hi-lo)/2
      if(mod(i,4)/=0.and.fh/=fl)then
        x=lo-fl/(fh-fl)*(hi-lo)
        if(x<=lo.or.x>=hi)x=lo+(hi-lo)/2
      endif
      if(x==lo.or.x==hi)exit
      jl=(momentum-rhog*x)/rhol
      fx=residual(alpha,p,jl,x,b,status)
      if(status/=CHEXAL_OK)return
      scale=max(abs(x),abs(jl),tiny(1.0_wp))
      if(abs(fx)<=16*epsilon(1.0_wp)*scale)then
        lo=x;hi=x;exit
      endif
      if(fx==0)then
        lo=x;hi=x;exit
      else if(sign(1.0_wp,fx)==sign(1.0_wp,fl))then
        lo=x;fl=fx
      else
        hi=x;fh=fx
      endif
    enddo
    x=lo+(hi-lo)/2;jl=(momentum-rhog*x)/rhol
    ug=x/alpha;ul=jl/(1-alpha)
    status=CHEXAL_OK
  end subroutine

  subroutine enumerate_cc_velocity(a,p,momentum,b,ul,ug,nroots,status)
    !! No arbitrary velocity cut-off: an analytic bound follows from eqs.3-36,
    !! 3-37 and the maximum C0 at B1=1/2 for positive gas Reynolds. See audit doc.
    real(wp),intent(in)::a,momentum
    type(parameters),intent(in)::p
    integer,intent(in)::b
    real(wp),intent(out)::ul,ug
    integer,intent(out)::nroots,status
    type(parameters)::pv
    real(wp)::direction,q,ratio,ko,power,l,cmax,e,c3max,vmax,c2,c4,c5,c7,den,lower,upper
    real(wp)::xx(0:320),ff(0:320),t,prev,fprev,x,fx,lo,hi,fl,fh,mid,fmid,jl,asel,rr(8)
    integer::i,j,k,n,st,sm,nn
    logical::had_previous
    ul=nan();ug=nan();nroots=0;status=CHEXAL_INVALID
    if(b==CHEXAL_BRANCH_AUTO)return
    pv=p;pv%theta=PI/2;direction=1;if(p%theta<0)direction=-1
    q=direction*momentum;ratio=p%rg/p%rl
    ko=.5_wp+.5_wp*ratio**.25_wp;power=2*(1+1.57_wp*ratio)
    l=min(1.15_wp*a**.45_wp,1.0_wp);cmax=l/(ko+(1-ko)*a**power)
    e=.5_wp*(.0381_wp/p%d)**.25_wp*350000.0_wp**.001_wp*.0025_wp**.0025_wp*exp(-.0025_wp)
    c3max=2*exp(20.0_wp**.4_wp)*(1+e)
    if(1/ratio<=18)then
      c2=.4757_wp*log(1/ratio)**.7_wp
    else
      c5=sqrt(150*ratio);c2=1
      if(c5<1)c2=1/(1-exp(-c5/(1-c5)))
    endif
    c7=(.09144_wp/p%d)**.6_wp;c4=1
    if(c7<1)c4=1/(1-exp(-c7/(1-c7)))
    vmax=1.41_wp*((p%rl-p%rg)*p%sigma*p%g/p%rl**2)**.25_wp*c2*c3max*c4*sqrt(1-a)
    den=1-a*cmax*(1-ratio)
    status=CHEXAL_SINGULAR
    if(den<=0)return
    lower=max(0.0_wp,q/p%rg)
    upper=a*(cmax*max(q/p%rl,0.0_wp)+vmax)/den
    status=CHEXAL_NO_ROOT
    if(upper<=lower)return
    n=256
    do i=0,n
      xx(i)=lower+(upper-lower)*.5_wp*(1-cos(PI*real(i,wp)/256))
    enddo
    do i=1,60
      n=n+1;xx(n)=lower+(upper-lower)*exp(-real(i,wp))
    enddo
    do i=1,n
      t=xx(i);j=i-1
      do while(j>=0)
        if(xx(j)<=t)exit
        xx(j+1)=xx(j);j=j-1
      enddo
      xx(j+1)=t
    enddo
    had_previous=.false.;rr=nan()
    do i=0,n
      x=xx(i)
      if(x<=lower)cycle
      jl=(q-p%rg*x)/p%rl
      fx=residual(a,pv,jl,x,b,st)
      if(st/=CHEXAL_OK)then
        ! Do not bridge a missing / forbidden CCFL interval.
        if(st==CHEXAL_DOMAIN.or.st==CHEXAL_CCFL_EXCEEDED)then
          had_previous=.false.;cycle
        endif
        ! A failed numerical evaluation is NOT evidence of an empty root set.
        status=st
        if(st==CHEXAL_NO_ROOT)status=CHEXAL_NUMERIC
        return
      endif
      if(.not.had_previous.and.fx==0)then
        call chexal_solve_alpha(jl,x,p%rl,p%rg,PI/2,p%d,p%g,p%ml,p%mg,p%sigma,asel,sm,branch=b)
        if(sm/=CHEXAL_OK.and.sm/=CHEXAL_DOMAIN.and.sm/=CHEXAL_NO_ROOT.and.sm/=CHEXAL_CCFL_EXCEEDED)then
          status=sm;return
        endif
        if(sm==CHEXAL_OK)then
          if(abs(asel-a)<=2e-9_wp*max(abs(a),tiny(1.0_wp)))then
            call add_root(x,rr,nroots,sm)
            if(sm/=CHEXAL_OK)then
              status=sm;return
            endif
          endif
        endif
      endif
      if(had_previous)then
        if(fx==0.or.sign(1.0_wp,fx)/=sign(1.0_wp,fprev))then
          lo=prev;hi=x;fl=fprev;fh=fx
          do k=1,120
            if(fl==0)then
              hi=lo;exit
            else if(fh==0)then
              lo=hi;exit
            endif
            mid=lo+(hi-lo)/2
            if(mid==lo.or.mid==hi)exit
            fmid=residual(a,pv,(q-p%rg*mid)/p%rl,mid,b,sm)
            if(sm/=CHEXAL_OK)then
              if(sm==CHEXAL_DOMAIN.or.sm==CHEXAL_CCFL_EXCEEDED)exit
              status=sm
              if(sm==CHEXAL_NO_ROOT)status=CHEXAL_NUMERIC
              return
            endif
            if(fmid==0)then
              lo=mid;hi=mid;exit
            else if(sign(1.0_wp,fmid)==sign(1.0_wp,fl))then
              lo=mid;fl=fmid
            else
              hi=mid;fh=fmid
            endif
          enddo
          mid=lo+(hi-lo)/2;jl=(q-p%rg*mid)/p%rl
          fmid=residual(a,pv,jl,mid,b,sm)
          if(sm/=CHEXAL_OK.and.sm/=CHEXAL_DOMAIN.and.sm/=CHEXAL_CCFL_EXCEEDED)then
            status=sm
            if(sm==CHEXAL_NO_ROOT)status=CHEXAL_NUMERIC
            return
          endif
          if(sm==CHEXAL_OK)then
            if(abs(fmid)<=256*epsilon(1.0_wp)*max(abs(mid),abs(jl),tiny(1.0_wp)))then
              call chexal_solve_alpha(jl,mid,p%rl,p%rg,PI/2,p%d,p%g,p%ml,p%mg,p%sigma,asel,sm,branch=b)
              if(sm/=CHEXAL_OK.and.sm/=CHEXAL_DOMAIN.and.sm/=CHEXAL_NO_ROOT.and.sm/=CHEXAL_CCFL_EXCEEDED)then
                status=sm;return
              endif
              if(sm==CHEXAL_OK)then
                if(abs(asel-a)<=2e-9_wp*max(abs(a),tiny(1.0_wp)))then
                  call add_root(mid,rr,nroots,sm)
                  if(sm/=CHEXAL_OK)then
                    status=sm;return
                  endif
                endif
              endif
            endif
          endif
        endif
      endif
      prev=x;fprev=fx;had_previous=.true.
    enddo
    if(nroots==0)then
      status=CHEXAL_NO_ROOT
    else if(nroots==1)then
      ug=direction*rr(1)/a;ul=direction*(q-p%rg*rr(1))/p%rl/(1-a);status=CHEXAL_OK
    else
      status=CHEXAL_BRANCH_REQUIRED
    endif
  end subroutine

  subroutine recover_countercurrent(a,p,momentum,b,ughint,ul,ug,status)
    !! Local Newton continuation inside the SAME vertical CC domain, followed
    !! by independent alpha-branch verification. It never searches a remote
    !! root after continuation fails or crosses the CCFL boundary.
    real(wp),intent(in)::a,momentum,ughint
    type(parameters),intent(in)::p
    integer,intent(in)::b
    real(wp),intent(out)::ul,ug
    integer,intent(out)::status
    type(parameters)::pv
    real(wp)::direction,q,x,bound,jl,f,fp,fm,h,der,step,trial,ft,anew,scale
    integer::i,k,st,stp,stm,n
    ul=nan();ug=nan();status=CHEXAL_DOMAIN
    if(abs(abs(p%theta)-PI/2)>ANGLE_TOL)return
    status=CHEXAL_BRANCH_REQUIRED
    if(b==CHEXAL_BRANCH_AUTO)return
    direction=1;if(p%theta<0)direction=-1
    pv=p;pv%theta=PI/2;q=direction*momentum
    x=direction*a*ughint;bound=max(0.0_wp,q/p%rg)
    status=CHEXAL_DOMAIN
    if(x<=bound)return
    do i=1,32
      jl=(q-p%rg*x)/p%rl
      f=residual(a,pv,jl,x,b,status)
      if(status/=CHEXAL_OK)return
      scale=max(abs(x),abs(jl),tiny(1.0_wp))
      if(abs(f)<=64*epsilon(1.0_wp)*scale)exit
      h=min(1e-5_wp*max(abs(x),1e-8_wp),0.1_wp*(x-bound))
      do k=1,20
        fp=residual(a,pv,(q-p%rg*(x+h))/p%rl,x+h,b,stp)
        fm=residual(a,pv,(q-p%rg*(x-h))/p%rl,x-h,b,stm)
        if(stp==CHEXAL_OK.and.stm==CHEXAL_OK)exit
        h=h/2
      enddo
      if(k>20.or.x+h==x.or.x-h==x)then
        status=CHEXAL_SINGULAR;return
      endif
      der=(fp-fm)/(2*h)
      if(abs(der)<=64*epsilon(1.0_wp))then
        status=CHEXAL_SINGULAR;return
      endif
      step=-f/der
      do k=1,30
        trial=x+step
        if(trial<=bound)then
          step=step/2;cycle
        endif
        ft=residual(a,pv,(q-p%rg*trial)/p%rl,trial,b,st)
        if(st==CHEXAL_OK)then
          if(abs(ft)<abs(f))exit
        endif
        step=step/2
      enddo
      if(k>30)then
        status=CHEXAL_NUMERIC;return
      endif
      x=trial
    enddo
    if(i>32)then
      status=CHEXAL_NUMERIC;return
    endif
    jl=(q-p%rg*x)/p%rl
    call chexal_solve_alpha(jl,x,p%rl,p%rg,PI/2,p%d,p%g,p%ml,p%mg,p%sigma,anew,status, &
      branch=b,nroots=n)
    if(status/=CHEXAL_OK)return
    if(abs(anew-a)>2e-9_wp*max(abs(a),tiny(1.0_wp)))then
      status=CHEXAL_BRANCH_REQUIRED;return
    endif
    ug=direction*x/a;ul=direction*jl/(1-a);status=CHEXAL_OK
  end subroutine
end module module_chexal_1992
