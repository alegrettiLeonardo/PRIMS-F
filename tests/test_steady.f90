program test_steady
  use mod_kinds,only:wp
  use module_steady_state
  use module_closure_law,only:BENDIKSEN,STRATIFIED,BHAGWAT,lei_fechamento_or_ndim_simp
  use mod_void_helpers,only:ffan
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  real(wp),parameter::pi=acos(-1.0_wp),d=.05254_wp,area=pi*d*d/4.0_wp
  real(wp),parameter::eps=4.6e-5_wp,g=9.81_wp,mul=.89e-3_wp,mug=1.85e-5_wp,sigma=.0728_wp
  real(wp),parameter::lp=57.4_wp,ca=3.018181937256615_wp,x=6.435_wp
  real(wp),parameter::lr=ca*sinh(x/ca),ztop=ca*(cosh(x/ca)-1.0_wp),beta=2*pi/180
  real(wp),parameter::wr=1000.0_wp,wu=1.0_wp,wpres=1e5_wp
  real(wp),parameter::cl=149.8_wp,cg=34.3_wp,rho0=.99705_wp,p0=1.01325_wp,tol=1e-9_wp
  real(wp),allocatable::s(:),p(:),rg(:),rl(:),a(:),ug(:),ul(:),theta(:),jl(:),jg(:),faces(:),pf(:)
  real(wp),allocatable::centers(:),nodes(:)
  real(wp)::alpha,gl,gg,expected,z,rhotop,err,plow,pjoin,al,ar,s0,pp,dps,h,hplus,hminus,source
  real(wp)::rhom,jmix,mu,re,ff,totalderiv,fast_alpha,scan_alpha,dfast,dscan,press
  real(wp)::low_root,high_root,ap,am,ap2,am2,residual_value
  real(wp),parameter::case_gl(4)=[.5922477_wp,.54538635_wp,.279174_wp,.5284365_wp]
  real(wp),parameter::case_gg(4)=[.8307813542792721_wp,1.396675899947762_wp, &
    .577934855150798_wp,.13244340430539123_wp]/wr
  integer::st,nroots,i,k,count
  count=0
  call malha_progressiva(5.0_wp,4.0_wp,4,3,1.2_wp,1.2_wp,centers,nodes,st)
  call check(st==0.and.size(centers)==7.and.size(nodes)==8,'mesh counts and success')
  call close(nodes(5),4.0_wp,1e-14_wp,'junction is exact face')
  call close(nodes(8),5.0_wp,1e-14_wp,'top is exact face')
  call close((nodes(3)-nodes(2))/(nodes(2)-nodes(1)),1.0_wp/1.2_wp,1e-13_wp,'pipe mesh ratio honored')
  call close((nodes(7)-nodes(6))/(nodes(6)-nodes(5)),1.2_wp,1e-13_wp,'riser mesh ratio honored')
  alpha=solve_alpha_from_superficials(.5_wp,0.0_wp,1.0_wp,.003_wp,.4_wp,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=BENDIKSEN,phase_hint=PHASE_LIQUID_ONLY)
  call check(st==0.and.alpha==0.0_wp,'declared liquid-only inventory has exact alpha')
  alpha=solve_alpha_from_superficials(0.0_wp,.5_wp,1.0_wp,.003_wp,.4_wp,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=BENDIKSEN,phase_hint=PHASE_GAS_ONLY)
  call check(st==0.and.alpha==1.0_wp,'declared gas-only inventory has exact alpha')
  alpha=solve_alpha_from_superficials(-10.0_wp,.01_wp,1.0_wp,.003_wp,pi/2,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=BENDIKSEN)
  call check(st/=0.and..not.ieee_is_finite(alpha),'no physical root returns failure, not clipping')
  alpha=solve_alpha_from_superficials(.5_wp,.1_wp,1.0_wp,.003_wp,.4_wp,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=99)
  call check(st/=0.and..not.ieee_is_finite(alpha),'unknown closure remains explicit')
  call test_bhagwat_roots()
  call test_bhagwat_base_trace()

  ! Zero throughput is compatible with a nonzero inventory. In a descending
  ! stratified pipe, stagnant gas must remain a physical interior closure root.
  alpha=solve_alpha_from_superficials(.53_wp,0.0_wp,.99714_wp,.00254_wp,-beta,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=STRATIFIED,nroots=nroots,alpha_hint=.56_wp)
  call check(st==0.and.alpha>.1_wp.and.alpha<.9_wp,'zero gas flow retains stagnant gas inventory')
  call lei_fechamento_or_ndim_simp(alpha,.99714_wp,.00254_wp,.53_wp/(1-alpha),0.0_wp, &
    -beta,d,area,eps,g,mul,mug,sigma,wu,wr,tol,residual_value,STRATIFIED,st)
  call check(st==0.and.abs(residual_value)<1e-9_wp,'stagnant gas root satisfies stratified closure')
  low_root=alpha;h=1e-7_wp
  ap=solve_alpha_from_superficials(.53_wp,h,.99714_wp,.00254_wp,-beta,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=STRATIFIED,alpha_hint=low_root)
  am=solve_alpha_from_superficials(.53_wp,-h,.99714_wp,.00254_wp,-beta,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=STRATIFIED,alpha_hint=low_root)
  call check(st==0.and.abs(ap-low_root)<1e-5_wp.and.abs(am-low_root)<1e-5_wp, &
    'inventory branch continuous through zero gas flux')
  alpha=solve_alpha_from_superficials(0.0_wp,.5_wp,1.0_wp,.003_wp,.4_wp,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=BENDIKSEN)
  call check(st==0.and.alpha>0.0_wp.and.alpha<1.0_wp,'zero liquid flow is not assumed liquid absence')
  ! Production endpoint-Brent agrees with the full branch audit across the
  ! pressure range used by the catenary cases; FD derivatives must also agree.
  do k=1,4
   gl=case_gl(k);gg=case_gg(k)
   do i=1,8
    press=1.8_wp+.45_wp*real(i,wp)
    call alpha_at_pressure(press,.false.,fast_alpha)
    call alpha_at_pressure(press,.true.,scan_alpha)
    call close(fast_alpha,scan_alpha,1e-11_wp,'fast root agrees with full logit scan')
    h=2e-5_wp
    call alpha_at_pressure(press+h,.false.,al)
    call alpha_at_pressure(press-h,.false.,ar)
    dfast=(al-ar)/(2*h)
    call alpha_at_pressure(press+h,.true.,al)
    call alpha_at_pressure(press-h,.true.,ar)
    dscan=(al-ar)/(2*h)
    call close(dfast,dscan,5e-8_wp,'fast and full-scan pressure derivatives agree')
   enddo
  enddo
  ! Countercurrent stratified roots: the default refuses ambiguity; an explicit
  ! trace hint selects a valid continuous branch without modifying the roots.
  alpha=multi_root(.1_wp,-1.0_wp)
  call check(st==STEADY_MULTIPLE_ROOTS.and.nroots==2.and..not.ieee_is_finite(alpha), &
    'ambiguous countercurrent roots remain an error without a hint')
  low_root=multi_root(.1_wp,0.0_wp)
  call check(st==0.and.nroots==2,'lower branch selected by limiting hint')
  high_root=multi_root(.1_wp,1.0_wp)
  call check(st==0.and.nroots==2.and.high_root>low_root+.1_wp,'upper branch remains distinct and admissible')
  call lei_fechamento_or_ndim_simp(low_root,1.0_wp,.003_wp,.1_wp/(1-low_root),-.01_wp/low_root, &
    -beta,d,area,eps,g,mul,mug,sigma,wu,wr,tol,residual_value,STRATIFIED,st)
  call check(st==0.and.abs(residual_value)<1e-9_wp,'selected lower branch satisfies physical closure')
  call lei_fechamento_or_ndim_simp(high_root,1.0_wp,.003_wp,.1_wp/(1-high_root),-.01_wp/high_root, &
    -beta,d,area,eps,g,mul,mug,sigma,wu,wr,tol,residual_value,STRATIFIED,st)
  call check(st==0.and.abs(residual_value)<1e-9_wp,'selected upper branch satisfies physical closure')
  h=2e-5_wp
  ap=multi_root(.1_wp+h,low_root);am=multi_root(.1_wp-h,low_root)
  ap2=multi_root(.1_wp+h/2,low_root);am2=multi_root(.1_wp-h/2,low_root)
  call close((ap-am)/(2*h),(ap2-am2)/h,2e-6_wp,'lower branch finite difference converges with continuation')
  ap=multi_root(.1_wp+h,high_root);am=multi_root(.1_wp-h,high_root)
  ap2=multi_root(.1_wp+h/2,high_root);am2=multi_root(.1_wp-h/2,high_root)
  call close((ap-am)/(2*h),(ap2-am2)/h,2e-6_wp,'upper branch finite difference converges with continuation')
  ! Independent exact hydrostatic solution for the compressible-liquid EOS.
  gl=0.0_wp;gg=0.0_wp
  call solve(12,10,1.04_wp,4.2_wp)
  call check(st==0,'hydrostatic initialization succeeds')
  call check(size(s)==22,'independent requested riser count used')
  call close(pf(size(pf)),4.2_wp,1e-14_wp,'separator pressure applied at top face')
  call check(abs(p(size(p))-4.2_wp)>1e-4_wp,'last center does not inherit face pressure')
  rhotop=rho0+(4.2_wp-p0)/cl**2
  err=0.0_wp
  do i=1,size(s)
    z=elevation(s(i))
    expected=p0+cl**2*(rhotop*exp(g*(ztop-z)/1498.0_wp**2)-rho0)
    err=max(err,abs(p(i)-expected))
  enddo
  call check(err<2e-8_wp,'compressible hydrostatics at all centers agrees with exact solution')
  call check(all(a==0.0_wp).and.all(ul==0.0_wp).and.all(ug==0.0_wp),'static pure-liquid state has exact limits')

  ! Two-phase case: independent phase mass invariants and junction traces.
  gl=.5284365_wp;gg=.13244340430539123_wp/wr
  call solve(8,8,1.02_wp,2.08_wp)
  call check(st==0,'two-phase catenary/pipeline initialization succeeds')
  call check(all(ieee_is_finite(p)).and.all(p>0.0_wp),'finite positive initial pressures')
  call check(maxval(abs((1-a)*rl*ul-gl))<1e-12_wp,'liquid mass flux preserved throughout')
  call check(maxval(abs(a*rg*ug-gg))<1e-12_wp,'gas mass flux preserved throughout')
  plow=pf(1);pjoin=pf(9)
  expected=pjoin/cg**2
  rhotop=rho0+(pjoin-p0)/cl**2
  al=solve_alpha_from_superficials(gl/rhotop,gg/expected,rhotop,expected,-beta,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=STRATIFIED,nroots=nroots)
  call check(st==0.and.nroots==1,'pipeline junction trace has a certified bracketed root')
  ar=solve_alpha_from_superficials(gl/rhotop,gg/expected,rhotop,expected,0.0_wp,d,area,eps,g,mul,mug, &
    sigma,wu,wr,tol,status=st,closure_id=BENDIKSEN)
  call check(st==0.and.abs(al-ar)>.01_wp,'same junction pressure and fluxes allow distinct alpha traces')
  call solve(16,16,1.0_wp,2.08_wp)
  call check(st==0,'refined initialization succeeds')
  call close(pf(1),plow,2e-7_wp,'initial pressure independent of output mesh')
  call close(pf(17),pjoin,2e-7_wp,'junction pressure independent of output mesh')

  ! Independent total derivative of H along the computed pressure tangent.
  ! This fails if the explicit catenary-angle contribution Htheta*theta' is absent.
  gl=.54538635_wp;gg=1.396675899947762_wp/wr
  s0=lp/lr+.35_wp;pp=3.0_wp
  dps=fun_dpds(s0,pp,gl,gg,lp,lr,ca,beta,d,area,eps,g,cl,cg,rho0,p0, &
    mul,mug,sigma,wpres,wu,wr,tol,status=st,closure_id_riser=BENDIKSEN)
  call check(st==0,'steady RHS evaluated with varying riser angle')
  h=2e-5_wp
  hplus=momentum(s0+h,pp+dps*h)
  hminus=momentum(s0-h,pp-dps*h)
  totalderiv=(hplus-hminus)/(2*h)
  expected=pp/cg**2;rhotop=rho0+(pp-p0)/cl**2
  alpha=solve_alpha_from_superficials(gl/rhotop,gg/expected,rhotop,expected,atan((s0*lr-lp)/ca), &
    d,area,eps,g,mul,mug,sigma,wu,wr,tol,status=st,closure_id=BENDIKSEN)
  rhom=(1-alpha)*rhotop+alpha*expected;jmix=gl/rhotop+gg/expected
  mu=(1-alpha)*mul+alpha*mug;re=wr*wu*d*rhom*abs(jmix)/mu
  ff=ffan(eps/d,re,st)
  source=-g*lr/wu**2*rhom*sin(atan((s0*lr-lp)/ca))-2*lr/d*ff*rhom*jmix*abs(jmix)
  call close(totalderiv,source,2e-7_wp,'full derivative of H equals gravity plus Fanning friction')
  write(*,'(a,i0,a)') 'PASS ',count,' steady checks'
contains
  subroutine test_bhagwat_base_trace()
    ! The catenary has theta>=0. At its lower face the angular derivative is
    ! the right trace, even when Bhagwat C4 changes sign for theta<0.
    real(wp)::sb,pr,rli,rgi,aa,jli,jgi,frsg,dps0,dpsnear,hs,h0,h1,h2,deriv,ssource
    real(wp)::rm0,j0,mu0,re0,ff0,base_pressure,inlet_pressure
    integer::code,nr,ii
    gl=case_gl(4);gg=case_gg(4);pr=3.0_wp;sb=lp/lr
    rli=rho0+(pr-p0)/cl**2;rgi=pr/cg**2;jli=gl/rli;jgi=gg/rgi
    aa=solve_alpha_from_superficials(jli,jgi,rli,rgi,0.0_wp,d,area,eps,g,mul,mug,sigma,wu,wr, &
      1e-12_wp,status=code,closure_id=BHAGWAT,nroots=nr)
    call check(code==0.and.nr==1.and.aa<.15_wp,'small-alpha Bhagwat root at the horizontal riser face')
    frsg=sqrt(rgi/(rli-rgi))*jgi*wu/sqrt(g*d)
    call check(frsg<.01_wp,'base regression lies in the low-Frsg C4 angular switch range')
    dps0=fun_dpds(sb,pr,gl,gg,lp,lr,ca,beta,d,area,eps,g,cl,cg,rho0,p0, &
      mul,mug,sigma,wpres,wu,wr,1e-11_wp,status=code,closure_id_riser=BHAGWAT,section_id=2)
    call check(code==0,'right-trace steady RHS exists at the lower riser face')
    rm0=(1-aa)*rli+aa*rgi;j0=jli+jgi;mu0=(1-aa)*mul+aa*mug
    re0=wr*wu*d*rm0*abs(j0)/mu0;ff0=ffan(eps/d,re0,code)
    call check(code==0,'base regression friction coefficient exists')
    ssource=-2*lr/d*ff0*rm0*j0*abs(j0)
    h0=momentum(sb,pr,BHAGWAT)
    hs=2e-5_wp
    do ii=1,2
      h1=momentum(sb+hs,pr+dps0*hs,BHAGWAT)
      h2=momentum(sb+2*hs,pr+2*dps0*hs,BHAGWAT)
      deriv=(4*(h1-h0)-(h2-h0))/(2*hs)
      call close(deriv,ssource,3e-7_wp,'one-sided dH/ds equals S at the base of the Bhagwat riser')
      hs=hs/2
    enddo
    dpsnear=fun_dpds(sb+1e-8_wp,pr+dps0*1e-8_wp,gl,gg,lp,lr,ca,beta,d,area,eps,g,cl,cg,rho0,p0, &
      mul,mug,sigma,wpres,wu,wr,1e-11_wp,status=code,closure_id_riser=BHAGWAT,section_id=2)
    call check(code==0,'steady RHS exists just above the riser face')
    call close(dpsnear,dps0,1e-6_wp,'base RHS agrees with its local right limit')
    ! The explicit trace tolerates geometric roundoff only, without allowing a
    ! physical continuation of the riser below the junction.
    dpsnear=fun_dpds(sb-spacing(sb),pr,gl,gg,lp,lr,ca,beta,d,area,eps,g,cl,cg,rho0,p0, &
      mul,mug,sigma,wpres,wu,wr,1e-11_wp,status=code,closure_id_riser=BHAGWAT,section_id=2)
    call check(code==0,'roundoff at the shared face does not select downhill Bhagwat C4')
    call close(dpsnear,dps0,1e-10_wp,'geometric roundoff preserves the right face trace')
    call solve(8,8,1.04_wp,2.08_wp,BHAGWAT)
    call check(st==0,'coarse Bhagwat steady profile reaches the exact lower face')
    base_pressure=pf(9);inlet_pressure=pf(1)
    call solve(16,16,1.0_wp,2.08_wp,BHAGWAT)
    call check(st==0,'refined Bhagwat steady profile reaches the exact lower face')
    call close(pf(17),base_pressure,2e-7_wp,'Bhagwat base-face pressure converges across output meshes')
    call close(pf(1),inlet_pressure,2e-7_wp,'Bhagwat inlet pressure converges across output meshes')
  end subroutine
  subroutine test_bhagwat_roots()
    ! The alpha-dependent correlation must satisfy its full nonlinear law at
    ! fixed superficial fluxes. Neither coefficients frozen at alpha=.5 nor
    ! the explicit Bendiksen formula is an admissible substitute.
    real(wp),parameter::ps_case(4)=[4.2_wp,2.06_wp,2.01_wp,2.08_wp]
    real(wp)::rli,rgi,jsl,jsg,aa,bb,cc,angle,pr,rr,delta,da1,da2,a0,scaled
    integer::case_i,angle_i,code,nr,ii
    do case_i=1,4
      pr=ps_case(case_i);rli=rho0+(pr-p0)/cl**2;rgi=pr/cg**2
      jsl=case_gl(case_i)/rli;jsg=case_gg(case_i)/rgi
      do angle_i=0,4
        angle=real(angle_i,wp)*pi/6
        if(angle_i==4)angle=1.319613718692701_wp
        aa=solve_alpha_from_superficials(jsl,jsg,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
          1e-11_wp,status=code,closure_id=BHAGWAT,nroots=nr)
        call check(code==0.and.nr==1.and.aa>0.0_wp.and.aa<1.0_wp,'Bhagwat full scan has one interior root')
        call lei_fechamento_or_ndim_simp(aa,rli,rgi,jsl/(1-aa),jsg/aa,angle,d,area,eps,g,mul,mug, &
          sigma,wu,wr,tol,rr,BHAGWAT,code)
        call check(code==0.and.abs(rr)<1e-10_wp,'Bhagwat alpha satisfies complete nonlinear closure')
      enddo
      gl=case_gl(case_i);gg=case_gg(case_i)
      call solve(4,4,1.0_wp,pr,BHAGWAT)
      call check(st==0.and.all(ieee_is_finite(p)),'D profile initialized with Bhagwat riser')
      call close(pf(size(pf)),pr,1e-13_wp,'Bhagwat separator pressure remains a face condition')
      call check(maxval(abs((1-a)*rl*ul-gl))<1e-12_wp,'Bhagwat steady liquid mass flux preserved')
      call check(maxval(abs(a*rg*ug-gg))<1e-12_wp,'Bhagwat steady gas mass flux preserved')
      do ii=5,8
        call lei_fechamento_or_ndim_simp(a(ii),rl(ii),rg(ii),ul(ii),ug(ii),theta(ii),d,area,eps,g,mul,mug, &
          sigma,wu,wr,tol,rr,BHAGWAT,code)
        call check(code==0.and.abs(rr)<1e-10_wp,'integrated riser retains Bhagwat closure at each output')
      enddo
    enddo
    ! The former ID2 Fr=3.5 barrier is not a coefficient switch in ID3. Check
    ! the alpha root and its derivative through that same total superficial flux.
    rli=.9970971476888178_wp;rgi=.001760530070787292_wp
    angle=1.319613718692701_wp;jsl=1.28491342508_wp;jsg=3.5_wp*sqrt(g*d)-jsl
    a0=solve_alpha_from_superficials(jsl,jsg,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
      1e-12_wp,status=code,closure_id=BHAGWAT,nroots=nr)
    call check(code==0.and.nr==1,'Bhagwat root at total Fr=3.5')
    delta=1e-5_wp
    do ii=1,2
      aa=solve_alpha_from_superficials(jsl,jsg+delta,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
        1e-12_wp,status=code,closure_id=BHAGWAT,alpha_hint=a0)
      call check(code==0,'Bhagwat root above former Froude threshold')
      bb=solve_alpha_from_superficials(jsl,jsg-delta,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
        1e-12_wp,status=code,closure_id=BHAGWAT,alpha_hint=a0)
      call check(code==0,'Bhagwat root below former Froude threshold')
      cc=(aa-bb)/(2*delta)
      if(ii==1)then
        da1=cc
      else
        da2=cc
      endif
      delta=delta/2
    enddo
    call close(da1,da2,1e-7_wp,'Bhagwat root derivative is stable through Fr=3.5')
    scaled=solve_alpha_from_superficials(jsl/2,jsg/2,rli*wr/333.0_wp,rgi*wr/333.0_wp,angle, &
      d,area,eps,g,mul,mug,sigma,2.0_wp,333.0_wp,1e-12_wp,status=code,closure_id=BHAGWAT)
    call check(code==0,'Bhagwat root with alternate nondimensional scales')
    call close(scaled,a0,1e-11_wp,'Bhagwat root invariant to nondimensional scales')
    aa=solve_alpha_from_superficials(.5_wp,-.1_wp,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
      tol,status=code,closure_id=BHAGWAT)
    call check(code/=0.and..not.ieee_is_finite(aa),'Bhagwat rejects unqualified countercurrent input')
    aa=solve_alpha_from_superficials(.5_wp,0.0_wp,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
      tol,status=code,closure_id=BHAGWAT)
    call check(code==STEADY_NO_ROOT.and..not.ieee_is_finite(aa),'zero gas flux alone does not declare Bhagwat phase absence')
    aa=solve_alpha_from_superficials(0.0_wp,0.0_wp,rli,rgi,angle,d,area,eps,0.0_wp,mul,mug,sigma,wu,wr, &
      tol,status=code,closure_id=BHAGWAT,phase_hint=PHASE_LIQUID_ONLY)
    call check(code==0.and.aa==0.0_wp,'declared single liquid phase bypasses two-phase Bhagwat domain')
    aa=solve_alpha_from_superficials(0.0_wp,.5_wp,rli,rgi,angle,d,area,eps,g,mul,mug,sigma,wu,wr, &
      tol,status=code,closure_id=BHAGWAT,phase_hint=PHASE_GAS_ONLY)
    call check(code==0.and.aa==1.0_wp,'declared single gas phase bypasses two-phase Bhagwat domain')
  end subroutine
  real(wp) function multi_root(liquid,hint) result(aa)
    real(wp),intent(in)::liquid,hint
    if(hint<0.0_wp)then
      aa=solve_alpha_from_superficials(liquid,-.01_wp,1.0_wp,.003_wp,-beta,d,area,eps,g,mul,mug, &
        sigma,wu,wr,1e-11_wp,status=st,closure_id=STRATIFIED,nroots=nroots)
    else
      aa=solve_alpha_from_superficials(liquid,-.01_wp,1.0_wp,.003_wp,-beta,d,area,eps,g,mul,mug, &
        sigma,wu,wr,1e-11_wp,status=st,closure_id=STRATIFIED,nroots=nroots,alpha_hint=hint)
    endif
  end function
  subroutine alpha_at_pressure(pressure,full_scan,aa)
    real(wp),intent(in)::pressure
    logical,intent(in)::full_scan
    real(wp),intent(out)::aa
    real(wp)::rli,rgi
    integer::code,nr
    rli=rho0+(pressure-p0)/cl**2;rgi=pressure/cg**2
    if(full_scan)then
      aa=solve_alpha_from_superficials(gl/rli,gg/rgi,rli,rgi, &
        -beta,d,area,eps,g,mul,mug,sigma,wu,wr,1e-11_wp,status=code,closure_id=STRATIFIED,nroots=nr)
      call check(code==0.and.nr==1,'full branch audit unique')
    else
      aa=solve_alpha_from_superficials(gl/rli,gg/rgi,rli,rgi, &
        -beta,d,area,eps,g,mul,mug,sigma,wu,wr,1e-11_wp,status=code,closure_id=STRATIFIED)
      call check(code==0,'fast branch converged')
    endif
  end subroutine
  subroutine solve(np,nr,ratio,ps,closure)
    integer,intent(in)::np,nr
    real(wp),intent(in)::ratio,ps
    integer,intent(in),optional::closure
    integer::model
    model=BENDIKSEN
    if(present(closure))model=closure
    call steady_state_ndim_simp(np,gl,gg,ps,lp,lr,ca,beta,d,area,eps,g,cl,cg,rho0,p0,mul,mug,sigma, &
      wpres,wu,wr,tol,s,p,rg,rl,a,ug,ul,theta,jl,jg,mesh_ratio=ratio,n_riser=nr,face_positions=faces, &
      closure_id_riser=model,status=st,face_pressures=pf)
  end subroutine
  real(wp) function elevation(ss) result(zz)
    real(wp),intent(in)::ss
    real(wp)::local
    if(ss*lr<lp)then
      zz=(lp-ss*lr)*sin(beta)
    else
      local=ss*lr-lp
      zz=sqrt(ca*ca+local*local)-ca
    endif
  end function
  real(wp) function momentum(ss,pressure,closure) result(hh)
    real(wp),intent(in)::ss,pressure
    integer,intent(in),optional::closure
    real(wp)::rli,rgi,aa,uli,ugi,th
    integer::code,model
    model=BENDIKSEN
    if(present(closure))model=closure
    rli=rho0+(pressure-p0)/cl**2;rgi=pressure/cg**2
    th=atan(lr*(ss-lp/lr)/ca)
    aa=solve_alpha_from_superficials(gl/rli,gg/rgi,rli,rgi,th,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
      status=code,closure_id=model)
    call check(code==0,'independent momentum evaluation admissible')
    uli=gl/(rli*(1-aa));ugi=gg/(rgi*aa)
    hh=gl*uli+gg*ugi+wpres/(wr*wu**2)*pressure
  end function
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    if(.not.ok)then
      write(*,'(a)')'FAIL '//label
      error stop 1
    endif
    count=count+1
  end subroutine
  subroutine close(actual,target,tolerance,label)
    real(wp),intent(in)::actual,target,tolerance
    character(*),intent(in)::label
    if(.not.ieee_is_finite(actual).or.abs(actual-target)>tolerance*max(1.0_wp,abs(target)))then
      write(*,'(a,3es22.12)')trim(label)//': ',actual,target,actual-target
      error stop 1
    endif
    count=count+1
  end subroutine
end program
