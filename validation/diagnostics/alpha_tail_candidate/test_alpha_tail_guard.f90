program test_alpha_tail_guard
  use mod_kinds,only:wp
  use, intrinsic::iso_fortran_env,only:int64
  use, intrinsic::ieee_arithmetic,only:ieee_is_finite
  use module_steady_state
  use module_closure_law,only:BHAGWAT,STRATIFIED,CLOSURE_DOMAIN,lei_fechamento_or_ndim_simp
  implicit none
  real(wp),parameter::pi=acos(-1.0_wp),d=.05254_wp,area=pi*d*d/4,eps=4.6e-5_wp,g=9.81_wp
  real(wp),parameter::mul=.00089_wp,mug=.0000185_wp,sigma=.0728_wp,wr=997.05_wp,wu=.53_wp,tol=1e-12_wp
  real(wp)::a,b,jl,jg,rl,rg,theta,expected,p,residual,low,high
  integer::id,i,code,nroots,unit,number,nfixtures,checks
  character(30)::mode
  checks=0
  call get_command_argument(1,mode)
  ! Byte-comparable roots/statuses along the old successful paths. The old and
  ! isolated modules are compiled against this identical caller at O0 and O3.
  do id=1,4
    do i=1,24
      jl=.02_wp+.04_wp*real(i,wp);jg=.005_wp+.03_wp*real(i,wp)
      theta=(real(i-1,wp)/23.0_wp)*(pi/2)
      if(id==STRATIFIED)theta=-2*pi/180
      a=solve_alpha_from_superficials(jl,jg,1.0_wp,.003_wp,theta,d,area,eps,g,mul,mug,sigma, &
        1.0_wp,1000.0_wp,tol,status=code,closure_id=id,nroots=nroots)
      write(*,'(a,4i6,z18.16)')'OLD ',id,i,code,nroots,transfer(a,0_int64)
    enddo
  enddo
  a=solve_alpha_from_superficials(.1_wp,-.01_wp,1.0_wp,.003_wp,-2*pi/180,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=STRATIFIED,nroots=nroots)
  call check(code==STEADY_MULTIPLE_ROOTS.and.nroots==2.and..not.ieee_is_finite(a),'multiple roots remain qualified')
  write(*,'(a,2i6,z18.16)')'MULTI ',code,nroots,transfer(a,0_int64)
  low=solve_alpha_from_superficials(.1_wp,-.01_wp,1.0_wp,.003_wp,-2*pi/180,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=STRATIFIED,nroots=nroots,alpha_hint=0.0_wp)
  call check(code==0.and.nroots==2,'lower multiple root hint retained')
  write(*,'(a,2i6,z18.16)')'LOWER ',code,nroots,transfer(low,0_int64)
  high=solve_alpha_from_superficials(.1_wp,-.01_wp,1.0_wp,.003_wp,-2*pi/180,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=STRATIFIED,nroots=nroots,alpha_hint=1.0_wp)
  call check(code==0.and.nroots==2.and.high>low,'upper multiple root hint retained')
  write(*,'(a,2i6,z18.16)')'UPPER ',code,nroots,transfer(high,0_int64)
  a=solve_alpha_from_superficials(.5_wp,0.0_wp,1.0_wp,.003_wp,0.0_wp,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=BHAGWAT,nroots=nroots,phase_hint=PHASE_LIQUID_ONLY)
  call check(code==0.and.a==0.0_wp.and.nroots==1,'liquid-only inventory limit retained')
  write(*,'(a,2i6,z18.16)')'LIQUID ',code,nroots,transfer(a,0_int64)
  a=solve_alpha_from_superficials(0.0_wp,.1_wp,1.0_wp,.003_wp,0.0_wp,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=BHAGWAT,nroots=nroots,phase_hint=PHASE_GAS_ONLY)
  call check(code==0.and.a==1.0_wp.and.nroots==1,'gas-only inventory limit retained')
  write(*,'(a,2i6,z18.16)')'GAS ',code,nroots,transfer(a,0_int64)
  if(trim(mode)=='emit')stop
  open(newunit=unit,file='tail_fixtures.dat',status='old',action='read')
  read(unit,*)nfixtures
  do i=1,nfixtures
    read(unit,*)number,p,jl,jg,theta,expected
    rl=(997.05_wp+(p-101325.0_wp)/1498.0_wp**2)/wr;rg=p/343.0_wp**2/wr
    a=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
      status=code,closure_id=BHAGWAT,nroots=nroots)
    call check(code==0.and.nroots==1.and.ieee_is_finite(a),'tiny positive flux returns one qualified root')
    call check(a>0.0_wp.and.a<1.0_wp,'root is positive without a phase floor')
    call check(abs(a/expected-1.0_wp)<5e-13_wp,'root agrees with independent 80-digit fixture')
    b=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
      status=code,closure_id=BHAGWAT,alpha_hint=expected)
    call check(code==0.and.transfer(a,0_int64)==transfer(b,0_int64),'nroots and hint do not change a unique tail root')
    call lei_fechamento_or_ndim_simp(a,rl,rg,jl/(1-a)/wu,jg/a/wu,theta,d,area,eps,g, &
      mul,mug,sigma,wu,wr,tol,residual,BHAGWAT,code)
    call check(code==0.and.abs(residual/(jg/wu))<5e-13_wp,'relative closure residual is small at tiny root')
    b=solve_alpha_from_superficials(jl/wu,-jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
      status=code,closure_id=BHAGWAT,nroots=nroots,alpha_hint=expected)
    call check(code==CLOSURE_DOMAIN.and.nroots==0.and..not.ieee_is_finite(b),'negative gas flux still rejected')
    write(*,'(a,i6,3es25.16)')'TAIL ',number,jg,a,a/expected-1
  enddo
  close(unit)
  ! Below the normal range, an ULP bound replaces an impossible fixed relative
  ! tolerance. The root must still come from a sign bracket, never from a floor.
  jl=.5295514577421165_wp;jg=1e-320_wp;theta=0.0_wp
  a=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
    status=code,closure_id=BHAGWAT,nroots=nroots)
  call check(code==0.and.nroots==1.and.a>0.0_wp,'subnormal positive root remains representable')
  expected=jg*1.1206300892060053278_wp
  call check(abs(a-expected)<=16*nearest(0.0_wp,1.0_wp),'subnormal root agrees within floating-point granularity')
  call lei_fechamento_or_ndim_simp(a,rl,rg,jl/(1-a)/wu,jg/a/wu,theta,d,area,eps,g, &
    mul,mug,sigma,wu,wr,tol,residual,BHAGWAT,code)
  call check(code==0.and.abs(residual)<=16*nearest(0.0_wp,1.0_wp),'subnormal closure residual bounded by ULP')
  b=solve_alpha_from_superficials(jl/wu,-jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
    status=code,closure_id=BHAGWAT,nroots=nroots)
  call check(code==CLOSURE_DOMAIN.and.nroots==0,'negative subnormal gas flux is rejected')
  jg=nearest(0.0_wp,1.0_wp)
  a=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
    status=code,closure_id=BHAGWAT,nroots=nroots)
  call check(code==0.and.nroots==1.and.a>=jg.and.a<=2*jg,'lowest representable root is bracketed')
  a=solve_alpha_from_superficials(2.0_wp/wu,jg/wu,rl,rg,theta,d,area,eps,g,mul,mug,sigma,wu,wr,tol, &
    status=code,closure_id=BHAGWAT,nroots=nroots)
  call check(code==STEADY_NO_ROOT.and.nroots==0.and..not.ieee_is_finite(a),'no clipping when a positive root lies below representability')
  a=solve_alpha_from_superficials(.5_wp,0.0_wp,1.0_wp,.003_wp,0.0_wp,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=BHAGWAT,nroots=nroots)
  call check(code==STEADY_NO_ROOT.and.nroots==0.and..not.ieee_is_finite(a),'zero flux is not a fabricated pure phase')
  a=solve_alpha_from_superficials(.5_wp,1e-16_wp,1.0_wp,.003_wp,-.1_wp,d,area,eps,g,mul,mug,sigma, &
    1.0_wp,1000.0_wp,tol,status=code,closure_id=BHAGWAT,nroots=nroots)
  call check(code==STEADY_NO_ROOT.and.nroots==0.and..not.ieee_is_finite(a),'fallback does not extend downhill correlation')
  write(*,'(a,i0,a)')'PASS ',checks,' alpha tail guard checks'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    checks=checks+1
    if(.not.ok)then
      write(*,'(a)')'FAIL: '//label
      error stop 1
    endif
  end subroutine
end program
