program test_chexal_1992_core
  use mod_kinds,only:wp
  use module_chexal_1992
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  implicit none
  real(wp),parameter::pi=acos(-1.0_wp),rl=998.0_wp,rg=1.2_wp,d=.05254_wp,g=9.81_wp
  real(wp),parameter::ml=.001_wp,mg=.000018_wp,sigma=.072_wp
  real(wp)::c,v,r,a,af,fs,roots(8),u,w,momentum,refa,refh,refl,x,y,h,fp,fm
  integer::st,n,checks,b,i,s1,s2
  checks=0
  call chexal_coefficients_si(.2_wp,rl,rg,.4_wp,.7_wp,pi/2,d,g,ml,mg,sigma,c,v,st)
  call ok(st==0,'vertical up status')
  call near(c,.83152878676784701845_wp,2e-14_wp,'vertical C0 independent 60d reference')
  call near(v,.54480865106966165970_wp,2e-14_wp,'vertical Vgj independent reference')
  call chexal_coefficients_si(.2_wp,rl,rg,.4_wp,.7_wp,0.0_wp,d,g,ml,mg,sigma,c,v,st)
  call near(c,1.0163058837066762876_wp,2e-14_wp,'horizontal C0 reference')
  call near(v,.54480865106966165970_wp,2e-14_wp,'horizontal Vgj remains positive')
  call chexal_coefficients_si(.6_wp,rl,rg,-.5_wp,-.2_wp,pi/2,d,g,ml,mg,sigma,c,v,st)
  call ok(st==0,'signed vertical downflow')
  call near(c,1.2510715870715231473_wp,2e-14_wp,'down C0 signed Reynolds reference')
  call near(v,.48100638255655177481_wp,2e-14_wp,'down Vgj C10 reference')
  call chexal_coefficients_si(.6_wp,rl,rg,.5_wp,.2_wp,-pi/2,d,g,ml,mg,sigma,c,v,st)
  call near(c,1.2510715870715231473_wp,2e-14_wp,'coordinate reversal C0')
  call near(v,-.48100638255655177481_wp,2e-14_wp,'coordinate reversal drift sign')
  call chexal_coefficients_si(.4_wp,rl,rg,-.1_wp,.1_wp,pi/2,d,g,ml,mg,sigma,c,v,st,branch=1)
  call near(c,1.1560561920771662692_wp,2e-14_wp,'CC C0 independent reference')
  call near(v,.40880128591713410621_wp,2e-14_wp,'CC Vgj independent reference')
  call chexal_residual_si(.4_wp,rl,rg,-.1_wp,.1_wp,pi/2,d,g,ml,mg,sigma,r,st,branch=1)
  call ok(st==0.and.ieee_is_finite(r),'j=0 residual finite')
  call near(r,.4_wp*v-.1_wp,2e-14_wp,'j=0 residual formula')
  call chexal_coefficients_si(.4_wp,rl,rg,-.1_wp,.1_wp,pi/2,d,g,ml,mg,sigma,c,v,st)
  call ok(st==CHEXAL_BRANCH_REQUIRED,'CC auto refuses branch')
  call chexal_coefficients_si(.4_wp,rl,rg,-.1_wp,.1_wp,76.473695_wp*pi/180,d,g,ml,mg,sigma,c,v,st,branch=1)
  call ok(st==CHEXAL_DOMAIN,'catenary is not vertical CC')
  call chexal_coefficients_si(.4_wp,rl,rg,-.1_wp,.1_wp,0.0_wp,d,g,ml,mg,sigma,c,v,st,branch=1)
  call ok(st==CHEXAL_DOMAIN,'horizontal CC rejected')
  call chexal_coefficients_si(.4_wp,rl,rg,-.1_wp,-.1_wp,0.0_wp,d,g,ml,mg,sigma,c,v,st)
  call ok(st==CHEXAL_DOMAIN,'horizontal negative flow rejected')
  call chexal_coefficients_si(.4_wp,rl,rg,.1_wp,-.1_wp,pi/2,d,g,ml,mg,sigma,c,v,st,branch=1)
  call ok(st==CHEXAL_DOMAIN,'gas-down liquid-up rejected')
  call chexal_coefficients_si(0.0_wp,rl,rg,.4_wp,0.0_wp,pi/2,d,g,ml,mg,sigma,c,v,st)
  call near(c,0.0_wp,1e-14_wp,'alpha0 C0 limit')
  call near(v,.62098758247803309918_wp,2e-14_wp,'alpha0 finite drift')
  call chexal_coefficients_si(1.0_wp,rl,rg,0.0_wp,.7_wp,pi/2,d,g,ml,mg,sigma,c,v,st)
  call near(c,1.0_wp,1e-14_wp,'alpha1 C0 limit')
  call near(v,0.0_wp,1e-14_wp,'alpha1 drift limit')
  call chexal_solve_alpha(.4_wp,.7_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,nroots=n)
  call ok(st==0.and.n==1,'co-current unique sampled root')
  call near(a,.42010691165400489012_wp,2e-14_wp,'co-current alpha reference')
  momentum=rl*.4_wp+rg*.7_wp
  call chexal_recover_velocities(a,rl,rg,momentum,pi/2,d,g,ml,mg,sigma,u,w,st)
  call ok(st==0,'co-current primitive recovery')
  call near((1-a)*u,.4_wp,2e-13_wp,'liquid flux round trip')
  call near(a*w,.7_wp,2e-13_wp,'gas flux round trip')
  do i=1,3
    x=-pi/real(2*i,wp)
    call chexal_solve_alpha(.4_wp,.7_wp,rl,rg,x,d,g,ml,mg,sigma,a,st,nroots=n)
    call ok(st==0,'descending cocurrent alpha')
    call chexal_recover_velocities(a,rl,rg,momentum,x,d,g,ml,mg,sigma,u,w,st, &
      ul_hint=.4_wp/(1-a),ug_hint=.7_wp/a)
    call ok(st==0,'descending cocurrent recovery')
    call near((1-a)*u,.4_wp,3e-12_wp,'descending liquid flux round trip')
    call near(a*w,.7_wp,3e-12_wp,'descending gas flux round trip')
  enddo
  call chexal_solve_alpha(.4_wp,1e-40_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st)
  call ok(st==0.and.a>0.and.a<1e-39_wp,'positive alpha tail retained')
  call chexal_solve_alpha(-.1_wp,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,branch=1)
  call ok(st==0,'j=0 high branch solve')
  call near(a,.96079675038909524439_wp,2e-12_wp,'j=0 high branch reference')
  call chexal_solve_alpha(-.1_wp,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,branch=-1)
  call ok(st==0,'j=0 low branch solve')
  call near(a,.21358229200180246032_wp,2e-12_wp,'j=0 low branch reference')
  call chexal_solve_alpha(-.1_wp,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,branch=99)
  call ok(st==CHEXAL_INVALID,'invalid branch rejected')
  call chexal_solve_alpha(.4_wp,.7_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,alpha_hint=-1.0_wp)
  call ok(st==CHEXAL_INVALID,'invalid hint rejected even for unique solution')
  call chexal_ccfl(.1_wp,rl,rg,d,g,ml,mg,sigma,fs,af,st)
  call ok(st==0,'CCFL solver status')
  call near(fs,-.44380008170984717617_wp,2e-10_wp,'CCFL jl high precision independent')
  call near(af,.34341477862448349592_wp,3e-8_wp,'CCFL alpha high precision independent')
  do b=-1,1,2
    call chexal_recover_velocities(af,rl,rg,rl*fs+rg*.1_wp,pi/2,d,g,ml,mg,sigma,u,w,st,branch=b)
    call ok(st==CHEXAL_SINGULAR,'CCFL singular candidate never falls back to cocurrent flow')
  enddo
  call chexal_residual_si(af,rl,rg,fs,.1_wp,pi/2,d,g,ml,mg,sigma,r,st,branch=1)
  call near(r,0.0_wp,1e-12_wp,'CCFL residual')
  h=1e-5_wp
  call chexal_residual_si(af+h,rl,rg,fs,.1_wp,pi/2,d,g,ml,mg,sigma,fp,st,branch=1)
  call chexal_residual_si(af-h,rl,rg,fs,.1_wp,pi/2,d,g,ml,mg,sigma,fm,st,branch=1)
  call near((fp-fm)/(2*h),0.0_wp,3e-8_wp,'CCFL tangency')
  call chexal_solve_alpha(1.01_wp*fs,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,branch=1)
  call ok(st==CHEXAL_CCFL_EXCEEDED,'beyond CCFL rejected')
  call chexal_coefficients_si(.4_wp,rl,rg,1.01_wp*fs,.1_wp,pi/2,d,g,ml,mg,sigma,c,v,st,branch=-1)
  call ok(st==CHEXAL_CCFL_EXCEEDED,'public LOW coefficients beyond CCFL rejected')
  call chexal_coefficients_si(.4_wp,rl,rg,1.01_wp*fs,.1_wp,pi/2,d,g,ml,mg,sigma,c,v,st,branch=1)
  call ok(st==CHEXAL_CCFL_EXCEEDED,'public HIGH coefficients beyond CCFL rejected')
  call chexal_coefficients_si(.4_wp,rl,rg,-1e-25_wp,.1_wp,pi/2,d,g,ml,mg,sigma,c,v,st,branch=-1)
  call ok(st==CHEXAL_OK,'active LOW Eq3-41 finite')
  call near(v,.68178808767336478791_wp,2e-14_wp,'active LOW Eq3-41 independent reference')
  call chexal_solve_alpha(fs,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,branch=1)
  call ok(st==CHEXAL_SINGULAR,'CCFL singular status')
  call chexal_solve_alpha(fs/2,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,nroots=n,roots=roots)
  call ok(st==CHEXAL_BRANCH_REQUIRED.and.n==2,'two physical branches need selection')
  call near(roots(1),.22021468253469897891_wp,2e-10_wp,'low alpha branch reference')
  call near(roots(2),.83984885953778943976_wp,2e-10_wp,'high alpha branch reference')
  do b=-1,1,2
    call chexal_solve_alpha(fs/2,.1_wp,rl,rg,pi/2,d,g,ml,mg,sigma,a,st,branch=b,nroots=n)
    call ok(st==0.and.n==2,'explicit CC branch selection')
    u=(fs/2)/(1-a);w=.1_wp/a;momentum=rl*fs/2+rg*.1_wp
    call chexal_recover_velocities(a,rl,rg,momentum,pi/2,d,g,ml,mg,sigma,x,y,st,branch=b,ul_hint=u,ug_hint=w)
    call ok(st==0,'CC inverse with same branch prior')
    call near(x,u,2e-10_wp,'CC liquid inverse')
    call near(y,w,2e-10_wp,'CC gas inverse')
    call chexal_recover_velocities(a,rl,rg,momentum,pi/2,d,g,ml,mg,sigma,x,y,st,branch=b, &
      ul_hint=u*.999_wp,ug_hint=w*1.001_wp)
    call ok(st==0,'CC continuation from perturbed prior')
    call near(x,u,2e-9_wp,'CC continued liquid velocity')
    call near(y,w,2e-9_wp,'CC continued gas velocity')
    call chexal_recover_velocities(a,rl,rg,momentum,pi/2,d,g,ml,mg,sigma,x,y,st,branch=-b,ul_hint=u,ug_hint=w)
    call ok(st/=0,'wrong physical alpha branch rejected')
    call chexal_recover_velocities(a,rl,rg,momentum,pi/2,d,g,ml,mg,sigma,x,y,st,branch=b)
    call ok(st==CHEXAL_OK,'unique CC inverse on explicit branch without hint')
    call near(x,u,2e-9_wp,'unique CC inverse liquid')
    call near(y,w,2e-9_wp,'unique CC inverse gas')
    call chexal_recover_velocities(a,rl,rg,momentum,pi/2,d,g,ml,mg,sigma,x,y,st)
    call ok(st==CHEXAL_BRANCH_REQUIRED,'AUTO refuses CC alternative without prior')
  enddo
  call chexal_coefficients_si(.2_wp,7*rl,7*rg,.4_wp,.7_wp,pi/2,d,g,7*ml,7*mg,7*sigma,c,v,st)
  call near(c,.83152878676784701845_wp,2e-14_wp,'mass-unit scaling C0')
  call near(v,.54480865106966165970_wp,2e-14_wp,'mass-unit scaling drift')
  call chexal_coefficients_si(.2_wp,rl,rg,1.2_wp,2.1_wp,pi/2,d,9*g,3*ml,3*mg,9*sigma,c,v,st)
  call near(c,.83152878676784701845_wp,2e-14_wp,'velocity similitude C0')
  call near(v,3*.54480865106966165970_wp,2e-14_wp,'velocity similitude drift')
  call chexal_coefficients_si(.2_wp,rl,rg,.4_wp,.7_wp,pi/2,d,g,ml,mg,sigma,c,v,st,signature=s1)
  call chexal_coefficients_si(.2_wp+1e-6_wp,rl,rg,.4_wp,.7_wp,pi/2,d,g,ml,mg,sigma,c,v,st,signature=s2)
  call ok(s1==s2,'same-piece signature for derivative perturbation')
  do i=4,8,2
    call chexal_ccfl(10.0_wp**(-i),rl,rg,d,g,ml,mg,sigma,fs,af,st)
    call ok(st==CHEXAL_NUMERIC,'unresolved tiny-gas CCFL is not accepted as a fold')
    call ok(.not.ieee_is_finite(fs).and..not.ieee_is_finite(af),'unresolved CCFL has no physical outputs')
  enddo
  write(*,'(a,i0,a)')'PASS ',checks,' Chexal 1992 core checks (not experimental validation)'
contains
  subroutine ok(condition,message)
    logical,intent(in)::condition
    character(*),intent(in)::message
    checks=checks+1
    if(.not.condition)then
      write(*,*)'FAIL ',trim(message),' status=',st
      error stop 1
    endif
  end subroutine
  subroutine near(value,expected,tolerance,message)
    real(wp),intent(in)::value,expected,tolerance
    character(*),intent(in)::message
    call ok(ieee_is_finite(value),trim(message)//' finite')
    if(abs(value-expected)>tolerance*max(1.0_wp,abs(expected)))then
      write(*,*)'actual, expected:',value,expected
    endif
    call ok(abs(value-expected)<=tolerance*max(1.0_wp,abs(expected)),message)
  end subroutine
end program test_chexal_1992_core
