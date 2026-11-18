program test_chexal_derivative_tail
  use mod_kinds,only:wp
  use module_closure_law,only:CHEXAL_1992,CLOSURE_DERIVATIVE, &
    CdUd_drift,Dlei_fechamento_or_ndim_simp,lei_fechamento_or_ndim_simp
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  real(wp),parameter::a=.4_wp,rl=998.0_wp,rg=2.3_wp,ug=.5_wp,theta=1.151159_wp
  real(wp),parameter::d=.05254_wp,area=acos(-1.0_wp)*d*d/4,g=9.81_wp,ml=.001_wp,mg=.000018_wp,sig=.072_wp
  real(wp),parameter::vel(5)=[2.2586e-8_wp,2.25e-8_wp,1e-10_wp,1e-14_wp,1e-100_wp]
  real(wp)::der(5),c0,vgj,expected,r,wr,wu,ul,critical
  integer::k,scale_id,status,checks
  checks=0
  do scale_id=1,3
    wr=1.0_wp;wu=1.0_wp
    if(scale_id==2)then
      wr=1000.0_wp;wu=.7_wp
    else if(scale_id==3)then
      wr=37.0_wp;wu=3.0_wp
    endif
    do k=1,size(vel)
      ul=vel(k)
      call lei_fechamento_or_ndim_simp(a,rl/wr,rg/wr,ul/wu,ug/wu,theta,d,area,0.0_wp,g,ml,mg,sig, &
        wu,wr,1e-12_wp,r,CHEXAL_1992,status)
      call check(status==0.and.ieee_is_finite(r),'base state is admissible, even below old derivative floor')
      call CdUd_drift(a,rl/wr,rg/wr,ul/wu,ug/wu,theta,d,area,0.0_wp,g,ml,mg,sig, &
        wu,wr,1e-12_wp,c0,vgj,CHEXAL_1992,status)
      call check(status==0,'SI coefficients')
      ! Independent local analytic derivative: Reg dominates; C3 is unsaturated.
      expected=a*(1-a)*(c0-vgj*rl*d/(ml*60000.0_wp))
      call derivatives(ul,ml)
      call check(status==0,'one-sided derivative remains on the qualified piece')
      call check(abs(der(4)-expected)<3e-10_wp,'derivative versus analytic local expression')
      if(k==1)call check(abs(der(4)-.18369089945766421933792045146128424325_wp)<3e-10_wp, &
        'independent 60-digit derivative at observed velocity scale')
      if(k==4)call check(abs(der(4)-.18369089829713939780027525161302797245_wp)<3e-10_wp, &
        'independent 60-digit derivative below centered-stencil range')
    enddo
    call derivatives(0.0_wp,ml)
    call check(status==CLOSURE_DERIVATIVE,'exact sign boundary has no invented bilateral derivative')
    critical=60000.0_wp*log(4.0_wp)*ml/(rl*d*(1-a))
    call derivatives(critical,ml)
    call check(status==CLOSURE_DERIVATIVE,'ordinary C3 kink stays rejected')
    critical=60000.0_wp*log(4.0_wp)*1e-9_wp/(rl*d*(1-a))
    call derivatives(critical,1e-9_wp)
    call check(status==CLOSURE_DERIVATIVE,'near-zero constitutive kink cannot masquerade as a sign boundary')
    ! The fixed 12.5% probe falsely rejected this smooth state although
    ! the old centered stencil succeeded. Reg still dominates at 0.95*ulcrit.
    critical=rg*a*.005_wp*ml/(rl*(1-a)*mg)
    ul=.95_wp*critical
    call CdUd_drift(a,rl/wr,rg/wr,ul/wu,.005_wp/wu,theta,d,area,0.0_wp,g,ml,mg,sig, &
      wu,wr,1e-12_wp,c0,vgj,CHEXAL_1992,status)
    call check(status==0,'smooth point near max-Reynolds switch has valid coefficients')
    expected=a*(1-a)*(c0-vgj*rl*d/(ml*60000.0_wp))
    call derivatives(ul,ml,.005_wp)
    call check(status==0,'adaptive bilateral guard preserves smooth near-switch state')
    call check(abs(der(4)-expected)<3e-10_wp,'near-switch derivative agrees with local analytic expression')
    call derivatives(critical,ml,.005_wp)
    call check(status==CLOSURE_DERIVATIVE,'exact low-speed max-Reynolds kink stays rejected')
  enddo
  write(*,'(a,i0,a)')'PASS Chexal velocity-tail derivatives: ',checks,' assertions'
contains
  subroutine derivatives(velocity,viscosity,gas_speed)
    real(wp),intent(in)::velocity,viscosity
    real(wp),intent(in),optional::gas_speed
    real(wp)::gas_velocity
    gas_velocity=ug
    if(present(gas_speed))gas_velocity=gas_speed
    call Dlei_fechamento_or_ndim_simp(a,rl/wr,rg/wr,velocity/wu,gas_velocity/wu,theta,d,area,0.0_wp,g,viscosity,mg,sig, &
      wu,wr,1e-12_wp,der(1),der(3),der(2),der(5),der(4),CHEXAL_1992,status)
  end subroutine
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    if(.not.ok)then
      write(*,*)'FAIL ',label,' status=',status,' derivative=',der(4),' expected=',expected
      error stop 1
    endif
    checks=checks+1
  end subroutine
end program
