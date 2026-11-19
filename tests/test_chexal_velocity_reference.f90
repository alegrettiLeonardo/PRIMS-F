program test_chexal_velocity_reference
  !! Independent references: tests/fixtures/chexal_velocity_derivatives.json.
  !! These partial derivatives are evaluated off the solution surface R=0.
  use mod_kinds,only:wp
  use module_closure_law,only:CHEXAL_1992,Dlei_fechamento_or_ndim_simp
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  real(wp),parameter::a=.4_wp,rl=998.0_wp,rg=2.3_wp,ug=.5_wp,theta=1.151159_wp
  real(wp),parameter::d=.05254_wp,g=9.81_wp,ml=.001_wp,mg=.000018_wp,sigma=.072_wp
  real(wp),parameter::pi=acos(-1.0_wp),area=pi*d*d/4
  real(wp),parameter::vel(2)=[2.2586e-8_wp,1e-14_wp]
  real(wp),parameter::reference(2)=[.18369089945766421933792045146128424325_wp, &
                                   .18369089829713939780027525161302797245_wp]
  real(wp),parameter::uscale(3)=[1.0_wp,3.0_wp,.2_wp],rscale(3)=[1.0_wp,1000.0_wp,10.0_wp]
  real(wp)::fa,fg,fl,fug,ful,ul,wu,wr
  integer::k,j,st,n
  n=0
  do k=1,3
    wu=uscale(k);wr=rscale(k)
    do j=1,2
      ul=vel(j)
      call derivative(ul)
      call check(st==0,'positive near-zero velocity derivative succeeds')
      call check(ieee_is_finite(ful),'liquid velocity derivative finite')
      call check(abs(ful-reference(j))<2e-9_wp,'independent 60-digit derivative and scale invariance')
    enddo
  enddo
  wu=1;wr=1
  call derivative(0.0_wp)
  call check(st/=0,'exact sign boundary retains derivative status')
  call derivative(-1e-14_wp)
  call check(st/=0,'inclined countercurrent state remains rejected')
  ul=(rg*a*ug*d/mg)*ml/(rl*(1-a)*d)
  call derivative(ul)
  call check(st/=0,'actual max-Reynolds kink remains rejected')
  write(*,'(a,i0,a)')'PASS ',n,' independent Chexal velocity derivative checks'
contains
  subroutine derivative(u)
    real(wp),intent(in)::u
    call Dlei_fechamento_or_ndim_simp(a,rl/wr,rg/wr,u/wu,ug/wu,theta,d,area,0.0_wp,g,ml,mg,sigma, &
      wu,wr,1e-12_wp,fa,fg,fl,fug,ful,CHEXAL_1992,st)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(*),intent(in)::msg
    n=n+1
    if(.not.ok)then
      write(*,*)'FAIL ',trim(msg),' status=',st
      error stop 1
    endif
  end subroutine
end program test_chexal_velocity_reference
