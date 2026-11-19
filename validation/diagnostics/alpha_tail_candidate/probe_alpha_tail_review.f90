program probe_alpha_tail_review
  use mod_kinds, only: wp
  use module_steady_state, only: solve_alpha_from_superficials
  use module_closure_law, only: BHAGWAT, lei_fechamento_or_ndim_simp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  real(wp), parameter :: pi=acos(-1.0_wp),diam=.05254_wp,area=pi*diam**2/4
  real(wp), parameter :: p=284053.9719206488_wp
  real(wp) :: jl,jg,theta,wr,wu,tol,rl,rg,alpha,hint0,hint1,residual
  integer :: label,ios,code,nroots,s0,s1,rs
  do
    read(*,*,iostat=ios)label,jl,jg,theta,wr,wu,tol
    if(ios/=0)exit
    rl=(997.05_wp+(p-101325.0_wp)/1498.0_wp**2)/wr
    rg=p/343.0_wp**2/wr
    alpha=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,diam,area,4.6e-5_wp,9.81_wp, &
      .00089_wp,1.85e-5_wp,.0728_wp,wu,wr,tol,status=code,closure_id=BHAGWAT,nroots=nroots)
    hint0=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,diam,area,4.6e-5_wp,9.81_wp, &
      .00089_wp,1.85e-5_wp,.0728_wp,wu,wr,tol,status=s0,closure_id=BHAGWAT,alpha_hint=0.0_wp)
    hint1=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,theta,diam,area,4.6e-5_wp,9.81_wp, &
      .00089_wp,1.85e-5_wp,.0728_wp,wu,wr,tol,status=s1,closure_id=BHAGWAT,alpha_hint=1.0_wp)
    residual=ieee_value(0.0_wp,ieee_quiet_nan);rs=-1
    if(code==0.and.ieee_is_finite(alpha))then
      call lei_fechamento_or_ndim_simp(alpha,rl,rg,jl/(1.0_wp-alpha)/wu,jg/alpha/wu,theta,diam, &
        area,4.6e-5_wp,9.81_wp,.00089_wp,1.85e-5_wp,.0728_wp,wu,wr,tol,residual,BHAGWAT,rs)
    endif
    write(*,'(6(i0,1x),4(es25.16e3,1x))')label,code,nroots,s0,s1,rs,alpha,hint0,hint1,residual
  enddo
end program
