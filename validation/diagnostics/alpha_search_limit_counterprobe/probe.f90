program probe_alpha_search_limit
  use mod_kinds, only: wp
  use module_closure_law, only: lei_fechamento_or_ndim_simp,BHAGWAT
  use module_steady_state, only: solve_alpha_from_superficials
  implicit none
  real(wp),parameter :: wu=0.53_wp,wrho=997.05_wp,d=0.05254_wp,pi=acos(-1.0_wp)
  real(wp),parameter :: area=pi*d*d/4,rough=4.6e-5_wp,g=9.81_wp,mul=8.9e-4_wp,mug=1.85e-5_wp
  real(wp),parameter :: sigma=0.0728_wp,tol=1.e-12_wp
  real(wp) :: th,p,jl,jg,aref,rl,rg,alpha,alpha_hint_result,residual
  integer :: i,unit,code,nroots,closure_status,hint_status
  open(newunit=unit,file='inputs.dat',status='old',action='read')
  do i=1,4
    read(unit,*)th,p,jl,jg,aref
    rl=(997.05_wp+(p-101325.0_wp)/1498.0_wp**2)/wrho
    rg=p/343.0_wp**2/wrho
    alpha=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,th,d,area,rough,g,mul,mug,sigma,wu,wrho,tol, &
      status=code,closure_id=BHAGWAT,nroots=nroots)
    call lei_fechamento_or_ndim_simp(aref,rl,rg,jl/(1-aref)/wu,jg/aref/wu,th,d,area,rough,g, &
      mul,mug,sigma,wu,wrho,tol,residual,closure_id=BHAGWAT,status=closure_status)
    alpha_hint_result=solve_alpha_from_superficials(jl/wu,jg/wu,rl,rg,th,d,area,rough,g,mul,mug,sigma,wu,wrho,tol, &
      status=hint_status,closure_id=BHAGWAT,alpha_hint=aref)
    write(*,'(3i6,es26.17,i6,es26.17,i6)')i,code,nroots,alpha,closure_status,residual*wu,hint_status
  enddo
  close(unit)
end program
