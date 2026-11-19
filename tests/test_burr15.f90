program test_burr15
  use mod_kinds, only: wp
  use module_closure_law, only: STRATIFIED, BENDIKSEN
  use module_config, only: sim_config_t
  use module_flow_regime, only: burr15_result_t, burr15_stability, BURR15_OK, BURR15_INVALID_STATE, &
    BURR15_PHASE_ABSENT, BURR15_DIRECTION, BURR15_ANGLE, BURR15_GEOMETRY, BURR15_NUMERIC
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use module_steady_state, only: solve_alpha_from_superficials
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none

  real(wp), parameter :: pi=acos(-1.0_wp), d=.05254_wp, area=pi*d*d/4, g=9.81_wp
  real(wp), parameter :: rl=998.0_wp, rg=2.3_wp, reference=.17927450395772561488_wp
  type(burr15_result_t) :: q, qp, qm
  type(sim_config_t) :: cfg
  type(pipeline_riser_simulator_t) :: sim
  real(wp) :: alpha, h, p, rln, rgn, ug, u(3), x(3), y(3), v(6), vp(6), vm(6)
  real(wp) :: f(3), fp(3), fm(3), s(3), aj(3,3), sj(3,3), work1(3,3), work2(3,3), source(3)
  real(wp) :: trace(6), lo, hi, gl, gg, updated(2), expected
  integer :: count, code, model, k, j

  count=0
  q=burr15_stability(.5_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_OK.and.q%unstable,'literal Burr 15 predicts instability for ug=.2 m/s')
  call near(q%hL_over_D,.5_wp,2e-14_wp,'half-full depth')
  call near(q%critical_ug2,reference**2,5e-15_wp,'independent half-full reference')
  call check(q%ug2<2*reference**2,'test discriminates from placing 1-hL/D inside the root')
  call check((.5_wp*.2_wp)**2<reference**2,'test discriminates real from superficial gas speed')
  call near(q%margin,.2_wp**2-reference**2,5e-15_wp,'dimensional squared-speed margin')
  q=burr15_stability(.5_wp,.1_wp,reference*(1-1e-7_wp),rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_OK.and..not.q%unstable,'below threshold')
  q=burr15_stability(.5_wp,.1_wp,reference*(1+1e-7_wp),rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_OK.and.q%unstable,'above threshold')
  q=burr15_stability(.5_wp,0.0_wp,0.0_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_OK.and..not.q%unstable,'zero velocity with both phases present')
  alpha=2.0_wp/3+sqrt(3.0_wp)/(4*pi)
  q=burr15_stability(alpha,.1_wp,.2_wp,rl,rg,-2*pi/180,d,area,g)
  expected=(.75_wp)**2*(1-rg/rl)*g*cos(2*pi/180)*(3*d/4)*alpha
  call near(q%hL_over_D,.25_wp,2e-14_wp,'analytic quarter-full depth')
  call near(q%critical_ug2,expected,5e-15_wp,'analytic wetted arc at quarter fill')
  qp=burr15_stability(alpha,.1_wp,.2_wp,2*rl,2*rg,-2*pi/180,4*d,16*area,g)
  call near(qp%critical_ug2,4*expected,2e-14_wp,'SI length scaling and density ratio')
  q=burr15_stability(.99_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,area,g)
  qp=burr15_stability(.9999_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_OK.and.qp%status==BURR15_OK,'approach liquid disappearance')
  call check(qp%hL_over_D<q%hL_over_D.and.qp%critical_ug2>q%critical_ug2,'shallow-liquid limit')
  q=burr15_stability(.01_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,area,g)
  qp=burr15_stability(.0001_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_OK.and.qp%status==BURR15_OK,'approach gas disappearance')
  call check(qp%hL_over_D>q%hL_over_D.and.qp%critical_ug2<q%critical_ug2,'shallow-gas limit')
  q=burr15_stability(0.0_wp,.1_wp,0.0_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_PHASE_ABSENT.and..not.ieee_is_finite(q%margin),'gas absent')
  q=burr15_stability(1.0_wp,0.0_wp,1.0_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_PHASE_ABSENT.and..not.ieee_is_finite(q%margin),'liquid absent')
  q=burr15_stability(.5_wp,.1_wp,-1.0_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_DIRECTION.and..not.q%unstable,'signed negative gas velocity is refused')
  q=burr15_stability(.5_wp,-.1_wp,1.0_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_DIRECTION.and..not.q%unstable,'negative liquid velocity is refused')
  q=burr15_stability(.5_wp,.1_wp,.2_wp,rl,rg,.01_wp,d,area,g)
  call check(q%status==BURR15_ANGLE,'uphill geometry is outside this detector')
  q=burr15_stability(.5_wp,.1_wp,.2_wp,rl,rg,-pi/2,d,area,g)
  call check(q%status==BURR15_ANGLE,'vertical geometry is outside this detector')
  q=burr15_stability(.5_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,2*area,g)
  call check(q%status==BURR15_GEOMETRY,'area inconsistent with circular diameter')
  q=burr15_stability(.5_wp,.1_wp,.2_wp,rg,rl,0.0_wp,d,area,g)
  call check(q%status==BURR15_INVALID_STATE,'inverted density order')
  q=burr15_stability(1.01_wp,.1_wp,.2_wp,rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_INVALID_STATE,'invalid alpha')
  q=burr15_stability(.5_wp,.1_wp,ieee_value(0.0_wp,ieee_quiet_nan),rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_INVALID_STATE,'nonfinite input')
  q=burr15_stability(.5_wp,.1_wp,huge(g),rl,rg,0.0_wp,d,area,g)
  call check(q%status==BURR15_NUMERIC,'unrepresentable squared speed')

  cfg%omega_u=.53_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp
  cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1;cfg%nP_l0=1
  cfg%D_H=d;cfg%AREA=area;cfg%eps=4.6e-5_wp;cfg%G=g;cfg%Lr=12.546253993197080_wp
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
  cfg%tol=1e-13_wp;cfg%tola=1e-13_wp;cfg%maxit=100
  cfg%N_pipe=2;cfg%num_dummy_cells=1;cfg%nCells=6
  cfg%closure_riser=BENDIKSEN;cfg%junction_model='MATCHED'
  call sim%init(cfg)
  p=2.2_wp
  rln=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2
  rgn=p/cfg%nCg**2

  ! A horizontal stratified candidate close to the threshold. Solve for its
  ! momentum without letting the detector select the constitutive equation.
  u(:2)=[.5_wp*rln,.5_wp*rgn]
  q=burr15_stability(.5_wp,0.0_wp,0.0_wp,rln*cfg%omega_rho,rgn*cfg%omega_rho,0.0_wp,d,area,g)
  ug=sqrt(q%critical_ug2)/cfg%omega_u
  lo=0.0_wp;hi=1.0_wp
  do k=1,44
    u(3)=(lo+hi)/2
    call sim%point_state(u,0.0_wp,v,f,s,aj,sj,code,.false.,STRATIFIED)
    call check(code==0,'stratified threshold bracket')
    if(v(6)<ug)then
      lo=u(3)
    else
      hi=u(3)
    endif
  enddo
  u(3)=(lo+hi)/2
  do model=STRATIFIED,BENDIKSEN
    sim%cfg%closure_pipe=model
    call sim%point_state(u,0.0_wp,v,f,s,aj,sj,code,closure_id=model)
    call check(code==0,'admissible primitive and analytic derivatives under fixed closure')
    call near(f(1)+f(2),u(3),2e-12_wp,'primitive round-trip momentum')
    call near((1-v(1))*v(2),u(1),2e-14_wp,'primitive round-trip liquid inventory')
    call near(v(1)*v(3),u(2),2e-14_wp,'primitive round-trip gas inventory')
    do j=1,3
      h=1e-6_wp*abs(u(j))
      x=u;y=u;x(j)=x(j)+h;y(j)=y(j)-h
      call sim%point_state(x,0.0_wp,vp,fp,source,work1,work2,code,.false.,model)
      call check(code==0,'same-closure plus perturbation')
      call sim%point_state(y,0.0_wp,vm,fm,source,work1,work2,code,.false.,model)
      call check(code==0,'same-closure minus perturbation')
      call check(maxval(abs(aj(:,j)-(fp-fm)/(2*h))/max(1.0_wp,abs(aj(:,j))))<1e-5_wp, &
        'flux derivative under unchanged closure')
      if(model==STRATIFIED.and.j==3)then
        qp=burr15_stability(vp(1),vp(5)*cfg%omega_u,vp(6)*cfg%omega_u, &
          vp(2)*cfg%omega_rho,vp(3)*cfg%omega_rho,0.0_wp,d,area,g)
        qm=burr15_stability(vm(1),vm(5)*cfg%omega_u,vm(6)*cfg%omega_u, &
          vm(2)*cfg%omega_rho,vm(3)*cfg%omega_rho,0.0_wp,d,area,g)
        call check(qp%status==BURR15_OK.and.qm%status==BURR15_OK,'accepted diagnostic candidates')
        call check(qp%unstable.and..not.qm%unstable,'diagnostic threshold crossed inside FD stencil')
      endif
    enddo
    call sim%face_pair(u,u,0.0_wp,0.0_wp,fm,fp,code,left_cell=2)
    call check(code==0,'same-state horizontal pipe face with mesh section')
    call check(maxval(abs(fm-f)/max(1.0_wp,abs(f)))<2e-12_wp,'equal-state face uses the cell closure')
    call check(maxval(abs(fp-fm))==0,'equal-state common face flux')

    ! A horizontal MATCHED junction, even with equal coordinate inclinations.
    gl=.6_wp*rln;gg=.75_wp*rgn
    call state_from_flux(model,x)
    call state_from_flux(BENDIKSEN,y)
    call sim%face_pair(x,y,0.0_wp,0.0_wp,fm,fp,code,junction_trace=trace,left_cell=3)
    call check(code==0,'MATCHED on horizontal pipeline/riser base')
    call check(all(ieee_is_finite(trace)),'MATCHED trace written')
    call near(trace(1),p,2e-10_wp,'MATCHED pressure continuity')
    call check(maxval(abs(fm(:2)-fp(:2)))==0,'MATCHED phase flux continuity')
    call near(fm(1),gl,2e-10_wp,'positive liquid phase flux')
    call near(fm(2),gg,2e-10_wp,'positive gas phase flux')
    call near(trace(6),gg,2e-10_wp,'positive unprojected gas trace')
    if(model==STRATIFIED)call check(abs(trace(4)-trace(5))>.01_wp,'matched void fractions may differ')
    h=1e-5_wp
    updated=((x(:2)-h*fm(:2)/.2_wp)-x(:2))*.2_wp + ((y(:2)+h*fp(:2)/.07_wp)-y(:2))*.07_wp
    call check(maxval(abs(updated))*cfg%omega_rho*area*cfg%Lr<1e-13_wp,'two-phase conservation at unequal widths')
  enddo
  write(*,'(a,i0,a)') 'PASS ',count,' Burr 15 and fixed pipeline checks'

contains
  subroutine state_from_flux(id, state)
    integer,intent(in)::id
    real(wp),intent(out)::state(3)
    real(wp)::a
    a=solve_alpha_from_superficials(gl/rln,gg/rgn,rln,rgn,0.0_wp,d,area,cfg%eps,g, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=code,closure_id=id)
    call check(code==0,'construct state with prescribed pressure and phase fluxes')
    state=[(1-a)*rln,a*rgn,gl+gg]
  end subroutine
  subroutine check(yes,label)
    logical,intent(in)::yes
    character(*),intent(in)::label
    count=count+1
    if(.not.yes)then
      write(*,'(a,i0,2a)') 'FAIL ',count,': ',label
      error stop 1
    endif
  end subroutine
  subroutine near(a,b,tol,label)
    real(wp),intent(in)::a,b,tol
    character(*),intent(in)::label
    call check(abs(a-b)<=tol,label)
  end subroutine
end program test_burr15
