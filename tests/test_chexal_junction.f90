program test_chexal_junction
  !! MATCHED phase balances/pressure and replay of a historical gas-negative
  !! numerical trace. No projection and no claim of a physical reverse flow.
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_steady_state,only:solve_alpha_from_superficials
  use module_boundary_conditions,only:boundary_flux
  use module_closure_law,only:lei_fechamento_or_ndim_simp
  use module_chexal_1992,only:CHEXAL_DOMAIN
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp),parameter::pi=acos(-1.0_wp),pipe_theta=-2*pi/180
  real(wp)::x(3),y(3),savedx(3),savedy(3),fm(3),fp(3),trace(6),rl,rg,p,pip,gl,gg,a
  real(wp)::vl(6),vr(6),fl(3),fr(3),src(3),jac(3,3),ds(3,3),p_left,p_right,balance(2),newx(3),newy(3)
  real(wp)::state(3),bc(3),f(3),dt,theta,negative_jg,phase_ul,phase_ug,residual
  integer::model,count,code,nroots,k
  logical::blocked
  count=0
  cfg%omega_u=.53_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp
  cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1;cfg%nP_l0=1
  cfg%D_H=.05254_wp;cfg%AREA=pi*cfg%D_H**2/4;cfg%eps=4.6e-5_wp
  cfg%G=9.81_wp;cfg%Lr=12.546253993197080_wp
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
  cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%quad_points=3;cfg%maxit=100
  cfg%closure_pipe=1;cfg%junction_model='MATCHED';cfg%flux_scheme='ROE'
  pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
  do model=5,7
    cfg%closure_pipe=1;cfg%closure_riser=model;sim%cfg=cfg
    p=2.2_wp;gl=.6_wp;gg=.002_wp
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    call state_from_fluxes(p,gl,gg,pipe_theta,1,x)
    call state_from_fluxes(p,gl,gg,0.0_wp,model,y)
    savedx=x;savedy=y
    call sim%point_state(x,pipe_theta,vl,fl,src,jac,ds,code)
    call check(code==0,'admissible stratified pipe state')
    call sim%point_state(y,0.0_wp,vr,fr,src,jac,ds,code)
    call check(code==0,'admissible Chexal base state')
    call sim%face_pair(x,y,pipe_theta,0.0_wp,fm,fp,code,blocking_active=blocked,junction_trace=trace)
    if(code/=0)write(*,*)'MATCHED code/model=',code,model
    call check(code==0.and..not.blocked,'Chexal MATCHED equilibrium without a barrier')
    call check(all(ieee_is_finite(trace)),'complete physical trace reported')
    call check(maxval(abs(fm(:2)-fp(:2)))==0,'both phases have one common junction flux')
    call near(fm(1),gl,2e-9_wp,'equilibrium liquid mass flux')
    call near(fm(2),gg,2e-9_wp,'equilibrium gas mass flux')
    call near(trace(1),p,2e-9_wp,'equilibrium common pressure')
    call check(abs(trace(4)-trace(5))>.01_wp,'different closure holdups remain distinct')
    p_left=(fm(3)-fm(1)**2/((1-trace(4))*rl)-fm(2)**2/(trace(4)*rg))/pip
    p_right=(fp(3)-fp(1)**2/((1-trace(5))*rl)-fp(2)**2/(trace(5)*rg))/pip
    call near(p_left,p_right,2e-9_wp,'pressure reconstructed independently from both momentum traces')
    call near(p_left,trace(1),2e-9_wp,'reconstructed pressure agrees with common trace')
    call check(maxval(abs(x-savedx))==0.and.maxval(abs(y-savedy))==0,'junction leaves cell inventories unchanged')
    do k=1,3
      dt=real(k,wp)*1e-6_wp
      newx=x-dt*fm/.2_wp;newy=y+dt*fp/.07_wp
      balance=((newx(:2)-x(:2))*.2_wp+(newy(:2)-y(:2))*.07_wp)*cfg%omega_rho*cfg%AREA*cfg%Lr
      call check(maxval(abs(balance))<2e-14_wp,'two-cell phase inventory telescopes at unequal widths')
    enddo

    ! Replay the full-step numerical trace in bg_ss1_stall/summary.json:
    ! p=298980.84692980215Pa, jl=0.5296144723724787m/s,
    ! jg=-6.800580025319073e-11m/s. The original run used Bhagwat; this
    ! replays its trace through each Chexal selector, not its transient solution.
    p=298980.84692980215_wp/cfg%omega_P
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    negative_jg=-6.800580025319073e-11_wp/cfg%omega_u
    a=solve_alpha_from_superficials(.5296144723724787_wp/cfg%omega_u,negative_jg,rl,rg,0.0_wp, &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho, &
      cfg%tol,status=code,closure_id=model,nroots=nroots)
    write(*,'(a,i0,a,i0,a,es15.7)')'HISTORICAL_NEGATIVE_TRACE model=',model,' status=',code,' jg_SI=',negative_jg*cfg%omega_u
    call check(code==CHEXAL_DOMAIN.and.nroots==0.and..not.ieee_is_finite(a),'horizontal reverse-gas trace is domain failure for every branch')
    call check(negative_jg<0,'negative trace is not clipped or replaced by phase absence')

    ! Boundary pure phases are inventory declarations and may not be inferred
    ! merely from a zero imposed phase flux.
    cfg%closure_pipe=model;cfg%nPs=2.2_wp;p=cfg%nPs;theta=.5_wp
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    state=[rl,0.0_wp,.3_wp*rl];cfg%mul=state(3);cfg%mug=0
    call boundary_flux(state,theta,cfg,.true.,f,code,bc)
    call check(code==0.and.f(2)==0.and.bc(2)==0,'Chexal inlet exact absent-gas inventory')
    call boundary_flux(state,theta,cfg,.false.,f,code,bc)
    call check(code==0.and.f(2)==0.and.bc(2)==0,'Chexal outlet exact absent-gas inventory')
    state=[0.0_wp,rg,.8_wp*rg];cfg%mul=0;cfg%mug=state(3)
    call boundary_flux(state,theta,cfg,.true.,f,code,bc)
    call check(code==0.and.f(1)==0.and.bc(1)==0,'Chexal inlet exact absent-liquid inventory')
    call boundary_flux(state,theta,cfg,.false.,f,code,bc)
    call check(code==0.and.f(1)==0.and.bc(1)==0,'Chexal outlet exact absent-liquid inventory')
    cfg%mul=.3_wp*rl;cfg%mug=0
    state=[.7_wp*rl,.3_wp*rg,.3_wp*rl]
    call boundary_flux(state,theta,cfg,.true.,f,code,bc)
    call check(code/=0,'zero gas flow with finite gas inventory is not a fabricated pure phase')
  enddo
  ! The explicit high/low source branch plus a unique compatible inverse
  ! determines the interior; outlet pressure perturbations continue that branch.
  do model=6,7
    cfg%closure_riser=model;cfg%closure_pipe=model;p=2.2_wp;theta=pi/2
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    gl=-.1_wp*rl/cfg%omega_u;gg=.1_wp*rg/cfg%omega_u
    call state_from_fluxes(p,gl,gg,theta,model,state)
    do k=-1,1
      cfg%nPs=p*(1+real(k,wp)*1e-4_wp)
      call boundary_flux(state,theta,cfg,.false.,f,code,bc)
      call check(code==0,'vertical CC outlet selects and continues explicit source branch')
      call check(f(1)<0.and.f(2)>0,'vertical outlet keeps falling liquid and rising gas')
      call near(bc(3),state(3),1e-13_wp,'pressure outlet retains total momentum')
      phase_ul=f(1)/bc(1);phase_ug=f(2)/bc(2)
      p_left=(f(3)-bc(1)*phase_ul**2-bc(2)*phase_ug**2)/pip
      call near(p_left,cfg%nPs,2e-10_wp,'vertical CC outlet imposed pressure')
      rl=cfg%nrho_l0+(cfg%nPs-cfg%nP_l0)/cfg%nCl**2;rg=cfg%nPs/cfg%nCg**2
      a=bc(2)/rg
      call lei_fechamento_or_ndim_simp(a,rl,rg,phase_ul,phase_ug,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
        cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,residual,model,code)
      call check(code==0.and.abs(residual*cfg%omega_u)<2e-9_wp,'outlet satisfies selected Chexal branch')
      if(k==0)then
        call near(f(1),gl,2e-9_wp,'vertical CC equal-pressure liquid round trip')
        call near(f(2),gg,2e-9_wp,'vertical CC equal-pressure gas round trip')
      endif
    enddo
  enddo
  write(*,'(a,i0,a)')'PASS ',count,' Chexal junction checks'
contains
  subroutine state_from_fluxes(p,gl,gg,theta,id,u)
    real(wp),intent(in)::p,gl,gg,theta
    integer,intent(in)::id
    real(wp),intent(out)::u(3)
    real(wp)::rl,rg,a
    integer::status
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    a=solve_alpha_from_superficials(gl/rl,gg/rg,rl,rg,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=status,closure_id=id)
    call check(status==0,'prescribed co-current junction state')
    u=[(1-a)*rl,a*rg,gl+gg]
  end subroutine
  subroutine check(yes,label)
    logical,intent(in)::yes
    character(*),intent(in)::label
    count=count+1
    if(.not.yes)then
      write(*,'(a,i0,a,i0)')'FAIL '//label//' status=',code,' model=',model
      error stop 1
    endif
  end subroutine
  subroutine near(x,y,tolerance,label)
    real(wp),intent(in)::x,y,tolerance
    character(*),intent(in)::label
    call check(ieee_is_finite(x).and.abs(x-y)<=tolerance*max(1.0_wp,abs(y)),label)
  end subroutine
end program
