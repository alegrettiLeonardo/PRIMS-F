module module_boundary_conditions
  !! Boundary conditions act on boundary faces and ghosts only. There is no
  !! inlet cache and no positivity repair of physical conservative states.
  use mod_kinds,only:wp,i4
  use module_config,only:sim_config_t
  use module_closure_law,only:is_chexal
  use module_find,only:compute_primitive_variables,pressure_fun_uj
  use module_steady_state,only:solve_alpha_from_superficials,PHASE_AUTO,PHASE_LIQUID_ONLY,PHASE_GAS_ONLY
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  public::apply_ghost_cells,apply_boundary_conditions,boundary_flux
contains
  subroutine apply_ghost_cells(U,nd,status)
    real(wp),intent(inout)::U(:,:)
    integer(i4),intent(in)::nd
    integer,intent(out),optional::status
    integer::n,k,first,last
    if(present(status))status=101
    n=size(U,1)
    if(size(U,2)/=3.or.nd<0.or.n<=2*nd)return
    first=nd+1;last=n-nd
    do k=1,nd
      U(first-k,:)=U(first,:)
      U(last+k,:)=U(last,:)
    enddo
    if(present(status))status=0
  end subroutine
  subroutine boundary_flux(Uinterior,theta,cfg,is_inlet,F,status,Ubc)
    real(wp),intent(in)::Uinterior(3),theta
    type(sim_config_t),intent(in)::cfg
    logical,intent(in)::is_inlet
    real(wp),intent(out)::F(3)
    integer,intent(out)::status
    real(wp),intent(out),optional::Ubc(3)
    real(wp)::p,rhol,rhog,alpha,alpha_interior,ul,ug,jl,jg,pip,state(3),pcheck
    real(wp)::ul_interior,ug_interior,ai,rli,rgi,pi
    integer::model,phase,inlet_phase,interior_status
    F=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(Ubc))Ubc=F
    status=101
    if(.not.all(ieee_is_finite(Uinterior)))return
    if(minval(Uinterior(1:2))<0.0_wp.or.sum(Uinterior(1:2))<=0.0_wp)return
    if(.not.all(ieee_is_finite([theta,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0, &
      cfg%omega_rho,cfg%omega_u,cfg%omega_P])))return
    if(min(cfg%nCg,cfg%nCl,cfg%omega_rho,cfg%omega_u,cfg%omega_P)<=0.0_wp)return
    p=pressure_fun_uj(Uinterior(1),Uinterior(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
    if(.not.ieee_is_finite(p))return
    if(p<=0.0_wp)return
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    if(is_inlet)then
      model=cfg%closure_pipe
      ! Inlet pressure follows the interior. Two independently prescribed
      ! phase mass fluxes supply the incoming boundary conditions.
      rhol=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2
      rhog=p/cfg%nCg**2
      if(min(rhol,rhog)<=0.0_wp)return
      alpha_interior=Uinterior(2)/rhog
      if(Uinterior(1)==0.0_wp)alpha_interior=1.0_wp
      if(alpha_interior<0.0_wp.or.alpha_interior>1.0_wp)return
      ! Zero imposed flow alone does not imply that a phase is absent.
      inlet_phase=PHASE_AUTO
      if(Uinterior(2)==0.0_wp.and.cfg%mug==0.0_wp)inlet_phase=PHASE_LIQUID_ONLY
      if(Uinterior(1)==0.0_wp.and.cfg%mul==0.0_wp)inlet_phase=PHASE_GAS_ONLY
      jl=cfg%mul/rhol;jg=cfg%mug/rhog
      alpha=solve_alpha_from_superficials(jl,jg,rhol,rhog,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
        cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=status,closure_id=model, &
        alpha_hint=alpha_interior,phase_hint=inlet_phase)
      if(status/=0)return
      ul=0.0_wp;ug=0.0_wp
      if(alpha<1.0_wp)ul=jl/(1.0_wp-alpha)
      if(alpha>0.0_wp)ug=jg/alpha
      state=[(1.0_wp-alpha)*rhol,alpha*rhog,cfg%mul+cfg%mug]
      F=[cfg%mul,cfg%mug,cfg%mul*ul+cfg%mug*ug+pip*p]
    else
      model=cfg%closure_riser
      rhog=p/cfg%nCg**2
      alpha=Uinterior(2)/rhog
      if(alpha<0.0_wp.or.alpha>1.0_wp)return
      if(is_chexal(model))then
        ! Physical interior velocities carry continuation information to Pout.
        ! They are recovered on the explicitly selected constitutive branch.
        call compute_primitive_variables(Uinterior(1),Uinterior(2),Uinterior(3),theta,cfg%nCg,cfg%nCl, &
          cfg%nrho_l0,cfg%nP_l0,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%njl, &
          cfg%omega_u,cfg%omega_rho,cfg%tol,cfg%tola,cfg%maxit,ai,rli,rgi,pi,ul_interior,ug_interior, &
          phase,status=interior_status,closure_id=model)
        if(interior_status/=0)then
          status=interior_status
          return
        endif
      endif
      ! Impose pressure at the face; retain interior void fraction and total
      ! momentum, then recover velocities under the selected closure at Pout.
      p=cfg%nPs
      if(.not.ieee_is_finite(p))return
      if(p<=0.0_wp)return
      rhol=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2
      rhog=p/cfg%nCg**2
      if(min(rhol,rhog)<=0.0_wp)return
      state=[(1.0_wp-alpha)*rhol,alpha*rhog,Uinterior(3)]
      if(is_chexal(model))then
        call compute_primitive_variables(state(1),state(2),state(3),theta,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0, &
          cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%njl,cfg%omega_u,cfg%omega_rho, &
          cfg%tol,cfg%tola,cfg%maxit,alpha,rhol,rhog,pcheck,ul,ug,phase,status=status,closure_id=model, &
          ul_prev=ul_interior,ug_prev=ug_interior)
      else
      call compute_primitive_variables(state(1),state(2),state(3),theta,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0, &
        cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%njl,cfg%omega_u,cfg%omega_rho, &
        cfg%tol,cfg%tola,cfg%maxit,alpha,rhol,rhog,pcheck,ul,ug,phase,status=status,closure_id=model)
      endif
      if(status/=0)return
      F=[state(1)*ul,state(2)*ug,state(1)*ul**2+state(2)*ug**2+pip*p]
    endif
    if(.not.all(ieee_is_finite(F)).or..not.all(ieee_is_finite(state)))then
      F=ieee_value(0.0_wp,ieee_quiet_nan);status=102;return
    endif
    if(present(Ubc))Ubc=state
    status=0
  end subroutine
  subroutine apply_boundary_conditions(U,inlet_idx,outlet_idx,theta_inlet,theta_outlet,cfg,status)
    real(wp),intent(inout)::U(:,:)
    integer(i4),intent(in)::inlet_idx,outlet_idx
    real(wp),intent(in)::theta_inlet,theta_outlet
    type(sim_config_t),intent(in)::cfg
    integer,intent(out),optional::status
    real(wp)::left(3),right(3),flux(3)
    integer::n,nd,k,code
    if(present(status))status=101
    n=size(U,1);nd=cfg%num_dummy_cells
    if(size(U,2)/=3.or.nd<0.or.n<=2*nd)return
    if(inlet_idx/=nd+1.or.outlet_idx/=n-nd)return
    if(nd==0)then
      if(present(status))status=0
      return
    endif
    call boundary_flux(U(inlet_idx,:),theta_inlet,cfg,.true.,flux,code,left)
    if(code/=0)then
      if(present(status))status=code
      return
    endif
    call boundary_flux(U(outlet_idx,:),theta_outlet,cfg,.false.,flux,code,right)
    if(code/=0)then
      if(present(status))status=code
      return
    endif
    ! Commit only when both boundary states have succeeded.
    do k=1,nd
      U(inlet_idx-k,:)=left
      U(outlet_idx+k,:)=right
    enddo
    if(present(status))status=0
  end subroutine
end module module_boundary_conditions
