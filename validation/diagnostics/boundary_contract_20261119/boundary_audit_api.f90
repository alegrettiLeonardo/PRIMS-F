module boundary_audit_api
  use iso_c_binding, only:c_int,c_double
  use mod_kinds, only:wp
  use module_config, only:sim_config_t
  use module_simulate_pipe_riser, only:pipeline_riser_simulator_t
  use module_boundary_conditions, only:boundary_flux
  use module_find, only:pressure_fun_uj
  use module_lapack_small, only:dgeev_3x3
  use offline_api, only:replay_cfg=>cfg,replay_sim=>sim
  implicit none
  type(sim_config_t),save :: cfg
  type(pipeline_riser_simulator_t),save :: sim
contains
  subroutine ba_setup(par,models) bind(C)
    real(c_double),intent(in)::par(20)
    integer(c_int),intent(in)::models(2)
    cfg=sim_config_t()
    cfg%omega_u=par(1);cfg%omega_rho=par(2);cfg%omega_P=par(3)
    cfg%omega_t=par(4);cfg%Lr=par(5);cfg%Lp=par(6)
    cfg%D_H=par(7);cfg%AREA=acos(-1.0_wp)*cfg%D_H**2/4
    cfg%eps=par(8);cfg%G=par(9);cfg%mu_l=par(10);cfg%mu_g=par(11);cfg%sigma=par(12)
    cfg%nCl=par(13);cfg%nCg=par(14);cfg%nrho_l0=par(15);cfg%nP_l0=par(16)
    cfg%nPs=par(17);cfg%mul=par(18);cfg%mug=par(19);cfg%njl=par(20)
    cfg%closure_pipe=models(1);cfg%closure_riser=models(2)
    cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%maxit=1000;cfg%num_dummy_cells=0
    cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%dissipation_mode='LAGGED'
    cfg%jac_dissipation_fd=.true.
    call sim%init(cfg)
  end subroutine
  subroutine ba_boundary(u,theta,side,replay,f,bc,pb,status) bind(C)
    real(c_double),intent(in)::u(3)
    real(c_double),value::theta
    integer(c_int),value::side,replay
    real(c_double),intent(out)::f(3),bc(3),pb
    integer(c_int),intent(out)::status
    type(sim_config_t)::c
    c=cfg
    if(replay/=0)c=replay_cfg
    call boundary_flux(u,theta,c,side==1,f,status,bc)
    pb=0
    if(status==0)pb=pressure_fun_uj(bc(1),bc(2),c%nCg,c%nCl,c%nrho_l0,c%nP_l0)
  end subroutine
  subroutine ba_point(u,theta,model,replay,v,f,s,a,ds,status) bind(C)
    real(c_double),intent(in)::u(3)
    real(c_double),value::theta
    integer(c_int),value::model,replay
    real(c_double),intent(out)::v(6),f(3),s(3),a(3,3),ds(3,3)
    integer(c_int),intent(out)::status
    if(replay==0)then
      call sim%point_state(u,theta,v,f,s,a,ds,status,closure_id=model)
    else
      call replay_sim%point_state(u,theta,v,f,s,a,ds,status,closure_id=model)
    endif
  end subroutine
  subroutine ba_assemble(n,np,u,old,theta,width,dt,jac,r,j,fin,fout,base,status,idx) bind(C)
    integer(c_int),value::n,np,jac
    real(c_double),intent(in)::u(n,3),old(n,3),theta(n),width(n)
    real(c_double),value::dt
    real(c_double),intent(out)::r(3,n),j(3*n,3*n),fin(3),fout(3),base(2,2)
    integer(c_int),intent(out)::status,idx
    real(wp)::a(3,3,n),b(3,3,n),c(3,3,n),norm
    integer::i,k
    sim%cfg=cfg;sim%cfg%nCells=n;sim%cfg%N_pipe=np
    sim%cfg%thetav=theta;sim%cfg%cell_width=width
    if(allocated(sim%lagged_dissipation))deallocate(sim%lagged_dissipation)
    call sim%prepare_dissipation(old,status)
    if(status/=0)return
    call sim%assemble(u,old,dt,r,a,b,c,norm,fin,fout,status,jac/=0,base,failure_index=idx)
    j=0
    if(status/=0)return
    do i=1,n
      k=3*(i-1)
      j(k+1:k+3,k+1:k+3)=b(:,:,i)
      if(i>1)j(k+1:k+3,k-2:k)=a(:,:,i)
      if(i<n)j(k+1:k+3,k+4:k+6)=c(:,:,i)
    enddo
  end subroutine
end module
