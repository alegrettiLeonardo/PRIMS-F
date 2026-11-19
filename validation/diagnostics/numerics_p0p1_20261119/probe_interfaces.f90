program probe_interfaces
  use mod_kinds, only: wp
  use module_config, only: sim_config_t
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use module_boundary_conditions, only: boundary_flux
  use module_lapack_small, only: dgeev_3x3
  use module_steady_state, only: solve_alpha_from_superficials
  implicit none
  integer, parameter :: nc=106, np=102, lo=3, hi=104, ib=53
  type(sim_config_t) :: cfg
  type(pipeline_riser_simulator_t) :: sim
  real(wp) :: old(nc,3), trial(nc,3), state(nc,3), dummy(3), t, dt, saved
  real(wp) :: ps, jl0, jg0, ml0, mg0, theta, width, trace(6), fm(3), fp(3), ubc(3)
  real(wp) :: rl, rg, alpha, mapped(3), v(6), f(3), s(3), a(3,3), ds(3,3)
  integer :: i, idx, iu, out, code, profile, side, reason, unused
  character(len=200) :: header
  character(len=64) :: tag
  open(newunit=iu,file='probe_parameters.dat',status='old')
  read(iu,*)cfg%omega_u,ps,jl0,jg0,ml0,mg0,unused
  close(iu)
  cfg%omega_rho=997.05_wp; cfg%omega_P=101325.0_wp
  cfg%Lr=12.546253993197080_wp; cfg%Lp=57.4_wp; cfg%omega_t=cfg%Lr/cfg%omega_u
  cfg%D_H=.05254_wp; cfg%AREA=acos(-1.0_wp)*cfg%D_H**2/4
  cfg%eps=4.6e-5_wp; cfg%G=9.81_wp
  cfg%mu_l=.89e-3_wp; cfg%mu_g=1.85e-5_wp; cfg%sigma=.0728_wp
  cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1; cfg%nP_l0=1; cfg%nPs=ps/cfg%omega_P
  cfg%mul=ml0/(cfg%omega_rho*cfg%omega_u); cfg%mug=mg0/(cfg%omega_rho*cfg%omega_u)
  cfg%njl=jl0/cfg%omega_u; cfg%njg=jg0/cfg%omega_u
  cfg%maxit=100; cfg%tol=1e-12_wp; cfg%tola=1e-13_wp
  cfg%tol_res=1e-8_wp; cfg%tol_inc=1e-8_wp
  cfg%closure_pipe=1; cfg%closure_riser=3; cfg%gas_eos='BAROTROPIC'
  cfg%junction_model='MATCHED'; cfg%flux_scheme='ROE'; cfg%dissipation_mode='LAGGED'
  cfg%nCells=nc; cfg%num_dummy_cells=2; cfg%N_pipe=51; cfg%quad_points=3
  allocate(cfg%cell_width(nc),cfg%thetav(nc))
  open(newunit=iu,file='last_rejected_state.dat',status='old')
  read(iu,'(a)')header
  read(iu,*)t,dt,saved,reason
  do i=lo,hi
    read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
    if(idx/=i)error stop 'snapshot order'
    cfg%thetav(i)=theta; cfg%cell_width(i)=width
  enddo
  close(iu)
  do i=1,2
    old(i,:)=old(lo,:); trial(i,:)=trial(lo,:)
    old(hi+i,:)=old(hi,:); trial(hi+i,:)=trial(hi,:)
    cfg%thetav(i)=cfg%thetav(lo); cfg%cell_width(i)=cfg%cell_width(lo)
    cfg%thetav(hi+i)=cfg%thetav(hi); cfg%cell_width(hi+i)=cfg%cell_width(hi)
  enddo
  call sim%init(cfg)
  call sim%prepare_dissipation(old,code)
  if(code/=0)error stop 'frozen dissipation refusal'
  open(newunit=out,file='state_operators.dat',status='replace')
  write(out,'(a)')'# blocks: tag; profile,index,closure,status; theta; U; v(alpha,rl,rg,p,ul,ug); F; A rows; eigenreal; eigenimag; eigenvector rows'
  do profile=1,2
    state=old
    if(profile==2)state=trial
    do side=1,2
      i=lo
      if(side==2)i=hi
      tag='inlet_interior'
      if(side==2)tag='outlet_interior'
      call dump_state(trim(tag),state(i,:),cfg%thetav(i),merge(1,3,side==1),profile,i)
      call boundary_flux(state(i,:),cfg%thetav(i),cfg,side==1,f,code,ubc)
      write(*,'(a,3i6)')'boundary status profile,side,status ',profile,side,code
      if(code/=0)cycle
      tag='inlet_boundary'
      if(side==2)tag='outlet_boundary'
      call dump_state(trim(tag),ubc,cfg%thetav(i),merge(1,3,side==1),profile,i)
    enddo
    call dump_state('junction_left_cell',state(ib,:),cfg%thetav(ib),1,profile,ib)
    call dump_state('junction_right_cell',state(ib+1,:),cfg%thetav(ib+1),3,profile,ib+1)
    call sim%point_state(state(ib+1,:),cfg%thetav(ib+1),v,f,s,a,ds,code,closure_id=3)
    if(code/=0)error stop 'right point'
    alpha=solve_alpha_from_superficials(f(1)/v(2),f(2)/v(3),v(2),v(3),cfg%thetav(ib), &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho, &
      cfg%tol,status=code,closure_id=1)
    if(code/=0)error stop 'mapped trace'
    mapped=[(1-alpha)*v(2),alpha*v(3),f(1)+f(2)]
    call dump_state('junction_mapped_right',mapped,cfg%thetav(ib),1,profile,ib+1)
    call sim%face_pair(state(ib,:),state(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code, &
      junction_trace=trace,dissipation_matrix=sim%lagged_dissipation(:,:,ib),left_cell=ib)
    write(*,'(a,2i6)')'junction status profile,status ',profile,code
    if(code/=0)cycle
    rl=cfg%nrho_l0+(trace(1)-cfg%nP_l0)/cfg%nCl**2; rg=trace(1)/cfg%nCg**2
    mapped=[(1-trace(4))*rl,trace(4)*rg,trace(2)+trace(3)]
    call dump_state('junction_left_trace',mapped,cfg%thetav(ib),1,profile,ib)
    mapped=[(1-trace(5))*rl,trace(5)*rg,trace(2)+trace(3)]
    call dump_state('junction_right_trace',mapped,0.0_wp,3,profile,ib+1)
  enddo
  close(out)
contains
  subroutine dump_state(name,u,angle,closure,which,cell)
    character(len=*),intent(in)::name
    real(wp),intent(in)::u(3),angle
    integer,intent(in)::closure,which,cell
    real(wp)::vv(6),ff(3),ss(3),aa(3,3),dd(3,3),wr(3),wi(3),vr(3,3)
    integer::status,j
    call sim%point_state(u,angle,vv,ff,ss,aa,dd,status,closure_id=closure)
    if(status==0)call dgeev_3x3(aa,wr,wi,vr,status)
    write(out,'(a)')name
    write(out,'(4i8)')which,cell,closure,status
    write(out,'(es26.17e3)')angle
    write(out,'(3es26.17e3)')u
    write(out,'(6es26.17e3)')vv
    write(out,'(3es26.17e3)')ff
    if(status/=0)then
      write(*,'(a,1x,a,3i6)')'state/eigen refusal',name,which,cell,status
      return
    endif
    do j=1,3
      write(out,'(3es26.17e3)')aa(j,:)
    enddo
    write(out,'(3es26.17e3)')wr
    write(out,'(3es26.17e3)')wi
    do j=1,3
      write(out,'(3es26.17e3)')vr(j,:)
    enddo
  end subroutine
end program
