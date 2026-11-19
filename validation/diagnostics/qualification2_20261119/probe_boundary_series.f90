program probe_boundary_series
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_boundary_conditions,only:boundary_flux
  use module_lapack_small,only:dgeev_3x3
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp)::t,angle,u(3),bc(3),f(3),v(6),s(3),a(3,3),ds(3,3),wr(3),wi(3),vr(3,3),pg(3)
  real(wp)::alpha,rl,rg,p,ul,ug,dummy,nan,den
  integer::i,j,unit_in,unit_out,ios,st,side,model,kind
  character(512)::input_name,output_name
  call get_command_argument(1,input_name)
  call get_command_argument(2,output_name)
  open(newunit=unit_in,file=trim(input_name),status='old')
  read(unit_in,*)cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l, &
    cfg%mu_g,cfg%sigma,cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0,cfg%nPs,cfg%mul,cfg%mug,cfg%njl
  cfg%Lr=12.546253993197080_wp;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp
  cfg%closure_pipe=1;cfg%closure_riser=3
  call sim%init(cfg)
  open(newunit=unit_out,file=trim(output_name),status='new')
  write(unit_out,'(a)')'time_s,side,kind,status,alpha,ul_m_s,ug_m_s,e1_m_s,e2_m_s,e3_m_s,'// &
    'i1_m_s,i2_m_s,i3_m_s,c11,c12,c13,c21,c22,c23,c31,c32,c33'
  nan=ieee_value(0.0_wp,ieee_quiet_nan)
  do
    read(unit_in,*,iostat=ios)t,side,angle,alpha,rl,rg,p,ul,ug
    if(ios<0)exit
    if(ios>0)error stop 'Malformed accepted boundary state'
    u=[(1-alpha)*rl,alpha*rg,(1-alpha)*rl*ul+alpha*rg*ug]
    model=cfg%closure_pipe
    if(side==2)model=cfg%closure_riser
    do kind=1,2
      bc=u
      st=0
      if(kind==2)call boundary_flux(u,angle,cfg,side==1,f,st,bc)
      wr=nan;wi=nan;a=nan;v=nan
      if(st==0)call sim%point_state(bc,angle,v,f,s,a,ds,st,closure_id=model)
      if(st==0)call dgeev_3x3(a,wr,wi,vr,st)
      if(st==0)then
        den=v(2)*v(1)/cfg%nCg**2+v(3)*(1-v(1))/cfg%nCl**2
        pg=[v(3),v(2),0.0_wp]/den
        ! Constraint derivatives projected on the three right eigenvectors.
        a(3,:)=pg
        do i=1,3
          a(i,:)=a(i,:)/norm2(a(i,:))
        enddo
        a=matmul(a,vr)
      endif
      write(unit_out,'(es25.16e3,",",3(i0,","),18(es25.16e3,:,","))')t,side,kind,st, &
        v(1),v(5:6)*cfg%omega_u,wr*cfg%omega_u,wi*cfg%omega_u,((a(i,j),j=1,3),i=1,3)
    enddo
  enddo
  close(unit_in);close(unit_out)
end program
