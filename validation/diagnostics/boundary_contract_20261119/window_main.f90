program boundary_window
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::u(:,:),s(:)
 real(wp)::duration
 integer::iu,n,i,k,st,maxsteps
 character(512)::inputfile,outdir
 call get_command_argument(1,inputfile);call get_command_argument(2,outdir)
 open(newunit=iu,file=trim(inputfile),status='old')
 read(iu,*)n,cfg%num_dummy_cells,cfg%N_pipe,cfg%closure_pipe,cfg%closure_riser,duration, &
   cfg%tol_res,cfg%tol_inc,cfg%newton_max,cfg%jac_dissipation_fd,maxsteps
 read(iu,*)cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%omega_t,cfg%Lr,cfg%Lp,cfg%D_H, &
   cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0, &
   cfg%nPs,cfg%mul,cfg%mug,cfg%njl
 read(iu,*)cfg%dt_start,cfg%dt_max,cfg%dt_min
 cfg%nCells=n;cfg%AREA=acos(-1.0_wp)*cfg%D_H**2/4
 cfg%out_dir=trim(outdir);cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%maxit=1000
 cfg%output_interval_s=.03125_wp;cfg%field_interval_s=.25_wp
 cfg%junction_model='MATCHED';cfg%dissipation_mode='LAGGED';cfg%flux_scheme='ROE'
 cfg%newton_globalization='LINE_SEARCH';cfg%gas_eos='BAROTROPIC'
 allocate(u(n,3),s(n),cfg%thetav(n),cfg%cell_width(n))
 do i=1,n
  read(iu,*)k,cfg%thetav(i),cfg%cell_width(i),u(i,:)
  if(k/=i)error stop 'bad checkpoint index'
 enddo
 close(iu);s=0
 call sim%init(cfg)
 call sim%simulate(u,s,cfg%thetav,duration/cfg%omega_t,max_iter=maxsteps,status=st)
 write(*,*)'SIMULATION_STATUS',st
end program
