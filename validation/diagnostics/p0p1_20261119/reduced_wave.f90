program reduced_wave
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::u(:,:)
 real(wp)::x,p,rho,base,c,velocity,step
 integer::n,phase,i,st
 character(512)::folder,arg
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)n
 call get_command_argument(3,arg);read(arg,*)phase
 cfg%omega_u=1;cfg%omega_rho=997.05_wp;cfg%omega_P=101325;cfg%Lr=10;cfg%omega_t=10
 cfg%Lp=5;cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=0
 cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
 cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=sqrt(287*293.15_wp)*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=206000/cfg%omega_P
 cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp;cfg%newton_max=25
 cfg%nCells=n+4;cfg%N_pipe=n/2;cfg%num_dummy_cells=2
 cfg%single_phase_mode='REDUCED';cfg%source_scheme='HYDROSTATIC';cfg%outlet_reservoir_phase=phase
 c=1498;if(phase==2)c=sqrt(287*293.15_wp)
 step=.2_wp*10/(n*c);cfg%dt_start=step/10;cfg%dt_max=step/10;cfg%dt_min=1e-10_wp
 cfg%out_dir=trim(folder);cfg%output_interval_s=step;cfg%field_interval_s=step;cfg%checkpoint_interval_s=1
 allocate(u(n+4,3),cfg%vns(n+4),cfg%thetav(n+4),cfg%cell_width(n+4),cfg%z_faces(0:n));u=0
 cfg%thetav=0;cfg%cell_width=1._wp/n;cfg%z_faces=0
 base=997.05_wp+(206000-101325)/1498._wp**2
 if(phase==2)base=206000/(287*293.15_wp)
 do i=1,n+4
  x=10*(i-2.5_wp)/n;cfg%vns(i)=x/10;p=206000+100*exp(-.5_wp*((x-2.5_wp)/.5_wp)**2)
  rho=997.05_wp+(p-101325)/1498._wp**2
  if(phase==2)rho=p/(287*293.15_wp)
  velocity=(p-206000)/(base*c)
  u(i,phase)=rho/cfg%omega_rho;u(i,3)=rho*velocity/cfg%omega_rho
 enddo
 call sim%init(cfg)
 call sim%simulate(u,cfg%vns,cfg%thetav,11/(c*cfg%omega_t),max_iter=100000,status=st)
 if(st/=0)stop 3
end program
