program restart_refinement
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::u(:,:)
 real(wp)::s,dtmax,duration,dummy(6),step_start
 integer::n,i,k,iu,st
 character(512)::folder,input,arg
 call get_command_argument(1,folder);call get_command_argument(2,input)
 call get_command_argument(3,arg);read(arg,*)n
 call get_command_argument(4,arg);read(arg,*)dtmax
 call get_command_argument(5,arg);read(arg,*)duration
 cfg%omega_u=1.16_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325
 cfg%Lr=12.54625399319708_wp;cfg%Lp=57.4_wp;cfg%CA=3.01818193717446_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
 cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%beta=2*acos(-1._wp)/180
 cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
 cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=206000/cfg%omega_P
 cfg%mul=997.05_wp*.547_wp/(cfg%omega_rho*cfg%omega_u)
 cfg%mug=1.396675899947762_wp/(cfg%omega_rho*cfg%omega_u);cfg%njl=.547_wp/cfg%omega_u;cfg%njg=1.16_wp/cfg%omega_u
 cfg%gas_eos='BAROTROPIC';cfg%closure_pipe=1;cfg%closure_riser=3
 cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp;cfg%newton_max=25
 cfg%dissipation_mode='LAGGED';cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%jac_dissipation_fd=.true.
 cfg%nCells=2*n+4;cfg%N_pipe=n;cfg%num_dummy_cells=2;cfg%quad_points=3
 cfg%dt_start=.001_wp/cfg%omega_t;cfg%dt_min=1e-7_wp/cfg%omega_t;cfg%dt_max=dtmax/cfg%omega_t
 cfg%out_dir=trim(folder);cfg%output_interval_s=.125_wp;cfg%field_interval_s=1;cfg%checkpoint_interval_s=30
 allocate(u(cfg%nCells,3),cfg%vns(cfg%nCells),cfg%thetav(cfg%nCells),cfg%cell_width(cfg%nCells),cfg%face_positions(0:2*n))
 do i=0,n
  cfg%face_positions(i)=cfg%Lp/cfg%Lr*i/n
  cfg%face_positions(n+i)=cfg%Lp/cfg%Lr+real(i,wp)/n
 enddo
 do i=1,2*n
  cfg%vns(i+2)=.5_wp*(cfg%face_positions(i-1)+cfg%face_positions(i))
  cfg%cell_width(i+2)=cfg%face_positions(i)-cfg%face_positions(i-1)
  cfg%thetav(i+2)=-cfg%beta
  if(i>n)cfg%thetav(i+2)=atan((cfg%vns(i+2)*cfg%Lr-cfg%Lp)/cfg%CA)
 enddo
 open(newunit=iu,file=trim(input),status='old')
 do i=3,2*n+2
  read(iu,*)u(i,:)
 enddo
 close(iu)
 do i=1,2
  u(i,:)=u(3,:);u(2*n+2+i,:)=u(2*n+2,:)
  cfg%vns(i)=cfg%vns(3)-(3-i)*cfg%cell_width(3);cfg%vns(2*n+2+i)=cfg%vns(2*n+2)+i*cfg%cell_width(2*n+2)
  cfg%cell_width(i)=cfg%cell_width(3);cfg%cell_width(2*n+2+i)=cfg%cell_width(2*n+2)
  cfg%thetav(i)=cfg%thetav(3);cfg%thetav(2*n+2+i)=cfg%thetav(2*n+2)
 enddo
 call get_command_argument(6,arg)
 if(len_trim(arg)>0)then
  read(arg,*)step_start
  cfg%dt_start=step_start/cfg%omega_t
 endif
 call sim%init(cfg)
 call sim%simulate(u,cfg%vns,cfg%thetav,duration/cfg%omega_t,max_iter=10000000,status=st)
 if(st/=0)stop 3
end program
