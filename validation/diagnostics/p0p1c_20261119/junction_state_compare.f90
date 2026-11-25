program junction_state_compare
 !! P0.1/P0.2 diagnostic at a saved terminal state (51+51 mesh, offline_api).
 !! For the accepted state (0) and the last refused Newton trial (1):
 !!  CELL  rows for the two cells adjacent to the junction:
 !!        alpha_g, rho_g, u_g, j_g, U2 (gas mass per volume), p, U3, closure,
 !!        Bhagwat drift velocity Ud of that state (physical velocity scale);
 !!  FACE  rows, with EXACTLY the same left/right cell states, for
 !!        Roe/CENTER (lagged dissipation of the run), Roe/STEADY,
 !!        E1/CENTER and E1/STEADY: status, p*, Gl*, Gg*, alpha_pipe*,
 !!        alpha_riser*, and the convective (centred) part of the Roe vector.
 !! All in SI. Nothing in the run is modified; E1 is evaluated, not run.
 !! Arguments: fixture folder, gas sound speed (m/s) of the run EOS, and
 !! optionally the junction_flux and junction_trace of the run that produced the state (its lagged
 !! dissipation and probe status are then evaluated with that operator).
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_closure_law,only:CdUd_bhagwat
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),cgas,fm(3),fp(3),tr(6),v(6),f(3),s(3),a(3,3),d(3,3),cd,ud,pipf
 integer(c_int)::st
 integer::which,k,i,code,cid
 character(512)::folder,arg
 character(16)::runflux,runtrace
 character(16)::flux(4)=['ROE         ','ROE         ','PHASE_UPWIND','PHASE_UPWIND']
 character(16)::trc(4)=['CENTER      ','STEADY      ','CENTER      ','STEADY      ']
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 runflux='ROE'
 runtrace='CENTER'
 if(command_argument_count()>2)call get_command_argument(3,runflux)
 if(command_argument_count()>3)call get_command_argument(4,runtrace)
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 ! probe_init evaluates the saved state with the default operator; its status
 ! is reported but only the operator of the run decides admissibility below.
 write(*,'(a,i6)')'PROBE_DEFAULT_STATUS|',st
 st=0
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P);cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp
 cfg%junction_flux=runflux;cfg%junction_trace=runtrace
 ! Geometry needed only by the opt-in STEADY trace (not stored in the fixture):
 ! faces/centres from the stored cell widths, catenary constant of the rig.
 cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
 if(.not.allocated(cfg%face_positions))allocate(cfg%face_positions(0:hi-lo+1))
 if(.not.allocated(cfg%vns))allocate(cfg%vns(nc))
 cfg%face_positions(0)=0
 do i=lo,hi
  cfg%face_positions(i-lo+1)=cfg%face_positions(i-lo)+cfg%cell_width(i)
  cfg%vns(i)=cfg%face_positions(i-lo)+.5_wp*cfg%cell_width(i)
 enddo
 do i=1,lo-1
  cfg%vns(i)=cfg%vns(lo)-(lo-i)*cfg%cell_width(lo);cfg%vns(hi+i)=cfg%vns(hi)+i*cfg%cell_width(hi)
 enddo
 write(*,'(a,es24.16,a,es24.16)')'GEOM|junction face s/Lr ',cfg%face_positions(cfg%N_pipe),' expected Lp/Lr ',cfg%Lp/cfg%Lr
 call sim%init(cfg);call sim%prepare_dissipation(old,code)
 write(*,'(a,2i6,a,i4)')'INIT|',st,code,'|junction left cell',ib
 if(st/=0.or.code/=0)stop 2
 pipf=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
 do which=0,1
  call probe_state(int(which,c_int),x,sc,meta);call unpack(x,uu)
  do i=ib,ib+1
   cid=merge(cfg%closure_pipe,cfg%closure_riser,i==ib)
   call sim%point_state(uu(i,:),cfg%thetav(i),v,f,s,a,d,code,.false.,cid)
   ud=-1
   if(code==0)then
    call CdUd_bhagwat(v(1),v(2),v(3),v(5),v(6),cfg%thetav(i),cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
      cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,cd,ud,k)
    if(k/=0)ud=-1
   endif
   write(*,'(a,i1,a,i4,a,i3,a,i4,a,10es16.8)')'CELL|',which,'|',i,'|',cid,'|',code,'|',v(1),v(3)*cfg%omega_rho, &
     v(6)*cfg%omega_u,v(1)*v(6)*cfg%omega_u,uu(i,2)*cfg%omega_rho,v(4)*cfg%omega_P,uu(i,3)*cfg%omega_rho*cfg%omega_u, &
     f(2)*cfg%omega_rho*cfg%omega_u,cfg%thetav(i),ud
  enddo
  do k=1,4
   cfg%junction_flux=flux(k);cfg%junction_trace=trc(k);call sim%init(cfg)
   if(k==1)then
    cfg%junction_flux=runflux;cfg%junction_trace=runtrace;call sim%init(cfg)
    call sim%prepare_dissipation(old,code)
    cfg%junction_flux=flux(k);cfg%junction_trace=trc(k);call sim%init(cfg)
    call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code, &
      junction_trace=tr,dissipation_matrix=sim%lagged_dissipation(:,:,ib),left_cell=ib)
   else
    call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code,junction_trace=tr,left_cell=ib)
   endif
   write(*,'(a,i1,a,a,a,a,a,i6,a,6es16.8)')'FACE|',which,'|',trim(flux(k)),'|',trim(trc(k)),'|',code,'|', &
     tr(1)*cfg%omega_P,tr(2)*cfg%omega_rho*cfg%omega_u,tr(3)*cfg%omega_rho*cfg%omega_u,tr(4),tr(5), &
     tr(6)*cfg%omega_rho*cfg%omega_u
  enddo
  cfg%junction_flux=runflux;cfg%junction_trace=runtrace;call sim%init(cfg);call sim%prepare_dissipation(old,code)
 enddo
end program
