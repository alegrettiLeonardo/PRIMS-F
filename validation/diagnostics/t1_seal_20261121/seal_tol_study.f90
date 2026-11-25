program seal_tol_study
  !! P1.3: face Newton tolerance of the seal solver. For one fixture state (accepted U),
  !! junction_seal_face_tol in {1e-10,1e-11,1e-12,1e-13}: status, configuration, flux
  !! (SI), difference to the tightest converged tolerance, CPU time per face (10 calls).
  !! Usage: seal_tol_study FIXTURE   (PRIMS_JUNCTION_TRANSMISSION=T1, PRIMS_JUNCTION_SEAL=LIQUID_BRIDGE)
  use offline_api
  implicit none
  character(512)::folder
  real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),fm(3,4),fp(3),tr(6),info(3),t0,t1,tols(4),ref(3)
  integer(c_int)::st
  integer::k,i,code(4)
  tols=[1e-10_wp,1e-11_wp,1e-12_wp,1e-13_wp]
  call get_command_argument(1,folder)
  call probe_init(trim(folder),len_trim(folder),1_c_int,st)
  if(st/=0)then
    write(*,'(a,i0)')'probe_init status ',st;stop 1
  endif
  call probe_state(0_c_int,x,sc,meta);call unpack(x,uu)
  do k=1,4
    cfg%junction_seal_face_tol=tols(k);call sim%init(cfg)
    call sim%prepare_dissipation(uu,code(k))
    call cpu_time(t0)
    do i=1,10
      call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm(:,k),fp,code(k),junction_trace=tr, &
        dissipation_matrix=sim%lagged_dissipation(:,:,ib),left_cell=ib,junction_seal_info=info)
    enddo
    call cpu_time(t1)
    fm(:,k)=fm(:,k)*[cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u**2]
    write(*,'(a,es8.1,a,i5,a,f3.0,a,3es22.14,a,f8.4)')'TOL|',tols(k),'|',code(k),'|',info(1),'|',fm(:,k),'|cpu_s_per_face ',(t1-t0)/10
  enddo
  ref=fm(:,4);if(code(4)/=0)ref=fm(:,3)
  do k=1,4
    if(code(k)==0)write(*,'(a,es8.1,a,3es11.3)')'DIFF|',tols(k),'|abs diff to tightest (Gl,Gg kg/m2/s; F3 Pa) ',abs(fm(:,k)-ref)
  enddo
end program
