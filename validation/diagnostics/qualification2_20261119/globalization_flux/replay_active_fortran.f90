program replay_active_fortran
  use offline_api
  implicit none
  character(len=1024) :: folder
  real(wp) :: u(nc,3),work(nc,3),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),dx(3,np)
  real(wp) :: aa(3,3,np),bb(3,3,np),cc(3,3,np),rr(3,np),norm,nt,fi(3),fo(3),lam,g0,gpred
  real(wp) :: total_r(2),balance(2),inventory_delta(2),mass_factor,step_seconds,final_margin
  integer :: status,kind,idx,it,ls,i,initial_status,location
  logical :: active,accepted
  character(len=32) :: stage
  call get_command_argument(1,folder)
  call probe_init(trim(folder),len_trim(folder),1,status)
  if(status/=0)stop 1
  u=trial
  call sim%assemble(u,old,dt,r,a,b,c,norm,fi,fo,status,.true.)
  if(status/=0.or.norm/=saved_norm.or.any(r/=saved_r))stop 2
  write(*,'(a,es25.16e3)') 'initial_norm ',norm
  do it=1,5
    call sim%assemble(u,old,dt,r,a,b,c,norm,fi,fo,status,.true.)
    if(status/=0)stop 3
    if(norm<=cfg%tol_res)exit
    call sim%ti%solve_block_tridiagonal_3x3(a,b,c,-r,dx,status)
    if(status/=0)stop 4
    work=u
    do i=lo,hi
      work(i,:)=u(i,:)+dx(:,i-lo+1)
    enddo
    call sim%assemble(work,old,dt,rr,aa,bb,cc,nt,fi,fo,initial_status,.false., &
      failure_index=location,failure_stage=stage)
    write(*,'(a,3(i0,1x),a)') 'original_full_proposal ',it,initial_status,location,trim(stage)
    call sim%constrained_direction(u,old,a,b,c,dx,active,status,kind,idx,g0,gpred)
    write(*,'(a,4(i0,1x),l1,2(es25.16e3,1x))') 'active_direction ',it,status,kind,idx,active,g0,gpred
    if(status/=0.or..not.active)stop 5
    accepted=.false.
    do ls=1,18
      lam=2.0_wp**(1-ls)
      work=u
      do i=lo,hi
        work(i,:)=u(i,:)+lam*dx(:,i-lo+1)
      enddo
      call sim%assemble(work,old,dt,rr,aa,bb,cc,nt,fi,fo,status,.false.)
      if(status/=0)cycle
      call sim%admissibility_margin(work,kind,idx,final_margin,status)
      if(status/=0.or.final_margin<0.0_wp)cycle
      if(nt<norm*(1-1e-4_wp*lam).or.nt<=cfg%tol_res)then
        accepted=.true.
        exit
      endif
    enddo
    write(*,'(a,2(i0,1x),l1,2(es25.16e3,1x))') 'accepted ',it,status,accepted,lam,nt
    if(.not.accepted)stop 6
    u=work
  enddo
  if(it>5)stop 7
  call sim%assemble(u,old,dt,r,a,b,c,norm,fi,fo,status,.false.)
  if(status/=0.or.norm>cfg%tol_res)stop 8
  if(minval(u(lo:hi,:2))<=0.0_wp)stop 9
  call sim%admissibility_margin(u,kind,idx,final_margin,status)
  if(status/=0.or.final_margin<0.0_wp)stop 10
  mass_factor=cfg%omega_rho*cfg%AREA*cfg%Lr
  step_seconds=dt*cfg%omega_t
  total_r=0.0_wp
  inventory_delta=0.0_wp
  do i=lo,hi
    total_r=total_r+r(:2,i-lo+1)*cfg%cell_width(i)*mass_factor
    inventory_delta=inventory_delta+(u(i,:2)-old(i,:2))*cfg%cell_width(i)*mass_factor
  enddo
  balance=inventory_delta-step_seconds*(fi(:2)-fo(:2))*cfg%omega_rho*cfg%omega_u*cfg%AREA
  write(*,'(a,es25.16e3)') 'final_norm ',norm
  write(*,'(a,es25.16e3)') 'final_margin ',final_margin
  write(*,'(a,2(es25.16e3,1x))') 'mass_min_nd ',minval(u(lo:hi,1)),minval(u(lo:hi,2))
  write(*,'(a,2(es25.16e3,1x))') 'phase_balance_kg ',balance
  write(*,'(a,2(es25.16e3,1x))') 'integral_residual_kg ',total_r
  write(*,'(a,2(es25.16e3,1x))') 'balance_identity_error_kg ',balance-total_r
  if(maxval(abs(balance-total_r))>1e-18_wp)stop 11
  write(*,'(a)') 'PASS frozen snapshot active-set replay'
end program replay_active_fortran
