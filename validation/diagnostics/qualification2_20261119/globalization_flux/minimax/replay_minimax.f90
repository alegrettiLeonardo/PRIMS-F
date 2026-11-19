program replay_minimax
  use offline_api
  implicit none
  character(len=1024) :: folder
  real(wp) :: r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),dx(3,np),original(3,np),work(nc,3)
  real(wp) :: norm,fi(3),fo(3),g0,predicted,margin,lam,nt,rr(3,np),aa(3,3,np),bb(3,3,np),cc(3,3,np)
  integer :: status,kind,idx,objective,ls,i,physical_status,margin_status
  logical :: active,passed
  call get_command_argument(1,folder)
  call probe_init(trim(folder),len_trim(folder),1,status)
  if(status/=0)stop 1
  call sim%assemble(trial,old,dt,r,a,b,c,norm,fi,fo,status,.true.)
  if(status/=0.or.norm/=saved_norm.or.any(r/=saved_r))stop 2
  call sim%ti%solve_block_tridiagonal_3x3(a,b,c,-r,original,status)
  if(status/=0)stop 3
  write(*,'(a,es25.16e3)') 'initial_norm ',norm
  do objective=1,2
    dx=original
    call sim%constrained_direction(trial,old,a,b,c,dx,active,status,kind,idx,g0,predicted,minimax=objective==2)
    write(*,'(a,4(i0,1x),l1,2(es25.16e3,1x))') 'direction ',objective,status,kind,idx,active,g0,predicted
    if(status/=0.or..not.active)cycle
    write(*,'(a,i0,3(es25.16e3,1x))') 'cell_delta ',objective,dx(:,idx-lo+1)
    passed=.false.
    do ls=1,18
      lam=2.0_wp**(1-ls)
      work=trial
      do i=lo,hi
        work(i,:)=trial(i,:)+lam*dx(:,i-lo+1)
      enddo
      call sim%assemble(work,old,dt,rr,aa,bb,cc,nt,fi,fo,physical_status,.false.)
      call sim%admissibility_margin(work,kind,idx,margin,margin_status)
      write(*,'(a,4(i0,1x),3(es25.16e3,1x))') 'proposal ',objective,ls,physical_status,margin_status,lam,nt,margin
      if(physical_status/=0.or.margin_status/=0)cycle
      if(margin<0.0_wp)cycle
      if(nt<=cfg%tol_res)then
        passed=.true.
        exit
      endif
    enddo
    write(*,'(a,i0,1x,l1)') 'original_residual_pass ',objective,passed
  enddo
end program replay_minimax
