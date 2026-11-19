program replay_frozen_terminal
  use offline_api
  use iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  character(len=1024) :: folder
  character(len=32) :: cell_argument
  character(len=32) :: stage
  real(wp) :: r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),dx(3,np),original(3,np),work(nc,3)
  real(wp) :: norm,fi(3),fo(3),g0,predicted,margin,lam,nt,rr(3,np),aa(3,3,np),bb(3,3,np),cc(3,3,np)
  integer :: status,kind,idx,objective,ls,i,physical_status,margin_status,location,unitcsv,focus
  logical :: active
  call get_command_argument(1,folder)
  call get_command_argument(2,cell_argument)
  read(cell_argument,*)focus
  if(focus<lo.or.focus>hi)stop 6
  call probe_init(trim(folder),len_trim(folder),1,status)
  if(status/=0)stop 1
  call sim%assemble(trial,old,dt,r,a,b,c,norm,fi,fo,status,.true.)
  if(status/=0)stop 2
  if(any(r/=saved_r))stop 3
  if(ieee_is_finite(saved_norm))then
    if(norm/=saved_norm)stop 4
  endif
  call sim%ti%solve_block_tridiagonal_3x3(a,b,c,-r,original,status)
  if(status/=0)stop 5
  write(*,'(a,es25.16e3)')'initial_norm ',norm
  open(newunit=unitcsv,file='proposal_cells_exact.csv',status='new')
  write(unitcsv,'(a)')'method,k,cell,lambda,theta,U1,U2,U3,theta_bits,U1_bits,U2_bits,U3_bits,'// &
    'physical_status,margin_status,margin_si,residual_norm,failure_index,failure_stage'
  call sim%admissibility_margin(old,2,focus,margin,margin_status)
  call write_row('OLD',-2,0.0_wp,old,0,margin_status,margin,0.0_wp,0,'NONE')
  call sim%admissibility_margin(trial,2,focus,margin,margin_status)
  call write_row('LINEARIZATION',-1,0.0_wp,trial,0,margin_status,margin,norm,0,'NONE')
  do objective=0,2
    dx=original
    kind=2;idx=focus;active=.true.;status=0
    if(objective>0)then
      call sim%constrained_direction(trial,old,a,b,c,dx,active,status,kind,idx,g0,predicted,minimax=objective==2)
      write(*,'(a,4(i0,1x),l1,2(es25.16e3,1x))')'direction ',objective,status,kind,idx,active,g0,predicted
    endif
    if(status/=0.or..not.active)cycle
    do ls=1,18
      lam=2.0_wp**(1-ls)
      work=trial
      do i=lo,hi
        work(i,:)=trial(i,:)+lam*dx(:,i-lo+1)
      enddo
      call sim%assemble(work,old,dt,rr,aa,bb,cc,nt,fi,fo,physical_status,.false., &
        failure_index=location,failure_stage=stage)
      call sim%admissibility_margin(work,2,focus,margin,margin_status)
      select case(objective)
      case(0)
        call write_row('NEWTON',ls,lam,work,physical_status,margin_status,margin,nt,location,stage)
      case(1)
        call write_row('L2',ls,lam,work,physical_status,margin_status,margin,nt,location,stage)
      case(2)
        call write_row('LINF',ls,lam,work,physical_status,margin_status,margin,nt,location,stage)
      end select
    enddo
  enddo
  close(unitcsv)
contains
  subroutine write_row(method,k,lambda,u,code,mcode,g,norm_value,index,phase)
    character(*),intent(in) :: method,phase
    integer,intent(in) :: k,code,mcode,index
    real(wp),intent(in) :: lambda,u(nc,3),g,norm_value
    integer :: j
    write(unitcsv,'(a,2(a,i0),5(a,es26.17e3),4(a,z16.16),2(a,i0),2(a,es26.17e3),a,i0,2a)') &
      trim(method),',',k,',',focus,',',lambda,',',cfg%thetav(focus),(',',u(focus,j),j=1,3), &
      ',',transfer(cfg%thetav(focus),0_int64),(',',transfer(u(focus,j),0_int64),j=1,3), &
      ',',code,',',mcode,',',g,',',norm_value,',',index,',',trim(phase)
  end subroutine write_row
end program replay_frozen_terminal
