! The fixture/API are scientific replay support, not a replacement residual.
include '../validation/diagnostics/qualification2_20261119/globalization_flux/offline_api.f90'
program test_active_replay
  use offline_api
  implicit none
  character(len=1024) :: rootdir,outputdir,folder
  character(len=32) :: case_name
  real(wp) :: u(nc,3),locations(nc),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np)
  real(wp) :: norm,fi(3),fo(3),fm(3),fp(3),margin,total_r(2),balance(2),delta_mass(2),factor
  integer :: mode,status,i,checks,case_index,margin_kind,margin_index
  call get_command_argument(1,rootdir)
  call get_command_argument(2,outputdir)
  checks=0
  do case_index=1,2
    case_name='D_SS1_active_terminal'
    margin_kind=1
    margin_index=ib
    if(case_index==2)then
      case_name='D_SS3_active_terminal'
      margin_kind=2
      margin_index=69
    endif
    folder=trim(rootdir)//'/validation/diagnostics/qualification2_20261119/globalization_flux/fixtures/'//trim(case_name)
    do mode=1,2
      call probe_init(trim(folder),len_trim(folder),1,status)
      call check(status==0,'frozen snapshot initialized')
      cfg%dt_start=dt
      cfg%dt_min=dt
      cfg%dt_max=dt
      cfg%newton_max=25
      cfg%newton_globalization='LINE_SEARCH'
      cfg%out_dir=trim(outputdir)//'/'//trim(case_name)//'_line'
      if(mode==2)then
        cfg%newton_globalization='ACTIVE_SET'
        cfg%out_dir=trim(outputdir)//'/'//trim(case_name)//'_active'
      endif
      call execute_command_line('mkdir -p '//trim(cfg%out_dir),exitstat=status)
      call check(status==0,'replay output directory')
      call sim%init(cfg)
      u=old
      locations=0.0_wp
      do i=lo,hi
        locations(i)=sum(cfg%cell_width(lo:i))-cfg%cell_width(i)/2
      enddo
      call sim%simulate(u,locations,cfg%thetav,dt,max_iter=1,status=status)
      if(mode==1)then
        call check(status==6000,'original Newton ray fails on the same accepted state and dt')
      else
        call check(status==0,'active direction completes the same conservative step')
        call sim%assemble(u,old,dt,r,a,b,c,norm,fi,fo,status,.false.)
        call check(status==0.and.norm<=cfg%tol_res,'original nonlinear residual tolerance')
        call check(minval(u(lo:hi,:2))>0.0_wp,'both phase masses remain positive')
        call sim%admissibility_margin(u,margin_kind,margin_index,margin,status)
        call check(status==0.and.margin>=0.0_wp,'actual selected branch margin is nonnegative')
        call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,status, &
          dissipation_matrix=sim%lagged_dissipation(:,:,ib),left_cell=ib)
        call check(status==0.and.all(fm(:2)==fp(:2)),'accepted phase fluxes coincide across MATCHED')
        factor=cfg%omega_rho*cfg%AREA*cfg%Lr
        total_r=0.0_wp
        delta_mass=0.0_wp
        do i=lo,hi
          total_r=total_r+r(:2,i-lo+1)*cfg%cell_width(i)*factor
          delta_mass=delta_mass+(u(i,:2)-old(i,:2))*cfg%cell_width(i)*factor
        enddo
        balance=delta_mass-dt*cfg%omega_t*(fi(:2)-fo(:2))*cfg%omega_rho*cfg%omega_u*cfg%AREA
        call check(maxval(abs(balance-total_r))<1e-18_wp,'integral mass balance equals actual residual')
      endif
    enddo
  enddo
  write(*,'(a,i0,a)') 'PASS ',checks,' active replay checks'
contains
  subroutine check(condition,message)
    logical,intent(in) :: condition
    character(*),intent(in) :: message
    if(.not.condition)then
      write(*,'(a)') 'FAIL '//message
      stop 1
    endif
    checks=checks+1
  end subroutine check
end program test_active_replay
