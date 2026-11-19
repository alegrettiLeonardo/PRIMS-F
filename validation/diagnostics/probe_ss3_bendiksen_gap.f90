program probe_ss3_bendiksen_gap
  ! An algebraic domain diagnostic, not a time integration or a fitted model.
  use mod_kinds, only: wp
  use module_find, only: compute_primitive_variables, pressure_fun_uj
  implicit none
  real(wp), parameter :: pi=acos(-1.0_wp), wrho=997.05_wp, wu=.48_wp, wp0=101325._wp
  real(wp), parameter :: dh=.05254_wp, grav=9.81_wp, cl=1498._wp, cg=343._wp
  real(wp) :: old(3),trial(3),r(3),theta,width,tt,dt,norm,un(3)
  real(wp) :: alpha,rhol,rhog,p,ul,ug,ncg,ncl,vref,c0l,c0h,udl,udh,al,ah,jcrit
  real(wp) :: mlow,mhigh,mlocal,jlow,jhigh,moments(7)
  integer :: iu,idx,reason,code,index,k,out
  character(160) :: header
  character(28), parameter :: labels(7)=[character(28):: &
    'accepted','below_low_edge','just_inside_gap','gap_midpoint','just_below_high_edge','above_high_edge','trial']
  ncg=cg*sqrt(wrho/wp0); ncl=cl*sqrt(wrho/wp0)
  open(newunit=iu,file='ss3_late_snapshot.dat',status='old')
  read(iu,'(a)')header
  read(iu,*)tt,dt,norm,reason
  do
    read(iu,*,end=99)idx,theta,width,old,trial,r
    if(idx==49)exit
  enddo
  close(iu)
  p=pressure_fun_uj(old(1),old(2),ncg,ncl,1._wp,1._wp)
  rhol=wrho*(1._wp+(p-1._wp)/ncl**2); rhog=wrho*p/ncg**2
  alpha=old(2)*wrho/rhog
  vref=sqrt(grav*dh);jcrit=3.5_wp*vref
  c0l=1.05_wp+.15_wp*sin(theta);c0h=1.2_wp
  udl=vref*(.35_wp*sin(theta)+.54_wp*cos(theta));udh=vref*.35_wp*sin(theta)
  al=rhol-alpha*(rhol-rhog)*c0l;ah=rhol-alpha*(rhol-rhog)*c0h
  mlow=al*jcrit-alpha*(rhol-rhog)*udl
  mhigh=ah*jcrit-alpha*(rhol-rhog)*udh
  moments=[old(3)*wrho*wu,mlow-1e-4_wp,mlow+1e-4_wp,(mlow+mhigh)/2, &
    mhigh-1e-4_wp,mhigh+1e-4_wp,trial(3)*wrho*wu]
  write(*,'(a,3es25.16)')'alpha, p_Pa, theta_rad ',alpha,p*wp0,theta
  write(*,'(a,3es25.16)')'m_low_SI, m_high_SI, gap_SI ',mlow,mhigh,mhigh-mlow
  open(newunit=out,file='ss3_bendiksen_gap_probe.csv',status='replace')
  write(out,'(a)')'state,momentum_SI,low_Fr,high_Fr,primitive_status,actual_Fr'
  do k=1,size(moments)
    un=old;un(3)=moments(k)/(wrho*wu)
    if(k==7)un=trial
    ! Recompute the two analytic candidates from the actual inventories.
    p=pressure_fun_uj(un(1),un(2),ncg,ncl,1._wp,1._wp)
    rhol=wrho*(1._wp+(p-1._wp)/ncl**2);rhog=wrho*p/ncg**2;alpha=un(2)*wrho/rhog
    al=rhol-alpha*(rhol-rhog)*c0l;ah=rhol-alpha*(rhol-rhog)*c0h
    mlocal=un(3)*wrho*wu
    jlow=(mlocal+alpha*(rhol-rhog)*udl)/al
    jhigh=(mlocal+alpha*(rhol-rhog)*udh)/ah
    call compute_primitive_variables(un(1),un(2),un(3),theta,ncg,ncl,1._wp,1._wp,dh,pi*dh**2/4, &
      4.6e-5_wp,grav,.89e-3_wp,1.85e-5_wp,.0728_wp,.28_wp/wu,wu,wrho,1e-13_wp,1e-14_wp,100, &
      alpha,rhol,rhog,p,ul,ug,index,status=code,closure_id=2)
    if(code==0)then
      write(out,'(a,",",3(es25.16,","),i0,",",es25.16)')trim(labels(k)),mlocal,jlow/vref,jhigh/vref, &
        code,((1-alpha)*ul+alpha*ug)*wu/vref
    else
      write(out,'(a,",",3(es25.16,","),i0,",null")')trim(labels(k)),mlocal,jlow/vref,jhigh/vref,code
    endif
    write(*,'(a,1x,i0)')trim(labels(k))//' primitive_status=',code
    if(k==4.and.code/=4)error stop 'mid-gap must return PRIMITIVE_NO_ROOT=4'
    if((k==1.or.k==2.or.k==6).and.code/=0)error stop 'valid exterior point must recover'
  enddo
  close(out)
  write(*,'(a)')'PASS: positive phase inventories do not guarantee a Bendiksen velocity root across Fr=3.5.'
  stop
99 error stop 'cell 49 absent'
end program probe_ss3_bendiksen_gap
