program test_active_set
  use mod_kinds, only: wp
  use module_active_set, only: active_set_direction
  use module_time_step, only: time_integrator_t
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(time_integrator_t) :: solver
  real(wp) :: a(3,3,2),b(3,3,2),c(3,3,2),scale(3,2),g(3,2),d(3,2),original(3,2)
  real(wp) :: residual(3,2),stationarity(3,2),mu,margin
  integer :: i,status,checks
  logical :: active
  checks=0
  a=0.0_wp;b=0.0_wp;c=0.0_wp
  do i=1,3
    b(i,i,:)=2.0_wp+i
    a(i,i,2)=-0.2_wp
    c(i,i,1)=0.1_wp
  enddo
  b(1,2,:)=-0.1_wp
  b(2,3,:)=0.3_wp
  scale=1.0_wp
  g=0.0_wp;g(1,1)=1.0_wp;g(3,2)=-0.7_wp
  original=-1.0_wp;original(3,2)=1.0_wp
  d=original; margin=0.2_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  call check(status==0.and.active,'active linear inequality')
  call check(abs(margin+sum(g*d))<1e-13_wp,'linear feasibility without clipping')
  residual(:,1)=matmul(b(:,:,1),d(:,1)-original(:,1))+matmul(c(:,:,1),d(:,2)-original(:,2))
  residual(:,2)=matmul(a(:,:,2),d(:,1)-original(:,1))+matmul(b(:,:,2),d(:,2)-original(:,2))
  stationarity(:,1)=matmul(transpose(b(:,:,1)),residual(:,1))+matmul(transpose(a(:,:,2)),residual(:,2))
  stationarity(:,2)=matmul(transpose(c(:,:,1)),residual(:,1))+matmul(transpose(b(:,:,2)),residual(:,2))
  mu=sum(stationarity*g)/sum(g*g)
  call check(mu>0.0_wp,'nonnegative active multiplier')
  call check(maxval(abs(stationarity-mu*g))<1e-12_wp,'independent KKT stationarity')
  d=-original
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  call check(status==0.and..not.active.and.all(d==-original),'admissible Newton direction unchanged')
  d=original
  call active_set_direction(a,b,c,scale,g,-margin,d,active,status)
  call check(status/=0.and..not.active.and.all(d==original),'infeasible linearization refused')
  scale(1,1)=0.0_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  call check(status/=0.and..not.active,'zero scale refused')
  scale=1.0_wp;g(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  call check(status/=0.and..not.active,'nonfinite gradient refused')
  a=0.0_wp;b=0.0_wp;c=0.0_wp
  do i=1,3
    b(i,i,:)=1.0_wp
  enddo
  b(1,1,1)=1e-200_wp
  g=0.0_wp;g(1,1)=1.0_wp;d=-1.0_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  call check(status==0.and.active,'large finite dual avoids squared-norm overflow')
  call check(abs(margin+d(1,1))<1e-13_wp,'near-singular representable correction')
  b(1,1,1)=0.0_wp;d=-1.0_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  call check(status/=0.and..not.active,'singular band matrix reports numerical status')
  write(*,'(a,i0,a)') 'PASS ',checks,' active_set checks'
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
end program test_active_set
