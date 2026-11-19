program test_active_set
  use mod_kinds, only: wp
  use module_active_set, only: active_set_direction
  use module_time_step, only: time_integrator_t
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(time_integrator_t) :: solver
  real(wp) :: a(3,3,2),b(3,3,2),c(3,3,2),scale(3,2),g(3,2),d(3,2),original(3,2)
  real(wp) :: residual(3,2),stationarity(3,2),mu,margin,dual_l1,required_residual,linf_norm,scaled_reference(3,2)
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
  a=0.0_wp;b=0.0_wp;c=0.0_wp;scale=1.0_wp
  do i=1,3
    b(i,i,:)=1.0_wp
  enddo
  b(1,1,1)=2.0_wp;b(2,2,1)=4.0_wp;b(3,3,1)=0.5_wp
  g=0.0_wp;g(:,1)=[1.0_wp,-1.0_wp,0.25_wp]
  original=0.0_wp;original(:,1)=[-1.0_wp,1.0_wp,-1.0_wp]
  margin=0.2_wp;d=original
  dual_l1=0.5_wp+0.25_wp+0.5_wp
  required_residual=(-margin-sum(g*original))/dual_l1
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status==0.and.active,'minimax direction available')
  residual(:,1)=matmul(b(:,:,1),d(:,1)-original(:,1))
  residual(:,2)=matmul(b(:,:,2),d(:,2)-original(:,2))
  linf_norm=maxval(abs(residual))
  call check(abs(linf_norm-required_residual)<1e-13_wp,'independent dual L1 lower bound attained')
  call check(abs(margin+sum(g*d))<1e-13_wp,'minimax linear margin')
  call check(all(d(:,2)==0.0_wp),'zero dual support remains unchanged')
  scaled_reference=d
  d=original
  call active_set_direction(a,b,c,scale,g,margin,d,active,status)
  residual(:,1)=matmul(b(:,:,1),d(:,1)-original(:,1))
  call check(maxval(abs(residual))>linf_norm,'minimax differs from L2 on unequal dual entries')
  d=original*1e100_wp;scale=1e100_wp;g=g/1e100_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status==0.and.active,'minimax accepts consistent units rescaling')
  call check(maxval(abs(d/1e100_wp-scaled_reference))<1e-13_wp,'minimax unit scaling invariance')
  scale=1.0_wp;g=0.0_wp;g(1,1)=1.0_wp;d=-1.0_wp
  b(1,1,1)=1e-200_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status==0.and.active.and.abs(d(1,1)+margin)<1e-13_wp,'minimax finite large dual')
  b(1,1,1)=0.0_wp;d=-1.0_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status/=0.and..not.active,'minimax singular solve refused')
  b(1,1,1)=2.0_wp;g=0.0_wp;g(1,1)=1.0_wp;d=-1.0_wp;scale=1.0_wp
  scale(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status/=0.and..not.active,'NaN scale returns status without trapping')
  scale=1.0_wp;d(1,1)=ieee_value(0.0_wp,ieee_quiet_nan)
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status/=0.and..not.active,'NaN direction returns status without trapping')
  d=-1.0_wp
  call active_set_direction(a,b,c,scale,g,ieee_value(0.0_wp,ieee_quiet_nan),d,active,status,minimax=.true.)
  call check(status/=0.and..not.active,'NaN margin returns status without trapping')
  scale=1e308_wp;g=0.0_wp;g(1,1)=1e-308_wp;d=-1e307_wp
  call active_set_direction(a,b,c,scale,g,0.0_wp,d,active,status,minimax=.true.)
  call check(status==0.and.active,'scaled coefficients avoid intermediate numerator overflow')
  call check(abs(d(1,1)/1e307_wp)<1e-13_wp,'extreme finite scaling retains feasible correction')
  scale=1e308_wp;g=1e308_wp;d=-1.0_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status/=0.and..not.active,'unrepresentable gradient scale refused before multiplication')
  scale=1e-308_wp;g=1.0_wp;d=-1e308_wp
  call active_set_direction(a,b,c,scale,g,margin,d,active,status,minimax=.true.)
  call check(status/=0.and..not.active,'unrepresentable scaled direction refused before division')
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
