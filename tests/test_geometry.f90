program test_geometry
  use mod_kinds,only:wp
  use module_catenary,only:solve_A,altura_catenaria,comprimento_arco,angulo_inclinacao,fungeo,Dfungeo
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  real(wp)::a,l,theta,h,d,v(2),ar
  integer::code,k,count
  count=0
  a=solve_A(9.886_wp,6.435_wp,status=code)
  call check(code==0,'catenary root status')
  call check(abs(a-3.018181937174462_wp)<1e-12_wp,'catenary parameter')
  l=comprimento_arco(6.435_wp,a)
  call check(abs(l-12.54625399319707_wp)<1e-11_wp,'arc length')
  call check(abs(altura_catenaria(6.435_wp,a)-9.886_wp)<1e-12_wp,'endpoint elevation')
  theta=angulo_inclinacao(l,a)
  call check(abs(theta*180/acos(-1.0_wp)-76.47369465598749_wp)<1e-11_wp,'top inclination')
  v=[a,l];h=1e-5_wp
  do k=0,10
    d=(fungeo(real(k,wp)/10+h,v)-fungeo(real(k,wp)/10-h,v))/(2*h)
    call check(abs(d-Dfungeo(real(k,wp)/10,v))<3e-9_wp,'inclination derivative')
  enddo
  do k=-6,6
    ar=solve_A(10.0_wp**k,1.0_wp,status=code)
    call check(code==0,'wide aspect ratio status')
    ! sinh form avoids cancellation for a shallow catenary.
    call check(abs(2*ar*sinh(0.5_wp/ar)**2/10.0_wp**k-1)<2e-13_wp,'wide aspect endpoint')
  enddo
  ar=solve_A(-1.0_wp,1.0_wp,status=code)
  call check(code/=0.and..not.ieee_is_finite(ar),'invalid geometry rejected')
  write(*,'(a,i0,a)') 'PASS geometry: ',count,' assertions'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then
      write(*,'(a)')'FAIL '//msg
      error stop 1
    endif
    count=count+1
  end subroutine
end program
