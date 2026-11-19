subroutine probe_cell(x,index,data,status) bind(C)
  use offline_api
  real(c_double),intent(in) :: x(nx)
  integer(c_int),value :: index
  real(c_double),intent(out) :: data(15)
  integer(c_int),intent(out) :: status
  real(wp) :: u(nc,3),v(6),f(3),s(3),a(3,3),ds(3,3)
  integer :: code
  call unpack(x,u)
  call sim%point_state(u(index,:),cfg%thetav(index),v,f,s,a,ds,code,closure_id=3)
  status=code
  data=ieee_value(0.0_wp,ieee_quiet_nan)
  if(code/=0)return
  data(1:6)=v*[1.0_wp,cfg%omega_rho,cfg%omega_rho,cfg%omega_P,cfg%omega_u,cfg%omega_u]
  data(7:9)=u(index,:)*[cfg%omega_rho,cfg%omega_rho,cfg%omega_rho*cfg%omega_u]
  data(10:15)=[cfg%thetav(index),cfg%D_H,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma]
end subroutine probe_cell
