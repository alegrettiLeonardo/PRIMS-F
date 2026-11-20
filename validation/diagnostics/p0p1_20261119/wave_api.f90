module wave_api
 use offline_api
contains
 subroutine wave_setup(folder,nchar,base,amat,primitives,angle,scales,status)bind(C)
  character(c_char),intent(in)::folder(*)
  integer(c_int),value::nchar
  real(c_double),intent(out)::base(3),amat(3,3),primitives(6),angle,scales(3)
  integer(c_int),intent(out)::status
  real(wp)::flux(3),source(3),ds(3,3)
  call probe_init(folder,nchar,1_c_int,status)
  if(status/=0)return
  base=old(80,:);angle=cfg%thetav(80);scales=[cfg%omega_u,cfg%omega_rho,cfg%omega_P]
  call sim%point_state(base,angle,primitives,flux,source,amat,ds,status,.true.,3)
 end subroutine
 subroutine wave_admissible(states,n,status)bind(C)
  integer(c_int),value::n
  real(c_double),intent(in)::states(n,3)
  integer(c_int),intent(out)::status(n)
  real(wp)::v(6),flux(3),source(3),aa(3,3),ds(3,3)
  integer::i
  do i=1,n
   call sim%point_state(states(i,:),cfg%thetav(80),v,flux,source,aa,ds,status(i),.false.,3)
  enddo
 end subroutine
end module
