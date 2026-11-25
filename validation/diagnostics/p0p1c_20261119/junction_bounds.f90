program junction_bounds
 !! Diagnostic only: per-phase junction flux bounds at a snapshot. Compares the
 !! MATCHED Roe face flux with the physical fluxes of the two traces (pipe cell,
 !! riser trace mapped with the same p and phase fluxes) and with the HLL flux
 !! reported by offline_api. Nothing is substituted in the solver.
 use offline_api
 use module_lapack_small,only:dgeev_3x3
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),d(48),cgas,s,dm(3,3),du(3),dup(3),wr(3),wi(3),ev(3,3),beta(3),evc(3,3),dp
 integer::piv(3),k,info
 integer(c_int)::st
 integer::w,cc
 character(512)::folder,arg
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P)
 call sim%init(cfg);call sim%prepare_dissipation(old,cc)
 s=cfg%omega_rho*cfg%omega_u
 do w=0,1
  call probe_state(int(w,c_int),x,sc,meta)
  call probe_interface(x,d,st)
  if(st/=0)then
   print *,'INTERFACE_STATUS',w,st;cycle
  endif
  print '(a,i2,a,f8.5,a,f8.5,a,f8.5)','STATE ',w,' alpha_pipe ',d(5),' alpha_riser ',d(11),' alpha_mapped_pipe_closure ',d(17)
  print '(a,3es14.6)','  Gl kg/m2s: left, mapped right, Roe face ',d(29)*s,d(32)*s,d(1)*s
  print '(a,3es14.6)','  Gg kg/m2s: left, mapped right, Roe face ',d(30)*s,d(33)*s,d(2)*s
  print '(a,2es14.6,a,2f9.4)','  HLL face Gl,Gg kg/m2s ',d(37)*s,d(38)*s,'  wave bounds sl,sr m/s ',d(35)*cfg%omega_u,d(36)*cfg%omega_u
  print '(a,l2,a,l2)','  Roe Gg inside [min,max] of trace fluxes:',d(2)>=min(d(30),d(33)).and.d(2)<=max(d(30),d(33)), &
    '   HLL Gg inside:',d(38)>=min(d(30),d(33)).and.d(38)<=max(d(30),d(33))
  ! Wave decomposition of the gas dissipation -0.5*|A|(U_mapped-U_left).
  dm=sim%lagged_dissipation(:,:,ib);du=d(26:28)-d(23:25)
  call dgeev_3x3(dm,wr,wi,ev,cc)
  if(cc/=0.or.any(wi/=0))then
   print *,'  EIGEN_FAILURE',cc;cycle
  endif
  evc=ev;beta=du;call dgesv(3,1,evc,3,piv,beta,3,info)
  print '(a,3es12.4)','  |lambda| of LAGGED |A| (nondim)            ',wr
  print '(a,4es13.5)','  gas dissipation by wave k=1..3, total (kg/m2s)',(-.5_wp*wr(k)*beta(k)*ev(2,k)*s,k=1,3),-.5_wp*dot_product(dm(2,:),du)*s
  dp=d(20)-d(8)
  dup=[(1-d(17))/cfg%nCl**2,d(17)/cfg%nCg**2,(1-d(17))*d(21)/cfg%nCl**2+d(17)*d(22)/cfg%nCg**2]*dp
  print '(a,es12.4,a,es13.5)','  p_mapped-p_left (Pa) ',dp*cfg%omega_P,'  gas dissipation of that pressure jump alone (kg/m2s) ',-.5_wp*dot_product(dm(2,:),dup)*s
 enddo
end program
