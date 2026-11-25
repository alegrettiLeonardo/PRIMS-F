program ss3_local_closures
 !! Standalone local constitutive study (diagnostic, NOT used by the solver).
 !! For a fixed conservative cell state U (alpha, p, rho fixed by U1,U2 and the
 !! EOS), the mixture momentum U3 links the signed velocities:
 !!   u_l = (U3 - U2 u_g)/U1.
 !! Each candidate closure Phi(alpha,rho_l,rho_g,u_l,u_g,theta)=0 is scanned in
 !! u_g over a co-current interval [0,U3/U2] AND its counter-current extension
 !! (u_g>U3/U2 => u_l<0). Sign changes are bracketed and bisected; closure
 !! domain refusals are counted, never skipped silently.
 !! Arguments: fixture folder, gas sound speed m/s, cell index, scan factor.
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_closure_law,only:lei_fechamento_or_ndim_simp
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),cgas,fac,p,rl,rg,alpha,ugmax,theta
 integer(c_int)::st
 integer::which,cell,ids(3),k,current
 character(512)::folder,arg
 character(12)::names(3)=['STRATIFIED  ','BENDIKSEN   ','BHAGWAT     ']
 ids=[1,2,3]
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call get_command_argument(3,arg);read(arg,*)cell;call get_command_argument(4,arg);read(arg,*)fac
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P)
 theta=cfg%thetav(cell)
 do which=0,1
  call probe_state(int(which,c_int),x,sc,meta);call unpack(x,uu)
  p=pressure_fun_uj(uu(cell,1),uu(cell,2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2;alpha=uu(cell,2)/rg
  ugmax=uu(cell,3)/uu(cell,2)
  write(*,'(a,i2,a,i4,a,f8.3,a,f9.6,a,f12.1,a,2f10.4,a,es12.4)')'STATE ',which,' cell',cell,' theta_deg',theta*180/acos(-1._wp), &
    ' alpha',alpha,' p_Pa',p*cfg%omega_P,' rho_l,g',rl*cfg%omega_rho,rg*cfg%omega_rho,' G_mix kg/m2s',uu(cell,3)*cfg%omega_rho*cfg%omega_u
  do k=1,3
   call scan(ids(k),names(k))
  enddo
 enddo
contains
 subroutine scan(id,name)
  integer,intent(in)::id
  character(*),intent(in)::name
  integer,parameter::m=6000
  real(wp)::ug,f,fprev,ugprev,a,b,fa,fb,c,fc,ul
  integer::i,code,nd,nroot,j,cprev
  current=id;nd=0;nroot=0;cprev=1;fprev=0;ugprev=0
  do i=0,m
   ug=fac*ugmax*(real(i,wp)/m)**2
   call phi(ug,f,code)
   if(code/=0)then
    nd=nd+1;cprev=code;cycle
   endif
   if(i>0.and.cprev==0.and.f*fprev<0)then
    a=ugprev;b=ug;fa=fprev;fb=f
    do j=1,200
     c=.5_wp*(a+b);call phi(c,fc,code)
     if(code/=0.or.c==a.or.c==b)exit
     if(fa*fc<=0)then
      b=c;fb=fc
     else
      a=c;fa=fc
     endif
    enddo
    nroot=nroot+1;ul=(uu(cell,3)-uu(cell,2)*c)/uu(cell,1)
    write(*,'(4x,a,a,i2,a,2f12.6,a,a)')name,' root',nroot,'  u_l,u_g (m/s)',ul*cfg%omega_u,c*cfg%omega_u, &
      '  branch ',merge('CO-CURRENT     ','COUNTER-CURRENT',ul>=0)
   endif
   fprev=f;ugprev=ug;cprev=code
  enddo
  write(*,'(4x,a,a,i3,a,i5,a,i5,a,f8.2,a)')name,' roots',nroot,'  domain refusals',nd,' of',m+1, &
    '  scan u_g up to',fac,' x U3/U2'
 end subroutine
 subroutine phi(ug,f,code)
  real(wp),intent(in)::ug
  real(wp),intent(out)::f
  integer,intent(out)::code
  real(wp)::ul
  ul=(uu(cell,3)-uu(cell,2)*ug)/uu(cell,1)
  call lei_fechamento_or_ndim_simp(alpha,rl,rg,ul,ug,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
    cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,f,closure_id=current,status=code)
  if(code==0.and..not.ieee_is_finite(f))code=-1
 end subroutine
end program
