program ss3_closure_map
 !! P1.1 standalone constitutive map at one cell (diagnostic, NOT used by the
 !! solver). Cell state U fixes alpha, p, rho_l, rho_g (EOS); the mixture
 !! momentum links the signed velocities, u_l = (U3 - U2 u_g)/U1, so u_l = 0
 !! at u_g = U3/U2. For each existing closure Phi_k(alpha,rho,u_l,u_g,theta):
 !!  MAP     Phi and dPhi/du_g on both sides of u_l = 0 (and u_g < 0);
 !!  ROOT    every sign change, bisected; derivative (simple vs multiple),
 !!          relative conditioning |u_g dPhi/du_g| / max|Phi|, and the
 !!          eigenvalues of the drift-flux flux Jacobian at the root
 !!          (hyperbolic iff all real);
 !!  TANGENT local minima of |Phi| without sign change (candidate double roots);
 !!  BRANCH  roots while the mixture momentum is scaled, U3 -> kappa U3,
 !!          showing where a branch reaches u_l = 0 and whether it continues.
 !! Domain refusals of a closure are counted and reported, never skipped.
 !! Arguments: fixture folder, gas sound speed m/s, cell index.
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_closure_law,only:lei_fechamento_or_ndim_simp
 use module_jacobian,only:DFduj
 use module_lapack_small,only:dgeev_3x3
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 integer,parameter::m=8000
 real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),cgas,p,rl,rg,alpha,theta,u1,u2,u3,u3ref,kappa
 real(wp)::ug(0:m),phi(0:m),ugmin,ugmax,scale
 integer::cell,k,i,code(0:m),id,ik,nroot
 character(512)::folder,arg
 character(12)::names(3)=['STRATIFIED  ','BENDIKSEN   ','BHAGWAT     ']
 integer(c_int)::st
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call get_command_argument(3,arg);read(arg,*)cell
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P)
 call probe_state(0_c_int,x,sc,meta);call unpack(x,uu)
 theta=cfg%thetav(cell);u1=uu(cell,1);u2=uu(cell,2);u3ref=uu(cell,3)
 p=pressure_fun_uj(u1,u2,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
 rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2;alpha=u2/rg
 write(*,'(a,i4,a,f8.3,a,f10.6,a,f12.1,a,2f10.4,a,es12.4,a,f10.6)')'STATE|cell',cell,'|theta_deg',theta*180/acos(-1._wp), &
   '|alpha',alpha,'|p_Pa',p*cfg%omega_P,'|rho_l,g',rl*cfg%omega_rho,rg*cfg%omega_rho,'|G_mix',u3ref*cfg%omega_rho*cfg%omega_u, &
   '|u_g(u_l=0) m/s',u3ref/u2*cfg%omega_u
 u3=u3ref
 ugmin=-1/cfg%omega_u;ugmax=5/cfg%omega_u
 do id=1,3
  do i=0,m
   ug(i)=ugmin+(ugmax-ugmin)*real(i,wp)/m
   call phi_eval(id,ug(i),phi(i),code(i))
  enddo
  scale=maxval(abs(phi),mask=code==0)
  write(*,'(a,a,a,i5,a,i5,a,es12.4)')'SUMMARY|',trim(names(id)),'|admissible',count(code==0),'|refused',count(code/=0),'|max|Phi|',scale
  do i=0,m,40
   write(*,'(a,a,a,2es14.6,a,i4,a,es14.6,a,es14.6)')'MAP|',trim(names(id)),'|u_l,u_g',ul_of(ug(i))*cfg%omega_u,ug(i)*cfg%omega_u, &
     '|',code(i),'|',phi(i),'|',dphi(id,ug(i))
  enddo
  nroot=0
  do i=1,m
   if(code(i)==0.and.code(i-1)==0)then
    if(phi(i)*phi(i-1)<0)then
     nroot=nroot+1;call root_report(id,ug(i-1),ug(i),phi(i-1),scale,nroot)
    endif
   endif
   if(i<m)then
    if(all(code(i-1:i+1)==0))then
     if(abs(phi(i))<abs(phi(i-1)).and.abs(phi(i))<abs(phi(i+1)).and. &
        phi(i)*phi(i-1)>0.and.phi(i)*phi(i+1)>0.and.abs(phi(i))<1e-3_wp*scale) &
       write(*,'(a,a,a,2es14.6,a,es12.4)')'TANGENT|',trim(names(id)),'|u_l,u_g',ul_of(ug(i))*cfg%omega_u,ug(i)*cfg%omega_u, &
         '|min|Phi|/scale',abs(phi(i))/scale
    endif
   endif
  enddo
  do i=1,m
   if((code(i)==0).neqv.(code(i-1)==0))write(*,'(a,a,a,2es14.6,a,i4,a,i4)')'DOMAIN_EDGE|',trim(names(id)),'|u_l,u_g', &
     ul_of(ug(i))*cfg%omega_u,ug(i)*cfg%omega_u,'|codes',code(i-1),'->',code(i)
  enddo
  write(*,'(a,a,a,i3)')'ROOTS|',trim(names(id)),'|',nroot
 enddo
 ! branch continuation in the mixture momentum
 do ik=0,40
  kappa=.5_wp+real(ik,wp)/40
  u3=kappa*u3ref
  do id=1,3
   do i=0,m
    ug(i)=ugmin+(ugmax-ugmin)*real(i,wp)/m
    call phi_eval(id,ug(i),phi(i),code(i))
   enddo
   do i=1,m
    if(code(i)==0.and.code(i-1)==0)then
     if(phi(i)*phi(i-1)<0)write(*,'(a,f6.3,a,a,a,2es14.6)')'BRANCH|kappa',kappa,'|',trim(names(id)),'|u_l,u_g', &
       ul_of(bisect(id,ug(i-1),ug(i),phi(i-1)))*cfg%omega_u,bisect(id,ug(i-1),ug(i),phi(i-1))*cfg%omega_u
    endif
   enddo
  enddo
 enddo
contains
 real(wp) function ul_of(g)
  real(wp),intent(in)::g
  ul_of=(u3-u2*g)/u1
 end function
 subroutine phi_eval(id,g,f,c)
  integer,intent(in)::id
  real(wp),intent(in)::g
  real(wp),intent(out)::f
  integer,intent(out)::c
  call lei_fechamento_or_ndim_simp(alpha,rl,rg,ul_of(g),g,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
    cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,f,closure_id=id,status=c)
  if(c==0.and..not.ieee_is_finite(f))c=-1
 end subroutine
 real(wp) function dphi(id,g)
  integer,intent(in)::id
  real(wp),intent(in)::g
  real(wp)::fp,fm,h
  integer::cp,cm
  h=1e-6_wp*max(abs(g),1e-3_wp)
  call phi_eval(id,g+h,fp,cp);call phi_eval(id,g-h,fm,cm)
  dphi=-huge(1._wp)
  if(cp==0.and.cm==0)dphi=(fp-fm)/(2*h)
 end function
 real(wp) function bisect(id,a0,b0,fa0)
  integer,intent(in)::id
  real(wp),intent(in)::a0,b0,fa0
  real(wp)::a,b,fa,c,fc
  integer::j,cc
  a=a0;b=b0;fa=fa0;c=.5_wp*(a+b)
  do j=1,200
   c=.5_wp*(a+b);call phi_eval(id,c,fc,cc)
   if(cc/=0.or.c==a.or.c==b)exit
   if(fa*fc<=0)then
    b=c
   else
    a=c;fa=fc
   endif
  enddo
  bisect=c
 end function
 subroutine root_report(id,a,b,fa,scale,n)
  integer,intent(in)::id,n
  real(wp),intent(in)::a,b,fa,scale
  real(wp)::g,d,jac(3,3),wr(3),wi(3),vr(3,3)
  integer::c
  g=bisect(id,a,b,fa);d=dphi(id,g)
  jac=DFduj(alpha,rl,rg,ul_of(g),g,cfg%nCl,cfg%nCg,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma, &
    cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%tol,status=c,closure_id=id)
  wr=0;wi=0
  if(c==0)call dgeev_3x3(jac,wr,wi,vr,c)
  write(*,'(a,a,a,i2,a,2es14.6,a,a,a,es12.4,a,es12.4,a,i4,a,3es12.4,a,3es12.4)')'ROOT|',trim(names(id)),'|',n,'|u_l,u_g', &
    ul_of(g)*cfg%omega_u,g*cfg%omega_u,'|',trim(merge('CO-CURRENT     ','COUNTER-CURRENT',ul_of(g)>=0)),'|dPhi/du_g',d, &
    '|cond',abs(g*d)/scale,'|jac_status',c,'|eig_re m/s',wr*cfg%omega_u,'|eig_im',wi*cfg%omega_u
 end subroutine
end program
