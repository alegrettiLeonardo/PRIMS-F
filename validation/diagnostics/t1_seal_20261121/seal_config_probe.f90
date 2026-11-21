program seal_config_probe
  !! P0.1/P0.2 configuration certificate for the T1 junction near the seal
  !! (DIAGNOSTIC, independent of the solver's t1_face; shares only the closures).
  !!
  !! Kinematic (2-family) integral curve parametrized by z = sqrt(alpha):
  !!   dU/dz = 2 z r2/(grad(alpha).r2)      (finite limit of r2/(grad alpha.r2) as
  !! alpha->0+ verified by the H4 map: lambda2 -> Vgj(0)), so the endpoint alpha=0
  !! is reached WITHOUT evaluating the closure at alpha=0 (field = 0 at z=0).
  !! Acoustic families: RK4 integral curves in p (8 steps).
  !!
  !!  G0   : H4 restricted to alpha_R*=0 (riser trace gas-free): solve Gl,F3 for
  !!         (p_L*, p_M); g0 = Fg(U_L*) (the pipe-side gas flux at the face).
  !!         g0 > 0 : the pipe pushes gas into a gas-free trace -> an H4 root with
  !!                  alpha_R* > 0 is sought (bisection in z, continuation);
  !!         g0 <= 0: no H4 root near alpha_R*=0+ ; LIQUID-BRIDGE configuration H7.
  !!  H4   : root of g(z) = Fg_L* - Fg_R* (inner Gl,F3 solves), bisection in z.
  !!  H7   : pipe 1-wave U_L1 = 1-curve(U_L;p_L1); pipe 2-front from U_L1 to a pure
  !!         liquid state U* (Rankine-Hugoniot, speed s = Fg(U_L1)/U2(U_L1) = u_g(U_L1));
  !!         U* = 2-curve endpoint alpha=0 of U_M = 3-curve(U_R;p_M). Unknowns
  !!         (p_L1,p_M); equations s[U1]=[F1], s[U3]=[F3]. Face flux F(U*) = (U3*,0,F3*).
  !!         At s=0 the H7 equations coincide with H4 at alpha_R*=0 with g0=0:
  !!         the face flux is continuous through the switch.
  !! Usage: seal_config_probe FIXTURE WHICH [NZ]
  use offline_api
  use module_lapack_small,only:dgeev_3x3
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  character(512)::folder,arg
  real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),xl(3),yr(3),thl,thr,vl(6),vr(6),f(3),lam(3),r(3,3)
  real(wp)::s(2),g0,o(24),pf,uf,rf,ff3,zlo,zhi,glo,ghi,zm,gm,ust(3,4)
  integer::which,code,k,nz
  call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)which
  nz=16
  if(command_argument_count()>2)then
    call get_command_argument(3,arg);read(arg,*)nz
  endif
  call probe_init(trim(folder),len_trim(folder),1_c_int,code)
  if(code/=0)then
    write(*,'(a,i0)')'probe_init status ',code;stop 1
  endif
  call probe_state(int(which,c_int),x,sc,meta);call unpack(x,uu)
  pf=cfg%omega_P;uf=cfg%omega_u;rf=cfg%omega_rho;ff3=rf*uf**2
  xl=uu(ib,:);yr=uu(ib+1,:);thl=cfg%thetav(ib);thr=cfg%thetav(ib+1)
  call eig(xl,thl,cfg%closure_pipe,lam,r,vl,f,code)
  call eig(yr,thr,cfg%closure_riser,lam,r,vr,f,code)
  write(*,'(a,es15.7,a,i0,a,i0)')'STATE|t=',meta(1),'|which=',which,'|nz=',nz
  ! G0
  s=[vl(4),vr(4)]
  call solve2(4,0.0_wp,s,o,ust,code)
  g0=o(7)
  write(*,'(a,i4,a,es14.6,a,2es16.8,a,3es15.7,a,es12.4)')'G0|',code,'|g0=Fg(U_L*) kg/m2/s ',g0*rf*uf,'|p_L* p_M ',s*pf, &
    '|Gl F3_L F3_R ',o(3)*rf*uf,o(5:6)*ff3,'|alpha_L* ',o(9)
  if(code/=0)stop
  if(g0>0)then
    ! H4 root: g(z)=Fg_L*-Fg_R*; g(0)=g0>0; find zhi with g<0
    zlo=0;glo=g0;zhi=sqrt(vr(1));s=[vl(4),vr(4)]
    do k=1,60
      call solve2(4,zhi,s,o,ust,code)
      if(code/=0)exit
      ghi=o(7)-o(8)
      if(ghi<0)exit
      zlo=zhi;glo=ghi;zhi=2*zhi
    enddo
    if(code/=0.or.ghi>=0)then
      write(*,'(a,i0)')'H4|NO_BRACKET|',code;stop
    endif
    do k=1,80
      zm=0.5_wp*(zlo+zhi)
      call solve2(4,zm,s,o,ust,code)
      if(code/=0)exit
      gm=o(7)-o(8)
      if(gm>0)then
        zlo=zm
      else
        zhi=zm
      endif
      if(zhi-zlo<=1e-15_wp*zhi)exit
    enddo
    write(*,'(a,i4,a,es14.6,a,3es15.7,a,2es16.8,a,9es11.3,a,i0)')'H4|',code,'|alpha_R* ',zm**2,'|Gl Gg F3 ',o(3)*rf*uf, &
      o(8)*rf*uf,o(6)*ff3,'|p_L* p_R* ',o(1)*pf,o(10)*pf,'|lam ',o(11:19)*uf,'|count ',nint(o(20))
  else
    s=[vl(4),vr(4)]
    call solve2(7,0.0_wp,s,o,ust,code)
    write(*,'(a,i4,a,es14.6,a,3es15.7,a,2es16.8,a,9es11.3,a,es12.4,a,es12.4)')'H7|',code,'|s m/s ',o(21)*uf,'|Gl 0 F3 ', &
      o(3)*rf*uf,0.0_wp,o(6)*ff3,'|p_L1 p* ',o(1)*pf,o(10)*pf,'|lam(U_L1) lam(U*) lam(U_M) ',o(11:19)*uf, &
      '|alpha_L1 ',o(9),'|u_g(U_L1) m/s ',o(22)*uf
  endif
contains
  subroutine eig(u,th,cid,lam,r,v,f,code)
    real(wp),intent(in)::u(3),th
    integer,intent(in)::cid
    real(wp),intent(out)::lam(3),r(3,3),v(6),f(3)
    integer,intent(out)::code
    real(wp)::q(3),a1(3,3),d1(3,3),wr(3),wi(3),vr(3,3)
    integer::idx(3),i,j,t
    lam=ieee_value(0.0_wp,ieee_quiet_nan);r=0
    call sim%point_state(u,th,v,f,q,a1,d1,code,need_jac=.true.,closure_id=cid)
    if(code/=0)return
    call dgeev_3x3(a1,wr,wi,vr,code)
    if(code/=0)return
    if(any(wi/=0.0_wp))then
      code=-7;return
    endif
    idx=[1,2,3]
    do i=1,2
      do j=i+1,3
        if(wr(idx(j))<wr(idx(i)))then
          t=idx(i);idx(i)=idx(j);idx(j)=t
        endif
      enddo
    enddo
    lam=wr(idx);r=vr(:,idx)
  end subroutine
  subroutine flux0(u,th,cid,f,v,code)
    !! flux only (valid at alpha=0: storage convention ug=0 is multiplied by U2=0)
    real(wp),intent(in)::u(3),th
    integer,intent(in)::cid
    real(wp),intent(out)::f(3),v(6)
    integer,intent(out)::code
    real(wp)::q(3),a1(3,3),d1(3,3)
    call sim%point_state(u,th,v,f,q,a1,d1,code,need_jac=.false.,closure_id=cid)
  end subroutine
  subroutine grads(v,gp,ga)
    real(wp),intent(in)::v(6)
    real(wp),intent(out)::gp(3),ga(3)
    real(wp)::den
    den=v(2)*v(1)/cfg%nCg**2+v(3)*(1-v(1))/cfg%nCl**2
    ga=[-v(1)/cfg%nCg**2,(1-v(1))/cfg%nCl**2,0.0_wp]/den
    gp=[v(3),v(2),0.0_wp]/den
  end subroutine
  subroutine fieldp(u,th,cid,fam,du,code)
    real(wp),intent(in)::u(3),th
    integer,intent(in)::cid,fam
    real(wp),intent(out)::du(3)
    integer,intent(out)::code
    real(wp)::lam(3),r(3,3),v(6),f(3),gp(3),ga(3)
    call eig(u,th,cid,lam,r,v,f,code)
    if(code/=0)return
    call grads(v,gp,ga)
    du=r(:,fam)/dot_product(gp,r(:,fam))
  end subroutine
  subroutine fieldz(u,z,th,cid,du,code)
    !! dU/dz along the 2-curve, z=sqrt(alpha); exactly 0 at z=0 (finite dU/dalpha)
    real(wp),intent(in)::u(3),z,th
    integer,intent(in)::cid
    real(wp),intent(out)::du(3)
    integer,intent(out)::code
    real(wp)::lam(3),r(3,3),v(6),f(3),gp(3),ga(3)
    du=0;code=0
    if(z==0.0_wp)return
    call eig(u,th,cid,lam,r,v,f,code)
    if(code/=0)return
    call grads(v,gp,ga)
    du=2*z*r(:,2)/dot_product(ga,r(:,2))
  end subroutine
  subroutine curvep(u0,th,cid,fam,p1,u1,code)
    real(wp),intent(in)::u0(3),th,p1
    integer,intent(in)::cid,fam
    real(wp),intent(out)::u1(3)
    integer,intent(out)::code
    real(wp)::v(6),f(3),k1(3),k2(3),k3(3),k4(3),h
    integer::i
    u1=u0
    call flux0(u0,th,cid,f,v,code);if(code/=0)return
    h=(p1-v(4))/8
    do i=1,8
      call fieldp(u1,th,cid,fam,k1,code);if(code/=0)return
      call fieldp(u1+h/2*k1,th,cid,fam,k2,code);if(code/=0)return
      call fieldp(u1+h/2*k2,th,cid,fam,k3,code);if(code/=0)return
      call fieldp(u1+h*k3,th,cid,fam,k4,code);if(code/=0)return
      u1=u1+h*(k1+2*k2+2*k3+k4)/6
    enddo
  end subroutine
  function onz(w,z) result(u)
    !! state on the 2-curve from (U1,U3) and alpha=z**2: U2 from the EOS
    !! (rho_l = U1/(1-alpha) gives p; U2 = alpha p/Cg^2). Exact identity, no projection.
    real(wp),intent(in)::w(2),z
    real(wp)::u(3),pp,a
    a=z*z
    pp=cfg%nP_l0+cfg%nCl**2*(w(1)/(1-a)-cfg%nrho_l0)
    u=[w(1),a*pp/cfg%nCg**2,w(2)]
  end function
  subroutine fz(w,z,th,cid,dw,code)
    real(wp),intent(in)::w(2),z,th
    integer,intent(in)::cid
    real(wp),intent(out)::dw(2)
    integer,intent(out)::code
    real(wp)::du(3)
    call fieldz(onz(w,z),z,th,cid,du,code)
    dw=du([1,3])
  end subroutine
  subroutine curvez(u0,th,cid,z1,u1,code)
    !! RK4 in z=sqrt(alpha) on (U1,U3); every stage state has U2 = z_stage^2 p/Cg^2 > 0
    real(wp),intent(in)::u0(3),th,z1
    integer,intent(in)::cid
    real(wp),intent(out)::u1(3)
    integer,intent(out)::code
    real(wp)::v(6),f(3),k1(2),k2(2),k3(2),k4(2),h,z,w(2)
    integer::i
    u1=u0
    call flux0(u0,th,cid,f,v,code);if(code/=0)return
    z=sqrt(v(1));h=(z1-z)/nz;w=u0([1,3])
    do i=1,nz
      call fz(w,z,th,cid,k1,code);if(code/=0)return
      call fz(w+h/2*k1,z+h/2,th,cid,k2,code);if(code/=0)return
      call fz(w+h/2*k2,z+h/2,th,cid,k3,code);if(code/=0)return
      call fz(w+h*k3,z+h,th,cid,k4,code);if(code/=0)return
      w=w+h*(k1+2*k2+2*k3+k4)/6;z=z+h
    enddo
    u1=onz(w,z1)
  end subroutine
  subroutine res2(cfgn,z,s,rr,o,ust,code)
    integer,intent(in)::cfgn
    real(wp),intent(in)::z,s(2)
    real(wp),intent(out)::rr(2),o(24),ust(3,4)
    integer,intent(out)::code
    real(wp)::ul1(3),um(3),urs(3),fl(3),fr(3),v(6),vv(6),sp
    o=0
    call curvep(xl,thl,cfg%closure_pipe,1,s(1),ul1,code);if(code/=0)return
    call curvep(yr,thr,cfg%closure_riser,3,s(2),um,code);if(code/=0)return
    call curvez(um,thr,cfg%closure_riser,z,urs,code);if(code/=0)return
    call flux0(ul1,thl,cfg%closure_pipe,fl,v,code);if(code/=0)return
    call flux0(urs,thr,cfg%closure_riser,fr,vv,code);if(code/=0)return
    if(cfgn==4)then
      rr=[fl(1)-fr(1),fl(3)-fr(3)]
      sp=0
    else
      sp=fl(2)/ul1(2)                       ! s = u_g(U_L1), gas Rankine-Hugoniot with U2*=0
      rr=[sp*(urs(1)-ul1(1))-(fr(1)-fl(1)),sp*(urs(3)-ul1(3))-(fr(3)-fl(3))]
    endif
    o(1:2)=s;o(3:4)=[fl(1),fr(1)];o(5:6)=[fl(3),fr(3)];o(7:8)=[fl(2),fr(2)];o(9)=v(1);o(10)=vv(4);o(21)=sp;o(22)=v(6)
    ust(:,1)=ul1;ust(:,2)=urs;ust(:,3)=um
  end subroutine
  subroutine solve2(cfgn,z,s,o,ust,code)
    integer,intent(in)::cfgn
    real(wp),intent(in)::z
    real(wp),intent(inout)::s(2)
    real(wp),intent(out)::o(24),ust(3,4)
    integer,intent(out)::code
    real(wp)::rr(2),rn(2),jm(2,2),det,ds(2),sn(2),h,lm,sc2(2),on(24),un(3,4),l1(3),l2(3),l3(3),rv(3,3),v(6),f(3)
    integer::it,j,k,c2
    call res2(cfgn,z,s,rr,o,ust,code);if(code/=0)return
    sc2=[max(abs(o(3)),1e-8_wp),max(abs(o(5)),1e-8_wp)]
    do it=1,60
      if(maxval(abs(rr)/sc2)<1e-14_wp)exit
      do j=1,2
        sn=s;h=1e-7_wp*max(abs(s(j)),1.0_wp);sn(j)=s(j)+h
        call res2(cfgn,z,sn,rn,on,un,code);if(code/=0)return
        jm(:,j)=(rn-rr)/h
      enddo
      det=jm(1,1)*jm(2,2)-jm(1,2)*jm(2,1)
      if(.not.(abs(det)>0.0_wp))then
        code=-20;return
      endif
      ds=-[jm(2,2)*rr(1)-jm(1,2)*rr(2),-jm(2,1)*rr(1)+jm(1,1)*rr(2)]/det
      lm=1
      do k=1,40
        sn=s+lm*ds
        call res2(cfgn,z,sn,rn,on,un,code)
        if(code==0)then
          if(maxval(abs(rn)/sc2)<maxval(abs(rr)/sc2))exit
        endif
        lm=lm/2
      enddo
      if(k>40)then
        if(maxval(abs(rr)/sc2)<1e-12_wp)exit   ! round-off floor of the closures
        code=-21;return
      endif
      s=sn;rr=rn;o=on;ust=un
    enddo
    if(maxval(abs(rr)/sc2)>=1e-12_wp)then
      code=-22;return
    endif
    call eig(ust(:,1),thl,cfg%closure_pipe,l1,rv,v,f,c2);o(11:13)=l1
    o(14:16)=ieee_value(0.0_wp,ieee_quiet_nan)
    if(ust(2,2)>0)then
      call eig(ust(:,2),thr,cfg%closure_riser,l2,rv,v,f,c2);o(14:16)=l2
    endif
    call eig(ust(:,3),thr,cfg%closure_riser,l3,rv,v,f,c2);o(17:19)=l3
    o(20)=count(l1<0)+count(o(14:16)>0)
  end subroutine
end program seal_config_probe
