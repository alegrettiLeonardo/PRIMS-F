program seal_face_probe
  !! P0.1 seal-limit certificate (DIAGNOSTIC; not used by the solver).
  !! Independent of the solver's t1_face: own wave curves (RK4 integral curves
  !! of the right eigenvectors of A=dF/dU, acoustic families parametrized by p,
  !! the kinematic family by ln(alpha), 64 steps), own inner Newton.  Shared with
  !! the solver: only the closures (sim%point_state) and the EOS.
  !!
  !! Usage: seal_face_probe FIXTURE WHICH [MODE]  (WHICH 0 = accepted U, 1 = refused trial)
  !!   MODE MAP (default): CELL + H4 + H5 maps;  MODE CONV: convergence of g(alpha_R*) in H4
  !!   with the number of RK4 steps (4..256) and with linearized acoustic waves (as the solver).
  !! Output (SI; '|' separated):
  !!  CELL|i|section|theta_deg|U1 U2 U3 (kg/m3, kg/m3, kg/m2/s)|p alpha rho_g rho_l u_g u_l j_g j_l C0 Vgj|F1 F2 F3|lam1 lam2 lam3|status
  !!  H4|alpha_R*|status|p_L* p_M (Pa)|Gl_L Gl_R F3_L F3_R (SI)|Fg_L Fg_R g=Fg_L-Fg_R|alpha_L* p_R*|lam(U_L*) lam(U_R*) lam(U_M)
  !!     configuration H4: pipe 1-wave (left) | stationary layer | riser 2- and 3-waves (right);
  !!     for a fixed alpha_R* the liquid-flux and F3 conditions are solved for (p_L*, p_M);
  !!     g(alpha_R*) is the remaining gas-flux mismatch. A root of g with the characteristic
  !!     count #{lam(U_L*)<0}+#{lam(U_R*)>0} = 3 is a T1 solution in H4.
  !!  H5|alpha_L*|... configuration H5: pipe 1- and 2-waves (left) | layer | riser 3-wave;
  !!     U_L1 = 1-curve(U_L;p_L1), U_L* = 2-curve(U_L1; ln alpha_L*), U_R* = 3-curve(U_R;p_M).
  !!  COUNT lines give the characteristic count at each converged map point.
  use offline_api
  use module_lapack_small,only:dgeev_3x3
  use module_chexal_1992,only:chexal_coefficients_si
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  character(512)::folder,arg
  real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),xl(3),yr(3),thl,thr,vl(6),vr(6),f(3),lam(3),r(3,3)
  real(wp)::a,amin,amax,s(2),out(20),pf,uf,rf,ff3,c0,vgj,jl,jg
  integer::which,i,code,k,n,cid,m
  integer::nst=64
  logical::lin_acoustic=.false.
  real(wp),parameter::alist(4)=[1e-3_wp,1e-5_wp,4.3569e-7_wp,1e-9_wp]
  character(16)::mode='MAP'
  call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)which
  if(command_argument_count()>2)call get_command_argument(3,mode)
  call probe_init(trim(folder),len_trim(folder),1_c_int,code)
  if(code/=0)then
    write(*,'(a,i0)')'probe_init status ',code;stop 1
  endif
  call probe_state(int(which,c_int),x,sc,meta);call unpack(x,uu)
  pf=cfg%omega_P;uf=cfg%omega_u;rf=cfg%omega_rho;ff3=rf*uf**2
  do i=ib-3,ib+4
    cid=cfg%closure_pipe;if(i>ib)cid=cfg%closure_riser
    call eig(uu(i,:),cfg%thetav(i),cid,lam,r,vl,f,code)
    c0=ieee_value(0.0_wp,ieee_quiet_nan);vgj=c0
    jl=(1-vl(1))*vl(5)*uf;jg=vl(1)*vl(6)*uf
    if(cid==5.and.code==0)call chexal_coefficients_si(vl(1),vl(2)*rf,vl(3)*rf,jl,jg,cfg%thetav(i),cfg%D_H,cfg%G, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,c0,vgj,k)
    write(*,'(a,i0,a,a,a,f8.3,a,3es14.6,a,10es14.6,a,3es14.6,a,3es12.4,a,i0)')'CELL|',i,'|',merge('pipe ','riser',i<=ib),'|', &
      cfg%thetav(i)*180/acos(-1.0_wp),'|',uu(i,1:2)*rf,uu(i,3)*rf*uf,'|',vl(4)*pf,vl(1),vl(3)*rf,vl(2)*rf,vl(6)*uf,vl(5)*uf, &
      jg,jl,c0,vgj,'|',f(1:2)*rf*uf,f(3)*ff3,'|',lam*uf,'|',code
  enddo
  xl=uu(ib,:);yr=uu(ib+1,:);thl=cfg%thetav(ib);thr=cfg%thetav(ib+1)
  call eig(xl,thl,cfg%closure_pipe,lam,r,vl,f,code)
  call eig(yr,thr,cfg%closure_riser,lam,r,vr,f,code)
  if(trim(mode)=='CONV')then
    do m=0,1
      lin_acoustic=m==1
      do k=2,8
        nst=2**k
        do i=1,4
          a=alist(i)
          s=[vl(4),vr(4)]
          call solve_map(4,a,s,out,code)
          write(*,'(a,l1,a,i4,a)',advance='no')'CONV|lin_acoustic=',lin_acoustic,'|nst=',nst,'|'
          call emit('H4',a,code,out)
        enddo
      enddo
    enddo
    stop
  endif
  ! H4 map over alpha_R*, from above the riser cell value down to 1e-14
  s=[vl(4),vr(4)]
  amax=min(0.9_wp,4*vr(1));amin=1e-14_wp;n=90
  do k=0,n
    a=exp(log(amax)+(log(amin)-log(amax))*k/real(n,wp))
    call solve_map(4,a,s,out,code)
    call emit('H4',a,code,out)
  enddo
  ! H5 map over alpha_L* (pipe kinematic wave), continuation from the pipe cell value
  s=[vl(4),vr(4)]
  do k=0,n
    a=vl(1)*exp(log(0.999_wp)*0+(log(1e-12_wp/vl(1)))*k/real(n,wp))
    call solve_map(5,a,s,out,code)
    call emit('H5',a,code,out)
  enddo
  s=[vl(4),vr(4)]
  do k=1,40
    a=vl(1)+(0.999_wp-vl(1))*k/40.0_wp
    call solve_map(5,a,s,out,code)
    call emit('H5',a,code,out)
  enddo
contains
  subroutine emit(tag,a,code,o)
    character(*),intent(in)::tag
    real(wp),intent(in)::a,o(20)
    integer,intent(in)::code
    write(*,'(2a,es13.5,a,i3,a,2es16.8,a,4es15.7,a,3es13.5,a,2es13.5,a,9es11.3,a,i0)')tag,'|',a,'|',code,'|',o(1:2)*pf,'|', &
      o(3:4)*rf*uf,o(5:6)*ff3,'|',o(7:8)*rf*uf,(o(7)-o(8))*rf*uf,'|',o(9),o(10)*pf,'|',o(11:19)*uf,'|',nint(o(20))
  end subroutine
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
  subroutine grads(v,gp,ga)
    real(wp),intent(in)::v(6)
    real(wp),intent(out)::gp(3),ga(3)
    real(wp)::den
    den=v(2)*v(1)/cfg%nCg**2+v(3)*(1-v(1))/cfg%nCl**2
    ga=[-v(1)/cfg%nCg**2,(1-v(1))/cfg%nCl**2,0.0_wp]/den
    gp=[v(3),v(2),0.0_wp]/den
  end subroutine
  subroutine field(u,th,cid,fam,du,code)
    !! fam 1,3: dU/dp ; fam 2: dU/d(ln alpha)
    real(wp),intent(in)::u(3),th
    integer,intent(in)::cid,fam
    real(wp),intent(out)::du(3)
    integer,intent(out)::code
    real(wp)::lam(3),r(3,3),v(6),f(3),gp(3),ga(3),dq
    call eig(u,th,cid,lam,r,v,f,code)
    if(code/=0)return
    call grads(v,gp,ga)
    if(fam==2)then
      dq=dot_product(ga,r(:,2))/v(1)
    else
      dq=dot_product(gp,r(:,fam))
    endif
    if(.not.(abs(dq)>0.0_wp))then
      code=-8;return
    endif
    du=r(:,fam)/dq
  end subroutine
  subroutine curve(u0,th,cid,fam,q1,u1,code)
    !! RK4 integral curve of family fam from u0 to p=q1 (fam 1,3) or ln alpha=q1 (fam 2)
    real(wp),intent(in)::u0(3),th,q1
    integer,intent(in)::cid,fam
    real(wp),intent(out)::u1(3)
    integer,intent(out)::code
    real(wp)::v(6),f(3),k1(3),k2(3),k3(3),k4(3),h,q0,lam(3),r(3,3),gp(3),ga(3)
    integer::i
    u1=u0
    call eig(u0,th,cid,lam,r,v,f,code)
    if(code/=0)return
    q0=merge(log(v(1)),v(4),fam==2)
    if(fam/=2.and.lin_acoustic)then
      ! the solver's acoustic waves: straight line along r_fam(u0), parametrized by p
      call grads(v,gp,ga)
      u1=u0+(q1-q0)*r(:,fam)/dot_product(gp,r(:,fam));return
    endif
    h=(q1-q0)/nst
    if(h==0.0_wp)return
    do i=1,nst
      call field(u1,th,cid,fam,k1,code);if(code/=0)return
      call field(u1+h/2*k1,th,cid,fam,k2,code);if(code/=0)return
      call field(u1+h/2*k2,th,cid,fam,k3,code);if(code/=0)return
      call field(u1+h*k3,th,cid,fam,k4,code);if(code/=0)return
      u1=u1+h*(k1+2*k2+2*k3+k4)/6
    enddo
  end subroutine
  subroutine states(cfgn,a,s,uls,urs,um,code)
    integer,intent(in)::cfgn
    real(wp),intent(in)::a,s(2)
    real(wp),intent(out)::uls(3),urs(3),um(3)
    integer,intent(out)::code
    real(wp)::ul1(3)
    if(cfgn==4)then
      call curve(xl,thl,cfg%closure_pipe,1,s(1),uls,code);if(code/=0)return
      call curve(yr,thr,cfg%closure_riser,3,s(2),um,code);if(code/=0)return
      call curve(um,thr,cfg%closure_riser,2,log(a),urs,code)
    else
      call curve(xl,thl,cfg%closure_pipe,1,s(1),ul1,code);if(code/=0)return
      call curve(ul1,thl,cfg%closure_pipe,2,log(a),uls,code);if(code/=0)return
      call curve(yr,thr,cfg%closure_riser,3,s(2),um,code);if(code/=0)return
      urs=um
    endif
  end subroutine
  subroutine res2(cfgn,a,s,rr,fl,fr,uls,urs,um,code)
    integer,intent(in)::cfgn
    real(wp),intent(in)::a,s(2)
    real(wp),intent(out)::rr(2),fl(3),fr(3),uls(3),urs(3),um(3)
    integer,intent(out)::code
    real(wp)::v(6),lam(3),r(3,3)
    call states(cfgn,a,s,uls,urs,um,code);if(code/=0)return
    call eig(uls,thl,cfg%closure_pipe,lam,r,v,fl,code);if(code/=0)return
    call eig(urs,thr,cfg%closure_riser,lam,r,v,fr,code);if(code/=0)return
    rr=[fl(1)-fr(1),fl(3)-fr(3)]
  end subroutine
  subroutine solve_map(cfgn,a,s,o,code)
    !! inner Newton on (p_L*, p_M) for Gl and F3 continuity at fixed alpha (continuation in s)
    integer,intent(in)::cfgn
    real(wp),intent(in)::a
    real(wp),intent(inout)::s(2)
    real(wp),intent(out)::o(20)
    integer,intent(out)::code
    real(wp)::rr(2),rn(2),jm(2,2),det,ds(2),sn(2),fl(3),fr(3),uls(3),urs(3),um(3),h,lam_,sc2(2),v(6),f(3)
    real(wp)::l1(3),l2(3),l3(3),r(3,3)
    integer::it,j,k,cnt
    o=ieee_value(0.0_wp,ieee_quiet_nan)
    call res2(cfgn,a,s,rr,fl,fr,uls,urs,um,code)
    if(code/=0)return
    sc2=[max(abs(fl(1)),1e-8_wp),max(abs(fl(3)),1e-8_wp)]
    do it=1,60
      if(maxval(abs(rr)/sc2)<1e-13_wp)exit
      do j=1,2
        sn=s;h=1e-7_wp*max(abs(s(j)),1.0_wp);sn(j)=s(j)+h
        call res2(cfgn,a,sn,rn,fl,fr,uls,urs,um,code)
        if(code/=0)return
        jm(:,j)=(rn-rr)/h
      enddo
      det=jm(1,1)*jm(2,2)-jm(1,2)*jm(2,1)
      if(.not.(abs(det)>0.0_wp))then
        code=-20;return
      endif
      ds=-[jm(2,2)*rr(1)-jm(1,2)*rr(2),-jm(2,1)*rr(1)+jm(1,1)*rr(2)]/det
      lam_=1
      do k=1,40
        sn=s+lam_*ds
        call res2(cfgn,a,sn,rn,fl,fr,uls,urs,um,code)
        if(code==0)then
          if(maxval(abs(rn)/sc2)<maxval(abs(rr)/sc2))exit
        endif
        lam_=lam_/2
      enddo
      if(code/=0.or.k>40)then
        code=-21;return
      endif
      s=sn;rr=rn
    enddo
    call res2(cfgn,a,s,rr,fl,fr,uls,urs,um,code)
    if(code/=0)return
    if(maxval(abs(rr)/sc2)>=1e-12_wp)then
      code=-22
    endif
    call eig(uls,thl,cfg%closure_pipe,l1,r,v,f,j);o(9)=v(1)
    call eig(urs,thr,cfg%closure_riser,l2,r,v,f,k);o(10)=v(4)
    call eig(um,thr,cfg%closure_riser,l3,r,v,f,k)
    cnt=count(l1<0)+count(l2>0)
    o(1:2)=s;o(3:4)=[fl(1),fr(1)];o(5:6)=[fl(3),fr(3)];o(7:8)=[fl(2),fr(2)]
    o(11:13)=l1;o(14:16)=l2;o(17:19)=l3;o(20)=cnt
  end subroutine
end program seal_face_probe
