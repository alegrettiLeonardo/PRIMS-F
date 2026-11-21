program junction_riemann_ref
 !! P0.2 INDEPENDENT REFERENCE for the junction transmission problem (diagnostic;
 !! not used by the solver). It does not reuse the Roe matrix, the E1 splitting
 !! or the junction mapping: only the physical flux F(U), its Jacobian A=dF/dU
 !! and the closures of each section (sim%point_state, same as the solver).
 !!
 !! Hypotheses (declared, not proven):
 !!  H1  Riemann limit: bounded sources (gravity, friction) do not act at t->0+
 !!      on the self-similar interface states; they are switched off here.
 !!  H2  The junction is a zero-length STATIONARY layer in which the closure
 !!      changes from the pipe law (theta_L) to the riser law (theta_R).
 !!      Transmission conditions (two variants):
 !!        T1 conservative layer:  Gl, Gg and F3 = sum G_k u_k + p continuous;
 !!        T2 PRIMS convention:    Gl, Gg and p continuous (point force).
 !!  H3  Waves are integral curves of the right eigenvectors (exact for
 !!      rarefactions, second order in strength for shocks/contacts). Valid only
 !!      for the weak waves present here; strong fronts are out of scope.
 !!  H4  Non-resonant wave configuration: lambda1<0 on the left, 0<lambda2<lambda3
 !!      on the right; checked at the solution, reported if violated.
 !! Unknowns s=(p_L*, ln alpha_R*, p_M): U_L* = 1-curve(U_L; p=p_L*);
 !! U_M = 3-curve(U_R; p=p_M); U_R* = 2-curve(U_M; ln alpha = s2).
 !! For identical closure/angle on both sides (partition, uniform) T1 reduces
 !! to the standard Riemann problem; there the reference must reproduce
 !! F(U,U)=F(U) and agree with Roe to second order in the jump.
 !! Riser alpha*->0: if Newton drives ln alpha_R* below ln(1e-14) the solution
 !! needs the closed limit alpha_R*=0 (not admitted by the Bhagwat inversion);
 !! this is reported as REF_LIMIT_ALPHA0 with the limiting gas flux.
 !! Arguments: fixture folder, c_gas m/s, mode (JUNCTION | TESTS).
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_boundary_trace,only:steady_junction_trace
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),cgas,pip,xl(3),yr(3),thl,thr,fm(3),fp(3),tr(6),dmat(3,3)
 real(wp)::ref1(12),ref2(12),eps,ul_(3),ur_(3),r(3,3),lam(3),fref(3),fa(3),fb(3)
 integer(c_int)::st
 integer::which,k,ic,code,it,istat1,istat2,cl,cr
 character(512)::folder,arg,mode
 logical::dbg=.false.
 character(16)::trc(2)=['CENTER          ','STEADY          ']
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 mode='JUNCTION'
 if(command_argument_count()>2)call get_command_argument(3,mode)
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
 call geometry()
 call sim%init(cfg)
 pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
 if(trim(mode)=='TESTS')then
  call probe_state(0_c_int,x,sc,meta);call unpack(x,uu)
  ! (a) uniform state, pipe closure, both sides identical: F(U,U)=F(U)
  call solve(uu(ib,:),uu(ib,:),cfg%thetav(ib),cfg%thetav(ib),cfg%closure_pipe,cfg%closure_pipe,1,ref1,istat1)
  call flux_of(uu(ib,:),cfg%thetav(ib),cfg%closure_pipe,fa,code)
  write(*,'(a,i4,a,3es14.6,a,3es12.4)')'TEST|UNIFORM_PIPE|',istat1,'|Fref',ref1(1:3),'|Fref-F(U) rel', &
    (ref1(1:3)-fa)/max(abs(fa),1e-12_wp)
  call solve(uu(ib+1,:),uu(ib+1,:),cfg%thetav(ib+1),cfg%thetav(ib+1),cfg%closure_riser,cfg%closure_riser,1,ref1,istat1)
  call flux_of(uu(ib+1,:),cfg%thetav(ib+1),cfg%closure_riser,fa,code)
  write(*,'(a,i4,a,3es14.6,a,3es12.4)')'TEST|UNIFORM_RISER|',istat1,'|Fref',ref1(1:3),'|Fref-F(U) rel', &
    (ref1(1:3)-fa)/max(abs(fa),1e-12_wp)
  dbg=.true.
  ! (b) artificial partition inside the pipe: real neighbouring pipe states,
  !     same closure and angle -> standard Riemann problem; compare with Roe.
  call solve(uu(ib-1,:),uu(ib,:),cfg%thetav(ib-1),cfg%thetav(ib),cfg%closure_pipe,cfg%closure_pipe,1,ref1,istat1)
  call sim%face_pair(uu(ib-1,:),uu(ib,:),cfg%thetav(ib-1),cfg%thetav(ib),fm,fp,code,left_cell=ib-1)
  write(*,'(a,i4,a,i4,a,3es14.6,a,3es14.6,a,3es12.4)')'TEST|PARTITION_PIPE|',istat1,'|',code,'|Fref',ref1(1:3),'|Froe',fm, &
    '|rel diff',(ref1(1:3)-fm)/max(abs(fm),1e-12_wp)
  dbg=.false.
  ! (c) small waves: U_R = U_L + eps r_k for each family, both pipe and riser.
  do cl=1,2
   if(cl==1)then
    ul_=uu(ib,:);thl=cfg%thetav(ib);cr=cfg%closure_pipe
   else
    ul_=uu(ib+2,:);thl=cfg%thetav(ib+2);cr=cfg%closure_riser
   endif
   call eig(ul_,thl,cr,lam,r,code)
   call flux_of(ul_,thl,cr,fa,code)
   do k=1,3
    do it=1,3
     eps=10._wp**(-1-it)
     ur_=ul_+eps*minval(abs(ul_)/max(abs(r(:,k)),1e-300_wp))*r(:,k)
     call flux_of(ur_,thl,cr,fb,code)
     call solve(ul_,ur_,thl,thl,cr,cr,1,ref1,istat1)
     fref=fa
     if(lam(k)<0)fref=fb
     write(*,'(a,i2,a,i2,a,es9.1,a,es11.3,a,i4,a,3es12.4)')'TEST|SMALL_WAVE|closure',cr,'|family',k,'|eps',eps,'|lambda m/s', &
       lam(k)*cfg%omega_u,'|',istat1,'|(Fref-Fupwind)/|F(UR)-F(UL)|',(ref1(1:3)-fref)/max(maxval(abs(fb-fa)),1e-300_wp)
    enddo
   enddo
  enddo
  ! (d) momentum jump of the junction operator for physically equivalent
  !     sides: the junction face with the right cell replaced by a pipe-closure
  !     state at the same angle is NOT a junction (match=.false.); the point
  !     force must vanish. Reported for the actual junction below instead.
  stop
 endif
 dbg=trim(mode)=='JUNCTION_DEBUG'
 do which=0,1
  call probe_state(int(which,c_int),x,sc,meta);call unpack(x,uu)
  do k=1,2
   xl=uu(ib,:);yr=uu(ib+1,:);thl=cfg%thetav(ib);thr=cfg%thetav(ib+1);ic=0
   if(k==2)then
    call steady_junction_trace(uu(ib,:),cfg%thetav(ib),cfg,.true.,xl,thl,ic)
    if(ic==0)call steady_junction_trace(uu(ib+1,:),cfg%thetav(ib+1),cfg,.false.,yr,thr,ic)
   endif
   if(ic/=0)then
    write(*,'(a,i1,a,a,a,i6)')'REF|',which,'|',trim(trc(k)),'|TRACE_FAIL|',ic;cycle
   endif
   ! The solver maps the riser trace with angle 0 at the face for the riser
   ! inversion; the reference keeps each trace with its own angle.
   call solve(xl,yr,thl,thr,cfg%closure_pipe,cfg%closure_riser,1,ref1,istat1)
   call solve(xl,yr,thl,thr,cfg%closure_pipe,cfg%closure_riser,2,ref2,istat2)
   write(*,'(a,i1,a,a,a,i4,a,12es14.6)')'REF|',which,'|',trim(trc(k)),'|T1|',istat1,'|',out_si(ref1)
   write(*,'(a,i1,a,a,a,i4,a,12es14.6)')'REF|',which,'|',trim(trc(k)),'|T2|',istat2,'|',out_si(ref2)
   ! Operators on the same cell states (their own trace handling).
   cfg%junction_trace=trc(k)
   cfg%junction_flux='ROE';call sim%init(cfg)
   call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code,junction_trace=tr,left_cell=ib, &
     used_dissipation=dmat)
   write(*,'(a,i1,a,a,a,i6,a,6es14.6)')'OP|',which,'|',trim(trc(k)),'|ROE|',code,'|',tr(1)*cfg%omega_P,tr(2:3)*cfg%omega_rho*cfg%omega_u, &
     tr(6)*cfg%omega_rho*cfg%omega_u,fm(3)*cfg%omega_rho*cfg%omega_u**2,fp(3)*cfg%omega_rho*cfg%omega_u**2
   cfg%junction_flux='PHASE_UPWIND';call sim%init(cfg)
   call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code,junction_trace=tr,left_cell=ib)
   write(*,'(a,i1,a,a,a,i6,a,6es14.6)')'OP|',which,'|',trim(trc(k)),'|E1|',code,'|',tr(1)*cfg%omega_P,fm(1:2)*cfg%omega_rho*cfg%omega_u, &
     tr(6)*cfg%omega_rho*cfg%omega_u,fm(3)*cfg%omega_rho*cfg%omega_u**2,fp(3)*cfg%omega_rho*cfg%omega_u**2
   cfg%junction_flux='ROE';cfg%junction_trace='CENTER';call sim%init(cfg)
  enddo
 enddo
contains
 function out_si(rf) result(o)
  !! [Gl*, Gg*, F3_L*, F3_R*, pL*, pR*, alphaL*, alphaR*, lam1(UL*), lam2(UR*), lam3(UM), newton its]
  real(wp),intent(in)::rf(12)
  real(wp)::o(12)
  o=rf
  o(1:2)=rf(1:2)*cfg%omega_rho*cfg%omega_u
  o(3:4)=rf(3:4)*cfg%omega_rho*cfg%omega_u**2
  o(5:6)=rf(5:6)*cfg%omega_P
  o(9:11)=rf(9:11)*cfg%omega_u
 end function
 subroutine geometry()
  integer::i
  if(.not.allocated(cfg%face_positions))allocate(cfg%face_positions(0:hi-lo+1))
  if(.not.allocated(cfg%vns))allocate(cfg%vns(nc))
  cfg%face_positions(0)=0
  do i=lo,hi
   cfg%face_positions(i-lo+1)=cfg%face_positions(i-lo)+cfg%cell_width(i)
   cfg%vns(i)=cfg%face_positions(i-lo)+.5_wp*cfg%cell_width(i)
  enddo
  do i=1,lo-1
   cfg%vns(i)=cfg%vns(lo)-(lo-i)*cfg%cell_width(lo);cfg%vns(hi+i)=cfg%vns(hi)+i*cfg%cell_width(hi)
  enddo
 end subroutine
 subroutine prim(u,th,cid,v,code)
  real(wp),intent(in)::u(3),th
  integer,intent(in)::cid
  real(wp),intent(out)::v(6)
  integer,intent(out)::code
  real(wp)::f(3),s(3),a(3,3),d(3,3)
  call sim%point_state(u,th,v,f,s,a,d,code,.false.,cid)
 end subroutine
 subroutine flux_of(u,th,cid,f,code)
  real(wp),intent(in)::u(3),th
  integer,intent(in)::cid
  real(wp),intent(out)::f(3)
  integer,intent(out)::code
  real(wp)::v(6),s(3),a(3,3),d(3,3)
  call sim%point_state(u,th,v,f,s,a,d,code,.false.,cid)
 end subroutine
 subroutine eig(u,th,cid,lam,r,code)
  !! Real eigenpairs of A=dF/dU sorted by speed; r(:,k) unnormalised.
  real(wp),intent(in)::u(3),th
  integer,intent(in)::cid
  real(wp),intent(out)::lam(3),r(3,3)
  integer,intent(out)::code
  real(wp)::v(6),f(3),s(3),a(3,3),d(3,3),wr(3),wi(3),vr(3,3)
  integer::idx(3),i,j,t
  lam=0;r=0
  call sim%point_state(u,th,v,f,s,a,d,code,.true.,cid)
  if(code/=0)return
  call dgeev_3x3(a,wr,wi,vr,code)
  if(code/=0)return
  if(any(wi/=0))then
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
 subroutine grads(u,th,cid,gp,ga,code)
  !! Gradients of p and alpha with respect to U (EOS identities, as in point_state).
  real(wp),intent(in)::u(3),th
  integer,intent(in)::cid
  real(wp),intent(out)::gp(3),ga(3)
  integer,intent(out)::code
  real(wp)::v(6),den
  call prim(u,th,cid,v,code)
  if(code/=0)return
  den=v(2)*v(1)/cfg%nCg**2+v(3)*(1-v(1))/cfg%nCl**2
  ga=[-v(1)/cfg%nCg**2,(1-v(1))/cfg%nCl**2,0.0_wp]/den
  gp=[v(3),v(2),0.0_wp]/den
 end subroutine
 subroutine field(u,th,cid,fam,kind,du,code)
  !! dU/dq along the integral curve of family fam; q=p (kind 1) or ln alpha (kind 2).
  real(wp),intent(in)::u(3),th
  integer,intent(in)::cid,fam,kind
  real(wp),intent(out)::du(3)
  integer,intent(out)::code
  real(wp)::lam(3),r(3,3),gp(3),ga(3),v(6),dq
  call eig(u,th,cid,lam,r,code)
  if(code/=0)return
  call grads(u,th,cid,gp,ga,code)
  if(code/=0)return
  if(kind==1)then
   dq=dot_product(gp,r(:,fam))
  else
   call prim(u,th,cid,v,code)
   dq=dot_product(ga,r(:,fam))/v(1)
  endif
  if(.not.(abs(dq)>0))then
   code=-8;return
  endif
  du=r(:,fam)/dq
 end subroutine
 subroutine curve(u0,th,cid,fam,kind,q1,u1,code)
  !! RK4 integration of the integral curve from u0 to parameter value q1.
  real(wp),intent(in)::u0(3),th,q1
  integer,intent(in)::cid,fam,kind
  real(wp),intent(out)::u1(3)
  integer,intent(out)::code
  integer,parameter::nst=48
  real(wp)::q0,h,k1(3),k2(3),k3(3),k4(3),v(6)
  integer::i
  u1=u0
  call prim(u0,th,cid,v,code)
  if(code/=0)return
  q0=v(4)
  if(kind==2)q0=log(v(1))
  h=(q1-q0)/nst
  if(h==0)return
  do i=1,nst
   call field(u1,th,cid,fam,kind,k1,code);if(code/=0)return
   call field(u1+h/2*k1,th,cid,fam,kind,k2,code);if(code/=0)return
   call field(u1+h/2*k2,th,cid,fam,kind,k3,code);if(code/=0)return
   call field(u1+h*k3,th,cid,fam,kind,k4,code);if(code/=0)return
   u1=u1+h*(k1+2*k2+2*k3+k4)/6
  enddo
 end subroutine
 subroutine resid(s,ul0,ur0,tl,tr_,cl_,cr_,tcond,res,ust,code)
  real(wp),intent(in)::s(3),ul0(3),ur0(3),tl,tr_
  integer,intent(in)::cl_,cr_,tcond
  real(wp),intent(out)::res(3),ust(3,3)
  integer,intent(out)::code
  real(wp)::uls(3),um(3),urs(3),fl_(3),fr_(3),vl_(6),vr_(6)
  res=huge(1._wp)
  call curve(ul0,tl,cl_,1,1,s(1),uls,code);if(code/=0)return
  call curve(ur0,tr_,cr_,3,1,s(3),um,code);if(code/=0)return
  call curve(um,tr_,cr_,2,2,s(2),urs,code);if(code/=0)return
  call flux_of(uls,tl,cl_,fl_,code);if(code/=0)return
  call flux_of(urs,tr_,cr_,fr_,code);if(code/=0)return
  res=fl_-fr_
  if(tcond==2)then
   call prim(uls,tl,cl_,vl_,code);if(code/=0)return
   call prim(urs,tr_,cr_,vr_,code);if(code/=0)return
   res(3)=pip*(vl_(4)-vr_(4))
  endif
  ust(:,1)=uls;ust(:,2)=urs;ust(:,3)=um
 end subroutine
 subroutine solve(ul0,ur0,tl,tr_,cl_,cr_,tcond,out,istat)
  real(wp),intent(in)::ul0(3),ur0(3),tl,tr_
  integer,intent(in)::cl_,cr_,tcond
  real(wp),intent(out)::out(12)
  integer,intent(out)::istat
  real(wp)::s(3),res(3),res2(3),jac(3,3),ust(3,3),ust2(3,3),ds(3),scl(3),vl_(6),vr_(6),fl_(3),fr_(3),h,lam(3),rr(3,3),lam2(3)
  real(wp)::jinv(3,3),det,lam3(3)
  integer::n,j,code
  out=ieee_value(0.0_wp,ieee_quiet_nan);istat=1
  call prim(ul0,tl,cl_,vl_,code);if(code/=0)return
  call prim(ur0,tr_,cr_,vr_,code);if(code/=0)return
  call flux_of(ul0,tl,cl_,fl_,code);if(code/=0)return
  s=[vl_(4),log(vr_(1)),vr_(4)]
  scl=[max(abs(fl_(1)),1e-8_wp),max(abs(fl_(2)),abs(fl_(1))*1e-6_wp,1e-14_wp),pip*vl_(4)]
  istat=2
  do n=1,40
   call resid(s,ul0,ur0,tl,tr_,cl_,cr_,tcond,res,ust,code)
   if(code/=0)then
    istat=100+abs(code);return
   endif
   if(maxval(abs(res)/scl)<1e-10_wp)exit
   if(dbg)write(*,'(a,i3,3es12.4,a,3es16.8)')'  NEWTON',n,res/scl,' s',s
   do j=1,3
    h=1e-7_wp*max(abs(s(j)),1e-3_wp)
    s(j)=s(j)+h
    call resid(s,ul0,ur0,tl,tr_,cl_,cr_,tcond,res2,ust2,code)
    s(j)=s(j)-h
    if(code/=0)then
     istat=200+abs(code);return
    endif
    jac(:,j)=(res2-res)/h
   enddo
   det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(2,3)*jac(3,2))-jac(1,2)*(jac(2,1)*jac(3,3)-jac(2,3)*jac(3,1)) &
      +jac(1,3)*(jac(2,1)*jac(3,2)-jac(2,2)*jac(3,1))
   if(det==0)then
    istat=3;return
   endif
   call inv3(jac,jinv)
   ds=-matmul(jinv,res)
   ! damped step: keep p positive, cap ln(alpha) move
   ds(2)=max(min(ds(2),3._wp),-3._wp)
   s=s+ds
   if(s(2)<log(1e-14_wp))then
    istat=50 ! REF_LIMIT_ALPHA0: riser trace alpha*->0 required
    s(2)=log(1e-14_wp)
    call resid(s,ul0,ur0,tl,tr_,cl_,cr_,tcond,res,ust,code)
    exit
   endif
  enddo
  if(istat/=50)then
   if(n>40)then
    istat=4;return
   endif
   istat=0
  endif
  call flux_of(ust(:,1),tl,cl_,fl_,code)
  call flux_of(ust(:,2),tr_,cr_,fr_,code)
  call prim(ust(:,1),tl,cl_,vl_,code)
  call prim(ust(:,2),tr_,cr_,vr_,code)
  call eig(ust(:,1),tl,cl_,lam,rr,code)
  call eig(ust(:,2),tr_,cr_,lam2,rr,code)
  call eig(ust(:,3),tr_,cr_,lam3,rr,code)
  out=[fl_(1),fl_(2),fl_(3),fr_(3),vl_(4),vr_(4),vl_(1),vr_(1),lam(1),lam2(2),lam3(3),real(n,wp)]
  if(istat==0.and..not.(lam(1)<0.and.lam2(2)>0.and.lam3(3)>0))istat=60 ! configuration violated (H4)
 end subroutine
 subroutine inv3(a,b)
  real(wp),intent(in)::a(3,3)
  real(wp),intent(out)::b(3,3)
  real(wp)::d
  d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2);b(1,2)=a(1,3)*a(3,2)-a(1,2)*a(3,3);b(1,3)=a(1,2)*a(2,3)-a(1,3)*a(2,2)
  b(2,1)=a(2,3)*a(3,1)-a(2,1)*a(3,3);b(2,2)=a(1,1)*a(3,3)-a(1,3)*a(3,1);b(2,3)=a(1,3)*a(2,1)-a(1,1)*a(2,3)
  b(3,1)=a(2,1)*a(3,2)-a(2,2)*a(3,1);b(3,2)=a(1,2)*a(3,1)-a(1,1)*a(3,2);b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
  b=b/d
 end subroutine
end program
