program seal_jac_check
  !! Round 20261122 (P0.1/P0.2): junction face of the GLOBAL operator on frozen accepted
  !! states (offline_api fixture: cells 53|54, stratified pipe | CHEXAL_1992 riser).
  !! Modes (argument 1), fixture folder (argument 2):
  !!  OPS     operators A (T1, curve in ln alpha with 8..128 RK4 steps), B (seal solver,
  !!          H4_ONLY) and C (LIQUID_BRIDGE), B/C with legacy adaptive steps (0) and fixed
  !!          8..128 steps; status, configuration, SI flux, pressure/alpha traces.
  !!  JAC     implicit Jacobian (IMPLICIT, fixed steps) against an INDEPENDENT reference:
  !!          central differences of the COMPLETE face solve (face tol 1e-12) at two steps
  !!          h and h/2 with Richardson extrapolation, every sample in the base
  !!          configuration; also the legacy FD face Jacobian; step study of the
  !!          elementary-map differences (junction_seal_ijac_eps).
  !!  SWITCH  pipe momentum U3_L as parameter: switch g0=0 by bisection; flux and
  !!          implicit Jacobian at offsets +-1e-4..1e-9 (lateral limits J+, J-).
  !!  RH      moving front of H7 (P0.5): with the converged (U_L1, U*, s) recorded by the
  !!          solver, independent fluxes (point_state for U_L1 with the pipe closure;
  !!          pure-liquid flux of U* from its curve pressure) and the jump residuals
  !!          s[U_k]-[F_k], k=1,2,3, in SI, relative to max(|F_k|,|s U_k|); gas flux
  !!          through the front Fg_L1 - s U2_L1 (zero) versus the ABSOLUTE pipe-side gas
  !!          flux Fg_L1 (not zero) and the fixed-face Gg* (zero). Pipe momentum U3_L
  !!          scaled by 1, 0.99, 0.97, 0.95, 0.9, 0.8 to sample H7 states.
  !!  NOISE   flux under relative perturbations 1e-15..1e-8 of U3_R (legacy rule and fixed
  !!          steps): discontinuities of the discrete operator and round-off level.
  !! Output scaling (SI): G in kg/(m2 s), F3 in Pa; Jacobians scaled as J_ik*U_k/Fs_i with
  !! Fs = (|F1|, max(|Fg|, G_g,ref), |F3|), G_g,ref = rho_g(Ps) jg0.
  use offline_api
  use module_simulate_pipe_riser,only:junction_seal_front
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  character(512)::folder,mode
  real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),ul(3),ur(3),thl,thr,fsi(3),ggref,t0,t1
  integer(c_int)::st
  integer::code,nsteps
  character(16)::envs
  call get_command_argument(1,mode);call get_command_argument(2,folder)
  nsteps=16;call get_environment_variable('PRIMS_SEAL_STEPS',envs,status=code);if(code==0)read(envs,*)nsteps
  call probe_init(trim(folder),len_trim(folder),1_c_int,st)
  if(st/=0)then
    write(*,'(a,i0)')'probe_init status ',st;stop 1
  endif
  call probe_state(0_c_int,x,sc,meta);call unpack(x,uu)
  ul=uu(ib,:);ur=uu(ib+1,:);thl=cfg%thetav(ib);thr=cfg%thetav(ib+1)
  fsi=[cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u**2]
  ggref=(cfg%nPs/cfg%nCg**2)*cfg%njg            ! nondimensional rho_g(Ps) jg0
  write(*,'(a,f12.6,a,3es14.6,a,3es14.6,a,i0)')'STATE|t_s ',meta(1),'|UL ',ul,'|UR ',ur,'|fixed curve steps ',nsteps
  select case(trim(mode))
  case('OPS');call ops()
  case('JAC');call jacs()
  case('SWITCH');call switch()
  case('NOISE');call noise()
  case('RH');call rh()
  case default
    stop 'mode: OPS | JAC | SWITCH | NOISE | RH'
  end select
contains
  subroutine setop(trans,seal,t1steps,ssteps,jmode,ftol,eps)
    character(*),intent(in)::trans,seal,jmode
    integer,intent(in)::t1steps,ssteps
    real(wp),intent(in)::ftol,eps
    cfg%junction_transmission=trans;cfg%junction_seal=seal;cfg%junction_t1_curve_steps=t1steps
    cfg%junction_seal_curve_steps=ssteps;cfg%junction_seal_jacobian=jmode;cfg%junction_seal_face_tol=ftol
    cfg%junction_seal_ijac_eps=eps
    call sim%init(cfg);call sim%prepare_dissipation(uu,code)
  end subroutine
  subroutine face(a,b,f,stt,tr,info,jl,jr)
    real(wp),intent(in)::a(3),b(3)
    real(wp),intent(out)::f(3),tr(6),info(3)
    integer,intent(out)::stt
    real(wp),intent(out),optional::jl(3,3),jr(3,3)
    real(wp)::fp(3)
    if(present(jl))then
      call sim%face_pair(a,b,thl,thr,f,fp,stt,jl,jr,junction_trace=tr,dissipation_matrix=sim%lagged_dissipation(:,:,ib), &
        left_cell=ib,junction_seal_info=info)
    else
      call sim%face_pair(a,b,thl,thr,f,fp,stt,junction_trace=tr,dissipation_matrix=sim%lagged_dissipation(:,:,ib), &
        left_cell=ib,junction_seal_info=info)
    endif
  end subroutine
  subroutine ops()
    integer::steps(6)=[0,8,16,32,64,128],k,m
    real(wp)::f(3),tr(6),info(3)
    character(16)::lab
    do k=2,6
      call setop('T1','OFF',steps(k),0,'FD',1e-11_wp,1e-6_wp)
      call cpu_time(t0);call face(ul,ur,f,code,tr,info);call cpu_time(t1)
      call out('A',steps(k),f,tr,info,t1-t0)
    enddo
    do m=1,2
      lab='H4_ONLY';if(m==2)lab='LIQUID_BRIDGE'
      do k=1,6
        call setop('T1',trim(lab),8,steps(k),'FD',1e-11_wp,1e-6_wp)
        call cpu_time(t0);call face(ul,ur,f,code,tr,info);call cpu_time(t1)
        call out(merge('B','C',m==1),steps(k),f,tr,info,t1-t0)
        if(steps(k)>=32)then   ! tighter face tolerance (P1.3 floor study: 1e-12)
          call setop('T1',trim(lab),8,steps(k),'FD',1e-12_wp,1e-6_wp)
          call face(ul,ur,f,code,tr,info)
          call out(merge('b','c',m==1),steps(k),f,tr,info,0.0_wp)
        endif
      enddo
    enddo
  end subroutine
  subroutine out(op,n,f,tr,info,cpu)
    character(*),intent(in)::op
    integer,intent(in)::n
    real(wp),intent(in)::f(3),tr(6),info(3),cpu
    write(*,'(a,a,a,i4,a,i6,a,f3.0,a,3es24.16,a,2es20.12,a,2es14.6,a,f8.4)')'OPS|',op,'|',n,'|',code,'|',info(1), &
      '|F ',f*fsi,'|pL*,pR* ',tr(1)*cfg%omega_P,tr(6)*cfg%omega_P,'|aL*,aR* ',tr(4),tr(5),'|cpu ',cpu
  end subroutine
  subroutine scaled(j,u,f,s)
    !! J_ik U_k / Fs_i (dimensionless, independent of the unit system)
    real(wp),intent(in)::j(3,3),u(3),f(3)
    real(wp),intent(out)::s(3,3)
    real(wp)::fs(3)
    integer::i
    fs=[abs(f(1)),max(abs(f(2)),ggref),abs(f(3))]
    do i=1,3
      s(i,:)=j(i,:)*abs(u)/fs(i)
    enddo
  end subroutine
  subroutine fdref(f0,cfg0,jref,ok)
    !! Richardson-extrapolated central differences of the complete face solve. Per
    !! column, the largest relative step h in {2e-4/4**k, k=0..8} for
    !! which the four samples (+-h, +-h/2) are admissible AND in the configuration of the
    !! base; the step used is reported (column hk).
    real(wp),intent(in)::f0(3),cfg0
    real(wp),intent(out)::jref(3,6)
    logical,intent(out)::ok
    real(wp)::w(6),wp_(6),wm(6),fp(3),fm(3),tr(6),info(3),d(3,2),h,hk(6)
    integer::k,m,s1,s2,lev
    logical::good
    w=[ul,ur];ok=.true.;hk=0;jref=ieee_value(0.0_wp,ieee_quiet_nan)
    do k=1,6
      do lev=0,8
        good=.true.
        do m=1,2
          h=2e-4_wp/4.0_wp**lev/2**(m-1)*abs(w(k))
          wp_=w;wm=w;wp_(k)=w(k)+h;wm(k)=w(k)-h
          call face(wp_(1:3),wp_(4:6),fp,s1,tr,info)
          if(s1/=0.or.nint(info(1))/=nint(cfg0))good=.false.
          call face(wm(1:3),wm(4:6),fm,s2,tr,info)
          if(s2/=0.or.nint(info(1))/=nint(cfg0))good=.false.
          if(.not.good)exit
          d(:,m)=(fp-fm)/(2*h)
        enddo
        if(good)then
          jref(:,k)=(4*d(:,2)-d(:,1))/3;hk(k)=2e-4_wp/4.0_wp**lev;exit
        endif
      enddo
      if(.not.good)ok=.false.
    enddo
    write(*,'(a,6es9.2)')'REF|relative step per column (0 = none admissible) ',hk
    if(.false.)write(*,*)f0
  end subroutine
  subroutine jacs()
    real(wp)::f(3),tr(6),info(3),jl(3,3),jr(3,3),jref(3,6),ji(3,6),jf(3,6),a(3,3),b(3,3),c0,es(5),ef
    real(wp)::epsv(5)=[1e-4_wp,1e-5_wp,1e-6_wp,1e-7_wp,1e-8_wp]
    logical::ok
    integer::k,stt
    ! reference: smooth operator (64 fixed steps), face tolerance 1e-12
    call setop('T1','LIQUID_BRIDGE',8,nsteps,'FD',1e-12_wp,1e-6_wp)
    call face(ul,ur,f,stt,tr,info)
    c0=info(1)
    write(*,'(a,i0,a,f3.0,a,3es22.14)')'BASE|status ',stt,'|config ',c0,'|F ',f*fsi
    if(stt/=0)return
    call cpu_time(t0);call fdref(f,c0,jref,ok);call cpu_time(t1)
    write(*,'(a,l1,a,f8.3)')'REF|richardson admissible same-configuration samples ',ok,'|cpu ',t1-t0
    do k=1,5
      call setop('T1','LIQUID_BRIDGE',8,nsteps,'IMPLICIT',1e-12_wp,epsv(k))
      call cpu_time(t0);call face(ul,ur,f,stt,tr,info,jl,jr);call cpu_time(t1)
      ji(:,1:3)=jl;ji(:,4:6)=jr
      call errs(ji,jref,es(k))
      write(*,'(a,es8.1,a,i0,a,es11.3,a,f8.4)')'IMPLICIT|eps ',epsv(k),'|status ',stt,'|max scaled err vs ref ',es(k),'|cpu ',t1-t0
      if(k==3)call show('J_implicit(1e-6)',ji)
    enddo
    call show('J_reference',jref)
    call setop('T1','LIQUID_BRIDGE',8,nsteps,'FD',1e-12_wp,1e-6_wp)
    call cpu_time(t0);call face(ul,ur,f,stt,tr,info,jl,jr);call cpu_time(t1)
    jf(:,1:3)=jl;jf(:,4:6)=jr;call errs(jf,jref,ef)
    write(*,'(a,i0,a,es11.3,a,f8.4)')'LEGACY_FD(fixed steps, tol 1e-12)|status ',stt,'|max scaled err vs ref ',ef,'|cpu ',t1-t0
    call setop('T1','LIQUID_BRIDGE',8,0,'FD',1e-11_wp,1e-6_wp)
    call cpu_time(t0);call face(ul,ur,f,stt,tr,info,jl,jr);call cpu_time(t1)
    jf(:,1:3)=jl;jf(:,4:6)=jr;call errs(jf,jref,ef)
    write(*,'(a,i0,a,es11.3,a,f8.4)')'LEGACY_FD(seal5 operator)|status ',stt,'|max scaled err vs ref ',ef,'|cpu ',t1-t0
    call setop('T1','LIQUID_BRIDGE',8,0,'IMPLICIT',1e-11_wp,1e-6_wp)
    call cpu_time(t0);call face(ul,ur,f,stt,tr,info,jl,jr);call cpu_time(t1)
    jf(:,1:3)=jl;jf(:,4:6)=jr;call errs(jf,jref,ef)
    write(*,'(a,i0,a,es11.3,a,f8.4)')'IMPLICIT(legacy steps, frozen)|status ',stt,'|max scaled err vs ref(fixed) ',ef,'|cpu ',t1-t0
    if(.false.)then
      a=0;b=0
    endif
  end subroutine
  subroutine errs(j,jr_,e)
    real(wp),intent(in)::j(3,6),jr_(3,6)
    real(wp),intent(out)::e
    real(wp)::s1(3,3),s2(3,3),t1_(3,3),t2_(3,3),f(3),tr(6),info(3)
    integer::stt
    call face(ul,ur,f,stt,tr,info)
    call scaled(j(:,1:3),ul,f,s1);call scaled(jr_(:,1:3),ul,f,s2)
    call scaled(j(:,4:6),ur,f,t1_);call scaled(jr_(:,4:6),ur,f,t2_)
    ! normalized by the largest scaled reference entry of each row
    e=max(rowerr(s1,s2,t1_,t2_,1),rowerr(s1,s2,t1_,t2_,2),rowerr(s1,s2,t1_,t2_,3))
  end subroutine
  real(wp) function rowerr(a,b,c,d,i) result(e)
    real(wp),intent(in)::a(3,3),b(3,3),c(3,3),d(3,3)
    integer,intent(in)::i
    real(wp)::nrm
    nrm=max(maxval(abs(b(i,:))),maxval(abs(d(i,:))),1e-300_wp)
    e=max(maxval(abs(a(i,:)-b(i,:))),maxval(abs(c(i,:)-d(i,:))))/nrm
    if(nrm<=1e-300_wp)e=max(maxval(abs(a(i,:))),maxval(abs(c(i,:))))
  end function
  subroutine show(lab,j)
    character(*),intent(in)::lab
    real(wp),intent(in)::j(3,6)
    integer::i
    do i=1,3
      write(*,'(a,a,i1,a,6es14.6)')'MAT|',lab,i,'|',j(i,:)
    enddo
  end subroutine
  subroutine switch()
    !! parameter: pipe U3 (momentum of cell 53)
    real(wp)::u(3),lo,hi,mid,c,f(3),fa(3),fb(3),tr(6),ia(3),ibb(3),jla(3,3),jra(3,3),jlb(3,3),jrb(3,3)
    real(wp)::ja(3,6),jb(3,6),off,d
    integer::k,stt,sa,sb
    call setop('T1','LIQUID_BRIDGE',8,nsteps,'IMPLICIT',1e-12_wp,1e-6_wp)
    call face(ul,ur,f,stt,tr,ia)
    write(*,'(a,i0,a,f3.0,a,es14.6)')'BASE|status ',stt,'|config ',ia(1),'|g0 ',ia(3)*fsi(2)
    lo=ul(3);hi=ul(3)
    do k=1,60
      if(nint(ia(1))==7)then
        hi=hi*1.0005_wp;u=ul;u(3)=hi
      else
        lo=lo/1.0005_wp;u=ul;u(3)=lo
      endif
      call face(u,ur,f,stt,tr,ibb)
      if(stt==0.and.nint(ibb(1))/=nint(ia(1)))exit
    enddo
    do k=1,200
      mid=0.5_wp*(lo+hi);u=ul;u(3)=mid
      call face(u,ur,f,stt,tr,ibb)
      if(stt/=0)exit
      if(nint(ibb(1))==4.and.ibb(3)>0)then
        hi=mid
      else
        lo=mid
      endif
      if(hi-lo<=4*epsilon(1.0_wp)*hi)exit
    enddo
    c=0.5_wp*(lo+hi)
    write(*,'(a,es24.16,a,i0)')'SWITCH|U3c (kg/m2/s) ',c*fsi(1),'|status ',stt
    do k=4,11
      off=10.0_wp**(-k)
      u=ul;u(3)=c*(1+off);call face(u,ur,fa,sa,tr,ia,jla,jra)
      u=ul;u(3)=c*(1-off);call face(u,ur,fb,sb,tr,ibb,jlb,jrb)
      ja(:,1:3)=jla;ja(:,4:6)=jra;jb(:,1:3)=jlb;jb(:,4:6)=jrb
      d=maxval(abs((ja-jb)))/max(maxval(abs(ja)),1e-300_wp)
      write(*,'(a,es8.1,a,2i5,a,2f3.0,a,3es12.4,a,es12.4,a,es12.4,a,es12.4)')'SIDE|offset ',off,'|status ',sa,sb, &
        '|config ',ia(1),ibb(1),'|F+ - F- (SI) ',(fa-fb)*fsi,'|Gg*(+) ',fa(2)*fsi(2),'|max|J+ - J-|/max|J+| ',d, &
        '|g0(+) ',ia(3)*fsi(2)
      if(k==6)then
        call show('J+(1e-6)',ja);call show('J-(1e-6)',jb)
      endif
    enddo
  end subroutine
  subroutine noise()
    real(wp)::u(3),f0(3),f(3),tr(6),info(3),d
    integer::k,m,stt,steps(4)
    steps=[0,nsteps,nsteps,-8]
    do m=1,4
      ! 0: seal5 operator; 64: fixed steps tol 1e-11; 64 (m=3): tol 1e-13; -8: operator A (T1)
      if(steps(m)<0)then
        call setop('T1','OFF',-steps(m),0,'FD',1e-11_wp,1e-6_wp)
      else
        call setop('T1','LIQUID_BRIDGE',8,steps(m),'FD',merge(1e-13_wp,1e-11_wp,m==3),1e-6_wp)
      endif
      call face(ul,ur,f0,stt,tr,info)
      write(*,'(a,i3,a,i0,a,f3.0,a,3es24.16)')'NOISE|steps ',steps(m),'|base status ',stt,'|config ',info(1),'|F ',f0*fsi
      do k=0,14
        d=10.0_wp**(-15+0.5_wp*k)
        u=ur;u(3)=ur(3)*(1+d)
        call face(ul,u,f,stt,tr,info)
        write(*,'(a,i3,a,es9.2,a,i0,a,3es13.5,a,3es13.5,a,es24.16,es24.16,es24.16)')'NOISE|steps ',steps(m),'|dU3R/U3R ',d,'|status ',stt, &
          '|(F-F0)/Fs ',(f-f0)/[abs(f0(1)),max(abs(f0(2)),ggref),abs(f0(3))],'|slope (dF/F)/(dU/U) ', &
          (f-f0)/[abs(f0(1)),max(abs(f0(2)),ggref),abs(f0(3))]/d,'|pL*,pR*,aR* ',tr(1),tr(6),tr(5)
      enddo
    enddo
  end subroutine
  subroutine rh()
    real(wp)::u(3),f(3),tr(6),info(3),s_,ul1(3),us(3),fl1(3),fs(3),fs2(3),v(6),sv(3),a(3,3),ds(3,3),fac(6),rk(3),nk(3)
    real(wp)::pip,usi(3)
    integer::k,stt,sp,sp2,cidl
    fac=[1.0_wp,0.99_wp,0.97_wp,0.95_wp,0.9_wp,0.8_wp]
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    usi=[cfg%omega_rho,cfg%omega_rho,cfg%omega_rho*cfg%omega_u]   ! U scale (SI)
    cidl=cfg%closure_pipe
    call setop('T1','LIQUID_BRIDGE',8,nsteps,'IMPLICIT',1e-11_wp,1e-5_wp)
    do k=1,6
      u=ul;u(3)=ul(3)*fac(k)
      call face(u,ur,f,stt,tr,info)
      if(stt/=0.or.nint(info(1))/=7)then
        write(*,'(a,f6.3,a,i0,a,f3.0)')'RH|U3L factor ',fac(k),'|status ',stt,'|config ',info(1);cycle
      endif
      s_=junction_seal_front(1);ul1=junction_seal_front(2:4);us=junction_seal_front(5:7)
      call sim%point_state(ul1,thl,v,fl1,sv,a,ds,sp,need_jac=.false.,closure_id=cidl)
      fs=[us(3),0.0_wp,us(3)**2/us(1)+pip*junction_seal_front(9)]
      call sim%point_state(us,thr,v,fs2,sv,a,ds,sp2,need_jac=.false.)
      rk=s_*(us-ul1)-(fs-fl1)
      nk=max(abs(fs),abs(fl1),abs(s_*us),abs(s_*ul1),1e-300_wp)
      write(*,'(a,f6.3,a,i0,a,es12.4,a,3es12.4,a,3es12.4,a,es12.4,a,es12.4,a,es12.4,a,i0,a,3es12.4)') &
        'RH|U3L factor ',fac(k),'|status ',stt,'|s (m/s) ',s_*cfg%omega_u,'|s[U]-[F] SI ',rk*fsi,'|rel ',rk/nk, &
        '|Fg_L1-s*U2_L1 ',(fl1(2)-s_*ul1(2))*fsi(2),'|Fg_L1 abs ',fl1(2)*fsi(2),'|Gg* face ',f(2)*fsi(2), &
        '|pt_state(U*) status ',sp2,'|F*(curve p)-F*(point_state) rel ',(fs-fs2)/max(abs(fs),1e-300_wp)
      write(*,'(a,3es14.6,a,3es14.6,a,2es14.6,a,es12.4)')'RH|U_L1 SI ',ul1*usi,'|U* SI ',us*usi, &
        '|p_L1,p* Pa ',junction_seal_front(8:9)*cfg%omega_P,'|face F* vs fs rel ',maxval(abs(f-fs)/max(abs(fs),1e-300_wp))
    enddo
  end subroutine
end program seal_jac_check
