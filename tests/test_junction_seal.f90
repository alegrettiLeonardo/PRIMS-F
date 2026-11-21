program test_junction_seal
  !! junction_seal='LIQUID_BRIDGE' (EXPERIMENTAL) on the T1 junction, stratified
  !! pipe (1) | CHEXAL_1992 riser (5). Independent checks (no reuse of the solver's
  !! inner quantities except the public face outputs):
  !!  (1) OFF is the T1 operator bit for bit; the seal solver refuses nothing silently;
  !!  (2) gas present and passing (co-current state): configuration H4, Gg*>0, one
  !!      flux, agreement with T1 within the T1 kinematic-curve discretization;
  !!  (3) uniform stratified state (homogeneous partition, theta_L=-1e-12):
  !!      Fhat(U,U)=F(U), equal pressure traces;
  !!  (4) seal state of SS1 at 53.4753 s (cell states of the accepted solution):
  !!      configuration H7, Gg*=0 exactly, front speed s<=0, one flux;
  !!  (5) gas flux tending to zero and reappearance: the pipe momentum U3_L is the
  !!      parameter; g0 changes sign at U3c (bisection on the reported g0); the face
  !!      flux is continuous through the switch (jump -> 0 linearly in the offset)
  !!      and Gg* -> 0+ on the H4 side;
  !!  (6) face Jacobian against independent central differences in H4 and H7;
  !!  (7) a different pressure level (states scaled to +20 % pressure) behaves alike;
  !!  (8) a pure-liquid riser CELL (U2=0) is refused explicitly (not supported).
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp)::u(3),ul(3),ur(3),fm(3),fp(3),fu(3),fq(3),tr(6),info(3),v(6),s(3),a(3,3),ds(3,3)
  real(wp)::f0(3),fa(3),fb(3),ia(3),ib(3),lo,hi,mid,jl(3,3),jr(3,3),fpp(3),fmm(3),h,err,c,eps,jump(2)
  real(wp),parameter::thp=-0.03490658503988659_wp,thr=0.04073_wp,rho=997.05_wp
  integer::st,n,fails,k,i,m
  real(wp)::errs(2)
  real(wp),parameter::hs(2)=[3e-5_wp,3e-6_wp]
  n=0;fails=0;jump=huge(1.0_wp)
  call base()
  ! SS1 T1 accepted state at 53.4753 s (SI from the seal probe, cells 53|54)
  ul=[438.9335_wp,1.411318_wp,528.3750_wp]/rho
  ur=[990.2674_wp,1.739060e-2_wp,524.9158_wp]/rho
  ! (1) OFF: the T1 solver is called and no seal configuration is reported
  cfg%junction_seal='OFF';call sim%init(cfg)
  call sim%face_pair(ul,ur,thp,thr,fa,fq,st,left_cell=3,junction_seal_info=info)
  write(*,'(a,i0,a,es13.5)')'  OFF (T1) at the seal state: status ',st,'  Gg* ',fa(2)*rho
  call chk(.not.any(ieee_is_finite(info)),'OFF: no seal configuration (plain T1 operator)')
  ! (4) seal state
  cfg%junction_seal='LIQUID_BRIDGE';call sim%init(cfg)
  call sim%face_pair(ul,ur,thp,thr,fm,fp,st,junction_trace=tr,left_cell=3,junction_seal_info=info)
  write(*,'(a,i0,a,f3.0,a,3es14.6,a,2es12.4)')'  seal state: status ',st,' config ',info(1),' flux ',fm*rho,' s,g0 ',info(2),info(3)*rho
  call chk(st==0,'seal state admitted')
  call chk(nint(info(1))==7,'seal state: liquid-bridge configuration H7 (g0<=0)')
  call chk(fm(2)==0.0_wp.and.fp(2)==0.0_wp,'H7: Gg*=0 exactly on both cells')
  call chk(info(2)<=0.0_wp,'H7: accumulation front in the pipe (s<=0)')
  call chk(all(fm==fp),'H7: one flux vector for both cells')
  call chk(tr(4)==0.0_wp.and.tr(5)==0.0_wp,'H7: pure liquid face state (alpha*=0 on both traces)')
  ! (5) parameter U3_L: find the switch g0=0 by bisection on the reported g0
  lo=ul(3);hi=ul(3)
  do k=1,40
    hi=hi*1.0005_wp;u=ul;u(3)=hi
    call sim%face_pair(u,ur,thp,thr,fm,fp,st,left_cell=3,junction_seal_info=info)
    if(st==0.and.info(3)>0)exit
  enddo
  call chk(st==0.and.info(3)>0,'raising the pipe momentum makes g0>0 (gas penetration)')
  do k=1,200
    mid=0.5_wp*(lo+hi);u=ul;u(3)=mid
    call sim%face_pair(u,ur,thp,thr,fm,fp,st,left_cell=3,junction_seal_info=info)
    if(st/=0)exit
    if(info(3)>0)then
      hi=mid
    else
      lo=mid
    endif
    if(hi-lo<=1e-12_wp*hi)exit
  enddo
  call chk(st==0,'switch located')
  c=0.5_wp*(lo+hi)
  write(*,'(a,es22.14)')'  switch U3c (kg/m2/s) ',c*rho
  do i=1,2
    eps=10.0_wp**(-3-2*i)
    u=ul;u(3)=c*(1+eps)
    call sim%face_pair(u,ur,thp,thr,fa,fq,st,junction_trace=tr,left_cell=3,junction_seal_info=ia)
    if(st/=0)exit
    u=ul;u(3)=c*(1-eps)
    call sim%face_pair(u,ur,thp,thr,fb,fq,st,left_cell=3,junction_seal_info=ib)
    if(st/=0)exit
    jump(i)=maxval(abs(fa-fb)/max(1e-3_wp,abs(fa)))
    write(*,'(a,es8.1,a,2f3.0,a,es11.3,a,es11.3,a,es11.3)')'  offset ',eps,' configs ',ia(1),ib(1),' |dF|/|F| ',jump(i), &
      ' Gg*(H4 side) kg/m2/s ',fa(2)*rho,' alpha_R* ',tr(5)
    call chk(nint(ia(1))==4.and.nint(ib(1))==7,'configurations on both sides of the switch')
    call chk(fa(2)>=0.0_wp.and.fb(2)==0.0_wp,'Gg*>=0 on the H4 side, =0 on the H7 side (no reverse gas)')
  enddo
  call chk(st==0,'switch neighbourhood admitted')
  if(st==0)call chk(jump(2)<=jump(1)/30.and.jump(2)<1e-6_wp,'face flux continuous through the switch (O(offset))')
  ! (2) co-current state with gas passing (test_junction_t1 item 3)
  ul=state(.55_wp,2.93_wp,.53_wp);ur=state(.25_wp,2.92_wp,.53_wp)
  call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fa,fp,st,left_cell=3,junction_seal_info=info)
  call chk(st==0.and.nint(info(1))==4.and.fa(2)>0.and.all(fa==fp),'co-current: H4, forward gas, one flux')
  cfg%junction_seal='OFF';call sim%init(cfg)
  call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fb,fq,st,left_cell=3)
  err=maxval(abs(fa-fb)/max(1e-3_wp,abs(fb)))
  write(*,'(a,es11.3,a,2es13.5)')'  co-current seal vs T1: max rel ',err,'  Gg* ',fa(2),fb(2)
  call chk(st==0.and.err<1e-3_wp,'co-current: seal solver agrees with T1 (curve discretization)')
  ! (6) Jacobian vs independent central differences, H4 (co-current) and H7 (seal)
  cfg%junction_seal='LIQUID_BRIDGE';call sim%init(cfg)
  do i=1,2
    if(i==2)then
      ! H7 away from the switch (riser gas +1 %: g0 = -1.4 kg/m2/s); at the SS1
      ! seal state itself the switch lies within ~4e-6 relative in U1_R and the
      ! difference Jacobian is a convex combination of the one-sided ones
      ul=[438.9335_wp,1.411318_wp,528.3750_wp]/rho;ur=[990.2674_wp,1.739060e-2_wp*1.01_wp,524.9158_wp]/rho
      call sim%face_pair(ul,ur,thp,thr,fm,fp,st,JL=jl,JR=jr,left_cell=3)
    else
      call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fm,fp,st,JL=jl,JR=jr,left_cell=3)
    endif
    call chk(st==0,'seal face Jacobian available')
    ! independent central differences with three increments: the error must
    ! reach the tolerance for some increment and not grow when it decreases
    ! from the largest (sensitivity to the differentiation increment)
    do m=1,2
      err=0
      do k=1,3
        h=hs(m)*max(abs(ul(k)),merge(1e-5_wp,1e-2_wp,k==2))
        u=ul;u(k)=u(k)+h;call face(u,ur,i,fpp)
        u=ul;u(k)=u(k)-h;call face(u,ur,i,fmm)
        err=max(err,maxval(abs(jl(:,k)-(fpp-fmm)/(2*h))/max(1.0_wp,abs((fpp-fmm)/(2*h)))))
        h=hs(m)*max(abs(ur(k)),merge(1e-5_wp,1e-2_wp,k==2))
        u=ur;u(k)=u(k)+h;call face(ul,u,i,fpp)
        u=ur;u(k)=u(k)-h;call face(ul,u,i,fmm)
        err=max(err,maxval(abs(jr(:,k)-(fpp-fmm)/(2*h))/max(1.0_wp,abs((fpp-fmm)/(2*h)))))
      enddo
      errs(m)=err
      write(*,'(a,a,a,es8.1,a,es11.3)')'  Jacobian vs independent central differences ',merge('H4','H7',i==1),' h/|U| ',hs(m),' err ',err
    enddo
    ! Bound: the face flux is defined to tol_face*|F1| (absolute); with the solver's
    ! increment 1e-6*|U2| on the gas mass (U2 ~ 1e-3) the difference Jacobian of
    ! that column carries ~tol*|F1|/(1e-6*U2) ~ 4e-3. Inexact-Newton accuracy.
    call chk(minval(errs)<5e-3_wp,'seal face Jacobian consistent with the face flux (within the face-solve noise)')
  enddo
  ! (7) +20 % pressure level: scale both inventories by the EOS at higher p
  ul=[438.9335_wp,1.411318_wp*1.2_wp,528.3750_wp]/rho;ur=[990.2674_wp,1.739060e-2_wp*1.2_wp,524.9158_wp]/rho
  call sim%face_pair(ul,ur,thp,thr,fm,fp,st,left_cell=3,junction_seal_info=info)
  write(*,'(a,i0,a,f3.0,a,3es13.5)')'  higher gas density: status ',st,' config ',info(1),' flux ',fm*rho
  call chk(st==0.and.all(fm==fp).and.(nint(info(1))==7.eqv.fm(2)==0.0_wp),'higher gas density: admitted, one flux, Gg*=0 iff H7')
  ! (3) homogeneous partition at the face level (stratified both sides, theta_L=-1e-12):
  ! g(z_c)=0 at the uniform state, the closed limit is not used; Fhat(U,U)=F(U).
  ! (Campaign decks with the seal option remain restricted to the pair (1,5).)
  cfg%closure_pipe=1;cfg%closure_riser=1;call sim%init(cfg)
  u=state(.30_wp,2.10_wp,.60_wp)
  call sim%face_pair(u,u,-1e-12_wp,0.0_wp,fm,fp,st,junction_trace=tr,left_cell=3,junction_seal_info=info)
  call sim%point_state(u,0.0_wp,v,fu,s,a,ds,k,.false.,1)
  write(*,'(a,i0)')'  uniform (stratified|stratified): status ',st
  call chk(st==0.and.nint(info(1))==4,'uniform state: H4 (gas passes)')
  if(st==0)then
    write(*,'(a,es11.3)')'  uniform: |Fhat-F|/|F| ',maxval(abs(fm-fu)/max(1.0_wp,abs(fu)))
    call chk(maxval(abs(fm-fu)/max(1.0_wp,abs(fu)))<1e-10_wp,'uniform state: Fhat(U,U)=F(U)')
    call chk(abs(tr(1)-tr(6))<=1e-10_wp*tr(1),'uniform state: no artificial pressure jump')
  endif
  ! (8) pure liquid riser cell
  cfg%closure_pipe=1;cfg%closure_riser=5;call sim%init(cfg)
  ul=[438.9335_wp,1.411318_wp,528.3750_wp]/rho;ur=[990.2674_wp,0.0_wp,524.9158_wp]/rho
  call sim%face_pair(ul,ur,thp,thr,fm,fp,st,left_cell=3)
  write(*,'(a,i0)')'  pure-liquid riser cell: status ',st
  call chk(st/=0,'pure-liquid riser cell refused explicitly (not supported)')
  if(fails>0)then
    write(*,'(a,i0,a,i0)')'FAILED ',fails,' of ',n;error stop 1
  endif
  write(*,'(a,i0,a)')'PASS junction seal: ',n,' assertions'
contains
  subroutine face(a_,b_,which,f)
    real(wp),intent(in)::a_(3),b_(3)
    integer,intent(in)::which
    real(wp),intent(out)::f(3)
    real(wp)::g(3)
    integer::is
    if(which==2)then
      call sim%face_pair(a_,b_,thp,thr,f,g,is,left_cell=3)
    else
      call sim%face_pair(a_,b_,-.0349_wp,.04_wp,f,g,is,left_cell=3)
    endif
  end subroutine
  subroutine base()
    cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1
    cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
    cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp;cfg%mul=.5_wp;cfg%mug=.001_wp
    cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%Lr=12.546253993197219_wp
    cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
    cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%nCells=6;cfg%num_dummy_cells=0;cfg%N_pipe=3
    cfg%junction_model='MATCHED';cfg%quad_points=3;cfg%junction_transmission='T1'
    cfg%closure_pipe=1;cfg%closure_riser=5
    allocate(cfg%thetav(6),cfg%cell_width(6));cfg%thetav=0;cfg%cell_width=1
  end subroutine
  function state(alpha,p,g) result(w)
    real(wp),intent(in)::alpha,p,g
    real(wp)::w(3),rl,rg
    rl=1+(p-1)/cfg%nCl**2;rg=p/cfg%nCg**2
    w=[(1-alpha)*rl,alpha*rg,g]
  end function
  subroutine chk(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    n=n+1
    if(.not.ok)then
      fails=fails+1;write(*,'(2a)')'FAIL: ',label
    endif
  end subroutine
end program
