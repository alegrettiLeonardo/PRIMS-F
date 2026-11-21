program test_junction_t1
  !! junction_transmission='T1' (EXPERIMENTAL): conservative transmission.
  !! (1) homogeneous partition: same closure on both sides of the MATCHED face
  !!     (stratified, theta_L=-1e-12 so the junction path is exercised): a
  !!     uniform state gives Fhat(U,U)=F(U), equal pressure traces and one flux;
  !! (2) small jump eps: T1 face flux equals the ordinary Roe face of the same
  !!     closure to second order in eps (ratio of errors at eps and eps/10 ~ 100);
  !! (3) heterogeneous stratified-Chexal junction: one flux vector for both
  !!     cells in all three components (no singular momentum source), status 0;
  !! (4) T2 (default) keeps its momentum jump F3_L* /= F3_R* at the same state;
  !! (5) face Jacobian blocks against independent central differences of the
  !!     face flux (same admissible states).
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp)::u(3),ul(3),ur(3),fm(3),fp(3),fu(3),fr_(3),tr(6),v(6),s(3),a(3,3),ds(3,3),e1,e2,jl(3,3),jr(3,3)
  real(wp)::fpp(3),fmm(3),fq(3),h,err,eps
  integer::st,n,fails,k
  n=0;fails=0
  call base()
  cfg%closure_pipe=1;cfg%closure_riser=1;cfg%junction_transmission='T1';call sim%init(cfg)
  u=state(.30_wp,2.10_wp,.60_wp)
  call sim%face_pair(u,u,-1e-12_wp,0.0_wp,fm,fp,st,junction_trace=tr,left_cell=3)
  call sim%point_state(u,0.0_wp,v,fu,s,a,ds,st,.false.,1)
  call chk(st==0,'uniform state admitted')
  call chk(maxval(abs(fm-fu)/max(1.0_wp,abs(fu)))<1e-10_wp,'T1 homogeneous partition: Fhat(U,U)=F(U)')
  call chk(all(fm==fp),'T1: one flux vector for both cells')
  call chk(abs(tr(1)-tr(6))<=1e-10_wp*tr(1),'T1 homogeneous: no artificial pressure jump p_L*=p_R*')
  ! (2) second-order agreement with the ordinary Roe face of the same closure
  do k=1,2
    eps=1e-3_wp/10**(k-1)
    ul=u;ur=u*(1+eps*[1.0_wp,-1.0_wp,1.0_wp])
    cfg%junction_transmission='T1';cfg%junction_model='MATCHED';call sim%init(cfg)
    call sim%face_pair(ul,ur,-1e-12_wp,0.0_wp,fm,fp,st,left_cell=3)
    call chk(st==0,'T1 small jump admitted')
    cfg%junction_model='CONSERVATIVE';call sim%init(cfg)
    call sim%face_pair(ul,ur,0.0_wp,0.0_wp,fr_,fq,st)
    err=maxval(abs(fm-fr_)/max(1.0_wp,abs(fr_)))
    if(k==1)e1=err
    if(k==2)e2=err
    write(*,'(a,es9.1,a,es11.3)')'  small jump eps=',eps,'  |T1-Roe|/|F| =',err
  enddo
  call chk(e2<=e1/30.and.e1<1e-5_wp,'T1 equals the ordinary Roe face to second order in the jump')
  cfg%junction_model='MATCHED'
  ! (3) heterogeneous stratified-Chexal junction
  cfg%closure_pipe=1;cfg%closure_riser=5;cfg%junction_transmission='T1';call sim%init(cfg)
  ul=state(.55_wp,2.93_wp,.53_wp);ur=state(.25_wp,2.92_wp,.53_wp)
  call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fm,fp,st,junction_trace=tr,left_cell=3)
  write(*,'(a,i0,a,3es13.5,a,2es13.5)')'  heterogeneous T1 status ',st,' flux ',fm,' p_L*,p_R* ',tr(1),tr(6)
  call chk(st==0,'heterogeneous stratified-Chexal T1 face admitted')
  call chk(all(fm==fp),'heterogeneous T1: identical flux in all three components')
  call chk(fm(2)>0,'heterogeneous T1: forward gas flux at this co-current state')
  ! (4) T2 keeps the point force at the same state
  cfg%junction_transmission='T2';call sim%init(cfg)
  call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fm,fp,st,left_cell=3)
  write(*,'(a,i0,a,3es13.5,a,es13.5)')'  T2 status ',st,' flux ',fm,' F3_R* ',fp(3)
  ! At this state the Roe/T2 face proposes a REVERSE gas flux on the
  ! horizontal riser trace: refused (1061 = Chexal domain), as in SS1, while
  ! T1 transmits a forward gas flux (item 3).
  call chk(st==1061.and.fm(2)<0,'T2 at the same state: reverse gas proposed and refused (1061)')
  ! A state where T2 is admissible keeps its momentum point force.
  ur=state(.40_wp,2.92_wp,.53_wp)
  call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fm,fp,st,left_cell=3)
  write(*,'(a,i0,a,3es13.5,a,es13.5)')'  T2 (alpha_R=.40) status ',st,' flux ',fm,' F3_R* ',fp(3)
  call chk(st==0.and.all(fm(1:2)==fp(1:2)).and.fm(3)/=fp(3),'T2 unchanged: common mass fluxes, momentum jump')
  ur=state(.25_wp,2.92_wp,.53_wp)
  ! (5) face Jacobian blocks versus independent central differences (T1)
  cfg%junction_transmission='T1';call sim%init(cfg)
  call sim%face_pair(ul,ur,-.0349_wp,.04_wp,fm,fp,st,JL=jl,JR=jr,left_cell=3)
  call chk(st==0,'T1 face Jacobian available')
  ! increments 30x the internal one (1e-6 relative): an independent stencil
  err=0
  do k=1,3
    h=3e-5_wp*max(abs(ul(k)),merge(1e-5_wp,1e-2_wp,k==2))
    u=ul;u(k)=u(k)+h;call sim%face_pair(u,ur,-.0349_wp,.04_wp,fpp,fq,st,left_cell=3)
    u=ul;u(k)=u(k)-h;call sim%face_pair(u,ur,-.0349_wp,.04_wp,fmm,fq,st,left_cell=3)
    err=max(err,maxval(abs(jl(:,k)-(fpp-fmm)/(2*h))/max(1.0_wp,abs((fpp-fmm)/(2*h)))))
    h=3e-5_wp*max(abs(ur(k)),merge(1e-5_wp,1e-2_wp,k==2))
    u=ur;u(k)=u(k)+h;call sim%face_pair(ul,u,-.0349_wp,.04_wp,fpp,fq,st,left_cell=3)
    u=ur;u(k)=u(k)-h;call sim%face_pair(ul,u,-.0349_wp,.04_wp,fmm,fq,st,left_cell=3)
    err=max(err,maxval(abs(jr(:,k)-(fpp-fmm)/(2*h))/max(1.0_wp,abs((fpp-fmm)/(2*h)))))
  enddo
  write(*,'(a,es11.3)')'  T1 face Jacobian vs independent central differences: ',err
  call chk(err<1e-4_wp,'T1 face Jacobian consistent with the face flux')
  if(fails>0)then
    write(*,'(a,i0,a,i0)')'FAILED ',fails,' of ',n;error stop 1
  endif
  write(*,'(a,i0,a)')'PASS junction T1: ',n,' assertions'
contains
  subroutine base()
    cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1
    cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
    cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp;cfg%mul=.5_wp;cfg%mug=.001_wp
    cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%Lr=12.546253993197219_wp
    cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
    cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%nCells=6;cfg%num_dummy_cells=0;cfg%N_pipe=3
    cfg%junction_model='MATCHED';cfg%quad_points=3
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
