module module_simulate_pipe_riser
  !! Conservative finite volumes, backward Euler and a pivoted band Newton solve.
  !! Every physical cell is an unknown. Phase mass is never clipped or reset.
  !! A matched junction has one phase flux and two momentum traces: its point
  !! force permits the pressure/phase-flux continuity prescribed by Burr.
  use mod_kinds, only: wp,i4
  use module_config, only: sim_config_t,first_physical_cell,last_physical_cell
  use module_find, only: compute_primitive_variables
  use module_flow_regime, only: burr15_result_t,burr15_stability,BURR15_NOT_EVALUATED
  use module_jacobian, only: DFduj
  use module_roe_solver, only: roe_solver_t,roe_absA_from_A
  use module_time_step, only: time_integrator_t
  use module_solver_utils, only: fluid_dynamics_t
  use module_boundary_conditions, only: boundary_flux,apply_ghost_cells
  use module_steady_state, only: solve_alpha_from_superficials
  use mod_void_helpers, only: ffan,Dffan
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only: output_unit,int64
  implicit none
  private
  public :: pipeline_riser_simulator_t
  type :: pipeline_riser_simulator_t
    type(sim_config_t) :: cfg
    type(roe_solver_t) :: roe_solver
    type(fluid_dynamics_t) :: fd
    type(time_integrator_t) :: ti
    real(wp) :: alpha_eps=1e-12_wp,min_dt=1e-12_wp,max_dt=1e-3_wp
    real(wp) :: dt_growth=1.2_wp,dt_shrink=0.5_wp
    integer(i4) :: newton_fast=3,newton_slow=8,max_backtrack=18,maxit_root=100
    real(wp),allocatable :: ul_prev(:),ug_prev(:)
    integer(i4),allocatable :: idx_prev(:)
    real(wp),allocatable :: lagged_dissipation(:,:,:)
  contains
    procedure :: init,ensure_prev_arrays,update_prev_arrays,process_cell,simulate
    procedure :: point_state,face_pair,assemble
    procedure :: prepare_dissipation
  end type
contains
  subroutine prepare_dissipation(this,U,status)
    !! Lag only the artificial-viscosity coefficients, consistently in the
    !! residual AND its derivative. Both physical fluxes remain implicit.
    !! This avoids a discontinuous numerical flux when a C0 empirical law's
    !! Jacobian crosses a Reynolds branch during Newton iterations.
    class(pipeline_riser_simulator_t),intent(inout)::this
    real(wp),intent(in)::U(:,:)
    integer,intent(out)::status
    real(wp)::fm(3),fp(3),dmat(3,3)
    integer::i,lo,hi
    status=0
    if(trim(this%cfg%dissipation_mode)/='LAGGED')then
      if(allocated(this%lagged_dissipation))deallocate(this%lagged_dissipation)
      return
    endif
    lo=first_physical_cell(this%cfg);hi=last_physical_cell(this%cfg)
    if(.not.allocated(this%lagged_dissipation))allocate(this%lagged_dissipation(3,3,size(U,1)))
    this%lagged_dissipation=0
    do i=lo,hi-1
      call this%face_pair(U(i,:),U(i+1,:),this%cfg%thetav(i),this%cfg%thetav(i+1),fm,fp,status,used_dissipation=dmat,left_cell=i)
      if(status/=0)return
      this%lagged_dissipation(:,:,i)=dmat
    enddo
  end subroutine
  subroutine init(this,cfg)
    class(pipeline_riser_simulator_t),intent(inout)::this
    type(sim_config_t),intent(in)::cfg
    this%cfg=cfg
    this%min_dt=cfg%dt_min;this%max_dt=cfg%dt_max
    this%maxit_root=cfg%maxit
    call this%ensure_prev_arrays(cfg%nCells)
  end subroutine

  subroutine ensure_prev_arrays(this,n_cells)
    class(pipeline_riser_simulator_t),intent(inout)::this
    integer(i4),intent(in)::n_cells
    if(allocated(this%ul_prev))then
      if(size(this%ul_prev)/=n_cells)deallocate(this%ul_prev,this%ug_prev,this%idx_prev)
    endif
    if(.not.allocated(this%ul_prev))then
      allocate(this%ul_prev(n_cells),this%ug_prev(n_cells),this%idx_prev(n_cells))
      this%ul_prev=0;this%ug_prev=0;this%idx_prev=0
    endif
  end subroutine

  subroutine update_prev_arrays(this,U,thetav)
    class(pipeline_riser_simulator_t),intent(inout)::this
    real(wp),intent(in)::U(:,:),thetav(:)
    real(wp)::v(6),f(3),s(3),a(3,3),ds(3,3)
    integer::i,code
    call this%ensure_prev_arrays(size(U,1))
    do i=1,size(U,1)
      call this%point_state(U(i,:),thetav(i),v,f,s,a,ds,code,.false.,cell_closure(this,thetav(i),i))
      if(code/=0)cycle
      this%ul_prev(i)=v(5);this%ug_prev(i)=v(6)
      this%idx_prev(i)=merge(0,1,cell_closure(this,thetav(i),i)==1)
    enddo
  end subroutine

  integer function cell_closure(this,theta,i) result(cid)
    class(pipeline_riser_simulator_t), intent(in) :: this
    real(wp), intent(in) :: theta
    integer, intent(in), optional :: i
    logical :: in_pipe
    ! The supplied mesh section resolves the shared horizontal angle.
    ! Retain the signed-angle convention for nonhorizontal operators.
    in_pipe=theta<0.0_wp
    if(theta==0.0_wp.and.present(i).and.this%cfg%N_pipe>0)in_pipe=i<=this%cfg%num_dummy_cells+this%cfg%N_pipe
    cid=this%cfg%closure_riser
    if(in_pipe)cid=this%cfg%closure_pipe
  end function cell_closure

  subroutine point_state(this,U,theta,v,F,S,A,dS,status,need_jac,closure_id)
    class(pipeline_riser_simulator_t),intent(in)::this
    real(wp),intent(in)::U(3),theta
    real(wp),intent(out)::v(6),F(3),S(3),A(3,3),dS(3,3)
    integer,intent(out)::status
    logical,intent(in),optional::need_jac
    integer,intent(in),optional::closure_id
    real(wp)::alpha,rl,rg,p,ul,ug,rhom,j,mu,re,ff,fr,fe,pip,grav,drag
    real(wp)::da(3),dp(3),dul(3),dug(3),dj(3),dr(3),dre(3),dm(3),den
    integer::cid,idx,code
    logical::want
    v=ieee_value(0.0_wp,ieee_quiet_nan);F=v(:3);S=F;A=0;dS=0
    status=1;want=.true.;if(present(need_jac))want=need_jac
    cid=cell_closure(this,theta)
    if(present(closure_id))cid=closure_id
    call compute_primitive_variables(U(1),U(2),U(3),theta,this%cfg%nCg,this%cfg%nCl, &
       this%cfg%nrho_l0,this%cfg%nP_l0,this%cfg%D_H,this%cfg%AREA,this%cfg%eps,this%cfg%G, &
       this%cfg%mu_l,this%cfg%mu_g,this%cfg%sigma,this%cfg%njl,this%cfg%omega_u,this%cfg%omega_rho, &
       this%cfg%tol,this%cfg%tola,this%maxit_root,alpha,rl,rg,p,ul,ug,idx,status=code,closure_id=cid)
    if(code/=0)then
      status=10+code;return
    endif
    v=[alpha,rl,rg,p,ul,ug]
    pip=this%cfg%omega_P/(this%cfg%omega_rho*this%cfg%omega_u**2)
    F=[U(1)*ul,U(2)*ug,U(1)*ul**2+U(2)*ug**2+pip*p]
    rhom=U(1)+U(2);j=(1-alpha)*ul+alpha*ug
    mu=this%cfg%mu_l*(1-alpha)+this%cfg%mu_g*alpha
    re=rhom*this%cfg%omega_rho*abs(j)*this%cfg%omega_u*this%cfg%D_H/mu
    grav=this%cfg%G*this%cfg%Lr/this%cfg%omega_u**2
    drag=2*this%cfg%Lr/this%cfg%D_H
    S=0
    if(re<2000.0_wp)then
      ! Analytic laminar extension through zero velocity, with finite derivative.
      S(3)=-rhom*grav*sin(theta)-32*this%cfg%Lr*mu*j/(this%cfg%omega_rho*this%cfg%omega_u*this%cfg%D_H**2)
    else
      ff=ffan(this%cfg%eps/this%cfg%D_H,re,code)
      if(code/=0)then
        status=20+code;return
      endif
      S(3)=-rhom*grav*sin(theta)-drag*rhom*ff*j*abs(j)
    endif
    if(want)then
      A=DFduj(alpha,rl,rg,ul,ug,this%cfg%nCl,this%cfg%nCg,theta,this%cfg%D_H,this%cfg%AREA, &
           this%cfg%eps,this%cfg%G,this%cfg%mu_l,this%cfg%mu_g,this%cfg%sigma,this%cfg%omega_u, &
           this%cfg%omega_rho,this%cfg%omega_P,this%cfg%tol,status=code,closure_id=cid)
      if(code/=0)then
        status=30+code;return
      endif
      den=rl*alpha/this%cfg%nCg**2+rg*(1-alpha)/this%cfg%nCl**2
      da=[-alpha/this%cfg%nCg**2,(1-alpha)/this%cfg%nCl**2,0.0_wp]/den
      dp=[rg,rl,0.0_wp]/den
      dul=A(1,:);dul(1)=dul(1)-ul;dul=dul/U(1)
      dug=A(2,:);dug(2)=dug(2)-ug;dug=dug/U(2)
      dj=(ug-ul)*da+(1-alpha)*dul+alpha*dug
      dr=[1.0_wp,1.0_wp,0.0_wp];dm=(this%cfg%mu_g-this%cfg%mu_l)*da
      if(re<2000.0_wp)then
        dS(3,:)=-grav*sin(theta)*dr-32*this%cfg%Lr/(this%cfg%omega_rho*this%cfg%omega_u*this%cfg%D_H**2)*(mu*dj+j*dm)
      else
        call Dffan(re,this%cfg%eps/this%cfg%D_H,fr,fe,code)
        if(code/=0)then
          status=40+code;return
        endif
        dre=re*(dr/rhom+dj/j-dm/mu)
        dS(3,:)=-grav*sin(theta)*dr-drag*(ff*j*abs(j)*dr+rhom*fr*j*abs(j)*dre+2*rhom*ff*abs(j)*dj)
      endif
    endif
    if(.not.all(ieee_is_finite(F)).or..not.all(ieee_is_finite(S)))return
    status=0
  end subroutine

  subroutine face_pair(this,UL,UR,thetaL,thetaR,FM,FP,status,JL,JR,KL,KR,speed,fd_one_sided,blocking_active,junction_trace,dissipation_matrix,used_dissipation,left_cell)
    use module_lapack_small,only:dgeev_3x3
    use module_matrix_roe,only:calculate_roe_matrix
    class(pipeline_riser_simulator_t),intent(in)::this
    real(wp),intent(in)::UL(3),UR(3),thetaL,thetaR
    real(wp),intent(out)::FM(3),FP(3)
    integer,intent(out)::status
    integer,intent(out),optional::fd_one_sided
    integer,intent(in),optional::left_cell
    logical,intent(out),optional::blocking_active
    real(wp),intent(out),optional::junction_trace(6)
    real(wp),intent(in),optional::dissipation_matrix(3,3)
    real(wp),intent(out),optional::used_dissipation(3,3)
    real(wp),intent(out),optional::JL(3,3),JR(3,3),KL(3,3),KR(3,3),speed
    real(wp)::vl(6),vr(6),fl(3),fr(3),sl(3),sr(3),al(3,3),ar(3,3),ds(3,3),ab(3,3),aa(3,3),amax
    real(wp)::up(3),um(3),fm1(3),fp1(3),fm2(3),fp2(3),h,lj(3,3),rj(3,3),lk(3,3),rk(3,3)
    integer::code,k,cp,cm,n_one_sided,id_left,id_right
    logical::jac,match,blocking_mode,base_blocked
    real(wp)::base_trace(6),override_speed,eig_real(3),eig_imag(3),eig_vec(3,3)
    if(present(used_dissipation))used_dissipation=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(dissipation_matrix))then
      ! The caller may freeze numerical viscosity at the beginning of a time
      ! step. This SAME matrix is then used in every base/FD/residual call;
      ! only freezing the Jacobian would leave a different nonlinear equation.
      status=72
      if(.not.all(ieee_is_finite(dissipation_matrix)))return
      call dgeev_3x3(dissipation_matrix,eig_real,eig_imag,eig_vec,code)
      if(code/=0)return
      override_speed=maxval(hypot(eig_real,eig_imag))
    endif
    n_one_sided=0
    if(present(fd_one_sided))fd_one_sided=0
    jac=present(JL).or.present(JR).or.present(KL).or.present(KR)
    blocking_mode=trim(this%cfg%junction_model)=='BLOCKING'
    id_left=cell_closure(this,thetaL);id_right=cell_closure(this,thetaR)
    match=thetaL<0.and.thetaR>=0
    if(present(left_cell).and.this%cfg%N_pipe>0)then
      id_left=cell_closure(this,thetaL,left_cell)
      id_right=cell_closure(this,thetaR,left_cell+1)
      if(thetaL==0.0_wp.and.thetaR>=0.0_wp)match=left_cell==this%cfg%num_dummy_cells+this%cfg%N_pipe
    endif
    match=match.and.(trim(this%cfg%junction_model)=='MATCHED'.or.blocking_mode)
    if(present(blocking_active))blocking_active=.false.
    if(present(junction_trace))junction_trace=ieee_value(0.0_wp,ieee_quiet_nan)
    call basic(UL,UR,FM,FP,status,al,ar,ab,amax,base_blocked,base_trace)
    if(status/=0)return
    if(present(speed))speed=amax
    if(present(used_dissipation))used_dissipation=ab
    if(present(blocking_active))blocking_active=base_blocked
    if(present(junction_trace))junction_trace=base_trace
    if(.not.jac)return
    lj=0.5_wp*(al+ab);rj=0.5_wp*(ar-ab);lk=lj;rk=rj
    if(this%cfg%jac_dissipation_fd.or.match)then
      do k=1,3
        call differentiate(1,k,lj(:,k),lk(:,k),code)
        if(code/=0)then
          status=70;return
        endif
        call differentiate(2,k,rj(:,k),rk(:,k),code)
        if(code/=0)then
          status=71;return
        endif
      enddo
    endif
    if(present(fd_one_sided))fd_one_sided=n_one_sided
    if(present(JL))JL=lj
    if(present(JR))JR=rj
    if(present(KL))KL=lk
    if(present(KR))KR=rk
  contains
    subroutine differentiate(side,component,dm,dp,istat)
      !! Differentiate the SAME admissible face operator. Near a model-domain
      !! boundary a centered perturbation can fail although the base is valid.
      !! Use a second-order one-sided difference then, with two admissible
      !! samples. If neither side supports that stencil, reduce h; no state,
      !! flux, pressure or phase mass is clipped to manufacture a derivative.
      integer,intent(in)::side,component
      real(wp),intent(out)::dm(3),dp(3)
      integer,intent(out)::istat
      real(wp)::xp(3),xm(3),yp(3),ym(3),fmp(3),fpp(3),fmm(3),fpm(3),fm2(3),fp2(3)
      real(wp)::qa(3,3),qb(3,3),qd(3,3),wave,step,base(3),step_min
      integer::ip,im,i2,it
      base=UL
      if(side==2)base=UR
      step=this%cfg%jac_dissipation_fd_eps*max(abs(base(component)),merge(1e-5_wp,1e-2_wp,component==2))
      if(component<=2)step=min(step,0.2_wp*base(component))
      step_min=32.0_wp*epsilon(1.0_wp)*max(abs(base(component)),1e-8_wp)
      istat=1
      do it=1,10
        if(step<=step_min)return
        xp=UL;xm=UL;yp=UR;ym=UR
        if(side==1)then
          xp(component)=UL(component)+step;xm(component)=UL(component)-step
        else
          yp(component)=UR(component)+step;ym(component)=UR(component)-step
        endif
        call basic(xp,yp,fmp,fpp,ip,qa,qb,qd,wave)
        call basic(xm,ym,fmm,fpm,im,qa,qb,qd,wave)
        if(ip==0.and.im==0)then
          dm=(fmp-fmm)/(2*step);dp=(fpp-fpm)/(2*step)
          istat=0;return
        else if(ip==0)then
          if(side==1)xp(component)=UL(component)+2*step
          if(side==2)yp(component)=UR(component)+2*step
          call basic(xp,yp,fm2,fp2,i2,qa,qb,qd,wave)
          if(i2==0)then
            dm=(-3*FM+4*fmp-fm2)/(2*step);dp=(-3*FP+4*fpp-fp2)/(2*step)
            n_one_sided=n_one_sided+1;istat=0;return
          endif
        else if(im==0)then
          if(side==1)xm(component)=UL(component)-2*step
          if(side==2)ym(component)=UR(component)-2*step
          call basic(xm,ym,fm2,fp2,i2,qa,qb,qd,wave)
          if(i2==0)then
            dm=(3*FM-4*fmm+fm2)/(2*step);dp=(3*FP-4*fpm+fp2)/(2*step)
            n_one_sided=n_one_sided+1;istat=0;return
          endif
        endif
        step=0.5_wp*step
      enddo
    end subroutine differentiate
    subroutine momentum_trace(gl,gg,rl,rg,alpha,pressure,momentum,istat)
      real(wp),intent(in)::gl,gg,rl,rg,alpha,pressure
      real(wp),intent(out)::momentum
      integer,intent(out)::istat
      momentum=ieee_value(0.0_wp,ieee_quiet_nan);istat=1
      if(.not.all(ieee_is_finite([gl,gg,rl,rg,alpha,pressure])))return
      if(min(rl,rg)<=0.or.alpha<0.or.alpha>1)return
      momentum=pressure
      if(gl/=0.0_wp)then
        if(alpha==1.0_wp)return
        momentum=momentum+gl*gl/((1-alpha)*rl)
      endif
      if(gg/=0.0_wp)then
        if(alpha==0.0_wp)return
        momentum=momentum+gg*gg/(alpha*rg)
      endif
      if(.not.ieee_is_finite(momentum))return
      istat=0
    end subroutine momentum_trace
    subroutine basic(x,y,fm_,fp_,istat,amat,bmat,dmat,sp,was_blocked,trace)
      real(wp),intent(in)::x(3),y(3)
      real(wp),intent(out)::fm_(3),fp_(3),amat(3,3),bmat(3,3),dmat(3,3),sp
      integer,intent(out)::istat
      logical,intent(out),optional::was_blocked
      real(wp),intent(out),optional::trace(6)
      real(wp)::mapped(3),vm(6),fs(3),ss(3),scratch(3,3),a0(3,3),pstar,rls,rgs,astl,astr
      real(wp)::pip,kin,jsl,jsg,ax(3,3),raw_gas_flux,path_theta,path_theta_right
      integer::ic,q
      logical::active
      active=.false.
      if(present(was_blocked))was_blocked=.false.
      if(present(trace))trace=ieee_value(0.0_wp,ieee_quiet_nan)
      call this%point_state(x,thetaL,vl,fl,sl,amat,scratch,istat,closure_id=id_left)
      if(istat/=0)return
      call this%point_state(y,thetaR,vr,fr,sr,bmat,scratch,istat,closure_id=id_right)
      if(istat/=0)return
      if(match)then
        ! Map the right trace to the pipe closure with the SAME p and phase fluxes.
        astl=solve_alpha_from_superficials(fr(1)/vr(2),fr(2)/vr(3),vr(2),vr(3),thetaL, &
          this%cfg%D_H,this%cfg%AREA,this%cfg%eps,this%cfg%G,this%cfg%mu_l,this%cfg%mu_g,this%cfg%sigma, &
          this%cfg%omega_u,this%cfg%omega_rho,this%cfg%tol,status=ic,closure_id=this%cfg%closure_pipe,alpha_hint=vl(1))
        if(ic/=0)then
          istat=800+ic;return
        endif
        mapped=[(1-astl)*vr(2),astl*vr(3),fr(1)+fr(2)]
        call this%point_state(mapped,thetaL,vm,fs,ss,bmat,scratch,istat,closure_id=id_left)
        if(istat/=0)return
        a0=0.5_wp*(amat+bmat)
      else
        mapped=y;fs=fr;a0=0.5_wp*(amat+bmat)
      endif
      if(present(dissipation_matrix))then
        dmat=dissipation_matrix;sp=override_speed
      else
        if(match.or.id_left==id_right)then
          path_theta=thetaL;path_theta_right=thetaR
          if(match)path_theta_right=thetaL
          ! Interpolate geometry along the conservative path on one closure.
          ! MATCHED already maps both traces to the same pipe inclination.
          ! This integral of F_U is a dissipative approximation: endpoint
          ! physical fluxes retain their angles; no exact A*dU secant is claimed.
          ! An override skips this work entirely: lagged viscosity is built
          ! once at the old state, then shared by residual and Jacobian.
          a0=calculate_roe_matrix(x,mapped,this%cfg%nCl,this%cfg%nCg,path_theta, &
             this%cfg%nrho_l0,this%cfg%nP_l0,this%cfg%D_H,this%cfg%AREA,this%cfg%eps,this%cfg%G, &
             this%cfg%mu_l,this%cfg%mu_g,this%cfg%sigma,this%cfg%njl,this%cfg%omega_u, &
             this%cfg%omega_rho,this%cfg%omega_P,this%cfg%tol,this%cfg%tola,this%cfg%quad_points, &
             ul_prev=vl(5),ug_prev=vl(6),status=ic,closure_id=id_left,theta_right=path_theta_right)
          if(ic/=0)then
            istat=1200+ic;return
          endif
        endif
        call roe_absA_from_A(a0,dmat,sp,istat)
        if(istat/=0)return
        if(trim(this%cfg%flux_scheme)=='RUSANOV')then
          dmat=0
          do q=1,3
            dmat(q,q)=sp
          enddo
        endif
      endif
      fm_=0.5_wp*(fl+fs)-0.5_wp*matmul(dmat,mapped-x)
      fp_=fm_
      if(match)then
        pip=this%cfg%omega_P/(this%cfg%omega_rho*this%cfg%omega_u**2)
        kin=0.5_wp*(fl(3)-pip*vl(4)+fs(3)-pip*vm(4))
        pstar=(fm_(3)-kin)/pip
        rls=this%cfg%nrho_l0+(pstar-this%cfg%nP_l0)/this%cfg%nCl**2
        rgs=pstar/this%cfg%nCg**2
        if(min(pstar,rls,rgs)<=0)then
          istat=89;return
        endif
        jsl=fm_(1)/rls;jsg=fm_(2)/rgs
        raw_gas_flux=fm_(2)
        astr=solve_alpha_from_superficials(jsl,jsg,rls,rgs,0.0_wp,this%cfg%D_H,this%cfg%AREA, &
          this%cfg%eps,this%cfg%G,this%cfg%mu_l,this%cfg%mu_g,this%cfg%sigma,this%cfg%omega_u, &
          this%cfg%omega_rho,this%cfg%tol,status=ic,closure_id=this%cfg%closure_riser)
        if(ic/=0)then
          ! EXPERIMENTAL unilateral gas-barrier model, opt-in only. The
          ! projection keeps the Roe pressure and liquid mass flux, and sets
          ! one COMMON gas flux to zero only when the upward liquid trace has
          ! no admissible drift state for the proposed reverse gas flux.
          ! This preserves phase balances and is continuous at activation; it
          ! is NOT a characteristic Riemann solution or the reference model's
          ! moving liquid-front equation. Admissible counterflow is retained.
          if(blocking_mode.and.ic==21.and.fm_(1)>0.0_wp.and.fm_(2)<0.0_wp.and. &
             (this%cfg%closure_riser==2.or.this%cfg%closure_riser==4))then
            fm_(2)=0.0_wp;fp_(2)=0.0_wp;jsg=0.0_wp;astr=0.0_wp
            active=.true.
          else
            istat=1000+ic;return
          endif
        endif
        astl=solve_alpha_from_superficials(jsl,jsg,rls,rgs,thetaL,this%cfg%D_H,this%cfg%AREA, &
          this%cfg%eps,this%cfg%G,this%cfg%mu_l,this%cfg%mu_g,this%cfg%sigma,this%cfg%omega_u, &
          this%cfg%omega_rho,this%cfg%tol,status=ic,closure_id=this%cfg%closure_pipe,alpha_hint=vl(1))
        if(ic/=0)then
          istat=900+ic;return
        endif
        call momentum_trace(fm_(1),fm_(2),rls,rgs,astl,pip*pstar,fm_(3),ic)
        if(ic/=0)then
          istat=1100+ic;return
        endif
        call momentum_trace(fm_(1),fm_(2),rls,rgs,astr,pip*pstar,fp_(3),ic)
        if(ic/=0)then
          istat=1100+ic;return
        endif
        if(present(was_blocked))was_blocked=active
        ! P*, Gl*, Gg*, alpha_pipe*, alpha_riser*, unprojected Roe Gg*.
        if(present(trace))trace=[pstar,fm_(1),fm_(2),astl,astr,raw_gas_flux]
      endif
      istat=0
    end subroutine
  end subroutine

  subroutine process_cell(this,i,U,delta_s,thetav,CFL,flux,source_avg,dt,dflux_dUL,dflux_dUR,dSrc_dU,status)
    class(pipeline_riser_simulator_t),intent(inout)::this
    integer(i4),intent(in)::i
    real(wp),intent(in)::U(:,:),delta_s(:),thetav(:),CFL
    real(wp),intent(out)::flux(3),source_avg(3),dt
    real(wp),intent(out),optional::dflux_dUL(3,3),dflux_dUR(3,3),dSrc_dU(3,3)
    integer,intent(out),optional::status
    real(wp)::fp(3),v(6),f(3),a(3,3),ds(3,3),jl(3,3),jr(3,3),sp
    integer::code
    call this%face_pair(U(i,:),U(i+1,:),thetav(i),thetav(i+1),flux,fp,code,jl,jr,speed=sp,left_cell=i)
    if(code==0)then
      call this%point_state(U(i,:),thetav(i),v,f,source_avg,a,ds,code,closure_id=cell_closure(this,thetav(i),i))
      dt=CFL*delta_s(i)/max(sp,tiny(1.0_wp))
      if(present(dflux_dUL))dflux_dUL=jl
      if(present(dflux_dUR))dflux_dUR=jr
      if(present(dSrc_dU))dSrc_dU=ds
    else
      dt=0
    endif
    if(present(status))status=code
  end subroutine

  subroutine assemble(this,U,Uold,dt,R,A,B,C,norm,fin,fout,status,with_jac,base_mass_flux)
    class(pipeline_riser_simulator_t),intent(in)::this
    real(wp),intent(in)::U(:,:),Uold(:,:),dt
    real(wp),intent(out)::R(:,:),A(:,:,:),B(:,:,:),C(:,:,:),norm,fin(3),fout(3)
    integer,intent(out)::status
    logical,intent(in)::with_jac
    ! Phase fluxes actually inserted in the residual: column 1 leaves the
    ! pipe, column 2 enters the riser. Keep both; never impose equality here.
    real(wp),intent(out),optional::base_mass_flux(2,2)
    integer::lo,hi,n,i,k,l,q,code,cp,cm
    real(wp)::v(6),f(3),s(3),fj(3,3),sj(3,3),fm(3),fp(3),jl(3,3),jr(3,3),kl(3,3),kr(3,3)
    real(wp)::up(3),um(3),plus(3),minus(3),h,jb(3,3),scale(3),ds
    lo=first_physical_cell(this%cfg);hi=last_physical_cell(this%cfg);n=hi-lo+1
    R=0;A=0;B=0;C=0;status=0;norm=huge(1.0_wp)
    if(present(base_mass_flux))base_mass_flux=ieee_value(0.0_wp,ieee_quiet_nan)
    do i=lo,hi
      l=i-lo+1;ds=this%cfg%cell_width(i)
      call this%point_state(U(i,:),this%cfg%thetav(i),v,f,s,fj,sj,code,with_jac,cell_closure(this,this%cfg%thetav(i),i))
      if(code/=0)then
        status=1000+code;return
      endif
      R(:,l)=U(i,:)-Uold(i,:)-dt*s
      if(with_jac)then
        B(:,:,l)=-dt*sj
        do k=1,3
          B(k,k,l)=B(k,k,l)+1
        enddo
      endif
    enddo
    do i=lo,hi-1
      l=i-lo+1
      if(with_jac.and.allocated(this%lagged_dissipation))then
        call this%face_pair(U(i,:),U(i+1,:),this%cfg%thetav(i),this%cfg%thetav(i+1),fm,fp,code,jl,jr,kl,kr, &
            dissipation_matrix=this%lagged_dissipation(:,:,i),left_cell=i)
      else if(with_jac)then
        call this%face_pair(U(i,:),U(i+1,:),this%cfg%thetav(i),this%cfg%thetav(i+1),fm,fp,code,jl,jr,kl,kr,left_cell=i)
      else if(allocated(this%lagged_dissipation))then
        call this%face_pair(U(i,:),U(i+1,:),this%cfg%thetav(i),this%cfg%thetav(i+1),fm,fp,code, &
            dissipation_matrix=this%lagged_dissipation(:,:,i),left_cell=i)
      else
        call this%face_pair(U(i,:),U(i+1,:),this%cfg%thetav(i),this%cfg%thetav(i+1),fm,fp,code,left_cell=i)
      endif
      if(code/=0)then
        status=2000+code;return
      endif
      if(present(base_mass_flux))then
        if(i==lo+this%cfg%N_pipe-1)then
          base_mass_flux(:,1)=fm(:2);base_mass_flux(:,2)=fp(:2)
        endif
      endif
      R(:,l)=R(:,l)+dt*fm/this%cfg%cell_width(i)
      R(:,l+1)=R(:,l+1)-dt*fp/this%cfg%cell_width(i+1)
      if(with_jac)then
        B(:,:,l)=B(:,:,l)+dt*jl/this%cfg%cell_width(i)
        C(:,:,l)=dt*jr/this%cfg%cell_width(i)
        A(:,:,l+1)=-dt*kl/this%cfg%cell_width(i+1)
        B(:,:,l+1)=B(:,:,l+1)-dt*kr/this%cfg%cell_width(i+1)
      endif
    enddo
    do q=1,2
      i=lo;l=1
      if(q==2)then
        i=hi;l=n
      endif
      call boundary_flux(U(i,:),this%cfg%thetav(i),this%cfg,q==1,f,code)
      if(code/=0)then
        status=3000+code;return
      endif
      if(q==1)then
        fin=f
      else
        fout=f
      endif
      ds=this%cfg%cell_width(i)
      if(q==1)ds=-ds
      R(:,l)=R(:,l)+dt*f/ds
      if(with_jac)then
        do k=1,3
          h=this%cfg%jac_dissipation_fd_eps*max(abs(U(i,k)),merge(1e-5_wp,1e-2_wp,k==2))
          if(k<=2)h=min(h,0.2_wp*U(i,k))
          up=U(i,:);um=up;up(k)=up(k)+h;um(k)=um(k)-h
          call boundary_flux(up,this%cfg%thetav(i),this%cfg,q==1,plus,cp)
          call boundary_flux(um,this%cfg%thetav(i),this%cfg,q==1,minus,cm)
          if(cp/=0.or.cm/=0)then
            status=3090;return
          endif
          jb(:,k)=(plus-minus)/(2*h)
        enddo
        B(:,:,l)=B(:,:,l)+dt*jb/ds
      endif
    enddo
    norm=0
    do i=lo,hi
      l=i-lo+1
      scale=[max(Uold(i,1),1e-3_wp),max(Uold(i,2),1e-7_wp),max(abs(Uold(i,3)),1.0_wp)]
      norm=max(norm,maxval(abs(R(:,l))/scale))
    enddo
    if(.not.ieee_is_finite(norm))status=4000
  end subroutine

  subroutine simulate(this,U,vns,thetav,T,max_iter,it_target,use_it_target,it_p,dt_fac_min,dt_fac_max, &
       dt_growth,dt_shrink,newton_fast,newton_slow,status)
    class(pipeline_riser_simulator_t),intent(inout)::this
    real(wp),intent(inout)::U(:,:)
    real(wp),intent(in)::vns(:),thetav(:),T
    integer(i4),intent(in),optional::max_iter,it_target,newton_fast,newton_slow
    logical,intent(in),optional::use_it_target
    real(wp),intent(in),optional::it_p,dt_fac_min,dt_fac_max,dt_growth,dt_shrink
    integer,intent(out),optional::status
    real(wp),allocatable::old(:,:),trial(:,:),work(:,:),r(:,:),rt(:,:),a(:,:,:),b(:,:,:),c(:,:,:),dx(:,:)
    real(wp),allocatable::az(:,:,:),bz(:,:,:),cz(:,:,:)
    real(wp)::tcur,dt,norm,nt,lambda,fin(3),fout(3),fi2(3),fo2(3),inventory0(2),inventory(2),in_sum(2),out_sum(2)
    real(wp)::next_output,next_field,last_print,wall0,wall1,remaining,dt_next
    integer::lo,hi,n,step,it,code,ls,info,maxsteps,rej,unitmon,unitfield,unitstat,i,k,lastcode,finstatus
    logical::accepted,converged,requested_full_jac,diagnose_burr15
    integer::unitburr15,unitsteps,unitmass,ib,initial_flux_status
    integer::jacobian_evaluations,residual_evaluations,linear_solves,backtrack_trials,newton_iterations,prepare_calls
    integer(int64)::clock_start,clock_now,clock_rate
    real(wp)::base_mass_flux(2,2),trial_base_mass_flux(2,2),accepted_base(2,2),accepted_in(2),accepted_out(2)
    real(wp)::base_sum(2,2),base_abs(2,2),in_abs(2),out_abs(2),massflow_scale,mass_scale
    real(wp)::mass0(2,3),mass_previous(2,3),mass_current(2,3),mass_residual(2,3),mass_residual_sum(2,3)
    real(wp)::step_balance(2,3),cpu_previous,cpu_now,elapsed_wall
    lo=first_physical_cell(this%cfg);hi=last_physical_cell(this%cfg);n=hi-lo+1
    if(present(status))status=1
    requested_full_jac=this%cfg%jac_dissipation_fd
    if(this%cfg%use_muscl)then
      write(*,'(a)') 'FV1 verified path requires use_muscl=.false.; MUSCL transient is not qualified.'
      return
    endif
    if(.not.allocated(this%cfg%cell_width))then
      write(*,'(a)') 'Exact finite-volume cell widths are required.'
      return
    endif
    allocate(old(size(U,1),3),trial(size(U,1),3),work(size(U,1),3),r(3,n),rt(3,n),a(3,3,n),b(3,3,n),c(3,3,n),dx(3,n), &
      az(3,3,n),bz(3,3,n),cz(3,3,n))
    maxsteps=10000000;if(present(max_iter))maxsteps=max_iter
    if(present(dt_growth))this%dt_growth=dt_growth
    if(present(dt_shrink))this%dt_shrink=dt_shrink
    dt=min(this%cfg%dt_start,this%max_dt)
    if(dt<=0)dt=min(0.01_wp/this%cfg%omega_t,this%max_dt)
    tcur=0;step=0;rej=0;in_sum=0;out_sum=0;lastcode=0;finstatus=0
    next_output=0;next_field=0;last_print=-huge(1.0_wp)
    ib=lo+this%cfg%N_pipe-1
    base_sum=0;base_abs=0;in_abs=0;out_abs=0;mass_residual=0;mass_residual_sum=0;step_balance=0
    accepted_in=ieee_value(0.0_wp,ieee_quiet_nan);accepted_out=accepted_in
    accepted_base=ieee_value(0.0_wp,ieee_quiet_nan);initial_flux_status=0
    massflow_scale=this%cfg%omega_rho*this%cfg%omega_u*this%cfg%AREA
    mass_scale=this%cfg%omega_rho*this%cfg%AREA*this%cfg%Lr
    jacobian_evaluations=0;residual_evaluations=0;linear_solves=0;backtrack_trials=0;newton_iterations=0;prepare_calls=0
    call cpu_time(wall0);cpu_previous=wall0
    call system_clock(clock_start,clock_rate)
    call mass_inventory(U,mass0);mass_previous=mass0
    call cpu_time(cpu_now)
    do k=1,2
      inventory0(k)=sum(U(lo:hi,k)*this%cfg%cell_width(lo:hi))*this%cfg%omega_rho*this%cfg%AREA*this%cfg%Lr
    enddo
    open(newunit=unitmon,file=trim(this%cfg%out_dir)//'/monitor.csv',status='replace')
    write(unitmon,'(a)') 'time_s,dt_s,p_in_Pa,p_pipe_end_Pa,p_riser_base_Pa,p_out_cell_Pa,p_separator_Pa,'// &
       'alpha_pipe_end,alpha_riser_base,alpha_out,jl_out_m_s,jg_out_m_s,mdot_l_out_kg_s,mdot_g_out_kg_s,'// &
       'Ml_kg,Mg_kg,Ml_pipe_kg,Mg_pipe_kg,Ml_riser_kg,Mg_riser_kg,'// &
       'Ml_in_integral_kg,Mg_in_integral_kg,Ml_out_integral_kg,Mg_out_integral_kg,balance_l_kg,balance_g_kg,newton_iterations,rejected_steps,'// &
       'mdot_l_base_kg_s,mdot_g_base_kg_s,p_junction_Pa,alpha_pipe_trace,alpha_riser_trace,mdot_g_unconstrained_kg_s,blocking_active,'// &
       'holdup_liquid_base,holdup_liquid_out,mdot_l_in_kg_s,mdot_g_in_kg_s,'// &
       'Ml_in_abs_integral_kg,Mg_in_abs_integral_kg,Ml_out_abs_integral_kg,Mg_out_abs_integral_kg,'// &
       'Ml_base_integral_kg,Mg_base_integral_kg,Ml_base_abs_integral_kg,Mg_base_abs_integral_kg,'// &
       'Ml_base_right_integral_kg,Mg_base_right_integral_kg,Ml_base_right_abs_integral_kg,Mg_base_right_abs_integral_kg'
    open(newunit=unitfield,file=trim(this%cfg%out_dir)//'/fields.csv',status='replace')
    write(unitfield,'(a)') 'time_s,s_m,theta_rad,alpha,holdup_liquid,p_Pa,ul_m_s,ug_m_s,rhol_kg_m3,rhog_kg_m3'
    diagnose_burr15=trim(this%cfg%pipe_regime_mode)=='AUTO_BURR15'
    if(diagnose_burr15)then
      open(newunit=unitburr15,file=trim(this%cfg%out_dir)//'/burr15.csv',status='replace')
      write(unitburr15,'(a)') 'time_s,s_m,closure_id,status,hL_over_D,ug2_m2_s2,ug2_critical_m2_s2,margin_m2_s2,unstable'
    endif
    open(newunit=unitsteps,file=trim(this%cfg%out_dir)//'/steps.csv',status='replace')
    write(unitsteps,'(a)') 'step,time_s,dt_s,newton_iterations,rejected_steps_total,residual_norm,'// &
      'newton_iterations_total,jacobian_evaluations_total,residual_evaluations_total,linear_solves_total,'// &
      'backtrack_trials_total,prepare_calls_total,cpu_total_s,cpu_since_previous_accept_s,wall_total_s'
    open(newunit=unitmass,file=trim(this%cfg%out_dir)//'/conservation.csv',status='replace')
    write(unitmass,'(a)') 'step,time_s,dt_s,residual_norm,flux_status,'// &
      'Ml_kg,Mg_kg,Ml_pipe_kg,Mg_pipe_kg,Ml_riser_kg,Mg_riser_kg,'// &
      'mdot_l_in_kg_s,mdot_g_in_kg_s,mdot_l_out_kg_s,mdot_g_out_kg_s,'// &
      'mdot_l_base_left_kg_s,mdot_g_base_left_kg_s,mdot_l_base_right_kg_s,mdot_g_base_right_kg_s,'// &
      'Ml_in_integral_kg,Mg_in_integral_kg,Ml_out_integral_kg,Mg_out_integral_kg,'// &
      'Ml_base_left_integral_kg,Mg_base_left_integral_kg,Ml_base_right_integral_kg,Mg_base_right_integral_kg,'// &
      'Ml_in_abs_integral_kg,Mg_in_abs_integral_kg,Ml_out_abs_integral_kg,Mg_out_abs_integral_kg,'// &
      'Ml_base_left_abs_integral_kg,Mg_base_left_abs_integral_kg,Ml_base_right_abs_integral_kg,Mg_base_right_abs_integral_kg,'// &
      'balance_l_kg,balance_g_kg,balance_l_pipe_kg,balance_g_pipe_kg,balance_l_riser_kg,balance_g_riser_kg,'// &
      'step_balance_l_kg,step_balance_g_kg,step_balance_l_pipe_kg,step_balance_g_pipe_kg,step_balance_l_riser_kg,step_balance_g_riser_kg,'// &
      'mass_residual_l_kg,mass_residual_g_kg,mass_residual_l_pipe_kg,mass_residual_g_pipe_kg,mass_residual_l_riser_kg,mass_residual_g_riser_kg,'// &
      'mass_residual_sum_l_kg,mass_residual_sum_g_kg,mass_residual_sum_l_pipe_kg,mass_residual_sum_g_pipe_kg,'// &
      'mass_residual_sum_l_riser_kg,mass_residual_sum_g_riser_kg'
    ! Initial fluxes are diagnostic only and do not contribute to an integral.
    ! An unavailable diagnostic must not alter the solver's acceptance path.
    call boundary_flux(U(lo,:),thetav(lo),this%cfg,.true.,fin,code)
    if(code==0)then
      accepted_in=fin(:2)
    else
      initial_flux_status=3000+code
    endif
    call boundary_flux(U(hi,:),thetav(hi),this%cfg,.false.,fout,code)
    if(code==0)accepted_out=fout(:2)
    if(code/=0)then
      initial_flux_status=3000+code
      finstatus=code
    else
      call write_monitor(0.0_wp,0)
      call write_field()
    endif
    call write_mass_record(0.0_wp,0.0_wp,initial_flux_status)
    do while(tcur<T.and.step<maxsteps.and.finstatus==0)
      remaining=T-tcur;dt=min(dt,remaining,this%max_dt)
      old=U;trial=U;accepted=.false.;converged=.false.
      prepare_calls=prepare_calls+1
      call this%prepare_dissipation(old,code)
      if(code/=0)then
        finstatus=8000+code
        exit
      endif
      do it=1,this%cfg%newton_max
        newton_iterations=newton_iterations+1
        ! A frozen |A| is an inexact Newton matrix. Upgrade to the complete
        ! flux derivative when its convergence slows; the residual never changes.
        if(it>=7)this%cfg%jac_dissipation_fd=.true.
        jacobian_evaluations=jacobian_evaluations+1;residual_evaluations=residual_evaluations+1
        call this%assemble(trial,old,dt,r,a,b,c,norm,fin,fout,code,.true.,base_mass_flux)
        if(code/=0)then
          lastcode=code;exit
        endif
        if(norm<=this%cfg%tol_res)then
          call capture_mass_residual(r)
          accepted=.true.;converged=.true.;exit
        endif
        linear_solves=linear_solves+1
        call this%ti%solve_block_tridiagonal_3x3(a,b,c,-r,dx,info)
        if(info/=0)then
          lastcode=5000+info;exit
        endif
        lambda=1
        do i=lo,hi
          do k=1,2
            if(dx(k,i-lo+1)<0)lambda=min(lambda,-0.90_wp*trial(i,k)/dx(k,i-lo+1))
          enddo
        enddo
        accepted=.false.
        do ls=1,this%max_backtrack
          work=trial
          do i=lo,hi
            work(i,:)=trial(i,:)+lambda*dx(:,i-lo+1)
          enddo
          residual_evaluations=residual_evaluations+1;backtrack_trials=backtrack_trials+1
          call this%assemble(work,old,dt,rt,az,bz,cz,nt,fi2,fo2,code,.false.,trial_base_mass_flux)
          if(code==0)then
            if(nt<norm*(1-1e-4_wp*lambda).or.nt<=this%cfg%tol_res)then
              accepted=.true.;exit
            endif
          else
            lastcode=code
          endif
          lambda=lambda*0.5_wp
        enddo
        if(.not.accepted)then
          if(.not.this%cfg%jac_dissipation_fd)then
            this%cfg%jac_dissipation_fd=.true.
            cycle
          endif
          lastcode=6000;exit
        endif
        trial=work
        if(nt<=this%cfg%tol_res)then
          norm=nt;fin=fi2;fout=fo2;base_mass_flux=trial_base_mass_flux
          call capture_mass_residual(rt)
          converged=.true.;exit
        endif
      enddo
      if(.not.converged)then
        call write_rejected_state()
        rej=rej+1;dt=dt*this%dt_shrink
        if(dt<this%min_dt)then
          write(*,'(a,es14.6,a,i0)') 'Time step exhausted at t(s)=',tcur*this%cfg%omega_t,' status=',lastcode
          finstatus=max(1,lastcode);exit
        endif
        cycle
      endif
      U=trial;call apply_ghost_cells(U,this%cfg%num_dummy_cells)
      this%cfg%jac_dissipation_fd=requested_full_jac
      tcur=tcur+dt;step=step+1
      accepted_in=fin(:2);accepted_out=fout(:2);accepted_base=base_mass_flux
      call mass_inventory(U,mass_current)
      step_balance(:,1)=mass_current(:,1)-mass_previous(:,1)- &
        dt*this%cfg%omega_t*massflow_scale*(accepted_in-accepted_out)
      step_balance(:,2)=mass_current(:,2)-mass_previous(:,2)- &
        dt*this%cfg%omega_t*massflow_scale*(accepted_in-accepted_base(:,1))
      step_balance(:,3)=mass_current(:,3)-mass_previous(:,3)- &
        dt*this%cfg%omega_t*massflow_scale*(accepted_base(:,2)-accepted_out)
      mass_previous=mass_current
      mass_residual_sum=mass_residual_sum+mass_residual
      base_sum=base_sum+dt*this%cfg%omega_t*massflow_scale*accepted_base
      base_abs=base_abs+dt*this%cfg%omega_t*massflow_scale*abs(accepted_base)
      in_abs=in_abs+dt*this%cfg%omega_t*massflow_scale*abs(accepted_in)
      out_abs=out_abs+dt*this%cfg%omega_t*massflow_scale*abs(accepted_out)
      in_sum=in_sum+dt*fin(:2)*this%cfg%omega_t*this%cfg%omega_rho*this%cfg%omega_u*this%cfg%AREA
      out_sum=out_sum+dt*fout(:2)*this%cfg%omega_t*this%cfg%omega_rho*this%cfg%omega_u*this%cfg%AREA
      if(tcur*this%cfg%omega_t>=next_output.or.tcur>=T)call write_monitor(dt*this%cfg%omega_t,it)
      if(tcur*this%cfg%omega_t>=next_field.or.tcur>=T)call write_field()
      call write_mass_record(dt*this%cfg%omega_t,norm,0)
      call cpu_time(cpu_now);call system_clock(clock_now)
      elapsed_wall=ieee_value(0.0_wp,ieee_quiet_nan)
      if(clock_rate>0)elapsed_wall=real(clock_now-clock_start,wp)/real(clock_rate,wp)
      write(unitsteps,'(i0,2(",",es24.16),2(",",i0),",",es24.16,6(",",i0),3(",",es24.16))') &
        step,tcur*this%cfg%omega_t,dt*this%cfg%omega_t,it,rej,norm,newton_iterations,jacobian_evaluations, &
        residual_evaluations,linear_solves,backtrack_trials,prepare_calls,cpu_now-wall0,cpu_now-cpu_previous,elapsed_wall
      cpu_previous=cpu_now
      call cpu_time(wall1)
      if(wall1-last_print>10.0_wp.or.tcur>=T)then
        write(output_unit,'(a,f12.4,a,i0,a,es11.3,a,i0,a,i0)') 't(s)=',tcur*this%cfg%omega_t, &
           ' step=',step,' dt(s)=',dt*this%cfg%omega_t,' Newton=',it,' rejects=',rej
        flush(output_unit);flush(unitmon);flush(unitsteps);flush(unitmass);last_print=wall1
      endif
      dt_next=dt
      if(it<=5)dt_next=dt*this%dt_growth
      if(it>=max(8,this%cfg%newton_max-4))dt_next=dt*0.8_wp
      dt=min(dt_next,this%max_dt)
    enddo
    if(tcur<T.and.finstatus==0)finstatus=7000
    call write_field()
    if(finstatus/=0)call write_monitor(0.0_wp,0)
    close(unitmon);close(unitfield);close(unitsteps);close(unitmass)
    if(diagnose_burr15)close(unitburr15)
    open(newunit=unitstat,file=trim(this%cfg%out_dir)//'/run_status.txt',status='replace')
    write(unitstat,'(a,i0)') 'status=',finstatus
    write(unitstat,'(a,es24.16)') 'reached_time_s=',tcur*this%cfg%omega_t
    write(unitstat,'(a,es24.16)') 'requested_time_s=',T*this%cfg%omega_t
    write(unitstat,'(a,i0)') 'accepted_steps=',step
    write(unitstat,'(a,i0)') 'rejected_steps=',rej
    write(unitstat,'(a,i0)') 'last_rejection_status=',lastcode
    call cpu_time(wall1)
    write(unitstat,'(a,es24.16)') 'cpu_seconds=',wall1-wall0
    call system_clock(clock_now)
    if(clock_rate>0)write(unitstat,'(a,es24.16)') 'wall_seconds=',real(clock_now-clock_start,wp)/real(clock_rate,wp)
    write(unitstat,'(a,i0)') 'newton_iterations_total=',newton_iterations
    write(unitstat,'(a,i0)') 'jacobian_evaluations_total=',jacobian_evaluations
    write(unitstat,'(a,i0)') 'residual_evaluations_total=',residual_evaluations
    write(unitstat,'(a,i0)') 'linear_solves_total=',linear_solves
    write(unitstat,'(a,i0)') 'backtrack_trials_total=',backtrack_trials
    write(unitstat,'(a,i0)') 'prepare_calls_total=',prepare_calls
    close(unitstat)
    if(present(status))status=finstatus
  contains
    subroutine mass_inventory(state,mass)
      real(wp),intent(in)::state(:,:)
      real(wp),intent(out)::mass(2,3)
      integer::phase
      do phase=1,2
        mass(phase,1)=sum(state(lo:hi,phase)*this%cfg%cell_width(lo:hi))*mass_scale
        mass(phase,2)=sum(state(lo:ib,phase)*this%cfg%cell_width(lo:ib))*mass_scale
        mass(phase,3)=sum(state(ib+1:hi,phase)*this%cfg%cell_width(ib+1:hi))*mass_scale
      enddo
    end subroutine
    subroutine capture_mass_residual(residual)
      real(wp),intent(in)::residual(:,:)
      integer::phase,np
      np=this%cfg%N_pipe
      ! Preserve the residual of the assembly that actually met tol_res.
      ! This is a signed defect in kg, not a corrected conservation balance.
      do phase=1,2
        mass_residual(phase,1)=sum(residual(phase,:)*this%cfg%cell_width(lo:hi))*mass_scale
        mass_residual(phase,2)=sum(residual(phase,1:np)*this%cfg%cell_width(lo:ib))*mass_scale
        mass_residual(phase,3)=sum(residual(phase,np+1:n)*this%cfg%cell_width(ib+1:hi))*mass_scale
      enddo
    end subroutine
    subroutine write_mass_record(dts,residual_norm,flux_status)
      real(wp),intent(in)::dts,residual_norm
      integer,intent(in)::flux_status
      real(wp)::mass(2,3),error(2,3)
      call mass_inventory(U,mass)
      error(:,1)=mass(:,1)-mass0(:,1)-in_sum+out_sum
      error(:,2)=mass(:,2)-mass0(:,2)-in_sum+base_sum(:,1)
      error(:,3)=mass(:,3)-mass0(:,3)-base_sum(:,2)+out_sum
      write(unitmass,'(i0,3(",",es24.16),",",i0,",",*(es24.16,:,","))') &
        step,tcur*this%cfg%omega_t,dts,residual_norm,flux_status,mass, &
        accepted_in*massflow_scale,accepted_out*massflow_scale,accepted_base*massflow_scale, &
        in_sum,out_sum,base_sum,in_abs,out_abs,base_abs,error,step_balance,mass_residual,mass_residual_sum
    end subroutine
    subroutine write_rejected_state()
      integer::ju,ii
      open(newunit=ju,file=trim(this%cfg%out_dir)//'/last_rejected_state.dat',status='replace')
      write(ju,'(a)') '# t_nd dt_nd residual code; then i,theta,width,Uold(3),Utrial(3),residual(3)'
      write(ju,'(3es25.16,i10)') tcur,dt,norm,lastcode
      do ii=lo,hi
        write(ju,'(i8,11es25.16)')ii,thetav(ii),this%cfg%cell_width(ii),old(ii,:),trial(ii,:),r(:,ii-lo+1)
      enddo
      close(ju)
    end subroutine
    subroutine write_monitor(dts,nit)
      real(wp),intent(in)::dts
      integer,intent(in)::nit
      real(wp)::pr(6,4),v(6),f(3),s(3),aj(3,3),sj(3,3),mpi(2),mr(2),err(2),conv,rlout,rgout
      real(wp)::baseflux(3),baseplus(3),trace(6),massflow_scale
      logical::blocked
      integer::ii,k,cc,indices(4),ib
      ib=lo+this%cfg%N_pipe-1;indices=[lo,ib,ib+1,hi]
      if(allocated(this%lagged_dissipation))then
        call this%face_pair(U(ib,:),U(ib+1,:),thetav(ib),thetav(ib+1),baseflux,baseplus,cc, &
          blocking_active=blocked,junction_trace=trace,dissipation_matrix=this%lagged_dissipation(:,:,ib),left_cell=ib)
      else
        call this%face_pair(U(ib,:),U(ib+1,:),thetav(ib),thetav(ib+1),baseflux,baseplus,cc, &
          blocking_active=blocked,junction_trace=trace,left_cell=ib)
      endif
      if(step==0)then
        if(cc==0)then
          accepted_base(:,1)=baseflux(:2);accepted_base(:,2)=baseplus(:2)
        else
          initial_flux_status=2000+cc
        endif
      else
        ! A later failed trial may have overwritten fout or the lagged cache.
        ! Report only the phase fluxes that advanced the accepted state.
        baseflux(:2)=accepted_base(:,1)
      endif
      massflow_scale=this%cfg%omega_rho*this%cfg%omega_u*this%cfg%AREA
      do k=1,4
        ii=indices(k)
        call this%point_state(U(ii,:),thetav(ii),v,f,s,aj,sj,cc,.false.,cell_closure(this,thetav(ii),ii))
        pr(:,k)=v
      enddo
      conv=this%cfg%omega_rho*this%cfg%AREA*this%cfg%Lr
      do k=1,2
        mpi(k)=sum(U(lo:ib,k)*this%cfg%cell_width(lo:ib))*conv
        mr(k)=sum(U(ib+1:hi,k)*this%cfg%cell_width(ib+1:hi))*conv
      enddo
      inventory=mpi+mr;err=inventory-inventory0-in_sum+out_sum
      rlout=this%cfg%nrho_l0+(this%cfg%nPs-this%cfg%nP_l0)/this%cfg%nCl**2
      rgout=this%cfg%nPs/this%cfg%nCg**2
      write(unitmon,'(*(es24.16,:,","))') tcur*this%cfg%omega_t,dts,pr(4,:)*this%cfg%omega_P, &
        this%cfg%nPs*this%cfg%omega_P,pr(1,2:4), &
        accepted_out(1)/rlout*this%cfg%omega_u,accepted_out(2)/rgout*this%cfg%omega_u, &
        accepted_out*this%cfg%omega_rho*this%cfg%omega_u*this%cfg%AREA, &
        inventory,mpi,mr,in_sum,out_sum,err,real(nit,wp),real(rej,wp),baseflux(:2)*massflow_scale, &
        trace(1)*this%cfg%omega_P,trace(4:5),trace(6)*massflow_scale,merge(1.0_wp,0.0_wp,blocked),1-pr(1,3),1-pr(1,4), &
        accepted_in*massflow_scale,in_abs,out_abs,base_sum(:,1),base_abs(:,1),base_sum(:,2),base_abs(:,2)
      next_output=tcur*this%cfg%omega_t+this%cfg%output_interval_s
    end subroutine
    subroutine write_field()
      real(wp)::v(6),f(3),s(3),aj(3,3),sj(3,3),nan
      type(burr15_result_t)::criterion
      integer::ii,cc,cid
      do ii=lo,hi
        cid=cell_closure(this,thetav(ii),ii)
        call this%point_state(U(ii,:),thetav(ii),v,f,s,aj,sj,cc,.false.,cid)
        write(unitfield,'(*(es24.16,:,","))') tcur*this%cfg%omega_t,vns(ii)*this%cfg%Lr,thetav(ii), &
          v(1),1-v(1),v(4)*this%cfg%omega_P,v(5:6)*this%cfg%omega_u,v(2:3)*this%cfg%omega_rho
        if(diagnose_burr15.and.ii<lo+this%cfg%N_pipe)then
          if(cc==0.and.cid==1)then
            criterion=burr15_stability(v(1),v(5)*this%cfg%omega_u,v(6)*this%cfg%omega_u, &
              v(2)*this%cfg%omega_rho,v(3)*this%cfg%omega_rho,thetav(ii),this%cfg%D_H,this%cfg%AREA,this%cfg%G)
          else
            nan=ieee_value(0.0_wp,ieee_quiet_nan)
            criterion=burr15_result_t(BURR15_NOT_EVALUATED,nan,nan,nan,nan,.false.)
          endif
          write(unitburr15,'(2(es24.16,","),2(i0,","),4(es24.16,","),i0)') &
            tcur*this%cfg%omega_t,vns(ii)*this%cfg%Lr,cid,criterion%status,criterion%hL_over_D, &
            criterion%ug2,criterion%critical_ug2,criterion%margin,merge(1,0,criterion%unstable)
        endif
      enddo
      next_field=tcur*this%cfg%omega_t+this%cfg%field_interval_s
    end subroutine
  end subroutine
end module module_simulate_pipe_riser
