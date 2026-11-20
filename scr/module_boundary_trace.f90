module module_boundary_trace
  !! Source-compatible, local steady extrapolation to an actual exterior face.
  !! It is NOT a nonreflecting characteristic or an exact nonlinear Riemann BC.
  use mod_kinds,only:wp
  use module_config,only:sim_config_t,first_physical_cell,last_physical_cell
  use module_find,only:compute_primitive_variables
  use module_steady_state,only:fun_dpds,state_at
  use module_jacobian,only:DFduj
  use module_lapack_small,only:dgeev_3x3
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  integer,parameter,public::BC_GEOMETRY_REQUIRED=130,BC_TRACE_FAILED=131,BC_DATA_REQUIRED=132
  integer,parameter,public::BC_SPECTRUM_INVALID=133,BC_PHASE_REQUIRED=134
  public::steady_boundary_trace,check_boundary_characteristics
contains
  subroutine steady_boundary_trace(u,theta,cfg,inlet,uf,angle,status)
    real(wp),intent(in)::u(3),theta
    type(sim_config_t),intent(in)::cfg
    logical,intent(in)::inlet
    real(wp),intent(out)::uf(3),angle
    integer,intent(out)::status
    real(wp)::v(6),alpha,rl,rg,p,ul,ug,jl,jg,kin,gl,gg,s,target,h,y,half,tmp,error,scale,leftover
    integer::i,model,section,phase,it,code
    status=BC_GEOMETRY_REQUIRED;uf=u;angle=theta
    if(.not.allocated(cfg%vns).or..not.allocated(cfg%face_positions).or.cfg%CA<=0)return
    i=last_physical_cell(cfg);model=cfg%closure_riser;section=2
    target=cfg%face_positions(ubound(cfg%face_positions,1))
    if(inlet)then
      i=first_physical_cell(cfg);model=cfg%closure_pipe;section=1
      target=cfg%face_positions(lbound(cfg%face_positions,1))
    endif
    call compute_primitive_variables(u(1),u(2),u(3),theta,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0, &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%njl,cfg%omega_u,cfg%omega_rho, &
      cfg%tol,cfg%tola,cfg%maxit,alpha,rl,rg,p,ul,ug,phase,status=code,closure_id=model)
    if(code/=0)then
      status=code;return
    endif
    if(alpha==0.or.alpha==1)then
      status=BC_PHASE_REQUIRED;return
    endif
    gl=u(1)*ul;gg=u(2)*ug
    if(gl==0.and.gg==0)then
      status=BC_DATA_REQUIRED;return ! zero fluxes do not determine mixed holdup
    endif
    s=cfg%vns(i);h=target-s
    do it=1,256
      leftover=target-s
      if(abs(leftover)<=32*epsilon(1._wp)*max(1._wp,abs(target)))exit
      if(abs(h)>abs(leftover))h=leftover
      if(abs(h)<64*epsilon(1._wp)*max(1._wp,abs(s)))then
        status=BC_TRACE_FAILED;return
      endif
      call rk(s,p,h,y,code)
      if(code==0)call rk(s,p,h/2,half,code)
      if(code==0)call rk(s+h/2,half,h/2,tmp,code)
      if(code/=0)then
        h=h/2;cycle
      endif
      error=abs(tmp-y)/15;scale=1e-11_wp*(1+max(abs(tmp),abs(p)))
      if(error<=scale)then
        s=s+h;p=tmp
      endif
      if(error==0)then
        h=2*h
      else
        h=h*min(2._wp,max(.2_wp,.9_wp*(scale/error)**.2_wp))
      endif
    enddo
    if(it>256)then
      status=BC_TRACE_FAILED;return
    endif
    angle=-cfg%beta
    if(.not.inlet)angle=atan((target*cfg%Lr-cfg%Lp)/cfg%CA)
    call state_at(p,angle,gl,gg,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%nCl,cfg%nCg,cfg%nrho_l0, &
      cfg%nP_l0,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,model, &
      rl,rg,alpha,ul,ug,jl,jg,kin,code)
    if(code/=0)then
      status=code;return
    endif
    uf=[(1-alpha)*rl,alpha*rg,gl+gg];status=0
  contains
    subroutine rhs(x,y,f,st)
      real(wp),intent(in)::x,y
      real(wp),intent(out)::f
      integer,intent(out)::st
      f=fun_dpds(x,y,gl,gg,cfg%Lp,cfg%Lr,cfg%CA,cfg%beta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
        cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0,cfg%mu_l,cfg%mu_g,cfg%sigma, &
        cfg%omega_P,cfg%omega_u,cfg%omega_rho,cfg%tol,st,cfg%closure_riser,section,cfg%closure_pipe)
    end subroutine
    subroutine rk(x,y,h,value,st)
      real(wp),intent(in)::x,y,h
      real(wp),intent(out)::value
      integer,intent(out)::st
      real(wp)::k(4)
      call rhs(x,y,k(1),st);if(st/=0)return
      call rhs(x+h/2,y+h*k(1)/2,k(2),st);if(st/=0)return
      call rhs(x+h/2,y+h*k(2)/2,k(3),st);if(st/=0)return
      call rhs(x+h,y+h*k(3),k(4),st);if(st/=0)return
      value=y+h*(k(1)+2*k(2)+2*k(3)+k(4))/6
      if(.not.ieee_is_finite(value).or.value<=0)st=BC_TRACE_FAILED
    end subroutine
  end subroutine

  subroutine check_boundary_characteristics(u,theta,cfg,inlet,status,ul_ref,ug_ref)
    real(wp),intent(in)::u(3),theta
    type(sim_config_t),intent(in)::cfg
    logical,intent(in)::inlet
    integer,intent(out)::status
    real(wp),intent(in),optional::ul_ref,ug_ref
    real(wp)::alpha,rl,rg,p,ul,ug,a(3,3),wr(3),wi(3),r(3,3)
    integer::model,phase,code,n
    model=cfg%closure_riser;if(inlet)model=cfg%closure_pipe
    call compute_primitive_variables(u(1),u(2),u(3),theta,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0, &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%njl,cfg%omega_u,cfg%omega_rho, &
      cfg%tol,cfg%tola,cfg%maxit,alpha,rl,rg,p,ul,ug,phase,status=code,closure_id=model,ul_prev=ul_ref,ug_prev=ug_ref)
    status=code;if(code/=0)return
    if(alpha==0.or.alpha==1)then
      status=BC_PHASE_REQUIRED;return
    endif
    a=DFduj(alpha,rl,rg,ul,ug,cfg%nCl,cfg%nCg,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l, &
      cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%tol,status=code,closure_id=model)
    status=BC_SPECTRUM_INVALID;if(code/=0)return
    call dgeev_3x3(a,wr,wi,r,code)
    if(code/=0.or.any(wi/=0).or.any(wr==0))return
    n=count(wr<0);if(inlet)n=count(wr>0)
    status=BC_DATA_REQUIRED
    if((inlet.and.n==2).or.(.not.inlet.and.n==1))status=0
  end subroutine
end module
