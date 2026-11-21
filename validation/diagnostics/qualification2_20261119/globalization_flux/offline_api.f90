module offline_api
  use iso_c_binding, only: c_int,c_double,c_char
  use mod_kinds, only: wp
  use module_config, only: sim_config_t
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use module_find, only: pressure_fun_uj
  use module_closure_law, only: CdUd_bhagwat
  use module_steady_state, only: solve_alpha_from_superficials
  use module_lapack_small, only: dgeev_3x3
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  integer,parameter :: nc=106,np=102,nx=306,lo=3,hi=104,ib=53,nmargin=53
  type(sim_config_t),save :: cfg
  type(pipeline_riser_simulator_t),save :: sim
  real(wp),save :: old(nc,3),trial(nc,3),saved_r(3,np),dt,t,saved_norm
contains
  subroutine probe_init(folder,nchar,mode,status) bind(C)
    character(c_char),intent(in) :: folder(*)
    integer(c_int),value :: nchar,mode
    integer(c_int),intent(out) :: status
    character(len=:),allocatable :: path
    character(len=200) :: header
    real(wp) :: ps,jl0,jg0,ml0,mg0,theta,width,dummy(3),wr(3),wi(3),vr(3,3),speed
    integer :: iu,i,idx,reason,unused,code
    allocate(character(nchar)::path)
    do i=1,nchar
      path(i:i)=folder(i)
    enddo
    cfg=sim_config_t()
    open(newunit=iu,file=path//'/probe_parameters.dat',status='old',iostat=code)
    if(code/=0)then
      status=90001;return
    endif
    read(iu,*)cfg%omega_u,ps,jl0,jg0,ml0,mg0,unused
    close(iu)
    cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp
    cfg%Lr=12.546253993197080_wp;cfg%Lp=57.4_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
    cfg%D_H=.05254_wp;cfg%AREA=acos(-1.0_wp)*cfg%D_H**2/4
    cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
    cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
    cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
    cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
    cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=ps/cfg%omega_P
    cfg%mul=ml0/(cfg%omega_rho*cfg%omega_u);cfg%mug=mg0/(cfg%omega_rho*cfg%omega_u)
    cfg%njl=jl0/cfg%omega_u;cfg%njg=jg0/cfg%omega_u
    cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%tol_res=1e-8_wp;cfg%tol_inc=1e-8_wp
    cfg%closure_pipe=1;cfg%closure_riser=3;cfg%gas_eos='BAROTROPIC'
    ! Riser closure of the fixture (default 3 = historical Bhagwat fixtures;
    ! PRIMS_PROBE_CLOSURE_RISER=5 for CHEXAL_1992 fixtures, 21/11/2026).
    call get_environment_variable('PRIMS_PROBE_CLOSURE_RISER',header,status=code)
    if(code==0)read(header,*)cfg%closure_riser
    call get_environment_variable('PRIMS_JUNCTION_TRANSMISSION',header,status=code)
    if(code==0)cfg%junction_transmission=trim(header)
    call get_environment_variable('PRIMS_JUNCTION_SEAL',header,status=code)
    if(code==0)cfg%junction_seal=trim(header)
    cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%dissipation_mode='LAGGED'
    cfg%jac_dissipation_fd=.true.
    cfg%nCells=nc;cfg%num_dummy_cells=2;cfg%N_pipe=51;cfg%quad_points=3
    allocate(cfg%cell_width(nc),cfg%thetav(nc))
    open(newunit=iu,file=path//'/last_rejected_state.dat',status='old')
    read(iu,'(a)')header
    read(iu,*)t,dt,saved_norm,reason
    do i=lo,hi
      read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
      if(idx/=i)then
        status=90002;close(iu);return
      endif
      cfg%thetav(i)=theta;cfg%cell_width(i)=width;saved_r(:,i-lo+1)=dummy
    enddo
    close(iu)
    do i=1,2
      old(i,:)=old(lo,:);trial(i,:)=trial(lo,:)
      old(hi+i,:)=old(hi,:);trial(hi+i,:)=trial(hi,:)
      cfg%thetav(i)=cfg%thetav(lo);cfg%cell_width(i)=cfg%cell_width(lo)
      cfg%thetav(hi+i)=cfg%thetav(hi);cfg%cell_width(hi+i)=cfg%cell_width(hi)
    enddo
    if(mode==2)cfg%flux_scheme='RUSANOV'
    call sim%init(cfg)
    call sim%prepare_dissipation(old,code)
    status=code
    if(code/=0)return
    if(mode==3)then
      call dgeev_3x3(sim%lagged_dissipation(:,:,ib),wr,wi,vr,code)
      status=code
      if(code/=0)return
      speed=maxval(hypot(wr,wi))
      sim%lagged_dissipation(:,:,ib)=0
      do i=1,3
        sim%lagged_dissipation(i,i,ib)=speed
      enddo
    endif
  end subroutine

  subroutine probe_state(which,x,scales,metadata) bind(C)
    integer(c_int),value :: which
    real(c_double),intent(out) :: x(nx),scales(nx),metadata(6)
    integer :: i,k
    do i=lo,hi
      k=3*(i-lo)
      x(k+1:k+3)=trial(i,:)
      if(which==0)x(k+1:k+3)=old(i,:)
      scales(k+1:k+3)=[max(old(i,1),1e-3_wp),max(old(i,2),1e-7_wp),max(abs(old(i,3)),1.0_wp)]
    enddo
    metadata=[t*cfg%omega_t,dt*cfg%omega_t,saved_norm,cfg%omega_u,cfg%omega_rho,cfg%omega_P]
  end subroutine

  subroutine unpack(x,u)
    real(wp),intent(in) :: x(nx)
    real(wp),intent(out) :: u(nc,3)
    integer :: i,k
    u=old
    do i=lo,hi
      k=3*(i-lo);u(i,:)=x(k+1:k+3)
    enddo
    u(1,:)=u(lo,:);u(2,:)=u(lo,:);u(105,:)=u(hi,:);u(106,:)=u(hi,:)
  end subroutine

  subroutine probe_eval(x,want_jac,res,jac,norm,fin,fout,status,index,stage) bind(C)
    real(c_double),intent(in) :: x(nx)
    integer(c_int),value :: want_jac
    real(c_double),intent(out) :: res(nx),jac(nx,nx),norm,fin(3),fout(3)
    integer(c_int),intent(out) :: status,index,stage
    real(wp) :: u(nc,3),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np)
    integer :: code,location,i,k
    character(len=32) :: phase
    call unpack(x,u)
    call sim%assemble(u,old,dt,r,a,b,c,norm,fin,fout,code,want_jac/=0, &
      failure_index=location,failure_stage=phase)
    status=code;index=location;stage=0
    select case(trim(phase))
    case('CELL');stage=1
    case('FACE');stage=2
    case('BOUNDARY');stage=3
    case('BOUNDARY_JAC');stage=4
    case('RESIDUAL');stage=5
    end select
    if(code/=0)then
      res=ieee_value(0.0_wp,ieee_quiet_nan)
      jac=res(1);norm=res(1);fin=res(1);fout=res(1)
      return
    endif
    res=reshape(r,[nx]);jac=0
    if(want_jac/=0)then
      do i=1,np
        k=3*(i-1)
        jac(k+1:k+3,k+1:k+3)=b(:,:,i)
        if(i>1)jac(k+1:k+3,k-2:k)=a(:,:,i)
        if(i<np)jac(k+1:k+3,k+4:k+6)=c(:,:,i)
      enddo
    endif
  end subroutine

  subroutine junction_raw(u,raw,pstar,vl,vr,vm,fl,fr,fs,mapped,status)
    real(wp),intent(in) :: u(nc,3)
    real(wp),intent(out) :: raw(3),pstar,vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),mapped(3)
    integer,intent(out) :: status
    real(wp) :: source(3),a(3,3),ds(3,3),alpha,pip,kin
    call sim%point_state(u(ib,:),cfg%thetav(ib),vl,fl,source,a,ds,status,closure_id=1)
    if(status/=0)return
    call sim%point_state(u(ib+1,:),cfg%thetav(ib+1),vr,fr,source,a,ds,status,closure_id=cfg%closure_riser)
    if(status/=0)return
    alpha=solve_alpha_from_superficials(fr(1)/vr(2),fr(2)/vr(3),vr(2),vr(3),cfg%thetav(ib), &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho, &
      cfg%tol,status=status,closure_id=1,alpha_hint=vl(1))
    if(status/=0)return
    mapped=[(1-alpha)*vr(2),alpha*vr(3),fr(1)+fr(2)]
    call sim%point_state(mapped,cfg%thetav(ib),vm,fs,source,a,ds,status,closure_id=1)
    if(status/=0)return
    raw=.5_wp*(fl+fs)-.5_wp*matmul(sim%lagged_dissipation(:,:,ib),mapped-u(ib,:))
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    kin=.5_wp*(fl(3)-pip*vl(4)+fs(3)-pip*vm(4))
    pstar=(raw(3)-kin)/pip
  end subroutine

  subroutine probe_margins(x,margin,status) bind(C)
    real(c_double),intent(in) :: x(nx)
    real(c_double),intent(out) :: margin(nmargin)
    integer(c_int),intent(out) :: status
    real(wp) :: u(nc,3),p,rl,rg,alpha,ug,cd,ud,raw(3),ps,vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),mapped(3)
    integer :: i,code
    margin=ieee_value(0.0_wp,ieee_quiet_nan)
    call unpack(x,u)
    call junction_raw(u,raw,ps,vl,vr,vm,fl,fr,fs,mapped,code)
    status=code
    if(code/=0)return
    margin(1)=raw(2)/(ps/cfg%nCg**2)*cfg%omega_u
    margin(2)=ps*cfg%omega_P
    do i=ib+1,hi
      p=pressure_fun_uj(u(i,1),u(i,2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
      rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2;alpha=u(i,2)/rg
      ug=u(i,3)/u(i,2)
      call CdUd_bhagwat(alpha,rl,rg,0.0_wp,ug,cfg%thetav(i),cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
        cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,cd,ud,code)
      if(code/=0)then
        status=code;return
      endif
      ! CdUd_bhagwat returns dimensional Ud despite nondimensional inputs.
      margin(2+i-ib)=(ug-cd*alpha*ug)*cfg%omega_u-ud
    enddo
  end subroutine

  subroutine probe_interface(x,data,status) bind(C)
    real(c_double),intent(in) :: x(nx)
    real(c_double),intent(out) :: data(48)
    integer(c_int),intent(out) :: status
    real(wp) :: u(nc,3),raw(3),ps,vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),mapped(3)
    real(wp) :: a(3,3),ds(3,3),source(3),v(6),f(3),wr(3),wi(3),ev(3,3),sl,sr,hll(3),middle(3)
    integer :: code
    data=ieee_value(0.0_wp,ieee_quiet_nan)
    call unpack(x,u)
    call junction_raw(u,raw,ps,vl,vr,vm,fl,fr,fs,mapped,code)
    status=code
    if(code/=0)return
    call sim%point_state(u(ib,:),cfg%thetav(ib),v,f,source,a,ds,code,closure_id=1)
    status=code
    if(code/=0)return
    call dgeev_3x3(a,wr,wi,ev,code)
    status=code
    if(code/=0)return
    if(any(wi/=0.0_wp))then
      status=90003;return
    endif
    sl=min(0.0_wp,minval(wr),vl(5),vl(6));sr=max(0.0_wp,maxval(wr),vl(5),vl(6))
    call sim%point_state(mapped,cfg%thetav(ib),v,f,source,a,ds,code,closure_id=1)
    status=code
    if(code/=0)return
    call dgeev_3x3(a,wr,wi,ev,code)
    status=code
    if(code/=0)return
    if(any(wi/=0.0_wp))then
      status=90003;return
    endif
    sl=min(sl,minval(wr),vm(5),vm(6));sr=max(sr,maxval(wr),vm(5),vm(6))
    hll=(sr*fl-sl*fs+sl*sr*(mapped-u(ib,:)))/(sr-sl)
    middle=(sr*mapped-sl*u(ib,:)-(fs-fl))/(sr-sl)
    data(1:3)=raw;data(4)=ps;data(5:10)=vl;data(11:16)=vr;data(17:22)=vm
    data(23:25)=u(ib,:);data(26:28)=mapped;data(29:31)=fl;data(32:34)=fs
    data(35)=sl;data(36)=sr;data(37:39)=hll;data(40:42)=middle
    data(43:45)=u(ib+1,:);data(46:48)=[cfg%omega_u,cfg%omega_rho,cfg%AREA]
  end subroutine
end module
