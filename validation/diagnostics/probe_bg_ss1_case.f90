program probe_bg_ss1_case
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_boundary_conditions,only:boundary_flux
  use module_steady_state,only:solve_alpha_from_superficials
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  implicit none
  integer::nc,np,nd,lo,hi,jface,nriser
  real(wp),parameter::pi=acos(-1.0_wp)
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp),allocatable::old(:,:),trial(:,:),state(:,:)
  real(wp)::t,dt,norm_saved,theta,width,dummy(3)
  real(wp)::v(6),f(3),s(3),a(3,3),ds(3,3),fm(3),fp(3),jl(3,3),jr(3,3),kl(3,3),kr(3,3)
  real(wp),allocatable::r(:,:),ab(:,:,:),bb(:,:,:),cb(:,:,:),dx(:,:),candidate(:,:)
  real(wp)::norm,fin(3),fout(3),smallalpha,lambda
  integer::i,idx,iu,reason,code,profile,facecode,boundcode,nonesided,k,info,pointfails
  character(200)::header
  ! Generated from the frozen run metadata by build_bg_ss1_case.py.
  open(newunit=iu,file='case_config.dat',status='old')
  read(iu,*)cfg%N_pipe,nriser,cfg%num_dummy_cells
  nd=cfg%num_dummy_cells;np=cfg%N_pipe+nriser;nc=np+2*nd;cfg%nCells=nc
  lo=nd+1;hi=nc-nd;jface=nd+cfg%N_pipe
  read(iu,*)cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%Lr,cfg%Lp
  cfg%omega_t=cfg%Lr/cfg%omega_u
  read(iu,*)cfg%D_H,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma
  cfg%AREA=pi*cfg%D_H**2/4
  read(iu,*)cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0,cfg%nPs
  read(iu,*)cfg%mul,cfg%mug,cfg%njl,cfg%njg
  read(iu,*)cfg%tol,cfg%tola,cfg%tol_res,cfg%tol_inc
  read(iu,*)cfg%closure_pipe,cfg%closure_riser,cfg%quad_points,cfg%maxit
  read(iu,*)cfg%jac_dissipation_fd
  read(iu,'(a)')cfg%junction_model
  read(iu,'(a)')cfg%flux_scheme
  read(iu,'(a)')cfg%dissipation_mode
  close(iu)
  allocate(old(nc,3),trial(nc,3),state(nc,3),candidate(nc,3))
  allocate(r(3,np),ab(3,3,np),bb(3,3,np),cb(3,3,np),dx(3,np))
  allocate(cfg%cell_width(nc),cfg%thetav(nc))
  open(newunit=iu,file='last_rejected_state.dat',status='old')
  read(iu,'(a)')header
  read(iu,*)t,dt,norm_saved,reason
  do i=lo,hi
    read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
    if(idx/=i)error stop 'snapshot cell order'
    cfg%thetav(i)=theta;cfg%cell_width(i)=width
  enddo
  close(iu)
  do i=1,nd
    old(i,:)=old(lo,:);old(hi+i,:)=old(hi,:)
    trial(i,:)=trial(lo,:);trial(hi+i,:)=trial(hi,:)
    cfg%thetav(i)=cfg%thetav(lo);cfg%thetav(hi+i)=cfg%thetav(hi)
    cfg%cell_width(i)=cfg%cell_width(lo);cfg%cell_width(hi+i)=cfg%cell_width(hi)
  enddo
  call sim%init(cfg)
  write(*,'(a,3es25.16,i8)')'snapshot time_s,dt_s,norm,status ',t*cfg%omega_t,dt*cfg%omega_t,norm_saved,reason
  call sim%prepare_dissipation(old,code)
  write(*,'(a,i8)')'prepare_dissipation status ',code
  if(code/=0)error stop 'old-state viscosity failed'
  open(newunit=iu,file='point_states.csv',status='replace')
  write(iu,'(a)')'profile,cell,status,theta_rad,alpha,p_Pa,ul_m_s,ug_m_s,jl_m_s,jg_m_s'
  do profile=1,2
    state=old
    if(profile==2)state=trial
    smallalpha=1
    do i=lo,hi
      call sim%point_state(state(i,:),cfg%thetav(i),v,f,s,a,ds,code)
      if(code/=0)write(*,'(a,3i8)')'point failure profile,cell,status ',profile,i,code
      if(code==0)smallalpha=min(smallalpha,v(1))
      write(iu,'(i2,a,i4,a,i6,7(a,es25.16))')profile,',',i,',',code,',',cfg%thetav(i),',',v(1), &
        ',',v(4)*cfg%omega_P,',',v(5)*cfg%omega_u,',',v(6)*cfg%omega_u, &
        ',',(1-v(1))*v(5)*cfg%omega_u,',',v(1)*v(6)*cfg%omega_u
    enddo
    write(*,'(a,i3,es25.16)')'minimum physical alpha profile ',profile,smallalpha
    do i=lo,hi-1
      call sim%face_pair(state(i,:),state(i+1,:),cfg%thetav(i),cfg%thetav(i+1),fm,fp,facecode, &
        dissipation_matrix=sim%lagged_dissipation(:,:,i))
      if(facecode/=0)write(*,'(a,3i8)')'face failure profile,leftcell,status ',profile,i,facecode
    enddo
    call boundary_flux(state(lo,:),cfg%thetav(lo),cfg,.true.,f,boundcode)
    write(*,'(a,2i8)')'inlet boundary profile,status ',profile,boundcode
    call boundary_flux(state(hi,:),cfg%thetav(hi),cfg,.false.,f,boundcode)
    write(*,'(a,2i8)')'outlet boundary profile,status ',profile,boundcode
    call junction_probe(state,profile)
    call sim%assemble(state,old,dt,r,ab,bb,cb,norm,fin,fout,code,.false.)
    write(*,'(a,i3,i8,es25.16)')'assembly residual profile,status,norm ',profile,code,norm
    call sim%assemble(state,old,dt,r,ab,bb,cb,norm,fin,fout,code,.true.)
    write(*,'(a,i3,i8,es25.16)')'assembly Jacobian profile,status,norm ',profile,code,norm
    if(profile==2.and.code==0)then
      call sim%ti%solve_block_tridiagonal_3x3(ab,bb,cb,-r,dx,info)
      write(*,'(a,i8,es25.16)')'Newton band solve status,max increment ',info,maxval(abs(dx))
      do k=0,4
        lambda=2.0_wp**(-k);candidate=state
        do i=lo,hi
          candidate(i,:)=state(i,:)+lambda*dx(:,i-lo+1)
        enddo
        pointfails=0
        do i=lo,hi
          call sim%point_state(candidate(i,:),cfg%thetav(i),v,f,s,a,ds,code)
          if(code/=0)pointfails=pointfails+1
        enddo
        call sim%face_pair(candidate(jface,:),candidate(jface+1,:),cfg%thetav(jface),cfg%thetav(jface+1),fm,fp,facecode, &
          dissipation_matrix=sim%lagged_dissipation(:,:,jface))
        write(*,'(a,es25.16,2i8)')'Newton candidate lambda,point failures,junction status ',lambda,pointfails,facecode
        call junction_probe(candidate,10+k)
      enddo
    endif
  enddo
  close(iu)
contains
  subroutine junction_probe(uu,iprofile)
    real(wp),intent(in)::uu(:,:)
    integer,intent(in)::iprofile
    real(wp)::vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),source(3),aa(3,3),ss(3,3)
    real(wp)::astl,astr,mapped(3),raw(3),pip,kin,pstar,rl,rg,jli,jga,df(3),convective(3),viscous(3)
    integer::st,nroots
    call sim%point_state(uu(jface,:),cfg%thetav(jface),vl,fl,source,aa,ss,st)
    call sim%point_state(uu(jface+1,:),cfg%thetav(jface+1),vr,fr,source,aa,ss,st)
    astl=solve_alpha_from_superficials(fr(1)/vr(2),fr(2)/vr(3),vr(2),vr(3),cfg%thetav(jface), &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol, &
      status=st,closure_id=1,alpha_hint=vl(1))
    write(*,'(a,i3,i8,es25.16)')'map right trace to pipe profile,status,alpha ',iprofile,st,astl
    if(st/=0)return
    mapped=[(1-astl)*vr(2),astl*vr(3),fr(1)+fr(2)]
    call sim%point_state(mapped,cfg%thetav(jface),vm,fs,source,aa,ss,st)
    convective=.5_wp*(fl+fs)
    viscous=.5_wp*matmul(sim%lagged_dissipation(:,:,jface),mapped-uu(jface,:))
    raw=convective-viscous
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    kin=.5_wp*(fl(3)-pip*vl(4)+fs(3)-pip*vm(4))
    pstar=(raw(3)-kin)/pip
    rl=cfg%nrho_l0+(pstar-cfg%nP_l0)/cfg%nCl**2;rg=pstar/cfg%nCg**2
    jli=raw(1)/rl;jga=raw(2)/rg
    write(*,'(a,i3,3es25.16)')'junction pstar_Pa,jl_m_s,jg_m_s profile ',iprofile,pstar*cfg%omega_P,jli*cfg%omega_u,jga*cfg%omega_u
    write(*,'(a,i3,3es25.16)')'junction gas convective,diffusive,net kg_s profile ',iprofile, &
      convective(2)*cfg%omega_rho*cfg%omega_u*cfg%AREA,viscous(2)*cfg%omega_rho*cfg%omega_u*cfg%AREA, &
      raw(2)*cfg%omega_rho*cfg%omega_u*cfg%AREA
    astr=solve_alpha_from_superficials(jli,jga,rl,rg,0.0_wp,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=st,closure_id=3,nroots=nroots)
    write(*,'(a,i3,2i8,es25.16)')'riser trace root profile,status,nroots,alpha ',iprofile,st,nroots,astr
    write(*,'(a,i3,6es25.16)')'neighbor velocities ulL,ugL,ulR,ugR,alphaL,alphaR profile ',iprofile, &
      vl(5)*cfg%omega_u,vl(6)*cfg%omega_u,vr(5)*cfg%omega_u,vr(6)*cfg%omega_u,vl(1),vr(1)
  end subroutine
end program
