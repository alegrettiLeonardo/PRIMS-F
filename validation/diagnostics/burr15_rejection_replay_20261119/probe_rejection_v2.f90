program probe_bg_ss1_stall
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_boundary_conditions,only:boundary_flux
  use module_steady_state,only:solve_alpha_from_superficials
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  implicit none
  integer,parameter::nc=106,np=102
  real(wp),parameter::pi=acos(-1.0_wp)
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp)::old(nc,3),trial(nc,3),state(nc,3),t,dt,norm_saved,theta,width,dummy(3),rsaved(3,np)
  real(wp)::v(6),f(3),s(3),a(3,3),ds(3,3),fm(3),fp(3),jl(3,3),jr(3,3),kl(3,3),kr(3,3)
  real(wp)::r(3,np),ab(3,3,np),bb(3,3,np),cb(3,3,np),norm,fin(3),fout(3),smallalpha,dx(3,np),candidate(nc,3),lambda,lambda0,nt,ps_si,jl_si,jg_si,mul_si,mug_si
  integer::i,idx,iu,reason,code,profile,facecode,boundcode,nonesided,k,info,pointfails,firstbad,firststatus,jacmode,firstface,firstfacestatus
  character(200)::header
  open(newunit=iu,file='probe_parameters.dat',status='old')
  read(iu,*)cfg%omega_u,ps_si,jl_si,jg_si,mul_si,mug_si,jacmode
  close(iu)
  cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp
  cfg%Lr=12.546253993197080_wp;cfg%Lp=57.4_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
  cfg%D_H=.05254_wp;cfg%AREA=pi*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
  cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
  cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=ps_si/cfg%omega_P
  cfg%mul=mul_si/(cfg%omega_rho*cfg%omega_u)
  cfg%mug=mug_si/(cfg%omega_rho*cfg%omega_u)
  cfg%njl=jl_si/cfg%omega_u;cfg%njg=jg_si/cfg%omega_u
  cfg%maxit=100;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp
  cfg%tol_res=1e-8_wp;cfg%tol_inc=1e-8_wp
  cfg%closure_pipe=1;cfg%closure_riser=3;cfg%jac_dissipation_fd=jacmode==1
  cfg%gas_eos='BAROTROPIC';cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%dissipation_mode='LAGGED'
  cfg%nCells=nc;cfg%num_dummy_cells=2;cfg%N_pipe=51;cfg%quad_points=3
  allocate(cfg%cell_width(nc),cfg%thetav(nc))
  open(newunit=iu,file='last_rejected_state.dat',status='old')
  read(iu,'(a)')header
  read(iu,*)t,dt,norm_saved,reason
  do i=3,104
    read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
    if(idx/=i)error stop 'snapshot cell order'
    cfg%thetav(i)=theta;cfg%cell_width(i)=width;rsaved(:,i-2)=dummy
  enddo
  close(iu)
  old(1,:)=old(3,:);old(2,:)=old(3,:);old(105,:)=old(104,:);old(106,:)=old(104,:)
  trial(1,:)=trial(3,:);trial(2,:)=trial(3,:);trial(105,:)=trial(104,:);trial(106,:)=trial(104,:)
  cfg%thetav(1:2)=cfg%thetav(3);cfg%thetav(105:106)=cfg%thetav(104)
  cfg%cell_width(1:2)=cfg%cell_width(3);cfg%cell_width(105:106)=cfg%cell_width(104)
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
    do i=3,104
      call sim%point_state(state(i,:),cfg%thetav(i),v,f,s,a,ds,code)
      if(code/=0)write(*,'(a,3i8)')'point failure profile,cell,status ',profile,i,code
      if(code==0)smallalpha=min(smallalpha,v(1))
      write(iu,'(i2,a,i4,a,i6,7(a,es25.16))')profile,',',i,',',code,',',cfg%thetav(i),',',v(1), &
        ',',v(4)*cfg%omega_P,',',v(5)*cfg%omega_u,',',v(6)*cfg%omega_u, &
        ',',(1-v(1))*v(5)*cfg%omega_u,',',v(1)*v(6)*cfg%omega_u
    enddo
    write(*,'(a,i3,es25.16)')'minimum physical alpha profile ',profile,smallalpha
    do i=3,103
      call sim%face_pair(state(i,:),state(i+1,:),cfg%thetav(i),cfg%thetav(i+1),fm,fp,facecode, &
        dissipation_matrix=sim%lagged_dissipation(:,:,i),left_cell=i)
      if(facecode/=0)write(*,'(a,3i8)')'face failure profile,leftcell,status ',profile,i,facecode
    enddo
    call boundary_flux(state(3,:),cfg%thetav(3),cfg,.true.,f,boundcode)
    write(*,'(a,2i8)')'inlet boundary profile,status ',profile,boundcode
    call boundary_flux(state(104,:),cfg%thetav(104),cfg,.false.,f,boundcode)
    write(*,'(a,2i8)')'outlet boundary profile,status ',profile,boundcode
    call junction_probe(state,profile)
    call sim%assemble(state,old,dt,r,ab,bb,cb,norm,fin,fout,code,.false.)
    write(*,'(a,i3,i8,es25.16)')'assembly residual profile,status,norm ',profile,code,norm
    call sim%assemble(state,old,dt,r,ab,bb,cb,norm,fin,fout,code,.true.)
    write(*,'(a,i3,i8,es25.16)')'assembly Jacobian profile,status,norm ',profile,code,norm
    if(profile==2.and.code==0)write(*,'(a,es25.16)')'residual difference from snapshot ',maxval(abs(r-rsaved))
    if(profile==2.and.code==0)then
      call sim%ti%solve_block_tridiagonal_3x3(ab,bb,cb,-r,dx,info)
      write(*,'(a,i8,es25.16)')'Newton band solve status,max increment ',info,maxval(abs(dx))
      lambda0=1
      do i=3,104
        do k=1,2
          if(dx(k,i-2)<0)lambda0=min(lambda0,-0.90_wp*state(i,k)/dx(k,i-2))
        enddo
      enddo
      do k=0,17
        lambda=lambda0*2.0_wp**(-k);candidate=state
        do i=3,104
          candidate(i,:)=state(i,:)+lambda*dx(:,i-2)
        enddo
        call sim%assemble(candidate,old,dt,r,ab,bb,cb,nt,fin,fout,code,.false.)
        write(*,'(a,i3,es25.16,i8,es25.16,l2)')'proposal k,lambda,status,norm,Armijo ',k,lambda,code,nt, &
          code==0.and.(nt<norm*(1-1e-4_wp*lambda).or.nt<=cfg%tol_res)
        pointfails=0;firstbad=0;firststatus=0
        do i=3,104
          call sim%point_state(candidate(i,:),cfg%thetav(i),v,f,s,a,ds,code)
          if(code/=0)then
            pointfails=pointfails+1
            if(firstbad==0)then
              firstbad=i;firststatus=code
            endif
          endif
        enddo
        call sim%face_pair(candidate(53,:),candidate(54,:),cfg%thetav(53),cfg%thetav(54),fm,fp,facecode, &
          dissipation_matrix=sim%lagged_dissipation(:,:,53),left_cell=53)
        write(*,'(a,es25.16,2i8)')'Newton candidate lambda,point failures,junction status ',lambda,pointfails,facecode
        write(*,'(a,2i8)')'first point failure cell,status ',firstbad,firststatus
        firstface=0;firstfacestatus=0
        if(pointfails==0)then
          do i=3,103
            call sim%face_pair(candidate(i,:),candidate(i+1,:),cfg%thetav(i),cfg%thetav(i+1),fm,fp,code, &
              dissipation_matrix=sim%lagged_dissipation(:,:,i),left_cell=i)
            if(code/=0)then
              firstface=i;firstfacestatus=code;exit
            endif
          enddo
          write(*,'(a,2i8)')'first face failure leftcell,status ',firstface,firstfacestatus
        endif
        call junction_probe(candidate,10+k)
      enddo
    endif
  enddo
  close(iu)
contains
  subroutine junction_probe(uu,iprofile)
    real(wp),intent(in)::uu(nc,3)
    integer,intent(in)::iprofile
    real(wp)::vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),source(3),aa(3,3),ss(3,3)
    real(wp)::astl,astr,mapped(3),raw(3),pip,kin,pstar,rl,rg,jli,jga,df(3),convective(3),viscous(3)
    integer::st,nroots
    call sim%point_state(uu(53,:),cfg%thetav(53),vl,fl,source,aa,ss,st)
    if(st/=0)then
      write(*,'(a,2i8)')'junction left cell status ',iprofile,st
      return
    endif
    call sim%point_state(uu(54,:),cfg%thetav(54),vr,fr,source,aa,ss,st)
    if(st/=0)then
      write(*,'(a,2i8)')'junction right cell status ',iprofile,st
      return
    endif
    astl=solve_alpha_from_superficials(fr(1)/vr(2),fr(2)/vr(3),vr(2),vr(3),cfg%thetav(53), &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol, &
      status=st,closure_id=1,alpha_hint=vl(1))
    write(*,'(a,i3,i8,es25.16)')'map right trace to pipe profile,status,alpha ',iprofile,st,astl
    if(st/=0)return
    mapped=[(1-astl)*vr(2),astl*vr(3),fr(1)+fr(2)]
    call sim%point_state(mapped,cfg%thetav(53),vm,fs,source,aa,ss,st)
    convective=.5_wp*(fl+fs)
    viscous=.5_wp*matmul(sim%lagged_dissipation(:,:,53),mapped-uu(53,:))
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
