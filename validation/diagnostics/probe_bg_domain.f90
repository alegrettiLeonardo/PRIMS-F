program probe_bg_domain
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_closure_law,only:CdUd_drift,lei_fechamento_or_ndim_simp
  use module_lapack_small,only:dgeev_3x3
  use module_find,only:pressure_fun_uj
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  integer::nc,np,ng,lo,hi,ib,physical_points_failed,eigen_failed,line_search_cell_failures,firstbad_common
  integer::min_status
  real(wp)::cl_si,cg_si,rho0_si,p0_si,ps_si,gl_si,gg_si,jl_si,jg_si,final_min_margin,nan
  real(wp)::candidate_limit,candidate_margin,base_limit,base_margin,base_min_ug,min_alpha,lin_error
  real(wp)::path_nodes(5),path_u(3),path_q(6),path_theta,path_limit,path_ug,path_margin
  real(wp)::path_min_margin
  real(wp)::interpolated_theta,interpolated_q(6),interpolated_limit,interpolated_margin
  integer::interpolated_status,line_search_junction_failures,junction_candidate_status
  integer::path_node_count,ju,jnode,path_status,prepare_status,first_face_failure,first_face_status
  real(wp)::new_q(6),physical_max_imag
  character(300)::config_file
  real(wp),parameter::pi=acos(-1.0_wp)
  real(wp),allocatable::old(:,:),trial(:,:),work(:,:),r(:,:),a(:,:,:),b(:,:,:),c(:,:,:),dx(:,:)
  real(wp)::norm,fin(3),fout(3),t,dt,norm_saved,dummy(3),theta,width
  real(wp)::v(6),f(3),src(3),jac(3,3),ds(3,3),eig(3),imag(3),vec(3,3),cd,ud,jmix,rem,mu
  real(wp)::fm(3),fp(3),trace(6),dmat(3,3),lambda,nt,mlimit,margin,uglimit
  real(wp),allocatable::vlo(:,:),linear_res(:,:),scaled(:,:)
  real(wp)::wave_imag,vel_min,dir(3),h,mp,mm,dmdt
  integer::iu,ic,i,k,code,idx,reason,st,es,ls,firstbad,mincell,badfaces
  character(300)::header,filename
  filename='last_rejected_state.dat';config_file='probe_config.dat'
  if(command_argument_count()>0)call get_command_argument(1,filename)
  if(command_argument_count()>1)call get_command_argument(2,config_file)
  open(newunit=iu,file=trim(config_file),status='old')
  read(iu,*)np,ng,cfg%N_pipe,cfg%closure_riser,cfg%quad_points,cfg%maxit
  read(iu,*)cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%Lr,cfg%Lp,cfg%D_H,cfg%eps,cfg%G
  read(iu,*)cfg%mu_l,cfg%mu_g,cfg%sigma,cl_si,cg_si,rho0_si,p0_si,ps_si
  read(iu,*)gl_si,gg_si,jl_si,jg_si,cfg%tol,cfg%tola,cfg%tol_res,cfg%tol_inc,cfg%jac_dissipation_fd_eps
  read(iu,*)cfg%flux_scheme,cfg%junction_model,cfg%dissipation_mode,cfg%gas_eos
  close(iu)
  cfg%omega_t=cfg%Lr/cfg%omega_u;cfg%AREA=pi*cfg%D_H**2/4
  cfg%nCl=cl_si*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=cg_si*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=rho0_si/cfg%omega_rho;cfg%nP_l0=p0_si/cfg%omega_P;cfg%nPs=ps_si/cfg%omega_P
  cfg%mul=gl_si/(cfg%omega_rho*cfg%omega_u);cfg%mug=gg_si/(cfg%omega_rho*cfg%omega_u)
  cfg%njl=jl_si/cfg%omega_u;cfg%njg=jg_si/cfg%omega_u
  cfg%closure_pipe=1;cfg%jac_dissipation_fd=.false.
  nc=np+2*ng;lo=ng+1;hi=ng+np;ib=lo+cfg%N_pipe-1
  cfg%nCells=nc;cfg%num_dummy_cells=ng
  allocate(cfg%cell_width(nc),cfg%thetav(nc))
  allocate(old(nc,3),trial(nc,3),work(nc,3),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),dx(3,np))
  allocate(vlo(6,nc),linear_res(3,np),scaled(3,np))
  nan=ieee_value(0._wp,ieee_quiet_nan);vlo=nan
  open(newunit=iu,file=trim(filename),status='old')
  read(iu,'(a)')header;read(iu,*)t,dt,norm_saved,reason
  do i=lo,hi
    read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
    cfg%thetav(i)=theta;cfg%cell_width(i)=width
  enddo
  close(iu)
  do i=1,ng
    old(i,:)=old(lo,:);trial(i,:)=trial(lo,:)
    old(hi+i,:)=old(hi,:);trial(hi+i,:)=trial(hi,:)
  enddo
  cfg%thetav(:ng)=cfg%thetav(lo);cfg%thetav(hi+1:)=cfg%thetav(hi)
  cfg%cell_width(:ng)=cfg%cell_width(lo);cfg%cell_width(hi+1:)=cfg%cell_width(hi)

  call sim%init(cfg)
  write(*,'(a,3es24.16,i8)')'snapshot time_s,dt_s,norm,status ',t*cfg%omega_t,dt*cfg%omega_t,norm_saved,reason
  vel_min=huge(1._wp);wave_imag=0;physical_points_failed=0;eigen_failed=0;mincell=0;min_alpha=1
  open(newunit=iu,file='point_states.csv',status='replace')
  write(iu,'(a)')'cell,theta,alpha,p_Pa,ul_m_s,ug_m_s,Re_mix,Re_BG,Cd,Ud_m_s,eig_min,eig_max,max_imag,m_min_ulzero_nd,m_margin_nd'
  do i=lo,hi
    call sim%point_state(trial(i,:),cfg%thetav(i),v,f,src,jac,ds,code)
    if(code/=0)then
      physical_points_failed=physical_points_failed+1
      write(*,*)'BASE POINT FAILURE',i,code;cycle
    endif
    vlo(:,i)=v;min_alpha=min(min_alpha,v(1))
    if(i>ib.and.v(5)<vel_min)then
      vel_min=v(5);mincell=i
    endif
    call dgeev_3x3(jac,eig,imag,vec,es)
    if(es/=0)then
      eigen_failed=eigen_failed+1;write(*,*)'EIG FAILURE',i,es
    endif
    wave_imag=max(wave_imag,maxval(abs(imag)))
    jmix=((1-v(1))*v(5)+v(1)*v(6))*cfg%omega_u
    mu=(1-v(1))*cfg%mu_l+v(1)*cfg%mu_g
    rem=(trial(i,1)+trial(i,2))*cfg%omega_rho*abs(jmix)*cfg%D_H/mu
    cd=0;ud=0;mlimit=0;margin=0
    if(i>ib)then
      call CdUd_drift(v(1),v(2),v(3),v(5),v(6),cfg%thetav(i),cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
                     cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,cd,ud,3,st)
      if(st/=0)write(*,*)'BG COEFF FAILURE',i,st
      call momentum_limit(v,cfg%thetav(i),mlimit,uglimit,st)
      if(st/=0)write(*,*)'LIMIT FAILURE',i,st
      margin=trial(i,3)-mlimit
    endif
    write(iu,'(i4,14(a,es24.16))')i,',',cfg%thetav(i),',',v(1),',',v(4)*cfg%omega_P, &
      ',',v(5)*cfg%omega_u,',',v(6)*cfg%omega_u,',',rem,',',v(2)*cfg%omega_rho*jmix*cfg%D_H/cfg%mu_l, &
      ',',cd,',',ud,',',minval(eig),',',maxval(eig),',',maxval(abs(imag)),',',mlimit,',',margin
  enddo
  close(iu)
  write(*,'(a,i4,2es24.16)')'MINIMUM LIQUID SPEED cell,SI,max spectrum imag ',mincell,vel_min*cfg%omega_u,wave_imag
  badfaces=0;path_node_count=0;path_min_margin=huge(1._wp);first_face_failure=0;first_face_status=0
  open(newunit=ju,file='roe_path_nodes.csv',status='replace')
  write(ju,'(a)')'left_cell,right_cell,path_fraction,theta_rad,point_status,alpha,p_Pa,U3_nd,limit_ulzero_nd,margin_nd,'// &
                'interpolated_theta_rad,interpolated_status,interpolated_limit_nd,interpolated_margin_nd'
  do i=lo,hi-1
    call sim%face_pair(trial(i,:),trial(i+1,:),cfg%thetav(i),cfg%thetav(i+1),fm,fp,st, &
                       junction_trace=trace,used_dissipation=dmat)
    if(st/=0)then
      badfaces=badfaces+1;write(*,*)'BASE FACE FAILURE',i,st
      if(first_face_failure==0)then
        first_face_failure=i;first_face_status=st
      endif
      ! Diagnose the actual three Gauss nodes plus endpoints for an interior
      ! riser face. The matched pipe/riser mapping is a separate operator.
      if(i>ib.and.cfg%quad_points==3)then
        path_nodes=[0._wp,.5_wp*(1._wp-.7745966692414834_wp),.5_wp, &
                        .5_wp*(1._wp+.7745966692414834_wp),1._wp]
        path_theta=.5_wp*(cfg%thetav(i)+cfg%thetav(i+1))
        do jnode=1,5
          path_u=trial(i,:)+path_nodes(jnode)*(trial(i+1,:)-trial(i,:))
          call sim%point_state(path_u,path_theta,path_q,f,src,jac,ds,path_status,.false.)
          path_q(4)=pressure_fun_uj(path_u(1),path_u(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
          path_q(2)=cfg%nrho_l0+(path_q(4)-cfg%nP_l0)/cfg%nCl**2
          path_q(3)=path_q(4)/cfg%nCg**2;path_q(1)=path_u(2)/path_q(3)
          call momentum_limit(path_q,path_theta,path_limit,path_ug,st)
          path_margin=path_u(3)-path_limit
          path_min_margin=min(path_min_margin,path_margin)
          path_node_count=path_node_count+1
          ! Diagnostic comparison only: production keeps theta fixed at its
          ! face mean. No alternative numerical operator is installed here.
          interpolated_theta=cfg%thetav(i)+path_nodes(jnode)*(cfg%thetav(i+1)-cfg%thetav(i))
          call sim%point_state(path_u,interpolated_theta,interpolated_q,f,src,jac,ds,interpolated_status,.false.)
          call momentum_limit(path_q,interpolated_theta,interpolated_limit,path_ug,st)
          interpolated_margin=path_u(3)-interpolated_limit
          write(ju,'(i6,a,i6,2(a,es26.16e3),a,i6,6(a,es26.16e3),a,i6,2(a,es26.16e3))') &
            i,',',i+1,',',path_nodes(jnode),',',path_theta,',',path_status, &
            ',',path_q(1),',',path_q(4)*cfg%omega_P,',',path_u(3),',',path_limit,',',path_margin, &
            ',',interpolated_theta,',',interpolated_status,',',interpolated_limit,',',interpolated_margin
        enddo
      endif
    endif
    if(i==ib)write(*,'(a,i5,6es24.16)')'JUNCTION status,trace ',st,trace
  enddo
  close(ju)
  write(*,'(a,i4)')'BASE FACE FAILURES ',badfaces
  call sim%prepare_dissipation(old,st)
  prepare_status=st
  write(*,'(a,i8)')'PREPARE LAGGED ',st
  if(st/=0)then
    call write_early_metrics()
    stop
  endif
  call sim%assemble(trial,old,dt,r,a,b,c,norm,fin,fout,code,.true.)
  write(*,'(a,i8,es24.16)')'BASE ASSEMBLE status,norm ',code,norm
  if(code/=0)stop
  call sim%ti%solve_block_tridiagonal_3x3(a,b,c,-r,dx,st)
  write(*,'(a,i8)')'NEWTON LINEAR SOLVE ',st
  linear_res=r
  do i=1,np
    linear_res(:,i)=linear_res(:,i)+matmul(b(:,:,i),dx(:,i))
    if(i>1)linear_res(:,i)=linear_res(:,i)+matmul(a(:,:,i),dx(:,i-1))
    if(i<np)linear_res(:,i)=linear_res(:,i)+matmul(c(:,:,i),dx(:,i+1))
  enddo
  lin_error=maxval(abs(linear_res))
  write(*,'(a,es24.16)')'LINEAR RESIDUAL MAX ',lin_error
  lambda=1;line_search_cell_failures=0;firstbad_common=0;line_search_junction_failures=0
  do ls=1,20
    work=trial
    do i=lo,hi
      work(i,:)=trial(i,:)+lambda*dx(:,i-lo+1)
    enddo
    call sim%assemble(work,old,dt,r,a,b,c,nt,fin,fout,code,.false.)
    firstbad=0;ic=0
    do i=lo,hi
      call sim%point_state(work(i,:),cfg%thetav(i),v,f,src,jac,ds,st,.false.)
      if(st/=0)then
        firstbad=i;ic=st;exit
      endif
    enddo
    if(firstbad>0)then
      line_search_cell_failures=line_search_cell_failures+1
      if(firstbad_common==0)firstbad_common=firstbad
      if(firstbad_common/=firstbad)firstbad_common=-1
    endif
    if(firstbad==0.and.code/=0)then
      call sim%face_pair(work(ib,:),work(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,junction_candidate_status, &
                        dissipation_matrix=sim%lagged_dissipation(:,:,ib))
      if(junction_candidate_status==1012)line_search_junction_failures=line_search_junction_failures+1
    endif
    write(*,'(a,i3,2es26.16e3,3i7)')'LINE SEARCH ls,lambda,norm,assembly,cell,status ',ls,lambda,nt,code,firstbad,ic
    lambda=lambda*.5_wp
  enddo
  i=mincell;v=vlo(:,i)
  call momentum_limit(v,cfg%thetav(i),mlimit,uglimit,st)
  write(*,'(a,2es24.16)')'LIMIT m_nd,ug_nd ',mlimit,uglimit
  write(*,'(a,3es24.16)')'CELL U ',trial(i,:)
  write(*,'(a,3es24.16)')'CELL Newton direction ',dx(:,i-lo+1)
  base_limit=mlimit;base_margin=trial(i,3)-mlimit;base_min_ug=v(6)*cfg%omega_u
  new_q=v
  work(i,:)=trial(i,:)+dx(:,i-lo+1)
  new_q(4)=pressure_fun_uj(work(i,1),work(i,2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  new_q(2)=cfg%nrho_l0+(new_q(4)-cfg%nP_l0)/cfg%nCl**2
  new_q(3)=new_q(4)/cfg%nCg**2;new_q(1)=work(i,2)/new_q(3)
  call momentum_limit(new_q,cfg%thetav(i),candidate_limit,uglimit,min_status)
  candidate_margin=work(i,3)-candidate_limit
  open(newunit=iu,file='metrics.txt',status='replace')
  write(iu,'(a,i0)')'physical_cells=',np
  write(iu,'(a,i0)')'pipe_cells=',cfg%N_pipe
  write(iu,'(a,i0)')'riser_cells=',np-cfg%N_pipe
  write(iu,'(a,i0)')'snapshot_rejection_status=',reason
  write(iu,'(a,es26.16e3)')'time_s=',t*cfg%omega_t
  write(iu,'(a,es26.16e3)')'dt_s=',dt*cfg%omega_t
  write(iu,'(a,es26.16e3)')'snapshot_norm=',norm_saved
  write(iu,'(a,es26.16e3)')'recomputed_norm=',norm
  write(iu,'(a,i0)')'point_failures=',physical_points_failed
  write(iu,'(a,i0)')'eigen_failures=',eigen_failed
  write(iu,'(a,es26.16e3)')'max_spectral_imag=',wave_imag
  write(iu,'(a,i0)')'face_failures=',badfaces
  write(iu,'(a,i0)')'first_failed_face=',first_face_failure
  write(iu,'(a,i0)')'first_failed_face_status=',first_face_status
  write(iu,'(a,i0)')'prepare_dissipation_status=',prepare_status
  write(iu,'(a,i0)')'path_nodes_examined=',path_node_count
  if(path_node_count>0)write(iu,'(a,es26.16e3)')'path_min_ulzero_margin_nd=',path_min_margin
  write(iu,'(a,es26.16e3)')'minimum_alpha=',min_alpha
  write(iu,'(a,i0)')'minimum_ul_cell_extended=',mincell
  write(iu,'(a,i0)')'minimum_ul_cell_physical=',mincell-ng
  write(iu,'(a,i0)')'minimum_ul_cell_riser=',mincell-ib
  write(iu,'(a,es26.16e3)')'critical_s_m=',(sum(cfg%cell_width(lo:mincell))-0.5_wp*cfg%cell_width(mincell))*cfg%Lr
  write(iu,'(a,es26.16e3)')'critical_theta_rad=',cfg%thetav(mincell)
  write(iu,'(a,es26.16e3)')'critical_alpha=',v(1)
  write(iu,'(a,es26.16e3)')'critical_pressure_Pa=',v(4)*cfg%omega_P
  write(iu,'(a,es26.16e3)')'critical_ul_m_s=',vel_min*cfg%omega_u
  write(iu,'(a,es26.16e3)')'critical_ug_m_s=',base_min_ug
  write(iu,'(a,es26.16e3)')'critical_momentum_nd=',trial(mincell,3)
  write(iu,'(a,es26.16e3)')'ulzero_momentum_nd=',base_limit
  write(iu,'(a,es26.16e3)')'ulzero_momentum_margin_nd=',base_margin
  write(iu,'(a,es26.16e3)')'newton_delta_momentum_nd=',dx(3,mincell-lo+1)
  write(iu,'(a,es26.16e3)')'linear_residual_max=',lin_error
  write(iu,'(a,i0)')'line_search_count=',20
  write(iu,'(a,i0)')'line_search_point_failures=',line_search_cell_failures
  write(iu,'(a,i0)')'line_search_junction_domain_failures=',line_search_junction_failures
  write(iu,'(a,i0)')'line_search_common_failure_cell=',firstbad_common
  write(iu,'(a,i0)')'candidate_limit_status=',min_status
  write(iu,'(a,es26.16e3)')'candidate_limit_momentum_nd=',candidate_limit
  write(iu,'(a,es26.16e3)')'candidate_limit_margin_nd=',candidate_margin
  close(iu)
contains
  subroutine write_early_metrics()
    open(newunit=iu,file='metrics.txt',status='replace')
    write(iu,'(a,i0)')'physical_cells=',np
    write(iu,'(a,i0)')'pipe_cells=',cfg%N_pipe
    write(iu,'(a,i0)')'riser_cells=',np-cfg%N_pipe
    write(iu,'(a,i0)')'snapshot_rejection_status=',reason
    write(iu,'(a,es26.16e3)')'time_s=',t*cfg%omega_t
    write(iu,'(a,es26.16e3)')'dt_s=',dt*cfg%omega_t
    write(iu,'(a,i0)')'point_failures=',physical_points_failed
    write(iu,'(a,i0)')'eigen_failures=',eigen_failed
    write(iu,'(a,es26.16e3)')'max_spectral_imag=',wave_imag
    write(iu,'(a,i0)')'face_failures=',badfaces
    write(iu,'(a,i0)')'first_failed_face=',first_face_failure
    write(iu,'(a,i0)')'first_failed_face_status=',first_face_status
    write(iu,'(a,i0)')'prepare_dissipation_status=',prepare_status
    write(iu,'(a,i0)')'path_nodes_examined=',path_node_count
    if(path_node_count>0)write(iu,'(a,es26.16e3)')'path_min_ulzero_margin_nd=',path_min_margin
    write(iu,'(a,es26.16e3)')'minimum_alpha=',min_alpha
    write(iu,'(a,i0)')'minimum_ul_cell_extended=',mincell
    write(iu,'(a,i0)')'minimum_ul_cell_physical=',mincell-ng
    write(iu,'(a,i0)')'minimum_ul_cell_riser=',mincell-ib
    write(iu,'(a,es26.16e3)')'critical_ul_m_s=',vel_min*cfg%omega_u
    close(iu)
  end subroutine write_early_metrics
  subroutine momentum_limit(v,theta,m,ug,status)
    real(wp),intent(in)::v(6),theta
    real(wp),intent(out)::m,ug
    integer,intent(out)::status
    real(wp)::lo,hi,mid,flo,fhi,fmid
    integer::it
    lo=1e-6_wp;hi=1
    call limit_value(v,theta,lo,flo,status)
    if(status/=0)return
    do it=1,40
      call limit_value(v,theta,hi,fhi,status)
      if(status/=0)return
      if(sign(1._wp,flo)/=sign(1._wp,fhi))exit
      hi=2*hi
    enddo
    if(it>40)then
      status=40;return
    endif
    do it=1,100
      mid=.5_wp*(lo+hi)
      call limit_value(v,theta,mid,fmid,status)
      if(status/=0)return
      if(abs(fmid)<1e-14_wp*v(1).or.hi-lo<16*epsilon(mid)*max(1._wp,abs(mid)))exit
      if(sign(1._wp,fmid)==sign(1._wp,flo))then
        lo=mid;flo=fmid
      else
        hi=mid
      endif
    enddo
    ug=mid;m=v(1)*v(3)*ug;status=0
  end subroutine momentum_limit
    subroutine limit_value(v,theta,x,f,status)
      real(wp),intent(in)::v(6),theta,x
      real(wp),intent(out)::f
      integer,intent(out)::status
      call lei_fechamento_or_ndim_simp(v(1),v(2),v(3),0._wp,x,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
                 cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,f,3,status)
    end subroutine limit_value
end program
