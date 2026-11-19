program probe_ss3_bg_stall
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_closure_law,only:CdUd_drift,lei_fechamento_or_ndim_simp
  use module_lapack_small,only:dgeev_3x3
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  integer,parameter::nc=52,np=48
  real(wp),parameter::pi=acos(-1.0_wp)
  real(wp)::old(nc,3),trial(nc,3),work(nc,3),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),dx(3,np)
  real(wp)::norm,fin(3),fout(3),t,dt,norm_saved,dummy(3),theta,width
  real(wp)::v(6),f(3),src(3),jac(3,3),ds(3,3),eig(3),imag(3),vec(3,3),cd,ud,jmix,rem,mu
  real(wp)::fm(3),fp(3),trace(6),dmat(3,3),lambda,nt,mlimit,margin,uglimit
  real(wp)::vlo(6,nc),wave_imag,vel_min,dir(3),h,mp,mm,dmdt,linear_res(3,np),scaled(3,np)
  integer::iu,ic,i,k,code,idx,reason,st,es,ls,firstbad,mincell,badfaces
  character(300)::header,filename
  filename='ss3_bg_stall_snapshot.dat'
  if(command_argument_count()>0)call get_command_argument(1,filename)
  cfg%omega_u=.48_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp
  cfg%Lr=12.546253993197080_wp;cfg%Lp=57.4_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
  cfg%D_H=.05254_wp;cfg%AREA=pi*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
  cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
  cfg%nCl=1498._wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=343._wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=201000._wp/cfg%omega_P
  cfg%mul=279.174_wp/(cfg%omega_rho*cfg%omega_u);cfg%mug=.577934855150798_wp/(cfg%omega_rho*cfg%omega_u)
  cfg%njl=.28_wp/cfg%omega_u;cfg%njg=.48_wp/cfg%omega_u
  cfg%maxit=100;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%tol_res=1e-8_wp;cfg%tol_inc=1e-8_wp
  cfg%closure_pipe=1;cfg%closure_riser=3;cfg%jac_dissipation_fd=.false.
  cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%dissipation_mode='LAGGED'
  cfg%nCells=nc;cfg%num_dummy_cells=2;cfg%N_pipe=24;cfg%quad_points=3
  allocate(cfg%cell_width(nc),cfg%thetav(nc))
  open(newunit=iu,file=trim(filename),status='old')
  read(iu,'(a)')header;read(iu,*)t,dt,norm_saved,reason
  do i=3,50
    read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
    cfg%thetav(i)=theta;cfg%cell_width(i)=width
  enddo
  close(iu)
  old(1,:)=old(3,:);old(2,:)=old(3,:);old(51,:)=old(50,:);old(52,:)=old(50,:)
  trial(1,:)=trial(3,:);trial(2,:)=trial(3,:);trial(51,:)=trial(50,:);trial(52,:)=trial(50,:)
  cfg%thetav(1:2)=cfg%thetav(3);cfg%thetav(51:52)=cfg%thetav(50)
  cfg%cell_width(1:2)=cfg%cell_width(3);cfg%cell_width(51:52)=cfg%cell_width(50)
  call sim%init(cfg)
  write(*,'(a,3es24.16,i8)')'snapshot time_s,dt_s,norm,status ',t*cfg%omega_t,dt*cfg%omega_t,norm_saved,reason
  vel_min=huge(1._wp);wave_imag=0
  open(newunit=iu,file='ss3_bg_stall_cells.csv',status='replace')
  write(iu,'(a)')'cell,theta,alpha,p_Pa,ul_m_s,ug_m_s,Re_mix,Re_BG,Cd,Ud_m_s,eig_min,eig_max,max_imag,m_min_ulzero_nd,m_margin_nd'
  do i=3,50
    call sim%point_state(trial(i,:),cfg%thetav(i),v,f,src,jac,ds,code)
    if(code/=0)then
      write(*,*)'BASE POINT FAILURE',i,code;cycle
    endif
    vlo(:,i)=v
    if(i>=27.and.v(5)<vel_min)then
      vel_min=v(5);mincell=i
    endif
    call dgeev_3x3(jac,eig,imag,vec,es)
    if(es/=0)write(*,*)'EIG FAILURE',i,es
    wave_imag=max(wave_imag,maxval(abs(imag)))
    jmix=((1-v(1))*v(5)+v(1)*v(6))*cfg%omega_u
    mu=(1-v(1))*cfg%mu_l+v(1)*cfg%mu_g
    rem=(trial(i,1)+trial(i,2))*cfg%omega_rho*abs(jmix)*cfg%D_H/mu
    cd=0;ud=0;mlimit=0;margin=0
    if(i>=27)then
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
  badfaces=0
  do i=3,49
    call sim%face_pair(trial(i,:),trial(i+1,:),cfg%thetav(i),cfg%thetav(i+1),fm,fp,st, &
                       junction_trace=trace,used_dissipation=dmat)
    if(st/=0)then
      badfaces=badfaces+1;write(*,*)'BASE FACE FAILURE',i,st
    endif
    if(i==26)write(*,'(a,i5,6es24.16)')'JUNCTION status,trace ',st,trace
  enddo
  write(*,'(a,i4)')'BASE FACE FAILURES ',badfaces
  call sim%prepare_dissipation(old,st)
  write(*,'(a,i8)')'PREPARE LAGGED ',st
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
  write(*,'(a,es24.16)')'LINEAR RESIDUAL MAX ',maxval(abs(linear_res))
  lambda=1
  do ls=1,20
    work=trial
    do i=3,50
      work(i,:)=trial(i,:)+lambda*dx(:,i-2)
    enddo
    call sim%assemble(work,old,dt,r,a,b,c,nt,fin,fout,code,.false.)
    firstbad=0;ic=0
    do i=3,50
      call sim%point_state(work(i,:),cfg%thetav(i),v,f,src,jac,ds,st,.false.)
      if(st/=0)then
        firstbad=i;ic=st;exit
      endif
    enddo
    write(*,'(a,i3,2es24.16,3i7)')'LINE SEARCH ls,lambda,norm,assembly,cell,status ',ls,lambda,nt,code,firstbad,ic
    lambda=lambda*.5_wp
  enddo
  i=mincell;v=vlo(:,i)
  call momentum_limit(v,cfg%thetav(i),mlimit,uglimit,st)
  write(*,'(a,2es24.16)')'LIMIT m_nd,ug_nd ',mlimit,uglimit
  write(*,'(a,3es24.16)')'CELL U ',trial(i,:)
  write(*,'(a,3es24.16)')'CELL Newton direction ',dx(:,i-2)
contains
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
