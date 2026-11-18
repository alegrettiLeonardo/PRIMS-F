program probe_ss3_stall
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 use mod_void_helpers,only:alpha2gamma
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 integer,parameter::nc=52,np=48
 real(wp),parameter::pi=acos(-1.0_wp)
 real(wp)::old(nc,3),trial(nc,3),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),rbase(3,np)
 real(wp)::norm,fin(3),fout(3),t,dt,norm_saved,dummy(3),theta,width
 real(wp)::v(6),f(3),src(3),jac(3,3),ds(3,3),vp(6),vm(6),fp(3),fm(3),sp(3),sm(3),tmp(3,3),tmp2(3,3)
 real(wp)::ga,gi,dl,dg,rel,reg,rem,fr,rm,j,mu,up(3),um(3),h,fd(3),ferror,serror,scale(3)
 real(wp)::eps_fd(3),bestf,bests,wave,diss(3,3),fmag,smag
 integer::iu,i,code,profile,k,ie,idx,reason,loc(2),imaxf,kmaxf,imaxs,kmaxs,st
 character(200)::header
 eps_fd=[1e-5_wp,1e-6_wp,1e-7_wp]
 cfg%omega_u=.48_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp
 cfg%Lr=12.546253993197093_wp;cfg%Lp=57.4_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
 cfg%D_H=.05254_wp;cfg%AREA=pi*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
 cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
 cfg%nCl=1498._wp*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%nCg=343._wp*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=201000._wp/cfg%omega_P
 cfg%mul=279.174_wp/(cfg%omega_rho*cfg%omega_u);cfg%mug=.577934855150798_wp/(cfg%omega_rho*cfg%omega_u)
 cfg%njl=.28_wp/cfg%omega_u;cfg%njg=.48_wp/cfg%omega_u
 cfg%maxit=100;cfg%tol_res=1e-8_wp;cfg%tol_inc=1e-8_wp
 cfg%closure_pipe=1;cfg%closure_riser=2;cfg%jac_dissipation_fd=.true.
 cfg%junction_model='MATCHED';cfg%flux_scheme='ROE';cfg%nCells=nc;cfg%num_dummy_cells=2;cfg%N_pipe=24
 allocate(cfg%cell_width(nc),cfg%thetav(nc))
 open(newunit=iu,file='ss3_stall_snapshot.dat',status='old')
 read(iu,'(a)')header
 read(iu,*)t,dt,norm_saved,reason
 do i=3,50
  read(iu,*)idx,theta,width,old(i,:),trial(i,:),dummy
  cfg%thetav(i)=theta;cfg%cell_width(i)=width
 enddo
 close(iu)
 old(1,:)=old(3,:);old(2,:)=old(3,:);old(51,:)=old(50,:);old(52,:)=old(50,:)
 trial(1,:)=trial(3,:);trial(2,:)=trial(3,:);trial(51,:)=trial(50,:);trial(52,:)=trial(50,:)
 cfg%thetav(1:2)=cfg%thetav(3);cfg%thetav(51:52)=cfg%thetav(50)
 cfg%cell_width(1:2)=cfg%cell_width(3);cfg%cell_width(51:52)=cfg%cell_width(50)
 write(*,'(a,3es20.10,i8)')'snapshot time_s,dt_s,norm,status ',t*cfg%omega_t,dt*cfg%omega_t,norm_saved,reason
 do profile=1,4
  cfg%tol=1e-11_wp;cfg%tola=1e-10_wp
  if(profile==2.or.profile==4)cfg%tol=1e-13_wp
  if(profile==3.or.profile==4)cfg%tola=1e-14_wp
  call sim%init(cfg)
  call sim%assemble(trial,old,dt,r,a,b,c,norm,fin,fout,code,.true.)
  write(*,'(a,i2,a,3es20.10,i8)')'profile ',profile,' tol,tola,norm,status ',cfg%tol,cfg%tola,norm,code
  if(code/=0)cycle
  if(profile==1)rbase=r
  loc=maxloc(abs(r))
  write(*,'(a,2i6,es20.10)')'max raw residual component,cell,value ',loc(1),loc(2)+2,r(loc(1),loc(2))
  write(*,'(a,es20.10)')'max residual difference vs baseline ',maxval(abs(r-rbase))
  if(profile==1)then
   open(newunit=iu,file='ss3_stall_regimes.csv',status='replace')
   write(iu,'(a)')'cell,theta,alpha,p_Pa,ul,ug,Re_l,Re_g,Re_mix,Fr_mix,r1,r2,r3'
   do i=3,50
    call sim%point_state(trial(i,:),cfg%thetav(i),v,f,src,jac,ds,code)
    if(code/=0)then
     write(*,*)'point error',i,code;cycle
    endif
    ga=alpha2gamma(v(1),1e-13_wp,st);gi=sin(pi*ga)/pi
    dl=(1-v(1))*cfg%D_H/ga;dg=v(1)*cfg%D_H/(1-ga+gi)
    rel=v(2)*cfg%omega_rho*abs(v(5))*cfg%omega_u*dl/cfg%mu_l
    reg=v(3)*cfg%omega_rho*abs(v(6))*cfg%omega_u*dg/cfg%mu_g
    rm=(1-v(1))*v(2)+v(1)*v(3);j=((1-v(1))*v(5)+v(1)*v(6))*cfg%omega_u
    mu=(1-v(1))*cfg%mu_l+v(1)*cfg%mu_g
    rem=rm*cfg%omega_rho*abs(j)*cfg%D_H/mu;fr=abs(j)/sqrt(cfg%G*cfg%D_H)
    write(iu,'(i4,12(a,es24.16))')i,',',cfg%thetav(i),',',v(1),',',v(4)*cfg%omega_P,',',v(5)*cfg%omega_u, &
      ',',v(6)*cfg%omega_u,',',rel,',',reg,',',rem,',',fr,',',r(1,i-2),',',r(2,i-2),',',r(3,i-2)
   enddo
   close(iu)
  endif
  do ie=1,3
   sim%cfg%jac_dissipation_fd_eps=eps_fd(ie)
   call sim%assemble(trial,old,dt,r,a,b,c,norm,fin,fout,code,.true.)
   write(*,'(a,es13.5,a,es20.10,i8)')'assemble FDstep ',eps_fd(ie),' norm,status ',norm,code
   bestf=0;bests=0;imaxf=0;kmaxf=0;imaxs=0;kmaxs=0
   do i=3,50
    call sim%point_state(trial(i,:),cfg%thetav(i),v,f,src,jac,ds,code)
    if(code/=0)cycle
    do k=1,3
     h=eps_fd(ie)*max(abs(trial(i,k)),merge(1e-5_wp,1e-2_wp,k==2))
     up=trial(i,:);um=up;up(k)=up(k)+h;um(k)=um(k)-h
     call sim%point_state(up,cfg%thetav(i),vp,fp,sp,tmp,tmp2,code,.false.)
     if(code/=0)cycle
     call sim%point_state(um,cfg%thetav(i),vm,fm,sm,tmp,tmp2,code,.false.)
     if(code/=0)cycle
     fd=(fp-fm)/(2*h);ferror=maxval(abs(fd-jac(:,k))/max(1._wp,abs(fd)))
     fd=(sp-sm)/(2*h);serror=maxval(abs(fd-ds(:,k))/max(1._wp,abs(fd)))
     if(ferror>bestf)then
      bestf=ferror;imaxf=i;kmaxf=k
     endif
     if(serror>bests)then
      bests=serror;imaxs=i;kmaxs=k
     endif
    enddo
   enddo
   write(*,'(a,es20.10,2i5)')'point flux Jacobian max relative error/cell/component ',bestf,imaxf,kmaxf
   write(*,'(a,es20.10,2i5)')'point source Jacobian max relative error/cell/component ',bests,imaxs,kmaxs
  enddo
 enddo
end program
