program replay_p0
 use offline_api
 implicit none
 real(wp)::x(nx),scale(nx),meta(6),rr(nx),jj(nx,nx),nn,fi(3),fo(3),dx(nx),amat(nx,nx),uu(nc,3)
 real(wp)::data(48),lambda,cgas,v(6),f(3),s(3),aj(3,3),sj(3,3),p,raw(3),ps,vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),mapped(3)
 real(wp)::rmat(3,np),fd(nx),jp(nx),rm(nx),jm(nx,nx),unusednorm,unusedfi(3),unusedfo(3),xp(nx),xm(nx),h,err
 integer(c_int)::st,idx,stage,code
 integer::i,k,which,ipiv(nx),info,largest,cc,selected,fluxmode
 character(512)::folder,arg
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 fluxmode=1;call get_command_argument(5,arg)
 if(len_trim(arg)>0)read(arg,*)fluxmode
 call probe_init(trim(folder),len_trim(folder),int(fluxmode,c_int),st)
 print *,'LEGACY_INIT_STATUS (uses 343 before explicit override)',st
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P);cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp
 call get_command_argument(3,arg)
 if(trim(arg)=='LOCAL_VOLUME')cfg%liquid_inlet_basis='LOCAL_VOLUME'
 call get_command_argument(4,arg)
 if(trim(arg)=='STEADY_TRACE')then
  cfg%boundary_method='STEADY_TRACE';cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
  allocate(cfg%vns(nc),cfg%face_positions(0:np));cfg%face_positions(0)=0
  do i=1,np
   cfg%face_positions(i)=cfg%face_positions(i-1)+cfg%cell_width(i+2)
   cfg%vns(i+2)=.5_wp*(cfg%face_positions(i-1)+cfg%face_positions(i))
  enddo
 endif
 call sim%init(cfg);call sim%prepare_dissipation(old,cc)
 print *,'ACTUAL_EOS_CG_M_S, PREPARE_STATUS',cgas,cc
 if(cc/=0)stop 2
 do which=0,1
  call probe_state(int(which,c_int),x,scale,meta)
  call probe_eval(x,1_c_int,rr,jj,nn,fi,fo,st,idx,stage)
  print *,'STATE old0/trial1,time_s,dt_s,status,index,stage,norm',which,meta(1:2),st,idx,stage,nn
  call unpack(x,uu)
  if(st/=0)then
   do i=lo,hi
    call sim%point_state(uu(i,:),cfg%thetav(i),v,f,s,aj,sj,cc,.false.,merge(1,3,i<=ib))
    if(cc/=0)print *,'POINT_FAILURE',i,cc,uu(i,:)*[cfg%omega_rho,cfg%omega_rho,cfg%omega_rho*cfg%omega_u]
   enddo
   cycle
  endif
  largest=maxloc(abs(rr),dim=1);selected=lo+(largest-1)/3
  print *,'MAX_RESIDUAL cell/component/value',selected,mod(largest-1,3)+1,rr(largest)
  do i=max(lo,selected-1),min(hi,selected+1)
   call sim%point_state(uu(i,:),cfg%thetav(i),v,f,s,aj,sj,cc,.false.,merge(1,3,i<=ib))
   print *,'SI_CELL i,closure,theta,status,U,alpha,rhol,rhog,p,ul,ug,phase_flux',i,merge(1,3,i<=ib),cfg%thetav(i),cc, &
     uu(i,:)*[cfg%omega_rho,cfg%omega_rho,cfg%omega_rho*cfg%omega_u], &
     v*[1._wp,cfg%omega_rho,cfg%omega_rho,cfg%omega_P,cfg%omega_u,cfg%omega_u],f(:2)*cfg%omega_rho*cfg%omega_u
  enddo
  call junction_raw(uu,raw,ps,vl,vr,vm,fl,fr,fs,mapped,cc)
  if(cc==0)print *,'JUNCTION status,conv_g,diss_g,raw_g_kg_m2s,p_Pa,ulug_left_right',cc, &
    .5_wp*(fl(2)+fs(2))*cfg%omega_rho*cfg%omega_u, &
    (raw(2)-.5_wp*(fl(2)+fs(2)))*cfg%omega_rho*cfg%omega_u,raw(2)*cfg%omega_rho*cfg%omega_u,ps*cfg%omega_P, &
    vl(5:6)*cfg%omega_u,vr(5:6)*cfg%omega_u
  err=0
  do k=1,3
   i=3*(selected-lo)+k;h=1e-5_wp*max(abs(x(i)),1e-8_wp);xp=x;xm=x;xp(i)=x(i)+h;xm(i)=x(i)-h
   call probe_eval(xp,0_c_int,fd,jm,unusednorm,unusedfi,unusedfo,st,idx,stage)
   if(st/=0)then
    print *,'FD_PLUS inadmissible',i,h,st,idx,stage;cycle
   endif
   call probe_eval(xm,0_c_int,rm,jm,unusednorm,unusedfi,unusedfo,st,idx,stage)
   if(st/=0)then
    print *,'FD_MINUS inadmissible',i,h,st,idx,stage;cycle
   endif
   err=max(err,maxval(abs((fd-rm)/(2*h)-jj(:,i))/max(1._wp,abs(jj(:,i)))))
  enddo
  print *,'LOCAL_JAC_ERROR same lagged operator',err
  amat=jj;dx=-rr
  call dgesv(nx,1,amat,nx,ipiv,dx,nx,info)
  print *,'LINEAR_STATUS, STEP_MAX',info,maxval(abs(dx))
  if(info/=0)cycle
  do k=0,18
   lambda=2._wp**(-k);xp=x+lambda*dx
   call probe_eval(xp,0_c_int,rr,jm,nn,fi,fo,st,idx,stage)
   print *,'PROPOSAL lambda,status,index,stage,residual',lambda,st,idx,stage,nn
  enddo
 enddo
end program
