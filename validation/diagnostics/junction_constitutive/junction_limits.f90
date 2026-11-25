program junction_limits
 !! P0.1 (diagnostic, NOT used by the solver): constitutive limits of the
 !! junction when the gas transport tends to zero, evaluated with the SAME
 !! closure routines the solver calls (module_closure_law, module_steady_state).
 !! Input: a saved terminal state (offline_api fixture, 51+51 mesh), the gas
 !! sound speed of the run EOS, and the junction_flux/junction_trace of the run.
 !! Records (SI):
 !!  CELL   two cells on each side of the junction: alpha, rho_g, u_g, j_g,
 !!         U_g=alpha rho_g, p, u_l, j_l, C0, Vgj (drift closures; for the
 !!         stratified pipe the apparent u_g/j and u_g-j are printed instead),
 !!         eigenvalues of the flux Jacobian.
 !!  FACE   for Roe/E1 x CENTER/STEADY on the SAME cell states: status, p*,
 !!         Gl*, Gg*, alpha*_pipe, alpha*_riser, F_conv=(F(x)+F(mapped))/2,
 !!         F_diss=-D(mapped-x)/2 (D = dissipation actually used), the raw
 !!         Roe vector F_conv+F_diss, and the transmitted vector (Gl*,Gg*,
 !!         momentum left, momentum right).
 !!  LIM1   alpha>0, u_g=0: Phi(alpha, u_l=j_l/(1-alpha), 0) on a logit grid;
 !!         sign changes = admissible stagnant-gas states.
 !!  LIM2   alpha->0 at fixed j_l: root u_g(alpha) of Phi (finite limit or not),
 !!         j_g=alpha u_g, and for drift closures C0, Vgj at that alpha.
 !!  LIM3   j_g -> 0+, 0, 0- at fixed j_l: production inversion
 !!         solve_alpha_from_superficials (status, number of roots, alpha, u_g).
 !! Pressure and j_l of each side are those of the adjacent cell. The riser
 !! inversion uses the face angle the solver uses (0 rad, catenary touchdown)
 !! and also the riser-cell angle.
 !! Arguments: fixture folder, c_gas (m/s), run junction_flux, run junction_trace.
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_closure_law,only:lei_fechamento_or_ndim_simp,CdUd_drift
 use module_steady_state,only:solve_alpha_from_superficials
 use module_boundary_trace,only:steady_junction_trace
 use module_jacobian,only:DFduj
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),cgas,fm(3),fp(3),tr(6),v(6),f(3),s(3),a(3,3),d(3,3),pipf
 real(wp)::wr(3),wi(3),vr(3,3),dmat(3,3),xl(3),yr(3),thl,thr,vl(6),vrr(6),fl(3),frr(3),vm(6),fs(3),mapped(3),astl
 real(wp)::cd,ud,jl,jg,fconv(3),fdiss(3),fu(3),cv(3)
 integer(c_int)::st
 integer::which,k,i,code,cid,ic,iside,lim_cell,lim_model
 real(wp)::lim_th,lim_jl
 character(512)::folder,arg
 character(16)::runflux,runtrace,runlimit='NONE'
 character(16)::flux(4)=['ROE         ','ROE         ','PHASE_UPWIND','PHASE_UPWIND']
 character(16)::trc(4)=['CENTER      ','STEADY      ','CENTER      ','STEADY      ']
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 runflux='ROE';runtrace='CENTER'
 if(command_argument_count()>2)call get_command_argument(3,runflux)
 if(command_argument_count()>3)call get_command_argument(4,runtrace)
 if(command_argument_count()>4)call get_command_argument(5,runlimit)
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P);cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp
 cfg%junction_flux=runflux;cfg%junction_trace=runtrace;cfg%junction_riser_limit=runlimit
 cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
 if(.not.allocated(cfg%face_positions))allocate(cfg%face_positions(0:hi-lo+1))
 if(.not.allocated(cfg%vns))allocate(cfg%vns(nc))
 cfg%face_positions(0)=0
 do i=lo,hi
  cfg%face_positions(i-lo+1)=cfg%face_positions(i-lo)+cfg%cell_width(i)
  cfg%vns(i)=cfg%face_positions(i-lo)+.5_wp*cfg%cell_width(i)
 enddo
 do i=1,lo-1
  cfg%vns(i)=cfg%vns(lo)-(lo-i)*cfg%cell_width(lo);cfg%vns(hi+i)=cfg%vns(hi)+i*cfg%cell_width(hi)
 enddo
 call sim%init(cfg);call sim%prepare_dissipation(old,code)
 if(code/=0)then
  write(*,'(a,i6)')'PREPARE_STATUS|',code;stop 2
 endif
 pipf=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
 write(*,'(a,a,a,f10.4,a,a,a,a)')'FIXTURE|',trim(folder),'|t_s',meta(1),'|run ',trim(runflux),'/',trim(runtrace)
 do which=0,1
  call probe_state(int(which,c_int),x,sc,meta);call unpack(x,uu)
  do i=ib-1,ib+2
   cid=merge(cfg%closure_pipe,cfg%closure_riser,i<=ib)
   call sim%point_state(uu(i,:),cfg%thetav(i),v,f,s,a,d,code,.true.,cid)
   wr=-999;wi=0;cd=-999;ud=-999
   if(code==0)then
    call dgeev_3x3(a,wr,wi,vr,ic)
    jl=(1-v(1))*v(5);jg=v(1)*v(6)
    if(cid/=1)then
     call CdUd_drift(v(1),v(2),v(3),v(5),v(6),cfg%thetav(i),cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
       cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,cd,ud,cid,ic)
     if(ic/=0)then;cd=-999;ud=-999;endif
    else
     ! stratified: apparent distribution parameter and drift (not closure inputs)
     cd=v(6)/(jl+jg);ud=(v(6)-(jl+jg))*cfg%omega_u
    endif
   endif
   write(*,'(a,i1,a,i4,a,i2,a,i4,a,f9.5,a,14es14.6)')'CELL|',which,'|',i,'|',cid,'|',code,'|',cfg%thetav(i)*180/acos(-1._wp),'|', &
     v(1),v(3)*cfg%omega_rho,v(6)*cfg%omega_u,v(1)*v(6)*cfg%omega_u,uu(i,2)*cfg%omega_rho,v(4)*cfg%omega_P, &
     v(5)*cfg%omega_u,(1-v(1))*v(5)*cfg%omega_u,cd,ud,wr*cfg%omega_u,wi(1)*cfg%omega_u
  enddo
  do k=1,4
   cfg%junction_flux=flux(k);cfg%junction_trace=trc(k);call sim%init(cfg)
   if(k==1)then
    cfg%junction_flux=runflux;cfg%junction_trace=runtrace;cfg%junction_riser_limit=runlimit;call sim%init(cfg)
    call sim%prepare_dissipation(old,code)
    cfg%junction_riser_limit='NONE'
    cfg%junction_flux=flux(k);cfg%junction_trace=trc(k);call sim%init(cfg)
    call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code, &
      junction_trace=tr,dissipation_matrix=sim%lagged_dissipation(:,:,ib),used_dissipation=dmat,left_cell=ib)
   else
    call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code,junction_trace=tr, &
      used_dissipation=dmat,left_cell=ib)
   endif
   ! Replicate the trace/mapping of the face operator to split F_conv/F_diss.
   fconv=-999;fdiss=-999;fu=-999
   xl=uu(ib,:);yr=uu(ib+1,:);thl=cfg%thetav(ib);thr=cfg%thetav(ib+1);ic=0
   if(trim(trc(k))=='STEADY')then
    call steady_junction_trace(uu(ib,:),cfg%thetav(ib),cfg,.true.,xl,thl,ic)
    if(ic==0)call steady_junction_trace(uu(ib+1,:),cfg%thetav(ib+1),cfg,.false.,yr,thr,ic)
   endif
   if(ic==0)call sim%point_state(xl,thl,vl,fl,s,a,d,ic,.false.,cfg%closure_pipe)
   if(ic==0)call sim%point_state(yr,thr,vrr,frr,s,a,d,ic,.false.,cfg%closure_riser)
   if(ic==0)then
    astl=solve_alpha_from_superficials(frr(1)/vrr(2),frr(2)/vrr(3),vrr(2),vrr(3),thl,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=ic,closure_id=cfg%closure_pipe,alpha_hint=vl(1))
   endif
   if(ic==0)then
    mapped=[(1-astl)*vrr(2),astl*vrr(3),frr(1)+frr(2)]
    call sim%point_state(mapped,thl,vm,fs,s,a,d,ic,.false.,cfg%closure_pipe)
   endif
   if(ic==0.and.all(ieee_is_finite(dmat)))then
    fconv=0.5_wp*(fl+fs);fdiss=-0.5_wp*matmul(dmat,mapped-xl);fu=fconv+fdiss
   endif
   cv=[cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u**2]
   write(*,'(a,i1,a,a,a,a,a,i6,a,5es14.6,a,3es14.6,a,3es14.6,a,3es14.6,a,4es14.6)')'FACE|',which,'|',trim(flux(k)),'|', &
     trim(trc(k)),'|',code,'|',tr(1)*cfg%omega_P,tr(2)*cv(1),tr(3)*cv(1),tr(4),tr(5), &
     '|Fconv',fconv*cv,'|Fdiss',fdiss*cv,'|Froe',fu*cv,'|Ftrans',fm(1)*cv(1),fm(2)*cv(1),fm(3)*cv(3),fp(3)*cv(3)
  enddo
  ! E2 (opt-in riser-trace limit) on the same states, Roe and E1, both traces.
  do k=1,4
   cfg%junction_flux=flux(k);cfg%junction_trace=trc(k);cfg%junction_riser_limit='ALPHA0';call sim%init(cfg)
   call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code,junction_trace=tr,left_cell=ib)
   write(*,'(a,i1,a,a,a,a,a,i6,a,5es14.6,a,es14.6)')'FACE_E2|',which,'|',trim(flux(k)),'|',trim(trc(k)),'|',code,'|', &
     tr(1)*cfg%omega_P,tr(2)*cfg%omega_rho*cfg%omega_u,tr(3)*cfg%omega_rho*cfg%omega_u,tr(4),tr(5),'|raw Gg',fm(2)*cfg%omega_rho*cfg%omega_u
   cfg%junction_riser_limit='NONE'
  enddo
  cfg%junction_flux=runflux;cfg%junction_trace=runtrace;cfg%junction_riser_limit=runlimit;call sim%init(cfg)
  call sim%prepare_dissipation(old,code)
  cfg%junction_riser_limit='NONE'
  if(which==0)then
   do iside=1,4
    call limits(iside)
   enddo
  endif
 enddo
contains
 subroutine limits(side)
  integer,intent(in)::side
  integer,parameter::m=4000
  real(wp)::p,rl,rg,jln,th,phi,phiprev,al,alprev,g,ugr,jgv,cdl,udl,aroot,ugroot
  integer::cell,model,j,c,cprev,nr,nsc,kk
  character(24)::lab
  select case(side)
  case(1);cell=ib;model=cfg%closure_pipe;th=cfg%thetav(ib);lab='PIPE_STRAT_cell'
  case(2);cell=ib+1;model=cfg%closure_riser;th=0;lab='RISER_BHAGWAT_face0'
  case(3);cell=ib+1;model=cfg%closure_riser;th=cfg%thetav(ib+1);lab='RISER_BHAGWAT_cell'
  case default;cell=ib+1;model=2;th=0;lab='RISER_BENDIKSEN_face0'
  end select
  p=pressure_fun_uj(uu(cell,1),uu(cell,2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
  call sim%point_state(uu(cell,:),cfg%thetav(cell),v,f,s,a,d,c,.false.,merge(cfg%closure_pipe,cfg%closure_riser,cell<=ib))
  jln=(1-v(1))*v(5)
  lim_cell=cell;lim_model=model;lim_th=th;lim_jl=jln
  write(*,'(a,a,a,i2,a,f9.4,a,f12.1,a,es14.6)')'LIMSIDE|',trim(lab),'|model',model,'|theta_deg',th*180/acos(-1._wp), &
    '|p_Pa',p*cfg%omega_P,'|j_l m/s',jln*cfg%omega_u
  ! LIM1: stagnant gas, alpha in (0,1)
  nsc=0;cprev=-1;phiprev=0;alprev=0
  do j=0,m
   al=1/(1+exp(-(-28+56*real(j,wp)/m)))
   call lei_fechamento_or_ndim_simp(al,rl,rg,jln/(1-al),0.0_wp,th,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
     cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,phi,closure_id=model,status=c)
   if(c==0.and.cprev==0)then
    if(phi*phiprev<0)then
     nsc=nsc+1
     write(*,'(a,a,a,2es14.6,a,es12.4)')'LIM1|',trim(lab),'|root alpha in',alprev,al,'|dPhi/dalpha sign',sign(1._wp,phi-phiprev)
    endif
   endif
   if(j==0.or.j==m.or.mod(j,500)==0)write(*,'(a,a,a,es12.4,a,i3,a,es14.6)')'LIM1MAP|',trim(lab),'|alpha',al,'|status',c,'|Phi',phi
   cprev=c;phiprev=phi;alprev=al
  enddo
  write(*,'(a,a,a,i3)')'LIM1|',trim(lab),'|stagnant-gas roots',nsc
  ! LIM2: alpha -> 0 at fixed j_l
  do kk=1,13
   al=10._wp**(-kk)
   call ug_root(al,ugroot,c)
   cdl=-999;udl=-999
   if(model/=1.and.c==0)then
    call CdUd_drift(al,rl,rg,jln/(1-al),ugroot,th,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma, &
      cfg%omega_u,cfg%omega_rho,cfg%tol,cdl,udl,model,nr)
   endif
   write(*,'(a,a,a,es10.2,a,i3,a,es16.8,a,es14.6,a,2es14.6)')'LIM2|',trim(lab),'|alpha',al,'|status',c,'|u_g m/s',ugroot*cfg%omega_u, &
     '|j_g m/s',al*ugroot*cfg%omega_u,'|C0,Vgj',cdl,udl
  enddo
  ! LIM3: production inversion at j_g -> 0+, 0, 0-
  do kk=-1,23
   if(kk<=10)then
    jgv=10._wp**(-kk-1)
   else if(kk==11)then
    jgv=0
   else
    jgv=-10._wp**(-(23-kk)-1)
   endif
   jgv=jgv/cfg%omega_u
   aroot=solve_alpha_from_superficials(jln,jgv,rl,rg,th,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma, &
     cfg%omega_u,cfg%omega_rho,cfg%tol,status=c,closure_id=model,nroots=nr)
   ugr=-999
   if(c==0)then
    if(aroot>0)ugr=jgv/aroot*cfg%omega_u
   endif
   write(*,'(a,a,a,es11.3,a,i3,a,i3,a,es16.8,a,es14.6)')'LIM3|',trim(lab),'|j_g m/s',jgv*cfg%omega_u,'|status',c,'|nroots',nr, &
     '|alpha',aroot,'|u_g m/s',ugr
  enddo
  g=0
 end subroutine
 subroutine ug_root(al,ugout,code)
  real(wp),intent(in)::al
  real(wp),intent(out)::ugout
  integer,intent(out)::code
  integer,parameter::m=20000
  real(wp)::p,rl,rg,jln,th,u0,u1,f0,f1,um,fmid
  integer::j,c0,c1,cm,it,nfound,cell,model
  cell=lim_cell;model=lim_model;th=lim_th
  p=pressure_fun_uj(uu(cell,1),uu(cell,2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
  jln=lim_jl
  code=21;ugout=-999;nfound=0
  u0=-5/cfg%omega_u
  call phi_u(al,rl,rg,jln,th,model,u0,f0,c0)
  do j=1,m
   u1=(-5+10*real(j,wp)/m)/cfg%omega_u
   call phi_u(al,rl,rg,jln,th,model,u1,f1,c1)
   if(c0/=0.or.c1/=0.or.nfound/=0)then
    u0=u1;f0=f1;c0=c1;cycle
   endif
   if(f0*f1<=0)then
    do it=1,200
     um=.5_wp*(u0+u1);call phi_u(al,rl,rg,jln,th,model,um,fmid,cm)
     if(cm/=0)exit
     if(f0*fmid<=0)then;u1=um;f1=fmid;else;u0=um;f0=fmid;endif
    enddo
    ugout=.5_wp*(u0+u1);code=0;nfound=1
   endif
   u0=u1;f0=f1;c0=c1
  enddo
 end subroutine
 subroutine phi_u(al,rl,rg,jln,th,model,ug,val,st_)
  real(wp),intent(in)::al,rl,rg,jln,th,ug
  integer,intent(in)::model
  real(wp),intent(out)::val
  integer,intent(out)::st_
  call lei_fechamento_or_ndim_simp(al,rl,rg,jln/(1-al),ug,th,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
    cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,val,closure_id=model,status=st_)
  if(st_==0)val=val/al   ! residual is proportional to alpha for drift closures; sign unchanged
  if(st_==0.and..not.ieee_is_finite(val))st_=-1
 end subroutine
end program
