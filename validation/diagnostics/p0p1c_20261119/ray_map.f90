program ray_map
 !! Diagnostic only. Maps the admissible part of the Newton ray at a terminal
 !! rejected state and the constitutive roots at the first inadmissibility.
 !! Separates (a) root-search failure, (b) proven absence of a root in the
 !! declared constitutive domain, (c) globalization limits (no admissible
 !! descent), without changing the solver, closures or tolerances.
 !! Arguments: fixture folder, gas sound speed (m/s) of the run's EOS.
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_closure_law,only:lei_fechamento_or_ndim_simp
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),rr(nx),jj(nx,nx),norm0,nn,fi(3),fo(3),dx(nx),amat(nx,nx),cgas
 real(wp)::lam,lam_ok,lam_bad,mid,xp(nx),best,lam_best,uu(nc,3),raw(3),ps,vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),mapped(3)
 integer(c_int)::st,idx,stage,st_ok,st_bad,idx_bad,stage_bad
 integer::i,k,info,ipiv(nx),cc,n_adm,cell
 character(512)::folder,arg
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P);cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp
 call sim%init(cfg);call sim%prepare_dissipation(old,cc)
 if(st/=0.or.cc/=0)then
  print *,'INIT_FAILURE',st,cc;stop 2
 endif
 call probe_state(1_c_int,x,sc,meta)
 call probe_eval(x,1_c_int,rr,jj,norm0,fi,fo,st,idx,stage)
 print '(a,2es24.16,3i6,es14.6)','TRIAL time_s dt_s status index stage norm ',meta(1:2),st,idx,stage,norm0
 if(st/=0)stop 3
 amat=jj;dx=-rr;call dgesv(nx,1,amat,nx,ipiv,dx,nx,info)
 print '(a,i4,es14.6)','NEWTON linear_status max|dx| ',info,maxval(abs(dx))
 if(info/=0)stop 4
 ! -------- admissibility and descent along the ray, geometric grid
 lam_ok=0;lam_bad=-1;best=norm0;lam_best=0;n_adm=0
 open(10,file=trim(folder)//'/ray_scan.csv',status='replace')
 write(10,'(a)')'lambda,status,index,stage,norm'
 do k=0,400
  lam=1e-9_wp*1.06_wp**k
  if(lam>1)lam=1
  xp=x+lam*dx;call probe_eval(xp,0_c_int,rr,jj,nn,fi,fo,st,idx,stage)
  write(10,'(es24.16,3(",",i0),",",es24.16)')lam,st,idx,stage,nn
  if(st==0)then
   n_adm=n_adm+1
   if(nn<best)then
    best=nn;lam_best=lam
   endif
   if(lam_bad<0)lam_ok=lam
  else if(lam_bad<0)then
   lam_bad=lam;st_bad=st;idx_bad=idx;stage_bad=stage
  endif
  if(lam==1)exit
 enddo
 close(10)
 print '(a,i5,2es14.6,es12.4)','RAY admissible_samples best_norm best_lambda norm_ratio ',n_adm,best,lam_best,best/norm0
 if(lam_bad<0)then
  print '(a)','RAY whole sampled ray admissible: limitation is descent/globalization, not domain'
  stop
 endif
 ! -------- bisection of the first domain boundary
 do k=1,200
  mid=lam_ok+.5_wp*(lam_bad-lam_ok)
  if(mid==lam_ok.or.mid==lam_bad)exit
  xp=x+mid*dx;call probe_eval(xp,0_c_int,rr,jj,nn,fi,fo,st,idx,stage)
  if(st==0)then
   lam_ok=mid;st_ok=st
  else
   lam_bad=mid;st_bad=st;idx_bad=idx;stage_bad=stage
  endif
 enddo
 print '(a,2es24.16,3i6)','BOUNDARY lambda_admissible lambda_refused status index stage ',lam_ok,lam_bad,st_bad,idx_bad,stage_bad
 xp=x+lam_ok*dx;call probe_eval(xp,0_c_int,rr,jj,nn,fi,fo,st,idx,stage)
 print '(a,es14.6,a,es14.6,a,l2)','BOUNDARY norm_at_edge ',nn,' norm_ratio ',nn/norm0,' Armijo(1e-4) ',nn<norm0*(1-1e-4_wp*lam_ok)
 ! -------- constitutive map at the failing location, both sides of the edge
 cell=idx_bad
 do k=0,2
  lam=0;if(k==1)lam=lam_ok;if(k==2)lam=lam_bad
  call unpack(x+lam*dx,uu)
  if(stage_bad==2.and.cell==ib)then
   call junction_raw(uu,raw,ps,vl,vr,vm,fl,fr,fs,mapped,cc)
   print '(a,es14.6,i6,4es16.8)','FACE lambda status raw_Gg_kg_m2s conv_Gg diss_Gg p*_Pa ',lam,cc, &
     raw(2)*cfg%omega_rho*cfg%omega_u,.5_wp*(fl(2)+fs(2))*cfg%omega_rho*cfg%omega_u, &
     (raw(2)-.5_wp*(fl(2)+fs(2)))*cfg%omega_rho*cfg%omega_u,ps*cfg%omega_P
  else if(cell>=lo.and.cell<=hi)then
   call root_map(uu(cell,:),cfg%thetav(cell),lam)
  endif
 enddo
contains
 subroutine root_map(u,theta,lam)
  !! Bhagwat co-current interval u_g in [0,u3/u2] (u_l>=0 by momentum).
  real(wp),intent(in)::u(3),theta,lam
  real(wp)::p,rl,rg,alpha,ug,ul,f,fprev,ugmax,fmin,ugmin,f0,f1
  integer::m,code,nsc,ndom
  p=pressure_fun_uj(u(1),u(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2;alpha=u(2)/rg
  if(u(3)<=0)then
   print '(a,es14.6)','CELL u3<=0: no co-current state with both inventories positive, lambda ',lam;return
  endif
  ugmax=u(3)/u(2);nsc=0;ndom=0;fmin=huge(1._wp);fprev=0
  do m=0,4000
   ! Dense near both ends: mapped Chebyshev-like abscissae.
   ug=ugmax*.5_wp*(1-cos(acos(-1._wp)*m/4000))
   ul=(u(3)-u(2)*ug)/u(1)
   call lei_fechamento_or_ndim_simp(alpha,rl,rg,ul,ug,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
     cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,f,closure_id=3,status=code)
   if(code/=0)then
    ndom=ndom+1;cycle
   endif
   if(m==0)f0=f
   if(m==4000)f1=f
   if(m>0.and.f*fprev<0)nsc=nsc+1
   if(abs(f)<fmin)then
    fmin=abs(f);ugmin=ug
   endif
   fprev=f
  enddo
  print '(a,es14.6,a,es13.5,a,f9.6,a,es13.5,a,i3,a,i5,a,2es13.5,a,es13.5,a,es13.5)','CELL lambda ',lam,' alpha ',alpha, &
    ' theta ',theta,' ug_max_m_s ',ugmax*cfg%omega_u,' sign_changes ',nsc,' closure_domain_errors ',ndom, &
    ' Phi(ug=0),Phi(ug_max) ',f0,f1,' min|Phi|/alpha ',fmin/alpha,' at ug_m_s ',ugmin*cfg%omega_u
 end subroutine
end program
