program certificate
 !! Deterministic blockage certificate at a saved terminal state (51+51 mesh,
 !! offline_api). Emits KEY|values records; make_certificate.py classifies.
 !! Nothing in the solver, closures or tolerances is changed.
 !! Arguments: fixture folder, gas sound speed (m/s) of the run EOS.
 use offline_api
 use module_find,only:pressure_fun_uj
 use module_closure_law,only:lei_fechamento_or_ndim_simp
 use module_steady_state,only:solve_alpha_from_superficials
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),r0(nx),rr(nx),jj(nx,nx),jw(nx,nx),norm0,nn,fi(3),fo(3),dx(nx),amat(nx,nx),cgas
 real(wp)::lam,lam_ok,lam_bad,mid,xp(nx),xm(nx),best,lam_best,uu(nc,3),margin(nmargin),anorm,rcond,work(4*nx)
 real(wp)::smid(nc),h,fd(nx),err,errmax
 integer(c_int)::st,idx,stage,st_bad,idx_bad,stage_bad,stm
 integer::i,k,info,ipiv(nx),iwork(nx),cc,n_adm,which,ord(5),col
 character(512)::folder,arg
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P);cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp
 call sim%init(cfg);call sim%prepare_dissipation(old,cc)
 write(*,'(a,2i6)')'INIT|',st,cc
 if(st/=0.or.cc/=0)stop 2
 smid(lo)=.5_wp*cfg%cell_width(lo)
 do i=lo+1,hi
  smid(i)=smid(i-1)+.5_wp*(cfg%cell_width(i-1)+cfg%cell_width(i))
 enddo
 smid=smid*cfg%Lr
 ! ---------------- accepted (old) and last trial states
 do which=0,1
  call probe_state(int(which,c_int),x,sc,meta)
  call probe_eval(x,1_c_int,rr,jj,nn,fi,fo,st,idx,stage)
  write(*,'(a,i2,a,2es24.16,a,3i6,a,es14.6)')'STATE|',which,'|',meta(1:2),'|',st,idx,stage,'|',nn
  if(which==1)then
   r0=rr;norm0=nn;jw=jj
  endif
 enddo
 if(st/=0)stop 3
 ! dominant residual components of the trial state (scaled as the solver norm)
 do k=1,5
  ord(k)=maxloc(abs(r0)/sc,dim=1,mask=[(all(i/=ord(1:k-1)),i=1,nx)])
  write(*,'(a,i2,a,i5,a,i2,a,es14.6)')'RESID|',k,'|',lo+(ord(k)-1)/3,'|',mod(ord(k)-1,3)+1,'|',r0(ord(k))/sc(ord(k))
 enddo
 ! global Jacobian: 1-norm condition estimate
 amat=jw;anorm=maxval(sum(abs(amat),dim=1))
 call dgetrf(nx,nx,amat,nx,ipiv,info)
 if(info==0)call dgecon('1',nx,amat,nx,anorm,rcond,work,iwork,info)
 write(*,'(a,i4,a,es14.6,a,es14.6)')'JACOBIAN|',info,'|',anorm,'|',merge(1/rcond,-1._wp,info==0.and.rcond>0)
 dx=-r0;amat=jw;call dgesv(nx,1,amat,nx,ipiv,dx,nx,info)
 write(*,'(a,i4,a,es14.6,a,i5,a,i2)')'NEWTON|',info,'|',maxval(abs(dx)/sc),'|',lo+(maxloc(abs(dx)/sc,dim=1)-1)/3,'|', &
   mod(maxloc(abs(dx)/sc,dim=1)-1,3)+1
 if(info/=0)stop 4
 ! admissible part of the Newton ray and first boundary
 lam_ok=0;lam_bad=-1;best=norm0;lam_best=0;n_adm=0
 do k=0,400
  lam=min(1._wp,1e-12_wp*1.075_wp**k)
  xp=x+lam*dx;call probe_eval(xp,0_c_int,rr,jj,nn,fi,fo,st,idx,stage)
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
 write(*,'(a,i5,a,es14.6,a,es14.6,a,es14.6)')'RAY|',n_adm,'|',best,'|',lam_best,'|',best/norm0
 if(lam_bad<0)then
  write(*,'(a)')'BOUNDARY|NONE'
 else
  do k=1,200
   mid=lam_ok+.5_wp*(lam_bad-lam_ok)
   if(mid==lam_ok.or.mid==lam_bad)exit
   xp=x+mid*dx;call probe_eval(xp,0_c_int,rr,jj,nn,fi,fo,st,idx,stage)
   if(st==0)then
    lam_ok=mid
   else
    lam_bad=mid;st_bad=st;idx_bad=idx;stage_bad=stage
   endif
  enddo
  xp=x+lam_ok*dx;call probe_eval(xp,0_c_int,rr,jj,nn,fi,fo,stm,idx,stage)
  write(*,'(a,2es24.16,a,3i6,a,es14.6)')'BOUNDARY|',lam_ok,lam_bad,'|',st_bad,idx_bad,stage_bad,'|',nn/norm0
  if(idx_bad>=lo.and.idx_bad<=hi)write(*,'(a,i5,a,f12.4,a,f9.4,a,f12.6)')'LOCATION|',idx_bad,'|',smid(idx_bad),'|', &
    cfg%thetav(idx_bad)*180/acos(-1._wp),'|',cfg%cell_width(idx_bad)*cfg%Lr
 endif
 ! FD consistency of Jacobian columns of the cells around the failure
 errmax=0
 do i=max(lo,idx_bad-1),min(hi,idx_bad+2)
  do k=1,3
   col=3*(i-lo)+k;h=1e-6_wp*max(abs(x(col)),sc(col)*1e-3_wp)
   xp=x;xm=x;xp(col)=x(col)+h;xm(col)=x(col)-h
   call probe_eval(xp,0_c_int,fd,jj,nn,fi,fo,st,idx,stage)
   call probe_eval(xm,0_c_int,rr,jj,nn,fi,fo,stm,idx,stage)
   if(st==0.and.stm==0)then
    err=maxval(abs((fd-rr)/(2*h)-jw(:,col)))/max(1._wp,maxval(abs(jw(:,col))))
   else if(st==0)then
    err=maxval(abs((fd-r0)/h-jw(:,col)))/max(1._wp,maxval(abs(jw(:,col))))
   else if(stm==0)then
    err=maxval(abs((r0-rr)/h-jw(:,col)))/max(1._wp,maxval(abs(jw(:,col))))
   else
    err=-1
   endif
   errmax=max(errmax,err)
   write(*,'(a,i5,a,i2,a,2i6,a,es14.6)')'JACFD|',i,'|',k,'|',st,stm,'|',err
  enddo
 enddo
 ! states of the cells around the failure, accepted / trial / edge / refused
 do k=0,3
  select case(k)
  case(0);call probe_state(0_c_int,xp,sc,meta)
  case(1);xp=x
  case(2);xp=x+lam_ok*dx
  case(3);xp=x+lam_bad*dx
  end select
  call unpack(xp,uu)
  do i=max(lo,idx_bad-1),min(hi,idx_bad+1)
   call cell_record(k,i,uu)
  enddo
  if(stage_bad==2.and.idx_bad==ib)call face_record(k,uu)
 enddo
 call probe_margins(x,margin,stm)
 write(*,'(a,i6,a,es14.6,a,es14.6,a,i5)')'MARGIN|',stm,'|',margin(1),'|',minval(margin(3:)),'|',ib+minloc(margin(3:),dim=1)
contains
 subroutine cell_record(tag,i,u)
  integer,intent(in)::tag,i
  real(wp),intent(in)::u(nc,3)
  real(wp)::v(6),f(3),s(3),a(3,3),ds(3,3),wr(3),wi(3),ev(3,3)
  integer::code,cid,c2
  cid=merge(1,3,i<=ib)
  call sim%point_state(u(i,:),cfg%thetav(i),v,f,s,a,ds,code,.true.,cid)
  wr=ieee_value(0._wp,ieee_quiet_nan);wi=wr
  if(code==0)call dgeev_3x3(a,wr,wi,ev,c2)
  write(*,'(a,i2,a,i5,a,i2,a,i5,a,3es16.8,a,6es16.8,a,6es14.6)')'CELL|',tag,'|',i,'|',cid,'|',code,'|', &
    u(i,:)*[cfg%omega_rho,cfg%omega_rho,cfg%omega_rho*cfg%omega_u],'|', &
    v*[1._wp,cfg%omega_rho,cfg%omega_rho,cfg%omega_P,cfg%omega_u,cfg%omega_u],'|',wr*cfg%omega_u,wi*cfg%omega_u
  if(code/=0.and.cid==3)call root_map(tag,i,u(i,:))
 end subroutine
 subroutine face_record(tag,u)
  integer,intent(in)::tag
  real(wp),intent(in)::u(nc,3)
  real(wp)::raw(3),ps,vl(6),vr(6),vm(6),fl(3),fr(3),fs(3),mapped(3),rl,rg,a,s
  integer::code
  call junction_raw(u,raw,ps,vl,vr,vm,fl,fr,fs,mapped,code)
  s=cfg%omega_rho*cfg%omega_u
  if(code/=0)then
   write(*,'(a,i2,a,i6)')'FACE|',tag,'|',code;return
  endif
  rl=cfg%nrho_l0+(ps-cfg%nP_l0)/cfg%nCl**2;rg=ps/cfg%nCg**2
  a=solve_alpha_from_superficials(raw(1)/rl,raw(2)/rg,rl,rg,0._wp,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
    cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=code,closure_id=3)
  write(*,'(a,i2,a,i6,a,9es15.7,a,es15.7)')'FACE|',tag,'|',code,'|',.5_wp*(fl(1:2)+fs(1:2))*s, &
    (raw(1:2)-.5_wp*(fl(1:2)+fs(1:2)))*s,raw(1:2)*s,fl(2)*s,fs(2)*s,ps*cfg%omega_P,'|',a
 end subroutine
 subroutine root_map(tag,i,u)
  integer,intent(in)::tag,i
  real(wp),intent(in)::u(3)
  real(wp)::p,rl,rg,alpha,ug,ul,f,fprev,ugmax,fmin,ugmin,f0,f1
  integer::m,code,nsc,ndom
  p=pressure_fun_uj(u(1),u(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2;alpha=u(2)/rg
  if(u(3)<=0)then
   write(*,'(a,i2,a,i5,a)')'ROOTMAP|',tag,'|',i,'|U3_NONPOSITIVE';return
  endif
  ugmax=u(3)/u(2);nsc=0;ndom=0;fmin=huge(1._wp);fprev=0;f0=ieee_value(0._wp,ieee_quiet_nan);f1=f0
  do m=0,4000
   ug=ugmax*.5_wp*(1-cos(acos(-1._wp)*m/4000))
   ul=(u(3)-u(2)*ug)/u(1)
   call lei_fechamento_or_ndim_simp(alpha,rl,rg,ul,ug,cfg%thetav(i),cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
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
  write(*,'(a,i2,a,i5,a,es14.6,a,i3,a,i5,a,2es14.6,a,2es14.6)')'ROOTMAP|',tag,'|',i,'|',alpha,'|',nsc,'|',ndom,'|',f0,f1, &
    '|',fmin/alpha,ugmin*cfg%omega_u
 end subroutine
end program
