program probe_ss3_branch
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_find,only:compute_primitive_variables
 use module_jacobian,only:DFduj
 use module_roe_solver,only:roe_absA_from_A
 use mod_void_helpers,only:alpha2gamma
 implicit none
 type(sim_config_t)::cfg
 integer,parameter::nc=52,np=48
 real(wp),parameter::pi=acos(-1.0_wp)
 real(wp)::old(nc,3),trial(nc,3),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np),rbase(3,np)
 real(wp)::norm,fin(3),fout(3),t,dt,norm_saved,dummy(3),theta,width
 real(wp)::v(6),f(3),src(3),jac(3,3),ds(3,3),vp(6),vm(6),fp(3),fm(3),sp(3),sm(3),tmp(3,3),tmp2(3,3)
 real(wp)::ga,gi,dl,dg,rel,reg,rem,fr,rm,j,mu,up(3),um(3),h,fd(3),ferror,serror,scale(3)
 real(wp)::eps_fd(3),bestf,bests,wave,diss(3,3),fmag,smag
 integer::iu,i,code,profile,k,ie,idx,reason,loc(2),imaxf,kmaxf,imaxs,kmaxs,st
 real(wp)::left(3),right(3),left0(3),lo,hi,mid,rlo,rhi,rtest,rlocal,momcrit,offset
 real(wp)::alocal(3,3),aright(3,3),amean(3,3),fface(3),fleft(3),fright(3),vlocal(6)
 real(wp)::fphysical0(3),rcrit,lambda,fnm
 integer::iteration,side,ko,uout
 character(40)::label
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
 cfg%tol=1e-13_wp;cfg%tola=1e-14_wp
 call get_command_argument(1,label)
 if(len_trim(label)==0)label='baseline'
 left0=trial(23,:);right=trial(24,:)
 lo=left0(3)-1e-4_wp;hi=left0(3)+1e-4_wp
 left=left0;left(3)=lo
 call decode(left,vlocal,fleft,alocal,rlo)
 left(3)=hi
 call decode(left,vlocal,fleft,alocal,rhi)
 if((rlo-2000._wp)*(rhi-2000._wp)>=0)then
  write(*,*)'critical momentum bracket failed',rlo,rhi
  stop 1
 endif
 do iteration=1,70
  mid=.5_wp*(lo+hi);left(3)=mid
  call decode(left,vlocal,fleft,alocal,rtest)
  if((rtest-2000._wp)*(rlo-2000._wp)>0)then
   lo=mid;rlo=rtest
  else
   hi=mid;rhi=rtest
  endif
 enddo
 momcrit=.5_wp*(lo+hi)
 write(*,'(a,es25.16)')'critical U3 for Re_g=2000 ',momcrit
 call decode(right,vlocal,fright,aright,rtest)
 open(newunit=uout,file='ss3_branch_'//trim(label)//'.csv',status='replace')
 write(uout,'(a)')'delta_U3,U3,Re_g,ul,ug,A11,A12,A13,A21,A22,A23,A31,A32,A33,wave_speed,Fg_physical,Fg_roe'
 do ko=5,13
  do side=-1,1,2
   offset=real(side,wp)*10._wp**(-ko);left=left0;left(3)=momcrit+offset
   call decode(left,vlocal,fleft,alocal,rlocal)
   amean=.5_wp*(alocal+aright)
   call roe_absA_from_A(amean,diss,lambda,code)
   if(code/=0)then
    write(*,*)'Roe failure',ko,side,code;cycle
   endif
   fface=.5_wp*(fleft+fright)-.5_wp*matmul(diss,right-left)
   write(uout,'(*(es24.16,:,","))')offset,left(3),rlocal,vlocal(5)*cfg%omega_u,vlocal(6)*cfg%omega_u, &
    alocal(1,:),alocal(2,:),alocal(3,:),lambda,fleft(2),fface(2)
  enddo
 enddo
 close(uout)
contains
 subroutine decode(u,v,f,aj,reg)
  real(wp),intent(in)::u(3)
  real(wp),intent(out)::v(6),f(3),aj(3,3),reg
  real(wp)::gamma,gammai,dg,pip
  integer::stat,index
  call compute_primitive_variables(u(1),u(2),u(3),-2*pi/180,cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0, &
   cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%njl,cfg%omega_u,cfg%omega_rho, &
   cfg%tol,cfg%tola,100,v(1),v(2),v(3),v(4),v(5),v(6),index,status=stat,closure_id=1)
  if(stat/=0)then
   write(*,*)'primitive failure',stat,u
   stop 2
  endif
  aj=DFduj(v(1),v(2),v(3),v(5),v(6),cfg%nCl,cfg%nCg,-2*pi/180,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
   cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%tol,status=stat,closure_id=1)
  if(stat/=0)stop 3
  pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
  f=[u(1)*v(5),u(2)*v(6),u(1)*v(5)**2+u(2)*v(6)**2+pip*v(4)]
  gamma=alpha2gamma(v(1),1e-13_wp,stat);gammai=sin(pi*gamma)/pi
  dg=v(1)*cfg%D_H/(1-gamma+gammai)
  reg=v(3)*cfg%omega_rho*abs(v(6))*cfg%omega_u*dg/cfg%mu_g
 end subroutine
end program
