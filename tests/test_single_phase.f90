program test_single_phase
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_single_phase
 implicit none
 integer,parameter::ncell=12,lo=3,hi=14
 type(sim_config_t)::cfg
 real(wp)::u(16,3),old(16,3),xp(16,3),xm(16,3),r(3,12),rp(3,12),rm(3,12),a(3,3,12),b(3,3,12),c(3,3,12)
 real(wp)::norm,nbad,fi(3),fo(3),base(2,2),x(3,3,12),y(3,3,12),z(3,3,12),dt,rhot,c2,zz,delta,avg,rnd,h,d(3),pred(3),f(3),mu,ub(3)
 integer::phase,ng,i,k,st,t,col,j,checks
 checks=0
 cfg%single_phase_mode='REDUCED';cfg%source_scheme='HYDROSTATIC'
 cfg%num_dummy_cells=2;cfg%nCells=16;cfg%N_pipe=6
 cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1;cfg%omega_t=10
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp
 cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=sqrt(287*293.15_wp)*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%Lr=10;cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp
 cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%mul=0;cfg%mug=0
 allocate(cfg%thetav(16),cfg%cell_width(16),cfg%z_faces(0:12))
 cfg%thetav=.7_wp;cfg%cell_width=1._wp/12
 do i=0,12
  cfg%z_faces(i)=cfg%Lr*sin(.7_wp)*i/12
 enddo
 dt=.001_wp
 do phase=1,2
  cfg%outlet_reservoir_phase=phase
  c2=1498._wp**2
  rhot=1+(cfg%nPs-1)/cfg%nCl**2
  if(phase==2)then
   c2=287*293.15_wp;rhot=cfg%nPs/cfg%nCg**2
  endif
  do ng=0,3
   cfg%G=4.905_wp*ng;u=0
   do i=lo,hi
    zz=.5_wp*(cfg%z_faces(i-lo)+cfg%z_faces(i-lo+1));delta=.5_wp*cfg%G*cfg%Lr*sin(.7_wp)/(12*c2)
    avg=1
    if(delta/=0)avg=sinh(delta)/delta
    u(i,phase)=rhot*exp(cfg%G*(cfg%z_faces(12)-zz)/c2)*avg
   enddo
   u(1,:)=u(lo,:);u(2,:)=u(lo,:);u(15,:)=u(hi,:);u(16,:)=u(hi,:);old=u
   call check(global_single_phase(u,cfg)==phase,'exact inventory reduction')
   cfg%source_scheme='HYDROSTATIC'
   call single_phase_assemble(u,old,dt,cfg,phase,r,a,b,c,norm,fi,fo,st,.false.,base)
   call check(st==0,'hydrostatic admissibility')
   rnd=64*epsilon(1._wp)*(c2/cfg%omega_u**2)*maxval(u(:,phase))*dt/minval(cfg%cell_width)
   call check(norm<=rnd,'hydrostatic residual bounded by conditioned roundoff')
   call check(all(r(3-phase,:)==0),'absent phase exactly zero')
   call check(maxval(abs(sum(r(:2,:),dim=2)/12-dt*(fo(:2)-fi(:2))))<1e-14_wp,'phase integral telescopes')
   call check(all(base(:,1)==base(:,2)),'common phase flux')
   cfg%source_scheme='POINT'
   call single_phase_assemble(u,old,dt,cfg,phase,r,x,y,z,nbad,fi,fo,st,.false.,base)
   call check(st==0,'unbalanced comparison admissible')
   if(ng>0)call check(nbad>100*norm,'hydrostatic reconstruction discriminates source error')
   write(*,'(a,2i4,3es20.10)')'HYDRO phase,g-index,WB,POINT,bound ',phase,ng,norm,nbad,rnd
  enddo
  cfg%G=9.81_wp;cfg%source_scheme='HYDROSTATIC';u(:,3)=.001_wp*u(:,phase);old=u
  call single_phase_assemble(u,old,dt,cfg,phase,r,a,b,c,norm,fi,fo,st,.true.,base)
  call check(st==0,'moving reduced Jacobian')
  do i=lo,hi,hi-lo
   do t=1,2
    k=phase;if(t==2)k=3
    h=2e-5_wp*max(abs(u(i,k)),1e-2_wp);xp=u;xm=u;xp(i,k)=u(i,k)+h;xm(i,k)=u(i,k)-h
    call single_phase_assemble(xp,old,dt,cfg,phase,rp,x,y,z,nbad,fi,fo,st,.false.,base)
    call check(st==0,'independent plus state')
    call single_phase_assemble(xm,old,dt,cfg,phase,rm,x,y,z,nbad,fi,fo,st,.false.,base)
    call check(st==0,'independent minus state')
    do j=max(lo,i-1),min(hi,i+1)
     col=j-lo+1;pred=b(:,k,col)
     if(j==i-1)pred=c(:,k,col)
     if(j==i+1)pred=a(:,k,col)
     d=(rp(:,col)-rm(:,col))/(2*h)
     call check(maxval(abs(d-pred)/max(1._wp,abs(pred)))<1e-5_wp,'independent complete Jacobian including boundary')
    enddo
   enddo
  enddo
  cfg%source_scheme='POINT'
  do i=1,16
   u(i,phase)=rhot;u(i,3)=-.1_wp*rhot
  enddo
  call single_phase_face(u,u,cfg,phase,7,f,st)
  call check(st==0.and.f(phase)<0,'negative signed flux is not inadmissible mass')
  cfg%outlet_reservoir_phase=0
  call single_phase_boundary(u,cfg,phase,.false.,f,st,ub)
  call check(st==PHASE_DATA_REQUIRED,'unknown external composition cannot be inferred from domain inventory')
  cfg%outlet_reservoir_phase=3-phase
  call single_phase_boundary(u,cfg,phase,.false.,f,st)
  call check(st==PHASE_TRANSITION_REQUIRED,'different reservoir phase needs conservative transition model')
  cfg%outlet_reservoir_phase=phase
  call single_phase_boundary(u,cfg,phase,.false.,f,st)
  call check(st==0.and.f(phase)<0,'subsonic pure backflow with explicitly declared reservoir phase')
  u(:,3)=-2*sqrt(c2)*u(:,phase)
  call single_phase_boundary(u,cfg,phase,.false.,f,st)
  call check(st==PHASE_DATA_REQUIRED,'supersonic incoming boundary requires more data')
 enddo
 u=0;u(:,1)=1+(cfg%nPs-1)/cfg%nCl**2;cfg%mug=.001_wp
 call single_phase_assemble(u,u,dt,cfg,1,r,a,b,c,norm,fi,fo,st,.false.,base)
 call check(st==PHASE_TRANSITION_REQUIRED,'new phase injection requires transition model, no velocity invented')
 write(*,'(a,i0,a)')'PASS single-phase: ',checks,' assertions'
contains
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label;error stop 1
  endif
  checks=checks+1
 end subroutine
end program
