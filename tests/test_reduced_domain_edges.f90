program test_reduced_domain_edges
 !! Regressions for the exact pure-phase reduction at the edge of its domain.
 !! A) Jacobian at the outlet backflow limit: the admissible equilibrium must
 !!    be differentiable from inside its domain; refusal codes are preserved.
 !! B) Inlet mass-flux search on the admissible subsonic, positive-pressure
 !!    branch; scale invariance; MASS x LOCAL_VOLUME equivalence; limits.
 !! Built with -ffpe-trap=invalid,zero,overflow by scripts/run_unit_tests.py.
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_single_phase
 implicit none
 integer,parameter::nc=16,lo=3,hi=14,n=12
 type(sim_config_t)::cfg
 real(wp)::u(nc,3),r(3,n),a(3,3,n),b(3,3,n),c(3,3,n),ar(3,3,n),br(3,3,n),cr(3,3,n)
 real(wp)::norm,fi(3),fo(3),base(2,2),rho0,err,rb,rb0,f(3),fr(3),ub(3),lam,kap,vb,rlv,ap,rp,rm,mp,mm,tg
 integer::phase,st,ns,k,checks
 real(wp),parameter::lams(3)=[1e-8_wp,1e-3_wp,1e6_wp],kaps(2)=[1e-3_wp,1e3_wp]
 real(wp),parameter::vbs(5)=[-.4_wp,-.05_wp,0._wp,.3_wp,.7_wp]
 character(len=12)::scheme
 checks=0
 ! ------------------------------------------------------------------ A
 do phase=1,2
  do k=1,2
   scheme='POINT';if(k==2)scheme='HYDROSTATIC'
   call setup(cfg);cfg%source_scheme=scheme
   rho0=cfg%nPs/cfg%nCg**2
   if(phase==1)rho0=cfg%nrho_l0+(cfg%nPs-cfg%nP_l0)/cfg%nCl**2
   u=0;u(:,phase)=rho0
   cfg%outlet_reservoir_phase=0
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.false.,base)
   write(*,'(a,i5,es12.4)')'A residual status, norm',st,norm
   call check(st==0.and.norm==0,'A equilibrium residual admissible and exact')
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.true.,base,ns)
   call check(st==0,'A equilibrium Jacobian exists without external reservoir data')
   call check(ns>0,'A outlet columns use the one-sided admissible stencil')
   ! Independent reference: same state, backflow of the same phase declared,
   ! so the residual is smooth across vb=0 and central differences apply.
   cfg%outlet_reservoir_phase=phase
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,ar,br,cr,norm,fi,fo,st,.true.,base,ns)
   call check(st==0.and.ns==0,'A reference Jacobian is central everywhere')
   err=max(maxval(abs(a-ar)/max(1._wp,abs(ar))),maxval(abs(b-br)/max(1._wp,abs(br))), &
     maxval(abs(c-cr)/max(1._wp,abs(cr))))
   write(*,'(a,i2,1x,a,es12.4)')'A phase/scheme one-sided vs central relative error',phase,trim(scheme),err
   call check(err<1e-7_wp,'A one-sided Jacobian equals smooth central reference')
   ! Physical backflow without data stays refused, with its own code.
   cfg%outlet_reservoir_phase=0
   u(hi:,3)=-1e-3_wp*rho0
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.false.,base)
   call check(st==PHASE_DATA_REQUIRED,'A effective backflow without reservoir composition refused (934)')
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.true.,base)
   call check(st==PHASE_DATA_REQUIRED,'A refusal reason is not converted to 935')
   cfg%outlet_reservoir_phase=3-phase
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.false.,base)
   call check(st==PHASE_TRANSITION_REQUIRED,'A different reservoir phase still requires transition (932)')
   ! A strictly interior state keeps central differences.
   cfg%outlet_reservoir_phase=0
   u(:,3)=1e-3_wp*rho0
   call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.true.,base,ns)
   call check(st==0.and.ns==0,'A interior outflow state uses central stencil only')
  enddo
 enddo
 ! ------------------------------------------------------------------ B (gas, c=1)
 call setup_unit(cfg)
 call inlet(2,1._wp,-.1_wp,-.1_wp,st,rb,f)
 call check(st==0.and.abs(rb-1)<=4*epsilon(1._wp),'B1 negative feed: uniform root found, no vacuum sampling')
 call check(f(2)==-.1_wp.and.abs(f(3)-(1*.01_wp+1))<=8*epsilon(1._wp),'B1 flux equals uniform state')
 call inlet(2,2._wp,.2_wp,.2_wp,st,rb,f)
 call check(st==0.and.abs(rb-2)<=8*epsilon(1._wp),'B2 positive feed: no log(0) under FPE traps')
 call inlet(2,1._wp,.3_wp,0._wp,st,rb,f)
 call check(st==0.and.abs(rb-exp(-.3_wp))<=4*epsilon(1._wp),'B zero feed: rarefaction root exp(-v/c)')
 call inlet(2,1._wp,0._wp,.5_wp,st,rb,f)
 vb=(rb-1)/sqrt(rb)
 call check(st==0.and.rb>1.and.abs(rb*vb-.5_wp)<=8*epsilon(1._wp).and.vb<1,'B compression root satisfies Hugoniot flux')
 ! Density scale invariance: (rho,m)->lambda*(rho,m) maps the root by lambda.
 call inlet(2,1._wp,.2_wp,.05_wp,st,rb0,f)
 call check(st==0,'B reference non-uniform root')
 do k=1,3
  lam=lams(k)
  call inlet(2,lam,.2_wp*lam,.05_wp*lam,st,rb,f)
  call check(st==0.and.abs(rb/lam-rb0)<=1e-14_wp*rb0,'B density-scale invariance')
  call inlet(2,lam,-.1_wp*lam,-.1_wp*lam,st,rb,f)
  call check(st==0.and.abs(rb/lam-1)<=8*epsilon(1._wp),'B scaled negative-feed uniform root')
 enddo
 ! Velocity scale invariance: (c,v,m)->kappa*(c,v,m) leaves the root.
 do k=1,2
  kap=kaps(k)
  cfg%nCg=kap
  call inlet(2,1._wp,.2_wp*kap,.05_wp*kap,st,rb,f)
  call check(st==0.and.abs(rb-rb0)<=1e-14_wp*rb0,'B velocity-scale invariance')
 enddo
 cfg%nCg=1
 ! MASS and LOCAL_VOLUME contracts describe the same physical state.
 do k=1,5
  vb=vbs(k)
  cfg%gas_inlet_basis='LOCAL_VOLUME';cfg%njg=vb
  call inlet(2,1._wp,.1_wp,0._wp,st,rlv,fr)
  call check(st==0,'B LOCAL_VOLUME state admissible')
  cfg%gas_inlet_basis='MASS';cfg%njg=0
  call inlet(2,1._wp,.1_wp,rlv*vb,st,rb,f)
  call check(st==0.and.abs(rb-rlv)<=1e-13_wp*rlv,'B MASS x LOCAL_VOLUME same boundary density')
  call check(maxval(abs(f-fr))<=1e-13_wp*max(1._wp,maxval(abs(fr))),'B MASS x LOCAL_VOLUME same flux')
 enddo
 ! Approach to the sonic limits of the admissible branch.
 ap=.8_wp;rp=(.5_wp*(ap+sqrt(ap*ap+4)))**2;mp=rp;rm=exp(-1.2_wp);mm=-rm
 call check(abs((sqrt(rp)-1/sqrt(rp))+.2_wp-1)<=4*epsilon(1._wp),'B closed-form compression sonic end')
 call inlet(2,1._wp,.2_wp,mp*(1-1e-10_wp),st,rb,f)
 vb=(rb-1)/sqrt(rb)+.2_wp
 call check(st==0.and.vb<1,'B just below compression sonic flux admissible')
 call inlet(2,1._wp,.2_wp,mp*(1+1e-10_wp),st,rb,f)
 call check(st==PHASE_DATA_REQUIRED,'B above compression sonic flux refused (934)')
 call inlet(2,1._wp,.2_wp,mp,st,rb,f)
 call check(st==PHASE_DATA_REQUIRED,'B exactly sonic flux refused')
 call inlet(2,1._wp,.2_wp,mm*(1-1e-10_wp),st,rb,f)
 call check(st==0.and.log(rb)+.2_wp>-1,'B just above rarefaction sonic flux admissible')
 call inlet(2,1._wp,.2_wp,mm*(1+1e-10_wp),st,rb,f)
 call check(st==PHASE_DATA_REQUIRED,'B below rarefaction sonic flux refused (934)')
 ! Unsupported interior states keep explicit codes.
 call inlet(2,1._wp,1._wp,.1_wp,st,rb,f)
 call check(st==PHASE_DATA_REQUIRED,'B sonic interior state needs more data')
 ! Liquid with the thermodynamic vacuum p=0 inside the rarefaction range.
 cfg%nCl=1;cfg%nrho_l0=1;cfg%nP_l0=.5_wp;cfg%nPs=2
 tg=.5_wp*log(.5_wp/.6_wp)
 call check(exp(-1._wp)*.6_wp<.5_wp,'B liquid case has the vacuum bound inside the rarefaction')
 call inlet(1,.6_wp,0._wp,tg*(1-1e-6_wp),st,rb,f)
 call check(st==0.and.rb>.5_wp,'B liquid root with positive pressure')
 call inlet(1,.6_wp,0._wp,tg*(1+1e-6_wp),st,rb,f)
 call check(st==PHASE_DATA_REQUIRED,'B liquid target below p=0 flux has no admissible root')
 call inlet(1,.45_wp,0._wp,0._wp,st,rb,f)
 call check(st==PHASE_INVALID,'B interior state at negative pressure invalid')
 write(*,'(a,i0,a)')'PASS reduced domain edges: ',checks,' assertions'
contains
 subroutine setup(cfg)
  type(sim_config_t),intent(inout)::cfg
  cfg=sim_config_t()
  cfg%single_phase_mode='REDUCED';cfg%source_scheme='POINT'
  cfg%num_dummy_cells=2;cfg%nCells=nc;cfg%N_pipe=6
  cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1;cfg%omega_t=10
  cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp
  cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=sqrt(287*293.15_wp)*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%Lr=10;cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%mul=0;cfg%mug=0;cfg%G=0
  allocate(cfg%thetav(nc),cfg%cell_width(nc),cfg%z_faces(0:n))
  cfg%thetav=.7_wp;cfg%cell_width=1._wp/n;cfg%z_faces=0
 end subroutine
 subroutine setup_unit(cfg)
  type(sim_config_t),intent(inout)::cfg
  call setup(cfg)
  cfg%omega_P=1;cfg%omega_rho=1;cfg%omega_u=1;cfg%nCg=1
 end subroutine
 subroutine inlet(ph,rho,mom,target,st,rb,f)
  integer,intent(in)::ph
  real(wp),intent(in)::rho,mom,target
  integer,intent(out)::st
  real(wp),intent(out)::rb,f(3)
  real(wp)::uu(nc,3),ubb(3)
  uu=0;uu(:,ph)=rho;uu(:,3)=mom
  cfg%mul=0;cfg%mug=0
  if(ph==1)cfg%mul=target
  if(ph==2)cfg%mug=target
  call single_phase_boundary(uu,cfg,ph,.true.,f,st,ubb)
  rb=ubb(ph)
 end subroutine
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label;error stop 1
  endif
  checks=checks+1
 end subroutine
end program
