program test_hydrostatic_outlet
 !! Discrete hydrostatic rest of the exact pure-phase reduction with G/=0:
 !! the separator face must see vb=0 exactly at the hydrostatic initial state
 !! (no arithmetic backflow), physical perturbations must give the signed
 !! outlet velocity of the legacy Riemann relations to full relative accuracy,
 !! and real backflow without declared reservoir composition stays refused.
 !! Grid: liquid/gas x 5 gravities x 4 meshes x 2 geometries (catenary riser,
 !! straight uphill pipe) x 2 nondimensional scales.
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_single_phase
 implicit none
 type(sim_config_t)::cfg
 real(wp),allocatable::u(:,:),r(:,:),a(:,:,:),b(:,:,:),c(:,:,:)
 real(wp)::norm,fi(3),fo(3),base(2,2),f(3),ub(3),rho0,vb_legacy,rhof,fl,fr,c2,rhob,v,eps_rel,g
 real(wp),parameter::gs(5)=[0._wp,4.905_wp,9.81_wp,19.62_wp,98.1_wp],omegas(2)=[1._wp,3.7_wp]
 integer,parameter::ns(4)=[4,8,16,64]
 integer::phase,ig,in,geo,io,st,n,top,checks,refused_at_rest,bad,k
 checks=0;refused_at_rest=0;bad=0
 do phase=1,2
  do ig=1,5
   do in=1,4
    do geo=1,2
     do io=1,2
      n=ns(in);g=gs(ig)
      call setup(n,geo,g,omegas(io))
      top=cfg%nCells-cfg%num_dummy_cells
      allocate(u(cfg%nCells,3),r(3,2*n),a(3,3,2*n),b(3,3,2*n),c(3,3,2*n))
      u=0;u(:,phase)=1;u(:,3)=0
      call initialize_single_phase_hydrostatic(cfg,phase,u);call ghosts()
      cfg%outlet_reservoir_phase=0
      call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.false.,base)
      if(st/=0)refused_at_rest=refused_at_rest+1
      call check(st==0,'hydrostatic rest admissible without reservoir data')
      call single_phase_boundary(u,cfg,phase,.false.,f,st,ub)
      call check(st==0.and.f(phase)==0.and.ub(3)==0,'separator face velocity exactly zero at discrete rest')
      call single_phase_assemble(u,u,1e-3_wp,cfg,phase,r,a,b,c,norm,fi,fo,st,.true.,base)
      call check(st==0,'complete Jacobian at rest')
      ! Physical perturbations of the top cell: compressed -> outflow; expanded -> backflow.
      do k=-1,1,2
       eps_rel=real(k,wp)
       call ghosts();u(top,phase)=u(top,phase)*(1+1e-6_wp*eps_rel);call ghosts()
       call coeff(phase,c2,rhob);call hydro_factors(cfg,phase,top,fl,fr,st)
       rhof=u(top,phase)*fr;v=u(top,3)/rhof
       if(rhob>rhof)then
        vb_legacy=v-sqrt(c2)*(rhob-rhof)/sqrt(rhob*rhof)
       else
        vb_legacy=v+sqrt(c2)*log(rhof/rhob)
       endif
       cfg%outlet_reservoir_phase=phase
       call single_phase_boundary(u,cfg,phase,.false.,f,st,ub)
       call check(st==0,'declared same-phase reservoir: signed outlet state exists')
       call check(sign(1._wp,ub(3))==eps_rel,'sign of outlet velocity follows the physical perturbation')
       if(abs(ub(3)/ub(phase)-vb_legacy)>1e-6_wp*abs(vb_legacy))bad=bad+1
       call check(abs(ub(3)/ub(phase)-vb_legacy)<=1e-6_wp*abs(vb_legacy),'deviation form equals legacy Riemann relation')
       cfg%outlet_reservoir_phase=0
       call single_phase_boundary(u,cfg,phase,.false.,f,st,ub)
       if(eps_rel<0)call check(st==PHASE_DATA_REQUIRED,'real backflow without composition refused (934)')
       if(eps_rel>0)call check(st==0,'physical outflow admissible')
       call initialize_single_phase_hydrostatic(cfg,phase,u);call ghosts()
      enddo
      deallocate(u,r,a,b,c)
     enddo
    enddo
   enddo
  enddo
 enddo
 write(*,'(a,i0,a,i0,a)')'rest states refused: ',refused_at_rest,'  deviation mismatches: ',bad,''
 write(*,'(a,i0,a)')'PASS hydrostatic outlet: ',checks,' assertions'
contains
 subroutine setup(n,geo,g,om)
  integer,intent(in)::n,geo
  real(wp),intent(in)::g,om
  integer::i
  cfg=sim_config_t()
  cfg%single_phase_mode='REDUCED';cfg%source_scheme='HYDROSTATIC'
  cfg%num_dummy_cells=2;cfg%N_pipe=n;cfg%nCells=2*n+4
  cfg%omega_rho=997.05_wp*om;cfg%omega_P=101325._wp*om**2;cfg%omega_u=1.16_wp*om;cfg%omega_t=10
  cfg%nrho_l0=997.05_wp/cfg%omega_rho;cfg%nP_l0=101325._wp/cfg%omega_P;cfg%nPs=206000._wp/cfg%omega_P
  cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=sqrt(287*293.15_wp)*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=g
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%mul=0;cfg%mug=0
  if(geo==1)then
   cfg%Lp=57.4_wp;cfg%Lr=12.54625399319708_wp;cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
  else
   cfg%Lp=40._wp;cfg%Lr=20._wp;cfg%CA=1e30_wp;cfg%beta=-acos(-1._wp)/6
  endif
  allocate(cfg%face_positions(0:2*n),cfg%cell_width(cfg%nCells),cfg%thetav(cfg%nCells))
  do i=0,n
   cfg%face_positions(i)=cfg%Lp/cfg%Lr*i/n
   cfg%face_positions(n+i)=cfg%Lp/cfg%Lr+real(i,wp)/n
  enddo
  do i=1,2*n
   cfg%cell_width(i+2)=cfg%face_positions(i)-cfg%face_positions(i-1)
  enddo
  cfg%cell_width(1:2)=cfg%cell_width(3);cfg%cell_width(2*n+3:)=cfg%cell_width(2*n+2);cfg%thetav=0
  call configure_pipe_riser_elevation(cfg)
 end subroutine
 subroutine ghosts()
  u(1,:)=u(3,:);u(2,:)=u(3,:);u(cfg%nCells-1,:)=u(cfg%nCells-2,:);u(cfg%nCells,:)=u(cfg%nCells-2,:)
 end subroutine
 subroutine coeff(phase,c2,rhob)
  integer,intent(in)::phase
  real(wp),intent(out)::c2,rhob
  real(wp)::pip
  pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
  c2=pip*cfg%nCg**2;rhob=cfg%nPs/cfg%nCg**2
  if(phase==1)then
   c2=pip*cfg%nCl**2;rhob=cfg%nrho_l0+(cfg%nPs-cfg%nP_l0)/cfg%nCl**2
  endif
 end subroutine
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label,' phase,G,n,geo,omega ',phase,g,n,geo,omegas(io)
   error stop 1
  endif
  checks=checks+1
 end subroutine
end program
