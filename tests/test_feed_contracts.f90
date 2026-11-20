program test_feed_contracts
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_feed_contracts
 use module_boundary_conditions,only:boundary_flux
 use module_find,only:pressure_fun_uj
 implicit none
 type(feed_contract_t)::feed
 type(sim_config_t)::cfg
 real(wp)::c,ci,rl,rg,p,u(3),bc(3),f(3),gl,gg,scale,epsp,rplus,rminus
 integer::st,n,k
 n=0
 call gas_eos_coefficient('ISOTHERMAL',287.0_wp,293.15_wp,343.0_wp,ci,st)
 call check(st==0,'isothermal EOS')
 call close(ci**2,287*293.15_wp,2e-16_wp,'isothermal dp/drho=RT')
 call gas_eos_coefficient('BAROTROPIC',287.0_wp,293.15_wp,343.0_wp,c,st)
 call check(st==0.and.c==343,'barotropic coefficient retained')
 call close(101300/ci**2,1.20403094823_wp,3e-12_wp,'independent standard ideal density')
 call close(101300/c**2,.861035792909_wp,2e-12_wp,'declared barotropic density')
 epsp=10;rplus=(101300+epsp)/ci**2;rminus=(101300-epsp)/ci**2
 call close((rplus-rminus)/(2*epsp),1/(287*293.15_wp),1e-15_wp,'EOS pressure derivative')
 call resolve_feed('LEGACY','STANDARD',.547_wp,1.16_wp,FLOW_UNSET,FLOW_UNSET,997.05_wp,206000._wp,c,287._wp,293.15_wp,101300._wp,feed,st)
 call check(st==0.and.feed%liquid_mode=='MASS'.and.feed%gas_mode=='MASS','resolved legacy mass modes')
 call close(feed%liquid_mass,545.38635_wp,1e-15_wp,'legacy liquid reference-density flux')
 call close(feed%gas_mass,1.396675899947762_wp,1e-15_wp,'standard gas conversion')
 call resolve_feed('LOCAL_VOLUME','MASS',.3_wp,.5_wp,.2_wp,.6_wp,997.05_wp,206000._wp,c,287._wp,293.15_wp,101300._wp,feed,st)
 call check(st==FEED_CONFLICT,'liquid volume/mass conflict explicit')
 call resolve_feed('MASS','LOCAL_VOLUME',.3_wp,.5_wp,.2_wp,.6_wp,997.05_wp,206000._wp,c,287._wp,293.15_wp,101300._wp,feed,st)
 call check(st==FEED_CONFLICT,'gas volume/mass conflict explicit')
 call resolve_feed('MASS','MASS',.3_wp,.5_wp,FLOW_UNSET,.6_wp,997.05_wp,206000._wp,c,287._wp,293.15_wp,101300._wp,feed,st)
 call check(st==FEED_REQUIRED,'explicit mass requires value')
 call resolve_feed('MASS','MASS',.3_wp,.5_wp,-.2_wp,-.6_wp,997.05_wp,206000._wp,c,287._wp,293.15_wp,101300._wp,feed,st)
 call check(st==0.and.feed%liquid_mass==-.2_wp.and.feed%gas_mass==-.6_wp,'signed mass input retained; not a domain certificate')
 cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4
 cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
 cfg%closure_pipe=2;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp
 do k=1,3
  cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=real(k,wp)
  if(k==3)then
   cfg%omega_rho=100.;cfg%omega_P=1000000.
  endif
  cfg%nCg=ci*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=997.05_wp/cfg%omega_rho;cfg%nP_l0=101325/cfg%omega_P
  p=230000/cfg%omega_P;rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
  cfg%njl=.3_wp/cfg%omega_u;cfg%njg=.5_wp/cfg%omega_u
  cfg%liquid_inlet_basis='LOCAL_VOLUME';cfg%gas_inlet_basis='MASS';cfg%mug=.6_wp/(cfg%omega_rho*cfg%omega_u)
  call inlet_mass_fluxes(p,cfg,gl,gg,st)
  call check(st==0,'local volume target')
  call close(gl/rl*cfg%omega_u,.3_wp,2e-15_wp,'liquid scale invariance')
  call close(gg*cfg%omega_rho*cfg%omega_u,.6_wp,2e-15_wp,'gas mass scale invariance')
  u=[.75_wp*rl,.25_wp*rg,.6_wp/cfg%omega_u]
  call boundary_flux(u,0._wp,cfg,.true.,f,st,bc)
  call check(st==0,'local-volume boundary primitive')
  call close(f(1)/rl*cfg%omega_u,.3_wp,1e-11_wp,'executed inlet liquid volume')
  call close(f(2)*cfg%omega_rho*cfg%omega_u,.6_wp,1e-11_wp,'executed inlet gas mass')
  cfg%gas_inlet_basis='LOCAL_VOLUME'
  call boundary_flux(u,0._wp,cfg,.true.,f,st,bc)
  call check(st==0,'gas local-volume boundary')
  call close(f(2)/rg*cfg%omega_u,.5_wp,1e-11_wp,'gas local volume scale invariance')
 enddo
 write(*,'(a,i0,a)')'PASS feed/EOS: ',n,' assertions'
contains
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label;error stop 1
  endif
  n=n+1
 end subroutine
 subroutine close(a,b,tol,label)
  real(wp),intent(in)::a,b,tol
  character(*),intent(in)::label
  call check(abs(a-b)<=tol*max(1._wp,abs(b)),label)
 end subroutine
end program
