program test_partition_invariance
 use mod_kinds,only:wp
 use module_config,only:sim_config_t
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp)::u(6,3),r(3,6),ref(3,6),a(3,3,6),b(3,3,6),c(3,3,6),norm,fi(3),fo(3),base(2,2),p,rl,rg,alpha
 integer::i,k,st,n
 cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1
 cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp;cfg%mul=.5_wp;cfg%mug=.001_wp
 cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%Lr=12.546253993197219_wp
 cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
 cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%nCells=6;cfg%num_dummy_cells=0
 cfg%closure_pipe=2;cfg%closure_riser=2;cfg%junction_model='MATCHED'
 allocate(cfg%thetav(6),cfg%cell_width(6));cfg%thetav=0;cfg%cell_width=1
 do i=1,6
  p=2.1_wp+.001_wp*(i-1);alpha=.25_wp+.01_wp*(i-1)
  rl=1+(p-1)/cfg%nCl**2;rg=p/cfg%nCg**2
  u(i,:)=[(1-alpha)*rl,alpha*rg,.6_wp+.003_wp*(i-1)]
 enddo
 n=0
 do k=0,4,2
  cfg%N_pipe=k;call sim%init(cfg);call sim%prepare_dissipation(u,st)
  call check(st==0,'prepare')
  call sim%assemble(u,u,.001_wp,r,a,b,c,norm,fi,fo,st,.false.,base)
  call check(st==0,'residual')
  if(k==0)ref=r
  call check(maxval(abs(r-ref))<1e-12_wp,'physical partition invariance including removal')
  call check(maxval(abs(sum(r(:2,:),dim=2)-.001_wp*(fo(:2)-fi(:2))))<1e-15_wp,'global phase conservation')
  if(k>0)call check(all(base(:,1)==base(:,2)),'one common phase flux')
 enddo
 write(*,'(a,i0,a)')'PASS partition: ',n,' assertions'
contains
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label;error stop 1
  endif
  n=n+1
 end subroutine
end program
