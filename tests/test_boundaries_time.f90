program test_boundaries_time
  use mod_kinds,only:wp,i4
  use module_config,only:sim_config_t
  use module_closure_law,only:CdUd_bendiksen
  use module_find,only:pressure_fun_uj
  use module_boundary_conditions,only:boundary_flux,apply_boundary_conditions,apply_ghost_cells
  use module_time_step,only:time_integrator_t
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(sim_config_t)::cfg
  type(time_integrator_t)::solver
  real(wp)::u(8,3),original(8,3),state(3),bc(3),f(3),f_first(3),f_second(3)
  real(wp)::p,alpha,rl,rg,ul,ug,cd,ud,theta,dt,resnorm
  real(wp)::a(3,3,3),b(3,3,3),c(3,3,3),d(3,3),x(3,3),expected(3,3)
  real(wp)::dense(9,9),rhs(9,1),matrix(3,3),inverse(3,3),eye(3,3)
  real(wp)::uu(6,3),un(6,3),ds(6),flux(6,3),source(6,3),r(6,3),balance(3)
  integer::st,k,i,j,row,col,piv(9),count
  interface
    subroutine dgesv(n,nrhs,aa,lda,ipiv,bb,ldb,info)
      import wp
      integer,intent(in)::n,nrhs,lda,ldb
      real(wp),intent(inout)::aa(lda,*),bb(ldb,*)
      integer,intent(out)::ipiv(*),info
    end subroutine
  end interface
  count=0
  cfg%num_dummy_cells=2;cfg%closure_pipe=2;cfg%closure_riser=2
  cfg%omega_u=1.0_wp;cfg%omega_rho=1000.0_wp;cfg%omega_P=1e5_wp
  cfg%nCg=sqrt(85000.0_wp)*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1.0_wp;cfg%nP_l0=1.0_wp;cfg%nPs=2.0_wp
  cfg%D_H=.0508_wp;cfg%AREA=acos(-1.0_wp)*cfg%D_H**2/4.0_wp
  cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
  cfg%mul=.5_wp;cfg%mug=.002_wp;cfg%tol=1e-12_wp;cfg%tola=1e-12_wp
  theta=.7_wp;p=2.3_wp;alpha=.2_wp;ul=1.0_wp
  rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
  call CdUd_bendiksen(alpha,rl,rg,ul,ul,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g, &
    cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,cd,ud,st)
  ug=(cd*(1-alpha)*ul+ud/cfg%omega_u)/(1-cd*alpha)
  state=[(1-alpha)*rl,alpha*rg,(1-alpha)*rl*ul+alpha*rg*ug]
  call boundary_flux(state,theta,cfg,.true.,f,st,bc)
  call check(st==0,'inlet flux status')
  call close(f(1),cfg%mul,1e-14_wp,'independent inlet liquid flux')
  call close(f(2),cfg%mug,1e-14_wp,'independent inlet gas flux')
  f_first=f
  p=pressure_fun_uj(bc(1),bc(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  call close(p,2.3_wp,1e-11_wp,'inlet pressure follows interior')
  cfg%mul=.4_wp;cfg%mug=.102_wp
  call boundary_flux(state,theta,cfg,.true.,f,st,bc)
  call check(st==0,'changed inlet flux status')
  call close(f(1),.4_wp,1e-14_wp,'same total but different liquid flux')
  call close(f(2),.102_wp,1e-14_wp,'same total but different gas flux')
  f_second=f
  cfg%mul=.5_wp;cfg%mug=.002_wp
  call boundary_flux(state,theta,cfg,.true.,f,st,bc)
  call check(maxval(abs(f-f_first))<1e-12_wp,'no cross-run inlet cache')
  call boundary_flux(state,theta,cfg,.false.,f,st,bc)
  call check(st==0,'outlet pressure state status')
  p=pressure_fun_uj(bc(1),bc(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  call close(p,cfg%nPs,1e-11_wp,'outlet exact prescribed pressure')
  call close(bc(3),state(3),1e-14_wp,'outlet keeps interior momentum')
  do i=1,8
    u(i,:)=state
  enddo
  original=u
  call apply_boundary_conditions(u,3,6,theta,theta,cfg,st)
  call check(st==0,'ghost boundary fill status')
  call check(all(u(3:6,:)==original(3:6,:)),'boundary fill never changes physical cells')
  u(4,1)=-1e-6_wp
  call apply_ghost_cells(u,2,st)
  call check(st==0.and.u(4,1)==-1e-6_wp,'ghost copying never repairs interior mass')
  ! A stagnant phase is still present. Pure limits require absent interior mass
  ! as well as zero prescribed flux; they are not inferred from the flux alone.
  cfg%closure_pipe=1;cfg%mug=0.0_wp
  call boundary_flux(state,-2.0_wp*acos(-1.0_wp)/180.0_wp,cfg,.true.,f,st,bc)
  call check(st==0.and.bc(2)>0.0_wp,'zero gas inlet flux retains stagnant gas inventory')
  call check(f(1)==cfg%mul.and.f(2)==0.0_wp,'stagnant gas keeps exact prescribed phase fluxes')
  state=[rl,0.0_wp,cfg%mul]
  call boundary_flux(state,-2.0_wp*acos(-1.0_wp)/180.0_wp,cfg,.true.,f,st,bc)
  call check(st==0.and.bc(2)==0.0_wp,'absent gas and zero gas flux use exact liquid-only limit')
  cfg%mul=0.0_wp;cfg%mug=.002_wp;state=[0.0_wp,rg,cfg%mug]
  call boundary_flux(state,-2.0_wp*acos(-1.0_wp)/180.0_wp,cfg,.true.,f,st,bc)
  call check(st==0.and.bc(1)==0.0_wp,'absent liquid and zero liquid flux use exact gas-only limit')
  cfg%closure_pipe=2;cfg%mul=.5_wp;cfg%mug=.002_wp
  matrix=0.0_wp;matrix(1,2)=1e8_wp;matrix(2,1)=1e-8_wp;matrix(3,3)=3.0_wp
  dt=solver%compute_time_step_lax_wendroff(.4_wp,.2_wp,matrix,status=st)
  call check(st==0,'spectral CFL status')
  call close(dt,.4_wp*.2_wp/3.0_wp,1e-12_wp,'CFL uses eigenvalues not maximum matrix entry')
  dt=solver%compute_time_step_lax_wendroff(.4_wp,.2_wp,matrix,acoustic_scale=.1_wp,status=st)
  call check(st/=0.and..not.ieee_is_finite(dt),'wave-speed underestimation rejected')
  matrix(1,2)=-1.0_wp;matrix(2,1)=1.0_wp
  dt=solver%compute_time_step_lax_wendroff(.4_wp,.2_wp,matrix,status=st)
  call check(st/=0.and..not.ieee_is_finite(dt),'nonhyperbolic acoustic CFL rejected')
  ! Full dense comparison of three coupled blocks.
  a=0.0_wp;b=0.0_wp;c=0.0_wp;dense=0.0_wp
  do k=1,3
    do i=1,3
      expected(i,k)=real(i+2*k,wp)/7.0_wp
      do j=1,3
        b(i,j,k)=.05_wp*real(i+j+k,wp)
        a(i,j,k)=-.03_wp*real(i+2*j,wp)
        c(i,j,k)=.02_wp*real(2*i+j,wp)
      enddo
      b(i,i,k)=b(i,i,k)+4.0_wp
    enddo
  enddo
  do k=1,3
    d(:,k)=matmul(b(:,:,k),expected(:,k))
    if(k>1)d(:,k)=d(:,k)+matmul(a(:,:,k),expected(:,k-1))
    if(k<3)d(:,k)=d(:,k)+matmul(c(:,:,k),expected(:,k+1))
    do i=1,3
      row=3*(k-1)+i
      do j=1,3
        col=3*(k-1)+j
        dense(row,col)=b(i,j,k)
        if(k>1)dense(row,col-3)=a(i,j,k)
        if(k<3)dense(row,col+3)=c(i,j,k)
      enddo
    enddo
  enddo
  rhs(:,1)=reshape(d,[9])
  call dgesv(9,1,dense,9,piv,rhs,9,st)
  call check(st==0,'independent full dense solve status')
  call solver%solve_block_tridiagonal_3x3(a,b,c,d,x,st)
  call check(st==0,'band solve status')
  call check(maxval(abs(reshape(x,[9])-rhs(:,1)))<1e-12_wp,'band solve agrees with full DGESV')
  call check(maxval(abs(x-expected))<1e-12_wp,'band solve recovers manufactured solution')
  ! Pivoting across blocks: leading diagonal block is singular, full matrix is not.
  a=0.0_wp;b=0.0_wp;c=0.0_wp
  eye=0.0_wp
  do i=1,3
    eye(i,i)=1.0_wp
  enddo
  c(:,:,1)=eye;a(:,:,2)=eye;b(:,:,3)=eye
  d(:,1)=expected(:,2);d(:,2)=expected(:,1);d(:,3)=expected(:,3)
  call solver%solve_block_tridiagonal_3x3(a,b,c,d,x,st)
  call check(st==0.and.maxval(abs(x-expected))<1e-12_wp,'cross-block pivoting without diagonal regularization')
  a=0.0_wp;b=0.0_wp;c=0.0_wp
  call solver%solve_block_tridiagonal_3x3(a,b,c,d,x,st)
  call check(st>0.and..not.all(ieee_is_finite(x)),'singular band matrix reported')
  call solver%invert_3x3(eye*1e-20_wp,inverse,st)
  call check(st==0,'small-scale nonsingular inverse accepted')
  call check(maxval(abs(inverse*1e-20_wp-eye))<1e-14_wp,'inverse has no epsilon regularization')
  uu=0.0_wp;un=0.0_wp;source=0.0_wp;ds=[1.0_wp,1.0_wp,2.0_wp,3.0_wp,4.0_wp,1.0_wp]
  do i=1,6
    flux(i,:)=[real(i,wp),real(i*i,wp),real(2*i,wp)]
  enddo
  call solver%compute_residual_vector(uu,un,ds,flux,source,2,5,.1_wp,r,resnorm,st)
  call check(st==0,'conservative residual status')
  balance=0.0_wp
  do i=2,5
    balance=balance+ds(i)*r(i,:)
  enddo
  call check(maxval(abs(balance-.1_wp*(flux(5,:)-flux(1,:))))<1e-13_wp,'residual telescopes to boundary fluxes')
  write(*,'(a,i0,a)')'PASS boundaries/time: ',count,' assertions'
contains
  subroutine check(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    if(.not.condition)then
      write(*,'(a)')'FAIL: '//label
      error stop 1
    endif
    count=count+1
  end subroutine
  subroutine close(value,target,tolerance,label)
    real(wp),intent(in)::value,target,tolerance
    character(*),intent(in)::label
    if(.not.all(ieee_is_finite([value,target])))then
      write(*,*)value,target
      call check(.false.,label)
    endif
    if(abs(value-target)>tolerance*max(1.0_wp,abs(value),abs(target)))write(*,*)label,value,target
    call check(abs(value-target)<=tolerance*max(1.0_wp,abs(value),abs(target)),label)
  end subroutine
end program test_boundaries_time
