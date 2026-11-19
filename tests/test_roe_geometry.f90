program test_roe_geometry
  !! Independent geometric chain rule and frozen SS3 domain regression.
  use mod_kinds, only: wp,i4
  use module_config, only: sim_config_t
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use module_steady_state, only: solve_alpha_from_superficials
  use module_matrix_roe, only: calculate_roe_matrix
  use module_roe_solver, only: roe_absA_from_A
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp)::left(3),right(3),x(3),fl(3),fr(3),fm(3),fp(3),fplus(3),fminus(3)
  real(wp)::avg(3,3),other(3,3),diss(3,3),used(3,3),geom(3),err,no_geom,t,theta,h,weight,speed
  real(wp)::theta_left,theta_right,v(6),src(3),jac(3,3),ds(3,3),jl(3,3),jr(3,3),kl(3,3),kr(3,3)
  real(wp)::up(3),um(3),fmp(3),fmm(3),fpp(3),fpm(3),derivative_error
  integer::count,code,j,k,cid
  count=0
  cfg%omega_u=.48_wp;cfg%omega_rho=997.05_wp;cfg%omega_P=101325.0_wp
  cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0=1;cfg%nP_l0=1
  cfg%D_H=.05254_wp;cfg%AREA=acos(-1.0_wp)*cfg%D_H**2/4
  cfg%eps=4.6e-5_wp;cfg%G=9.81_wp;cfg%Lr=12.546253993197080_wp
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp
  cfg%tol=1e-12_wp;cfg%tola=1e-13_wp;cfg%quad_points=5
  cfg%closure_pipe=1;cfg%junction_model='MATCHED';cfg%flux_scheme='ROE'
  cfg%jac_dissipation_fd=.true.;cfg%jac_dissipation_fd_eps=1e-6_wp
  sim%cfg=cfg
  theta_left=.4_wp;theta_right=.45_wp
  do cid=2,3
    sim%cfg%closure_riser=cid
    call state_from_fluxes(2.2_wp,.5_wp,.001_wp,theta_left,left)
    call state_from_fluxes(2.21_wp,.505_wp,.00101_wp,theta_right,right)
    avg=path_matrix(left,right,theta_left,theta_right,5,code)
    call check(code==0,'augmented path on a fixed closure')
    other=path_matrix(right,left,theta_right,theta_left,5,code)
    call check(code==0.and.maxval(abs(other-avg)/max(1.0_wp,abs(avg)))<1e-10_wp,'path reversal invariance')
    other=path_matrix(left,right,theta_left,theta_right,5,code,2.0_wp,7.0_wp)
    call check(code==0.and.maxval(abs(other-avg)/max(1.0_wp,abs(avg)))<1e-11_wp,'augmented path reparameterization')
    call physical_flux(left,theta_left,fl)
    call physical_flux(right,theta_right,fr)
    ! Composite Simpson uses independent flux evaluations at fixed U and
    ! centered theta differences, rather than the production Jacobian.
    geom=0;h=1e-5_wp
    do j=0,64
      t=real(j,wp)/64;theta=theta_left+t*(theta_right-theta_left)
      x=left+t*(right-left)
      call physical_flux(x,theta+h,fplus)
      call physical_flux(x,theta-h,fminus)
      weight=2; if(mod(j,2)==1)weight=4
      if(j==0.or.j==64)weight=1
      geom=geom+weight*(fplus-fminus)/(2*h)
    enddo
    geom=geom*(theta_right-theta_left)/(3*64)
    err=maxval(abs(matmul(avg,right-left)+geom-(fr-fl))/max(1.0_wp,abs(fr-fl)))
    no_geom=maxval(abs(matmul(avg,right-left)-(fr-fl))/max(1.0_wp,abs(fr-fl)))
    write(*,'(a,i0,2(a,es14.6))')'GEOMETRY closure=',cid,' chain_error=',err,' omitted_term_error=',no_geom
    call check(err<2e-7_wp,'augmented chain rule includes integral F_theta')
    call check(no_geom>100*max(err,1e-10_wp),'ordinary A*dU secant is not asserted for variable geometry')
    call roe_absA_from_A(avg,diss,speed,code)
    call check(code==0.and.all(ieee_is_finite(diss)),'path matrix has a real usable eigenstructure')
    call sim%face_pair(left,right,theta_left,theta_right,fm,fp,code,jl,jr,kl,kr,used_dissipation=used)
    call check(code==0,'variable-angle face and derivatives')
    call check(maxval(abs(fm-fp))==0,'one common conservative face flux')
    call check(maxval(abs(used-diss)/max(1.0_wp,abs(diss)))<1e-11_wp,'face uses augmented path dissipation')
    call check(maxval(abs(fm-(.5_wp*(fl+fr)-.5_wp*matmul(diss,right-left))))<1e-11_wp,'endpoint physical fluxes retain geometry')
    do k=1,3
      h=3e-6_wp*max(abs(left(k)),merge(1e-5_wp,1e-2_wp,k==2))
      up=left;um=left;up(k)=up(k)+h;um(k)=um(k)-h
      call sim%face_pair(up,right,theta_left,theta_right,fmp,fpp,code)
      call check(code==0,'variable-angle left plus perturbation')
      call sim%face_pair(um,right,theta_left,theta_right,fmm,fpm,code)
      call check(code==0,'variable-angle left minus perturbation')
      derivative_error=maxval(abs(jl(:,k)-(fmp-fmm)/(2*h))/max(1.0_wp,abs(jl(:,k))))
      call check(derivative_error<2e-4_wp,'variable-angle complete face derivative')
    enddo
    ! Omitted optional theta_right preserves the original constant-angle API.
    avg=path_matrix(left,right,theta_left,theta_left,5,code)
    other=calculate_roe_matrix(left,right,cfg%nCl,cfg%nCg,theta_left,cfg%nrho_l0,cfg%nP_l0, &
      cfg%D_H,cfg%AREA,cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,0.0_wp, &
      cfg%omega_u,cfg%omega_rho,cfg%omega_P,cfg%tol,cfg%tola,5_i4,status=code,closure_id=cid)
    call check(code==0.and.maxval(abs(avg-other))==0,'constant-angle API exactly retained')
  enddo

  ! Frozen baseline05ad799 snapshot: SS3 N51+51 dtmax .25, t=97.1927147578s.
  ! Source roe_path_nodes.csv has endpoints alpha, p_SI, U3_nd. EOS/scales from
  ! its run_metadata.txt; these constants preserve the diagnostic, not new data.
  sim%cfg%closure_riser=3;sim%cfg%quad_points=3
  call snapshot_state(.62292403044157918_wp,228011.27491391142_wp,.0015401562717948822_wp,left)
  call snapshot_state(.62668767630261180_wp,227260.90373115471_wp,.0015167684259870570_wp,right)
  theta_left=.98521318670228453_wp;theta_right=1.0092051213023341_wp
  call sim%point_state(left,theta_left,v,fl,src,jac,ds,code)
  call check(code==0,'frozen SS3 left cell is admissible at its angle')
  call sim%point_state(right,theta_right,v,fr,src,jac,ds,code)
  call check(code==0.and.v(5)>0,'frozen SS3 right liquid speed remains positive')
  avg=path_matrix(left,right,.5_wp*(theta_left+theta_right),.5_wp*(theta_left+theta_right),3,code)
  call check(code==14.and..not.all(ieee_is_finite(avg)),'frozen SS3 reproduces old mean-angle primitive failure')
  avg=path_matrix(left,right,theta_left,theta_right,3,code)
  call check(code==0,'frozen SS3 augmented path is admitted without clipping')
  call roe_absA_from_A(avg,diss,speed,code)
  call check(code==0.and.all(ieee_is_finite(diss)),'frozen SS3 augmented path retains real spectrum')
  call sim%face_pair(left,right,theta_left,theta_right,fm,fp,code)
  call check(code==0.and.maxval(abs(fm-fp))==0,'frozen SS3 conservative face succeeds')
  avg=path_matrix(left,right,theta_left,ieee_value(0.0_wp,ieee_quiet_nan),3,code)
  call check(code/=0.and..not.all(ieee_is_finite(avg)),'nonfinite geometry fails explicitly')
  write(*,'(a,i0,a)')'PASS ',count,' Roe geometry checks'
contains
  subroutine check(ok,label)
    logical,intent(in)::ok
    character(*),intent(in)::label
    count=count+1
    if(.not.ok)then
      write(*,'(a)')'FAIL: '//label
      error stop 1
    endif
  end subroutine
  subroutine physical_flux(u,theta,f)
    real(wp),intent(in)::u(3),theta
    real(wp),intent(out)::f(3)
    real(wp)::v(6),s(3),a(3,3),ds(3,3)
    integer::status
    call sim%point_state(u,theta,v,f,s,a,ds,status,.false.)
    call check(status==0,'independent physical flux status')
  end subroutine
  subroutine state_from_fluxes(p,gl,gg,theta,u)
    real(wp),intent(in)::p,gl,gg,theta
    real(wp),intent(out)::u(3)
    real(wp)::rl,rg,a
    integer::status
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    a=solve_alpha_from_superficials(gl/rl,gg/rg,rl,rg,theta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status,sim%cfg%closure_riser)
    call check(status==0,'prescribed phase flux state')
    u=[(1-a)*rl,a*rg,gl+gg]
  end subroutine
  subroutine snapshot_state(a,p_si,momentum,u)
    real(wp),intent(in)::a,p_si,momentum
    real(wp),intent(out)::u(3)
    real(wp)::p,rl,rg
    p=p_si/cfg%omega_P
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    u=[(1-a)*rl,a*rg,momentum]
  end subroutine
  function path_matrix(x,y,thl,thr,n,status,start,finish) result(a)
    real(wp),intent(in)::x(3),y(3),thl,thr
    integer,intent(in)::n
    integer,intent(out)::status
    real(wp),intent(in),optional::start,finish
    real(wp)::a(3,3)
    a=calculate_roe_matrix(x,y,cfg%nCl,cfg%nCg,thl,cfg%nrho_l0,cfg%nP_l0,cfg%D_H,cfg%AREA, &
      cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,0.0_wp,cfg%omega_u,cfg%omega_rho,cfg%omega_P, &
      cfg%tol,cfg%tola,n,a=start,b=finish,status=status,closure_id=sim%cfg%closure_riser,theta_right=thr)
  end function
end program
