program test_fv_operators
  !! Operator-level verification; no long transient and no slug-validation claim.
  use mod_kinds,only:wp,i4
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use module_steady_state,only:solve_alpha_from_superficials
  use module_closure_law,only:CdUd_bendiksen,lei_fechamento_or_ndim_simp
  use module_matrix_roe,only:calculate_roe_matrix
  use module_roe_solver,only:roe_absA_from_A
  use mod_void_helpers,only:ffan
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  real(wp),parameter::pi=acos(-1.0_wp)
  real(wp)::u(3),v(6),f(3),s(3),aj(3,3),sj(3,3),theta,alpha,p,rl,rg,ul,ug,cd,ud
  real(wp)::re,mu,jmix,rhom,expected,ff,pip,grav,ds,error,dt,conv,balance(2),weighted(2),start,finish
  real(wp)::diss(3,3),diss_used(3,3),path_a(3,3),path_abs(3,3),path_speed
  real(wp)::left(3),right(3),fm(3),fp(3),fl(3),fr(3),vl(6),vr(6),jl(3,3),jr(3,3),kl(3,3),kr(3,3)
  real(wp)::state(6,3),old(6,3),saved(6,3),r(3,4),a(3,3,4),b(3,3,4),c(3,3,4),norm,fin(3),fout(3)
  real(wp)::rt(3,4),at(3,3,4),bt(3,3,4),ct(3,3,4),n2,fi2(3),fo2(3),dR(3,4),col(3,4)
  real(wp)::up(3),um(3),vp(6),vm(6),fxp(3),fxm(3),sp(3),sm(3),tmpj(3,3),tmps(3,3),h
  integer::st,st2,count,k,m,ii,l
  count=0
  call configure()
  call sim%init(cfg)
  ! Laminar, zero-speed and turbulent source terms, independently evaluated.
  do k=1,4
    sim%cfg%G=cfg%G
    theta=.6_wp;alpha=.25_wp;p=2.0_wp
    if(k<=2)sim%cfg%G=0.0_wp
    select case(k)
    case(1)
      ul=.002_wp
    case(2)
      ul=0.0_wp
    case(3)
      ul=2.0_wp
    case(4)
      ul=-.6_wp
    end select
    call drift_state(alpha,p,ul,theta,u)
    call sim%point_state(u,theta,v,f,s,aj,sj,st)
    call check(st==0,'point state admitted')
    rhom=u(1)+u(2);jmix=(1-v(1))*v(5)+v(1)*v(6)
    mu=sim%cfg%mu_l*(1-v(1))+sim%cfg%mu_g*v(1)
    re=rhom*cfg%omega_rho*abs(jmix)*cfg%omega_u*cfg%D_H/mu
    grav=sim%cfg%G*cfg%Lr/cfg%omega_u**2
    if(re<2000.0_wp)then
      expected=-rhom*grav*sin(theta)-32*cfg%Lr*mu*jmix/(cfg%omega_rho*cfg%omega_u*cfg%D_H**2)
    else
      ff=ffan(cfg%eps/cfg%D_H,re,st)
      expected=-rhom*grav*sin(theta)-2*cfg%Lr/cfg%D_H*ff*rhom*jmix*abs(jmix)
    endif
    call close(s(3),expected,1e-12_wp,'Fanning source uses 2L/D, or its exact laminar limit')
    call check(maxval(abs(s(:2)))==0.0_wp,'phase sources exactly zero')
    do m=1,3
      h=3e-6_wp*max(abs(u(m)),merge(1e-5_wp,1e-2_wp,m==2))
      up=u;um=u;up(m)=up(m)+h;um(m)=um(m)-h
      call sim%point_state(up,theta,vp,fxp,sp,tmpj,tmps,st,.false.)
      call check(st==0,'source plus state')
      call sim%point_state(um,theta,vm,fxm,sm,tmpj,tmps,st,.false.)
      call check(st==0,'source minus state')
      call close(sj(3,m),(sp(3)-sm(3))/(2*h),5e-5_wp,'source Jacobian versus independent finite difference')
    enddo
  enddo
  sim%cfg%G=cfg%G

  ! Complete face derivatives include the variation of the dissipation matrix.
  call drift_state(.25_wp,2.0_wp,.5_wp,.6_wp,left)
  call drift_state(.255_wp,2.002_wp,.52_wp,.6_wp,right)
  call face_derivative_check(left,right,.6_wp,.6_wp)
  call sim%face_pair(left,right,.6_wp,.6_wp,fl,fr,st,used_dissipation=diss)
  call check(st==0,'capture base dissipation for lagged step')
  path_a=calculate_roe_matrix(left,right,cfg%nCl,cfg%nCg,.6_wp,cfg%nrho_l0,cfg%nP_l0,cfg%D_H,cfg%AREA, &
      cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,0.0_wp,cfg%omega_u,cfg%omega_rho,cfg%omega_P, &
      cfg%tol,cfg%tola,cfg%quad_points,status=st,closure_id=2)
  call check(st==0,'independent configured Roe quadrature')
  call roe_absA_from_A(path_a,path_abs,path_speed,st)
  call check(st==0.and.maxval(abs(path_abs-diss))<1e-12_wp,'face dissipation honors configured Roe quadrature')
  call sim%face_pair(left,right,.6_wp,.6_wp,fm,fp,st,dissipation_matrix=diss,used_dissipation=diss_used)
  call check(st==0.and.maxval(abs(diss-diss_used))==0.0_wp,'dissipation override used exactly')
  call check(maxval(abs(fm-fl))<1e-12_wp.and.maxval(abs(fp-fr))<1e-12_wp,'cached base matrix preserves base flux')
  sim%cfg%jac_dissipation_fd=.false.
  call face_derivative_check(left,right,.6_wp,.6_wp,diss)
  sim%cfg%jac_dissipation_fd=.true.

  ! A prescribed junction equilibrium has equal P, Gl and Gg, with distinct
  ! holdups under stratified and drift closures. Both momentum traces matter.
  call flux_state(2.0_wp,cfg%mul,cfg%mug,-.02_wp,1,left)
  call flux_state(2.0_wp,cfg%mul,cfg%mug,0.0_wp,2,right)
  call sim%point_state(left,-.02_wp,vl,fl,s,aj,sj,st)
  call check(st==0,'matched left pipe trace')
  call sim%point_state(right,0.0_wp,vr,fr,s,aj,sj,st)
  call check(st==0,'matched right riser trace')
  call check(abs(vl(1)-vr(1))>.01_wp,'matched junction retains distinct holdups')
  call cpu_time(start)
  call sim%face_pair(left,right,-.02_wp,0.0_wp,fm,fp,st,jl,jr,kl,kr)
  call cpu_time(finish)
  write(*,'(a,es12.4)')'Matched face with full Jacobian CPU seconds: ',finish-start
  call check(st==0,'matched equilibrium face status')
  call check(maxval(abs(fm(:2)-fp(:2)))==0.0_wp,'one common phase flux at junction')
  call check(maxval(abs(fm-fl)/max(1.0_wp,abs(fl)))<2e-9_wp,'left matched momentum and phase flux')
  call check(maxval(abs(fp-fr)/max(1.0_wp,abs(fr)))<2e-9_wp,'right matched momentum and phase flux')
  pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
  call close((fm(3)-fm(1)*vl(5)-fm(2)*vl(6))/pip,2.0_wp,1e-10_wp,'matched left pressure')
  call close((fp(3)-fp(1)*vr(5)-fp(2)*vr(6))/pip,2.0_wp,1e-10_wp,'matched right pressure')
  call check(abs(fm(3)-fp(3))>1e-4_wp,'junction force is not erased by forcing momentum flux equal')
  call face_derivative_check(left,right,-.02_wp,0.0_wp)

  ! The nonlinear Bhagwat transcription uses the same recovery and operator
  ! contracts, with its co-current domain enforced at every perturbed state.
  sim%cfg%closure_riser=3
  call flux_state(2.0_wp,cfg%mul,cfg%mug,.6_wp,3,left)
  call sim%point_state(left,.6_wp,vl,fl,s,aj,sj,st)
  call check(st==0,'Bhagwat finite-volume point state')
  call close(fl(1),cfg%mul,1e-10_wp,'Bhagwat prescribed liquid mass flux recovered')
  call close(fl(2),cfg%mug,1e-10_wp,'Bhagwat prescribed gas mass flux recovered')
  do m=1,3
    h=3e-6_wp*max(abs(left(m)),merge(1e-5_wp,1e-2_wp,m==2))
    up=left;um=left;up(m)=up(m)+h;um(m)=um(m)-h
    call sim%point_state(up,.6_wp,vp,fxp,sp,tmpj,tmps,st,.false.)
    call check(st==0,'Bhagwat source plus state')
    call sim%point_state(um,.6_wp,vm,fxm,sm,tmpj,tmps,st,.false.)
    call check(st==0,'Bhagwat source minus state')
    call close(sj(3,m),(sp(3)-sm(3))/(2*h),5e-5_wp,'Bhagwat source Jacobian versus conservative differences')
  enddo
  call flux_state(2.002_wp,1.002_wp*cfg%mul,.997_wp*cfg%mug,.6_wp,3,right)
  call face_derivative_check(left,right,.6_wp,.6_wp)
  call flux_state(2.0_wp,cfg%mul,cfg%mug,-.02_wp,1,left)
  call flux_state(2.0_wp,cfg%mul,cfg%mug,0.0_wp,3,right)
  call sim%point_state(left,-.02_wp,vl,fl,s,aj,sj,st)
  call check(st==0,'Bhagwat matched left pipe trace')
  call sim%point_state(right,0.0_wp,vr,fr,s,aj,sj,st)
  call check(st==0,'Bhagwat matched right riser trace')
  call sim%face_pair(left,right,-.02_wp,0.0_wp,fm,fp,st)
  call check(st==0,'Bhagwat matched equilibrium face')
  call check(maxval(abs(fm(:2)-fp(:2)))==0.0_wp,'Bhagwat matched common phase fluxes')
  call check(maxval(abs(fm-fl)/max(1.0_wp,abs(fl)))<2e-9_wp,'Bhagwat matched left momentum trace')
  call check(maxval(abs(fp-fr)/max(1.0_wp,abs(fr)))<2e-9_wp,'Bhagwat matched right momentum trace')
  call close((fm(3)-fm(1)*vl(5)-fm(2)*vl(6))/pip,2.0_wp,1e-10_wp,'Bhagwat matched left pressure')
  call close((fp(3)-fp(1)*vr(5)-fp(2)*vr(6))/pip,2.0_wp,1e-10_wp,'Bhagwat matched right pressure')
  call face_derivative_check(left,right,-.02_wp,0.0_wp)
  sim%cfg%closure_riser=2

  ! Four physical cells, one ghost each side; nonuniform widths. The weighted
  ! residual must equal accumulation plus boundary mass transport, in kg.
  do ii=2,5
    l=ii-1;p=2.02_wp-.006_wp*real(l-1,wp)
    if(ii<=3)then
      call flux_state(p,cfg%mul,cfg%mug,cfg%thetav(ii),1,state(ii,:))
    else
      call flux_state(p,cfg%mul,cfg%mug,cfg%thetav(ii),2,state(ii,:))
    endif
  enddo
  state(1,:)=state(2,:);state(6,:)=state(5,:)
  old=state
  old(2:5,1)=old(2:5,1)*.999_wp
  old(2:5,2)=old(2:5,2)*1.001_wp
  saved=state
  dt=1e-5_wp
  call sim%assemble(state,old,dt,r,a,b,c,norm,fin,fout,st,.true.)
  call check(st==0,'full FV assembly status')
  call check(maxval(abs(state-saved))==0.0_wp,'assembly does not mutate physical or ghost states')
  conv=cfg%omega_rho*cfg%AREA*cfg%Lr
  do k=1,2
    weighted(k)=sum(r(k,:)*cfg%cell_width(2:5))*conv
    balance(k)=sum((state(2:5,k)-old(2:5,k))*cfg%cell_width(2:5))*conv+ &
      dt*cfg%omega_t*cfg%omega_rho*cfg%omega_u*cfg%AREA*(fout(k)-fin(k))
    call close(weighted(k),balance(k),1e-12_wp,'weighted phase residual telescopes in physical kg')
  enddo
  do ii=2,5
    do m=1,3
      h=3e-6_wp*max(abs(state(ii,m)),merge(1e-5_wp,1e-2_wp,m==2))
      state=saved;state(ii,m)=state(ii,m)+h
      call sim%assemble(state,old,dt,rt,at,bt,ct,n2,fi2,fo2,st,.false.)
      call check(st==0,'assembled Jacobian plus state')
      dR=rt
      state=saved;state(ii,m)=state(ii,m)-h
      call sim%assemble(state,old,dt,rt,at,bt,ct,n2,fi2,fo2,st,.false.)
      call check(st==0,'assembled Jacobian minus state')
      dR=(dR-rt)/(2*h)
      l=ii-1;col=0
      col(:,l)=b(:,m,l)
      if(l>1)col(:,l-1)=c(:,m,l-1)
      if(l<4)col(:,l+1)=a(:,m,l+1)
      error=maxval(abs(col-dR)/max(1.0_wp,abs(dR)))
      if(error>3e-4_wp)write(*,*)'FV Jacobian error at cell/component:',ii,m,error
      call check(error<3e-4_wp,'assembled block Jacobian equals full residual derivative')
    enddo
  enddo
  state=saved
  call sim%assemble(state,old,0.0_wp,r,a,b,c,norm,fin,fout,st,.true.)
  call check(st==0,'zero-time assembly status')
  do ii=2,5
    l=ii-1
    call check(maxval(abs(r(:,l)-(state(ii,:)-old(ii,:))))<1e-15_wp,'every physical cell retains accumulation equation')
  enddo
  call check(maxval(abs(a))==0.0_wp.and.maxval(abs(c))==0.0_wp,'zero-time neighbor blocks vanish')
  do l=1,4
    do m=1,3
      call close(b(m,m,l),1.0_wp,1e-15_wp,'all physical diagonals contain identity')
    enddo
  enddo

  ! Exact stationary limit: uniform pressure, G=0 and zero phase velocities.
  sim%cfg%G=0.0_wp;sim%cfg%mul=0;sim%cfg%mug=0;sim%cfg%thetav=.3_wp
  call drift_state(.25_wp,sim%cfg%nPs,0.0_wp,.3_wp,u)
  do ii=1,6
    state(ii,:)=u
  enddo
  old=state
  call sim%assemble(state,old,1e-3_wp,r,a,b,c,norm,fin,fout,st,.true.)
  call check(st==0.and.maxval(abs(r))<1e-12_wp,'exact uniform static FV equilibrium')
  call ss1_junction_regression()
  write(*,'(a,i0,a)')'PASS finite-volume operators: ',count,' assertions'
contains
  subroutine configure()
    cfg%omega_u=.5_wp;cfg%omega_rho=1000.0_wp;cfg%omega_P=1e5_wp
    cfg%Lr=7.7_wp;cfg%Lp=7.7_wp;cfg%omega_t=cfg%Lr/cfg%omega_u
    cfg%D_H=.0508_wp;cfg%AREA=pi*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
    cfg%mu_l=.89e-3_wp;cfg%mu_g=1.85e-5_wp;cfg%sigma=.0728_wp
    cfg%nCl=150.0_wp;cfg%nCg=30.0_wp;cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2
    cfg%tol=1e-12_wp;cfg%tola=1e-12_wp;cfg%maxit=100
    cfg%closure_pipe=1;cfg%closure_riser=2;cfg%jac_dissipation_fd=.true.;cfg%jac_dissipation_fd_eps=1e-6_wp
    cfg%junction_model='MATCHED';cfg%nCells=6;cfg%num_dummy_cells=1;cfg%N_pipe=2
    allocate(cfg%cell_width(6),cfg%thetav(6))
    cfg%cell_width=[.3_wp,.3_wp,.5_wp,.4_wp,.8_wp,.8_wp]
    cfg%thetav=[-.02_wp,-.02_wp,-.02_wp,0.0_wp,.3_wp,.3_wp]
    cfg%mul=.6_wp*(cfg%nrho_l0+(2-cfg%nP_l0)/cfg%nCl**2)
    cfg%mug=2.0_wp/cfg%nCg**2
  end subroutine
  subroutine drift_state(a,p,u_liquid,theta,U)
    real(wp),intent(in)::a,p,u_liquid,theta
    real(wp),intent(out)::U(3)
    real(wp)::rl,rg,c,d,ug
    integer::code
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    call CdUd_bendiksen(a,rl,rg,u_liquid,u_liquid,theta,cfg%D_H,cfg%AREA,cfg%eps,sim%cfg%G,cfg%mu_l, &
      cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,c,d,code)
    call check(code==0,'drift test-state coefficient status')
    ug=(c*(1-a)*u_liquid+d/cfg%omega_u)/(1-c*a)
    U=[(1-a)*rl,a*rg,(1-a)*rl*u_liquid+a*rg*ug]
  end subroutine
  subroutine flux_state(p,gl,gg,theta,model,U)
    real(wp),intent(in)::p,gl,gg,theta
    integer,intent(in)::model
    real(wp),intent(out)::U(3)
    real(wp)::rl,rg,al
    integer::code
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    al=solve_alpha_from_superficials(gl/rl,gg/rg,rl,rg,theta,cfg%D_H,cfg%AREA,cfg%eps,sim%cfg%G, &
      cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,code,model)
    if(code/=0)write(*,*)'flux state failure:',code,p,gl,gg,theta
    call check(code==0,'phase-flux test state status')
    U=[(1-al)*rl,al*rg,gl+gg]
  end subroutine
  subroutine face_derivative_check(x,y,thl,thr,fixedD)
    real(wp),intent(in)::x(3),y(3),thl,thr
    real(wp),intent(in),optional::fixedD(3,3)
    real(wp)::fm(3),fp(3),jl(3,3),jr(3,3),kl(3,3),kr(3,3),xp(3),xm(3),yp(3),ym(3)
    real(wp)::fm1(3),fm2(3),fp1(3),fp2(3),h,actual(6),expected(6),error
    integer::j,code
    call sim%face_pair(x,y,thl,thr,fm,fp,code,jl,jr,kl,kr,dissipation_matrix=fixedD)
    if(code/=0)write(*,*)'face status:',code,thl,thr
    call check(code==0,'face Jacobian status')
    do j=1,3
      h=3e-6_wp*max(abs(x(j)),merge(1e-5_wp,1e-2_wp,j==2))
      xp=x;xm=x;xp(j)=xp(j)+h;xm(j)=xm(j)-h
      call sim%face_pair(xp,y,thl,thr,fm1,fp1,code,dissipation_matrix=fixedD)
      call check(code==0,'left face plus state')
      call sim%face_pair(xm,y,thl,thr,fm2,fp2,code,dissipation_matrix=fixedD)
      call check(code==0,'left face minus state')
      actual=[jl(:,j),kl(:,j)];expected=[fm1-fm2,fp1-fp2]/(2*h)
      error=maxval(abs(actual-expected)/max(1.0_wp,abs(expected)))
      if(error>2e-4_wp)write(*,*)'left face derivative error:',j,error
      call check(error<2e-4_wp,'complete left face derivatives')
      h=3e-6_wp*max(abs(y(j)),merge(1e-5_wp,1e-2_wp,j==2))
      yp=y;ym=y;yp(j)=yp(j)+h;ym(j)=ym(j)-h
      call sim%face_pair(x,yp,thl,thr,fm1,fp1,code,dissipation_matrix=fixedD)
      call check(code==0,'right face plus state')
      call sim%face_pair(x,ym,thl,thr,fm2,fp2,code,dissipation_matrix=fixedD)
      call check(code==0,'right face minus state')
      actual=[jr(:,j),kr(:,j)];expected=[fm1-fm2,fp1-fp2]/(2*h)
      error=maxval(abs(actual-expected)/max(1.0_wp,abs(expected)))
      if(error>2e-4_wp)write(*,*)'right face derivative error:',j,error
      call check(error<2e-4_wp,'complete right face derivatives')
    enddo
  end subroutine
  subroutine ss1_junction_regression()
    !! Accepted interface states immediately before the former 2071 failure,
    !! pilot_D_SS1 at 39.352838046571165 s. Constants make the regression
    !! independent of optional/large transient output files.
    real(wp)::x(3),y(3),fm(3),fp(3),aa(3,3),bb(3,3),cc(3,3),dd(3,3)
    real(wp)::yp(3),fm1(3),fp1(3),fm2(3),fp2(3),h,fdm(3),fdp(3),trace(6)
    real(wp)::rl,rg,ul,closure_res,pip,xnew(3),ynew(3),dt_local,phase_balance(2),savedx(3),savedy(3),bound_lo,bound_hi
    integer::code,n_one_sided,k,it
    real(wp)::cachedD(3,3)
    logical::blocked
    cfg%Lr=12.546253993197093_wp;cfg%Lp=57.4_wp
    cfg%omega_P=101325.0_wp;cfg%omega_rho=997.05_wp;cfg%omega_u=.53_wp
    cfg%omega_t=cfg%Lr/cfg%omega_u
    cfg%nCl=1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
    cfg%nCg=343.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
    cfg%nrho_l0=1.0_wp;cfg%nP_l0=1.0_wp
    cfg%D_H=.05254_wp;cfg%AREA=pi*cfg%D_H**2/4
    cfg%G=9.81_wp;cfg%tol=1e-11_wp;cfg%tola=1e-10_wp
    cfg%jac_dissipation_fd=.false.;cfg%junction_model='MATCHED'
    call sim%init(cfg)
    x=[4.40242700906511131e-01_wp,1.42728296420726068e-03_wp,1.00050605101821510e+00_wp]
    y=[9.89113033765173899e-01_wp,2.79928483699783725e-05_wp,9.84310916863125107e-01_wp]
    call sim%face_pair(x,y,-.03490658503988659_wp,.08638636790679209_wp,fm,fp,code)
    call check(code==0,'frozen SS1 snapshot remains admissible with configured Roe quadrature')
    ! Locate the nearby model-domain boundary for the selected quadrature.
    ! The former endpoint-average matrix placed this boundary at a slightly
    ! different cell mass. This constructs a numerical derivative test; it
    ! does not change any production state or tune a physical closure.
    bound_lo=y(1);bound_hi=y(1)*(1+1e-4_wp);yp=y;yp(1)=bound_hi
    call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm,fp,code)
    if(code/=1021)write(*,*)'SS1 boundary bracket status:',code,bound_lo,bound_hi
    call check(code==1021,'nearby SS1 upper bracket fails only the riser star domain')
    do it=1,30
      yp=y;yp(1)=.5_wp*(bound_lo+bound_hi)
      call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm,fp,code)
      if(code==0)then
        bound_lo=yp(1)
      else if(code==1021)then
        bound_hi=yp(1)
      else
        call check(.false.,'unexpected failure while locating SS1 model-domain boundary')
      endif
    enddo
    y(1)=bound_lo-.25_wp*cfg%jac_dissipation_fd_eps*bound_lo
    call sim%face_pair(x,y,-.03490658503988659_wp,.08638636790679209_wp,fm,fp,code,aa,bb,cc,dd, &
                       fd_one_sided=n_one_sided)
    if(code/=0)write(*,*)'SS1 recovered branch regression status:',code
    call check(code==0,'SS1 accepted state and its face derivatives remain admissible')
    call check(all(ieee_is_finite([fm,fp])),'SS1 junction flux is finite')
    call check(maxval(abs(fm(:2)-fp(:2)))==0.0_wp,'SS1 branch continuation preserves common phase flux')
    call check(n_one_sided>0,'SS1 uses reported domain-respecting one-sided columns')
    call check(all(ieee_is_finite(aa)).and.all(ieee_is_finite(bb)).and. &
               all(ieee_is_finite(cc)).and.all(ieee_is_finite(dd)),'SS1 complete face Jacobian is finite')
    h=3e-6_wp*y(1);yp=y;yp(1)=y(1)+h
    call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm1,fp1,code)
    call check(code/=0,'SS1 inadmissible two-phase star remains an explicit failure')
    yp=y;yp(1)=y(1)-h
    call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm1,fp1,code)
    call check(code==0,'SS1 first backward sample admissible')
    yp=y;yp(1)=y(1)-2*h
    call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm2,fp2,code)
    call check(code==0,'SS1 second backward sample admissible')
    fdm=(3*fm-4*fm1+fm2)/(2*h);fdp=(3*fp-4*fp1+fp2)/(2*h)
    call check(maxval(abs(bb(:,1)-fdm)/max(1.0_wp,abs(fdm)))<3e-4_wp,'SS1 one-sided left trace derivative convergence')
    call check(maxval(abs(dd(:,1)-fdp)/max(1.0_wp,abs(fdp)))<3e-4_wp,'SS1 one-sided right trace derivative convergence')

    ! BLOCKING is an explicitly selected experimental unilateral interface
    ! constraint, not a claim that the reference moving-front model was solved.
    ! Inactive states must reproduce MATCHED exactly.
    sim%cfg%junction_model='BLOCKING'
    call sim%face_pair(x,y,-.03490658503988659_wp,.08638636790679209_wp,fm1,fp1,code, &
                      blocking_active=blocked,junction_trace=trace)
    call check(code==0.and..not.blocked,'inactive blocking operator is reported')
    call check(maxval(abs(fm1-fm))<1e-12_wp.and.maxval(abs(fp1-fp))<1e-12_wp,'inactive BLOCKING equals MATCHED')
    yp=y;yp(1)=y(1)*(1+1e-4_wp);savedx=x;savedy=yp
    call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm1,fp1,code,aa,bb,cc,dd, &
                      blocking_active=blocked,junction_trace=trace,used_dissipation=cachedD)
    call check(code==0.and.blocked,'inadmissible reverse gas proposal activates explicit barrier')
    call check(trace(6)<0.0_wp.and.trace(3)==0.0_wp,'raw gas proposal retained in diagnostic trace')
    call check(fm1(2)==0.0_wp.and.fp1(2)==0.0_wp,'blocked interface uses one zero gas flux')
    call check(fm1(1)==fp1(1).and.fm1(1)>0.0_wp,'blocked liquid transport remains common and positive')
    call check(trace(4)>0.0_wp.and.trace(4)<1.0_wp.and.trace(5)==0.0_wp,'blocked pipe keeps stagnant gas and riser trace is liquid')
    call check(maxval(abs(x-savedx))==0.0_wp.and.maxval(abs(yp-savedy))==0.0_wp,'blocking changes no conservative cell inventory')
    rl=cfg%nrho_l0+(trace(1)-cfg%nP_l0)/cfg%nCl**2;rg=trace(1)/cfg%nCg**2
    ul=trace(2)/((1-trace(4))*rl)
    call lei_fechamento_or_ndim_simp(trace(4),rl,rg,ul,0.0_wp,-.03490658503988659_wp,cfg%D_H,cfg%AREA, &
      cfg%eps,cfg%G,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,closure_res,1,code)
    call check(code==0.and.abs(closure_res)<1e-9_wp,'blocked pipe trace satisfies stratified closure with stagnant gas')
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    call close((fm1(3)-fm1(1)*ul)/pip,trace(1),1e-12_wp,'blocked pipe pressure trace')
    call close((fp1(3)-fp1(1)**2/rl)/pip,trace(1),1e-12_wp,'blocked pure-liquid riser pressure trace')
    do k=1,3
      dt_local=real(k,wp)*1e-7_wp
      xnew=x-dt_local*fm1/.2_wp;ynew=yp+dt_local*fp1/.04_wp
      phase_balance=((xnew(:2)-x(:2))*.2_wp+(ynew(:2)-yp(:2))*.04_wp)* &
                     cfg%omega_rho*cfg%AREA*cfg%Lr
      call check(maxval(abs(phase_balance))<1e-13_wp,'blocked interface mass telescopes at different time increments')
    enddo
    call face_derivative_check(x,yp,-.03490658503988659_wp,.08638636790679209_wp,cachedD)
    call face_derivative_check(x,yp,-.03490658503988659_wp,.08638636790679209_wp)
    sim%cfg%junction_model='MATCHED'
    call sim%face_pair(x,yp,-.03490658503988659_wp,.08638636790679209_wp,fm2,fp2,code)
    call check(code/=0,'strict MATCHED remains a distinct model and rejects inadmissible trace')

  end subroutine
  subroutine check(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    if(.not.condition)then
      write(*,'(a)')'FAIL: '//label
      error stop 1
    endif
    count=count+1
  end subroutine
  subroutine close(x,y,eps,label)
    real(wp),intent(in)::x,y,eps
    character(*),intent(in)::label
    if(.not.all(ieee_is_finite([x,y])))call check(.false.,label)
    if(abs(x-y)>eps*max(1.0_wp,abs(x),abs(y)))write(*,*)label,x,y,abs(x-y)
    call check(abs(x-y)<=eps*max(1.0_wp,abs(x),abs(y)),label)
  end subroutine
end program test_fv_operators
