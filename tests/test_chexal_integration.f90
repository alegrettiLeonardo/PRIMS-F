program test_chexal_integration
  !! Numerical integration tests. These are not experimental validation data.
  use mod_kinds,only:wp,i4
  use module_closure_law,only:CHEXAL_1992,CHEXAL_1992_HIGH,CHEXAL_1992_LOW, &
    lei_fechamento_or_ndim_simp,Dlei_fechamento_or_ndim_simp
  use module_find,only:compute_primitive_variables
  use module_steady_state,only:solve_alpha_from_superficials,PHASE_LIQUID_ONLY,PHASE_GAS_ONLY
  use module_jacobian,only:DFduj
  use module_roe_solver,only:roe_eigenstructure
  use module_boundary_conditions,only:boundary_flux
  use module_config,only:sim_config_t
  use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  real(wp),parameter::pi=acos(-1.0_wp),diam=.05254_wp,area=pi*diam**2/4
  real(wp),parameter::rough=4.6e-5_wp,g=9.81_wp,mul=.89e-3_wp,mug=1.85e-5_wp,sigma=.0728_wp,tol=1e-12_wp
  real(wp)::wu,wr,wpres,cl,cg,rl0,p0,p,rl,rg,alpha,ul,ug,theta,jl,jg,rr
  real(wp)::state(3),up(3),um(3),q(6),qp(6),qm(6),fp(3),fm(3),jac(3,3),jfd(3,3)
  real(wp)::eig(3),imag(3),right(3,3),left(3,3),h,bc(3),f(3),base_alpha,base_ul,base_ug
  real(wp)::v(6),src(3),sj(3,3),face_left(3),face_right(3),der(5)
  integer::count,k,m,model,st,nr,idx
  logical::ok
  type(sim_config_t)::cfg
  type(pipeline_riser_simulator_t)::sim
  count=0
  do k=1,3
    call scales(k)
    do m=1,3
      theta=real(m-1,wp)*pi/4
      model=CHEXAL_1992
      jl=.4_wp/wu;jg=.2_wp/wu
      call build_state(jl,jg,theta,model,state)
      call recover(state,theta,model,q,st)
      call check(st==0,'Chexal U to primitive succeeds')
      call near(q(1),alpha,2e-9_wp,'alpha round trip')
      call near(q(5)*wu,ul*wu,2e-8_wp,'liquid velocity round trip SI')
      call near(q(6)*wu,ug*wu,2e-8_wp,'gas velocity round trip SI')
      call near(state(1)*q(5)+state(2)*q(6),state(3),1e-11_wp,'momentum identity')
      if(m==2)then
        if(k==1)then
          base_alpha=alpha;base_ul=ul*wu;base_ug=ug*wu
        else
          call near(alpha,base_alpha,2e-10_wp,'alpha independent of nondimensional scales')
          call near(ul*wu,base_ul,2e-10_wp,'liquid SI scale invariance')
          call near(ug*wu,base_ug,2e-10_wp,'gas SI scale invariance')
        endif
      endif
      call lei_fechamento_or_ndim_simp(q(1),q(2),q(3),q(5),q(6),theta,diam,area,rough,g,mul,mug,sigma, &
        wu,wr,tol,rr,model,st)
      call check(st==0,'Chexal wrapper residual status')
      call near(rr*wu,0.0_wp,2e-10_wp,'R = alpha*(C0*j+Vgj)-jg')
      jac=DFduj(q(1),q(2),q(3),q(5),q(6),cl,cg,theta,diam,area,rough,g,mul,mug,sigma, &
        wu,wr,wpres,tol,st,model)
      call check(st==0,'same branch flux Jacobian')
      do idx=1,3
        h=2e-6_wp*abs(state(idx))
        up=state;um=state;up(idx)=up(idx)+h;um(idx)=um(idx)-h
        call recover(up,theta,model,qp,st)
        call check(st==0,'Jacobian plus state remains admissible')
        call recover(um,theta,model,qm,st)
        call check(st==0,'Jacobian minus state remains admissible')
        fp=flux(up,qp);fm=flux(um,qm)
        jfd(:,idx)=(fp-fm)/(2*h)
      enddo
      call check(maxval(abs(jac-jfd)/(1.0_wp+abs(jfd)))<8e-5_wp,'Jacobian versus independent U flux differences')
      call roe_eigenstructure(jac,eig,imag,right,left,ok,st)
      call check(ok.and.st==0.and.all(imag==0.0_wp),'co-current fixture hyperbolicity without spectral repair')
    enddo
  enddo
  call scales(1)
  theta=.6_wp;model=CHEXAL_1992
  call build_state(.4_wp/wu,.2_wp/wu,theta,model,state)
  call configure()
  call boundary_flux(state,theta,cfg,.true.,f,st,bc)
  call check(st==0,'Chexal inlet')
  call near(f(1),cfg%mul,1e-13_wp,'inlet liquid mass conservation')
  call near(f(2),cfg%mug,1e-13_wp,'inlet gas mass conservation')
  call boundary_flux(state,theta,cfg,.false.,f,st,bc)
  call check(st==0,'Chexal outlet')
  call recover(bc,theta,model,q,st)
  call check(st==0,'outlet conservative state recovery')
  call near(q(4),cfg%nPs,2e-10_wp,'outlet pressure')
  call sim%init(cfg)
  call sim%point_state(state,theta,v,f,src,jac,sj,st)
  call check(st==0,'Chexal source/flux/Jacobian operator')
  call check(all(src(1:2)==0.0_wp),'no phase mass sources')
  call sim%face_pair(state,state,theta,theta,face_left,face_right,st)
  call check(st==0,'Chexal Roe face')
  call check(maxval(abs(face_left-face_right))==0.0_wp,'single conservative face flux')
  call check(maxval(abs(face_left-f)/(1.0_wp+abs(f)))<2e-10_wp,'consistent equal-state physical flux')

  ! Exact inventory limits are separate from two-phase spectral derivatives.
  state=[rl,0.0_wp,rl*.3_wp/wu]
  call recover(state,theta,model,q,st)
  call check(st==0.and.q(1)==0.0_wp.and.q(6)==0.0_wp,'absent gas exact inventory limit')
  call near(q(5)*wu,.3_wp,1e-12_wp,'pure liquid speed')
  state=[0.0_wp,rg,rg*.8_wp/wu]
  call recover(state,theta,model,q,st)
  call check(st==0.and.q(1)==1.0_wp.and.q(5)==0.0_wp,'absent liquid exact inventory limit')
  call near(q(6)*wu,.8_wp,1e-12_wp,'pure gas speed')
  alpha=solve_alpha_from_superficials(.3_wp,0.0_wp,rl,rg,theta,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,tol,st,model,phase_hint=PHASE_LIQUID_ONLY)
  call check(st==0.and.alpha==0.0_wp,'explicit liquid phase hint')
  alpha=solve_alpha_from_superficials(0.0_wp,.8_wp,rl,rg,theta,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,tol,st,model,phase_hint=PHASE_GAS_ONLY)
  call check(st==0.and.alpha==1.0_wp,'explicit gas phase hint')
  jac=DFduj(1.0_wp,rl,rg,0.0_wp,.8_wp/wu,cl,cg,theta,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,wpres,tol,st,model)
  call check(st/=0,'two-phase spectrum at absent phase is explicitly undefined')
  call lei_fechamento_or_ndim_simp(.4_wp,rl,rg,-.1_wp,.2_wp,.6_wp,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,tol,rr,CHEXAL_1992_HIGH,st)
  call check(st/=0,'inclined countercurrent remains excluded')
  call lei_fechamento_or_ndim_simp(.4_wp,rl,rg,.1_wp,-.2_wp,0.0_wp,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,tol,rr,CHEXAL_1992,st)
  call check(st/=0,'negative horizontal MATCHED gas trace remains excluded')
  alpha=solve_alpha_from_superficials(.3_wp,0.0_wp,rl,rg,theta,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,tol,st,model)
  call check(st/=0.and..not.ieee_is_finite(alpha),'zero flow without inventory never silently creates pure phase')
  alpha=solve_alpha_from_superficials(.3_wp,.2_wp,rl,rg,theta,diam,2*area,rough,g,mul,mug,sigma, &
    wu,wr,tol,st,model)
  call check(st/=0,'steady rejects inconsistent area/diameter')
  call compute_primitive_variables(state(1),state(2),state(3),theta,cg,cl,rl0,p0,diam,2*area,rough,g,mul,mug,sigma, &
    0.0_wp,wu,wr,tol,tol,100_i4,q(1),q(2),q(3),q(4),q(5),q(6),idx,status=st,closure_id=model)
  call check(st/=0,'primitive rejects inconsistent area/diameter even at pure-phase limit')
  ! C3=max(.5,2 exp(-Rel/60000)) has an actual constitutive kink.
  ul=60000.0_wp*log(4.0_wp)*mul/(rl*wr*diam)/(1.0_wp-.2_wp)/wu
  call Dlei_fechamento_or_ndim_simp(.2_wp,rl,rg,ul,.1_wp/(.2_wp*wu),pi/2,diam,area,rough,g,mul,mug,sigma, &
    wu,wr,tol,der(1),der(3),der(2),der(5),der(4),CHEXAL_1992,st)
  call check(st/=0,'constitutive kink cannot produce a smooth artificial Jacobian')
  call test_countercurrent()
  write(*,'(a,i0,a)')'PASS Chexal integration: ',count,' assertions'
contains
  subroutine test_countercurrent()
    ! Independent 60-digit fixtures: Chexal1992, air/water, at fixed jl/jg.
    real(wp),parameter::cc_jl(2)=[-.221900040854923588_wp,-.1_wp]
    real(wp),parameter::expected_low(2)=[.220214682534698979_wp,.213582292001802460_wp]
    real(wp),parameter::expected_high(2)=[.839848859537789440_wp,.960796750389095244_wp]
    real(wp)::ac,uc(3),qq(6),qqp(6),qqm(6),plus(3),minus(3),diff(3,3),matrix(3,3)
    real(wp)::liquid,gas,hh,expect,rf,rt(3,3),lt(3,3),ev(3),iv(3)
    integer::scale_id,branch_id,case_id,jj,regime,cc_status,roots_count,orientation
    logical::hyperbolic
    do scale_id=1,3
      call scales(scale_id)
      p=1.2_wp*343.0_wp**2/wpres;p0=p;rl0=998.0_wp/wr
      rl=rl0;rg=1.2_wp/wr
      do orientation=-1,1,2
      theta=orientation*pi/2
      do branch_id=CHEXAL_1992_HIGH,CHEXAL_1992_LOW
        do case_id=1,2
          ac=solve_alpha_from_superficials(orientation*cc_jl(case_id)/wu,orientation*.1_wp/wu,rl,rg,theta,diam,area,rough,g, &
            .001_wp,.000018_wp,.072_wp,wu,wr,tol,cc_status,branch_id,roots_count)
          call check(cc_status==0,'explicit vertical countercurrent alpha branch')
          expect=expected_high(case_id)
          if(branch_id==CHEXAL_1992_LOW)expect=expected_low(case_id)
          call near(ac,expect,3e-9_wp,'independent countercurrent alpha reference')
          liquid=orientation*cc_jl(case_id)/(1-ac)/wu;gas=orientation*.1_wp/ac/wu
          uc=[(1-ac)*rl,ac*rg,orientation*(rl*cc_jl(case_id)/wu+rg*.1_wp/wu)]
          call cc_recover(uc,qq,cc_status,branch_id,liquid,gas)
          call check(cc_status==0,'countercurrent primitive continuation')
          call near(qq(5)*wu,liquid*wu,2e-8_wp,'CC liquid round trip')
          call near(qq(6)*wu,gas*wu,2e-8_wp,'CC gas round trip')
          if(scale_id==1)then
            call compute_primitive_variables(uc(1),uc(2),uc(3),theta,cg,cl,rl0,p0,diam,area,rough,g, &
              .001_wp,.000018_wp,.072_wp,0.0_wp,wu,wr,tol,tol,100_i4, &
              qq(1),qq(2),qq(3),qq(4),qq(5),qq(6),regime,status=cc_status,closure_id=branch_id)
            call check(cc_status==0,'unique declared CC branch recoverable without fabricated hint')
            call near(qq(6)*wu,gas*wu,2e-8_wp,'unhinted recovery stays on declared CC root')
          endif
          call lei_fechamento_or_ndim_simp(ac,rl,rg,liquid,gas,theta,diam,area,rough,g, &
            .001_wp,.000018_wp,.072_wp,wu,wr,tol,rf,branch_id,cc_status)
          call check(cc_status==0,'CC residual including total volume flux zero')
          call near(rf*wu,0.0_wp,3e-10_wp,'CC residual reference')
          matrix=DFduj(ac,rl,rg,liquid,gas,cl,cg,theta,diam,area,rough,g, &
            .001_wp,.000018_wp,.072_wp,wu,wr,wpres,tol,cc_status,branch_id)
          call check(cc_status==0,'countercurrent same branch Jacobian')
          do jj=1,3
            hh=2e-6_wp*abs(uc(jj));plus=uc;minus=uc
            plus(jj)=plus(jj)+hh;minus(jj)=minus(jj)-hh
            call cc_recover(plus,qqp,cc_status,branch_id,liquid,gas)
            call check(cc_status==0,'countercurrent Jacobian plus continuation')
            call cc_recover(minus,qqm,cc_status,branch_id,liquid,gas)
            call check(cc_status==0,'countercurrent Jacobian minus continuation')
            diff(:,jj)=(flux(plus,qqp)-flux(minus,qqm))/(2*hh)
          enddo
          call check(maxval(abs(matrix-diff)/(1+abs(diff)))<2e-4_wp,'CC Jacobian versus same branch U perturbations')
          call roe_eigenstructure(matrix,ev,iv,rt,lt,hyperbolic,cc_status)
          if(any(iv/=0.0_wp))then
            call check(.not.hyperbolic.and.cc_status/=0,'complex CC spectrum rejected without repair')
          else
            call check(hyperbolic.and.cc_status==0,'real CC spectrum verified')
          endif
          write(*,'(a,3i3,a,3es18.9)')'CC spectrum scale/branch/case ',scale_id,branch_id,case_id,' imag=',iv
        enddo
      enddo
      enddo
    enddo
  end subroutine
    subroutine cc_recover(u,qq,status,branch_id,liquid,gas)
      real(wp),intent(in)::u(3)
      real(wp),intent(out)::qq(6)
      integer,intent(out)::status
      integer,intent(in)::branch_id
      real(wp),intent(in)::liquid,gas
      integer(i4)::regime
      call compute_primitive_variables(u(1),u(2),u(3),theta,cg,cl,rl0,p0,diam,area,rough,g, &
        .001_wp,.000018_wp,.072_wp,0.0_wp,wu,wr,tol,tol,100_i4, &
        qq(1),qq(2),qq(3),qq(4),qq(5),qq(6),regime,status=status,closure_id=branch_id, &
        ul_prev=liquid,ug_prev=gas)
    end subroutine
  subroutine scales(k)
    integer,intent(in)::k
    select case(k)
    case(1)
      wu=.7_wp;wr=1000.0_wp;wpres=1e5_wp
    case(2)
      wu=1.0_wp;wr=1.0_wp;wpres=1.0_wp
    case(3)
      wu=3.0_wp;wr=37.0_wp;wpres=7e4_wp
    end select
    cl=1498.0_wp*sqrt(wr/wpres);cg=343.0_wp*sqrt(wr/wpres)
    rl0=997.05_wp/wr;p0=101300.0_wp/wpres;p=2e5_wp/wpres
    rl=rl0+(p-p0)/cl**2;rg=p/cg**2
  end subroutine
  subroutine build_state(jl,jg,theta,model,u)
    real(wp),intent(in)::jl,jg,theta
    integer,intent(in)::model
    real(wp),intent(out)::u(3)
    alpha=solve_alpha_from_superficials(jl,jg,rl,rg,theta,diam,area,rough,g,mul,mug,sigma, &
      wu,wr,tol,st,model,nr)
    call check(st==0.and.nr==1,'unique co-current root')
    ul=jl/(1.0_wp-alpha);ug=jg/alpha
    u=[(1.0_wp-alpha)*rl,alpha*rg,rl*jl+rg*jg]
  end subroutine
  subroutine recover(u,theta,model,q,status)
    real(wp),intent(in)::u(3),theta
    integer,intent(in)::model
    real(wp),intent(out)::q(6)
    integer,intent(out)::status
    integer(i4)::regime
    call compute_primitive_variables(u(1),u(2),u(3),theta,cg,cl,rl0,p0,diam,area,rough,g,mul,mug,sigma, &
      0.0_wp,wu,wr,tol,tol,100_i4,q(1),q(2),q(3),q(4),q(5),q(6),regime,status=status,closure_id=model)
  end subroutine
  function flux(u,q) result(f)
    real(wp),intent(in)::u(3),q(6)
    real(wp)::f(3)
    f=[u(1)*q(5),u(2)*q(6),u(1)*q(5)**2+u(2)*q(6)**2+wpres/(wr*wu**2)*q(4)]
  end function
  subroutine configure()
    cfg%closure_pipe=model;cfg%closure_riser=model
    cfg%omega_u=wu;cfg%omega_rho=wr;cfg%omega_P=wpres
    cfg%nCl=cl;cfg%nCg=cg;cfg%nrho_l0=rl0;cfg%nP_l0=p0;cfg%nPs=p
    cfg%D_H=diam;cfg%AREA=area;cfg%eps=rough;cfg%G=g
    cfg%mu_l=mul;cfg%mu_g=mug;cfg%sigma=sigma;cfg%Lr=12.546254_wp
    cfg%mul=rl*.4_wp/wu;cfg%mug=rg*.2_wp/wu;cfg%njl=.4_wp/wu
    cfg%tol=tol;cfg%tola=tol;cfg%maxit=100;cfg%quad_points=3
    cfg%junction_model='MATCHED';cfg%flux_scheme='ROE'
  end subroutine
  subroutine check(yes,label)
    logical,intent(in)::yes
    character(*),intent(in)::label
    if(.not.yes)then
      write(*,'(a,i0)')'FAIL '//label//' status=',st
      error stop 1
    endif
    count=count+1
  end subroutine
  subroutine near(x,y,eps,label)
    real(wp),intent(in)::x,y,eps
    character(*),intent(in)::label
    call check(ieee_is_finite(x).and.abs(x-y)<=eps*max(1.0_wp,abs(y)),label)
  end subroutine
end program
