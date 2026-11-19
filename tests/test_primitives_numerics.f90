program test_primitives_numerics
  use mod_kinds, only: wp,i4
  use module_find
  use module_closure_law, only: lei_fechamento_or_ndim_simp,CdUd_bendiksen
  use module_jacobian, only: DFduj
  use module_matrix_roe, only: calculate_roe_matrix
  use module_flow_regime, only: h_from_alphap,alphap_from_h
  use module_lapack_small, only: dgeev_3x3,inv3x3
  use module_roe_solver, only: roe_solver_t,roe_absA_from_A,roe_eigenstructure
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(wp),parameter :: pi=acos(-1.0_wp),diam=.0508_wp,area=pi*diam**2/4.0_wp
  real(wp),parameter :: rough=4.6e-5_wp,grav=9.81_wp,mul=.89e-3_wp,mug=1.85e-5_wp,sig=.0728_wp
  real(wp),parameter :: tol=1e-12_wp,cl_si=1500.0_wp,cg_si=300.0_wp,p0_si=1e5_wp,rho0_si=1000.0_wp
  real(wp) :: wu,wrho,wpres,cl,cg,rho0,p0,u(3),q(6),q2(6),f(3),f2(3),u2(3),up(3),um(3)
  real(wp) :: jac(3,3),jfd(3,3),roe(3,3),roe2(3,3),a(3,3),r(3,3),l(3,3),ident(3,3),ainv(3,3)
  real(wp) :: eig(3),imag(3),absa(3,3),amax,cd,ud,alpha,p,rl,rg,ul,ug,th,h,val,scale(3)
  real(wp) :: aa(7),pp(3),fractions(5),oldq(6),oldj(3,3),physj(3,3),rs(3),cs(3)
  integer(i4) :: st,idx,count,k,j,m,cid
  logical :: ok
  type(roe_solver_t) :: solver
  count=0
  call scales(.5_wp,1000.0_wp,1e5_wp)
  fractions=[0.0_wp,.2_wp,.5_wp,.8_wp,1.0_wp]
  do k=1,size(fractions)
    h=h_from_alphap(fractions(k),1e-13_wp)
    call close(alphap_from_h(h),fractions(k),2e-12_wp,'decreasing gas-fraction geometry round trip')
  end do
  call close(h_from_alphap(0.0_wp,tol),1.0_wp,1e-14_wp,'zero gas has full liquid level')
  call close(h_from_alphap(1.0_wp,tol),0.0_wp,1e-14_wp,'all gas has zero liquid level')
  aa=[1e-12_wp,1e-8_wp,.01_wp,.25_wp,.9_wp,1.0_wp-1e-8_wp,1.0_wp-1e-12_wp]
  pp=[1e-4_wp,2.0_wp,1e3_wp]
  do k=1,size(aa)
    do j=1,size(pp)
      alpha=aa(k);p=pp(j)
      rl=rho0+(p-p0)/(cl*cl);rg=p/(cg*cg)
      u=[(1-alpha)*rl,alpha*rg,0.0_wp]
      call close(pressure_fun_uj(u(1),u(2),cg,cl,rho0,p0),p,5e-8_wp,'stable EOS quadratic')
      call close(alpha_fun_uj(u(1),u(2),cg,cl,rho0,p0,tol),alpha,1e-10_wp,'EOS gas mass retained')
    end do
  end do
  u=[-.1_wp,.001_wp,1.0_wp]
  call recover(u,.6_wp,2,q,st,idx)
  call check(st/=0.and..not.ieee_is_finite(q(1)),'negative phase mass is rejected without clipping')

  ! Analytically generated, admissible drift states, including reverse liquid.
  do k=1,3
    alpha=.25_wp;p=2.0_wp;th=.6_wp
    rl=rho0+(p-p0)/(cl*cl);rg=p/(cg*cg)
    select case(k)
    case(1)
      ul=.4_wp
    case(2)
      ul=-.6_wp
    case default
      ul=20.0_wp
    end select
    call CdUd_bendiksen(alpha,rl,rg,ul,ul,th,diam,area,rough,grav,mul,mug,sig,wu,wrho,tol,cd,ud,st)
    call check(st==0,'Bendiksen coefficients')
    ug=(cd*(1-alpha)*ul+ud/wu)/(1-cd*alpha)
    u=[(1-alpha)*rl,alpha*rg,(1-alpha)*rl*ul+alpha*rg*ug]
    call recover(u,th,2,q,st,idx)
    call check(st==0.and.idx==1,'drift recovery status separate from regime')
    call close(q(5),ul,1e-10_wp,'liquid velocity round trip')
    call close(q(6),ug,1e-10_wp,'gas velocity round trip')
    call invariants(u,q,th,2)
    call derivative_check(u,th,2)
  end do

  ! Bhagwat transcription (Azevedo/Balino/Burr 2015): nonlinear coefficients
  ! depend on alpha and velocities. Generate states with ul fixed and solve
  ! for ug, then independently recover velocities with total momentum fixed.
  do k=1,4
    fractions(:4)=[.01_wp,.25_wp,.6_wp,.8_wp]
    alpha=fractions(k)
    p=2.0_wp;th=.35_wp*real(k-1,wp)
    ul=.04_wp*10.0_wp**real(mod(k-1,3),wp)
    call bg_state(alpha,p,ul,th,u,ug)
    call recover(u,th,3,q,st,idx)
    call check(st==0.and.idx==1,'Bhagwat conservative inversion status')
    call close(q(5),ul,1e-9_wp,'Bhagwat fixed-liquid generated round trip')
    call close(q(6),ug,1e-9_wp,'Bhagwat gas velocity round trip')
    call invariants(u,q,th,3)
    call derivative_check(u,th,3)
  end do
  call bg_state(.25_wp,2.0_wp,0.0_wp,.6_wp,u,ug)
  call recover(u,.6_wp,3,q,st,idx)
  call check(st==0,'Bhagwat stationary-liquid boundary is admissible')
  call close(q(5),0.0_wp,1e-9_wp,'stationary-liquid boundary is not forced negative')
  call close(q(6),ug,1e-9_wp,'stationary-liquid boundary gas speed')
  u(3)=-abs(u(3))
  call recover(u,.6_wp,3,q,st,idx)
  call check(st/=0.and..not.ieee_is_finite(q(5)),'Bhagwat unsupported countercurrent state fails explicitly')

  ! Positive trace gas remains a solved phase, even when the full admissible
  ! gas-speed endpoint from total momentum is billions of times larger.
  aa(:3)=[1e-10_wp,1e-8_wp,1e-5_wp]
  do k=1,3
    call bg_state(aa(k),2.0_wp,.4_wp,.6_wp,u,ug)
    call recover(u,.6_wp,3,q,st,idx)
    call check(st==0,'Bhagwat small positive gas fraction recovery')
    call close(q(1)/aa(k),1.0_wp,2e-10_wp,'trace gas mass is retained')
    call close(q(5),.4_wp,1e-9_wp,'trace gas liquid velocity round trip')
    call close(q(6),ug,1e-9_wp,'trace gas velocity round trip')
    call invariants(u,q,.6_wp,3)
    call bg_value(q(1),q(2),q(3),q(5),q(6),.6_wp,val,st)
    call check(st==0.and.abs(val/q(1))<1e-9_wp,'trace gas normalized closure residual')
  enddo
  call bg_state(.4_wp,2.0_wp,.4_wp,.6_wp,u,ug)
  call scale_check(u,.6_wp,3)
  rl=rho0+(2.0_wp-p0)/(cl*cl);rg=2.0_wp/(cg*cg)
  jac=DFduj(0.0_wp,rl,rg,.4_wp,ug,cl,cg,.6_wp,diam,area,rough,grav,mul,mug,sig, &
             wu,wrho,wpres,tol,st,3)
  call check(st/=0.and..not.all(ieee_is_finite(jac)),'exact zero-gas two-phase Jacobian is explicit domain failure')
  jac=DFduj(1.0_wp,rl,rg,.4_wp,ug,cl,cg,.6_wp,diam,area,rough,grav,mul,mug,sig, &
             wu,wrho,wpres,tol,st,3)
  call check(st/=0.and..not.all(ieee_is_finite(jac)),'exact zero-liquid two-phase Jacobian is explicit domain failure')

  call bg_state(.25_wp,2.0_wp,.4_wp,.6_wp,u,ug)
  call bg_state(.251_wp,2.002_wp,.402_wp,.6_wp,u2,ug)
  roe=calculate_roe_matrix(u,u2,cl,cg,.6_wp,rho0,p0,diam,area,rough,grav,mul,mug,sig,0.0_wp,wu,wrho,wpres, &
                         tol,tol,5_i4,status=st,closure_id=3)
  call check(st==0,'Bhagwat Roe path recovery')
  call flux(u,.6_wp,3,f,st)
  call check(st==0,'Bhagwat left physical flux')
  call flux(u2,.6_wp,3,f2,st)
  call check(st==0,'Bhagwat right physical flux')
  call check(maxval(abs(matmul(roe,u2-u)-(f2-f))/max(1.0_wp,abs(f2-f)))<1e-7_wp,'Bhagwat smooth Roe secant identity')
  call roe_absA_from_A(roe,absa,amax,st)
  call check(st==0.and.all(ieee_is_finite(absa)),'Bhagwat representative Roe path is hyperbolic')
  call bhagwat_domain_boundary_regression()

  ! The pipe is stratified by default; countercurrent roots are admissible.
  alpha=.3_wp;p=2.0_wp;th=-.02_wp
  rl=rho0+(p-p0)/(cl*cl);rg=p/(cg*cg)
  u=[(1-alpha)*rl,alpha*rg,2.0_wp]
  call recover(u,th,1,q,st,idx)
  call check(st==0.and.idx==0,'valid stratified index zero')
  call compute_primitive_variables(u(1),u(2),u(3),th,cg,cl,rho0,p0,diam,area,rough,grav,mul,mug,sig, &
      0.0_wp,wu,wrho,tol,tol,100_i4,q2(1),q2(2),q2(3),q2(4),q2(5),q2(6),idx,status=st)
  call check(st==0.and.idx==0,'default pipe selects fixed stratified closure')
  call check(maxval(abs(q-q2))<1e-10_wp,'default and explicit closure agree')
  call set_stability_params(enable_switch=.true.)
  call compute_primitive_variables(u(1),u(2),u(3),th,cg,cl,rho0,p0,diam,area,rough,grav,mul,mug,sig, &
      0.0_wp,wu,wrho,tol,tol,100_i4,q2(1),q2(2),q2(3),q2(4),q2(5),q2(6),idx,status=st)
  call check(st==PRIMITIVE_UNSUPPORTED_SWITCH,'unqualified regime switch is explicit unsupported')
  call set_stability_params(enable_switch=.false.)
  call invariants(u,q,th,1)
  call derivative_check(u,th,1)
  oldq=q
  oldj=DFduj(q(1),q(2),q(3),q(5),q(6),cl,cg,th,diam,area,rough,grav,mul,mug,sig,wu,wrho,wpres,tol,st,1)
  rs=[wrho*wu,wrho*wu,wrho*wu*wu];cs=[wrho,wrho,wrho*wu]
  do j=1,3
    oldj(:,j)=oldj(:,j)*rs/cs(j)
  end do
  ! Change all arbitrary scales while retaining the identical SI state.
  u2=u*[wrho,wrho,wrho*wu]
  oldq=oldq*[1.0_wp,wrho,wrho,wpres,wu,wu]
  call scales(2.0_wp,1.0_wp,1e6_wp)
  u2=u2/[wrho,wrho,wrho*wu]
  call recover(u2,th,1,q2,st,idx)
  call check(st==0,'recovery under changed scales')
  q2=q2*[1.0_wp,wrho,wrho,wpres,wu,wu]
  call check(maxval(abs(q2-oldq)/max(abs(oldq),1.0_wp))<1e-8_wp,'physical primitives invariant under scale changes')
  call recover(u2,th,1,q,st,idx)
  physj=DFduj(q(1),q(2),q(3),q(5),q(6),cl,cg,th,diam,area,rough,grav,mul,mug,sig,wu,wrho,wpres,tol,st,1)
  rs=[wrho*wu,wrho*wu,wrho*wu*wu];cs=[wrho,wrho,wrho*wu]
  do j=1,3
    physj(:,j)=physj(:,j)*rs/cs(j)
  end do
  call check(maxval(abs(physj-oldj)/max(1.0_wp,abs(oldj)))<2e-5_wp,'physical flux Jacobian scale covariance')
  call scales(.5_wp,1000.0_wp,1e5_wp)

  ! Roe quadrature is independent of its path parameterization.
  th=.6_wp;alpha=.25_wp;p=2.0_wp
  rl=rho0+(p-p0)/(cl*cl);rg=p/(cg*cg)
  u=[(1-alpha)*rl,alpha*rg,.5_wp]
  u2=u*[1.001_wp,.999_wp,1.002_wp]
  roe=calculate_roe_matrix(u,u2,cl,cg,th,rho0,p0,diam,area,rough,grav,mul,mug,sig,0.0_wp,wu,wrho,wpres, &
                         tol,tol,5_i4,status=st,closure_id=2)
  call check(st==0,'Roe quadrature status')
  roe2=calculate_roe_matrix(u,u2,cl,cg,th,rho0,p0,diam,area,rough,grav,mul,mug,sig,0.0_wp,wu,wrho,wpres, &
                         tol,tol,5_i4,a=2.0_wp,b=7.0_wp,status=st,closure_id=2)
  call check(st==0.and.maxval(abs(roe2-roe))<1e-8_wp,'normalized Roe quadrature interval')
  call flux(u,th,2,f,st)
  call flux(u2,th,2,f2,st)
  call check(maxval(abs(matmul(roe,u2-u)-(f2-f))/max(1.0_wp,abs(f2-f)))<1e-7_wp,'Roe path secant identity on smooth branch')
  call recover(u,th,2,q,st,idx)
  jac=DFduj(q(1),q(2),q(3),q(5),q(6),cl,cg,th,diam,area,rough,grav,mul,mug,sig,wu,wrho,wpres,tol,st,2)
  roe=calculate_roe_matrix(u,u,cl,cg,th,rho0,p0,diam,area,rough,grav,mul,mug,sig,0.0_wp,wu,wrho,wpres, &
                         tol,tol,3_i4,status=st,closure_id=2)
  call check(st==0.and.maxval(abs(roe-jac))/maxval(abs(jac))<1e-12_wp,'equal-state Roe consistency')
  call solver%roe_riemann_solver(f,f,u,u,roe,f2,st)
  call check(st==0.and.maxval(abs(f2-f))<1e-12_wp,'uniform-state numerical flux')

  ident=0.0_wp
  do j=1,3
    ident(j,j)=1.0_wp
  end do
  call roe_eigenstructure(ident,eig,imag,r,l,ok,st)
  call check(ok.and.st==0,'repeated identity eigenvalue has full basis')
  call check(maxval(abs(matmul(r,l)-ident))<1e-14_wp,'repeated-eigenvalue inverse')
  a=ident*1e8_wp;a(1,1)=a(1,1)+1;a(2,2)=a(2,2)+2;a(3,3)=a(3,3)+3
  call dgeev_3x3(a,eig,imag,r,st)
  call check(st==0.and.maxval(abs(imag))==0.0_wp,'translated real spectrum remains real')
  call check(maxval(abs(eig-[1e8_wp+1,1e8_wp+2,1e8_wp+3]))<1e-7_wp,'translated eigenvalue gaps retained')
  a=reshape([1.5_wp,-.5_wp,0.0_wp,-.5_wp,1.5_wp,0.0_wp,0.0_wp,0.0_wp,3.0_wp],[3,3])+1e8_wp*ident
  call roe_eigenstructure(a,eig,imag,r,l,ok,st)
  call check(ok.and.st==0,'translated nondiagonal symmetric spectrum')
  call check(maxval(abs(eig-[1e8_wp+1,1e8_wp+2,1e8_wp+3]))<1e-6_wp,'translated nondiagonal eigenvalue gaps')
  a=0.0_wp;a(1,1)=-100.0_wp;a(2,2)=.2_wp;a(3,3)=100.0_wp
  call roe_absA_from_A(a,absa,amax,st)
  call check(st==0.and.abs(absa(2,2)-.2_wp)<1e-14_wp,'slow material mode avoids acoustic-scale dissipation')
  ainv=inv3x3(1e-5_wp*ident,st)
  call check(st==0.and.maxval(abs(1e-5_wp*ainv-ident))<1e-14_wp,'inverse not rejected by determinant magnitude')
  a=0.0_wp
  call roe_absA_from_A(a,absa,amax,st)
  call check(st==0.and.maxval(abs(absa))==0.0_wp.and.amax==0.0_wp,'zero operator has zero dissipation')
  a(1,2)=-1.0_wp;a(2,1)=1.0_wp;a(3,3)=2.0_wp
  call roe_absA_from_A(a,absa,amax,st)
  call check(st==11.and..not.all(ieee_is_finite(absa)),'complex spectrum explicitly rejected')
  call recover(u,th,99,q,st,idx)
  call check(st/=0.and..not.ieee_is_finite(q(5)),'invalid closure identifier has explicit failure')
  write(*,'(a,i0,a)') 'PASS primitives and numerics: ',count,' assertions'
contains
  subroutine scales(us,rs,ps)
    real(wp),intent(in)::us,rs,ps
    wu=us;wrho=rs;wpres=ps
    cl=cl_si*sqrt(wrho/wpres);cg=cg_si*sqrt(wrho/wpres)
    rho0=rho0_si/wrho;p0=p0_si/wpres
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
    if(abs(x-y)>eps*max(1.0_wp,abs(x),abs(y)))write(*,*)label,x,y
    call check(abs(x-y)<=eps*max(1.0_wp,abs(x),abs(y)),label)
  end subroutine
  subroutine recover(uu,theta,model,qq,istat,regime)
    real(wp),intent(in)::uu(3),theta
    integer(i4),intent(in)::model
    real(wp),intent(out)::qq(6)
    integer(i4),intent(out)::istat,regime
    call compute_primitive_variables(uu(1),uu(2),uu(3),theta,cg,cl,rho0,p0,diam,area,rough,grav,mul,mug,sig, &
      0.0_wp,wu,wrho,tol,tol,100_i4,qq(1),qq(2),qq(3),qq(4),qq(5),qq(6),regime,status=istat,closure_id=model)
  end subroutine
  subroutine invariants(uu,qq,theta,model)
    real(wp),intent(in)::uu(3),qq(6),theta
    integer(i4),intent(in)::model
    real(wp)::rebuilt(3),resid
    integer::istat
    rebuilt=[(1-qq(1))*qq(2),qq(1)*qq(3),(1-qq(1))*qq(2)*qq(5)+qq(1)*qq(3)*qq(6)]
    call check(maxval(abs(rebuilt-uu)/max(1e-12_wp,abs(uu)))<1e-10_wp,'phase masses and momentum round trip')
    call lei_fechamento_or_ndim_simp(qq(1),qq(2),qq(3),qq(5),qq(6),theta,diam,area,rough,grav,mul,mug,sig, &
                                    wu,wrho,tol,resid,model,istat)
    call check(istat==0.and.abs(resid)<1e-10_wp,'same closure residual after recovery')
  end subroutine
  subroutine flux(uu,theta,model,ff,istat)
    real(wp),intent(in)::uu(3),theta
    integer(i4),intent(in)::model
    real(wp),intent(out)::ff(3)
    integer(i4),intent(out)::istat
    real(wp)::qq(6)
    integer::regime
    call recover(uu,theta,model,qq,istat,regime)
    if(istat/=0)return
    ff=[uu(1)*qq(5),uu(2)*qq(6),uu(1)*qq(5)**2+uu(2)*qq(6)**2+wpres/(wrho*wu**2)*qq(4)]
  end subroutine
  subroutine bg_state(a,p,ul_known,theta,uu,ug_known)
    real(wp),intent(in)::a,p,ul_known,theta
    real(wp),intent(out)::uu(3),ug_known
    real(wp)::rl,rg,lo,hi,mid,flo,fhi,fmid
    integer::istat,it
    rl=rho0+(p-p0)/(cl*cl);rg=p/(cg*cg)
    lo=0.0_wp;hi=max(1.0_wp,2.0_wp*ul_known)
    if(ul_known==0.0_wp)lo=1e-5_wp
    call bg_value(a,rl,rg,ul_known,lo,theta,flo,istat)
    call check(istat==0,'Bhagwat generator lower velocity endpoint')
    do it=1,30
      call bg_value(a,rl,rg,ul_known,hi,theta,fhi,istat)
      call check(istat==0,'Bhagwat generator positive velocity endpoint')
      if(sign(1.0_wp,flo)/=sign(1.0_wp,fhi))exit
      hi=2*hi
    enddo
    call check(it<=30,'Bhagwat independent fixed-liquid root bracket')
    do it=1,80
      mid=.5_wp*(lo+hi)
      call bg_value(a,rl,rg,ul_known,mid,theta,fmid,istat)
      if(istat/=0)error stop 'Bhagwat generator residual'
      if(fmid==0.0_wp.or.hi-lo<16*epsilon(1.0_wp)*max(1.0_wp,abs(mid)))exit
      if(sign(1.0_wp,fmid)==sign(1.0_wp,flo))then
        lo=mid;flo=fmid
      else
        hi=mid
      endif
    enddo
    ug_known=mid
    uu=[(1-a)*rl,a*rg,(1-a)*rl*ul_known+a*rg*ug_known]
  end subroutine bg_state
  subroutine bg_value(a,rl,rg,vl,vg,theta,val,status)
    real(wp),intent(in)::a,rl,rg,vl,vg,theta
    real(wp),intent(out)::val
    integer,intent(out)::status
    call lei_fechamento_or_ndim_simp(a,rl,rg,vl,vg,theta,diam,area,rough,grav,mul,mug,sig, &
                                    wu,wrho,tol,val,3,status)
  end subroutine bg_value
  subroutine scale_check(uu,theta,model)
    real(wp),intent(in)::uu(3),theta
    integer,intent(in)::model
    real(wp)::physical_u(3),physical_q(6),changed_q(6),scaled_u(3),qq(6)
    real(wp)::physical_jac(3,3),changed_jac(3,3),row_scales(3),column_scales(3)
    real(wp)::saved_wu,saved_wrho,saved_wpres
    integer::istat,regime,jj
    saved_wu=wu;saved_wrho=wrho;saved_wpres=wpres
    call recover(uu,theta,model,qq,istat,regime)
    call check(istat==0,'scale covariance base recovery')
    physical_u=uu*[wrho,wrho,wrho*wu]
    physical_q=qq*[1.0_wp,wrho,wrho,wpres,wu,wu]
    physical_jac=DFduj(qq(1),qq(2),qq(3),qq(5),qq(6),cl,cg,theta,diam,area,rough,grav,mul,mug,sig, &
                         wu,wrho,wpres,tol,istat,model)
    call check(istat==0,'scale covariance base Jacobian')
    row_scales=[wrho*wu,wrho*wu,wrho*wu*wu];column_scales=[wrho,wrho,wrho*wu]
    do jj=1,3
      physical_jac(:,jj)=physical_jac(:,jj)*row_scales/column_scales(jj)
    enddo
    call scales(2.0_wp,1.0_wp,1e6_wp)
    scaled_u=physical_u/[wrho,wrho,wrho*wu]
    call recover(scaled_u,theta,model,qq,istat,regime)
    call check(istat==0,'scale covariance changed recovery')
    changed_q=qq*[1.0_wp,wrho,wrho,wpres,wu,wu]
    call check(maxval(abs(changed_q-physical_q)/max(1.0_wp,abs(physical_q)))<1e-9_wp, &
               'Bhagwat physical primitives invariant under three scale changes')
    changed_jac=DFduj(qq(1),qq(2),qq(3),qq(5),qq(6),cl,cg,theta,diam,area,rough,grav,mul,mug,sig, &
                        wu,wrho,wpres,tol,istat,model)
    call check(istat==0,'scale covariance changed Jacobian')
    row_scales=[wrho*wu,wrho*wu,wrho*wu*wu];column_scales=[wrho,wrho,wrho*wu]
    do jj=1,3
      changed_jac(:,jj)=changed_jac(:,jj)*row_scales/column_scales(jj)
    enddo
    call check(maxval(abs(changed_jac-physical_jac)/max(1.0_wp,abs(physical_jac)))<1e-7_wp, &
               'Bhagwat physical flux Jacobian scale covariance')
    call scales(saved_wu,saved_wrho,saved_wpres)
  end subroutine scale_check
  subroutine bhagwat_domain_boundary_regression()
    ! Frozen SS3 pilot, cell41, t=94.48170685397338s. This is a domain
    ! regression, not an expectation that the limited closure finishes SS3.
    real(wp),parameter::rho_scale=997.05_wp,u_scale=.48_wp,p_scale=101325._wp
    real(wp),parameter::d=.05254_wp,theta=1.1918637472270275_wp
    real(wp)::uu(3),base(3),step(3),qq(6),lambda,cln,cgn,pressure,jj(3,3)
    integer::istat,regime,k
    cln=1498._wp*sqrt(rho_scale/p_scale);cgn=343._wp*sqrt(rho_scale/p_scale)
    base=[.47423216860070572_wp,.00099880267138943864_wp,.0010366845397290843_wp]
    step=[3.1352300711581445e-12_wp,-1.7734242497486271e-12_wp,-1.0961662827351605e-8_wp]
    call compute_primitive_variables(base(1),base(2),base(3),theta,cgn,cln,1._wp,1._wp,d,pi*d*d/4, &
      rough,grav,mul,mug,sig,0._wp,u_scale,rho_scale,1e-12_wp,1e-13_wp,100_i4, &
      qq(1),qq(2),qq(3),qq(4),qq(5),qq(6),regime,status=istat,closure_id=3)
    call check(istat==0,'SS3 cocurrent boundary base is admissible')
    call check(qq(5)>=0._wp.and.qq(5)*u_scale<1e-10_wp,'SS3 liquid reaches zero speed with gas present')
    call check(qq(1)>.5_wp.and.qq(1)<.6_wp.and.qq(6)>0._wp,'SS3 domain boundary is not phase disappearance')
    jj=DFduj(qq(1),qq(2),qq(3),qq(5),qq(6),cln,cgn,theta,d,pi*d*d/4,rough,grav,mul,mug,sig, &
                 u_scale,rho_scale,p_scale,1e-12_wp,istat,3)
    call check(istat==0.and.all(ieee_is_finite(jj)),'SS3 boundary Jacobian stays finite on admitted side')
    do k=1,2
      lambda=1._wp
      if(k==2)lambda=1e-4_wp
      uu=base+lambda*step
      pressure=pressure_fun_uj(uu(1),uu(2),cgn,cln,1._wp,1._wp)
      call check(minval(uu(:2))>0._wp.and.ieee_is_finite(pressure).and.pressure>0._wp, &
                 'SS3 rejected Newton trial still has positive phase masses and pressure')
      call compute_primitive_variables(uu(1),uu(2),uu(3),theta,cgn,cln,1._wp,1._wp,d,pi*d*d/4, &
        rough,grav,mul,mug,sig,0._wp,u_scale,rho_scale,1e-12_wp,1e-13_wp,100_i4, &
        qq(1),qq(2),qq(3),qq(4),qq(5),qq(6),regime,status=istat,closure_id=3)
      call check(istat==PRIMITIVE_NO_PHYSICAL_ROOT,'positive lower bound excludes a cocurrent SS3 root')
      call check(.not.ieee_is_finite(qq(5)).and..not.ieee_is_finite(qq(6)), &
                 'SS3 unsupported reversal does not fabricate phase velocities')
    enddo
    ! Full Newton proposal at extended cell 69, t=97.1318523387 s,
    ! frozen 8a85e37 run. The independent 80-digit certificate bounds the
    ! entire q=[0,1] residual below by +5.984975257e-7 m/s.
    uu=[.361748887683525366_wp,.00124300510327740479_wp,.00161735202739144292_wp]
    call compute_primitive_variables(uu(1),uu(2),uu(3),.901238667226846091_wp,cgn,cln,1._wp,1._wp, &
      d,pi*d*d/4,rough,grav,mul,mug,sig,0._wp,u_scale,rho_scale,1e-12_wp,1e-13_wp,100_i4, &
      qq(1),qq(2),qq(3),qq(4),qq(5),qq(6),regime,status=istat,closure_id=3)
    call check(istat==PRIMITIVE_NO_PHYSICAL_ROOT,'cell 69 full proposal excluded over the momentum interval')
    call check(abs(qq(1)-.638271663623163248_wp)<1e-12_wp,'cell 69 has positive inventories in both phases')
  end subroutine bhagwat_domain_boundary_regression
  subroutine derivative_check(uu,theta,model)
    real(wp),intent(in)::uu(3),theta
    integer(i4),intent(in)::model
    real(wp)::qq(6),aj(3,3),fd(3,3),plus(3),minus(3),fp(3),fm(3),hh
    integer(i4)::istat,regime,jj
    call recover(uu,theta,model,qq,istat,regime)
    call check(istat==0,'Jacobian base recovery')
    aj=DFduj(qq(1),qq(2),qq(3),qq(5),qq(6),cl,cg,theta,diam,area,rough,grav,mul,mug,sig, &
               wu,wrho,wpres,tol,status=istat,closure_id=model)
    call check(istat==0,'implicit Jacobian status')
    do jj=1,3
      hh=1e-5_wp*max(abs(uu(jj)),1e-5_wp)
      plus=uu;minus=uu;plus(jj)=plus(jj)+hh;minus(jj)=minus(jj)-hh
      call flux(plus,theta,model,fp,istat)
      call check(istat==0,'plus Jacobian recovery')
      call flux(minus,theta,model,fm,istat)
      call check(istat==0,'minus Jacobian recovery')
      fd(:,jj)=(fp-fm)/(2*hh)
    end do
    if(maxval(abs(aj-fd)/max(1.0_wp,abs(fd)))>3e-5_wp)then
      write(*,*)'Jacobian analytic:',aj
      write(*,*)'Jacobian finite difference:',fd
    endif
    call check(maxval(abs(aj-fd)/max(1.0_wp,abs(fd)))<3e-5_wp,'flux Jacobian versus conservative finite differences')
    call check(maxval(abs(aj(1,:)+aj(2,:)-[0.0_wp,0.0_wp,1.0_wp]))<1e-10_wp,'d(F1+F2)/dU=(0,0,1)')
  end subroutine
end program test_primitives_numerics
