program test_junction_phase_upwind
 !! EXPERIMENTAL junction vector E1 (junction_flux='PHASE_UPWIND', opt-in).
 !! A. The vector alone, on synthetic traces: consistency for equal traces,
 !!    reflection antisymmetry, positive/negative/zero phase velocities
 !!    (continuity through zero, upwind side selection), finiteness.
 !! B. The MATCHED junction of an independently integrated stationary profile:
 !!    transparency when the riser trace is the riser-closure image of the
 !!    pipe trace (artificial partition), stationary recovery with
 !!    STEADY traces, one mass flux per phase, exact telescoping of the
 !!    assembled residual, and the assembled Jacobian against central FD.
 !! The default 'ROE' must remain the operator unless E1 is requested.
 use mod_kinds,only:wp
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t,junction_phase_upwind
 use module_config,only:sim_config_t
 use module_steady_state,only:EstadoEstacionario_ndim_simp,solve_alpha_from_superficials
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::s(:),p(:),rg(:),rl(:),alpha(:),ug(:),ul(:),th(:),jl(:),jg(:),faces(:),pf(:),u(:,:)
 real(wp),allocatable::res(:,:),rplus(:,:),rminus(:,:),az(:,:,:),bz(:,:,:),cz(:,:,:),daz(:,:,:),dbz(:,:,:),dcz(:,:,:),direction(:,:),work(:,:)
 real(wp)::va(6),vb(6),ua(3),ub(3),gl,gg,ps,gl2,gg2,ps2,pip,fm(3),fp(3),tr(6),tr2(6),norm,nt,fin(3),fout(3),h
 real(wp)::prediction(3),observed(3),jac_error,pj,y(3),e(3),delta
 integer::st,n,checks,i,jj,ib,model,k
 checks=0
 ! ---------------- A. vector alone ----------------
 cfg=sim_config_t()
 cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1.16_wp
 cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=290.06_wp*sqrt(cfg%omega_rho/cfg%omega_P)
 pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
 call trace(va,ua,.56_wp,2.78_wp,1.04_wp,.046_wp);call trace(vb,ub,.56_wp,2.78_wp,1.04_wp,.046_wp)
 call junction_phase_upwind(cfg,va,ua,vb,ub,gl,gg,ps,st)
 e=[abs(gl-ua(1)*va(5))/abs(ua(1)*va(5)),abs(gg-ua(2)*va(6))/abs(ua(2)*va(6)),abs(ps-va(4))/va(4)]
 write(*,'(a,3es11.3)')'A1 equal traces, rel err Gl,Gg,p* ',e
 call check(st==0.and.maxval(e)<=4*epsilon(1._wp),'consistency: equal traces give the physical vector')
 do k=1,3
  ! k=1: both gas velocities positive with a 73 Pa centre pressure jump (SS1/OSC scale)
  ! k=2: mirrored (velocities reversed, sides swapped): fluxes change sign exactly
  ! k=3: both negative
  call trace(va,ua,.56_wp,2.7807_wp,1.20_wp,.053_wp);call trace(vb,ub,.30_wp,2.7800_wp,1.30_wp,.40_wp)
  if(k==1)then
   call junction_phase_upwind(cfg,va,ua,vb,ub,gl,gg,ps,st)
   write(*,'(a,2es12.4,a,es12.4)')'A2 (+,+) Gl*,Gg* ',gl,gg,'  left gas flux ',ua(2)*va(6)
   call check(st==0.and.gg>0.and.gl>0,'positive phase velocities: positive phase fluxes')
   call check(abs(gg-ua(2)*va(6))<=abs(gg-ub(2)*vb(6)),'co-directional: upstream (left) trace flux plus acoustic shift')
   call check(gg<=max(ua(2)*va(6),ub(2)*vb(6))*1.05_wp,'no amplification beyond the trace fluxes')
  else if(k==2)then
   call junction_phase_upwind(cfg,va,ua,vb,ub,gl,gg,ps,st)
   va(5:6)=-va(5:6);vb(5:6)=-vb(5:6);ua(3)=-ua(3);ub(3)=-ub(3)
   call junction_phase_upwind(cfg,vb,ub,va,ua,gl2,gg2,ps2,st)
   write(*,'(a,2es12.4,a,2es12.4)')'A3 reflection Gl,Gg ',gl,gg,' mirrored ',gl2,gg2
   call check(st==0.and.abs(gl+gl2)<=1e-14_wp*abs(gl).and.abs(gg+gg2)<=1e-14_wp*abs(gg).and.abs(ps-ps2)<=1e-14_wp*ps, &
     'reflection antisymmetry of the mass fluxes, symmetric p*')
  else
   va(5:6)=-va(5:6);vb(5:6)=-vb(5:6);ua(3)=-ua(3);ub(3)=-ub(3)
   call junction_phase_upwind(cfg,va,ua,vb,ub,gl,gg,ps,st)
   write(*,'(a,2es12.4)')'A4 (-,-) Gl*,Gg* ',gl,gg
   call check(st==0.and.gg<0.and.gl<0,'negative phase velocities: negative phase fluxes')
   call check(abs(gg-ub(2)*vb(6))<abs(gg-ua(2)*va(6)),'co-directional negative: upstream (right) trace flux')
  endif
 enddo
 ! zero gas velocity at the face approached from both signs: continuous, no clipping
 call trace(va,ua,.56_wp,2.78_wp,1.2_wp,0._wp);call trace(vb,ub,.30_wp,2.78_wp,1.2_wp,0._wp)
 do k=-1,1
  delta=real(k,wp)*1e-9_wp
  va(6)=delta;vb(6)=delta;ua(3)=ua(1)*va(5)+ua(2)*va(6);ub(3)=ub(1)*vb(5)+ub(2)*vb(6)
  call junction_phase_upwind(cfg,va,ua,vb,ub,gl,gg,ps,st)
  write(*,'(a,es10.2,a,es12.4)')'A5 u_g=',delta,'  Gg* ',gg
  call check(st==0.and.abs(gg)<=max(ua(2),ub(2))*(abs(delta)+1e-15_wp),'gas flux -> 0 continuously as u_g -> 0')
  if(k/=0)call check(sign(1._wp,gg)==sign(1._wp,delta),'sign of gas flux follows the gas velocity')
 enddo
 ! ---------------- B. MATCHED junction ----------------
 do model=2,3
  call setup(model)
  ib=cfg%num_dummy_cells+cfg%N_pipe;pj=pf(cfg%N_pipe+1)
  call check(trim(cfg%junction_flux)=='ROE','Roe is the default junction vector')
  ! B1 artificial partition: riser trace = riser-closure image of the pipe trace
  cfg%junction_flux='PHASE_UPWIND';call sim%init(cfg)
  y=riser_image(u(ib,:),cfg%thetav(ib),cfg%thetav(ib+1),model)
  call sim%face_pair(u(ib,:),y,cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
  call check(st==0,'E1 admissible on a transparent junction')
  e=[abs(tr(2)-u(ib,1)*ulp(u(ib,:)))/abs(tr(2)),abs(tr(3)-u(ib,2)*ugp(u(ib,:)))/abs(tr(3)),0._wp]
  write(*,'(a,i2,a,2es11.3)')'B1 closure_riser',model,' transparent junction rel err Gl*,Gg* ',e(:2)
  call check(maxval(e)<1e-8_wp,'transparent junction: phase fluxes of the common state')
  ! B2 stationary profile, CENTER and STEADY traces, Roe vs E1
  do k=1,2
   cfg%junction_trace=merge('CENTER','STEADY',k==1)
   cfg%junction_flux='ROE';call sim%init(cfg)
   call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr2,left_cell=ib)
   call check(st==0,'Roe junction admissible')
   cfg%junction_flux='PHASE_UPWIND';call sim%init(cfg)
   call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
   call check(st==0,'E1 junction admissible')
   call check(all(fm(:2)==fp(:2)),'one mass flux per phase shared by both cells')
   write(*,'(a,a,a,3es11.3,a,3es11.3)')'B2 ',trim(cfg%junction_trace),' rel err p*,Gl*,Gg*  Roe ', &
     abs(tr2(1)-pj)/pj,abs(tr2(2)-cfg%mul)/cfg%mul,abs(tr2(3)-cfg%mug)/cfg%mug,'  E1 ', &
     abs(tr(1)-pj)/pj,abs(tr(2)-cfg%mul)/cfg%mul,abs(tr(3)-cfg%mug)/cfg%mug
   if(k==2)call check(max(abs(tr(1)-pj)/pj,abs(tr(2)-cfg%mul)/cfg%mul,abs(tr(3)-cfg%mug)/cfg%mug)<1e-7_wp, &
     'E1 with steady traces recovers the stationary junction vector')
  enddo
  cfg%junction_trace='CENTER'
  ! B3 conservation and Jacobian of the assembled operator with E1
  allocate(res(3,n),rplus(3,n),rminus(3,n),az(3,3,n),bz(3,3,n),cz(3,3,n),daz(3,3,n),dbz(3,3,n),dcz(3,3,n), &
    direction(n+4,3),work(n+4,3))
  cfg%junction_flux='PHASE_UPWIND';call sim%init(cfg);call sim%prepare_dissipation(u,st)
  call check(st==0,'lagged dissipation with E1')
  call sim%assemble(u,u,.001_wp,res,az,bz,cz,norm,fin,fout,st,.true.)
  call check(st==0,'E1 residual and complete Jacobian')
  do i=1,2
   nt=sum(res(i,:)*cfg%cell_width(3:n+2))-.001_wp*(fout(i)-fin(i))
   write(*,'(a,i1,a,es12.4)')'B3 phase ',i,' telescoping defect ',nt
   call check(abs(nt)<1e-14_wp,'phase residual telescopes: junction flux cancels exactly')
  enddo
  direction=0
  direction(ib,:)=u(ib,:)*[.37_wp,-.23_wp,.11_wp];direction(ib+1,:)=u(ib+1,:)*[-.31_wp,.19_wp,.43_wp]
  direction(ib-1,:)=u(ib-1,:)*[.05_wp,.07_wp,-.03_wp]
  h=1e-6_wp;work=u+h*direction
  call sim%assemble(work,u,.001_wp,rplus,daz,dbz,dcz,nt,fin,fout,st,.false.)
  call check(st==0,'plus direction admissible')
  work=u-h*direction
  call sim%assemble(work,u,.001_wp,rminus,daz,dbz,dcz,nt,fin,fout,st,.false.)
  call check(st==0,'minus direction admissible')
  jac_error=0
  do jj=1,n
   prediction=matmul(bz(:,:,jj),direction(jj+2,:))
   if(jj>1)prediction=prediction+matmul(az(:,:,jj),direction(jj+1,:))
   if(jj<n)prediction=prediction+matmul(cz(:,:,jj),direction(jj+3,:))
   observed=(rplus(:,jj)-rminus(:,jj))/(2*h)
   jac_error=max(jac_error,maxval(abs(observed-prediction)/max(1e-3_wp,abs(prediction))))
  enddo
  write(*,'(a,es12.4)')'B3 E1 Jacobian directional error ',jac_error
  call check(jac_error<1e-5_wp,'assembled Jacobian consistent with the E1 residual')
  deallocate(res,rplus,rminus,az,bz,cz,daz,dbz,dcz,direction,work)
 enddo
 write(*,'(a,i0,a)')'PASS junction phase upwind (E1): ',checks,' assertions'
contains
 subroutine trace(v,uu,alpha_,p_,ul_,ug_)
  !! Synthetic trace from (alpha, p [nondim], u_l, u_g [nondim]).
  real(wp),intent(out)::v(6),uu(3)
  real(wp),intent(in)::alpha_,p_,ul_,ug_
  real(wp)::rl_,rg_
  rl_=1+(p_-1)/cfg%nCl**2;rg_=p_/cfg%nCg**2
  v=[alpha_,rl_,rg_,p_,ul_,ug_]
  uu=[(1-alpha_)*rl_,alpha_*rg_,(1-alpha_)*rl_*ul_+alpha_*rg_*ug_]
 end subroutine
 real(wp) function ulp(x)
  real(wp),intent(in)::x(3)
  real(wp)::v(6),f(3),sx(3),a(3,3),d(3,3)
  integer::c
  call sim%point_state(x,cfg%thetav(ib),v,f,sx,a,d,c,.false.,1)
  ulp=v(5)
 end function
 real(wp) function ugp(x)
  real(wp),intent(in)::x(3)
  real(wp)::v(6),f(3),sx(3),a(3,3),d(3,3)
  integer::c
  call sim%point_state(x,cfg%thetav(ib),v,f,sx,a,d,c,.false.,1)
  ugp=v(6)
 end function
 function riser_image(x,thl,thr,model) result(yy)
  real(wp),intent(in)::x(3),thl,thr
  integer,intent(in)::model
  real(wp)::yy(3),v(6),f(3),sx(3),a(3,3),d(3,3),aa
  integer::c
  call sim%point_state(x,thl,v,f,sx,a,d,c,.false.,1)
  aa=solve_alpha_from_superficials(f(1)/v(2),f(2)/v(3),v(2),v(3),thr,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
    cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,1e-13_wp,status=c,closure_id=model)
  if(c/=0)print *,'riser image status',c,' jl,jg',f(1)/v(2),f(2)/v(3),' theta_r',thr
  call check(c==0,'riser-closure image of the pipe trace exists')
  yy=[(1-aa)*v(2),aa*v(3),f(1)+f(2)]
 end function
 subroutine setup(model)
  integer,intent(in)::model
  cfg=sim_config_t()
  cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1;cfg%omega_t=12.546253993197219_wp
  cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp;cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%Lp=57.4_wp;cfg%Lr=cfg%omega_t;cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
  cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp
  cfg%closure_pipe=1;cfg%closure_riser=model;cfg%mul=.5_wp;cfg%mug=.001_wp;cfg%njl=.5_wp
  cfg%jac_dissipation_fd=.true.
  if(allocated(s))deallocate(s,p,rg,rl,alpha,ug,ul,th,jl,jg,faces,pf,u)
  call EstadoEstacionario_ndim_simp(9,cfg%mul,cfg%mug,cfg%nPs,cfg%Lp,cfg%Lr,cfg%CA,cfg%beta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
    cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_P,cfg%omega_u,cfg%omega_rho,cfg%tol, &
    s,p,rg,rl,alpha,ug,ul,th,jl,jg,n_riser=9,closure_id_riser=model,closure_id_pipe=1,status=st,face_positions=faces,face_pressures=pf)
  call check(st==0,'independent global stationary integration')
  n=size(s);cfg%nCells=n+4;cfg%N_pipe=9;cfg%num_dummy_cells=2
  allocate(u(n+4,3),cfg%vns(n+4),cfg%thetav(n+4),cfg%cell_width(n+4),cfg%face_positions(0:n))
  cfg%face_positions=faces
  do i=1,n
   u(i+2,:)=[(1-alpha(i))*rl(i),alpha(i)*rg(i),(1-alpha(i))*rl(i)*ul(i)+alpha(i)*rg(i)*ug(i)]
   cfg%vns(i+2)=s(i);cfg%thetav(i+2)=th(i);cfg%cell_width(i+2)=faces(i+1)-faces(i)
  enddo
  do i=1,2
   u(i,:)=u(3,:);u(n+2+i,:)=u(n+2,:);cfg%thetav(i)=th(1);cfg%thetav(n+2+i)=th(n)
   cfg%cell_width(i)=cfg%cell_width(3);cfg%cell_width(n+2+i)=cfg%cell_width(n+2)
   cfg%vns(i)=s(1)-(3-i)*cfg%cell_width(3);cfg%vns(n+2+i)=s(n)+i*cfg%cell_width(n+2)
  enddo
  cfg%U0=u
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
