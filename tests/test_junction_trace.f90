program test_junction_trace
 !! Well-balanced junction option (junction_trace='STEADY') on an independently
 !! integrated stationary pipeline-riser profile. The default 'CENTER' is the
 !! legacy operator and must stay unchanged; the opt-in trace must recover the
 !! stationary junction pressure and phase fluxes, keep one flux per phase and
 !! a complete Jacobian consistent with the residual.
 use mod_kinds,only:wp
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 use module_config,only:sim_config_t
 use module_steady_state,only:EstadoEstacionario_ndim_simp
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::s(:),p(:),rg(:),rl(:),alpha(:),ug(:),ul(:),th(:),jl(:),jg(:),faces(:),pf(:),u(:,:)
 real(wp),allocatable::res(:,:),rc(:,:),rplus(:,:),rminus(:,:),az(:,:,:),bz(:,:,:),cz(:,:,:),daz(:,:,:),dbz(:,:,:),dcz(:,:,:),direction(:,:),work(:,:)
 real(wp)::fm(3),fp(3),tr(6),ec(3),es(3),norm,nt,fin(3),fout(3),h,prediction(3),observed(3),jac_error,pj,gasc,gass
 integer::st,n,checks,i,jj,ib,model
 checks=0
 do model=2,3
  call setup(model)
  ib=cfg%num_dummy_cells+cfg%N_pipe;pj=pf(cfg%N_pipe+1)
  call check(trim(cfg%junction_trace)=='CENTER','legacy junction operator is the default')
  call sim%init(cfg)
  call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
  call check(st==0,'legacy MATCHED junction admissible on stationary profile')
  ec=[abs(tr(1)-pj)/pj,abs(tr(2)-cfg%mul)/cfg%mul,abs(tr(3)-cfg%mug)/cfg%mug]
  cfg%junction_trace='STEADY';call sim%init(cfg)
  call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
  call check(st==0,'steady junction trace admissible')
  es=[abs(tr(1)-pj)/pj,abs(tr(2)-cfg%mul)/cfg%mul,abs(tr(3)-cfg%mug)/cfg%mug]
  write(*,'(a,i2,a,3es12.4,a,3es12.4)')'closure_riser',model,' CENTER rel err p*,Gl*,Gg* ',ec,'  STEADY ',es
  call check(maxval(es)<1e-7_wp,'stationary junction pressure and phase fluxes recovered')
  call check(es(3)<ec(3)/100,'spurious junction gas-flux correction removed')
  call check(all(fm(:2)==fp(:2)),'one mass flux per phase shared by both cells')
  ! Residual of the stationary profile: the gas mass residual at the junction
  ! cells must not be dominated by the centre-to-centre pressure difference.
  allocate(res(3,n),rc(3,n),rplus(3,n),rminus(3,n),az(3,3,n),bz(3,3,n),cz(3,3,n),daz(3,3,n),dbz(3,3,n),dcz(3,3,n), &
    direction(n+4,3),work(n+4,3))
  cfg%junction_trace='CENTER';call sim%init(cfg);call sim%prepare_dissipation(u,st)
  call sim%assemble(u,u,.001_wp,rc,az,bz,cz,norm,fin,fout,st,.false.)
  call check(st==0,'legacy stationary residual')
  cfg%junction_trace='STEADY';call sim%init(cfg);call sim%prepare_dissipation(u,st)
  call check(st==0,'lagged dissipation with steady junction')
  call sim%assemble(u,u,.001_wp,res,az,bz,cz,norm,fin,fout,st,.true.)
  call check(st==0,'steady-junction residual and complete Jacobian')
  gasc=maxval(abs(rc(2,ib-2:ib-1)));gass=maxval(abs(res(2,ib-2:ib-1)))
  write(*,'(a,2es12.4)')'  junction-cell gas residual CENTER, STEADY ',gasc,gass
  call check(gass<gasc,'junction gas residual reduced on the stationary profile')
  nt=sum(res(2,:)*cfg%cell_width(3:n+2))-.001_wp*(fout(2)-fin(2))
  write(*,'(a,es12.4)')'  gas telescoping defect (residual-weighted minus boundary fluxes) ',nt
  call check(abs(nt)<1e-14_wp,'gas residual telescopes: junction flux cancels exactly')
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
  write(*,'(a,es12.4)')'  steady-junction Jacobian directional error ',jac_error
  call check(jac_error<1e-5_wp,'assembled Jacobian includes the junction traces')
  deallocate(res,rc,rplus,rminus,az,bz,cz,daz,dbz,dcz,direction,work)
 enddo
 write(*,'(a,i0,a)')'PASS junction trace: ',checks,' assertions'
contains
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
