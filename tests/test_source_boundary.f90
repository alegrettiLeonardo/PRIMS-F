program test_source_boundary
 use mod_kinds,only:wp
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 use module_config,only:sim_config_t
 use module_steady_state,only:EstadoEstacionario_ndim_simp
 use module_boundary_conditions,only:boundary_flux
 use module_boundary_trace,only:check_boundary_characteristics,BC_DATA_REQUIRED
 use module_find,only:pressure_fun_uj
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::res(:,:),rplus(:,:),rminus(:,:),az(:,:,:),bz(:,:,:),cz(:,:,:),daz(:,:,:),dbz(:,:,:),dcz(:,:,:),direction(:,:),work(:,:)
 real(wp)::resnorm,nt,fin(3),fout(3),prediction(3),observed(3),jac_error
 integer::jj,kk
 real(wp),allocatable::s(:),p(:),rg(:),rl(:),alpha(:),ug(:),ul(:),th(:),jl(:),jg(:),faces(:),pf(:),u(:,:)
 real(wp)::f(3),bc(3),legacy_in_error,trace_in_error,legacy_out_error,trace_out_error,up(3),um(3),fp(3),fm(3),fhp(3),fhm(3),h,j2(3),j4(3),rln,rgn
 integer::st,n,checks,i,k,side,ih
 checks=0;cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1;cfg%omega_t=12.546253993197219_wp
 cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp;cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%Lp=57.4_wp;cfg%Lr=cfg%omega_t;cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
 cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
 cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp
 cfg%closure_pipe=1;cfg%closure_riser=2;cfg%mul=.5_wp;cfg%mug=.001_wp;cfg%njl=.5_wp
 call EstadoEstacionario_ndim_simp(9,cfg%mul,cfg%mug,cfg%nPs,cfg%Lp,cfg%Lr,cfg%CA,cfg%beta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
   cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_P,cfg%omega_u,cfg%omega_rho,cfg%tol, &
   s,p,rg,rl,alpha,ug,ul,th,jl,jg,n_riser=9,closure_id_riser=2,closure_id_pipe=1,status=st,face_positions=faces,face_pressures=pf)
 call check(st==0,'independent global stationary integration')
 n=size(s);cfg%nCells=n+4;cfg%N_pipe=9;cfg%num_dummy_cells=2
 allocate(u(n+4,3),cfg%vns(n+4),cfg%thetav(n+4),cfg%cell_width(n+4),cfg%face_positions(0:n))
 cfg%face_positions=faces
 do i=1,n
  u(i+2,:)=[(1-alpha(i))*rl(i),alpha(i)*rg(i),(1-alpha(i))*rl(i)*ul(i)+alpha(i)*rg(i)*ug(i)]
  cfg%vns(i+2)=s(i);cfg%thetav(i+2)=th(i);cfg%cell_width(i+2)=faces(i+1)-faces(i)
 enddo
 call boundary_flux(u(3,:),th(1),cfg,.true.,f,st,bc);call check(st==0,'legacy inlet')
 legacy_in_error=abs(pressure_fun_uj(bc(1),bc(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)-pf(1))
 call boundary_flux(u(n+2,:),th(n),cfg,.false.,f,st,bc);call check(st==0,'legacy outlet')
 legacy_out_error=maxval(abs(f(:2)-[cfg%mul,cfg%mug]))
 cfg%boundary_method='STEADY_TRACE'
 call boundary_flux(u(3,:),th(1),cfg,.true.,f,st,bc);call check(st==0,'source-compatible inlet')
 trace_in_error=abs(pressure_fun_uj(bc(1),bc(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)-pf(1))
 call check(trace_in_error<5e-9_wp,'face pressure agrees with independently integrated stationary face')
 call check(trace_in_error<legacy_in_error/100,'center/face source error removed')
 call check(f(1)==cfg%mul.and.f(2)==cfg%mug,'same prescribed inlet mass components')
 call boundary_flux(u(n+2,:),th(n),cfg,.false.,f,st,bc);call check(st==0,'source-compatible outlet')
 trace_out_error=maxval(abs(f(:2)-[cfg%mul,cfg%mug]))
 call check(trace_out_error<5e-9_wp,'stationary phase outflows recovered at physical face')
 call check(trace_out_error<legacy_out_error/100,'stationary phase split improved')
 call check(abs(pressure_fun_uj(bc(1),bc(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)-cfg%nPs)<1e-11_wp,'absolute pressure imposed at face')
 do side=1,2
  i=3;if(side==2)i=n+2
  do k=1,3
   do ih=1,4
   h=10._wp**(-ih-2)*max(abs(u(i,k)),1e-4_wp);up=u(i,:);um=up;up(k)=up(k)+h;um(k)=um(k)-h
   call boundary_flux(up,cfg%thetav(i),cfg,side==1,fp,st);call check(st==0,'FD full plus')
   call boundary_flux(um,cfg%thetav(i),cfg,side==1,fm,st);call check(st==0,'FD full minus')
   j2=(fp-fm)/(2*h);up=u(i,:);um=up;up(k)=up(k)+h/2;um(k)=um(k)-h/2
   call boundary_flux(up,cfg%thetav(i),cfg,side==1,fhp,st);call check(st==0,'FD half plus')
   call boundary_flux(um,cfg%thetav(i),cfg,side==1,fhm,st);call check(st==0,'FD half minus')
   j4=(4*(fhp-fhm)/h-j2)/3
   write(*,'(a,2i3,4es20.10)')'FD side,column,h,error,j2(3),j4(3) ',side,k,h,maxval(abs(j4-j2)/max(1._wp,abs(j4))),j2(3),j4(3)
   if(ih>=3)call check(maxval(abs(j4-j2)/max(1._wp,abs(j4)))<1e-5_wp,'source-compatible complete derivative on resolved h stencil')
   enddo
  enddo
 enddo
 allocate(res(3,n),rplus(3,n),rminus(3,n),az(3,3,n),bz(3,3,n),cz(3,3,n),daz(3,3,n),dbz(3,3,n),dcz(3,3,n),direction(n+4,3),work(n+4,3))
 call sim%init(cfg);call sim%prepare_dissipation(u,st);call check(st==0,'full residual frozen dissipation')
 call sim%assemble(u,u,.001_wp,res,az,bz,cz,resnorm,fin,fout,st,.true.)
 call check(st==0,'full source-compatible residual and Jacobian')
 direction=0;direction(3,:)=u(3,:)*[.37_wp,-.23_wp,.11_wp];direction(n+2,:)=u(n+2,:)*[-.31_wp,.19_wp,.43_wp]
 h=1e-5_wp;work=u+h*direction
 call sim%assemble(work,u,.001_wp,rplus,daz,dbz,dcz,nt,fin,fout,st,.false.)
 call check(st==0,'full residual plus, boundary recomputed')
 work=u-h*direction
 call sim%assemble(work,u,.001_wp,rminus,daz,dbz,dcz,nt,fin,fout,st,.false.)
 call check(st==0,'full residual minus, boundary recomputed')
 jac_error=0
 do jj=1,n
  prediction=matmul(bz(:,:,jj),direction(jj+2,:))
  if(jj>1)prediction=prediction+matmul(az(:,:,jj),direction(jj+1,:))
  if(jj<n)prediction=prediction+matmul(cz(:,:,jj),direction(jj+3,:))
  observed=(rplus(:,jj)-rminus(:,jj))/(2*h)
  jac_error=max(jac_error,maxval(abs(observed-prediction)/max(1._wp,abs(prediction))))
 enddo
 call check(jac_error<1e-5_wp,'assembled complete Jacobian agrees with independent boundary direction')
 write(*,'(a,es20.10)')'FULL_SOURCE_JAC_ERROR ',jac_error
 rln=1+(cfg%nPs-1)/cfg%nCl**2;rgn=cfg%nPs/cfg%nCg**2
 up=[.75_wp*rln,.25_wp*rgn,-.6_wp];cfg%closure_riser=2
 call check_boundary_characteristics(up,0._wp,cfg,.false.,st)
 call check(st==BC_DATA_REQUIRED,'reversed incoming material mode needs reservoir data')
 cfg%boundary_method='LEGACY'
 call boundary_flux(up,0._wp,cfg,.false.,f,st)
 call check(st==BC_DATA_REQUIRED,'actual outlet does not certify backflow without reservoir data')
 write(*,'(a,4es20.10)')'SOURCE_BOUNDARY inlet old/new, phase out old/new ',legacy_in_error,trace_in_error,legacy_out_error,trace_out_error
 write(*,'(a,i0,a)')'PASS source boundary: ',checks,' assertions'
contains
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label;error stop 1
  endif
  checks=checks+1
 end subroutine
end program
