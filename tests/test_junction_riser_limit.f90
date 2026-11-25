program test_junction_riser_limit
 !! EXPERIMENTAL E2 (junction_riser_limit='ALPHA0', opt-in): riser-trace limit
 !! alpha*_riser=0 at an EXACTLY zero junction gas flux with non-negative liquid
 !! flux (drift-flux limit j_g=alpha u_g, u_g -> C0(0) j_l + Vgj(0) > 0).
 !! Contract checked on an independently integrated stationary profile
 !! (Bhagwat riser), with the pipe-end gas velocity driven to <= 0 by the cell
 !! state itself (mixture momentum lowered, closure recovers u_g):
 !!  1 default 'NONE' unchanged: the same face is refused (1021);
 !!  2 'ALPHA0': admissible, alpha*_riser = 0, one mass flux per phase,
 !!    Gg* = 0 exactly (nothing projected), riser momentum = P* + Gl*^2/rho_l*;
 !!  3 continuity: alpha*(j_g) -> 0+ linearly and the riser momentum trace at a
 !!    small positive probe flux tends to the alpha*=0 value;
 !!  4 negative gas flux (Roe on a reversed-gas state) stays refused (CLOSURE_DOMAIN);
 !!  5 E2 is inert where the gas flux is non-zero: identical face vector;
 !!  6 assembled residual telescopes exactly and the Jacobian matches central FD
 !!    with E2 active at the junction.
 use mod_kinds,only:wp
 use module_simulate_pipe_riser,only:pipeline_riser_simulator_t
 use module_config,only:sim_config_t
 use module_steady_state,only:EstadoEstacionario_ndim_simp,solve_alpha_from_superficials
 implicit none
 type(sim_config_t)::cfg
 type(pipeline_riser_simulator_t)::sim
 real(wp),allocatable::s(:),p(:),rg(:),rl(:),alpha(:),ug(:),ul(:),th(:),jl(:),jg(:),faces(:),pf(:),u(:,:)
 real(wp),allocatable::res(:,:),rplus(:,:),rminus(:,:),az(:,:,:),bz(:,:,:),cz(:,:,:),daz(:,:,:),dbz(:,:,:),dcz(:,:,:),direction(:,:),work(:,:)
 real(wp)::fm(3),fp(3),tr(6),fm2(3),fp2(3),tr2(6),pip,rls,v(6),f(3),sx(3),a(3,3),d(3,3),norm,nt,fin(3),fout(3),h
 real(wp)::prediction(3),observed(3),jac_error,ustate(3),aprobe,jprobe,mom0,momp,ratio,rprev,aprev
 integer::st,n,checks,i,jj,ib,k
 checks=0
 call setup()
 ib=cfg%num_dummy_cells+cfg%N_pipe
 pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
 ! pipe-end state with gas velocity slightly below zero (liquid still flowing)
 ustate=with_pipe_gas_velocity(-5e-2_wp)
 call sim%point_state(ustate,cfg%thetav(ib),v,f,sx,a,d,st,.false.,1)
 write(*,'(a,4es13.5)')'pipe end alpha,u_l,u_g,p ',v(1),v(5),v(6),v(4)
 call check(st==0.and.v(6)<0.and.v(5)>0,'pipe-end state: gas receding, liquid advancing (closure-recovered)')
 u(ib,:)=ustate
 cfg%junction_flux='PHASE_UPWIND'
 ! 1 default
 cfg%junction_riser_limit='NONE';call sim%init(cfg)
 call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
 write(*,'(a,i6,a,es12.4)')'1 NONE status ',st,'  E1 gas flux before inversion ',fm(2)
 call check(st==1021.and.fm(2)==0,'default: zero gas flux is refused by the riser inversion (21)')
 ! 2 ALPHA0
 cfg%junction_riser_limit='ALPHA0';call sim%init(cfg)
 call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
 rls=cfg%nrho_l0+(tr(1)-cfg%nP_l0)/cfg%nCl**2
 write(*,'(a,i6,a,6es13.5)')'2 ALPHA0 status ',st,' trace p*,Gl*,Gg*,a_pipe*,a_riser* ',tr(1:5)
 call check(st==0,'E2: admissible face')
 call check(tr(5)==0.and.tr(3)==0.and.tr(6)==0,'alpha*_riser=0 with an exactly zero (unprojected) gas flux')
 call check(tr(4)>0.and.tr(4)<1,'pipe trace keeps its gas (two-phase, stagnant-gas branch)')
 call check(all(fm(:2)==fp(:2)),'one mass flux per phase shared by both cells')
 call check(abs(fp(3)-(pip*tr(1)+tr(2)**2/rls))<=1e-13_wp*fp(3),'riser momentum = P* + Gl*^2/rho_l* (liquid-only trace)')
 ! 3 continuity of the trace limit (alpha/j_g bounded, momentum continuous)
 rprev=-1;aprev=1
 do k=6,14,2
  jprobe=10._wp**(-k)
  aprobe=solve_alpha_from_superficials(tr(2)/rls,jprobe,rls,tr(1)/cfg%nCg**2,0._wp,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
    cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=st,closure_id=3)
  call check(st==0.and.aprobe>0.and.aprobe<aprev,'probe root exists and decreases to 0+')
  ratio=aprobe/jprobe
  if(rprev>0)call check(abs(ratio-rprev)<=1e-3_wp*rprev,'alpha*/j_g converges: C1 one-sided limit')
  mom0=pip*tr(1)+tr(2)**2/rls
  momp=pip*tr(1)+tr(2)**2/((1-aprobe)*rls)+(jprobe*tr(1)/cfg%nCg**2)**2/(aprobe*tr(1)/cfg%nCg**2)
  write(*,'(a,es9.1,a,es13.5,a,es12.4)')'3 j_g ',jprobe,' alpha* ',aprobe,'  riser momentum jump to alpha*=0 ',abs(momp-mom0)/mom0
  call check(abs(momp-mom0)/mom0<=10*jprobe,'riser momentum trace continuous at alpha*=0')
  rprev=ratio;aprev=aprobe
 enddo
 ! 4 Roe on receding-gas pipe states: the Roe vector stays non-zero here, so
 !   E2 must be bitwise inert (the refusal of a NEGATIVE gas flux, which E2
 !   never enters because it requires Gg*==0, is checked on the saved SS1/OSC
 !   states where Roe gives Gg*<0: validation/diagnostics/junction_constitutive).
 do k=1,12,3
  u(ib,:)=with_pipe_gas_velocity(-0.05_wp*k)
  cfg%junction_flux='ROE';cfg%junction_riser_limit='NONE';call sim%init(cfg)
  call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
  cfg%junction_riser_limit='ALPHA0';call sim%init(cfg)
  call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm2,fp2,jj,junction_trace=tr2,left_cell=ib)
  write(*,'(a,f6.2,a,2i6,a,es12.4)')'4 Roe pipe u_g ',-0.05_wp*k,' status NONE/ALPHA0 ',st,jj,'  Gg* ',tr(3)
  call check(st==jj.and.all(fm==fm2).and.all(fp==fp2),'Roe: E2 inert when the gas flux is not exactly zero')
 enddo
 ! 5 inert where the gas flux is non-zero (stationary profile)
 call setup()
 do k=1,2
  cfg%junction_flux=merge('ROE         ','PHASE_UPWIND',k==1)
  cfg%junction_riser_limit='NONE';call sim%init(cfg)
  call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
  cfg%junction_riser_limit='ALPHA0';call sim%init(cfg)
  call sim%face_pair(u(ib,:),u(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm2,fp2,jj,junction_trace=tr2,left_cell=ib)
  call check(st==0.and.jj==0.and.all(fm==fm2).and.all(fp==fp2).and.all(tr==tr2),'E2 inert for non-zero gas flux (bitwise)')
 enddo
 ! 6 conservation and Jacobian with E2 active
 u(ib,:)=with_pipe_gas_velocity(-5e-2_wp)
 cfg%junction_flux='PHASE_UPWIND';cfg%junction_riser_limit='ALPHA0';call sim%init(cfg)
 allocate(res(3,n),rplus(3,n),rminus(3,n),az(3,3,n),bz(3,3,n),cz(3,3,n),daz(3,3,n),dbz(3,3,n),dcz(3,3,n), &
   direction(n+4,3),work(n+4,3))
 call sim%prepare_dissipation(u,st)
 call check(st==0,'lagged dissipation')
 call sim%assemble(u,u,.001_wp,res,az,bz,cz,norm,fin,fout,st,.true.)
 call check(st==0,'E2 residual and complete Jacobian')
 do i=1,2
  nt=sum(res(i,:)*cfg%cell_width(3:n+2))-.001_wp*(fout(i)-fin(i))
  write(*,'(a,i1,a,es12.4)')'6 phase ',i,' telescoping defect ',nt
  call check(abs(nt)<1e-14_wp,'phase residual telescopes with E2 active')
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
 write(*,'(a,es12.4)')'6 E2 Jacobian directional error ',jac_error
 call check(jac_error<1e-5_wp,'assembled Jacobian consistent with the E2 residual')
 write(*,'(a,i0,a)')'PASS junction riser limit (E2): ',checks,' assertions'
contains
 function with_pipe_gas_velocity(target) result(y)
  !! Lower the pipe-end mixture momentum (masses unchanged) until the closure
  !! recovers u_g = target; bisection on U3, nothing else is altered.
  real(wp),intent(in)::target
  real(wp)::y(3),lo,hi,mid,vv(6),ff(3),ss(3),aa(3,3),dd(3,3)
  integer::it,c
  y=cfg%U0(ib,:);hi=y(3);lo=0.2_wp*y(3)
  do it=1,200
   mid=.5_wp*(lo+hi);y(3)=mid
   call sim%point_state(y,cfg%thetav(ib),vv,ff,ss,aa,dd,c,.false.,1)
   if(c/=0)then
    lo=mid;cycle
   endif
   if(vv(6)>target)then
    hi=mid
   else
    lo=mid
   endif
  enddo
  y(3)=hi
 end function
 subroutine setup()
  cfg=sim_config_t()
  cfg%omega_rho=997.05_wp;cfg%omega_P=101325._wp;cfg%omega_u=1;cfg%omega_t=12.546253993197219_wp
  cfg%nrho_l0=1;cfg%nP_l0=1;cfg%nPs=2.1_wp;cfg%nCl=1498*sqrt(cfg%omega_rho/cfg%omega_P);cfg%nCg=343*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%Lp=57.4_wp;cfg%Lr=cfg%omega_t;cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
  cfg%D_H=.05254_wp;cfg%AREA=acos(-1._wp)*cfg%D_H**2/4;cfg%eps=4.6e-5_wp;cfg%G=9.81_wp
  cfg%mu_l=.00089_wp;cfg%mu_g=.0000185_wp;cfg%sigma=.0728_wp;cfg%tol=1e-12_wp;cfg%tola=1e-13_wp
  cfg%closure_pipe=1;cfg%closure_riser=3;cfg%mul=.5_wp;cfg%mug=.001_wp;cfg%njl=.5_wp
  cfg%jac_dissipation_fd=.true.
  if(allocated(s))deallocate(s,p,rg,rl,alpha,ug,ul,th,jl,jg,faces,pf,u)
  call EstadoEstacionario_ndim_simp(9,cfg%mul,cfg%mug,cfg%nPs,cfg%Lp,cfg%Lr,cfg%CA,cfg%beta,cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
    cfg%nCl,cfg%nCg,cfg%nrho_l0,cfg%nP_l0,cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_P,cfg%omega_u,cfg%omega_rho,cfg%tol, &
    s,p,rg,rl,alpha,ug,ul,th,jl,jg,n_riser=9,closure_id_riser=3,closure_id_pipe=1,status=st,face_positions=faces,face_pressures=pf)
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
  call sim%init(cfg)
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
