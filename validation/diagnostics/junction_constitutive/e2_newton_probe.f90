program e2_newton_probe
 !! Gate P0.3 diagnostic (not used by the solver): reproduce, offline, the
 !! Newton direction of the last refused step of an E2 replay and walk along it.
 !! For the saved iterate (state 1 = trial) and the saved old state: assemble
 !! residual and FULL Jacobian with the run operator (E1 + ALPHA0), solve
 !! J d = -R (LAPACK dgesv), then for lambda = 1, 1e-1, ..., 1e-8 evaluate the
 !! residual at U + lambda d and report the junction quantities: E1 internals
 !! (pipe/riser trace velocities, acoustic shift dj), the transmitted Gl*, Gg*
 !! and the assembly status. Also reports the size of d in the junction cells.
 !! Arguments: fixture folder, c_gas m/s, junction_flux, junction_trace, riser limit.
 use offline_api
 use module_simulate_pipe_riser,only:junction_phase_upwind
 use module_steady_state,only:solve_alpha_from_superficials
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),res(nx),jac(nx,nx),norm,fin(3),fout(3),d(nx),cgas,lam,xt(nx),uu(nc,3)
 real(wp)::fm(3),fp(3),tr(6),va(6),vb(6),f(3),s(3),a(3,3),dd(3,3),gl,gg,ps,ja,jb,pip
 integer(c_int)::st,idx,stage
 integer::ipiv(nx),info,k,code,i,nr
 real(wp)::fr(3),astl,mapped(3),vm(6),za,zb,dj,compl
 character(512)::folder,arg
 character(16)::runflux,runtrace,runlimit
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call get_command_argument(3,runflux);call get_command_argument(4,runtrace);call get_command_argument(5,runlimit)
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P);cfg%tol_res=1e-10_wp;cfg%tol_inc=1e-10_wp
 cfg%junction_flux=runflux;cfg%junction_trace=runtrace;cfg%junction_riser_limit=runlimit
 call sim%init(cfg);call sim%prepare_dissipation(old,code)
 write(*,'(a,i6)')'PREPARE|',code
 pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
 call probe_state(1_c_int,x,sc,meta)
 call probe_eval(x,1_c_int,res,jac,norm,fin,fout,st,idx,stage)
 write(*,'(a,i6,a,es12.4,a,es12.4)')'ITERATE|status',st,'|norm',norm,'|dt_s',meta(2)
 if(st/=0)stop 3
 d=-res
 call dgesv(nx,1,jac,nx,ipiv,d,nx,info)
 write(*,'(a,i4,a,es12.4)')'SOLVE|info',info,'|max|d|',maxval(abs(d))
 do i=ib-1,ib+2
  k=3*(i-lo)
  write(*,'(a,i4,a,3es13.5,a,3es13.5)')'DIR|cell',i,'|U',x(k+1:k+3),'|d',d(k+1:k+3)
 enddo
 ! Replicate the MATCHED mapping of the riser trace to the pipe closure (CENTER)
 call unpack(x,uu)
 call sim%point_state(uu(ib,:),cfg%thetav(ib),va,f,s,a,dd,code,.false.,cfg%closure_pipe)
 call sim%point_state(uu(ib+1,:),cfg%thetav(ib+1),vb,fr,s,a,dd,code,.false.,cfg%closure_riser)
 astl=solve_alpha_from_superficials(fr(1)/vb(2),fr(2)/vb(3),vb(2),vb(3),cfg%thetav(ib),cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
   cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=code,closure_id=cfg%closure_pipe,alpha_hint=va(1))
 mapped=[(1-astl)*vb(2),astl*vb(3),fr(1)+fr(2)]
 call sim%point_state(mapped,cfg%thetav(ib),vm,f,s,a,dd,code,.false.,cfg%closure_pipe)
 compl=va(1)/(va(3)*pip*cfg%nCg**2)+(1-va(1))/(va(2)*pip*cfg%nCl**2);za=sqrt((uu(ib,1)+uu(ib,2))/compl)
 compl=vm(1)/(vm(3)*pip*cfg%nCg**2)+(1-vm(1))/(vm(2)*pip*cfg%nCl**2);zb=sqrt((mapped(1)+mapped(2))/compl)
 dj=-pip*(vm(4)-va(4))/(za+zb)
 write(*,'(a,es14.6,a,es14.6,a,es14.6,a,es14.6,a,es14.6)')'MAPPED|alpha_map',astl,'|ug_map m/s',vm(6)*cfg%omega_u, &
   '|dj m/s',dj*cfg%omega_u,'|ug_map+dj',(vm(6)+dj)*cfg%omega_u,'|ug_pipe+dj',(va(6)+dj)*cfg%omega_u
 call solve_count(fr(1)/vb(2),fr(2)/vb(3),vb(2),vb(3))
 do k=0,8
  lam=10._wp**(-k)
  xt=x+lam*d
  call probe_eval(xt,0_c_int,res,jac,norm,fin,fout,st,idx,stage)
  call unpack(xt,uu)
  call sim%point_state(uu(ib,:),cfg%thetav(ib),va,f,s,a,dd,code,.false.,cfg%closure_pipe)
  call sim%point_state(uu(ib+1,:),cfg%thetav(ib+1),vb,f,s,a,dd,i,.false.,cfg%closure_riser)
  gl=-999;gg=-999;ps=-999;ja=-999;jb=-999
  if(code==0.and.i==0)then
   call junction_phase_upwind(cfg,va,uu(ib,:),vb,uu(ib+1,:),gl,gg,ps,i)
   ja=(1-va(1))*va(5)+va(1)*va(6);jb=(1-vb(1))*vb(5)+vb(1)*vb(6)
  endif
  call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,i,junction_trace=tr,left_cell=ib)
  write(*,'(a,es9.1,a,i6,a,i4,a,i2,a,es11.3,a,i6,a,4es12.4,a,2es12.4,a,es12.4)')'WALK|lambda',lam,'|assembly',st,'|index',idx, &
    '|stage',stage,'|norm',norm,'|face',i,'|ug_pipe,ug_riser,alpha_riser,dP_Pa',va(6)*cfg%omega_u,vb(6)*cfg%omega_u,vb(1), &
    pip*(vb(4)-va(4))*cfg%omega_rho*cfg%omega_u**2,'|E1 Gl,Gg',gl*cfg%omega_rho*cfg%omega_u,gg*cfg%omega_rho*cfg%omega_u, &
    '|ul_riser',vb(5)*cfg%omega_u
 enddo
contains
 subroutine solve_count(jl_,jg_,rl_,rg_)
  real(wp),intent(in)::jl_,jg_,rl_,rg_
  real(wp)::aa
  integer::c
  aa=solve_alpha_from_superficials(jl_,jg_,rl_,rg_,cfg%thetav(ib),cfg%D_H,cfg%AREA,cfg%eps,cfg%G, &
    cfg%mu_l,cfg%mu_g,cfg%sigma,cfg%omega_u,cfg%omega_rho,cfg%tol,status=c,closure_id=cfg%closure_pipe,nroots=nr)
  write(*,'(a,i4,a,i3,a,es12.4)')'MAPPING_ROOTS|status',c,'|nroots',nr,'|j_g m/s',jg_*cfg%omega_u
 end subroutine
end program
