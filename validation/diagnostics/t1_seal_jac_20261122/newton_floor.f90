program newton_floor
  !! Round 20261122 (P0.3): offline Newton on the GLOBAL residual of one time step from an
  !! accepted state (fixture built from checkpoint.dat by ckpt_fixture.py: Uold = trial =
  !! accepted U, dt = next attempted step). Full Newton steps x <- x - J^{-1} R with the
  !! solver's own assemble (junction_seal_curve_steps=16, IMPLICIT face Jacobian, face tol
  !! 1e-11), dense LU (dgesv). Per iteration: solver norm, and the five largest scaled
  !! entries |R_k|/scale_k with cell, component and phase fractions: locates the floor.
  !! Argument 1: fixture folder. Env as seal_jac_check (closure 5, T1, LIQUID_BRIDGE).
  use offline_api
  implicit none
  character(512)::folder
  real(wp)::x(nx),sc(nx),meta(6),res(nx),norm,fin(3),fout(3),dx(nx),q(nx),uu(nc,3)
  real(wp),allocatable::jac(:,:)
  integer(c_int)::st,idx,stg
  integer::it,k,m,info,ipiv(nx),cell,comp,top(5),code
  character(16)::envs
  allocate(jac(nx,nx))
  call get_command_argument(1,folder)
  call probe_init(trim(folder),len_trim(folder),1_c_int,st)
  if(st/=0)then
    write(*,'(a,i0)')'probe_init status ',st;stop 1
  endif
  cfg%junction_seal_curve_steps=16;cfg%junction_seal_jacobian='IMPLICIT';cfg%junction_seal_face_tol=1e-11_wp
  call get_environment_variable('PRIMS_SEAL_STEPS',envs,status=code);if(code==0)read(envs,*)cfg%junction_seal_curve_steps
  call sim%init(cfg);call sim%prepare_dissipation(old,code)
  if(code/=0)stop 'prepare_dissipation'
  call probe_state(0_c_int,x,sc,meta)
  write(*,'(a,f12.6,a,es12.4)')'STATE|t_s ',meta(1),'|dt_s ',meta(2)
  do it=0,10
    call probe_eval(x,1_c_int,res,jac,norm,fin,fout,st,idx,stg)
    if(st/=0)then
      write(*,'(a,i0,a,i0,a,i0)')'ITER ',it,'|assemble status ',st,'|index ',idx;exit
    endif
    q=abs(res)/sc
    do m=1,5
      top(m)=maxloc(q,1);q(top(m))=-1
    enddo
    q=abs(res)/sc
    write(*,'(a,i0,a,es12.4,a,es12.4)')'ITER ',it,'|solver norm ',norm,'|max|R|/scale ',maxval(q)
    call unpack(x,uu)
    do m=1,5
      k=top(m);cell=(k-1)/3+lo;comp=mod(k-1,3)+1
      write(*,'(a,i4,a,i1,a,es11.3,a,es11.3,a,es11.3,a,3es13.5)')'   cell ',cell,'|comp ',comp,'|scaled ',q(k), &
        '|R ',res(k),'|scale ',sc(k),'|U ',uu(cell,:)
    enddo
    if(norm<1e-13_wp)exit
    dx=-res
    call dgesv(nx,1,jac,nx,ipiv,dx,nx,info)
    if(info/=0)then
      write(*,'(a,i0)')'dgesv info ',info;exit
    endif
    write(*,'(a,es12.4)')'   max|dx|/scale ',maxval(abs(dx)/sc)
    x=x+dx
  enddo
end program newton_floor
