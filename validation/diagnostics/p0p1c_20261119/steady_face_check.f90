program steady_face_check
 !! Diagnostic only: junction face with junction_trace='CENTER' and 'STEADY'
 !! at a saved state (old/trial), against the physical gas fluxes of the two
 !! adjacent cells. Geometry of the 51+51 benchmark mesh is rebuilt as in
 !! replay_p0 (STEADY_TRACE branch). Arguments: folder, gas sound speed m/s.
 use offline_api
 implicit none
 real(wp)::x(nx),sc(nx),meta(6),cgas,s,fm(3),fp(3),tr(6),uu(nc,3),v(6),f(3),src(3),a(3,3),ds(3,3),gl,gr
 integer(c_int)::st
 integer::i,w,m,cc
 character(512)::folder,arg
 character(6)::mode(2)=['CENTER','STEADY']
 call get_command_argument(1,folder);call get_command_argument(2,arg);read(arg,*)cgas
 call probe_init(trim(folder),len_trim(folder),1_c_int,st)
 cfg%nCg=cgas*sqrt(cfg%omega_rho/cfg%omega_P)
 cfg%CA=3.01818193717446_wp;cfg%beta=2*acos(-1._wp)/180
 allocate(cfg%vns(nc),cfg%face_positions(0:np));cfg%face_positions(0)=0;cfg%vns=0
 do i=1,np
  cfg%face_positions(i)=cfg%face_positions(i-1)+cfg%cell_width(i+2)
  cfg%vns(i+2)=.5_wp*(cfg%face_positions(i-1)+cfg%face_positions(i))
 enddo
 s=cfg%omega_rho*cfg%omega_u
 do w=0,1
  call probe_state(int(w,c_int),x,sc,meta);call unpack(x,uu)
  call sim%point_state(uu(ib,:),cfg%thetav(ib),v,f,src,a,ds,cc,.false.,1);gl=f(2)
  call sim%point_state(uu(ib+1,:),cfg%thetav(ib+1),v,f,src,a,ds,cc,.false.,3);gr=f(2)
  do m=1,2
   cfg%junction_trace=mode(m);call sim%init(cfg)
   call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,st,junction_trace=tr,left_cell=ib)
   print '(a,i2,1x,a,a,i6,a,2es13.5,a,es13.5,a,es14.6)','STATE',w,mode(m),' status',st,' Gg cells L,R',gl*s,gr*s, &
     ' Gg face ',tr(3)*s,' p* Pa ',tr(1)*cfg%omega_P
  enddo
 enddo
end program
