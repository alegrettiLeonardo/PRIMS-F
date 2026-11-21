program t1_face_check
  !! P0-B/C: junction face of the GLOBAL operator with junction_transmission
  !! T2 and T1 on the SAME accepted conservative state (offline_api fixture,
  !! riser closure from PRIMS_PROBE_CLOSURE_RISER), lagged dissipation built at
  !! that state as in the solver. Reports (SI):
  !!   FACE|tr|status|Gl* Gg* F3_L* F3_R* p_L* p_R* alpha_pipe* alpha_riser*
  !!   JUMP|tr|fm-fp (3 components; T1 must be 0 in all three)
  !!   BAL |tr|k| sum_i w_i R_i,k - [sum_i w_i (U-Uold)_i,k - dt sum_i w_i S_i,k + dt (Fout-Fin)_k]
  !!        which equals dt * sum_faces (fm - fp)_k (the singular source of the scheme):
  !!        0 for T1 in every component; dt*(F3_L*-F3_R*) for the T2 momentum.
  !! Usage: t1_face_check FIXTURE C_GAS
  use offline_api
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  character(512)::folder,arg
  real(wp)::x(nx),sc(nx),meta(6),uu(nc,3),uo(nc,3),fm(3),fp(3),tr6(6),r(3,np),a(3,3,np),b(3,3,np),c(3,3,np)
  real(wp)::norm,fin(3),fout(3),lhs(3),acc(3),src(3),v(6),f(3),s(3),aj(3,3),sj(3,3),dtn,wscale(3),fscale(3)
  integer(c_int)::st
  integer::k,i,code,cid
  character(2)::trs(2)=['T2','T1']
  call get_command_argument(1,folder)
  call probe_init(trim(folder),len_trim(folder),1_c_int,st)
  if(st/=0)then
    write(*,'(a,i0)')'probe_init status ',st;stop 1
  endif
  call probe_state(0_c_int,x,sc,meta);call unpack(x,uu)
  uo=uu;uo(:,1)=uu(:,1)*(1-1e-3_wp);uo(:,2)=uu(:,2)*(1+1e-3_wp)
  dtn=0.05_wp/cfg%omega_t
  fscale=[cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u,cfg%omega_rho*cfg%omega_u**2]
  do k=1,2
    cfg%junction_transmission=trs(k);call sim%init(cfg)
    call sim%prepare_dissipation(uu,code)
    if(code/=0)then
      write(*,'(3a,i0)')'FACE|',trs(k),'|prepare_dissipation ',code;cycle
    endif
    call sim%face_pair(uu(ib,:),uu(ib+1,:),cfg%thetav(ib),cfg%thetav(ib+1),fm,fp,code,junction_trace=tr6, &
      dissipation_matrix=sim%lagged_dissipation(:,:,ib),left_cell=ib)
    if(trs(k)=='T2')then
      write(*,'(3a,i6,a,8es15.7)')'FACE|',trs(k),'|',code,'|',fm(1:2)*fscale(1:2),fm(3)*fscale(3),fp(3)*fscale(3), &
        tr6(1)*cfg%omega_P,tr6(1)*cfg%omega_P,tr6(4:5)
    else
      write(*,'(3a,i6,a,8es15.7)')'FACE|',trs(k),'|',code,'|',fm(1:2)*fscale(1:2),fm(3)*fscale(3),fp(3)*fscale(3), &
        tr6(1)*cfg%omega_P,tr6(6)*cfg%omega_P,tr6(4:5)
    endif
    write(*,'(3a,3es12.4)')'JUMP|',trs(k),'|',(fm-fp)*fscale
    call sim%assemble(uu,uo,dtn,r,a,b,c,norm,fin,fout,code,.false.)
    if(code/=0)then
      write(*,'(3a,i0)')'BAL |',trs(k),'|assemble ',code;cycle
    endif
    lhs=0;acc=0;src=0;wscale=0
    do i=lo,hi
      cid=cfg%closure_pipe;if(i>ib)cid=cfg%closure_riser
      call sim%point_state(uu(i,:),cfg%thetav(i),v,f,s,aj,sj,code,.false.,cid)
      lhs=lhs+cfg%cell_width(i)*r(:,i-lo+1)
      acc=acc+cfg%cell_width(i)*(uu(i,:)-uo(i,:))
      src=src+cfg%cell_width(i)*s
      wscale=wscale+cfg%cell_width(i)*abs(uu(i,:)-uo(i,:))+dtn*cfg%cell_width(i)*abs(s)
    enddo
    write(*,'(3a,3es13.4,a,3es13.4,a,3es13.4)')'BAL |',trs(k),'|identity residual ',lhs-(acc-dtn*src+dtn*(fout-fin)), &
      ' | dt*(fm-fp) ',dtn*(fm-fp),' | scale ',wscale+dtn*(abs(fout)+abs(fin))
  enddo
end program t1_face_check
