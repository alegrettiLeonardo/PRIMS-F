module audit_observer
 use mod_kinds,only:wp
 use module_config,only:sim_config_t,first_physical_cell,last_physical_cell
 use module_boundary_conditions,only:boundary_flux
 use module_find,only:pressure_fun_uj
 use iso_c_binding,only:c_int,c_char,c_null_char
 implicit none
 integer,save::iu=0,ju=0
 interface
  function atomic_rename(old,new)bind(C,name='rename')result(st)
   import c_int,c_char
   character(c_char),intent(in)::old(*),new(*)
   integer(c_int)::st
  end function
 end interface
contains
 subroutine audit_record(cfg,u,step,t,dt,fin,fout)
  type(sim_config_t),intent(in)::cfg
  real(wp),intent(in)::u(:,:),t,dt,fin(3),fout(3)
  integer,intent(in)::step
  integer::l,h,stl,str,ck,i
  real(wp)::fb(3),lb(3),rb(3),ff(3),pl,pr,rl,rg,sc,mismatch
  l=first_physical_cell(cfg);h=last_physical_cell(cfg)
  if(iu==0)then
   open(newunit=iu,file=trim(cfg%out_dir)//'/boundary_contract.csv',status='new')
   write(iu,'(a)')'step,time_s,dt_s,inlet_status,outlet_status,mdot_l_eff_kg_s,mdot_g_eff_kg_s,'// &
     'mdot_l_Qtarget_kg_s,mdot_l_mass_target_kg_s,mdot_g_target_kg_s,pB_Pa,Pt_Pa,p_in_face_Pa,'// &
     'jl_eff_m_s,jl_target_m_s,rhol_in_kg_m3,rhog_in_kg_m3,rm_l_kg_s,rm_g_kg_s,rp_Pa,'// &
     'flux_recompute_max_nd'
   open(newunit=ju,file=trim(cfg%out_dir)//'/boundary_states.csv',status='new')
   write(ju,'(a)')'step,time_s,theta_in_rad,theta_out_rad,UiL1,UiL2,UiL3,UiR1,UiR2,UiR3,'// &
     'UBL1,UBL2,UBL3,UBR1,UBR2,UBR3'
  endif
  call boundary_flux(u(l,:),cfg%thetav(l),cfg,.true.,fb,stl,lb)
  call boundary_flux(u(h,:),cfg%thetav(h),cfg,.false.,ff,str,rb)
  if(stl/=0.or.str/=0)error stop 'observer: accepted boundary reconstruction failed'
  mismatch=max(maxval(abs(fb-fin)),maxval(abs(ff-fout)))
  pl=pressure_fun_uj(lb(1),lb(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  pr=pressure_fun_uj(rb(1),rb(2),cfg%nCg,cfg%nCl,cfg%nrho_l0,cfg%nP_l0)
  rl=cfg%nrho_l0+(pl-cfg%nP_l0)/cfg%nCl**2;rg=pl/cfg%nCg**2
  sc=cfg%AREA*cfg%omega_rho*cfg%omega_u
  write(iu,'(i0,2(",",es25.16e3),2(",",i0),16(",",es25.16e3))')step,t,dt,stl,str, &
   sc*fin(1:2),sc*rl*cfg%njl,sc*cfg%mul,sc*cfg%mug,pr*cfg%omega_P,cfg%nPs*cfg%omega_P, &
   pl*cfg%omega_P,fin(1)/rl*cfg%omega_u,cfg%njl*cfg%omega_u,rl*cfg%omega_rho,rg*cfg%omega_rho, &
   sc*(fin(1)-rl*cfg%njl),sc*(fin(2)-cfg%mug),(pr-cfg%nPs)*cfg%omega_P,mismatch
  flush(iu)
  write(ju,'(i0,15(",",es25.16e3))')step,t,cfg%thetav(l),cfg%thetav(h),u(l,:),u(h,:),lb,rb
  flush(ju)
  ! Full conservative checkpoint, binary64 round-trip decimal; C/POSIX rename
  ! atomically replaces the previous accepted snapshot, after closing .next.
  open(newunit=ck,file=trim(cfg%out_dir)//'/audit_checkpoint.next',status='replace')
  write(ck,'(i0,2(1x,es25.16e3))')step,t,dt
  do i=1,size(u,1)
   write(ck,'(i0,5(1x,es25.16e3))')i,cfg%thetav(i),cfg%cell_width(i),u(i,:)
  enddo
  close(ck)
  stl=atomic_rename(trim(cfg%out_dir)//'/audit_checkpoint.next'//c_null_char, &
    trim(cfg%out_dir)//'/audit_checkpoint.dat'//c_null_char)
  if(stl/=0)error stop 'observer: checkpoint rename failed'
 end subroutine
end module
