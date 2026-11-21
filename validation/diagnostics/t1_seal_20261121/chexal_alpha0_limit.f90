program chexal_alpha0_limit
  !! P0.2: CHEXAL_1992 (air-water) as alpha -> 0+ at the SS1 riser-base state
  !! (theta = 2.334 deg, jl = 0.5264 m/s, rho_l = 997.137, rho_g = 2.524 kg/m3,
  !! D = 0.05254 m). For each alpha: the drift relation u_g = C0 (jl+jg) + Vgj with
  !! jg = alpha u_g is solved by fixed point; prints C0, Vgj, u_g, jg, rho_g jg,
  !! the local exponents d ln C0/d ln alpha and d ln jg/d ln alpha, and the
  !! velocity recovery of chexal_recover_velocities from the conservative momentum,
  !! including alpha = 0 exactly (storage convention u_g = 0, flux alpha*u_g = 0).
  use mod_kinds,only:wp
  use module_chexal_1992,only:chexal_coefficients_si,chexal_recover_velocities
  implicit none
  real(wp),parameter::rl=997.137_wp,rg=2.524_wp,jl=0.5264_wp,d=0.05254_wp,g=9.81_wp,mul=0.00089_wp,mug=0.0000185_wp, &
    sig=0.0728_wp,th=2.334_wp*acos(-1.0_wp)/180
  real(wp)::a,ug,jg,c0,v,c0o,jgo,ao,ul,ug2,mom
  integer::k,it,st,st2
  ao=-1
  do k=1,12
    a=10.0_wp**(-k)
    ug=0.35_wp
    do it=1,200
      jg=a*ug
      call chexal_coefficients_si(a,rl,rg,jl,jg,th,d,g,mul,mug,sig,c0,v,st)
      if(st/=0)exit
      ug=c0*(jl+jg)+v
    enddo
    jg=a*ug
    mom=rl*(1-a)*(jl/(1-a))+rg*jg
    call chexal_recover_velocities(a,rl,rg,mom,th,d,g,mul,mug,sig,ul,ug2,st2)
    if(ao>0)then
      write(*,'(a,es8.1,a,i3,a,2es14.6,a,es14.6,a,2es12.4,a,2f8.4,a,es12.4,i4)')'ALPHA|',a,'|',st,'|C0,Vgj ',c0,v,'|u_g ',ug, &
        '|jg,rho_g jg ',jg,rg*jg,'|dlnC0/dlna,dlnjg/dlna ',log(c0/c0o)/log(a/ao),log(jg/jgo)/log(a/ao),'|recovered u_g ',ug2,st2
    else
      write(*,'(a,es8.1,a,i3,a,2es14.6,a,es14.6,a,2es12.4)')'ALPHA|',a,'|',st,'|C0,Vgj ',c0,v,'|u_g ',ug,'|jg,rho_g jg ',jg,rg*jg
    endif
    c0o=c0;jgo=jg;ao=a
  enddo
  call chexal_coefficients_si(0.0_wp,rl,rg,jl,0.0_wp,th,d,g,mul,mug,sig,c0,v,st)
  write(*,'(a,i3,a,2es14.6)')'ALPHA|0 exactly|',st,'|C0,Vgj (limit coefficients at jg=0) ',c0,v
  call chexal_recover_velocities(0.0_wp,rl,rg,rl*jl,th,d,g,mul,mug,sig,ul,ug2,st2)
  write(*,'(a,i3,a,2es14.6)')'RECOVER|alpha=0|',st2,'|u_l, u_g (storage convention) ',ul,ug2
  ! jg -> 0- with liquid upward: outside the correlation's domain at this inclination
  call chexal_coefficients_si(1e-4_wp,rl,rg,jl,-1e-6_wp,th,d,g,mul,mug,sig,c0,v,st)
  write(*,'(a,i3)')'JG_NEGATIVE|status (61 = CHEXAL_DOMAIN) ',st
end program
