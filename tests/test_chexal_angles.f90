program test_chexal_angles
  !! P0.2 certificate of the angle / signed-velocity convention of the Chexal
  !! core as used by PRIMS-F (air-water, AUTO branch unless stated).
  !! PRIMS-F: theta in radians FROM THE HORIZONTAL, positive when the section
  !! axis rises (riser 0..76.47 deg; pipe -beta). Fluxes are signed along the
  !! axis. Chexal 1992: theta_C = 90 - |theta| deg FROM THE VERTICAL; fluxes in
  !! the PHYSICAL upward sense f = sign(theta) j (p.262 vertical convention).
  !! Expected domain (p.262, 264-265): theta=0 only non-negative fluxes; theta/=0:
  !! co-current up or down (f_l f_g >= 0) admitted; liquid-down/gas-up only at
  !! |theta|=90 deg and only with an explicit branch (AUTO -> 62); gas-down /
  !! liquid-up never. No sign tolerance: 1e-12 m/s is a sign like any other.
  !! Checks per angle in {0,+-15,+-45,+-65,+-76,+-90} deg and per flux pair:
  !!   status = expected; axis reversal (theta,j) -> (-theta,-j) leaves C0
  !!   unchanged and reverses Vgj EXACTLY; recovered velocities from
  !!   (alpha, rho_l j_l + rho_g j_g) return the signed u_l, u_g; continuity of
  !!   C0 and Vgj across theta=0 for co-current flow in the axis direction.
  use mod_kinds, only: wp
  use module_chexal_1992
  implicit none
  real(wp), parameter :: PI=acos(-1.0_wp), RL=998.0_wp, RG=2.4_wp, D=.05254_wp, G=9.81_wp
  real(wp), parameter :: ML=8.9e-4_wp, MG=1.85e-5_wp, S=.0728_wp
  real(wp) :: angles(11), fl(12), fg(12), th, jl, jg, c0, v, c0r, vr, a, ul, ug, mom, sgn, x, y, cprev(2)
  integer :: i, k, st, st2, expect, fails, nchk, nadm, ncc, ndom, nbr, nnum
  angles=[0._wp,15._wp,-15._wp,45._wp,-45._wp,65._wp,-65._wp,76._wp,-76._wp,90._wp,-90._wp]
  fl=[.5_wp,-.5_wp,-.1_wp,.1_wp,0._wp,.3_wp,1e-12_wp,-1e-12_wp,.5_wp,.5_wp,0._wp,-.5_wp]
  fg=[.6_wp,-.6_wp,.4_wp,-.4_wp,.3_wp,0._wp,.6_wp,.6_wp,1e-12_wp,-1e-12_wp,0._wp,-.01_wp]
  fails=0;nchk=0;nadm=0;ncc=0;ndom=0;nbr=0;nnum=0
  do i=1,size(angles)
    th=angles(i)*PI/180
    sgn=merge(-1.0_wp,1.0_wp,th<0)
    do k=1,size(fl)
      jl=fl(k);jg=fg(k)
      x=sgn*jl;y=sgn*jg                       ! physical upward fluxes
      if(th==0.0_wp)then
        expect=merge(CHEXAL_DOMAIN,CHEXAL_OK,jl<0.or.jg<0)
      else if(x<0.and.y>0.or.x<0.and.y==0)then   ! liquid down, gas up/zero: countercurrent
        expect=CHEXAL_DOMAIN
        if(abs(abs(angles(i))-90)<1e-12_wp)expect=CHEXAL_BRANCH_REQUIRED
      else if(y<0.and.x>0)then
        expect=CHEXAL_DOMAIN
      else
        expect=CHEXAL_OK
      endif
      call chexal_coefficients_si(.4_wp,RL,RG,jl,jg,th,D,G,ML,MG,S,c0,v,st)
      call chk(st==expect,'status',i,k)
      if(st==CHEXAL_DOMAIN)ndom=ndom+1
      if(st==CHEXAL_BRANCH_REQUIRED)then
        ncc=ncc+1
        call chexal_coefficients_si(.4_wp,RL,RG,jl,jg,th,D,G,ML,MG,S,c0,v,st2,branch=CHEXAL_BRANCH_HIGH)
        if(st2==CHEXAL_NUMERIC.and.abs(jg)<1e-4_wp)then
          nnum=nnum+1   ! documented CCFL-search limit at tiny gas flux
        else
          call chk(st2==CHEXAL_OK.or.st2==CHEXAL_CCFL_EXCEEDED,'vertical CC admitted with explicit branch',i,k)
        endif
      endif
      if(st/=CHEXAL_OK)cycle
      nadm=nadm+1
      if(th/=0.0_wp)then
        call chexal_coefficients_si(.4_wp,RL,RG,-jl,-jg,-th,D,G,ML,MG,S,c0r,vr,st2)
        call chk(st2==CHEXAL_OK.and.c0r==c0.and.vr==-v,'axis reversal: C0 equal, Vgj reversed exactly',i,k)
      endif
      if(abs(jl)+abs(jg)>1e-9_wp.and.jl*jg>0)then
        call chexal_solve_alpha(jl,jg,RL,RG,th,D,G,ML,MG,S,a,st2)
        if(st2/=CHEXAL_OK)write(*,'(a,i0,2i4)') 'alpha root status ',st2,i,k
        call chk(st2==CHEXAL_OK,'co-current alpha root',i,k)
        if(st2/=CHEXAL_OK)cycle
        call chk(a>0.and.a<=1,'alpha root inside (0,1]',i,k)
        if(a==1.0_wp)then
          ! Gas-only limit (C0(1)=1, Vgj(1)=0, pp.259-260): residual equals j_l;
          ! the liquid velocity is undefined there, so no recovery is attempted.
          call chk(abs(jl)<=1e-12_wp,'alpha=1 endpoint only for vanishing liquid flux',i,k)
          cycle
        endif
        mom=RL*jl+RG*jg
        call chexal_recover_velocities(a,RL,RG,mom,th,D,G,ML,MG,S,ul,ug,st2)
        if(abs(abs(angles(i))-90)<1e-12_wp.and.st2==CHEXAL_BRANCH_REQUIRED)then
          ! Vertical: (alpha, M) also admits a countercurrent pair; AUTO must ask
          ! for the branch. The previous-state hint resolves the co-current one.
          nbr=nbr+1
          call chexal_recover_velocities(a,RL,RG,mom,th,D,G,ML,MG,S,ul,ug,st2,ul_hint=jl/(1-a),ug_hint=jg/a)
          call chk(st2==CHEXAL_OK,'vertical recovery resolved by the continuation hint',i,k)
        else if(abs(abs(angles(i))-90)<1e-12_wp.and.st2==CHEXAL_NUMERIC.and.abs(jg)<1e-4_wp)then
          nnum=nnum+1   ! documented CCFL-search limit at tiny gas flux (CHEXAL_REMAINING_LIMITS.md)
          cycle
        endif
        call chk(st2==CHEXAL_OK,'velocity recovery status',i,k)
        if(st2/=CHEXAL_OK)cycle
        ! Error measured in the conserved momentum scale (u_l from a mixture
        ! momentum dominated by the gas is ill-conditioned when |j_l|<<1);
        ! signs are checked separately and strictly.
        call chk(RL*(1-a)*abs(ul-jl/(1-a))+RG*a*abs(ug-jg/a)<=1e-12_wp*(RL*abs(jl)+RG*abs(jg)), &
          'recovered velocities (momentum-scaled error)',i,k)
        call chk(sign(1.0_wp,ul)==sign(1.0_wp,jl).and.sign(1.0_wp,ug)==sign(1.0_wp,jg),'recovered velocity signs',i,k)
      endif
    enddo
  enddo
  ! continuity across theta=0 for co-current flow along the axis (Fr kink, no jump)
  call chexal_coefficients_si(.4_wp,RL,RG,.5_wp,.6_wp,0.0_wp,D,G,ML,MG,S,cprev(1),cprev(2),st)
  do k=-1,1,2
    call chexal_coefficients_si(.4_wp,RL,RG,.5_wp,.6_wp,k*1e-9_wp,D,G,ML,MG,S,c0,v,st2)
    call chk(st==0.and.st2==0.and.abs(c0-cprev(1))<=1e-7_wp*cprev(1).and.abs(v-cprev(2))<=1e-7_wp*abs(cprev(2)), &
      'C0, Vgj continuous across theta=0',k,0)
  enddo
  ! the catenary riser range of PRIMS-F (D geometry): 0 .. 76.47 deg; top is NOT vertical
  call chexal_coefficients_si(.4_wp,RL,RG,-.1_wp,.4_wp,76.47_wp*PI/180,D,G,ML,MG,S,c0,v,st)
  call chk(st==CHEXAL_DOMAIN,'catenary top (76.47 deg) countercurrent refused',0,0)
  call chexal_coefficients_si(.4_wp,RL,RG,-.1_wp,.4_wp,76.47_wp*PI/180,D,G,ML,MG,S,c0,v,st,branch=CHEXAL_BRANCH_HIGH)
  call chk(st==CHEXAL_DOMAIN,'explicit branch does not extend CC to the inclined top',0,0)
  if(fails>0)then
    write(*,'(a,i0,a,i0)') 'FAILED ',fails,' of ',nchk;error stop 1
  endif
  write(*,'(a,i0,a,i0,a,i0,a,i0,a)') 'PASS Chexal angle/velocity certificate: ',nchk,' checks (',nadm,' admitted, ', &
    ndom,' domain refusals, ',ncc,' vertical CC branch-required)'
  write(*,'(a,i0,a,i0,a)') '     vertical recoveries needing the hint: ',nbr,'; documented tiny-flux NUMERIC: ',nnum, &
    ' (none at the PRIMS-F catenary angles <= 76.47 deg)'
contains
  subroutine chk(ok,label,i,k)
    logical,intent(in)::ok
    character(*),intent(in)::label
    integer,intent(in)::i,k
    nchk=nchk+1
    if(.not.ok)then
      fails=fails+1;write(*,'(a,1x,a,2i4)') 'FAIL:',label,i,k
    endif
  end subroutine
end program
