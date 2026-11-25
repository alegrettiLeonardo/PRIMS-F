program test_chexal_domain_contract
  !! P0.1 executable domain contract of the CHEXAL_1992 AIR_WATER path.
  !! Accepted: states inside the air-water data base of the source (Table I-I,
  !! Table 5-1: p 1-6.8 bar -> rho_g/rho_l <= 8.1e-3 at 293 K; D 0.010-0.305 m),
  !! co-current in any orientation, vertical countercurrent with a branch.
  !! Refused with CHEXAL_UNQUALIFIED (68): density ratio or diameter outside the
  !! base, in coefficients, residual, alpha root and velocity recovery.
  !! Refused with CHEXAL_DOMAIN (61): inclined countercurrent, horizontal negative
  !! flux, gas down / liquid up (p.262-265). The STEAM_WATER verification mode is
  !! not subject to the air-water contract.
  use mod_kinds, only: wp
  use module_chexal_1992
  implicit none
  real(wp), parameter :: PI=acos(-1.0_wp), RL=998.0_wp, G=9.81_wp, ML=8.9e-4_wp, MG=1.85e-5_wp, S=.0728_wp
  real(wp) :: c0, v, a, ul, ug, r
  integer :: st, fails, n
  fails=0; n=0
  ! --- accepted: PRIMS-F riser envelope (rho_g 1.2..6 kg/m3, D=0.05254, 0..76.47 deg)
  call co(RL, 1.2_wp, .05254_wp, 0.0_wp, .5_wp, .6_wp, CHEXAL_OK, 'PRIMS-F low pressure, horizontal base')
  call co(RL, 6.0_wp, .05254_wp, 76.47_wp, .5_wp, .6_wp, CHEXAL_OK, 'PRIMS-F 5 bar, catenary top')
  call co(RL, 8.0_wp, .05254_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_OK, 'ratio 8.0e-3 inside the base')
  call co(RL, 2.4_wp, 0.010_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_OK, 'D at the lower edge of the base')
  call co(RL, 2.4_wp, 0.305_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_OK, 'D at the upper edge of the base')
  call co(RL, 2.4_wp, .05254_wp, -30.0_wp, .5_wp, .6_wp, CHEXAL_OK, 'co-current along a descending axis')
  ! --- refused: outside the air-water data base
  call co(RL, 8.2_wp, .05254_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_UNQUALIFIED, 'rho_g/rho_l 8.2e-3 > base')
  call co(RL, 40.0_wp, .05254_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_UNQUALIFIED, 'high-pressure ratio (EPRI 1000 psia class)')
  call co(RL, 2.4_wp, 0.0099_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_UNQUALIFIED, 'D below the base')
  call co(RL, 2.4_wp, 0.31_wp, 45.0_wp, .5_wp, .6_wp, CHEXAL_UNQUALIFIED, 'D above the base')
  call chexal_solve_alpha(.5_wp, .6_wp, RL, 40.0_wp, PI/4, .05254_wp, G, ML, MG, S, a, st)
  call chk(st == CHEXAL_UNQUALIFIED, 'alpha root refused outside the base')
  call chexal_residual_si(.3_wp, RL, 40.0_wp, .5_wp, .6_wp, PI/4, .05254_wp, G, ML, MG, S, r, st)
  call chk(st == CHEXAL_UNQUALIFIED, 'residual refused outside the base')
  call chexal_recover_velocities(.3_wp, RL, 40.0_wp, 500.0_wp, PI/4, .05254_wp, G, ML, MG, S, ul, ug, st)
  call chk(st == CHEXAL_UNQUALIFIED, 'velocity recovery refused outside the base')
  call chexal_recover_velocities(.3_wp, RL, 2.4_wp, 500.0_wp, PI/4, 0.5_wp, G, ML, MG, S, ul, ug, st)
  call chk(st == CHEXAL_UNQUALIFIED, 'velocity recovery refused for D outside the base')
  ! --- refused: outside the source orientation domain (unchanged)
  call co(RL, 2.4_wp, .05254_wp, 65.0_wp, -.1_wp, .4_wp, CHEXAL_DOMAIN, 'inclined countercurrent (liquid down)')
  call co(RL, 2.4_wp, .05254_wp, 0.0_wp, .5_wp, -.1_wp, CHEXAL_DOMAIN, 'horizontal reverse gas (junction trace)')
  call co(RL, 2.4_wp, .05254_wp, 45.0_wp, .5_wp, -.1_wp, CHEXAL_DOMAIN, 'gas down / liquid up')
  call co(RL, 2.4_wp, .05254_wp, 90.0_wp, -.1_wp, .4_wp, CHEXAL_BRANCH_REQUIRED, 'vertical countercurrent needs a branch')
  ! --- the steam-water verification mode is outside the air-water contract
  call chexal_coefficients_si(.5_wp, 740.0_wp, 36.0_wp, .5_wp, .6_wp, PI/2, .05254_wp, G, 9.5e-5_wp, 1.9e-5_wp, .0178_wp, c0, v, st, &
    fluid=CHEXAL_FLUID_STEAM_WATER, pressure=1000.0_wp, pcrit=3208.2_wp)
  call chk(st == CHEXAL_OK, 'steam-water verification mode not subject to the air-water base')
  if (fails > 0) then
    write(*,'(a,i0,a,i0)') 'FAILED ', fails, ' of ', n; error stop 1
  endif
  write(*,'(a,i0,a)') 'PASS Chexal air-water domain contract: ', n, ' checks'
contains
  subroutine co(rl_, rg_, d_, deg, jl, jg, expect, label)
    real(wp), intent(in) :: rl_, rg_, d_, deg, jl, jg
    integer, intent(in) :: expect
    character(*), intent(in) :: label
    call chexal_coefficients_si(.4_wp, rl_, rg_, jl, jg, deg*PI/180, d_, G, ML, MG, S, c0, v, st)
    if (st /= expect) write(*,'(a,i0,a,i0)') '  status ', st, ' expected ', expect
    call chk(st == expect, label)
  end subroutine
  subroutine chk(ok, label)
    logical, intent(in) :: ok
    character(*), intent(in) :: label
    n = n + 1
    if (.not. ok) then
      fails = fails + 1; write(*,'(2a)') 'FAIL: ', label
    endif
  end subroutine
end program
