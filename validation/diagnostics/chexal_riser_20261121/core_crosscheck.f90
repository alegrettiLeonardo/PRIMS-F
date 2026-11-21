program core_crosscheck
  !! Reads lines "alpha rho_l rho_g theta jl jg branch" (branch 1 = HIGH, else AUTO) (SI, theta rad from horizontal)
  !! and prints the Fortran core C0, Vgj and status (AUTO branch). Diagnostic.
  use mod_kinds, only: wp
  use module_chexal_1992, only: chexal_coefficients_si, CHEXAL_BRANCH_AUTO, CHEXAL_BRANCH_HIGH
  implicit none
  real(wp) :: a, rl, rg, th, jl, jg, c0, v
  integer :: st, ios, br
  do
    read(*, *, iostat=ios) a, rl, rg, th, jl, jg, br
    if (ios /= 0) exit
    call chexal_coefficients_si(a, rl, rg, jl, jg, th, 0.05254_wp, 9.81_wp, 0.89e-3_wp, 1.85e-5_wp, 0.0728_wp, c0, v, st, &
      branch=merge(CHEXAL_BRANCH_HIGH, CHEXAL_BRANCH_AUTO, br == 1))
    write(*, '(i4,2es25.16)') st, c0, v
  enddo
end program
