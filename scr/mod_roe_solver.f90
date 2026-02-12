module mod_roe_solver
  !! Thin compatibility wrapper around `module_roe_solver`.
  !!
  !! The heavy lifting (adaptive characteristic MUSCL + Roe/Rusanov flux)
  !! lives in `module_roe_solver`. Here we expose a single entry point
  !! `roe_flux` that computes the Roe numerical flux given already
  !! reconstructed interface states and the Roe matrix.

  use mod_kinds,        only: wp
  use module_roe_solver, only: roe_solver_t
  implicit none
  private
  public :: roe_solver_t, roe_flux

contains

  subroutine roe_flux(solver, F_L, F_R, U_L, U_R, RoeA, flux)
    class(roe_solver_t), intent(in) :: solver
    real(wp), intent(in) :: F_L(3), F_R(3)
    real(wp), intent(in) :: U_L(3), U_R(3)
    real(wp), intent(in) :: RoeA(3,3)
    real(wp), intent(out) :: flux(3)

    call solver%roe_riemann_solver(F_L, F_R, U_L, U_R, RoeA, flux)
  end subroutine roe_flux

end module mod_roe_solver
