module module_solver_utils
  !! Fortran translation of solver_utils.py (core fluid-dynamics helpers).
  use mod_kinds, only: wp, i4
  use mod_void_helpers, only: ffan
  implicit none
  private

  public :: fluid_dynamics_t

  type :: fluid_dynamics_t
  contains
    procedure :: adimensionalizar
    procedure :: compute_conservative_variables
    procedure :: compute_flux_vector
    procedure :: compute_mixture_parameters
    procedure :: compute_source_term_vector
  end type fluid_dynamics_t

contains

  subroutine adimensionalizar(this, Cg, Cl, P_l0, rho_l0, Ps, Pb, jl, jg, mul, mug, &
                              nCg, nCl, nP_l0, nrho_l0, nPs, nPb, njl, njg, nmul, nmug, &
                              omega_c, omega_P, omega_rho, omega_u)
    class(fluid_dynamics_t), intent(in) :: this
    real(wp), intent(in) :: Cg, Cl, P_l0, rho_l0, Ps, Pb, jl, jg, mul, mug
    real(wp), intent(out) :: nCg, nCl, nP_l0, nrho_l0, nPs, nPb, njl, njg, nmul, nmug
    real(wp), intent(out) :: omega_c, omega_P, omega_rho, omega_u
    real(wp) :: omega_u_local

    omega_P   = P_l0
    omega_rho = rho_l0
    omega_c   = sqrt(omega_rho/omega_P)
    omega_u_local = max(jl, jg)
    if (omega_u_local <= 0.0_wp) omega_u_local = 1.0_wp
    omega_u = omega_u_local

    nCg     = Cg * omega_c
    nCl     = Cl * omega_c
    nP_l0   = P_l0 / omega_P
    nrho_l0 = rho_l0 / omega_rho
    nPs     = Ps / omega_P
    nPb     = Pb / omega_P
    njl     = jl / omega_u
    njg     = jg / omega_u
    nmul    = mul / (omega_rho * omega_u)
    nmug    = mug / (omega_rho * omega_u)
  end subroutine adimensionalizar


  subroutine compute_conservative_variables(this, alpha, rho_l, ul, rho_g, ug, U)
    class(fluid_dynamics_t), intent(in) :: this
    real(wp), intent(in) :: alpha, rho_l, ul, rho_g, ug
    real(wp), intent(out) :: U(3)
    U(1) = (1.0_wp - alpha) * rho_l
    U(2) = alpha * rho_g
    U(3) = (1.0_wp - alpha) * rho_l * ul + alpha * rho_g * ug
  end subroutine compute_conservative_variables


  subroutine compute_flux_vector(this, alpha, rho_l, rho_g, u_l, u_g, P, F)
    class(fluid_dynamics_t), intent(in) :: this
    real(wp), intent(in) :: alpha, rho_l, rho_g, u_l, u_g, P
    real(wp), intent(out) :: F(3)
    F(1) = (1.0_wp - alpha) * rho_l * u_l
    F(2) = alpha * rho_g * u_g
    F(3) = (1.0_wp - alpha) * rho_l * (u_l*u_l) + alpha * rho_g * (u_g*u_g) + P
  end subroutine compute_flux_vector


  subroutine compute_mixture_parameters(this, alpha_i, rho_l_i, rho_g_i, u_l_i, u_g_i, D, edia, mu_g, mu_l, &
                                        rho_m, f_m, j, Re_m)
    class(fluid_dynamics_t), intent(in) :: this
    real(wp), intent(in) :: alpha_i, rho_l_i, rho_g_i, u_l_i, u_g_i
    real(wp), intent(in) :: D, edia, mu_g, mu_l
    real(wp), intent(out) :: rho_m, f_m, j, Re_m
    real(wp) :: mu_m

    mu_m  = mu_g * alpha_i + (1.0_wp - alpha_i) * mu_l
    rho_m = alpha_i * rho_g_i + (1.0_wp - alpha_i) * rho_l_i
    j     = alpha_i * u_g_i + (1.0_wp - alpha_i) * u_l_i

    if (abs(mu_m) <= 1.0e-30_wp) then
      Re_m = 0.0_wp
      f_m  = 0.0_wp
    else
      Re_m = (rho_m * abs(j) * D) / mu_m
      f_m  = ffan(edia / D, Re_m)
    end if
  end subroutine compute_mixture_parameters


  subroutine compute_source_term_vector(this, alpha_i, rho_m_i, theta_i, f_m_i, j_i, g, D, S)
    class(fluid_dynamics_t), intent(in) :: this
    real(wp), intent(in) :: alpha_i, rho_m_i, theta_i, f_m_i, j_i, g, D
    real(wp), intent(out) :: S(3)
    S(1) = 0.0_wp
    S(2) = 0.0_wp
    S(3) = -rho_m_i * g * sin(theta_i) - (D/2.0_wp) * rho_m_i * f_m_i * abs(j_i) * j_i
  end subroutine compute_source_term_vector

end module module_solver_utils
