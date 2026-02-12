module module_boundary_conditions
  !! Boundary conditions and ghost-cell filling for PRIMS-F.
  !!
  !! IMPORTANT INDEXING CONVENTION (1-based):
  !!   - total cells include ghosts: N = size(U,1)
  !!   - left ghosts:   1 .. nd
  !!   - physical:      nd+1 .. N-nd
  !!   - right ghosts:  N-nd+1 .. N
  !!
  !! This module is written defensively to avoid out-of-bounds access
  !! (a common source of SIGSEGV when porting from 0-based Python).

  use mod_kinds,     only: wp, i4
  use module_config, only: sim_config_t
  use module_find,   only: compute_primitive_variables
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: apply_ghost_cells
  public :: apply_boundary_conditions

  ! Cache a reference inlet state so SS1/SS2/SS3 can keep a strong inlet
  ! Dirichlet BC (the solver excludes the inlet physical cell from Newton).
  logical, save :: inlet_cached = .false.
  real(wp), save :: u1_in_ref = 0.0_wp
  real(wp), save :: u2_in_ref = 0.0_wp

contains

  subroutine apply_ghost_cells(U, nd)
    !! Simple and SAFE ghost fill: constant extrapolation from nearest physical cells.
    !! This is robust and avoids any access to index 0/-1 or N+1.
    real(wp), intent(inout) :: U(:,:)
    integer(i4), intent(in) :: nd
    integer(i4) :: N, iL, iR, k, idx
    real(wp) :: bcmin

    N = size(U,1)
    if (N < 1) return
    if (nd < 1) return
    if (N <= 2*nd) return

    iL = nd + 1
    iR = N - nd
    bcmin = 1.0e-12_wp

    ! Left ghosts
    do k = 1, nd
      idx = iL - k
      if (idx < 1) exit
      U(idx,:) = U(iL,:)
    end do

    ! Right ghosts
    do k = 1, nd
      idx = iR + k
      if (idx > N) exit
      U(idx,:) = U(iR,:)
    end do

    ! Positivity clamp (helps avoid NaNs inside primitive recovery)
    U(:,1) = max(U(:,1), bcmin)
    U(:,2) = max(U(:,2), bcmin)
  end subroutine apply_ghost_cells


  subroutine apply_boundary_conditions(U, inlet_idx, outlet_idx, theta_inlet, theta_outlet, cfg)
    !! Applies:
    !!   - Inlet: strong mass-flux BC (cfg%mul, cfg%mug) by fixing the inlet
    !!            physical cell conservative state (and copying to left ghosts).
    !!   - Outlet: separator pressure BC (cfg%nPs) via right ghost cells.
    !!
    !! NOTES:
    !!   1) The solver residual at the outlet uses U(outlet_idx+1,:) (ghost).
    !!   2) The solver may exclude inlet_idx from Newton for non-OSC regimes;
    !!      therefore we keep inlet_idx fixed.
    real(wp), intent(inout)        :: U(:,:)
    integer(i4), intent(in)        :: inlet_idx, outlet_idx
    real(wp), intent(in)           :: theta_inlet, theta_outlet
    type(sim_config_t), intent(in) :: cfg

    integer(i4) :: N, nd, k, idx
    real(wp) :: bcmin
    real(wp) :: P_bc, rhog_bc, rhol_bc, alpha_out
    real(wp) :: ul_out, ug_out, P_out, rhol_out, rhog_out
    integer(i4) :: phase_out
    real(wp) :: u1_bc, u2_bc, u3_bc
    real(wp), parameter :: rho_eps = 1.0e-16_wp
    real(wp), parameter :: alpha_eps = 1.0e-12_wp

    N  = size(U,1)
    nd = cfg%num_dummy_cells
    if (N < 1) return
    if (nd < 0) return

    bcmin = max(cfg%bc_min_value, 1.0e-16_wp)

    ! -----------------------------
    ! Inlet strong BC (fix inlet cell)
    ! -----------------------------
    if (inlet_idx >= 1 .and. inlet_idx <= N) then
      if (.not. inlet_cached) then
        u1_in_ref = max(U(inlet_idx,1), bcmin)
        u2_in_ref = max(U(inlet_idx,2), bcmin)
        inlet_cached = .true.
      end if

      U(inlet_idx,1) = max(u1_in_ref, bcmin)
      U(inlet_idx,2) = max(u2_in_ref, bcmin)
      ! Enforce inlet phase mass fluxes: m_l = u1*ul, m_g = u2*ug.
      ! By construction, mixture momentum becomes m_l + m_g.
      U(inlet_idx,3) = cfg%mul + cfg%mug

      ! Fill left ghost cells as copies of inlet state
      do k = 1, nd
        idx = inlet_idx - k
        if (idx < 1) exit
        U(idx,:) = U(inlet_idx,:)
      end do
    end if

    ! -----------------------------
    ! Outlet pressure BC via right ghosts
    ! -----------------------------
    if (outlet_idx >= 1 .and. outlet_idx <= N) then
      ! Recover outlet primitives (to get alpha and velocities)
      call compute_primitive_variables( &
        U(outlet_idx,1), U(outlet_idx,2), U(outlet_idx,3), theta_outlet, &
        cfg%nCg, cfg%nCl, cfg%nrho_l0, cfg%nP_l0, &
        cfg%D_H, cfg%AREA, cfg%eps, cfg%G, &
        cfg%mu_l, cfg%mu_g, cfg%sigma, cfg%njl, &
        cfg%omega_u, cfg%omega_rho, &
        cfg%tol, cfg%tola, cfg%maxit, &
        alpha_out, rhol_out, rhog_out, P_out, ul_out, ug_out, phase_out )

      if (.not. ieee_is_finite(alpha_out)) alpha_out = 0.5_wp
      alpha_out = max(min(alpha_out, 1.0_wp - alpha_eps), alpha_eps)
      if (.not. ieee_is_finite(ul_out)) ul_out = 0.0_wp
      if (.not. ieee_is_finite(ug_out)) ug_out = 0.0_wp

      P_bc    = max(cfg%nPs, 1.0e-12_wp)
      rhog_bc = max(P_bc / max(cfg%nCg*cfg%nCg, 1.0e-30_wp), rho_eps)
      rhol_bc = max(cfg%nrho_l0 + (P_bc - cfg%nP_l0) / max(cfg%nCl*cfg%nCl, 1.0e-30_wp), rho_eps)

      u1_bc = max((1.0_wp - alpha_out) * rhol_bc, bcmin)
      u2_bc = max(alpha_out * rhog_bc, bcmin)
      u3_bc = u1_bc * ul_out + u2_bc * ug_out

      ! Fill right ghost cells
      do k = 1, nd
        idx = outlet_idx + k
        if (idx > N) exit
        U(idx,1) = u1_bc
        U(idx,2) = u2_bc
        U(idx,3) = u3_bc
      end do
    end if

    ! Final safety clamp
    U(:,1) = max(U(:,1), bcmin)
    U(:,2) = max(U(:,2), bcmin)
  end subroutine apply_boundary_conditions

end module module_boundary_conditions
