module module_simulate_pipe_riser
  !! Fortran translation (core skeleton) of simulatePipeRiser.py.
  !!
  !! This module provides a TYPE `pipeline_riser_simulator_t` mirroring the
  !! Python class `PipelineRiserSimulator`.
  !!
  !! Implicit upgrade (Euler Backward):
  !! - Time advance uses a Newton-Raphson loop on the implicit residual.
  !! - The linear Newton correction system is solved with a 3x3 block
  !!   tridiagonal Thomas algorithm.
  !! - The global Jacobian is approximated using frozen numerical-flux
  !!   linearization (Roe/Rusanov), yielding a block tridiagonal matrix.

  use mod_kinds,               only: wp, i4
  use module_config,           only: sim_config_t, first_physical_cell, last_physical_cell
  use module_roe_solver, only: roe_solver_t, roe_absA_from_A
  use module_time_step,        only: time_integrator_t
  use module_solver_utils,     only: fluid_dynamics_t
  use module_matrix_roe,       only: calculate_roe_matrix
  use module_lapack_small,     only: dgeev_3x3
  use module_jacobian,         only: DFduj
  use module_find,             only: compute_primitive_variables
  use module_boundary_conditions, only: apply_boundary_conditions
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use omp_lib, only: omp_get_max_threads

  implicit none
  private

  public :: pipeline_riser_simulator_t

  type :: pipeline_riser_simulator_t
    type(sim_config_t) :: cfg
    type(roe_solver_t) :: roe_solver
    type(fluid_dynamics_t) :: fd
    type(time_integrator_t) :: ti
    real(wp) :: alpha_eps   = 1.0e-12_wp
    real(wp) :: min_dt      = 1.0e-12_wp
    real(wp) :: max_dt      = 1.0e-3_wp
    real(wp) :: dt_growth   = 1.2_wp
    real(wp) :: dt_shrink   = 0.5_wp
    integer(i4) :: newton_fast   = 3
    integer(i4) :: newton_slow   = 6
    integer(i4) :: max_backtrack = 20
    integer(i4) :: maxit_root    = 1000

    real(wp), allocatable :: ul_prev(:)
    real(wp), allocatable :: ug_prev(:)
    integer(i4), allocatable :: idx_prev(:)
  contains
    procedure :: init
    procedure :: ensure_prev_arrays
    procedure :: update_prev_arrays
    procedure :: process_cell
    procedure :: simulate
  end type pipeline_riser_simulator_t

  logical, save :: banner_printed = .false.

contains


  pure function to_lower(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: out
    integer :: i, c
    out = str
    do i = 1, len(str)
      c = iachar(out(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) out(i:i) = achar(c + 32)
    end do
  end function to_lower

  subroutine init(this, cfg)
    class(pipeline_riser_simulator_t), intent(inout) :: this
    type(sim_config_t), intent(in) :: cfg


    integer :: iu
    character(len=256) :: fvisit

    this%cfg = cfg

    call execute_command_line('mkdir -p ' // trim(cfg%out_dir))
    if (cfg%save_every > 0) then
      fvisit = trim(cfg%out_dir)//'/results_'//trim(cfg%regime)//'.visit'
      open(newunit=iu, file=trim(fvisit), status='replace', action='write', form='formatted')
      write(iu,'(a)') '!NBLOCKS 1'
      close(iu)
    end if

    ! --- time-step bounds from configuration ---
    this%min_dt = max(cfg%dt_min, 1.0e-12_wp)
    if (cfg%dt_max > 0.0_wp) then
      this%max_dt = cfg%dt_max
    else
      this%max_dt = max(cfg%T / 10.0_wp, 10.0_wp*cfg%dt_min)
    end if

    this%maxit_root = cfg%maxit

    ! ------------------------------------------------------------------
    ! CRITICAL (segfault fix):
    !   process_cell() always passes idx_prev/ul_prev/ug_prev into
    !   compute_primitive_variables() as optional warm-start/hysteresis
    !   inputs. If these arrays are not allocated yet, any access like
    !     this%idx_prev(i)
    !   will dereference an unallocated allocatable -> SIGSEGV.
    !
    !   The dt-estimate loop at the start of simulate() calls process_cell()
    !   before update_prev_arrays() runs, so we must guarantee allocation
    !   here.
    ! ------------------------------------------------------------------
    if (cfg%nCells > 0) then
      call this%ensure_prev_arrays(cfg%nCells)
    end if
  end subroutine init


  subroutine ensure_prev_arrays(this, n_cells)
    class(pipeline_riser_simulator_t), intent(inout) :: this
    integer(i4), intent(in) :: n_cells

    if (.not. allocated(this%ul_prev) .or. size(this%ul_prev) /= n_cells) then
      if (allocated(this%ul_prev)) deallocate(this%ul_prev)
      if (allocated(this%ug_prev)) deallocate(this%ug_prev)
      if (allocated(this%idx_prev)) deallocate(this%idx_prev)
      allocate(this%ul_prev(n_cells), this%ug_prev(n_cells), this%idx_prev(n_cells))
      this%ul_prev = 0.0_wp
      this%ug_prev = 0.0_wp
      this%idx_prev = 0
    end if
  end subroutine ensure_prev_arrays


  subroutine update_prev_arrays(this, U, thetav)
    !! Updates ul_prev/ug_prev in ALL cells (including ghosts), as in Python.
    class(pipeline_riser_simulator_t), intent(inout) :: this
    real(wp), intent(in) :: U(:,:)
    real(wp), intent(in) :: thetav(:)

    integer(i4) :: k, n
    real(wp) :: alpha, rhol, rhog, P, ul, ug
    integer(i4) :: idx

    n = size(U,1)
    call this%ensure_prev_arrays(n)

    do k = 1, n
      call compute_primitive_variables( &
        U(k,1), U(k,2), U(k,3), thetav(k), &
        this%cfg%nCg, this%cfg%nCl, &
        this%cfg%nrho_l0, this%cfg%nP_l0, &
        this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, &
        this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
        this%cfg%omega_u, this%cfg%omega_rho, &
        this%cfg%tol, this%cfg%tola, this%maxit_root, &
        alpha, rhol, rhog, P, ul, ug, idx, &
        idx_prev=this%idx_prev(k), ul_prev=this%ul_prev(k), ug_prev=this%ug_prev(k) )

      if (ieee_is_finite(ul)) this%ul_prev(k) = ul
      if (ieee_is_finite(ug)) this%ug_prev(k) = ug
      this%idx_prev(k) = idx
    end do
  end subroutine update_prev_arrays


  subroutine process_cell(this, i, U, delta_s, thetav, CFL, flux, source_avg, dt, dflux_dUL, dflux_dUR, dSrc_dU)
    !! Numerical flux/source at face i | i+1.
    !! Optionally returns frozen linearization of the numerical flux:
    !!   dflux_dUL = d(flux_i)/d(U_i)
    !!   dflux_dUR = d(flux_i)/d(U_{i+1})

    class(pipeline_riser_simulator_t), intent(inout) :: this
    integer(i4), intent(in) :: i
    real(wp), intent(in) :: U(:,:)        ! (N,3)
    real(wp), intent(in) :: delta_s(:)    ! (N-1)
    real(wp), intent(in) :: thetav(:)     ! (N)
    real(wp), intent(in) :: CFL
    real(wp), intent(out) :: flux(3)
    real(wp), intent(out) :: source_avg(3)
    real(wp), intent(out) :: dt
    real(wp), intent(out), optional :: dflux_dUL(3,3)
    real(wp), intent(out), optional :: dflux_dUR(3,3)
    real(wp), intent(out), optional :: dSrc_dU(3,3)

    integer(i4) :: nd, n_cells
    integer(i4) :: idx_pb, buffer_cells
    real(wp) :: theta_L, theta_R, theta_m
    real(wp) :: U_stencil(5,3)
    real(wp) :: U_L(3), U_R(3)
    real(wp) :: alpha_L, rho_l_L, rho_g_L, P_L, u_l_L, u_g_L
    real(wp) :: alpha_R, rho_l_R, rho_g_R, P_R, u_l_R, u_g_R
    integer(i4) :: idx
    real(wp) :: F_L(3), F_R(3)
    real(wp) :: RoeA(3,3), RoeA_char(3,3)
    real(wp) :: rho_m_L, f_m_L, j_L, Re_L
    real(wp) :: rho_m_R, f_m_R, j_R, Re_R
    real(wp) :: S_L(3), S_R(3)
    real(wp) :: P_scale
    logical :: interface, near_interface
    real(wp) :: amax
    real(wp) :: wr(3), wi(3), VR(3,3)
    integer :: info, k
    integer(i4) :: info_absA
    real(wp) :: mag
    real(wp) :: ds
    real(wp) :: u_material_nd

    ! Jacobian helpers
    real(wp) :: DF_L(3,3), DF_R(3,3)
    real(wp) :: absA(3,3), invVR(3,3), Dabs(3,3), I3(3,3)
    real(wp), parameter :: imag_tol = 1.0e-10_wp

    nd = this%cfg%num_dummy_cells
    n_cells = size(U,1)

    flux = 0.0_wp
    source_avg = 0.0_wp
    dt = this%min_dt

    if (present(dflux_dUL)) dflux_dUL = 0.0_wp
    if (present(dflux_dUR)) dflux_dUR = 0.0_wp
    if (present(dSrc_dU))   dSrc_dU   = 0.0_wp
    if (i < 1) return
    if (i+1 > n_cells) return

    theta_L = thetav(i)
    theta_R = thetav(i+1)
    ! ------------------------------------------------------------------
    ! Reconstruction
    !   - If cfg%use_muscl=.false. -> always first-order
    !   - If stencil not available (near boundaries) -> first-order
    !   - Keep first-order near the pipe-riser base (idx_pb +/- 1)
    ! ------------------------------------------------------------------

    ! Python: idx_pb = num_dummy_cells + N_pipe - 1 (0-based)
    ! Fortran (1-based): base-of-riser face ~ nd + N_pipe
    idx_pb = nd + this%cfg%N_pipe

    if (.not. this%cfg%use_muscl) then
      U_L = U(i,   :)
      U_R = U(i+1, :)
    else if ( (i-2 < 1) .or. (i+2 > n_cells) .or. (nd < 2) ) then
      ! Not enough neighbors for 5-point stencil
      U_L = U(i,   :)
      U_R = U(i+1, :)
    else if (abs(i - idx_pb) <= 1) then
      ! Keep first-order near pipe-riser transition.
      U_L = U(i,   :)
      U_R = U(i+1, :)
    else
      U_stencil = U(i-2:i+2, :)

      ! MUSCL in characteristic space (frozen for Jacobian linearization).
      theta_m = 0.5_wp*(theta_L + theta_R)
      RoeA_char = calculate_roe_matrix( &
          U(i,:), U(i+1,:), this%cfg%nCl, this%cfg%nCg, theta_m, this%cfg%nrho_l0, this%cfg%nP_l0, &
          this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
          this%cfg%omega_u, this%cfg%omega_rho, this%cfg%omega_P, this%cfg%tol, this%cfg%tola, this%cfg%quad_points, &
          ul_prev=0.5_wp*(this%ul_prev(i)+this%ul_prev(i+1)), &
          ug_prev=0.5_wp*(this%ug_prev(i)+this%ug_prev(i+1)) )

      call this%roe_solver%muscl_reconstruction_characteristic(U_stencil, RoeA_char, U_L, U_R)
    end if
    ! Convert to primitives.
    call compute_primitive_variables( &
      U_L(1), U_L(2), U_L(3), theta_L, &
      this%cfg%nCg, this%cfg%nCl, this%cfg%nrho_l0, this%cfg%nP_l0, &
      this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, &
      this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
      this%cfg%omega_u, this%cfg%omega_rho, &
      this%cfg%tol, this%cfg%tola, this%maxit_root, &
      alpha_L, rho_l_L, rho_g_L, P_L, u_l_L, u_g_L, idx, &
      idx_prev=this%idx_prev(i), ul_prev=this%ul_prev(i), ug_prev=this%ug_prev(i) )

    call compute_primitive_variables( &
      U_R(1), U_R(2), U_R(3), theta_R, &
      this%cfg%nCg, this%cfg%nCl, this%cfg%nrho_l0, this%cfg%nP_l0, &
      this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, &
      this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
      this%cfg%omega_u, this%cfg%omega_rho, &
      this%cfg%tol, this%cfg%tola, this%maxit_root, &
      alpha_R, rho_l_R, rho_g_R, P_R, u_l_R, u_g_R, idx, &
      idx_prev=this%idx_prev(i+1), ul_prev=this%ul_prev(i+1), ug_prev=this%ug_prev(i+1) )

    ! Clip alpha.
    alpha_L = max(this%alpha_eps, min(1.0_wp-this%alpha_eps, alpha_L))
    alpha_R = max(this%alpha_eps, min(1.0_wp-this%alpha_eps, alpha_R))

    ! Mixture (superficial) velocities in nondimensional form (used for dt estimate)
    j_L = (1.0_wp - alpha_L) * u_l_L + alpha_L * u_g_L
    j_R = (1.0_wp - alpha_R) * u_l_R + alpha_R * u_g_R

    ! Physical flux vectors (pressure scaling kept consistent with existing fd%compute_flux_vector usage).
    P_scale = this%cfg%omega_P / (this%cfg%omega_rho * (this%cfg%omega_u**2))
    call this%fd%compute_flux_vector(alpha_L, rho_l_L, rho_g_L, u_l_L, u_g_L, P_L*P_scale, F_L)
    call this%fd%compute_flux_vector(alpha_R, rho_l_R, rho_g_R, u_l_R, u_g_R, P_R*P_scale, F_R)

    theta_m = 0.5_wp*(theta_L + theta_R)
    RoeA = calculate_roe_matrix( &
        U_L, U_R, this%cfg%nCl, this%cfg%nCg, theta_m, this%cfg%nrho_l0, this%cfg%nP_l0, &
        this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
        this%cfg%omega_u, this%cfg%omega_rho, this%cfg%omega_P, this%cfg%tol, this%cfg%tola, this%cfg%quad_points, &
        ul_prev=0.5_wp*(this%ul_prev(i)+this%ul_prev(i+1)), &
        ug_prev=0.5_wp*(this%ug_prev(i)+this%ug_prev(i+1)) )
    ! ------------------------------------------------------------------
    ! Source term at CELL i (not face-average).
    ! This keeps the Newton Jacobian consistent with the residual term S(U_i).
    ! Also optionally computes a numerical Jacobian dS/dU at cell i.
    ! ------------------------------------------------------------------
    call compute_source_from_state(U(i,:), thetav(i), source_avg)

    if (present(dSrc_dU)) then
      call compute_source_jacobian_fd(U(i,:), thetav(i), source_avg, dSrc_dU)
    end if

    ! Interface switching: use Rusanov near pipe-riser or stratified/non-stratified transitions.
    interface = ((theta_L <= 0.0_wp .and. theta_R > 0.0_wp) .or. (theta_R <= 0.0_wp .and. theta_L > 0.0_wp))
    buffer_cells = this%cfg%interface_buffer_cells
    near_interface = interface .or. (abs(i - idx_pb) <= buffer_cells)

    ! Eigen info for dt and |A| (frozen)
    call dgeev_3x3(RoeA, wr, wi, VR, info)
    amax = 1.0_wp
    if (info == 0) then
      amax = 0.0_wp
      do k = 1, 3
        mag = sqrt(wr(k)*wr(k) + wi(k)*wi(k))
        if (ieee_is_finite(mag)) amax = max(amax, abs(mag))
      end do
      if (amax <= 0.0_wp) amax = 1.0_wp
    end if

    if (near_interface) then
      call this%roe_solver%rusanov_flux(F_L, F_R, U_L, U_R, amax, flux)
    else
      call this%roe_solver%roe_riemann_solver(F_L, F_R, U_L, U_R, RoeA, flux)
    end if
    ds = delta_s(min(i, size(delta_s)))
    ! In the nondimensional system, superficial velocities are already scaled.
    ! Using an extra division by omega_u artificially reduces dt.
    u_material_nd = max(abs(j_L), abs(j_R))
    dt = abs(this%ti%compute_time_step_lax_wendroff(CFL, ds, RoeA, mode='material', u_material=u_material_nd))
    if (dt < this%min_dt) dt = this%min_dt
    ! Frozen linearization of the numerical flux (for implicit Jacobian assembly).
    if (present(dflux_dUL) .or. present(dflux_dUR)) then
      if (this%cfg%jac_dissipation_fd) then
        ! High-quality approximate Jacobian:
        !   d( F̂(U_L,U_R) ) / dU_L , d( F̂(U_L,U_R) ) / dU_R
        ! computed by FD of the *same* numerical flux used (Roe/Rusanov),
        ! thus capturing the sensitivity of the Roe dissipation.
        call compute_flux_jacobian_fd(this, U_L, U_R, theta_L, theta_R, near_interface, flux, DF_L, DF_R)
        if (present(dflux_dUL)) dflux_dUL = DF_L
        if (present(dflux_dUR)) dflux_dUR = DF_R
      else
        ! Frozen Roe linearization (cheap, but inconsistent when Roe matrix varies strongly).
        I3 = 0.0_wp
        I3(1,1) = 1.0_wp; I3(2,2) = 1.0_wp; I3(3,3) = 1.0_wp

        DF_L = DFduj(alpha_L, rho_l_L, rho_g_L, u_l_L, u_g_L, this%cfg%nCl, this%cfg%nCg, theta_L, &
                    this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, &
                    this%cfg%omega_u, this%cfg%omega_rho, this%cfg%omega_P, this%cfg%tol)

        DF_R = DFduj(alpha_R, rho_l_R, rho_g_R, u_l_R, u_g_R, this%cfg%nCl, this%cfg%nCg, theta_R, &
                    this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, &
                    this%cfg%omega_u, this%cfg%omega_rho, this%cfg%omega_P, this%cfg%tol)
        absA = 0.0_wp

        ! Build |A| consistent with the Roe flux (entropy-fixed |lambda|).
        call roe_absA_from_A(RoeA, absA, amax, info_absA)

        if (near_interface .or. info /= 0 .or. info_absA /= 0) then
          ! Robust fallback (Rusanov) using an infinity-norm estimate of spectral radius
          amax = max(amax, maxval(sum(abs(RoeA), dim=2)))
          amax = max(amax, 1.0e-6_wp)
          absA = amax * I3
        else
          if (.not. all(ieee_is_finite(absA))) then
            amax = max(amax, maxval(sum(abs(RoeA), dim=2)))
            amax = max(amax, 1.0e-6_wp)
            absA = amax * I3
          end if
        end if

        if (present(dflux_dUL)) dflux_dUL = 0.5_wp*DF_L + 0.5_wp*absA
        if (present(dflux_dUR)) dflux_dUR = 0.5_wp*DF_R - 0.5_wp*absA
      end if
    end if



    contains

    subroutine compute_source_from_state(Uvec, theta, Svec)
      real(wp), intent(in)  :: Uvec(3)
      real(wp), intent(in)  :: theta
      real(wp), intent(out) :: Svec(3)

      real(wp) :: alpha, rhol, rhog, P, ul, ug
      integer(i4) :: idx_loc
      real(wp) :: rho_m, f_m, j_m, Re_m

      Svec = 0.0_wp

      call compute_primitive_variables( &
        Uvec(1), Uvec(2), Uvec(3), theta, &
        this%cfg%nCg, this%cfg%nCl, &
        this%cfg%nrho_l0, this%cfg%nP_l0, &
        this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, &
        this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
        this%cfg%omega_u, this%cfg%omega_rho, &
        this%cfg%tol, this%cfg%tola, this%maxit_root, &
        alpha, rhol, rhog, P, ul, ug, idx_loc, &
        idx_prev=this%idx_prev(i), ul_prev=this%ul_prev(i), ug_prev=this%ug_prev(i) )

      if (.not. (ieee_is_finite(alpha) .and. ieee_is_finite(rhol) .and. ieee_is_finite(rhog) .and. &
                ieee_is_finite(ul) .and. ieee_is_finite(ug))) then
        Svec = 0.0_wp
        return
      end if

      alpha = max(this%alpha_eps, min(1.0_wp-this%alpha_eps, alpha))
      rhol  = max(rhol, 1.0e-16_wp)
      rhog  = max(rhog, 1.0e-16_wp)

      call this%fd%compute_mixture_parameters( &
        alpha, &
        rhol*this%cfg%omega_rho, rhog*this%cfg%omega_rho, &
        ul*this%cfg%omega_u,     ug*this%cfg%omega_u, &
        this%cfg%D_H, this%cfg%eps, this%cfg%mu_g, this%cfg%mu_l, &
        rho_m, f_m, j_m, Re_m )

      call this%fd%compute_source_term_vector( &
        alpha, rho_m/this%cfg%omega_rho, theta, f_m, j_m/this%cfg%omega_u, &
        this%cfg%G * this%cfg%Lr / (this%cfg%omega_u**2), &
        this%cfg%Lr / this%cfg%D_H, &
        Svec )

      if (.not. all(ieee_is_finite(Svec))) Svec = 0.0_wp
    end subroutine compute_source_from_state

  subroutine eval_numerical_flux(this, U_Lx, U_Rx, theta_Lx, theta_Rx, near_interface_x, flux_x, ok)
    class(pipeline_riser_simulator_t), intent(inout) :: this
    real(wp), intent(in)  :: U_Lx(3), U_Rx(3)
    real(wp), intent(in)  :: theta_Lx, theta_Rx
    logical, intent(in)   :: near_interface_x
    real(wp), intent(out) :: flux_x(3)
    logical, intent(out)  :: ok

    real(wp) :: alpha_Lx, rho_l_Lx, rho_g_Lx, P_Lx, u_l_Lx, u_g_Lx
    real(wp) :: alpha_Rx, rho_l_Rx, rho_g_Rx, P_Rx, u_l_Rx, u_g_Rx
    integer(i4) :: idxL, idxR

    real(wp) :: F_Lx(3), F_Rx(3)
    real(wp) :: theta_m
    real(wp) :: RoeA(3,3)
    real(wp) :: wr(3), wi(3), VR(3,3)
    real(wp) :: amax
    integer(i4) :: info
    real(wp) :: P_scale
    real(wp) :: ul_guess_Lx, ug_guess_Lx, ul_guess_Rx, ug_guess_Rx
    real(wp), parameter :: mass_eps = 1.0e-12_wp

    ok = .false.
    flux_x = 0.0_wp

    ! Consistent pressure scaling with process_cell()
    P_scale = this%cfg%omega_P/(this%cfg%omega_rho*this%cfg%omega_u*this%cfg%omega_u)

    ! Provide a reasonable guess for ul_prev/ug_prev (helps root solver stability)
    ul_guess_Lx = U_Lx(3) / max(U_Lx(1)+U_Lx(2), mass_eps)
    ug_guess_Lx = ul_guess_Lx
    ul_guess_Rx = U_Rx(3) / max(U_Rx(1)+U_Rx(2), mass_eps)
    ug_guess_Rx = ul_guess_Rx

    call compute_primitive_variables( U_Lx(1), U_Lx(2), U_Lx(3), theta_Lx, &
        this%cfg%nCg, this%cfg%nCl, this%cfg%nrho_l0, this%cfg%nP_l0, &
        this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
        this%cfg%omega_u, this%cfg%omega_rho, this%cfg%tol, this%cfg%tola, this%maxit_root, &
        alpha_Lx, rho_l_Lx, rho_g_Lx, P_Lx, u_l_Lx, u_g_Lx, idxL, &
        ul_prev=ul_guess_Lx, ug_prev=ug_guess_Lx )

    call compute_primitive_variables( U_Rx(1), U_Rx(2), U_Rx(3), theta_Rx, &
        this%cfg%nCg, this%cfg%nCl, this%cfg%nrho_l0, this%cfg%nP_l0, &
        this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
        this%cfg%omega_u, this%cfg%omega_rho, this%cfg%tol, this%cfg%tola, this%maxit_root, &
        alpha_Rx, rho_l_Rx, rho_g_Rx, P_Rx, u_l_Rx, u_g_Rx, idxR, &
        ul_prev=ul_guess_Rx, ug_prev=ug_guess_Rx )

    if (idxL /= 1_i4 .or. idxR /= 1_i4) return

    call this%fd%compute_flux_vector(alpha_Lx, rho_l_Lx, rho_g_Lx, u_l_Lx, u_g_Lx, P_Lx*P_scale, F_Lx)
    call this%fd%compute_flux_vector(alpha_Rx, rho_l_Rx, rho_g_Rx, u_l_Rx, u_g_Rx, P_Rx*P_scale, F_Rx)

    theta_m = 0.5_wp*(theta_Lx + theta_Rx)

    RoeA = calculate_roe_matrix( U_Lx, U_Rx, this%cfg%nCl, this%cfg%nCg, theta_m, &
        this%cfg%nrho_l0, this%cfg%nP_l0, this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, &
        this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, this%cfg%omega_u, this%cfg%omega_rho, this%cfg%omega_P, &
        this%cfg%tol, this%cfg%tola, this%cfg%quad_points, &
        ul_prev=0.5_wp*(u_l_Lx+u_l_Rx), ug_prev=0.5_wp*(u_g_Lx+u_g_Rx) )

    call dgeev_3x3(RoeA, wr, wi, VR, info)
    if (info /= 0_i4) return
    amax = maxval(abs(wr))

    if (near_interface_x) then
      call this%roe_solver%rusanov_flux(F_Lx, F_Rx, U_Lx, U_Rx, amax, flux_x)
    else
      call this%roe_solver%roe_riemann_solver(F_Lx, F_Rx, U_Lx, U_Rx, RoeA, flux_x)
    end if

    ok = .true.
  end subroutine eval_numerical_flux


subroutine compute_flux_jacobian_fd(this, U_L0, U_R0, theta_L0, theta_R0, near_interface0, flux0, dFdUL, dFdUR)
  class(pipeline_riser_simulator_t), intent(inout) :: this
  real(wp), intent(in) :: U_L0(3), U_R0(3)
  real(wp), intent(in) :: theta_L0, theta_R0
  logical, intent(in)  :: near_interface0
  real(wp), intent(in) :: flux0(3)
  real(wp), intent(out) :: dFdUL(3,3), dFdUR(3,3)

  real(wp), parameter :: mass_eps = 1.0e-12_wp
  real(wp) :: epsU, h, htry
  real(wp) :: ULp(3), ULm(3), URp(3), URm(3)
  real(wp) :: Fp(3), Fm(3)
  logical :: okp, okm
  integer(i4) :: k, tries

  epsU = 1.0e-6_wp  ! FD step for flux Jacobian


  dFdUL = 0.0_wp
  dFdUR = 0.0_wp

  ! dF/dU_L
  do k = 1, 3
    h = epsU*max(1.0_wp, abs(U_L0(k)))
    h = max(h, 1.0e-12_wp)

    tries = 0
    do
      tries = tries + 1
      htry = h

      ULp = U_L0; ULm = U_L0
      ULp(k) = U_L0(k) + htry
      ULm(k) = U_L0(k) - htry

      okp = .true.; okm = .true.

      ! Keep positivity for the mass-like components.
      if (k <= 2) then
        if (ULp(k) < mass_eps) okp = .false.
        if (ULm(k) < mass_eps) okm = .false.
      end if

      if (okp) call eval_numerical_flux(this, ULp, U_R0, theta_L0, theta_R0, near_interface0, Fp, okp)
      if (okm) call eval_numerical_flux(this, ULm, U_R0, theta_L0, theta_R0, near_interface0, Fm, okm)

      if (okp .and. okm) then
        dFdUL(:,k) = (Fp - Fm)/(2.0_wp*htry)
        exit
      else if (okp .and. .not. okm) then
        dFdUL(:,k) = (Fp - flux0)/htry
        exit
      else if (.not. okp .and. okm) then
        dFdUL(:,k) = (flux0 - Fm)/htry
        exit
      else
        if (tries >= 2) then
          dFdUL(:,k) = 0.0_wp
          exit
        end if
        h = 0.5_wp*h
      end if
    end do
  end do

  ! dF/dU_R
  do k = 1, 3
    h = epsU*max(1.0_wp, abs(U_R0(k)))
    h = max(h, 1.0e-12_wp)

    tries = 0
    do
      tries = tries + 1
      htry = h

      URp = U_R0; URm = U_R0
      URp(k) = U_R0(k) + htry
      URm(k) = U_R0(k) - htry

      okp = .true.; okm = .true.

      if (k <= 2) then
        if (URp(k) < mass_eps) okp = .false.
        if (URm(k) < mass_eps) okm = .false.
      end if

      if (okp) call eval_numerical_flux(this, U_L0, URp, theta_L0, theta_R0, near_interface0, Fp, okp)
      if (okm) call eval_numerical_flux(this, U_L0, URm, theta_L0, theta_R0, near_interface0, Fm, okm)

      if (okp .and. okm) then
        dFdUR(:,k) = (Fp - Fm)/(2.0_wp*htry)
        exit
      else if (okp .and. .not. okm) then
        dFdUR(:,k) = (Fp - flux0)/htry
        exit
      else if (.not. okp .and. okm) then
        dFdUR(:,k) = (flux0 - Fm)/htry
        exit
      else
        if (tries >= 2) then
          dFdUR(:,k) = 0.0_wp
          exit
        end if
        h = 0.5_wp*h
      end if
    end do
  end do

  if (.not. all(ieee_is_finite(dFdUL))) dFdUL = 0.0_wp
  if (.not. all(ieee_is_finite(dFdUR))) dFdUR = 0.0_wp
end subroutine compute_flux_jacobian_fd

subroutine compute_source_jacobian_fd(U0, theta, S0, dSdU)
      !! Finite-difference Jacobian of source term at a single cell:
      !!   dSdU(:,m) = dS/dU_m
      !!
      !! Robust strategy:
      !! - Use centered differences when possible.
      !! - Fall back to one-sided (forward/backward) if the negative perturbation
      !!   would violate positivity (for mass-like conservative components) or
      !!   if the source evaluation fails (non-finite).
      !! - Retry once with a smaller step if needed.
      real(wp), intent(in)  :: U0(3)
      real(wp), intent(in)  :: theta
      real(wp), intent(in)  :: S0(3)
      real(wp), intent(out) :: dSdU(3,3)

      real(wp) :: Up(3), Um(3), Sp(3), Sm(3)
      real(wp) :: h, epsU, hmin
      integer  :: m, attempt
      logical  :: okp, okm, use_central

      dSdU = 0.0_wp

      ! Perturbation size: relative O(1e-7) with a small absolute floor.
      epsU = 1.0e-7_wp
      hmin = 1.0e-10_wp

      do m = 1, 3
        h = max(epsU * max(1.0_wp, abs(U0(m))), hmin)

        do attempt = 1, 2
          Up = U0
          Um = U0

          Up(m) = U0(m) + h

          ! Decide if centered FD is allowed.
          use_central = .true.

          ! For conservative "mass" components, do not allow negative perturbations.
          if (m == 1 .or. m == 2) then
            if (U0(m) <= h) then
              use_central = .false.
            else
              Um(m) = U0(m) - h
            end if
          else
            Um(m) = U0(m) - h
          end if

          call compute_source_from_state(Up, theta, Sp)
          okp = all(ieee_is_finite(Sp))

          if (use_central) then
            call compute_source_from_state(Um, theta, Sm)
            okm = all(ieee_is_finite(Sm))
          else
            okm = .false.
          end if

          if (okp .and. okm) then
            dSdU(:,m) = (Sp - Sm) / (2.0_wp*h)
            exit
          else if (okp) then
            dSdU(:,m) = (Sp - S0) / h
            exit
          else if (okm) then
            dSdU(:,m) = (S0 - Sm) / h
            exit
          else
            ! Retry with smaller h
            h = 0.1_wp * h
          end if
        end do

        if (.not. all(ieee_is_finite(dSdU(:,m)))) dSdU(:,m) = 0.0_wp
      end do

    end subroutine compute_source_jacobian_fd
  end subroutine process_cell


  subroutine print_implicit_banner(cfg)
    type(sim_config_t), intent(in) :: cfg
    integer :: nth

    if (banner_printed) return
    banner_printed = .true.

    nth = omp_get_max_threads()

    write(*,'(a)') ' _______________________________________________________________________________'
    write(*,'(a)') '|                                                                               |'
    write(*,'(a)') '|   PRIMS-F: PIPELINE RISER IMPLICIT MULTIPHASE SOLVER - FORTRAN                |'
    write(*,'(a)') '|   Implicit Time Integration: Backward Euler / Newton-Raphson                  |'
    write(*,'(a)') '|_______________________________________________________________________________|'
    write(*,'(a)') ''
    write(*,'(a,i0)') ' > OpenMP Parallel Threads:             ', nth
    write(*,'(a)') ''
    write(*,'(a)') ' INITIAL CONDITIONS & MODEL PARAMETERS:'
    write(*,'(a)') ' ------------------------------------------------------------'
    write(*,'(a,i0)') '  Number of Nodes:                             ', cfg%nCells
    write(*,'(a,1pe12.4)') '  Inlet Pressure (Ps) [nd]:               ', cfg%nPs
    write(*,'(a,1pe12.4)') '  Outlet Pressure (Pb) [nd]:              ', cfg%nPb
    write(*,'(a,1pe12.4)') '  Liquid superficial vel (jl) [nd]:       ', cfg%njl
    write(*,'(a,1pe12.4)') '  Gas superficial vel (jg) [nd]:          ', cfg%njg
    write(*,'(a,1pe12.4)') '  Termination Time [nd]:                  ', cfg%T
    write(*,'(a)') ' ------------------------------------------------------------'
    write(*,'(a)') ''

    write(*,'(a)') ' Iteration    Time       Time Step   Newton   Residual RMS    Increment RMS'
    write(*,'(a)') ' ____________________________________________________________________________'

  end subroutine print_implicit_banner


  subroutine print_step_line(step, time, dt, newt, resn, incn, status)
    integer(i4), intent(in) :: step
    real(wp), intent(in) :: time, dt
    integer(i4), intent(in) :: newt
    real(wp), intent(in) :: resn, incn
    character(len=*), intent(in) :: status

    write(*,'(i9,1x,1pe11.4,1x,1pe11.4,1x,i6,1x,1pe13.6,1x,1pe13.6,2x,a)') step, time, dt, newt, resn, incn, trim(status)
  end subroutine print_step_line


  subroutine print_newton_iter(newt, res_rms, inc_rms)
    integer(i4), intent(in) :: newt
    real(wp),    intent(in) :: res_rms, inc_rms

    write(*,'(A, I0, A, 2(A,ES12.4))') '           --- Newton Iteration ', newt, ' ---', &
                                      '  Resid_rms: ', res_rms, '  Incr_rms: ', inc_rms
  end subroutine print_newton_iter


    subroutine simulate(this, U, vns, thetav, T, max_iter, &
                       it_target, use_it_target, it_p, dt_fac_min, dt_fac_max, &
                       newton_fast, newton_slow, dt_growth, dt_shrink)
      class(pipeline_riser_simulator_t), intent(inout) :: this
      real(wp), intent(inout) :: U(:,:)
      real(wp), intent(in)    :: vns(:)
      real(wp), intent(in)    :: thetav(:)
      real(wp), intent(in)    :: T
      integer(i4), intent(in), optional :: max_iter
      integer(i4), intent(in), optional :: it_target, newton_fast, newton_slow
      logical,     intent(in), optional :: use_it_target
      real(wp),    intent(in), optional :: it_p, dt_fac_min, dt_fac_max
      real(wp),    intent(in), optional :: dt_growth, dt_shrink

      integer(i4) :: N, step, ntime
      integer(i4) :: inlet_idx, outlet_idx, first_face, last_face
      integer(i4) :: i_first, i_last
      integer(i4) :: nblk
      integer(i4) :: i, j, k, jj
      integer(i4) :: newt, newt_max, info
      integer(i4) :: retries
      logical     :: converged
      real(wp)    :: time, CFL
      real(wp)    :: dt_step, dt_cfl
      real(wp)    :: res_norm, inc_norm, res_old, res_base
      integer(i4) :: ls
      integer(i4), parameter :: ls_max = 10
      real(wp) :: lambda, res_trial, lambda_pos
      real(wp), parameter :: ls_shrink = 0.5_wp
      real(wp), parameter :: armijo_c = 1.0e-6_wp
      real(wp), parameter :: mass_eps = 1.0e-12_wp
      logical :: ls_ok
      real(wp)    :: res_rms, inc_rms
      integer(i4) :: ndof_total
      real(wp)    :: flux(3), src(3), dt
      real(wp)    :: dFdL(3,3), dFdR(3,3), dSdU(3,3)

      ! --- dt-control (can be overridden by optional args from main)
      logical     :: use_it_target_loc
      logical     :: retried_with_fd
      integer(i4) :: it_target_loc, newton_fast_loc, newton_slow_loc
      real(wp)    :: it_p_loc, dt_fac_min_loc, dt_fac_max_loc
      real(wp)    :: dt_growth_loc, dt_shrink_loc
      real(wp)    :: U_rms, dt_fac

      real(wp), allocatable :: delta_s(:)
      real(wp), allocatable :: flux_faces(:,:), source_faces(:,:), dt_faces(:)
      real(wp), allocatable :: dFdUL_faces(:,:,:), dFdUR_faces(:,:,:)
      real(wp), allocatable :: dSdU_cells(:,:,:)
      real(wp), allocatable :: U_n(:,:), U_guess(:,:), U_trial(:,:), R(:,:)
      real(wp), allocatable :: Ablk(:,:,:), Bblk(:,:,:), Cblk(:,:,:)
      real(wp), allocatable :: rhs(:,:), dU(:,:)

      real(wp), parameter :: I3(3,3) = reshape([ &
        1.0_wp,0.0_wp,0.0_wp, &
        0.0_wp,1.0_wp,0.0_wp, &
        0.0_wp,0.0_wp,1.0_wp], [3,3])
      N = size(U,1)

      ! Defensive: guarantee prev-arrays match the actual state size.
      call this%ensure_prev_arrays(N)

      ! --------------------------------------------------------------
      ! Indices (1-based)
      !   - physical cells are [nd+1 .. N-nd]
      !   - residual needs flux_faces(inlet_idx-1) .. flux_faces(outlet_idx)
      ! --------------------------------------------------------------
      inlet_idx  = first_physical_cell(this%cfg)
      outlet_idx = last_physical_cell(this%cfg)

      ! Physical unknowns solved implicitly:
      !   - If inlet is strong Dirichlet (non-OSC), exclude the inlet physical cell.
      !   - Otherwise include it.
      if (index(to_lower(this%cfg%regime), 'osc') > 0) then
        i_first = inlet_idx
      else
        i_first = inlet_idx + 1
      end if
      i_last  = outlet_idx

      ! Faces needed for residual/Jacobian
      first_face = inlet_idx - 1
      last_face  = outlet_idx
      ! Number of interior unknowns solved implicitly (cells i_first .. i_last)
      nblk = i_last - i_first + 1
      ndof_total = 3*nblk
      if (ndof_total < 1) ndof_total = 1

      allocate(delta_s(N))
      allocate(flux_faces(N,3), source_faces(N,3), dt_faces(N))
      allocate(dFdUL_faces(N,3,3), dFdUR_faces(N,3,3))
      allocate(dSdU_cells(N,3,3))
      allocate(U_n(N,3), U_guess(N,3), U_trial(N,3), R(N,3))
      allocate(Ablk(3,3,nblk), Bblk(3,3,nblk), Cblk(3,3,nblk))
      allocate(rhs(3,nblk), dU(3,nblk))

      ! cell sizes (faces), in nondimensional s
      delta_s = 0.0_wp
      do i = 1, N-1
        delta_s(i) = max(abs(vns(i+1) - vns(i)), 1.0e-12_wp)
      end do
      delta_s(N) = delta_s(N-1)

      time = 0.0_wp
      step = 0

      ! max number of time steps (safety). Normally max_iter is provided from main.
      ntime = int(T / max(this%min_dt, 1.0e-16_wp)) + 1
      if (present(max_iter)) ntime = max_iter

      newt_max = max(1_i4, this%cfg%newton_max)

      ! Defaults for adaptive dt logic (can be overridden by optional args)
      use_it_target_loc = .false.
      it_target_loc     = 6
      it_p_loc          = 0.7_wp
      dt_fac_min_loc    = 0.5_wp
      dt_fac_max_loc    = 2.0_wp
      newton_fast_loc   = this%newton_fast
      newton_slow_loc   = this%newton_slow
      dt_growth_loc     = this%dt_growth
      dt_shrink_loc     = this%dt_shrink

      if (present(use_it_target)) use_it_target_loc = use_it_target
      if (present(it_target))     it_target_loc     = it_target
      if (present(it_p))          it_p_loc          = it_p
      if (present(dt_fac_min))    dt_fac_min_loc    = dt_fac_min
      if (present(dt_fac_max))    dt_fac_max_loc    = dt_fac_max
      if (present(newton_fast))   newton_fast_loc   = newton_fast
      if (present(newton_slow))   newton_slow_loc   = newton_slow
      if (present(dt_growth))     dt_growth_loc     = dt_growth
      if (present(dt_shrink))     dt_shrink_loc     = dt_shrink

      it_target_loc  = max(1, it_target_loc)
      it_p_loc       = max(0.1_wp, it_p_loc)
      dt_fac_min_loc = max(0.05_wp, dt_fac_min_loc)
      dt_fac_max_loc = max(dt_fac_min_loc, dt_fac_max_loc)
      dt_growth_loc  = max(1.0_wp, dt_growth_loc)
      dt_shrink_loc  = min(0.99_wp, max(0.05_wp, dt_shrink_loc))

      ! start with a reasonable explicit-like seed, but allow backtracking below min_dt if needed
          ! initial dt: default to 1% of total non-dimensional time (can be overridden by cfg%dt_start)
    if (this%cfg%dt_start > 0.0_wp) then
      dt_step = this%cfg%dt_start
    else
      dt_step = 0.01_wp * T
    end if
    dt_step = min(max(dt_step, this%min_dt), this%max_dt)

      call print_implicit_banner(this%cfg)

       call this%update_prev_arrays(U, thetav)

      ! Write initial snapshot (step 0)
      if (this%cfg%save_every > 0) then
        call write_output_snapshot(this, U, vns, thetav, step, time)
      end if

      do while (time < T .and. step < ntime)

        U_n = U
        retries = 0
        retried_with_fd = .false.
        attempt: do
    CFL = this%cfg%CFL
    if (CFL <= 0.0_wp) CFL = this%cfg%cfl_target
          ! ------------------------------------------------------------------
          ! CFL estimate (now based on MATERIAL speed inside process_cell)
          ! ------------------------------------------------------------------
          dt_faces = this%max_dt
          !$omp parallel do default(none) private(j,flux,src,dt,dFdL,dFdR,dSdU) &
          !$omp shared(this,U_n,delta_s,thetav,CFL,dt_faces,first_face,inlet_idx,last_face)
          do j = first_face, last_face
            call this%process_cell(j, U_n, delta_s, thetav, CFL, flux, src, dt, dFdL, dFdR)
            dt_faces(j) = dt
          end do
          !$omp end parallel do

          ! NOTE:
          !   dt_faces can be artificially tiny in ghost regions because vns/delta_s
          !   may be duplicated there (ds ~ 0 -> clamped). Do NOT let those ghost
          !   faces limit the global dt in an implicit solve.
          !   Use only faces adjacent to the unknown/physical cells.
          dt_cfl = minval(dt_faces(inlet_idx:outlet_idx-1))
          if (.not. ieee_is_finite(dt_cfl)) dt_cfl = this%max_dt

          ! Enforce bounds and remaining time
          dt_step = max(dt_step, this%min_dt)
          dt_step = min(dt_step, this%max_dt)

          ! ------------------------------------------------------------------
          ! IMPORTANT (implicit BE):
          !   Do NOT hard-limit dt_step by a CFL estimate.
          !   A too-small dt_cfl (e.g. from any pathological face) can freeze the
          !   whole simulation at dt_min and trigger premature backtracking failure.
          !
          !   Here we use dt_cfl only as a *soft* cap when it looks meaningful.
          ! ------------------------------------------------------------------
          if (ieee_is_finite(dt_cfl)) then
            if (dt_cfl > 10.0_wp*this%min_dt) then
              dt_step = min(dt_step, max(1.0_wp, this%cfg%cfl_target) * dt_cfl)
            end if
          end if
          if (time + dt_step > T) dt_step = max(T - time, 0.0_wp)
          if (dt_step <= 0.0_wp) exit

          ! ------------------------------------------------------------------
          ! Newton-Raphson for Backward Euler step
          ! ------------------------------------------------------------------
          U_guess = U_n
          call apply_boundary_conditions(U_guess, inlet_idx, outlet_idx, thetav(inlet_idx), thetav(outlet_idx), this%cfg)
          call this%update_prev_arrays(U_guess, thetav)
          converged = .false.

          ! Residual/increment norms (RMS per DOF recommended)
          res_norm = 0.0_wp
          inc_norm = 0.0_wp
          res_rms  = 0.0_wp
          inc_rms  = 0.0_wp

          ! print line at beginning of step (LS-DYNA style)
          call print_step_line(step, time, dt_step, 0, res_rms, 1.0_wp, 'converging')

          res_old  = huge(1.0_wp)
          inc_rms  = 0.0_wp
          do newt = 1, newt_max

            ! NOTE: ul_prev/ug_prev are frozen during Newton (updated only after accepting a step)

            ! Fluxes, sources and Roe-derivatives at faces using current U_guess
            flux_faces   = 0.0_wp
            source_faces = 0.0_wp
            dFdUL_faces  = 0.0_wp
            dFdUR_faces  = 0.0_wp
            dSdU_cells  = 0.0_wp
            !$omp parallel do default(none) private(j,flux,src,dt,dFdL,dFdR,dSdU) &
            !$omp shared(this,U_guess,delta_s,thetav,CFL,flux_faces,source_faces,dFdUL_faces,dFdUR_faces,dSdU_cells,first_face,inlet_idx,last_face)
            do j = first_face, last_face
              call this%process_cell(j, U_guess, delta_s, thetav, CFL, flux, src, dt, dFdL, dFdR, dSdU)
              flux_faces(j,:)   = flux(:)
              source_faces(j,:) = src(:)
              dFdUL_faces(j,:,:) = dFdL(:,:)
              dFdUR_faces(j,:,:) = dFdR(:,:)
              dSdU_cells(j,:,:)  = dSdU(:,:)
            end do
            !$omp end parallel do

            call this%ti%compute_residual_vector(U_guess, U_n, delta_s, flux_faces, source_faces, i_first, i_last, dt_step, R, res_norm)

            ! RMS per DOF (recommended): ||x||_2 / sqrt(ndof)
            res_rms = res_norm / sqrt(real(ndof_total, wp))
            ! Keep a copy of current residual RMS for line-search
            res_base = res_rms

            ! Do not hard-abort on non-monotone residuals; a line-search will handle it.
            ! (Nonlinear problems often have occasional increases even when converging.)
            ! Mixed convergence criterion (recommended)
            if (res_rms <= this%cfg%tol_res) then
              converged = .true.
              exit
            end if

            ! ----------------------------------------------------------------
                        ! Build block-tridiagonal Jacobian for implicit unknowns (cells i_first .. i_last)
                        Ablk = 0.0_wp ; Bblk = 0.0_wp ; Cblk = 0.0_wp ; rhs = 0.0_wp
                        do j = i_first, i_last
                          k = j - i_first + 1
                          rhs(:,k) = -R(j,:)

                          ! Base diagonal: Identity
                          Bblk(:,:,k) = 0.0_wp
                          Bblk(1,1,k) = 1.0_wp
                          Bblk(2,2,k) = 1.0_wp
                          Bblk(3,3,k) = 1.0_wp

                          ! Coupling to left unknown cell (U_{j-1}) only if j>i_first
                          if (j > i_first) then
                            Ablk(:,:,k) = -(dt_step / delta_s(j)) * dFdUL_faces(j-1,:,:)
                          else
                            Ablk(:,:,k) = 0.0_wp
                          end if

                          ! Coupling to right unknown cell (U_{j+1}) only if j<i_last
                          if (j < i_last) then
                            Cblk(:,:,k) =  (dt_step / delta_s(j)) * dFdUR_faces(j,:,:)
                          else
                            Cblk(:,:,k) = 0.0_wp
                          end if

                          ! Flux contributions to diagonal (always include both faces; ghost/fixed neighbors are treated as constants)
                          Bblk(:,:,k) = Bblk(:,:,k) + (dt_step / delta_s(j)) * ( dFdUL_faces(j,:,:) - dFdUR_faces(j-1,:,:) )

                          ! Source Jacobian contribution (SIGN: -dt*dS/dU  (residual: U-U_n - dt*rhs))
                          Bblk(:,:,k) = Bblk(:,:,k) - dt_step * dSdU_cells(j,:,:)

                        end do

            call this%ti%solve_block_tridiagonal_3x3(Ablk, Bblk, Cblk, rhs, dU, info)
            if (info /= 0) then
              write(*,'(a,i0)') ' [ERROR] Block Thomas failed, info=', info
              converged = .false.
              exit
            end if

            inc_norm = sqrt(sum(dU**2))
            inc_rms = inc_norm / sqrt(real(ndof_total, wp))
            ! Relative increment criterion: inc_rms <= tol_inc * max(1, U_rms)
            U_rms = sqrt(sum(U_guess(i_first:i_last,:)**2) / real(ndof_total, wp))
            if (res_rms <= this%cfg%tol_res .and. inc_rms <= this%cfg%tol_inc * max(1.0_wp, U_rms)) then
              converged = .true.
              exit
            end if

            write(*,'(A, I0, A, 2(A,ES12.4))') '           --- Newton Iteration ', newt, ' ---', '  Resid_rms: ', res_rms, '  Incr_rms: ', inc_rms
            ! --- Line-search / damping on Newton step to preserve admissibility and reduce residual
            ls_ok = .false.
            lambda = 1.0_wp
    ! Start from a positivity-safe lambda (fraction-to-the-boundary for u1,u2)
            lambda_pos = 1.0_wp
            do jj = 1, nblk
              j = i_first + jj - 1
              if (dU(1,jj) < 0.0_wp) lambda_pos = min(lambda_pos, 0.99_wp * (U_guess(j,1) - mass_eps) / (-dU(1,jj)))
              if (dU(2,jj) < 0.0_wp) lambda_pos = min(lambda_pos, 0.99_wp * (U_guess(j,2) - mass_eps) / (-dU(2,jj)))
            end do
            if (.not. ieee_is_finite(lambda_pos)) lambda_pos = 1.0_wp
            lambda_pos = max(min(lambda_pos, 1.0_wp), 1.0e-12_wp)
            lambda = min(lambda, lambda_pos)

            do ls = 1, ls_max
              U_trial = U_guess
              U_trial(i_first:i_last,:) = U_guess(i_first:i_last,:) + lambda * transpose(dU)

              ! Quick admissibility: keep u1,u2 strictly positive
              if (minval(U_trial(i_first:i_last,1)) <= mass_eps .or. minval(U_trial(i_first:i_last,2)) <= mass_eps) then
                lambda = lambda * ls_shrink
                cycle
              end if

              call apply_boundary_conditions(U_trial, inlet_idx, outlet_idx, thetav(inlet_idx), thetav(outlet_idx), this%cfg)

              ! Recompute residual at trial state (no Jacobians needed here)
              flux_faces   = 0.0_wp
              source_faces = 0.0_wp
              !$omp parallel do default(none) private(j,flux,src,dt) &
              !$omp shared(this,U_trial,delta_s,thetav,CFL,flux_faces,source_faces,first_face,last_face)
              do j = first_face, last_face
                call this%process_cell(j, U_trial, delta_s, thetav, CFL, flux, src, dt)
                flux_faces(j,:)   = flux(:)
                source_faces(j,:) = src(:)
              end do
              !$omp end parallel do

              call this%ti%compute_residual_vector(U_trial, U_n, delta_s, flux_faces, source_faces, i_first, i_last, dt_step, R, res_norm)
              res_trial = (res_norm / sqrt(real(ndof_total, wp)))

              if (ieee_is_finite(res_trial)) then
                if (res_trial <= (1.0_wp - armijo_c*lambda) * res_base .or. res_trial < res_base) then
                  ls_ok = .true.
                  exit
                end if
              end if

              lambda = lambda * ls_shrink
            end do

            if (.not. ls_ok) then
              converged = .false.
              exit
            end if

            U_guess = U_trial
            res_old = res_trial
          end do

          if (converged) then
            U = U_guess
            call print_step_line(step, time, dt_step, newt, res_rms, inc_rms, 'CONVERGED')
            time = time + dt_step
            step = step + 1



            ! Save output snapshot every cfg%save_every accepted steps
            if (this%cfg%save_every > 0) then
              if (mod(step, this%cfg%save_every) == 0) then
                call write_output_snapshot(this, U, vns, thetav, step, time)
              end if
            end if
            call this%update_prev_arrays(U, thetav)

            ! Adaptive dt control
            if (use_it_target_loc) then
              dt_fac = (real(it_target_loc, wp) / real(max(newt,1), wp))**it_p_loc
              dt_fac = min(max(dt_fac, dt_fac_min_loc), dt_fac_max_loc)
              dt_step = min(max(dt_step * dt_fac, this%min_dt), this%max_dt)
            else
              if (newt <= newton_fast_loc) then
                dt_step = min(dt_step * dt_growth_loc, this%max_dt)
              else if (newt >= newton_slow_loc) then
                dt_step = max(dt_step * dt_shrink_loc, this%min_dt)
              end if
            end if

            exit attempt

          else
            call print_step_line(step, time, dt_step, newt, res_rms, inc_rms, 'FAILED')

            ! If analytic/frozen-dissipation Jacobian struggles, retry with full FD Jacobian of numerical flux
            if (.not. this%cfg%jac_dissipation_fd .and. .not. retried_with_fd) then
              retried_with_fd = .true.
              this%cfg%jac_dissipation_fd = .true.
              write(*,'(a)') '  [JACOBIAN] Newton failed: retrying with jac_dissipation_fd=.true. (full FD Jacobian of numerical flux).'
              cycle attempt
            end if

            retries = retries + 1
            dt_step = dt_step * dt_shrink_loc
            write(*,'(a,i0,a,1pe12.4)') '  [BACKTRACK] tentativa ', retries, ' -> novo dt = ', dt_step

            if (dt_step < this%min_dt) error stop 'Erro: dt muito pequeno e sem convergencia.'
            if (retries >= this%max_backtrack) error stop 'Erro: muitas tentativas de backtracking.'

            cycle attempt
          end if
        end do attempt

      end do

      ! Ensure the final state is also written to the VisIt series
      if (this%cfg%save_every > 0) then
        if (step == 0 .or. mod(step, this%cfg%save_every) /= 0) then
          call write_output_snapshot(this, U, vns, thetav, step, time)
        end if
      end if

      deallocate(delta_s, flux_faces, source_faces, dt_faces)
      deallocate(dFdUL_faces, dFdUR_faces, dSdU_cells)
      deallocate(U_n, U_guess, U_trial, R)
      deallocate(Ablk, Bblk, Cblk, rhs, dU)

    end subroutine simulate

  !> Writes a snapshot of the solution every N accepted time steps.

  !> Writes a snapshot of the solution every N accepted time steps (VTK PolyData for VisIt).
  subroutine write_output_snapshot(this, U, vns, thetav, step, time)
    class(pipeline_riser_simulator_t), intent(in) :: this
    real(wp), intent(in) :: U(:,:), vns(:), thetav(:)
    integer(i4), intent(in) :: step
    real(wp), intent(in) :: time

    integer(i4) :: i, k, inlet_idx, outlet_idx, npts, idx
    real(wp)    :: alpha, rhol, rhog, P, ul, ug
    real(wp)    :: ds, thm, s_dim
    character(len=256) :: fname, fname_base, fvisit
    integer :: iu

    real(wp), allocatable :: x(:), y(:), z(:), s(:)
    real(wp), allocatable :: alpha_v(:), rhol_v(:), rhog_v(:), P_v(:), ul_v(:), ug_v(:), theta_v(:)
    real(wp), allocatable :: u1_v(:), u2_v(:), u3_v(:)

    inlet_idx  = first_physical_cell(this%cfg)
    outlet_idx = last_physical_cell(this%cfg)
    npts = outlet_idx - inlet_idx + 1
    if (npts <= 1) return

    allocate(x(npts), y(npts), z(npts), s(npts))
    allocate(alpha_v(npts), rhol_v(npts), rhog_v(npts), P_v(npts), ul_v(npts), ug_v(npts), theta_v(npts))
    allocate(u1_v(npts), u2_v(npts), u3_v(npts))

    ! Geometry reconstruction (x,z) from theta(s): dx/ds = cos(theta), dz/ds = sin(theta)
    x(1) = 0.0_wp
    y(1) = 0.0_wp
    z(1) = 0.0_wp

    do k = 1, npts
      i = inlet_idx + k - 1

      call compute_primitive_variables( &
        U(i,1), U(i,2), U(i,3), thetav(i), &
        this%cfg%nCg, this%cfg%nCl, &
        this%cfg%nrho_l0, this%cfg%nP_l0, &
        this%cfg%D_H, this%cfg%AREA, this%cfg%eps, this%cfg%G, &
        this%cfg%mu_l, this%cfg%mu_g, this%cfg%sigma, this%cfg%njl, &
        this%cfg%omega_u, this%cfg%omega_rho, &
        this%cfg%tol, this%cfg%tola, this%maxit_root, &
        alpha, rhol, rhog, P, ul, ug, idx, &
        idx_prev=this%idx_prev(i), ul_prev=this%ul_prev(i), ug_prev=this%ug_prev(i) )

      if (idx < 0_i4 .or. (.not. ieee_is_finite(alpha)) .or. (.not. ieee_is_finite(P)) ) then
        alpha = max(this%alpha_eps, min(1.0_wp - this%alpha_eps, U(i,1)))
        rhol  = this%cfg%nrho_l0
        rhog  = max(1.0e-12_wp, this%cfg%nPb/(this%cfg%nCg*this%cfg%nCg))
        P     = this%cfg%nP_l0
        ul    = 0.0_wp
        ug    = 0.0_wp
      end if

      s_dim = vns(i) * this%cfg%Lr
      s(k)  = s_dim

      if (k > 1) then
        ds  = (vns(i) - vns(i-1)) * this%cfg%Lr
        thm = 0.5_wp * (thetav(i) + thetav(i-1))
        x(k) = x(k-1) + ds * cos(thm)
        y(k) = 0.0_wp
        z(k) = z(k-1) + ds * sin(thm)
      end if

      alpha_v(k) = alpha
      rhol_v(k)  = rhol * this%cfg%omega_rho
      rhog_v(k)  = rhog * this%cfg%omega_rho
      P_v(k)     = P    * this%cfg%omega_P
      ul_v(k)    = ul   * this%cfg%omega_u
      ug_v(k)    = ug   * this%cfg%omega_u
      theta_v(k) = thetav(i)

      u1_v(k) = U(i,1)
      u2_v(k) = U(i,2)
      u3_v(k) = U(i,3)
    end do

    write(fname_base,'(a,a,a,i6.6,a)') 'results_', trim(this%cfg%regime), '_it', step, '.vtk'
    fname = trim(this%cfg%out_dir)//'/'//trim(fname_base)

    open(newunit=iu, file=trim(fname), status='replace', action='write', form='formatted')

    write(iu,'(a)') '# vtk DataFile Version 2.0'
    write(iu,'(a)') 'PRIMS-F pipeline-riser snapshot'
    write(iu,'(a)') 'ASCII'
    write(iu,'(a)') 'DATASET POLYDATA'
    write(iu,'(a,i0,a)') 'POINTS ', npts, ' double'
    do k = 1, npts
      write(iu,'(3(1x,es23.15))') x(k), y(k), z(k)
    end do

    write(iu,'(a,i0,1x,i0)') 'LINES ', 1, npts+1
    write(iu,'(i0)', advance='no') npts
    do k = 0, npts-1
      write(iu,'(1x,i0)', advance='no') k
    end do
    write(iu,*)

    write(iu,'(a,i0)') 'POINT_DATA ', npts
    call write_scalar(iu, 's',    s)
    call write_scalar(iu, 'alpha',  alpha_v)
    call write_scalar(iu, 'rho_l',  rhol_v)
    call write_scalar(iu, 'rho_g',  rhog_v)
    call write_scalar(iu, 'P',   P_v)
    call write_scalar(iu, 'u_l',    ul_v)
    call write_scalar(iu, 'u_g',    ug_v)
    call write_scalar(iu, 'theta',  theta_v)
    call write_scalar(iu, 'u1',     u1_v)
    call write_scalar(iu, 'u2',     u2_v)
    call write_scalar(iu, 'u3',     u3_v)

    close(iu)

    if (this%cfg%save_every > 0) then
      fvisit = trim(this%cfg%out_dir)//'/results_'//trim(this%cfg%regime)//'.visit'
      open(newunit=iu, file=trim(fvisit), status='old', action='write', position='append', form='formatted')
      write(iu,'(a)') trim(fname_base)
      close(iu)
    end if

    write(*,'(a,1x,a,1x,a,i0,a,1x,es12.4)') ' [OUTPUT] wrote', trim(fname), 'at step', step, ', time=', time

    deallocate(x,y,z,s, alpha_v,rhol_v,rhog_v,P_v,ul_v,ug_v,theta_v, u1_v,u2_v,u3_v)

  contains

    subroutine write_scalar(iu, name, arr)
      integer, intent(in) :: iu
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: arr(:)
      integer :: kk
      write(iu,'(a,1x,a,1x,a)') 'SCALARS', trim(name), 'double 1'
      write(iu,'(a)') 'LOOKUP_TABLE default'
      do kk = 1, size(arr)
        if (.not. ieee_is_finite(arr(kk))) then
          write(iu,'(es23.15)') 0.0_wp
        else
          write(iu,'(es23.15)') arr(kk)
        end if
      end do
    end subroutine write_scalar

  end subroutine write_output_snapshot

end module module_simulate_pipe_riser