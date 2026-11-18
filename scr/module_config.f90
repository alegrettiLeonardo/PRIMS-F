module module_config
  !! Shared configuration type replacing the Python dict `config`.
  !!
  !! Indexing/ghost cells:
  !!   total cells include ghosts: N = cfg%nCells
  !!   left ghosts:   1..nd
  !!   physical:      nd+1..N-nd
  !!   right ghosts:  N-nd+1..N
  use mod_kinds, only: wp, i4
  implicit none
  private

  integer, parameter, public :: STR_LEN = 64
  integer, parameter, public :: PATH_LEN = 256

  type, public :: sim_config_t
    ! --- tolerances
    real(wp) :: tol      = 1.0e-15_wp   ! legacy (root-finding etc.)
    real(wp) :: tola     = 1.0e-13_wp   ! legacy (steady-state etc.)
    real(wp) :: tol_res  = 1.0e-15_wp   ! Newton residual RMS tolerance
    real(wp) :: tol_inc  = 1.0e-15_wp   ! Newton increment RMS tolerance

    ! --- nonlinear solve controls (Python: maxit)
    integer(i4) :: maxit = 1000
    integer(i4) :: newton_max = 10  ! max Newton iterations per time step

    ! --- scenario/IO
    character(len=STR_LEN) :: regime = 'default'
    character(len=PATH_LEN) :: out_dir = 'output'
    integer(i4) :: max_workers = 1
    integer(i4) :: save_every = 100

    ! --- mesh
    integer(i4) :: num_dummy_cells = 2
    integer(i4) :: nCells = 0
    integer(i4) :: N_pipe = 0          ! number of pipe cells in the extended vector (Python: nodes-1)
    integer(i4) :: closure_pipe = 1
    integer(i4) :: closure_riser = 2
    character(len=16) :: flux_scheme = 'ROE'
    character(len=16) :: gas_eos = 'ISOTHERMAL'
    character(len=16) :: junction_model = 'MATCHED'
    character(len=20) :: dissipation_mode = 'LAGGED'
    real(wp) :: output_interval_s = 1.0_wp
    real(wp) :: field_interval_s = 30.0_wp
    real(wp) :: perturbation_amplitude = 0.01_wp

    ! --- geometry/physics
    real(wp) :: Lr = 0.0_wp
    real(wp) :: Lp = 0.0_wp
    real(wp) :: D_H = 0.0_wp
    real(wp) :: AREA = 0.0_wp
    real(wp) :: eps = 0.0_wp
    real(wp) :: G = 0.0_wp
    real(wp) :: beta = 0.0_wp

    ! --- dimensional material props used in closure/source
    real(wp) :: mu_g = 0.0_wp
    real(wp) :: mu_l = 0.0_wp
    real(wp) :: sigma = 0.0_wp

    ! --- nondimensional sound speeds and reference states
    real(wp) :: nCg = 0.0_wp
    real(wp) :: nCl = 0.0_wp
    real(wp) :: nrho_l0 = 0.0_wp
    real(wp) :: nP_l0   = 0.0_wp
    real(wp) :: nPs     = 0.0_wp
    real(wp) :: nPb     = 0.0_wp

    ! --- nondimensional inlet mass flows
    real(wp) :: mul = 0.0_wp
    real(wp) :: mug = 0.0_wp

    ! --- nondimensional superficial velocities (for logging / BC reference)
    ! These are returned by adimensionalizar() in module_solver_utils.
    real(wp) :: njl = 0.0_wp
    real(wp) :: njg = 0.0_wp

    ! --- scaling factors (Python: omega_*)
    real(wp) :: omega_u   = 1.0_wp
    real(wp) :: omega_rho = 1.0_wp
    real(wp) :: omega_P   = 1.0_wp
    real(wp) :: omega_t = 1.0_wp     ! time scale [s] (dim -> nd: t_nd = t_dim/omega_t)

    ! --- time integration
    real(wp) :: CFL = 0.95_wp
    real(wp) :: T   = 0.0_wp
    real(wp) :: dt_min    = 1.0e-9_wp   ! hard floor (machine/robustness)
    real(wp) :: dt_max    = -1.0_wp     ! computed from T (suggested: T/10)
    real(wp) :: dt_start  = -1.0_wp     ! initial guess (suggested: 0.01*T)
    real(wp) :: cfl_target = 5.0_wp     ! default implicit CFL if CFL<=0
    integer(i4) :: iter_cfl_switch = 900000

    ! --- Jacobian quality / Newton robustness
    logical :: jac_dissipation_fd = .false.   ! include FD d|A|/dU * dU correction in Roe flux Jacobian
    real(wp) :: jac_dissipation_fd_eps = 1.0e-6_wp
    logical :: debug_dfduj_check = .false.
    integer(i4) :: debug_dfduj_cell = 0
    real(wp) :: debug_dfduj_eps = 1.0e-6_wp


    ! --- Roe/quad
    integer(i4) :: quad_points = 3
    integer(i4) :: interface_buffer_cells = 2


    ! --- reconstruction control
    logical :: use_muscl = .false.
    ! --- boundary condition safety
    real(wp) :: bc_min_value = 1.0e-12_wp

    ! --- initial conditions (extended with ghosts)
    real(wp), allocatable :: U0(:,:)     ! (N,3)
    real(wp), allocatable :: vns(:)      ! (N)
    real(wp), allocatable :: thetav(:)   ! (N)
    real(wp), allocatable :: cell_width(:) ! (N), dimensionless TRUE FV widths
    real(wp), allocatable :: face_positions(:) ! physical faces, 0..nphysical
  end type sim_config_t

  public :: first_physical_cell, last_physical_cell

contains

  pure integer(i4) function first_physical_cell(cfg) result(k)
    type(sim_config_t), intent(in) :: cfg
    k = cfg%num_dummy_cells + 1
  end function first_physical_cell

  pure integer(i4) function last_physical_cell(cfg) result(k)
    type(sim_config_t), intent(in) :: cfg
    k = cfg%nCells - cfg%num_dummy_cells
  end function last_physical_cell

end module module_config
