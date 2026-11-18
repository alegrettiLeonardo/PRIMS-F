program main
  !! Driver similar to Python main.py, using module_steady_state for U0.

  use mod_kinds,               only: wp, i4
  use module_config,           only: sim_config_t
  use module_solver_utils,     only: fluid_dynamics_t
  use module_catenary,         only: solve_A
  use module_steady_state,     only: EstadoEstacionario_ndim_simp
  use module_boundary_conditions, only: apply_ghost_cells
  use module_find,             only: compute_primitive_variables, set_stability_params
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  ! ------------------------- Namelist parameters -------------------------
  integer(i4)        :: nodes, num_dummy_cells
  real(wp)           :: X, Z, Lp
  real(wp)           :: beta_deg
  real(wp)           :: rho_l0, P_l0
  real(wp)           :: Cl, Rg
  real(wp)           :: mu_l, mu_g
  ! Inlet mass fluxes (per area) [kg/(m^2.s)]
  ! NOTE: In the original Python code:
  !   - mul, mug  = inlet *mass fluxes* (rho*j)
  !   - MUL, MUG  = viscosities
  ! Fortran is case-insensitive, so viscosities are named mu_l/mu_g here.
  ! If mul/mug are negative, liquid mass flux uses rho_l0*jl and gas mass
  ! flux uses the explicitly selected STANDARD or SEPARATOR reference.
  real(wp)           :: mul, mug
  real(wp)           :: eps, D, sigma
  real(wp)           :: CFL, tol, tola, tol_res, tol_inc, G
  ! --- implicit time-step controls (non-dimensional) ---
  real(wp)           :: dt_min, dt_max, dt_start, cfl_target
  ! --- adaptive dt control based on Newton iterations ---
  integer(i4)        :: it_target, newton_fast, newton_slow
  logical            :: use_it_target
  real(wp)           :: it_p, dt_fac_min, dt_fac_max, dt_growth, dt_shrink
  ! --- Jacobian quality / debug ---
  logical            :: jac_dissipation_fd, debug_dfduj_check
  real(wp)           :: jac_dissipation_fd_eps, debug_dfduj_eps
  integer(i4)        :: debug_dfduj_cell

  real(wp)           :: Ps, jl, jg, Tgas, time_dim
  ! Backward-compatible aliases (some inputs use different names)
  real(wp)           :: mesh_ratio
  real(wp)           :: T_scenario, time_phys
  real(wp)           :: gamma_g
  character(len=16) :: gas_eos, gas_flow_basis, flux_scheme, junction_model
  character(len=20) :: dissipation_mode
  real(wp) :: Cg_fixed, P_standard, T_standard, dt_min_s, dt_max_s, dt_start_s
  real(wp) :: perturbation_amplitude, output_interval_s, field_interval_s
  integer(i4) :: n_riser, closure_riser, sim_status, steady_status
  integer(i4)        :: quad_points, max_iter, max_workers, save_every, newton_max
  character(len=32)  :: regime

  ! (Optional) SS1 control parameters used in some Python decks.
  ! Parsed for input compatibility (may be unused depending on your build).
  real(wp)           :: pb_dp_on, pb_dp_off

  ! --- Pipeline stability switch hysteresis (CriteriaStabFlow) ---
  real(wp)           :: stab_si_on, stab_si_off
  logical            :: stab_use_hysteresis


  character(len=256) :: out_dir
  logical           :: use_muscl
  namelist /PARAMS/ nodes, num_dummy_cells, X, Z, Lp, beta_deg, D, eps, mesh_ratio, &
                    rho_l0, P_l0, Cl, Rg, mu_l, mu_g, mul, mug, sigma, G, &
                    CFL, tol, tola, tol_res, tol_inc, dt_min, dt_max, dt_start, cfl_target, &
                    it_target, use_it_target, it_p, dt_fac_min, dt_fac_max, dt_growth, dt_shrink, newton_fast, newton_slow, quad_points, Ps, jl, jg, Tgas, time_dim, T_scenario, time_phys, &
                    gamma_g, max_iter, max_workers, save_every, newton_max, regime, use_muscl, out_dir, jac_dissipation_fd, jac_dissipation_fd_eps, debug_dfduj_check, debug_dfduj_cell, debug_dfduj_eps, &
                    pb_dp_on, pb_dp_off, &
                    stab_si_on, stab_si_off, stab_use_hysteresis, &
                    gas_eos, gas_flow_basis, Cg_fixed, P_standard, T_standard, n_riser, closure_riser, &
                    dt_min_s, dt_max_s, dt_start_s, perturbation_amplitude, output_interval_s, field_interval_s, flux_scheme, junction_model, dissipation_mode
  ! ------------------------------ Locals ---------------------------------
  type(sim_config_t)            :: cfg
  type(pipeline_riser_simulator_t) :: sim
  type(fluid_dynamics_t)        :: fd

  integer(i4) :: i, nd, Nss, N, ios
  real(wp)    :: CA, Lr, AREA
  real(wp)    :: Cg
  real(wp)    :: Pb, rhog, rhol
  real(wp)    :: ml_flux, mg_flux
  real(wp)           :: ntime, omega_t
  real(wp)    :: beta_rad

  ! Adimensionalization outputs (match Python solver_utils.adimensionalizar)
  real(wp) :: nCg, nCl, nP_l0, nrho_l0
  real(wp) :: nPs, nPb, njl, njg
  real(wp) :: nmul, nmug
  real(wp) :: omega_c, omega_P, omega_rho, omega_u

  ! Steady-state vectors (dimensionless)
  real(wp), allocatable :: vns(:), vnp(:)
  real(wp), allocatable :: physical_faces(:)
  real(wp), allocatable :: nrhogv(:), nrholv(:)
  real(wp), allocatable :: alphav(:), nugv(:), nulv(:)
  real(wp), allocatable :: thetav(:), njlv(:), njgv(:)

  ! Extended arrays (with ghost cells)
  real(wp), allocatable :: U(:,:), vns_ext(:), thetav_ext(:)

  ! Output primitives
  real(wp) :: alpha_i, nrhol_i, nrhog_i, P_nd, ul_nd, ug_nd
  integer(i4) :: phase
  real(wp) :: rho_l_dim, rho_g_dim, ul_dim, ug_dim, P_dim, s_dim
  character(len=128) :: outname
  character(len=256) :: infile
  character(len=256) :: iomsg

  ! -------------------------- Defaults (safe) ----------------------------
  nodes           = 41
  num_dummy_cells = 2
  X = 6.435_wp
  Z = 9.886_wp
  Lp = 57.4_wp
  beta_deg = 2.0_wp
  rho_l0 = 997.05_wp
  P_l0   = 101325.0_wp
  Cl = 1498.0_wp
  Rg = 287.053_wp
  mu_l = 0.89e-3_wp
  mu_g = 1.85e-5_wp
  ! Optional inlet mass fluxes (kg/m^2/s). If <=0, they will be computed from jl/jg.
  mul  = -1.0_wp
  mug  = -1.0_wp
  eps  = 4.6e-5_wp
  D    = (2.54e-2_wp)*2.0_wp
  sigma = 7.28e-2_wp
  ! For implicit BE/Newton, we can use a much larger "effective CFL".
  CFL = 50.0_wp
  tol  = 1.0e-10_wp
  tola = 1.0e-8_wp
  ! If tol_res / tol_inc are not provided in the namelist, they will be set = tol after reading.
  tol_res = -1.0_wp
  tol_inc = -1.0_wp
  ! --- dt bounds (non-dimensional). If <=0 in input, defaults are used later.
  dt_min = -1.0_wp
  dt_max = -1.0_wp
  dt_start = -1.0_wp
  cfl_target = -1.0_wp

  ! --- adaptive dt control defaults ---
  use_it_target = .true.
  it_target     = 7
  it_p          = 0.7_wp
  dt_fac_min    = 0.5_wp
  dt_fac_max    = 2.0_wp
  dt_growth     = 1.2_wp
  dt_shrink     = 0.5_wp
  newton_fast   = 5
  newton_slow   = 10

  ! --- Jacobian quality / debug defaults ---
  jac_dissipation_fd      = .false.
  jac_dissipation_fd_eps  = 1.0e-6_wp
  debug_dfduj_check       = .false.
  debug_dfduj_cell        = 0
  debug_dfduj_eps         = 1.0e-6_wp

  G    = 9.81_wp
  Ps   = 2.06e5_wp
  jl   = 0.84_wp
  jg   = 0.46_wp
  Tgas = 311.95_wp
  time_dim = 20.0_wp
  mesh_ratio = 1.0_wp
  T_scenario = -1.0_wp
  time_phys  = -1.0_wp
  gamma_g = 1.4_wp

  ! Optional SS1 thresholds (parsed for compatibility; may be unused)
  pb_dp_on  = 1.05_wp
  pb_dp_off = 0.98_wp

  ! Pipeline stability (hysteresis)
  stab_si_on  = 1.0_wp
  stab_si_off = 0.9_wp
  stab_use_hysteresis = .true.
  quad_points = 3
  max_iter = 10000000
  newton_max = 18
  max_workers = 1
    save_every = 10
  regime = 'OSC'
  out_dir = 'output'

  use_muscl = .false.
  n_riser = 0
  closure_riser = 2
  gas_eos = 'ISOTHERMAL'
  gas_flow_basis = 'STANDARD'
  Cg_fixed = 343.0_wp
  P_standard = 101300.0_wp
  T_standard = 293.15_wp
  dt_min_s = -1.0_wp
  dt_max_s = 0.20_wp
  dt_start_s = 0.01_wp
  perturbation_amplitude = 0.01_wp
  output_interval_s = 1.0_wp
  field_interval_s = 30.0_wp
  flux_scheme = 'ROE'
  junction_model = 'MATCHED'
  dissipation_mode = 'LAGGED'
  ! ------------------------------ Read -----------------------------------
  infile = 'input_params.nml'
  if (command_argument_count() >= 1) then
    call get_command_argument(1, infile)
  end if

  write(output_unit,'(a,1x,a)') '[INFO] Reading input file:', trim(infile)
  flush(output_unit)

  open(unit=10, file=trim(infile), status='old', action='read', iostat=ios, iomsg=iomsg)
  if (ios /= 0) then
    write(error_unit,'(a,1x,a)') '[ERROR] Could not open input file. file=', trim(infile)
    write(error_unit,'(a,1x,a)') '  iomsg:', trim(iomsg)
    flush(error_unit)
    error stop 1
  else
    read(10, nml=PARAMS, iostat=ios, iomsg=iomsg)
    if (ios /= 0) then
      write(error_unit,'(a,i0)') '[ERROR] Failed reading namelist PARAMS. iostat=', ios
      write(error_unit,'(a,1x,a)') '  file=', trim(infile)
      write(error_unit,'(a,1x,a)') '  iomsg:', trim(iomsg)
      stop 1
    end if

    ! mul/mug are inlet mass fluxes (rho*j). If provided (>0), they will override
    ! the internally computed fluxes from jl/jg (see below).

    ! Accept alternative names used by some input decks
    if (T_scenario > 0.0_wp) Tgas = T_scenario
    if (time_phys  > 0.0_wp) time_dim = time_phys
    close(10)

    write(output_unit,'(a,i0,2a)') '[INFO] Input parsed: nodes=', nodes, '  regime=', trim(regime)
    flush(output_unit)

    ! --- If tol_res/tol_inc were not explicitly provided, follow tol ---
    if (tol_res <= 0.0_wp) tol_res = tol
    if (tol_inc <= 0.0_wp) tol_inc = tol

    ! Configure pipeline stability hysteresis thresholds (CriteriaStabFlow)
    call set_stability_params(si_on_in=stab_si_on, si_off_in=stab_si_off, enable_hysteresis=stab_use_hysteresis)

  end if

  nd = num_dummy_cells
  if (n_riser <= 0) n_riser = max(2, nodes/2)
  if (nodes<2 .or. nd<1 .or. min(X,Z,Lp,D,Cl,Rg,Tgas,P_standard,T_standard,Ps,rho_l0,P_l0)<=0.0_wp) &
    error stop 'Invalid geometry, mesh or material domain'
  if (beta_deg<0.0_wp .or. beta_deg>=90.0_wp) error stop 'beta_deg is a positive downhill angle'
  if (trim(flux_scheme)/='ROE'.and.trim(flux_scheme)/='RUSANOV') error stop 'Unknown flux_scheme'
  if (trim(junction_model)/='MATCHED'.and.trim(junction_model)/='BLOCKING'.and. &
      trim(junction_model)/='CONSERVATIVE') error stop 'Unknown junction_model'
  if (trim(dissipation_mode)/='LAGGED'.and.trim(dissipation_mode)/='FULLY_IMPLICIT') error stop 'Unknown dissipation_mode'
  if (mu_l<=0.0_wp .or. mu_g<=0.0_wp .or. eps<0.0_wp .or. eps/D>=1.0_wp) error stop 'Invalid viscosities/roughness'
  beta_rad = beta_deg * (acos(-1.0_wp)/180.0_wp)

  ! ---------------------- Geometry & properties --------------------------
  AREA = acos(-1.0_wp) * (D*D) / 4.0_wp

  CA = solve_A(Z, X)
  Lr = CA * sinh(X / CA)

  ! The derivative of an ISOTHERMAL ideal-gas EOS is R*T, not gamma*R*T.
  select case(trim(gas_eos))
  case('ISOTHERMAL')
    Cg = sqrt(Rg*Tgas)
  case('BAROTROPIC')
    Cg = Cg_fixed
  case default
    error stop 'gas_eos must be ISOTHERMAL or BAROTROPIC'
  end select

  ! bottom pressure and densities for initial fluxes
  Pb   = Ps + rho_l0 * G * Z
  rhog = Pb / (Cg*Cg)
  rhol = rho_l0 + (Pb - P_l0) / (Cl*Cl)
  ! Default inlet mass fluxes per area (kg/m^2/s): mul=rho_l*jl, mug=rho_g*jg
  ml_flux  = rho_l0 * jl
  select case(trim(gas_flow_basis))
  case('STANDARD')
    mg_flux = (P_standard/(Rg*T_standard))*jg
  case('SEPARATOR')
    mg_flux = (Ps/(Cg*Cg))*jg
  case default
    error stop 'gas_flow_basis must be STANDARD or SEPARATOR; alternatively supply mug explicitly'
  end select

  ! If the user provided mul/mug in the namelist, use them as overrides.
  ! (Useful when inlet mass flux is specified directly, instead of jl/jg.)
  if (mul >= 0.0_wp) ml_flux = mul
  if (mug >= 0.0_wp) mg_flux = mug

  ! ----------------------- Adimensionalization ---------------------------
  call fd%adimensionalizar(Cg, Cl, P_l0, rho_l0, Ps, Pb, jl, jg, ml_flux, mg_flux, &
                           nCg, nCl, nP_l0, nrho_l0, nPs, nPb, njl, njg, nmul, nmug, &
                           omega_c, omega_P, omega_rho, omega_u)

  ! target nondimensional final time
  omega_t = Lr / max(omega_u, 1.0e-16_wp)
  ntime   = time_dim / omega_t
  ! ------------------------ Steady-state (IC) ----------------------------
  call EstadoEstacionario_ndim_simp(nodes, nmul, nmug, nPs, Lp, Lr, CA, beta_rad, D, AREA, eps, G, &
                                    nCl, nCg, nrho_l0, nP_l0, mu_l, mu_g, sigma, omega_P, omega_u, omega_rho, tol, &
                                    vns, vnp, nrhogv, nrholv, alphav, nugv, nulv, thetav, njlv, njgv, &
                                    mesh_ratio=mesh_ratio, n_riser=n_riser, face_positions=physical_faces, &
                                    closure_id_riser=closure_riser, status=steady_status)
  if (steady_status/=0) then
    write(error_unit,'(a,i0)') 'Steady initialization failed, status=',steady_status
    error stop 2
  end if

  call print_steady_state_summary(nodes, vns, vnp, alphav, nrholv, nrhogv, njlv, njgv, omega_u, omega_rho, omega_P, Lr, Lp)

  Nss = size(vns)
  N   = Nss + 2*nd

  allocate(U(N,3), vns_ext(N), thetav_ext(N))
  U = 0.0_wp
  vns_ext = 0.0_wp
  thetav_ext = 0.0_wp

  do i = 1, Nss
    call fd%compute_conservative_variables(alphav(i), nrholv(i), nulv(i), nrhogv(i), nugv(i), U(nd+i,:))
    vns_ext(nd+i)   = vns(i)
    thetav_ext(nd+i)= thetav(i)
  end do

  ! Fill ghost geometry (coordinates and inclination)
  call fill_ghost_geometry(vns_ext, thetav_ext, nd)

  ! Fill ghost cells (simple linear extrapolation) before perturbation
  call apply_ghost_cells(U, nd)

  ! ------------------------- Initial perturbation -------------------------
  ! The SAME localized momentum perturbation is used for every regime label.
  ! It changes neither phase inventory; no time-periodic forcing is imposed.
  do i=nd+1,N-nd
    U(i,3)=U(i,3)*(1.0_wp+perturbation_amplitude*exp(-((vns_ext(i)-Lp/Lr)/0.15_wp)**2))
  end do
  call apply_ghost_cells(U,nd)

  ! ---------------------------- Config -----------------------------------
  cfg%regime         = trim(adjustl(regime))
  cfg%num_dummy_cells= nd
  cfg%nCells         = N
  cfg%N_pipe         = nodes
  cfg%closure_riser  = closure_riser
  cfg%flux_scheme    = flux_scheme
  cfg%junction_model = junction_model
  cfg%dissipation_mode = dissipation_mode
  cfg%gas_eos        = gas_eos
  cfg%output_interval_s = output_interval_s
  cfg%field_interval_s = field_interval_s
  cfg%perturbation_amplitude = perturbation_amplitude

  cfg%Lp   = Lp
  cfg%Lr   = Lr
  cfg%beta = beta_rad

  cfg%D_H  = D
  cfg%AREA = AREA
  cfg%eps  = eps
  cfg%G    = G
  cfg%mu_l = mu_l
  cfg%mu_g = mu_g
  cfg%sigma= sigma

  cfg%CFL  = CFL
  cfg%tol  = tol
  cfg%tola = tola
  cfg%tol_res = tol_res
  cfg%tol_inc = tol_inc
  cfg%newton_max = newton_max
  cfg%quad_points = quad_points
  cfg%use_muscl   = use_muscl
  cfg%max_workers = max_workers

  cfg%jac_dissipation_fd = jac_dissipation_fd
  cfg%jac_dissipation_fd_eps = jac_dissipation_fd_eps
  cfg%debug_dfduj_check = debug_dfduj_check
  cfg%debug_dfduj_cell  = debug_dfduj_cell
  cfg%debug_dfduj_eps   = debug_dfduj_eps

  cfg%save_every = save_every
  cfg%nCg      = nCg
  cfg%nCl      = nCl
  cfg%nP_l0    = nP_l0
  cfg%nrho_l0  = nrho_l0
  cfg%nPs      = nPs
  cfg%nPb      = nPb
  cfg%mul      = nmul
  cfg%mug      = nmug
  cfg%njl      = njl
  cfg%njg      = njg
  cfg%omega_P  = omega_P
  cfg%omega_rho= omega_rho
  cfg%omega_u  = omega_u
  cfg%omega_t = omega_t
  cfg%T       = ntime

  ! --- output directory ---
  cfg%out_dir = trim(out_dir)
  call execute_command_line('mkdir -p -- ' // shell_quote(trim(cfg%out_dir)),exitstat=ios)
  if(ios/=0)error stop 'Unable to create output directory'

  ! --- implicit time-step bounds (non-dimensional) ---
  ! If the user provides dt_min/dt_max/dt_start/cfl_target in the namelist, use them.
  ! Otherwise, fall back to safe defaults based on total non-dimensional time cfg%T.
  if (dt_min_s > 0.0_wp) then
    cfg%dt_min = dt_min_s / omega_t
  else if (dt_min > 0.0_wp) then
    cfg%dt_min = dt_min
  else
    cfg%dt_min = 1.0e-9_wp
  end if

  if (dt_max > 0.0_wp) then
    cfg%dt_max = dt_max
  else
    cfg%dt_max = dt_max_s / omega_t
  end if

  if (dt_start > 0.0_wp) then
    cfg%dt_start = dt_start
  else
    cfg%dt_start = dt_start_s / omega_t
  end if

  if (cfl_target > 0.0_wp) then
    cfg%cfl_target = cfl_target
  else
    cfg%cfl_target = 5.0_wp
  end if
  if (allocated(cfg%U0))     deallocate(cfg%U0)
  if (allocated(cfg%vns))    deallocate(cfg%vns)
  if (allocated(cfg%thetav)) deallocate(cfg%thetav)

  allocate(cfg%U0(N,3), cfg%vns(N), cfg%thetav(N))
  cfg%U0     = U
  cfg%vns    = vns_ext
  cfg%thetav = thetav_ext
  allocate(cfg%cell_width(N), cfg%face_positions(0:Nss))
  cfg%face_positions=physical_faces
  do i=1,Nss
    cfg%cell_width(nd+i)=cfg%face_positions(i)-cfg%face_positions(i-1)
  end do
  cfg%cell_width(:nd)=cfg%cell_width(nd+1)
  cfg%cell_width(N-nd+1:)=cfg%cell_width(N-nd)
  open(newunit=ios,file=trim(out_dir)//'/run_metadata.txt',status='replace')
  write(ios,nml=PARAMS)
  write(ios,'(a,8(es24.16,1x))') 'Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u: ', &
       Lr,CA,Cg,ml_flux,mg_flux,omega_P,omega_rho,omega_u
  write(ios,'(a)') 'FV1 + backward Euler; mass fluxes kg/m2/s; pressure Pa; time s. Regime is a label only.'
  close(ios)

  ! ---------------------------- Simulate ---------------------------------
  call sim%init(cfg)
  call sim%simulate(U, vns_ext, thetav_ext, cfg%T, max_iter=max_iter, &
                   it_target=it_target, use_it_target=use_it_target, it_p=it_p, &
                   dt_fac_min=dt_fac_min, dt_fac_max=dt_fac_max, &
                   dt_growth=dt_growth, dt_shrink=dt_shrink, &
                   newton_fast=newton_fast, newton_slow=newton_slow, status=sim_status)
  if (sim_status/=0) then
    write(error_unit,'(a,i0)') 'Transient did not reach requested final time; status=',sim_status
    error stop 3
  end if

  ! ----------------------------- Output ----------------------------------
  outname = trim(cfg%out_dir)//'/results_'//trim(adjustl(regime))//'.dat'
  open(unit=20, file=trim(outname), status='replace', action='write')
  write(20,'(a)') '# s_nd  s_m   alpha   P_Pa   ul_mps   ug_mps'

  do i = nd+1, N-nd
    call compute_primitive_variables( &
      U(i,1), U(i,2), U(i,3), thetav_ext(i), &
      nCg, nCl, nrho_l0, nP_l0, &
      D, AREA, eps, G, mu_l, mu_g, sigma, njl, &
      omega_u, omega_rho, tol, tola, max_iter, &
      alpha_i, nrhol_i, nrhog_i, P_nd, ul_nd, ug_nd, phase)

    rho_l_dim = nrhol_i * omega_rho
    rho_g_dim = nrhog_i * omega_rho
    ul_dim    = ul_nd * omega_u
    ug_dim    = ug_nd * omega_u

    ! barotropic pressure (same as Python steady_at)
    P_dim = P_l0 + (Cl*Cl) * (rho_l_dim - rho_l0)

    s_dim = vns_ext(i) * Lr

    write(20,'(6(1x,es16.8))') vns_ext(i), s_dim, alpha_i, P_dim, ul_dim, ug_dim
  end do

  close(20)
  write(*,'(a,1x,a)') 'Done. Wrote', trim(outname)

contains
  function shell_quote(s) result(quoted)
    character(len=*),intent(in)::s
    character(len=:),allocatable::quoted
    integer::k
    quoted=achar(39)
    do k=1,len(s)
      if(s(k:k)==achar(39))then
        quoted=quoted//achar(39)//achar(92)//achar(39)//achar(39)
      else
        quoted=quoted//s(k:k)
      endif
    enddo
    quoted=quoted//achar(39)
  end function

  subroutine print_steady_state_summary(nodes, vns, vnp, alphav, nrholv, nrhogv, njlv, njgv, omega_u, omega_rho, omega_P, Lr, Lp)
    use mod_kinds, only: wp, i4
    integer(i4), intent(in) :: nodes
    real(wp), intent(in)    :: vns(:), vnp(:), alphav(:), nrholv(:), nrhogv(:), njlv(:), njgv(:)
    real(wp), intent(in)    :: omega_u, omega_rho, omega_P, Lr, Lp
    integer(i4) :: n, i_in, i_base, i_out, i
    real(wp)    :: s_base_nd
    character(len=5) :: tag

    n = size(vns)
    i_in   = 1
    i_out  = n
    ! "BASE" = juncao pipeline-riser. Para malhas uniformes em s_nd = s/Lr:
    ! s_base_nd ~ Lp/Lr. Usa o ponto mais proximo.
    s_base_nd = Lp / max(Lr, tiny(1.0_wp))
    i_base    = int(minloc(abs(vns - s_base_nd), dim=1), kind=i4)

    write(*,'(a)') ''
    write(*,'(a)') ' >>> Estado estacionario usado como condicao inicial (resumo) <<<'
    write(*,'(a)') '  (valores dimensionais recuperados via omegas)'
    write(*,'(a)') '  idx   tag    s_nd         s_m[m]       alpha       P[bar]      jl[m/s]     jg[m/s]     rho_l      rho_g'

    do i = 1, n
      tag = '     '
      if (i == i_in) then
        tag = 'INLET'
      else if (i == i_base) then
        tag = 'BASE '
      else if (i == i_out) then
        tag = 'OUT  '
      end if
      call print_point(tag, i, vns, vnp, alphav, nrholv, nrhogv, njlv, njgv, omega_u, omega_rho, omega_P, Lr)
    end do

  end subroutine print_steady_state_summary

  subroutine print_point(tag, idx, vns, vnp, alphav, nrholv, nrhogv, njlv, njgv, omega_u, omega_rho, omega_P, Lr)
    use mod_kinds, only: wp, i4
    character(len=*), intent(in) :: tag
    integer(i4), intent(in) :: idx
    real(wp), intent(in)    :: vns(:), vnp(:), alphav(:), nrholv(:), nrhogv(:), njlv(:), njgv(:)
    real(wp), intent(in)    :: omega_u, omega_rho, omega_P, Lr
    real(wp) :: s_nd, s_m, alpha, P_bar, jl_dim, jg_dim, rhol_dim, rhog_dim

    s_nd     = vns(idx)
    s_m      = vns(idx) * Lr
    alpha    = alphav(idx)
    P_bar    = (vnp(idx) * omega_P) / 1.0e5_wp
    jl_dim   = njlv(idx) * omega_u
    jg_dim   = njgv(idx) * omega_u
    rhol_dim = nrholv(idx) * omega_rho
    rhog_dim = nrhogv(idx) * omega_rho

    write(*,'(i5,1x,a5,1x,es11.3,1x,es11.3,1x,es10.3,1x,es10.3,1x,es10.3,1x,es10.3,1x,es10.3,1x,es10.3)') &
         idx, tag, s_nd, s_m, alpha, P_bar, jl_dim, jg_dim, rhol_dim, rhog_dim

  end subroutine print_point


  subroutine fill_ghost_geometry(vns, theta, nd)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(wp), intent(inout) :: vns(:)
    real(wp), intent(inout) :: theta(:)
    integer, intent(in) :: nd
    integer :: N, k, iL, iR
    real(wp) :: dxL, dxR

    N = size(vns)
    if (N < 2*nd + 2) return

    iL = nd + 1
    dxL = vns(iL+1) - vns(iL)
    if (.not. ieee_is_finite(dxL) .or. abs(dxL) < 1.0e-12_wp) dxL = 1.0e-3_wp

    do k = 1, nd
      vns(iL-k) = vns(iL) - real(k,wp)*dxL
      theta(iL-k) = theta(iL)
    end do

    iR = N - nd
    dxR = vns(iR) - vns(iR-1)
    if (.not. ieee_is_finite(dxR) .or. abs(dxR) < 1.0e-12_wp) dxR = dxL

    do k = 1, nd
      vns(iR+k) = vns(iR) + real(k,wp)*dxR
      theta(iR+k) = theta(iR)
    end do

  end subroutine fill_ghost_geometry

end program main
