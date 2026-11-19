program main
  !! Driver similar to Python main.py, using module_steady_state for U0.

  use mod_kinds,               only: wp, i4
  use module_feed_contracts,only:feed_contract_t,FLOW_UNSET,gas_eos_coefficient,resolve_feed
  use module_config,           only: sim_config_t
  use module_solver_utils,     only: fluid_dynamics_t
  use module_single_phase,only:configure_pipe_riser_elevation,global_single_phase,initialize_single_phase_hydrostatic
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
  character(len=16)::boundary_method
  integer::outlet_reservoir_phase
  character(len=16)::single_phase_mode,source_scheme
  character(len=16)::liquid_flow_basis
  type(feed_contract_t)::feed
  integer::feed_status,initialization_iteration
  real(wp)::target_gl,target_gg,inlet_contract_error,initial_pin
  real(wp),allocatable::initial_face_pressures(:)
  character(len=20) :: dissipation_mode
  character(len=16) :: newton_globalization
  character(len=32) :: pipe_regime_mode
  real(wp)::checkpoint_interval_s
  real(wp) :: Cg_fixed, P_standard, T_standard, dt_min_s, dt_max_s, dt_start_s
  real(wp) :: perturbation_amplitude, output_interval_s, field_interval_s
  integer(i4) :: n_riser, closure_riser, closure_pipe, sim_status, steady_status, output_closure, output_status
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
                    gas_eos, gas_flow_basis, liquid_flow_basis, single_phase_mode, source_scheme, outlet_reservoir_phase, boundary_method, Cg_fixed, P_standard, T_standard, n_riser, closure_riser, closure_pipe, pipe_regime_mode, &
                    checkpoint_interval_s, dt_min_s, dt_max_s, dt_start_s, perturbation_amplitude, output_interval_s, field_interval_s, flux_scheme, junction_model, dissipation_mode, newton_globalization
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
  mul  = FLOW_UNSET
  mug  = FLOW_UNSET
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
  ! Sentinel distinguishes an explicitly supplied drift closure from a default.
  closure_pipe = -1
  pipe_regime_mode = 'FIXED_STRATIFIED'
  gas_eos = 'ISOTHERMAL'
  gas_flow_basis = 'STANDARD'
  liquid_flow_basis = 'LEGACY'
  single_phase_mode='STRICT';source_scheme='POINT';boundary_method='LEGACY'
  outlet_reservoir_phase=0
  checkpoint_interval_s=60
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
  newton_globalization = 'LINE_SEARCH'
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

  select case(trim(pipe_regime_mode))
  case('FIXED_STRATIFIED','AUTO_BURR15')
    if(closure_pipe/=-1.and.closure_pipe/=1) &
      error stop 'FIXED_STRATIFIED and diagnostic AUTO_BURR15 require closure_pipe=1 or omission'
    closure_pipe=1
  case('FIXED_DRIFT')
    if(closure_pipe<2.or.closure_pipe>7) &
      error stop 'FIXED_DRIFT requires explicit closure_pipe in 2..7'
  case default
    error stop 'Unknown pipe_regime_mode: use FIXED_STRATIFIED, FIXED_DRIFT or AUTO_BURR15'
  end select
  if(trim(pipe_regime_mode)=='AUTO_BURR15') &
    write(output_unit,'(a)') '[MODEL] AUTO_BURR15 is diagnostic only; pipe operator remains stratified'
  if(closure_riser<1.or.closure_riser>7)error stop 'Unsupported riser closure ID: use 1..7'
  select case(closure_riser)
  case(5)
    write(output_unit,'(a)') '[MODEL] CHEXAL_1992 air-water; branch selection required for countercurrent ambiguity'
  case(6)
    write(output_unit,'(a)') '[MODEL] CHEXAL_1992 air-water; explicitly selected high-alpha branch'
  case(7)
    write(output_unit,'(a)') '[MODEL] CHEXAL_1992 air-water; explicitly selected low-alpha branch'
  end select
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
  if(trim(newton_globalization)/='LINE_SEARCH'.and.trim(newton_globalization)/='ACTIVE_SET') &
    error stop 'Unknown newton_globalization'
  beta_rad = beta_deg * (acos(-1.0_wp)/180.0_wp)

  ! ---------------------- Geometry & properties --------------------------
  AREA = acos(-1.0_wp) * (D*D) / 4.0_wp

  CA = solve_A(Z, X)
  Lr = CA * sinh(X / CA)

  ! The derivative of an ISOTHERMAL ideal-gas EOS is R*T, not gamma*R*T.
  call gas_eos_coefficient(gas_eos,Rg,Tgas,Cg_fixed,Cg,feed_status)
  if(feed_status/=0)error stop 'Invalid gas EOS: ISOTHERMAL uses dp/drho=R*T; BAROTROPIC uses Cg_fixed^2'
  Pb=Ps+rho_l0*G*Z
  rhog=Pb/Cg**2;rhol=rho_l0+(Pb-P_l0)/Cl**2
  call resolve_feed(liquid_flow_basis,gas_flow_basis,jl,jg,mul,mug,rho_l0,Ps,Cg,Rg,T_standard,P_standard,feed,feed_status)
  if(feed_status/=0)then
    write(error_unit,'(a,i0)')'Invalid/missing/conflicting inlet contract, status=',feed_status
    error stop 3
  endif
  ml_flux=feed%liquid_mass;mg_flux=feed%gas_mass
  write(output_unit,'(4a)')'[INLET] resolved liquid=',trim(feed%liquid_mode),' gas=',trim(feed%gas_mode)
  ! ----------------------- Adimensionalization ---------------------------
  call fd%adimensionalizar(Cg, Cl, P_l0, rho_l0, Ps, Pb, jl, jg, ml_flux, mg_flux, &
                           nCg, nCl, nP_l0, nrho_l0, nPs, nPb, njl, njg, nmul, nmug, &
                           omega_c, omega_P, omega_rho, omega_u)

  ! target nondimensional final time
  omega_t = Lr / max(omega_u, 1.0e-16_wp)
  ntime   = time_dim / omega_t
  ! ------------------------ Steady-state (IC) ----------------------------
  ! For a local volumetric contract the steady mass flux is unknown until
  ! the inlet FACE pressure is known. Iterate the actual steady ODE, with
  ! fixed EOS/geometry; each iterate conserves one mass flux per phase.
  do initialization_iteration=1,50
  call EstadoEstacionario_ndim_simp(nodes, nmul, nmug, nPs, Lp, Lr, CA, beta_rad, D, AREA, eps, G, &
                                    nCl, nCg, nrho_l0, nP_l0, mu_l, mu_g, sigma, omega_P, omega_u, omega_rho, tol, &
                                    vns, vnp, nrhogv, nrholv, alphav, nugv, nulv, thetav, njlv, njgv, &
                                    mesh_ratio=mesh_ratio, n_riser=n_riser, face_positions=physical_faces, &
                                    closure_id_riser=closure_riser, status=steady_status, closure_id_pipe=closure_pipe, &
                                    face_pressures=initial_face_pressures)
  if (steady_status/=0) then
    write(error_unit,'(a,i0)') 'Steady initialization failed, status=',steady_status
    error stop 2
  end if

    initial_pin=initial_face_pressures(1)
    target_gl=nmul;target_gg=nmug
    if(feed%liquid_mode=='LOCAL_VOLUME')target_gl=(nrho_l0+(initial_pin-nP_l0)/nCl**2)*njl
    if(feed%gas_mode=='LOCAL_VOLUME')target_gg=initial_pin/nCg**2*njg
    inlet_contract_error=0
    if(target_gl/=nmul)inlet_contract_error=abs(target_gl-nmul)/max(abs(target_gl),abs(nmul))
    if(target_gg/=nmug)inlet_contract_error=max(inlet_contract_error,abs(target_gg-nmug)/max(abs(target_gg),abs(nmug)))
    if(inlet_contract_error<=1e-11_wp)exit
    nmul=target_gl;nmug=target_gg
  enddo
  if(initialization_iteration>50)error stop 'Local-volume steady inlet fixed point did not converge; no global root exclusion claimed'
  ml_flux=nmul*omega_rho*omega_u;mg_flux=nmug*omega_rho*omega_u

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
  cfg%closure_pipe   = closure_pipe
  cfg%pipe_regime_mode = trim(pipe_regime_mode)
  cfg%flux_scheme    = flux_scheme
  cfg%junction_model = junction_model
  cfg%dissipation_mode = dissipation_mode
  cfg%newton_globalization = newton_globalization
  cfg%gas_eos        = gas_eos
  cfg%liquid_inlet_basis=feed%liquid_mode;cfg%gas_inlet_basis=feed%gas_mode
  cfg%output_interval_s = output_interval_s
  cfg%field_interval_s = field_interval_s
  cfg%perturbation_amplitude = perturbation_amplitude

  cfg%checkpoint_interval_s=checkpoint_interval_s
  cfg%boundary_method=boundary_method
  cfg%single_phase_mode=single_phase_mode;cfg%source_scheme=source_scheme;cfg%CA=CA
  cfg%outlet_reservoir_phase=outlet_reservoir_phase
  if(outlet_reservoir_phase<0.or.outlet_reservoir_phase>2)error stop 'Unknown outlet_reservoir_phase'
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
  if(trim(boundary_method)/='LEGACY'.and.trim(boundary_method)/='STEADY_TRACE')error stop 'Unknown boundary_method'
  if(trim(single_phase_mode)/='STRICT'.and.trim(single_phase_mode)/='REDUCED')error stop 'Unknown single_phase_mode'
  if(trim(source_scheme)/='POINT'.and.trim(source_scheme)/='HYDROSTATIC')error stop 'Unknown source_scheme'
  call configure_pipe_riser_elevation(cfg)
  if(source_scheme=='HYDROSTATIC')then
    if(global_single_phase(U,cfg)==0)error stop 'HYDROSTATIC currently requires a declared global single phase; no mixed-flow extension'
    if(cfg%mul==0.and.cfg%mug==0)then
      call initialize_single_phase_hydrostatic(cfg,global_single_phase(U,cfg),U)
      call apply_ghost_cells(U,nd);cfg%U0=U
    endif
  endif
  open(newunit=ios,file=trim(out_dir)//'/run_metadata.txt',status='replace')
  write(ios,nml=PARAMS)
  write(ios,'(a,8(es24.16,1x))') 'Lr,CA,Cg,Gl,Gg,omega_P,omega_rho,omega_u: ', &
       Lr,CA,Cg,ml_flux,mg_flux,omega_P,omega_rho,omega_u
  write(ios,'(a)') 'FV1 + backward Euler; mass fluxes kg/m2/s; pressure Pa; time s. Regime is a label only.'
  close(ios)
  open(newunit=ios,file=trim(out_dir)//'/inlet_contract.txt',status='replace')
  write(ios,'(a,a)')'liquid_resolved=',trim(feed%liquid_mode)
  write(ios,'(a,a)')'gas_resolved=',trim(feed%gas_mode)
  write(ios,'(a,es25.16)')'gas_dp_drho_m2_s2=',Cg**2
  write(ios,'(a,es25.16)')'rho_standard_ideal_kg_m3=',P_standard/(Rg*T_standard)
  write(ios,'(a,es25.16)')'rho_at_Pstandard_EOS_kg_m3=',P_standard/Cg**2
  write(ios,'(a,es25.16)')'initial_inlet_face_pressure_Pa=',initial_pin*omega_P
  write(ios,'(a,es25.16)')'initial_liquid_mass_flux_kg_m2_s=',ml_flux
  write(ios,'(a,es25.16)')'initial_gas_mass_flux_kg_m2_s=',mg_flux
  write(ios,'(a,es25.16)')'steady_inlet_contract_relative_error=',inlet_contract_error
  write(ios,'(a,i0)')'initialization_iterations=',initialization_iteration
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
    ! This legacy output used implicit Bendiksen in every riser. Recover with
    ! the simulated section closure, including a horizontal fixed-drift pipe.
    ! Only this post-processing file changes; monitor/fields use the operator.
    output_closure=cfg%closure_riser
    if(i<=nd+cfg%N_pipe)output_closure=cfg%closure_pipe
    call compute_primitive_variables( &
      U(i,1), U(i,2), U(i,3), thetav_ext(i), &
      nCg, nCl, nrho_l0, nP_l0, &
      D, AREA, eps, G, mu_l, mu_g, sigma, njl, &
      omega_u, omega_rho, tol, tola, max_iter, &
      alpha_i, nrhol_i, nrhog_i, P_nd, ul_nd, ug_nd, phase,status=output_status,closure_id=output_closure)
    if(output_status/=0)then
      write(error_unit,'(a,i0,a,i0)') 'Final output recovery failed at cell ',i,', status=',output_status
      error stop 4
    endif

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
