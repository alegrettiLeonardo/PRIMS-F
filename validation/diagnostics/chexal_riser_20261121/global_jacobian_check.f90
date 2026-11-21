program global_jacobian_check
  !! Coupled-domain check of the assembled Newton Jacobian at an ACCEPTED solver
  !! state (pipe stratified + riser closure CLOSURE, MATCHED junction, Roe,
  !! LAGGED dissipation), against central differences of the FULL residual
  !! (independent of the analytic/Richardson derivative code), for three
  !! increments h = c*scale, c in {1e-4,1e-5,1e-6}. Diagnostic only.
  !! Also: junction face status, common phase flux, and the riser-side trace
  !! alpha* checked against the Chexal residual at the junction (theta=0).
  !! Usage: global_jacobian_check PARAMS STATE CLOSURE [JUNCTION [OLD_STATE]]  (JUNCTION: face/trace only)
  use mod_kinds, only: wp
  use module_config, only: sim_config_t
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use module_chexal_1992, only: chexal_residual_si
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(sim_config_t) :: cfg
  type(pipeline_riser_simulator_t) :: sim
  integer, parameter :: nd=2
  character(len=512) :: fparams, fstate, arg, mode
  character(len=16) :: eos
  real(wp) :: omega_u, ps, mul0, mug0, cg_dim, lr, dt_s, jl0, jg0
  real(wp), allocatable :: u(:,:), uold(:,:), saved(:,:), r(:,:), a(:,:,:), b(:,:,:), c(:,:,:)
  real(wp), allocatable :: rp(:,:), rm(:,:), ta(:,:,:), tb(:,:,:), tc(:,:,:), dR(:,:), col(:,:)
  real(wp) :: norm, fin(3), fout(3), h, cs(3), err(3), best, trace(6), fm(3), fp(3), res, worst(3)
  integer :: np, nc, lo, hi, ib, i, k, l, m, iu, st, closure, nbad, q, ncol(3)
  real(wp) :: theta, width, v(3)
  call get_command_argument(1, fparams); call get_command_argument(2, fstate)
  call get_command_argument(3, arg); read(arg, *) closure
  call get_command_argument(4, mode)
  open(newunit=iu, file=trim(fparams), status='old')
  read(iu, *) omega_u, ps, jl0, jg0, mul0, mug0, cg_dim, lr, dt_s, np
  read(iu, *) eos
  close(iu)
  nc = np + 2*nd; lo = nd + 1; hi = nd + np; ib = nd + np/2
  cfg = sim_config_t()
  cfg%omega_u = omega_u; cfg%omega_rho = 997.05_wp; cfg%omega_P = 101325.0_wp
  cfg%Lr = lr; cfg%Lp = 57.4_wp; cfg%omega_t = cfg%Lr/cfg%omega_u
  cfg%D_H = .05254_wp; cfg%AREA = acos(-1.0_wp)*cfg%D_H**2/4
  cfg%eps = 4.6e-5_wp; cfg%G = 9.81_wp
  cfg%mu_l = .89e-3_wp; cfg%mu_g = 1.85e-5_wp; cfg%sigma = .0728_wp
  cfg%nCl = 1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nCg = cg_dim*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0 = 1; cfg%nP_l0 = 1; cfg%nPs = ps/cfg%omega_P
  cfg%mul = mul0/(cfg%omega_rho*cfg%omega_u); cfg%mug = mug0/(cfg%omega_rho*cfg%omega_u)
  cfg%njl = jl0/cfg%omega_u; cfg%njg = jg0/cfg%omega_u
  cfg%tol = 1e-12_wp; cfg%tola = 1e-13_wp; cfg%tol_res = 1e-10_wp; cfg%tol_inc = 1e-8_wp
  cfg%closure_pipe = 1; cfg%closure_riser = closure; cfg%gas_eos = trim(eos)
  cfg%junction_model = 'MATCHED'; cfg%flux_scheme = 'ROE'; cfg%dissipation_mode = 'LAGGED'
  cfg%jac_dissipation_fd = .false.
  cfg%nCells = nc; cfg%num_dummy_cells = nd; cfg%N_pipe = np/2; cfg%quad_points = 3
  allocate(cfg%cell_width(nc), cfg%thetav(nc), u(nc,3))
  open(newunit=iu, file=trim(fstate), status='old')
  do i = lo, hi
    read(iu, *) theta, width, v
    cfg%thetav(i) = theta; cfg%cell_width(i) = width; u(i,:) = v
  enddo
  close(iu)
  do i = 1, nd
    u(i,:) = u(lo,:); u(hi+i,:) = u(hi,:)
    cfg%thetav(i) = cfg%thetav(lo); cfg%cell_width(i) = cfg%cell_width(lo)
    cfg%thetav(hi+i) = cfg%thetav(hi); cfg%cell_width(hi+i) = cfg%cell_width(hi)
  enddo
  call sim%init(cfg)
  ! Old state: phase inventories displaced by 1e-3 relative (a genuine time step).
  allocate(uold(nc,3), saved(nc,3))
  uold = u; uold(:,1) = u(:,1)*(1 - 1e-3_wp); uold(:,2) = u(:,2)*(1 + 1e-3_wp)
  if (trim(mode) == 'JUNCTION') then
    ! As the solver: LAGGED dissipation frozen at the ACCEPTED state (5th
    ! argument, a state file; default = the evaluated state itself).
    uold = u
    call get_command_argument(5, arg)
    if (len_trim(arg) > 0) then
      open(newunit=iu, file=trim(arg), status='old')
      do i = lo, hi
        read(iu, *) theta, width, v
        uold(i,:) = v
      enddo
      close(iu)
      call ghosts(uold)
    endif
  endif
  call sim%prepare_dissipation(uold, st)
  write(*,'(a,i0)') 'prepare_dissipation status=', st
  if (st /= 0) stop 1
  if (trim(mode) == 'JUNCTION') then
    call sim%face_pair(u(ib,:), u(ib+1,:), cfg%thetav(ib), cfg%thetav(ib+1), fm, fp, st, junction_trace=trace, &
                       dissipation_matrix=sim%lagged_dissipation(:,:,ib), left_cell=ib)
    write(*,'(a,i0,a,2es13.5)') 'JUNCTION status=', st, ' Gl*,Gg* (nondim)=', fm(1:2)
    write(*,'(a,6es14.6)') 'JUNCTION trace [p*,Gl*,Gg*,alpha_pipe*,alpha_riser*,rawGg*]=', trace
    stop
  endif
  allocate(r(3,np), a(3,3,np), b(3,3,np), c(3,3,np), rp(3,np), rm(3,np), ta(3,3,np), tb(3,3,np), tc(3,3,np), dR(3,np), col(3,np))
  call sim%assemble(u, uold, dt_s/cfg%omega_t, r, a, b, c, norm, fin, fout, st, .true.)
  write(*,'(a,i0,a,es12.4)') 'assemble status=', st, ' norm=', norm
  if (st /= 0) stop 2
  ! Junction face: one common phase flux, riser trace on the selected closure.
  call sim%face_pair(u(ib,:), u(ib+1,:), cfg%thetav(ib), cfg%thetav(ib+1), fm, fp, st, junction_trace=trace, &
                     dissipation_matrix=sim%lagged_dissipation(:,:,ib), left_cell=ib)
  write(*,'(a,i0,a,2es12.4,a,es10.2)') 'JUNCTION status=', st, ' Gl*,Gg*=', fm(1:2), ' |fm-fp|(1:2)=', maxval(abs(fm(1:2)-fp(1:2)))
  write(*,'(a,6es14.6)') 'JUNCTION trace [p*,Gl*,Gg*,alpha_pipe*,alpha_riser*,rawGg*]=', trace
  if (st == 0 .and. closure >= 5) then
    block
      real(wp) :: pp, rl, rg, al, jl, jg
      pp = trace(1); rl = cfg%nrho_l0 + (pp - cfg%nP_l0)/cfg%nCl**2; rg = pp/cfg%nCg**2; al = trace(5)
      jl = trace(2)/rl*cfg%omega_u; jg = trace(3)/rg*cfg%omega_u
      call chexal_residual_si(al, rl*cfg%omega_rho, rg*cfg%omega_rho, jl, jg, 0.0_wp, cfg%D_H, cfg%G, cfg%mu_l, cfg%mu_g, &
                              cfg%sigma, res, st)
      write(*,'(a,i0,a,es12.4,a,es12.4)') 'JUNCTION riser trace Chexal(theta=0) status=', st, ' residual m/s=', res, &
        ' relative=', abs(res)/max(abs(jg), 1e-300_wp)
    end block
  endif
  if (trim(mode) == 'JUNCTION') stop
  saved = u
  worst = 0; ncol = 0; nbad = 0
  write(*,'(a)') 'COLUMN cell comp section err(1e-4) err(1e-5) err(1e-6) best'
  do i = lo, hi
    do m = 1, 3
      cs = [1e-4_wp, 1e-5_wp, 1e-6_wp]
      do q = 1, 3
        h = cs(q)*max(abs(saved(i,m)), merge(1e-5_wp, 1e-2_wp, m == 2))
        u = saved; u(i,m) = u(i,m) + h; call ghosts(u)
        call sim%assemble(u, uold, dt_s/cfg%omega_t, rp, ta, tb, tc, norm, fin, fout, st, .false.)
        if (st /= 0) then
          err(q) = huge(1.0_wp); cycle
        endif
        u = saved; u(i,m) = u(i,m) - h; call ghosts(u)
        call sim%assemble(u, uold, dt_s/cfg%omega_t, rm, ta, tb, tc, norm, fin, fout, st, .false.)
        if (st /= 0) then
          err(q) = huge(1.0_wp); cycle
        endif
        dR = (rp - rm)/(2*h)
        l = i - nd; col = 0
        col(:,l) = b(:,m,l)
        if (l > 1) col(:,l-1) = c(:,m,l-1)
        if (l < np) col(:,l+1) = a(:,m,l+1)
        err(q) = maxval(abs(col - dR)/max(1.0_wp, abs(dR)))
      enddo
      best = minval(err)
      k = merge(1, 3, i < ib); if (i == ib .or. i == ib+1) k = 2
      worst(k) = max(worst(k), best); ncol(k) = ncol(k) + 1
      if (best > 3e-4_wp) nbad = nbad + 1
      write(*,'(a,2i5,i3,4es11.3)') 'COLUMN', i, m, k, err, best
    enddo
  enddo
  u = saved
  write(*,'(a,3es11.3,a,3i5,a,i0)') 'SUMMARY best-h max error pipe/junction/riser=', worst, ' columns=', ncol, ' columns_above_3e-4=', nbad
contains
  subroutine ghosts(x)
    real(wp), intent(inout) :: x(:,:)
    integer :: j
    do j = 1, nd
      x(j,:) = x(lo,:); x(hi+j,:) = x(hi,:)
    enddo
  end subroutine
end program global_jacobian_check
