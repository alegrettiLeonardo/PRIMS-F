program boundary_jacobian_probe
  !! P0.4: point flux Jacobian A = dF/dU of the riser operator (CHEXAL_1992,
  !! closure 5) near constitutive boundaries, against independent differences
  !! of F(U) (point_state without Jacobian). Stencils never cross into an
  !! inadmissible state: central if U+-h are both admitted, else second-order
  !! one-sided on the admitted side, else none (reported). Increments
  !! h = c*max(|U_m|,1e-14), c in {1e-2,...,1e-7}; central, forward and backward
  !! whenever admissible; a kink is flagged when the two one-sided slopes keep differing (>1e-4) at the two smallest h; error per column
  !! |A(:,m) - D_h F(:,m)| / max(1,|D_h F|), best over h and its trend.
  !! Eigenvalues of A (dgeev_3x3): real spectrum and |lambda|max/|lambda|min.
  !! States are made CONSISTENT with Chexal: for given (alpha, j_l, theta) the
  !! gas flux solves the drift relation (bisection on chexal_residual_si).
  !! Usage: boundary_jacobian_probe PARAMS (same params file as global_jacobian_check)
  use mod_kinds, only: wp
  use module_config, only: sim_config_t
  use module_simulate_pipe_riser, only: pipeline_riser_simulator_t
  use module_chexal_1992, only: chexal_residual_si, CHEXAL_OK
  use module_lapack_small, only: dgeev_3x3
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  type(sim_config_t) :: cfg
  type(pipeline_riser_simulator_t) :: sim
  character(len=512) :: fparams
  character(len=16) :: eos
  real(wp) :: omega_u, ps, mul0, mug0, cg_dim, lr, dt_s, jl0, jg0, pa
  integer :: np, iu, k
  real(wp), parameter :: PI=acos(-1.0_wp), DEG=PI/180
  real(wp) :: rel_l, jl_kink, a_l
  call get_command_argument(1, fparams)
  open(newunit=iu, file=trim(fparams), status='old')
  read(iu, *) omega_u, ps, jl0, jg0, mul0, mug0, cg_dim, lr, dt_s, np
  read(iu, *) eos
  close(iu)
  cfg = sim_config_t()
  cfg%omega_u = omega_u; cfg%omega_rho = 997.05_wp; cfg%omega_P = 101325.0_wp
  cfg%Lr = lr; cfg%Lp = 57.4_wp; cfg%omega_t = cfg%Lr/cfg%omega_u
  cfg%D_H = .05254_wp; cfg%AREA = PI*cfg%D_H**2/4; cfg%eps = 4.6e-5_wp; cfg%G = 9.81_wp
  cfg%mu_l = .89e-3_wp; cfg%mu_g = 1.85e-5_wp; cfg%sigma = .0728_wp
  cfg%nCl = 1498.0_wp*sqrt(cfg%omega_rho/cfg%omega_P); cfg%nCg = cg_dim*sqrt(cfg%omega_rho/cfg%omega_P)
  cfg%nrho_l0 = 1; cfg%nP_l0 = 1; cfg%nPs = ps/cfg%omega_P
  cfg%mul = mul0/(cfg%omega_rho*cfg%omega_u); cfg%mug = mug0/(cfg%omega_rho*cfg%omega_u)
  cfg%njl = jl0/cfg%omega_u; cfg%njg = jg0/cfg%omega_u
  cfg%tol = 1e-12_wp; cfg%tola = 1e-13_wp; cfg%closure_pipe = 1; cfg%closure_riser = 5; cfg%gas_eos = trim(eos)
  cfg%junction_model = 'MATCHED'; cfg%flux_scheme = 'ROE'; cfg%dissipation_mode = 'LAGGED'
  cfg%nCells = 6; cfg%num_dummy_cells = 2; cfg%N_pipe = 1; cfg%quad_points = 3
  allocate(cfg%cell_width(6), cfg%thetav(6)); cfg%cell_width = 0.1_wp; cfg%thetav = 0.5_wp
  call sim%init(cfg)
  pa = 2.07e5_wp
  write(*,'(a)') 'STATE family label alpha jl jg theta_deg | status | col best_err trend stencil | eig real? lam_min lam_max ratio'
  ! A: alpha -> 0
  do k = 2, 6, 2
    call probe('A_alpha_to_0', 10.0_wp**(-k), 0.5_wp, 45*DEG)
  enddo
  ! B: alpha -> 1
  call probe('B_alpha_to_1', 0.99_wp, 0.01_wp, 45*DEG)
  call probe('B_alpha_to_1', 0.9999_wp, 0.01_wp, 45*DEG)
  ! C: stagnant-liquid boundary at the catenary top (OSC/SS3 blocking states)
  do k = 2, 9
    call probe('C_jl_to_0+_76deg', 0.5773_wp, 10.0_wp**(-k), 76.08_wp*DEG)
  enddo
  call probe('C_jl_to_0+_68deg', 0.4951_wp, 1e-8_wp, 68.72_wp*DEG)
  ! D: vanishing gas at the horizontal riser base (junction-trace type)
  call probe('D_jg_to_0+_0deg', 1e-3_wp, 1.0_wp, 0.0_wp)
  call probe('D_jg_to_0+_0deg', 1e-6_wp, 1.0_wp, 0.0_wp)
  ! E: correlation pieces. Re_l = 83178 (B1 = 0.8 cap and C3 = 0.5 floor coincide)
  rel_l = 60000*log(4.0_wp)
  jl_kink = rel_l*cfg%mu_l/(cfg%omega_rho*cfg%D_H)   ! approximate (rho_l ~ rho_l0)
  call probe('E_Rel_kink-1e-3', 0.5_wp, jl_kink*(1-1e-3_wp), 45*DEG)
  call probe('E_Rel_kink-1e-7', 0.5_wp, jl_kink*(1-1e-7_wp), 45*DEG)
  call probe('E_Rel_kink+1e-7', 0.5_wp, jl_kink*(1+1e-7_wp), 45*DEG)
  a_l = (1/1.15_wp)**(1/0.45_wp)                        ! L = 1 cap of Eq.3-7
  call probe('E_Lcap-1e-3', a_l*(1-1e-3_wp), 0.5_wp, 45*DEG)
  call probe('E_Lcap-1e-8', a_l*(1-1e-8_wp), 0.5_wp, 45*DEG)
  call probe('E_Lcap+1e-8', a_l*(1+1e-8_wp), 0.5_wp, 45*DEG)
  call probe('E_Lcap_exact', a_l, 0.5_wp, 45*DEG)
  call probe('E_Lcap+1e-12', a_l*(1+1e-12_wp), 0.5_wp, 45*DEG)
  ! F: angle limits of the riser
  call probe('F_theta_1e-9rad', 0.3_wp, 0.5_wp, 1e-9_wp)
  call probe('F_theta_top_76.47', 0.3_wp, 0.5_wp, 76.47_wp*DEG)
contains
  subroutine probe(label, alpha, jl, theta)
    character(*), intent(in) :: label
    real(wp), intent(in) :: alpha, jl, theta
    real(wp) :: rl, rg, jg, u(3), v(6), f(3), s(3), a(3,3), ds(3,3), fp(3), fm(3), f2(3), best, c(6)
    real(wp) :: ec(6), ef(6), eb(6), gap(6), colf(3), colb(3)
    logical :: kink
    real(wp) :: wr(3), wi(3), vr(3,3), h, lo, hi, rlo, rmid
    integer :: st, sp, sm, s2, m, q, info, ist
    character(len=8) :: sten
    character(len=36) :: trend
    rl = cfg%omega_rho + (pa - cfg%omega_P)/1498.0_wp**2
    rg = pa/cg_dim**2
    ! consistent gas flux: R(jg) = alpha(C0(jl+jg)+Vgj) - jg, bisection on [0, 50] m/s
    lo = 0; hi = 50
    call chexal_residual_si(alpha, rl, rg, jl, lo, theta, cfg%D_H, cfg%G, cfg%mu_l, cfg%mu_g, cfg%sigma, rlo, ist)
    if (ist /= CHEXAL_OK) then
      write(*,'(a,1x,a,es10.2,a,i0)') 'STATE', label, alpha, ' no consistent state, status ', ist; return
    endif
    do q = 1, 200
      jg = 0.5_wp*(lo + hi)
      call chexal_residual_si(alpha, rl, rg, jl, jg, theta, cfg%D_H, cfg%G, cfg%mu_l, cfg%mu_g, cfg%sigma, rmid, ist)
      if (ist /= CHEXAL_OK) exit
      if (sign(1.0_wp, rmid) == sign(1.0_wp, rlo)) then
        lo = jg; rlo = rmid
      else
        hi = jg
      endif
    enddo
    u = [(1-alpha)*rl/cfg%omega_rho, alpha*rg/cfg%omega_rho, (rl*jl + rg*jg)/(cfg%omega_rho*cfg%omega_u)]
    call sim%point_state(u, theta, v, f, s, a, ds, st, need_jac=.true., closure_id=5)
    write(*,'(a,1x,a,4es12.4,a,i0)') 'STATE', label, alpha, jl, jg, theta/DEG, ' | status ', st
    if (st /= 0) return
    c = [1e-2_wp, 1e-3_wp, 1e-4_wp, 1e-5_wp, 1e-6_wp, 1e-7_wp]
    do m = 1, 3
      ec = huge(1.0_wp); ef = huge(1.0_wp); eb = huge(1.0_wp); gap = 0
      do q = 1, 6
        h = c(q)*max(abs(u(m)), 1e-14_wp)
        call evalf(u, theta, m, h, fp, sp); call evalf(u, theta, m, -h, fm, sm)
        if (sp == 0 .and. sm == 0) ec(q) = maxval(abs(a(:,m) - (fp - fm)/(2*h))/max(1.0_wp, abs((fp - fm)/(2*h))))
        if (sp == 0) then
          call evalf(u, theta, m, 2*h, f2, s2)
          if (s2 == 0) then
            colf = (-3*f + 4*fp - f2)/(2*h); ef(q) = maxval(abs(a(:,m) - colf)/max(1.0_wp, abs(colf)))
          endif
        endif
        if (sm == 0) then
          call evalf(u, theta, m, -2*h, f2, s2)
          if (s2 == 0) then
            colb = (3*f - 4*fm + f2)/(2*h); eb(q) = maxval(abs(a(:,m) - colb)/max(1.0_wp, abs(colb)))
          endif
        endif
        if (ef(q) < huge(1.0_wp) .and. eb(q) < huge(1.0_wp)) gap(q) = maxval(abs(colf - colb)/max(1.0_wp, abs(colf)))
      enddo
      best = min(minval(ec), minval(ef), minval(eb))
      sten = 'central'
      if (minval(ef) == best) sten = 'forward'
      if (minval(eb) == best) sten = 'backward'
      if (minval(ec) == best) sten = 'central'
      ! a kink inside the stencil: one-sided slopes keep differing as h decreases
      kink = gap(6) > 1e-4_wp .and. gap(5) > 1e-4_wp
      write(trend,'(6f6.1)') log10(max(min(ec, ef, eb), 1e-300_wp))
      write(*,'(a,i2,es11.3,1x,a,1x,a,a,l1,a,es9.2,a,es9.2)') '   col', m, best, trend, sten, ' kink=', kink, &
        ' best_fwd=', minval(ef), ' best_bwd=', minval(eb)
    enddo
    call dgeev_3x3(a, wr, wi, vr, info)
    write(*,'(a,i0,a,l1,3es12.4)') '   eig info=', info, ' real=', all(wi == 0.0_wp), minval(abs(wr)), maxval(abs(wr)), &
      maxval(abs(wr))/max(minval(abs(wr)), tiny(1.0_wp))
  end subroutine
  subroutine evalf(u, theta, m, h, fo, st)
    real(wp), intent(in) :: u(3), theta, h
    integer, intent(in) :: m
    real(wp), intent(out) :: fo(3)
    integer, intent(out) :: st
    real(wp) :: w(3), v(6), s(3), a(3,3), ds(3,3)
    w = u; w(m) = w(m) + h
    call sim%point_state(w, theta, v, fo, s, a, ds, st, need_jac=.false., closure_id=5)
  end subroutine
end program boundary_jacobian_probe
