module module_steady_state
  !! Steady-state utilities (pressure ODE along pipeline+riser) - port of steadyState.py
  !!
  !! Provides:
  !!   - malha_progressiva: progressive mesh (cell centers + nodal points)
  !!   - fun_dpds: dp/ds right-hand side (dimensionless)
  !!   - steady_state_ndim_simp: compute steady state along the system on a progressive mesh

  use mod_kinds, only: wp, i4
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  use mod_rootfinding, only: brent_root, bisect_root
  use mod_void_helpers, only: ffan
  use module_catenary, only: fun_or_geo
  use module_closure_law, only: lei_fechamento_or_ndim_simp, Dlei_fechamento_or_ndim_simp
  implicit none
  private

  public :: malha_progressiva
  public :: fun_dpds
  public :: steady_state_ndim_simp
  public :: EstadoEstacionario_ndim_simp

contains

  !----------------------------------------
  ! Python-compatible wrapper name
  !----------------------------------------
  subroutine EstadoEstacionario_ndim_simp(N, mul, mug, Ps, Lp, Lr, CA, beta, DH, AREA, EPS, G, Cl, Cg, rho_l0, P_l0, &
                                         mu_l_dyn, mu_g_dyn, sigma, w_p, w_u, w_rho, tol, &
                                         vns, vnp, nrhogv, nrholv, alphav, nugv, nulv, thetav, njlv, njgv)
    integer(i4), intent(in) :: N
    real(wp),    intent(in) :: mul, mug, Ps
    real(wp),    intent(in) :: Lp, Lr, CA, beta
    real(wp),    intent(in) :: DH, AREA, EPS, G
    real(wp),    intent(in) :: Cl, Cg, rho_l0, P_l0
    real(wp),    intent(in) :: mu_l_dyn, mu_g_dyn, sigma
    real(wp),    intent(in) :: w_p, w_u, w_rho, tol

    real(wp), allocatable, intent(out) :: vns(:), vnp(:)
    real(wp), allocatable, intent(out) :: nrhogv(:), nrholv(:)
    real(wp), allocatable, intent(out) :: alphav(:), nugv(:), nulv(:)
    real(wp), allocatable, intent(out) :: thetav(:)
    real(wp), allocatable, intent(out) :: njlv(:), njgv(:)

    call steady_state_ndim_simp(N, mul, mug, Ps, Lp, Lr, CA, beta, DH, AREA, EPS, G, Cl, Cg, rho_l0, P_l0, &
                               mu_l_dyn, mu_g_dyn, sigma, w_p, w_u, w_rho, tol, &
                               vns, vnp, nrhogv, nrholv, alphav, nugv, nulv, thetav, njlv, njgv)
  end subroutine EstadoEstacionario_ndim_simp

  subroutine malha_progressiva(L_total, s_riser, N_before, N_after, ratio_before, ratio_after, &
                               malha_centros, malha_nodal)
    !! Progressive mesh generator (ported from steadyState.py)
    !!
    !! Before s_riser: decreasing spacing (towards s_riser)
    !! After  s_riser: increasing spacing
    real(wp), intent(in) :: L_total, s_riser
    integer(i4), intent(in) :: N_before, N_after
    real(wp), intent(in) :: ratio_before, ratio_after
    real(wp), allocatable, intent(out) :: malha_centros(:)
    real(wp), allocatable, intent(out) :: malha_nodal(:)

    real(wp) :: L_before, L_after
    real(wp) :: r_b, r_a, d_b, d_a
    integer(i4) :: n_b, n_a
    integer(i4) :: i, n_nodes
    real(wp), allocatable :: delta_before(:), delta_after(:)

    n_b = max(1_i4, N_before)
    n_a = max(1_i4, N_after)
    r_b = ratio_before
    r_a = ratio_after

    L_before = s_riser
    L_after  = L_total - s_riser
    if (L_before < 0.0_wp) L_before = 0.0_wp
    if (L_after  < 0.0_wp) L_after  = 0.0_wp

    allocate(delta_before(n_b), delta_after(n_a))

    ! --- before s_riser (decreasing towards the junction)
    if (abs(r_b - 1.0_wp) < 1.0e-14_wp) then
      d_b = L_before / real(n_b, wp)
      do i = 1, n_b
        delta_before(i) = d_b
      end do
    else
      d_b = L_before * (1.0_wp - r_b) / (1.0_wp - r_b**real(n_b, wp))
      do i = 1, n_b
        ! reversed geometric progression: r_b^(n_b-1 .. 0)
        delta_before(i) = d_b * r_b**real(n_b - i, wp)
      end do
    end if

    ! --- after s_riser (increasing away from the junction)
    if (abs(r_a - 1.0_wp) < 1.0e-14_wp) then
      d_a = L_after / real(n_a, wp)
      do i = 1, n_a
        delta_after(i) = d_a
      end do
    else
      d_a = L_after * (1.0_wp - r_a) / (1.0_wp - r_a**real(n_a, wp))
      do i = 1, n_a
        delta_after(i) = d_a * r_a**real(i - 1, wp)
      end do
    end if

    ! nodal mesh: [0 .. s_riser] with n_b cells + [s_riser .. L_total] with n_a cells
    n_nodes = (n_b + 1) + n_a
    allocate(malha_nodal(n_nodes))
    malha_nodal(1) = 0.0_wp
    do i = 1, n_b
      malha_nodal(i+1) = malha_nodal(i) + delta_before(i)
    end do
    ! force exact junction
    malha_nodal(n_b+1) = s_riser
    do i = 1, n_a
      malha_nodal(n_b+1+i) = s_riser + sum(delta_after(1:i))
    end do
    ! force exact end
    malha_nodal(n_nodes) = L_total

    allocate(malha_centros(n_nodes-1))
    do i = 1, n_nodes-1
      malha_centros(i) = 0.5_wp * (malha_nodal(i) + malha_nodal(i+1))
    end do

    deallocate(delta_before, delta_after)
  end subroutine malha_progressiva


  real(wp) function fun_dpds(s, p, mul, mug, Lp, Lr, CA, beta, DH, AREA, EPS, G, &
                             Cl, Cg, rho_l0, P_l0, mu_l_dyn, mu_g_dyn, sigma, w_p, w_u, w_rho, tol) result(dpds)
    !! RHS for dp/ds along pipeline-riser (dimensionless).
    !! Direct port of steadyState.py:fun_dpds (without printing).
    real(wp), intent(in) :: s, p
    real(wp), intent(in) :: mul, mug
    real(wp), intent(in) :: Lp, Lr, CA, beta
    real(wp), intent(in) :: DH, AREA, EPS, G
    real(wp), intent(in) :: Cl, Cg, rho_l0, P_l0
    real(wp), intent(in) :: mu_l_dyn, mu_g_dyn, sigma
    real(wp), intent(in) :: w_p, w_u, w_rho, tol

    real(wp) :: PIP, PIG, PIF
    real(wp) :: nrhog, nrhol, njl, njg
    real(wp) :: theta, alpha
    real(wp) :: nug, nul
    real(wp) :: dfda, dfdrhog, dfdrhol, dfdug, dfdul
    real(wp) :: A11, A12, A21, A22
    real(wp) :: B1, B2
    real(wp) :: auxa, auxb, auxc
    real(wp) :: rhom, jt, dmu, Rem, fm

    ! Nondimensional groups
    PIP = w_p / (w_rho * (w_u*w_u))
    PIG = G * Lr / (w_u*w_u)
    PIF = Lr / DH

    nrhog = p / (Cg*Cg)
    nrhol = rho_l0 + (p - P_l0) / (Cl*Cl)

    njl = mul / nrhol
    njg = mug / nrhog

    theta = fun_or_geo(s * Lr, Lp, beta, CA)

    alpha = solve_alpha_from_superficials(njl, njg, nrhol, nrhog, theta, DH, AREA, EPS, G, mu_l_dyn, mu_g_dyn, sigma, w_u, w_rho, tol)

    nug = njg / max(alpha, 1.0e-14_wp)
    nul = njl / max(1.0_wp - alpha, 1.0e-14_wp)

    call Dlei_fechamento_or_ndim_simp(alpha, nrhol, nrhog, nul, nug, theta, DH, AREA, EPS, G, mu_l_dyn, mu_g_dyn, sigma, w_u, w_rho, tol, &
                                     dfda, dfdrhol, dfdrhog, dfdul, dfdug)

    ! A matrix entries (i,j = 1,2)
    A11 = (1.0_wp - alpha) * nrhol + nrhol * nul * dfdul / dfda
    A12 = nrhol * nul * dfdug / dfda
    A21 = -nrhog * nug * dfdug / dfda
    A22 = alpha * nrhog - nrhog * nug * dfdug / dfda

    ! B vector entries
    B1 = -((1.0_wp - alpha) * nul / (Cl*Cl) + nrhol * nul * (dfdrhol / (Cl*Cl) + dfdrhog / (Cg*Cg)) / dfda)
    B2 = -(alpha * nug / (Cg*Cg) - nrhog * nug * (dfdrhol / (Cl*Cl) + dfdrhog / (Cg*Cg)) / dfda)

    ! Aux terms that multiply dp/ds
    auxa = (B1 * A22 - B2 * A12) / (A22 * A11 - A21 * A12)
    auxb = (B2 * A11 - B1 * A21) / (A22 * A11 - A21 * A12)
    auxc = mul * auxa + mug * auxb + PIP

    rhom = (1.0_wp - alpha) * nrhol + alpha * nrhog
    auxa = -PIG * sin(theta) * rhom

    jt = njl + njg
    dmu = mu_g_dyn / mu_l_dyn
    Rem = (w_rho * w_u * DH / mu_l_dyn) * rhom * abs(jt) / ((1.0_wp - alpha) + dmu * alpha)

    fm = ffan(EPS / DH, Rem)
    auxb = -PIF * rhom * fm * jt * abs(jt) / 2.0_wp

    dpds = (auxa + auxb) / auxc
  end function fun_dpds


  subroutine steady_state_ndim_simp(N, mul, mug, Ps, Lp, Lr, CA, beta, DH, AREA, EPS, G, &
                                   Cl, Cg, rho_l0, P_l0, mu_l_dyn, mu_g_dyn, sigma, w_p, w_u, w_rho, tol, &
                                   vns, vnp, nrhogv, nrholv, alphav, nugv, nulv, thetav, njlv, njgv)
    !! Steady-state solver (nondimensional) on a progressive mesh.
    !!
    !! Port of steadyState.py:EstadoEstacionario_ndim_simp (centers-based mesh).
    integer(i4), intent(in) :: N
    real(wp), intent(in) :: mul, mug, Ps
    real(wp), intent(in) :: Lp, Lr, CA, beta
    real(wp), intent(in) :: DH, AREA, EPS, G
    real(wp), intent(in) :: Cl, Cg, rho_l0, P_l0
    real(wp), intent(in) :: mu_l_dyn, mu_g_dyn, sigma
    real(wp), intent(in) :: w_p, w_u, w_rho, tol

    real(wp), allocatable, intent(out) :: vns(:), vnp(:)
    real(wp), allocatable, intent(out) :: nrhogv(:), nrholv(:)
    real(wp), allocatable, intent(out) :: alphav(:), nugv(:), nulv(:)
    real(wp), allocatable, intent(out) :: thetav(:), njlv(:), njgv(:)

    real(wp), allocatable :: s_cent(:), s_nodes(:)
    real(wp), allocatable :: tspan(:), pspan(:)
    integer(i4) :: N_before, N_after, N_total
    real(wp) :: L_total_nd, s_riser
    integer(i4) :: i
    real(wp) :: nrhog, nrhol, njl, njg, theta, alpha

    ! Mesh parameters as in the python implementation
    N_before = N
    N_after  = max(1_i4, N/2)
    L_total_nd = (Lp + Lr) / Lr
    s_riser    = Lp / Lr

    call malha_progressiva(L_total_nd, s_riser, N_before, N_after, 1.04_wp, 1.04_wp, s_cent, s_nodes)
    N_total = size(s_cent)

    allocate(vns(N_total), vnp(N_total))
    allocate(nrhogv(N_total), nrholv(N_total))
    allocate(alphav(N_total), nugv(N_total), nulv(N_total))
    allocate(thetav(N_total), njlv(N_total), njgv(N_total))

    vns = s_cent

    ! integrate from top to bottom (reverse centers)
    allocate(tspan(N_total), pspan(N_total))
    do i = 1, N_total
      tspan(i) = s_cent(N_total - i + 1)
    end do

    call integrate_rk45(tspan, Ps, pspan, 1.0e-12_wp, tol)

    ! reverse back to increasing order
    do i = 1, N_total
      vnp(i) = pspan(N_total - i + 1)
    end do

    do i = 1, N_total
      nrhog = vnp(i) / (Cg*Cg)
      nrhol = rho_l0 + (vnp(i) - P_l0) / (Cl*Cl)
      nrhogv(i) = nrhog
      nrholv(i) = nrhol

      njg = mug / nrhog
      njl = mul / nrhol
      njlv(i) = njl
      njgv(i) = njg

      theta = fun_or_geo(vns(i) * Lr, Lp, beta, CA)
      thetav(i) = theta

      alpha = solve_alpha_from_superficials(njl, njg, nrhol, nrhog, theta, DH, AREA, EPS, G, mu_l_dyn, mu_g_dyn, sigma, w_u, w_rho, tol)
      alphav(i) = alpha
      nugv(i) = njg / max(alpha, 1.0e-14_wp)
      nulv(i) = njl / max(1.0_wp - alpha, 1.0e-14_wp)
    end do

    deallocate(s_cent, s_nodes, tspan, pspan)

  contains

    subroutine integrate_rk45(s_eval, y0, y_eval, rtol, atol)
      !! Adaptive Dormand-Prince RK45 integration for a scalar ODE.
      !! s_eval is monotonic (here: decreasing). y_eval is returned at s_eval points.
      real(wp), intent(in) :: s_eval(:)
      real(wp), intent(in) :: y0
      real(wp), intent(out) :: y_eval(:)
      real(wp), intent(in) :: rtol, atol

      integer(i4) :: n, i
      real(wp) :: s, y, s_tgt
      real(wp) :: h, dir
      real(wp) :: k1, k2, k3, k4, k5, k6, k7
      real(wp) :: y5, y4, err, sc, fac
      integer(i4) :: it

      n = size(s_eval)
      if (size(y_eval) /= n) error stop 'integrate_rk45: size mismatch'

      y_eval(1) = y0

      do i = 1, n-1
        s = s_eval(i)
        y = y_eval(i)
        s_tgt = s_eval(i+1)

        dir = sign(1.0_wp, s_tgt - s)  ! expected -1
        h = (s_tgt - s)
        if (abs(h) > 0.0_wp) then
          h = h / 10.0_wp
        else
          h = dir * 1.0e-6_wp
        end if

        it = 0
        do while (dir*(s - s_tgt) < 0.0_wp)
          it = it + 1
          if (it > 200000) error stop 'integrate_rk45: too many steps'

          if (abs(h) > abs(s_tgt - s)) h = s_tgt - s

          k1 = rhs(s, y)
          k2 = rhs(s + h*(1.0_wp/5.0_wp), y + h*(1.0_wp/5.0_wp)*k1)
          k3 = rhs(s + h*(3.0_wp/10.0_wp), y + h*((3.0_wp/40.0_wp)*k1 + (9.0_wp/40.0_wp)*k2))
          k4 = rhs(s + h*(4.0_wp/5.0_wp), y + h*((44.0_wp/45.0_wp)*k1 + (-56.0_wp/15.0_wp)*k2 + (32.0_wp/9.0_wp)*k3))
          k5 = rhs(s + h*(8.0_wp/9.0_wp), y + h*((19372.0_wp/6561.0_wp)*k1 + (-25360.0_wp/2187.0_wp)*k2 + (64448.0_wp/6561.0_wp)*k3 + (-212.0_wp/729.0_wp)*k4))
          k6 = rhs(s + h,                y + h*((9017.0_wp/3168.0_wp)*k1 + (-355.0_wp/33.0_wp)*k2 + (46732.0_wp/5247.0_wp)*k3 + (49.0_wp/176.0_wp)*k4 + (-5103.0_wp/18656.0_wp)*k5))

          y5 = y + h*((35.0_wp/384.0_wp)*k1 + (500.0_wp/1113.0_wp)*k3 + (125.0_wp/192.0_wp)*k4 + (-2187.0_wp/6784.0_wp)*k5 + (11.0_wp/84.0_wp)*k6)
          k7 = rhs(s + h, y5)

          y4 = y + h*((5179.0_wp/57600.0_wp)*k1 + (7571.0_wp/16695.0_wp)*k3 + (393.0_wp/640.0_wp)*k4 + (-92097.0_wp/339200.0_wp)*k5 + (187.0_wp/2100.0_wp)*k6 + (1.0_wp/40.0_wp)*k7)

          err = abs(y5 - y4)
          sc  = atol + rtol*max(abs(y), abs(y5))

          if (.not. ieee_is_finite(err) .or. .not. ieee_is_finite(sc)) then
            ! Reject and reduce step
            h = h * 0.5_wp
            cycle
          end if

          if (err <= sc) then
            ! accept
            s = s + h
            y = y5
          end if

          if (abs(err) <= 100.0_wp*epsilon(1.0_wp)) then
            fac = 5.0_wp
          else
            fac = 0.9_wp * (sc / err)**0.2_wp
            fac = min(5.0_wp, max(0.2_wp, fac))
          end if

          if (err <= sc) then
            h = h * fac
          else
            ! error estimate is O(h^5) -> use exponent 1/5
            h = h * max(0.2_wp, 0.9_wp*(sc/err)**0.2_wp)
          end if

          if (abs(h) < 1.0e-14_wp) error stop 'integrate_rk45: step underflow'
        end do

        y_eval(i+1) = y
      end do

    end subroutine integrate_rk45


    real(wp) function rhs(sloc, ploc) result(f)
      real(wp), intent(in) :: sloc, ploc
      f = fun_dpds(sloc, ploc, mul, mug, Lp, Lr, CA, beta, DH, AREA, EPS, G, Cl, Cg, rho_l0, P_l0, &
                   mu_l_dyn, mu_g_dyn, sigma, w_p, w_u, w_rho, tol)
    end function rhs

  end subroutine steady_state_ndim_simp


  real(wp) function solve_alpha_from_superficials(njl, njg, nrhol, nrhog, theta, DH, AREA, EPS, G, mu_l_dyn, mu_g_dyn, sigma, w_u, w_rho, tol) result(alpha)
    !! Solve for void fraction alpha from superficial velocities (njl, njg)
    !! using the same closure selector as lei_fechamento_or_ndim_simp.
    real(wp), intent(in) :: njl, njg, nrhol, nrhog, theta
    real(wp), intent(in) :: DH, AREA, EPS, G, mu_l_dyn, mu_g_dyn, sigma, w_u, w_rho, tol

    real(wp), parameter :: tiny_a = 1.0e-14_wp
    integer(i4), parameter :: nscan = 200_i4
    real(wp) :: a_lo, a_hi, a, b, fa, fb
    real(wp) :: best_a, best_f
    integer(i4) :: k
    integer(i4) :: info

    ! Python reference uses alpha_min = sqrt(tol) and alpha_max = 1 - sqrt(tol)
    a_lo = max(1.0e-12_wp, sqrt(max(tol, tiny_a)))
    a_hi = 1.0_wp - a_lo

    best_a = 0.5_wp
    best_f = huge(1.0_wp)

    ! scan for a sign change bracket
    fa = alpha_residual(a_lo)
    best_a = a_lo
    best_f = abs(fa)
    b = a_lo
    fb = fa

    ! Scan uniformly including both endpoints (matches np.linspace in Python)
    do k = 0, nscan-1
      if (nscan > 1) then
        a = a_lo + (a_hi - a_lo) * real(k, wp) / real(nscan-1, wp)
      else
        a = a_lo
      end if
      fa = alpha_residual(a)
      if (ieee_is_finite(fa)) then
        if (abs(fa) < best_f) then
          best_f = abs(fa)
          best_a = a
        end if
        if (ieee_is_finite(fb)) then
          if (abs(fa) <= 100.0_wp*epsilon(1.0_wp)) then
            alpha = a
            return
          end if
          if (abs(fb) <= 100.0_wp*epsilon(1.0_wp)) then
            alpha = b
            return
          end if
          if (fa*fb < 0.0_wp) then
            ! Bracket found -> use bisection (matches Python root_scalar(method='bisect'))
            block
              real(wp)   :: root
              integer(i4):: info
              root = bisect_root(alpha_residual, b, a, xtol=tol, maxit=200_i4, info=info)
              if (info == 0_i4) then
                alpha = min(a_hi, max(a_lo, root))
                return
              end if
            end block
          end if
        end if
      end if
      b  = a
      fb = fa
    end do

    ! fallback: best point from scan
    alpha = min(a_hi, max(a_lo, best_a))
    return

  contains
    real(wp) function alpha_residual(a) result(val)
      real(wp), intent(in) :: a
      real(wp) :: ul, ug
      ! Reject only invalid void fractions (outside (0,1)).
      if (a <= 0.0_wp .or. a >= 1.0_wp) then
        val = huge(1.0_wp)
        return
      end if
      ul = njl / max(1.0_wp - a, 1.0e-14_wp)
      ug = njg / max(a, 1.0e-14_wp)
      call lei_fechamento_or_ndim_simp(a, nrhol, nrhog, ul, ug, theta, DH, AREA, EPS, G, mu_l_dyn, mu_g_dyn, sigma, w_u, w_rho, tol, val)
      if (ieee_is_nan(val)) val = huge(1.0_wp)
    end function alpha_residual
  end function solve_alpha_from_superficials

end module module_steady_state
