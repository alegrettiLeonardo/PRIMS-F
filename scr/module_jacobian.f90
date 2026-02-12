module module_jacobian
  !! Fortran translation of jacobianMatrix.py
  !!
  !! Provides DFduj: Jacobian of flux vector F wrt conservative variables U=(u1,u2,u3)
  !! returning a 3x3 matrix.
  use mod_kinds, only: wp
  use module_closure_law, only: Dlei_fechamento_or_ndim_simp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: DFduj

contains

  pure real(wp) function sdiv(a, b) result(r)
    !! Safe divide with sign-preserving epsilon regularization.
    !! This avoids discontinuous Jacobian jumps when b crosses zero.
    real(wp), intent(in) :: a, b
    real(wp), parameter :: eps = 1.0e-12_wp
    real(wp) :: denom
    denom = sign(1.0_wp, b) * max(abs(b), eps)
    r = a / denom
  end function sdiv


  function DFduj(alpha, rhol, rhog, ul, ug, Cl, Cg, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, w_P, tol) result(MJF)
    !! Jacobian (3x3) of F with respect to u_j.
    real(wp), intent(in) :: alpha, rhol, rhog, ul, ug
    real(wp), intent(in) :: Cl, Cg, theta, DH, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in) :: w_u, w_rho, w_P, tol
    real(wp) :: MJF(3,3)

    real(wp) :: dfda, dfdrhog, dfdrhol, dfdug, dfdul
    real(wp) :: Cg2, Cl2, PIP, t2, t7, a1, den, tA
    real(wp) :: num00, num012
    real(wp) :: t6, t8, inv_t8, inv_den, t35
    real(wp) :: ul_over_t8, termA, termB, termC, termD
    real(wp) :: t10, den21, inv_d21, inv_t7
    real(wp) :: ug2, a1_sq, rhol2, t49, t62
    real(wp) :: term21A, term21B, term21C, term21D
    real(wp) :: t2_22, den22
    integer :: i, j

    call Dlei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, dfda, dfdrhog, dfdrhol, dfdug, dfdul)

    Cg2  = Cg*Cg
    Cl2  = Cl*Cl
    PIP  = w_P / (w_rho * (w_u*w_u))
    t2   = rhog * Cg2
    t7   = alpha * (-rhol*Cl2 + t2) - t2
    a1   = alpha - 1.0_wp
    den  = rhog*alpha*dfdul + dfdug*a1*rhol
    tA   = (-dfda*alpha*Cl2 + (dfdrhol*Cg2 + dfdrhog*Cl2) * rhog)

    MJF = 0.0_wp

    ! ---------------- F1 ----------------
    num00  = alpha * (dfdul*t7*ul - a1*tA*rhol) * rhog
    MJF(1,1) = sdiv(1.0_wp, t7) * sdiv(num00, den)

    num012 = rhol * a1 * (-dfdug*ug*t7 + alpha*rhog*(dfda*a1*Cg2 - (dfdrhol*Cg2 + dfdrhog*Cl2)*rhol))
    MJF(1,2) = sdiv(1.0_wp, t7) * sdiv(num012, den)

    MJF(1,3) = -sdiv(dfdug * rhol * (1.0_wp - alpha), den)

    ! ---------------- F2 ----------------
    MJF(2,1) = -sdiv(1.0_wp, t7) * sdiv(num00,  den)
    MJF(2,2) = -sdiv(1.0_wp, t7) * sdiv(num012, den)
    MJF(2,3) =  sdiv(dfdul * alpha * rhog, den)

    ! ---------------- F3 ----------------
    t6     = -rhol*Cl2 + t2
    t8     = alpha*t6 - t2
    inv_t8 = sdiv(1.0_wp, t8)
    inv_den= sdiv(1.0_wp, den)
    t35    = -dfda*alpha*Cl2 + (dfdrhol*Cg2 + dfdrhog*Cl2) * rhog

    ul_over_t8 = ul*inv_t8

    termA = ug * (2.0_wp*inv_t8*inv_den*alpha*a1*t35*rhol*rhog &
                 -2.0_wp*ul_over_t8*inv_den*dfdul*t8*(alpha*rhog))

    termB = (ul*ul) * inv_t8 * inv_den * (-(a1*rhol)*dfdug*t8 &
                 - (alpha*rhog)*dfdul * (-(alpha*t6) + t2))

    termC = -2.0_wp*ul_over_t8*inv_den*a1*t35*rhol*(alpha*rhog)

    termD = inv_t8*inv_den * ( -dfdug*rhol*a1*Cg2*Cl2*rhog*PIP &
                 - (rhog*rhog)*alpha*PIP*dfdul*Cl2*Cg2 )

    MJF(3,1) = termA + termB + termC + termD

    t10     = dfdul*alpha
    den21   = rhog*t10 + dfdug*(a1*rhol)
    inv_d21 = sdiv(1.0_wp, den21)
    inv_t7  = sdiv(1.0_wp, t7)
    ug2     = ug*ug
    a1_sq   = a1*a1
    rhol2   = rhol*rhol
    t49     = 2.0_wp*dfda*a1_sq*rhol*Cg2 - 2.0_wp*(dfdrhol*Cg2 + dfdrhog*Cl2) * (a1*rhol2)
    t62     = Cl2*Cg2

    term21A = ug2*inv_t7*inv_d21 * (dfdug*rhol*a1*t7 - t10*t7*rhog)
    term21B = ug * (-2.0_wp*ul*inv_t7*inv_d21*dfdug*t7*(a1*rhol) &
                    -inv_t7*inv_d21*alpha*t49*rhog)
    term21C = ul * inv_t7*inv_d21 * alpha*t49*rhog
    term21D = inv_t7*inv_d21 * (-PIP*alpha*rhog*rhol*dfdul*t62 &
                                -dfdug*t62*PIP*(a1*rhol2))

    MJF(3,2) = term21A + term21B + term21C + term21D

    t2_22  = a1*rhol
    den22  = rhog*alpha*dfdul + dfdug*t2_22
    MJF(3,3) = sdiv(2.0_wp*dfdul*ug*alpha*rhog + 2.0_wp*ul*dfdug*t2_22, den22)

    ! Clean any NaN/Inf (Python: nan_to_num)
    do i = 1, 3
      do j = 1, 3
        if (.not. ieee_is_finite(MJF(i,j))) MJF(i,j) = 0.0_wp
      end do
    end do

  end function DFduj

end module module_jacobian
