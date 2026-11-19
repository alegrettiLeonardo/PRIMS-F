module module_roe_solver
  !! Roe solver + MUSCL reconstruction (characteristic space) for the 3-equation
  !! two-phase drift-flux model.
  !!
  !! Key features:
  !! - Adaptive limiter: MC in smooth regions, Superbee near steep gradients
  !!   (hysteresis shock sensor), applied to characteristic variables.
  !! - Roe flux: F = 0.5(F_L+F_R) - 0.5*|Ahat|*(U_R-U_L)
  !!   with the exact matrix absolute value (no acoustic-scale contact fix).
  !!
  !! Robustness:
  !! - A failed or non-real eigenstructure is reported explicitly. Choosing a
  !!   dissipative alternative is the caller's decision, recorded separately
  !!   from a valid Roe flux; complex spectra are not erased.

  use mod_kinds,           only: wp, i4
  use module_lapack_small, only: dgeev_3x3, inv3x3, matvec3
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private

  public :: roe_solver_t

  public :: roe_absA_from_A
  ! Expose Roe eigenstructure so other modules (e.g., dt estimate, diagnostics)
  ! can reuse the *same* robust eigen handling (real-only, inversion checks, etc.).
  public :: roe_eigenstructure

  integer(i4), parameter :: LIMITER_MC       = 1
  integer(i4), parameter :: LIMITER_SUPERBEE = 2

  type :: roe_solver_t
    real(wp) :: shock_threshold_on  = 0.20_wp
    real(wp) :: shock_threshold_off = 0.10_wp
    real(wp) :: shock_eps           = 1.0e-12_wp
    logical  :: prev_shock_state    = .false.
    logical  :: has_prev_shock_state = .false.
  contains
    procedure :: select_adaptive_limiter
    procedure :: muscl_reconstruction_adaptive
    procedure :: muscl_reconstruction_characteristic
    procedure :: roe_riemann_solver
    procedure :: rusanov_flux
  end type roe_solver_t

contains

  pure subroutine sort_eigenpairs_3(wr, wi, R)
    !! Sort eigenpairs by increasing real eigenvalue.
    !! This is important for characteristic MUSCL: if eigenvectors reorder
    !! from face to face, the limiter sees “scrambled” characteristic variables.
    real(wp), intent(inout) :: wr(3), wi(3)
    real(wp), intent(inout) :: R(3,3)  ! columns are right eigenvectors
    integer :: p(3)
    integer :: i, j, kmin, tmp
    real(wp) :: wr0(3), wi0(3)
    real(wp) :: R0(3,3)

    p = [1, 2, 3]
    do i = 1, 2
      kmin = i
      do j = i+1, 3
        if (wr(p(j)) < wr(p(kmin))) kmin = j
      end do
      if (kmin /= i) then
        tmp   = p(i)
        p(i)  = p(kmin)
        p(kmin) = tmp
      end if
    end do

    wr0 = wr
    wi0 = wi
    R0  = R
    do i = 1, 3
      wr(i)   = wr0(p(i))
      wi(i)   = wi0(p(i))
      R(:, i) = R0(:, p(i))
    end do
  end subroutine sort_eigenpairs_3

  subroutine roe_eigenstructure(A,wr,wi,R,L,ok,status)
    real(wp), intent(in) :: A(3,3)
    real(wp), intent(out) :: wr(3),wi(3),R(3,3),L(3,3)
    logical, intent(out) :: ok
    integer(i4), intent(out), optional :: status
    integer :: info,j
    real(wp) :: residual(3,3),scale,condition
    ok=.false.
    L=ieee_value(0.0_wp,ieee_quiet_nan)
    if (present(status)) status=10
    call dgeev_3x3(A,wr,wi,R,info)
    if (info/=0) return
    if (present(status)) status=11
    ! LAPACK stores a complex pair as real and imaginary vector columns.
    ! Such columns must never be treated as independent real eigenvectors,
    ! even if the imaginary part is small compared with an arbitrary scale.
    if (any(wi/=0.0_wp)) return
    call sort_eigenpairs_3(wr,wi,R)
    if (present(status)) status=12
    L=inv3x3(R,info)
    if (info/=0) return
    condition=maxval(sum(abs(R),dim=2))*maxval(sum(abs(L),dim=2))
    if (condition>1.0_wp/sqrt(epsilon(1.0_wp))) return
    residual=matmul(A,R)
    do j=1,3
      residual(:,j)=residual(:,j)-wr(j)*R(:,j)
    end do
    scale=max(1.0_wp,maxval(abs(A)))*max(1.0_wp,maxval(abs(R)))
    if (present(status)) status=13
    if (maxval(abs(residual))>1024.0_wp*epsilon(1.0_wp)*scale) return
    ok=.true.
    if (present(status)) status=0
  end subroutine roe_eigenstructure

  pure subroutine matvec3_strided(A, Ust, irow, y)
    !! y = A * Ust(irow,:)
    !! Ust(irow,:) is strided in column-major storage; compute explicitly.
    real(wp), intent(in) :: A(3,3)
    real(wp), intent(in) :: Ust(5,3)
    integer(i4), intent(in) :: irow
    real(wp), intent(out) :: y(3)
    integer(i4) :: k
    do k = 1, 3
      y(k) = A(k,1)*Ust(irow,1) + A(k,2)*Ust(irow,2) + A(k,3)*Ust(irow,3)
    end do
  end subroutine matvec3_strided

  pure real(wp) function minmod2(a, b)
    real(wp), intent(in) :: a, b
    if (a*b <= 0.0_wp) then
      minmod2 = 0.0_wp
    else
      minmod2 = sign(1.0_wp, a) * min(abs(a), abs(b))
    end if
  end function minmod2

  pure real(wp) function minmod3(a, b, c)
    real(wp), intent(in) :: a, b, c
    if (a*b <= 0.0_wp .or. a*c <= 0.0_wp) then
      minmod3 = 0.0_wp
    else
      minmod3 = sign(1.0_wp, a) * min(abs(a), min(abs(b), abs(c)))
    end if
  end function minmod3

  pure real(wp) function superbee(a, b)
    real(wp), intent(in) :: a, b
    real(wp) :: s1, s2
    s1 = minmod2(2.0_wp*a, b)
    s2 = minmod2(a, 2.0_wp*b)
    if (abs(s1) >= abs(s2)) then
      superbee = s1
    else
      superbee = s2
    end if
  end function superbee

  pure real(wp) function mc_limiter(a, b)
    real(wp), intent(in) :: a, b
    mc_limiter = minmod3(0.5_wp*(a+b), 2.0_wp*a, 2.0_wp*b)
  end function mc_limiter

  pure real(wp) function apply_limiter(lim_id, a, b)
    integer(i4), intent(in) :: lim_id
    real(wp), intent(in) :: a, b
    select case (lim_id)
    case (LIMITER_SUPERBEE)
      apply_limiter = superbee(a, b)
    case (LIMITER_MC)
      apply_limiter = mc_limiter(a, b)
    case default
      apply_limiter = mc_limiter(a, b)
    end select
  end function apply_limiter


  integer(i4) function select_adaptive_limiter(this, U_stencil) result(lim_id)
    class(roe_solver_t), intent(in) :: this
    real(wp), intent(in) :: U_stencil(5,3)

    real(wp) :: U_left(3), U_cent(3), U_right(3)
    real(wp) :: dL(3), dR(3)
    real(wp) :: denom(3)
    real(wp) :: relL(3), relR(3)
    real(wp) :: sensor
    logical  :: is_shock

    U_left  = U_stencil(2,:)
    U_cent  = U_stencil(3,:)
    U_right = U_stencil(4,:)

    dL = U_cent  - U_left
    dR = U_right - U_cent

    denom = max(abs(U_cent), this%shock_eps)
    relL  = abs(dL) / denom
    relR  = abs(dR) / denom

    sensor = max(maxval(relL), maxval(relR))

    ! Stateless shock sensor (thread-safe):
    ! - use SUPERBEE only when the local sensor exceeds the ON threshold.
    is_shock = (sensor > this%shock_threshold_on)

    if (is_shock) then
      lim_id = LIMITER_SUPERBEE
    else
      lim_id = LIMITER_MC
    end if
  end function select_adaptive_limiter
  subroutine build_characteristic_basis(A, wr, wi, R, L, ok)
    real(wp), intent(in)  :: A(3,3)
    real(wp), intent(out) :: wr(3), wi(3)
    real(wp), intent(out) :: R(3,3)
    real(wp), intent(out) :: L(3,3)
    logical, intent(out)  :: ok

    call roe_eigenstructure(A, wr, wi, R, L, ok)
  end subroutine build_characteristic_basis


  subroutine muscl_reconstruction_adaptive(this, U_stencil, U_L, U_R)
    !! Component-wise MUSCL (fallback).
    class(roe_solver_t), intent(inout) :: this
    real(wp), intent(in) :: U_stencil(5,3)
    real(wp), intent(out) :: U_L(3), U_R(3)

    integer(i4) :: lim
    integer(i4) :: m
    real(wp) :: dL, dR, slope_i, slope_ip1

    lim = this%select_adaptive_limiter(U_stencil)

    do m = 1, 3
      dL = U_stencil(3,m) - U_stencil(2,m)
      dR = U_stencil(4,m) - U_stencil(3,m)
      slope_i = apply_limiter(lim, dL, dR)

      dL = U_stencil(4,m) - U_stencil(3,m)
      dR = U_stencil(5,m) - U_stencil(4,m)
      slope_ip1 = apply_limiter(lim, dL, dR)

      U_L(m) = U_stencil(3,m) + 0.5_wp*slope_i
      U_R(m) = U_stencil(4,m) - 0.5_wp*slope_ip1
    end do
  end subroutine muscl_reconstruction_adaptive


  subroutine muscl_reconstruction_characteristic(this, U_stencil, RoeA, U_L, U_R)
    !! MUSCL in characteristic space using the eigenbasis of RoeA.
    class(roe_solver_t), intent(inout) :: this
    real(wp), intent(in)  :: U_stencil(5,3)
    real(wp), intent(in)  :: RoeA(3,3)
    real(wp), intent(out) :: U_L(3), U_R(3)

    integer(i4) :: lim, m
    real(wp) :: wr(3), wi(3)
    real(wp) :: R(3,3), Lmat(3,3)
    logical :: ok

    real(wp) :: w2(3), w3(3), w4(3), w5(3)
    ! NOTE: Avoid name clash with eigenvalue arrays wr/wi (Fortran is case-insensitive).
    real(wp) :: wL_face(3), wR_face(3)
    real(wp) :: dL, dR, slope_i, slope_ip1

    lim = this%select_adaptive_limiter(U_stencil)

    call build_characteristic_basis(RoeA, wr, wi, R, Lmat, ok)
    if (.not. ok) then
      call this%muscl_reconstruction_adaptive(U_stencil, U_L, U_R)
      return
    end if

    ! Explicit 3x3 mat-vec products (avoid temporaries from MATMUL on strided slices).
    call matvec3_strided(Lmat, U_stencil, 2, w2)
    call matvec3_strided(Lmat, U_stencil, 3, w3)
    call matvec3_strided(Lmat, U_stencil, 4, w4)
    call matvec3_strided(Lmat, U_stencil, 5, w5)

    do m = 1, 3
      dL = w3(m) - w2(m)
      dR = w4(m) - w3(m)
      slope_i = apply_limiter(lim, dL, dR)

      dL = w4(m) - w3(m)
      dR = w5(m) - w4(m)
      slope_ip1 = apply_limiter(lim, dL, dR)

      wL_face(m) = w3(m) + 0.5_wp*slope_i
      wR_face(m) = w4(m) - 0.5_wp*slope_ip1
    end do

    call matvec3(R, wL_face, U_L)
    call matvec3(R, wR_face, U_R)

    if (.not. all(ieee_is_finite(U_L)) .or. .not. all(ieee_is_finite(U_R))) then
      call this%muscl_reconstruction_adaptive(U_stencil, U_L, U_R)
    end if
  end subroutine muscl_reconstruction_characteristic


  subroutine rusanov_flux(this,F_L,F_R,U_L,U_R,amax,flux,status)
    class(roe_solver_t), intent(in) :: this
    real(wp), intent(in) :: F_L(3),F_R(3),U_L(3),U_R(3),amax
    real(wp), intent(out) :: flux(3)
    integer(i4), intent(out), optional :: status
    flux=ieee_value(0.0_wp,ieee_quiet_nan)
    if (present(status)) status=1
    if (.not. all(ieee_is_finite([F_L,F_R,U_L,U_R,amax]))) return
    if (amax<0.0_wp) return
    flux=0.5_wp*(F_L+F_R)-0.5_wp*amax*(U_R-U_L)
    if (present(status)) status=0
  end subroutine rusanov_flux

  subroutine roe_riemann_solver(this,F_L,F_R,U_L,U_R,Roe_matrix,flux,status)
    class(roe_solver_t), intent(in) :: this
    real(wp), intent(in) :: F_L(3),F_R(3),U_L(3),U_R(3),Roe_matrix(3,3)
    real(wp), intent(out) :: flux(3)
    integer(i4), intent(out), optional :: status
    real(wp) :: absA(3,3),amax
    integer(i4) :: code
    flux=ieee_value(0.0_wp,ieee_quiet_nan)
    if (present(status)) status=1
    if (.not. all(ieee_is_finite([F_L,F_R,U_L,U_R]))) return
    call roe_absA_from_A(Roe_matrix,absA,amax,code)
    if (present(status)) status=code
    if (code/=0) return
    flux=0.5_wp*(F_L+F_R)-0.5_wp*matmul(absA,U_R-U_L)
    if (.not. all(ieee_is_finite(flux))) then
      if (present(status)) status=14
    end if
  end subroutine roe_riemann_solver

  subroutine roe_absA_from_A(A,absA,amax,info)
    !! Exact matrix absolute value, shared by flux and frozen dissipation.
    !! An acoustic-scale Harten fix on all modes would artificially diffuse
    !! the slow material wave. A transonic fix requires left/right acoustic
    !! information not present in this API; none is silently imposed here.
    real(wp), intent(in) :: A(3,3)
    real(wp), intent(out) :: absA(3,3),amax
    integer(i4), intent(out) :: info
    real(wp) :: wr(3),wi(3),R(3,3),L(3,3),scaled_R(3,3),absLam(3)
    integer :: j
    logical :: ok
    absA=ieee_value(0.0_wp,ieee_quiet_nan)
    amax=ieee_value(0.0_wp,ieee_quiet_nan)
    call roe_eigenstructure(A,wr,wi,R,L,ok,info)
    if (.not. ok) return
    absLam=abs(wr)
    do j=1,3
      scaled_R(:,j)=R(:,j)*absLam(j)
    end do
    absA=matmul(scaled_R,L)
    amax=maxval(absLam)
    if (.not. all(ieee_is_finite(absA))) info=14
  end subroutine roe_absA_from_A
end module module_roe_solver
