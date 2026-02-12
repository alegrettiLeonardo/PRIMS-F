module module_lapack_small
  !! Small linear algebra utilities for 3x3 real matrices.
  !!
  !! Provides:
  !!   - dgeev_3x3: eigenvalues and right eigenvectors (pure Fortran, no LAPACK)
  !!   - inv3x3: explicit inverse of a 3x3 matrix
  !!   - matvec3: y = A*x
  !!
  !! Rationale:
  !!   Some Linux distributions ship LAPACK/BLAS stacks that may internally
  !!   use additional threads or ILP64 builds. In PRIMS-F, eigen decompositions
  !!   are called frequently and sometimes inside OpenMP regions; depending on
  !!   the BLAS/LAPACK stack, this can lead to sporadic SIGSEGV.
  !!   This module removes that dependency for the 3x3 case.
  use mod_kinds, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  public :: dgeev_3x3, inv3x3, matvec3

contains

  subroutine dgeev_3x3(Ain, wr, wi, VR, info)
    !! Eigen-decomposition for a 3x3 real matrix.
    !! Returns right eigenvectors in VR (columns) and eigenvalues (wr + i*wi).
    !!
    !! NOTE:
    !! - For the hyperbolic Roe matrix we expect 3 real eigenvalues.
    !! - If complex eigenvalues are detected (or numerical issues arise),
    !!   info is set nonzero and outputs are zeroed.
    real(wp), intent(in)  :: Ain(3,3)
    real(wp), intent(out) :: wr(3), wi(3)
    real(wp), intent(out) :: VR(3,3)
    integer, intent(out)  :: info

    real(wp) :: A(3,3)
    real(wp) :: c2, c1, c0
    real(wp) :: a2, b2, c2p
    real(wp) :: p, q, disc
    real(wp) :: r, phi
    real(wp) :: roots(3)
    real(wp) :: B(3,3)
    real(wp) :: v(3), nrm
    real(wp), parameter :: pi = 3.14159265358979323846264338327950288_wp
    integer :: i

    A = Ain

    ! Characteristic polynomial: det(lambda I - A) = lambda^3 - c2*lambda^2 + c1*lambda - c0
    c2 = A(1,1) + A(2,2) + A(3,3)
    c1 = (A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3)) - (A(1,2)*A(2,1) + A(1,3)*A(3,1) + A(2,3)*A(3,2))
    c0 = det3(A)

    ! Convert to monic cubic: x^3 + a2*x^2 + b2*x + c2p = 0
    a2  = -c2
    b2  =  c1
    c2p = -c0

    ! Depressed cubic: x = t - a2/3
    p = b2 - a2*a2/3.0_wp
    q = 2.0_wp*a2*a2*a2/27.0_wp - a2*b2/3.0_wp + c2p

    disc = (0.5_wp*q)*(0.5_wp*q) + (p/3.0_wp)*(p/3.0_wp)*(p/3.0_wp)

    wi = 0.0_wp
    VR = 0.0_wp
    wr = 0.0_wp
    info = 0

    if (.not. ieee_is_finite(disc) .or. .not. ieee_is_finite(p) .or. .not. ieee_is_finite(q)) then
      info = 1
      return
    end if

    if (disc > 1.0e-22_wp) then
      ! One real root (others complex): not expected for Roe matrix
      roots = 0.0_wp
      roots(1) = cbrt_real(-0.5_wp*q + sqrt(disc)) + cbrt_real(-0.5_wp*q - sqrt(disc)) - a2/3.0_wp
      info = 2
      wr(1) = roots(1)
      wi(1) = 0.0_wp
      VR(:,1) = 0.0_wp
      return
    else
      ! Three real roots (disc <= ~0)
      if (p >= 0.0_wp) then
        ! Degenerate case (p ~ 0) -> triple root
        roots = -a2/3.0_wp
      else
        r = sqrt(-p/3.0_wp)
        if (r <= 0.0_wp .or. .not. ieee_is_finite(r)) then
          info = 3
          return
        end if
        phi = acos( clamp( (-0.5_wp*q) / (r*r*r), -1.0_wp, 1.0_wp ) )
        roots(1) = 2.0_wp*r*cos(phi/3.0_wp)           - a2/3.0_wp
        roots(2) = 2.0_wp*r*cos((phi+2.0_wp*pi)/3.0_wp) - a2/3.0_wp
        roots(3) = 2.0_wp*r*cos((phi+4.0_wp*pi)/3.0_wp) - a2/3.0_wp
      end if
    end if

    ! Fill eigenvalues
    wr = roots

    ! Build right eigenvectors
    do i = 1, 3
      B = A
      B(1,1) = B(1,1) - wr(i)
      B(2,2) = B(2,2) - wr(i)
      B(3,3) = B(3,3) - wr(i)

      call nullvec_from_rows(B, v, nrm)
      if (nrm <= 1.0e-28_wp .or. .not. ieee_is_finite(nrm)) then
        info = 4
        VR = 0.0_wp
        wr = 0.0_wp
        wi = 0.0_wp
        return
      end if
      VR(:,i) = v / nrm
    end do

    ! Sort eigenpairs by eigenvalue (ascending) for reproducibility
    call sort_eigenpairs(wr, wi, VR)

    if (.not. all(ieee_is_finite(wr)) .or. .not. all(ieee_is_finite(VR))) then
      info = 5
      VR = 0.0_wp
      wr = 0.0_wp
      wi = 0.0_wp
      return
    end if

  end subroutine dgeev_3x3


  function inv3x3(A, info) result(Ainv)
    !! Explicit inverse of a 3x3 matrix.
    !! info = 0 success, info != 0 singular/ill-conditioned.
    real(wp), intent(in) :: A(3,3)
    integer, intent(out) :: info
    real(wp) :: Ainv(3,3)

    real(wp) :: det

    Ainv(1,1) =  A(2,2)*A(3,3) - A(2,3)*A(3,2)
    Ainv(1,2) = -A(1,2)*A(3,3) + A(1,3)*A(3,2)
    Ainv(1,3) =  A(1,2)*A(2,3) - A(1,3)*A(2,2)

    Ainv(2,1) = -A(2,1)*A(3,3) + A(2,3)*A(3,1)
    Ainv(2,2) =  A(1,1)*A(3,3) - A(1,3)*A(3,1)
    Ainv(2,3) = -A(1,1)*A(2,3) + A(1,3)*A(2,1)

    Ainv(3,1) =  A(2,1)*A(3,2) - A(2,2)*A(3,1)
    Ainv(3,2) = -A(1,1)*A(3,2) + A(1,2)*A(3,1)
    Ainv(3,3) =  A(1,1)*A(2,2) - A(1,2)*A(2,1)

    det = A(1,1)*Ainv(1,1) + A(1,2)*Ainv(2,1) + A(1,3)*Ainv(3,1)

    if (abs(det) <= 1.0e-30_wp .or. .not. ieee_is_finite(det)) then
      info = 1
      Ainv = 0.0_wp
      return
    end if

    Ainv = Ainv / det
    info = 0
  end function inv3x3


  subroutine matvec3(A, x, y)
    real(wp), intent(in)  :: A(3,3)
    real(wp), intent(in)  :: x(3)
    real(wp), intent(out) :: y(3)

    y(1) = A(1,1)*x(1) + A(1,2)*x(2) + A(1,3)*x(3)
    y(2) = A(2,1)*x(1) + A(2,2)*x(2) + A(2,3)*x(3)
    y(3) = A(3,1)*x(1) + A(3,2)*x(2) + A(3,3)*x(3)
  end subroutine matvec3


  pure real(wp) function det3(A) result(d)
    real(wp), intent(in) :: A(3,3)
    d = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - &
        A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + &
        A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
  end function det3


  pure real(wp) function cbrt_real(x) result(y)
    real(wp), intent(in) :: x
    if (x >= 0.0_wp) then
      y = abs(x)**(1.0_wp/3.0_wp)
    else
      y = -abs(x)**(1.0_wp/3.0_wp)
    end if
  end function cbrt_real


  pure real(wp) function clamp(x, lo, hi) result(y)
    real(wp), intent(in) :: x, lo, hi
    y = min(hi, max(lo, x))
  end function clamp


  pure subroutine cross3(a, b, c)
    real(wp), intent(in) :: a(3), b(3)
    real(wp), intent(out) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end subroutine cross3


  pure real(wp) function norm2_3(v) result(n)
    real(wp), intent(in) :: v(3)
    n = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
  end function norm2_3


  subroutine nullvec_from_rows(B, v, nrm)
    !! Find a (non-normalized) null vector for B using cross products of rows.
    real(wp), intent(in) :: B(3,3)
    real(wp), intent(out) :: v(3)
    real(wp), intent(out) :: nrm

    real(wp) :: r1(3), r2(3), r3(3)
    real(wp) :: v12(3), v13(3), v23(3)
    real(wp) :: n12, n13, n23

    r1 = B(1,:)
    r2 = B(2,:)
    r3 = B(3,:)

    call cross3(r1, r2, v12)
    call cross3(r1, r3, v13)
    call cross3(r2, r3, v23)

    n12 = norm2_3(v12)
    n13 = norm2_3(v13)
    n23 = norm2_3(v23)

    if (n12 >= n13 .and. n12 >= n23) then
      v = v12
      nrm = n12
    elseif (n13 >= n12 .and. n13 >= n23) then
      v = v13
      nrm = n13
    else
      v = v23
      nrm = n23
    end if

    ! Degenerate: pick a safe direction
    if (nrm <= 1.0e-28_wp) then
      v = [1.0_wp, 0.0_wp, 0.0_wp]
      nrm = 1.0_wp
    end if
  end subroutine nullvec_from_rows


  subroutine sort_eigenpairs(wr, wi, VR)
    real(wp), intent(inout) :: wr(3), wi(3)
    real(wp), intent(inout) :: VR(3,3)

    integer :: i, j
    real(wp) :: tmp
    real(wp) :: col(3)

    do i = 1, 2
      do j = i+1, 3
        if (wr(j) < wr(i)) then
          tmp = wr(i); wr(i) = wr(j); wr(j) = tmp
          tmp = wi(i); wi(i) = wi(j); wi(j) = tmp
          col = VR(:,i); VR(:,i) = VR(:,j); VR(:,j) = col
        end if
      end do
    end do
  end subroutine sort_eigenpairs

end module module_lapack_small
