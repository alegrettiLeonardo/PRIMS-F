module module_lapack_small
  !! Small matrices use reference LAPACK algorithms; no characteristic
  !! polynomial is formed, so translated and repeated spectra are retained.
  use mod_kinds, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: dgeev_3x3, inv3x3, matvec3
  interface
    subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
      import wp
      character(len=1), intent(in) :: jobvl, jobvr
      integer, intent(in) :: n, lda, ldvl, ldvr, lwork
      real(wp), intent(inout) :: a(lda,*), work(*)
      real(wp), intent(out) :: wr(*), wi(*), vl(ldvl,*), vr(ldvr,*)
      integer, intent(out) :: info
    end subroutine dgeev
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import wp
      integer, intent(in) :: n, nrhs, lda, ldb
      real(wp), intent(inout) :: a(lda,*), b(ldb,*)
      integer, intent(out) :: ipiv(*), info
    end subroutine dgesv
  end interface
contains
  subroutine dgeev_3x3(Ain, wr, wi, VR, info)
    !! Complex pairs retain LAPACK's real/imaginary eigenvector convention.
    !! info=0 means convergence, not hyperbolicity.
    real(wp), intent(in) :: Ain(3,3)
    real(wp), intent(out) :: wr(3), wi(3), VR(3,3)
    integer, intent(out) :: info
    real(wp) :: a(3,3), vl(1,1), work(128), nan
    nan = ieee_value(0.0_wp, ieee_quiet_nan)
    wr = nan; wi = nan; VR = nan
    info = -100
    if (.not. all(ieee_is_finite(Ain))) return
    a = Ain
    call dgeev('N','V',3,a,3,wr,wi,vl,1,VR,3,work,size(work),info)
    if (info /= 0) return
    if (.not. all(ieee_is_finite(wr)) .or. .not. all(ieee_is_finite(wi)) .or. &
        .not. all(ieee_is_finite(VR))) info = -101
  end subroutine dgeev_3x3

  function inv3x3(A, info) result(Ainv)
    !! Pivoted solve A X=I. Scale alone is not a singularity criterion.
    real(wp), intent(in) :: A(3,3)
    integer, intent(out) :: info
    real(wp) :: Ainv(3,3), work(3,3), identity(3,3), residual, scale
    integer :: pivots(3), k
    Ainv = ieee_value(0.0_wp, ieee_quiet_nan)
    info = -100
    if (.not. all(ieee_is_finite(A))) return
    identity = 0.0_wp
    do k=1,3
      identity(k,k)=1.0_wp
    end do
    work=A; Ainv=identity
    call dgesv(3,3,work,3,pivots,Ainv,3,info)
    if (info /= 0) then
      Ainv=ieee_value(0.0_wp, ieee_quiet_nan)
      return
    end if
    if (.not. all(ieee_is_finite(Ainv))) then
      info=-101
      return
    end if
    scale=maxval(sum(abs(A),dim=2))*maxval(sum(abs(Ainv),dim=2))
    residual=maxval(abs(matmul(A,Ainv)-identity))
    if (residual > 128.0_wp*epsilon(1.0_wp)*max(1.0_wp,scale)) info=-102
  end function inv3x3

  pure subroutine matvec3(A,x,y)
    real(wp), intent(in) :: A(3,3),x(3)
    real(wp), intent(out) :: y(3)
    y=matmul(A,x)
  end subroutine matvec3
end module module_lapack_small
