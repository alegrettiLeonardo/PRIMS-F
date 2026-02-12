program debug_closure_fd
  !! Finite-difference check for closure Jacobians:
  !! compares Dlei_fechamento_or_ndim_simp (analytic) vs centered FD.
  !!
  !! Compile (example):
  !!   gfortran -std=f2008 -O0 -g -Jtmp -Itmp \
  !!     scr/mod_kinds.f90 scr/mod_void_helpers.f90 scr/module_closure_law.f90 \
  !!     debug_closure_fd.f90 -o bin/debug_closure_fd.out
  !!
  !! NOTE: adjust the compile line to match your project structure.

  use mod_kinds,          only: wp
  use module_closure_law, only: lei_fechamento_or_ndim_simp, Dlei_fechamento_or_ndim_simp
  implicit none

  real(wp) :: alpha, rhol, rhog, ul, ug
  real(wp) :: theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol
  real(wp) :: f0, fp, fm, h
  real(wp) :: dfda, dfdrhog, dfdrhol, dfdug, dfdul
  real(wp) :: fd_dfda, fd_dfdrhog, fd_dfdrhol, fd_dfdug, fd_dfdul

  call case_default(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol)

  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, f0)
  call Dlei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, &
                                   dfda, dfdrhog, dfdrhol, dfdug, dfdul)

  ! alpha
  h = 1.0e-6_wp * max(abs(alpha), 1.0_wp)
  call lei_fechamento_or_ndim_simp(alpha+h, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fp)
  call lei_fechamento_or_ndim_simp(alpha-h, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fm)
  fd_dfda = (fp - fm) / (2.0_wp*h)

  ! rhog
  h = 1.0e-6_wp * max(abs(rhog), 1.0_wp)
  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog+h, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fp)
  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog-h, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fm)
  fd_dfdrhog = (fp - fm) / (2.0_wp*h)

  ! rhol
  h = 1.0e-6_wp * max(abs(rhol), 1.0_wp)
  call lei_fechamento_or_ndim_simp(alpha, rhol+h, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fp)
  call lei_fechamento_or_ndim_simp(alpha, rhol-h, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fm)
  fd_dfdrhol = (fp - fm) / (2.0_wp*h)

  ! ug
  h = 1.0e-6_wp * max(abs(ug), 1.0_wp)
  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug+h, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fp)
  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug-h, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fm)
  fd_dfdug = (fp - fm) / (2.0_wp*h)

  ! ul
  h = 1.0e-6_wp * max(abs(ul), 1.0_wp)
  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul+h, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fp)
  call lei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul-h, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, fm)
  fd_dfdul = (fp - fm) / (2.0_wp*h)

  write(*,'(a)') '============================================================='
  write(*,'(a)') ' Closure Jacobian check: Dlei_fechamento_or_ndim_simp'
  write(*,'(a)') ' (analytic vs centered finite-difference)'
  write(*,'(a)') '============================================================='
  write(*,'(a,1pe12.4)') ' f0 = ', f0

  call print_cmp('df/dalpha ', dfda,   fd_dfda)
  call print_cmp('df/drhog  ', dfdrhog,fd_dfdrhog)
  call print_cmp('df/drhol  ', dfdrhol,fd_dfdrhol)
  call print_cmp('df/dug    ', dfdug,  fd_dfdug)
  call print_cmp('df/dul    ', dfdul,  fd_dfdul)

contains

  subroutine case_default(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol)
    real(wp), intent(out) :: alpha, rhol, rhog, ul, ug
    real(wp), intent(out) :: theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol

    ! ---- Example values (edit as needed) ----
    alpha = 0.54_wp
    rhol  = 997.0_wp
    rhog  = 3.9_wp
    ul    = 1.0_wp
    ug    = 3.0_wp

    theta = -2.0_wp * acos(-1.0_wp) / 180.0_wp   ! radians (negative => stratified branch)
    DH    = 0.0508_wp
    AREA  = 0.25_wp * acos(-1.0_wp) * DH*DH
    EPS   = 4.6e-5_wp
    G     = 9.81_wp
    MUL   = 0.89e-3_wp
    MUG   = 1.85e-5_wp
    sigma = 7.28e-2_wp

    ! nondimensionalization scales (edit to match your omegas if needed)
    w_u   = 1.0_wp
    w_rho = 1.0_wp

    tol = 1.0e-10_wp
  end subroutine case_default

  subroutine print_cmp(lbl, a, b)
    character(len=*), intent(in) :: lbl
    real(wp), intent(in) :: a, b
    real(wp) :: re
    re = relerr(a,b)
    write(*,'(a,2x,a,1pe12.4,2x,a,1pe12.4,2x,a,1pe12.4)') trim(lbl), 'ana=', a, 'fd=', b, 'rel=', re
  end subroutine print_cmp

  pure real(wp) function relerr(a,b) result(r)
    real(wp), intent(in) :: a, b
    real(wp) :: den
    den = max(1.0_wp, abs(a), abs(b))
    r = abs(a-b) / den
  end function relerr

end program debug_closure_fd
