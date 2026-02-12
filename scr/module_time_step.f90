module module_time_step
  use mod_kinds,            only: wp, i4
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  private
  public :: time_integrator_t

  type :: time_integrator_t
  contains
    procedure :: compute_time_step_lax_wendroff
    procedure :: solver_ivp
    procedure :: compute_residual_vector
    procedure :: solve_block_tridiagonal_3x3
    procedure :: invert_3x3
  end type time_integrator_t

contains

  real(wp) function compute_time_step_lax_wendroff(this, CFL, delta_x, Roe_matrix, mode, acoustic_scale, u_material) result(dt)
    !! CFL-based time-step estimate.
    !! - mode='acoustic' (default): dt ~ CFL*dx/(acoustic_scale*max|A|)
    !! - mode='material' : dt ~ CFL*dx/max(|u_material|, tiny)
    !! For implicit stepping, users can select large CFL and/or small acoustic_scale.
    class(time_integrator_t), intent(in) :: this
    real(wp), intent(in) :: CFL
    real(wp), intent(in) :: delta_x
    real(wp), intent(in) :: Roe_matrix(3,3)
    character(len=*), intent(in), optional :: mode
    real(wp), intent(in), optional :: acoustic_scale
    real(wp), intent(in), optional :: u_material

    real(wp) :: amax, scale
    logical  :: use_material

    use_material = .false.
    if (present(mode)) then
      select case (trim(mode))
      case('material','MATERIAL','Material')
        use_material = .true.
      case default
        use_material = .false.
      end select
    end if

    if (use_material) then
      if (present(u_material)) then
        amax = abs(u_material)
      else
        amax = maxval(abs(Roe_matrix))
      end if
    else
      amax = maxval(abs(Roe_matrix))
      if (present(acoustic_scale)) then
        scale = max(acoustic_scale, 0.0_wp)
        amax  = amax * scale
      end if
    end if

    amax = max(amax, 1.0e-12_wp)
    dt   = CFL * abs(delta_x) / amax

    if (.not. ieee_is_finite(dt)) dt = 0.0_wp

    call dummy_use(this)
  end function compute_time_step_lax_wendroff


  subroutine solver_ivp(this, U, numerical_flux_values, source_term_values, delta_s, delta_t, tol, U_out)
    class(time_integrator_t), intent(in) :: this
    real(wp), intent(in) :: U(:,:)
    real(wp), intent(in) :: numerical_flux_values(:,:)
    real(wp), intent(in) :: source_term_values(:,:)
    real(wp), intent(in) :: delta_s(:)
    real(wp), intent(in) :: delta_t
    real(wp), intent(in) :: tol
    real(wp), intent(out):: U_out(:,:)

    integer :: i
    real(wp) :: ds

    U_out = U
    do i = 2, size(U,1)-1
      ds = max(abs(delta_s(i)), 1.0e-12_wp)
      U_out(i,:) = U(i,:) - (delta_t/ds) * &
                   (numerical_flux_values(i,:) - numerical_flux_values(i-1,:)) &
                   + delta_t * source_term_values(i,:)
    end do

    call dummy_use(this)
    call dummy_use_real(tol)
  end subroutine solver_ivp


  subroutine compute_residual_vector(this, U, U_n, delta_s, flux_faces, source_faces, i_first, i_last, dt, R, res_norm)
    class(time_integrator_t), intent(in) :: this
    real(wp), intent(in)  :: U(:,:)
    real(wp), intent(in)  :: U_n(:,:)
    real(wp), intent(in)  :: delta_s(:)
    real(wp), intent(in)  :: flux_faces(:,:)
    real(wp), intent(in)  :: source_faces(:,:)
    integer, intent(in)   :: i_first, i_last
    real(wp), intent(in)  :: dt
    real(wp), intent(out) :: R(:,:)
    real(wp), intent(out) :: res_norm

    integer :: i
    real(wp) :: ds

    R = 0.0_wp
    res_norm = 0.0_wp

    if (i_first < 2 .or. i_last < i_first) then
      res_norm = huge(1.0_wp)
      return
    end if

    do i = i_first, i_last
      ds = max(abs(delta_s(i)), 1.0e-12_wp)

      R(i,:) = U(i,:) - U_n(i,:) - dt * ( &
               - (flux_faces(i,:) - flux_faces(i-1,:)) / ds + source_faces(i,:) )

      if (.not. all(ieee_is_finite(R(i,:)))) then
        write(*,'(a,i0)') '[ERROR] Non-finite residual at cell i=', i
        write(*,'(a,3es25.16)') '  R(i,:)     =', R(i,1), R(i,2), R(i,3)
        write(*,'(a,3es25.16)') '  U(i,:)     =', U(i,1), U(i,2), U(i,3)
        write(*,'(a,3es25.16)') '  U_n(i,:)   =', U_n(i,1), U_n(i,2), U_n(i,3)
        write(*,'(a,3es25.16)') '  fluxR      =', flux_faces(i,1), flux_faces(i,2), flux_faces(i,3)
        write(*,'(a,3es25.16)') '  fluxL      =', flux_faces(i-1,1), flux_faces(i-1,2), flux_faces(i-1,3)
        write(*,'(a,3es25.16)') '  source(i,:) =', source_faces(i,1), source_faces(i,2), source_faces(i,3)
        write(*,'(a,es25.16,a,es25.16)') '  dt=', dt, '  ds=', ds
        res_norm = huge(1.0_wp)
        call dummy_use(this)
        return
      end if

      res_norm = res_norm + sum(R(i,:)*R(i,:))
    end do

    res_norm = sqrt(res_norm)
    if (.not. ieee_is_finite(res_norm)) res_norm = huge(1.0_wp)

    call dummy_use(this)
  end subroutine compute_residual_vector


  subroutine solve_block_tridiagonal_3x3(this, A, B, C, D, X, info)
    !! Solve a block tridiagonal system with 3x3 blocks using a Thomas-like algorithm.
    !!
    !! System form (for i=1..n):
    !!   A_i * x_{i-1} + B_i * x_i + C_i * x_{i+1} = D_i
    !! with A_1 unused and C_n unused.
    !!
    !! Dummy arguments are assumed-shape to satisfy Fortran rules and avoid
    !! nonstandard mixes like (3,3,:).
    class(time_integrator_t), intent(in) :: this
    real(wp), intent(in)  :: A(:,:,:), B(:,:,:), C(:,:,:)
    real(wp), intent(in)  :: D(:,:)
    real(wp), intent(out) :: X(:,:)
    integer(i4), intent(out) :: info

    integer(i4) :: n, i, info_loc
    real(wp), allocatable :: Cp(:,:,:), Dp(:,:), M(:,:)

    info = 0_i4

    ! Basic shape checks (avoid silent out-of-bounds / wrong calls)
    if (size(D,1) /= 3 .or. size(X,1) /= 3) then
      info = -2_i4
      call dummy_use(this)
      return
    end if
    if (size(A,1) /= 3 .or. size(A,2) /= 3 .or. size(B,1) /= 3 .or. size(B,2) /= 3 .or. size(C,1) /= 3 .or. size(C,2) /= 3) then
      info = -1_i4
      call dummy_use(this)
      return
    end if

    n = int(size(D,2), kind=i4)

    if (n < 1_i4) then
      info = -4_i4
      call dummy_use(this)
      return
    end if

    if (size(A,3) < n .or. size(B,3) < n .or. size(C,3) < n .or. size(X,2) < n) then
      info = -3_i4
      call dummy_use(this)
      return
    end if

    allocate(Cp(3,3,n), Dp(3,n), M(3,3))

    ! Forward elimination
    call solve_3x3_mat(B(:,:,1), C(:,:,1), Cp(:,:,1), info_loc)
    if (info_loc /= 0_i4) then
      info = 1_i4
      deallocate(Cp, Dp, M)
      call dummy_use(this)
      return
    end if

    call solve_3x3_vec(B(:,:,1), D(:,1), Dp(:,1), info_loc)
    if (info_loc /= 0_i4) then
      info = 1_i4
      deallocate(Cp, Dp, M)
      call dummy_use(this)
      return
    end if

    do i = 2_i4, n
      M = B(:,:,i) - matmul(A(:,:,i), Cp(:,:,i-1_i4))

      call solve_3x3_mat(M, C(:,:,i), Cp(:,:,i), info_loc)
      if (info_loc /= 0_i4) then
        info = i
        deallocate(Cp, Dp, M)
        call dummy_use(this)
        return
      end if

      call solve_3x3_vec(M, D(:,i) - matmul(A(:,:,i), Dp(:,i-1_i4)), Dp(:,i), info_loc)
      if (info_loc /= 0_i4) then
        info = i
        deallocate(Cp, Dp, M)
        call dummy_use(this)
        return
      end if
    end do

    ! Back substitution
    X(:,n) = Dp(:,n)
    do i = n-1_i4, 1_i4, -1_i4
      X(:,i) = Dp(:,i) - matmul(Cp(:,:,i), X(:,i+1_i4))
    end do

    deallocate(Cp, Dp, M)
    call dummy_use(this)
  end subroutine solve_block_tridiagonal_3x3


  subroutine invert_3x3(this, A, Ainv, info)
    !! Kept for compatibility; prefers solve_3x3_* regularized path.
    class(time_integrator_t), intent(in) :: this
    real(wp), intent(in)  :: A(3,3)
    real(wp), intent(out) :: Ainv(3,3)
    integer(i4), intent(out)  :: info

    real(wp) :: det
    real(wp), parameter :: eps_reg = 1.0e-12_wp
    real(wp) :: Areg(3,3)

    Areg = A
    det = det3(Areg)
    if (.not. ieee_is_finite(det) .or. abs(det) < eps_reg) then
      Areg(1,1)=Areg(1,1)+eps_reg
      Areg(2,2)=Areg(2,2)+eps_reg
      Areg(3,3)=Areg(3,3)+eps_reg
      det = det3(Areg)
    end if

    if (.not. ieee_is_finite(det) .or. abs(det) < eps_reg) then
      info = 1_i4
      Ainv = 0.0_wp
      call dummy_use(this)
      return
    end if

    Ainv(1,1) =  (Areg(2,2)*Areg(3,3) - Areg(2,3)*Areg(3,2)) / det
    Ainv(1,2) = -(Areg(1,2)*Areg(3,3) - Areg(1,3)*Areg(3,2)) / det
    Ainv(1,3) =  (Areg(1,2)*Areg(2,3) - Areg(1,3)*Areg(2,2)) / det
    Ainv(2,1) = -(Areg(2,1)*Areg(3,3) - Areg(2,3)*Areg(3,1)) / det
    Ainv(2,2) =  (Areg(1,1)*Areg(3,3) - Areg(1,3)*Areg(3,1)) / det
    Ainv(2,3) = -(Areg(1,1)*Areg(2,3) - Areg(1,3)*Areg(2,1)) / det
    Ainv(3,1) =  (Areg(2,1)*Areg(3,2) - Areg(2,2)*Areg(3,1)) / det
    Ainv(3,2) = -(Areg(1,1)*Areg(3,2) - Areg(1,2)*Areg(3,1)) / det
    Ainv(3,3) =  (Areg(1,1)*Areg(2,2) - Areg(1,2)*Areg(2,1)) / det

    info = 0_i4
    call dummy_use(this)
  end subroutine invert_3x3


  subroutine solve_3x3_mat(A, B, X, info)
    real(wp), intent(in)  :: A(3,3), B(3,3)
    real(wp), intent(out) :: X(3,3)
    integer(i4), intent(out)  :: info

    real(wp), parameter :: eps_reg = 1.0e-12_wp
    real(wp) :: Areg(3,3), invA(3,3), det

    Areg = A
    det  = det3(Areg)

    if (.not. ieee_is_finite(det) .or. abs(det) < eps_reg) then
      Areg(1,1)=Areg(1,1)+eps_reg
      Areg(2,2)=Areg(2,2)+eps_reg
      Areg(3,3)=Areg(3,3)+eps_reg
      det = det3(Areg)
    end if

    if (.not. ieee_is_finite(det) .or. abs(det) < eps_reg) then
      info = 1_i4
      X    = 0.0_wp
      return
    end if

    invA(1,1) =  (Areg(2,2)*Areg(3,3) - Areg(2,3)*Areg(3,2)) / det
    invA(1,2) = -(Areg(1,2)*Areg(3,3) - Areg(1,3)*Areg(3,2)) / det
    invA(1,3) =  (Areg(1,2)*Areg(2,3) - Areg(1,3)*Areg(2,2)) / det
    invA(2,1) = -(Areg(2,1)*Areg(3,3) - Areg(2,3)*Areg(3,1)) / det
    invA(2,2) =  (Areg(1,1)*Areg(3,3) - Areg(1,3)*Areg(3,1)) / det
    invA(2,3) = -(Areg(1,1)*Areg(2,3) - Areg(1,3)*Areg(2,1)) / det
    invA(3,1) =  (Areg(2,1)*Areg(3,2) - Areg(2,2)*Areg(3,1)) / det
    invA(3,2) = -(Areg(1,1)*Areg(3,2) - Areg(1,2)*Areg(3,1)) / det
    invA(3,3) =  (Areg(1,1)*Areg(2,2) - Areg(1,2)*Areg(2,1)) / det

    X = matmul(invA, B)
    info = 0_i4
  end subroutine solve_3x3_mat


  subroutine solve_3x3_vec(A, b, x, info)
    real(wp), intent(in)  :: A(3,3), b(3)
    real(wp), intent(out) :: x(3)
    integer(i4), intent(out)  :: info

    real(wp), parameter :: eps_reg = 1.0e-12_wp
    real(wp) :: Areg(3,3), invA(3,3), det

    Areg = A
    det  = det3(Areg)

    if (.not. ieee_is_finite(det) .or. abs(det) < eps_reg) then
      Areg(1,1)=Areg(1,1)+eps_reg
      Areg(2,2)=Areg(2,2)+eps_reg
      Areg(3,3)=Areg(3,3)+eps_reg
      det = det3(Areg)
    end if

    if (.not. ieee_is_finite(det) .or. abs(det) < eps_reg) then
      info = 1_i4
      x    = 0.0_wp
      return
    end if

    invA(1,1) =  (Areg(2,2)*Areg(3,3) - Areg(2,3)*Areg(3,2)) / det
    invA(1,2) = -(Areg(1,2)*Areg(3,3) - Areg(1,3)*Areg(3,2)) / det
    invA(1,3) =  (Areg(1,2)*Areg(2,3) - Areg(1,3)*Areg(2,2)) / det
    invA(2,1) = -(Areg(2,1)*Areg(3,3) - Areg(2,3)*Areg(3,1)) / det
    invA(2,2) =  (Areg(1,1)*Areg(3,3) - Areg(1,3)*Areg(3,1)) / det
    invA(2,3) = -(Areg(1,1)*Areg(2,3) - Areg(1,3)*Areg(2,1)) / det
    invA(3,1) =  (Areg(2,1)*Areg(3,2) - Areg(2,2)*Areg(3,1)) / det
    invA(3,2) = -(Areg(1,1)*Areg(3,2) - Areg(1,2)*Areg(3,1)) / det
    invA(3,3) =  (Areg(1,1)*Areg(2,2) - Areg(1,2)*Areg(2,1)) / det

    x = matmul(invA, b)
    info = 0_i4
  end subroutine solve_3x3_vec


  pure real(wp) function det3(A) result(det)
    real(wp), intent(in) :: A(3,3)
    det = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
  end function det3


  subroutine dummy_use(this)
    !! Utility to mark polymorphic dummy arguments as "used" without I/O.
    class(time_integrator_t), intent(in) :: this
    select type (this)
    type is (time_integrator_t)
      continue
    class default
      continue
    end select
  end subroutine dummy_use


  subroutine dummy_use_real(x)
    real(wp), intent(in) :: x
    if (.false.) write(*,*) x
  end subroutine dummy_use_real

end module module_time_step
