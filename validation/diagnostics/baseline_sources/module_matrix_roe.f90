module module_matrix_roe
  !! Fortran translation of matrixRoe.py (Roe matrix via Gauss-Legendre quadrature).
  !!
  !! This module mirrors:
  !!   calculate_roe_matrix(U_i, U_ip1, ... , num_points, a=0, b=1, ul_prev, ug_prev)
  !!
  !! To keep the core algorithm faithful, the quadrature and interpolation
  !! are implemented directly. Caching is omitted (you can add it later).
  use mod_kinds, only: wp, i4
  use module_find, only: compute_primitive_variables
  use module_jacobian, only: DFduj
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private

  public :: calculate_roe_matrix

contains

  subroutine gauss_legendre_points(n, a, b, x, w)
    integer(i4), intent(in) :: n
    real(wp), intent(in) :: a, b
    real(wp), intent(out) :: x(n), w(n)
    real(wp) :: half, mid

    half = 0.5_wp*(b-a)
    mid  = 0.5_wp*(a+b)

    select case(n)
    case(1)
      x(1) = 0.0_wp
      w(1) = 2.0_wp
    case(2)
      x = [ -0.5773502691896257_wp, 0.5773502691896257_wp ]
      w = [  1.0_wp, 1.0_wp ]
    case(3)
      x = [ -0.7745966692414834_wp, 0.0_wp, 0.7745966692414834_wp ]
      w = [  0.5555555555555556_wp, 0.8888888888888888_wp, 0.5555555555555556_wp ]
    case(4)
      x = [ -0.8611363115940526_wp, -0.3399810435848563_wp, 0.3399810435848563_wp, 0.8611363115940526_wp ]
      w = [  0.3478548451374539_wp,  0.6521451548625461_wp, 0.6521451548625461_wp, 0.3478548451374539_wp ]
    case(5)
      x = [ -0.9061798459386640_wp, -0.5384693101056831_wp, 0.0_wp, 0.5384693101056831_wp, 0.9061798459386640_wp ]
      w = [  0.2369268850561891_wp,  0.4786286704993665_wp, 0.5688888888888889_wp, 0.4786286704993665_wp, 0.2369268850561891_wp ]
    case default
      x=ieee_value(0.0_wp,ieee_quiet_nan)
      w=x
      return
    end select

    ! transform [-1,1] -> [a,b]
    x = mid + half * x
    w = half * w
  end subroutine gauss_legendre_points


  function calculate_roe_matrix(U_i,U_ip1,Cl,Cg,theta,rho_l0,P_l0,D,AREA,EPS,G,MUL,MUG,sigma,njl_ref, &
                                w_u,w_rho,w_P,tol,tola,num_points,a,b,ul_prev,ug_prev,status,fixed_regime,closure_id) result(RoeA)
    !! Quadrature of dF/dU on the conservative segment. The closure branch
    !! is selected once and held fixed at every node. This is an approximate
    !! path integral, not an exact algebraic Roe average for a general closure.
    real(wp), intent(in) :: U_i(3),U_ip1(3),Cl,Cg,theta,rho_l0,P_l0
    real(wp), intent(in) :: D,AREA,EPS,G,MUL,MUG,sigma,njl_ref,w_u,w_rho,w_P,tol,tola
    integer(i4), intent(in) :: num_points
    real(wp), intent(in), optional :: a,b,ul_prev,ug_prev
    integer(i4), intent(out), optional :: status
    integer(i4), intent(in), optional :: fixed_regime,closure_id
    real(wp) :: RoeA(3,3),aa,bb,xi(5),wi(5),t,Ux(3),partial(3,3),Jmat(3,3)
    real(wp) :: alpha_i,rho_l_i,rho_g_i,P_i,u_l_i,u_g_i,ul_seed,ug_seed
    integer(i4) :: idx,j,cid,code
    RoeA=ieee_value(0.0_wp,ieee_quiet_nan)
    if (present(status)) status=1
    if (num_points<1 .or. num_points>5) return
    if (.not. all(ieee_is_finite(U_i)) .or. .not. all(ieee_is_finite(U_ip1))) return
    aa=0.0_wp; bb=1.0_wp
    if (present(a)) aa=a
    if (present(b)) bb=b
    if (.not. all(ieee_is_finite([aa,bb]))) return
    if (bb<=aa) return
    cid=2
    if (theta<0.0_wp) cid=1
    if (present(fixed_regime)) then
      select case(fixed_regime)
      case(0)
        cid=1
      case(1)
        cid=2
      case default
        return
      end select
    end if
    if (present(closure_id)) cid=closure_id
    call gauss_legendre_points(num_points,aa,bb,xi(:num_points),wi(:num_points))
    partial=0.0_wp
    ul_seed=0.0_wp; ug_seed=0.0_wp
    if (present(ul_prev)) ul_seed=ul_prev
    if (present(ug_prev)) ug_seed=ug_prev
    do j=1,num_points
      t=(xi(j)-aa)/(bb-aa)
      Ux=U_i+t*(U_ip1-U_i)
      call compute_primitive_variables(Ux(1),Ux(2),Ux(3),theta,Cg,Cl,rho_l0,P_l0,D,AREA,EPS,G, &
            MUL,MUG,sigma,njl_ref,w_u,w_rho,tol,tola,100_i4,alpha_i,rho_l_i,rho_g_i,P_i,u_l_i,u_g_i,idx, &
            ul_prev=ul_seed,ug_prev=ug_seed,status=code,closure_id=cid)
      if (code/=0) then
        if (present(status)) status=10+code
        return
      end if
      ul_seed=u_l_i; ug_seed=u_g_i
      Jmat=DFduj(alpha_i,rho_l_i,rho_g_i,u_l_i,u_g_i,Cl,Cg,theta,D,AREA,EPS,G,MUL,MUG,sigma, &
                 w_u,w_rho,w_P,tol,status=code,closure_id=cid)
      if (code/=0) then
        if (present(status)) status=100+code
        return
      end if
      ! Normalized path weights: arbitrary [a,b] parameterizations must give
      ! the same average as [0,1], including when the interval length is not 1.
      partial=partial+wi(j)/(bb-aa)*Jmat
    end do
    if (.not. all(ieee_is_finite(partial))) return
    RoeA=partial
    if (present(status)) status=0
  end function calculate_roe_matrix
end module module_matrix_roe
