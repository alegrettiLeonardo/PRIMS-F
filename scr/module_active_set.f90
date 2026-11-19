module module_active_set
  !! A single active linear inequality in a scaled Newton L2 or infinity-norm model.
  !! The nonlinear residual and all physical acceptance tests stay with caller.
  use mod_kinds, only: wp
  use module_time_step, only: time_integrator_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: active_set_direction
contains
  subroutine active_set_direction(a,b,c,scale,gradient,margin,direction,active,status,minimax,solve_count)
    real(wp), intent(in) :: a(:,:,:),b(:,:,:),c(:,:,:),scale(:,:),gradient(:,:),margin
    real(wp), intent(inout) :: direction(:,:)
    logical, intent(out) :: active
    logical, intent(in), optional :: minimax
    integer, intent(out) :: status
    integer, intent(out), optional :: solve_count
    type(time_integrator_t) :: solver
    real(wp), allocatable :: sa(:,:,:),sb(:,:,:),sc(:,:,:),ta(:,:,:),tb(:,:,:),tc(:,:,:)
    real(wp), allocatable :: normal(:,:),dual(:,:),correction(:,:)
    real(wp) :: margin_scale,predicted,denominator,multiplier,dual_scale,term_scale,term_sum,increment
    integer :: n,i,j,k
    logical :: use_minimax

    if(present(solve_count))solve_count=0
    use_minimax=.false.
    if(present(minimax))use_minimax=minimax
    active=.false.
    status=-100
    n=size(scale,2)
    if(size(scale,1)/=3.or.n<1)return
    if(any(shape(direction)/=shape(scale)).or.any(shape(gradient)/=shape(scale)))return
    if(any(shape(a)/=[3,3,n]).or.any(shape(b)/=[3,3,n]).or.any(shape(c)/=[3,3,n]))return
    if(.not.all(ieee_is_finite(scale)))return
    if(any(scale<=0.0_wp))return
    if(.not.all(ieee_is_finite(gradient)))return
    if(.not.ieee_is_finite(margin))return
    if(.not.all(ieee_is_finite(direction)))return
    if(margin<0.0_wp)return
    if(.not.all(ieee_is_finite(a)))return
    if(.not.all(ieee_is_finite(b)))return
    if(.not.all(ieee_is_finite(c)))return
    do k=1,n
      do i=1,3
        if(abs(gradient(i,k))>1.0_wp)then
          if(scale(i,k)>huge(1.0_wp)/abs(gradient(i,k)))then
            status=-105
            return
          endif
        endif
        if(scale(i,k)<1.0_wp)then
          if(abs(direction(i,k))>huge(1.0_wp)*scale(i,k))then
            status=-105
            return
          endif
        endif
      enddo
    enddo
    margin_scale=max(abs(margin),maxval(abs(gradient*scale)))
    if(margin_scale==0.0_wp)then
      status=0
      return
    endif
    allocate(normal(3,n),dual(3,n),correction(3,n))
    normal=gradient*scale/margin_scale
    correction=normal*(direction/scale)
    term_scale=maxval(abs(correction))
    predicted=margin/margin_scale
    if(term_scale>0.0_wp)then
      term_sum=sum(correction/term_scale)
      if(term_scale>1.0_wp)then
        if(abs(term_sum)>huge(1.0_wp)/term_scale)then
          status=-105
          return
        endif
      endif
      predicted=predicted+term_sum*term_scale
    endif
    status=0
    if(predicted>=0.0_wp)return

    allocate(sa(3,3,n),sb(3,3,n),sc(3,3,n),ta(3,3,n),tb(3,3,n),tc(3,3,n))
    sa=0.0_wp
    sc=0.0_wp
    do k=1,n
      do j=1,3
        do i=1,3
          sb(i,j,k)=rescale_coefficient(b(i,j,k),scale(j,k),scale(i,k),status)
          if(status/=0)return
          if(k>1)then
            sa(i,j,k)=rescale_coefficient(a(i,j,k),scale(j,k-1),scale(i,k),status)
            if(status/=0)return
          endif
          if(k<n)then
            sc(i,j,k)=rescale_coefficient(c(i,j,k),scale(j,k+1),scale(i,k),status)
            if(status/=0)return
          endif
        enddo
      enddo
    enddo
    ta=0.0_wp
    tc=0.0_wp
    do k=1,n
      tb(:,:,k)=transpose(sb(:,:,k))
      if(k>1)ta(:,:,k)=transpose(sc(:,:,k-1))
      if(k<n)tc(:,:,k)=transpose(sa(:,:,k+1))
    enddo
    ! A = W J S, b = W R, z0 = S^-1 d_Newton, c = grad(g) S/gscale.
    ! q=A^-T c; z=z0+A^-1 q*(-g/gscale-c.z0)/(q.q).
    ! This solves the KKT equations without forming A^T A.
    if(present(solve_count))solve_count=1
    call solver%solve_block_tridiagonal_3x3(ta,tb,tc,normal,dual,status)
    if(status/=0)return
    ! Normalize before the dot product: a finite near-singular solve can
    ! produce q~1e200 while its squared norm is not representable.
    dual_scale=maxval(abs(dual))
    if(.not.ieee_is_finite(dual_scale))then
      status=-101
      return
    endif
    if(dual_scale==0.0_wp)then
      status=-101
      return
    endif
    dual=dual/dual_scale
    denominator=sum(dual*dual)
    if(dual_scale<1.0_wp)then
      if(abs(predicted)>huge(1.0_wp)*dual_scale)then
        status=-103
        return
      endif
    endif
    if(use_minimax)then
      ! Holder: q.r <= ||q||_1 ||r||_infinity. Equality is attained by
      ! r_i=delta*sign(q_i)/||q||_1, with zero entries where q_i=0.
      ! Normalize q first so both its L1 and squared norms are finite.
      multiplier=(-predicted/dual_scale)/sum(abs(dual))
      if(multiplier==0.0_wp)then
        status=-104
        return
      endif
      where(dual/=0.0_wp)
        dual=sign(multiplier,dual)
      elsewhere
        dual=0.0_wp
      endwhere
    else
      multiplier=(-predicted/dual_scale)/denominator
      if(multiplier==0.0_wp)then
        status=-104
        return
      endif
      dual=dual*multiplier
    endif
    if(present(solve_count))solve_count=2
    call solver%solve_block_tridiagonal_3x3(sa,sb,sc,dual,correction,status)
    if(status/=0)return
    do k=1,n
      do i=1,3
        increment=rescale_coefficient(correction(i,k),scale(i,k),1.0_wp,status)
        if(status/=0)return
        if((direction(i,k)>0.0_wp).eqv.(increment>0.0_wp))then
          if(abs(direction(i,k))>huge(1.0_wp)-abs(increment))then
            status=-107
            return
          endif
        endif
        correction(i,k)=direction(i,k)+increment
      enddo
    enddo
    direction=correction
    active=.true.
  end subroutine active_set_direction
  function rescale_coefficient(value,column_scale,row_scale,status) result(result_value)
    !! Evaluate value*column_scale/row_scale without an overflowing numerator.
    real(wp),intent(in) :: value,column_scale,row_scale
    integer,intent(out) :: status
    real(wp) :: result_value,mantissa
    integer :: power
    status=0
    result_value=0.0_wp
    if(value==0.0_wp)return
    mantissa=fraction(value)*fraction(column_scale)/fraction(row_scale)
    power=exponent(value)+exponent(column_scale)-exponent(row_scale)+exponent(mantissa)
    mantissa=fraction(mantissa)
    if(power>maxexponent(1.0_wp))then
      status=-106
      return
    endif
    if(power==maxexponent(1.0_wp))then
      if(abs(mantissa)>fraction(huge(1.0_wp)))then
        status=-106
        return
      endif
    endif
    if(power<minexponent(1.0_wp)-digits(1.0_wp))then
      status=-106
      return
    endif
    result_value=scale(mantissa,power)
    if(result_value==0.0_wp)status=-106
  end function rescale_coefficient
end module module_active_set
