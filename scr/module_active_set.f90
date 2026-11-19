module module_active_set
  !! A single active linear inequality in a scaled Newton least-squares model.
  !! The nonlinear residual and all physical acceptance tests stay with caller.
  use mod_kinds, only: wp
  use module_time_step, only: time_integrator_t
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: active_set_direction
contains
  subroutine active_set_direction(a,b,c,scale,gradient,margin,direction,active,status)
    real(wp), intent(in) :: a(:,:,:),b(:,:,:),c(:,:,:),scale(:,:),gradient(:,:),margin
    real(wp), intent(inout) :: direction(:,:)
    logical, intent(out) :: active
    integer, intent(out) :: status
    type(time_integrator_t) :: solver
    real(wp), allocatable :: sa(:,:,:),sb(:,:,:),sc(:,:,:),ta(:,:,:),tb(:,:,:),tc(:,:,:)
    real(wp), allocatable :: normal(:,:),dual(:,:),correction(:,:)
    real(wp) :: margin_scale,predicted,denominator,multiplier,dual_scale
    integer :: n,i,j,k

    active=.false.
    status=-100
    n=size(scale,2)
    if(size(scale,1)/=3.or.n<1)return
    if(any(shape(direction)/=shape(scale)).or.any(shape(gradient)/=shape(scale)))return
    if(any(shape(a)/=[3,3,n]).or.any(shape(b)/=[3,3,n]).or.any(shape(c)/=[3,3,n]))return
    if(.not.all(ieee_is_finite(scale)).or.any(scale<=0.0_wp))return
    if(.not.all(ieee_is_finite(gradient)).or..not.ieee_is_finite(margin))return
    if(.not.all(ieee_is_finite(direction)).or.margin<0.0_wp)return
    margin_scale=max(abs(margin),maxval(abs(gradient*scale)))
    if(margin_scale<=tiny(1.0_wp))then
      status=0
      return
    endif
    allocate(normal(3,n),dual(3,n),correction(3,n))
    normal=gradient*scale/margin_scale
    predicted=margin/margin_scale+sum(normal*direction/scale)
    status=0
    if(predicted>=0.0_wp)return

    allocate(sa(3,3,n),sb(3,3,n),sc(3,3,n),ta(3,3,n),tb(3,3,n),tc(3,3,n))
    sa=0.0_wp
    sc=0.0_wp
    do k=1,n
      do j=1,3
        do i=1,3
          sb(i,j,k)=b(i,j,k)*scale(j,k)/scale(i,k)
          if(k>1)sa(i,j,k)=a(i,j,k)*scale(j,k-1)/scale(i,k)
          if(k<n)sc(i,j,k)=c(i,j,k)*scale(j,k+1)/scale(i,k)
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
    call solver%solve_block_tridiagonal_3x3(ta,tb,tc,normal,dual,status)
    if(status/=0)return
    ! Normalize before the dot product: a finite near-singular solve can
    ! produce q~1e200 while its squared norm is not representable.
    dual_scale=maxval(abs(dual))
    if(.not.ieee_is_finite(dual_scale).or.dual_scale<=tiny(1.0_wp))then
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
    multiplier=(-predicted/dual_scale)/denominator
    call solver%solve_block_tridiagonal_3x3(sa,sb,sc,dual*multiplier,correction,status)
    if(status/=0)return
    correction=direction+scale*correction
    if(.not.all(ieee_is_finite(correction)))then
      status=-102
      return
    endif
    direction=correction
    active=.true.
  end subroutine active_set_direction
end module module_active_set
