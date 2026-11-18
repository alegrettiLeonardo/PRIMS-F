module module_time_step
  use mod_kinds,only:wp,i4
  use module_lapack_small,only:dgeev_3x3,inv3x3
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  public::time_integrator_t
  type::time_integrator_t
  contains
    procedure::compute_time_step_lax_wendroff
    procedure::solver_ivp
    procedure::compute_residual_vector
    procedure::solve_block_tridiagonal_3x3
    procedure::invert_3x3
  end type
  interface
    subroutine dgbsv(n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
      import wp
      integer,intent(in)::n,kl,ku,nrhs,ldab,ldb
      real(wp),intent(inout)::ab(ldab,*),b(ldb,*)
      integer,intent(out)::ipiv(*),info
    end subroutine
  end interface
contains
  real(wp) function compute_time_step_lax_wendroff(this,CFL,delta_x,Roe_matrix,mode,acoustic_scale,u_material,status) result(dt)
    class(time_integrator_t),intent(in)::this
    real(wp),intent(in)::CFL,delta_x,Roe_matrix(3,3)
    character(*),intent(in),optional::mode
    real(wp),intent(in),optional::acoustic_scale,u_material
    integer,intent(out),optional::status
    real(wp)::wr(3),wi(3),vr(3,3),amax,scale
    integer::info
    dt=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(status))status=-100
    if(.not.all(ieee_is_finite([CFL,delta_x])).or.CFL<=0.0_wp.or.delta_x<=0.0_wp)return
    call dgeev_3x3(Roe_matrix,wr,wi,vr,info)
    if(info/=0)then
      if(present(status))status=info
      return
    endif
    amax=maxval(abs(wr))
    if(maxval(abs(wi))>1e-10_wp*max(1.0_wp,amax))then
      if(present(status))status=-103
      return
    endif
    if(present(acoustic_scale))then
      ! Downscaling wave speeds would cease to be a spectral CFL estimate.
      if(.not.ieee_is_finite(acoustic_scale).or.acoustic_scale<1.0_wp)return
      amax=amax*acoustic_scale
    endif
    if(present(mode))then
      select case(trim(mode))
      case('acoustic','ACOUSTIC','Acoustic')
      case('material','MATERIAL','Material')
        ! Material speed may tighten the estimate, never remove acoustic waves.
        if(.not.present(u_material))return
        if(.not.ieee_is_finite(u_material))return
        amax=max(amax,abs(u_material))
      case default
        return
      end select
    endif
    if(amax==0.0_wp)then
      dt=huge(1.0_wp)
    else
      dt=CFL*delta_x/amax
    endif
    if(.not.ieee_is_finite(dt).or.dt<=0.0_wp)then
      dt=ieee_value(0.0_wp,ieee_quiet_nan);return
    endif
    if(present(status))status=0
  end function
  subroutine solver_ivp(this,U,numerical_flux_values,source_term_values,delta_s,delta_t,tol,U_out,status)
    class(time_integrator_t),intent(in)::this
    real(wp),intent(in)::U(:,:),numerical_flux_values(:,:),source_term_values(:,:),delta_s(:),delta_t,tol
    real(wp),intent(out)::U_out(:,:)
    integer,intent(out),optional::status
    integer::i,n
    U_out=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(status))status=-100
    n=size(U,1)
    if(size(U,2)/=3.or.any(shape(U_out)/=shape(U)).or.n<3)return
    if(size(delta_s)<n.or.size(numerical_flux_values,1)<n-1.or.size(numerical_flux_values,2)/=3)return
    if(any(shape(source_term_values)/=shape(U)))return
    if(.not.ieee_is_finite(delta_t).or.delta_t<0.0_wp)return
    if(.not.all(ieee_is_finite(delta_s(2:n-1))).or.minval(delta_s(2:n-1))<=0.0_wp)return
    U_out=U
    do i=2,n-1
      U_out(i,:)=U(i,:)-delta_t*(numerical_flux_values(i,:)-numerical_flux_values(i-1,:))/delta_s(i)+ &
        delta_t*source_term_values(i,:)
    enddo
    if(.not.all(ieee_is_finite(U_out)))return
    if(present(status))status=0
  end subroutine
  subroutine compute_residual_vector(this,U,U_n,delta_s,flux_faces,source_faces,i_first,i_last,dt,R,res_norm,status)
    class(time_integrator_t),intent(in)::this
    real(wp),intent(in)::U(:,:),U_n(:,:),delta_s(:),flux_faces(:,:),source_faces(:,:),dt
    integer,intent(in)::i_first,i_last
    real(wp),intent(out)::R(:,:),res_norm
    integer,intent(out),optional::status
    integer::i
    R=ieee_value(0.0_wp,ieee_quiet_nan);res_norm=ieee_value(0.0_wp,ieee_quiet_nan)
    if(present(status))status=-100
    if(size(U,2)/=3.or.any(shape(U_n)/=shape(U)).or.any(shape(R)/=shape(U)))return
    if(i_first<2.or.i_last<i_first.or.i_last>size(U,1))return
    if(size(delta_s)<i_last.or.size(flux_faces,1)<i_last.or.size(flux_faces,2)/=3)return
    if(size(source_faces,1)<i_last.or.size(source_faces,2)/=3)return
    if(.not.ieee_is_finite(dt).or.dt<0.0_wp)return
    if(.not.all(ieee_is_finite(delta_s(i_first:i_last))).or.minval(delta_s(i_first:i_last))<=0.0_wp)return
    R=0.0_wp
    do i=i_first,i_last
      R(i,:)=U(i,:)-U_n(i,:)+dt*(flux_faces(i,:)-flux_faces(i-1,:))/delta_s(i)-dt*source_faces(i,:)
    enddo
    if(.not.all(ieee_is_finite(R)))return
    res_norm=norm2(R(i_first:i_last,:))
    if(.not.ieee_is_finite(res_norm))return
    if(present(status))status=0
  end subroutine
  subroutine solve_block_tridiagonal_3x3(this,A,B,C,D,X,info)
    !! Pivoted LAPACK band solve. Full 3x3 off-diagonal blocks imply KL=KU=5.
    !! AB(KL+KU+1+i-j,j)=matrix(i,j); LDAB=2*KL+KU+1=16.
    !! No diagonal epsilon or per-block nonsingularity assumption is introduced.
    class(time_integrator_t),intent(in)::this
    real(wp),intent(in)::A(:,:,:),B(:,:,:),C(:,:,:),D(:,:)
    real(wp),intent(out)::X(:,:)
    integer(i4),intent(out)::info
    integer,parameter::kl=5,ku=5,ldab=16,diag=11
    real(wp),allocatable::ab(:,:),rhs(:,:)
    integer,allocatable::ipiv(:)
    integer::nb,n,k,i,j,row,col
    X=ieee_value(0.0_wp,ieee_quiet_nan);info=-100
    if(size(D,1)/=3.or.any(shape(X)/=shape(D)))return
    if(size(A,1)/=3.or.size(A,2)/=3.or.any(shape(A)/=shape(B)).or.any(shape(B)/=shape(C)))return
    nb=size(D,2)
    if(nb<1.or.size(B,3)/=nb)return
    if(.not.all(ieee_is_finite(A)).or..not.all(ieee_is_finite(B)).or. &
      .not.all(ieee_is_finite(C)).or..not.all(ieee_is_finite(D)))return
    n=3*nb
    allocate(ab(ldab,n),rhs(n,1),ipiv(n))
    ab=0.0_wp
    do k=1,nb
      do j=1,3
        col=3*(k-1)+j
        do i=1,3
          row=3*(k-1)+i
          ab(diag+row-col,col)=B(i,j,k)
          if(k>1)ab(diag+row-(col-3),col-3)=A(i,j,k)
          if(k<nb)ab(diag+row-(col+3),col+3)=C(i,j,k)
        enddo
      enddo
    enddo
    rhs(:,1)=reshape(D,[n])
    call dgbsv(n,kl,ku,1,ab,ldab,ipiv,rhs,n,info)
    if(info/=0)return
    if(.not.all(ieee_is_finite(rhs)))then
      info=-101;return
    endif
    X=reshape(rhs(:,1),shape(X))
  end subroutine
  subroutine invert_3x3(this,A,Ainv,info)
    class(time_integrator_t),intent(in)::this
    real(wp),intent(in)::A(3,3)
    real(wp),intent(out)::Ainv(3,3)
    integer(i4),intent(out)::info
    Ainv=inv3x3(A,info)
  end subroutine
end module module_time_step
