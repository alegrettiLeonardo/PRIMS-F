module module_find
  !! Conservative recovery with exact mass/momentum identities and explicit
  !! failures. The regime index (0=stratified, 1=drift flux) is NOT a status.
  !! EOS: rho_g=P/Cg**2, rho_l=rho_l0+(P-P_l0)/Cl**2.
  use mod_kinds, only: wp, i4
  use module_closure_law, only: lei_fechamento_or_ndim_simp, CdUd_bendiksen
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: compute_primitive_variables, pressure_fun_uj, alpha_fun_uj, set_stability_params
  public :: PRIMITIVE_OK, PRIMITIVE_DOMAIN, PRIMITIVE_EOS, PRIMITIVE_CLOSURE
  public :: PRIMITIVE_NO_ROOT, PRIMITIVE_SINGULAR, PRIMITIVE_UNSUPPORTED_SWITCH
  integer(i4), parameter :: PRIMITIVE_OK=0, PRIMITIVE_DOMAIN=1, PRIMITIVE_EOS=2
  integer(i4), parameter :: PRIMITIVE_CLOSURE=3, PRIMITIVE_NO_ROOT=4, PRIMITIVE_SINGULAR=5
  integer(i4), parameter :: PRIMITIVE_UNSUPPORTED_SWITCH=6
  logical, save :: requested_unqualified_switch=.false.
contains
  subroutine set_stability_params(si_on_in,si_off_in,enable_hysteresis,enable_switch)
    !! Retained source compatibility. The reference model fixes the pipe
    !! closure. An unqualified empirical switching criterion is not activated
    !! by changing its hysteresis thresholds. Explicit enable_switch returns
    !! status 6 at recovery until a traceable criterion is supplied.
    real(wp), intent(in), optional :: si_on_in,si_off_in
    logical, intent(in), optional :: enable_hysteresis,enable_switch
    if (present(enable_switch)) requested_unqualified_switch=enable_switch
  end subroutine set_stability_params

  pure real(wp) function pressure_fun_uj(u1,u2,Cg,Cl,rho_l0,P_l0) result(P)
    real(wp), intent(in) :: u1,u2,Cg,Cl,rho_l0,P_l0
    real(wp) :: a,b,q,disc,y
    P=ieee_value(0.0_wp,ieee_quiet_nan)
    if (.not. all(ieee_is_finite([u1,u2,Cg,Cl,rho_l0,P_l0]))) return
    if (min(u1,u2)<0.0_wp .or. u1+u2<=0.0_wp .or. min(Cg,Cl,rho_l0)<=0.0_wp) return
    a=rho_l0-P_l0/(Cl*Cl)
    if (a<=0.0_wp) return
    q=u2*(Cg/Cl)**2
    b=a-u1-q
    disc=hypot(b,2.0_wp*sqrt(a)*sqrt(q))
    if (b>=0.0_wp) then
      if (b+disc<=0.0_wp) return
      y=2.0_wp*a*(q/(b+disc))
    else
      y=0.5_wp*(disc-b)
    end if
    if (y<=0.0_wp) return
    P=y*Cl*Cl
  end function pressure_fun_uj

  pure real(wp) function alpha_fun_uj(u1,u2,Cg,Cl,rho_l0,P_l0,tol) result(alpha)
    !! tol is retained for compatibility; a small positive gas mass is never
    !! silently replaced by a gas-free state.
    real(wp), intent(in) :: u1,u2,Cg,Cl,rho_l0,P_l0,tol
    real(wp) :: P
    alpha=ieee_value(0.0_wp,ieee_quiet_nan)
    P=pressure_fun_uj(u1,u2,Cg,Cl,rho_l0,P_l0)
    if (.not. ieee_is_finite(P)) return
    if (u2==0.0_wp) then
      alpha=0.0_wp
    else if (u1==0.0_wp) then
      alpha=1.0_wp
    else
      alpha=u2/(P/(Cg*Cg))
      if (alpha<=0.0_wp .or. alpha>=1.0_wp) alpha=ieee_value(0.0_wp,ieee_quiet_nan)
    end if
  end function alpha_fun_uj

  subroutine compute_primitive_variables(u1,u2,u3,theta,Cg,Cl,rho_l0,P_l0,D,AREA,EPS,G,MUL,MUG,sigma,njl_ref, &
                  w_u,w_rho,tol,tola,maxit,alpha,rhol,rhog,P,ul,ug,index,idx_prev,ul_prev,ug_prev, &
                  status,fixed_regime,closure_id)
    real(wp), intent(in) :: u1,u2,u3,theta,Cg,Cl,rho_l0,P_l0,D,AREA,EPS,G,MUL,MUG,sigma,njl_ref
    real(wp), intent(in) :: w_u,w_rho,tol,tola
    integer(i4), intent(in) :: maxit
    real(wp), intent(out) :: alpha,rhol,rhog,P,ul,ug
    integer(i4), intent(out) :: index
    integer(i4), intent(in), optional :: idx_prev,fixed_regime,closure_id
    real(wp), intent(in), optional :: ul_prev,ug_prev
    integer(i4), intent(out), optional :: status
    real(wp) :: nan,x,fx,xlo,xhi,flo,fhi,width,h,fp,fm,deriv,xnew,dx,xtol,ftol,scale
    real(wp) :: c0,ud,det,ds,ul_try,ug_try,fr,vol,ul_seed,ug_seed,velocity_upper,hplus,hminus,endpoint_tol
    integer(i4) :: cid,code,it,expansion
    logical :: solve_gas,bracketed
    nan=ieee_value(0.0_wp,ieee_quiet_nan)
    alpha=nan; rhol=nan; rhog=nan; P=nan; ul=nan; ug=nan
    index=-1
    if (present(status)) status=PRIMITIVE_DOMAIN
    if (.not. all(ieee_is_finite([u1,u2,u3,theta,Cg,Cl,rho_l0,P_l0,w_u,w_rho,tol,tola]))) return
    if (min(u1,u2)<0.0_wp .or. u1+u2<=0.0_wp .or. min(Cg,Cl,w_u,w_rho)<=0.0_wp) return
    if (maxit<=0) return
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
    if (cid<1 .or. cid>4) return
    index=1
    if (cid==1) index=0
    if (requested_unqualified_switch .and. theta<0.0_wp .and. &
        .not. present(fixed_regime) .and. .not. present(closure_id)) then
      if (present(status)) status=PRIMITIVE_UNSUPPORTED_SWITCH
      return
    end if
    if (present(status)) status=PRIMITIVE_EOS
    P=pressure_fun_uj(u1,u2,Cg,Cl,rho_l0,P_l0)
    if (.not. ieee_is_finite(P)) return
    rhog=P/(Cg*Cg)
    rhol=rho_l0+(P-P_l0)/(Cl*Cl)
    if (.not. all(ieee_is_finite([rhog,rhol]))) return
    if (min(rhol,rhog)<=0.0_wp) return
    if (u2==0.0_wp) then
      alpha=0.0_wp; ul=u3/u1; ug=0.0_wp
      if (present(status)) status=PRIMITIVE_OK
      return
    else if (u1==0.0_wp) then
      alpha=1.0_wp; ul=0.0_wp; ug=u3/u2
      if (present(status)) status=PRIMITIVE_OK
      return
    end if
    alpha=u2/rhog
    vol=u1/rhol+alpha
    if (alpha<=0.0_wp .or. alpha>=1.0_wp) return
    if (abs(vol-1.0_wp)>1024.0_wp*epsilon(1.0_wp)) return
    if (present(status)) status=PRIMITIVE_CLOSURE
    ul_seed=u3/(u1+u2); ug_seed=ul_seed
    if (present(ul_prev)) then
      if (ieee_is_finite(ul_prev)) ul_seed=ul_prev
    end if
    if (present(ug_prev)) then
      if (ieee_is_finite(ug_prev)) ug_seed=ug_prev
    end if

    ! Bendiksen is linear in velocities within each Froude branch. Recheck
    ! the branch AFTER solving, including negative liquid velocity. No scan.
    if (cid==2) then
      do it=1,4
        call CdUd_bendiksen(alpha,rhol,rhog,ul_seed,ug_seed,theta,D,AREA,EPS,G,MUL,MUG,sigma, &
                           w_u,w_rho,tol,c0,ud,status=code)
        if (code/=0) return
        det=u1*(1.0_wp-c0*alpha)+u2*c0*(1.0_wp-alpha)
        ds=abs(u1*(1.0_wp-c0*alpha))+abs(u2*c0*(1.0_wp-alpha))
        if (abs(det)<=64.0_wp*epsilon(1.0_wp)*ds) then
          if (present(status)) status=PRIMITIVE_SINGULAR
          return
        end if
        ul_try=((1.0_wp-c0*alpha)*u3-u2*ud/w_u)/det
        ug_try=(u1*ud/w_u+c0*(1.0_wp-alpha)*u3)/det
        call lei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul_try,ug_try,theta,D,AREA,EPS,G,MUL,MUG,sigma, &
                                        w_u,w_rho,tol,fr,closure_id=cid,status=code)
        if (code/=0) return
        if (abs(fr)<=128.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(ul_try),abs(ug_try))) then
          ul=ul_try; ug=ug_try
          if (present(status)) status=PRIMITIVE_OK
          return
        end if
        ul_seed=ul_try; ug_seed=ug_try
      end do
      ! A discontinuous Froude branch can have no admissible algebraic root;
      ! the general solver below must satisfy the residual, not just a jump.
    end if

    solve_gas=(u1>=u2)
    if (cid==3) then
      ! The verified Bhagwat transcription is co-current only. Bound the
      ! SOLVER variable, not U: momentum then gives the other speed >=0.
      if (u3<=0.0_wp) return
      if (solve_gas) then
        velocity_upper=u3/u2
      else
        velocity_upper=u3/u1
      end if
    end if
    if (solve_gas) then
      x=ug_seed
    else
      x=ul_seed
    end if
    if (cid==3) x=min(max(x,0.0_wp),velocity_upper)
    call residual(x,fx,code)
    if (code/=0) return
    if (fx==0.0_wp) then
      call accept_velocity(x)
      return
    end if
    width=max(0.25_wp,0.5_wp*abs(x))
    if (G>=0.0_wp .and. D>0.0_wp) width=max(width,sqrt(G*D)/w_u)
    bracketed=.false.
    do expansion=1,40
      if (cid==3) then
        ! Keep the lower-mass velocity as the independent unknown to avoid
        ! cancellation in its recovery. Expand only as far as required to
        ! bracket the root: the full co-current endpoint may be enormous
        ! when one phase has a very small positive mass.
        xlo=max(0.0_wp,x-width); xhi=min(velocity_upper,x+width)
      else
        xlo=x-width; xhi=x+width
      end if
      call residual(xlo,flo,code)
      if (code/=0) return
      call residual(xhi,fhi,code)
      if (code/=0) return
      if (cid==3) then
        endpoint_tol=128.0_wp*epsilon(1.0_wp)*alpha*max(1.0_wp,abs(xlo),abs(xhi),u3/max(u1,u2))
        if (abs(flo)<=endpoint_tol) then
          call accept_velocity(xlo)
          return
        end if
        if (abs(fhi)<=endpoint_tol) then
          call accept_velocity(xhi)
          return
        end if
      end if
      if (flo==0.0_wp) then
        call accept_velocity(xlo)
        return
      end if
      if (fhi==0.0_wp) then
        call accept_velocity(xhi)
        return
      end if
      if (sign(1.0_wp,flo)/=sign(1.0_wp,fhi)) then
        bracketed=.true.
        exit
      end if
      if (cid==3) then
        if (xlo==0.0_wp .and. xhi==velocity_upper) exit
      end if
      width=2.0_wp*width
    end do
    if (present(status)) status=PRIMITIVE_NO_ROOT
    if (.not. bracketed) return
    scale=max(abs(flo),abs(fhi),tiny(1.0_wp))
    ftol=max(128.0_wp*epsilon(1.0_wp),min(abs(tol),1.0e-11_wp))*scale
    do it=1,min(maxit,100_i4)
      if (fx==0.0_wp) then
        call accept_velocity(x)
        return
      end if
      if (sign(1.0_wp,fx)==sign(1.0_wp,flo)) then
        xlo=x; flo=fx
      else
        xhi=x; fhi=fx
      end if
      xtol=max(128.0_wp*epsilon(1.0_wp),min(abs(tola),1.0e-10_wp))*(1.0_wp+abs(x))
      h=epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(1.0_wp,abs(x))
      hplus=h; hminus=h
      if (cid==3) then
        hplus=min(h,velocity_upper-x); hminus=min(h,x)
      end if
      fp=fx; fm=fx
      if (hplus>0.0_wp) then
        call residual(x+hplus,fp,code)
        if (code/=0) return
      end if
      if (hminus>0.0_wp) then
        call residual(x-hminus,fm,code)
        if (code/=0) return
      end if
      if (hplus+hminus<=0.0_wp) return
      deriv=(fp-fm)/(hplus+hminus)
      xnew=0.5_wp*(xlo+xhi)
      if (abs(deriv)>tiny(1.0_wp)) then
        dx=fx/deriv
        if (abs(dx)<=xtol .and. abs(fx)<=ftol) then
          call accept_velocity(x)
          return
        end if
        if (x-dx>xlo .and. x-dx<xhi) xnew=x-dx
      end if
      if (xhi-xlo<=xtol) then
        if (abs(fx)<=ftol) call accept_velocity(x)
        return
      end if
      x=xnew
      call residual(x,fx,code)
      if (code/=0) return
    end do
  contains
    subroutine residual(v,f,istat)
      real(wp), intent(in) :: v
      real(wp), intent(out) :: f
      integer(i4), intent(out) :: istat
      real(wp) :: vl,vg
      if (solve_gas) then
        vg=v; vl=(u3-u2*v)/u1
      else
        vl=v; vg=(u3-u1*v)/u2
      end if
      if (cid==3) then
        ! The upper endpoint means the other phase is stationary exactly;
        ! avoid subtractive roundoff inventing a negative boundary speed.
        if (v==velocity_upper) then
          if (solve_gas) then
            vl=0.0_wp
          else
            vg=0.0_wp
          end if
        end if
      end if
      call lei_fechamento_or_ndim_simp(alpha,rhol,rhog,vl,vg,theta,D,AREA,EPS,G,MUL,MUG,sigma, &
                                      w_u,w_rho,tol,f,closure_id=cid,status=istat)
      if (.not. ieee_is_finite(f)) istat=PRIMITIVE_CLOSURE
    end subroutine residual
    subroutine accept_velocity(v)
      real(wp), intent(in) :: v
      if (solve_gas) then
        ug=v; ul=(u3-u2*v)/u1
      else
        ul=v; ug=(u3-u1*v)/u2
      end if
      if (cid==3) then
        if (v==velocity_upper) then
          if (solve_gas) then
            ul=0.0_wp
          else
            ug=0.0_wp
          end if
        end if
      end if
      if (.not. all(ieee_is_finite([ul,ug]))) then
        ul=nan; ug=nan
        if (present(status)) status=PRIMITIVE_NO_ROOT
      else
        if (present(status)) status=PRIMITIVE_OK
      end if
    end subroutine accept_velocity
  end subroutine compute_primitive_variables
end module module_find
