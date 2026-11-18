module module_jacobian
  !! Flux derivative on ONE selected closure branch. The EOS equations
  !! determine d(alpha,P)/dU; momentum and the differential of that same
  !! closure determine d(ul,ug)/dU. No denominator clipping is applied.
  use mod_kinds, only: wp, i4
  use module_closure_law, only: Dlei_fechamento_or_ndim_simp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: DFduj
contains
  function DFduj(alpha,rhol,rhog,ul,ug,Cl,Cg,theta,DH,AREA,EPS,G,MUL,MUG,sigma, &
                w_u,w_rho,w_P,tol,status,closure_id) result(MJF)
    real(wp), intent(in) :: alpha,rhol,rhog,ul,ug,Cl,Cg,theta,DH,AREA,EPS,G,MUL,MUG,sigma
    real(wp), intent(in) :: w_u,w_rho,w_P,tol
    integer(i4), intent(out), optional :: status
    integer(i4), intent(in), optional :: closure_id
    real(wp) :: MJF(3,3)
    real(wp) :: fa,frg,frl,fug,ful,fp,ml,mg,dl,dg,eosdet,vdet,vscale,pip
    real(wp) :: da(3),dp(3),dul,dug,rhs_m,rhs_f,dml,dmg
    integer(i4) :: code,j
    MJF=ieee_value(0.0_wp,ieee_quiet_nan)
    if (present(status)) status=1
    if (.not. all(ieee_is_finite([alpha,rhol,rhog,ul,ug,Cl,Cg,w_u,w_rho,w_P]))) return
    if (alpha<=0.0_wp .or. alpha>=1.0_wp .or. min(rhol,rhog,Cl,Cg,w_u,w_rho,w_P)<=0.0_wp) return
    call Dlei_fechamento_or_ndim_simp(alpha,rhol,rhog,ul,ug,theta,DH,AREA,EPS,G,MUL,MUG,sigma, &
                                    w_u,w_rho,tol,fa,frg,frl,fug,ful,closure_id=closure_id,status=code)
    if (code/=0) then
      if (present(status)) status=10+code
      return
    end if
    if (.not. all(ieee_is_finite([fa,frg,frl,fug,ful]))) return
    ml=(1.0_wp-alpha)*rhol; mg=alpha*rhog
    dl=(1.0_wp-alpha)/(Cl*Cl); dg=alpha/(Cg*Cg)
    eosdet=rhol*dg+rhog*dl
    vdet=ml*fug-mg*ful
    vscale=abs(ml*fug)+abs(mg*ful)
    if (eosdet<=tiny(1.0_wp) .or. abs(vdet)<=64.0_wp*epsilon(1.0_wp)*vscale) then
      if (present(status)) status=2
      return
    end if
    da=[-dg,dl,0.0_wp]/eosdet
    dp=[rhog,rhol,0.0_wp]/eosdet
    fp=frl/(Cl*Cl)+frg/(Cg*Cg)
    pip=w_P/(w_rho*w_u*w_u)
    do j=1,3
      dml=0.0_wp; dmg=0.0_wp; rhs_m=0.0_wp
      if (j==1) dml=1.0_wp
      if (j==2) dmg=1.0_wp
      if (j==3) rhs_m=1.0_wp
      rhs_m=rhs_m-ul*dml-ug*dmg
      rhs_f=-fa*da(j)-fp*dp(j)
      dul=(rhs_m*fug-mg*rhs_f)/vdet
      dug=(ml*rhs_f-ful*rhs_m)/vdet
      MJF(1,j)=ul*dml+ml*dul
      MJF(2,j)=ug*dmg+mg*dug
      MJF(3,j)=ul*ul*dml+ug*ug*dmg+2.0_wp*(ml*ul*dul+mg*ug*dug)+pip*dp(j)
    end do
    if (.not. all(ieee_is_finite(MJF))) then
      if (present(status)) status=3
      return
    end if
    if (present(status)) status=0
  end function DFduj
end module module_jacobian
