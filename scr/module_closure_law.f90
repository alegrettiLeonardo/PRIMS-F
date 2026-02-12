module module_closure_law
  !! Fortran translation of the core closure laws in closureLaw.py
  !!
  !! Implements:
  !!   - CdUd_swananda / DCdUd_swananda
  !!   - drift_flux_swananda_ndim and derivatives
  !!   - RelEquilLocalPipe_comp and derivatives
  !!   - lei_fechamento_or_ndim_simp and derivatives
  !!
  !! Notes:
  !!   * The derivative d(gamma)/d(alpha) in DRelEquiLocal_pipe_comp is computed
  !!     with the same centered finite-difference used in the Python code.
  use mod_kinds, only: wp
  use mod_void_helpers, only: ffan, Dffan, alpha2gamma
  use, intrinsic :: iso_fortran_env, only: i4 => int32
  implicit none
  private

  public :: CdUd_swananda
  public :: DCdUd_swananda
  public :: drift_flux_swananda_ndim
  public :: Ddrift_flux_swananda_ndim
  public :: RelEquilLocalPipe_comp
  public :: DRelEquiLocal_pipe_comp
  public :: lei_fechamento_or_ndim_simp
  public :: Dlei_fechamento_or_ndim_simp

contains

  pure real(wp) function clamp(x, xmin, xmax) result(y)
    real(wp), intent(in) :: x, xmin, xmax
    y = max(xmin, min(xmax, x))
  end function clamp

  pure real(wp) function safe_sqrt(x) result(y)
    real(wp), intent(in) :: x
    y = sqrt(max(x, 1.0e-30_wp))
  end function safe_sqrt

  pure real(wp) function sign0(x) result(s)
    !! Sign with sign(0)=0
    real(wp), intent(in) :: x
    if (abs(x) < 1.0e-30_wp) then
      s = 0.0_wp
    else
      s = sign(1.0_wp, x)
    end if
  end function sign0

  subroutine CdUd_swananda(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, Cd, Ud)
    real(wp), intent(in)  :: alpha, rhol, rhog, ul, ug
    real(wp), intent(in)  :: theta, DH, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in)  :: w_u, w_rho, tol
    real(wp), intent(out) :: Cd, Ud

    real(wp), parameter :: ALPHA_EPS=1.0e-12_wp, MU_EPS=1.0e-30_wp, JT_EPS=1.0e-12_wp
    real(wp), parameter :: RE_EPS=1.0e-6_wp, F_EPS=1.0e-8_wp
    integer(i4) :: it

    real(wp) :: a, alp, jl, jg, jt, jt_safe
    real(wp) :: Fr, Retp, Retp_abs, beta, x, La
    real(wp) :: df, ftp, b, f_new, arg
    real(wp) :: C1, C2, C3, C4, C01, C0a, C0b, aux
    real(wp) :: c_safe, muL_s, muG_s, den_rho, rhog_s

    a = 1.8_wp
    alp = max(ALPHA_EPS, min(1.0_wp-ALPHA_EPS, alpha))

    muL_s = max(MUL, MU_EPS)
    muG_s = max(MUG, MU_EPS)
    rhog_s = max(rhog, 1.0e-12_wp)

    jl = rhol*(1.0_wp-alp)*ul
    jg = rhog*alp*ug
    jt = jl + jg
    jt_safe = sign0(jt) * max(abs(jt), JT_EPS)

    c_safe = max(abs(cos(theta)), 1.0e-12_wp)
    den_rho = max(rhol-rhog, 1.0e-12_wp*rhol)

    Fr = (w_u/(w_rho*jt_safe)) * sqrt((w_rho*rhol*rhog)/den_rho) * sqrt(a/G*DH*c_safe)

    Retp = (w_rho*DH/muL_s) * (jl + jg*(muL_s/muG_s))
    Retp_abs = max(abs(Retp), RE_EPS)

    ! friction factor: laminar fallback + bounded Colebrook
    if (Retp_abs < 2300.0_wp) then
      ftp = 16.0_wp/Retp_abs
    else
      ftp = 0.01_wp
      do it = 1, 50
        arg = EPS/DH/3.7_wp + 1.256_wp/(Retp_abs*sqrt(max(ftp,F_EPS)))
        arg = max(arg, 1.0e-30_wp)
        b = -4.0_wp*log10(arg)
        f_new = 1.0_wp/(b*b)
        df = abs(f_new - ftp)
        ftp = f_new
        if (df <= tol) exit
      end do
    end if
    ftp = max(ftp, F_EPS)

    beta = max(0.0001_wp, min(10.0_wp, (muG_s/muL_s)/max(ftp,F_EPS)))

    x = (w_rho/w_u) * (jl/(w_rho*rhog_s))
    x = min(x, 20.0_wp)

    La = (w_rho*sigma*DH)/(muL_s*w_u) * sqrt((w_rho*rhol)/(w_rho*rhog_s))

    C1 = 1.0_wp/(1.0_wp + (alp/(1.0_wp-alp))**a)
    C2 = beta/((1.0_wp + 75.0_wp*C1)*sqrt(max(Fr,1.0e-30_wp)))
    C3 = 1.0_wp + 0.7_wp*(1.0_wp-C1)*(1.0_wp-exp(-2.2_wp*sqrt(max(beta,1.0e-30_wp))))
    C4 = (x + (1.0_wp-exp(-18.0_wp*alp)))/max(C3, 1.0e-30_wp)

    C01 = (1.0_wp + C2*C4) / (1.0_wp + C2)

    aux = 1.0_wp + 0.55_wp*(1.0_wp-exp(-0.02_wp*La))*(sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp)

    C0a = C01
    C0b = 1.0_wp + (1.0_wp-C1)*aux

    Cd = max(C0a, C0b)

    Ud = 0.2_wp*sqrt(G*DH)*(sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp)
    Ud = Ud*sqrt(max((rhol-rhog)/max(rhol,1.0e-12_wp), 0.0_wp))*sqrt(c_safe)
    Ud = w_u*Ud
  end subroutine CdUd_swananda

  subroutine DCdUd_swananda(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, &
                           dC0da, dC0drhol, dC0drhog, dC0dul, dC0dug, dUdda, dUddrhol, dUddrhog, dUddul, dUddug)
    !! Python: DCdUd_swananda(...)
    !! Returns analytic derivatives of C0 and Ud.
    real(wp), intent(in)  :: alpha, rhol, rhog, ul, ug
    real(wp), intent(in)  :: theta, DH, AREA, EPS, G, MUL, MUG, sigma
    real(wp), intent(in)  :: w_u, w_rho, tol
    real(wp), intent(out) :: dC0da, dC0drhol, dC0drhog, dC0dul, dC0dug
    real(wp), intent(out) :: dUdda, dUddrhol, dUddrhog, dUddul, dUddug

    real(wp) :: a, alp
    real(wp) :: jl, jg, jt
    real(wp) :: Fr, Retp
    real(wp) :: df, ftp, b, f_new
    real(wp) :: beta, x, La

    real(wp) :: djl_da, djl_drhol, djl_dul
    real(wp) :: djg_da, djg_drhog, djg_dug
    real(wp) :: djt_da, djt_drhol, djt_drhog, djt_dul, djt_dug

    real(wp) :: dFr_da, dFr_drhol, dFr_drhog, dFr_dul, dFr_dug
    real(wp) :: dRetp_da, dRetp_drhol, dRetp_drhog, dRetp_dul, dRetp_dug
    real(wp) :: dftp_dRetp

    real(wp) :: dbeta_da, dbeta_drhol, dbeta_drhog, dbeta_dul, dbeta_dug
    real(wp) :: dx_da, dx_drhol, dx_drhog, dx_dul, dx_dug
    real(wp) :: dLa_da, dLa_drhol, dLa_drhog, dLa_dul, dLa_dug

    real(wp) :: C1, C2, C3, C4, C01, aux, C0a, C0b, Cd
    real(wp) :: dC1_da
    real(wp) :: dC2_da, dC2_drhol, dC2_drhog, dC2_dul, dC2_dug
    real(wp) :: dC3_da
    real(wp) :: dC4_da, dC4_drhol, dC4_drhog, dC4_dul, dC4_dug
    real(wp) :: dC01_da, dC01_drhol, dC01_drhog, dC01_dul, dC01_dug
    real(wp) :: daux_da, daux_drhol, daux_drhog
    real(wp) :: dC0b_da, dC0b_drhol, dC0b_drhog

    real(wp) :: Ud
    real(wp) :: dUd1_da
    real(wp) :: dUd2_drhol, dUd2_drhog

    real(wp) :: c_safe, dc_safe

    ! --- local helpers (match Python safe guards)
    real(wp), parameter :: TOL_A = 1.0e-12_wp

    a = 1.8_wp
    alp = clamp(alpha, TOL_A, 1.0_wp - TOL_A)

    jl = rhol*(1.0_wp - alp)*ul
    jg = rhog*alp*ug
    jt = jl + jg

    djl_da   = -rhol*ul
    djl_drhol= (1.0_wp - alp)*ul
    djl_dul  = rhol*(1.0_wp - alp)

    djg_da   = rhog*ug
    djg_drhog= alp*ug
    djg_dug  = rhog*alp

    djt_da    = djl_da + djg_da
    djt_drhol = djl_drhol
    djt_drhog = djg_drhog
    djt_dul   = djl_dul
    djt_dug   = djg_dug

    c_safe = max(abs(cos(theta)), 1.0e-12_wp)
    dc_safe = 0.0_wp ! theta is treated as constant in this derivative set

    Fr = (w_u/(w_rho*jt))*sqrt((w_rho*rhol*rhog)/(rhol-rhog))*sqrt(a/G*DH*c_safe)

    dFr_da    = Fr*(-djt_da/jt)
    dFr_dul   = Fr*(-djt_dul/jt)
    dFr_dug   = Fr*(-djt_dug/jt)
    dFr_drhol = Fr*(0.5_wp*(1.0_wp/rhol + 1.0_wp/(rhol-rhog)) - djt_drhol/jt)
    dFr_drhog = Fr*(0.5_wp*(1.0_wp/rhog + 1.0_wp/(rhog-rhol)) - djt_drhog/jt)

    Retp = (w_rho*DH/MUL)*(jl + jg*(MUL/MUG))

    dRetp_da    = (w_rho*DH/MUL)*(djl_da + (MUL/MUG)*djg_da)
    dRetp_dul   = (w_rho*DH/MUL)*(djl_dul)
    dRetp_dug   = (w_rho*DH/MUL)*((MUL/MUG)*djg_dug)
    dRetp_drhol = (w_rho*DH/MUL)*(djl_drhol)
    dRetp_drhog = (w_rho*DH/MUL)*((MUL/MUG)*djg_drhog)

    df  = 1.0_wp
    ftp = 1.0_wp
    do while (df > tol)
      b     = -4.0_wp*log10(EPS/DH/3.7_wp + 1.256_wp/(Retp*sqrt(ftp)))
      f_new = 1.0_wp/(b*b)
      df    = abs(f_new - ftp)
      ftp   = f_new
    end do

    ! dftp/dRetp (Python literal)
    dftp_dRetp = -0.20096e11_wp*safe_sqrt(ftp)/(Retp*(safe_sqrt(ftp)*Retp*EPS/DH + 0.46472e10_wp))

    beta = (MUL/MUG) / max(ftp, 1.0e-30_wp) * (MUG/MUL)**2
    beta = clamp(beta, 0.0001_wp, 10.0_wp)

    ! dbeta = -(MUG/MUL)/ftp^2 * dftp
    dbeta_da    = -((MUG/MUL)/(ftp*ftp)) * dftp_dRetp * dRetp_da
    dbeta_dul   = -((MUG/MUL)/(ftp*ftp)) * dftp_dRetp * dRetp_dul
    dbeta_dug   = -((MUG/MUL)/(ftp*ftp)) * dftp_dRetp * dRetp_dug
    dbeta_drhol = -((MUG/MUL)/(ftp*ftp)) * dftp_dRetp * dRetp_drhol
    dbeta_drhog = -((MUG/MUL)/(ftp*ftp)) * dftp_dRetp * dRetp_drhog

    x = (w_rho/w_u)*(jl/(w_rho*rhog))
    x = min(x, 20.0_wp)

    dx_da    = -(w_rho/w_u)*(rhol/(w_rho*rhog))*ul
    dx_drhol =  (w_rho/w_u)*((1.0_wp-alp)/(w_rho*rhog))*ul
    dx_drhog = -(w_rho/w_u)*(jl/(w_rho*rhog*rhog))
    dx_dul   =  (w_rho/w_u)*(rhol*(1.0_wp-alp)/(w_rho*rhog))
    dx_dug   = 0.0_wp

    La = (w_rho*sigma*DH)/(MUL*w_u)*sqrt((w_rho*rhol)/(w_rho*rhog))

    dLa_da    = 0.0_wp
    dLa_dul   = 0.0_wp
    dLa_dug   = 0.0_wp
    dLa_drhol = 0.5_wp*La/rhol
    dLa_drhog = -0.5_wp*La/rhog

    C1 = 1.0_wp/(1.0_wp + (alp/(1.0_wp-alp))**a)
    dC1_da = (-a*(alp/(1.0_wp-alp))**a*(alp - 1.0_wp))/(alp*((alp/(1.0_wp-alp))**a + 1.0_wp)**2)

    C2 = beta/((1.0_wp + 75.0_wp*C1)*safe_sqrt(Fr))

    dC2_da    = dbeta_da/((1.0_wp+75.0_wp*C1)*safe_sqrt(Fr)) + beta*( (-75.0_wp*dC1_da)/((1.0_wp+75.0_wp*C1)**2*safe_sqrt(Fr)) ) + beta*( (-0.5_wp*dFr_da)/( (1.0_wp+75.0_wp*C1)*Fr*safe_sqrt(Fr) ) )
    dC2_dul   = dbeta_dul/((1.0_wp+75.0_wp*C1)*safe_sqrt(Fr)) + beta*( (-0.5_wp*dFr_dul)/( (1.0_wp+75.0_wp*C1)*Fr*safe_sqrt(Fr) ) )
    dC2_dug   = dbeta_dug/((1.0_wp+75.0_wp*C1)*safe_sqrt(Fr)) + beta*( (-0.5_wp*dFr_dug)/( (1.0_wp+75.0_wp*C1)*Fr*safe_sqrt(Fr) ) )
    dC2_drhol = dbeta_drhol/((1.0_wp+75.0_wp*C1)*safe_sqrt(Fr)) + beta*( (-0.5_wp*dFr_drhol)/( (1.0_wp+75.0_wp*C1)*Fr*safe_sqrt(Fr) ) )
    dC2_drhog = dbeta_drhog/((1.0_wp+75.0_wp*C1)*safe_sqrt(Fr)) + beta*( (-0.5_wp*dFr_drhog)/( (1.0_wp+75.0_wp*C1)*Fr*safe_sqrt(Fr) ) )

    C3 = 1.0_wp + 0.7_wp*(1.0_wp-C1)*(1.0_wp-exp(-2.2_wp*safe_sqrt(beta)))
    dC3_da = 0.7_wp*(-dC1_da)*(1.0_wp-exp(-2.2_wp*safe_sqrt(beta)))

    C4 = (x + (1.0_wp-exp(-18.0_wp*alp)))/C3

    dC4_da    = (dx_da + 18.0_wp*exp(-18.0_wp*alp))/C3 + (x + 1.0_wp-exp(-18.0_wp*alp))*(-dC3_da)/(C3*C3)
    dC4_dul   = (dx_dul)/C3
    dC4_dug   = (dx_dug)/C3
    dC4_drhol = (dx_drhol)/C3
    dC4_drhog = (dx_drhog)/C3

    C01 = (1.0_wp + C2*C4) / (1.0_wp + C2)

    dC01_da    = ( (dC2_da*C4 + C2*dC4_da)*(1.0_wp+C2) - (1.0_wp + C2*C4)*dC2_da ) / ((1.0_wp+C2)**2)
    dC01_dul   = ( (dC2_dul*C4 + C2*dC4_dul)*(1.0_wp+C2) - (1.0_wp + C2*C4)*dC2_dul ) / ((1.0_wp+C2)**2)
    dC01_dug   = ( (dC2_dug*C4 + C2*dC4_dug)*(1.0_wp+C2) - (1.0_wp + C2*C4)*dC2_dug ) / ((1.0_wp+C2)**2)
    dC01_drhol = ( (dC2_drhol*C4 + C2*dC4_drhol)*(1.0_wp+C2) - (1.0_wp + C2*C4)*dC2_drhol ) / ((1.0_wp+C2)**2)
    dC01_drhog = ( (dC2_drhog*C4 + C2*dC4_drhog)*(1.0_wp+C2) - (1.0_wp + C2*C4)*dC2_drhog ) / ((1.0_wp+C2)**2)

    aux = 1.0_wp + 0.55_wp*(1.0_wp-exp(-0.02_wp*La))*(safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp)

    daux_da    = 0.55_wp*(1.0_wp-exp(-0.02_wp*La))*( -9.0_wp*exp(-18.0_wp*alp) / safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) )
    daux_drhol = 0.55_wp*( (0.02_wp*exp(-0.02_wp*La)*dLa_drhol) )*(safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp)
    daux_drhog = 0.55_wp*( (0.02_wp*exp(-0.02_wp*La)*dLa_drhog) )*(safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp)

    C0a = C01
    C0b = 1.0_wp + (1.0_wp-C1)*aux

    dC0b_da    = (-dC1_da)*aux + (1.0_wp-C1)*daux_da
    dC0b_drhol = (1.0_wp-C1)*daux_drhol
    dC0b_drhog = (1.0_wp-C1)*daux_drhog

    Cd = max(C0a, C0b)

    if (C0a >= C0b) then
      dC0da    = dC01_da
      dC0drhol = dC01_drhol
      dC0drhog = dC01_drhog
      dC0dul   = dC01_dul
      dC0dug   = dC01_dug
    else
      dC0da    = dC0b_da
      dC0drhol = dC0b_drhol
      dC0drhog = dC0b_drhog
      dC0dul   = 0.0_wp
      dC0dug   = 0.0_wp
    end if

    Ud = 0.2_wp*safe_sqrt(G*DH)*(safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp)
    Ud = Ud*safe_sqrt((rhol - rhog)/rhol)*safe_sqrt(c_safe)
    Ud = w_u*Ud

    dUd1_da = 0.2_wp*safe_sqrt(G*DH)*( -9.0_wp*exp(-18.0_wp*alp)/safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) )

    ! term2 = sqrt((rhol-rhog)/rhol)
    dUd2_drhol = 0.5_wp*safe_sqrt((rhol-rhog)/rhol)*( (1.0_wp/(rhol-rhog) - 1.0_wp/rhol) )
    dUd2_drhog = -0.5_wp*safe_sqrt((rhol-rhog)/rhol)*( 1.0_wp/(rhol-rhog) )

    dUdda    = w_u*(dUd1_da*safe_sqrt((rhol-rhog)/rhol)*safe_sqrt(c_safe))
    dUddrhol = w_u*((0.2_wp*safe_sqrt(G*DH)*(safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp))*dUd2_drhol*safe_sqrt(c_safe))
    dUddrhog = w_u*((0.2_wp*safe_sqrt(G*DH)*(safe_sqrt(1.0_wp + exp(-18.0_wp*alp)) - 1.0_wp))*dUd2_drhog*safe_sqrt(c_safe))
    dUddul   = 0.0_wp
    dUddug   = 0.0_wp

  end subroutine DCdUd_swananda


  subroutine drift_flux_swananda_ndim(ul, ug, rhol, rhog, ALPHA, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, value)
    !! Python: drift_flux_swananda_ndim(...)
    real(wp), intent(in)  :: ul, ug, rhol, rhog, ALPHA, theta
    real(wp), intent(in)  :: DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol
    real(wp), intent(out) :: value

    real(wp) :: Cd, Ud
    real(wp) :: njl, njg

    call CdUd_swananda(ALPHA, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, Cd, Ud)

    njl = (1.0_wp - ALPHA)*ul
    njg = ALPHA*ug

    value = -njg + ALPHA*( Cd*(njl + njg) + Ud/w_u )
  end subroutine drift_flux_swananda_ndim


  subroutine Ddrift_flux_swananda_ndim(ul, ug, rhol, rhog, ALPHA, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, &
                                      dfda, dfdrhog, dfdrhol, dfdug, dfdul)
    !! Python: Ddrift_flux_swananda_ndim(...)
    real(wp), intent(in)  :: ul, ug, rhol, rhog, ALPHA, theta
    real(wp), intent(in)  :: DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol
    real(wp), intent(out) :: dfda, dfdrhog, dfdrhol, dfdug, dfdul

    real(wp) :: Cd, Ud
    real(wp) :: dC0da, dC0drhol, dC0drhog, dC0dul, dC0dug
    real(wp) :: dUdda, dUddrhol, dUddrhog, dUddul, dUddug
    real(wp) :: njl, njg, jt

    njl = (1.0_wp - ALPHA)*ul
    njg = ALPHA*ug
    jt  = njl + njg

    call CdUd_swananda(ALPHA, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, Cd, Ud)
    call DCdUd_swananda(ALPHA, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, &
                        dC0da, dC0drhol, dC0drhog, dC0dul, dC0dug, dUdda, dUddrhol, dUddrhog, dUddul, dUddug)

    ! f = -njg + alpha*(Cd*jt + Ud/w_u)
    dfda    = -ug + (Cd*jt + Ud/w_u) + ALPHA*( dC0da*jt + Cd*(ug - ul) + dUdda/w_u )
    dfdrhog = ALPHA*( dC0drhog*jt + dUddrhog/w_u )
    dfdrhol = ALPHA*( dC0drhol*jt + dUddrhol/w_u )
    dfdug   = -ALPHA + ALPHA*( dC0dug*jt + Cd*ALPHA + dUddug/w_u )
    dfdul   = ALPHA*( dC0dul*jt + Cd*(1.0_wp - ALPHA) + dUddul/w_u )
  end subroutine Ddrift_flux_swananda_ndim


  subroutine RelEquilLocalPipe_comp(ul, ug, rhol, rhog, ALPHA, BETA, DH, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, value)
    real(wp), intent(in)  :: ul, ug, rhol, rhog, ALPHA, BETA
    real(wp), intent(in)  :: DH, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol
    real(wp), intent(out) :: value

    real(wp), parameter :: ALPHA_EPS=1.0e-12_wp, MU_EPS=1.0e-30_wp, DU_EPS=1.0e-12_wp
    real(wp) :: alp, gamma
    real(wp) :: RESD, Rew, f
    real(wp) :: DL, DG, vL, vG
    real(wp) :: t1, t2, du
    real(wp) :: muL_s, muG_s

    alp = max(ALPHA_EPS, min(1.0_wp-ALPHA_EPS, ALPHA))

    du = ug - ul
    if (abs(du) < DU_EPS) then
      ! no relative slip -> shear term vanishes, keep only hydrostatic part
      t1 = alp*(1.0_wp-alp)*(rhol-rhog)*sin(BETA)*DH
      t1 = t1 / (4.0_wp*w_rho*(w_u**2))
      value = t1
      return
    end if

    muL_s = max(MUL, MU_EPS)
    muG_s = max(MUG, MU_EPS)

    RESD = (w_rho/w_u) * DH * abs(du)
    Rew  = (w_rho*DH*w_u) / muL_s

    gamma = alpha2gamma(alp, tol)
    gamma = max(ALPHA_EPS, min(1.0_wp-ALPHA_EPS, gamma))

    ! IMPORTANT: ffan in your python is ffan(epsD, Re)
    f = ffan(EPS/DH, max(Rew*RESD, 1.0e-6_wp))

    DL = (muL_s*w_u)/(w_rho*DH)
    DG = (muG_s*w_u)/(w_rho*DH)

    vL = DL*(0.045_wp + 0.9_wp*alp)*((0.37_wp + 0.6_wp*gamma)/max(gamma**2, 1.0e-24_wp))*f
    vG = DG*(0.045_wp + 0.9_wp*alp)*((0.37_wp + 0.6_wp*gamma)/max((1.0_wp-gamma)**2, 1.0e-24_wp))*f

    t1 = alp*(1.0_wp-alp)*(rhol-rhog)*sin(BETA)*DH
    t1 = t1 / (4.0_wp*w_rho*(w_u**2))

    t2 = (vL*rhol*(1.0_wp-alp) + vG*rhog*alp)*RESD*du
    t2 = t2 / (DH*DH)

    value = t1 + t2
  end subroutine RelEquilLocalPipe_comp


  subroutine DRelEquiLocal_pipe_comp(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, &
                                    dfda, dfdrhog, dfdrhol, dfdug, dfdul)
    !! Python: DRelEquiLocal_pipe_comp(alpha, rhol, rhog, ul, ug, theta, ...)
    real(wp), intent(in)  :: alpha, rhol, rhog, ul, ug
    real(wp), intent(in)  :: theta, DH, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol
    real(wp), intent(out) :: dfda, dfdrhog, dfdrhol, dfdug, dfdul

    real(wp) :: gamma_p, gamma_m, dgamma_dalpha
    real(wp) :: delta
    real(wp) :: PID, RESD, Rew, Rem
    real(wp) :: DL, DG
    real(wp) :: f, dfdRem, dfdEDIA
    real(wp) :: vp, vm

    real(wp) :: B1, B2
    real(wp) :: C1L, C1G
    real(wp) :: C2L, C2G
    real(wp) :: C3L, C3G
    real(wp) :: C4L, C4G

    real(wp) :: dB1da
    real(wp) :: dC1Lda, dC1Gda
    real(wp) :: dC2Lda, dC2Gda
    real(wp) :: dC3Lda, dC3Gda
    real(wp) :: dC4Lda, dC4Gda

    real(wp) :: dvLda, dvGda
    real(wp) :: dvLdug, dvLdul
    real(wp) :: dvGdug, dvGdul
    real(wp) :: dvLdrhol, dvGdrhog

    real(wp) :: dt1da, dt1drhol, dt1drhog
    real(wp) :: dt2da, dt2drhol, dt2drhog
    real(wp) :: dt2dul, dt2dug

    delta = 1.0e-5_wp
    gamma_p = alpha2gamma(clamp(alpha + delta, 1.0e-12_wp, 1.0_wp-1.0e-12_wp), tol)
    gamma_m = alpha2gamma(clamp(alpha - delta, 1.0e-12_wp, 1.0_wp-1.0e-12_wp), tol)
    dgamma_dalpha = (gamma_p - gamma_m)/(2.0_wp*delta)

    PID = (DH*G)/(4.0_wp*w_u*w_u)

    RESD = (w_rho/w_u) * DH * abs(ug - ul)
    Rew  = (w_rho*DH*w_u) / MUL
    Rem  = Rew*RESD

    DL = (MUL*w_u)/(w_rho*DH)
    DG = (MUG*w_u)/(w_rho*DH)

    f = ffan(EPS/DH, Rem)
    call Dffan(Rem, EPS/DH, dfdRem, dfdEDIA)

    ! dv/df
    vp = (0.37_wp + 0.6_wp*gamma_p)
    vm = (0.37_wp + 0.6_wp*gamma_m)

    ! Common factor B1 = (0.045 + 0.9 alpha)
    B1 = 0.045_wp + 0.9_wp*alpha
    dB1da = 0.9_wp

    C1L = (0.37_wp + 0.6_wp*alpha2gamma(alpha, tol))
    C1G = C1L

    ! C2 terms (gamma powers)
    C2L = (alpha2gamma(alpha, tol))**(-2)
    C2G = (1.0_wp - alpha2gamma(alpha, tol))**(-2)

    ! dC1/dalpha = 0.6 dgamma/dalpha
    dC1Lda = 0.6_wp*dgamma_dalpha
    dC1Gda = dC1Lda

    ! dC2L/dalpha = -2 gamma^-3 dgamma
    dC2Lda = -2.0_wp*(alpha2gamma(alpha, tol))**(-3) * dgamma_dalpha
    ! dC2G/dalpha = -2 (1-gamma)^-3 * d(1-gamma) = +2 (1-gamma)^-3 dgamma
    dC2Gda =  2.0_wp*(1.0_wp - alpha2gamma(alpha, tol))**(-3) * dgamma_dalpha

    ! C3 = f (depend via Rem)
    C3L = f
    C3G = f

    ! dC3/du = dfdRem * dRem/du
    B2 = (w_rho/w_u)*DH
    C4L = (w_rho/w_u)*DH * sign0(ug - ul)
    C4G = C4L

    dC3Lda = 0.0_wp
    dC3Gda = 0.0_wp

    dC4Lda = 0.0_wp
    dC4Gda = 0.0_wp

    ! dRem/dug = Rew * B2 * sign
    dvLdug = DL * (B1*C1L*C2L) * (dfdRem * Rew * B2 * sign0(ug - ul))
    dvGdug = DG * (B1*C1G*C2G) * (dfdRem * Rew * B2 * sign0(ug - ul))

    dvLdul = -DL * (B1*C1L*C2L) * (dfdRem * Rew * B2 * sign0(ug - ul))
    dvGdul = -DG * (B1*C1G*C2G) * (dfdRem * Rew * B2 * sign0(ug - ul))

    dvLdrhol = 0.0_wp
    dvGdrhog = 0.0_wp

    ! dv/da via B1, C1, C2 (f independent of alpha except via Rem; Rem depends on |ug-ul| only)
    dvLda = DL*( dB1da*C1L*C2L*C3L + B1*dC1Lda*C2L*C3L + B1*C1L*dC2Lda*C3L )
    dvGda = DG*( dB1da*C1G*C2G*C3G + B1*dC1Gda*C2G*C3G + B1*C1G*dC2Gda*C3G )

    ! t1 = PID*alpha*(1-alpha)*(rhol-rhog)*sin(theta)
    dt1da    = PID*(1.0_wp-2.0_wp*alpha)*(rhol-rhog)*sin(theta)
    dt1drhol = PID*alpha*(1.0_wp-alpha)*sin(theta)
    dt1drhog = -PID*alpha*(1.0_wp-alpha)*sin(theta)

    ! t2 = (vL*rhol*(1-alpha) + vG*rhog*alpha) * RESD*(ug-ul) / DH^2
    dt2dul = 0.0_wp
    dt2dug = 0.0_wp

    ! RESD*(ug-ul) = (w_rho/w_u)*DH*|ug-ul|*(ug-ul)
    ! Derivative: d/dx [|x|*x] = 2|x|
    ! -> d/dug = (w_rho/w_u)*DH * 2*abs(ug-ul)
    dt2dug = ( (dvLdug*rhol*(1.0_wp-alpha) + dvGdug*rhog*alpha) ) * (RESD*(ug-ul))/(DH*DH) + &
             ( ( (DL*(B1*C1L*C2L*C3L))*rhol*(1.0_wp-alpha) + (DG*(B1*C1G*C2G*C3G))*rhog*alpha ) ) * ((w_rho/w_u)*DH*2.0_wp*abs(ug-ul))/(DH*DH)

    dt2dul = ( (dvLdul*rhol*(1.0_wp-alpha) + dvGdul*rhog*alpha) ) * (RESD*(ug-ul))/(DH*DH) - &
             ( ( (DL*(B1*C1L*C2L*C3L))*rhol*(1.0_wp-alpha) + (DG*(B1*C1G*C2G*C3G))*rhog*alpha ) ) * ((w_rho/w_u)*DH*2.0_wp*abs(ug-ul))/(DH*DH)

    dt2drhol = (dvLdrhol*rhol*(1.0_wp-alpha) + (DL*(B1*C1L*C2L*C3L))*(1.0_wp-alpha)) * (RESD*(ug-ul))/(DH*DH)
    dt2drhog = (dvGdrhog*rhog*alpha + (DG*(B1*C1G*C2G*C3G))*alpha) * (RESD*(ug-ul))/(DH*DH)

    dt2da = (dvLda*rhol*(1.0_wp-alpha) - (DL*(B1*C1L*C2L*C3L))*rhol + &
             dvGda*rhog*alpha + (DG*(B1*C1G*C2G*C3G))*rhog) * (RESD*(ug-ul))/(DH*DH)

    dfda    = dt1da + dt2da
    dfdrhol = dt1drhol + dt2drhol
    dfdrhog = dt1drhog + dt2drhog
    dfdul   = dt2dul
    dfdug   = dt2dug

  end subroutine DRelEquiLocal_pipe_comp


  subroutine lei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, value)
    !! Python: lei_fechamento_or_ndim_simp(...)
    real(wp), intent(in)  :: alpha, rhol, rhog, ul, ug
    real(wp), intent(in)  :: theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol
    real(wp), intent(out) :: value

    if (theta < 0.0_wp) then
      call RelEquilLocalPipe_comp(ul, ug, rhol, rhog, alpha, -theta, DH, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, value)
    else
      call drift_flux_swananda_ndim(ul, ug, rhol, rhog, alpha, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, value)
    end if
  end subroutine lei_fechamento_or_ndim_simp


  subroutine Dlei_fechamento_or_ndim_simp(alpha, rhol, rhog, ul, ug, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, &
                                         dfda, dfdrhog, dfdrhol, dfdug, dfdul)
    !! Dlei_fechamento_or_ndim_simp(...)
    real(wp), intent(in)  :: alpha, rhol, rhog, ul, ug
    real(wp), intent(in)  :: theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol
    real(wp), intent(out) :: dfda, dfdrhog, dfdrhol, dfdug, dfdul

    if (theta < 0.0_wp) then
      call DRelEquiLocal_pipe_comp(alpha, rhol, rhog, ul, ug, -theta, DH, AREA, EPS, G, MUL, MUG, w_u, w_rho, tol, dfda, dfdrhog, dfdrhol, dfdug, dfdul)
    else
      call Ddrift_flux_swananda_ndim(ul, ug, rhol, rhog, alpha, theta, DH, AREA, EPS, G, MUL, MUG, sigma, w_u, w_rho, tol, dfda, dfdrhog, dfdrhol, dfdug, dfdul)
    end if
  end subroutine Dlei_fechamento_or_ndim_simp

end module module_closure_law
