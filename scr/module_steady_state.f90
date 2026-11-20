module module_steady_state
  !! Steady mass fluxes Gl,Gg and total momentum flux
  !! H(P,theta)=Gl*ul+Gg*ug+PIP*P obey dH/ds=S3.
  !! Pressure integration includes explicit theta(s) dependence of the closure.
  !! Inputs P,rho,u and global s/Lr are nondimensional; geometry/properties are SI.
  !! At the pipeline/riser junction pressure is matched, with separate closure
  !! traces, as specified by the reference's coupling assumption.
  use mod_kinds, only: wp,i4
  use mod_void_helpers, only: ffan
  use module_closure_law, only: lei_fechamento_or_ndim_simp,CdUd_drift, &
    STRATIFIED,BENDIKSEN,BHAGWAT,MONTGOMERY,CLOSURE_UNSUPPORTED,is_chexal,chexal_branch_id
  use module_chexal_1992, only: chexal_solve_alpha
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  integer,parameter,public :: STEADY_OK=0,STEADY_INVALID=20,STEADY_NO_ROOT=21
  integer,parameter,public :: STEADY_MULTIPLE_ROOTS=22,STEADY_SINGULAR=23
  integer,parameter,public :: STEADY_INTEGRATION=24,STEADY_PHASE_REQUIRED=25
  integer,parameter,public :: PHASE_AUTO=0,PHASE_LIQUID_ONLY=1,PHASE_GAS_ONLY=2
  public :: malha_progressiva,fun_dpds,steady_state_ndim_simp,EstadoEstacionario_ndim_simp
  public :: solve_alpha_from_superficials,state_at ! same physical state for boundary stationary continuation
contains
  pure real(wp) function nan_value() result(v)
    v=ieee_value(0.0_wp,ieee_quiet_nan)
  end function

  subroutine malha_progressiva(L_total,s_riser,N_before,N_after,ratio_before,ratio_after, &
      malha_centros,malha_nodal,status)
    real(wp),intent(in)::L_total,s_riser,ratio_before,ratio_after
    integer(i4),intent(in)::N_before,N_after
    real(wp),allocatable,intent(out)::malha_centros(:),malha_nodal(:)
    integer,intent(out),optional::status
    real(wp),allocatable::d(:)
    real(wp)::exponent,shift
    integer::i,nb,na,nt
    if(present(status))status=STEADY_INVALID
    if(.not.all(ieee_is_finite([L_total,s_riser,ratio_before,ratio_after])))return
    if(L_total<=0.0_wp.or.s_riser<0.0_wp.or.s_riser>L_total)return
    if(ratio_before<=0.0_wp.or.ratio_after<=0.0_wp)return
    nb=N_before;na=N_after;nt=nb+na
    if(nb<0.or.na<0.or.nt<1)return
    if((s_riser>0.0_wp.and.nb<1).or.(s_riser==0.0_wp.and.nb/=0))return
    if((L_total>s_riser.and.na<1).or.(L_total==s_riser.and.na/=0))return
    allocate(d(nt))
    if(nb>0)then
      shift=max(0.0_wp,real(nb-1,wp)*log(ratio_before))
      do i=1,nb
        exponent=real(nb-i,wp)*log(ratio_before)-shift
        d(i)=exp(exponent)
      enddo
      d(1:nb)=s_riser*d(1:nb)/sum(d(1:nb))
    endif
    if(na>0)then
      shift=max(0.0_wp,real(na-1,wp)*log(ratio_after))
      do i=1,na
        exponent=real(i-1,wp)*log(ratio_after)-shift
        d(nb+i)=exp(exponent)
      enddo
      d(nb+1:nt)=(L_total-s_riser)*d(nb+1:nt)/sum(d(nb+1:nt))
    endif
    if(any(d<=0.0_wp).or..not.all(ieee_is_finite(d)))return
    allocate(malha_nodal(nt+1),malha_centros(nt))
    malha_nodal(1)=0.0_wp
    do i=1,nt
      malha_nodal(i+1)=malha_nodal(i)+d(i)
    enddo
    if(nb>0)malha_nodal(nb+1)=s_riser
    malha_nodal(nt+1)=L_total
    if(any(malha_nodal(2:)<=malha_nodal(:nt)))then
      deallocate(malha_nodal,malha_centros)
      return
    endif
    malha_centros=0.5_wp*(malha_nodal(:nt)+malha_nodal(2:))
    if(present(status))status=STEADY_OK
  end subroutine

  real(wp) function solve_alpha_from_superficials(njl,njg,nrhol,nrhog,theta,DH,AREA,EPS,G, &
      mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol,status,closure_id,nroots,alpha_hint,phase_hint) result(alpha)
    !! Explicit pure-phase inventory limits are exact; zero flux alone does not
    !! imply phase absence. For stratified flow and alpha-dependent Bhagwat, scan a
    !! logit grid and bisect every detected sign-changing bracket. Multiple roots
    !! are reported unless the caller supplies alpha_hint for explicit branch
    !! continuation (closest logit root). There is no best-point fallback or
    !! clipping of a nonphysical root. nroots requests the full diagnostic scan.
    real(wp),intent(in)::njl,njg,nrhol,nrhog,theta,DH,AREA,EPS,G
    real(wp),intent(in)::mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol
    integer,intent(out),optional::status,nroots
    integer,intent(in),optional::closure_id,phase_hint
    real(wp),intent(in),optional::alpha_hint
    integer,parameter::nscan=160
    real(wp)::a,b,fa,fb,mid,cd,ud,den,xtol,x,rootlast,re_l_min
    real(wp)::roots(nscan+2),hint_logit,distance,best_distance
    integer::id,i,ns,ierr,selected
    logical::previous_valid,tied
    alpha=nan_value()
    if(present(status))status=STEADY_INVALID
    if(present(nroots))nroots=0
    if(present(alpha_hint))then
      if(.not.ieee_is_finite(alpha_hint))return
      if(alpha_hint<0.0_wp.or.alpha_hint>1.0_wp)return
    endif
    if(.not.all(ieee_is_finite([njl,njg,nrhol,nrhog,theta,DH,AREA,EPS,G, &
      mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol])))return
    if(nrhol<=0.0_wp.or.nrhog<=0.0_wp.or.DH<=0.0_wp.or.AREA<=0.0_wp)return
    if(EPS<0.0_wp.or.G<0.0_wp.or.mu_l_dyn<=0.0_wp.or.mu_g_dyn<=0.0_wp)return
    if(sigma<0.0_wp.or.w_u<=0.0_wp.or.w_rho<=0.0_wp.or.tol<=0.0_wp)return
    id=BENDIKSEN
    if(theta<0.0_wp)id=STRATIFIED
    if(present(closure_id))id=closure_id
    if(id/=STRATIFIED.and.id/=BENDIKSEN.and.id/=BHAGWAT.and.id/=MONTGOMERY.and..not.is_chexal(id))then
      if(present(status))status=CLOSURE_UNSUPPORTED
      return
    endif
    if(is_chexal(id))then
      if(abs(AREA-acos(-1.0_wp)*DH**2/4.0_wp)>1e-6_wp*AREA)return
    endif
    ! A zero PHASE FLUX does not establish that the phase inventory is absent.
    ! Pure phases must be declared by a caller that knows the inventory.
    if(present(phase_hint))then
      select case(phase_hint)
      case(PHASE_LIQUID_ONLY)
        if(njg/=0.0_wp)return
        alpha=0.0_wp
        if(present(status))status=0
        if(present(nroots))nroots=1
        return
      case(PHASE_GAS_ONLY)
        if(njl/=0.0_wp)return
        alpha=1.0_wp
        if(present(status))status=0
        if(present(nroots))nroots=1
        return
      case(PHASE_AUTO)
        continue
      case default
        return
      end select
    endif
    if(is_chexal(id))then
      call chexal_solve_alpha(njl*w_u,njg*w_u,nrhol*w_rho,nrhog*w_rho,theta,DH,G, &
        mu_l_dyn,mu_g_dyn,sigma,alpha,ierr,branch=chexal_branch_id(id),nroots=ns,alpha_hint=alpha_hint)
      ! An endpoint root of the algebraic equation is not evidence that its
      ! phase inventory is absent. Exact inventory hints were handled above.
      if(ierr==0.and.(alpha==0.0_wp.or.alpha==1.0_wp))then
        ierr=STEADY_PHASE_REQUIRED
        alpha=nan_value()
      endif
      if(present(status))status=ierr
      if(present(nroots))nroots=ns
      return
    endif
    if(id==BHAGWAT)then
      ! Qualify the global two-phase input domain through the closure itself.
      ! This is not an explicit alpha solution: Bhagwat's C0 and Ud depend on
      ! alpha, and every residual below recomputes them at the candidate alpha.
      ! Pure phases above are explicit inventory limits, outside this two-phase
      ! correlation. In particular, zero flow alone never selects alpha=0 or 1.
      call CdUd_drift(0.5_wp,nrhol,nrhog,2.0_wp*njl,2.0_wp*njg,theta,DH,AREA,EPS,G, &
        mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol,cd,ud,id,ierr)
      if(ierr/=0)then
        if(present(status))status=ierr
        return
      endif
    endif
    if(njl==0.0_wp.and.njg==0.0_wp)then
      if(G==0.0_wp.or.(id==STRATIFIED.and.sin(theta)==0.0_wp))then
        ! Continuum of stationary roots: the known inventory/hint is required.
        if(present(nroots))nroots=-1
        if(present(status))status=STEADY_MULTIPLE_ROOTS
        if(present(alpha_hint))then
          alpha=alpha_hint
          if(present(status))status=0
        endif
        return
      endif
    endif
    if(id==BENDIKSEN.or.id==MONTGOMERY)then
      ! At fixed superficial velocities coefficients do not depend on alpha.
      call CdUd_drift(0.5_wp,nrhol,nrhog,2.0_wp*njl,2.0_wp*njg,theta,DH,AREA,EPS,G, &
        mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol,cd,ud,id,ierr)
      if(ierr/=0)then
        if(present(status))status=ierr
        return
      endif
      den=cd*(njl+njg)+ud/w_u
      if(den==0.0_wp)then
        if(present(status))status=STEADY_NO_ROOT
        return
      endif
      a=njg/den
      if(.not.ieee_is_finite(a))return
      if(a<0.0_wp.or.a>1.0_wp.or.(a==0.0_wp.and.njg/=0.0_wp).or.(a==1.0_wp.and.njl/=0.0_wp))then
        if(present(status))status=STEADY_NO_ROOT
        return
      endif
      alpha=a
      if(present(status))status=0
      if(present(nroots))nroots=1
      return
    endif
    xtol=max(8.0_wp*epsilon(1.0_wp),min(1.0e-12_wp,0.01_wp*tol))
    ! Positive co-current downhill flow, with no interfacial-velocity transition,
    ! uses the regular single-root branch. Endpoint bracketing + safeguarded Brent
    ! is the fast production path. Requesting nroots always enables the full scan
    ! for branch audits; countercurrent/low-Re states always use that scan too.
    re_l_min=nrhol*w_rho*njl*w_u*DH/mu_l_dyn
    if(id==STRATIFIED.and.njl>0.0_wp.and.njg>0.0_wp.and.theta<=0.0_wp.and. &
      re_l_min>=2200.0_wp.and..not.present(nroots))then
      a=1.0_wp/(1.0_wp+exp(28.0_wp));b=1.0_wp-a
      call residual(a,fa,ierr)
      if(ierr==0)then
        call residual(b,fb,ierr)
        if(ierr==0)then
          if((fa>0.0_wp.and.fb<0.0_wp).or.(fa<0.0_wp.and.fb>0.0_wp))then
            call root_bracket(a,b,fa,fb,mid,ierr)
            if(ierr==0)then
              alpha=mid
              if(present(status))status=0
              return
            endif
          endif
        endif
      endif
    endif
    ns=0;rootlast=-1.0_wp;previous_valid=.false.
    b=0.0_wp;fb=0.0_wp
    ! |logit(alpha)|<=28 covers alpha down to6.9e-13. No finite grid proves
    ! absence of a tangential root; such degeneracy is not silently accepted.
    do i=0,nscan
      x=-28.0_wp+56.0_wp*real(i,wp)/real(nscan,wp)
      a=1.0_wp/(1.0_wp+exp(-x))
      call residual(a,fa,ierr)
      if(ierr/=0.or..not.ieee_is_finite(fa))then
        previous_valid=.false.
        cycle
      endif
      if(fa==0.0_wp)then
        call record(a)
      else if(previous_valid)then
        if((fa>0.0_wp.and.fb<0.0_wp).or.(fa<0.0_wp.and.fb>0.0_wp))then
          call root_bracket(b,a,fb,fa,mid,ierr)
          if(ierr==0)call record(mid)
        endif
      endif
      b=a;fb=fa;previous_valid=.true.
    enddo
    if(ns==0.and.id==BHAGWAT.and.theta>=0.0_wp.and.njl>=0.0_wp.and.njg>0.0_wp)then
      ! Only extend a failed scan. Existing successes, multiple-root selection,
      ! pure-phase declarations and the qualified countercurrent exclusions are
      ! unchanged. The regular co-current riser tail is strictly monotone; the
      ! helper checks a positive C0 bound before using that property.
      a=1.0_wp/(1.0_wp+exp(28.0_wp))
      call lower_riser_tail(a,mid,ierr)
      if(ierr==0)call record(mid)
    endif
    if(present(nroots))nroots=ns
    if(ns==1)then
      alpha=rootlast
      if(present(status))status=0
    else if(ns>1)then
      if(present(status))status=STEADY_MULTIPLE_ROOTS
      if(present(alpha_hint))then
        ! Continuation is an explicit caller choice, not a root-clipping rule.
        ! Endpoint hints select the corresponding limit in logit distance.
        selected=0;tied=.false.
        if(alpha_hint==0.0_wp)then
          selected=1
        else if(alpha_hint==1.0_wp)then
          selected=ns
        else
          hint_logit=log(alpha_hint/(1.0_wp-alpha_hint))
          best_distance=huge(1.0_wp)
          do i=1,ns
            distance=abs(log(roots(i)/(1.0_wp-roots(i)))-hint_logit)
            if(distance<best_distance-1e-12_wp)then
              selected=i;best_distance=distance;tied=.false.
            else if(abs(distance-best_distance)<=1e-12_wp)then
              tied=.true.
            endif
          enddo
        endif
        if(selected>0.and..not.tied)then
          alpha=roots(selected)
          if(present(status))status=0
        endif
      endif
    else
      if(present(status))status=STEADY_NO_ROOT
    endif
  contains
    subroutine lower_riser_tail(upper,rr,code)
      !! Geometric bracketing below the old alpha minimum, then relative
      !! bisection. A representable bracket is required; no alpha floor is
      !! substituted for a root. This path never changes a previously found
      !! root and does not extend the closure to negative phase velocities.
      real(wp),intent(in)::upper
      real(wp),intent(out)::rr
      integer,intent(out)::code
      real(wp)::left,right,fleft,fright,midpoint,fmid,smallest,cd_tail,ud_tail,cd_margin
      integer::step,iteration,local_status
      rr=nan_value();code=STEADY_NO_ROOT
      call CdUd_drift(upper,nrhol,nrhog,njl/(1.0_wp-upper),njg/upper,theta,DH,AREA,EPS,G, &
        mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol,cd_tail,ud_tail,id,local_status)
      if(local_status/=0)return
      ! At fixed superficial fluxes and 0<=theta<=pi/2, 0<=dC0/da<=ln(2)/5.
      ! This sufficient bound ensures C0>0 throughout 0<a<=upper; with Ud>=0,
      ! d[a(C0*j+Ud)]/da>0 for a<2/3. The lower tail then has at most one root.
      cd_margin=0.2_wp*log(2.0_wp)*upper+32.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(cd_tail))
      if(cd_tail<=cd_margin.or.ud_tail<0.0_wp)return
      right=upper
      call residual(right,fright,local_status)
      if(local_status/=0.or..not.ieee_is_finite(fright))return
      if(fright<=0.0_wp)return
      smallest=nearest(0.0_wp,1.0_wp)
      do step=1,maxexponent(1.0_wp)-minexponent(1.0_wp)+digits(1.0_wp)+8
        if(right<=smallest)return
        left=0.5_wp*right
        if(left<=0.0_wp.or.left>=right)return
        ! Avoid an overflow before calling a closure that takes phase speeds.
        if(left<njg/huge(1.0_wp))return
        call residual(left,fleft,local_status)
        if(local_status/=0.or..not.ieee_is_finite(fleft))return
        if(fleft==0.0_wp)then
          rr=left;code=0;return
        endif
        if(fleft<0.0_wp)then
          do iteration=1,digits(1.0_wp)+8
            midpoint=left+0.5_wp*(right-left)
            if(midpoint==left.or.midpoint==right.or.right-left<=xtol*right)then
              rr=left
              if(abs(fright)<abs(fleft))rr=right
              code=0;return
            endif
            call residual(midpoint,fmid,local_status)
            if(local_status/=0.or..not.ieee_is_finite(fmid))return
            if(fmid==0.0_wp)then
              rr=midpoint;code=0;return
            else if(fmid<0.0_wp)then
              left=midpoint;fleft=fmid
            else
              right=midpoint;fright=fmid
            endif
          enddo
          return
        endif
        right=left;fright=fleft
      enddo
    end subroutine
    subroutine residual(aa,value,code)
      real(wp),intent(in)::aa
      real(wp),intent(out)::value
      integer,intent(out)::code
      call lei_fechamento_or_ndim_simp(aa,nrhol,nrhog,njl/(1.0_wp-aa),njg/aa,theta,DH,AREA,EPS,G, &
        mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol,value,id,code)
    end subroutine
    subroutine root_bracket(left,right,fleft,fright,rr,code)
      !! Brent-Dekker, called directly (no procedure pointer or executable stack).
      real(wp),intent(in)::left,right,fleft,fright
      real(wp),intent(out)::rr
      integer,intent(out)::code
      real(wp)::aa,bb,cc,ff_a,ff_b,ff_c,dd,ee,midpoint,tol_b,ss,pp,qq,ratio
      real(wp)::old_b,old_fb
      integer::iteration
      aa=left;bb=right;cc=bb
      ff_a=fleft;ff_b=fright;ff_c=ff_b
      dd=bb-aa;ee=dd;rr=nan_value();code=STEADY_NO_ROOT
      do iteration=1,100
        if((ff_b>0.0_wp.and.ff_c>0.0_wp).or.(ff_b<0.0_wp.and.ff_c<0.0_wp))then
          cc=aa;ff_c=ff_a;dd=bb-aa;ee=dd
        endif
        if(abs(ff_c)<abs(ff_b))then
          old_b=bb;old_fb=ff_b
          aa=old_b;ff_a=old_fb;bb=cc;ff_b=ff_c;cc=old_b;ff_c=old_fb
        endif
        tol_b=2.0_wp*epsilon(1.0_wp)*abs(bb)+0.5_wp*xtol
        midpoint=0.5_wp*(cc-bb)
        if(abs(midpoint)<=tol_b.or.ff_b==0.0_wp)then
          rr=bb;code=0;return
        endif
        if(abs(ee)>=tol_b.and.abs(ff_a)>abs(ff_b))then
          ss=ff_b/ff_a
          if(aa==cc)then
            pp=2.0_wp*midpoint*ss;qq=1.0_wp-ss
          else
            qq=ff_a/ff_c;ratio=ff_b/ff_c
            pp=ss*(2.0_wp*midpoint*qq*(qq-ratio)-(bb-aa)*(ratio-1.0_wp))
            qq=(qq-1.0_wp)*(ratio-1.0_wp)*(ss-1.0_wp)
          endif
          if(pp>0.0_wp)qq=-qq
          pp=abs(pp)
          if(2.0_wp*pp<min(3.0_wp*midpoint*qq-abs(tol_b*qq),abs(ee*qq)))then
            ee=dd;dd=pp/qq
          else
            dd=midpoint;ee=dd
          endif
        else
          dd=midpoint;ee=dd
        endif
        aa=bb;ff_a=ff_b
        if(abs(dd)>tol_b)then
          bb=bb+dd
        else
          bb=bb+sign(tol_b,midpoint)
        endif
        call residual(bb,ff_b,code)
        if(code/=0)return
        if(.not.ieee_is_finite(ff_b))then
          code=STEADY_NO_ROOT;return
        endif
      enddo
      code=STEADY_NO_ROOT
    end subroutine
    subroutine record(rr)
      real(wp),intent(in)::rr
      if(ns>0)then
        if(abs(rr-rootlast)<=8.0_wp*xtol)return
      endif
      ns=ns+1;rootlast=rr;roots(ns)=rr
    end subroutine
  end function

  subroutine state_at(p,theta,gl,gg,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
      w_u,w_rho,tol,id,rl,rg,alpha,ul,ug,jl,jg,kin,status,phase_hint)
    real(wp),intent(in)::p,theta,gl,gg,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0
    real(wp),intent(in)::mu_l_dyn,mu_g_dyn,sigma,w_u,w_rho,tol
    integer,intent(in)::id
    integer,intent(in),optional::phase_hint
    real(wp),intent(out)::rl,rg,alpha,ul,ug,jl,jg,kin
    integer,intent(out)::status
    rl=nan_value();rg=rl;alpha=rl;ul=rl;ug=rl;jl=rl;jg=rl;kin=rl
    status=STEADY_INVALID
    if(.not.all(ieee_is_finite([p,Cl,Cg,rho_l0,P_l0,gl,gg])))return
    if(p<=0.0_wp.or.Cl<=0.0_wp.or.Cg<=0.0_wp.or.rho_l0<=0.0_wp)return
    rg=p/(Cg*Cg);rl=rho_l0+(p-P_l0)/(Cl*Cl)
    if(rl<=0.0_wp.or.rg<=0.0_wp)return
    jl=gl/rl;jg=gg/rg
    if(present(phase_hint))then
      if(phase_hint/=PHASE_AUTO)then
        alpha=solve_alpha_from_superficials(jl,jg,rl,rg,theta,DH,AREA,EPS,G,mu_l_dyn,mu_g_dyn, &
          sigma,w_u,w_rho,tol,status,id,phase_hint=phase_hint)
        if(status/=0)return
        ul=0;ug=0
        if(alpha==0)ul=jl
        if(alpha==1)ug=jg
        kin=gl*ul+gg*ug
        if(.not.ieee_is_finite(kin))status=STEADY_INVALID
        return
      endif
    endif
    if(gl==0.0_wp.and.gg==0.0_wp)then
      ! Defined initialization scenario: stationary liquid-filled hydrostatics.
      ! This is a declared inventory choice, not inferred for a flowing phase.
      alpha=solve_alpha_from_superficials(jl,jg,rl,rg,theta,DH,AREA,EPS,G,mu_l_dyn,mu_g_dyn, &
        sigma,w_u,w_rho,tol,status,id,phase_hint=PHASE_LIQUID_ONLY)
    else
      alpha=solve_alpha_from_superficials(jl,jg,rl,rg,theta,DH,AREA,EPS,G,mu_l_dyn,mu_g_dyn, &
        sigma,w_u,w_rho,tol,status,id)
    endif
    if(status/=0)return
    if(alpha==0.0_wp)then
      ul=jl;ug=0.0_wp
    else if(alpha==1.0_wp)then
      ul=0.0_wp;ug=jg
    else
      ul=jl/(1.0_wp-alpha);ug=jg/alpha
    endif
    kin=gl*ul+gg*ug
    if(.not.ieee_is_finite(kin))status=STEADY_INVALID
  end subroutine

  real(wp) function fun_dpds(s,p,mul,mug,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0, &
      mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol,status,closure_id_riser,section_id,closure_id_pipe,phase_hint) result(dpds)
    !! Optional section_id forces the pipe/riser trace at the shared junction;
    !! 1=pipe and2=riser. It never changes the closure's model definition.
    real(wp),intent(in)::s,p,mul,mug,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0
    real(wp),intent(in)::mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol
    integer,intent(out),optional::status
    integer,intent(in),optional::closure_id_riser,section_id,closure_id_pipe,phase_hint
    real(wp)::theta,theta_s,sr,pip,pig,pif,rl,rg,alpha,ul,ug,jl,jg,kin
    real(wp)::rl1,rg1,a1,ul1,ug1,jl1,jg1,kp,km,hp,ht,delta,h_theta,source
    real(wp)::rm,j,mu,re,f,fric,scale,geometry_tol,theta_other
    integer::cid,sec,ierr
    dpds=nan_value()
    if(present(status))status=STEADY_INVALID
    if(.not.all(ieee_is_finite([s,p,Lp,Lr,CA,beta,w_p,w_u,w_rho,DH,AREA,EPS,G, &
      Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma,tol])))return
    if(Lr<=0.0_wp.or.Lp<0.0_wp.or.CA<=0.0_wp.or.w_p<=0.0_wp.or.w_u<=0.0_wp.or.w_rho<=0.0_wp)return
    if(DH<=0.0_wp.or.G<0.0_wp.or.p<=0.0_wp.or.tol<=0.0_wp)return
    sec=2
    if(s<Lp/Lr)sec=1
    if(present(section_id))sec=section_id
    if(sec/=1.and.sec/=2)return
    cid=BENDIKSEN
    if(present(closure_id_riser))cid=closure_id_riser
    if(sec==1)then
      theta=-beta;theta_s=0.0_wp;cid=STRATIFIED
      if(present(closure_id_pipe))cid=closure_id_pipe
    else
      ! Compute the coordinate relative to the exact shared mesh face first.
      ! A sub-ulp roundoff excursion must not select the downhill C4 trace of
      ! Bhagwat at a geometrically horizontal riser boundary. A real excursion
      ! outside the riser remains invalid; this is not a physical angle clip.
      sr=Lr*(s-Lp/Lr)
      geometry_tol=64.0_wp*epsilon(1.0_wp)*max(Lp,Lr,abs(s*Lr),CA)
      if(abs(sr)<=geometry_tol)sr=0.0_wp
      if(sr<0.0_wp)return
      theta=atan(sr/CA)
      theta_s=Lr*CA/(CA*CA+sr*sr)
    endif
    pip=w_p/(w_rho*w_u*w_u);pig=G*Lr/(w_u*w_u);pif=Lr/DH
    call state_at(p,theta,mul,mug,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
      w_u,w_rho,tol,cid,rl,rg,alpha,ul,ug,jl,jg,kin,ierr,phase_hint=phase_hint)
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    ! Separate pressure term analytically to avoid cancellation of PIP*P.
    delta=min(0.25_wp*p,epsilon(1.0_wp)**(1.0_wp/3.0_wp)*max(abs(p),1.0_wp))
    call state_at(p+delta,theta,mul,mug,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
      w_u,w_rho,tol,cid,rl1,rg1,a1,ul1,ug1,jl1,jg1,kp,ierr,phase_hint=phase_hint)
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    call state_at(p-delta,theta,mul,mug,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
      w_u,w_rho,tol,cid,rl1,rg1,a1,ul1,ug1,jl1,jg1,km,ierr,phase_hint=phase_hint)
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    hp=pip+(kp-km)/(2.0_wp*delta)
    scale=max(abs(pip),abs((kp-km)/(2.0_wp*delta)),1.0_wp)
    if(.not.ieee_is_finite(hp).or.abs(hp)<1.0e-10_wp*scale)then
      if(present(status))status=STEADY_SINGULAR
      return
    endif
    ht=0.0_wp
    if(theta_s/=0.0_wp.and.mul/=0.0_wp.and.mug/=0.0_wp)then
      h_theta=epsilon(1.0_wp)**(1.0_wp/3.0_wp)
      ! Htheta is the derivative of the RISER trace. Its geometric domain is
      ! 0<=theta<pi/2. At the lower face, a centered stencil would cross the real
      ! Bhagwat C4 jump at theta=0 and inject an O(1/h) term into the ODE. Use
      ! second-order one-sided derivatives at angular boundaries, not smoothing.
      if(theta<h_theta)then
        theta_other=theta+h_theta
      else if(theta+h_theta>0.5_wp*acos(-1.0_wp))then
        theta_other=theta-h_theta
      else
        theta_other=theta+h_theta
      endif
      call state_at(p,theta_other,mul,mug,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
        w_u,w_rho,tol,cid,rl1,rg1,a1,ul1,ug1,jl1,jg1,kp,ierr,phase_hint=phase_hint)
      if(ierr/=0)then
        if(present(status))status=ierr
        return
      endif
      if(theta<h_theta)then
        theta_other=theta+2.0_wp*h_theta
      else if(theta+h_theta>0.5_wp*acos(-1.0_wp))then
        theta_other=theta-2.0_wp*h_theta
      else
        theta_other=theta-h_theta
      endif
      call state_at(p,theta_other,mul,mug,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
        w_u,w_rho,tol,cid,rl1,rg1,a1,ul1,ug1,jl1,jg1,km,ierr,phase_hint=phase_hint)
      if(ierr/=0)then
        if(present(status))status=ierr
        return
      endif
      if(theta<h_theta)then
        ht=(4.0_wp*(kp-kin)-(km-kin))/(2.0_wp*h_theta)
      else if(theta+h_theta>0.5_wp*acos(-1.0_wp))then
        ht=(4.0_wp*(kin-kp)-(kin-km))/(2.0_wp*h_theta)
      else
        ht=(kp-km)/(2.0_wp*h_theta)
      endif
    endif
    rm=(1.0_wp-alpha)*rl+alpha*rg;j=jl+jg
    mu=(1.0_wp-alpha)*mu_l_dyn+alpha*mu_g_dyn
    fric=0.0_wp
    if(j/=0.0_wp)then
      re=w_rho*w_u*DH*rm*abs(j)/mu
      if(re<2000.0_wp)then
        ! Analytic finite signed limit of -2*(Lr/D)*f_F*rm*j*abs(j).
        fric=-32.0_wp*Lr*mu*j/(DH*DH*w_rho*w_u)
      else
        f=ffan(EPS/DH,re,ierr)
        if(ierr/=0)then
          if(present(status))status=ierr
          return
        endif
        fric=-2.0_wp*pif*f*rm*j*abs(j)
      endif
    endif
    source=-pig*rm*sin(theta)+fric
    dpds=(source-ht*theta_s)/hp
    if(.not.ieee_is_finite(dpds))then
      dpds=nan_value()
      if(present(status))status=STEADY_INVALID
    else
      if(present(status))status=STEADY_OK
    endif
  end function

  subroutine steady_state_ndim_simp(N,mul,mug,Ps,Lp,Lr,CA,beta,DH,AREA,EPS,G, &
      Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol, &
      vns,vnp,nrhogv,nrholv,alphav,nugv,nulv,thetav,njlv,njgv, &
      mesh_ratio,n_riser,face_positions,closure_id_riser,status,face_pressures,closure_id_pipe,phase_hint)
    integer(i4),intent(in)::N
    real(wp),intent(in)::mul,mug,Ps,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0
    real(wp),intent(in)::mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol
    real(wp),allocatable,intent(out)::vns(:),vnp(:),nrhogv(:),nrholv(:),alphav(:),nugv(:),nulv(:)
    real(wp),allocatable,intent(out)::thetav(:),njlv(:),njgv(:)
    real(wp),intent(in),optional::mesh_ratio
    integer,intent(in),optional::n_riser,closure_id_riser,closure_id_pipe,phase_hint
    real(wp),allocatable,intent(out),optional::face_positions(:),face_pressures(:)
    integer,intent(out),optional::status
    real(wp),allocatable::faces(:),pf(:)
    real(wp)::ratio,s,p,theta,kin,rtol
    integer::nb,na,nt,i,sec,cid,iriser,ipipe,ierr
    if(present(status))status=STEADY_INVALID
    if(.not.all(ieee_is_finite([Ps,Lp,Lr,CA,tol])))return
    if(Ps<=0.0_wp.or.Lp<0.0_wp.or.Lr<=0.0_wp.or.CA<=0.0_wp.or.tol<=0.0_wp)return
    nb=N
    if(Lp==0.0_wp)nb=0
    na=max(1,N/2)
    if(present(n_riser))na=n_riser
    ratio=1.0_wp
    if(present(mesh_ratio))ratio=mesh_ratio
    iriser=BENDIKSEN
    if(present(closure_id_riser))iriser=closure_id_riser
    ipipe=STRATIFIED
    if(present(closure_id_pipe))ipipe=closure_id_pipe
    call malha_progressiva((Lp+Lr)/Lr,Lp/Lr,nb,na,ratio,ratio,vns,faces,ierr)
    if(ierr/=0)then
      if(present(status))status=ierr
      return
    endif
    nt=size(vns)
    allocate(vnp(nt),nrhogv(nt),nrholv(nt),alphav(nt),nugv(nt),nulv(nt),thetav(nt),njlv(nt),njgv(nt),pf(nt+1))
    vnp=nan_value();nrhogv=vnp;nrholv=vnp;alphav=vnp;nugv=vnp;nulv=vnp;thetav=vnp;njlv=vnp;njgv=vnp;pf=nan_value()
    if(present(face_positions))face_positions=faces
    ! The boundary is at the TOP FACE, never at the last cell center.
    s=faces(nt+1);p=Ps;pf(nt+1)=Ps
    rtol=max(tol,1.0e-11_wp) ! FD pressure derivatives do not support an arbitrarily tiny ODE tolerance.
    do i=nt,1,-1
      sec=2;cid=iriser
      if(i<=nb)then
        sec=1;cid=ipipe
      endif
      call advance(vns(i),sec,s,p,ierr)
      if(ierr/=0)exit
      vnp(i)=p
      if(sec==1)then
        theta=-beta
      else
        theta=atan((vns(i)*Lr-Lp)/CA)
      endif
      thetav(i)=theta
      call state_at(p,theta,mul,mug,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0,mu_l_dyn,mu_g_dyn,sigma, &
        w_u,w_rho,tol,cid,nrholv(i),nrhogv(i),alphav(i),nulv(i),nugv(i),njlv(i),njgv(i),kin,ierr,phase_hint=phase_hint)
      if(ierr/=0)exit
      call advance(faces(i),sec,s,p,ierr)
      if(ierr/=0)exit
      pf(i)=p
    enddo
    if(present(face_pressures))face_pressures=pf
    if(present(status))status=ierr
  contains
    subroutine advance(target,section,ss,pp,code)
      real(wp),intent(in)::target
      integer,intent(in)::section
      real(wp),intent(inout)::ss,pp
      integer,intent(out)::code
      real(wp)::h,remaining,y5,y4,err,sc,fac,k(7)
      integer::it,ic
      code=0
      if(target==ss)return
      if(target>ss)then
        code=STEADY_INVALID;return
      endif
      h=target-ss
      do it=1,100000
        remaining=target-ss
        if(abs(remaining)<=16.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(ss)))then
          ss=target;return
        endif
        if(abs(h)>abs(remaining))h=remaining
        if(abs(h)<64.0_wp*epsilon(1.0_wp)*max(1.0_wp,abs(ss)))then
          code=STEADY_INTEGRATION;return
        endif
        call rhs(ss,pp,section,k(1),ic)
        if(ic/=0)then
          code=ic;return
        endif
        call rhs(ss+h/5.0_wp,pp+h*k(1)/5.0_wp,section,k(2),ic)
        if(ic==0)call rhs(ss+3.0_wp*h/10.0_wp,pp+h*(3.0_wp*k(1)/40.0_wp+9.0_wp*k(2)/40.0_wp), &
          section,k(3),ic)
        if(ic==0)call rhs(ss+4.0_wp*h/5.0_wp,pp+h*(44.0_wp*k(1)/45.0_wp-56.0_wp*k(2)/15.0_wp+ &
          32.0_wp*k(3)/9.0_wp),section,k(4),ic)
        if(ic==0)call rhs(ss+8.0_wp*h/9.0_wp,pp+h*(19372.0_wp*k(1)/6561.0_wp-25360.0_wp*k(2)/2187.0_wp+ &
          64448.0_wp*k(3)/6561.0_wp-212.0_wp*k(4)/729.0_wp),section,k(5),ic)
        if(ic==0)call rhs(ss+h,pp+h*(9017.0_wp*k(1)/3168.0_wp-355.0_wp*k(2)/33.0_wp+ &
          46732.0_wp*k(3)/5247.0_wp+49.0_wp*k(4)/176.0_wp-5103.0_wp*k(5)/18656.0_wp),section,k(6),ic)
        if(ic==0)then
          y5=pp+h*(35.0_wp*k(1)/384.0_wp+500.0_wp*k(3)/1113.0_wp+125.0_wp*k(4)/192.0_wp- &
            2187.0_wp*k(5)/6784.0_wp+11.0_wp*k(6)/84.0_wp)
          call rhs(ss+h,y5,section,k(7),ic)
        endif
        if(ic/=0)then
          h=0.5_wp*h
          cycle
        endif
        y4=pp+h*(5179.0_wp*k(1)/57600.0_wp+7571.0_wp*k(3)/16695.0_wp+393.0_wp*k(4)/640.0_wp- &
          92097.0_wp*k(5)/339200.0_wp+187.0_wp*k(6)/2100.0_wp+k(7)/40.0_wp)
        err=abs(y5-y4);sc=rtol*(1.0_wp+max(abs(pp),abs(y5)))
        if(.not.ieee_is_finite(err).or.y5<=0.0_wp)then
          h=0.5_wp*h
          cycle
        endif
        if(err<=sc)then
          ss=ss+h;pp=y5
        endif
        if(err<=tiny(1.0_wp))then
          fac=2.0_wp
        else
          fac=min(2.0_wp,max(0.2_wp,0.9_wp*(sc/err)**0.2_wp))
        endif
        h=h*fac
      enddo
      code=STEADY_INTEGRATION
    end subroutine
    subroutine rhs(ss,pp,section,value,code)
      real(wp),intent(in)::ss,pp
      integer,intent(in)::section
      real(wp),intent(out)::value
      integer,intent(out)::code
      value=fun_dpds(ss,pp,mul,mug,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0, &
        mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol,code,iriser,section,ipipe,phase_hint)
    end subroutine
  end subroutine

  subroutine EstadoEstacionario_ndim_simp(N,mul,mug,Ps,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0, &
      mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol,vns,vnp,nrhogv,nrholv,alphav,nugv,nulv,thetav,njlv,njgv, &
      mesh_ratio,n_riser,face_positions,closure_id_riser,status,face_pressures,closure_id_pipe,phase_hint)
    integer(i4),intent(in)::N
    real(wp),intent(in)::mul,mug,Ps,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0
    real(wp),intent(in)::mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol
    real(wp),allocatable,intent(out)::vns(:),vnp(:),nrhogv(:),nrholv(:),alphav(:),nugv(:),nulv(:)
    real(wp),allocatable,intent(out)::thetav(:),njlv(:),njgv(:)
    real(wp),intent(in),optional::mesh_ratio
    integer,intent(in),optional::n_riser,closure_id_riser,closure_id_pipe,phase_hint
    real(wp),allocatable,intent(out),optional::face_positions(:),face_pressures(:)
    integer,intent(out),optional::status
    call steady_state_ndim_simp(N,mul,mug,Ps,Lp,Lr,CA,beta,DH,AREA,EPS,G,Cl,Cg,rho_l0,P_l0, &
      mu_l_dyn,mu_g_dyn,sigma,w_p,w_u,w_rho,tol,vns,vnp,nrhogv,nrholv,alphav,nugv,nulv,thetav,njlv,njgv, &
      mesh_ratio,n_riser,face_positions,closure_id_riser,status,face_pressures,closure_id_pipe,phase_hint)
  end subroutine
end module module_steady_state
