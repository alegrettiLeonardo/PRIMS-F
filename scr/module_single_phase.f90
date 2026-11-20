module module_single_phase
  !! Exact absent-inventory reduction: TWO active conservative equations.
  !! The inactive phase stays exactly zero; it is never assigned a wave speed.
  use mod_kinds,only:wp
  use module_config,only:sim_config_t,first_physical_cell,last_physical_cell
  use mod_void_helpers,only:ffan
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  private
  integer,parameter,public::PHASE_TRANSITION_REQUIRED=932,PHASE_GEOMETRY_REQUIRED=933
  integer,parameter,public::PHASE_DATA_REQUIRED=934,PHASE_INVALID=935
  real(wp),parameter::qx(8)=[-.9602898564975363_wp,-.7966664774136267_wp,-.5255324099163290_wp, &
    -.1834346424956498_wp,.1834346424956498_wp,.5255324099163290_wp,.7966664774136267_wp,.9602898564975363_wp]
  real(wp),parameter::qw(8)=[.1012285362903763_wp,.2223810344533745_wp,.3137066458778873_wp, &
    .3626837833783620_wp,.3626837833783620_wp,.3137066458778873_wp,.2223810344533745_wp,.1012285362903763_wp]
  public::global_single_phase,single_phase_assemble,single_phase_face,single_phase_boundary,hydro_factors
  public::configure_pipe_riser_elevation,initialize_single_phase_hydrostatic
contains
  subroutine configure_pipe_riser_elevation(cfg)
    type(sim_config_t),intent(inout)::cfg
    integer::i,j,n
    real(wp)::left,right,s
    n=size(cfg%face_positions)-1
    allocate(cfg%z_faces(0:n),cfg%z_quadrature(8,n))
    do i=0,n
      cfg%z_faces(i)=elevation(cfg%face_positions(i)*cfg%Lr)
    enddo
    do i=1,n
      left=cfg%face_positions(i-1)*cfg%Lr;right=cfg%face_positions(i)*cfg%Lr
      do j=1,8
        s=.5_wp*(left+right)+.5_wp*(right-left)*qx(j)
        cfg%z_quadrature(j,i)=elevation(s)
      enddo
    enddo
  contains
    real(wp) function elevation(s)result(z)
      real(wp),intent(in)::s
      real(wp)::arc
      z=-min(s,cfg%Lp)*sin(cfg%beta)
      if(s>cfg%Lp)then
        arc=s-cfg%Lp
        z=z+arc**2/(sqrt(cfg%CA**2+arc**2)+cfg%CA)
      endif
    end function
  end subroutine
  subroutine initialize_single_phase_hydrostatic(cfg,phase,u)
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase
    real(wp),intent(inout)::u(:,:)
    real(wp)::c2,p0,rho0,mu,rhot,zref
    integer::i,lo,hi,cell
    call coefficients(cfg,phase,c2,p0,rho0,mu)
    rhot=cfg%nPs/cfg%nCg**2
    if(phase==1)rhot=cfg%nrho_l0+(cfg%nPs-cfg%nP_l0)/cfg%nCl**2
    zref=cfg%z_faces(ubound(cfg%z_faces,1));lo=first_physical_cell(cfg);hi=last_physical_cell(cfg)
    do i=lo,hi
      cell=i-cfg%num_dummy_cells
      u(i,:)=0
      ! Mean normalized by the stored weights: constants are integrated exactly.
      u(i,phase)=rhot*sum(qw*exp(-cfg%G*(cfg%z_quadrature(:,cell)-zref)/(cfg%omega_u**2*c2)))/sum(qw)
    enddo
  end subroutine
  integer function global_single_phase(u,cfg)result(phase)
    real(wp),intent(in)::u(:,:)
    type(sim_config_t),intent(in)::cfg
    integer::lo,hi
    phase=0
    if(trim(cfg%single_phase_mode)/='REDUCED')return
    lo=first_physical_cell(cfg);hi=last_physical_cell(cfg)
    if(.not.all(ieee_is_finite(u(lo:hi,:))))return
    if(all(u(lo:hi,2)==0._wp).and.all(u(lo:hi,1)>0._wp))phase=1
    if(all(u(lo:hi,1)==0._wp).and.all(u(lo:hi,2)>0._wp))phase=2
  end function
  subroutine coefficients(cfg,phase,c2,p0,rho0,mu)
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase
    real(wp),intent(out)::c2,p0,rho0,mu
    real(wp)::pip
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    if(phase==1)then
      c2=pip*cfg%nCl**2;p0=cfg%nP_l0;rho0=cfg%nrho_l0;mu=cfg%mu_l
    else
      c2=pip*cfg%nCg**2;p0=0;rho0=0;mu=cfg%mu_g
    endif
  end subroutine
  real(wp) function pressure(rho,cfg,phase)result(p)
    real(wp),intent(in)::rho
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase
    p=rho*cfg%nCg**2
    if(phase==1)p=cfg%nP_l0+(rho-cfg%nrho_l0)*cfg%nCl**2
  end function
  subroutine hydro_factors(cfg,phase,i,fl,fr,status)
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase,i
    real(wp),intent(out)::fl,fr
    integer,intent(out)::status
    real(wp)::c2,p0,rho0,mu,zl,zr,zm,z(8),ex(8),average
    integer::j,cell
    fl=1;fr=1;status=0
    if(trim(cfg%source_scheme)=='POINT')return
    if(trim(cfg%source_scheme)/='HYDROSTATIC')then
      status=PHASE_INVALID;return
    endif
    call coefficients(cfg,phase,c2,p0,rho0,mu)
    cell=i-cfg%num_dummy_cells
    if(allocated(cfg%z_faces))then
      if(cell<1.or.cell>size(cfg%z_faces)-1)then
        status=PHASE_GEOMETRY_REQUIRED;return
      endif
      zl=cfg%z_faces(cell-1);zr=cfg%z_faces(cell)
    else
      ! A constant inclination in one cell defines its elevation difference.
      zl=-.5_wp*cfg%cell_width(i)*cfg%Lr*sin(cfg%thetav(i));zr=-zl
    endif
    zm=.5_wp*(zl+zr)
    do j=1,8
      z(j)=zm+.5_wp*(zr-zl)*qx(j)
    enddo
    if(allocated(cfg%z_quadrature))z=cfg%z_quadrature(:,cell)
    ex=-cfg%G*(z-zm)/(cfg%omega_u**2*c2)
    if(maxval(abs(ex))>log(huge(1._wp))-2)then
      status=PHASE_INVALID;return
    endif
    ! Normalize by the stored weights (0.5*sum(qw)=1+2.2e-16 in binary64):
    ! G=0 then gives fl=fr=1 exactly, never a roundoff-signed outlet velocity.
    average=sum(qw*exp(ex))/sum(qw)
    fl=exp(-cfg%G*(zl-zm)/(cfg%omega_u**2*c2))/average
    fr=exp(-cfg%G*(zr-zm)/(cfg%omega_u**2*c2))/average
    if(.not.all(ieee_is_finite([fl,fr])).or.min(fl,fr)<=0)status=PHASE_INVALID
  end subroutine
  subroutine physical_flux(w,cfg,phase,f,status)
    real(wp),intent(in)::w(2)
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase
    real(wp),intent(out)::f(2)
    integer,intent(out)::status
    real(wp)::p,pip
    status=PHASE_INVALID;f=ieee_value(0._wp,ieee_quiet_nan)
    if(.not.all(ieee_is_finite(w)))return
    if(w(1)<=0)return
    p=pressure(w(1),cfg,phase)
    if(p<=0)return
    pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
    f=[w(2),w(2)**2/w(1)+pip*p];status=0
  end subroutine
  subroutine single_phase_face(u,old,cfg,phase,i,f,status)
    real(wp),intent(in)::u(:,:),old(:,:)
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase,i
    real(wp),intent(out)::f(3)
    integer,intent(out)::status
    real(wp)::wl(2),wr(2),ol(2),orr(2),lf,rf,f1(2),f2(2),ab(2,2),aa(2,2),v,c,c2,p0,rho0,mu
    call hydro_factors(cfg,phase,i,lf,rf,status);if(status/=0)return
    wl=[u(i,phase)*rf,u(i,3)];ol=[old(i,phase)*rf,old(i,3)]
    call hydro_factors(cfg,phase,i+1,lf,rf,status);if(status/=0)return
    wr=[u(i+1,phase)*lf,u(i+1,3)];orr=[old(i+1,phase)*lf,old(i+1,3)]
    call physical_flux(wl,cfg,phase,f1,status);if(status/=0)return
    call physical_flux(wr,cfg,phase,f2,status);if(status/=0)return
    if(trim(cfg%dissipation_mode)/='LAGGED')then
      ol=wl;orr=wr
    endif
    if(min(ol(1),orr(1))<=0)then
      status=PHASE_INVALID;return
    endif
    v=(ol(2)/sqrt(ol(1))+orr(2)/sqrt(orr(1)))/(sqrt(ol(1))+sqrt(orr(1)))
    call coefficients(cfg,phase,c2,p0,rho0,mu);c=sqrt(c2)
    aa(1,:)=[0._wp,1._wp];aa(2,:)=[c2-v*v,2*v]
    ! Exact |A| of the physical TWO-wave barotropic system, no third eigenvalue.
    if(trim(cfg%flux_scheme)=='RUSANOV')then
      ab=0;ab(1,1)=abs(v)+c;ab(2,2)=ab(1,1)
    else if(v>=c)then
      ab=aa
    else if(v<=-c)then
      ab=-aa
    else
      ab=(v/c)*aa;ab(1,1)=ab(1,1)+c-v*v/c;ab(2,2)=ab(2,2)+c-v*v/c
    endif
    f1=.5_wp*(f1+f2)-.5_wp*matmul(ab,wr-wl)
    f=0;f(phase)=f1(1);f(3)=f1(2);status=0
  end subroutine

  subroutine single_phase_boundary(u,cfg,phase,inlet,f,status,ub)
    real(wp),intent(in)::u(:,:)
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase
    logical,intent(in)::inlet
    real(wp),intent(out)::f(3)
    integer,intent(out)::status
    real(wp),intent(out),optional::ub(3)
    real(wp)::fl,fr,rho,v,c2,c,p0,rho0,mu,rhob,vb,target,left,right,mid,gleft,gright,gm,val(2),tmp(3),lower
    integer::i,j,st
    f=ieee_value(0._wp,ieee_quiet_nan)
    if(present(ub))ub=f
    i=last_physical_cell(cfg);if(inlet)i=first_physical_cell(cfg)
    status=PHASE_INVALID
    if(.not.all(ieee_is_finite(u(i,:))))return
    call coefficients(cfg,phase,c2,p0,rho0,mu);c=sqrt(c2)
    call hydro_factors(cfg,phase,i,fl,fr,status);if(status/=0)return
    rho=u(i,phase)*fr;if(inlet)rho=u(i,phase)*fl
    v=u(i,3)/rho
    if(abs(v)>=c)then
      status=PHASE_DATA_REQUIRED;return
    endif
    if(inlet)then
      if((phase==1.and.((cfg%gas_inlet_basis=='MASS'.and.cfg%mug/=0).or. &
        (cfg%gas_inlet_basis=='LOCAL_VOLUME'.and.cfg%njg/=0))).or. &
        (phase==2.and.((cfg%liquid_inlet_basis=='MASS'.and.cfg%mul/=0).or. &
        (cfg%liquid_inlet_basis=='LOCAL_VOLUME'.and.cfg%njl/=0))))then
        status=PHASE_TRANSITION_REQUIRED;return
      endif
      if((phase==1.and.cfg%liquid_inlet_basis=='LOCAL_VOLUME').or. &
        (phase==2.and.cfg%gas_inlet_basis=='LOCAL_VOLUME'))then
        vb=cfg%njl;if(phase==2)vb=cfg%njg
        val(1)=(vb-v)/c
        if(vb>v)then
          rhob=rho*(.5_wp*(val(1)+sqrt(val(1)**2+4)))**2
        else
          rhob=rho*exp(val(1))
        endif
      else
        target=cfg%mul;if(phase==2)target=cfg%mug
        ! Admissible branch of the boundary state: subsonic |vb|<c and p>0.
        ! On it m(r)=r*vb(r) is strictly increasing (dm/dr=vb+c on the
        ! rarefaction curve, >=vb+c on the compression curve), so a root exists
        ! iff m(left)<target<m(right). The sonic ends are closed form:
        !   rarefaction vb=-c: r=rho*exp(-1-v/c), r/rho in (exp(-2),1);
        !   compression vb=+c: sqrt(r/rho)-sqrt(rho/r)=(c-v)/c.
        ! The thermodynamic vacuum p=0 is an open bound, never a mass floor.
        lower=0
        if(phase==1)lower=cfg%nrho_l0-cfg%nP_l0/cfg%nCl**2
        if(rho<=lower)then
          status=PHASE_INVALID;return
        endif
        val(1)=(c-v)/c
        right=rho*(.5_wp*(val(1)+sqrt(val(1)**2+4)))**2
        left=max(rho*exp(-1-v/c),lower)
        gleft=mass_residual(left);gright=mass_residual(right)
        if(.not.ieee_is_finite(gleft).or..not.ieee_is_finite(gright))then
          status=PHASE_INVALID;return
        endif
        if(gleft>=0.or.gright<=0)then
          ! Proven absence of a subsonic, positive-pressure boundary state.
          status=PHASE_DATA_REQUIRED;return
        endif
        gm=gleft;mid=left
        do j=1,200
          mid=left+.5_wp*(right-left);gm=mass_residual(mid)
          if(mid==left.or.mid==right.or.gm==0)exit
          if(gm<0)then
            left=mid;gleft=gm
          else
            right=mid;gright=gm
          endif
        enddo
        if(gm/=0)then
          ! Bracket exhausted at adjacent floating-point numbers.
          mid=left
          if(abs(gright)<abs(gleft))mid=right
        endif
        rhob=mid;vb=target/rhob
      endif
    else
      rhob=cfg%nPs/cfg%nCg**2
      if(phase==1)rhob=cfg%nrho_l0+(cfg%nPs-cfg%nP_l0)/cfg%nCl**2
      if(rhob>rho)then
        vb=v-c*(rhob-rho)/sqrt(rhob*rho)
      else
        vb=v+c*log(rho/rhob)
      endif
    endif
    if(abs(vb)>=c)then
      status=PHASE_DATA_REQUIRED;return
    endif
    if(.not.inlet.and.vb<0)then
      if(cfg%outlet_reservoir_phase==0)then
        status=PHASE_DATA_REQUIRED;return
      else if(cfg%outlet_reservoir_phase/=phase)then
        status=PHASE_TRANSITION_REQUIRED;return
      endif
    endif
    call physical_flux([rhob,rhob*vb],cfg,phase,val,status);if(status/=0)return
    if(inlet)then
      ! Use the imposed mass component exactly; velocity is its same state.
      if((phase==1.and.cfg%liquid_inlet_basis=='MASS').or.(phase==2.and.cfg%gas_inlet_basis=='MASS'))then
        val(1)=target
      endif
    endif
    f=0;f(phase)=val(1);f(3)=val(2)
    tmp=0;tmp(phase)=rhob;tmp(3)=rhob*vb
    if(present(ub))ub=tmp
  contains
    real(wp) function mass_residual(r)result(value)
      real(wp),intent(in)::r
      real(wp)::vv
      if(r>rho)then
        vv=v+c*(r-rho)/sqrt(r*rho)
      else
        vv=v+c*log(r/rho)
      endif
      value=r*vv-target
    end function
  end subroutine

  subroutine single_phase_assemble(u,old,dt,cfg,phase,r,a,b,c,norm,fin,fout,status,jac,base)
    real(wp),intent(in)::u(:,:),old(:,:),dt
    type(sim_config_t),intent(in)::cfg
    integer,intent(in)::phase
    real(wp),intent(out)::r(:,:),a(:,:,:),b(:,:,:),c(:,:,:),norm,fin(3),fout(3),base(2,2)
    integer,intent(out)::status
    logical,intent(in)::jac
    real(wp)::xp(size(u,1),3),xm(size(u,1),3),rp(size(r,1),size(r,2)),rm(size(r,1),size(r,2))
    real(wp)::fi(3),fo(3),bb(2,2),h,sc(3),d(3),err
    integer::lo,hi,i,j,k,l,t,cp,cm,idx(2)
    lo=first_physical_cell(cfg);hi=last_physical_cell(cfg);idx=[phase,3]
    a=0;b=0;c=0;norm=huge(1._wp)
    r=ieee_value(0._wp,ieee_quiet_nan);fin=r(:,1);fout=fin;base=fin(1)
    status=PHASE_INVALID
    if(.not.all(ieee_is_finite(u(lo:hi,:))).or..not.all(ieee_is_finite(old(lo:hi,:))))return
    if(global_single_phase(old,cfg)/=phase)then
      status=PHASE_TRANSITION_REQUIRED;return
    endif
    call residual(u,r,fin,fout,base,status);if(status/=0)return
    norm=0
    do i=lo,hi
      sc=[max(old(i,1),1e-3_wp),max(old(i,2),1e-7_wp),max(abs(old(i,3)),1._wp)]
      norm=max(norm,maxval(abs(r(:,i-lo+1))/sc))
    enddo
    if(.not.jac)return
    do i=lo,hi
      ! An exact inactive constraint, not a fabricated physical wave/Jacobian.
      b(3-phase,3-phase,i-lo+1)=1
      do t=1,2
        k=idx(t);h=1e-5_wp*max(abs(u(i,k)),1e-2_wp)
        if(k==phase)h=min(h,.2_wp*u(i,k))
        xp=u;xm=u;xp(i,k)=u(i,k)+h;xm(i,k)=u(i,k)-h
        call residual(xp,rp,fi,fo,bb,cp);call residual(xm,rm,fi,fo,bb,cm)
        if(cp/=0.or.cm/=0)then
          status=PHASE_INVALID;return
        endif
        do j=max(lo,i-1),min(hi,i+1)
          l=j-lo+1;d=(rp(:,l)-rm(:,l))/(2*h)
          if(j==i)b(:,k,l)=d
          if(j==i-1)c(:,k,l)=d
          if(j==i+1)a(:,k,l)=d
        enddo
      enddo
    enddo
  contains
    subroutine residual(x,rr,fi,fo,fb,st)
      real(wp),intent(in)::x(:,:)
      real(wp),intent(out)::rr(:,:),fi(3),fo(3),fb(2,2)
      integer,intent(out)::st
      real(wp)::f(3),fl,fr,rho,v,pip,p0,rho0,mu,c2,gravity,friction,re,ff,pL,pR
      integer::ii,ll,code
      st=PHASE_INVALID;rr=0;fb=ieee_value(0._wp,ieee_quiet_nan)
      call coefficients(cfg,phase,c2,p0,rho0,mu);pip=cfg%omega_P/(cfg%omega_rho*cfg%omega_u**2)
      do ii=lo,hi
        ll=ii-lo+1;rho=x(ii,phase)
        if(rho<=0.or.x(ii,3-phase)/=0.or.pressure(rho,cfg,phase)<=0)then
          st=PHASE_INVALID;return
        endif
        v=x(ii,3)/rho;rr(:,ll)=x(ii,:)-old(ii,:)
        call hydro_factors(cfg,phase,ii,fl,fr,st);if(st/=0)return
        if(trim(cfg%source_scheme)=='HYDROSTATIC')then
          pL=pressure(rho*fl,cfg,phase);pR=pressure(rho*fr,cfg,phase)
          gravity=pip*(pR-pL)/cfg%cell_width(ii)
        else
          gravity=-rho*cfg%G*cfg%Lr/cfg%omega_u**2*sin(cfg%thetav(ii))
        endif
        re=rho*cfg%omega_rho*abs(v)*cfg%omega_u*cfg%D_H/mu
        friction=-32*cfg%Lr*mu*v/(cfg%omega_rho*cfg%omega_u*cfg%D_H**2)
        if(re>=2000)then
          ff=ffan(cfg%eps/cfg%D_H,re,code)
          if(code/=0)then
            st=PHASE_INVALID;return
          endif
          friction=-2*cfg%Lr/cfg%D_H*ff*rho*v*abs(v)
        endif
        rr(3,ll)=rr(3,ll)-dt*(gravity+friction)
      enddo
      do ii=lo,hi-1
        call single_phase_face(x,old,cfg,phase,ii,f,st);if(st/=0)return
        ll=ii-lo+1
        rr(:,ll)=rr(:,ll)+dt*f/cfg%cell_width(ii)
        rr(:,ll+1)=rr(:,ll+1)-dt*f/cfg%cell_width(ii+1)
        if(ii==lo+cfg%N_pipe-1)then
          fb(:,1)=f(:2);fb(:,2)=f(:2)
        endif
      enddo
      call single_phase_boundary(x,cfg,phase,.true.,fi,st);if(st/=0)return
      call single_phase_boundary(x,cfg,phase,.false.,fo,st);if(st/=0)return
      rr(:,1)=rr(:,1)-dt*fi/cfg%cell_width(lo)
      rr(:,hi-lo+1)=rr(:,hi-lo+1)+dt*fo/cfg%cell_width(hi)
      st=0
    end subroutine
  end subroutine
end module
