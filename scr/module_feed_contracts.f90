module module_feed_contracts
  !! SI input contracts. No thermodynamic gamma is used in an isothermal EOS.
  use mod_kinds,only:wp
  use module_config,only:sim_config_t
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  real(wp),parameter,public::FLOW_UNSET=-huge(1.0_wp)
  integer,parameter,public::FEED_INVALID=120,FEED_CONFLICT=121,FEED_REQUIRED=122
  type,public::feed_contract_t
    character(16)::liquid_mode='MASS',gas_mode='MASS'
    real(wp)::liquid_mass=0,gas_mass=0,liquid_j=0,gas_j=0
  end type
  public::gas_eos_coefficient,resolve_feed,inlet_mass_fluxes
contains
  subroutine gas_eos_coefficient(eos,R,T,c_fixed,c,status)
    character(*),intent(in)::eos
    real(wp),intent(in)::R,T,c_fixed
    real(wp),intent(out)::c
    integer,intent(out)::status
    status=FEED_INVALID;c=0
    select case(trim(eos))
    case('ISOTHERMAL')
      if(.not.all(ieee_is_finite([R,T])))return
      if(min(R,T)<=0)return
      c=sqrt(R*T)
    case('BAROTROPIC')
      if(.not.ieee_is_finite(c_fixed))return
      if(c_fixed<=0)return
      c=c_fixed
    case default
      return
    end select
    status=0
  end subroutine

  subroutine resolve_feed(lb,gb,jl,jg,ml,mg,rho_l0,ps,cg,R,Tstd,Pstd,feed,status)
    character(*),intent(in)::lb,gb
    real(wp),intent(in)::jl,jg,ml,mg,rho_l0,ps,cg,R,Tstd,Pstd
    type(feed_contract_t),intent(out)::feed
    integer,intent(out)::status
    status=FEED_INVALID
    if(.not.all(ieee_is_finite([jl,jg,ml,mg,rho_l0,ps,cg,R,Tstd,Pstd])))return
    if(min(rho_l0,ps,cg,R,Tstd,Pstd)<=0)return
    feed%liquid_j=jl;feed%gas_j=jg
    select case(trim(lb))
    case('LEGACY')
      feed%liquid_mass=rho_l0*jl
      if(ml>=0)feed%liquid_mass=ml
    case('MASS')
      if(ml==FLOW_UNSET)then
        status=FEED_REQUIRED;return
      endif
      feed%liquid_mass=ml
    case('LOCAL_VOLUME')
      if(ml/=FLOW_UNSET)then
        status=FEED_CONFLICT;return
      endif
      feed%liquid_mode='LOCAL_VOLUME';feed%liquid_mass=rho_l0*jl ! IC guess only
    case default
      return
    end select
    select case(trim(gb))
    case('STANDARD')
      feed%gas_mass=Pstd/(R*Tstd)*jg
      if(mg>=0)feed%gas_mass=mg ! explicit historical mass override
    case('SEPARATOR')
      feed%gas_mass=ps/(cg*cg)*jg
      if(mg>=0)feed%gas_mass=mg ! explicit historical mass override
    case('MASS')
      if(mg==FLOW_UNSET)then
        status=FEED_REQUIRED;return
      endif
      feed%gas_mass=mg
    case('LOCAL_VOLUME')
      if(mg/=FLOW_UNSET)then
        status=FEED_CONFLICT;return
      endif
      feed%gas_mode='LOCAL_VOLUME';feed%gas_mass=ps/(cg*cg)*jg ! IC guess only
    case default
      return
    end select
    status=0
  end subroutine

  subroutine inlet_mass_fluxes(p,cfg,gl,gg,status)
    real(wp),intent(in)::p
    type(sim_config_t),intent(in)::cfg
    real(wp),intent(out)::gl,gg
    integer,intent(out)::status
    real(wp)::rl,rg
    gl=cfg%mul;gg=cfg%mug;status=FEED_INVALID
    if(.not.all(ieee_is_finite([p,cfg%nCl,cfg%nCg,cfg%mul,cfg%mug,cfg%njl,cfg%njg])))return
    if(min(p,cfg%nCl,cfg%nCg)<=0)return
    rl=cfg%nrho_l0+(p-cfg%nP_l0)/cfg%nCl**2;rg=p/cfg%nCg**2
    if(min(rl,rg)<=0)return
    select case(trim(cfg%liquid_inlet_basis))
    case('MASS')
    case('LOCAL_VOLUME');gl=rl*cfg%njl
    case default;return
    end select
    select case(trim(cfg%gas_inlet_basis))
    case('MASS')
    case('LOCAL_VOLUME');gg=rg*cfg%njg
    case default;return
    end select
    status=0
  end subroutine
end module
