program test_chexal_np3989
  !! EPRI NP-3989-SR Rev.1 (1986), Appendix A, Table A-1 (PDF p.100), row by row,
  !! with the Fortran core in its STEAM_WATER verification mode (Eq.3-5/3-6).
  !! Columns: case | version | properties | reference | Fortran | error | cause | status.
  !! Properties are declared (the report prints none): classical saturated steam
  !! tables (14.7 psia vf .01672, vg 26.80; 1000 psia vf .02159, vg .4456 ft3/lb),
  !! sigma IAPWS 1994 at Tsat, viscosities from the PRINTED Re (so Re_f, Re_g are
  !! the printed ones), p_crit 3208.2 psia. The independent transcription is
  !! validation/diagnostics/chexal_riser_20261121/np3989_audit.py (V1992 values
  !! below, full precision). Causes are those demonstrated in that audit:
  !!   rows 1,4,8 reproduced (property band); row 5 reproduced only by the 1986
  !!   downflow equations (version revision); row 3 violates its own drift
  !!   identity beyond rounding (defective); rows 2,6,7 unexplained -> FAIL.
  !! The test fails if the Fortran core departs from the independent 1992
  !! transcription or if any row changes classification.
  use mod_kinds, only: wp
  use module_chexal_1992
  implicit none
  real(wp), parameter :: FT=0.3048_wp, LB=16.01846337_wp, G=9.80665_wp, PC=3208.2_wp
  real(wp) :: tab(9,8), py(3,8), rf, rg, sig, d, jf, jg, ml, mg, a, c0, v, e(3), etab(3)
  integer :: k, st, nr, fails, n
  character(len=40) :: cause(8)
  character(len=34) :: status_txt(8)
  real(wp), parameter :: tolband(3)=[5e-4_wp, 1.2e-2_wp, 5e-4_wp]   ! C0, Vgj, alpha: property band of the audit
  ! case P D jf jg Ref Reg | C0 Vgj alpha (table)
  tab(:,1)=[14.7_wp,.05_wp,5._wp,10._wp,78684._wp,2297._wp,1.2037_wp,.5979_wp,.5361_wp]
  tab(:,2)=[1000._wp,.05_wp,5._wp,10._wp,184219._wp,87257._wp,1.1116_wp,.1410_wp,.5947_wp]
  tab(:,3)=[1000._wp,1._wp,5._wp,10._wp,3684867._wp,1745362._wp,1.1119_wp,.2234_wp,.5914_wp]
  tab(:,4)=[14.7_wp,1._wp,5._wp,10._wp,1573883._wp,45947._wp,1.1922_wp,.9054_wp,.5323_wp]
  tab(:,5)=[14.7_wp,1._wp,-5._wp,-10._wp,-1573883._wp,-45947._wp,1.0400_wp,5.1172_wp,.9538_wp]
  tab(:,6)=[1000._wp,1._wp,-5._wp,-10._wp,-3684867._wp,-1745362._wp,1.0498_wp,4.8444_wp,.9171_wp]
  tab(:,7)=[1000._wp,.05_wp,-5._wp,-10._wp,-184219._wp,-87257._wp,1.2948_wp,1.8775_wp,.5700_wp]
  tab(:,8)=[14.7_wp,.05_wp,-5._wp,-10._wp,-78684._wp,-2297._wp,1.3036_wp,4.6590_wp,.6714_wp]
  py(:,1)=[1.20375272217602824_wp,0.598411418028930941_wp,0.536057872467044394_wp]
  py(:,2)=[1.11066056634767985_wp,0.137348124928460880_wp,0.595335311363152542_wp]
  py(:,3)=[1.11089345166989406_wp,0.226265922356566840_wp,0.592077960274875803_wp]
  py(:,4)=[1.19221237529317081_wp,0.906286270196161814_wp,0.532212935703436418_wp]
  py(:,5)=[1.58692842440878845_wp,11.3908858200494443_wp,0.805604393450245260_wp]
  py(:,6)=[1.09206903156341362_wp,4.63743079328305452_wp,0.851527301228673261_wp]
  py(:,7)=[1.29337630791095504_wp,1.81925320443212857_wp,0.568783196073148423_wp]
  py(:,8)=[1.30365557489340911_wp,4.65842701591130126_wp,0.671302835878019755_wp]
  cause=[character(len=40)::'reproduced (both versions)','unexplained, not a version difference', &
    'printed row violates drift identity','reproduced (both versions)','1992 downflow revision; 1986 reproduces', &
    'unexplained, not a version difference','unexplained, not a version difference','reproduced (both versions)']
  fails=0
  write(*,'(a)') 'case | version | properties | reference C0/Vgj/a | Fortran C0/Vgj/a | max rel err (ref) | cause | status'
  do k=1,8
    if(tab(1,k)<100)then
      rf=LB/.01672_wp;rg=LB/26.80_wp;sig=sigma(373.12_wp)
    else
      rf=LB/.02159_wp;rg=LB/.4456_wp;sig=sigma(557.93_wp)
    endif
    d=tab(2,k)*FT;jf=tab(3,k)*FT;jg=tab(4,k)*FT
    ml=rf*abs(jf)*d/abs(tab(5,k));mg=rg*abs(jg)*d/abs(tab(6,k))
    call chexal_solve_alpha(jf,jg,rf,rg,acos(-1.0_wp)/2,d,G,ml,mg,sig,a,st,nroots=nr, &
      fluid=CHEXAL_FLUID_STEAM_WATER,pressure=tab(1,k),pcrit=PC)
    call chk(st==CHEXAL_OK.and.nr==1,'single steam-water root, case',k)
    call chexal_coefficients_si(a,rf,rg,jf,jg,acos(-1.0_wp)/2,d,G,ml,mg,sig,c0,v,st, &
      fluid=CHEXAL_FLUID_STEAM_WATER,pressure=tab(1,k),pcrit=PC)
    call chk(st==CHEXAL_OK,'coefficients status, case',k)
    e=abs([c0,v/FT,a]-py(:,k))/abs(py(:,k))
    call chk(maxval(e)<1e-11_wp,'Fortran = independent 1992 transcription, case',k)
    etab=abs([c0,v/FT,a]-tab(7:9,k))/tab(7:9,k)
    status_txt(k)='FAIL'
    if(all(etab<=tolband))status_txt(k)='PASS'
    if(k==5)status_txt(k)='PASS (1986 version; 1992 differs)'
    if(k==3)status_txt(k)='FAIL (reference row inconsistent)'
    write(*,'(i2,a,f6.1,a,3f9.4,a,3f9.4,a,es9.2,a,a,a,a)') k,' | 1992 | steam-water P=',tab(1,k),' psia |', &
      tab(7:9,k),' |',c0,v/FT,a,' | ',maxval(etab),' | ',trim(cause(k)),' | ',trim(status_txt(k))
    ! classification must stay what the audit demonstrated
    select case(k)
    case(1,4,8);call chk(all(etab<=tolband),'reproduced row stays inside the property band, case',k)
    case(5);call chk(etab(2)>0.5_wp,'1992 downflow revision differs from the 1986 row, case',k)
    case default;call chk(.not.all(etab<=tolband),'1000 psia row stays outside the band (FAIL recorded), case',k)
    end select
  enddo
  ! default (air-water) path unchanged and invalid fluid inputs refused
  call chexal_coefficients_si(.5_wp,998._wp,1.2_wp,.5_wp,.6_wp,.9_wp,.05254_wp,9.81_wp,1e-3_wp,1.8e-5_wp,.072_wp,c0,v,st)
  call chexal_coefficients_si(.5_wp,998._wp,1.2_wp,.5_wp,.6_wp,.9_wp,.05254_wp,9.81_wp,1e-3_wp,1.8e-5_wp,.072_wp,e(1),e(2),n, &
    fluid=CHEXAL_FLUID_AIR_WATER)
  call chk(st==0.and.n==0.and.c0==e(1).and.v==e(2),'explicit AIR_WATER equals default bit for bit',0)
  call chexal_coefficients_si(.5_wp,998._wp,1.2_wp,.5_wp,.6_wp,.9_wp,.05254_wp,9.81_wp,1e-3_wp,1.8e-5_wp,.072_wp,c0,v,st, &
    fluid=CHEXAL_FLUID_STEAM_WATER)
  call chk(st==CHEXAL_INVALID,'steam-water without pressure is invalid',0)
  call chexal_coefficients_si(.5_wp,998._wp,1.2_wp,.5_wp,.6_wp,.9_wp,.05254_wp,9.81_wp,1e-3_wp,1.8e-5_wp,.072_wp,c0,v,st, &
    fluid=CHEXAL_FLUID_STEAM_WATER,pressure=PC,pcrit=PC)
  call chk(st==CHEXAL_INVALID,'steam-water at p>=p_crit is invalid',0)
  call chexal_coefficients_si(.5_wp,998._wp,1.2_wp,.5_wp,.6_wp,.9_wp,.05254_wp,9.81_wp,1e-3_wp,1.8e-5_wp,.072_wp,c0,v,st, &
    fluid=7)
  call chk(st==CHEXAL_INVALID,'unknown fluid is invalid',0)
  if(fails>0)then
    write(*,'(a,i0)') 'FAILED checks: ',fails;error stop 1
  endif
  write(*,'(a)') 'PASS NP-3989 Table A-1 row audit (8 rows; classification fixed; not experimental validation)'
contains
  real(wp) function sigma(T)
    real(wp),intent(in)::T
    real(wp)::tau
    tau=1-T/647.096_wp
    sigma=0.2358_wp*tau**1.256_wp*(1-0.625_wp*tau)
  end function
  subroutine chk(ok,label,k)
    logical,intent(in)::ok
    character(*),intent(in)::label
    integer,intent(in)::k
    if(.not.ok)then
      fails=fails+1;write(*,'(a,1x,a,1x,i0)') 'FAIL:',label,k
    endif
  end subroutine
end program
