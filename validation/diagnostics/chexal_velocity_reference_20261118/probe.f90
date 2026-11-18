program probe_velocity_derivative
  use mod_kinds,only:wp
  use module_closure_law,only:CHEXAL_1992,Dlei_fechamento_or_ndim_simp
  implicit none
  real(wp),parameter::a=.4_wp,rl=998.0_wp,rg=2.3_wp,ug=.5_wp,theta=1.151159_wp
  real(wp),parameter::d=.05254_wp,g=9.81_wp,ml=.001_wp,mg=.000018_wp,sigma=.072_wp
  real(wp),parameter::pi=acos(-1.0_wp),area=pi*d*d/4
  real(wp),parameter::vel(2)=[2.2586e-8_wp,1e-14_wp]
  real(wp),parameter::uscale(3)=[1.0_wp,3.0_wp,.2_wp],rscale(3)=[1.0_wp,1000.0_wp,10.0_wp]
  real(wp)::fa,fg,fl,fug,ful,wu,wr
  integer::k,j,st
  write(*,'(a)')'ul_si,wu,wr,status,dR_dul'
  do k=1,3
    wu=uscale(k);wr=rscale(k)
    do j=1,2
      call Dlei_fechamento_or_ndim_simp(a,rl/wr,rg/wr,vel(j)/wu,ug/wu,theta,d,area,0.0_wp,g,ml,mg,sigma, &
        wu,wr,1e-12_wp,fa,fg,fl,fug,ful,CHEXAL_1992,st)
      if(st==0)then
        write(*,'(es24.16,a,es24.16,a,es24.16,a,i0,a,es24.16)')vel(j),',',wu,',',wr,',',st,',',ful
      else
        write(*,'(es24.16,a,es24.16,a,es24.16,a,i0,a)')vel(j),',',wu,',',wr,',',st,',UNDEFINED'
      endif
    enddo
  enddo
end program probe_velocity_derivative
