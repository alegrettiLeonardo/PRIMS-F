program test_initial_inventory
 use mod_kinds,only:wp
 use module_steady_state
 implicit none
 real(wp)::rho_scale,p_scale,u_scale,cl,cg,p,rl,rg,a,ul,ug,jl,jg,kin,gl,gg,f,expected,rho,velocity,c_si,mu
 integer::phase,sgn,st,checks
 checks=0;rho_scale=997.05_wp;p_scale=101325;u_scale=1
 cl=1498*sqrt(rho_scale/p_scale);cg=sqrt(287*293.15_wp)*sqrt(rho_scale/p_scale);p=2.2_wp
 do phase=1,2
  do sgn=-1,1,2
   gl=0;gg=0
   if(phase==1)gl=sgn*.0001_wp
   if(phase==2)gg=sgn*.0001_wp
   call state_at(p,-.1_wp,gl,gg,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp,cl,cg,1._wp,1._wp, &
      .00089_wp,.0000185_wp,.0728_wp,u_scale,rho_scale,1e-12_wp,1,rl,rg,a,ul,ug,jl,jg,kin,st,phase_hint=phase)
   call check(st==0,'explicit inventory has exact pure state for either velocity sign')
   call check(a==real(phase-1,wp),'exact void endpoint belongs to declared inventory')
   rho=rl*rho_scale;velocity=ul*u_scale;c_si=1498;mu=.00089_wp
   if(phase==2)then
    rho=rg*rho_scale;velocity=ug*u_scale;c_si=sqrt(287*293.15_wp);mu=.0000185_wp
   endif
   call check(abs(kin-(gl*ul+gg*ug))<1e-18_wp,'same state supplies kinetic momentum')
   f=fun_dpds(.1_wp,p,gl,gg,5._wp,10._wp,3._wp,.1_wp,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp, &
       cl,cg,1._wp,1._wp,.00089_wp,.0000185_wp,.0728_wp,p_scale,u_scale,rho_scale,1e-12_wp,st,3,1,1,phase)
   expected=10/p_scale*(-rho*9.81_wp*sin(-.1_wp)-32*mu*velocity/.05254_wp**2)/(1-velocity**2/c_si**2)
   call check(st==0.and.abs(f-expected)<1e-10_wp*max(abs(expected),1e-10_wp),'pure stationary momentum agrees with independent SI reduced equation')
  enddo
 enddo
 call state_at(p,0._wp,0._wp,.0001_wp,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp,cl,cg,1._wp,1._wp, &
   .00089_wp,.0000185_wp,.0728_wp,u_scale,rho_scale,1e-12_wp,1,rl,rg,a,ul,ug,jl,jg,kin,st)
 if(st==0)call check(a<1,'zero liquid flux does not infer liquid absence')
 call state_at(p,0._wp,.0001_wp,.0001_wp,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp,cl,cg,1._wp,1._wp, &
   .00089_wp,.0000185_wp,.0728_wp,u_scale,rho_scale,1e-12_wp,1,rl,rg,a,ul,ug,jl,jg,kin,st,phase_hint=PHASE_GAS_ONLY)
 call check(st==STEADY_INVALID,'declared gas inventory cannot accept nonzero liquid feed')
 write(*,'(a,i0,a)')'PASS initial inventory: ',checks,' assertions'
contains
 subroutine check(ok,label)
  logical,intent(in)::ok
  character(*),intent(in)::label
  if(.not.ok)then
   print *,'FAIL ',label;error stop 1
  endif
  checks=checks+1
 end subroutine
end program
