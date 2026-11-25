program test_junction_vanishing_gas
 !! Contract of the riser face trace when the junction gas flux vanishes
 !! (Bhagwat riser closure at the catenary base, theta=0, ascending liquid):
 !!  * j_g -> 0+ : admissible roots, alpha -> 0+ continuously, alpha/j_g bounded
 !!    (drift-flux limit alpha ~ j_g/(C0 j + Ud));
 !!  * j_g = 0   : refused; zero flux alone never selects the pure-liquid trace;
 !!  * j_g < 0   : CLOSURE_DOMAIN (12), outside the co-current closure.
 !! This is the documented local limit of the present formulation: the
 !! mixed->pure passage of a face trace is not representable.
 use mod_kinds,only:wp
 use module_steady_state,only:solve_alpha_from_superficials
 use module_closure_law,only:CLOSURE_DOMAIN
 implicit none
 real(wp)::orho,ou,rl,rg,jl,a,aprev,ratio,rprev,jg
 integer::st,k,checks
 checks=0;orho=997.05_wp;ou=1.16_wp
 rl=(997.05_wp+(270600._wp-101325._wp)/1498._wp**2)/orho;rg=(270600._wp/343._wp**2)/orho
 jl=.55_wp/ou;aprev=1;rprev=-1
 do k=4,14,2
  jg=10._wp**(-k)/ou
  a=solve_alpha_from_superficials(jl,jg,rl,rg,0._wp,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp, &
    .89e-3_wp,1.85e-5_wp,.0728_wp,ou,orho,1e-12_wp,status=st,closure_id=3)
  ratio=a/jg
  write(*,'(a,es10.2,a,i4,a,es14.6,a,es14.6)')'j_g(m/s)',10._wp**(-k),' status',st,' alpha',a,' alpha/j_g',ratio
  call check(st==0.and.a>0.and.a<aprev,'admissible root, alpha decreasing to 0+')
  if(rprev>0)call check(abs(ratio-rprev)<=1e-3_wp*rprev,'alpha/j_g converges (continuous drift-flux limit)')
  aprev=a;rprev=ratio
 enddo
 a=solve_alpha_from_superficials(jl,0._wp,rl,rg,0._wp,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp, &
   .89e-3_wp,1.85e-5_wp,.0728_wp,ou,orho,1e-12_wp,status=st,closure_id=3)
 write(*,'(a,i4)')'j_g=0 status',st
 call check(st/=0,'exact zero gas flux does not select a pure-liquid trace')
 a=solve_alpha_from_superficials(jl,-1e-12_wp,rl,rg,0._wp,.05254_wp,acos(-1._wp)*.05254_wp**2/4,4.6e-5_wp,9.81_wp, &
   .89e-3_wp,1.85e-5_wp,.0728_wp,ou,orho,1e-12_wp,status=st,closure_id=3)
 write(*,'(a,i4)')'j_g<0 status',st
 call check(st==CLOSURE_DOMAIN,'negative junction gas flux under ascending liquid: CLOSURE_DOMAIN')
 write(*,'(a,i0,a)')'PASS junction vanishing gas: ',checks,' assertions'
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
