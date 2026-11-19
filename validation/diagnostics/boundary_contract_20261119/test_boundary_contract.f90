program test_boundary_contract
 use mod_kinds,only:wp
 use boundary_audit_api,only:ba_setup,ba_boundary,ba_assemble
 implicit none
 real(wp)::par(20),u(6,3),old(6,3),theta(6),width(6),r1(3,6),r2(3,6),jac(18,18),fi(3),fo(3),base(2,2)
 real(wp)::bc(3),f(3),p,rl,rg,err,alpha,pi
 integer::models(2),i,st,idx
 character(32)::mode
 call get_command_argument(1,mode)
 par=[1.0_wp,997.05_wp,101325.0_wp,12.546253993197219_wp,12.546253993197219_wp,57.4_wp, &
 .05254_wp,4.6e-5_wp,9.81_wp,.00089_wp,.0000185_wp,.0728_wp,1498*sqrt(997.05_wp/101325), &
 343*sqrt(997.05_wp/101325),1.0_wp,1.0_wp,2.1_wp,.5_wp,.001_wp,.5_wp]
 models=2;call ba_setup(par,models)
 do i=1,6
  pi=2.1_wp+.001_wp*(i-1);alpha=.25_wp+.01_wp*(i-1)
  rl=1+(pi-1)/par(13)**2;rg=pi/par(14)**2
  u(i,:)=[(1-alpha)*rl,alpha*rg,.6_wp+.003_wp*(i-1)]
 enddo
 if(trim(mode)=='contract')then
  call ba_boundary(u(1,:),0.0_wp,1,0,f,bc,p,st)
  if(st/=0)error stop 'inlet domain'
  rl=1+(p-1)/par(13)**2
  err=abs(f(1)/rl/par(20)-1)
  write(*,'(a,es24.16)')'relative prescribed-Q error=',err
  if(err>1e-10_wp)error stop 'FAIL prescribed liquid volumetric flux'
 elseif(trim(mode)=='partition')then
  old=u;theta=0;width=1
  call ba_assemble(6,2,u,old,theta,width,.001_wp,0,r1,jac,fi,fo,base,st,idx)
  if(st/=0)error stop 'first partition inadmissible'
  call ba_assemble(6,4,u,old,theta,width,.001_wp,0,r2,jac,fi,fo,base,st,idx)
  if(st/=0)error stop 'second partition inadmissible'
  err=maxval(abs(r1-r2));write(*,'(a,es24.16)')'partition residual difference=',err
  if(err>1e-12_wp)error stop 'FAIL identical-model artificial partition'
 else
  error stop 'expected contract or partition'
 endif
 write(*,'(a)')'PASS '//trim(mode)
end program
