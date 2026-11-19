program probe_ccfl_limits
 use mod_kinds,only:wp
 use module_chexal_1992
 implicit none
 real(wp),parameter::rl=998.0_wp,rg=1.2_wp,d=.05254_wp,g=9.81_wp,ml=.001_wp,mg=.000018_wp,s=.072_wp
 real(wp),parameter::gas(6)=[.1_wp,.01_wp,1e-4_wp,1e-6_wp,1e-8_wp,1e-12_wp]
 real(wp)::jl,a,r,c,v,scale,ul,ug
 integer::i,st,k
 write(*,'(a)')'jg,status,jl_CCFL,alpha_CCFL,R,normalized_R'
 do i=1,6
   call chexal_ccfl(gas(i),rl,rg,d,g,ml,mg,s,jl,a,st)
   if(st==0)then
     call chexal_residual_si(a,rl,rg,jl,gas(i),acos(-1.0_wp)/2,d,g,ml,mg,s,r,k,branch=1,jl_ccfl=jl)
     call chexal_coefficients_si(a,rl,rg,jl,gas(i),acos(-1.0_wp)/2,d,g,ml,mg,s,c,v,k,branch=1,jl_ccfl=jl)
     scale=max(abs(a*c*(jl+gas(i))),abs(a*v),abs(gas(i)))
     write(*,'(es24.16,a,i0,4(a,es24.16))')gas(i),',',st,',',jl,',',a,',',r,',',r/scale
   else
     write(*,'(es24.16,a,i0)')gas(i),',',st
   endif
 enddo
 call chexal_ccfl(.1_wp,rl,rg,d,g,ml,mg,s,jl,a,st)
 do i=-1,1,2
   call chexal_recover_velocities(a,rl,rg,rl*jl+rg*.1_wp,acos(-1.0_wp)/2,d,g,ml,mg,s,ul,ug,st,branch=i)
   write(*,*)'recover_at_CCFL branch,status=',i,st
   if(st==0)write(*,*)'recovered jl,jg=',(1-a)*ul,a*ug
 enddo
end program
