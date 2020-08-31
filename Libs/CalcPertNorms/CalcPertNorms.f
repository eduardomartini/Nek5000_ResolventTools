      subroutine CalcPertNorms
      include 'SIZE'
      include 'TSTEP'           ! ISTEP
      include 'INPUT'           ! PARAM
      include 'SOLN'            ! V[XYZ]
      include 'GEOM'            ! xm1,ym1
      include 'MASS'
      integer ntot , ipert
      real  pertnorm(lpert) 
      ntot = nx1*ny1*nz1*nelt
      

        if (nio==0 .and. istep ==0) then !clear/create file on the first loop
          open(unit=88, file='norm_perturbations.dat',
     $      status="replace", action="write")
          close(88)
        endif   

      do ipert = 1,npert
        pertnorm(ipert)=0 
        do i=1,ntot
          pertnorm(ipert)= pertnorm(ipert) +  bm1(i,1,1,1)*
     $       (VXP(i,ipert)**2.0+VYP(i,ipert)**2.0)
          if (ndim>2) then
            pertnorm(ipert)= pertnorm(ipert) +  bm1(i,1,1,1)*
     $          VZP(i,pert)**2.0
          end if
        enddo
        pertnorm(ipert) = glsum(pertnorm(ipert),1)
        pertnorm(ipert) = sqrt(pertnorm(ipert))
      enddo
      if(nio==0) then
        write(*,*) 'Perturbation magnitude ' , pertnorm
        open(unit=88, file='norm_perturbations.dat',
     $            status='old',action='write',position="append")
        write(88,*) time , pertnorm       
        close(88)    
      endif
      end subroutine

