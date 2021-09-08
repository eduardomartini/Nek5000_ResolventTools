
      include 'FFT.f'      
      include 'SavePerturbations.f'
      include 'Load_FLD_To.f'
      include 'FLD2Force.f'
      include 'MultiHarm.f'


c========================================================
c     SSRM Module
c
c     Contains function to be added to .usr file in order to use the
c     SSRM.
c
c     Subroutines
c           SSRM_Init    : Initialize the module.
c                   Should be called in userdat. 
c           SSRM_IC      : Sets initial conditions. 
c                   Should be called in useric. 
c                  read (used in exit conditions).
c           SSRM_userchk : Sets up initialization, force and FT updates.
c                   Should be called in userchk    
c           SSRM_Force   : Creates external forcing terms.
c                   Should be called in userf  
c
c========================================================

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-------------------------- SSRM Module --------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      module SSRM
      integer currIter,iFFtSteps,nHarmonics,iTransTime,iFFT_sampling
      real forceNorm


      contains
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-------------------------- SSRM Functions -----------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine SSRM_Init
      use FFT
      use SavePerturb
      use MultiHarm
      implicit none
      include 'SIZE'
      include 'TSTEP'           ! ISTEP
      include 'INPUT'           ! PARAM
      include 'SOLN'            ! V[XYZ]
      include 'GEOM'            ! xm1,ym1
      include 'ADJOINT'
      include 'MASS'
      real resParams(5)
      real,external::glsum
      integer ntot,ntot2,i,j , ipert
      character(len=100) refFile 
      logical file_exists 
      real tmpFreq(2000)

      integer,save ::  pertIO

      real B(lx1,ly1,lz1,lelt,3),C(lx1,ly1,lz1,lelt,3)
      COMMON /SSRM_BC/ B,C

      
      open (unit=98, file='params.input',
     $                status='old',action='read')
      read(98,*) resParams ! currIter,transTime,maxf,nHarmonics
      close(98)
      currIter      = int(resParams(1))
      iTransTime    = int(resParams(2)/dt)
      nHarmonics    = int(resParams(4))
      iFFtSteps     = int(   1/(resParams(3)*dt) )*(nHarmonics-1)
      forceNorm     = resParams(5)
      iFFT_sampling = iFFtSteps
      
      if (NIO==0) write(*,*) 'Looking for FFT sampling freq...'
      do iFFT_sampling=int(1/(resParams(3)*dt) )+1,1,-1
        if (NIO==0) write(*,*) iFFT_sampling,
     $   mod(iFFtSteps,iFFT_sampling),1.0/(iFFT_sampling*dt),
     $   resParams(3)
        if (  ( mod(iFFtSteps,iFFT_sampling)==0   )  .and.
     $        ( (1.0/(iFFT_sampling*dt))>resParams(3)*2 ) ) exit
      enddo
      
      if (NIO==0) then
      write(*,*) 'dt',dt
      write(*,*) 'steps ',iFFtSteps,iFFT_sampling
      write(*,*) 'times ',iFFtSteps*dt,iFFT_sampling*dt
      write(*,*) 'freqs low:',1/(iFFtSteps*dt),
     $ ' samp :',1/(iFFT_sampling*dt)
      write(*,*) 'tar freqs',resParams(3)/(nHarmonics-1),resParams(3)
      write(*,*) 'npoints' , 1.0*iFFtSteps/iFFT_sampling
      endif





      end subroutine 

c-------------------------- useric call --------------------------------
      subroutine SSRM_IC(ux,uy,uz) 
      include 'SIZE'
      include 'TOTAL'
      
      real ux,uy,uz

      real B(lx1,ly1,lz1,lelt,3),C(lx1,ly1,lz1,lelt,3)
      COMMON /SSRM_BC/ B,C

      if (currIter<0 .and. jp > 0) then
        call RANDOM_NUMBER(ux)
        call RANDOM_NUMBER(uy)
        call RANDOM_NUMBER(uz)
        ux = ux*2-1
        uy = uy*2-1
        uz = uz*2-1
      else
        ux  =0.0 
        uy  =0.0
        uz  =0.0
      endif
      end subroutine

c-------------------------- userchk call  -----------------------------
      subroutine SSRM_userchk
      use FFT
      use SavePerturb
      use MultiHarm
      implicit none
      include 'SIZE'
      include 'TSTEP'           ! ISTEP
      include 'INPUT'           ! PARAM
      include 'SOLN'            ! V[XYZ]
      include 'GEOM'            ! xm1,ym1
      include 'ADJOINT'
      include 'MASS'
      real,external::glsum
      integer ntot,ntot2,i,j , ipert
      real pertNorm
      character(len=100) refFile 
      logical file_exists 
      real tmpFreq(2000)

      integer,save ::  pertIO

      real B(lx1,ly1,lz1,lelt,3),C(lx1,ly1,lz1,lelt,3)
      COMMON /SSRM_BC/ B,C


      ntot = nx1*ny1*nz1*nelv
      ntot2= nx2*ny2*nz2*nelv
                      
      if (istep==0) then


      ntot = nx1*ny1*nz1*nelv
      ntot2= nx2*ny2*nz2*nelv
                      
      if (istep==0) then
        if (NIO==0) write(*,*) "Initializing Routines"
        !Remove end run conditions from reaFile
        FINTIM=1e9
        NSTEPS=1e9
        param(10:11)=1e9

        ! change IO as to not save baseflow
        pertIO = IOSTEP
        IOSTEP = 0
        param(14) = 0
        param(15) = 0

        !Save Mass Matrix for posprocessing
        ifxyo = .true.
        call outpost2(bm1,bm1,bm1,bm2,bm1,0,'bm1')

        ! Setup B and C matrices 
        ! If files exists, read them, otherwise initialize to 1.
        refFile='B.fld';
        INQUIRE(FILE=refFile, EXIST=file_exists)
        if (file_exists) then
            if (NIO==0) write(*,*) 'Reading B matrix from file...'
            refFile='B.fld'
            call Load_FLD_To(refFile,
     $        B(:,:,:,:,1),
     $        B(:,:,:,:,2),
     $        B(:,:,:,:,3))
            if (NIO==0) write(*,*) '    B matrix read...'
        else
            if (NIO==0) write(*,*) 'Setting B as identity'
            B=1
        endif
        INQUIRE(FILE="C.fld", EXIST=file_exists)
        if (file_exists) then
            if (NIO==0) write(*,*) 'Reading C matrix from file...'
            refFile='C.fld'
            call Load_FLD_To(refFile,
     $           C(:,:,:,:,1),
     $           C(:,:,:,:,2),
     $           C(:,:,:,:,3))
            if (NIO==0) write(*,*) '   C matrix read...'
        else
            if (NIO==0) write(*,*) 'Setting B as identity'
            C=1
        endif
        if (NIO==0) write(*,*) 'Matrices B and C set ! '

        ! Set as direct or adjoint run
        if ( modulo(currIter,2)  == 0) then
          ifadj = .false.
        else
          ifadj = .true.
        endif

        ! initialize target frequencies as the first harmonics
        do i=1,nHarmonics
          tmpFreq(i)=(i-1)*1.0 /(iFFtSteps*dt) !first harmonic
        enddo
        call InitFFts_Array(tmpFreq,nHarmonics) ! FFT module initialization

        ! initialize harmonic shapes 
        if (currIter==0) then
          ! Random initialization
          call MultiHarm_RandomInitShapes(tmpFreq,nHarmonics)
          !Pre-multiply by B remove zero gain components
          call MultiHarm_Multi(B(:,:,:,:,1),B(:,:,:,:,2),B(:,:,:,:,3))
            

        else if (currIter>0) then
          ! From File
          refFile = 'harmCompList.txt'
          call MultiHarm_LoadShapes(refFile,.false.)
        endif 
        call MultiHarm_Normalize(forceNorm)
        call MultiHarm_OutPost()


        ! FFT module setups
        if (.not. ifadj) then
          call FFts_SetOffSet(itransTime*dt)
          call MultiHarm_SetOffSet(itransTime*dt)
        else
          call FFts_SetOffSet(iTransTime*dt+iFFtSteps*dt*1.0)
          call MultiHarm_SetOffSet(iTransTime*dt+iFFtSteps*dt*1.0)
        endif    

        call MultiHarm_FFTConvMinus( .not. ifadj )   ! q(t) =  Real( (\sum a exp(-i ometa t))  = \sum ar cos(wt)+ai sin(wt) )
        call FFTs_FFTConvMinus(   .not. ifadj )       ! a    = \int q(t)exp(iwt) dt = \int q(t)(cos(wt)+i sin(wt)) dt
        call FFts_SetNormTime(iFFtSteps*dt*1.0)
        call FFts_SetTrapzIntRule(.false.)


        ! Apply C/B matrix to the output for adj/dir run
        if (ifadj) then
            call MultiHarm_Multi(C(:,:,:,:,1),C(:,:,:,:,2),C(:,:,:,:,3))
        else
            call MultiHarm_Multi(B(:,:,:,:,1),B(:,:,:,:,2),B(:,:,:,:,3))
        endif
        
      endif        
      endif

      ! calculate current norm and save to disk
      if (NIO==0) write(*,*) 'Computing run Norm...'
      
      pertNorm = 0
      do i=1,ntot 
        pertNorm = pertNorm + ( vxp(i,1)**2 + vyp(i,1)**2
     $                         )*bm1(i,1,1,1)
        if (ndim==3) pertNorm = pertNorm + 
     $                 vzp(i,1)**2*bm1(i,1,1,1)
      enddo
      pertNorm = sqrt(glsum(pertNorm,1))
     
      ! IO Operations   
      if (pertIO/=0 .and. mod(ISTEP,abs(pertIO))==0) 
     $            call SavePerturbationsNow('i')

      if ( (pertIO/=0 .and. mod(ISTEP,abs(pertIO) )==0 ) .or.  
     $     (pertIO==0 .and. mod(ISTEP,   1000    )==0 )  ) then
       if (NIO==0) then   
          if (istep == 0) then
            open (unit=99, file='runNorm.txt',
     $          status='replace',action='write')
          else
            open (unit=99, file='runNorm.txt',
     $      position="append" ,status='old',action='write')
         endif
        
          write(99,*) time , pertNorm , j 
          close(99)
        endif ! NIO==0
      endif ! ISTEP

      !update ffts every output step
      if (NIO==0) write(*,*) 'TIME:',iTransTime*dt,iFFtSteps

      if (istep>iTransTime .and. currIter >= 0 )  then
        if (NIO==0) write(*,*) 'Updating FFts...'
        if (mod(istep-iTransTime-1,iFFT_sampling )==0) then
              call UpdateFFts(dt*iFFT_sampling)
        endif
        
        ! End of run!
        if (istep==iTransTime+iFFtSteps) then
            ! Apply B/C matrix to the output for adj/dir run
            if (ifadj) then
              call FFt_Multi(1,1,B(:,:,:,:,1),B(:,:,:,:,2),B(:,:,:,:,3))
            else
              call FFt_Multi(1,1,C(:,:,:,:,1),C(:,:,:,:,2),C(:,:,:,:,3))
            endif

            ! Write data
            call OutPostFFts(1,1)
            if (NIO==0) write(*,*) 'Completed Run, exitting...'
            if (NIO==0) write(*,*) currIter
            call exitt
        endif
      endif
        
      if (currIter<0 .and. istep>1 .and.  pertNorm < 1e-6) then
        if (NIO==0) write(*,*) 'Completed Transient time run, exitting.'
        if (NIO==0) write(*,*) pertNorm
        call exitt 
      endif    


      return
      end


      end module

      subroutine SSRM_Force(ix,iy,iz,e,time,ffx,ffy,ffz)
        use MultiHarm
        real ffx,ffy,ffz,time
        integer ix,iy,iz,e

        if (currIter<0) then
            ffx=0
            ffy=0
            ffz=0
        else 

            ffx = MultiHarm_GetSignal(ix,iy,iz,e,1,time)
            ffy = MultiHarm_GetSignal(ix,iy,iz,e,2,time)
            ffz = MultiHarm_GetSignal(ix,iy,iz,e,3,time)
        endif


      end subroutine




