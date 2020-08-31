      module MultiHarm
      implicit none
      integer,parameter :: nMaxHarm = 50
      integer,save      :: nHarm 
      integer,save      :: FourSignCon = 1 
      real, save        :: freqs(nMaxHarm),offSetTime=0;
      real,parameter    :: MH_pi = 3.1415926535897932384626433832795 


      contains

      function  MultiHarm_GetSignal(ix,iy,iz,e,ui,time)
        include 'MultiHarm_DEF'
        integer ix,iy,iz,e,ui
        real MultiHarm_GetSignal
        real arg(nMaxHarm),time,factor(nMaxHarm)
        real,save :: ccos(nMaxHarm),csin(nMaxHarm)
        real,save :: prevTime = -1e99



        ! compute sin and cos and saves it until time changes
        if (time /= prevTime ) then
          factor=2.
          factor(1)=1
C           arg(1:nHarm) = 2*MH_pi*freqs(1:nHarm)*(time-offSetTime)
C           ccos(1:nHarm) = cos(arg(1:nHarm))*factor(1:nHarm)
C           csin(1:nHarm) = sin(arg(1:nHarm))*factor(1:nHarm)*FourSignCon
          ccos(1:nHarm) = cos(2*MH_pi*freqs(1:nHarm)*(time-offSetTime))
     $                     *factor(1:nHarm)
          csin(1:nHarm) = sin(2*MH_pi*freqs(1:nHarm)*(time-offSetTime))
     $                     *factor(1:nHarm)*FourSignCon
          prevTime = time 
        endif 

        if (ui<=ndim) then
          MultiHarm_GetSignal = sum(
     $    cosShapes(ix,iy,iz,e,ui,1:nHarm)*ccos(1:nHarm)+
     $    sinShapes(ix,iy,iz,e,ui,1:nHarm)*csin(1:nHarm) )
        else
          MultiHarm_GetSignal = 0
        end if


      endfunction

      subroutine  MultiHarm_LoadShapes(harmParamsFile,normalize_in)
        include 'MultiHarm_DEF'
        include 'TOTAL'
        integer,parameter :: fileUnit = 88
        real dummy(lx1,ly1,lz1,lelt),timebak
        real norm
        character(len=100):: harmParamsFile
        character(len=100):: currFilename
        logical, intent(in),optional :: normalize_in
        logical normalize
        integer i
        real,external :: glsum

        if ( present(normalize_in))  then
          normalize=normalize_in
        else 
          normalize = .false.
        endif

        open (unit=fileUnit, file=harmParamsFile,
     $         status='old',action='read')
        timebak = time 

        ! Reads total number of harm components
        read(fileUnit,*)  nHarm
        ! Reads real (cos) parts
        do i=1,nHarm
          read(fileUnit,'(A)')  currFilename
          if(NIO==0) then
            write(*,*) 'MultiHarm: Real Part file ' , currFilename
          endif 
          if (ndim==2) then
            call Load_FLD_Time_To_wTime(currFilename,
     $        cosShapes(:,:,:,:,1,i),
     $        cosShapes(:,:,:,:,2,i),
     $        dummy,dummy,dummy)
          else
            call Load_FLD_Time_To_wTime(currFilename,
     $        cosShapes(:,:,:,:,1,i),
     $        cosShapes(:,:,:,:,2,i),
     $        cosShapes(:,:,:,:,3,i),
     $        dummy,dummy)     
          end if
          freqs(i) = time
        enddo

        ! Reads imag (sin) parts
        do i=1,nHarm
          read(fileUnit,'(A)')  currFilename
          if(NIO==0) then
            write(*,*) 'MultiHarm: Imag Part file ' , currFilename
          endif
          if (ndim==2) then
            call Load_FLD_Time_To_wTime(currFilename,
     $        sinShapes(:,:,:,:,1,i),
     $        sinShapes(:,:,:,:,2,i),
     $        dummy,dummy,dummy)
          else
            call Load_FLD_Time_To_wTime(currFilename,
     $        sinShapes(:,:,:,:,1,i),
     $        sinShapes(:,:,:,:,2,i),
     $        sinShapes(:,:,:,:,ndim,i),
     $        dummy,dummy)   
          end if  
          if (freqs(i) /= time) then
            write(*,*) 'Error loading multiHarm Shapes: '
     $     ,'Missmatch between real and imaginary part frequencies... ', 
     $          freqs(i) , time, '.  Aborting!'
          call exitt
          end if
        enddo

        time = timebak 

        if (NIO==0) then
        do i=1,nHarm
          write(*,*) 'MultiHarm freqs :',i,freqs(i)
        enddo
        endif
        close(fileUnit)
      endsubroutine


      subroutine  MultiHarm_RandomInitShapes(frequencies,nfreqs)
        include 'MultiHarm_DEF'
        include 'TOTAL'
        integer nfreqs
        real frequencies(nfreqs),norm
        integer i,j,k
        real xx,yy,zz
        real, external :: glsum

        nHarm = nfreqs
        freqs(1:nfreqs) = frequencies(1:nfreqs)

        if (NIO==0) write(*,*) 'Random initialization of Harm comp '
C #IF .false.
C         ! Uses fortran pseudo random 
C         write(*,*) 'Building pseudo-Random start condition'
C         do i=1,nx1*ny1*nz1*nelt
C           xx=xm1(i,1,1,1)
C           yy=ym1(i,1,1,1)
C           zz=zm1(i,1,1,1)
C         do j=1,ndim
C           do k=1,nfreqs
C             cosShapes(:,:,:,:,j,k)=cos( 
C      $    j*3.e4*(xm1*ym1)**2-1.5e3*xm1*ym1 +85.e3*xm1+zm1+1e3*zm1**2.)
C             sinShapes(:,:,:,:,j,k)=cos( 
C      $    j*5.e4*(xm1*ym1)**2-6.5e3*xm1*ym1 +35.e3*xm1+zm1+4e3*zm1**2.)
C           enddo
C         enddo

C         enddo
C #ELSE
        ! Uses fortran built-in random number generator
        write(*,*) 'Building Random start condition with Fortran rand'
        call RANDOM_NUMBER(cosShapes)
        call RANDOM_NUMBER(sinShapes)
        cosShapes=cosShapes*2.-1.
        sinShapes=sinShapes*2.-1.
        do j=1,ndim
          do k=1,nfreqs
             call dsavg(cosShapes(:,:,:,:,j,k))
             call dsavg(sinShapes(:,:,:,:,j,k))
          enddo
        enddo
C #ENDIF  
        if (NIO==0) write(*,*) 'Done! '
        if (freqs(1)==0) then        
          sinShapes(:,:,:,:,:,1)=0
        endif
      endsubroutine

      subroutine MultiHarm_Normalize(targetNorm)
        include 'MultiHarm_DEF'
        include 'TOTAL'
        real targetNorm
        real, external :: glsum
        integer k
        real norm
        if (NIO==0) write(*,*) 'Normalizing Harmonic Components'
        do k=1,nHarm
          norm = sum( (cosShapes(:,:,:,:,1,k)**2. + 
     $                 cosShapes(:,:,:,:,2,k)**2. + 
     $                 sinShapes(:,:,:,:,1,k)**2. + 
     $                 sinShapes(:,:,:,:,2,k)**2.    )*bm1)

          if (ndim==3) then
            norm = norm+ sum( (cosShapes(:,:,:,:,ndim,k)**2. + 
     $                         sinShapes(:,:,:,:,ndim,k)**2. )*bm1) 
          endif

          norm=sqrt(glsum(norm,1))
          cosShapes(:,:,:,:,:,k)=cosShapes(:,:,:,:,:,k)*targetNorm/norm
          sinShapes(:,:,:,:,:,k)=sinShapes(:,:,:,:,:,k)*targetNorm/norm
        enddo
        if (freqs(1)==0) then        
          sinShapes(:,:,:,:,:,1)=0
        endif
      endsubroutine


      subroutine  MultiHarm_RandomInitShapesFromFile(freqsFile)
        integer :: fileUnit = 88
        real currFreqs(nMaxHarm)
        integer currNHarm
        logical endOfFile
        integer Reason
        character(len=100) freqsFile
        open (unit=fileUnit, file=freqsFile,
     $                     status='old',action='read')
        currNHarm = 1
        do while (.not. endOfFile)
          if (currNHarm>nMaxHarm) then
            write(*,*) 'Error initializing multiHarm from file' ,
     $      ' increase nMaxHarm! Aborting'
            call exitt
          end if

          read(fileUnit,*,IOSTAT=Reason) currFreqs(currNHarm)
          IF (Reason > 0)  THEN
            write(*,*) "Error reading ", freqsFile 
            endOfFile = .true.
          else if (Reason < 0) THEN
            endOfFile = .true.
          ELSE
            currNHarm = currNHarm+1
          END IF
        end do
        close(fileUnit)
        call MultiHarm_RandomInitShapes(currFreqs,currNHarm)
      endsubroutine

      subroutine MultiHarm_FFTConvMinus(sign)
        logical sign
        if (sign)  then
          FourSignCon = -1
        else
          FourSignCon = +1
        endif
      endsubroutine



      subroutine MultiHarm_OutPost()
        include 'MultiHarm_DEF'
        include 'TOTAL'
        CHARACTER(LEN=3)  :: prefix
        integer ntot,i,j,k 
        real timesave
        logical ifxyo_bak
      !------ Writes Output Files --------------------
        ifxyo_bak = ifxyo
        timesave = time
        
c       Output perturbation states
        do k= 1,nHarm
          time = freqs(k)
          ifxyo = (k==1)
          write(prefix,'(A2,I0.1)') 'fc', 1
          call outpost2(cosShapes(1,1,1,1,1,k)
     $                 ,cosShapes(1,1,1,1,2,k)
     $                 ,cosShapes(1,1,1,1,max(2,ndim),k)
     $                 ,pr,tp,0,prefix)



         write(prefix,'(A2,I0.1)') 'fs', 1
          call outpost2(sinShapes(1,1,1,1,1,k)
     $                 ,sinShapes(1,1,1,1,2,k)
     $                 ,sinShapes(1,1,1,1,max(2,ndim),k)
     $                 ,pr,tp,0,prefix)
        end do
        time = timesave
        ifxyo= ifxyo_bak
        return
      end subroutine      

      subroutine MultiHarm_SetOffSet(offSettime_in)
        real offsetTime_in
        offsetTime =offSettime_in
        return
      endsubroutine

      end module MultiHarm
