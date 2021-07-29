      module FFT
      implicit none
      integer, parameter:: nmaxFreq       = 051
      integer           :: nFreqs
      real              :: frequencies(nmaxFreq)
      real              :: lastFreqTime 
      real,save         :: offsetTime = 0
      integer,save      :: FourSignCon = 1 
      real,save         :: normalizatinTime=1;
      logical, save     :: trapzIntRule       = .true.
      real,parameter    :: FFT_pi = 3.1415926535897932384626433832795 


      contains 
c========================================================
c     InitFFts
c
c     Initializes FFT module, reads frequencis file and
c     initializes variables
c     freqsFile : file containing frequencies to be computed
c     list : array of n,7 for returning data
c     maxListSize : size of n (lines of list)
c     nLines : number of lines read from file
c========================================================
      subroutine InitFFts_Array(frequencies_in,nFreqs_in)
      include 'SIZE'
      include 'TOTAL'
      include 'FFT_DEF'
      integer nFreqs_in,i
      real frequencies_in(nFreqs_in)
      vx_FFTs_c(:,:) = 0
      vx_FFTs_s(:,:) = 0
      vy_FFTs_c(:,:) = 0
      vy_FFTs_s(:,:) = 0
      vz_FFTs_c(:,:) = 0
      vz_FFTs_s(:,:) = 0

      vxp_FFTs_c(:,:,:) = 0
      vxp_FFTs_s(:,:,:) = 0
      vyp_FFTs_c(:,:,:) = 0
      vyp_FFTs_s(:,:,:) = 0
      vzp_FFTs_c(:,:,:) = 0
      vzp_FFTs_s(:,:,:) = 0
      lastFreqTime =  -1e99
      if (nFreqs_in>nmaxFreq) then
        write(*,*) 'FFT Init Failed: Increase nMaxFreq ! Aborting'
        call exitt
      endif
      nFreqs = nFreqs_in
      frequencies(1:nFreqs) = frequencies_in(1:nFreqs)

      if(NIO==0) then
            do i=1,nFreqs
            write(*,'(I0.2,A,f8.5)') i,
     $           " FFT Freqs : " , frequencies(i)
            enddo
      end if

      end subroutine


      subroutine InitFFts(freqsFile)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FFT_DEF'

      character(LEN=100) freqsFile
      logical endOfFile  
      integer reason,i,lineCount

c    ---------------- Read Input Sensors List   ------------------------------
      if(NIO==0) write(*,*) "Shapes: Frequencies File :" ,
     $ freqsFile 

      lineCount = 0
      endOfFile = .false.
      open (unit=98, file=freqsFile,
     $                     status='old',action='read')
      do while (.not. endOfFile)
            read(98,*,IOSTAT=Reason) frequencies(lineCount+1)
            IF (Reason > 0)  THEN
              write(*,*) "Error reading ", freqsFile 
              endOfFile = .true.
            else if (Reason < 0) THEN
              endOfFile = .true.
            ELSE
             lineCount = lineCount+1
            END IF
      end do
      close(98)
      nFreqs = lineCount

      if(NIO==0) then
            do i=1,lineCount
            write(*,'(I0.2,A,f8.5)') i,
     $           " FFT Freqs : " , frequencies(i)
            enddo
      end if
      vx_FFTs_c(:,:) = 0
      vx_FFTs_s(:,:) = 0
      vy_FFTs_c(:,:) = 0
      vy_FFTs_s(:,:) = 0
      vz_FFTs_c(:,:) = 0
      vz_FFTs_s(:,:) = 0

      vxp_FFTs_c(:,:,:) = 0
      vxp_FFTs_s(:,:,:) = 0
      vyp_FFTs_c(:,:,:) = 0
      vyp_FFTs_s(:,:,:) = 0
      vzp_FFTs_c(:,:,:) = 0
      vzp_FFTs_s(:,:,:) = 0
      lastFreqTime =  -1e99

      return
      end subroutine 

c========================================================
c     UpdateFFts
c
c     Updates the integral for FFT calculations
c========================================================


      subroutine UpdateFFts(deltaT,jpStart_in, jpEnd_in)
      include 'SIZE'
      include 'TOTAL'
      include 'FFT_DEF'
      integer,intent(in),optional :: jpStart_in,jpEnd_in
      integer :: jpStart,jpEnd
      real deltaT,weight
      logical,save :: firstCall=.true.
      real, external ::glmax

      integer ntot,i,j,k 
      real weigth ,c(nMaxFreq) ,s(nMaxFreq) 
      ntot = nx1*ny1*nz1*nelt
      jpStart=0
      jpEnd  =npert
      
      if(present(jpStart_in)) jpStart=jpStart_in
      if(present(jpEnd_in  )) jpEnd  =jpEnd_in
c    ---------------- Read Input Sensors List   ------------------------------
      if(NIO==0) write(*,*) "Updating FFTs ..."


      ! Weight for trapesoidal rule. equals 0.5 to boundary terms and 1 otherwise
      ! For now forward integration rule
      if (trapzIntRule) then
        if (firstCall) then
          weight = .5    
          firstCall=.false.
        else
          weight = 1.    
        endif
      else
        weight = 1.    
        firstCall=.false.
      end if

      if (lastep==1)  then
        weight = 0.5
        call OutPostFFts()
        return
      endif
      weight = weight/normalizatinTime

      if (NIO==0) write(*,*) 'FFT Integration dt : ',
     $                    deltaT,time

      c(1:nFreqs) =cos(2*FFT_pi*frequencies(1:nFreqs)*(time-offsetTime))
      s(1:nFreqs) =sin(2*FFT_pi*frequencies(1:nFreqs)*(time-offsetTime))
     $                        *FourSignCon
      if (jpStart==0) then
        do i= 1,ntot
          ! Base flow (jp==0)
            ! Cosine transform
          vx_FFTs_c(i,:)= vx_FFTs_c(i,:) + Vx(i,1,1,1)*c*deltaT*weight
          vy_FFTs_c(i,:)= vy_FFTs_c(i,:) + Vy(i,1,1,1)*c*deltaT*weight
             ! Sine transform
          vx_FFTs_s(i,:)= vx_FFTs_s(i,:) + Vx(i,1,1,1)*s*deltaT*weight
          vy_FFTs_s(i,:)= vy_FFTs_s(i,:) + Vy(i,1,1,1)*s*deltaT*weight
            ! if 3d
          if ( ndim == 3) then
            vz_FFTs_c(i,:)= vz_FFTs_c(i,:) + Vz(i,1,1,1)*c*deltaT*weight
            vz_FFTs_s(i,:)= vz_FFTs_s(i,:) + Vz(i,1,1,1)*s*deltaT*weight
          end if
        end do ! i
      endif ! jpstart

      do k= 1,nFreqs
        ! Perturbations (jp>0)
        i=max(1,jpStart)
        j=jpEnd
        ! Cosine transform
        vxp_FFTs_c(:,i:j,k)= vxp_FFTs_c(:,i:j,k) +
     $        VXP(:,i:j)*c(k)*deltaT*weight
        vyp_FFTs_c(:,i:j,k)= vyp_FFTs_c(:,i:j,k) +
     $        VYP(:,i:j)*c(k)*deltaT*weight
         ! Sine transform
        vxp_FFTs_s(:,i:j,k)= vxp_FFTs_s(:,i:j,k) +
     $        VXP(:,i:j)*s(k)*deltaT*weight
        vyp_FFTs_s(:,i:j,k)= vyp_FFTs_s(:,i:j,k) +
     $        VYP(:,i:j)*s(k)*deltaT*weight
        ! if 3d
        if ( ndim == 3) then
          vzp_FFTs_c(:,i:j,k)= vzp_FFTs_c(:,i:j,k) +
     $        VZP(:,i:j)*c(k)*deltaT*weight
          vzp_FFTs_s(:,i:j,k)= vzp_FFTs_s(:,i:j,k) +
     $        VZP(:,i:j)*s(k)*deltaT*weight
        end if
      enddo ! k
      return
      end subroutine
      

      subroutine OutPostFFts( iOutputStart_, iOutputEnd_)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FFT_DEF'
      integer,intent(in),optional :: iOutputStart_
      integer,intent(in),optional :: iOutputEnd_
      integer iOutputStart,iOutputEnd
      CHARACTER(LEN=3)  :: prefix
      integer ntot,i,j,k 
      real timesave, freqNorms(0:lpert,nmaxFreq)
      logical ifxyo_bak
      real, external :: glsum
      !default values for start
      iOutputStart=0
      iOutputEnd  =npert
      
      if(present(iOutputStart_)) iOutputStart=iOutputStart_
      if(present(iOutputEnd_  )) iOutputEnd  =iOutputend_
      
c     ------ Writes Output Files --------------------
      ifxyo_bak = ifxyo

      timesave = time
      if (NIO==0)  then
          write(*,*) "Outposting FFts, for ", 
     $             npert, "perturbations and ", nFreqs, 'frequencies.'
      end if
        
c       Output perturbation states
      ifxyo_bak = ifxyo
      do k= 1,nFreqs
        time = frequencies(k)
        ifxyo = (k==1)
        if (iOutputStart==0) then
          write(prefix,'(A1,I0.2)') 'c', 0
          call outpost2(vx_FFTs_c(1,k)
     $                 ,vy_FFTs_c(1,k)
     $                 ,vz_FFTs_c(1,k)
     $                 ,pr,tp,0,prefix)
          write(prefix,'(A1,I0.2)') 's', 0
          call outpost2(vx_FFTs_s(1,k)
     $                 ,vy_FFTs_s(1,k)
     $                 ,vz_FFTs_s(1,k)
     $                 ,pr,tp,0,prefix)
        endif 

        do j= max(iOutputStart,1),iOutputEnd
         write(prefix,'(A1,I0.2)') 'c', j
          call outpost2(vxp_FFTs_c(1,j,k)
     $                 ,vyp_FFTs_c(1,j,k)
     $                 ,vzp_FFTs_c(1,j,k)
     $                 ,pr,tp,0,prefix)
         write(prefix,'(A1,I0.2)') 's', j
          call outpost2(vxp_FFTs_s(1,j,k)
     $                 ,vyp_FFTs_s(1,j,k)
     $                 ,vzp_FFTs_s(1,j,k)
     $                 ,pr,tp,0,prefix)
        enddo
      end do
      
      if (NIO==0) open (unit=80, file='freqNorm.txt',
     $           status='replace',action='write')
      ntot = nx1*ny1*nz1*nelt
      do k= 1,nFreqs
        j=0
        freqNorms(j,k)=0
        do i=1,ntot
           freqNorms(j,k)=freqNorms(j,k)+
     $      (vx_FFTs_c(i,k)**2.+vy_FFTs_c(i,k)**2.+
     $       vx_FFTs_s(i,k)**2.+vy_FFTs_s(i,k)**2.
     $       )*bm1(i,1,1,1)
           if (ndim==3)
     $     freqNorms(j,k)=freqNorms(j,k)+
     $      (vz_FFTs_c(i,k)**2.+vz_FFTs_c(i,k)**2.
     $       )*bm1(i,1,1,1)
        enddo
        freqNorms(j,k) = glsum(freqNorms(j,k),1)
        freqNorms(j,k) = sqrt(freqNorms(j,k))
      
        do j= 1,npert
         freqNorms(j,k)=0
         do i=1,ntot
          freqNorms(j,k)=freqNorms(j,k)+
     $      (vxp_FFTs_c(i,j,k)**2.+vyp_FFTs_c(i,j,k)**2.+
     $       vxp_FFTs_s(i,j,k)**2.+vyp_FFTs_s(i,j,k)**2.
     $       )*bm1(i,1,1,1)
          if (ndim==3)
     $     freqNorms(j,k)=freqNorms(j,k)+
     $      (vzp_FFTs_c(i,j,k)**2.+vzp_FFTs_c(i,j,k)**2.
     $       )*bm1(i,1,1,1)
          enddo
          freqNorms(j,k) = glsum(freqNorms(j,k),1)
          freqNorms(j,k) = sqrt(freqNorms(j,k))
        enddo
        if (NIO==0) write(80,*) frequencies(k),
     $            freqNorms(iOutputStart:iOutputEnd,k)
      enddo      
      if (NIO==0) close(80) 
      time = timesave
      ifxyo= ifxyo_bak
      return
      end subroutine      
      


      
      subroutine FFt_Multi(jpStart_in, jpEnd_in,Mx,My,Mz)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'FFT_DEF'

      real Mx(lx1,ly1,lz1,lelt)
      real My(lx1,ly1,lz1,lelt)
      real Mz(lx1,ly1,lz1,lelt)
      integer i,j,k,ntot,jpStart,jpEnd
      integer,intent(in),optional  :: jpStart_in,jpEnd_in

      ntot = nx1*ny1*nz1*nelt

      jpStart = 0 
      jpEnd = npert

      if(present(jpStart_in)) jpStart=jpStart_in
      if(present(jpEnd_in  )) jpEnd  =jpEnd_in
      
      if (jpStart==0) then
        do i= 1,ntot
          ! Base flow (jp==0)
            ! Cosine transform
          vx_FFTs_c(i,:)= vx_FFTs_c(i,:)*Mx(i,1,1,1)
          vy_FFTs_c(i,:)= vy_FFTs_c(i,:)*My(i,1,1,1)
             ! Sine transform
          vx_FFTs_s(i,:)= vx_FFTs_s(i,:)*Mx(i,1,1,1)
          vy_FFTs_s(i,:)= vy_FFTs_s(i,:)*My(i,1,1,1)
            ! if 3d
          if ( ndim == 3) then
              vz_FFTs_c(i,:)= vx_FFTs_c(i,:)*Mz(i,1,1,1)
              vz_FFTs_s(i,:)= vz_FFTs_s(i,:)*Mz(i,1,1,1)
          end if
        end do ! i
      endif ! jpstart

      ! Perturbations (jp>0)
      do k= 1,nFreqs
        do i=1,ntot
            do j=max(1,jpStart),jpEnd
                ! Cosine transform
                vxp_FFTs_c(i,j,k)= vxp_FFTs_c(i,j,k)*Mx(i,1,1,1)
                vyp_FFTs_c(i,j,k)= vyp_FFTs_c(i,j,k)*Mx(i,1,1,1)
                 ! Sine transform
                vxp_FFTs_s(i,j,k)= vxp_FFTs_s(i,j,k)*My(i,1,1,1)
                vyp_FFTs_s(i,j,k)= vyp_FFTs_s(i,j,k)*My(i,1,1,1)
                ! if 3d
                if ( ndim == 3) then
                    vzp_FFTs_c(i,j,k)= vzp_FFTs_c(i,j,k)*Mz(i,1,1,1)
                    vzp_FFTs_s(i,j,k)= vzp_FFTs_s(i,j,k)*Mz(i,1,1,1)
                end if
            enddo ! k
        enddo
      enddo
      end subroutine
      
      subroutine ApplyFft(fileListtxt,freqsTxt,deltaT)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer iipert,Reason,ifile
      integer lastRead_iFile
      integer i,j,k,imod
      integer ntot
      real deltaT
      real bm1bak(lx1,ly1,lz1,lelt)
      real tmpV(lx1*ly1*lz1*lelt),tmpPt(lx2*ly2*lz2*lelt)
      real tmpTp(lx1*ly1*lz1*lelt*ldimt)
      character(len=100) fileListtxt,currFilename,freqsTxt      
      logical allFilesRead
      allFilesRead =.false.
      lastRead_iFile = -1
      ntot = nx1*ny1*nz1*nelt
      
      
      
      !backup integration weigths (they are lost during file loading
      bm1bak=bm1
      
      call InitFFts(freqsTxt)

c     Cicle files and get Ft 
      open(unit=80, file=fileListtxt,status='old',action='read')      
      ifile=0  ! Filter position
      do while (.true.)      
        ifile=ifile+1
c        Read one file into filter memory
        if (.not. allFilesRead) then
         read(80,'(A)',IOSTAT=Reason) currFilename
         if (NIO==0) write(*,*) "FFT : File to load ",
     $                          currFilename,'--',Reason
         if (Reason > 0)  THEN
           write(*,*) "FFT : Error reading file lsit "
         else if (Reason < 0) THEN
           if (NIO==0) 
     $       write(*,*) "FFT : Closing file ", 80, Reason
             close(80)           
             exit
         end if ! Reason
        endif !allFilesRead
        call load_fld(currFilename)
c        call Load_FLD_Time_To_wTime(currFilename,
c     $      vxp,vyp,vzp,prp,tp)
        call UpdateFFts(deltaT) 
      enddo 
      !restore integration weigths (they are lost during file loading
      
      bm1=bm1bak
           
      call OutPostFFts()
      
      return
      end subroutine

      subroutine FFTs_FFTConvMinus(sign)
        logical sign
        if (sign)  then
          FourSignCon = -1
        else
          FourSignCon = +1
        endif
      endsubroutine

      subroutine FFts_SetNormTime(normTime_in)
        real normTime_in
        normalizatinTime = normTime_in
      endsubroutine

      subroutine FFts_SetOffSet(offSettime_in)
        real offsetTime_in
        offsetTime =offSettime_in
        return
      endsubroutine

      subroutine FFts_SetTrapzIntRule(trapzIntRule_in)
        logical trapzIntRule_in
        trapzIntRule =trapzIntRule_in
        return
      endsubroutine

      end module FFT
