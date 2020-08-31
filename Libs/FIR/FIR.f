      module FIR
      
      contains 
      
      subroutine ApplyFIR(fileListtxt,fileCoefstxt,deltat,extTol)
      use FFT
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'FIR_DEF'  
      real, intent(in),optional :: extTol
      integer iipert,Reason,ifile
      integer lastRead_iFile
      integer i,j,k,imod
      integer ntot
      real deltat,filerTlag
      real tmpV(lx1*ly1*lz1*lelt),tmpPt(lx2*ly2*lz2*lelt)
      real tmpTp(lx1*ly1*lz1*lelt*ldimt)
      real norm,absTol, relTol
      character(len=100) fileListtxt,currFilename,fileCoefstxt      
      character(len=100) freqList
      logical allFilesRead
      logical doFFT
      logical firstSave
      logical firstFile
      real,external :: glsum

      firstSave=.true. 
      doFFT = .true.
      allFilesRead =.false.
      lastRead_iFile = -1
      ntot = nx1*ny1*nz1*nelt
      absTol = 0
      if (present(extTol)) then
            relTol = extTol
      else
            relTol = 1e-6
      endif
      freqList = 'freqList.input' 
      if (doFFT) call InitFFTs(freqList)
      
c     Initialize Filter Memory---------------      
      if (NIO==0) write(*,*) "FIR : INITIALIZING ..."
      filterMem(:,:,:,:,:,:)=0
      firCoefs(:)=0
            
c     Read Filter Coefs----------------------      
      i=1
      if (NIO==0) write(*,*) "FIR : Reading Filter Coefs ..."
      open (unit=80, file=fileCoefstxt,status='old',action='read')      
      do while (.true.)
         read(80,*,IOSTAT=Reason) firCoefs(i)
         if (Reason > 0)  THEN
           write(*,*) "FIR : Error reading FIR COEFS ", i  
         else if (Reason < 0) THEN
            if (NIO==0) write(*,*) "FIR : Closing file ", 80, Reason
            close(80)
            exit
         else
          if (NIO==0) write(*,*) "FIR : Coef ", i , "=" , firCoefs(i) 
          i=i+1
         end if ! Reason
      enddo
      filterOrder = i-1
      filerTlag = deltaT*(filterOrder/2.)
      if (NIO==0) write(*,*) "FIR : Filter Order ", filterOrder , 
     $  ', time lat ', filerTlag

c     Cicle files and apply FIR filter
      open(unit=80, file=fileListtxt,status='old',action='read')      
      ifile=0  ! Filter position
      do while (.true.)      
        ifile=ifile+1
        if(NIO==0) write(*,*) 'IFILE 1 ' ,  ifile
        
c        Read one file into filter memory
        if (.not. allFilesRead) then
         read(80,'(A)',IOSTAT=Reason) currFilename
         if (NIO==0) write(*,*) "FIR : File to load ",
     $                          currFilename,'--',Reason
         if (Reason > 0)  THEN
           write(*,*) "FIR : Error reading file lsit "
         else if (Reason < 0) THEN
           if (NIO==0) 
     $       write(*,*) "FIR : Closing file ", 80, Reason
             allFilesRead=.true.
             lastRead_iFile =ifile-1
             close(80)           
         end if ! Reason
        endif !allFilesRead
        
        ! If all files have been read, pad with zeros, otherwize read data from file
        imod=modulo(ifile-1,filterOrder)+1
        if(NIO==0) write(*,*) 'IFILE 2 ' ,  ifile
        if (allFilesRead) then
          if (NIO==0) 
     $     write(*,*) 'Setting to zero ' , imod, ifile,filterOrder
          filterMem(:,:,:,:,:,imod) = 0
          time = time + deltaT
        else
          if (NIO==0) 
     $     write(*,*) 'Reading ' , currFilename, imod, ifile,filterOrder
          if (ldim==2) then                  
            call Load_FLD_Time_To_wTime(currFilename,
     $        filterMem(:,:,:,:,1,imod),
     $        filterMem(:,:,:,:,2,imod),
     $        tmpV,tmpPt,tmpTp)
          else
            call Load_FLD_Time_To_wTime(currFilename,
     $        filterMem(:,:,:,:,1,imod),
     $        filterMem(:,:,:,:,2,imod),
     $        filterMem(:,:,:,:,3,imod),
     $        tmpPt,tmpTp)
          endif !ndim
          time  = time - filerTlag ! discount FIR lag
        endif !allFilesRead
        if(NIO==0) write(*,*) 'IFILE 3 ' ,  ifile

c      Compute filter and dump it
        Vx(:,:,:,:)=0
        Vy(:,:,:,:)=0
        if (ndim>2) Vz(:,:,:,:)=0
        do j=1,filterOrder
         imod=modulo(ifile-1-(j-1),filterOrder)+1
         Vx(:,:,:,:)=Vx(:,:,:,:)+filterMem(:,:,:,:,1,imod)*firCoefs(j)
         Vy(:,:,:,:)=Vy(:,:,:,:)+filterMem(:,:,:,:,2,imod)*firCoefs(j)
         if (ndim>2) then
         Vz(:,:,:,:)=Vz(:,:,:,:)+filterMem(:,:,:,:,3,imod)*firCoefs(j)
         endif
        enddo

c       Compute filtered data FFT
        if (doFFT) call UpdateFFts(deltat)

        if (NIO==0) write(*,*) 'FIR : Saving Filtered data'
        ifvo = .true.
        ifpo = .false. 
        ifto = .false. 
        ifpso =.false.  
        
        if (ndim>2) then
         norm  = sum( (Vx**2+Vy**2      )*bm1)
        else
         norm  = sum( (Vx**2+Vy**2+Vz**2)*bm1)
        endif
        
        norm = sqrt(glsum(norm,1))
        absTol=max(absTol,norm*relTol)
        
        if (NIO==0) write(*,*) 'Curr filtered Norm : ',
     $   time , norm , firstSave , absTol
        if (norm > absTol) then
              ifxyo = firstSave
              call outpost(Vx,Vy,Vz,tmpPt,tmpTp,'fir')
              firstSave=.false.
        endif
        
        if (NIO==0) write(*,*) 'Loging Norms '
        if(NIO==0) then
        if (ifile == 1) then
          open (unit=98, file='runNormFiltered.txt',
     $             status='replace',action='write')
          open (unit=99, file='pertNorm.txt',
     $             status='replace',action='write')
        else
          open (unit=98, file='runNormFiltered.txt',
     $        position="append" ,status='old',action='write')
          open (unit=99, file='pertNorm.txt',
     $        position="append" ,status='old',action='write')
        endif
        if (norm > absTol) then
              write(98,*) time , norm , 1
              write(99,*) time , norm 
        else
              write(98,*) time , norm , 0
        endif
              
        close(98)
        endif
        
        if (lastRead_iFile+filterOrder==ifile .and. allFilesRead) then
           IF(NIO==0) write(*,*) 'FIR : Finishing ApplyFIR subroutine.'
     $                          , doFFT
           if (doFFT) then
             IF(NIO==0) write(*,*) 'FIR : Outpostting FFTs ...'
             call OutPostFFts(0,0)
             IF(NIO==0) write(*,*) 'FIR : Outpostting FFTs done'
           endif
           return
        endif
      enddo
      return
      end subroutine
      
      end module FIR
