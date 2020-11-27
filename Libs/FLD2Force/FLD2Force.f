

c========================================================
c     FLD2Force Module 
c
c     Loads files from disk and uses them as external force.
c     For now assumes iteration between direct and adjoint solver,
c     and thus there is a time inversion between them.
c     Forces for times in between two files are obtained by a linear interpolation
c           (higher interpolation order needs to be implemented at some point)
c
c     Subroutines
c           FLD2ForceInit              : Initialize the module.
c           FLD2Force_AreAllFilesRead  : Check if all files in list have been 
c                  read (used in exit conditions).
c           FLD2Force_AreAllFilesReadJ :Check if all files for the j-th perturbation
c                  field in list have been read (used in exit conditions)
c           FLD2ForceUpdate            : Update forces
c           FLD2Force_GetF             : Get force component (previously computed by FLD2ForceUpdate
c
c========================================================
    
      Module FLD2FORCE
      integer, parameter :: maxOrder = 6
      integer,parameter ::  fileUnit = 40
      integer ::  lastIOSTEP(maxOrder)
      logical :: allFilesRead , forceActive
      real ::  fldTimes(maxOrder)
      real ::  fldRefTime
      
      contains
      
c========================================================
c     Initialize module 
c
c     Open list of files and setup variables
c========================================================
            
      subroutine FLD2ForceInit(currFilename)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'FLD2Force_DEF'

      character(len=100) currFilename
      logical,save :: FirstCall = .true. 

      integer i,j,k, Reason
            
      if ( FirstCall ) then ! Check Previously Initialization
        !Initialization Sequence
        if (NIO==0) write(*,*) "Initializing FLD2Force with file",
     $                 currFilename
        FirstCall    = .false. 
        allFilesRead =.false.
        forceActive = .true.
        fldRefTime=-1e99
        ! initialize fields
        fldTimes(:)=-1e99
        fldU(:,:,:,:,:,:) = 0
        fldT(:,:,:,:,:,:) = 0
        
        !open files for reading
        currFilename=trim(currFilename)
        open (unit=fileUnit, file=currFilename,status='old',
     $                                  action='read')
      else ! Debug Message
        if (NIO==0) write(*,*) "FLD2Force reviously Initialized...", 
     $   "Skipping"
      endif !First Call
      return
      end subroutine
     
     
     
      
c========================================================
c     FLD2Force_AreAllFilesRead 
c
c     Check if all files in list have been read. 
c     ( fileUnit <0)
c========================================================

      function FLD2Force_StillActive()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'FLD2Force_DEF'
      logical FLD2Force_StillActive
      FLD2Force_StillActive=forceActive
      return
      end function
      
      
c========================================================
c     FLD2ForceUpdate 
c
c     Check if  the current simulation time is in between the 
c     the time of the loaded files. 
c========================================================
      subroutine FLD2ForceUpdate()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'FLD2Force_DEF'      
      
      integer Reason

      character(len=100) currFilename
      integer ntot,ntot2
      logical, save :: firstCall = .true.
      
      ntot = lx1*ly1*lz1*nelt
      ntot2= lx2*ly2*lz2*nelt

      if ( .not. allFilesRead  ) then 
        !read data if necessary. For interpolation to be in the middle of the domain, read until file maxOrder/2.
        do while (time>=fldTimes(maxOrder/2+1) .or. 
     $            fldTimes(1)<-.5e99) ! if simulation time is larger than the last loaded file.
         read(fileUnit,'(A)',IOSTAT=Reason) currFilename
         if (NIO==0) ! debug messages
     $      write(*,*) "DEBUG : fileUnit to load ", fileUnit,
     $      currFilename , Reason , 
     $      time,fldRefTime,fldTimes
         if (Reason > 0)  THEN
          write(*,*) "DEBUG : Error reading fileUnit ", fileUnit
          call exitt
         else if (Reason < 0) THEN
          if (NIO==0) write(*,*) "DEBUG : Closing file ", 
     $          fileUnit, Reason,currFilename
          close(fileUnit)
          allFilesRead =.true.
          return
         end if ! Reason
          
         if (NIO==0) ! debug messages
     $      write(*,*) "DEBUG : Cycling fields ",maxorder,ntot
          
         !Cicle data and files 
         do k=1,maxorder-1
           fldTimes(k) = fldTimes(k+1)
           fldU(:,:,:,:,:,k) = fldU(:,:,:,:,:,k+1)
           fldT(:,:,:,:,:,k) = fldT(:,:,:,:,:,k+1)
           fldP(:,:,:,:,k)   = fldP(:,:,:,:,k+1) 
         enddo
         
         if (NIO==0) write(*,*) "DEBUG : Loading file" ! debug messages
         tmpTime = time   ! backput current simulation time
        ! Load new file into the end of the arrays
         call Load_FLD_Time_To_wTime(currFilename,
     $    fldU(1,1,1,1,1,maxorder),fldU(1,1,1,1,2,maxorder),
     $    fldU(1,1,1,1,3,maxorder),fldP (1,1,1,1,maxorder),
     $    fldT(1,1,1,1,1,maxorder) )

         if (NIO==0) write(*,*) 'Done'
         fldTimes(maxorder) = -time ! Time reversal, using first time-stamp as reference
     
       ! Save first timestamp as reference time
         if (firstCall) then 
          fldRefTime = ftime 
          time = -Time      ! set initial time as the time of the first fld file read
          firstCall=.false.
         else
          time = tmpTime    ! restore simulation time
         end if

       enddo ! while
      endif
     
      forceActive = fldTimes(maxorder)+dt/2>time
      return
      end subroutine
      
      

c========================================================
c     FLD2Force_GetF 
c
c     Linearly interpolate between the two files read and 
c     return force values.
c========================================================
      subroutine FLD2Force_GetF(ix,iy,iz,eg,
     $                      fffx,fffy,fffz)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'FLD2Force_DEF'      
      
      integer,intent(in) :: ix,iy,iz,eg
      integer e,ntot    
      real,save :: prevTime = -1e99
      real :: xx(maxorder),yy(maxorder,3), xx0
      real,intent(out) ::  fffx,fffy,fffz

       if (prevTime /= time .and. NIO==0 )then  
          write(*,*) 'FLD2Force_GetF Computing Forces :', forceActive
          write(*,*) 'FLD2Force_GetF Times :', fldTimes(1:maxorder)
          write(*,*) 'Time  :', time
          prevTime = time
        endif
      
      if (.not. forceActive) then 
          fffx = 0.0
          fffy = 0.0
          fffz = 0.0
          return ! no extra data to be used
      endif

      e = gllel(eg)
      ! interpolate to get current forcing
      ! switch conditions for different interpolation schemes.
       xx0 = time
       xx = fldTimes(1:maxorder)
       yy(:,1) = fldU(ix,iy,iz,e,1,1:maxorder)
       yy(:,2) = fldU(ix,iy,iz,e,2,1:maxorder)
       yy(:,3) = fldU(ix,iy,iz,e,3,1:maxorder)

      if (maxorder==2) then  ! linear interpolation
        call LinInterp(xx,yy(:,1),xx0,fffx)
        call LinInterp(xx,yy(:,2),xx0,fffy)
        if (ndim>2) then
        call LinInterp(xx,yy(:,3),xx0,fffz)
        endif
      else if (maxorder==4) then  ! cubic interpolation
        call CubicInterp(xx,yy(:,1),xx0,fffx)
        call CubicInterp(xx,yy(:,2),xx0,fffy)
        if (ndim>2) then
        call CubicInterp(xx,yy(:,3),xx0,fffz)
        endif
      else if (maxorder==6) then  ! cubic interpolation
        call FifthOrderInterp(xx,yy(:,1),xx0,fffx)
        call FifthOrderInterp(xx,yy(:,2),xx0,fffy)
        if (ndim>2) then
        call FifthOrderInterp(xx,yy(:,3),xx0,fffz)
        endif
      else
        write(*,*) 'ABORT: Invalid Polinomial Order (',
     $               maxOrder,')'
        call exitt 
      endif
      
      return
      end subroutine
      
c========================================================
c     Interpolation Routines
c
c     Linear, cubic and 5-th order interpolations (C0,C1,C2)
c========================================================
      subroutine LinInterp(x,y,xout,yout)
      real :: x(2),y(2),xout,yout
      real :: deltaT1,deltaT2
      deltaT1 = x(2)-xout
      deltaT2 = xout-x(1) 
      yout = ( y(1)*deltaT1 + y(2)*deltaT2 )
     $                / (deltaT1+deltaT2)
      return
      end subroutine      
      
      subroutine CubicInterp(x,y,xout,yout)
      real :: x(4),y(4),xout,yout ,h 
      real :: y1,y2,x1,x2,dx1,dx2,pl0,pl1,pr0,pr1
      
      h  =  x(3)-x(2);
      
      ! reference for coefficients 
      ! https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference
      ! Wikipedia, I known....
      if (xout>=x(2) .and. xout<=x(3) ) then
        x1 = x(2)
        x2 = x(3)
        
        y1 = y(2)
        y2 = y(3)

        y1p =(y(3)-y(1))/(2*h); ! centered finite differences
        y2p =(y(4)-y(2))/(2*h); ! centered finite differences
      elseif ( xout<=x(2) ) then
        x1 = x(1)
        x2 = x(2)
        
        y1 = y(1)
        y2 = y(2)

        y1p =(y(2)-y(1))/(h); ! forward finite differences
        y2p =(y(3)-y(1))/(2*h); ! centered finite differences
      elseif ( xout>=x(3)  ) then
        x1 = x(3)
        x2 = x(4)
        
        y1 = y(3)
        y2 = y(4)

        y1p =(y(4)-y(2))/(2*h); ! centered finite differences
        y2p =(y(4)-y(3))/(h); ! backward finite differences
      endif

      if ( xout>x(4) .or. xout<x(1) ) then
       write(*,*) 'Warning!!! Interpolation out of limits!!! ',
     $             xout,' out of ' ,x(:)
      endif
      
      ! polynomials based on dofs on valus and
      ! derivatives around the center points
      dx1 = (xout-x1)/h;
      pl0  = 1     - 3*dx1**2 + 2*dx1**3 ;  ! value of 1 at the left
      pl1  =    dx1 - 2*dx1**2 +   dx1**3 ; ! derivative of value 1 at the left
      dx2 = (x2-xout)/h;
      pr0  =   1     - 3*dx2**2 + 2*dx2**3 ; ! value of 1 at the right
      pr1  =-(  dx2  - 2*dx2**2 +   dx2**3); !derivative of value 1 at the left

      yout = y1*pl0 + y2*pr0 + y1p*pl1*h + y2p*pr1*h; !builds up polynomial and get final value

      return
      end subroutine      
       
      subroutine FifthOrderInterp(x,y,xout,yout)
      implicit none
      real :: x(6),y(6),xout,yout ,h 
      real :: dx1,dx2,pl0,pl1,pl2,pr0,pr1,pr2
      real :: xa,xb
      real :: ya,yap,yapp
      real :: yb,ybp,ybpp
      
      
      h  =  x(4)-x(3);
      ! reference for coefficients 
      ! https://en.wikipedia.org/wiki/Finite_difference_coefficient#Central_finite_difference
      ! Wikipedia, I known....
      if (xout>=x(3) .and. xout<=x(4) ) then
       !interpolate between x2 and x3 using numerical derivatives on each limit
       xa = x(3)
       xb = x(4)
       
       ya = y(3)
       yb = y(4)

       yap  = (-y(5) +8*y(4)         -8*y(2)+y(1) )/(12*h   ); ! 4-th order centered finite differences
       yapp = (-y(5)+16*y(4)-30*y(3)+16*y(2)-y(1) )/(12*h**2); ! 4-th order centered finite differences
      
       ybp  = (-y(6) +8*y(5)         -8*y(3)+y(2) )/(12*h   );! 4-th order centered finite differences
       ybpp = (-y(6)+16*y(5)-30*y(4)+16*y(3)-y(2) )/(12*h**2); ! 4-th order centered finite differences
      
      elseif (xout>=x(2) .and. xout<=x(3) ) then
       !interpolate between x2 and x3 using numerical derivatives on each limit
       xa = x(2)
       xb = x(3)
       
       ya = y(2)
       yb = y(3)

       yap  = ( y(3)       -y(1) )/(2*h    ); ! 2-th order centered finite differences
       yapp = ( y(3)-2*y(2)+y(1) )/(   h**2); ! 2-th order centered finite differences

       ybp  = (-y(5) +8*y(4)         -8*y(2)+y(1) )/(12*h   ); ! 4-th order centered finite differences
       ybpp = (-y(5)+16*y(4)-30*y(3)+16*y(2)-y(1) )/(12*h**2); ! 4-th order centered finite differences
      elseif ( xout<=x(2) ) then
       !interpolate between x2 and x3 using numerical derivatives on each limit
       xa = x(1)
       xb = x(2)
       
       ya = y(1)
       yb = y(2)

       yap  = (y(2)-y(1) )/(h); ! 1-th order centered finite differences
       yapp =  0; 

       ybp  = ( y(3)       -y(1) )/(2*h   ); ! 2-th order centered finite differences
       ybpp = ( y(3)-2*y(2)+y(1) )/(  h**2); ! 2-th order centered finite differences

      elseif (xout>=x(4) .and. xout<=x(5) ) then
       !interpolate between x2 and x3 using numerical derivatives on each limit
       xa = x(4)
       xb = x(5)
       
       ya = y(4)
       yb = y(5)

       yap  = (-y(6) +8*y(5)         -8*y(3)+y(2) )/(12*h   );! 4-th order centered finite differences
       yapp = (-y(6)+16*y(5)-30*y(4)+16*y(3)-y(2) )/(12*h**2); ! 4-th order centered finite differences

       ybp  = ( y(6)       -y(4) )/(2*h   ); ! 2-th order centered finite differences
       ybpp = ( y(6)-2*y(5)+y(4) )/(  h**2); ! 4-th order centered finite differences
      elseif (xout>=x(5)  ) then
       !interpolate between x2 and x3 using numerical derivatives on each limit
       xa = x(5)
       xb = x(6)
       
       ya = y(5)
       yb = y(6)

       yap  = ( y(6)       -y(4) )/(2*h   ); ! 2-th order centered finite differences
       yapp = ( y(6)-2*y(5)+y(4) )/(  h**2); ! 2-th order centered finite differences
     
       ybp  = ( y(6)-  y(5)      )/(  h   ); ! 1-th order centered finite differences
       ybpp =  0; 

      endif

      if ( xout>x(6) .or. xout<x(1) ) then
       write(*,*) 'Warning!!! Interpolation out of limits!!! ',
     $             xout,' out of ' ,x(:)
      endif
      ! polynomials based on dofs on values and
      ! derivatives around the center points

      dx1 = (xout-xa)/h;
      pl0  = 1         -10  *dx1**3+15  *dx1**4-6  *dx1**5 ;
      pl1  =     dx1   - 6  *dx1**3+ 8  *dx1**4-3  *dx1**5 ;
      pl2  =  .5*dx1**2- 1.5*dx1**3+ 1.5*dx1**4- .5*dx1**5 ;

      dx2 = (xb-xout)/h;
      pr0  = 1         -10  *dx2**3+15  *dx2**4-6  *dx2**5 ;
      pr1  =     dx2   - 6  *dx2**3+ 8  *dx2**4-3  *dx2**5 ;
      pr2  =  .5*dx2**2- 1.5*dx2**3+ 1.5*dx2**4- .5*dx2**5 ;

      yout = ya*pl0 + yap*pl1*h + yapp*pl2*h**2 +
     $       yb*pr0 - ybp*pr1*h + ybpp*pr2*h**2 ;  

      return
      end subroutine      
      
      end module
