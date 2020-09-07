c========================================================
c     Control Module 
c
c     Contains Routines for using a convolution of sensor
c     readings and a Kernel for control.
c
c     Subroutines
c           ControlInit           : Initialize the module.
c           ControlUpdateReadings: Initialize the module.
c           Control_CalcActuations: Initialize the module.
c           Control_LinInterp: Linear Interpolation routine.
c
c========================================================

      module Control
      implicit none
      integer, parameter :: maxInputSensors   = 10    ! defines max number of sensors
      integer, parameter :: maxActuators      = 10    ! defines max number of actuators
      integer, parameter :: maxControlLawLines= 10000 ! defines max number of points in the control law
      integer, parameter :: SensorHistoryLen  = 10000 ! defines max sensor history size
      integer, parameter :: fileNameLength  = 40      ! max filenime sizes

      

      real,save:: InputSensorHistory(SensorHistoryLen,maxInputSensors) ! array containing sensor histories
      real,save:: timeHistory(SensorHistoryLen)                        ! array containing time values for sensor history

      ! Array contaning control kernel
      real,save:: ControlLaw(maxControlLawLines,maxInputSensors,
     $                maxActuators)
     
      real,save:: control_dt , actuatorForce(maxActuators) ! delta t used in the Kernel / computed acuations

      integer,save:: nControlLawLines, nSensors,nActs  ! number of lines used in the kernel/ number of sensors and actuators used
      
      contains 
      
c========================================================
c     ControlInit
c
c     Reads parameter file, creates output files, performs
c     sensors readings calculations, and outputs readings.
c========================================================

      subroutine ControlInit(controlLawParamFile)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      character(LEN=40) controlLawFile
      character(LEN=40) controlLawParamFile
      integer i , j ,k  ,Reason ,tmp
      logical endOfFile

c     Initi Histories to Zero
      do i=1,SensorHistoryLen
         TimeHistory(i) = -1e99
           do j=1,maxInputSensors
               InputSensorHistory(i,j)=0
           enddo ! j
      enddo ! i
c     Read Control Laws

      if (NIO==0) write(*,*) "AdjEst: Reading Control Law "
      
      
c     Reads control law file, which has
c           "#Sensors" "#actuators" "# kernel points (control lines)" "time discretization"
c           Kernel for actuator 1
c           Kernel for actuator 2
c           ...

      open (unit=98, file=controlLawParamFile,
     $               status='old',action='read')
     
      endOfFile = .false.
      i = 0
      read(98,*,IOSTAT=Reason) nSensors, nActs  ,
     $                          nControlLawLines , control_dt
      do while (.not. endOfFile)
        read(98,*,IOSTAT=Reason) controlLawParamFile
        IF (Reason > 0)  THEN
          if (NIO==0) write(*,*) "Error reading ", controlLawFile
          endOfFile = .true.
        else if (Reason < 0) THEN
          endOfFile = .true.
        else
          i = i+1
          open (unit=99, file=controlLawParamFile,
     $                   status='old',action='read')
          do j=1,nControlLawLines
            read(99,*) ControlLaw(j,1:nSensors,i)
          enddo
          close(99)
        endif
        
        if (NIO==0) then
            write(*,*) 'Control Law ', i , ' (first 5 lines)'
          do j=1,5
            write(*,*) ControlLaw(j,1:nSensors,i)
            enddo
        endif
      end do
      close(98)
      
      call ControlPrintReadingsHist() 
      
      return
      end subroutine
      
      
      
c========================================================
c     ControlUpdateReadings
c
c     Updates the sensor and time histories with current values.
c========================================================
      subroutine ControlUpdateReadings(
     $              currtime,readings)
c      integer  nSensors
      integer i,j
      real currtime, readings(nSensors)
      ! loop time history and sensors
      if (currtime - TimeHistory(2) >control_dt/3) then
        do j=SensorHistoryLen,2,-1
          TimeHistory(j) = TimeHistory(j-1)
            do i=1,nSensors
                 InputSensorHistory(j  ,i) =
     $           InputSensorHistory(j-1,i)
            enddo
          enddo
      endif
      

      TimeHistory(1) = currtime
      do i=1,nSensors
         InputSensorHistory(1,i) = readings(i)
      enddo
      

      return
      end subroutine
      


c========================================================
c     Control_CalcActuations
c
c     Compute actuations from the kernels.
c     Kernel is interpolated on the time history, and 
c     intergrals computed using a trapeizodal rule     
c========================================================
      subroutine Control_CalcActuations
      include 'SIZE'
      include 'TOTAL'
      real interp_hist(maxControlLawLines,maxInputSensors)
      real tau(maxControlLawLines), weight
      integer i,j,k
c      integer,external :: NIO

      ! Compute time lag (tau)
      do i=1,nControlLawLines
          tau(i)= time-(i-1)*control_dt
      enddo

      ! Interpolate kernels 
      do i=1,nSensors
        call Control_LinInterp(
     $   SensorHistoryLen, timeHistory , InputSensorHistory(:,i) ,
     $   maxControlLawLines, tau,interp_hist(:,i))
      enddo


      ! Compute integral using the Trapeizodal rule 
      do k=1,nActs
        actuatorForce(k)=0
         do j=1,nControlLawLines
           do i=1,nSensors
            ! weights for trapeizodal rule.
            if (j==1 .or. j==nControlLawLines ) then
                  weight = 0.5
            else      
                  weight = 1.0
            endif
            
            if(isnan(interp_hist(j,i))) interp_hist(j,i) = 0 ! GAMBETA!!! CORRIGIR INTERPOLACAO
            actuatorForce(k) = actuatorForce(k)
     $          +  interp_hist(j,i)*ControlLaw(j,i,k)*control_dt*weight
         enddo ! i, input
        enddo ! j, time
      enddo ! k control actuator

      
      return
      end subroutine

c========================================================
c      Control_Actuation
c
c     Returns the i-th acuation, previously computed with  Control_CalcActuation.
c     Just an interface function to the modules "actuatorForce" array
c========================================================
      function Control_Actuation(i)
      real Control_Actuation
      integer i
      Control_Actuation = actuatorForce(i)
      return
      end function


c========================================================
c      Control_LinInterp
c
c     Linear interpolation function. Literally found it online.
c========================================================
      subroutine Control_LinInterp ( nd, xd, yd, ni, xi ,yi)

      !*****************************************************************************80
      !! Control_LinInterp evaluates the piecewise linear interpolant.
      !  Discussion:
      !    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
      !    linear function which interpolates the data (XD(I),YD(I)) for I = 1
      !    to ND.
      !  Licensing:
      !    This code is distributed under the GNU LGPL license.
      !  Modified:
      !    22 September 2012
      !  Author:
      !    John Burkardt
      !  Parameters:
      !    Input, integer ( kind = 4 ) ND, the number of data points.
      !    ND must be at least 1.
      !    Input, real ( kind = 8 ) XD(ND), the data points.
      !    Input, real ( kind = 8 ) YD(ND), the data values.
      !    Input, integer ( kind = 4 ) NI, the number of interpolation points.
      !    Input, real ( kind = 8 ) XI(NI), the interpolation points.
      !    Output, real ( kind = 8 ) YI(NI), the interpolated values.
        implicit none

        integer ( kind = 4 ) nd
        integer ( kind = 4 ) ni

        integer ( kind = 4 ) i
        integer ( kind = 4 ) k
        real ( kind = 8 ) t
        real ( kind = 8 ) xd(nd)
        real ( kind = 8 ) yd(nd)
        real ( kind = 8 ) xi(ni)
        real ( kind = 8 ) yi(ni)
        logical flipD, flipI

        flipD = xd(1)>xd(nd)
        flipI = xi(1)>xi(nd)

        if (flipD) then
          xd = xd(nd:1:-1)
          yd = yd(nd:1:-1)
        endif
        if (flipI) then
          xi = xi(ni:1:-1)
          yi = yi(ni:1:-1)
        endif


        yi(1:ni) = 0.0D+00

        if ( nd == 1 ) then
          yi(1:ni) = yd(1)
          return
        end if

        do i = 1, ni
          if ( xi(i) < xd(1) ) then
c            t = ( xi(i) - xd(1) ) / ( xd(2) - xd(1) )
c            yi(i) = ( 1.0D+00 - t ) * yd(1) + t * yd(2)
            yi(i) = 0
          else if ( xd(nd) < xi(i) ) then
c            t = ( xi(i) - xd(nd-1) ) / ( xd(nd) - xd(nd-1) )
c            yi(i) = ( 1.0D+00 - t ) * yd(nd-1) + t * yd(nd)
            yi(i) = 0
          else
            do k = 2, nd
              if ( xd(k-1) <= xi(i) .and. xi(i) <= xd(k) ) then
                t = ( xi(i) - xd(k-1) ) / ( xd(k) - xd(k-1) )
                yi(i) = ( 1.0D+00 - t ) * yd(k-1) + t * yd(k)
                exit
              end if
            end do
          end if
        end do

        if (flipD) then
          xd = xd(nd:1:-1)
          yd = yd(nd:1:-1)
        endif
        if (flipI) then
          xi = xi(ni:1:-1)
          yi = yi(ni:1:-1)
        endif

      return
      end
      
      subroutine ControlPrintReadingsHist()
      include 'SIZE'
      include 'TOTAL'      
      integer i,j 
      if (NIO==0) then
        write(*,*) 'First 10 lines of sensor history ...'
        do i=1,10
          write(*,*) TimeHistory(i),
     $          InputSensorHistory(i,1:nSensors)               
        enddo
      endif
      
      return
      end subroutine


      
      
      subroutine ControlPrintActuation()
      include 'SIZE'
      include 'TOTAL'    
c      integer i,j , nSensors
      if (NIO==0) then
          write(*,*) 'Current Actuations ', time ,
     $          actuatorForce(1:nActs)               
      endif
      
      return
      end subroutine
      
      
      
      subroutine SaveControlHist()
      include 'SIZE'
      include 'TOTAL' 
      logical,save :: firstHistSave  = .true.  
      
      if (NIO==0) then
        if (firstHistSave) then
          open (unit=98, file='Control_actHist.txt',
     $            status='replace',action='write')
          open (unit=99, file='Control_sensHist.txt',
     $            status='replace',action='write')
          firstHistSave = .false.
        else
          open (unit=98, file='Control_actHist.txt',
     $         position="append" ,status='old',action='write')
          open (unit=99, file='Control_sensHist.txt',
     $         position="append" ,status='old',action='write')
        endif
        
        write(98,*) time, actuatorForce(1:nActs) 
        write(99,*) TimeHistory(1),
     $              InputSensorHistory(1,1:nSensors) 
        close(98)
        close(99)
        
        
      endif
      
      return
      end subroutine 


      end module Control
      
      
      
