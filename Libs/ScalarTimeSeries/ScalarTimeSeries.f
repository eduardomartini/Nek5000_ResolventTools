      module ScalarTimeSeries
      implicit none

      integer, parameter :: maxNcolumns = 100
      integer            :: ncolumns
      integer, parameter :: interpOrder = 6
      integer, parameter :: fileUnit = 41
      real               :: t_data(interpOrder)
      real               :: data(interpOrder,maxNcolumns)
      logical            :: flagFileEnded,flagInvertTime,flagSetTime

      contains
c========================================================
c     Initialize module 
c========================================================
      subroutine ScalarTimeSeries_Init(currFilename,ncolumns_in,
     $                                        finvT,fsetT)
        include 'SIZE'
        include 'TOTAL'
        character(len=100) currFilename
        integer i,j,k, Reason,ncolumns_in
        logical finvT,fsetT

        flagInvertTime=finvT
        flagSetTime   = fsetT

        ncolumns=ncolumns_in 

        if (NIO==0) write(*,*) "Initializing ScalarTimeSeries ",
     $                 "with file",currFilename
        flagFileEnded =.false.
        t_data        = -1e99
        data          =  0
        !open files for reading
        currFilename=trim(currFilename)
        open (unit=fileUnit,file=currFilename,
     $      status='old',action='read')
        
        call ScalarTimeSeries_Update(time)
        return
      end subroutine
            
c========================================================
c     Check if time series is over 
c========================================================
      function ScalarTimeSeries_StillActive()
        logical ScalarTimeSeries_StillActive
        ScalarTimeSeries_StillActive= .not. flagFileEnded
      end function
c========================================================
c     Get current Value
c========================================================
      subroutine ScalarTimeSeries_Get(t_curr,range,values)
        real   ,intent(in)  :: t_curr
        integer,intent(in)  :: range(2)
        real   ,intent(out) :: values(1:range(2)-range(1)+1)
        real   ,save        :: t_prev = -1e99
        integer             :: i,j

        if (flagFileEnded) then
          values=0
          return
        endif

        !check if data needs updating
        if (t_curr /= t_prev ) then
          call ScalarTimeSeries_Update(t_curr)
          t_prev = t_curr
        endif
        ! Interpolate data
        do i=0,range(2)-range(1)
          j=range(1)+j
          if (interpOrder == 2) then
              call LinInterp        ( t_data(:),data(:,j),
     $                                       t_curr,values(i+1))
          elseif (interpOrder == 4) then
              call CubicInterp      ( t_data(:),data(:,j),
     $                                       t_curr,values(i+1))
          elseif (interpOrder == 6) then
              call FifthOrderInterp ( t_data(:),data(:,j),
     $                                       t_curr,values(i+1))
          else
            write(*,*) 'Invalid Interpolation Order!!! Exiting'  
            call exitt
          endif
        enddo
      end subroutine

c========================================================
c     Update Values
c========================================================
      subroutine ScalarTimeSeries_Update(t_curr_in)
        include 'SIZE'
        include 'TOTAL'      
        real,intent(in)  :: t_curr_in
        real             :: t_curr
        integer          :: i,j,reason
        t_curr=t_curr_in

        !check if data needs updating
        if (.not. flagFileEnded) then
          ! update conditions, t < middle of the data and file not ended
          ! Also updates if the end of the data <t, for initialization
          do while (t_curr>=t_data(interpOrder/2+1) .or. 
     $              t_data(1)<-1e98 .and. 
     $        .not. flagFileEnded)
             t_data   (1:interpOrder-1)  =t_data   (2:interpOrder)
             data     (1:interpOrder-1,:)=data     (2:interpOrder,:) 


             read(fileUnit,*,IOSTAT=Reason) 
     $                 t_data(interpOrder),data(interpOrder,1:ncolumns)

            if (reason > 0)  THEN
              write(*,*) "ScalarTimeSeries_Update : ",
     $                    " Error reading fileUnit ", fileUnit
              call exitt
            else if (Reason < 0) then
             if (NIO==0) write(*,*) "ScalarTimeSeries_Update : ",
     $              " Closing file ",  fileUnit, Reason
              flagFileEnded =.true.
              close(fileUnit)
            endif
            if (flagInvertTime) 
     $              t_data(interpOrder)=-t_data(interpOrder)
            if (flagSetTime) then
                time = t_data(interpOrder)
                t_curr = time
                flagSetTime=.false.
            endif 
         enddo
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
      real :: y1,y2,x1,x2,dx1,dx2,pl0,pl1,pr0,pr1,y1p,y2p
      
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