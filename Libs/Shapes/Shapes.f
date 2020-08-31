      module Shapes
      implicit none
      integer, parameter :: shapesNParams  = 9
      integer, parameter :: nmaxShapes     = 230
      integer, parameter :: shapeTypeLen   = 3      
      integer nShapes
      real shapeParams(nmaxShapes,shapesNParams)
      character(len=shapeTypeLen), 
     $         dimension(nmaxShapes):: shapeType


      contains
c      character shapeType(shapesNParams,3)
      
c========================================================
c     Shape module 
c           (use of fortran modules still to be implemented)
c     Reads a shape file and define domain shapes.
c     For gaussian sensors shapes.txt lines read
c     GGG xc yc zc sigx sigy sigz sx sy sz
c     with corresponding shape on the X,Y,Z direction being
c           shape_x(x,y,z)=A sx exp( -(x-xc)^2/(2*sigx)-(y-yc)^2/(2*sigy)-(z-zc)^2/(2*sigz)).
c           shape_y(x,y,z)=A sy exp( -(x-xc)^2/(2*sigx)-(y-yc)^2/(2*sigy)-(z-zc)^2/(2*sigz)).
c           shape_z(x,y,z)=A sz exp( -(x-xc)^2/(2*sigx)-(y-yc)^2/(2*sigy)-(z-zc)^2/(2*sigz)).
c           with a being A normalization term.
c========================================================



c========================================================
c     ReadShapeFile
c
c     Reads sensor parameters from file
c     shapeFile : file containing sensor information on
c       the form of a matrix of nxnMaxParams.
c     fileNameLength : lenght of shapeFile
c     list : array of n,7 for returning data
c     maxListSize : size of n (lines of list)
c     nLines : number of lines read from file
c========================================================
      subroutine ReadShapesFile(shapeFile)
C       include 'Shapes_DEF'
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      integer fileNameLength
      character(LEN=100) shapeFile
      integer lineCount ,  Reason , nLines , i
      logical endOfFile
c    ---------------- Read Input Sensors List   ------------------------------
      if(NIO==0) write(*,*) "Shapes: Reading File :" ,
     $ shapeFile 

      lineCount = 0
      endOfFile = .false.
      open (unit=98, file=shapeFile,
     $                     status='old',action='read')
      do while (.not. endOfFile)
            read(98,*,IOSTAT=Reason) shapeType(lineCount+1),
     $                               shapeParams(lineCount+1,:)

            IF (Reason > 0)  THEN
              write(*,*) "Error reading ", shapeFile 
              endOfFile = .true.
            else if (Reason < 0) THEN
              endOfFile = .true.
            ELSE
             lineCount = lineCount+1
            END IF
      end do
      close(98)
      write(*,*) "Left Loop" , lineCount
      nShapes = lineCount

      if(NIO==0) then
            do i=1,lineCount
            write(*,'(I0.2,A,A3,A,*(f8.2))') i, " ",
     $                      shapeType(i)," : " , shapeParams(i,:)
            enddo
      end if
      
      return
      end subroutine

      
c========================================================
c     exportShapes
c
c     Save shapes to disk as fld files. For visualisation/ debug
c========================================================
      subroutine exportShapes()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C       include 'Shapes_DEF'
      integer ntot , i , ishape
      real sx(nx1,ny1,nz1,nelt)
      real sy(nx1,ny1,nz1,nelt)
      real sz(nx1,ny1,nz1,nelt)
      real sp(nx1,ny1,nz1,nelt)

      
      real xx,yy,zz
      ntot = nx1*ny1*nz1*nelt
      
c      write(*,*) nx1,ny1,nz1,nelt,nShapes,ntot

      do ishape = 1,nShapes
        do i = 1,ntot
          xx = xm1(i,1,1,1)
          yy = ym1(i,1,1,1)
          zz = zm1(i,1,1,1)
          call getShape(ishape,xx,yy,zz,sx(i,1,1,1),
     $                        sy(i,1,1,1),sz(i,1,1,1))
          sp(i,1,1,1) = 0
        enddo
        call outpost2(sx,sy,sz,sp,sp,0,'sha')
c        write(*,*) "done"
      enddo
      return
      end subroutine
      
      
c========================================================
c     getShape
c
c     Get shape value (sx,sy,sz) at the position xx,yy,zz
c========================================================
      subroutine getShape(ishape,xx,yy,zz,sx,sy,sz)
      include 'SIZE'
      include 'TOTAL'
c      include 'NEKUSE'
C       include 'Shapes_DEF'
      real xyz(3),xx,yy,zz,sx,sy,sz
      real tmp , x_tmp , L ,sensShape
      character(len=3) currShapeType
      integer i ,ishape
      
      currShapeType = shapeType(ishape)

      sensShape = 1.0
      PI=4.D0*DATAN(1.D0)
      
      !Fourier - Cheb - Fourrier, with dirletrich - dirletrich - periodic
      if (currShapeType(1:1) == 'F' .and. 
     $    currShapeType(2:2) == 'C' .and.
     $    currShapeType(3:3) == 'F'      ) then

      ! Sensor parameter : xc , xL , xn , yL , yN , zk
          if (abs(xx-shapeParams(ishape,1))<shapeParams(ishape,2)) then
            x_tmp = (xx-shapeParams(ishape,1))/shapeParams(ishape,2) + 1
            if (shapeParams(ishape,3) > 0) then
              sensShape = sensShape*
     $              sin(0.5*pi*shapeParams(ishape,3)*x_tmp)
            else
              sensShape = sensShape*
     $              cos(0.5*pi*shapeParams(ishape,3)*x_tmp)
            end if
                          
          else 
            sensShape=0          
          end if
          L = shapeParams(ishape,4)
c          x_tmp = (yy-L)/(yy+L)
c          sensShape =sensShape*cos(shapeParams(ishape,5)*acos(x_tmp))
          if (yy<L) then
            if (shapeParams(ishape,5) > 0) then
             sensShape = sensShape*sin(pi*shapeParams(ishape,5)*yy/L)
           else
             sensShape = sensShape*cos(pi*shapeParams(ishape,5)*yy/L)
           end if
          else  
             sensShape = 0
          end if
          
          if (shapeParams(ishape,6) >= 0) then
             sensShape = sensShape*cos(2.*pi*zz*shapeParams(ishape,6))
          else
             sensShape = sensShape*sin(2.*pi*zz*shapeParams(ishape,6))          
          end if
          sx = sensShape*shapeParams(ishape,7)
          sy = sensShape*shapeParams(ishape,8)
          sz = sensShape*shapeParams(ishape,9)
          return
      end if
      
		
      xyz(1) = xx
      xyz(2) = yy
      xyz(3) = zz
      do i=1,ndim
        if (currShapeType(i:i) == 'g') then
            sensShape = sensShape*
     $      exp( -(xyz(i)-shapeParams(ishape,i))**2. / 
     $              (2.*shapeParams(ishape,i+3) **2.))
        elseif (currShapeType(i:i) == 'G') then
            sensShape = sensShape*
     $      exp( -(xyz(i)-shapeParams(ishape,i))**2./
     $              (2.*shapeParams(ishape,i+3) **2.)) /
     $      (sqrt(2.*pi)*shapeParams(ishape,i+3))
        elseif (currShapeType(i:i) == 't') then
            if (shapeParams(ishape,i) >= 0) then
             sensShape = sensShape*
     $       cos(2.*pi*xyz(i)*shapeParams(ishape,i))
            else
             sensShape = sensShape*
     $       sin(2.*pi*xyz(i)*shapeParams(ishape,i))
            end if
        elseif (currShapeType(i:i) == 's') then
            sensShape = sensShape*xyz(i)*
     $    (1-tanh(  (xyz(i)-shapeParams(ishape,i  ))/
     $                      shapeParams(ishape,i+3)  ))/2.0
        elseif (currShapeType(i:i) == 'r') then
            sensShape = sensShape*0.5*(1-
     $       tanh(0.1*(abs(xyz(i)-shapeParams(ishape,i  )) -
     $                            shapeParams(ishape,i+3) )))
        end if

      end do
      sx = sensShape*shapeParams(ishape,7)
      sy = sensShape*shapeParams(ishape,8)
      sz = sensShape*shapeParams(ishape,9)
        

      return
      end subroutine


c========================================================
c     ProjectOnShapes
c
c     Project flow perturbation (inner produc) onto the shapes 
c========================================================
      function ProjectOnShapes(jjp)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C       include 'Shapes_DEF'
      real  :: ProjectOnShapes(nShapes) 
      real xx,yy , sx,sy,sz
      integer ishape,i,ntot,jjp
      real, external:: glsum
      ntot = nx1*ny1*nz1*nelt
      ProjectOnShapes(:) = 0

      if (jjp==0) then
        !routine for reading non-linear run sensors

        do i = 1,ntot
          do ishape = 1,nShapes
            call getShape(ishape,xm1(i,1,1,1),ym1(i,1,1,1),
     $                           zm1(i,1,1,1),sx,sy,sz)
            ProjectOnShapes(ishape) = ProjectOnShapes(ishape)  +
     $         (Vx(i,1,1,1)*sx + Vy(i,1,1,1)*sy )*bm1(i,1,1,1)
            if (ndim>2) ProjectOnShapes(ishape) = 
     $         ProjectOnShapes(ishape)+Vz(i,1,1,1)*sz*bm1(i,1,1,1)
          enddo
        enddo
      else
        !routine for reading linearized run sensors
        do i = 1,ntot
          do ishape = 1,nShapes
            call getShape(ishape,xm1(i,1,1,1),ym1(i,1,1,1),
     $                         zm1(i,1,1,1),sx,sy,sz)
            ProjectOnShapes(ishape) = ProjectOnShapes(ishape)  +
     $       (Vxp(i,jjp)*sx + Vyp(i,jjp)*sy )*bm1(i,1,1,1)
            if (ndim>2) ProjectOnShapes(ishape) = 
     $        ProjectOnShapes(ishape)+Vzp(i,jjp)*sz*bm1(i,1,1,1)
          enddo
        enddo
      endif

      do ishape = 1,nShapes
          ProjectOnShapes(ishape) = glsum(ProjectOnShapes(ishape),1)
      enddo


      return
      end function
      
c========================================================
c     SaveProjUp2Shapes
c
c     Project main flow (inner produc) onto the 
c     shapes and save to file
c========================================================
      subroutine SaveProjUp2Shapes(jp1,jp2)
      ! Saves projections on velocity fields from jp1 to jp2 to file
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
C       include 'Shapes_DEF'
      logical , save :: firstCall = .true.
      character(len=100) fileName
      integer i,j,jp1,jp2

      real projections(nShapes,jp1:jp2) 
      do i=jp1,jp2
        projections(:,i) = ProjectOnShapes(i)
      enddo

      if (NIO==0 ) then
       do i = 1,npert
          write(fileName,'(A,I0.2,A)') 'ProjShapes_',i,'.dat'
          if (firstCall ) then
            open(unit=88, file=fileName,
     $        status="replace", action="write")
            write(88,'(A)', advance="no") 'time' 
            do j = 1,nShapes
              write(88,'(A,I0.2)',advance="no") '   ProjOnShape_' , i
            enddo ! header loop
            write(88,*) ' ' 
          else
            open(unit=88, file=fileName,
     $       status='old',action='write',position="append")
          endif !firs tcall
c          write(*,*) time , nShapes , projections(1:nShapes,i)
          write(88,*) time , projections(1:nShapes,i)
          close(88)
        enddo ! perturbations loop
        endif  ! NIO 
      firstCall = .false.
      return
      end subroutine

c========================================================
c     Get Number of Shapes shape
c
c========================================================
      function GetNumShapes()
        integer GetNumShapes
        GetNumShapes=nShapes
        return 
      end function





      end module
