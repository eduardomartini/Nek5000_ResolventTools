c========================================================
c     Load_FLT_TO Module 
c
c     Loads data from snapshots and saves to arrays. 
c     It backups simulation fields, uses uses neks "loadfld"
c     uses neks "loadfld" subroutine to laod the snapshot
c     into the simulation fields, copy them to the desired
c     arrays and restore the simulation status.
c
c     Subroutines
c           Load_FLD_To            : Load velocity fields of a file into the target arrays.
c           Load_FLD_Time_To_wTime : Load velocity fields of a file into the target arrays. Does not restore simulation time
c
c========================================================

      subroutine Load_FLD_To(fldFile,tarVx,tarVy,tarVz)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      CHARACTER(LEN=100) :: fldFile
      real  tarVx(lx1*ly1*lz1*lelt)  ,
     $     tarVy(lx1*ly1*lz1*lelt) ,  tarVz(lx1*ly1*lz1*lelt)  
      real dummyT(lx1*ly1*lz1*lelt*ldimt)
      real dummyP(lx2*ly2*lz2*lelt)
      
      prevTime = time
      call  Load_FLD_Time_To_wTime(fldFile,tarVx,tarVy,tarVz,
     $                                     dummyP,dummyT)
      time = prevTime
      return
      end subroutine
      
      
      subroutine Load_FLD_Time_To_wTime(fldFile,tarVx,tarVy,
     $                                  tarVz,tarP,tarT)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      INCLUDE 'RESTART'
      
      integer n,n2,i ,prevjp
      CHARACTER(LEN=100) :: fldFile
      real,intent(out) :: tarVx(lx1,ly1,lz1,lelt)  ,
     $     tarVy(lx1,ly1,lz1,lelt),  tarVz(lx1,ly1,lz1,lelt),  
     $     tarP(lx2,ly2,lz2,lelt) ,  tarT(lx1,ly1,lz1,lelt,ldimt)  

      real tmp_mesh(lx1,ly1,lz1,lelt,ldim)
      real tmp_mesh2(lx2,ly2,lz2,lelt,ldim)
      real tmp_vp(lx1*ly1*lz1*lelt,lpert,ldim)
      real tmp_v(lx1,ly1,lz1,lelt,ldim) !,tmp_vp(lx1,ly1,lz1,lelt,lpert,4)
      real tmp_p(lx2,ly2,lz2,lelt)  
      real tmp_t(lx1,ly1,lz1,lelt,ldimt)
      real tmp_pp(lx2*ly2*lz2*lelt,lpert) 
      real tmp_tp(lx1*ly1*lz1*lelt,ldimt,lpert)


c     Backup Flow
      prevjp = jp

      ! backup mesh
      tmp_mesh(:,:,:,:,1) = xm1
      tmp_mesh(:,:,:,:,2) = ym1
      if (ndim==3)  tmp_mesh(:,:,:,:,ndim) = zm1
      tmp_mesh2(:,:,:,:,1) = xm2
      tmp_mesh2(:,:,:,:,2) = ym2
      if (ndim==3)  tmp_mesh2(:,:,:,:,ndim) = zm2
      
      ! backup velocities
      tmp_v(:,:,:,:,1)=VX
      tmp_v(:,:,:,:,2)=VY
      if (ndim == 3)  tmp_v(:,:,:,:,ndim) = VZ
      
      ! backup pressure and scalars
      tmp_p = Pr
      tmp_t = T
      ! backup perturbations
      tmp_vp(:,:,1) = VXP
      tmp_vp(:,:,2) = VYP
      if (ndim ==3 ) tmp_vp(:,:,ndim) = VZP
      tmp_pp = Prp
      tmp_tp = Tp

c     read file 
      call load_fld(fldFile)

c     copy fields to target
      tarVx = Vx
      tarVy = Vy
      if (ndim == 3)  tarVz=Vz
      tarP = Pr
      tarT = T

c     Restore Velocities
      VX = tmp_v(:,:,:,:,1)
      VY = tmp_v(:,:,:,:,2)
      if (ndim == 3 )  Vz = tmp_v(:,:,:,:,ndim)
      Pr = tmp_p
      T = tmp_t

      VXP=tmp_vp(:,:,1)
      VYP=tmp_vp(:,:,2)
      if (ndim==3) VZP=tmp_vp(:,:,ndim)
      Prp = tmp_pp
      Tp  = tmp_tp

      xm1 = tmp_mesh(:,:,:,:,1) 
      ym1 = tmp_mesh(:,:,:,:,2) 
      if (ndim==3)  zm1 = tmp_mesh(:,:,:,:,ndim) 
      xm2 = tmp_mesh2(:,:,:,:,1) 
      ym2 = tmp_mesh2(:,:,:,:,2) 
      if (ndim==3) zm2 = tmp_mesh2(:,:,:,:,ndim) 
      
      param(59) = 1. ! Force nek5 to  recognize element deformation.
      jp = prevjp

      return
      end subroutine
