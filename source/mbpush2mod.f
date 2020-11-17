!-----------------------------------------------------------------------
!
      module mbpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mbpush2lib.f
! mbpush2mod.f contains multi-tasking interface procedures to process
!              particles with magnetic fields:
!              defines module mbpush2d
! djpost => imgjpost2 deposits current density, with various
!           interpolations and optimizations.
!           calls MGJPOST2, MGSJPOST2, MGSJPOST2X, MGJPOST2L,
!           MGSJPOST2L, MGSJPOST2XL, MGJPOST2C, MGJPOST22, MGSJPOST22,
!           MGSJPOST22X, MGJPOST22L, MGSJPOST22L, MGSJPOST22XL, or 
!           MGJPOST22C
! push => imgbpush2 push particles with magnetic field and 2 component
!         electric field, with various interpolations and optimizations.
!         calls MGBPUSH2, MGSBPUSH2, MGBPUSH2L, MGSBPUSH2L, MGBPUSH2C,
!         MGBPUSH22, MGSBPUSH22, MGBPUSH22L, MGSBPUSH22L, or MGBPUSH22C
! push3 => imgbpush23 push particles with magnetic field and 3 component
!          electric field, with various interpolations and optimizations
!          calls MGBPUSH23, MGSBPUSH23, MGBPUSH23L, MGSBPUSH23L,
!          MGBPUSH23C, MGBPUSH22, MGSBPUSH22, MGBPUSH22L, MGSBPUSH22L,
!          or MGBPUSH22C
! djpostgl => imdjpost2gl deposits current density, using gridless
!             method.
!             calls MDJPOST2GL
! push3gl => imbpush23gl push particles with 2-1/2d electromagnetic
!            fields, using gridless method.
!            calls MBPUSH23GL
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 28, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use bpush2d, only: wtimer, retard
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: djpost, push, push3, retard, djpostgl, push3gl
!
! buffer data for current deposit
      real, dimension(:,:,:,:), allocatable :: cup
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MDJPOST2(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv,&
     &cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: cux, cuy, cuz
         real, dimension(nxv,ny,3,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,&
     &cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipb&
     &c,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: cu
         real, dimension(3,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy,cup&
     &,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxvy) :: cu
         real, dimension(3*nxvy,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ip&
     &bc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDJPOST2L(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv&
     &,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: cux, cuy, cuz
         real, dimension(nxv,ny,3,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ip&
     &bc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy,cu&
     &p,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxvy) :: cu
         real, dimension(3*nxvy,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,i&
     &pbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ip&
     &bc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: cu
         real, dimension(2,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST22X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,i&
     &pbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxyv) :: cu
         real, dimension(2*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipb&
     &c,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,i&
     &pbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSJPOST22XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,&
     &ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGJPOST2C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGJPOST22C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipb&
     &c,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MBPUSH2(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MBPUSH2L(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MBPUSH2CQ(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH2CQ(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nx&
     &v,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MBPUSH2CL(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH2CL(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nx&
     &v,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nxv,ny) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nxv,ny) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGBPUSH22C(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nxv,ny) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxvh,n&
     &yv,ipbc,cup,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc, nmt, ierr
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: cu
         complex, dimension(nxvh) :: sctx
         real, dimension(3,2*nxvh,nyv,nmt) :: cup
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop&
     &,nx,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: fxy, bxy
         complex, dimension(nxvh) :: sctx
         complex, dimension(nxvh,nmt) :: sctxp
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure imgjpost2
      end interface
!
      interface djpostgl
         module procedure imdjpost2gl
      end interface
!
      interface push
         module procedure imgbpush2
!        module procedure imgbpush2cq
      end interface
!
      interface push3
         module procedure imgbpush23
      end interface
!
      interface push3gl
         module procedure imbpush23gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imgjpost2(part,cu,nop,qm,dt,tdjpost,nx,ny,ipbc,inord&
     &er,djopt)
! multi-tasking current deposit
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, opt, ierr
         integer :: nnxyv
         real :: tj
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(cu,2); nyv = size(cu,3)
         nxyv = nxv*nyv; nnxyv = size(cu,1)*nxyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of buffer has changed
         if (szbuf < nnxyv) then
            if (szbuf /= 0) deallocate(cup)
! allocate buffer
            allocate(cup(size(cu,1),nxv,nyv,ntasks))
            szbuf = nnxyv
         endif
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSJPOST22XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nx&
     &yv,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc,cup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGJPOST22C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipb&
     &c,cup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSJPOST22X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc,cup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc,cup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGJPOST2C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &,cup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,&
     &ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ip&
     &bc,cup,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine imgjpost2
!
         subroutine imgbpush2(part,fxy,bxy,nop,qbm,dt,dtc,ek,tpush,nx,ny&
     &,ipbc,inorder,popt)
! multi-tasking particle push with 2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, opt, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3)
         nxyv = nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif

            else if (order==CUBIC) then
               call MGBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)

            else
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgbpush2
!
         subroutine imgbpush2cq(part,fxy,bxy,nop,qbm,dt,ek,tpush,nx,ny,i&
     &pbc,inorder)
! multi-tasking particle push with 2d electromagnetic fields,
! with correction to Boris Mover
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3)
         nxyv = nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tp,ltime,-1)
         if (order==LINEAR) then
            call MGBPUSH2CL(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,ekp,idtask,nmt,ierr)
         else
            call MGBPUSH2CQ(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,ekp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgbpush2cq
!
         subroutine imgbpush23(part,fxy,bxy,nop,qbm,dt,dtc,ek,tpush,nx,n&
     &y,ipbc,inorder,popt)
! multi-tasking particle push with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, opt, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3)
         nxyv = nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgbpush23
!
         subroutine imdjpost2gl(part,cu,nop,qm,dt,nx,ny,ipbc,tdjpost)
! multi-tasking deposit current using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qm, dt, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: idimp, nxv, nxvh, nyv, nnxyv, nmt, ltime, ierr
         real :: tj
         complex, dimension(size(cu,2)/2) :: sctx
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
         complex, dimension(size(cu,2)/2,ntasks) :: sctxp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(cu,2); nxvh = nxv/2; nyv = size(cu,3)
         nnxyv = size(cu,1)*nxv*nyv
         nmt = ntasks
! check if size of buffer has changed
         if (szbuf < nnxyv) then
            if (szbuf /= 0) deallocate(cup)
! allocate buffer
            allocate(cup(size(cu,1),nxv,nyv,ntasks))
            szbuf = nnxyv
         endif
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(cu,1))
         case (3)
            call MDJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxvh,nyv,&
     &ipbc,cup,sctxp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine imdjpost2gl
!
         subroutine imbpush23gl(part,fxy,bxy,nop,qbm,dt,dtc,ek,nx,ny,ipb&
     &c,tpush)
! multi-tasking particle push with 2-1/2d electromagnetic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, dtc, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxvh, nyv, nmt, ltime, ierr
         real :: tp
         complex, dimension(size(fxy,2)/2) :: sctx
         complex, dimension(size(fxy,2)/2,ntasks) :: sctxp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         nmt = ntasks
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(bxy,1))
         case (3)
            call MBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imbpush23gl
!
      end module mbpush2d
