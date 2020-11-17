!-----------------------------------------------------------------------
!
      module mrbpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mrbpush2lib.f
! mrbpush2mod.f contains multi-tasking interface procedures to process
!              relativistic particles with magnetic fields:
!              defines module mrbpush2d
! rdjpost => imgrjpost2 deposits relativistic current density, with
!            various interpolations and optimizations.
!            calls MGRJPOST2, MGSRJPOST2, MGSRJPOST2X, MGRJPOST2L,
!            MGSRJPOST2L, MGSRJPOST2XL, MGRJPOST2C, MGRJPOST22,
!            MGSRJPOST22, MGSRJPOST22X, MGRJPOST22L, MGSRJPOST22L,
!            MGSRJPOST22XL, or MGRJPOST22C
! rpush => imgrpush2 push relativistic particles with 2 component
!          electric field, with various interpolations and optimizations
!          calls MGRPUSH2, MGSRPUSH2, MGRPUSH2L, MGSRPUSH2L,
!          or MGRPUSH2C
! rpush => imgrbpush2 push relativistic particles with magnetic field
!          and 2 component electric field, with various interpolations
!          and optimizations.
!          calls MGRBPUSH2, MGSRBPUSH2, MGRBPUSH2L, MGSRBPUSH2L,
!          MGRBPUSH2C, MGRBPUSH22, MGSRBPUSH22, MGRBPUSH22L,
!          MGSRBPUSH22L, or MGRBPUSH22C
! rpush3 => imgrbpush23 push relativistic particles with magnetic field
!           and 3 component electric field, with various interpolations
!           and optimizations.
!           calls MGRBPUSH23, MGSRBPUSH23, MGRBPUSH23L, MGSRBPUSH23L,
!           MGRBPUSH23C, MGRBPUSH22, MGSRBPUSH22, MGRBPUSH22L,
!           MGSRBPUSH22L, or MGRBPUSH22C
! rpushzf => imrpush2zf push 2d relativistic particles with no forces.
!            calls MRPUSH2ZF
! rpush3zf => imrpush23zf push 2-1/2d relativistic particles with no
!             forces.
!             calls MRPUSH23ZF, or MRPUSH2ZF
! rdjpostgl => imrdjpost2gl deposits relativistic current density, using
!              gridless method.
!              calls MRDJPOST2GL
! rpushgl => imrpush2gl push relativistic particles with 2d electrostatic
!            fields, using gridless method.
!            calls MRPUSH2GL
! rpushglx => imrpush2glx push relativistic particles with 2d
!             electrostatic fields, using optimized gridless method.
!             calls MRPUSH2GLX, or MRPUSH2GL
! rpush3gl => imrbpush23gl push relativistic particles with 2-1/2d
!             electromagnetic fields, using gridless method.
!             calls MRBPUSH23GL
! rgcjpost => imgrcjpost2 deposits time-centered relativistic current
!             density with 2d electrostatic fields.
!             calls MGRCJPOST2, MGRCJPOST2L, or MGRCJPOST2C
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 28, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use rbpush2d, only: wtimer, retard, icptov2
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: rdjpost, rpush, rpush3, retard, icptov2
      public :: rpushzf, rpush3zf
      public :: rpushgl, rpushglx, rdjpostgl, rpush3gl, rgcjpost
!
! buffer data for current deposit
      real, dimension(:,:,:,:), allocatable :: cup
      integer :: szbuf = 0
      save
! define interface to original Fortran77 procedures
!
      interface
         subroutine MGRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: cu
         real, dimension(3,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST2X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST2XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nx&
     &yv,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: cu
         real, dimension(2,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nx&
     &yv,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxyv) :: cu
         real, dimension(2*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv&
     &,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nx&
     &yv,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &xyv,ipbc,cup,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         real, dimension(3*nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv&
     &,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nxv,ny) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nxv,ny) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRPUSH2C(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,ny) :: fxy
         real, dimension(3,nxv,ny) :: bxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRBPUSH22C(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nxv,ny) :: bz
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,idt&
     &ask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc, nmt, ierr
         real  :: dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRPUSH23ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,id&
     &task,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc, nmt, ierr
         real  :: dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nx&
     &vh,nyv,ipbc,cup,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc, nmt, ierr
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: cu
         complex, dimension(nxvh) :: sctx
         real, dimension(3,2*nxvh,nyv,nmt) :: cup
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRPUSH2GL(part,fxy,sctx,qbm,dt,ci,ek,idimp,nop,nx,ny&
     &,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc, nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh) :: sctx
         complex, dimension(nxvh,nmt) :: sctxp
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,no&
     &p,nx,ny,nxvh,nyv,nnp,ipbc,sctxpp,exypp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nnp, ipbc, nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh,nnp) :: sctxp
         double precision, dimension(2,nnp) :: exyp
         complex, dimension(nxvh,nnp,nmt) :: sctxpp
         double precision, dimension(2,nnp,nmt) :: exypp
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp&
     &,nop,nx,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc, nmt, ierr
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: fxy, bxy
         complex, dimension(nxvh) :: sctx
         complex, dimension(nxvh,nmt) :: sctxp
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRCJPOST2(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,ny&
     &v,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv, nmt, ierr
         real  :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRCJPOST2L(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,n&
     &yv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv, nmt, ierr
         real  :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRCJPOST2C(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,n&
     &yv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv, nmt, ierr
         real  :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdjpost
         module procedure imgrjpost2
      end interface
!
      interface rdjpostgl
         module procedure imrdjpost2gl
      end interface
!
      interface rpush
         module procedure imgrpush2
         module procedure imgrbpush2
      end interface
!
      interface rpush3
         module procedure imgrbpush23
      end interface
!
      interface rpushzf
         module procedure imrpush2zf
      end interface
!
      interface rpush3zf
         module procedure imrpush23zf
      end interface
!
      interface rpushgl
         module procedure imrpush2gl
      end interface
!
      interface rpushglx
         module procedure imrpush2glx
      end interface
!
      interface rpush3gl
         module procedure imrbpush23gl
      end interface
!
      interface rgcjpost
         module procedure imgrcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imgrpush2(part,fxy,nop,qbm,dt,ci,ek,tpush,nx,ny,ipbc&
     &,inorder,popt)
! multi-tasking relativistic particle push with 2d electrostatic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
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
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MGSRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv&
     &,nxyv,ipbc,ekp,idtask,nmt,ierr)
            else
               call MGRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nyv,ipbc,ekp,idtask,nmt,ierr)
            endif
         else if (order==CUBIC) then
            call MGRPUSH2C(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv&
     &,ipbc,ekp,idtask,nmt,ierr)
         else
            if (opt==LOOKAHEAD) then
               call MGSRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc,ekp,idtask,nmt,ierr)
            else
               call MGRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,ekp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgrpush2
!
         subroutine imgrjpost2(part,cu,nop,qm,dt,ci,tdjpost,nx,ny,ipbc,i&
     &norder,djopt)
! multi-tasking, relativisitc current deposit
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
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
                  call MGSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv&
     &,nxyv,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nx&
     &v,nxyv,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nyv,ipbc,cup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRJPOST22C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv&
     &,ipbc,cup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nxyv,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv&
     &,nxyv,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &yv,ipbc,cup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nxyv,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSRJPOST2XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv&
     &,nxyv,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &yv,ipbc,cup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRJPOST2C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc,cup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &xyv,ipbc,cup,idtask,nmt,ierr)
               else if (opt==VECTOR) then
                  call MGSRJPOST2X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nxyv,ipbc,cup,idtask,nmt,ierr)
               else
                  call MGRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,ny&
     &v,ipbc,cup,idtask,nmt,ierr)
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
         end subroutine imgrjpost2
!
         subroutine imgrbpush2(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,tpush,n&
     &x,ny,ipbc,inorder,popt)
! multi-tasking, relativistic particle push with 2d electromagnetic
! fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
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
                  call MGSRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
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
         end subroutine imgrbpush2
!
         subroutine imgrbpush23(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,tpush,&
     &nx,ny,ipbc,inorder,popt)
! multi-tasking, relativistic particle push with 2-1/2d electromagnetic
! fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
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
                  call MGSRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
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
         end subroutine imgrbpush23
!
         subroutine imrpush2zf(part,nop,dt,ci,ek,tpush,nx,ny,ipbc)
! multi-tasking relativistic 2d particle push with no forces
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, nmt, ltime, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nmt = ntasks
! initialize timer
         call wtimer(tp,ltime,-1)
         call MRPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,idtask,nm&
     &t,ierr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imrpush2zf
!
         subroutine imrpush23zf(part,nop,dt,ci,ek,tpush,nx,ny,ipbc,ndim)
! multi-tasking relativistic particle push with no forces
         implicit none
         integer :: nop, nx, ny, ipbc, ndim
         real :: dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, nmt, ltime, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nmt = ntasks
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(ndim)
         case (2)
            call MRPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,idtask&
     &,nmt,ierr)
         case (3)
            call MRPUSH23ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,idtas&
     &k,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imrpush23zf
!
         subroutine imrdjpost2gl(part,cu,nop,qm,dt,ci,nx,ny,ipbc,tdjpost&
     &)
! multi-tasking deposit current using gridless method
! for relativistic particles
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qm, dt, ci, tdjpost
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
            call MRDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nxvh,&
     &nyv,ipbc,cup,sctxp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine imrdjpost2gl
!
         subroutine imrpush2gl(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpus&
     &h)
! multi-tasking relativistic particle push with 2d electrostatic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
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
         call MRPUSH2GL(part,fxy,sctx,qbm,dt,ci,ek,idimp,nop,nx,ny,nxvh,&
     &nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imrpush2gl
!
         subroutine imrpush2glx(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpu&
     &sh,popt)
! multi-tasking relativistic particle push with 2d electrostatic fields
! using gridless method, optimized version
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ci, ek, tpush
         integer, optional :: popt
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer, parameter :: nnp = 128
         integer :: idimp, nxvh, nyv, nmt, opt, ltime, ierr
         real :: tp
         complex, dimension(size(fxy,2)/2,nnp) :: sctxp
         double precision, dimension(2,nnp) :: exyp
         complex, dimension(size(fxy,2)/2,nnp,ntasks) :: sctxpp
         double precision, dimension(2,nnp,ntasks) :: exypp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         nmt = ntasks
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         if (opt==LOOKAHEAD) then
            call MRPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,nop,n&
     &x,ny,nxvh,nyv,nnp,ipbc,sctxpp,exypp,ekp,idtask,nmt,ierr)
         else
            call MRPUSH2GL(part,fxy,sctxp,qbm,dt,ci,ek,idimp,nop,nx,ny,n&
     &xvh,nyv,ipbc,sctxpp,ekp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imrpush2glx
!
         subroutine imrbpush23gl(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,nx,ny&
     &,ipbc,tpush)
! multi-tasking particle push with 2-1/2d electromagnetic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, dtc, ci, ek, tpush
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
            call MRBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imrbpush23gl
!
         subroutine imgrcjpost2(part,fxy,cu,nop,qm,qbm,dt,ci,tdcjpost,in&
     &order)
! multi-tasking current density deposit for relativistic particles
! with 2d electrostatic fields
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, cu
! local data
         integer :: idimp, nxv, nyv, nnxyv, nmt, ltime, order, ierr
         real :: tdc
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3)
         nnxyv = size(cu,1)*nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! check if size of buffer has changed
         if (szbuf < nnxyv) then
            if (szbuf /= 0) deallocate(cup)
! allocate buffer
            allocate(cup(size(cu,1),nxv,nyv,ntasks))
            szbuf = nnxyv
         endif
! initialize timer
         call wtimer(tdc,ltime,-1)
         if (order==LINEAR) then
            call MGRCJPOST2L(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv,&
     &cup,idtask,nmt,ierr)
         else if (order==CUBIC) then
            call MGRCJPOST2C(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv,&
     &cup,idtask,nmt,ierr)
         else
            call MGRCJPOST2(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv,c&
     &up,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgrcjpost2
!
      end module mrbpush2d
