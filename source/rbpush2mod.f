!-----------------------------------------------------------------------
!
      module rbpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library rbpush2lib.f
! rbpush2mod.f contains interface procedures to process relativistic
!              particles with magnetic fields:
!              defines module rbpush2d
! rdjpost => igrjpost2 deposits relativistic current density, with
!            various interpolations and optimizations.
!            calls GRJPOST2, GSRJPOST2, GSRJPOST2X, GRJPOST2L,
!            GSRJPOST2L, GSRJPOST2XL, GRJPOST2C, GRJPOST22, GSRJPOST22,
!            GSRJPOST22X, GRJPOST22L, GSRJPOST22L, GSRJPOST22XL, or
!            GRJPOST22C
! rpush => igrpush2 push relativistic particles with 2 component
!          electric field, with various interpolations and optimizations
!          calls GRPUSH2, GSRPUSH2, GRPUSH2L, GSRPUSH2L, or GRPUSH2C
! rpush => igrbpush2 push relativistic particles with magnetic field and
!          2 component electric field, with various interpolations and
!          optimizations.
!          calls GRBPUSH2, GSRBPUSH2, GRBPUSH2L, GSRBPUSH2L, GRBPUSH2C,
!          GRBPUSH22, GSRBPUSH22, GRBPUSH22L, GSRBPUSH22L, or GRBPUSH22C
! rpush3 => igrbpush23 push relativistic particles with magnetic field
!           and 3 component electric field, with various interpolations
!           and optimizations.
!           calls GRBPUSH23, GSRBPUSH23, GRBPUSH23L, GSRBPUSH23L,
!           GRBPUSH23C, GRBPUSH22, GSRBPUSH22, GRBPUSH22L, GSRBPUSH22L,
!           or GRBPUSH22C
! rpushzf => irpush2zf push 2d relativistic particles with no forces.
!            calls RPUSH2ZF
! rpush3zf => irpush23zf push 2-1/2d relativistic particles with no
!             forces.
!             calls RPUSH23ZF, or RPUSH2ZF
! retard => irretard2 retard relativistic particle position a half
!           time-step.
!           calls RRETARD2, or RRETARD22
! icptov2 converts momentum to velocity for relativistic particles.
!         calls CPTOV2, or CPTOV22
! rdjpostgl => irdjpost2gl deposits relativistic current density, using
!              gridless method.
!              calls RDJPOST2GL
! rpushgl => irpush2gl push relativistic particles with 2d electrostatic
!            fields, using gridless method.
!            calls RPUSH2GL
! rpushglx => irpush2glx push relativistic particles with 2d
!             electrostatic fields, using optimized gridless method.
!             calls RPUSH2GLX, or RPUSH2GL
! rpush3gl => irbpush23gl push relativistic particles with 2-1/2d
!             electromagnetic fields, using gridless method.
!             calls RBPUSH23GL
! rgcjpost => igrcjpost2 deposits time-centered relativistic current
!             density with 2d electrostatic fields.
!             calls GRCJPOST2, GRCJPOST2L, or GRCJPOST2C
!*****
! ordjpost => calls iogrjpost2 deposits OSIRIS charge conserving
!             current with various interpolation schemes. Calls PSPLIT2
!             before calling OGRJPOST2, 2L, 2C to deposit the current 
!             to the virtual current array 
!*****
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: september 24, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, rdjpost, rpush, rpush3, retard, icptov2
      public :: ordjpost, iogrjpost2
      public :: rpushzf, rpush3zf
      public :: rdjpostgl, rpushgl, rpushglx, rpush3gl, rgcjpost
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ip&
     &bc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,&
     &ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST2X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST2XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nx&
     &yv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GRJPOST2C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GRJPOST22C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GSRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GSRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GRPUSH2C(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBPUSH22C(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine RRETARD2(part,dtc,ci,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc
         real :: dtc, ci
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine RRETARD22(part,dtc,ci,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc
         real :: dtc, ci
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine CPTOV2(part,ci,nop,idimp)
         integer :: nop, idimp
         real :: ci
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine CPTOV22(part,ci,nop,idimp)
         integer :: nop, idimp
         real :: ci
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine RPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc
         real :: dt, ci, ek
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine RPUSH23ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc
         real :: dt, ci, ek
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine RDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nxv&
     &h,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: cu
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine RPUSH2GL(part,fxy,sctx,qbm,dt,ci,ek,idimp,nop,nx,ny,&
     &nxvh,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine RPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxvh,nyv,npp,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, npp, ipbc
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh,npp) :: sctxp
         double precision, dimension(2,npp) :: exyp
         end subroutine
      end interface
      interface
         subroutine RBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,&
     &nop,nx,ny,nxvh,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt, dtc, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: fxy, bxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GRCJPOST2(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv&
     &)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, cu
         end subroutine
      end interface
      interface
         subroutine GRCJPOST2L(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,ny&
     &v)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, cu
         end subroutine
      end interface
      interface
         subroutine GRCJPOST2C(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,ny&
     &v)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, cu
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdjpost
         module procedure igrjpost2
      end interface
!
      interface ordjpost
         module procedure iogrjpost2
      end interface
!
      interface rdjpostgl
         module procedure irdjpost2gl
      end interface
!
      interface rpush
         module procedure igrpush2
         module procedure igrbpush2
      end interface
!
      interface rpush3
         module procedure igrbpush23
      end interface
!
      interface rpushzf
         module procedure irpush2zf
      end interface
!
      interface rpush3zf
         module procedure irpush23zf
      end interface
!
      interface retard
         module procedure irretard2
      end interface
!
      interface rpushgl
         module procedure irpush2gl
      end interface
!
      interface rpushglx
         module procedure irpush2glx
      end interface
!
      interface rpush3gl
         module procedure irbpush23gl
      end interface
!
      interface rgcjpost
         module procedure igrcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igrpush2(part,fxy,nop,qbm,dt,ci,ek,tpush,nx,ny,ipbc,&
     &inorder,popt)
! push relativistic particles with 2d electrostatic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call GSRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc)
            else
               call GRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
            endif
         else if (order==CUBIC) then
            call GRPUSH2C(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc)
         else
            if (opt==LOOKAHEAD) then
               call GSRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc)
            else
               call GRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
            endif
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrpush2
!
         subroutine igrjpost2(part,cu,nop,qm,dt,ci,tdjpost,nx,ny,ipbc,in&
     &order,djopt)
! deposit relativistic current
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tj
         idimp = size(part,1)
         nxv = size(cu,2); nyv = size(cu,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nxyv,ipbc)
               else if (opt==VECTOR) then
                  call GSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv&
     &,nxyv,ipbc)
               else
                  call GRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &yv,ipbc)
               endif
            else if (order==CUBIC) then
               call GRJPOST22C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
     &ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &xyv,ipbc)
               else if (opt==VECTOR) then
                  call GSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nxyv,ipbc)
               else
                  call GRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,ny&
     &v,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &xyv,ipbc)
               else if (opt==VECTOR) then
                  call GSRJPOST2XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
     &nxyv,ipbc)
               else
                  call GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,ny&
     &v,ipbc)
               endif
            else if (order==CUBIC) then
               call GRJPOST2C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nx&
     &yv,ipbc)
               else if (opt==VECTOR) then
                  call GSRJPOST2X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
     &xyv,ipbc)
               else
                  call GRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv&
     &,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine igrjpost2
!
         subroutine iogrjpost2(part,vpart,cu,vcu,nop,qm,dth,ci,tdjpost&
        &,nx,ny,ipbc,inorder,djopt,dt)
! deposit relativistic charge conserving current following OSIRIS scheme
         implicit none
         integer :: nop, nx, ny, ipbc!, inorder
         integer, optional :: djopt, inorder
         real :: qm, dth, dt, ci, tdjpost
         real, dimension(:,:), pointer :: part,vpart
         real, dimension(:,:,:), pointer :: cu, vcu
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime, nsplit
         integer :: j, cj
         real :: tj
         !real, dimension(size(cu,2)) :: cuTy
         !real, dimension(size(cu,3)) :: cuTx
         idimp = size(part,1)
         nxv = size(cu,2); nyv = size(cu,3); nxyv = nxv*nyv
         !order = inorder
         order = QUADRATIC
         if (present(inorder)) order = inorder
         !order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,ltime,-1)

! split the particles
         !print*, 'order1 =', order
         call PSPLIT2(part,vpart,nop,dt,nsplit,ci,idimp,nx,ny, &
        &ipbc,order)
         !print*, 'order3 =', order

         select case(size(cu,1))
         case (2)
!            if (order==LINEAR) then
!               if (opt==LOOKAHEAD) then
!                  call GSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
!     &nxyv,ipbc)
!               else if (opt==VECTOR) then
!                  call GSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv&
!     &,nxyv,ipbc)
!               else
!                  call GRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
!     &yv,ipbc)
!               endif
!            else if (order==CUBIC) then
!               call GRJPOST22C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,&
!     &ipbc)
!            else
!               if (opt==LOOKAHEAD) then
!                  call GSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,n&
!     &xyv,ipbc)
!               else if (opt==VECTOR) then
!                  call GSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,&
!     &nxyv,ipbc)
!               else
!                  call GRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,ny&
!     &v,ipbc)
!               endif
!            endif
         case (3)
            if (order==LINEAR) then
                 call OGRJPOST2L(part,vpart,vcu,qm,dth,ci,nsplit,nop,idimp&
       &,nx,ny,nxv,nyv,ipbc,dt)
               else if(order==CUBIC) then
		   call OGRJPOST2C(part,vpart,vcu,qm,dth,ci,nsplit,nop,idimp&
       &,nx,ny,nxv,nyv,ipbc,dt)
               else if (order==2) then
                 call OGRJPOST2(part,vpart,vcu,qm,dth,ci,nsplit,nop,idimp&
       &,nx,ny,nxv,nyv,ipbc,dt)
               endif
         end select
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine iogrjpost2
!
         subroutine igrbpush2(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,tpush,nx&
     &,ny,ipbc,inorder,popt)
! push relativistic particles with 2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3); nxyv = nxv*nyv
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
                  call GSRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GRBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GRBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrbpush2
!
         subroutine igrbpush23(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,tpush,n&
     &x,ny,ipbc,inorder,popt)
! push relativistic particles with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3); nxyv = nxv*nyv
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
                  call GSRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GRBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH22(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GRBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ipbc)
               else
                  call GRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrbpush23
!
         subroutine irretard2(part,nop,dtc,ci,nx,ny,ipbc,ndim)
! retards relativistic particle positions half time-step
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: ndim
         real :: dtc, ci
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, nd
         idimp = size(part,1)
         nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call RRETARD22(part,dtc,ci,idimp,nop,nx,ny,ipbc)
         case (3)
            call RRETARD2(part,dtc,ci,idimp,nop,nx,ny,ipbc)
         end select
         end subroutine irretard2
!
         subroutine icptov2(part,nop,ci,ndim)
! convert momentum to velocity for relativistic particles
         implicit none
         integer :: nop
         integer, optional :: ndim
         real :: ci
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, nd
         idimp = size(part,1); nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call CPTOV22(part,ci,nop,idimp)
         case (3)
            call CPTOV2(part,ci,nop,idimp)
         end select
         end subroutine icptov2
!
         subroutine irpush2zf(part,nop,dt,ci,ek,tpush,nx,ny,ipbc)
! push 2d relativistic particles with no forces
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, ltime
         real :: tp
         idimp = size(part,1)
! initialize timer
         call wtimer(tp,ltime,-1)
         call RPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine irpush2zf
!
         subroutine irpush23zf(part,nop,dt,ci,ek,tpush,nx,ny,ipbc,ndim)
! push relativistic particles with no forces
         implicit none
         integer :: nop, nx, ny, ipbc, ndim
         real :: dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, ltime
         real :: tp
         idimp = size(part,1)
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(ndim)
         case (2)
            call RPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
         case (3)
            call RPUSH23ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine irpush23zf
!
         subroutine irdjpost2gl(part,cu,nop,qm,dt,ci,nx,ny,ipbc,tdjpost)
! deposit current using gridless method for relativistic particles
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qm, dt, ci, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tj
         complex, dimension(size(cu,2)/2) :: sctx
         idimp = size(part,1)
         nxvh = size(cu,2)/2; nyv = size(cu,3)
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(cu,1))
         case (3)
            call RDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nxvh,n&
     &yv,ipbc)
         end select
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine irdjpost2gl
!
         subroutine irpush2gl(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpush&
     &)
! push relativistic particles with 2d electrostatic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tp
         complex, dimension(size(fxy,2)/2) :: sctx
         idimp = size(part,1)
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
! initialize timer
         call wtimer(tp,ltime,-1)
         call RPUSH2GL(part,fxy,sctx,qbm,dt,ci,ek,idimp,nop,nx,ny,nxvh,n&
     &yv,ipbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine irpush2gl
!
         subroutine irpush2glx(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpus&
     &h,popt)
! push relativistic particles with 2d electrostatic fields
! using gridless method, optimized method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ci, ek, tpush
         integer, optional :: popt
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer, parameter :: npp = 32
         integer :: idimp, nxvh, nyv, opt, ltime
         real :: tp
         complex, dimension(size(fxy,2)/2,npp) :: sctxp
         double precision, dimension(2,npp) :: exyp
         idimp = size(part,1)
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         if (opt==LOOKAHEAD) then
            call RPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,nop,nx&
     &,ny,nxvh,nyv,npp,ipbc)
         else
            call RPUSH2GL(part,fxy,sctxp,qbm,dt,ci,ek,idimp,nop,nx,ny,nx&
     &vh,nyv,ipbc)
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine irpush2glx
!
         subroutine irbpush23gl(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,nx,ny,&
     &ipbc,tpush)
! push relativistic particles with 2-1/2d electromagnetic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tp
         complex, dimension(size(fxy,2)/2) :: sctx
         idimp = size(part,1)
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(bxy,1))
         case (3)
            call RBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,nop&
     &,nx,ny,nxvh,nyv,ipbc)
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine irbpush23gl
!
         subroutine igrcjpost2(part,fxy,cu,nop,qm,qbm,dt,ci,tdcjpost,ino&
     &rder)
! deposit relativistic current density with 2d electrostatic fields
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, cu
! local data
         integer :: idimp, nxv, nyv, order, ltime
         real :: tdc
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tdc,ltime,-1)
         if (order==LINEAR) then
            call GRCJPOST2L(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv)
         else if (order==CUBIC) then
            call GRCJPOST2C(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv)
         else
            call GRCJPOST2(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv)
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igrcjpost2
!
      end module rbpush2d
