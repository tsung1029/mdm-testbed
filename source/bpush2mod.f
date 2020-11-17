!-----------------------------------------------------------------------
!
      module bpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library bpush2lib.f
! bpush2mod.f contains interface procedures to process particles with
!             magnetic fields:
!             defines module bpush2d
! djpost => igjpost2 deposits current density, with various
!           interpolations and optimizations.
!           calls GJPOST2, GSJPOST2, GSJPOST2X, GJPOST2L, GSJPOST2L,
!           GSJPOST2XL, GJPOST2C, GJPOST22, GSJPOST22, GSJPOST22X,
!           GJPOST22L, GSJPOST22L, GSJPOST22XL, or GJPOST22C
! push => igbpush2 push particles with magnetic field and 2 component
!         electric field, with various interpolations and optimizations.
!         calls GBPUSH2, GSBPUSH2, GBPUSH2L, GSBPUSH2L, GBPUSH2C,
!         GBPUSH22, GSBPUSH22, GBPUSH22L, GSBPUSH22L, or GBPUSH22C
! push3 => igbpush23 push particles with magnetic field and 3 component
!          electric field, with various interpolations and optimizations
!          calls GBPUSH23, GSBPUSH23, GBPUSH23L, GSBPUSH23L, GBPUSH23C,
!          GBPUSH22, GSBPUSH22, GBPUSH22L, GSBPUSH22L, or GBPUSH22C
! retard => iretard2 retard particle position a half time-step.
!           calls RETARD2
! djpostgl => idjpost2gl deposits current density, using gridless method
!             calls DJPOST2GL
! push3gl => ibpush23gl push particles with 2-1/2d electromagnetic
!            fields, using gridless method.
!            calls BPUSH23GL
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: november 28, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: BPUSH2, BPUSH2L
      public :: wtimer
      public :: djpost, push, push3, retard, djpostgl, push3gl
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine DJPOST2(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv)
         implicit none
         integer :: nop, idimp, nx, ny, nxv
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: cux, cuy, cuz
         end subroutine
      end interface
      interface
         subroutine GJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc&
     &)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine SJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxvy) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipb&
     &c)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine DJPOST2L(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv)
         implicit none
         integer :: nop, idimp, nx, ny, nxv
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: cux, cuy, cuz
         end subroutine
      end interface
      interface
         subroutine GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipb&
     &c)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine SJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxvy) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ip&
     &bc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipb&
     &c)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST22X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ip&
     &bc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ip&
     &bc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GSJPOST22XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,i&
     &pbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GJPOST2C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine GJPOST22C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine BPUSH2(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny&
     &,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         end subroutine
      end interface
      interface
         subroutine GBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,n&
     &xv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine BPUSH2L(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         end subroutine
      end interface
      interface
         subroutine GBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine BPUSH2CQ(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         end subroutine
      end interface
      interface
         subroutine GBPUSH2CQ(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv&
     &,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine BPUSH2CL(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy, bx, by, bz
         end subroutine
      end interface
      interface
         subroutine GBPUSH2CL(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv&
     &,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy
         real, dimension(3,nxyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,n&
     &xv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy
         real, dimension(3,nxv,nyv) :: bxy
         end subroutine
      end interface
      interface
         subroutine GBPUSH22C(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,&
     &nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine RETARD2(part,dtc,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc
         real :: dtc
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine DJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxvh,ny&
     &v,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: cu
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine BPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,&
     &nx,ny,nxvh,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt, dtc, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: fxy, bxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface djpost
         module procedure igjpost2
      end interface
!
      interface djpostgl
         module procedure idjpost2gl
      end interface
!
      interface push
         module procedure igbpush2
!        module procedure igbpush2cq
      end interface
!
      interface push3
         module procedure igbpush23
      end interface
!
      interface push3gl
         module procedure ibpush23gl
      end interface
!
      interface retard
         module procedure iretard2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igjpost2(part,cu,nop,qm,dt,tdjpost,nx,ny,ipbc,inorde&
     &r,djopt)
! deposit current
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, djopt
         real :: qm, dt, tdjpost
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
                  call GSJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc)
               else if (opt==VECTOR) then
                  call GSJPOST22XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxy&
     &v,ipbc)
               else
                  call GJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,i&
     &pbc)
               endif
            else if (order==CUBIC) then
               call GJPOST22C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &)
            else
               if (opt==LOOKAHEAD) then
                  call GSJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,&
     &ipbc)
               else if (opt==VECTOR) then
                  call GSJPOST22X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc)
               else
                  call GJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ip&
     &bc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,&
     &ipbc)
               else if (opt==VECTOR) then
                  call GSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv&
     &,ipbc)
               else
                  call GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ip&
     &bc)
               endif
            else if (order==CUBIC) then
               call GJPOST2C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,i&
     &pbc)
               else if (opt==VECTOR) then
                  call GSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,&
     &ipbc)
               else
                  call GJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipb&
     &c)
               endif
            endif
         end select
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine igjpost2
!
         subroutine igbpush2(part,fxy,bxy,nop,qbm,dt,dtc,ek,tpush,nx,ny,&
     &ipbc,inorder,popt)
! push particles with 2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
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
                  call GSBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
                  call GBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igbpush2
!
         subroutine igbpush2cq(part,fxy,bxy,nop,qbm,dt,ek,tpush,nx,ny,ip&
     &bc,inorder)
! push particles with 2d electromagnetic fields,
! with correction to Boris Mover
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tp,ltime,-1)
         if (order==LINEAR) then
            call GBPUSH2CL(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
         else
            call GBPUSH2CQ(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igbpush2cq
!
         subroutine igbpush23(part,fxy,bxy,nop,qbm,dt,dtc,ek,tpush,nx,ny&
     &,ipbc,inorder,popt)
! push particles with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ek, tpush
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
                  call GSBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH22L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GBPUSH22C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)

            else
               if (opt==LOOKAHEAD) then
                  call GSBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH22(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
               endif
            else if (order==CUBIC) then
               call GBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
            else
               if (opt==LOOKAHEAD) then
                  call GSBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc)
               else
                  call GBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igbpush23
!
         subroutine iretard2(part,nop,dtc,nx,ny,ipbc,ndim)
! retards particle positions half time-step
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: ndim
         real :: dtc
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call RETARD2(part,dtc,idimp,nop,nx,ny,ipbc)
         end subroutine iretard2
!
         subroutine idjpost2gl(part,cu,nop,qm,dt,nx,ny,ipbc,tdjpost)
! deposit current using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qm, dt, tdjpost
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
            call DJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxvh,nyv,i&
     &pbc)
         end select
! record time
         call wtimer(tj,ltime)
         tdjpost = tdjpost + tj
         end subroutine idjpost2gl
!
         subroutine ibpush23gl(part,fxy,bxy,nop,qbm,dt,dtc,ek,nx,ny,ipbc&
     &,tpush)
! push particles with 2-1/2d electromagnetic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, dtc, ek, tpush
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
            call BPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,nx,&
     &ny,nxvh,nyv,ipbc)
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ibpush23gl
!
      end module bpush2d
