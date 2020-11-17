!-----------------------------------------------------------------------
!
      module rdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library rdpush2lib.f
! rdpush2mod.f contains interface procedures to process relativistic
!              particles with darwin electric and magnetic fields:
!              defines module rdpush2d
! rdmjpost => igrmjpost2 deposits relativistic momentum flux, with
!             various interpolations and optimizations.
!             calls GRMJPOST2, GSRMJPOST2, GRMJPOST2L, GSRMJPOST2L,
!             GRMJPOST2C, GRMJPOST22, GSRMJPOST22, GRMJPOST22L,
!             GSRMJPOST22L, or GRMJPOST22C
! rdcjpost => igrdcjpost2 deposits relativistic momentum flux,
!             acceleration density, and current density, with various
!             interpolations and optimizations.
!             calls GRDCJPOST2, GSRDCJPOST2, GRDCJPOST2L, GSRDCJPOST2L,
!             GRDCJPOST2C, GRDCJPOST22, GSRDCJPOST22, GRDCJPOST22L,
!             GSRDCJPOST22L, or GRDCJPOST22C
! rdmjpostgl => igrmjpost2gl deposits relativistic momentum flux, using
!               gridless method.
!               calls GRMJPOST2GL
! rdcjpostgl => igrdcjpost2gl deposits relativistic momentum flux,
!               acceleration density, and current density, using
!               gridless method.
!               calls GRDCJPOST2GL
! written by viktor k. decyk, ucla
! copyright 2006, regents of the university of california
! update: august 14, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use diag2d, only: wtimer
      use rbpush2d, only: rdjpost, rpush3, icptov2
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, rdmjpost, rdcjpost, rdmjpostgl, rdcjpostgl
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GRMJPOST2C(part,amu,qm,ci,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRMJPOST22C(part,amu,qm,ci,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GRMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxvh&
     &,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GRDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,&
     &ci,idimp,nop,nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxvh,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdmjpost
         module procedure igrmjpost2
      end interface
!
      interface rdmjpostgl
         module procedure igrmjpost2gl
      end interface
!
      interface rdcjpost
         module procedure igrdcjpost2
      end interface
!
      interface rdcjpostgl
         module procedure igrdcjpost2gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igrmjpost2(part,amu,nop,qm,ci,tdcjpost,inorder,djopt&
     &)
! deposit momentum flux with 2-1/2d electromagnetic fields
! and with relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tdc
         idimp = size(part,1)
         nxv = size(amu,2); nyv = size(amu,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nxyv)
               else
                  call GRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GRMJPOST22C(part,amu,qm,ci,nop,idimp,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nxyv)
               else
                  call GRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nyv)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nxyv)
               else
                  call GRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GRMJPOST2C(part,amu,qm,ci,nop,idimp,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nxyv)
               else
                  call GRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igrmjpost2
!
         subroutine igrdcjpost2(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,ci&
     &,tdcjpost,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields and with relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy, cu, dcu, amu
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tdc
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,c&
     &i,idimp,nop,nxv,nxyv)
               else
                  call GRDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci&
     &,idimp,nop,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GRDCJPOST22C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSRDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci&
     &,idimp,nop,nxv,nxyv)
               else
                  call GRDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,&
     &idimp,nop,nxv,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci&
     &,idimp,nop,nxv,nxyv)
               else
                  call GRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,&
     &idimp,nop,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GRDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,&
     &idimp,nop,nxv,nxyv)
               else
                  call GRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,i&
     &dimp,nop,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igrdcjpost2
!
         subroutine igrmjpost2gl(part,amu,nop,qm,ci,nx,ny,tdcjpost)
! deposit momentum flux using gridless method
! with 2-1/2d electromagnetic fields and with relativistic particles
         implicit none
         integer :: nop, nx, ny
         real :: qm, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tdc
         complex, dimension(size(amu,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(amu,2)/2; nyv = size(amu,3)
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(amu,1))
         case (4)
            call GRMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxvh,ny&
     &v)
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igrmjpost2gl
!
         subroutine igrdcjpost2gl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,&
     &ci,nx,ny,tdcjpost)
! deposit momentum flux, acceleration density, and curent density
! with 2-1/2d electromagnetic fields, using gridless method
! and with relativistic particles
         implicit none
         integer :: nop, nx, ny
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy, cu, dcu, amu
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tdc
         complex, dimension(size(fxy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(bxy,1))
         case (3)
            call GRDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,ci,&
     &idimp,nop,nx,ny,nxvh,nyv)
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igrdcjpost2gl
!
      end module rdpush2d
