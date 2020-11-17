!-----------------------------------------------------------------------
!
      module dpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library dpush2lib.f
! dpush2mod.f contains interface procedures to process particles with
!             darwin electric and magnetic fields:
!             defines module dpush2d
! dmjpost => igmjpost2 deposits momentum flux, with various
!            interpolations and optimizations.
!            calls GMJPOST2, GSMJPOST2, GMJPOST2L, GSMJPOST2L,
!            GMJPOST2C, GMJPOST22, GSMJPOST22, GMJPOST22L, GSMJPOST22L,
!            or GMJPOST22C
! dcjpost => igdcjpost2 deposits momentum flux, acceleration density,
!            and current density, with various interpolations and
!            optimizations.
!            calls GDCJPOST2, GSDCJPOST2, GDCJPOST2L, GSDCJPOST2L,
!            GDCJPOST2C, GDCJPOST22, GSDCJPOST22, GDCJPOST22L,
!            GSDCJPOST22L, or GDCJPOST22C
! dmjpostgl => igmjpost2gl deposits momentum flux, using gridless method
!              calls GMJPOST2GL
! dcjpostgl => igdcjpost2gl deposits momentum flux, acceleration
!              density, and current density, using gridless method.
!              calls GDCJPOST2GL
! written by viktor k. decyk, ucla
! copyright 2006, regents of the university of california
! update: august 14, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, dmjpost, dcjpost, dmjpostgl, dcjpostgl
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GMJPOST2(part,amu,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSMJPOST2(part,amu,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,no&
     &p,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GMJPOST22(part,amu,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSMJPOST22(part,amu,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,no&
     &p,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GMJPOST2L(part,amu,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSMJPOST2L(part,amu,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GMJPOST22L(part,amu,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GSMJPOST22L(part,amu,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GSDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GMJPOST2C(part,amu,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GMJPOST22C(part,amu,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         end subroutine
      end interface
      interface
         subroutine GDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         end subroutine
      end interface
      interface
         subroutine GMJPOST2GL(part,amu,sctx,qm,nop,idimp,nx,ny,nxvh,nyv&
     &)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,2*nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,i&
     &dimp,nop,nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,2*nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dmjpost
         module procedure igmjpost2
      end interface
!
      interface dmjpostgl
         module procedure igmjpost2gl
      end interface
!
      interface dcjpost
         module procedure igdcjpost2
      end interface
!
      interface dcjpostgl
         module procedure igdcjpost2gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igmjpost2(part,amu,nop,qm,tdcjpost,inorder,djopt)
! deposit momentum flux with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, tdcjpost
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
                  call GSMJPOST22L(part,amu,qm,nop,idimp,nxv,nxyv)
               else
                  call GMJPOST22L(part,amu,qm,nop,idimp,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GMJPOST22C(part,amu,qm,nop,idimp,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSMJPOST22(part,amu,qm,nop,idimp,nxv,nxyv)
               else
                  call GMJPOST22(part,amu,qm,nop,idimp,nxv,nyv)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSMJPOST2L(part,amu,qm,nop,idimp,nxv,nxyv)
               else
                  call GMJPOST2L(part,amu,qm,nop,idimp,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GMJPOST2C(part,amu,qm,nop,idimp,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSMJPOST2(part,amu,qm,nop,idimp,nxv,nxyv)
               else
                  call GMJPOST2(part,amu,qm,nop,idimp,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igmjpost2
!
         subroutine igdcjpost2(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,tdc&
     &jpost,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, tdcjpost
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
                  call GSDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,id&
     &imp,nop,nxv,nxyv)
               else
                  call GDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idi&
     &mp,nop,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GDCJPOST22C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idi&
     &mp,nop,nxv,nxyv)
               else
                  call GDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idim&
     &p,nop,nxv,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idi&
     &mp,nop,nxv,nxyv)
               else
                  call GDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idim&
     &p,nop,nxv,nyv)
               endif
            else if (order==CUBIC) then
               call GDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv)
            else
               if (opt==LOOKAHEAD) then
                  call GSDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idim&
     &p,nop,nxv,nxyv)
               else
                  call GDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp&
     &,nop,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igdcjpost2
!
         subroutine igmjpost2gl(part,amu,nop,qm,nx,ny,tdcjpost)
! deposit momentum flux using gridless method
! with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny
         real :: qm, tdcjpost
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
            call GMJPOST2GL(part,amu,sctx,qm,nop,idimp,nx,ny,nxvh,nyv)
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igmjpost2gl
!
         subroutine igdcjpost2gl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,n&
     &x,ny,tdcjpost)
! deposit momentum flux, acceleration density, and curent density
! with 2-1/2d electromagnetic fields, using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, qbm, dt, tdcjpost
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
            call GDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,idim&
     &p,nop,nx,ny,nxvh,nyv)
         end select
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igdcjpost2gl
!
      end module dpush2d
