!-----------------------------------------------------------------------
!
      module mrdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mrdpush2lib.f
! mrdpush2mod.f contains multi-tasking interface procedures to process
!               relativistic particles with darwin electric and magnetic
!               fields:
!               defines module mrdpush2d
! rdmjpost => imgrmjpost2 deposits relativistic momentum flux, with
!             various interpolations and optimizations.
!             calls MGRMJPOST2, MGSRMJPOST2, MGRMJPOST2L, MGSRMJPOST2L,
!             MGRMJPOST2C, MGRMJPOST22, MGSRMJPOST22, MGRMJPOST22L,
!             MGSRMJPOST22L, or MGRMJPOST22C
! rdcjpost => imgrdcjpost2 deposits relativistic momentum flux,
!             acceleration density, and current density, with various
!             interpolations and optimizations.
!             calls MGRDCJPOST2, MGSRDCJPOST2, MGRDCJPOST2L,
!             MGSRDCJPOST2L,MGRDCJPOST2C, MGRDCJPOST22, MGSRDCJPOST22,
!             MGRDCJPOST22L, MGSRDCJPOST22L, or MGRDCJPOST22C
! rdmjpostgl => imgrmjpost2gl deposits relativistic momentum flux, using
!               gridless method.
!               calls MGRMJPOST2GL
! rdcjpostgl => imgrdcjpost2gl deposits relativistic momentum flux,
!               acceleration density, and current density, using
!               gridless method.
!               calls MGRDCJPOST2GL
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 30, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use rdpush2d, only: wtimer
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, rdmjpost, rdcjpost, rdmjpostgl, rdcjpostgl
!
! buffer data for  momentum flux, acceleration, and current deposit
      real, dimension(:,:,:,:), allocatable :: cup, dcup, amup
      integer :: szbufc = 0, szbufd = 0, szbufa = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MGRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idt&
     &ask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv,  nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,i&
     &dtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,id&
     &task,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,&
     &idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,id&
     &task,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv,  nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         real, dimension(2,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,&
     &idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         real, dimension(2,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,i&
     &dtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         real, dimension(2,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup&
     &,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         real, dimension(2,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(3,nxv,nyv,nmt) :: cup, dcup
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         real, dimension(3,nxyv,nmt) :: cup, dcup
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(3,nxv,nyv,nmt) :: cup, dcup
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,i&
     &dimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         real, dimension(3,nxyv,nmt) :: cup, dcup
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         real, dimension(2,nxv,nyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         real, dimension(2,nxyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         real, dimension(2,nxv,nyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,i&
     &dimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         real, dimension(2,nxyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRMJPOST2C(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,id&
     &task,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv,  nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRMJPOST22C(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,i&
     &dtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv,  nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         real, dimension(2,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(3,nxv,nyv,nmt) :: cup, dcup
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         real, dimension(2,nxv,nyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxv&
     &h,nyv,amup,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nmt, ierr
         real :: qm, ci
         real, dimension(idimp,nop) :: part
         real, dimension(4,2*nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         real, dimension(4,2*nxvh,nyv,nmt) :: amup
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt&
     &,ci,idimp,nop,nx,ny,nxvh,nyv,cup,dcup,amup,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,2*nxvh,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,2*nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         real, dimension(3,2*nxvh,nyv,nmt) :: cup, dcup
         real, dimension(4,2*nxvh,nyv,nmt) :: amup
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rdmjpost
         module procedure imgrmjpost2
      end interface
!
      interface rdmjpostgl
         module procedure imgrmjpost2gl
      end interface
!
      interface rdcjpost
         module procedure imgrdcjpost2
      end interface
!
      interface rdcjpostgl
         module procedure imgrdcjpost2gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imgrmjpost2(part,amu,nop,qm,ci,tdcjpost,inorder,djop&
     &t)
! multi-tasking momentum flux deposit
! for relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, opt, ierr
         integer :: nnxyv
         real :: tdc
!        real, dimension(size(amu,1),size(amu,2),size(amu,3),ntasks) :: &
!    &amup
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(amu,2); nyv = size(amu,3)
         nxyv = nxv*nyv; nnxyv = size(amu,1)*nxyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of momentum flux buffer has changed
         if (szbufa < nnxyv) then
            if (szbufa /= 0) deallocate(amup)
! allocate buffer
            allocate(amup(size(amu,1),nxv,nyv,ntasks))
            szbufa = nnxyv
         endif
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nxyv,a&
     &mup,idtask,nmt,ierr)
               else
                  call MGRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nyv,amu&
     &p,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRMJPOST22C(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,i&
     &dtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nxyv,am&
     &up,idtask,nmt,ierr)
               else
                  call MGRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nyv,amup&
     &,idtask,nmt,ierr)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nxyv,am&
     &up,idtask,nmt,ierr)
               else
                  call MGRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nyv,amup&
     &,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRMJPOST2C(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,id&
     &task,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nxyv,amu&
     &p,idtask,nmt,ierr)
               else
                  call MGRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,&
     &idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgrmjpost2
!
         subroutine imgrdcjpost2(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,c&
     &i,tdcjpost,inorder,djopt)
! multi-tasking momentum flux, acceleration, and current deposit
! for relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy, cu, dcu, amu
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, opt, ierr
         integer :: nnxyv
         real :: tdc
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
!        real, dimension(size(dcu,1),size(dcu,2),size(dcu,3),ntasks) :: &
!    &dcup
!        real, dimension(size(amu,1),size(amu,2),size(amu,3),ntasks) :: &
!    &amup
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fxy,2); nyv = size(fxy,3)
         nxyv = nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of current buffer has changed
         nnxyv = size(cu,1)*nxyv
         if (szbufc < nnxyv) then
            if (szbufc /= 0) deallocate(cup)
            allocate(cup(size(cu,1),nxv,nyv,ntasks))
            szbufc = nnxyv
         endif
! check if size of acceleration buffer has changed
         nnxyv = size(dcu,1)*nxyv
         if (szbufd < nnxyv) then
            if (szbufd /= 0) deallocate(dcup)
            allocate(dcup(size(dcu,1),nxv,nyv,ntasks))
            szbufd = nnxyv
         endif
! check if size of momentum flux buffer has changed
         nnxyv = size(amu,1)*nxyv
         if (szbufa < nnxyv) then
            if (szbufa /= 0) deallocate(amup)
            allocate(amup(size(amu,1),nxv,nyv,ntasks))
            szbufa = nnxyv
         endif
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,&
     &ci,idimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGRDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,c&
     &i,idimp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRDCJPOST22C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,i&
     &dimp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,c&
     &i,idimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGRDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci&
     &,idimp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,c&
     &i,idimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci&
     &,idimp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGRDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci&
     &,idimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,&
     &idimp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgrdcjpost2
!
         subroutine imgrmjpost2gl(part,amu,nop,qm,ci,nx,ny,tdcjpost)
! multi-tasking deposit momentum flux using gridless method
! for relativistic particles
         implicit none
         integer :: nop, nx, ny
         real :: qm, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: idimp, nxv, nxvh, nyv, nnxyv, nmt, ltime, ierr
         real :: tdc
         complex, dimension(size(amu,2)/2) :: sctx
!        real, dimension(size(amu,1),size(amu,2),size(amu,3),ntasks) :: &
!    &amup
         complex, dimension(size(amu,2)/2,ntasks) :: sctxp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(amu,2); nxvh = nxv/2; nyv = size(amu,3)
         nnxyv = size(amu,1)*nxv*nyv
         nmt = ntasks
! check if size of momentum flux buffer has changed
         if (szbufa < nnxyv) then
            if (szbufa /= 0) deallocate(amup)
! allocate buffer
            allocate(amup(size(amu,1),nxv,nyv,ntasks))
            szbufa = nnxyv
         endif
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(amu,1))
         case (4)
            call MGRMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxvh,n&
     &yv,amup,sctxp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgrmjpost2gl
!
         subroutine imgrdcjpost2gl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt&
     &,ci,nx,ny,tdcjpost)
! multi-tasking momentum flux, acceleration, and current deposit
! using gridless method for relativistic particles
         implicit none
         integer :: nop, nx, ny
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy, cu, dcu, amu
! local data
         integer :: idimp, nxv, nxvh, nyv, nxyv, nnxyv, nmt, ltime, ierr
         real :: tdc
         complex, dimension(size(amu,2)/2) :: sctx
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
!        real, dimension(size(dcu,1),size(dcu,2),size(dcu,3),ntasks) :: &
!    &dcup
!        real, dimension(size(amu,1),size(amu,2),size(amu,3),ntasks) :: &
!    &amup
         complex, dimension(size(amu,2)/2,ntasks) :: sctxp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(fxy,2); nxvh = nxv/2; nyv = size(fxy,3)
         nxyv = nxv*nyv
         nmt = ntasks
! check if size of current buffer has changed
         nnxyv = size(cu,1)*nxyv
         if (szbufc < nnxyv) then
            if (szbufc /= 0) deallocate(cup)
            allocate(cup(size(cu,1),nxv,nyv,ntasks))
            szbufc = nnxyv
         endif
! check if size of acceleration buffer has changed
         nnxyv = size(dcu,1)*nxyv
         if (szbufd < nnxyv) then
            if (szbufd /= 0) deallocate(dcup)
            allocate(dcup(size(dcu,1),nxv,nyv,ntasks))
            szbufd = nnxyv
         endif
! check if size of momentum flux buffer has changed
         nnxyv = size(amu,1)*nxyv
         if (szbufa < nnxyv) then
            if (szbufa /= 0) deallocate(amup)
            allocate(amup(size(amu,1),nxv,nyv,ntasks))
            szbufa = nnxyv
         endif
! initialize timer
         call wtimer(tdc,ltime,-1)
         select case(size(bxy,1))
         case (3)
            call MGRDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,ci&
     &,idimp,nop,nx,ny,nxvh,nyv,cup,dcup,amup,sctxp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgrdcjpost2gl
!
      end module mrdpush2d
