!-----------------------------------------------------------------------
!
      module mdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mdpush2lib.f
! mdpush2mod.f contains multi-tasking interface procedures to process
!              particles with darwin electric and magnetic fields:
!              defines module mdpush2d
! dmjpost => imgmjpost2 deposits momentum flux, with various
!            interpolations and optimizations.
!            calls MGMJPOST2, MGSMJPOST2, MGMJPOST2L, MGSMJPOST2L,
!            MGMJPOST2C, MGMJPOST22, MGSMJPOST22, MGMJPOST22L, 
!            MGSMJPOST22L, or MGMJPOST22C
! dcjpost => imgdcjpost2 deposits momentum flux, acceleration density,
!            and current density, with various interpolations and
!            optimizations.
!            calls MGDCJPOST2, MGSDCJPOST2, MGDCJPOST2L, MGSDCJPOST2L,
!            MGDCJPOST2C, MGDCJPOST22, MGSDCJPOST22, MGDCJPOST22L,
!            MGSDCJPOST22L, or MGDCJPOST22C
! dmjpostgl => imgmjpost2gl deposits momentum flux, using gridless method
!              calls MGMJPOST2GL
! dcjpostgl => imgdcjpost2gl deposits momentum flux, acceleration
!              density, and current density, using gridless method.
!              calls MGDCJPOST2GL
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 30, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use dpush2d, only: wtimer
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: dmjpost, dcjpost, dmjpostgl, dcjpostgl
!
! buffer data for  momentum flux, acceleration, and current deposit
      real, dimension(:,:,:,:), allocatable :: cup, dcup, amup
      integer :: szbufc = 0, szbufd = 0, szbufa = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MGMJPOST2(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,&
     &nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv,  nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSMJPOST2(part,amu,qm,nop,idimp,nxv,nxyv,amup,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGMJPOST2L(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask&
     &,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSMJPOST2L(part,amu,qm,nop,idimp,nxv,nxyv,amup,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxyv) :: amu
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGMJPOST22(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask&
     &,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv,  nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         real, dimension(2,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSMJPOST22(part,amu,qm,nop,idimp,nxv,nxyv,amup,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         real, dimension(2,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGMJPOST22L(part,amu,qm,nop,idimp,nxv,nyv,amup,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         real, dimension(2,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSMJPOST22L(part,amu,qm,nop,idimp,nxv,nxyv,amup,idt&
     &ask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: amu
         real, dimension(2,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGMJPOST2C(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask&
     &,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGMJPOST22C(part,amu,qm,nop,idimp,nxv,nyv,amup,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: amu
         real, dimension(2,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(3,nxv,nyv,nmt) :: cup, dcup
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         real, dimension(3,nxyv,nmt) :: cup, dcup
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(3,nxv,nyv,nmt) :: cup, dcup
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp&
     &,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxyv) :: amu
         real, dimension(3,nxyv,nmt) :: cup, dcup
         real, dimension(4,nxyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         real, dimension(2,nxv,nyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         real, dimension(2,nxyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         real, dimension(2,nxv,nyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp&
     &,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, cu, dcu, amu
         real, dimension(nxyv) :: bz
         real, dimension(2,nxyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, bxy, cu, dcu
         real, dimension(4,nxv,nyv) :: amu
         real, dimension(3,nxv,nyv,nmt) :: cup, dcup
         real, dimension(4,nxv,nyv,nmt) :: amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
       interface
         subroutine MGDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu, dcu, amu
         real, dimension(nxv,nyv) :: bz
         real, dimension(2,nxv,nyv,nmt) :: cup, dcup, amup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGMJPOST2GL(part,amu,sctx,qm,nop,idimp,nx,ny,nxvh,ny&
     &v,amup,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(4,2*nxvh,nyv) :: amu
         complex, dimension(nxvh) :: sctx
         real, dimension(4,2*nxvh,nyv,nmt) :: amup
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,&
     &idimp,nop,nx,ny,nxvh,nyv,cup,dcup,amup,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nmt, ierr
         real :: qm, qbm, dt
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
      interface dmjpost
         module procedure imgmjpost2
      end interface
!
      interface dmjpostgl
         module procedure imgmjpost2gl
      end interface
!
      interface dcjpost
         module procedure imgdcjpost2
      end interface
!
      interface dcjpostgl
         module procedure imgdcjpost2gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imgmjpost2(part,amu,nop,qm,tdcjpost,inorder,djopt)
! multi-tasking momentum flux deposit
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, tdcjpost
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
                  call MGSMJPOST22L(part,amu,qm,nop,idimp,nxv,nxyv,amup,&
     &idtask,nmt,ierr)
               else
                  call MGMJPOST22L(part,amu,qm,nop,idimp,nxv,nyv,amup,id&
     &task,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGMJPOST22C(part,amu,qm,nop,idimp,nxv,nyv,amup,idtas&
     &k,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSMJPOST22(part,amu,qm,nop,idimp,nxv,nxyv,amup,i&
     &dtask,nmt,ierr)
               else
                  call MGMJPOST22(part,amu,qm,nop,idimp,nxv,nyv,amup,idt&
     &ask,nmt,ierr)
               endif
            endif
         case (4)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSMJPOST2L(part,amu,qm,nop,idimp,nxv,nxyv,amup,i&
     &dtask,nmt,ierr)
               else
                  call MGMJPOST2L(part,amu,qm,nop,idimp,nxv,nyv,amup,idt&
     &ask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGMJPOST2C(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask&
     &,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSMJPOST2(part,amu,qm,nop,idimp,nxv,nxyv,amup,id&
     &task,nmt,ierr)
               else
                  call MGMJPOST2(part,amu,qm,nop,idimp,nxv,nyv,amup,idta&
     &sk,nmt,ierr)
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
         end subroutine imgmjpost2
!
         subroutine imgdcjpost2(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,td&
     &cjpost,inorder,djopt)
! multi-tasking momentum flux, acceleration, and current deposit
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, tdcjpost
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
                  call MGSDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,i&
     &dimp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGDCJPOST22L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,id&
     &imp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGDCJPOST22C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp&
     &,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,id&
     &imp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGDCJPOST22(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idi&
     &mp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,id&
     &imp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idi&
     &mp,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
               endif
            else if (order==CUBIC) then
               call MGDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
            else
               if (opt==LOOKAHEAD) then
                  call MGSDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idi&
     &mp,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
               else
                  call MGDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idim&
     &p,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
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
         end subroutine imgdcjpost2
!
         subroutine imgmjpost2gl(part,amu,nop,qm,nx,ny,tdcjpost)
! multi-tasking deposit momentum flux using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, tdcjpost
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
            call MGMJPOST2GL(part,amu,sctx,qm,nop,idimp,nx,ny,nxvh,nyv,a&
     &mup,sctxp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgmjpost2gl
!
         subroutine imgdcjpost2gl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,&
     &nx,ny,tdcjpost)
! multi-tasking momentum flux, acceleration, and current deposit
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, qbm, dt, tdcjpost
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
            call MGDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,idi&
     &mp,nop,nx,ny,nxvh,nyv,cup,dcup,amup,sctxp,idtask,nmt,ierr)
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgdcjpost2gl
!
      end module mdpush2d
