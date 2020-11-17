!-----------------------------------------------------------------------
!
      module mpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mpush2lib.f
! mpush2mod.f contains multi-tasking interface procedures to process
!             particles:
!             defines module mpush2d
! dpost => imgpost2 deposits charge density, with various interpolations
!          and optimizations.
!          calls MGPOST2, MGSPOST2, MGSPOST2X, MGPOST2L, MGSPOST2L,
!          MGSPOST2XL, or MGPOST2C
! push => imgpush2 push particles, with various interpolations and
!         optimizations.
!         calls MGPUSH2, MGSPUSH2, MGPUSH2L, MGSPUSH2L, or MGPUSH2C
! sortp => imsortp2y sorts particles by y grid using memory-conserving
!          bin sort, with various interpolations.
!          calls MSORTP2Y, or MSORTP2YL
! sortp => imdsortp2y sorts particles by y grid using optimized bin sort
!          with various interpolations.
!          calls MDSORTP2Y, or MDSORTP2YL
! pushzf => impush2zf, push particles with no forces.
!           calls MPUSH2ZF
! dpostgl => imdpost2gl deposits charge density using gridless method.
!            calls MDPOST2GL
! dpostglx => imdpost2glx deposits charge density using gridless method.
!             calls MDPOST2GL or MDPOST2GLX
! pushgl => impush2gl push particles using optimzed gridless method.
!           calls MPUSH2GL
! pushglx => impush2glx push particles using optimized gridless method.
!            calls MPUSH2GL or MPUSH2GLX
! gcjpost => imgcjpost2 deposits time-centered current density with
!            2d electrostatic fields.
!            calls MGCJPOST2, MGCJPOST2L, or MGCJPOST2C
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: january 8, 2010
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR
      use push2d, only: wtimer, rmove, initmomt2, premoment2, primoment2
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: dpost, push, sortp, rmove, dpostgl, pushzf, pushgl
      public :: dpostglx, pushglx, gcjpost
      public :: initmomt2, premoment2, primoment2
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MDPOST2(part,q,qm,nop,idimp,nx,ny,nxv,qp,idtask,nmt,&
     &ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: q
         real, dimension(nxv,ny,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPOST2(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,ie&
     &rr)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,nyv) :: q
         real, dimension(nxv,nyv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST2(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt,&
     &ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         real, dimension(nxyv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSPOST2X(part,q,qm,nop,idimp,nx,ny,nxv,nxvy,qp,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxvy) :: q
         real, dimension(nxvy,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST2X(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt&
     &,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         real, dimension(nxyv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDPOST2L(part,q,qm,nop,idimp,nx,ny,nxv,qp,idtask,nmt&
     &,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: q
         real, dimension(nxv,ny,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPOST2L(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,i&
     &err)
         implicit none
         integer :: nop, idimp, nxv, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,nyv) :: q
         real, dimension(nxv,nyv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST2L(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt&
     &,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         real, dimension(nxyv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSPOST2XL(part,q,qm,nop,idimp,nx,ny,nxv,nxvy,qp,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxvy) :: q
         real, dimension(nxvy,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPOST2XL(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nm&
     &t,ierr)
         implicit none
         integer :: nop, idimp, nxv, nxyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         real, dimension(nxyv,nmt) :: qp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPUSH2(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ekp,&
     &idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPUSH2L(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ekp&
     &,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,ny) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real  :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSORTP2Y(part,pt,ip,npic,idimp,nop,ny1,npicp,idtask,&
     &nmt,ierr)
         implicit none
         integer :: idimp, nop, ny1, nmt, ierr
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(ny1) :: npic
         integer, dimension(ny1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MSORTP2YL(part,pt,ip,npic,idimp,nop,ny1,npicp,idtask&
     &,nmt,ierr)
         implicit none
         integer :: idimp, nop, ny1, nmt, ierr
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(ny1) :: npic
         integer, dimension(ny1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDSORTP2Y(parta,partb,npic,idimp,nop,ny1,npicp,idtas&
     &k,nmt,ierr)
         implicit none
         integer :: idimp, nop, ny1, nmt, ierr
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(ny1) :: npic
         integer, dimension(ny1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDSORTP2YL(parta,partb,npic,idimp,nop,ny1,npicp,idta&
     &sk,nmt,ierr)
         implicit none
         integer :: idimp, nop, ny1, nmt, ierr
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(ny1) :: npic
         integer, dimension(ny1,nmt) :: npicp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPUSH2ZF(part,dt,ek,idimp,nop,nx,ny,ipbc,ekp,idtask,&
     &nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc, nmt, ierr
         real  :: dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxvh,nyv,qp&
     &,sctxp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh,nyv) :: q
         complex, dimension(nxvh) :: sctx
         real, dimension(2*nxvh,nyv,nmt) :: qp
         complex, dimension(nxvh,nmt) :: sctxp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MDPOST2GLX(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,&
     &nnp,qp,sctxpp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nnp, nmt, ierr
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh,nyv) :: q
         complex, dimension(nxvh,nnp) :: sctxp
         real, dimension(2*nxvh,nyv,nmt) :: qp
         complex, dimension(nxvh,nnp,nmt) :: sctxpp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxv&
     &h,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh) :: sctx
         complex, dimension(nxvh,nmt) :: sctxp
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxvh,nyv,nnp,ipbc,sctxpp,exypp,ekp,idtask,nmt,ierr)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, nnp, ipbc, nmt, ierr
         real :: qbm, dt, ek
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
         subroutine MGCJPOST2(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cu&
     &p,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv, nmt, ierr
         real  :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGCJPOST2L(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,c&
     &up,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv, nmt, ierr
         real  :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGCJPOST2C(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,c&
     &up,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv, nmt, ierr
         real  :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, cu
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure imgpost2
      end interface
!
      interface dpostgl
         module procedure imdpost2gl
      end interface
!
      interface dpostglx
         module procedure imdpost2glx
      end interface
!
      interface push
         module procedure imgpush2
      end interface
!
      interface pushzf
         module procedure impush2zf
      end interface
!
      interface pushgl
         module procedure impush2gl
      end interface
!
      interface pushglx
         module procedure impush2glx
      end interface
!
      interface sortp
         module procedure imsortp2y
         module procedure imdsortp2y
      end interface
!
      interface gcjpost
         module procedure imgcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imgpost2(part,q,nop,qm,tdpost,inorder,dopt)
! multi-tasking charge deposit
         implicit none
         integer :: nop
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, ltime, order, opt, ierr
         real :: td
!        real, dimension(size(q,1),size(q,2),ntasks) :: qp
         real, dimension(:,:,:), allocatable, save :: qp
         integer, save :: szbuf = 0
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(q,1); nyv = size(q,2)
         nxyv = nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! check if size of buffer has changed
         if (szbuf < nxyv) then
            if (szbuf /= 0) deallocate(qp)
! allocate buffer
            allocate(qp(nxv,nyv,ntasks))
            szbuf = nxyv
         endif
! initialize timer
         call wtimer(td,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call MGSPOST2L(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt&
     &,ierr)
            else if (opt==VECTOR) then
               call MGSPOST2XL(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nm&
     &t,ierr)
            else
               call MGPOST2L(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,i&
     &err)
            endif
         else if (order==CUBIC) then
            call MGPOST2C(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,ierr&
     &)
         else
            if (opt==LOOKAHEAD) then
               call MGSPOST2(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt,&
     &ierr)
            else if (opt==VECTOR) then
               call MGSPOST2X(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt&
     &,ierr)
            else
               call MGPOST2(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,ie&
     &rr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine imgpost2
!
         subroutine imgpush2(part,fxy,nop,qbm,dt,ek,tpush,nx,ny,ipbc,ino&
     &rder,popt)
! multi-tasking particle push with 2d electrostatic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
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
               call MGSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc,ekp,idtask,nmt,ierr)
            else
               call MGPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc,ekp,idtask,nmt,ierr)
            endif
         else if (order==CUBIC) then
            call MGPUSH2C(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipb&
     &c,ekp,idtask,nmt,ierr)
         else
            if (opt==LOOKAHEAD) then
               call MGSPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc,ekp,idtask,nmt,ierr)
            else
               call MGPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc,ekp,idtask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgpush2
!
         subroutine imsortp2y(part,pt,ip,nop,npic,tsort,inorder)
! multi-tasking particle sort by y grid using memory-conserving bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: pt
         integer, dimension(:), pointer :: ip, npic
! local data
         integer, dimension(size(npic),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         integer :: idimp, ny1, nmt, ltime, order, ierr
         real :: ts
         idimp = size(part,1)
         ny1 = size(npic)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call MSORTP2YL(part,pt,ip,npic,idimp,nop,ny1,npicp,idtask,nm&
     &t,ierr)
         else
            call MSORTP2Y(part,pt,ip,npic,idimp,nop,ny1,npicp,idtask,nmt&
     &,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine imsortp2y
!
         subroutine imdsortp2y(parta,partb,nop,npic,tsort,inorder)
! multi-tasking particle sort by y grid using optimized bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npic
         real, dimension(:,:), pointer :: part
! local data
         integer, dimension(size(npic),ntasks) :: npicp
         integer, dimension(ntasks) :: idtask
         integer :: idimp, ny1, nmt, ltime, order, ierr
         real :: ts
         idimp = size(parta,1)
         ny1 = size(npic)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call MDSORTP2YL(parta,partb,npic,idimp,nop,ny1,npicp,idtask,&
     &nmt,ierr)
         else
            call MDSORTP2Y(parta,partb,npic,idimp,nop,ny1,npicp,idtask,n&
     &mt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine imdsortp2y
!
         subroutine impush2zf(part,nop,dt,ek,tpush,nx,ny,ipbc)
! multi-tasking particle push with no forces
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: dt, ek, tpush
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
         call MPUSH2ZF(part,dt,ek,idimp,nop,nx,ny,ipbc,ekp,idtask,nmt,ie&
     &rr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine impush2zf
!
         subroutine imdpost2gl(part,q,nop,qm,nx,ny,tdpost)
! multi-tasking deposit charge using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         integer :: idimp, nxvh, nyv, nmt, ltime, ierr
         real :: td
         complex, dimension(size(q,1)/2) :: sctx
         real, dimension(size(q,1),size(q,2),ntasks) :: qp
         complex, dimension(size(q,1)/2,ntasks) :: sctxp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxvh = size(q,1)/2; nyv = size(q,2)
         nmt = ntasks
! initialize timer
         call wtimer(td,ltime,-1)
         call MDPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxvh,nyv,qp,sctxp&
     &,idtask,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine imdpost2gl
!
         subroutine imdpost2glx(part,q,nop,qm,nx,ny,tdpost,dopt)
! multi-tasking deposit charge using gridless method
! optimized version
         implicit none
         integer :: nop, nx, ny
         real :: qm, tdpost
         integer, optional :: dopt
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         integer, parameter :: nnp = 128
         integer :: idimp, nxvh, nyv, nmt, opt, ltime, ierr
         real :: td
         complex, dimension(size(q,1)/2,nnp) :: sctxp
         real, dimension(size(q,1),size(q,2),ntasks) :: qp
         complex, dimension(size(q,1)/2,nnp,ntasks) :: sctxpp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxvh = size(q,1)/2; nyv = size(q,2)
         nmt = ntasks
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         call wtimer(td,ltime,-1)
         if (opt==LOOKAHEAD) then
            call MDPOST2GLX(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,nnp&
     &,qp,sctxpp,idtask,nmt,ierr)
         else
            call MDPOST2GL(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,qp,s&
     &ctxpp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine imdpost2glx
!
         subroutine impush2gl(part,fxy,nop,qbm,dt,ek,nx,ny,ipbc,tpush)
! multi-tasking particle push with 2d electrostatic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ek, tpush
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
         call MPUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxvh,nyv,&
     &ipbc,sctxp,ekp,idtask,nmt,ierr)
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine impush2gl
!
         subroutine impush2glx(part,fxy,nop,qbm,dt,ek,nx,ny,ipbc,tpush,p&
     &opt)
! multi-tasking particle push with 2d electrostatic fields
! using gridless method, optimized version
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ek, tpush
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
            call MPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ek,idimp,nop,nx,ny&
     &,nxvh,nyv,nnp,ipbc,sctxpp,exypp,ekp,idtask,nmt,ierr)
         else
            call MPUSH2GL(part,fxy,sctxp,qbm,dt,ek,idimp,nop,nx,ny,nxvh,&
     &nyv,ipbc,sctxpp,ekp,idtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine impush2glx
!
         subroutine imgcjpost2(part,fxy,cu,nop,qm,qbm,dt,tdcjpost,inorde&
     &r)
! multi-tasking current density deposit with 2d electrostatic fields
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qm, qbm, dt, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, cu
! local data
         integer :: idimp, nxv, nyv, nmt, ltime, order, ierr, nnxyv
         real :: tdc
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
         real, dimension(:,:,:,:), allocatable, save :: cup
         integer, save :: szbuf = 0
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
            call MGCJPOST2L(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cup,&
     &idtask,nmt,ierr)
         else if (order==CUBIC) then
             call MGCJPOST2C(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cup&
     &,idtask,nmt,ierr)
         else
            call MGCJPOST2(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cup,i&
     &dtask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine imgcjpost2
!
      end module mpush2d
