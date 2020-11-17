!-----------------------------------------------------------------------
!
      module push2d
!
! Fortran90 interface to 2d PIC Fortran77 library push2lib.f
! push2mod.f contains interface procedures to process particles:
!            defines module push2d
! dpost => igpost2 deposits charge density, with various interpolations
!          and optimizations.
!          calls GPOST2, GSPOST2, GSPOST2X, GPOST2L, GSPOST2L, GSPOST2XL
!          or GPOST2C
! push => igpush2 push particles, with various interpolations and
!         optimizations.
!         calls GPUSH2, GSPUSH2, GPUSH2L, GSPUSH2L, or GPUSH2C
! sortp => isortp2y sorts particles by y grid using memory-conserving
!          bin sort, with various interpolations.
!          calls SORTP2Y, or SORTP2YL
! sortp => idsortp2y sorts particles by y grid using optimized bin sort
!          with various interpolations.
!          calls DSORTP2Y, or DSORTP2YL
! rmove => irmove2, removes particles which would normally be reflected.
!          calls RMOVE2
! pushzf => ipush2zf, push particles with no forces.
!           calls PUSH2ZF
! dpostgl => idpost2gl deposits charge density using gridless method.
!            calls DPOST2GL
! dpostglx => idpost2gl deposits charge density using gridless method.
!             calls DPOST2GL or DPOST2GLX
! pushgl => ipush2gl push particles using optimzed gridless method.
!           calls PUSH2GL
! pushglx => ipush2glx push particles using optimized gridless method.
!            calls PUSH2GL or PUSH2GLX
! gcjpost => igcjpost2 deposits time-centered current density with
!            2d electrostatic fields.
!            calls GCJPOST2, GCJPOST2L, or GCJPOST2C
! initmomt2 calculates initial momentum, for 2 or 2-1/2d code.
! premoment2 print outs electron and field momentum, calculates total
!            momentum for 2 or 2-1/2d code.
! primoment2 prints out ion momentum, adds total momentum, for 2
!            or 2-1/2d code
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
      public :: DPOST2, SPOST2X, DPOST2L, SPOST2XL, PUSH2, PUSH2L
      public :: wtimer
      public :: dpost, push, sortp, sort2p, rmove, pushzf
      public :: dpostgl, pushgl, dpostglx, pushglx, gcjpost
      public :: initmomt2, premoment2, primoment2
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine DPOST2(part,q,qm,nop,idimp,nx,ny,nxv)
         implicit none
         integer :: nop, idimp, nx, ny, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: q
         end subroutine
      end interface
      interface
         subroutine GPOST2(part,q,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST2(part,q,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         end subroutine
      end interface
      interface
         subroutine SPOST2X(part,q,qm,nop,idimp,nx,ny,nxv,nxvy)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxvy) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST2X(part,q,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         end subroutine
      end interface
      interface
         subroutine DPOST2L(part,q,qm,nop,idimp,nx,ny,nxv)
         implicit none
         integer :: nop, idimp, nx, ny, nxv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: q
         end subroutine
      end interface
      interface
         subroutine GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST2L(part,q,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         end subroutine
      end interface
      interface
         subroutine SPOST2XL(part,q,qm,nop,idimp,nx,ny,nxv,nxvy)
         implicit none
         integer :: nop, idimp, nx, ny, nxv, nxvy
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxvy) :: q
         end subroutine
      end interface
      interface
         subroutine GSPOST2XL(part,q,qm,nop,idimp,nxv,nxyv)
         implicit none
         integer :: nop, idimp, nxv, nxyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxyv) :: q
         end subroutine
      end interface
      interface
         subroutine GPOST2C(part,q,qm,nop,idimp,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,nyv) :: q
         end subroutine
      end interface
      interface
         subroutine PUSH2(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy
         end subroutine
      end interface
      interface
         subroutine GPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ip&
     &bc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GSPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,&
     &ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine PUSH2L(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: fx, fy
         end subroutine
      end interface
      interface
         subroutine GPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine GPUSH2C(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy
         end subroutine
      end interface
      interface
         subroutine SORTP2Y(part,pt,ip,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(ny1) :: npic
         end subroutine
      end interface
      interface
         subroutine SORTP2YL(part,pt,ip,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(ny1) :: npic
         end subroutine
      end interface
      interface
         subroutine DSORTP2Y(parta,partb,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(ny1) :: npic
         end subroutine
      end interface
      interface
         subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: parta, partb
         integer, dimension(ny1) :: npic
         end subroutine
      end interface
      interface
         subroutine IPSORTP2Y(part,pt2,ip2,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: part
         real, dimension(idimp,2) :: pt2
         integer, dimension(2,nop) :: ip2
         integer, dimension(ny1) :: npic
         end subroutine
      end interface
      interface
         subroutine IPSORTP2YL(part,pt2,ip2,npic,idimp,nop,ny1)
         implicit none
         integer :: idimp, nop, ny1
         real, dimension(idimp,nop) :: part
         real, dimension(idimp,2) :: pt2
         integer, dimension(2,nop) :: ip2
         integer, dimension(ny1) :: npic
         end subroutine
      end interface
      interface
         subroutine SORTP2(part,pt,ip,npic,idimp,nop,nx1,nxy1)
         implicit none
         integer :: idimp, nop, nx1, nxy1
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(nxy1) :: npic
         end subroutine
      end interface
      interface
         subroutine SORTP2L(part,pt,ip,npic,idimp,nop,nx1,nxy1)
         implicit none
         integer :: idimp, nop, nx1, nxy1
         real, dimension(idimp,nop) :: part
         real, dimension(nop) :: pt
         integer, dimension(nop) :: ip
         integer, dimension(nxy1) :: npic
         end subroutine
      end interface
      interface
         subroutine IPSORTP2(part,pt2,ip2,npic,idimp,nop,nx1,nxy1)
         implicit none
         integer :: idimp, nop, nx1, nxy1
         real, dimension(idimp,nop) :: part
         real, dimension(idimp,2) :: pt2
         integer, dimension(2,nop) :: ip2
         integer, dimension(nxy1) :: npic
         end subroutine
      end interface
      interface
         subroutine IPSORTP2L(part,pt2,ip2,npic,idimp,nop,nx1,nxy1)
         implicit none
         integer :: idimp, nop, nx1, nxy1
         real, dimension(idimp,nop) :: part
         real, dimension(idimp,2) :: pt2
         integer, dimension(2,nop) :: ip2
         integer, dimension(nxy1) :: npic
         end subroutine
      end interface
      interface
         subroutine RMOVE2(part,ihole,nx,ny,idimp,nop,ntmax,ipbc)
         implicit none
         integer :: nx, ny, idimp, nop, ntmax, ipbc
         real, dimension(idimp,nop) :: part
         integer, dimension(ntmax) :: ihole
         end subroutine
      end interface
      interface
         subroutine PUSH2ZF(part,dt,ek,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, ipbc
         real :: dt, ek
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine DPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh,nyv) :: q
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine DPOST2GLX(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,n&
     &pp)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, npp
         real :: qm
         real, dimension(idimp,nop) :: part
         real, dimension(2*nxvh,nyv) :: q
         complex, dimension(nxvh,npp) :: sctxp
         end subroutine
      end interface
      interface
         subroutine PUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxvh&
     &,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine PUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxvh,nyv,npp,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, npp, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,2*nxvh,nyv) :: fxy
         complex, dimension(nxvh,npp) :: sctxp
         double precision, dimension(2,npp) :: exyp
         end subroutine
      end interface
      interface
         subroutine GCJPOST2(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, cu
         end subroutine
      end interface
      interface
         subroutine GCJPOST2L(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, cu
         end subroutine
      end interface
      interface
         subroutine GCJPOST2C(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
         implicit none
         integer :: nop, idimp, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, cu
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface dpost
         module procedure igpost2
      end interface
!
      interface dpostgl
         module procedure idpost2gl
      end interface
!
      interface dpostglx
         module procedure idpost2glx
      end interface
!
      interface push
         module procedure igpush2
      end interface
!
      interface pushzf
         module procedure ipush2zf
      end interface
!
      interface pushgl
         module procedure ipush2gl
      end interface
!
      interface pushglx
         module procedure ipush2glx
      end interface
!
      interface sortp
         module procedure isortp2y
         module procedure idsortp2y
         module procedure iipsortp2y
      end interface
!
      interface sort2p
         module procedure isortp2
         module procedure iipsortp2
      end interface
!
      interface rmove
         module procedure irmove2
      end interface
!
      interface gcjpost
         module procedure igcjpost2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine igpost2(part,q,nop,qm,tdpost,inorder,dopt)
! deposit charge
         implicit none
         integer :: nop
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: td
         idimp = size(part,1)
         nxv = size(q,1); nyv = size(q,2); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         !print*, 'here2'
         call wtimer(td,ltime,-1)
         if (order==LINEAR) then
            if (opt==LOOKAHEAD) then
               call GSPOST2L(part,q,qm,nop,idimp,nxv,nxyv)
            else if (opt==VECTOR) then
               call GSPOST2XL(part,q,qm,nop,idimp,nxv,nxyv)
            else
               call GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
            endif
         else if (order==CUBIC) then
            call GPOST2C(part,q,qm,nop,idimp,nxv,nyv)
         else
            if (opt==LOOKAHEAD) then
               call GSPOST2(part,q,qm,nop,idimp,nxv,nxyv)
            else if (opt==VECTOR) then
               call GSPOST2X(part,q,qm,nop,idimp,nxv,nxyv)
            else
               call GPOST2(part,q,qm,nop,idimp,nxv,nyv)
            endif
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine igpost2
!
         subroutine igpush2(part,fxy,nop,qbm,dt,ek,tpush,nx,ny,ipbc,inor&
     &der,popt)
! push particles with 2d electrostatic fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
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
               call GSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc)
            else
               call GPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
            endif
         else if (order==CUBIC) then
            call GPUSH2C(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipbc&
     &)
         else
            if (opt==LOOKAHEAD) then
               call GSPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,&
     &ipbc)
            else
               call GPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ip&
     &bc)
            endif
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igpush2
!
         subroutine isortp2y(part,pt,ip,nop,npic,tsort,inorder)
! sort particles by y grid using memory-conserving bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: pt
         integer, dimension(:), pointer :: ip, npic
! local data
         integer :: idimp, ny1, order, ltime
         real :: ts
         idimp = size(part,1)
         ny1 = size(npic)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call SORTP2YL(part,pt,ip,npic,idimp,nop,ny1)
         else
            call SORTP2Y(part,pt,ip,npic,idimp,nop,ny1)
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine isortp2y
!
         subroutine idsortp2y(parta,partb,nop,npic,tsort,inorder)
! sort particles by y grid using optimized bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: parta, partb
         integer, dimension(:), pointer :: npic
! local data
         integer :: idimp, ny1, order, ltime
         real :: ts
         real, dimension(:,:), pointer :: part
         idimp = size(parta,1)
         ny1 = size(npic)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
         else
             call DSORTP2Y(parta,partb,npic,idimp,nop,ny1)
         endif
         part => parta
         parta => partb
         partb => part
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine idsortp2y
!
         subroutine iipsortp2y(part,pt2,ip2,nop,npic,tsort,inorder)
! sort particles by y grid using in place bin sort
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: tsort
         real, dimension(:,:), pointer :: part, pt2
         integer, dimension(:,:), pointer :: ip2
         integer, dimension(:), pointer :: npic
! local data
         integer :: idimp, ny1, order, ltime
         real :: ts
         idimp = size(part,1)
         ny1 = size(npic)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call IPSORTP2YL(part,pt2,ip2,npic,idimp,nop,ny1)
         else
            call IPSORTP2Y(part,pt2,ip2,npic,idimp,nop,ny1)
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine iipsortp2y
!
         subroutine isortp2(part,pt,ip,nop,npic,tsort,nx1,inorder)
! sort particles by x, y grid using memory-conserving bin sort
         implicit none
         integer :: nop, nx1
         real :: tsort
         integer, optional :: inorder
         real, dimension(:,:), pointer :: part
         real, dimension(:), pointer :: pt
         integer, dimension(:), pointer :: ip, npic
! local data
         integer :: idimp, nxy1, order, ltime
         real :: ts
         idimp = size(part,1)
         nxy1 = size(npic)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call SORTP2L(part,pt,ip,npic,idimp,nop,nx1,nxy1)
         else
            call SORTP2(part,pt,ip,npic,idimp,nop,nx1,nxy1)
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine isortp2
!
         subroutine iipsortp2(part,pt2,ip2,nop,npic,tsort,nx1,inorder)
! sort particles by x, y grid using in place bin sort
         implicit none
         integer :: nop, nx1
         real :: tsort
         integer, optional :: inorder
         real, dimension(:,:), pointer :: part, pt2
         integer, dimension(:,:), pointer :: ip2
         integer, dimension(:), pointer :: npic
! local data
         integer :: idimp, nxy1, order, ltime
         real :: ts
         idimp = size(part,1)
         nxy1 = size(npic)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(ts,ltime,-1)
         if ((order==LINEAR).or.(order==CUBIC)) then
            call IPSORTP2L(part,pt2,ip2,npic,idimp,nop,nx1,nxy1)
         else
            call IPSORTP2(part,pt2,ip2,npic,idimp,nop,nx1,nxy1)
         endif
! record time
         call wtimer(ts,ltime)
         tsort = tsort + ts
         end subroutine iipsortp2
!
         subroutine irmove2(part,nop,nx,ny,ntmax,ipbc)
! removes particles which would normally be reflected
         implicit none
         integer :: nop, nx, ny, ntmax, ipbc
         real, dimension(:,:), pointer :: part
! local data
         integer, dimension(ntmax) :: ihole
         integer :: idimp
         idimp = size(part,1)
         call RMOVE2(part,ihole,nx,ny,idimp,nop,ntmax,ipbc)
         end subroutine irmove2
!
         subroutine ipush2zf(part,nop,dt,ek,tpush,nx,ny,ipbc)
! push particles with no forces
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: dt, ek, tpush
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, ltime
         real :: tp
         idimp = size(part,1)
! initialize timer
         call wtimer(tp,ltime,-1)
         call PUSH2ZF(part,dt,ek,idimp,nop,nx,ny,ipbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ipush2zf
!
         subroutine idpost2gl(part,q,nop,qm,nx,ny,tdpost)
! deposit charge using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: td
         complex, dimension(size(q,1)/2) :: sctx
         idimp = size(part,1)
         nxvh = size(q,1)/2; nyv = size(q,2)
! initialize timer
         call wtimer(td,ltime,-1)
         call DPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxvh,nyv)
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine idpost2gl
!
         subroutine idpost2glx(part,q,nop,qm,nx,ny,tdpost,dopt)
! deposit charge using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, tdpost
         integer, optional :: dopt
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         integer, parameter :: npp = 32
         integer :: idimp, nxvh, nyv, opt, ltime
         real :: td
         complex, dimension(size(q,1)/2,npp) :: sctxp
         idimp = size(part,1)
         nxvh = size(q,1)/2; nyv = size(q,2)
         opt = STANDARD
         if (present(dopt)) opt = dopt
! initialize timer
         call wtimer(td,ltime,-1)
         if (opt==LOOKAHEAD) then
           call DPOST2GLX(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,npp)
         else
            call DPOST2GL(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv)
         endif
! record time
         call wtimer(td,ltime)
         tdpost = tdpost + td
         end subroutine idpost2glx
!
         subroutine ipush2gl(part,fxy,nop,qbm,dt,ek,nx,ny,ipbc,tpush)
! push particles with 2d electrostatic fields using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ek, tpush
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
         call PUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxvh,nyv,i&
     &pbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ipush2gl
!
         subroutine ipush2glx(part,fxy,nop,qbm,dt,ek,nx,ny,ipbc,tpush,po&
     &pt)
! push particles with 2d electrostatic fields using gridless method
! optimized method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ek, tpush
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
            call PUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ek,idimp,nop,nx,ny,&
     &nxvh,nyv,npp,ipbc)
         else
            call PUSH2GL(part,fxy,sctxp,qbm,dt,ek,idimp,nop,nx,ny,nxvh,n&
     &yv,ipbc)
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ipush2glx
!
         subroutine igcjpost2(part,fxy,cu,nop,qm,qbm,dt,tdcjpost,inorder&
     &)
! deposit current density with 2d electrostatic fields
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qm, qbm, dt, tdcjpost
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
            call GCJPOST2L(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
         else if (order==CUBIC) then
            call GCJPOST2C(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
         else
            call GCJPOST2(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
         endif
! record time
         call wtimer(tdc,ltime)
         tdcjpost = tdcjpost + tdc
         end subroutine igcjpost2
!
         subroutine initmomt2(part,nop,px,py,pz,ndim)
! calculate initial momentum, for 2 or 2-1/2d code
         integer :: nop
         integer, optional :: ndim
         real :: px, py, pz
         real, dimension(:,:), pointer :: part
! local data
         integer :: j, nd
         double precision :: sum1, sum2, sum3
         nd = 3
         if (present(ndim)) nd = ndim
! calculate momentum at t=t-dt/2
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         select case(nd)
         case (2)
            do j = 1, nop
            sum1 = sum1 + part(3,j)
            sum2 = sum2 + part(4,j)
            enddo
            px = sum1
            py = sum2
            pz = 0.0
         case (3)
            do j = 1, nop
            sum1 = sum1 + part(3,j)
            sum2 = sum2 + part(4,j)
            sum3 = sum3 + part(5,j)
            enddo
            px = sum1
            py = sum2
            pz = sum3
         end select
         end subroutine initmomt2
!
         subroutine premoment2(part,itime,nop,iunit,px,py,pz,sx,sy,sz,wx&
     &,wy,wz,ndim,nprint)
! print out electron and field momentum, calculate total momentum
! for 2 or 2-1/2d code
         integer :: itime, nop, iunit
         integer, optional :: ndim, nprint
         real :: px, py, pz, sx, sy, sz, wx, wy, wz
         real, dimension(:,:), pointer :: part
! local data
         integer :: j, nd, np
         double precision :: sum1, sum2, sum3
  991    format (' T = ',i7)
  994    format (' electron momentum = ',3e14.7)
  996    format (' field momentum = ',3e14.7)
         nd = 3
         if (present(ndim)) nd = ndim
         np = 1
         if (present(nprint)) np = nprint
         if (np < 0) return
         if (np==0) then
            px = 0.0; py = 0.0; pz = 0.0
         else if (np==1) then
            write (iunit,991) itime
         endif
! calculate and print electron momentum
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         select case(nd)
         case (2)
            do j = 1, nop
            sum1 = sum1 + part(3,j)
            sum2 = sum2 + part(4,j)
            enddo
            px = 0.5*(px + sum1)
            py = 0.5*(py + sum2)
            pz = 0.0
         case (3)
            do j = 1, nop
            sum1 = sum1 + part(3,j)
            sum2 = sum2 + part(4,j)
            sum3 = sum3 + part(5,j)
            enddo
            px = 0.5*(px + sum1)
            py = 0.5*(py + sum2)
            pz = 0.5*(pz + sum3)
         end select
         if (np==1) then
! print electron momentum
            write (iunit,994) px, py, pz
! print field momentum
            write (iunit,996) sx, sy, sz
! calculate total momentum
            wx = px + sx
            wy = py + sy
            wz = pz + sz
         endif
         px = sum1
         py = sum2
         pz = sum3
         end subroutine premoment2
!
         subroutine primoment2(parti,nopi,iunit,rmass,px,py,pz,wx,wy,wz,&
     &ndim,nprint)
! print out ion momentum, adds total momentum, for 2 or 2-1/2d code
         integer :: nopi, iunit
         integer, optional :: ndim, nprint
         real :: rmass, px, py, pz, wx, wy, wz
         real, dimension(:,:), pointer :: parti
! local data
         integer :: j, nd, np
         real :: at1
         double precision :: sum1, sum2, sum3
  995    format (' ion momentum = ',3e14.7)
         nd = 3
         if (present(ndim)) nd = ndim
         np = 1
         if (present(nprint)) np = nprint
         if (np < 0) return
         at1 = 0.5*rmass
         if (np==0) then
            px = 0.0; py = 0.0; pz = 0.0
         endif
! calculate and print ion momentum
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         select case(nd)
         case (2)
            do j = 1, nopi
            sum1 = sum1 + parti(3,j)
            sum2 = sum2 + parti(4,j)
            enddo
            px = at1*(px + sum1)
            py = at1*(py + sum2)
            pz = 0.0
         case (3)
            do j = 1, nopi
            sum1 = sum1 + parti(3,j)
            sum2 = sum2 + parti(4,j)
            sum3 = sum3 + parti(5,j)
            enddo
            px = at1*(px + sum1)
            py = at1*(py + sum2)
            pz = at1*(pz + sum3)
         end select
         if (np==1) then
! print ion momentum
            write (iunit,995) px, py, pz
! add to total momentum
            wx = wx + px
            wy = wy + py
            wz = wz + pz
         endif
         px = sum1
         py = sum2
         pz = sum3
         end subroutine primoment2
!
      end module push2d
