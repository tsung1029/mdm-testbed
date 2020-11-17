!-----------------------------------------------------------------------
!
      module hdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library hdpush2lib.f
! written by viktor k. decyk, ucla
! copyright 2007, regents of the university of california
! update: october 27, 2007
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, cr8pcm, dhjpost, hpush, hpushxh, hpushx
      public :: precmoment2, pricmoment2
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine CR8PCM23(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nx&
     &v,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,n&
     &xv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ip&
     &bc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,&
     &ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine CR8PCM22(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,n&
     &xv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,&
     &nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny&
     &,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ip&
     &bc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,&
     &ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine CR8PCM23L(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,n&
     &xv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,&
     &nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine CR8PCM22L(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,&
     &nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop&
     &,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,nxyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine CR8PCM23GL(part,axy,sctx,qbm,sx,sy,sz,idimp,nop,nx,n&
     &y,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qbm, sx, sy, sz
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: axy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GHJPOST2GL(part,fxy,axy,daxy,cu,sctx,qm,qbm,dt,idimp&
     &,nop,nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: fxy, axy, cu
         complex, dimension(5,nxvh,nyv) :: daxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GHPUSH23GL(part,fxy,axy,daxy,sctx,qbm,dt,ek,idimp,no&
     &p,nx,ny,nxvh,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: fxy, axy
         complex, dimension(5,nxvh,nyv) :: daxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GHPUSHXH23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nx&
     &vh,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: axy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GHPUSHX23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxv&
     &h,nyv,ipbc)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv, ipbc
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: axy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface cr8pcm
         module procedure icr8pcm2
         module procedure icr8pcm2gl
      end interface
!
      interface dhjpost
         module procedure ighjpost2
         module procedure ighjpost2gl
      end interface
!
      interface hpush
         module procedure ighpush23
         module procedure ighpush23gl
      end interface
!
      interface hpushxh
         module procedure ighpushxh23
         module procedure ighpushxh23gl
      end interface
!
      interface hpushx
         module procedure ighpushx23
         module procedure ighpushx23gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine icr8pcm2(part,axy,nop,qbm,inorder)
! calculate canonical momentum for particles with 2-1/2d darwin fields
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, order
         real :: sx, sy, sz
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               call CR8PCM22L(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv)
            else
               call CR8PCM22(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv)
            endif
         case (3)
            if (order==LINEAR) then
               call CR8PCM23L(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv)
            else
               call CR8PCM23(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv)
            endif
         end select
         end subroutine icr8pcm2
!
         subroutine ighjpost2(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,tdhjpos&
     &t,inorder,djopt)
! deposit current density with 2-1/2d darwin fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, tdhjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, axy, daxy, cu
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tj
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv)
               else
                  call GHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,n&
     &op,nxv,nxyv)
               else
                  call GHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,no&
     &p,nxv,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,n&
     &op,nxv,nxyv)
               else
                  call GHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,no&
     &p,nxv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,no&
     &p,nxv,nxyv)
               else
                  call GHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop&
     &,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tj,ltime)
         tdhjpost = tdhjpost + tj
         end subroutine ighjpost2
!
         subroutine ighpush23(part,fxy,axy,daxy,nop,qbm,dt,ek,tpush,nx,n&
     &y,ipbc,inorder,popt)
! update particle canonical momentum with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, axy, daxy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc)
               else
                  call GHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
               else
                  call GHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc)
               else
                  call GHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,n&
     &x,ny,nxv,nxyv,ipbc)
               else
                  call GHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxv,nyv,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ighpush23
!
         subroutine ighpushxh23(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc,ino&
     &rder,popt)
! update particle position to t+dt/2 with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc)
               else
                  call GHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc)
               else
                  call GHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc)
               else
                  call GHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc)
               else
                  call GHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ighpushxh23
!
         subroutine ighpushx23(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc,inor&
     &der,popt)
! update particle position to t+dt with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, opt, ltime
         real :: tp
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(popt)) opt = popt
! initialize timer
         call wtimer(tp,ltime,-1)
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc)
               else
                  call GHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nx&
     &yv,ipbc)
               else
                  call GHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv&
     &,ipbc)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc)
               else
                  call GHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nx&
     &yv,ipbc)
               else
                  call GHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv&
     &,ipbc)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ighpushx23
!
         subroutine icr8pcm2gl(part,axy,nop,nx,ny,qbm)
! calculate canonical momentum for particles with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qbm
         real, dimension(:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxvh, nyv
         real :: sx, sy, sz
         complex, dimension(size(axy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(axy,2); nyv = size(axy,3)
         call CR8PCM23GL(part,axy,sctx,qbm,sx,sy,sz,idimp,nop,nx,ny,nxvh&
     &,nyv)
         end subroutine icr8pcm2gl
!
         subroutine ighjpost2gl(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,tdhjp&
     &ost,nx,ny)
! deposit current density with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, qbm, dt, tdhjpost
         real, dimension(:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: fxy, axy, daxy, cu
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tj
         complex, dimension(size(axy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(axy,2); nyv = size(axy,3)
! initialize timer
         call wtimer(tj,ltime,-1)
         call GHJPOST2GL(part,fxy,axy,daxy,cu,sctx,qm,qbm,dt,idimp,nop,n&
     &x,ny,nxvh,nyv)
! record time
         call wtimer(tj,ltime)
         tdhjpost = tdhjpost + tj
         end subroutine ighjpost2gl
!
         subroutine ighpush23gl(part,fxy,axy,daxy,nop,qbm,dt,ek,tpush,nx&
     &,ny,ipbc)
! update particle canonical momentum with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, ek, tpush
         real, dimension(:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: fxy, axy, daxy
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tp
         complex, dimension(size(axy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(axy,2); nyv = size(axy,3)
! initialize timer
         call wtimer(tp,ltime,-1)
         call GHPUSH23GL(part,fxy,axy,daxy,sctx,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxvh,nyv,ipbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ighpush23gl
!
         subroutine ighpushxh23gl(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc)
! update particle position to t+dt/2 with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, tpush
         real, dimension(:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tp
         complex, dimension(size(axy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(axy,2); nyv = size(axy,3)
! initialize timer
         call wtimer(tp,ltime,-1)
         call GHPUSHXH23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxvh,nyv&
     &,ipbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ighpushxh23gl
!
         subroutine ighpushx23gl(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc)
! update particle position to t+dt with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc
         real :: qbm, dt, tpush
         real, dimension(:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxvh, nyv, ltime
         real :: tp
         complex, dimension(size(axy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(axy,2); nyv = size(axy,3)
! initialize timer
         call wtimer(tp,ltime,-1)
         call GHPUSHX23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxvh,nyv,&
     &ipbc)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine ighpushx23gl
!
         subroutine precmoment2(part,itime,nop,iunit,wx,wy,wz,ndim)
! print out electron and field momentum, calculate total momentum
! for 2 or 2-1/2d code
         integer :: itime, nop, iunit
         integer, optional :: ndim
         real :: wx, wy, wz
         real, dimension(:,:), pointer :: part
! local data
         integer :: j, nd = 3
         real :: px, py, pz, sx, sy, sz
         double precision :: sum1, sum2, sum3, sum4, sum5, sum6
  991    format (' T = ',i7)
  994    format (' electron momentum = ',3e14.7)
  996    format (' electron field momentum = ',3e14.7)
         if (present(ndim)) nd = ndim
         write (iunit,991) itime
! calculate and print electron momentum
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         sum4 = 0.0d0
         sum5 = 0.0d0
         sum6 = 0.0d0
         select case(nd)
         case (2)
            do j = 1, nop
            sum1 = sum1 + part(3,j)
            sum2 = sum2 + part(4,j)
            sum3 = sum3 + part(5,j)
            sum4 = sum4 + part(6,j)
            enddo
            px = sum1
            py = sum2
            wx = sum3
            wy = sum4
            pz = 0.0
            wz = 0.0
         case (3)
            do j = 1, nop
            sum1 = sum1 + part(3,j)
            sum2 = sum2 + part(4,j)
            sum3 = sum3 + part(5,j)
            sum4 = sum4 + part(6,j)
            sum5 = sum5 + part(7,j)
            sum6 = sum6 + part(8,j)
            enddo
            px = sum1
            py = sum2
            pz = sum3
            wx = sum4
            wy = sum5
            wz = sum6
         end select
         write (iunit,994) px, py, pz
! print field momentum
         sx = wx - px
         sy = wy - py
         sz = wz - pz
         write (iunit,996) sx, sy, sz
         end subroutine precmoment2
!
         subroutine pricmoment2(parti,nopi,iunit,rmass,wx,wy,wz,ndim)
! print out ion momentum, adds total momentum, for 2 or 2-1/2d code
         integer :: nopi, iunit
         integer, optional :: ndim
         real :: rmass, wx, wy, wz
         real, dimension(:,:), pointer :: parti
! local data
         integer :: j, nd = 3
         real :: px, py, pz, sx, sy, sz
         double precision :: sum1, sum2, sum3, sum4, sum5, sum6
  995    format (' ion momentum = ',3e14.7)
  996    format (' ion field momentum = ',3e14.7)
         if (present(ndim)) nd = ndim
! calculate and print ion momentum
         sum1 = 0.0d0
         sum2 = 0.0d0
         sum3 = 0.0d0
         select case(nd)
         case (2)
            do j = 1, nopi
            sum1 = sum1 + parti(3,j)
            sum2 = sum2 + parti(4,j)
            sum3 = sum3 + parti(5,j)
            sum4 = sum4 + parti(6,j)
            enddo
            px = rmass*sum1
            py = rmass*sum2
            sx = rmass*sum3
            sy = rmass*sum4
            pz = 0.0
            sz = 0.0
         case (3)
            do j = 1, nopi
            sum1 = sum1 + parti(3,j)
            sum2 = sum2 + parti(4,j)
            sum3 = sum3 + parti(5,j)
            sum4 = sum4 + parti(6,j)
            sum5 = sum5 + parti(7,j)
            sum6 = sum6 + parti(8,j)
            enddo
            px = rmass*sum1
            py = rmass*sum2
            pz = rmass*sum3
            sx = rmass*sum4
            sy = rmass*sum5
            sz = rmass*sum6
         end select
         write (iunit,995) px, py, pz
! add total momentum
         wx = wx + sx
         wy = wy + sy
         wz = wz + sz
! print ion field momentum
         sx = sx - px
         sy = sy - py
         sz = sz - pz
         write (iunit,996) sx, sy, sz
         end subroutine pricmoment2
!
      end module hdpush2d
