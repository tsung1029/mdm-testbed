!-----------------------------------------------------------------------
!
      module rhdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library rhdpush2lib.f
! written by viktor k. decyk, ucla
! copyright 2007, regents of the university of california
! update: september 12, 2007
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: wtimer, rcr8pcm, rdhjpost, rhpush, rhpushxh, rhpushx
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine RCR8PCM23(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,nyv&
     &)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, ci, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,no&
     &p,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,n&
     &op,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &yv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,ny&
     &v)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &xyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine RCR8PCM22(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, ci, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,n&
     &op,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,&
     &nop,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,n&
     &x,ny,nxv,nyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &yv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,ny&
     &v)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &xyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine RCR8PCM23L(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,ny&
     &v)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, ci, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,n&
     &op,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,&
     &nop,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &yv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine RCR8PCM22L(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qbm, ci, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,&
     &nop,nxv,nyv)
         implicit none
         integer :: idimp, nop, nxv, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp&
     &,nop,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nxv, nxyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxv,nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &yv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         end subroutine
      end interface
      interface
         subroutine GSRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nxyv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         end subroutine
      end interface
      interface
         subroutine RCR8PCM23GL(part,axy,sctx,qbm,ci,sx,sy,sz,idimp,nop,&
     &nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qbm, ci, sx, sy, sz
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: axy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine RGHJPOST2GL(part,fxy,axy,daxy,cu,sctx,qm,qbm,dt,ci,i&
     &dimp,nop,nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: fxy, axy, cu
         complex, dimension(5,nxvh,nyv) :: daxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GRHPUSH23GL(part,fxy,axy,daxy,sctx,qbm,dt,ci,ek,idim&
     &p,nop,nx,ny,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: fxy, axy
         complex, dimension(5,nxvh,nyv) :: daxy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GRHPUSHXH23GL(part,axy,sctx,qbm,dt,ci,idimp,nop,nx,n&
     &y,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: axy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
      interface
         subroutine GRHPUSHX23GL(part,axy,sctx,qbm,dt,ci,idimp,nop,nx,ny&
     &,nxvh,nyv)
         implicit none
         integer :: nop, idimp, nx, ny, nxvh, nyv
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         complex, dimension(3,nxvh,nyv) :: axy
         complex, dimension(nxvh) :: sctx
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rcr8pcm
         module procedure ircr8pcm2
         module procedure ircr8pcm2gl
      end interface
!
      interface rdhjpost
         module procedure igrhjpost2
         module procedure igrhjpost2gl
      end interface
!
      interface rhpush
         module procedure igrhpush23
         module procedure igrhpush23gl
      end interface
!
      interface rhpushxh
         module procedure igrhpushxh23
         module procedure igrhpushxh23gl
      end interface
!
      interface rhpushx
         module procedure igrhpushx23
         module procedure igrhpushx23gl
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ircr8pcm2(part,axy,nop,qbm,ci,inorder)
! calculate canonical momentum for particles with 2-1/2d darwin fields
! for relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qbm, ci
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
               call RCR8PCM22L(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv)
            else
               call RCR8PCM22(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv)
            endif
         case (3)
            if (order==LINEAR) then
               call RCR8PCM23L(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,ny&
     &v)
            else
               call RCR8PCM23(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,nyv&
     &)
            endif
         end select
         end subroutine ircr8pcm2
!
         subroutine igrhjpost2(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,ci,tdh&
     &jpost,inorder,djopt)
! deposit current density with 2-1/2d darwin fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdhjpost
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
                  call GSRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv)
               else
                  call GRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nxyv)
               else
                  call GRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nxyv)
               else
                  call GRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nxyv)
               else
                  call GRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp&
     &,nop,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tj,ltime)
         tdhjpost = tdhjpost + tj
         end subroutine igrhjpost2
!
         subroutine igrhpush23(part,fxy,axy,daxy,nop,qbm,dt,ci,ek,tpush,&
     &nx,ny,inorder,popt)
! update particle canonical momentum with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
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
                  call GSRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv)
               else
                  call GRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv)
               else
                  call GRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv)
               else
                  call GRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,n&
     &op,nx,ny,nxv,nxyv)
               else
                  call GRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,no&
     &p,nx,ny,nxv,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrhpush23
!
         subroutine igrhpushxh23(part,axy,nop,qbm,dt,ci,tpush,nx,ny,inor&
     &der,popt)
! update particle position to t+dt/2 with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, tpush
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
                  call GSRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nxyv)
               else
                  call GRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nxyv)
               else
                  call GRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nxyv)
               else
                  call GRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nxyv)
               else
                  call GRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrhpushxh23
!
         subroutine igrhpushx23(part,axy,nop,qbm,dt,ci,tpush,nx,ny,inord&
     &er,popt)
! update particle position to t+dt with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, tpush
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
                  call GSRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nxyv)
               else
                  call GRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nxyv)
               else
                  call GRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nyv)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call GSRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nxyv)
               else
                  call GRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nyv)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call GSRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nxyv)
               else
                  call GRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nyv)
               endif
            endif
         end select
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrhpushx23
!
         subroutine ircr8pcm2gl(part,axy,nop,nx,ny,qbm,ci)
! calculate canonical momentum for particles with 2-1/2d darwin fields
! using gridless method for relativistic particles
         implicit none
         integer :: nop, nx, ny
         real :: qbm, ci
         real, dimension(:,:), pointer :: part
         complex, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxvh, nyv
         real :: sx, sy, sz
         complex, dimension(size(axy,2)) :: sctx
         idimp = size(part,1)
         nxvh = size(axy,2); nyv = size(axy,3)
         call RCR8PCM23GL(part,axy,sctx,qbm,ci,sx,sy,sz,idimp,nop,nx,ny,&
     &nxvh,nyv)
         end subroutine ircr8pcm2gl
!
         subroutine igrhjpost2gl(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,ci,t&
     &dhjpost,nx,ny)
! deposit current density with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qm, qbm, dt, ci, tdhjpost
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
         call GRHJPOST2GL(part,fxy,axy,daxy,cu,sctx,qm,qbm,dt,ci,idimp,n&
     &op,nx,ny,nxvh,nyv)
! record time
         call wtimer(tj,ltime)
         tdhjpost = tdhjpost + tj
         end subroutine igrhjpost2gl
!
         subroutine igrhpush23gl(part,fxy,axy,daxy,nop,qbm,dt,ci,ek,tpus&
     &h,nx,ny)
! update particle canonical momentum with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qbm, dt, ci, ek, tpush
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
         call GRHPUSH23GL(part,fxy,axy,daxy,sctx,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxvh,nyv)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrhpush23gl
!
         subroutine igrhpushxh23gl(part,axy,nop,qbm,dt,ci,tpush,nx,ny)
! update particle position to t+dt/2 with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qbm, dt, ci, tpush
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
         call GRHPUSHXH23GL(part,axy,sctx,qbm,dt,ci,idimp,nop,nx,ny,nxvh&
     &,nyv)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrhpushxh23gl
!
         subroutine igrhpushx23gl(part,axy,nop,qbm,dt,ci,tpush,nx,ny)
! update particle position to t+dt with 2-1/2d darwin fields
! using gridless method
         implicit none
         integer :: nop, nx, ny
         real :: qbm, dt, ci, tpush
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
         call GRHPUSHX23GL(part,axy,sctx,qbm,dt,ci,idimp,nop,nx,ny,nxvh,&
     &nyv)
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine igrhpushx23gl
!
      end module rhdpush2d
