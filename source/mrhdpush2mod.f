!-----------------------------------------------------------------------
!
      module mrhdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mrhdpush2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 19, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use rhdpush2d, only: wtimer
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: rcr8pcm, rdhjpost, rhpush, rhpushxh, rhpushx
!
! buffer data for current deposit
      real, dimension(:,:,:,:), allocatable :: cup
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MRCR8PCM23(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,ny&
     &v,sp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, ci, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         real, dimension(3,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRCR8PCM23L(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,n&
     &yv,sp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, ci, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         real, dimension(3,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRCR8PCM22(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv,s&
     &p,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, ci, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         real, dimension(2,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MRCR8PCM22L(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv,&
     &sp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, ci, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         real, dimension(2,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,n&
     &op,nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,&
     &nop,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         real, dimension(3,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,&
     &nop,nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp&
     &,nop,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         real, dimension(3,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,&
     &nop,nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp&
     &,nop,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         real, dimension(2,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp&
     &,nop,nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         real, dimension(2,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop&
     &,nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,no&
     &p,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &yv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n&
     &yv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,&
     &nyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv&
     &,nxyv,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv,  nmt, ierr
         real :: qbm, dt, ci
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface rcr8pcm
         module procedure imrcr8pcm2
      end interface
!
      interface rdhjpost
         module procedure imgrhjpost2
      end interface
!
      interface rhpush
         module procedure imgrhpush23
      end interface
!
      interface rhpushxh
         module procedure imgrhpushxh23
      end interface
!
      interface rhpushx
         module procedure imgrhpushx23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imrcr8pcm2(part,axy,nop,qbm,ci,inorder)
! multi-tasking canonical momentum calculation with 2-1/2d darwin fields
! for relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, order, ierr
         real :: ci, sx, sy, sz
         real, dimension(size(axy,1),ntasks) :: sp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3)
         nxyv = nxv*nyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               call MRCR8PCM22L(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv,&
     &sp,idtask,nmt,ierr)
            else
               call MRCR8PCM22(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv,s&
     &p,idtask,nmt,ierr)
            endif
         case (3)
            if (order==LINEAR) then
               call MRCR8PCM23L(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,n&
     &yv,sp,idtask,nmt,ierr)
            else
               call MRCR8PCM23(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,ny&
     &v,sp,idtask,nmt,ierr)
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
         end subroutine imrcr8pcm2
!
         subroutine imgrhjpost2(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,ci,td&
     &hjpost,inorder,djopt)
! multi-tasking deposit current deposit with 2-1/2d darwin fields
! for relativistic particles
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdhjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, axy, daxy, cu
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, order, opt, ltime, ierr
         integer :: nnxyv
         real :: tj
!        real, dimension(size(cu,1),size(cu,2),size(cu,3),ntasks) :: cup
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         nnxyv = size(cu,1)*nxyv
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
         opt = STANDARD
         if (present(djopt)) opt = djopt
! check if size of buffer has changed
         if (szbuf < nnxyv) then
            if (szbuf /= 0) deallocate(cup)
! allocate buffer
            allocate(cup(size(cu,1),nxv,nyv,ntasks))
            szbuf = nnxyv
         endif
! initialize timer
         call wtimer(tj,ltime,-1)
         select case(size(axy,1))
         case (2)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,i&
     &dimp,nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,id&
     &imp,nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idi&
     &mp,nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idim&
     &p,nop,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tj,ltime)
         tdhjpost = tdhjpost + tj
         end subroutine imgrhjpost2
!
         subroutine imgrhpush23(part,fxy,axy,daxy,nop,qbm,dt,ci,ek,tpush&
     &,nx,ny,inorder,popt)
! multi-tasking canonical momentum update with 2-1/2d darwin fields
! for relativistic particles
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, axy, daxy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, order, opt, ltime, ierr
         real :: tp
         real, dimension(ntasks) :: ekp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         nmt = ntasks
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
                  call MGSRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp&
     &,nop,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
               else
                  call MGRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,&
     &nop,nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
               else
                  call MGRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp&
     &,nop,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
               else
                  call MGRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,&
     &nop,nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,&
     &nop,nx,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
               else
                  call MGRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,n&
     &op,nx,ny,nxv,nyv,ekp,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgrhpush23
!
         subroutine imgrhpushxh23(part,axy,nop,qbm,dt,ci,tpush,nx,ny,ino&
     &rder,popt)
! multi-tasking canonical momentum update of particle position to t+dt/2
! with 2-1/2d darwin fields for relativistic particles
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, order, opt, ltime, ierr
         real :: tp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         nmt = ntasks
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
                  call MGSRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny&
     &,nxv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nyv,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nyv,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny&
     &,nxv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nyv,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nyv,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgrhpushxh23
!
         subroutine imgrhpushx23(part,axy,nop,qbm,dt,ci,tpush,nx,ny,inor&
     &der,popt)
! multi-tasking canonical momentum update of particle position to t+dt
! with 2-1/2d darwin fields for relativistic particles
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, nxyv, order, nmt, opt, ltime, ierr
         real :: tp
         integer, dimension(ntasks) :: idtask
         idimp = size(part,1)
         nxv = size(axy,2); nyv = size(axy,3); nxyv = nxv*nyv
         nmt = ntasks
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
                  call MGSRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nyv,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nyv,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,&
     &nxv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nyv,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,n&
     &xv,nxyv,idtask,nmt,ierr)
               else
                  call MGRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nx&
     &v,nyv,idtask,nmt,ierr)
               endif
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tp,ltime)
         tpush = tpush + tp
         end subroutine imgrhpushx23
!
      end module mrhdpush2d
