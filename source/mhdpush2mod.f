!-----------------------------------------------------------------------
!
      module mhdpush2d
!
! Fortran90 interface to 2d PIC Fortran77 library mhdpush2lib.f
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: november 19, 2009
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      use hdpush2d, only: wtimer, precmoment2, pricmoment2
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: cr8pcm, dhjpost, hpush, hpushxh, hpushx
      public :: precmoment2, pricmoment2
!
! buffer data for current deposit
      real, dimension(:,:,:,:), allocatable :: cup
      integer :: szbuf = 0
      save
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MCR8PCM23(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv,sp&
     &,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         real, dimension(3,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MCR8PCM23L(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv,s&
     &p,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, sx, sy, sz
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         real, dimension(3,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MCR8PCM22(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv,sp,id&
     &task,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         real, dimension(2,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MCR8PCM22L(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv,sp,i&
     &dtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qbm, sx, sy
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         real, dimension(2,nmt) :: sp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,n&
     &xv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,&
     &nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         real, dimension(3,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,&
     &nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy, cu
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(3,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop&
     &,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy, cu
         real, dimension(5,nxyv) :: daxy
         real, dimension(3,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,&
     &nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop&
     &,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         real, dimension(2,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop&
     &,nxv,nyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy, cu
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(2,nxv,nyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,no&
     &p,nxv,nxyv,cup,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nxv, nxyv,  nmt, ierr
         real :: qm, qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy, cu
         real, dimension(3,nxyv) :: daxy
         real, dimension(2,nxyv,nmt) :: cup
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: fxy, axy
         real, dimension(5,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: fxy, axy
         real, dimension(5,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,n&
     &y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: fxy, axy
         real, dimension(3,nxv,nyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx&
     &,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt, ek
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: fxy, axy
         real, dimension(3,nxyv) :: daxy
         real, dimension(nmt) :: ekp
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv&
     &,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nx&
     &yv,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv&
     &,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nx&
     &yv,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(3,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv&
     &,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,&
     &ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxv,nyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
      interface
         subroutine MGSHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxy&
     &v,ipbc,idtask,nmt,ierr)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nxyv, ipbc, nmt, ierr
         real :: qbm, dt
         real, dimension(idimp,nop) :: part
         real, dimension(2,nxyv) :: axy
         integer, dimension(nmt) :: idtask
         end subroutine
      end interface
!
! define generic interface to Fortran90 library
!
      interface cr8pcm
         module procedure imcr8pcm2
      end interface
!
      interface dhjpost
         module procedure imghjpost2
      end interface
!
      interface hpush
         module procedure imghpush23
      end interface
!
      interface hpushxh
         module procedure imghpushxh23
      end interface
!
      interface hpushx
         module procedure imghpushx23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imcr8pcm2(part,axy,nop,qbm,inorder)
! multi-tasking canonical momentum calculation with 2-1/2d darwin fields
         implicit none
         integer :: nop
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: idimp, nxv, nyv, nxyv, nmt, order, ierr
         real :: sx, sy, sz
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
               call MCR8PCM22L(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv,sp,i&
     &dtask,nmt,ierr)
            else
               call MCR8PCM22(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv,sp,id&
     &task,nmt,ierr)
            endif
         case (3)
            if (order==LINEAR) then
               call MCR8PCM23L(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv,s&
     &p,idtask,nmt,ierr)
            else
               call MCR8PCM23(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv,sp&
     &,idtask,nmt,ierr)
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
         end subroutine imcr8pcm2
!
         subroutine imghjpost2(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,tdhjpo&
     &st,inorder,djopt)
! multi-tasking deposit current deposit with 2-1/2d darwin fields
         implicit none
         integer :: nop
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, tdhjpost
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
                  call MGSHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp&
     &,nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,&
     &nop,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,&
     &nop,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,n&
     &op,nxv,nyv,cup,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,n&
     &op,nxv,nxyv,cup,idtask,nmt,ierr)
               else
                  call MGHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,no&
     &p,nxv,nyv,cup,idtask,nmt,ierr)
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
         end subroutine imghjpost2
!
         subroutine imghpush23(part,fxy,axy,daxy,nop,qbm,dt,ek,tpush,nx,&
     &ny,ipbc,inorder,popt)
! multi-tasking canonical momentum update with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, ek, tpush
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
                  call MGSHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop&
     &,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop&
     &,nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,&
     &nx,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,&
     &nx,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
               else
                  call MGHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,n&
     &x,ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
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
         end subroutine imghpush23
!
         subroutine imghpushxh23(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc,in&
     &order,popt)
! update particle position to t+dt/2 with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, tpush
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
                  call MGSHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv&
     &,nxyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nyv,ipbc,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv&
     &,nxyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nyv,ipbc,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,idtask,nmt,ierr)
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
         end subroutine imghpushxh23
!
         subroutine imghpushx23(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc,ino&
     &rder,popt)
! update particle position to t+dt with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder, popt
         real :: qbm, dt, tpush
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
                  call MGSHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc,idtask,nmt,ierr)
               endif
            endif
         case (3)
            if (order==LINEAR) then
               if (opt==LOOKAHEAD) then
                  call MGSHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,&
     &nxyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc,idtask,nmt,ierr)
               endif
            else
               if (opt==LOOKAHEAD) then
                  call MGSHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,n&
     &xyv,ipbc,idtask,nmt,ierr)
               else
                  call MGHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc,idtask,nmt,ierr)
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
         end subroutine imghpushx23
!
      end module mhdpush2d
