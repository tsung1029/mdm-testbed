!-----------------------------------------------------------------------
!
      module cfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library cfield2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: october 17, 2007
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: zdbl, poisc2_init, poisc3_init, poisc
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine ZDBL2C(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(2,nx2v,ny2) :: cu2
         end subroutine
      end interface
      interface
         subroutine ZDBL2D(q,q2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,ny2) :: q2
         end subroutine
      end interface
      interface
         subroutine ZDBL2B(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx2v,ny2) :: cu2
         end subroutine
      end interface
      interface
         subroutine FORMC2(ffg,f,fpotc,mixup2,sct2,affp,ar,indx1,indy1,n&
     &x1d,ny1d,nxv,ny2v,nxhy2,nxyh2)
         implicit none
         integer :: indx1, indy1, nx1d, ny1d, nxv, ny2v, nxhy2, nxyh2
         real :: ar, affp
         real, dimension(4,nx1d,ny1d) :: ffg
         real, dimension(2*nxv,ny2v) :: f
         integer, dimension(nxhy2) :: mixup2
         complex, dimension(nxyh2) :: sct2
         real, external :: fpotc
         end subroutine
      end interface
      interface
         subroutine POISC2(q,fx,fy,isign,ffg,we,nx,ny,nx2v,ny2d,nx1d,ny1&
     &d)
         implicit none
         integer :: isign, nx, ny, nx2v, ny2d, nx1d, ny1d
         real :: we
         real, dimension(nx2v,ny2d) :: q, fx, fy
         real, dimension(4,nx1d,ny1d) :: ffg
         end subroutine
      end interface
      interface
         subroutine POISC22(q,fxy,ffg,we,nx,ny,nxv,ny2d,nx1d,ny1d)
         implicit none
         integer :: nx, ny, nxv, ny2d, nx1d, ny1d
         real :: we
         real, dimension(2*nxv,ny2d) :: q
         real, dimension(2,2*nxv,ny2d) :: fxy
         real, dimension(4,nx1d,ny1d) :: ffg
         end subroutine
      end interface
      interface
         function POTC3(r,affp,ari,ifun)
         implicit none
         integer :: ifun
         real :: POTC3, r, affp, ari
         end function
      end interface
      interface
         function POTC2(r,affp,ari,ifun)
         implicit none
         integer :: ifun
         real :: POTC2, r, affp, ari
         end function
      end interface 
!
! define generic interfaces to Fortran90 library
!
      interface zdbl
         module procedure izdbl2c
         module procedure izdbl2d
      end interface
!
      interface poisc2_init
         module procedure ipoisc2init
      end interface
!
      interface poisc3_init
         module procedure ipoisc3init
      end interface
!
      interface poisc
         module procedure ipoisc2
         module procedure ispoisc2
         module procedure ipoisc22
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine izdbl2c(cu,cu2,nx,ny,inorder)
! double array in each dimension for 2d vector data, zeroing copies
! for open boundary conditions
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
         real, dimension(:,:,:), pointer :: cu2
! local data
         integer :: nxv, nyv, nx2v, ny2, order
         nxv = size(cu,2);  nyv = size(cu,3)
         nx2v = size(cu2,2);  ny2 = size(cu2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call ZDBL2B(cu(1,1,1),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call ZDBL2B(cu(1,2,2),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call ZDBL2C(cu(1,1,1),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call ZDBL2C(cu(1,2,2),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         endif
         end subroutine izdbl2c
!
         subroutine izdbl2d(q,q2,nx,ny,inorder)
! double array in each dimension for 2d scalar data, zeroing copies
! for open boundary conditions
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
         real, dimension(:,:), pointer :: q2
! local data
         integer :: nxv, nyv, nx2v, ny2, order
         nxv = size(q,1);  nyv = size(q,2)
         nx2v = size(q2,1);  ny2 = size(q2,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call ZDBL2D(q(1,1),q2,nx,ny,nxv,nyv,nx2v,ny2)
         else
            call ZDBL2D(q(2,2),q2,nx,ny,nxv,nyv,nx2v,ny2)
         endif
         end subroutine izdbl2d
!
         subroutine ipoisc2(q,fx,ffg,we,nx,ny)
! poisson solver for 2d potential, open boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q, fx
         real, dimension(:,:,:), pointer :: ffg
! local data
         integer :: isign = 1, nx2v, ny2d, nx1d, ny1d
         real, dimension(1,1) :: fy
         nx2v = size(q,1); ny2d = size(q,2)
         nx1d = size(ffg,2); ny1d = size(ffg,3)
         call POISC2(q,fx,fy,isign,ffg,we,nx,ny,nx2v,ny2d,nx1d,ny1d)
         end subroutine ipoisc2
!
         subroutine ispoisc2(q,fy,ffg,nx,ny)
! smoother for 2d scalar field, open boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, fy
         real, dimension(:,:,:), pointer :: ffg
! local data
         integer :: isign = 2, nx2v, ny2d, nx1d, ny1d
         real :: we
         real, dimension(1,1) :: fx
         nx2v = size(q,1); ny2d = size(q,2)
         nx1d = size(ffg,2); ny1d = size(ffg,3)
         call POISC2(q,fx,fy,isign,ffg,we,nx,ny,nx2v,ny2d,nx1d,ny1d)
         end subroutine ispoisc2
!
         subroutine ipoisc2init(ffg,f,mixup2,sct2,ar,affp,indx,indy)
! initialize 2d poisson solver, open boundary conditions
         implicit none
         integer :: indx, indy
         real :: ar, affp
         real, dimension(:,:,:), pointer :: ffg
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup2
         complex, dimension(:), pointer :: sct2
! local data
         real, external :: POTC2
         integer :: indx1, indy1, nx1d, ny1d, nxv, ny2v, nxhy2, nxyh2
         indx1 = indx + 1; indy1 = indy + 1
         nx1d = size(ffg,2); ny1d = size(ffg,3)
         nxv = size(f,1)/2; ny2v = size(f,2)
         nxhy2 = size(mixup2); nxyh2 = size(sct2)
         call FORMC2(ffg,f,POTC2,mixup2,sct2,affp,ar,indx1,indy1,nx1d,ny&
     &1d,nxv,ny2v,nxhy2,nxyh2)
         end subroutine ipoisc2init
!
         subroutine ipoisc3init(ffg,f,mixup2,sct2,ar,affp,indx,indy)
! initialize 2d poisson solver with 3d fields, open boundary conditions
         implicit none
         integer :: indx, indy
         real :: ar, affp
         real, dimension(:,:,:), pointer :: ffg
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup2
         complex, dimension(:), pointer :: sct2
! local data
         real, external :: POTC3
         integer :: indx1, indy1, nx1d, ny1d, nxv, ny2v, nxhy2, nxyh2
         indx1 = indx + 1; indy1 = indy + 1
         nx1d = size(ffg,2); ny1d = size(ffg,3)
         nxv = size(f,1)/2; ny2v = size(f,2)
         nxhy2 = size(mixup2); nxyh2 = size(sct2)
         call FORMC2(ffg,f,POTC3,mixup2,sct2,affp,ar,indx1,indy1,nx1d,ny&
     &1d,nxv,ny2v,nxhy2,nxyh2)
         end subroutine ipoisc3init
!
         subroutine ipoisc22(q,fxy,ffg,we,nx,ny)
! poisson solver for 2d electric field, open boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         real, dimension(:,:,:), pointer :: ffg
! local data
         integer :: nxv, ny2d, nx1d, ny1d
         nxv = size(q,1)/2; ny2d = size(q,2)
         nx1d = size(ffg,2); ny1d = size(ffg,3)
         if (size(fxy,1)==2) then
            call POISC22(q,fxy,ffg,we,nx,ny,nxv,ny2d,nx1d,ny1d)
!        else if (size(fxy,1)==3) then
!           call POISC23(q,fxy,ffg,we,nx,ny,nxv,ny2d,nx1d,ny1d)
         endif
         end subroutine ipoisc22
!
      end module cfield2d
