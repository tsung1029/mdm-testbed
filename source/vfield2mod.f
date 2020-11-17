!-----------------------------------------------------------------------
!
      module vfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library vfield2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: october 29, 2005
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: ncguard, bndryv, poisb_init, poisb
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine NCGUARD2(fxy,bv,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         real, dimension(ny,6) :: bv
         end subroutine
      end interface
      interface
         subroutine NDGUARD2(q,bv,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         real, dimension(ny,6) :: bv
         end subroutine
      end interface
      interface
         subroutine NCGUARD2L(fxy,bv,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         real, dimension(ny,6) :: bv
         end subroutine
      end interface
      interface
         subroutine NDGUARD2L(q,bv,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         real, dimension(ny,6) :: bv
         end subroutine
      end interface
      interface
         subroutine BNDRYV2(q,ffc,bv,nx,ny,nxv,nxd,nyd,nyhd)
         implicit none
         integer :: nx, ny, nxv, nxd, nyd, nyhd
         real :: q
         complex, dimension(nxd/2,nyhd) :: ffc
         real, dimension(nyd,4) :: bv
         end subroutine
      end interface
      interface
         subroutine POISB2(fx,fy,isign,ffb,bv,bcd,mixup,sct,t,we,affp,in&
     &dx,ny,nxv,nxd,nxhd,nyd,nyhd)
         implicit none
         integer :: isign, indx, ny, nxv, nxd, nxhd, nyd, nyhd
         real ::  we, affp
         real :: fx, fy
         complex, dimension(nxd/2,nyhd) :: ffb
         real, dimension(nyd,6) :: bv
         real, dimension(nyhd) :: bcd
         integer, dimension(nxhd):: mixup
         complex, dimension(nxhd) :: sct
         real, dimension(nxd) :: t
         end subroutine
      end interface
      interface
         subroutine POISB22(fxy,isign,ffb,bv,bcd,mixup,sct,t,we,affp,ind&
     &x,ny,nxvh,nxhd,nyhd)
         implicit none
         integer :: isign, indx, ny, nxvh, nxhd, nyhd
         real :: we, affp
         real :: fxy
         complex, dimension(nxhd,nyhd) :: ffb
         real, dimension(2*nyhd,6) :: bv
         real, dimension(nyhd) :: bcd
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: sct, t
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface ncguard
         module procedure incguard2
         module procedure indguard2
      end interface
!
      interface bndryv
         module procedure ibndryv2
      end interface
!
       interface poisb_init
         module procedure ipoisb22init
      end interface
!
      interface poisb
         module procedure ipoisb2
         module procedure ipoisb22
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine incguard2(fxy,bv,nx,ny,inorder)
! copy guard cells for 2d vector data from bounary values
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         real, dimension(:,:), pointer :: bv
! local data
         integer :: nxe, nye, order
         nxe = size(fxy,2); nye = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call NCGUARD2L(fxy,bv,nx,ny,nxe,nye)
         else
            call NCGUARD2(fxy,bv,nx,ny,nxe,nye)
         endif
         end subroutine incguard2
!
         subroutine indguard2(q,bv,nx,ny,inorder)
! copy guard cells for 2d scalar data from bounary values
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
         real, dimension(:,:), pointer :: bv
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call NDGUARD2L(q,bv,nx,ny,nxe,nye)
         else
            call NDGUARD2(q,bv,nx,ny,nxe,nye)
         endif
         end subroutine indguard2
!
         subroutine ibndryv2(q,ffc,bv,nx,ny,inorder)
! calculates boundary values of electric field of periodic solution of
! poisson's equation
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, bv
         complex, dimension(:,:), pointer :: ffc
         integer, optional :: inorder
! local data
         integer :: nxv, nxd, nyd, nyhd, order
         nxv = size(q,1)
         nxd = 2*size(ffc,1); nyhd = size(ffc,2); nyd = size(bv,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call BNDRYV2(q(1,1),ffc,bv,nx,ny,nxv,nxd,nyd,nyhd)
         else
            call BNDRYV2(q(2,2),ffc,bv,nx,ny,nxv,nxd,nyd,nyhd)
         endif
         end subroutine ibndryv2
!
         subroutine ipoisb2(fx,ffb,bv,bcd,we,affp,indx,ny,inorder)
! finds corrections to 2d poisson's equation for potential
! with vacuum boundary conditions
         implicit none
         integer :: indx, ny
         integer, optional :: inorder
         real ::  we, affp
         real, dimension(:,:), pointer :: fx, bv
         complex, dimension(:,:), pointer :: ffb
         real, dimension(:), pointer :: bcd
! local data
         integer :: isign = 1, nxv, nxd, nxhd, nyd, nyhd, order
         real, dimension(1,1) :: fy, t
         integer, dimension(1) :: mixup
         complex, dimension(1) :: sct
         nxv = size(fx,1); nxd = 2*size(ffb,1); nyhd = size(ffb,2)
         nxhd = size(mixup); nyd = size(bv,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call POISB2(fx(1,1),fy(1,1),isign,ffb,bv,bcd,mixup,sct,t,we,&
     &affp,indx,ny,nxv,nxd,nxhd,nyd,nyhd)
         else
            call POISB2(fx(2,2),fy(1,1),isign,ffb,bv,bcd,mixup,sct,t,we,&
     &affp,indx,ny,nxv,nxd,nxhd,nyd,nyhd)
         endif
         end subroutine ipoisb2
!
         subroutine ipoisb22init(ffb,bcd,affp,indx,ny)
! initialize 2d electric field solver, vacuum boundaries
         implicit none
         integer :: indx, ny
         real :: affp
         complex, dimension(:,:), pointer :: ffb
         real, dimension(:), pointer :: bcd
! local data
         integer :: isign = 0, nxvh = 1, nxhd, nyhd
         real :: we, fxy
         integer, dimension(size(ffb,1)) :: mixup
         complex, dimension(size(ffb,1)) :: sct, t
         real, dimension(1,1) :: bv
         nxhd = size(ffb,1); nyhd = size(ffb,2)
         call POISB22(fxy,isign,ffb,bv,bcd,mixup,sct,t,we,affp,indx,ny,n&
     &xvh,nxhd,nyhd)
         end subroutine ipoisb22init
!
         subroutine ipoisb22(fxy,ffb,bv,bcd,we,affp,indx,ny,inorder)
! finds corrections to 2d poisson's equation for force/charge
! with vacuum boundary conditions
         implicit none
         integer :: indx, ny
         integer, optional :: inorder
         real :: we, affp
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffb
         real, dimension(:,:) :: bv
         real, dimension(:) :: bcd
! local data
         integer :: isign = -1, nxvh, nxhd, nyhd, order
         complex, dimension(1,1) :: t
         integer, dimension(1) :: mixup
         complex, dimension(1) :: sct
         nxvh = size(fxy,2)/2; nxhd = size(ffb,1); nyhd = size(ffb,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call POISB22(fxy(1,1,1),isign,ffb,bv,bcd,mixup,sct,t,we,affp&
     &,indx,ny,nxvh,nxhd,nyhd)
         else
            call POISB22(fxy(1,2,2),isign,ffb,bv,bcd,mixup,sct,t,we,affp&
     &,indx,ny,nxvh,nxhd,nyhd)
         endif
         end subroutine ipoisb22
!
      end module vfield2d
