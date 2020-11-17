!-----------------------------------------------------------------------
!
      module npfield2d
!
! Fortran90 interface to 2d PIC Fortran77 libraries field2lib.f,
! dfield2lib.f, bfield2lib.f, cfield2lib, nfield2lib, mfield2lib,
! hfield2lib
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: jyly 11, 2011
!
      use globals, only: LINEAR, QUADRATIC
      use field2d
      use dfield2d
      use bfield2d
      use cfield2d
      use nfield2d
      use mfield2d
      use hfield2d
      implicit none
      public
      private :: isguard2p, iscguard2p, ismcguard2p, isfcguard2p
      private :: iaguard2p, iacguard2p, iamcguard2p, iafcguard2p
      private :: icguard2p, idguard2p
!
! define generic interfaces to Fortran90 library
!
      interface sguardp
         module procedure isguard2p
         module procedure iscguard2p
         module procedure ismcguard2p
      end interface
!
      interface sfguardp
         module procedure isfcguard2p
      end interface
!
      interface aguardp
         module procedure iaguard2p
         module procedure iacguard2p
      end interface
! 
      interface amguardp
         module procedure iamcguard2p
      end interface
!
      interface afguardp
         module procedure iafcguard2p
      end interface
!
      interface cguardp
         module procedure icguard2p
         module procedure idguard2p
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine iscguard2p(cu,xj0,yj0,zj0,nx,ny,ipbc,inorder)
! initialize 2d vector field
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: ngx = 1, ngy = 1, nxe, nye, order
         nxe = size(cu,2); nye = size(cu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (3)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SCGUARD2L(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
               else
                  call SCGUARD2(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSCGUARD2L(cu,xj0,yj0,zj0,nx,ny,ngx,ngy,nxe,nye)
               else
                  call LSCGUARD2(cu,xj0,yj0,zj0,nx,ny,ngx,ngy,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MSCGUARD2L(cu,xj0,yj0,zj0,nx,ny,ngx,nxe,nye)
               else
                  call MSCGUARD2(cu,xj0,yj0,zj0,nx,ny,ngx,nxe,nye)
               endif
            endif
         case (2)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SCGUARD22L(cu,xj0,yj0,nx,ny,nxe,nye)
               else
                  call SCGUARD22(cu,xj0,yj0,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSCGUARD22L(cu,xj0,yj0,nx,ny,ngx,ngy,nxe,nye)
               else
                  call LSCGUARD22(cu,xj0,yj0,nx,ny,ngx,ngy,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MSCGUARD22L(cu,xj0,yj0,nx,ny,ngx,nxe,nye)
               else
                  call MSCGUARD22(cu,xj0,yj0,nx,ny,ngx,nxe,nye)
               endif
            endif
         case (1)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SGUARD2L(cu,xj0,nx,ny,nxe,nye)
               else
                  call SGUARD2(cu,xj0,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSGUARD2L(cu,xj0,nx,ny,ngx,ngy,nxe,nye)
               else
                  call LSGUARD2(cu,xj0,nx,ny,ngx,ngy,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MSGUARD2L(cu,xj0,nx,ny,ngx,nxe,nye)
               else
                  call MSGUARD2(cu,xj0,nx,ny,ngx,nxe,nye)
               endif
            endif
         end select
         end subroutine iscguard2p
!
         subroutine isguard2p(q,qi0,nx,ny,ipbc,inorder)
! initialize non-uniform 2d scalar field
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:), pointer :: q
! local data
         integer :: ngx = 1, ngy = 1, nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call SGUARD2L(q,qi0,nx,ny,nxe,nye)
            else
               call SGUARD2(q,qi0,nx,ny,nxe,nye)
            endif
         else if (ipbc==2) then
            if (order==LINEAR) then
               call LSGUARD2L(q,qi0,nx,ny,ngx,ngy,nxe,nye)
            else
               call LSGUARD2(q,qi0,nx,ny,ngx,ngy,nxe,nye)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call MSGUARD2L(q,qi0,nx,ny,ngx,nxe,nye)
            else
               call MSGUARD2(q,qi0,nx,ny,ngx,nxe,nye)
            endif
         endif
         end subroutine isguard2p
!
         subroutine isfcguard2p(cus,cu,ecu,q2m0,nx,ny,ipbc,inorder)
! initialize 2d vector field with scaled field
! copy cu field to ecu in some cases
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: q2m0
         real, dimension(:,:,:), pointer :: cu, cus, ecu
! local data
         integer :: nxe, nye, order
         real :: zero = 0.0
         nxe = size(cus,2); nye = size(cus,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cus,1))
         case (2)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
               else
                  call SCFGUARD22(cus,cu,q2m0,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
               else
                  call SCGUARD22(cus,zero,zero,nx,ny,nxe,nye)
                  ecu = cu
               endif
!           else if (ipbc==3) then
            endif
         case (3)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
               else
                  call SCFGUARD2(cus,cu,q2m0,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
               else
                  call SCGUARD2(cus,zero,zero,zero,nx,ny,nxe,nye)
                  ecu = cu
               endif
!           else if (ipbc==3) then
            endif
         end select
         end subroutine isfcguard2p
!
         subroutine ismcguard2p(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ipbc,ino&
     &rder)
! initialize 2d tensor field
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: x2y2m0, xym0, zxm0, zym0
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: ngx = 1, ngy = 1, nxe, nye, order
         nxe = size(amu,2); nye = size(amu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(amu,1))
         case (4)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,ny&
     &e)
               else
                  call SMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye&
     &)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ngx,n&
     &gy,nxe,nye)
               else
                  call LSMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ngx,ng&
     &y,nxe,nye)
               endif
            endif
         case (2)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call SMCGUARD22L(amu,x2y2m0,xym0,nx,ny,nxe,nye)
               else
                  call SMCGUARD22(amu,x2y2m0,xym0,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==LINEAR) then
                  call LSMCGUARD22L(amu,x2y2m0,xym0,nx,ny,ngx,ngy,nxe,ny&
     &e)
               else
                  call LSMCGUARD22(amu,x2y2m0,xym0,nx,ny,ngx,ngy,nxe,nye&
     &)
               endif
            endif
         end select
         end subroutine ismcguard2p
!
         subroutine iacguard2p(cu,nx,ny,ipbc,inorder)
! add guard cells for 2d vector data
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nye, order
         nxe = size(cu,2); nye = size(cu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (3)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call ACGUARD2L(cu,nx,ny,nxe,nye)
               else
                  call ACGUARD2(cu,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LACGUARD2(cu(1,2,2),nx-2,ny-2,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MACGUARD2L(cu,nx,ny,nxe,nye)
               else
                  call MACGUARD2(cu(1,2,1),nx-2,ny,nxe,nye)
               endif
            endif
          case (2)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call ACGUARD22L(cu,nx,ny,nxe,nye)
               else
                  call ACGUARD22(cu,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LACGUARD22(cu(1,2,2),nx-2,ny-2,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MACGUARD22L(cu,nx,ny,nxe,nye)
               else
                  call MACGUARD22(cu(1,2,1),nx-2,ny,nxe,nye)
               endif
            endif
         case (1)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call AGUARD2L(cu,nx,ny,nxe,nye)
               else
                  call AGUARD2(cu,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LAGUARD2(cu(1,2,2),nx-2,ny-2,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MAGUARD2L(cu,nx,ny,nxe,nye)
               else
                  call MAGUARD2(cu(1,2,1),nx-2,ny,nxe,nye)
               endif
            endif
         end select
         end subroutine iacguard2p
!
         subroutine iaguard2p(q,nx,ny,ipbc,inorder)
! add guard cells for 2d scalar data
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call AGUARD2L(q,nx,ny,nxe,nye)
            else
               call AGUARD2(q,nx,ny,nxe,nye)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call LAGUARD2(q(2,2),nx-2,ny-2,nxe,nye)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call MAGUARD2L(q,nx,ny,nxe,nye)
            else
               call MAGUARD2(q(2,1),nx-2,ny,nxe,nye)
            endif
         endif
         end subroutine iaguard2p
!
         subroutine iamcguard2p(amu,nx,ny,ipbc,inorder)
! add guard cells for 2d tensor data
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: ndim, nxe, nye, order
         ndim = size(amu,1)
         nxe = size(amu,2); nye = size(amu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call AMCGUARD2L(amu,nx,ny,nxe,nye,ndim)
            else
               call AMCGUARD2(amu,nx,ny,nxe,nye,ndim)
            endif
         else if (ipbc==2) then
            if (order /= LINEAR) then
               call LAMCGUARD2(amu(1,2,2),nx-2,ny-2,nxe,nye,ndim)
            endif
         endif
         end subroutine iamcguard2p
!
         subroutine iafcguard2p(cus,ecu,q2m0,nx,ny,ipbc,inorder)
! increment 2d vector field with scaled field
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real :: q2m0
         real, dimension(:,:,:), pointer :: cus, ecu
! local data
         integer :: nxe, nye, order
         nxe = size(cus,2); nye = size(cus,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cus,1))
         case (2)
            if (ipbc==2) then
               if (order /= LINEAR) then
                  call LACFGUARD22(cus,ecu,q2m0,nx,ny,nxe,nye)
               endif
!           else if (ipbc==3) then
            endif
         case (3)
            if (ipbc==2) then
               if (order /= LINEAR) then
                  call LACFGUARD2(cus,ecu,q2m0,nx,ny,nxe,nye)
               endif
!           else if (ipbc==3) then
            endif
         end select
         end subroutine iafcguard2p
!
         subroutine icguard2p(fxy,nx,ny,ipbc,inorder)
! copy guard cells for 2d vector data
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: ndim, nxe, nye, order
         nxe = size(fxy,2); nye = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call DGUARD2L(fxy,nx,ny,nxe,nye)
               else
                  call DGUARD2(fxy,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LDGUARD2(fxy,nx,ny,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MDGUARD2L(fxy,nx,ny,nxe,nye)
               else
                  call MDGUARD2(fxy,nx,ny,nxe,nye)
               endif
            endif
         case (2)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call CGUARD2L(fxy,nx,ny,nxe,nye)
               else
                  call CGUARD2(fxy,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LCGUARD2(fxy,nx,ny,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MCGUARD2L(fxy,nx,ny,nxe,nye)
               else
                  call MCGUARD2(fxy,nx,ny,nxe,nye)
               endif
            endif
         case (3)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call BGUARD2L(fxy,nx,ny,nxe,nye)
               else
                  call BGUARD2(fxy,nx,ny,nxe,nye)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LBGUARD2(fxy,nx,ny,nxe,nye)
               endif
            else if (ipbc==3) then
               if (order==LINEAR) then
                  call MBGUARD2L(fxy,nx,ny,nxe,nye)
               else
                  call MBGUARD2(fxy,nx,ny,nxe,nye)
               endif
            endif
         case (4:)
            ndim = size(fxy,1)
            if (ipbc==1) then
               if (order==LINEAR) then
                  call DCGUARD2L(fxy,nx,ny,nxe,nye,ndim)
               else
                  call DCGUARD2(fxy,nx,ny,nxe,nye,ndim)
               endif
            else if (ipbc==2) then
               if (order==QUADRATIC) then
                  call LDCGUARD2(fxy,nx,ny,nxe,nye,ndim)
               endif
!           else if (ipbc==3) then
            endif
         end select
         end subroutine icguard2p
!
         subroutine idguard2p(q,nx,ny,ipbc,inorder)
! copy guard cells for 2d scalar data
         implicit none
         integer :: nx, ny, ipbc
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (ipbc==1) then
            if (order==LINEAR) then
               call DGUARD2L(q,nx,ny,nxe,nye)
            else
               call DGUARD2(q,nx,ny,nxe,nye)
            endif
         else if (ipbc==2) then
            if (order==QUADRATIC) then
               call LDGUARD2(q,nx,ny,nxe,nye)
            endif
         else if (ipbc==3) then
            if (order==LINEAR) then
               call MDGUARD2L(q,nx,ny,nxe,nye)
            else
               call MDGUARD2(q,nx,ny,nxe,nye)
            endif
         endif
         end subroutine idguard2p
!
      end module npfield2d
