!-----------------------------------------------------------------------
!
      module bfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library bfield2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: april 22, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: MCGUARD2, MDGUARD2, MBGUARD2, MSCGUARD2, MSCGUARD22
      public :: MSGUARD2, MACGUARD2, MACGUARD22, MAGUARD2, MCGUARD2L
      public :: MDGUARD2L, MBGUARD2L, MSCGUARD2L, MSCGUARD22L, MSGUARD2L
      public :: maguard, mcguard
      public :: sglsin, hafsgl, poism_init, poismx, poism, cuperpm
      public :: bpoism, ibpoism, maxwelm, cmfieldm, emfieldm, cpfieldm
      public :: avpotm
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MCGUARD2(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         end subroutine
      end interface
      interface
         subroutine MDGUARD2(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine MBGUARD2(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: bxy
         end subroutine
      end interface
      interface
         subroutine MSCGUARD2(cu,xj0,yj0,zj0,nx,ny,ngx,nxe,nye)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ny, ngx, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine MSCGUARD22(cu,xj0,yj0,nx,ny,ngx,nxe,nye)
         implicit none
         real :: xj0, yj0
         integer :: nx, ny, ngx, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine MSGUARD2(q,qi0,nx,ny,ngx,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, ngx, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine MACGUARD2(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: cu
         end subroutine
      end interface
      interface
         subroutine MACGUARD22(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: cu
         end subroutine
      end interface
      interface
         subroutine MAGUARD2(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: q
         end subroutine
      end interface
      interface
         subroutine MCGUARD2L(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         end subroutine
      end interface
      interface
         subroutine MDGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine MBGUARD2L(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: fxy
         end subroutine
      end interface
      interface
         subroutine MSCGUARD2L(cu,xj0,yj0,zj0,nx,ny,ngx,nxe,nye)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ny, ngx, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine MSCGUARD22L(cu,xj0,yj0,nx,ny,ngx,nxe,nye)
         implicit none
         real :: xj0, yj0
         integer :: nx, ny, ngx, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine MSGUARD2L(q,qi0,nx,ny,ngx,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, ngx, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine MACGUARD2L(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine MACGUARD22L(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine MAGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine SGLSIN2C(cu,cu1,nx,ny,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v
!        real, dimension(* :: cu
         real :: cu
         real, dimension(2,nx2v,nyv) :: cu1
         end subroutine
      end interface
      interface
         subroutine SGLSIN2D(q,q1,nx,ny,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,nyv) :: q1
         end subroutine
      end interface
      interface
         subroutine SGLSIN2B(cu,cu1,nx,ny,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx2v,nyv) :: cu1
         end subroutine
      end interface
      interface
         subroutine HAFSGL2C(fxy,fxy1,nx,ny,nxe,nye,nx2v)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v
!        real, dimension(*) :: fxy
         real :: fxy
         real, dimension(2,nx2v,nye) :: fxy1
         end subroutine
      end interface
      interface
         subroutine HAFSGL2D(q,q1,nx,ny,nxe,nye,nx2v)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,nye) :: q1
         end subroutine
      end interface
      interface
         subroutine HAFSGL2B(bxy,bxy1,nx,ny,nxe,nye,nx2v)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v
!        real, dimension(*) :: bxy
         real :: bxy
         real, dimension(3,nx2v,nye) :: bxy1
         end subroutine
      end interface
      interface
         subroutine POISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,n&
     &yv,nyhd)
         implicit none
         integer :: isign, nx, ny, nx2v, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(nx2v,nyv) :: q, fx, fy
         complex, dimension(nx2v/2,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISM2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nye&
     &h,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fx, fy
         real :: q, fx, fy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv&
     &,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxv,nyv) :: q
         real, dimension(2,2*nxv,nyv) :: fxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISM22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nyeh&
     &,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISMX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv&
     &,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxv,nyv) :: q
         real, dimension(3,2*nxv,nyv) :: fxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISM23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nyeh&
     &,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine CUPERPMX2(cu,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(3,2*nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPM2(cu,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPMX22(cu,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(2,2*nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPM22(cu,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine BPOISMX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nx&
     &v,nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(3,2*nxv,nyv) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine BPOISM23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxe&
     &,nyeh,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy
         real :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine BPOISMX22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,ny&
     &,nxv,nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(2,2*nxv,nyv) :: cu, bxy
         real, dimension(2*nxv,nyv) :: bz
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface 
         subroutine BPOISM22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,ny,&
     &nxe,nyeh,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy, bz
         real :: cu, bxy, bz
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine IBPOISMX23(cu,bxy,ffb,ci,wm,nx,ny,nxv,nyv,nyhd)
         implicit none
         integer :: nx, ny, nxv, nyv, nyhd
         real :: ci, wm
         real, dimension(3,2*nxv,nyv) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine IBPOISM23(cu,bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,nyhd)
         implicit none
         integer :: nx, ny, nxe, nyeh, nxv, nyhd
         real :: ci, wm
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(3,nxe/2,2*nyeh) :: bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine MAXWELMX2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxv,nyv,n&
     &yhd)
         implicit none
         integer :: nx, ny, nxv, nyv, nyhd
         real :: ci, dt, wf, wm
         real, dimension(3,2*nxv,nyv) :: exy, bxy, cu
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine MAXWELM2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxe,nyeh,n&
     &xv,nyhd)
         implicit none
         integer :: nx, ny, nxe, nyeh, nxv, nyhd
         real :: ci, dt, wf, wm
         complex, dimension(3,nxe/2,2*nyeh) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine DMFIELDM2(q1,q,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(2*nxv,nyv) :: q1
         real, dimension(nxe,2*nyeh) :: q
         end subroutine
      end interface
      interface
         subroutine CMFIELDM2(cu1,cu,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(3,2*nxv,nyv) :: cu1
         real, dimension(3,nxe,2*nyeh) :: cu
         end subroutine
      end interface
      interface
         subroutine EMFIELDM2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,n&
     &yhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nxe, nyeh, nyhd
         real, dimension(3,2*nxv,nyv) :: fxy
         complex, dimension(3,nxe/2,2*nyeh) :: exy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine PMFIELDM2(pot1,pot,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(2*nxv,nyv) :: pot1
         real, dimension(nxe,2*nyeh) :: pot
         end subroutine
      end interface
      interface
         subroutine CPFIELDM2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(3,2*nxv,nyv) :: fxy
         real, dimension(3,nxe,2*nyeh) :: exy
         end subroutine
      end interface
      interface
         subroutine AVPOTMX23(bxy,axy,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(3,2*nxv,nyv) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine AVPOTM23(bxy,axy,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
         complex, dimension(3,nxe/2,2*nyeh) :: bxy
!        real, dimension(* :: axy
         real :: axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!      
      interface maguard
         module procedure imaguard2
      end interface
!
      interface mcguard
         module procedure imcguard2
         module procedure imdguard2
      end interface
!
      interface sglsin
         module procedure isglsin2c
         module procedure isglsin2d
      end interface
!
      interface hafsgl
         module procedure ihafsgl2c
         module procedure ihafsgl2d
      end interface
!
       interface poism_init
         module procedure ipoism22init
      end interface
!
      interface poismx
         module procedure ipoismx2
         module procedure ispoismx2
         module procedure ipoismx22
      end interface
!
      interface poism
         module procedure ipoism2
         module procedure ispoism2
         module procedure ipoism23
      end interface
!
      interface cuperpm
         module procedure icuperpm2
      end interface
!
      interface bpoism
         module procedure jbpoism23
      end interface
!
      interface ibpoism
         module procedure jibpoism23
      end interface
!
      interface maxwelm
         module procedure imaxwelm2
      end interface
!
      interface cmfieldm
         module procedure icmfieldm2
         module procedure idmfieldm2
      end interface
!
      interface emfieldm
         module procedure iemfieldm2
      end interface
!
      interface cpfieldm
         module procedure icpfieldm2
      end interface
!
      interface avpotm
         module procedure iavpotm23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imaguard2(q,nx,ny,inorder)
! add guard cells for 2d scalar data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call MAGUARD2L(q,nx,ny,nxe,nye)
         else
            call MAGUARD2(q(2,1),nx-2,ny,nxe,nye)
         endif
         end subroutine imaguard2
!
         subroutine imcguard2(fxy,nx,ny,inorder)
! copy guard cells for 2d vector data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: nxe, nye, order
         nxe = size(fxy,2); nye = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(fxy,1)==2) then
            if (order==LINEAR) then
               call MCGUARD2L(fxy,nx,ny,nxe,nye)
            else
               call MCGUARD2(fxy,nx,ny,nxe,nye)
            endif
         else if (size(fxy,1)==3) then
            if (order==LINEAR) then
               call MBGUARD2L(fxy,nx,ny,nxe,nye)
            else
               call MBGUARD2(fxy,nx,ny,nxe,nye)
            endif
         endif
         end subroutine imcguard2
!
         subroutine imdguard2(q,nx,ny,inorder)
! copy guard cells for 2d scalar data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call MDGUARD2L(q,nx,ny,nxe,nye)
         else
            call MDGUARD2(q,nx,ny,nxe,nye)
         endif
         end subroutine imdguard2
!
         subroutine isglsin2c(cu,cu1,nx,ny,inorder)
! double array in x dimension for 2d vector data
! for mixed periodic/dirichlet boundary condition
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
         real, dimension(:,:,:), pointer :: cu1
! local data
         integer :: nxv, nyv, nx2v, order
         nxv = size(cu,2);  nyv = size(cu,3)
         nx2v = size(cu1,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call SGLSIN2B(cu(1,1,1),cu1,nx,ny,nxv,nyv,nx2v)
            else
               call SGLSIN2B(cu(1,2,2),cu1,nx,ny,nxv,nyv,nx2v)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call SGLSIN2C(cu(1,1,1),cu1,nx,ny,nxv,nyv,nx2v)
            else
               call SGLSIN2C(cu(1,2,2),cu1,nx,ny,nxv,nyv,nx2v)
            endif
         endif
         end subroutine isglsin2c
!
         subroutine isglsin2d(q,q1,nx,ny,inorder)
! double array in x dimension for 2d scalar data
! for mixed periodic/dirichlet boundary conditions
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
         real, dimension(:,:), pointer :: q1
! local data
         integer :: nxv, nyv, nx2v, order
         nxv = size(q,1);  nyv = size(q,2)
         nx2v = size(q1,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call SGLSIN2D(q(1,1),q1,nx,ny,nxv,nyv,nx2v)
         else
            call SGLSIN2D(q(2,2),q1,nx,ny,nxv,nyv,nx2v)
         endif
         end subroutine isglsin2d
!
         subroutine ihafsgl2c(fxy,fxy1,nx,ny,inorder)
! copy from double to normal array in x dimension for 2d vector data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         real, dimension(:,:,:), pointer :: fxy1
! local data
         integer :: nxe, nye, nx2v, order
         nxe = size(fxy,2);  nye = size(fxy,3)
         nx2v = size(fxy1,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(fxy,1)==2) then
            if (order==LINEAR) then
               call HAFSGL2C(fxy(1,1,1),fxy1,nx,ny,nxe,nye,nx2v)
            else
               call HAFSGL2C(fxy(1,2,2),fxy1,nx,ny,nxe,nye,nx2v)
            endif
         else if (size(fxy,1)==3) then
            if (order==LINEAR) then
               call HAFSGL2B(fxy(1,1,1),fxy1,nx,ny,nxe,nye,nx2v)
            else
               call HAFSGL2B(fxy(1,2,2),fxy1,nx,ny,nxe,nye,nx2v)
            endif
         endif
         end subroutine ihafsgl2c
!
         subroutine ihafsgl2d(q,q1,nx,ny,inorder)
! copy from double to normal array in x dimension for 2d scalar data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
         real, dimension(:,:), pointer :: q1
! local data
         integer :: nxe, nye, nx2v, order
         nxe = size(q,1);  nye = size(q,2)
         nx2v = size(q1,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call HAFSGL2D(q(1,1),q1,nx,ny,nxe,nye,nx2v)
         else
            call HAFSGL2D(q(2,2),q1,nx,ny,nxe,nye,nx2v)
         endif
         end subroutine ihafsgl2d
!
         subroutine ipoismx2(q,fx,ffb,we,nx,ny)
! poisson solver for 2d potential, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 1, nx2v, nyv, nyhd
         real :: ax, ay, affp
         real, dimension(1,1) :: fy
         nx2v = size(q,1); nyv = size(q,2); nyhd = size(ffb,2)
         call POISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,nyv,nyh&
     &d)
         end subroutine ipoismx2
!
         subroutine ispoismx2(q,fy,ffb,nx,ny)
! smoother for 2d scalar field, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 2, nx2v, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(1,1) :: fx
         nx2v = size(q,1); nyv = size(q,2); nyhd = size(ffb,2)
         call POISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,nyv,nyh&
     &d)
         end subroutine ispoismx2
!
         subroutine ipoism22init(ffb,ax,ay,affp,nx,ny)
! initialize 2d electric field solver,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: ax, ay, affp
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 0, nxv, nyv = 1, nyhd
         real :: we
         real, dimension(1,1) :: q
         real, dimension(2,1,1) :: fxy
         nxv = size(ffb,1); nyhd = size(ffb,2)
         call POISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,nyhd)
         end subroutine ipoism22init
!
         subroutine ipoismx22(q,fxy,ffb,we,nx,ny)
! poisson solver for 2d electric field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = -1, nxv, nyv, nyhd
         real :: ax, ay, affp
         nxv = size(q,1)/2; nyv = size(q,2); nyhd = size(ffb,2)
         if (size(fxy,1)==2) then
            call POISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,ny&
     &hd)
         else if (size(fxy,1)==3) then
            call POISMX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,ny&
     &hd)
         endif
         end subroutine ipoismx22
!
         subroutine ipoism2(q,fx,ffb,we,nx,ny,order)
! poisson solver for 2d potential, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 1, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp
         real, dimension(2,2) :: fy
         nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffb,1)
         nyhd = size(ffb,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISM2(q(1,1),fx(1,1),fy(1,1),isign,ffb,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
         else
            call POISM2(q(2,2),fx(2,2),fy(2,2),isign,ffb,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
         endif
         end subroutine ipoism2
!
         subroutine ispoism2(q,fy,ffb,nx,ny,order)
! smoother for 2d scalar field, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 2, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp, we
         real, dimension(2,2) :: fx
         nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffb,1)
         nyhd = size(ffb,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISM2(q(1,1),fx(1,1),fy(1,1),isign,ffb,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
         else
            call POISM2(q(2,2),fx(2,2),fy(2,2),isign,ffb,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
         endif
         end subroutine ispoism2
!
         subroutine ipoism23(q,fxy,ffb,we,nx,ny,order)
! poisson solver for 2-1/2d electric field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = -1, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp
         nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffb,1)
         nyhd = size(ffb,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (size(fxy,1)==2) then
            if (inorder==LINEAR) then
               call POISM22(q(1,1),fxy(1,1,1),isign,ffb,ax,ay,affp,we,nx&
     &,ny,nxe,nyeh,nxv,nyhd)
            else
               call POISM22(q(2,2),fxy(1,2,2),isign,ffb,ax,ay,affp,we,nx&
     &,ny,nxe,nyeh,nxv,nyhd)
            endif
         else if (size(fxy,1)==3) then
            if (inorder==LINEAR) then
               call POISM23(q(1,1),fxy(1,1,1),isign,ffb,ax,ay,affp,we,nx&
     &,ny,nxe,nyeh,nxv,nyhd)
            else
               call POISM23(q(2,2),fxy(1,2,2),isign,ffb,ax,ay,affp,we,nx&
     &,ny,nxe,nyeh,nxv,nyhd)
            endif
         endif
         end subroutine ipoism23
!
         subroutine icuperpm2(cu,nx,ny,order)
! calculates transverse part of 2d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nyeh, inorder
         nxe = size(cu,2); nyeh = size(cu,3)/2
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (size(cu,1)==2) then
            if (inorder==LINEAR) then
               call CUPERPM22(cu(1,1,1),nx,ny,nxe,nyeh)
            else
               call CUPERPM22(cu(1,2,2),nx,ny,nxe,nyeh)
            endif
         else if (size(cu,1)==3) then
            if (inorder==LINEAR) then
               call CUPERPM2(cu(1,1,1),nx,ny,nxe,nyeh)
            else
               call CUPERPM2(cu(1,2,2),nx,ny,nxe,nyeh)
            endif
         endif
         end subroutine icuperpm2
!
         subroutine jbpoism23(cu,bxy,ffb,ci,wm,nx,ny,order)
! calculates static vector potential for 2d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 1, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp
         nxe = size(cu,2); nyeh = size(cu,3)/2
         nxv = size(ffb,1); nyhd = size(ffb,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call BPOISM23(cu(1,1,1),bxy(1,1,1),isign,ffb,ax,ay,affp,ci,w&
     &m,nx,ny,nxe,nyeh,nxv,nyhd)
         else
            call BPOISM23(cu(1,2,2),bxy(1,2,2),isign,ffb,ax,ay,affp,ci,w&
     &m,nx,ny,nxe,nyeh,nxv,nyhd)
         endif
         end subroutine jbpoism23
!
!                  call ibpois(cu,bxyz,ffc,ci,wm,nx,ny,inorder)
         subroutine jibpoism23(cu,bxy,ffb,ci,wm,nx,ny,order)
! calculates static magnetic field for periodic 2d vector field
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: bxy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: nxe, nyeh, nxv, nyhd, inorder
         nxe = size(cu,2); nyeh = size(cu,3)/2
         nxv = size(ffb,1); nyhd = size(ffb,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call IBPOISM23(cu(1,1,1),bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,ny&
     &hd)
         else
            call IBPOISM23(cu(1,2,2),bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,ny&
     &hd)
         endif
         end subroutine jibpoism23
!
         subroutine imaxwelm2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,order)
! calculates maxwell's equation for 2d vector field,
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, dt, wf, wm
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: nxe, nyeh, nxv, nyhd, inorder
         nxe = size(cu,2); nyeh = size(cu,3)/2; nxv = size(ffb,1)
         nyhd = size(ffb,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call MAXWELM2(exy,bxy,cu(1,1,1),ffb,ci,dt,wf,wm,nx,ny,nxe,ny&
     &eh,nxv,nyhd)
         else
            call MAXWELM2(exy,bxy,cu(1,2,2),ffb,ci,dt,wf,wm,nx,ny,nxe,ny&
     &eh,nxv,nyhd)
         endif
         end subroutine imaxwelm2    
!
         subroutine idmfieldm2(q1,q,nx,ny)
! copies from double to normal array in y dimension for 2d scalar data
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q1, q
! local data
         integer :: nxv, nyv, nxe, nyeh
         nxv = size(q1,1)/2; nyv = size(q1,2)
         nxe = size(q,1); nyeh = size(q,2)/2
         call DMFIELDM2(q1,q,nx,ny,nxv,nyv,nxe,nyeh)
         end subroutine idmfieldm2
!
         subroutine icmfieldm2(cu1,cu,nx,ny)
! copies from double to normal array in y dimension for 2d vector data
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu1, cu
! local data
         integer :: nxv, nyv, nxe, nyeh
         nxv = size(cu1,2)/2; nyv = size(cu1,3)
         nxe = size(cu,2); nyeh = size(cu,3)/2
         call CMFIELDM2(cu1,cu,nx,ny,nxv,nyv,nxe,nyeh)
         end subroutine icmfieldm2
!
         subroutine iemfieldm2(fxy,exy,ffb,isign,nx,ny)
! combines and smooths 2d vector fields,
! mixed conducting/periodic boundaries
         implicit none
         integer :: isign, nx, ny
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: nxv, nyv, nxe, nyeh, nyhd
         nxv = size(fxy,2)/2; nyv = size(fxy,3)
         nxe = 2*size(exy,2); nyeh = size(exy,3)/2
         nyhd = size(ffb,2)
         call EMFIELDM2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,nyhd)
         end subroutine iemfieldm2
!
         subroutine icpfieldm2(fxy,exy,nx,ny)
! combines 2d electric fields, mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: fxy, exy
! local data
         integer :: nxv, nyv, nxe, nyeh
         nxv = size(fxy,2)/2; nyv = size(fxy,3)
         nxe = size(exy,2); nyeh = size(exy,3)/2
         call CPFIELDM2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
         end subroutine icpfieldm2
!
         subroutine iavpotm23(bxy,axy,nx,ny,order)
! calculates 2d vector potential from magnetic field
! mixed conducting/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         complex, dimension(:,:,:), pointer :: bxy
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: nxe, nyeh, inorder
         nxe = size(axy,2); nyeh = size(axy,3)/2
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call AVPOTM23(bxy,axy(1,1,1),nx,ny,nxe,nyeh)
         else
            call AVPOTM23(bxy,axy(1,2,2),nx,ny,nxe,nyeh)
         endif
         end subroutine iavpotm23
!
      end module bfield2d
