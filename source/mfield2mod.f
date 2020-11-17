!-----------------------------------------------------------------------
!
      module mfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library mfield2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: april 21, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: sglcos, poismn_init, poismn, cuperpmn, bpoismn, ibpoismn
      public :: maxwelmn, cmfieldmn, emfieldmn, cpfieldmn, avpotmn
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine SGLCOS2C(cu,cu1,nx,ny,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v
!        real, dimension(* :: cu
         real :: cu
         real, dimension(2,nx2v,nyv) :: cu1
         end subroutine
      end interface
      interface
         subroutine SGLCOS2D(q,q1,nx,ny,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,nyv) :: q1
         end subroutine
      end interface
      interface
         subroutine SGLCOS2B(cu,cu1,nx,ny,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx2v,nyv) :: cu1
         end subroutine
      end interface
      interface
         subroutine POISMNX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,&
     &nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nx2v, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(nx2v,nyv) :: q, fx, fy
         complex, dimension(nx2v/2,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISMNX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,ny&
     &v,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxv,nyv) :: q
         real, dimension(2,2*nxv,nyv) :: fxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine POISMNX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,ny&
     &v,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nxv,nyv) :: q
         real, dimension(3,2*nxv,nyv) :: fxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine CUPERPMNX2(cu,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(3,2*nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPMN2(cu,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
         real, dimension(3,nxe,2*nyeh) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPMNX22(cu,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(2,2*nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine BPOISMNX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,n&
     &xv,nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(3,2*nxv,nyv) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine BPOISMN23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nx&
     &e,nyeh,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(3,nxe,2*nyeh) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine BPOISMNX22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,n&
     &y,nxv,nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(2,2*nxv,nyv) :: cu, bxy
         real, dimension(2*nxv,nyv) :: bz
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine IBPOISMNX23(cu,bxy,ffb,ci,wm,nx,ny,nxv,nyv,nyhd)
         implicit none
         integer :: nx, ny, nxv, nyv, nyhd
         real :: ci, wm
         real, dimension(3,2*nxv,nyv) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine IBPOISMN23(cu,bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,nyhd)
         implicit none
         integer :: nx, ny, nxe, nyeh, nxv, nyhd
         real :: ci, wm
         real, dimension(3,nxe,2*nyeh) :: cu
         complex, dimension(3,nxe/2,2*nyeh) :: bxy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine MAXWELMNX2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxv,nyv,&
     &nyhd)
         implicit none
         integer :: nx, ny, nxv, nyv, nyhd
         real :: ci, dt, wf, wm
         real, dimension(3,2*nxv,nyv) :: exy, bxy, cu
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine MAXWELMN2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxe,nyeh,&
     &nxv,nyhd)
         implicit none
         integer :: nx, ny, nxe, nyeh, nxv, nyhd
         real :: ci, dt, wf, wm
         complex, dimension(3,nxe/2,2*nyeh) :: exy, bxy
         real, dimension(3,nxe,2*nyeh) :: cu
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine CMFIELDMN2(cu1,cu,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(3,2*nxv,nyv) :: cu1
         real, dimension(3,nxe,2*nyeh) :: cu
         end subroutine
      end interface
      interface
         subroutine EMFIELDMN2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,&
     &nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nxe, nyeh, nyhd
         real, dimension(3,2*nxv,nyv) :: fxy
         complex, dimension(3,nxe/2,2*nyeh) :: exy
         complex, dimension(nxv,nyhd) :: ffb
         end subroutine
      end interface
      interface
         subroutine CPFIELDMN2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(3,2*nxv,nyv) :: fxy
         real, dimension(3,nxe,2*nyeh) :: exy
         end subroutine
      end interface
      interface
         subroutine AVPOTMNX23(bxy,axy,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(3,2*nxv,nyv) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine AVPOTMN23(bxy,axy,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
         complex, dimension(3,nxe/2,2*nyeh) :: bxy
         real, dimension(3,nxe,2*nyeh) :: axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface sglcos
         module procedure isglcos2c
         module procedure isglcos2d
      end interface
!
       interface poismn_init
         module procedure ipoismn22init
      end interface
!
      interface poismn
         module procedure ipoismn2
         module procedure ispoismn2
         module procedure ipoismn22
      end interface
!
!     interface poism3
!        module procedure ipoism23
!     end interface
!
      interface cuperpmn
         module procedure icuperpmn2
      end interface
!
      interface bpoismn
         module procedure jbpoismn23
      end interface
!
      interface ibpoismn
         module procedure jibpoismn23
      end interface
!
      interface maxwelmn
         module procedure imaxwelmn2
      end interface
!
      interface cmfieldmn
         module procedure icmfieldmn2
      end interface
!
      interface emfieldmn
         module procedure iemfieldmn2
      end interface
!
      interface cpfieldmn
         module procedure icpfieldmn2
      end interface
!
      interface avpotmn
         module procedure iavpotmn23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine isglcos2c(cu,cu1,nx,ny,inorder)
! double array in x dimension for 2d vector data
! for mixed periodic/neumann boundary condition
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
               call SGLCOS2B(cu(1,1,1),cu1,nx,ny,nxv,nyv,nx2v)
            else
               call SGLCOS2B(cu(1,2,2),cu1,nx,ny,nxv,nyv,nx2v)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call SGLCOS2C(cu(1,1,1),cu1,nx,ny,nxv,nyv,nx2v)
            else
               call SGLCOS2C(cu(1,2,2),cu1,nx,ny,nxv,nyv,nx2v)
            endif
         endif
         end subroutine isglcos2c
!
         subroutine isglcos2d(q,q1,nx,ny,inorder)
! double array in x dimension for 2d scalar data
! for mixed periodic/neumann boundary conditions
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
            call SGLCOS2D(q(1,1),q1,nx,ny,nxv,nyv,nx2v)
         else
            call SGLCOS2D(q(2,2),q1,nx,ny,nxv,nyv,nx2v)
         endif
         end subroutine isglcos2d
!
         subroutine ipoismn2(q,fx,ffb,we,nx,ny)
! poisson solver for 2d potential, mixed neumann/periodic boundaries
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
         call POISMNX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,nyv,ny&
     &hd)
         end subroutine ipoismn2
!
         subroutine ispoismn2(q,fy,ffb,nx,ny)
! smoother for 2d scalar field, mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 2, nx2v, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(1,1) :: fx
         nx2v = size(q,1); nyv = size(q,2); nyhd = size(ffb,2)
         call POISMNX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,nyv,ny&
     &hd)
         end subroutine ispoismn2
!
         subroutine ipoismn22init(ffb,ax,ay,affp,nx,ny)
! initialize 2d electric field solver,
! mixed neumann/periodic boundaries
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
         call POISMNX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,nyhd&
     &)
         end subroutine ipoismn22init
!
         subroutine ipoismn22(q,fxy,ffb,we,nx,ny)
! poisson solver for 2d electric field,
! mixed neumann/periodic boundaries
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
            call POISMNX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,n&
     &yhd)
         else if (size(fxy,1)==3) then
            call POISMNX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,n&
     &yhd)
         endif
         end subroutine ipoismn22
!
!        subroutine ipoism23(q,fxy,ffb,we,nx,ny)
! poisson solver for 2-1/2d electric field,
! mixed neumann/periodic boundaries
!        implicit none
!        integer :: nx, ny
!        real :: we
!        real, dimension(:,:), pointer :: q
!        real, dimension(:,:,:), pointer :: fxy
!        complex, dimension(:,:), pointer :: ffb
! local data
!        integer :: isign = -1, nxe, nyeh, nxv, nyhd
!        real :: ax, ay, affp
!        nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffb,1)
!        nyhd = size(ffb,2)
!        call POISM23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nyeh,nxv,n&
!    &yhd)
!        end subroutine ipoism23
!
         subroutine icuperpmn2(cu,nx,ny)
! calculates transverse part of 2d vector field,
! mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nyeh
         nxe = size(cu,2); nyeh = size(cu,3)/2
         call CUPERPMN2(cu,nx,ny,nxe,nyeh)
         end subroutine icuperpmn2
!
         subroutine jbpoismn23(cu,bxy,ffb,ci,wm,nx,ny)
! calculates static vector potential for 2d vector field,
! mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: isign = 1, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp
         nxe = size(cu,2); nyeh = size(cu,3)/2
         nxv = size(ffb,1); nyhd = size(ffb,2)
         call BPOISMN23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxe,nyeh&
     &,nxv,nyhd)
         end subroutine jbpoismn23
!
         subroutine jibpoismn23(cu,bxy,ffb,ci,wm,nx,ny)
! calculates static magnetic field for periodic 2d vector field
! mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: bxy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: nxe, nyeh, nxv, nyhd
         nxe = size(cu,2); nyeh = size(cu,3)/2
         nxv = size(ffb,1); nyhd = size(ffb,2)
         call IBPOISMN23(cu,bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,nyhd)
         end subroutine jibpoismn23
!
         subroutine imaxwelmn2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny)
! calculates maxwell's equation for 2d vector field,
! mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: ci, dt, wf, wm
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: nxe, nyeh, nxv, nyhd
         nxe = size(cu,2); nyeh = size(cu,3)/2; nxv = size(ffb,1)
         nyhd = size(ffb,2)
         call MAXWELMN2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxe,nyeh,nxv,ny&
     &hd)
         end subroutine imaxwelmn2    
!
         subroutine icmfieldmn2(cu1,cu,nx,ny)
! copies from double to normal array in y dimension for 2d vector data
! mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu1, cu
! local data
         integer :: nxv, nyv, nxe, nyeh
         nxv = size(cu1,2)/2; nyv = size(cu1,3)
         nxe = size(cu,2); nyeh = size(cu,3)/2
         call CMFIELDMN2(cu1,cu,nx,ny,nxv,nyv,nxe,nyeh)
         end subroutine icmfieldmn2
!
         subroutine iemfieldmn2(fxy,exy,ffb,isign,nx,ny)
! combines and smooths 2d vector fields,
! mixed neumann/periodic boundaries
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
         call EMFIELDMN2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,nyhd)
         end subroutine iemfieldmn2
!
         subroutine icpfieldmn2(fxy,exy,nx,ny)
! combines 2d electric fields, mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: fxy, exy
! local data
         integer :: nxv, nyv, nxe, nyeh
         nxv = size(fxy,2)/2; nyv = size(fxy,3)
         nxe = size(exy,2); nyeh = size(exy,3)/2
         call CPFIELDMN2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
         end subroutine icpfieldmn2
!
         subroutine iavpotmn23(bxy,axy,nx,ny)
! calculates 2d vector potential from magnetic field
! mixed neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         complex, dimension(:,:,:), pointer :: bxy
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: nxe, nyeh
         nxe = size(axy,2); nyeh = size(axy,3)/2
         call AVPOTMN23(bxy,axy,nx,ny,nxe,nyeh)
         end subroutine iavpotmn23
!
      end module mfield2d
