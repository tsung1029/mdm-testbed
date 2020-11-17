!-----------------------------------------------------------------------
!
      module nfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library nfield2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: april 21, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: sglvcs
      public :: dblcos, poisn_init, poisn, cuperpn, bpoisn, ibpoisn
      public :: maxweln, cmfieldn, emfieldn, cpfieldn, avpotn
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine DBLCOS2C(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(2,nx2v,ny2) :: cu2
         end subroutine
      end interface
      interface
         subroutine DBLCOS2D(q,q2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,ny2) :: q2
         end subroutine
      end interface
      interface
         subroutine DBLCOS2B(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx2v,ny2) :: cu2
         end subroutine
      end interface
      interface
         subroutine SGLVCS2C(cu,cu1,nx,ny,ncsx,ncsy,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, ncsx, ncsy, nxv, nyv, nx2v
!        real, dimension(* :: cu
         real :: cu
         real, dimension(2,nx2v,nyv) :: cu1
         end subroutine
      end interface
      interface
         subroutine SGLVCS2D(q,q1,nx,ny,ncs,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, ncs, nxv, nyv, nx2v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,nyv) :: q1
         end subroutine
      end interface
      interface
         subroutine SGLVCS2B(cu,cu1,nx,ny,ncsx,ncsy,ncsz,nxv,nyv,nx2v)
         implicit none
         integer :: nx, ny, ncsx, ncsy, ncsz, nxv, nyv, nx2v
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx2v,nyv) :: cu1
         end subroutine
      end interface
      interface
         subroutine POISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,n&
     &y2d)
         implicit none
         integer :: isign, nx, ny, nx2v, ny2d
         real :: ax, ay, affp, we
         real, dimension(nx2v,ny2d) :: q, fx, fy
         complex, dimension(nx2v/2,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine POISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2&
     &d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, we
         real, dimension(2*nxv,ny2d) :: q
         real, dimension(2,2*nxv,ny2d) :: fxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine POISNX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2&
     &d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, we
         real, dimension(2*nxv,ny2d) :: q
         real, dimension(3,2*nxv,ny2d) :: fxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine CUPERPNX2(cu,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real, dimension(3,2*nxv,ny2d) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPN2(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPNX22(cu,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real, dimension(2,2*nxv,ny2d) :: cu
         end subroutine
      end interface
      interface
         subroutine BPOISNX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nx&
     &v,ny2d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, ci, wm
         real, dimension(3,2*nxv,ny2d) :: cu, bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine BPOISN23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxe&
     &,nye,nxv)
         implicit none
         integer :: isign, nx, ny, nxe, nye, nxv
         real :: ax, ay, affp, ci, wm
         real, dimension(3,nxe,nye) :: cu, bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine BPOISNX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
     &,nxv,ny2d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, ci, wm
         real, dimension(2,2*nxv,ny2d) :: cu, bxy
         real, dimension(2*nxv,ny2d) :: bz
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine IBPOISNX23(cu,bxy,ffd,ci,wm,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real :: ci, wm
         real, dimension(3,2*nxv,ny2d) :: cu, bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine IBPOISN23(cu,bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
         implicit none
         integer :: nx, ny, nxe, nye, nxv
         real :: ci, wm
         real, dimension(3,nxe,nye) :: cu
         complex, dimension(3,nxe/2,nye) :: bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine MAXWELNX2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real :: ci, dt, wf, wm
         real, dimension(3,2*nxv,ny2d) :: exy, bxy, cu
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine MAXWELN2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxe,nye,nx&
     &v)
         implicit none
         integer :: nx, ny, nxe, nye, nxv
         real :: ci, dt, wf, wm
         complex, dimension(3,nxe/2,nye) :: exy, bxy
         real, dimension(3,nxe,nye) :: cu
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine CMFIELDN2(cu2,cu,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: nx, ny, nxv, ny2d, nxe, nye
         real, dimension(3,2*nxv,ny2d) :: cu2
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine EMFIELDN2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d, nxe, nye
         real, dimension(3,2*nxv,ny2d) :: fxy
         complex, dimension(3,nxe/2,nye) :: exy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine CPFIELDN2(fxy,exy,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: nx, ny, nxv, ny2d, nxe, nye
         real, dimension(3,2*nxv,ny2d) :: fxy
         real, dimension(3,nxe,nye) :: exy
         end subroutine
      end interface
      interface
         subroutine AVPOTNX23(bxy,axy,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         complex, dimension(3,nxv,ny2d) :: bxy, axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface dblcos
         module procedure idblcos2c
         module procedure idblcos2d
      end interface
!
      interface sglvcs
         module procedure isglvcs2c
         module procedure isglvcs2d
      end interface
!
      interface poisn_init
         module procedure ipoisn22init
      end interface
!
      interface poisn
         module procedure ipoisn2
         module procedure ispoisn2
         module procedure ipoisn22
      end interface
!
!     interface poisd3
!        module procedure ipoisd23
!     end interface
!
      interface cuperpn
         module procedure icuperpn2
      end interface
!
      interface bpoisn
         module procedure jbpoisn23
      end interface
!
      interface ibpoisn
         module procedure jibpoisn23
      end interface
!
      interface maxweln
         module procedure imaxweln2
      end interface
!
      interface cmfieldn
         module procedure icmfieldn2
!        module procedure idmfieldd2
      end interface
!
      interface emfieldn
         module procedure iemfieldn2
      end interface
!
      interface cpfieldn
         module procedure icpfieldn2
      end interface
!
      interface avpotn
         module procedure iavpotn23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine idblcos2c(cu,cu2,nx,ny,inorder)
! double array in each dimension for 2d vector data
! for neumann boundary conditions
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
               call DBLCOS2B(cu(1,1,1),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call DBLCOS2B(cu(1,2,2),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call DBLCOS2C(cu(1,1,1),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call DBLCOS2C(cu(1,2,2),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         endif
         end subroutine idblcos2c
!
         subroutine idblcos2d(q,q2,nx,ny,inorder)
! double array in each dimension for 2d scalar data
! for neumann boundary conditions
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
            call DBLCOS2D(q(1,1),q2,nx,ny,nxv,nyv,nx2v,ny2)
         else
            call DBLCOS2D(q(2,2),q2,nx,ny,nxv,nyv,nx2v,ny2)
         endif
         end subroutine idblcos2d
!
         subroutine isglvcs2c(cu,cu1,nx,ny,ncsx,ncsy,ncsz,inorder)
! double array in x dimension for 2d vector data
! for various boundary condition
         implicit none
         integer :: nx, ny, ncsx, ncsy, ncsz
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
               call SGLVCS2B(cu(1,1,1),cu1,nx,ny,ncsx,ncsy,ncsz,nxv,nyv,&
     &nx2v)
            else
               call SGLVCS2B(cu(1,2,2),cu1,nx,ny,ncsx,ncsy,ncsz,nxv,nyv,&
     &nx2v)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call SGLVCS2C(cu(1,1,1),cu1,nx,ny,ncsx,ncsy,nxv,nyv,nx2v)
            else
               call SGLVCS2C(cu(1,2,2),cu1,nx,ny,ncsx,ncsy,nxv,nyv,nx2v)
            endif
         endif
         end subroutine isglvcs2c
!
         subroutine isglvcs2d(q,q1,nx,ny,ncs,inorder)
! double array in x dimension for 2d scalar data
! for various boundary condition
         implicit none
         integer :: nx, ny, ncs
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
            call SGLVCS2D(q(1,1),q1,nx,ny,ncs,nxv,nyv,nx2v)
         else
            call SGLVCS2D(q(2,2),q1,nx,ny,ncs,nxv,nyv,nx2v)
         endif
         end subroutine isglvcs2d
!
         subroutine ipoisn2(q,fx,ffd,we,nx,ny)
! poisson solver for 2d potential, neumann boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 1, nx2v, ny2d
         real :: ax, ay, affp
         real, dimension(1,1) :: fy
         nx2v = size(q,1); ny2d = size(q,2)
         call POISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,ny2d)
         end subroutine ipoisn2
!
         subroutine ispoisn2(q,fy,ffd,nx,ny)
! smoother for 2d scalar field, neumann boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 2, nx2v, ny2d
         real :: ax, ay, affp, we
         real, dimension(1,1) :: fx
         nx2v = size(q,1); ny2d = size(q,2)
         call POISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,ny2d)
         end subroutine ispoisn2
!
         subroutine ipoisn22init(ffd,ax,ay,affp,nx,ny)
! initialize 2d electric field solver, neumann boundaries
         implicit none
         integer :: nx, ny
         real :: ax, ay, affp
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 0, nxv, ny2d = 1
         real :: we
         real, dimension(1,1) :: q
         real, dimension(2,1,1) :: fxy
         nxv = size(ffd,1)
         call POISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
         end subroutine ipoisn22init
!
         subroutine ipoisn22(q,fxy,ffd,we,nx,ny)
! poisson solver for 2d electric field, neumann boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = -1, nxv, ny2d
         real :: ax, ay, affp
         nxv = size(q,1)/2; ny2d = size(q,2)
         if (size(fxy,1)==2) then
            call POISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
         else if (size(fxy,1)==3) then
            call POISNX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
         endif
         end subroutine ipoisn22
!
!        subroutine ipoisd23(q,fxy,ffd,we,nx,ny)
! poisson solver for 2-1/2d electric field, neumann boundaries
!        implicit none
!        integer :: nx, ny
!        real :: we
!        real, dimension(:,:), pointer :: q
!        real, dimension(:,:,:), pointer :: fxy
!        complex, dimension(:,:), pointer :: ffd
! local data
!        integer :: isign = -1, nxe, nye, nxv
!        real :: ax, ay, affp
!        nxe = size(q,1); nye = size(q,2); nxv = size(ffd,1)
!        call POISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye,nxv)
!        end subroutine ipoisd23
!
         subroutine icuperpn2(cu,nx,ny)
! calculates transverse part of 2d vector field, neumann boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nye
         nxe = size(cu,2); nye = size(cu,3)
         call CUPERPN2(cu,nx,ny,nxe,nye)
         end subroutine icuperpn2
!
         subroutine jbpoisn23(cu,bxy,ffd,ci,wm,nx,ny)
! calculates static vector potential for 2d vector field,
! neumann boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 1, nxe, nye, nxv
         real :: ax, ay, affp
         nxe = size(cu,2); nye = size(cu,3)
         nxv = size(ffd,1)
         call BPOISN23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxe,nye,n&
     &xv)
         end subroutine jbpoisn23
!
         subroutine jibpoisn23(cu,bxy,ffd,ci,wm,nx,ny)
! calculates static magnetic field for periodic 2d vector field
! neumann boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: bxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxe, nye, nxv
         nxe = size(cu,2); nye = size(cu,3)
         nxv = size(ffd,1)
         call IBPOISN23(cu,bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
         end subroutine jibpoisn23
!
         subroutine imaxweln2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny)
! calculates maxwell's equation for 2d vector field,
! neumann boundaries
         implicit none
         integer :: nx, ny
         real :: ci, dt, wf, wm
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxe, nye, nxv
         nxe = size(cu,2); nye = size(cu,3); nxv = size(ffd,1)
         call MAXWELN2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxe,nye,nxv)
         end subroutine imaxweln2    
!
         subroutine icmfieldn2(cu2,cu,nx,ny)
! copies from double to normal array in y dimension for 2d vector data
! neumann boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu2, cu
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(cu2,2)/2; ny2d = size(cu2,3)
         nxe = size(cu,2); nye = size(cu,3)
         call CMFIELDN2(cu2,cu,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine icmfieldn2
!
!        subroutine idmfieldd2(q2,q,nx,ny)
! copies from double to normal array in y dimension for 2d scalar data
!        implicit none
!        integer :: nx, ny
!        real, dimension(:,:), pointer :: q2, q
! local data
!        integer :: nxv, ny2d, nxe, nye
!        nxv = size(q2,1)/2; ny2d = size(q2,2)
!        nxe = size(q,1); nye = size(q,2)
!        call DMFIELDD2(q2,q,nx,ny,nxv,ny2d,nxe,nye)
!        end subroutine idmfieldd2
!
         subroutine iemfieldn2(fxy,exy,ffd,isign,nx,ny)
! combines and smooths 2d vector fields, neumann boundaries
         implicit none
         integer :: isign, nx, ny
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(fxy,2)/2; ny2d = size(fxy,3)
         nxe = 2*size(exy,2); nye = size(exy,3)
         call EMFIELDN2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine iemfieldn2
!
         subroutine icpfieldn2(fxy,exy,nx,ny)
! combines 2d electric fields, neumann boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: fxy, exy
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(fxy,2)/2; ny2d = size(fxy,3)
         nxe = size(exy,2); nye = size(exy,3)
         call CPFIELDN2(fxy,exy,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine icpfieldn2
!
         subroutine iavpotn23(bxy,axy,nx,ny)
! calculates 2d vector potential from magnetic field
! neumann boundaries
         implicit none
         integer :: nx, ny
         complex, dimension(:,:,:), pointer :: bxy
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: nxe, nye
         nxe = size(axy,2); nye = size(axy,3)
         call AVPOTN23(bxy,axy,nx,ny,nxe,nye)
         end subroutine iavpotn23
!
      end module nfield2d
