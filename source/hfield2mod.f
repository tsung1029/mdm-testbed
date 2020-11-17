!-----------------------------------------------------------------------
!
      module hfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library hfield2lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: may 2, 2006
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: sgldsin, poismd_init, poismdx, poismd, cuperpmd
      public :: bpoismd, ibpoismd, maxwelmd, cmfieldmd, emfieldmd
      public :: cpfieldmd, avpotmd
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine SGLDSIN2C(cu,cu3,nx,ny,nxv,nyv,nx4v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx4v
!        real, dimension(* :: cu
         real :: cu
         real, dimension(2,nx4v,nyv) :: cu3
         end subroutine
      end interface
      interface
         subroutine SGLDSIN2D(q,q3,nx,ny,nxv,nyv,nx4v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx4v
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx4v,nyv) :: q3
         end subroutine
      end interface
      interface
         subroutine SGLDSIN2B(cu,cu3,nx,ny,nxv,nyv,nx4v)
         implicit none
         integer :: nx, ny, nxv, nyv, nx4v
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx4v,nyv) :: cu3
         end subroutine
      end interface
      interface
         subroutine POISMDX2(q,fx,fy,isign,ffh,ax,ay,affp,we,nx,ny,nx2v,&
     &nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nx2v, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(2*nx2v,nyv) :: q, fx, fy
         complex, dimension(nx2v/2,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine POISMD2(q,fx,fy,isign,ffh,ax,ay,affp,we,nx,ny,nxe,ny&
     &eh,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fx, fy
         real :: q, fx, fy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine POISMDX22(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxv,ny&
     &v,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(4*nxv,nyv) :: q
         real, dimension(2,4*nxv,nyv) :: fxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine POISMD22(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxe,nye&
     &h,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine POISMDX23(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxv,ny&
     &v,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(4*nxv,nyv) :: q
         real, dimension(3,4*nxv,nyv) :: fxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine POISMD23(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxe,nye&
     &h,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine CUPERPMDX2(cu,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(3,4*nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPMD2(cu,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPMDX22(cu,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(2,4*nxv,nyv) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPMD22(cu,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine BPOISMDX23(cu,bxy,isign,ffh,ax,ay,affp,ci,wm,nx,ny,n&
     &xv,nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(3,4*nxv,nyv) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine BPOISMD23(cu,bxy,isign,ffh,ax,ay,affp,ci,wm,nx,ny,nx&
     &e,nyeh,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy
         real :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine BPOISMDX22(cu,bxy,bz,isign,ffh,ax,ay,affp,ci,wm,nx,n&
     &y,nxv,nyv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nyhd
         real :: ax, ay, affp, ci, wm
         real, dimension(2,4*nxv,nyv) :: cu, bxy
         real, dimension(2*nxv,nyv) :: bz
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine BPOISMD22(cu,bxy,bz,isign,ffh,ax,ay,affp,ci,wm,nx,ny&
     &,nxe,nyeh,nxv,nyhd)
         implicit none
         integer :: isign, nx, ny, nxe, nyeh, nxv, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy, bz
         real :: cu, bxy, bz
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine IBPOISMDX23(cu,bxy,ffh,ci,wm,nx,ny,nxv,nyv,nyhd)
         implicit none
         integer :: nx, ny, nxv, nyv, nyhd
         real :: ci, wm
         real, dimension(3,4*nxv,nyv) :: cu, bxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine IBPOISMD23(cu,bxy,ffh,ci,wm,nx,ny,nxe,nyeh,nxv,nyhd)
         implicit none
         integer :: nx, ny, nxe, nyeh, nxv, nyhd
         real :: ci, wm
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(3,nxe,nyeh) :: bxy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine MAXWELMDX2(exy,bxy,cu,ffh,ci,dt,wf,wm,nx,ny,nxv,nyv,&
     &nyhd)
         implicit none
         integer :: nx, ny, nxv, nyv, nyhd
         real :: ci, dt, wf, wm
         real, dimension(3,4*nxv,nyv) :: exy, bxy, cu
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine MAXWELMD2(exy,bxy,cu,ffh,ci,dt,wf,wm,nx,ny,nxe,nyeh,&
     &nxv,nyhd)
         implicit none
         integer :: nx, ny, nxe, nyeh, nxv, nyhd
         real :: ci, dt, wf, wm
         complex, dimension(3,nxe,nyeh) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine DMFIELDMD2(q3,q,nx,ny,nx2v,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nx2v, nyv, nxe, nyeh
         real, dimension(2*nx2v,nyv) :: q3
         real, dimension(nxe,2*nyeh) :: q
         end subroutine
      end interface
      interface
         subroutine CMFIELDMD2(cu3,cu,nx,ny,nx2v,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nx2v, nyv, nxe, nyeh
         real, dimension(3,2*nx2v,nyv) :: cu3
         real, dimension(3,nxe,2*nyeh) :: cu
         end subroutine
      end interface
      interface
         subroutine EMFIELDMD2(fxy,exy,ffh,isign,nx,ny,nxv,nyv,nxe,nyeh,&
     &nyhd)
         implicit none
         integer :: isign, nx, ny, nxv, nyv, nxe, nyeh, nyhd
         real, dimension(3,4*nxv,nyv) :: fxy
         complex, dimension(3,nxe/2,2*nyeh) :: exy
         complex, dimension(nxv,nyhd) :: ffh
         end subroutine
      end interface
      interface
         subroutine PMFIELDMD2(pot3,pot,nx,ny,nx2v,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nx2v, nyv, nxe, nyeh
         real, dimension(2*nx2v,nyv) :: pot3
         real, dimension(nxe,2*nyeh) :: pot
         end subroutine
      end interface
      interface
         subroutine CPFIELDMD2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxv, nyv, nxe, nyeh
         real, dimension(3,4*nxv,nyv) :: fxy
         real, dimension(3,nxe,2*nyeh) :: exy
         end subroutine
      end interface
      interface
         subroutine AVPOTMDX23(bxy,axy,nx,ny,nxv,nyv)
         implicit none
         integer :: nx, ny, nxv, nyv
         real, dimension(3,4*nxv,nyv) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine AVPOTMD23(bxy,axy,nx,ny,nxe,nyeh)
         implicit none
         integer :: nx, ny, nxe, nyeh
         complex, dimension(3,nxe,nyeh) :: bxy
!        real, dimension(* :: axy
         real :: axy
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface sgldsin
         module procedure isgldsin2c
         module procedure isgldsin2d
      end interface
!
       interface poismd_init
         module procedure ipoismd22init
      end interface
!
      interface poismdx
         module procedure ipoismdx2
         module procedure ispoismdx2
         module procedure ipoismdx22
      end interface
!
      interface poismd
         module procedure ipoismd2
         module procedure ispoismd2
         module procedure ipoismd23
      end interface
!
      interface cuperpmd
         module procedure icuperpmd2
      end interface
!
      interface bpoismd
         module procedure jbpoismd23
      end interface
!
      interface ibpoismd
         module procedure jibpoismd23
      end interface
!
      interface maxwelmd
         module procedure imaxwelmd2
      end interface
!
      interface cmfieldmd
         module procedure icmfieldmd2
         module procedure idmfieldmd2
      end interface
!
      interface emfieldmd
         module procedure iemfieldmd2
      end interface
!
      interface cpfieldmd
         module procedure icpfieldmd2
      end interface
!
      interface avpotmd
         module procedure iavpotmd23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine isgldsin2c(cu,cu3,nx,ny,inorder)
! quadruple array in x dimension for 2d vector data
! for mixed periodic/dirichlet-neumann boundary condition
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
         real, dimension(:,:,:), pointer :: cu3
! local data
         integer :: nxv, nyv, nx4v, order
         nxv = size(cu,2);  nyv = size(cu,3)
         nx4v = size(cu3,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(cu,1)==3) then
            if (order==LINEAR) then
               call SGLDSIN2B(cu(1,1,1),cu3,nx,ny,nxv,nyv,nx4v)
            else
               call SGLDSIN2B(cu(1,2,2),cu3,nx,ny,nxv,nyv,nx4v)
            endif
         else if (size(cu,1)==2) then
            if (order==LINEAR) then
               call SGLDSIN2C(cu(1,1,1),cu3,nx,ny,nxv,nyv,nx4v)
            else
               call SGLDSIN2C(cu(1,2,2),cu3,nx,ny,nxv,nyv,nx4v)
            endif
         endif
         end subroutine isgldsin2c
!
         subroutine isgldsin2d(q,q3,nx,ny,inorder)
! quadruple array in x dimension for 2d scalar data
! for mixed periodic/dirichlet-neumann boundary conditions
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
         real, dimension(:,:), pointer :: q3
! local data
         integer :: nxv, nyv, nx4v, order
         nxv = size(q,1);  nyv = size(q,2)
         nx4v = size(q3,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call SGLDSIN2D(q(1,1),q3,nx,ny,nxv,nyv,nx4v)
         else
            call SGLDSIN2D(q(2,2),q3,nx,ny,nxv,nyv,nx4v)
         endif
         end subroutine isgldsin2d
!
         subroutine ipoismdx2(q,fx,ffh,we,nx,ny)
! poisson solver for 2d potential, mixed dirichlet-neumann/periodic
! boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = 1, nx2v, nyv, nyhd
         real :: ax, ay, affp
         real, dimension(1,1) :: fy
         nx2v = size(q,1)/2; nyv = size(q,2); nyhd = size(ffh,2)
         call POISMDX2(q,fx,fy,isign,ffh,ax,ay,affp,we,nx,ny,nx2v,nyv,ny&
     &hd)
         end subroutine ipoismdx2
!
         subroutine ispoismdx2(q,fy,ffh,nx,ny)
! smoother for 2d scalar field, mixed dirichlet-neumann/periodic
! boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = 2, nx2v, nyv, nyhd
         real :: ax, ay, affp, we
         real, dimension(1,1) :: fx
         nx2v = size(q,1)/2; nyv = size(q,2); nyhd = size(ffh,2)
         call POISMDX2(q,fx,fy,isign,ffh,ax,ay,affp,we,nx,ny,nx2v,nyv,ny&
     &hd)
         end subroutine ispoismdx2
!
         subroutine ipoismd22init(ffh,ax,ay,affp,nx,ny)
! initialize 2d electric field solver,
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: ax, ay, affp
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = 0, nxv, nyv = 1, nyhd
         real :: we
         real, dimension(1,1) :: q
         real, dimension(2,1,1) :: fxy
         nxv = size(ffh,1); nyhd = size(ffh,2)
         call POISMDX22(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxv,nyv,nyhd&
     &)
         end subroutine ipoismd22init
!
         subroutine ipoismdx22(q,fxy,ffh,we,nx,ny)
! poisson solver for 2d electric field,
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = -1, nxv, nyv, nyhd
         real :: ax, ay, affp
         nxv = size(q,1)/4; nyv = size(q,2); nyhd = size(ffh,2)
         if (size(fxy,1)==2) then
            call POISMDX22(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxv,nyv,n&
     &yhd)
         else if (size(fxy,1)==3) then
            call POISMDX23(q,fxy,isign,ffh,ax,ay,affp,we,nx,ny,nxv,nyv,n&
     &yhd)
         endif
         end subroutine ipoismdx22
!
         subroutine ipoismd2(q,fx,ffh,we,nx,ny,order)
! poisson solver for 2d potential, mixed dirichlet-neumann/periodic
! boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = 1, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp
         real, dimension(2,2) :: fy
         nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffh,1)
         nyhd = size(ffh,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISMD2(q(1,1),fx(1,1),fy(1,1),isign,ffh,ax,ay,affp,we,&
     &nx,ny,nxe,nyeh,nxv,nyhd)
         else
            call POISMD2(q(2,2),fx(2,2),fy(2,2),isign,ffh,ax,ay,affp,we,&
     &nx,ny,nxe,nyeh,nxv,nyhd)
         endif
         end subroutine ipoismd2
!
         subroutine ispoismd2(q,fy,ffh,nx,ny,order)
! smoother for 2d scalar field, mixed dirichlet-neumann/periodic
! boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = 2, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp, we
         real, dimension(2,2) :: fx
         nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffh,1)
         nyhd = size(ffh,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISMD2(q(1,1),fx(1,1),fy(1,1),isign,ffh,ax,ay,affp,we,&
     &nx,ny,nxe,nyeh,nxv,nyhd)
         else
            call POISMD2(q(2,2),fx(2,2),fy(2,2),isign,ffh,ax,ay,affp,we,&
     &nx,ny,nxe,nyeh,nxv,nyhd)
         endif
         end subroutine ispoismd2
!
         subroutine ipoismd23(q,fxy,ffh,we,nx,ny,order)
! poisson solver for 2-1/2d electric field,
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = -1, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp
         nxe = size(q,1); nyeh = size(q,2)/2; nxv = size(ffh,1)
         nyhd = size(ffh,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (size(fxy,1)==2) then
            if (inorder==LINEAR) then
               call POISMD22(q(1,1),fxy(1,1,1),isign,ffh,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
            else
               call POISMD22(q(2,2),fxy(1,2,2),isign,ffh,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
            endif
         else if (size(fxy,1)==3) then
            if (inorder==LINEAR) then
               call POISMD23(q(1,1),fxy(1,1,1),isign,ffh,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
            else
               call POISMD23(q(2,2),fxy(1,2,2),isign,ffh,ax,ay,affp,we,n&
     &x,ny,nxe,nyeh,nxv,nyhd)
            endif
         endif
         end subroutine ipoismd23
!
         subroutine icuperpmd2(cu,nx,ny,order)
! calculates transverse part of 2d vector field,
! mixed dirichlet-neumann/periodic boundaries
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
               call CUPERPMD22(cu(1,1,1),nx,ny,nxe,nyeh)
            else
               call CUPERPMD22(cu(1,2,2),nx,ny,nxe,nyeh)
            endif
         else if (size(cu,1)==3) then
            if (inorder==LINEAR) then
               call CUPERPMD2(cu(1,1,1),nx,ny,nxe,nyeh)
            else
               call CUPERPMD2(cu(1,2,2),nx,ny,nxe,nyeh)
            endif
         endif
         end subroutine icuperpmd2
!
         subroutine jbpoismd23(cu,bxy,ffh,ci,wm,nx,ny,order)
! calculates static vector potential for 2d vector field,
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: isign = 1, nxe, nyeh, nxv, nyhd, inorder
         real :: ax, ay, affp
         nxe = size(cu,2); nyeh = size(cu,3)/2
         nxv = size(ffh,1); nyhd = size(ffh,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call BPOISMD23(cu(1,1,1),bxy(1,1,1),isign,ffh,ax,ay,affp,ci,&
     &wm,nx,ny,nxe,nyeh,nxv,nyhd)
         else
            call BPOISMD23(cu(1,2,2),bxy(1,2,2),isign,ffh,ax,ay,affp,ci,&
     &wm,nx,ny,nxe,nyeh,nxv,nyhd)
         endif
         end subroutine jbpoismd23
!
         subroutine jibpoismd23(cu,bxy,ffh,ci,wm,nx,ny,order)
! calculates static magnetic field for periodic 2d vector field
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: bxy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: nxe, nyeh, nxv, nyhd, inorder
         nxe = size(cu,2); nyeh = size(cu,3)/2
         nxv = size(ffh,1); nyhd = size(ffh,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call IBPOISMD23(cu(1,1,1),bxy,ffh,ci,wm,nx,ny,nxe,nyeh,nxv,n&
     &yhd)
         else
            call IBPOISMD23(cu(1,2,2),bxy,ffh,ci,wm,nx,ny,nxe,nyeh,nxv,n&
     &yhd)
         endif
         end subroutine jibpoismd23
!
         subroutine imaxwelmd2(exy,bxy,cu,ffh,ci,dt,wf,wm,nx,ny,order)
! calculates maxwell's equation for 2d vector field,
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, dt, wf, wm
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: nxe, nyeh, nxv, nyhd, inorder
         nxe = size(cu,2); nyeh = size(cu,3)/2; nxv = size(ffh,1)
         nyhd = size(ffh,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call MAXWELMD2(exy,bxy,cu(1,1,1),ffh,ci,dt,wf,wm,nx,ny,nxe,n&
     &yeh,nxv,nyhd)
         else
            call MAXWELMD2(exy,bxy,cu(1,2,2),ffh,ci,dt,wf,wm,nx,ny,nxe,n&
     &yeh,nxv,nyhd)
         endif
         end subroutine imaxwelmd2    
!
         subroutine idmfieldmd2(q3,q,nx,ny)
! copies from quadruple to normal array in y dimension
! for 2d scalar data with mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q3, q
! local data
         integer :: nx2v, nyv, nxe, nyeh
         nx2v = size(q3,1)/2; nyv = size(q3,2)
         nxe = size(q,1); nyeh = size(q,2)/2
         call DMFIELDMD2(q3,q,nx,ny,nx2v,nyv,nxe,nyeh)
         end subroutine idmfieldmd2
!
         subroutine icmfieldmd2(cu3,cu,nx,ny)
! copies from quadruple to normal array in y dimension
! for 2d vector data with mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu3, cu
! local data
         integer :: nx2v, nyv, nxe, nyeh
         nx2v = size(cu3,2)/2; nyv = size(cu3,3)
         nxe = size(cu,2); nyeh = size(cu,3)/2
         call CMFIELDMD2(cu3,cu,nx,ny,nx2v,nyv,nxe,nyeh)
         end subroutine icmfieldmd2
!
         subroutine iemfieldmd2(fxy,exy,ffh,isign,nx,ny)
! combines and smooths 2d vector fields,
! mixed dirichlet-neumann/periodic boundaries
         implicit none
         integer :: isign, nx, ny
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffh
! local data
         integer :: nxv, nyv, nxe, nyeh, nyhd
         nxv = size(fxy,2)/4; nyv = size(fxy,3)
         nxe = 2*size(exy,2); nyeh = size(exy,3)/2
         nyhd = size(ffh,2)
         call EMFIELDMD2(fxy,exy,ffh,isign,nx,ny,nxv,nyv,nxe,nyeh,nyhd)
         end subroutine iemfieldmd2
!
         subroutine icpfieldmd2(fxy,exy,nx,ny)
! combines 2d electric fields, mixed dirichlet-neumann/periodic
! boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: fxy, exy
! local data
         integer :: nxv, nyv, nxe, nyeh
         nxv = size(fxy,2)/4; nyv = size(fxy,3)
         nxe = size(exy,2); nyeh = size(exy,3)/2
         call CPFIELDMD2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
         end subroutine icpfieldmd2
!
         subroutine iavpotmd23(bxy,axy,nx,ny,order)
! calculates 2d vector potential from magnetic field
! mixed dirichlet-neumann/periodic boundaries
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
            call AVPOTMD23(bxy,axy(1,1,1),nx,ny,nxe,nyeh)
         else
            call AVPOTMD23(bxy,axy(1,2,2),nx,ny,nxe,nyeh)
         endif
         end subroutine iavpotmd23
!
      end module hfield2d
