!-----------------------------------------------------------------------
!
      module dfield2d
!
! Fortran90 interface to 2d PIC Fortran77 library dfield2lib.f
! dfield2mod.f contains procedures to manage guard cells and solve fields
!              in fourier space with dirichlet boundary conditions:
!              defines module dfield2d
! laguard => ilaguard2 add guard cells for non-periodic scalar array
!            to replace quadratic with linear interpolation at the edges
!            calls LAGUARD2
! lcguard => ilcguard2 copy guard cells for scalar or 2 and 3 component
!            non-periodic vector arrays to replace quadratic with linear
!            interpolation at the edges.
!            calls LDGUARD2, LCGUARD2, LBGUARD2, or LDCGUARD2
! lcguard => ildguard2 copy guard cells for scalar arrays with various
!            interpolations.
!            calls LDGUARD2
! dblsin => idblsin2c double array for 2d vector data for dirichlet
!           boundary conditions.
!           calls DBLSIN2B or DBLSIN2C
! dblsin => idblsin2d double array for 2d scalar data for dirichlet
!           boundary conditions.
!           calls DBLSIN2D
! hafdbl => ihafdbl2c extract data from doubled array for 2d vector data
!           calls HAFDBL2D, HAFDBL2C, or HAFDBL2B
! hafdbl => ihafdbl2d extract data from doubled array for 2d scalar data
!           calls HAFDBL2D
! poisd_init => ipoisd22init initializes tables for field solvers, with
!               dirichlet boundary conditions.
!               calls POISDX22
! poisdx => ipoisdx2 solves 2d poisson equation for potential, with
!           dirichlet boundary conditions, using ffts.
!           calls POISDX2
! poisdx => ispoisdx2 smoother for 2d periodic scalar field, with
!           dirichlet boundary conditions, using ffts.
!           calls POISDX2
! poisdx => ipoisdx23 solves 2-1/2d poisson equation for electric force.
!           with dirichlet boundary conditions, using ffts.
!           calls POISDX22, or POISDX23
! poisd => ipoisd2 solves 2d poisson equation for potential, with
!           dirichlet boundary conditions, using sine/cosine transforms.
!           calls POISD2
! poisd => ispoisd2 smoother for 2d periodic scalar field, with
!           dirichlet boundary conditions, using sine/cosine transforms.
!           calls POISD2
! poisd => ipoisd23 solves 2-1/2d poisson equation for electric force.
!           with dirichlet boundary conditions, using sine/cosine
!           transforms.
!           calls POISD22, or POISD23
! cuperpd => icuperpd2 calculates the transverse part of periodic 2-1/2d
!           vector field with dirichlet boundary conditions, using
!           sine/cosine transforms.
!           calls CUPERPD22, or CUPERPD2
! bpoisd => jbpoisd23 solves 2-1/2d vector poisson equation for magnetic
!           force with dirichlet boundary conditions, using sine/cosine
!           transforms.
!           calls BPOISD22, or BPOISD23
! bpoisd => sbpoisd23 smoother for 2-1/2d non-periodic vector field with
!           dirichlet boundary conditions, using sine/cosine transforms.
!           calls BPOISD22, or BPOISD23
! apoisd => iapoisd23 solves 2-1/2d vector poisson equation for vector
!           potential with dirichlet boundary conditions, using
!           sine/cosine transforms.
!           calls BPOISD22, or BPOISD23
! ibpoisd => jibpoisd23 solves vector poisson equation for magnetic
!            field with dirichlet boundary conditions, using sine/cosine
!            transforms.
!            calls IBPOISD23
! maxweld => imaxweld2 solves maxwell equation for electric and magnetic
!            fields with dirichlet boundary conditions, using
!            sine/cosine transforms.
!            calls MAXWELD2
! cmfieldd => icmfieldd2 copies vector data from doubled fft to
!             sine/cosine format.
!             calls CMFIELDD2
! cmfieldd => idmfieldd2 copies scalar data from doubled fft to
!             sine/cosine format.
!             calls DMFIELDD2
! emfieldd => iemfieldd2 calculates periodic electric and magnetic
!             forces from fields given by maxwell and poisson equations
!             with dirichlet boundary conditions.
!             calls EMFIELDD2
! cpfieldd => icpfieldd2 copies electric in sine/cosine to doubled fft
!             format.
!             calls CPFIELDD2
! avpotd => iavpotd23 calculates vector potential from magnetic field
!           with dirichlet boundary conditions, using sine/cosine
!           transforms.
!           calls AVPOTD23
! poyntdx => idpoyntdx2 calculates the momentum in the darwin field,
!            with dirichlet boundary conditions, using ffts.
!            calls DPOYNTDX2
! dblsinm => idblsin2m double array for 2d tensor data for dirichlet
!            boundary conditions.
!            calls DBLSIN2M
! dcuperpdx => idcuperpdx23 calculate transverse derivative of 2-1/2d
!              current density from momentum flux, with dirichlet
!              boundary conditions, using ffts.
!              calls DCUPERPDX23
! adcuperpdx => iadcuperpdx23 calculate transverse derivative of 2-1/2d
!               current density from momentum flux and acceleration
!               density, with dirichlet boundary conditions, using ffts.
!               calls ADCUPERPDX23
! epoisd_init => jepoisdx23 initializes tables for darwin field solver,
!                with dirichlet boundary conditions.
!                calls EPOISDX23
! epoisdx => jepoisdx23 solves 2-1/2d vector poisson equation for
!            transverse electric force with dirichlet boundary
!            conditions, using ffts.
!            calls EPOISDX23
! iepoisdx => jiepoisdx23 solves 2-1/2d vector poisson equation
!             for transverse electric field without smoothing, with
!             dirichlet boundary conditions, using ffts.
!             calls EPOISDX23
! dapoisdx => idapoisdx23 calculates derivatives of smoothed 2-1/2d
!             vector potential with dirichlet boundary conditions, using
!             ffts.
!             calls APOISDX23
! sapoisdx => isapoisdx23 calculates smoothed 2-1/2d vector potential
!             with dirichlet boundary conditions, using ffts.
!             calls APOISDX23
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: july 21, 2010
!
      use globals, only: LINEAR, QUADRATIC
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: LCGUARD2, LDGUARD2, LBGUARD2, LDCGUARD2
      public :: LSGUARD2, LAGUARD2, LSGUARD2L
      public :: LSCGUARD2, LSCGUARD22, LSCGUARD2L, LSCGUARD22L
      public :: LACGUARD2, LACGUARD22
      public :: LACFGUARD2, LACFGUARD22, LSCFGUARD2L, LSCFGUARD22L
      public :: LSMCGUARD2, LSMCGUARD22, LSMCGUARD2L, LSMCGUARD22L
      public :: LAMCGUARD2
      public :: laguard, lcguard
      public :: dblsin, hafdbl, poisd_init, poisdx, poisd, cuperpd
      public :: bpoisd, apoisd, ibpoisd, maxweld, cmfieldd, emfieldd
      public :: cpfieldd, avpotd, dblsinm, poyntdx, dcuperpdx
      public :: adcuperpdx, epoisd_init, epoisdx, iepoisdx, dapoisdx
      public :: sapoisdx
! debug
      public :: icuperpdx2
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine LCGUARD2(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         end subroutine
      end interface
      interface
         subroutine LDGUARD2(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine LBGUARD2(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: bxy
         end subroutine
      end interface
      interface
         subroutine LSCGUARD2(cu,xj0,yj0,zj0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine LSCGUARD22(cu,xj0,yj0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: xj0, yj0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine LSGUARD2(q,qi0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine LACGUARD2(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: cu
         end subroutine
      end interface
      interface
         subroutine LACGUARD22(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: cu
         end subroutine
      end interface
      interface
         subroutine LAGUARD2(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: q
         end subroutine
      end interface
      interface
         subroutine LSCGUARD2L(cu,xj0,yj0,zj0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: xj0,yj0,zj0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine LSCGUARD22L(cu,xj0,yj0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: xj0,yj0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine LSGUARD2L(q,qi0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine DBLSIN2C(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(2,nx2v,ny2) :: cu2
         end subroutine
      end interface
      interface
         subroutine DBLSIN2D(q,q2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,ny2) :: q2
         end subroutine
      end interface
      interface
         subroutine DBLSIN2B(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: cu
         real :: cu
         real, dimension(3,nx2v,ny2) :: cu2
         end subroutine
      end interface
      interface
         subroutine HAFDBL2C(fxy,fxy2,nx,ny,nxe,nye,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v, ny2
!        real, dimension(*) :: fxy
         real :: fxy
         real, dimension(2,nx2v,ny2) :: fxy2
         end subroutine
      end interface
      interface
         subroutine HAFDBL2D(q,q2,nx,ny,nxe,nye,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v, ny2
!        real, dimension(*) :: q
         real :: q
         real, dimension(nx2v,ny2) :: q2
         end subroutine
      end interface
      interface
         subroutine HAFDBL2B(bxy,bxy2,nx,ny,nxe,nye,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v, ny2
!        real, dimension(*) :: bxy
         real :: bxy
         real, dimension(3,nx2v,ny2) :: bxy2
         end subroutine
      end interface
      interface
         subroutine POISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,n&
     &y2d)
         implicit none
         integer :: isign, nx, ny, nx2v, ny2d
         real :: ax, ay, affp, we
         real, dimension(nx2v,ny2d) :: q, fx, fy
         complex, dimension(nx2v/2,ny) :: ffd
         end subroutine
      end interface
       interface
         subroutine POISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye&
     &,nx2v)
         implicit none
         integer :: isign, nx, ny, nxe, nye, nx2v
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fx, fy
         real :: q, fx, fy
         complex, dimension(nx2v/2,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine POISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2&
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
         subroutine POISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye,&
     &nxv)
         implicit none
         integer :: isign, nx, ny, nxe, nye, nxv
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine POISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2&
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
         subroutine POISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye,&
     &nxv)
         implicit none
         integer :: isign, nx, ny, nxe, nye, nxv
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine CUPERPDX2(cu,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real, dimension(3,2*nxv,ny2d) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPD2(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPDX22(cu,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real, dimension(2,2*nxv,ny2d) :: cu
         end subroutine
      end interface
      interface
         subroutine CUPERPD22(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine BPOISDX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nx&
     &v,ny2d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, ci, wm
         real, dimension(3,2*nxv,ny2d) :: cu, bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine BPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxe&
     &,nye,nxv)
         implicit none
         integer :: isign, nx, ny, nxe, nye, nxv
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy
         real :: cu, bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine BPOISDX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
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
         subroutine BPOISD22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,&
     &nxe,nye,nxv)
         implicit none
         integer :: isign, nx, ny, nxe, nye, nxv
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy, bz
         real :: cu, bxy, bz
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine IBPOISDX23(cu,bxy,ffd,ci,wm,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real :: ci, wm
         real, dimension(3,2*nxv,ny2d) :: cu, bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine IBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
         implicit none
         integer :: nx, ny, nxe, nye, nxv
         real :: ci, wm
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(3,nxe/2,nye) :: bxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine MAXWELDX2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real :: ci, dt, wf, wm
         real, dimension(3,2*nxv,ny2d) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine MAXWELD2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxe,nye,nx&
     &v)
         implicit none
         integer :: nx, ny, nxe, nye, nxv
         real :: ci, dt, wf, wm
         complex, dimension(3,nxe/2,nye) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine DMFIELDD2(q2,q,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: nx, ny, nxv, ny2d, nxe, nye
         real, dimension(2*nxv,ny2d) :: q2
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine CMFIELDD2(cu2,cu,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: nx, ny, nxv, ny2d, nxe, nye
         real, dimension(3,2*nxv,ny2d) :: cu2
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine EMFIELDD2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d, nxe, nye
         real, dimension(3,2*nxv,ny2d) :: fxy
         complex, dimension(3,nxe/2,nye) :: exy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine PMFIELDD2(pot2,pot,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: nx, ny, nxv, ny2d, nxe, nye
         real, dimension(2*nxv,ny2d) :: pot2
         real, dimension(nxe,nye) :: pot
         end subroutine
      end interface
      interface
         subroutine CPFIELDD2(fxy,exy,nx,ny,nxv,ny2d,nxe,nye)
         implicit none
         integer :: nx, ny, nxv, ny2d, nxe, nye
         real, dimension(3,2*nxv,ny2d) :: fxy
         real, dimension(3,nxe,nye) :: exy
         end subroutine
      end interface
      interface
         subroutine AVPOTDX23(bxy,axy,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         complex, dimension(3,nxv,ny2d) :: bxy, axy
         end subroutine
      end interface
      interface
         subroutine AVPOTD23(bxy,axy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         complex, dimension(3,nxe/2,nye) :: bxy
!        real, dimension(*) :: axy
         real :: axy
         end subroutine
      end interface
      interface
         subroutine DPOYNTDX2(q,cu,ffd,ci,sx,sy,sz,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real :: ci, sx, sy, sz
         real, dimension(2*nxv,ny2d) :: q
         real, dimension(3,2*nxv,ny2d) :: cu
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine LACFGUARD2(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine LACFGUARD22(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine LSCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine LSCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine LSMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ngx,ngy,n&
     &xe,nye)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(4,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine LSMCGUARD22(amu,x2y2m0,xym0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: x2y2m0, xym0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(2,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine LAMCGUARD2(amu,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
!        real, dimension(*) :: amu
         real :: amu
         end subroutine
      end interface
      interface
         subroutine LSMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ngx,ngy,&
     &nxe,nye)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(4,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine LSMCGUARD22L(amu,x2y2m0,xym0,nx,ny,ngx,ngy,nxe,nye)
         implicit none
         real :: x2y2m0, xym0
         integer :: nx, ny, ngx, ngy, nxe, nye
         real, dimension(2,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine DBLSIN2M(amu,amu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: amu
         real :: amu
         real, dimension(4,nx2v,ny2) :: amu2
         end subroutine
      end interface
      interface
         subroutine DBLSIN22M(amu,amu2,nx,ny,nxv,nyv,nx2v,ny2)
         implicit none
         integer :: nx, ny, nxv, nyv, nx2v, ny2
!        real, dimension(*) :: amu
         real :: amu
         real, dimension(2,nx2v,ny2) :: amu2
         end subroutine
      end interface
      interface
         subroutine DCUPERPDX23(dcu,amu,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real, dimension(3,2*nxv,ny2d) :: dcu
         real, dimension(4,2*nxv,ny2d) :: amu
         end subroutine
      end interface
      interface
         subroutine ADCUPERPDX23(dcu,amu,nx,ny,nxv,ny2d)
         implicit none
         integer :: nx, ny, nxv, ny2d
         real, dimension(3,2*nxv,ny2d) :: dcu
         real, dimension(4,2*nxv,ny2d) :: amu
         end subroutine
      end interface
      interface
         subroutine EPOISDX23(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,&
     &ny,nxv,ny2d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, wp0, ci, wf
         real, dimension(3,2*nxv,ny2d) :: dcu, exy
         complex, dimension(nxv,ny) :: fff
         end subroutine
      end interface
      interface
         subroutine HAFDBL2N(daxy,daxy2,nx,ny,nxe,nye,nx2v,ny2,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, nx2v, ny2, ndim
!        real, dimension(*) :: daxy
         real :: daxy
         real, dimension(ndim,nx2v,ny2) :: daxy2
         end subroutine
      end interface
      interface
         subroutine LDCGUARD2(daxy,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: daxy
         end subroutine
      end interface
      interface
         subroutine APOISDX23(cu,daxy,axy,isign,ffd,ax,ay,affp,ci,wm,nx,&
     &ny,nxv,ny2d)
         implicit none
         integer :: isign, nx, ny, nxv, ny2d
         real :: ax, ay, affp, ci, wm
         real, dimension(3,2*nxv,ny2d) :: cu, axy
         real, dimension(5,2*nxv,ny2d) :: daxy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!      
      interface laguard
         module procedure ilaguard2
      end interface
!
      interface lcguard
         module procedure ilcguard2
         module procedure ildguard2
      end interface
!
      interface dblsin
         module procedure idblsin2c
         module procedure idblsin2d
      end interface
!
      interface hafdbl
         module procedure ihafdbl2c
         module procedure ihafdbl2d
      end interface
!
      interface poisd_init
         module procedure ipoisd22init
      end interface
!
      interface poisdx
         module procedure ipoisdx2
         module procedure ispoisdx2
         module procedure ipoisdx23
      end interface
!
      interface poisd
         module procedure ipoisd2
         module procedure ispoisd2
         module procedure ipoisd23
      end interface
!
      interface cuperpd
         module procedure icuperpd2
      end interface
!
      interface bpoisd
         module procedure jbpoisd23
         module procedure sbpoisd23
      end interface
!
      interface apoisd
         module procedure iapoisd23
      end interface
!
      interface ibpoisd
         module procedure jibpoisd23
      end interface
!
      interface maxweld
         module procedure imaxweld2
      end interface
!
      interface cmfieldd
         module procedure icmfieldd2
         module procedure idmfieldd2
      end interface
!
      interface emfieldd
         module procedure iemfieldd2
      end interface
!
      interface cpfieldd
         module procedure icpfieldd2
      end interface
!
      interface avpotd
         module procedure iavpotd23
      end interface
!
      interface poyntdx
         module procedure idpoyntdx2
      end interface
!
      interface dblsinm
         module procedure idblsin2m
      end interface
!
      interface dcuperpdx
         module procedure idcuperpdx23
      end interface
!
      interface adcuperpdx
         module procedure iadcuperpdx23
      end interface
!
       interface epoisd_init
         module procedure iepoisdx23init
      end interface
!
      interface epoisdx
         module procedure jepoisdx23
      end interface
!
      interface iepoisdx
         module procedure jiepoisdx23
      end interface
!
      interface dapoisdx
         module procedure idapoisdx23
      end interface
!
      interface sapoisdx
         module procedure isapoisdx23
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ilaguard2(q,nx,ny,inorder)
! add guard cells for 2d scalar data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call LAGUARD2(q(2,2),nx-2,ny-2,nxe,nye)
         endif
         end subroutine ilaguard2
!
         subroutine ilcguard2(fxy,nx,ny,inorder)
! copy guard cells for 2d vector data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: ndim, nxe, nye, order
         nxe = size(fxy,2); nye = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order==QUADRATIC) then
               call LDGUARD2(fxy,nx,ny,nxe,nye)
            endif
         case (2)
            if (order==QUADRATIC) then
               call LCGUARD2(fxy,nx,ny,nxe,nye)
            endif
         case (3)
            if (order==QUADRATIC) then
               call LBGUARD2(fxy,nx,ny,nxe,nye)
            endif
         case (4:)
            ndim = size(fxy,1)
            if (order==QUADRATIC) then
               call LDCGUARD2(fxy,nx,ny,nxe,nye,ndim)
            endif
         end select
         end subroutine ilcguard2
!
         subroutine ildguard2(q,nx,ny,inorder)
! copy guard cells for 2d scalar data
! disable quadratic interpolation at edge
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order /= LINEAR) then
            call LDGUARD2(q,nx,ny,nxe,nye)
         endif
         end subroutine ildguard2
!
         subroutine idblsin2c(cu,cu2,nx,ny,inorder)
! double array in each dimension for 2d vector data
! for dirichlet boundary conditions
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
         select case(size(cu,1))
         case (3)
            if (order==LINEAR) then
               call DBLSIN2B(cu(1,1,1),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call DBLSIN2B(cu(1,2,2),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         case (2)
            if (order==LINEAR) then
               call DBLSIN2C(cu(1,1,1),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call DBLSIN2C(cu(1,2,2),cu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         end select
         end subroutine idblsin2c
!
         subroutine idblsin2d(q,q2,nx,ny,inorder)
! double array in each dimension for 2d scalar data
! for dirichlet boundary conditions
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
            call DBLSIN2D(q(1,1),q2,nx,ny,nxv,nyv,nx2v,ny2)
         else
            call DBLSIN2D(q(2,2),q2,nx,ny,nxv,nyv,nx2v,ny2)
         endif
         end subroutine idblsin2d
!
         subroutine ihafdbl2c(fxy,fxy2,nx,ny,inorder)
! copy from double to normal array in each dimension for 2d vector data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         real, dimension(:,:,:), pointer :: fxy2
! local data
         integer :: ndim, nxe, nye, nx2v, ny2, order
         ndim = size(fxy,1)
         nxe = size(fxy,2);  nye = size(fxy,3)
         nx2v = size(fxy2,2);  ny2 = size(fxy2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(fxy,1))
         case (1)
            if (order==LINEAR) then
               call HAFDBL2D(fxy(1,1,1),fxy2,nx,ny,nxe,nye,nx2v,ny2)
            else
               call HAFDBL2D(fxy(1,2,2),fxy2,nx,ny,nxe,nye,nx2v,ny2)
            endif
         case (2)
            if (order==LINEAR) then
               call HAFDBL2C(fxy(1,1,1),fxy2,nx,ny,nxe,nye,nx2v,ny2)
            else
               call HAFDBL2C(fxy(1,2,2),fxy2,nx,ny,nxe,nye,nx2v,ny2)
            endif
         case (3)
            if (order==LINEAR) then
               call HAFDBL2B(fxy(1,1,1),fxy2,nx,ny,nxe,nye,nx2v,ny2)
            else
               call HAFDBL2B(fxy(1,2,2),fxy2,nx,ny,nxe,nye,nx2v,ny2)
            endif
         case (4:)
            if (order==LINEAR) then
               call HAFDBL2N(fxy(1,1,1),fxy2,nx,ny,nxe,nye,nx2v,ny2,ndim&
     &)
            else
               call HAFDBL2N(fxy(1,2,2),fxy2,nx,ny,nxe,nye,nx2v,ny2,ndim&
     &)
            endif
         end select
         end subroutine ihafdbl2c
!
         subroutine ihafdbl2d(q,q2,nx,ny,inorder)
! copy from double to normal array in each dimension for 2d scalar data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q
         real, dimension(:,:), pointer :: q2
! local data
         integer :: nxe, nye, nx2v, ny2, order
         nxe = size(q,1);  nye = size(q,2)
         nx2v = size(q2,1);  ny2 = size(q2,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call HAFDBL2D(q(1,1),q2,nx,ny,nxe,nye,nx2v,ny2)
         else
            call HAFDBL2D(q(2,2),q2,nx,ny,nxe,nye,nx2v,ny2)
         endif
         end subroutine ihafdbl2d
!
         subroutine ipoisdx2(q,fx,ffd,we,nx,ny)
! poisson solver for 2d potential, conducting boundaries
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
         call POISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,ny2d)
         end subroutine ipoisdx2
!
         subroutine ispoisdx2(q,fy,ffd,nx,ny)
! smoother for 2d scalar field, conducting boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 2, nx2v, ny2d
         real :: ax, ay, affp, we
         real, dimension(1,1) :: fx
         nx2v = size(q,1); ny2d = size(q,2)
         call POISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,ny2d)
         end subroutine ispoisdx2
!
         subroutine ipoisd22init(ffd,ax,ay,affp,nx,ny)
! initialize 2d electric field solver, conducting boundaries
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
         call POISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
         end subroutine ipoisd22init
!
         subroutine ipoisdx23(q,fxy,ffd,we,nx,ny)
! poisson solver for 2d electric field, conducting boundaries
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
         select case(size(fxy,1))
         case (2)
            call POISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
         case (3)
            call POISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
         end select
         end subroutine ipoisdx23
!
         subroutine ipoisd2(q,fx,ffd,we,nx,ny,order)
! poisson solver for 2d potential, conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 1, nxe, nye, nxv, nx2v, inorder
         real :: ax, ay, affp
         real, dimension(2,2) :: fy
         nxe = size(q,1); nye = size(q,2); nxv = size(ffd,1)
         nx2v = 2*nxv
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISD2(q(1,1),fx(1,1),fy(1,1),isign,ffd,ax,ay,affp,we,n&
     &x,ny,nxe,nye,nx2v)
         else
            call POISD2(q(2,2),fx(2,2),fy(2,2),isign,ffd,ax,ay,affp,we,n&
     &x,ny,nxe,nye,nx2v)
         endif
         end subroutine ipoisd2
!
         subroutine ispoisd2(q,fy,ffd,nx,ny,order)
! smoother for 2d scalar field, conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 2, nxe, nye, nxv, nx2v, inorder
         real :: ax, ay, affp, we
         real, dimension(2,2) :: fx
         nxe = size(q,1); nye = size(q,2); nxv = size(ffd,1)
         nx2v = 2*nxv
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call POISD2(q(1,1),fx(1,1),fy(1,1),isign,ffd,ax,ay,affp,we,n&
     &x,ny,nxe,nye,nx2v)
         else
            call POISD2(q(2,2),fx(2,2),fy(2,2),isign,ffd,ax,ay,affp,we,n&
     &x,ny,nxe,nye,nx2v)
         endif
         end subroutine ispoisd2
!
         subroutine ipoisd23(q,fxy,ffd,we,nx,ny,order)
! poisson solver for 2-1/2d electric field, conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: we
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = -1, nxe, nye, nxv, inorder
         real :: ax, ay, affp
         nxe = size(q,1); nye = size(q,2); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(fxy,1))
         case (2)
            if (inorder==LINEAR) then
               call POISD22(q(1,1),fxy(1,1,1),isign,ffd,ax,ay,affp,we,nx&
     &,ny,nxe,nye,nxv)
            else
               call POISD22(q(2,2),fxy(1,2,2),isign,ffd,ax,ay,affp,we,nx&
     &,ny,nxe,nye,nxv)
            endif
         case (3)
            if (inorder==LINEAR) then
               call POISD23(q(1,1),fxy(1,1,1),isign,ffd,ax,ay,affp,we,nx&
     &,ny,nxe,nye,nxv)
            else
               call POISD23(q(2,2),fxy(1,2,2),isign,ffd,ax,ay,affp,we,nx&
     &,ny,nxe,nye,nxv)
            endif
         end select
         end subroutine ipoisd23
!
         subroutine icuperpdx2(cu,nx,ny)
! calculates transverse part of 2d vector field, conducting boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxv, ny2d
         nxv = size(cu,2)/2; ny2d = size(cu,3)
         select case(size(cu,1))
         case (2)
            call CUPERPDX22(cu,nx,ny,nxv,ny2d)
         case (3)
            call CUPERPDX2(cu,nx,ny,nxv,ny2d)
         end select
         end subroutine icuperpdx2
!
         subroutine icuperpd2(cu,nx,ny,order)
! calculates transverse part of 2d vector field, conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nye, inorder
         nxe = size(cu,2); nye = size(cu,3)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call CUPERPD22(cu(1,1,1),nx,ny,nxe,nye)
            else
               call CUPERPD22(cu(1,2,2),nx,ny,nxe,nye)
            endif
         case (3)
            if (inorder==LINEAR) then
               call CUPERPD2(cu(1,1,1),nx,ny,nxe,nye)
            else
               call CUPERPD2(cu(1,2,2),nx,ny,nxe,nye)
            endif
         end select
         end subroutine icuperpd2
!
         subroutine jbpoisd23(cu,bxy,ffd,ci,wm,nx,ny,order)
! calculates static magnetic field for 2d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = -1, nxe, nye, nxv, inorder
         real :: ax, ay, affp, bxy0
         nxe = size(cu,2); nye = size(cu,3); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call BPOISD22(cu(1,1,1),bxy0,bxy(1,1,1),isign,ffd,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxe,nye,nxv)
            else
               call BPOISD22(cu(1,2,2),bxy0,bxy(1,2,2),isign,ffd,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxe,nye,nxv)
            endif
         case (3)
            if (inorder==LINEAR) then
               call BPOISD23(cu(1,1,1),bxy(1,1,1),isign,ffd,ax,ay,affp,c&
     &i,wm,nx,ny,nxe,nye,nxv)
            else
               call BPOISD23(cu(1,2,2),bxy(1,2,2),isign,ffd,ax,ay,affp,c&
     &i,wm,nx,ny,nxe,nye,nxv)
            endif
         end select
         end subroutine jbpoisd23
!
         subroutine sbpoisd23(cu,bxy,ffd,nx,ny,order)
! smoother for 2d vector field, conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 2, nxe, nye, nxv, inorder
         real :: ax, ay, affp, bz0, ci, wm
         nxe = size(cu,2); nye = size(cu,3); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call BPOISD22(cu(1,1,1),bxy(1,1,1),bz0,isign,ffd,ax,ay,af&
     &fp,ci,wm,nx,ny,nxe,nye,nxv)
            else
               call BPOISD22(cu(1,2,2),bxy(1,2,2),bz0,isign,ffd,ax,ay,af&
     &fp,ci,wm,nx,ny,nxe,nye,nxv)
            endif
         case (3)
            if (inorder==LINEAR) then
               call BPOISD23(cu(1,1,1),bxy(1,1,1),isign,ffd,ax,ay,affp,c&
     &i,wm,nx,ny,nxe,nye,nxv)
            else
               call BPOISD23(cu(1,2,2),bxy(1,2,2),isign,ffd,ax,ay,affp,c&
     &i,wm,nx,ny,nxe,nye,nxv)
            endif
         end select
         end subroutine sbpoisd23
!
         subroutine iapoisd23(cu,axy,ffd,ci,wm,nx,ny,order)
! calculates static vector potential for 2d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, axy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 1, nxe, nye, nxv, inorder
         real :: ax, ay, affp, bz0
         nxe = size(cu,2); nye = size(cu,3); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         select case(size(cu,1))
         case (2)
            if (inorder==LINEAR) then
               call BPOISD22(cu(1,1,1),axy(1,1,1),bz0,isign,ffd,ax,ay,af&
     &fp,ci,wm,nx,ny,nxe,nye,nxv)
            else
               call BPOISD22(cu(1,2,2),axy(1,2,2),bz0,isign,ffd,ax,ay,af&
     &fp,ci,wm,nx,ny,nxe,nye,nxv)
            endif
         case (3)
            if (inorder==LINEAR) then
               call BPOISD23(cu(1,1,1),axy(1,1,1),isign,ffd,ax,ay,affp,c&
     &i,wm,nx,ny,nxe,nye,nxv)
            else
               call BPOISD23(cu(1,2,2),axy(1,2,2),isign,ffd,ax,ay,affp,c&
     &i,wm,nx,ny,nxe,nye,nxv)
            endif
         end select
         end subroutine iapoisd23
!
         subroutine jibpoisd23(cu,bxy,ffd,ci,wm,nx,ny,order)
! calculates static magnetic field for periodic 2d vector field
! conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: bxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxe, nye, nxv, inorder
         nxe = size(cu,2); nye = size(cu,3); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call IBPOISD23(cu(1,1,1),bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
         else
            call IBPOISD23(cu(1,2,2),bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
         endif
         end subroutine jibpoisd23
!
         subroutine imaxweld2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,order)
! calculates maxwell's equation for 2d vector field,
! conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         real :: ci, dt, wf, wm
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxe, nye, nxv, inorder
         nxe = size(cu,2); nye = size(cu,3); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call MAXWELD2(exy,bxy,cu(1,1,1),ffd,ci,dt,wf,wm,nx,ny,nxe,ny&
     &e,nxv)
         else
            call MAXWELD2(exy,bxy,cu(1,2,2),ffd,ci,dt,wf,wm,nx,ny,nxe,ny&
     &e,nxv)
         endif
         end subroutine imaxweld2    
!
         subroutine icmfieldd2(cu2,cu,nx,ny)
! copies from double to normal array in y dimension for 2d vector data
! conducting boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: cu2, cu
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(cu2,2)/2; ny2d = size(cu2,3)
         nxe = size(cu,2); nye = size(cu,3)
         call CMFIELDD2(cu2,cu,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine icmfieldd2
!
         subroutine idmfieldd2(q2,q,nx,ny)
! copies from double to normal array in y dimension for 2d scalar data
         implicit none
         integer :: nx, ny
         real, dimension(:,:), pointer :: q2, q
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(q2,1)/2; ny2d = size(q2,2)
         nxe = size(q,1); nye = size(q,2)
         call DMFIELDD2(q2,q,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine idmfieldd2
!
         subroutine iemfieldd2(fxy,exy,ffd,isign,nx,ny)
! combines and smooths 2d vector fields, conducting boundaries
         implicit none
         integer :: isign, nx, ny
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(fxy,2)/2; ny2d = size(fxy,3)
         nxe = 2*size(exy,2); nye = size(exy,3)
         call EMFIELDD2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine iemfieldd2
!
         subroutine icpfieldd2(fxy,exy,nx,ny)
! combines 2d electric fields, conducting boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: fxy, exy
! local data
         integer :: nxv, ny2d, nxe, nye
         nxv = size(fxy,2)/2; ny2d = size(fxy,3)
         nxe = size(exy,2); nye = size(exy,3)
         call CPFIELDD2(fxy,exy,nx,ny,nxv,ny2d,nxe,nye)
         end subroutine icpfieldd2
!
         subroutine iavpotd23(bxy,axy,nx,ny,order)
! calculates 2d vector potential from magnetic field
! conducting boundaries
         implicit none
         integer :: nx, ny
         integer, optional :: order
         complex, dimension(:,:,:), pointer :: bxy
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: nxe, nye, inorder
         nxe = size(axy,2); nye = size(axy,3)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call AVPOTD23(bxy,axy(1,1,1),nx,ny,nxe,nye)
         else
            call AVPOTD23(bxy,axy(1,2,2),nx,ny,nxe,nye)
         endif
         end subroutine iavpotd23
!
         subroutine idpoyntdx2(q,cu,ffd,ci,sx,sy,sz,nx,ny)
! calculates the momentum in the darwin field, conducting boundaries
         implicit none
         integer :: nx, ny
         real :: ci, sx, sy, sz
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxv, ny2d
         nxv = size(q,1)/2; ny2d = size(q,2)
         select case(size(cu,1))
!        case (2)
!           call DPOYNTDX22(q,cu,ffd,ci,sx,sy,sz,nx,ny,nxv,ny2d)
!           sz = 0.0
         case (3)
            call DPOYNTDX2(q,cu,ffd,ci,sx,sy,sz,nx,ny,nxv,ny2d)
         end select
         end subroutine idpoyntdx2
!
         subroutine idblsin2m(amu,amu2,nx,ny,inorder)
! double array in each dimension for 2d tensor data
! for dirichlet boundary conditions
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: amu
         real, dimension(:,:,:), pointer :: amu2
! local data
         integer :: nxv, nyv, nx2v, ny2, order
         nxv = size(amu,2);  nyv = size(amu,3)
         nx2v = size(amu2,2);  ny2 = size(amu2,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (size(amu,1)==2) then
            if (order==LINEAR) then
               call DBLSIN22M(amu(1,1,1),amu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call DBLSIN22M(amu(1,2,2),amu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         else if (size(amu,1)==4) then
            if (order==LINEAR) then
               call DBLSIN2M(amu(1,1,1),amu2,nx,ny,nxv,nyv,nx2v,ny2)
            else
               call DBLSIN2M(amu(1,2,2),amu2,nx,ny,nxv,nyv,nx2v,ny2)
            endif
         endif
         end subroutine idblsin2m
!
         subroutine idcuperpdx23(dcu,amu,nx,ny)
! calculates the transverse part of periodic 2d vector field
! from momentum flux tensor, conducting boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: dcu, amu
! local data
         integer :: nxv, ny2d
         nxv = size(dcu,2)/2; ny2d = size(dcu,3)
         select case(size(dcu,1))
!        case (2)
!           call DCUPERPDX22(dcu,amu,nx,ny,nxv,ny2d)
         case (3)
            call DCUPERPDX23(dcu,amu,nx,ny,nxv,ny2d)
         end select
         end subroutine idcuperpdx23
!
         subroutine iadcuperpdx23(dcu,amu,nx,ny)
! calculates the transverse part of periodic 2d vector field
! from acceleration vector and momentum flux tensor
! conducting boundaries
         implicit none
         integer :: nx, ny
         real, dimension(:,:,:), pointer :: dcu, amu
! local data
         integer :: nxv, ny2d
         nxv = size(dcu,2)/2; ny2d = size(dcu,3)
         select case(size(dcu,1))
!        case (2)
!           call ADCUPERPDX22(dcu,amu,nx,ny,nxv,ny2d)
         case (3)
            call ADCUPERPDX23(dcu,amu,nx,ny,nxv,ny2d)
         end select
         end subroutine iadcuperpdx23
!
         subroutine iepoisdx23init(fff,ax,ay,affp,wp0,ci,nx,ny)
! initialize 2d periodic transverse electric field solver
         implicit none
         integer :: nx, ny
         real :: ax, ay, affp, wp0, ci
         complex, dimension(:,:), pointer :: fff
! local data
         integer :: isign = 0, nxv, ny2d = 1
         real :: wf
         real, dimension(1,1,1) :: dcu, exy
         nxv = size(fff,1)
         call EPOISDX23(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,ny,nxv&
     &,ny2d)
         end subroutine iepoisdx23init
!
         subroutine jepoisdx23(dcu,exy,fff,ci,wf,nx,ny)
! calculates transverse electric field for periodic 2d vector field
! conducting boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wf
         real, dimension(:,:,:), pointer :: dcu, exy
         complex, dimension(:,:), pointer :: fff
! local data
         integer :: isign = -1, nxv, ny2d
         real :: ax, ay, affp, wp0
         nxv = size(dcu,2)/2; ny2d = size(dcu,3)
         select case(size(dcu,1))
!        case (2)
!           call EPOISDX22(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,ny,&
!    &nxv,ny2d)
         case (3)
            call EPOISDX23(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,ny,&
     &nxv,ny2d)
         end select
         end subroutine jepoisdx23
!
         subroutine jiepoisdx23(dcu,exy,fff,ci,wf,nx,ny)
! calculates transverse electric field for periodic 2d vector field
! without smoothing, conducting boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wf
         real, dimension(:,:,:), pointer :: dcu, exy
         complex, dimension(:,:), pointer :: fff
! local data
         integer :: isign = 1, nxv, ny2d
         real :: ax, ay, affp, wp0
         nxv = size(dcu,2)/2; ny2d = size(dcu,3)
         select case(size(dcu,1))
!        case (2)
!           call EPOISDX22(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,ny,&
!    &nxv,ny2d)
         case (3)
            call EPOISDX23(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,ny,&
     &nxv,ny2d)
         end select
         end subroutine jiepoisdx23
!
         subroutine idapoisdx23(cu,daxy,ffd,ci,wm,nx,ny)
! calculates derivatives of smoothed vector potential
! for conducting boundaries
         implicit none
         integer :: nx, ny
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, daxy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = -1, nxv, ny2d
         real :: ax, ay, affp
         real, dimension(1,1,1) :: axy0
         nxv = size(cu,2)/2; ny2d = size(cu,3)
         select case(size(cu,1))
!        case (2)
!           call APOISDX22(cu,daxy,axy0,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
!    &,nxv,ny2d)
         case (3)
            call APOISDX23(cu,daxy,axy0,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
     &,nxv,ny2d)
         end select
         end subroutine idapoisdx23
!
         subroutine isapoisdx23(cu,axy,ffd,ci,nx,ny)
! calculates smoothed vector potential, for conducting boundaries
         implicit none
         integer :: nx, ny
         real :: ci
         real, dimension(:,:,:), pointer :: cu, axy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: isign = 1, nxv, ny2d
         real :: ax, ay, affp, wm
         real, dimension(1,1,1) :: daxy0
         nxv = size(cu,2)/2; ny2d = size(cu,3)
         select case(size(cu,1))
!        case (2)
!           call APOISDX22(cu,daxy0,axy,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
!    &,nxv,ny2d)
         case (3)
            call APOISDX23(cu,daxy0,axy,isign,ffd,ax,ay,affp,ci,wm,nx,ny&
     &,nxv,ny2d)
         end select
         end subroutine isapoisdx23
!
      end module dfield2d
