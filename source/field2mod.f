!-----------------------------------------------------------------------
!
      module field2d
!
! Fortran90 interface to 2d PIC Fortran77 library field2lib.f
! field2mod.f contains procedures to manage guard cells and solve fields
!             in fourier space:
!             defines module field2d
! cguard => icguard2 copy guard cells for scalar or 2 and 3 component
!           vector arrays with various interpolations.
!           calls CGUARD2, BGUARD2, DGUARD2, CGUARD2L, BGUARD2L,
!           DGUARD2L, or CGUARD2C
! cguard => idguard2 copy guard cells for scalar arrays with various
!           interpolations.
!           calls DGUARD2, DGUARD2L, or DGUARD2C
! sguard => isguard2 initialize field for scalar array with various
!           interpolations.
!           calls SGUARD2, SGUARD2L, or SGUARD2C
! sguard => iscguard2 initialize field for 2 or 3 component vector array
!           with various interpolations.
!           calls SCGUARD2, SCGUARD22, SCGUARD2L, SCGUARD2L, SCGUARD2C,
!           or SCGUARD22C
! sguard => isfguard2 initialize 2 or 3 component field with scaled
!           vector array with various interpolations.
!           calls SCFGUARD2, SCFGUARD22, SCFGUARD2L, or SCFGUARD22L
! sguard => ismcguard2 initialize 2d tensor field  with various
!           interpolations.
!           calls SMCGUARD2, SMCGUARD22, SMCGUARD2L, SMCGUARD22L,
!           SMCGUARD2C,  or SMCGUARD22C
! aguard => iaguard2 add guard cells for scalar array with various
!           interpolations.
!           calls AGUARD2,  AGUARD2L, or AGUARD2C
! aguard => iacguard2 add guard cells for 2 or 3 component vector array
!           with various interpolations.
!           calls ACGUARD2, ACGUARD22, ACGUARD2L, ACGUARD22L, or
!           ACGUARD2C
! amcguard => iamcguard2 add guard cells for 2d tensor field with
!             various interpolations.
!             calls AMCGUARD2, AMCGUARD2LL, or ACGUARD2C
! pois_init => ipois22init initializes tables for field solvers.
!              calls POIS22
! pois => ipois2 solves poisson equation for electric force, potential,
!         or smoothing.
!         calls POISP2
! pois => ipois22 solves 2d poisson equation for electric force.
!         calls POIS22
! spois => ispois2 smoother for 2d periodic scalar field.
!          calls POISP2
! pois3 => ipois23 solves 2-1/2d poisson equation for electric force.
!          calls POIS22, or POIS23
! cuperp => icuperp2 calculates the transverse part of periodic 2-1/2d
!           vector field.
!           calls CUPERP2, or CUPERP22
! bpois => jbpois23 solves 2-1/2d vector poisson equation for magnetic
!          force.
!          calls BPOIS23, or BPOIS22
! sbpois => sbpois23 smoother for 2-1/2d periodic vector field.
!           calls BPOIS23, or BPOIS22
! apois => iapois23 solves 2-1/2d vector poisson equation for vector
!          potential.
!          calls BPOIS23, or BPOIS22
! ibpois => iibpois23 solves vector poisson equation for magnetic field.
!           calls IBPOIS23
! maxwel => imaxwel2 solves maxwell equation for electric and magnetic
!           fields.
!           calls MAXWEL2
! emfield => iemfield2 calculates periodic electric and magnetic forces
!            from fields given by maxwell and poisson equations.
!            calls EMFIELD2
! emfieldr => iemfieldr2 calculates electric and magnetic forces from
!             fields given by maxwell and poisson equations for
!             sine-cosine transforms.
!             calls EMFIELDR2
! emfieldc => iemfieldc2 calculates electric and magnetic forces from
!             fields given by maxwell and poisson equations for
!             sine-cosine/periodic transforms.
!             calls EMFIELDC2
! avpot => iavpot23 calculates vector potential from magnetic field.
!          calls AVPOT23
! avrpot => iavrpot23 calculates radiative vector potential from
!           magnetic field and current.
!           calls AVRPOT23
! gtmodes => igtmodes2 extracts selected fourier components from
!            potential array.
!            calls GTMODES2
! gtmodes => igtvmodes2 extracts selected fourier components from vector
!            potential array.
!            calls GTVMODES2
! ptmodes => iptmodes2 places selected fourier components into potential
!            array.
!            calls PTMODES2
! ptmodes => iptvmodes2 places selected fourier components into vector
!            potential array.
!            calls PTVMODES2
! poynt => ipoynt2 calculates the momentum in the electromagnetic field.
!          calls POYNT2
! poynt => idpoynt2 calculates the momentum in the darwin field.
!          calls DPOYNT2, or DPOYNT22
! dcuperp => idcuperp23 calculate transverse derivative of 2-1/2d
!            current density from momentum flux.
!            calls DCUPERP23, or DCUPERP22
! adcuperp => iadcuperp23 calculate transverse derivative of 2-1/2d
!             current density from momentum flux and acceleration
!             density.
!             calls ADCUPERP23, or ADCUPERP22
! epois_init => iepois23init initializes tables for darwin field solver.
!               calls EPOIS23
! epois => iepois23 solves 2-1/2d vector vector poisson equation for
!          transverse electric force.
!          calls EPOIS23, or EPOIS22
! iepois => iiepois23 solves 2-1/2d vector vector poisson equation for
!           transverse electric field.
!           calls EPOIS23, or EPOIS22
! dapois => idapois23 calculates derivatives of smoothed 2-1/2d vector
!           potential.
!           calls APOIS23, or APOIS22
! sapois => isapois23 calculates smoothed 2-1/2d vector potential.
!           calls APOIS23, or APOIS22
! addqei => iaddqei2 adds electron and ion densities.
!           calls ADDQEI2
! addqei => iaddqei2x adds electron and ion densities, and calculates
!           maximum and minimum plasma frequency.
!           calls ADDQEI2X
! baddext => ibaddext2 adds constant to magnetic field in real space for
!            2-1/2d code.
!            calls BADDEXT2, or BADDEXT22
! imoment => iimoment2 calculates ion momentum from integral of qi*fxy.
!            calls IMOMENT2, or IMOMENT22
! addfields => iaddvrfield2 calculates a = b + c for real vector fields.
!              calls ADDVRFIELD2

!****
! cushift => calls icushift to shift current array by half cell in fourier
!            space. calls CUSHIFT2.
! parforce2 => calls iparforce to calculate parallel part of fxyz for
!               debugging. calls PARFORCE2.
! modechop2 => calls imodeshop2 for zeroing out high k modes for
!               debugging. calls MODECHOP2
! avcguard => calls iavcguagd2 for adding the virtual current guard
!             cells. calls AVCGUARD2, AVCGUARD2L, or AVCGUARD2L
!             dependent on weighting scheme
! fldmod_init => calls ifldmodinit2 to shift and modify fields at the
!                 initial time step. calls FLDMODINIT2.
! fldshift => calls ifldshift2 to shift the fields between integer grid
!             and their respective yee positions. calls FLDSHIFT2.
!****
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: november 15, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC
      public :: DGUARD2, CGUARD2, BGUARD2, DGUARD2L, CGUARD2L, BGUARD2L
      public :: DCGUARD2, DCGUARD2L
      public :: SGUARD2, SGUARD2L, SGUARD2C, AGUARD2, AGUARD2L, AGUARD2C
      public :: SCGUARD2, SCGUARD22, SCGUARD2L, SCGUARD22L
      public :: ACGUARD2, ACGUARD22, ACGUARD2L, ACGUARD22L
      public :: SMCGUARD2, SMCGUARD22, SMCGUARD2L, SMCGUARD22L
      public :: DGUARD2C, CGUARD2C, SCGUARD2C, SCGUARD22C, ACGUARD2C
      public :: AMCGUARD2, AMCGUARD2L
      public :: cguard, sguard, aguard, acguard
      public :: pois_init, pois, spois, pois3, cuperp, bpois, sbpois
      public :: ibpois, maxwel, emfield, emfieldr, emfieldc
      public :: apois, avpot, avrpot, gtmodes, ptmodes, poynt
      public :: amcguard, dcuperp, adcuperp, epois_init, epois, iepois
      public :: dapois, sapois, addqei, baddext
      public :: imoment, idivf2, igradf2, icurlf2, icurlf22, ilaplace2
      public :: ivccopy, addfields, ilsremfield2,ilsremfield3,ilsremfield4
      public :: isavcuz2, imaxz2
      public :: icushift2, cushift, iparforce2, parforce, imodechop2, modechop
      public :: avcguard, iavcguard2, fldmod_init, ifldmodinit2, fldshift
      public :: ifldshift2
      public :: maxwelbyee, maxweleyee, maxwelb, maxwele
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine CGUARD2(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         end subroutine
      end interface 
      interface
         subroutine DGUARD2(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine BGUARD2(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: bxy
         end subroutine
      end interface
      interface
         subroutine SCGUARD2(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine SCGUARD22(cu,xj0,yj0,nx,ny,nxe,nye)
         implicit none
         real :: xj0, yj0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine SGUARD2(q,qi0,nx,ny,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine ACGUARD2(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine ACGUARD22(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine AGUARD2(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine CGUARD2L(fxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: fxy
         end subroutine
      end interface
      interface
         subroutine DGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: bxy
         end subroutine
      end interface
      interface
         subroutine SCGUARD2L(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine SCGUARD22L(cu,xj0,yj0,nx,ny,nxe,nye)
         implicit none
         real :: xj0, yj0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine SGUARD2L(q,qi0,nx,ny,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine ACGUARD22L(cu,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine AGUARD2L(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine CGUARD2C(fxy,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: fxy
         end subroutine
      end interface
      interface
         subroutine DGUARD2C(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
        subroutine SCGUARD2C(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
         implicit none
         real :: xj0, yj0, zj0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine SCGUARD22C(cu,xj0,yj0,nx,ny,nxe,nye)
         implicit none
         real :: xj0, yj0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine SGUARD2C(q,qi0,nx,ny,nxe,nye)
         implicit none
         real :: qi0
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine ACGUARD2C(cu,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: cu
         end subroutine
      end interface
      interface
         subroutine AGUARD2C(q,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real, dimension(nxe,nye) :: q
         end subroutine
      end interface
      interface
         subroutine POISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,nxv,nyh&
     &d)
         implicit none
         integer :: isign, nx, ny, nxv, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fx, fy
         real :: q, fx, fy
         complex, dimension(nxv/2,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine POIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,&
     &nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, we
!        real, dimension(*) :: q, fxy
         real :: q, fxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,&
     &nxhd,nyhd,yee)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd, yee
         real :: ax, ay, affp, we
!        real, dimension :: q, fxy
         real :: q, fxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine DIVF2(f,df,nx,ny,ndim,nxvh,nyv)
         implicit none
         integer :: nx, ny, ndim, nxvh, nyv
!        real, dimension(*) :: f, df
         real :: f, df
         end subroutine
      end interface
      interface
         subroutine GRADF2(df,f,nx,ny,ndim,nxvh,nyv)
         implicit none
         integer :: nx, ny, ndim, nxvh, nyv
!        real, dimension(*) :: f, df
         real :: f, df
         end subroutine
      end interface
      interface
         subroutine CURLF2(f,g,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: f, g
         real :: f, g
         end subroutine
      end interface
      interface
         subroutine CURLF22(f,g,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: f, g
         real :: f, g
         end subroutine
      end interface
      interface
         subroutine LAPLACE23(f,g,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: f, g
         real :: f, g
         end subroutine
      end interface
      interface
         subroutine LAPLACE22(f,g,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: f, g
         real :: f, g
         end subroutine
      end interface
      interface
         subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
!      interface
!         subroutine CUSHIFT2(cu,nx,ny,nxvh,nyv)
!         implicit none
!         integer :: nx, ny, nxvh, nyv
!         real :: cu
!         end subroutine
!      end interface
      interface
         subroutine CUPERP22(cu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: cu
         real :: cu
         end subroutine
      end interface
      interface
         subroutine BPOIS23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,nxvh&
     &,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy
         real :: cu, bxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine BPOIS22(cu,bxy,bz,isign,ffc,ax,ay,affp,ci,wm,nx,ny,n&
     &xvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, bxy, bz
         real :: cu, bxy, bz
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine IBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd,yee)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd,yee
         real :: ci, wm
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(3,nxvh,nyv) :: bxy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine MAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nx&
     &hd,nyhd,parmax,yee)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd, parmax, yee
         real :: ci, dt, wf, wm
         complex, dimension(3,nxvh,nyv) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!*****
      interface
         subroutine MAXWELBYEE2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nx&
     &hd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, dt, wf, wm
         complex, dimension(3,nxvh,nyv) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWELEYEE2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nx&
     &hd,nyhd,parmax)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd, parmax
         real :: ci, dt, wf, wm
         complex, dimension(3,nxvh,nyv) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!
      interface
         subroutine MAXWELB2(exy,bxy,ci,dt,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, dt, wm
         complex, dimension(3,nxvh,nyv) :: exy, bxy
         end subroutine
      end interface
!
      interface
         subroutine MAXWELE2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nx&
     &hd,nyhd,parmax)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd, parmax
         real :: ci, dt, wf, wm
         complex, dimension(3,nxvh,nyv) :: exy, bxy
!        real, dimension(*) :: cu
         real :: cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
!*****
      interface
         subroutine EMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real, dimension(*) :: fxy
         real :: fxy
         complex, dimension(3,nxvh,nyv) :: exy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine LSREMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real, dimension(*) :: fxy
         real :: fxy
         complex, dimension(3,nxvh,nyv) :: exy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine BDEMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real, dimension(*) :: fxy
         real :: fxy
         complex, dimension(3,nxvh,nyv) :: exy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine EMFIELDR2(fxy,exy,ffd,isign,nx,ny,nxv,nxe,nye)
         implicit none
         integer :: isign, nx, ny, nxv, nxe, nye
!        real, dimension(*) :: fxy
         real :: fxy
         complex, dimension(3,nxe/2,nye) :: exy
         complex, dimension(nxv,ny) :: ffd
         end subroutine
      end interface
      interface
         subroutine EMFIELDC2(fxy,exy,ffb,isign,nx,ny,nxv,nxe,nyeh)
         implicit none
         integer :: isign, nx, ny, nxv, nxe, nyeh
!        real, dimension(*) :: fxy
         real :: fxy
         complex, dimension(3,nxe,nyeh) :: exy
         complex, dimension(nxv,ny/2) :: ffb
         end subroutine
      end interface
      interface
         subroutine AVPOT23(bxy,axy,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
         complex, dimension(3,nxvh,nyv) :: bxy
!        real, dimension(*) :: axy
         real :: axy
         end subroutine
      end interface
      interface
         subroutine AVRPOT23(axy,bxy,ffc,ci,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci
         complex, dimension(3,nxvh,nyv) :: bxy
!        real, dimension(*) :: axy
         real :: axy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine GTMODES2(pot,pott,nx,ny,it,modesx,modesy,nxe,nye,nt2&
     &,modesxd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, nxe, nye
         integer :: nt2, modesxd, modesyd
!        real, dimension(*) :: pot
         real :: pot
         complex, dimension(nt2/2,modesxd,modesyd) :: pott
         end subroutine
      end interface
      interface
         subroutine PTMODES2(pot,pott,nx,ny,it,modesx,modesy,nxe,nye,nt2&
     &,modesxd,modesyd)  
         implicit none
         integer :: nx, ny, it, modesx, modesy, nxe, nye
         integer :: nt2, modesxd, modesyd
!        real, dimension(*) :: pot
         real :: pot
         complex, dimension(nt2/2,modesxd,modesyd) :: pott
         end subroutine
      end interface
      interface
         subroutine GTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,nxv&
     &h,nyv,nt,modesxd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, ndim, nxvh, nyv
         integer :: nt, modesxd, modesyd
!        complex, dimension(*) :: vpot
         real :: vpot
         complex, dimension(nt,ndim,modesxd,modesyd) :: vpott
         end subroutine
      end interface
      interface
         subroutine PTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,nxv&
     &h,nyv,nt,modesxd,modesyd)
         implicit none
         integer :: nx, ny, it, modesx, modesy, ndim, nxvh, nyv
         integer :: nt, modesxd, modesyd
!        complex, dimension(*) :: vpot
         real :: vpot
         complex, dimension(nt,ndim,modesxd,modesyd) :: vpott
         end subroutine
      end interface
      interface
         subroutine POYNT2(q,exy,bxy,ffc,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,ny&
     &hd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: sx, sy, sz
         complex, dimension(3,nxvh,nyv) :: exy, bxy
!        real, dimension(*) :: q
         real :: q
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine DPOYNT2(q,cu,ffc,ci,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,nyh&
     &d)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, sx, sy, sz
!        real, dimension(*) :: q, cu
         real :: q, cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine DPOYNT22(q,cu,ffc,ci,sx,sy,nx,ny,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ci, sx, sy
!        real, dimension(*) :: q, cu
         real :: q, cu
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine SCFGUARD2(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine SCFGUARD22(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine SMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine SMCGUARD22(amu,x2y2m0,xym0,nx,ny,nxe,nye)
         implicit none
         real :: x2y2m0, xym0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine AMCGUARD2(amu,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine SCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(3,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine SCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
         implicit none
         real :: q2m0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: cus, cu
         end subroutine
      end interface
      interface
         subroutine SMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine SMCGUARD22L(amu,x2y2m0,xym0,nx,ny,nxe,nye)
         implicit none
         real :: x2y2m0, xym0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine AMCGUARD2L(amu,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine SMCGUARD2C(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
         implicit none
         real :: x2y2m0, xym0, zxm0, zym0
         integer :: nx, ny, nxe, nye
         real, dimension(4,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine SMCGUARD22C(amu,x2y2m0,xym0,nx,ny,nxe,nye)
         implicit none
         real :: x2y2m0, xym0
         integer :: nx, ny, nxe, nye
         real, dimension(2,nxe,nye) :: amu
         end subroutine
      end interface
      interface
         subroutine DCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine DCUPERP22(dcu,amu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine ADCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
      interface
         subroutine ADCUPERP22(dcu,amu,nx,ny,nxvh,nyv)
         implicit none
         integer :: nx, ny, nxvh, nyv
!        real, dimension(*) :: dcu, amu
         real :: dcu, amu
         end subroutine
      end interface
!     interface
!        subroutine EPOIS23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny&
!    &,nxvh,nyv,nxhd,nyhd)
!        implicit none
!        integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real :: ax, ay, affp, wp0, ci, wf
!        real, dimension(*) :: dcu, exy
!        real :: dcu, exy
!        complex, dimension(nxhd,nyhd) :: ffe
!        end subroutine
!     end interface
!     interface
!        subroutine EPOIS22(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny&
!    &,nxvh,nyv,nxhd,nyhd)
!        implicit none
!        integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
!        real :: ax, ay, affp, wp0, ci, wf
!        real, dimension(*) :: dcu, exy
!        real :: dcu, exy
!        complex, dimension(nxhd,nyhd) :: ffe
!        end subroutine
!     end interface
      interface
         subroutine DCGUARD2(daxy,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: daxy
         end subroutine
      end interface
      interface
         subroutine DCGUARD2L(daxy,nx,ny,nxe,nye,ndim)
         implicit none
         integer :: nx, ny, nxe, nye, ndim
         real, dimension(ndim,nxe,nye) :: daxy
         end subroutine
      end interface
      interface
         subroutine APOIS23(cu,daxy,axy,isign,ffc,ax,ay,affp,ci,wm,nx,ny&
     &,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, daxy, axy
         real :: cu, daxy, axy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine APOIS22(cu,daxy,axy,isign,ffc,ax,ay,affp,ci,wm,nx,ny&
     &,nxvh,nyv,nxhd,nyhd)
         implicit none
         integer :: isign, nx, ny, nxvh, nyv, nxhd, nyhd
         real :: ax, ay, affp, ci, wm
!        real, dimension(*) :: cu, daxy, axy
         real :: cu, daxy, axy
         complex, dimension(nxhd,nyhd) :: ffc
         end subroutine
      end interface
      interface
         subroutine ADDQEI2(qe,qi,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
!        real, dimension(*,*) :: qe, qi
         real :: qe, qi
         end subroutine
      end interface
      interface
         subroutine ADDQEI2X(qe,qi,qbme,qbmi,wpmax,wpmin,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: qbme, qbmi, wpmax, wpmin
!        real, dimension(*,*) :: qe, qi
         real :: qe, qi
         end subroutine
      end interface
      interface
         subroutine BADDEXT2(bxy,omx,omy,omz,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: omx, omy, omz
!        real, dimension(*,*,*) :: bxy
         real :: bxy
         end subroutine
      end interface
      interface
         subroutine BADDEXT22(bz,omz,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: omz
!        real, dimension(*,*) :: bz
         real :: bz
         end subroutine
      end interface
      interface
         subroutine IMOMENT2(qi,fxy,pxi,pyi,pzi,dt,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: pxi, pyi, pzi, dt
!        real, dimension(*) :: qi, fxy
         real :: qi, fxy
         end subroutine
      end interface
      interface
         subroutine IMOMENT22(qi,fxy,pxi,pyi,dt,nx,ny,nxe,nye)
         implicit none
         integer :: nx, ny, nxe, nye
         real :: pxi, pyi, dt
!        real, dimension(*) :: qi, fxy
         real :: qi, fxy
         end subroutine
      end interface
      interface
         subroutine VCCOPY2(f,g,nx,ny,ndim,nxv,nyv)
         implicit none
         integer nx, ny, ndim, nxv, nyv
         complex g
!        dimension f(*)
         real f
         dimension g(ndim,nxv,nyv)
         end subroutine
      end interface
      interface
         subroutine ADDVRFIELD2(a,b,c,ndim,nxe,nye)
         implicit none
         integer :: ndim, nxe, nye
         real, dimension(ndim,nxe,nye) :: a, b, c
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface cguard
         module procedure icguard2
         module procedure idguard2
      end interface
!
      interface bguard
         module procedure ibguard2
      end interface
!
      interface sguard
         module procedure isguard2
         module procedure iscguard2
         module procedure isfguard2
         module procedure ismcguard2
      end interface
!
      interface aguard
         module procedure iaguard2
         module procedure iacguard2
      end interface
!
      interface acguard
         module procedure iacguard2
      end interface
!
      interface avcguard
         module procedure iavcguard2
      end interface
!
       interface pois_init
         module procedure ipois22init
      end interface
!
      interface pois
         module procedure ipois2
         module procedure ipois22
      end interface
!
      interface spois
         module procedure ispois2
      end interface
!
      interface pois3
         module procedure ipois23
      end interface
!
      interface cuperp
         module procedure icuperp2
      end interface
!
      interface modechop
         module procedure imodechop2
      end interface
 
      interface parforce
         module procedure iparforce2
      end interface
!
      interface cushift
         module procedure icushift2
      end interface
!
      interface bpois
         module procedure jbpois23
      end interface
!
      interface sbpois
         module procedure sbpois23
      end interface
!
      interface apois
         module procedure iapois23
      end interface
!
      interface ibpois
         module procedure iibpois23
      end interface
!
      interface maxwel
         module procedure imaxwel2
      end interface
!*****
      interface maxwelbyee
         module procedure imaxwelbyee2
      end interface
!
      interface maxweleyee
         module procedure imaxweleyee2
      end interface
      interface maxwelb
         module procedure imaxwelb2
      end interface
!
      interface maxwele
         module procedure imaxwele2
      end interface
!*****
      interface emfield
         module procedure iemfield2
      end interface
!
      interface emfieldr
         module procedure iemfieldr2
      end interface
!
      interface emfieldc
         module procedure iemfieldc2
      end interface
!
      interface avpot
         module procedure iavpot23
      end interface
!
      interface avrpot
         module procedure iavrpot23
      end interface
!
      interface gtmodes
         module procedure igtmodes2
         module procedure igtvmodes2
      end interface
!
      interface ptmodes
         module procedure iptmodes2
         module procedure iptvmodes2
      end interface
!
      interface poynt
         module procedure ipoynt2
         module procedure idpoynt2
      end interface
!
      interface amcguard
         module procedure iamcguard2
      end interface
!
      interface dcuperp
         module procedure idcuperp23
      end interface
!
      interface adcuperp
         module procedure iadcuperp23
      end interface
!
      interface fldshift
         module procedure ifldshift2
      end interface
!
      interface fldmod_init
         module procedure ifldmodinit2
      end interface
!
       interface epois_init
         module procedure iepois23init
      end interface
!
      interface epois
         module procedure iepois23
      end interface
!
      interface iepois
         module procedure iiepois23
      end interface
!
      interface dapois
         module procedure idapois23
      end interface
!
      interface sapois
         module procedure isapois23
      end interface
!
      interface addqei
         module procedure iaddqei2
         module procedure iaddqei2x
      end interface
!
      interface baddext
         module procedure ibaddext2
      end interface
!
      interface imoment
         module procedure iimoment2
      end interface
!
      interface ivccopy
         module procedure ivccopy2
      end interface
!
      interface addfields
         module procedure iaddvrfield2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine icguard2(fxy,nx,ny,inorder)
! copy guard cells for periodic 2d vector data, for N component vectors
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: nxe, nye, ndim, order
         nxe = size(fxy,2); nye = size(fxy,3)
         ndim = size(fxy,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==CUBIC) then
            call CGUARD2C(fxy,nx,ny,nxe,nye,ndim)
            return
         endif
         select case(size(fxy,1))
         case (1)
               !print*, 'here1'
            if (order==LINEAR) then
               call DGUARD2L(fxy,nx,ny,nxe,nye)
            else
               call DGUARD2(fxy,nx,ny,nxe,nye)
            endif
         case (2)
               !print*, 'here2'
            if (order==LINEAR) then
               call CGUARD2L(fxy,nx,ny,nxe,nye)
            else
               call CGUARD2(fxy,nx,ny,nxe,nye)
            endif
         case (3)
               !print*, 'here3'
            if (order==LINEAR) then
               !print*, 'here'
               call BGUARD2L(fxy,nx,ny,nxe,nye)
            else
               !print*, 'there'
               call BGUARD2(fxy,nx,ny,nxe,nye)
            endif
         case (4:)
               !print*, 'here4'
            if (order==LINEAR) then
               call DCGUARD2L(fxy,nx,ny,nxe,nye,ndim)
            else
               call DCGUARD2(fxy,nx,ny,nxe,nye,ndim)
            endif
         end select
         end subroutine icguard2
!
         subroutine idguard2(q,nx,ny,inorder)
! copy guard cells for periodic 2d scalar data
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
            call DGUARD2L(q,nx,ny,nxe,nye)
         else if (order==CUBIC) then
            call DGUARD2C(q,nx,ny,nxe,nye)
         else
            call DGUARD2(q,nx,ny,nxe,nye)
         endif
         end subroutine idguard2
!
         subroutine ibguard2(bxy,nx,ny,inorder)
! copy guard cells for periodic 2d vector data, for 3 component vectors
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: bxy
! local data
         integer :: nxe, nye, order
         nxe = size(bxy,2); nye = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call BGUARD2L(bxy,nx,ny,nxe,nye)
         else
            call BGUARD2(bxy,nx,ny,nxe,nye)
         endif
         end subroutine ibguard2
!
!                 call sguard(cu,zero,zero,zero,nx,ny,inorder)
         subroutine iscguard2(cu,xj0,yj0,zj0,nx,ny,inorder)
! initialize periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: xj0, yj0, zj0
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nye, order
         nxe = size(cu,2); nye = size(cu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call SCGUARD22L(cu,xj0,yj0,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call SCGUARD22C(cu,xj0,yj0,nx,ny,nxe,nye)
            else
               call SCGUARD22(cu,xj0,yj0,nx,ny,nxe,nye)
            endif
         case (3)
            if (order==LINEAR) then
               call SCGUARD2L(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call SCGUARD2C(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
            else
               call SCGUARD2(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
            endif
         end select
         end subroutine iscguard2
!
         subroutine isguard2(q,qi0,nx,ny,inorder)
! initialize periodic 2d scalar field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: qi0
         real, dimension(:,:), pointer :: q
! local data
         integer :: nxe, nye, order
         nxe = size(q,1); nye = size(q,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call SGUARD2L(q,qi0,nx,ny,nxe,nye)
         else if (order==CUBIC) then
            call SGUARD2C(q,qi0,nx,ny,nxe,nye)
         else
            call SGUARD2(q,qi0,nx,ny,nxe,nye)
         endif
         end subroutine isguard2
!
         subroutine iacguard2(cu,nx,ny,inorder)
! add guard cells for periodic 2d vector data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxe, nye, ndim, order
         nxe = size(cu,2); nye = size(cu,3)
         ndim = size(cu,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==CUBIC) then
            call ACGUARD2C(cu,nx,ny,nxe,nye,ndim)
            return
         endif
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call ACGUARD22L(cu,nx,ny,nxe,nye)
            else
               call ACGUARD22(cu,nx,ny,nxe,nye)
            endif
         case (3)
            if (order==LINEAR) then
               call ACGUARD2L(cu,nx,ny,nxe,nye)
            else
               call ACGUARD2(cu,nx,ny,nxe,nye)
            endif
         end select
         end subroutine iacguard2
!
         subroutine iavcguard2(cu,vcu,nx,ny,inorder)
! add guard cells for the virtual curre
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu, vcu
! local data
         integer :: nxe, nye, ndim, order
         nxe = size(cu,2); nye = size(cu,3)
         ndim = size(cu,1)
         order = QUADRATIC
         if (present(inorder)) order = inorder
!         if (order==CUBIC) then
!            call ACGUARD2C(cu,nx,ny,nxe,nye,ndim)
!            return
!         endif
!         select case(size(cu,1))
!         case (2)
!            if (order==LINEAR) then
!               call ACGUARD22L(cu,nx,ny,nxe,nye)
!            else
!               call ACGUARD22(cu,nx,ny,nxe,nye)
!            endif
!         case (3)
            if (order==LINEAR) then
               call AVCGUARD2L(cu,vcu,nx,ny,nxe,nye)
            else if(order==QUADRATIC) then
               call AVCGUARD2(cu,vcu,nx,ny,nxe,nye)
            else if(order==CUBIC) then
               call AVCGUARD2C(cu,vcu,nx,ny,nxe,nye)
            endif
!         end select
         end subroutine iavcguard2
!
         subroutine iaguard2(q,nx,ny,inorder)
! add guard cells for periodic 2d scalar data
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
            call AGUARD2L(q,nx,ny,nxe,nye)
         else if (order==CUBIC) then
            call AGUARD2C(q,nx,ny,nxe,nye)
         else
            call AGUARD2(q,nx,ny,nxe,nye)
         endif
         end subroutine iaguard2
!
         subroutine ipois2(q,fx,ffc,we,nx,ny,inorder)
! poisson solver for periodic 2d potential
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: we
         real, dimension(:,:), pointer :: q, fx
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = 1, nxv, nyhd, order
         real :: ax, ay, affp
         real, dimension(1,1) :: fy
         nxv = size(q,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call POISP2(q(1,1),fx(1,1),fy(1,1),isign,ffc,ax,ay,affp,we,n&
     &x,ny,nxv,nyhd)
         else if (order==CUBIC) then
            call POISP2(q(3,3),fx(3,3),fy(1,1),isign,ffc,ax,ay,affp,we,n&
     &x,ny,nxv,nyhd)
         else
            call POISP2(q(2,2),fx(2,2),fy(1,1),isign,ffc,ax,ay,affp,we,n&
     &x,ny,nxv,nyhd)
         endif
         end subroutine ipois2
!
         subroutine ispois2(q,fy,ffc,nx,ny,inorder)
! smoother for periodic 2d scalar field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: q, fy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = 2, nxv, nyhd, order
         real :: ax, ay, affp, we
         real, dimension(1,1) :: fx
         nxv = size(q,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call POISP2(q(1,1),fx(1,1),fy(1,1),isign,ffc,ax,ay,affp,we,n&
     &x,ny,nxv,nyhd)
         else if (order==CUBIC) then
            call POISP2(q(3,3),fx(1,1),fy(3,3),isign,ffc,ax,ay,affp,we,n&
     &x,ny,nxv,nyhd)
         else
            call POISP2(q(2,2),fx(1,1),fy(2,2),isign,ffc,ax,ay,affp,we,n&
     &x,ny,nxv,nyhd)
         endif
         end subroutine ispois2
!
         subroutine ipois22init(ffc,ax,ay,affp,nx,ny,smoothtype,smoothfactor,&
      &smoothmask)
! initialize 2d periodic electric field solver
         implicit none
         integer :: nx, ny, smoothtype
         real :: ax, ay, affp, smoothfactor, smoothmask
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = 0, nxvh = 1, nyv = 1, nxhd, nyhd
         real :: we, q, fxy
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         if (smoothtype == 1) then
            call POIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         elseif (smoothtype == 2) then
            call POIS22cf(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nxhd,n&
     &yhd,smoothfactor,smoothmask)
         else
            print *, 'smooth type error'
         endif
         end subroutine ipois22init
!
         subroutine ipois22(q,fxy,ffc,we,tfield,nx,ny,inorder,yee)
! poisson solver for periodic 2d electric field
         implicit none
         integer :: nx, ny, yee
         integer, optional :: inorder
         real :: we, tfield
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = -1, nxv, nxvh, nyv, nxhd, nyhd, ltime, order
         real :: ax, ay, affp, tf
         nxv = size(q,1); nxvh = nxv/2; nyv = size(q,2)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         if (order==LINEAR) then
            call POIS22(q(1,1),fxy(1,1,1),isign,ffc,ax,ay,affp,we,nx,ny,&
     &nxvh,nyv,nxhd,nyhd)
         else if (order==CUBIC) then
            call POIS22(q(3,3),fxy(1,3,3),isign,ffc,ax,ay,affp,we,nx,ny,&
     &nxvh,nyv,nxhd,nyhd)
         else
            call POIS22(q(2,2),fxy(1,2,2),isign,ffc,ax,ay,affp,we,nx,ny,&
     &nxvh,nyv,nxhd,nyhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine ipois22
!
         subroutine ipois23(q,fxy,ffc,we,tfield,nx,ny,inorder,yee)
! poisson solver for periodic 2-1/2d electric field
         implicit none
         integer :: nx, ny, yee
         integer, optional :: inorder
         real :: we, tfield
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = -1, nxv, nxvh, nyv, nxhd, nyhd, ltime, order
         real :: ax, ay, affp, tf
         nxv = size(q,1); nxvh = nxv/2; nyv = size(q,2)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(fxy,1))
         case (2)
            if (order==LINEAR) then
               call POIS22(q(1,1),fxy(1,1,1),isign,ffc,ax,ay,affp,we,nx,&
     &ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call POIS22(q(3,3),fxy(1,3,3),isign,ffc,ax,ay,affp,we,nx,&
     &ny,nxvh,nyv,nxhd,nyhd)
            else
               call POIS22(q(2,2),fxy(1,2,2),isign,ffc,ax,ay,affp,we,nx,&
     &ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            print*, 'using POIS23 solver for parallel field'
            if (order==LINEAR) then
               call POIS23(q(1,1),fxy(1,1,1),isign,ffc,ax,ay,affp,we,nx,&
     &ny,nxvh,nyv,nxhd,nyhd,yee)
            else if (order==CUBIC) then
               call POIS23(q(3,3),fxy(1,3,3),isign,ffc,ax,ay,affp,we,nx,&
     &ny,nxvh,nyv,nxhd,nyhd,yee)
            else
               call POIS23(q(2,2),fxy(1,2,2),isign,ffc,ax,ay,affp,we,nx,&
     &ny,nxvh,nyv,nxhd,nyhd,yee)
            endif
         end select
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine ipois23
!
         subroutine idivf2(f,df,nx,ny,inorder)
! calculates the divergence of periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
         real, dimension(:,:), pointer :: df
! local data
         integer :: ndim, nxvh, nyv, order
         ndim = size(f,1)
         nxvh = size(f,2)/2; nyv = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DIVF2(f(1,1,1),df(1,1),nx,ny,ndim,nxvh,nyv)
         else if (order==CUBIC) then
            call DIVF2(f(1,3,3),df(3,3),nx,ny,ndim,nxvh,nyv)
         else
            call DIVF2(f(1,2,2),df(2,2),nx,ny,ndim,nxvh,nyv)
         endif
         end subroutine idivf2
!
         subroutine igradf2(df,f,nx,ny,inorder)
! calculates the gradient of periodic 2d scalar field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: df
         real, dimension(:,:,:), pointer :: f
! local data
         integer :: ndim, nxvh, nyv, order
         ndim = size(f,1)
         nxvh = size(df,1)/2; nyv = size(df,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call GRADF2(df(1,1),f(1,1,1),nx,ny,ndim,nxvh,nyv)
         else if (order==CUBIC) then
            call GRADF2(df(3,3),f(1,3,3),nx,ny,ndim,nxvh,nyv)
         else
            call GRADF2(df(2,2),f(1,2,2),nx,ny,ndim,nxvh,nyv)
         endif
         end subroutine igradf2
!
         subroutine icurlf2(f,g,nx,ny,inorder)
! calculates the curl of periodic 2-1/2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f, g
! local data
         integer :: nxvh, nyv, order
         nxvh = size(f,2)/2; nyv = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call CURLF2(f(1,1,1),g(1,1,1),nx,ny,nxvh,nyv)
         else if (order==CUBIC) then
            call CURLF2(f(1,3,3),g(1,3,3),nx,ny,nxvh,nyv)
         else
            call CURLF2(f(1,2,2),g(1,2,2),nx,ny,nxvh,nyv)
         endif
         end subroutine icurlf2
!
         subroutine icurlf22(f,g,nx,ny,inorder)
! calculates the curl of periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f
         real, dimension(:,:), pointer :: g
! local data
         integer :: nxvh, nyv, order
         nxvh = size(f,2)/2; nyv = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call CURLF22(f(1,1,1),g(1,1),nx,ny,nxvh,nyv)
         else if (order==CUBIC) then
            call CURLF22(f(1,3,3),g(3,3),nx,ny,nxvh,nyv)
         else
            call CURLF22(f(1,2,2),g(2,2),nx,ny,nxvh,nyv)
         endif
         end subroutine icurlf22
!
         subroutine ilaplace2(f,g,nx,ny,inorder)
! calculates the laplacian of periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: f, g
! local data
         integer :: nxvh, nyv, order
         nxvh = size(f,2)/2; nyv = size(f,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(f,1))
         case (3)
            if (order==LINEAR) then
               call LAPLACE23(f(1,1,1),g(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call LAPLACE23(f(1,3,3),g(1,3,3),nx,ny,nxvh,nyv)
            else
               call LAPLACE23(f(1,2,2),g(1,2,2),nx,ny,nxvh,nyv)
            endif
         case (2)
            if (order==LINEAR) then
               call LAPLACE22(f(1,1,1),g(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call LAPLACE22(f(1,3,3),g(1,3,3),nx,ny,nxvh,nyv)
            else
               call LAPLACE22(f(1,2,2),g(1,2,2),nx,ny,nxvh,nyv)
            endif
         end select
         end subroutine ilaplace2
!
         subroutine icuperp2(cu,nx,ny,inorder)
! calculates the transverse part of periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxvh, nyv, order
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call CUPERP22(cu(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call CUPERP22(cu(1,3,3),nx,ny,nxvh,nyv)
            else
               call CUPERP22(cu(1,2,2),nx,ny,nxvh,nyv)
            endif
         case (3)
            if (order==LINEAR) then
               call CUPERP2(cu(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call CUPERP2(cu(1,3,3),nx,ny,nxvh,nyv)
            else
               call CUPERP2(cu(1,2,2),nx,ny,nxvh,nyv)
            endif
         end select
         end subroutine icuperp2
!
         subroutine imodechop2(fxy,kxcutoff,kycutoff,nx,ny,inorder)
! chops off high k for a periodic 2d vector field
! for debugging
         implicit none
         integer :: nx, ny, kxcutoff, kycutoff
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: nxvh, nyv, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
            if (order==LINEAR) then
               call MODECHOP2(fxy(1,1,1),kxcutoff,kycutoff,nx,ny,nxvh,nyv)
            else if (order==QUADRATIC) then
              call MODECHOP2(fxy(1,2,2),kxcutoff,kycutoff,nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
              call MODECHOP2(fxy(1,3,3),kxcutoff,kycutoff,nx,ny,nxvh,nyv)
            endif
         end subroutine imodechop2
!
         subroutine ifldshift2(fxy,exy,bxy,isign,nx,ny)
         implicit none
         integer :: nx, ny, isign
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy, bxy
         integer :: nxvh, nyv, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
              call FLDSHIFT2(exy,bxy,isign,nx,ny,nxvh,nyv)
         end subroutine ifldshift2
!
         subroutine ifldmodinit2(fxy,exy,bxy,nx,ny,inorder)
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy, bxy
! local data
         integer :: nxvh, nyv, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
            if (order==LINEAR) then
               call FLDMODINIT2(fxy(1,1,1),exy,bxy,nx,ny,nxvh,nyv)
            else if (order==QUADRATIC) then
              call FLDMODINIT2(fxy(1,2,2),exy,bxy,nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
              call FLDMODINIT2(fxy(1,3,3),exy,bxy,nx,ny,nxvh,nyv)
            endif
         end subroutine ifldmodinit2
!
         subroutine iparforce2(fxy,fp,nx,ny,inorder)
! calculates the parallel part of periodic 2d vector field
! for debugging
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy, fp
! local data
         integer :: nxvh, nyv, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
!         select case(size(cu,1))
!         case (2)
!            if (order==LINEAR) then
!               call CUPERP22(cu(1,1,1),nx,ny,nxvh,nyv)
!            else if (order==CUBIC) then
!               call CUPERP22(cu(1,3,3),nx,ny,nxvh,nyv)
!            else
!               call CUPERP22(cu(1,2,2),nx,ny,nxvh,nyv)
!            endif
!         case (3)
            if (order==LINEAR) then
               call PARFORCE2(fxy(1,1,1),fp(1,1,1),nx,ny,nxvh,nyv)
            else if (order==QUADRATIC) then
              call PARFORCE2(fxy(1,2,2),fp(1,2,2),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
              call PARFORCE2(fxy(1,3,3),fp(1,3,3),nx,ny,nxvh,nyv)
            endif
!         end select
         end subroutine iparforce2
!   
         subroutine icushift2(cu,nx,ny,inorder)
! calculates the transverse part of periodic 2d vector field
         implicit none
         integer :: nx, ny, order!, yee
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: nxvh, nyv
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         !print*, 'order =', order
         if (order==LINEAR) then
           call CUSHIFT2(cu(1,1,1),nx,ny,nxvh,nyv)
         else if (order==QUADRATIC) then
           call CUSHIFT2(cu(1,2,2),nx,ny,nxvh,nyv)
         else if(order==CUBIC) then
           call CUSHIFT2(cu(1,3,3),nx,ny,nxvh,nyv)
         end if
         end subroutine icushift2
!
         subroutine jbpois23(cu,bxy,ffc,ci,wm,tfield,nx,ny,inorder)
! calculates static magnetic field for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, wm, tfield
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = -1, nxvh, nyv, nxhd, nyhd, ltime, order
         real :: ax, ay, affp, bxy0, tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call BPOIS22(cu(1,1,1),bxy0,bxy(1,1,1),isign,ffc,ax,ay,af&
     &fp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call BPOIS22(cu(1,3,3),bxy0,bxy(1,3,3),isign,ffc,ax,ay,af&
     &fp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call BPOIS22(cu(1,2,2),bxy0,bxy(1,2,2),isign,ffc,ax,ay,af&
     &fp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call BPOIS23(cu(1,1,1),bxy(1,1,1),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call BPOIS23(cu(1,3,3),bxy(1,3,3),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call BPOIS23(cu(1,2,2),bxy(1,2,2),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         end select
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine jbpois23
!
         subroutine sbpois23(cu,bxy,ffc,nx,ny,inorder)
! smoother for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu, bxy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = 2, nxvh, nyv, nxhd, nyhd, order
         real :: ax, ay, affp, bz0, ci, wm
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call BPOIS22(cu(1,1,1),bxy(1,1,1),bz0,isign,ffc,ax,ay,aff&
     &p,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call BPOIS22(cu(1,3,3),bxy(1,3,3),bz0,isign,ffc,ax,ay,aff&
     &p,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call BPOIS22(cu(1,2,2),bxy(1,2,2),bz0,isign,ffc,ax,ay,aff&
     &p,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call BPOIS23(cu(1,1,1),bxy(1,1,1),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call BPOIS23(cu(1,3,3),bxy(1,3,3),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call BPOIS23(cu(1,2,2),bxy(1,2,2),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         end select
         end subroutine sbpois23
!
         subroutine iapois23(cu,axy,ffc,ci,wm,nx,ny,inorder)
! calculates static vector potential for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu, axy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = 1, nxvh, nyv, nxhd, nyhd, order
         real :: ax, ay, affp, bz0
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call BPOIS22(cu(1,1,1),axy(1,1,1),bz0,isign,ffc,ax,ay,aff&
     &p,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call BPOIS22(cu(1,3,3),axy(1,3,3),bz0,isign,ffc,ax,ay,aff&
     &p,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call BPOIS22(cu(1,2,2),axy(1,2,2),bz0,isign,ffc,ax,ay,aff&
     &p,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call BPOIS23(cu(1,1,1),axy(1,1,1),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call BPOIS23(cu(1,3,3),axy(1,3,3),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call BPOIS23(cu(1,2,2),axy(1,2,2),isign,ffc,ax,ay,affp,ci&
     &,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         end select
         end subroutine iapois23
!
         subroutine iibpois23(cu,bxy,ffc,ci,wm,nx,ny,inorder,yee)
! calculates static magnetic field for periodic 2d vector field
         implicit none
         integer :: nx, ny, yee
         integer, optional :: inorder
         real :: ci, wm
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:,:), pointer :: bxy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call IBPOIS23(cu(1,1,1),bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,ny&
     &hd,yee)
         else if (order==CUBIC) then
            call IBPOIS23(cu(1,3,3),bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,ny&
     &hd,yee)
         else
            call IBPOIS23(cu(1,2,2),bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,ny&
     &hd,yee)
         endif
         end subroutine iibpois23
!
         subroutine imaxwel2(exy,bxy,cu,ffc,ci,dt,wf,wm,tfield,nx,ny,ino&
     &rder,solvertype,parmax,yee)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny, parmax, yee
         integer, optional :: inorder, solvertype
         real :: ci, dt, wf, wm, tfield
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, ltime, order
         real :: tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         if (solvertype == 1) then
          if(parmax==1) then
          print*, 'using MAXWEL2 solver for total field'
          else
          print *, 'using MAXWEL2 solver for perp. field'
          endif
          if (order==LINEAR) then
            call MAXWEL2(exy,bxy,cu(1,1,1),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax,yee)
          else if (order==CUBIC) then
            call MAXWEL2(exy,bxy,cu(1,3,3),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax,yee)
          else
            call MAXWEL2(exy,bxy,cu(1,2,2),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax,yee)
          endif
         elseif (solvertype == 2) then
          print *, 'using MAXWEL2PS solver'
          if (order==LINEAR) then
            call MAXWEL2PS(exy,bxy,cu(1,1,1),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd)
          else if (order==CUBIC) then
            call MAXWEL2PS(exy,bxy,cu(1,3,3),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd)
          else
            call MAXWEL2PS(exy,bxy,cu(1,2,2),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd)
          endif
         endif
          
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine imaxwel2
!*****
         subroutine imaxwelbyee2(exy,bxy,cu,ffc,ci,dt,wf,wm,tfield,nx,ny,ino&
     &rder)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, dt, wf, wm, tfield
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, ltime, order
         real :: tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
          print *, 'using MAXWELBYEE2 solver'
          if (order==LINEAR) then
            call MAXWELBYEE2(exy,bxy,cu(1,1,1),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd)
          else if (order==CUBIC) then
            call MAXWELBYEE2(exy,bxy,cu(1,3,3),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd)
          else
            call MAXWELBYEE2(exy,bxy,cu(1,2,2),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd)
          endif
          
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine imaxwelbyee2
!
         subroutine imaxweleyee2(exy,bxy,cu,ffc,ci,dt,wf,wm,tfield,nx,ny,ino&
     &rder,parmax)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny, parmax
         integer, optional :: inorder
         real :: ci, dt, wf, wm, tfield
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, ltime, order
         real :: tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
          print *, 'using MAXWELEYEE2 solver'
          if (order==LINEAR) then
            call MAXWELEYEE2(exy,bxy,cu(1,1,1),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax)
          else if (order==CUBIC) then
            call MAXWELEYEE2(exy,bxy,cu(1,3,3),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax)
          else
            call MAXWELEYEE2(exy,bxy,cu(1,2,2),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax)
          endif
          
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine imaxweleyee2
!
         subroutine imaxwelb2(exy,bxy,ci,dt,wm,tfield,nx,ny,nxe,nye,nxhd,&
      &nyhd,inorder)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny, nxe, nye, nxhd, nyhd
         integer, optional :: inorder
         real :: ci, dt, wm, tfield
         complex, dimension(:,:,:), pointer :: exy, bxy
! local data
         integer :: nxvh, nyv, ltime, order
         real :: tf
         nxvh = nxe/2; nyv = nye
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
          print *, 'using MAXWELB2 solver'
          if (order==LINEAR) then
            call MAXWELB2(exy,bxy,ci,dt,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
          else if (order==CUBIC) then
            call MAXWELB2(exy,bxy,ci,dt,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
          else
            call MAXWELB2(exy,bxy,ci,dt,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
          endif
          
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine imaxwelb2
!
         subroutine imaxwele2(exy,bxy,cu,ffc,ci,dt,wf,wm,tfield,nx,ny,ino&
     &rder,parmax)
! calculates maxwell's equation for periodic 2d vector field
         implicit none
         integer :: nx, ny, parmax
         integer, optional :: inorder
         real :: ci, dt, wf, wm, tfield
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, ltime, order
         real :: tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
          print *, 'using MAXWELE2 solver'
          if (order==LINEAR) then
            call MAXWELE2(exy,bxy,cu(1,1,1),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax)
          else if (order==CUBIC) then
            call MAXWELE2(exy,bxy,cu(1,3,3),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax)
          else
            call MAXWELE2(exy,bxy,cu(1,2,2),ffc,ci,dt,wf,wm,nx,ny,nxvh,ny&
            &v,nxhd,nyhd,parmax)
          endif
          
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine imaxwele2
!*****
!                call emfield(fxyze,exyz,ffc,isign,nx,ny,inorder)
         subroutine iemfield2(fxy,exy,ffc,isign,nx,ny,inorder)
! combines and smooths periodic 2d vector fields
         implicit none
         integer :: isign, nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call EMFIELD2(fxy(1,1,1),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else if (order==CUBIC) then
            call EMFIELD2(fxy(1,3,3),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else
            call EMFIELD2(fxy(1,2,2),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         endif
         end subroutine iemfield2
         
!        add laser
         subroutine ilsremfield2(fxy,exy,ffc,isign,nx,ny,inorder)
! combines and smooths periodic 2d vector fields
         implicit none
         integer :: isign, nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call LSREMFIELD2(fxy(1,1,1),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else if (order==CUBIC) then
            call LSREMFIELD2(fxy(1,3,3),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else
            call LSREMFIELD2(fxy(1,2,2),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         endif
         end subroutine ilsremfield2
         
!        exyz to lsrfxyze (fourier to fourier)
         subroutine ilsremfield4(fxy,exy,ffc,isign,nx,ny,inorder)
! combines and smooths periodic 2d vector fields
         implicit none
         integer :: isign, nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call LSREMFIELD4(fxy(1,1,1),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else if (order==CUBIC) then
            call LSREMFIELD4(fxy(1,3,3),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else
            call LSREMFIELD4(fxy(1,2,2),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         endif
         end subroutine ilsremfield4
         
         subroutine isavcuz2(cu,cu0,inorder)
! save current at ntime == 0
         implicit none
         integer :: inorder
         real, dimension(:,:,:), pointer :: cu
         real, dimension(:), pointer :: cu0
! local data
         integer :: nxvh,nyv
         nxvh = size(cu,2); nyv = size(cu,3)
         if (inorder == LINEAR) then
            call SAVCUZ2(cu(1,1,1),cu0,nxvh,nyv)
         elseif (inorder == CUBIC) then
            call SAVCUZ2(cu(1,3,3),cu0,nxvh,nyv)
         else
            call SAVCUZ2(cu(1,2,2),cu0,nxvh,nyv)
            !print *, 'quad save current'
         endif
         end subroutine isavcuz2
         
         subroutine imaxz2(exy,cu,cu0,ffc,dt,wf,nx,ny,inorder)
!        solve for k = 0 mode
         implicit none
         complex, dimension(:,:,:), pointer :: exy
         real, dimension(:,:,:), pointer :: cu
         real, dimension(:), pointer :: cu0
         complex, dimension(:,:), pointer :: ffc
         real :: dt, wf
         integer :: nx, ny, inorder
         
         integer :: nxeh, nye, nxh, nyh
         nxeh = size(cu,2); nye = size(cu,3)
         nxh = size(ffc,1); nyh = size(ffc,2)
         if (inorder == LINEAR) then
            call MAXZ2(exy,cu(1,1,1),cu0,ffc,dt,wf,nx,ny,nxeh,nye,nxh,nyh)
         elseif (inorder == CUBIC) then
            call MAXZ2(exy,cu(1,3,3),cu0,ffc,dt,wf,nx,ny,nxeh,nye,nxh,nyh)
         else
            call MAXZ2(exy,cu(1,2,2),cu0,ffc,dt,wf,nx,ny,nxeh,nye,nxh,nyh)
         endif
         
         end subroutine imaxz2
         
!        set zero boundary
         subroutine ilsremfield3(fxy,exy,ffc,isign,nx,ny,inorder)
! combines and smooths periodic 2d vector fields
         implicit none
         integer :: isign, nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(fxy,2)/2; nyv = size(fxy,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call LSREMFIELD3(fxy(1,1,1),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else if (order==CUBIC) then
            call LSREMFIELD3(fxy(1,3,3),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         else
            call LSREMFIELD3(fxy(1,2,2),exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,n&
     &yhd)
         endif
         end subroutine ilsremfield3
!
         subroutine iemfieldr2(fxy,exy,ffd,isign,nx,ny,order)
! combines and smooths real 2d vector fields
         implicit none
         integer :: isign, nx, ny
         integer, optional :: order
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffd
! local data
         integer :: nxe, nye, nxv, inorder
         nxe = size(fxy,2); nye = size(fxy,3); nxv = size(ffd,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call EMFIELDR2(fxy(1,1,1),exy,ffd,isign,nx,ny,nxv,nxe,nye)
         else if (inorder==CUBIC) then
            call EMFIELDR2(fxy(1,3,3),exy,ffd,isign,nx,ny,nxv,nxe,nye)
         else
            call EMFIELDR2(fxy(1,2,2),exy,ffd,isign,nx,ny,nxv,nxe,nye)
         endif
         end subroutine iemfieldr2
!
         subroutine iemfieldc2(fxy,exy,ffb,isign,nx,ny,order)
! combines and smooths real 2d vector fields
         implicit none
         integer :: isign, nx, ny
         integer, optional :: order
         real, dimension(:,:,:), pointer :: fxy
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffb
! local data
         integer :: nxe, nyeh, nxv, inorder
         nxe = size(fxy,2); nyeh = size(fxy,3)/2; nxv = size(ffb,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call EMFIELDC2(fxy(1,1,1),exy,ffb,isign,nx,ny,nxv,nxe,nyeh)
         else if (inorder==CUBIC) then
            call EMFIELDC2(fxy(1,3,3),exy,ffb,isign,nx,ny,nxv,nxe,nyeh)
         else
            call EMFIELDC2(fxy(1,2,2),exy,ffb,isign,nx,ny,nxv,nxe,nyeh)
         endif
         end subroutine iemfieldc2
!
         subroutine iavpot23(bxy,axy,nx,ny,inorder)
! calculates periodic 2d vector potential from magnetic field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         complex, dimension(:,:,:), pointer :: bxy
         real, dimension(:,:,:), pointer :: axy
! local data
         integer :: nxvh, nyv, order
         nxvh = size(bxy,2); nyv = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call AVPOT23(bxy,axy(1,1,1),nx,ny,nxvh,nyv)
         else if (order==CUBIC) then
            call AVPOT23(bxy,axy(1,3,3),nx,ny,nxvh,nyv)
         else
            call AVPOT23(bxy,axy(1,2,2),nx,ny,nxvh,nyv)
         endif
         end subroutine iavpot23
!
         subroutine iavrpot23(axy,bxy,ffc,ci,nx,ny,inorder)
! calculates periodic 2d radiative vector potential from magnetic field
! and current
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci
         complex, dimension(:,:,:), pointer :: bxy
         real, dimension(:,:,:), pointer :: axy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(axy,2)/2; nyv = size(axy,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call AVRPOT23(axy(1,1,1),bxy,ffc,ci,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
         else if (order==CUBIC) then
            call AVRPOT23(axy(1,3,3),bxy,ffc,ci,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
         else
            call AVRPOT23(axy(1,2,2),bxy,ffc,ci,nx,ny,nxvh,nyv,nxhd,nyhd&
     &)
         endif
         end subroutine iavrpot23
!
         subroutine igtmodes2(pot,pott,nx,ny,modesx,modesy,order)
! extracts lowest order modes from periodic 2d scalar field
         implicit none
         integer :: nx, ny, modesx, modesy
         integer, optional :: order
         real, dimension(:,:), pointer :: pot
         complex, dimension(:,:), pointer :: pott
! local data
         integer :: nxe, nye, it, nt2, modesxd, modesyd, inorder
         nxe = size(pot,1); nye = size(pot,2)
         it = 1; nt2 = 2
         modesxd = size(pott,1); modesyd = size(pott,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call GTMODES2(pot(1,1),pott,nx,ny,it,modesx,modesy,nxe,nye,n&
     &t2,modesxd,modesyd)
         else if (inorder==CUBIC) then
            call GTMODES2(pot(3,3),pott,nx,ny,it,modesx,modesy,nxe,nye,n&
     &t2,modesxd,modesyd)
         else
            call GTMODES2(pot(2,2),pott,nx,ny,it,modesx,modesy,nxe,nye,n&
     &t2,modesxd,modesyd)
         endif
         end subroutine igtmodes2
!
         subroutine iptmodes2(pot,pott,nx,ny,modesx,modesy,order)
! extracts lowest order modes to periodic 2d scalar field
         implicit none
         integer :: nx, ny, modesx, modesy
         integer, optional :: order
         real, dimension(:,:), pointer :: pot
         complex, dimension(:,:), pointer :: pott
! local data
         integer :: nxe, nye, it, nt2, modesxd, modesyd, inorder
         nxe = size(pot,1); nye = size(pot,2)
         it = 1; nt2 = 2
         modesxd = size(pott,1); modesyd = size(pott,2)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call PTMODES2(pot(1,1),pott,nx,ny,it,modesx,modesy,nxe,nye,n&
     &t2,modesxd,modesyd) 
         else if (inorder==CUBIC) then
            call PTMODES2(pot(3,3),pott,nx,ny,it,modesx,modesy,nxe,nye,n&
     &t2,modesxd,modesyd) 
         else
            call PTMODES2(pot(2,2),pott,nx,ny,it,modesx,modesy,nxe,nye,n&
     &t2,modesxd,modesyd) 
         endif
         end subroutine iptmodes2
!
         subroutine igtvmodes2(vpot,vpott,nx,ny,modesx,modesy,order)
! extracts lowest order modes from periodic 2d vector field
         implicit none
         integer :: nx, ny, modesx, modesy
         integer, optional :: order
         real, dimension(:,:,:), pointer :: vpot
         complex, dimension(:,:,:), pointer :: vpott
! local data
         integer :: ndim, nxvh, nyv, it, nt, modesxd, modesyd, inorder
         ndim = size(vpot,1); nxvh = size(vpot,2)/2; nyv = size(vpot,3)
         it = 1; nt = 1
         modesxd = size(vpott,2); modesyd = size(vpott,3)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call GTVMODES2(vpot(1,1,1),vpott,nx,ny,it,modesx,modesy,ndim&
     &,nxvh,nyv,nt,modesxd,modesyd)
         else if (inorder==CUBIC) then
            call GTVMODES2(vpot(1,3,3),vpott,nx,ny,it,modesx,modesy,ndim&
     &,nxvh,nyv,nt,modesxd,modesyd)
         else
            call GTVMODES2(vpot(1,2,2),vpott,nx,ny,it,modesx,modesy,ndim&
     &,nxvh,nyv,nt,modesxd,modesyd)
         endif
         end subroutine igtvmodes2
!
         subroutine iptvmodes2(vpot,vpott,nx,ny,modesx,modesy,order)
! extracts lowest order modes to periodic 2d vector field
         implicit none
         integer :: nx, ny, modesx, modesy
         integer, optional :: order
         real, dimension(:,:,:), pointer :: vpot
         complex, dimension(:,:,:), pointer :: vpott
! local data
         integer :: ndim, nxvh, nyv, it, nt, modesxd, modesyd, inorder
         ndim = size(vpot,1); nxvh = size(vpot,2)/2; nyv = size(vpot,3)
         it = 1; nt = 1
         modesxd = size(vpott,2); modesyd = size(vpott,3)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (inorder==LINEAR) then
            call PTVMODES2(vpot(1,1,1),vpott,nx,ny,it,modesx,modesy,ndim&
     &,nxvh,nyv,nt,modesxd,modesyd)
         else if (inorder==CUBIC) then
            call PTVMODES2(vpot(1,3,3),vpott,nx,ny,it,modesx,modesy,ndim&
     &,nxvh,nyv,nt,modesxd,modesyd) 
         else
            call PTVMODES2(vpot(1,2,2),vpott,nx,ny,it,modesx,modesy,ndim&
     &,nxvh,nyv,nt,modesxd,modesyd) 
         endif
         end subroutine iptvmodes2
!
         subroutine ipoynt2(q,exy,bxy,ffc,sx,sy,sz,nx,ny,inorder)
! calculates the momentum in the electromagnetic field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: sx, sy, sz
         complex, dimension(:,:,:), pointer :: exy, bxy
         real, dimension(:,:), pointer :: q
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(q,1)/2; nyv = size(q,2)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call POYNT2(q(1,1),exy,bxy,ffc,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,&
     &nyhd)
         else if (order==CUBIC) then
            call POYNT2(q(3,3),exy,bxy,ffc,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,&
     &nyhd)
         else
            call POYNT2(q(2,2),exy,bxy,ffc,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,&
     &nyhd)
         endif
         end subroutine ipoynt2
!
         subroutine idpoynt2(q,cu,ffc,ci,sx,sy,sz,nx,ny,inorder)
! calculates the momentum in the darwin field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, sx, sy, sz
         real, dimension(:,:), pointer :: q
         real, dimension(:,:,:), pointer :: cu
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: nxvh, nyv, nxhd, nyhd, order
         nxvh = size(q,1)/2; nyv = size(q,2)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call DPOYNT22(q(1,1),cu(1,1,1),ffc,ci,sx,sy,nx,ny,nxvh,ny&
     &v,nxhd,nyhd)
            else if (order==CUBIC) then
               call DPOYNT22(q(3,3),cu(1,3,3),ffc,ci,sx,sy,nx,ny,nxvh,ny&
     &v,nxhd,nyhd)
            else
               call DPOYNT22(q(2,2),cu(1,2,2),ffc,ci,sx,sy,nx,ny,nxvh,ny&
     &v,nxhd,nyhd)
            endif
            sz = 0.0
         case (3)
            if (order==LINEAR) then
               call DPOYNT2(q(1,1),cu(1,1,1),ffc,ci,sx,sy,sz,nx,ny,nxvh,&
     &nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call DPOYNT2(q(3,3),cu(1,3,3),ffc,ci,sx,sy,sz,nx,ny,nxvh,&
     &nyv,nxhd,nyhd)
            else
               call DPOYNT2(q(2,2),cu(1,2,2),ffc,ci,sx,sy,sz,nx,ny,nxvh,&
     &nyv,nxhd,nyhd)
            endif
         end select
         end subroutine idpoynt2
!
         subroutine isfguard2(cus,cu,q2m0,nx,ny,inorder)
! initialize periodic 2d vector field with scaled field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: q2m0
         real, dimension(:,:,:), pointer :: cu, cus
! local data
         integer :: nxe, nye, order
         nxe = size(cus,2); nye = size(cus,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(cus,1))
         case (2)
            if (order==LINEAR) then
               call SCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
            else
               call SCFGUARD22(cus,cu,q2m0,nx,ny,nxe,nye)
            endif
         case (3)
            if (order==LINEAR) then
               call SCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
            else
               call SCFGUARD2(cus,cu,q2m0,nx,ny,nxe,nye)
            endif
         end select
         end subroutine isfguard2
!
         subroutine ismcguard2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,inorder)
! initialize periodic 2d tensor field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: x2y2m0, xym0, zxm0, zym0
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: nxe, nye, order
         nxe = size(amu,2); nye = size(amu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(amu,1))
         case (2)
            if (order==LINEAR) then
               call SMCGUARD22L(amu,x2y2m0,xym0,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call SMCGUARD22C(amu,x2y2m0,xym0,nx,ny,nxe,nye)
            else
               call SMCGUARD22(amu,x2y2m0,xym0,nx,ny,nxe,nye)
            endif
         case (4)
            if (order==LINEAR) then
               call SMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call SMCGUARD2C(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
            else
               call SMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
            endif
         end select
         end subroutine ismcguard2
!
         subroutine iamcguard2(amu,nx,ny,inorder)
! add guard cells for periodic 2d tensor data
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: amu
! local data
         integer :: ndim, nxe, nye, order
         ndim = size(amu,1); nxe = size(amu,2); nye = size(amu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call AMCGUARD2L(amu,nx,ny,nxe,nye,ndim)
         else if (order==CUBIC) then
            call ACGUARD2C(amu,nx,ny,nxe,nye,ndim)
         else
            call AMCGUARD2(amu,nx,ny,nxe,nye,ndim)
         endif
         end subroutine iamcguard2
!
         subroutine idcuperp23(dcu,amu,nx,ny,inorder)
! calculates the transverse part of periodic 2d vector field
! from momentum flux tensor
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: dcu, amu
! local data
         integer :: nxvh, nyv, order
         nxvh = size(dcu,2)/2; nyv = size(dcu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call DCUPERP22(dcu(1,1,1),amu(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call DCUPERP22(dcu(1,3,3),amu(1,3,3),nx,ny,nxvh,nyv)
            else
               call DCUPERP22(dcu(1,2,2),amu(1,2,2),nx,ny,nxvh,nyv)
            endif
         case (3)
            if (order==LINEAR) then
               call DCUPERP23(dcu(1,1,1),amu(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call DCUPERP23(dcu(1,3,3),amu(1,3,3),nx,ny,nxvh,nyv)
            else
               call DCUPERP23(dcu(1,2,2),amu(1,2,2),nx,ny,nxvh,nyv)
            endif
         end select
         end subroutine idcuperp23
!
         subroutine iadcuperp23(dcu,amu,nx,ny,inorder)
! calculates the transverse part of periodic 2d vector field
! from acceleration vector and momentum flux tensor
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: dcu, amu
! local data
         integer :: nxvh, nyv, order
         nxvh = size(dcu,2)/2; nyv = size(dcu,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call ADCUPERP22(dcu(1,1,1),amu(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call ADCUPERP22(dcu(1,3,3),amu(1,3,3),nx,ny,nxvh,nyv)
            else
               call ADCUPERP22(dcu(1,2,2),amu(1,2,2),nx,ny,nxvh,nyv)
            endif
         case (3)
            if (order==LINEAR) then
               call ADCUPERP23(dcu(1,1,1),amu(1,1,1),nx,ny,nxvh,nyv)
            else if (order==CUBIC) then
               call ADCUPERP23(dcu(1,3,3),amu(1,3,3),nx,ny,nxvh,nyv)
            else
               call ADCUPERP23(dcu(1,2,2),amu(1,2,2),nx,ny,nxvh,nyv)
            endif
         end select
         end subroutine iadcuperp23
!
         subroutine iepois23init(ffe,ax,ay,affp,wp0,ci,nx,ny)
! initialize 2d periodic transverse electric field solver
         implicit none
         integer :: nx, ny
         real :: ax, ay, affp, wp0, ci
         complex, dimension(:,:), pointer :: ffe
! local data
         integer :: isign = 0, nxvh = 1, nyv = 1, nxhd, nyhd
         real :: wf, dcu, exy
         nxhd = size(ffe,1); nyhd = size(ffe,2)
         call EPOIS23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,nxvh,&
     &nyv,nxhd,nyhd)
         end subroutine iepois23init
!
         subroutine iepois23(dcu,exy,ffe,ci,wf,tfield,nx,ny,inorder)
! calculates transverse electric field for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, wf, tfield
         real, dimension(:,:,:), pointer :: dcu, exy
         complex, dimension(:,:), pointer :: ffe
! local data
         integer :: isign = -1, nxvh, nyv, nxhd, nyhd, ltime, order
         real :: ax, ay, affp, wp0, tf
         nxvh = size(dcu,2)/2; nyv = size(dcu,3)
         nxhd = size(ffe,1); nyhd = size(ffe,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call EPOIS22(dcu(1,1,1),exy(1,1,1),isign,ffe,ax,ay,affp,w&
     &p0,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call EPOIS22(dcu(1,3,3),exy(1,3,3),isign,ffe,ax,ay,affp,w&
     &p0,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call EPOIS22(dcu(1,2,2),exy(1,2,2),isign,ffe,ax,ay,affp,w&
     &p0,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call EPOIS23(dcu(1,1,1),exy(1,1,1),isign,ffe,ax,ay,affp,w&
     &p0,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call EPOIS23(dcu(1,3,3),exy(1,3,3),isign,ffe,ax,ay,affp,w&
     &p0,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call EPOIS23(dcu(1,2,2),exy(1,2,2),isign,ffe,ax,ay,affp,w&
     &p0,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         end select
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine iepois23
!
         subroutine iiepois23(dcu,exy,ffe,ci,wf,nx,ny,inorder,yee)
! calculates transverse electric field for periodic 2d vector field
! without smoothing
         implicit none
         integer :: nx, ny, yee
         integer, optional :: inorder
         real :: ci, wf
         real, dimension(:,:,:), pointer :: dcu
         complex, dimension(:,:,:), pointer :: exy
         complex, dimension(:,:), pointer :: ffe
! local data
         integer :: isign = 1, nxvh, nyv, nxhd, nyhd, order
         real :: ax, ay, affp, wp0
         nxvh = size(dcu,2)/2; nyv = size(dcu,3)
         nxhd = size(ffe,1); nyhd = size(ffe,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(dcu,1))
         case (2)
            if (order==LINEAR) then
               call EPOIS22(dcu(1,1,1),exy,isign,ffe,ax,ay,affp,wp0,ci,w&
     &f,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call EPOIS22(dcu(1,3,3),exy,isign,ffe,ax,ay,affp,wp0,ci,w&
     &f,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call EPOIS22(dcu(1,2,2),exy,isign,ffe,ax,ay,affp,wp0,ci,w&
     &f,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call EPOIS23(dcu(1,1,1),exy,isign,ffe,ax,ay,affp,wp0,ci,w&
     &f,nx,ny,nxvh,nyv,nxhd,nyhd,yee)
            else if (order==CUBIC) then
               call EPOIS23(dcu(1,3,3),exy,isign,ffe,ax,ay,affp,wp0,ci,w&
     &f,nx,ny,nxvh,nyv,nxhd,nyhd,yee)
            else
               call EPOIS23(dcu(1,2,2),exy,isign,ffe,ax,ay,affp,wp0,ci,w&
     &f,nx,ny,nxvh,nyv,nxhd,nyhd,yee)
            endif
         end select
         end subroutine iiepois23
!
         subroutine idapois23(cu,daxy,ffc,ci,wm,tfield,nx,ny,inorder)
! calculates derivatives of smoothed vector potential
! for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, wm, tfield
         real, dimension(:,:,:), pointer :: cu, daxy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = -1, nxvh, nyv, nxhd, nyhd, ltime, order
         real :: ax, ay, affp, axy0, tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call APOIS22(cu(1,1,1),daxy(1,1,1),axy0,isign,ffc,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call APOIS22(cu(1,3,3),daxy(1,3,3),axy0,isign,ffc,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call APOIS22(cu(1,2,2),daxy(1,2,2),axy0,isign,ffc,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call APOIS23(cu(1,1,1),daxy(1,1,1),axy0,isign,ffc,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call APOIS23(cu(1,3,3),daxy(1,3,3),axy0,isign,ffc,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call APOIS23(cu(1,2,2),daxy(1,2,2),axy0,isign,ffc,ax,ay,a&
     &ffp,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         end select
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine idapois23
!
         subroutine isapois23(cu,axy,ffc,ci,tfield,nx,ny,inorder)
! calculates smoothed vector potential for periodic 2d vector field
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: ci, wf, tfield
         real, dimension(:,:,:), pointer :: cu, axy
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: isign = 1, nxvh, nyv, nxhd, nyhd, ltime, order
         real :: ax, ay, affp, daxy0, tf
         nxvh = size(cu,2)/2; nyv = size(cu,3)
         nxhd = size(ffc,1); nyhd = size(ffc,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(cu,1))
         case (2)
            if (order==LINEAR) then
               call APOIS22(cu(1,1,1),daxy0,axy(1,1,1),isign,ffc,ax,ay,a&
     &ffp,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call APOIS22(cu(1,3,3),daxy0,axy(1,3,3),isign,ffc,ax,ay,a&
     &ffp,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call APOIS22(cu(1,2,2),daxy0,axy(1,2,2),isign,ffc,ax,ay,a&
     &ffp,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         case (3)
            if (order==LINEAR) then
               call APOIS23(cu(1,1,1),daxy0,axy(1,1,1),isign,ffc,ax,ay,a&
     &ffp,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else if (order==CUBIC) then
               call APOIS23(cu(1,3,3),daxy0,axy(1,3,3),isign,ffc,ax,ay,a&
     &ffp,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            else
               call APOIS23(cu(1,2,2),daxy0,axy(1,2,2),isign,ffc,ax,ay,a&
     &ffp,ci,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
            endif
         end select
         call wtimer(tf,ltime)
         tfield = tfield + tf
         end subroutine isapois23
!
         subroutine iaddqei2(qe,qi,nx,ny,inorder)
! adds electron and ion densities
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real, dimension(:,:), pointer :: qe, qi
! local data
         integer :: nxe, nye, order
         nxe = size(qe,1); nye = size(qe,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call ADDQEI2(qe(1,1),qi(1,1),nx,ny,nxe,nye)
         else if (order==CUBIC) then
            call ADDQEI2(qe(3,3),qi(3,3),nx,ny,nxe,nye)
         else
            call ADDQEI2(qe(2,2),qi(2,2),nx,ny,nxe,nye)
         endif
         end subroutine iaddqei2
!
         subroutine iaddqei2x(qe,qi,qbme,qbmi,wpmax,wpmin,nx,ny,inorder)
! adds electron and ion densities, and calculates maximum and minimum
! plasma frequency
         implicit none
         integer :: nx, ny
         integer, optional :: inorder
         real :: qbme, qbmi, wpmax, wpmin
         real, dimension(:,:), pointer :: qe, qi
! local data
         integer :: nxe, nye, order
         nxe = size(qe,1); nye = size(qe,2)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call ADDQEI2X(qe(1,1),qi(1,1),qbme,qbmi,wpmax,wpmin,nx,ny,nx&
     &e,nye)
         else if (order==CUBIC) then
            call ADDQEI2X(qe(3,3),qi(3,3),qbme,qbmi,wpmax,wpmin,nx,ny,nx&
     &e,nye)
         else
            call ADDQEI2X(qe(2,2),qi(2,2),qbme,qbmi,wpmax,wpmin,nx,ny,nx&
     &e,nye)
         endif
         end subroutine iaddqei2x
!
         subroutine ibaddext2(bxy,omx,omy,omz,nx,ny,inorder)
! adds constant to magnetic field
         implicit none
         integer :: nx, ny
         real :: omx, omy, omz
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: bxy
! local data
         integer :: nxe, nye, order
         nxe = size(bxy,2); nye = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call BADDEXT22(bxy(1,1,1),omz,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call BADDEXT22(bxy(1,3,3),omz,nx,ny,nxe,nye)
            else
               call BADDEXT22(bxy(1,2,2),omz,nx,ny,nxe,nye)
            endif
         case (3)
            if (order==LINEAR) then
               call BADDEXT2(bxy(1,1,1),omx,omy,omz,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call BADDEXT2(bxy(1,3,3),omx,omy,omz,nx,ny,nxe,nye)
            else
               call BADDEXT2(bxy(1,2,2),omx,omy,omz,nx,ny,nxe,nye)
            endif
         end select
         end subroutine ibaddext2
!
         subroutine iimoment2(qi,fxy,iunit,px,py,pz,dt,wx,wy,wz,nx,ny,in&
     &order)
! calculate ion momentum from integral of qi*fxy,
! and prints it, and adds it to total momentum, for 2 or 2-1/2d code
         implicit none
         integer :: nx, ny, iunit
         real :: px, py, pz, dt, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:), pointer :: qi
         real, dimension(:,:,:), pointer :: fxy
! local data
         integer :: nxe, nye, order
  995    format (' ion momentum = ',3e14.7)
         nxe = size(fxy,2); nye = size(fxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
! calculate and print ion momentum
         select case(size(fxy,1))
         case (2)
            if (order==LINEAR) then
               call IMOMENT22(qi(1,1),fxy(1,1,1),px,py,dt,nx,ny,nxe,nye)
            else if (order==CUBIC) then
               call IMOMENT22(qi(3,3),fxy(1,3,3),px,py,dt,nx,ny,nxe,nye)
            else
               call IMOMENT22(qi(2,2),fxy(1,2,2),px,py,dt,nx,ny,nxe,nye)
            endif
            pz = 0.0
         case (3)
            if (order==LINEAR) then
               call IMOMENT2(qi(1,1),fxy(1,1,1),px,py,pz,dt,nx,ny,nxe,ny&
     &e)
            else if (order==CUBIC) then
               call IMOMENT2(qi(3,3),fxy(1,3,3),px,py,pz,dt,nx,ny,nxe,ny&
     &e)
            else
               call IMOMENT2(qi(2,2),fxy(1,2,2),px,py,pz,dt,nx,ny,nxe,ny&
     &e)
            endif
         end select
         write (iunit,995) px, py, pz
! add to total momentum
         wx = wx + px
         wy = wy + py
         wz = wz + pz
         end subroutine iimoment2
!
         subroutine ivccopy2(f,g,nx,ny,inorder)
! copies complex vector array elements from f to g
         implicit none
         integer nx, ny
         integer, optional :: inorder
         real, dimension(:,:,:) :: f
         complex, dimension(:,:,:) :: g
! local data
         integer :: ndim, nxv, nyv, order
         ndim = size(g,1); nxv = size(g,2); nyv = size(g,3)
         if ((nx > nxv) .or. (ny > nyv)) return
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call VCCOPY2(f(1,1,1),g,nx,ny,ndim,nxv,nyv)
         else if (order==CUBIC) then
            call VCCOPY2(f(1,3,3),g,nx,ny,ndim,nxv,nyv)
         else
            call VCCOPY2(f(1,2,2),g,nx,ny,ndim,nxv,nyv)
         endif
         end subroutine ivccopy2
!
         subroutine iaddvrfield2(a,b,c)
! calculate a = b + c for real vector fields
         implicit none
         real, dimension(:,:,:), pointer :: a, b, c
! local data
         integer :: ndim, nxe, nye
         ndim = size(a,1); nxe = size(a,2); nye = size(a,3)
         call ADDVRFIELD2(a,b,c,ndim,nxe,nye)
         end subroutine iaddvrfield2
!
      end module field2d
