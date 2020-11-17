c-----------------------------------------------------------------------
c 2d PIC library for solving field equations
c field2lib.f contains procedures to manage guard cells and solve
c             fields equations in fourier space:
c CGUARD2 copy guard cells for 2 component vector array, quadratic
c         interpolation.
c BGUARD2 copy guard cells for 3 component vector array, quadratic
c         interpolation.
c DGUARD2 copy guard cells for scalar array, quadratic interpolation.
c SCGUARD2 initialize field for 3 component vector array, quadratic
c          interpolation.
c SCGUARD22 initialize field for 2 component vector array, quadratic
c           interpolation.
c SGUARD2 initialize field for scalar array, quadratic interpolation.
c ACGUARD2 add guard cells for 3 component vector array, quadratic
c          interpolation.
c ACGUARD22 add guard cells for 2 component vector array, quadratic
c           interpolation.
c AGUARD2 add guard cells for scalar array, quadratic interpolation. 
c CGUARD2L copy guard cells for 2 component vector array, linear
c          interpolation.
c BGUARD2L copy guard cells for 3 component vector array, linear
c          interpolation.
c DGUARD2L copy guard cells for scalar array, linear interpolation.
c SCGUARD2L initialize field for 3 component vector array, linear
c           interpolation.
c SCGUARD22L initialize field for 2 component vector array, linear
c            interpolation.
c SGUARD2L initialize field for scalar array, linear interpolation.
c ACGUARD2L add guard cells for 3 component vector array, linear
c           interpolation.
c ACGUARD22L add guard cells for 2 component vector array, linear
c            interpolation.
c AGUARD2L add guard cells for scalar array, linear interpolation.
c CGUARD2C copy guard cells for n component vector array, cpois23232323ubic
c          interpolation.
c DGUARD2C copy guard cells for scalar array, cubic interpolation.
c SCGUARD2C initialize field for 3 component vector array, cubic
c           interpolation.
c SCGUARD22C initialize field for 2 component vector array, cubic
c            interpolation.
c SGUARD2C initialize field for scalar array, cubic interpolation.
c ACGUARD2C add guard cells for n component vector array, cubic
c           interpolation.
c AGUARD2C add guard cells for scalar array, cubic interpolation.
c POISP2 solve 2d poisson equation for electric force, potential, or
c        smoothing.
c POIS22 solve 2d poisson equation for electric force.
c POIS23 solve 2-1/2d poisson equation for electric force.
c POISP2T solve 2d poisson equation for electric force, potential, or
c         smoothing for transposed data.
c POIS23T solve 2-1/2d poisson equation for electric force force
c         for transposed data.
c DIVF2 calculates 2d divergence of n component vector in fourier space.
c GRADF2 calculates 2d gradient of scalar field in fourier space.
c CURLF2 calculates 2d divergence of 3 component vector in fourier space
c CURLF22 calculates 2d divergence of 2 component vector in fourier
c         space.
c LAPLACE23 solves 2d vector laplacian of 3 component vector in fourier
c           space.
c LAPLACE22 solves 2d vector laplacian of 2 component vector in fourier
c           space.
c CUPERP2 calculates 2d transverse current of 3 component vector in
c         fourier space.
c CUPERP22 calculates 2d transverse current of 2 component vector in
c         fourier space.
c CUPERP2T calculates 2d transverse current of 3 component vector in
c          fourier space for transposed data.
c BPOIS23 solve 2-1/2d vector poisson equation for magnetic force,
c         vector potential, or smoothing.
c BPOIS22 solve 2d vector poisson equation for magnetic force, vector
c         potential, or smoothing.
c IBPOIS23 solve 2-1/2d vector poisson equation for magnetic field.
c IBPOIS23T solve 2-1/2d vector poisson equation for magnetic field
c           for transposed data.
c MAXWEL2 solve 2d maxwell equation for electric and magnetic fields.
c MAXWEL2T solve 2d maxwell equation for electric and magnetic fields
c          for transposed data.
c EMFIELD2 combines and smooths 2d periodic electric or magnetic forces.
c EMFIELDR2 combines and smooths 2d real electric or magnetic forces for
c           sine-cosine transforms.
c EMFIELDC2 combines and smooths 2d real electric or magnetic forces for
c           sine-cosine/periodic transforms.
c EMFIELD2T combines and smooths 2d periodic electric or magnetic forces
c           for transposed data.
c AVPOT23 calculate 2-1/2d vector potential from magnetic field.
c AVRPOT23 calculate 2-1/2d radiative part of the vector potential.
c GTMODES2 extracts selected 2d fourier components from potential array.
c PTMODES2 places selected 2d fourier components into potential array.
c GTVMODES2 extracts selected 2d fourier components from vector
c           potential array.
c PTVMODES2 places selected 2d fourier components into vector potential
c           array.
c POYNT2 calculate poynting electromagnetic flux.
c DPOYNT2 calculate electromagnetic flux in 2-1/2d Darwin field.
c DPOYNT22 calculate electromagnetic flux in 2d Darwin field.
c SCFGUARD2 initialize 3 component field with scaled vector array,
c           quadratic interpolation.
c SCFGUARD22 initialize 2 component field with scaled vector array,
c            quadratic interpolation.
c SMCGUARD2 initialize field for 4 component tensor array, quadratic
c           interpolation.
c SMCGUARD22 initialize field for 2 component tensor array, quadratic
c            interpolation.
c AMCGUARD2 add guard cells for 4 component tensor array, quadratic
c           interpolation.
c SCFGUARD2L initialize 3 component field with scaled vector array,
c            linear interpolation.
c SCFGUARD22L initialize 2 component field with scaled vector array,
c             linear interpolation.
c SMCGUARD2L initialize field for 4 component tensor array, linear
c            interpolation.
c SMCGUARD22L initialize field for 2 component tensor array, linear
c             interpolation.
c AMCGUARD2L add guard cells for 4 component tensor array, linear
c            interpolation.
c SMCGUARD2C initialize field for 4 component tensor array, cubic
c            interpolation.
c SMCGUARD22C initialize field for 2 component tensor array, cubic
c             interpolation.
c DCUPERP23 calculate 2-1/d transverse derivative of current density
c           from momentum flux.
c DCUPERP22 calculate 2d transverse derivative of current density from
c           momentum flux.
c ADCUPERP23 calculate 2-1/2d transverse derivative of current density
c            from momentum flux and acceleration density.
c ADCUPERP22 calculate 2d transverse derivative of current density
c            from momentum flux and acceleration density.
c EPOIS23 solve vector poisson equation for 2-1/2d transverse electric
c         field or force.
c EPOIS22 solve vector poisson equation for 2d transverse electric field
c         or force.
c APOIS23 solves 2-1/2d poisson equation for smoothed derivative of
c         vector potential or smoothed vector potential.
c APOIS22 solves 2d poisson equation for smoothed derivative of
c         vector potential or smoothed vector potential.
c ADDQEI2 adds electron and ion densities.
c ADDQEI2X adds electron and ion densities and calculates maximum and
c          minimum plasma frequency.
c BADDEXT2 adds constant to magnetic field in real space for 2-1/2d code
c BADDEXT22 adds constant to magnetic field in real space for 2d code.
c IMOMENT2 calculates ion momentum for 2-1/2d code from qi*fxy.
c IMOMENT22 calculates ion momentum for 2d code from qi*fxy.
c POIS23GL solves 2-1/2d poisson's equation for electric force for
c          gridless code.
c ADDVRFIELD2 calculates a = b + c for real vector fields
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: march 1, 2012
!*****
!AGUARD2C fixed bug for only one cell wide in y
!AVCGUARD2, 2L, 2C adds guard cells for quadratic, linear and cubic
!                  particle shapes
!FLDSHIFT2 shifts fields between integer and yee mesh positions
!FLDMODINIT initializes fields by shifting them and
!MODECHOP2 chops off high k modes for debugging.
!IBPOIS23 added option for yee cell.
!MAXWEL2 added option for yee cell.
!PARFORCE2 claculates parallel component of electric force, for debugging.
!EMFIELD2 added block to add fxyz to exyz 
!****
c-----------------------------------------------------------------------
      subroutine CGUARD2(fxy,nx,ny,nxe,nye)
c replicate extended periodic field
c quadratic interpolation
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer i, j, k
      do 20 k = 1, ny
      do 10 i = 1, 2
      fxy(i,1,k+1) = fxy(i,nx+1,k+1)
      fxy(i,nx+2,k+1) = fxy(i,2,k+1)
      fxy(i,nx+3,k+1) = fxy(i,3,k+1)
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 2
      fxy(i,j+1,1) = fxy(i,j+1,ny+1)
      fxy(i,j+1,ny+2) = fxy(i,j+1,2)
      fxy(i,j+1,ny+3) = fxy(i,j+1,3)
   30 continue
   40 continue
      do 50 i = 1, 2
      fxy(i,1,1) = fxy(i,nx+1,ny+1)
      fxy(i,nx+2,1) = fxy(i,2,ny+1)
      fxy(i,nx+3,1) = fxy(i,3,ny+1)
      fxy(i,1,ny+2) = fxy(i,nx+1,2)
      fxy(i,nx+2,ny+2) = fxy(i,2,2)
      fxy(i,nx+3,ny+2) = fxy(i,3,2)
      fxy(i,1,ny+3) = fxy(i,nx+1,3)
      fxy(i,nx+2,ny+3) = fxy(i,2,3)
      fxy(i,nx+3,ny+3) = fxy(i,3,3)
   50 continue
      return
      end
      subroutine DGUARD2(q,nx,ny,nxe,nye)
c replicate extended periodic scalar field
c quadratic interpolation
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
      do 10 k = 1, ny
      q(1,k+1) = q(nx+1,k+1)
      q(nx+2,k+1) = q(2,k+1)
      q(nx+3,k+1) = q(3,k+1)
   10 continue
      do 20 j = 1, nx
      q(j+1,1) = q(j+1,ny+1)
      q(j+1,ny+2) = q(j+1,2)
      q(j+1,ny+3) = q(j+1,3)
   20 continue
      q(1,1) = q(nx+1,ny+1)
      q(nx+2,1) = q(2,ny+1)
      q(nx+3,1) = q(3,ny+1)
      q(1,ny+2) = q(nx+1,2)
      q(nx+2,ny+2) = q(2,2)
      q(nx+3,ny+2) = q(3,2)
      q(1,ny+3) = q(nx+1,3)
      q(nx+2,ny+3) = q(2,3)
      q(nx+3,ny+3) = q(3,3)
      return
      end
      subroutine BGUARD2(bxy,nx,ny,nxe,nye)
c replicate extended periodic vector field
c quadratic interpolation
      implicit none
      real bxy
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer i, j, k
      do 20 k = 1, ny
      do 10 i = 1, 3
      bxy(i,1,k+1) = bxy(i,nx+1,k+1)
      bxy(i,nx+2,k+1) = bxy(i,2,k+1)
      bxy(i,nx+3,k+1) = bxy(i,3,k+1)
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 3
      bxy(i,j+1,1) = bxy(i,j+1,ny+1)
      bxy(i,j+1,ny+2) = bxy(i,j+1,2)
      bxy(i,j+1,ny+3) = bxy(i,j+1,3)
   30 continue
   40 continue
      do 50 i = 1, 3
      bxy(i,1,1) = bxy(i,nx+1,ny+1)
      bxy(i,nx+2,1) = bxy(i,2,ny+1)
      bxy(i,nx+3,1) = bxy(i,3,ny+1)
      bxy(i,1,ny+2) = bxy(i,nx+1,2)
      bxy(i,nx+2,ny+2) = bxy(i,2,2)
      bxy(i,nx+3,ny+2) = bxy(i,3,2)
      bxy(i,1,ny+3) = bxy(i,nx+1,3)
      bxy(i,nx+2,ny+3) = bxy(i,2,3)
      bxy(i,nx+3,ny+3) = bxy(i,3,3)
   50 continue
      return
      end
      subroutine SCGUARD2(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
c initialize extended periodic field
c quadratic interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      cu(1,j+1,k+1) = xj0
      cu(2,j+1,k+1) = yj0
      cu(3,j+1,k+1) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+1) = 0.
      cu(i,nx+2,k+1) = 0.
      cu(i,nx+3,k+1) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 3
      cu(i,j+1,1) = 0.
      cu(i,j+1,ny+2) = 0.
      cu(i,j+1,ny+3) = 0.
   40 continue
   50 continue
      do 60 i = 1, 3
      cu(i,1,1) = 0.
      cu(i,nx+2,1) = 0.
      cu(i,nx+3,1) = 0.
      cu(i,1,ny+2) = 0.
      cu(i,nx+2,ny+2) = 0.
      cu(i,nx+3,ny+2) = 0.
      cu(i,1,ny+3) = 0.
      cu(i,nx+2,ny+3) = 0.
      cu(i,nx+3,ny+3) = 0.
   60 continue
      return
      end
      subroutine SCGUARD22(cu,xj0,yj0,nx,ny,nxe,nye)
c initialize extended periodic field
c quadratic interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      cu(1,j+1,k+1) = xj0
      cu(2,j+1,k+1) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+1) = 0.
      cu(i,nx+2,k+1) = 0.
      cu(i,nx+3,k+1) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 2
      cu(i,j+1,1) = 0.
      cu(i,j+1,ny+2) = 0.
      cu(i,j+1,ny+3) = 0.
   40 continue
   50 continue
      do 60 i = 1, 2
      cu(i,1,1) = 0.
      cu(i,nx+2,1) = 0.
      cu(i,nx+3,1) = 0.
      cu(i,1,ny+2) = 0.
      cu(i,nx+2,ny+2) = 0.
      cu(i,nx+3,ny+2) = 0.
      cu(i,1,ny+3) = 0.
      cu(i,nx+2,ny+3) = 0.
      cu(i,nx+3,ny+3) = 0.
   60 continue
      return
      end
      subroutine SGUARD2(q,qi0,nx,ny,nxe,nye)
c initialize extended periodic scalar field
c quadratic interpolation
      implicit none
      real q, qi0
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c initialize extended field, with zero in the edges
      do 20 k = 1, ny
      do 10 j = 1, nx
      q(j+1,k+1) = qi0
   10 continue
      q(1,k+1) = 0.
      q(nx+2,k+1) = 0.
      q(nx+3,k+1) = 0.
   20 continue
      do 30 j = 1, nx
      q(j+1,1) = 0.
      q(j+1,ny+2) = 0.
      q(j+1,ny+3) = 0.
   30 continue
      q(1,1) = 0.
      q(nx+2,1) = 0.
      q(nx+3,1) = 0.
      q(1,ny+2) = 0.
      q(nx+2,ny+2) = 0.
      q(nx+3,ny+2) = 0.
      q(1,ny+3) = 0.
      q(nx+2,ny+3) = 0.
      q(nx+3,ny+3) = 0.
      return
      end
      subroutine ACGUARD2(cu,nx,ny,nxe,nye)
c accumulate extended periodic 3 component vector field
c quadratic interpolation
c handles special case of ny=1 properly
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, 3
      cu(i,2,k+1) = cu(i,2,k+1) + cu(i,nx+2,k+1)
      cu(i,3,k+1) = cu(i,3,k+1) + cu(i,nx+3,k+1)
      cu(i,nx+1,k+1) = cu(i,nx+1,k+1) + cu(i,1,k+1)
      cu(i,nx+2,k+1) = 0.0
      cu(i,nx+3,k+1) = 0.0
      cu(i,1,k+1) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx+3
      do 30 i = 1, 3
      cu(i,j,3) = cu(i,j,3) + cu(i,j,ny+3)
      cu(i,j,2) = cu(i,j,2) + cu(i,j,ny+2)
      cu(i,j,ny+1) = cu(i,j,ny+1) + cu(i,j,1)
      cu(i,j,ny+3) = 0.0
      cu(i,j,ny+2) = 0.0
      cu(i,j,1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 3
      cu(i,2,3) = cu(i,2,3) + cu(i,nx+2,3)
      cu(i,3,3) = cu(i,3,3) + cu(i,nx+3,3)
      cu(i,nx+1,3) = cu(i,nx+1,3) + cu(i,1,3)
      cu(i,nx+2,3) = 0.0
      cu(i,nx+3,3) = 0.0
      cu(i,1,3) = 0.0
      cu(i,2,2) = cu(i,2,2) + cu(i,nx+2,2)
      cu(i,3,2) = cu(i,3,2) + cu(i,nx+3,2)
      cu(i,nx+1,2) = cu(i,nx+1,2) + cu(i,1,2)
      cu(i,nx+2,2) = 0.0
      cu(i,nx+3,2) = 0.0
      cu(i,1,2) = 0.0
      cu(i,2,ny+1) = cu(i,2,ny+1) + cu(i,nx+2,ny+1)
      cu(i,3,ny+1) = cu(i,3,ny+1) + cu(i,nx+3,ny+1)
      cu(i,nx+1,ny+1) = cu(i,nx+1,ny+1) + cu(i,1,ny+1)
      cu(i,nx+2,ny+1) = 0.0
      cu(i,nx+3,ny+1) = 0.0
      cu(i,1,ny+1) = 0.0
   50 continue
      return
      end
      subroutine ACGUARD22(cu,nx,ny,nxe,nye)
c accumulate extended periodic 2 component vector field
c quadratic interpolation
c handles special case of ny=1 properly
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, 2
      cu(i,2,k+1) = cu(i,2,k+1) + cu(i,nx+2,k+1)
      cu(i,3,k+1) = cu(i,3,k+1) + cu(i,nx+3,k+1)
      cu(i,nx+1,k+1) = cu(i,nx+1,k+1) + cu(i,1,k+1)
      cu(i,nx+2,k+1) = 0.0
      cu(i,nx+3,k+1) = 0.0
      cu(i,1,k+1) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx+3
      do 30 i = 1, 2
      cu(i,j,3) = cu(i,j,3) + cu(i,j,ny+3)
      cu(i,j,2) = cu(i,j,2) + cu(i,j,ny+2)
      cu(i,j,ny+1) = cu(i,j,ny+1) + cu(i,j,1)
      cu(i,j,ny+3) = 0.0
      cu(i,j,ny+2) = 0.0
      cu(i,j,1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 2
      cu(i,2,3) = cu(i,2,3) + cu(i,nx+2,3)
      cu(i,3,3) = cu(i,3,3) + cu(i,nx+3,3)
      cu(i,nx+1,3) = cu(i,nx+1,3) + cu(i,1,3)
      cu(i,nx+2,3) = 0.0
      cu(i,nx+3,3) = 0.0
      cu(i,1,3) = 0.0
      cu(i,2,2) = cu(i,2,2) + cu(i,nx+2,2)
      cu(i,3,2) = cu(i,3,2) + cu(i,nx+3,2)
      cu(i,nx+1,2) = cu(i,nx+1,2) + cu(i,1,2)
      cu(i,nx+2,2) = 0.0
      cu(i,nx+3,2) = 0.0
      cu(i,1,2) = 0.0
      cu(i,2,ny+1) = cu(i,2,ny+1) + cu(i,nx+2,ny+1)
      cu(i,3,ny+1) = cu(i,3,ny+1) + cu(i,nx+3,ny+1)
      cu(i,nx+1,ny+1) = cu(i,nx+1,ny+1) + cu(i,1,ny+1)
      cu(i,nx+2,ny+1) = 0.0
      cu(i,nx+3,ny+1) = 0.0
      cu(i,1,ny+1) = 0.0
   50 continue
      return
      end
      subroutine AGUARD2(q,nx,ny,nxe,nye)
c accumulate extended periodic scalar field
c quadratic interpolation
c handles special case of ny=1 properly
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
      integer j, k
c accumulate edges of extended field
      do 10 k = 1, ny
      q(2,k+1) = q(2,k+1) + q(nx+2,k+1)
      q(3,k+1) = q(3,k+1) + q(nx+3,k+1)
      q(nx+1,k+1) = q(nx+1,k+1) + q(1,k+1)
      q(nx+2,k+1) = 0.0
      q(nx+3,k+1) = 0.0
      q(1,k+1) = 0.0
   10 continue
      do 20 j = 1, nx+3
      q(j,3) = q(j,3) + q(j,ny+3)
      q(j,2) = q(j,2) + q(j,ny+2)
      q(j,ny+1) = q(j,ny+1) + q(j,1)
      q(j,ny+3) = 0.0
      q(j,ny+2) = 0.0
      q(j,1) = 0.0
   20 continue
      q(2,3) = q(2,3) + q(nx+2,3)
      q(3,3) = q(3,3) + q(nx+3,3)
      q(nx+1,3) = q(nx+1,3) + q(1,3)
      q(nx+2,3) = 0.0
      q(nx+3,3) = 0.0
      q(1,3) = 0.0
      q(2,2) = q(2,2) + q(nx+2,2)
      q(3,2) = q(3,2) + q(nx+3,2)
      q(nx+1,2) = q(nx+1,2) + q(1,2)
      q(nx+2,2) = 0.0
      q(nx+3,2) = 0.0
      q(1,2) = 0.0
      q(2,ny+1) = q(2,ny+1) + q(nx+2,ny+1)
      q(3,ny+1) = q(3,ny+1) + q(nx+3,ny+1)
      q(nx+1,ny+1) = q(nx+1,ny+1) + q(1,ny+1)
      q(nx+2,ny+1) = 0.0
      q(nx+3,ny+1) = 0.0
      q(1,ny+1) = 0.0
      return
      end
      subroutine CGUARD2L(fxy,nx,ny,nxe,nye)
c replicate extended periodic field
c linear interpolation
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer j, k
      do 10 k = 1, ny
      fxy(1,nx+1,k) = fxy(1,1,k)
      fxy(2,nx+1,k) = fxy(2,1,k)
   10 continue
      do 20 j = 1, nx
      fxy(1,j,ny+1) = fxy(1,j,1)
      fxy(2,j,ny+1) = fxy(2,j,1)
   20 continue
      fxy(1,nx+1,ny+1) = fxy(1,1,1)
      fxy(2,nx+1,ny+1) = fxy(2,1,1)
      return
      end
      subroutine DGUARD2L(q,nx,ny,nxe,nye)
c replicate extended scalar periodic field
c linear interpolation
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
      do 10 k = 1, ny
      q(nx+1,k) = q(1,k)
   10 continue
      do 20 j = 1, nx
      q(j,ny+1) = q(j,1)
   20 continue
      q(nx+1,ny+1) = q(1,1)
      return
      end
      subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
c replicate extended vector periodic field
c linear interpolation
      implicit none
      real bxy
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer j, k
      do 10 k = 1, ny
      !bxy(1,nx+1,k) = bxy(1,1,k)
      !bxy(2,nx+1,k) = bxy(2,1,k)
      !bxy(3,nx+1,k) = bxy(3,1,k)
   10 continue
      do 20 j = 1, nx
      bxy(1,j,ny+1) = bxy(1,j,1)
      bxy(2,j,ny+1) = bxy(2,j,1)
      bxy(3,j,ny+1) = bxy(3,j,1)
   20 continue
      bxy(1,nx+1,ny+1) = bxy(1,1,1)
      bxy(2,nx+1,ny+1) = bxy(2,1,1)
      bxy(3,nx+1,ny+1) = bxy(3,1,1)
      return
      end
      subroutine SCGUARD2L(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
c initialize extended periodic field
c linear interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      cu(1,j,k) = xj0
      cu(2,j,k) = yj0
      cu(3,j,k) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,nx+1,k) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 3
      cu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 i = 1, 3
      cu(i,nx+1,ny+1) = 0.
   60 continue
      return
      end
      subroutine SCGUARD22L(cu,xj0,yj0,nx,ny,nxe,nye)
c initialize extended periodic field
c linear interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      cu(1,j,k) = xj0
      cu(2,j,k) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,nx+1,k) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 2
      cu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 i = 1, 2
      cu(i,nx+1,ny+1) = 0.
   60 continue
      return
      end
      subroutine SGUARD2L(q,qi0,nx,ny,nxe,nye)
c initialize extended periodic scalar field
c linear interpolation
      implicit none
      real q, qi0
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c initialize extended field, with zero in the edges
      do 20 k = 1, ny
      do 10 j = 1, nx
      q(j,k) = qi0
   10 continue
      q(nx+1,k) = 0.
   20 continue
      do 30 j = 1, nx
      q(j,ny+1) = 0.
   30 continue
      q(nx+1,ny+1) = 0.
      return
      end
      subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
c accumulate extended periodic field
c linear interpolation
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, 3
      cu(i,1,k) = cu(i,1,k) + cu(i,nx+1,k)
      cu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 3
      cu(i,j,1) = cu(i,j,1) + cu(i,j,ny+1)
      cu(i,j,ny+1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 3
      cu(i,1,1) = cu(i,1,1) + cu(i,nx+1,ny+1)
      cu(i,nx+1,ny+1) = 0.0
   50 continue
      return
      end
      subroutine ACGUARD22L(cu,nx,ny,nxe,nye)
c accumulate extended periodic field
c linear interpolation
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, 2
      cu(i,1,k) = cu(i,1,k) + cu(i,nx+1,k)
      cu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 2
      cu(i,j,1) = cu(i,j,1) + cu(i,j,ny+1)
      cu(i,j,ny+1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 2
      cu(i,1,1) = cu(i,1,1) + cu(i,nx+1,ny+1)
      cu(i,nx+1,ny+1) = 0.0
   50 continue
      return
      end
      subroutine AVCGUARD2L(cu,vcu,nx,ny,nxe,nye)
c accumulate extended periodic current
c linear interpolation
      implicit none
      real cu, vcu
      integer nx, ny, nxe, nye, j
      dimension cu(3,nxe,nye), vcu(3,nx+3,ny+3)
!for 1-D case
      if(ny==1) then
        do 10 j = 1, nx+3
           vcu(1,j,2) = sum(vcu(1,j,:))
           vcu(2,j,2) = sum(vcu(2,j,:))
           vcu(3,j,2) = sum(vcu(3,j,:))
   10 continue
!j_x
      cu(1,1:nx,1) = cu(1,1:nx,1) + vcu(1,2:nx+1,2)!a
      cu(1,nx,1) = cu(1,nx,1) + vcu(1,1,2)!b
      cu(1,1,1) = cu(1,1,1) + vcu(1,nx+2,2)!i
!j_y
      cu(2,1:nx,1) = cu(2,1:nx,1) + vcu(2,2:nx+1,2)!a
      cu(2,nx,1) = cu(2,nx,1) + vcu(2,1,2)!d
      cu(2,1:2,1) = cu(2,1:2,1) + vcu(2,nx+2:nx+3,2)!b
!j_z
      cu(3,1:nx,1) = cu(3,1:nx,1) + vcu(3,2:nx+1,2)!a
      cu(3,nx,1) = cu(3,nx,1) + vcu(3,1,2)!d
      cu(3,1:2,1) = cu(3,1:2,1) + vcu(3,nx+2:nx+3,2)!b
      else
!for regular 2D case     
!For j_x
      cu(1,1:nx,1:ny) = cu(1,1:nx,1:ny) + vcu(1,2:nx+1,2:ny+1)!a
      cu(1,nx,1:ny) = cu(1,nx,1:ny) + vcu(1,1,2:ny+1)!b
      cu(1,1:nx,1:2) = cu(1,1:nx,1:2) + vcu(1,2:nx+1,ny+2:ny+3)!c
      cu(1,nx,1:2) = cu(1,nx,1:2) + vcu(1,1,ny+2:ny+3)!d
      cu(1,1:nx,ny) = cu(1,1:nx,ny) + vcu(1,2:nx+1,1)!e
      cu(1,nx,ny) = cu(1,nx,ny) + vcu(1,1,1)!f
      cu(1,1,1:2) = cu(1,1,1:2) + vcu(1,nx+2,ny+2:ny+3)!h
      cu(1,1,1:ny) = cu(1,1,1:ny) + vcu(1,nx+2,2:ny+1)!i
      cu(1,1,ny) = cu(1,1,ny) + vcu(1,nx+2,1)!l
!for j_y
      cu(2,1:nx,1:ny) = cu(2,1:nx,1:ny) + vcu(2,2:nx+1,2:ny+1)!a
      cu(2,1:nx,ny) = cu(2,1:nx,ny) + vcu(2,2:nx+1,1)!b
      cu(2,1:2,1:ny) = cu(2,1:2,1:ny) + vcu(2,nx+2:nx+3,2:ny+1)!c
      cu(2,1:2,ny) = cu(2,1:2,ny) + vcu(2,nx+2:nx+3,1)!d
      cu(2,nx,1:ny) = cu(2,nx,1:ny) + vcu(2,1,2:ny+1)!e
      cu(2,nx,ny) = cu(2,nx,ny) + vcu(2,1,1)!f
      cu(2,1:2,1) = cu(2,1:2,1) + vcu(2,nx+2:nx+3,ny+2)!h
      cu(2,1:nx,1) = cu(2,1:nx,1) + vcu(2,2:nx+1,ny+2)!i
      cu(2,nx,1) = cu(2,nx,1) + vcu(2,1,ny+2)!l 
!for j_z
      cu(3,1:nx,1:ny) = cu(3,1:nx,1:ny) + vcu(3,2:nx+1,2:ny+1)!a
      cu(3,1:2,1:ny) = cu(3,1:2,1:ny) + vcu(3,nx+2:nx+3,2:ny+1)!b
      !cu(3,1,1:ny) = cu(3,1,1:ny) + vcu(3,nx+2,2:ny+1)!c
      cu(3,nx,1:ny) = cu(3,nx,1:ny) + vcu(3,1,2:ny+1)!d
      cu(3,1:nx,1:2) = cu(3,1:nx,1:2) + vcu(3,2:nx+1,ny+2:ny+3)!e
      cu(3,1:nx,ny) = cu(3,1:nx,ny) + vcu(3,2:nx+1,1)!g
      cu(3,1:2,1:2) = cu(3,1:2,1:2) + vcu(3,nx+2:nx+3,ny+2:ny+3)!h
      cu(3,nx,1:2) = cu(3,nx,1:2) + vcu(3,1,ny+2:ny+3)!k, double split at corner
      cu(3,1:2,ny) = cu(3,1:2,ny) + vcu(3,nx+2:nx+3,1)!l, DS at corner
      cu(3,nx,ny) = cu(3,nx,ny) + vcu(3,1,1)!o, DS at corner
      endif
      return
      end
      subroutine AVCGUARD2(cu,vcu,nx,ny,nxe,nye)
c accumulate extended periodic current
c quadratic interpolation
      implicit none
      real cu, vcu
      integer nx, ny, nxe, nye, j
      dimension cu(3,nxe,nye), vcu(3,nx+4,ny+4)
      if(ny==1) then
!For 1D case
        do 10 j = 1, nx+4
           vcu(1,j,3) = sum(vcu(1,j,:))
           vcu(2,j,3) = sum(vcu(2,j,:))
           vcu(3,j,3) = sum(vcu(3,j,:))
   10 continue
!j_x
      cu(1,2:nx+1,2) = cu(1,2:nx+1,2) + vcu(1,3:nx+2,3)!a
      cu(1,nx:nx+1,2) = cu(1,nx:nx+1,2) + vcu(1,1:2,3)!b
      cu(1,2,2) = cu(1,2,2) + vcu(1,nx+3,3)!i
!j_y
      cu(2,2:nx+1,2) = cu(2,2:nx+1,2) + vcu(2,3:nx+2,3)!a
      cu(2,nx:nx+1,2) = cu(2,nx:nx+1,2) + vcu(2,1:2,3)!e
      cu(2,1:2,2) = cu(2,1:2,2) + vcu(2,nx+3:nx+4,3)!b
!j_z
      cu(3,2:nx+1,2) = cu(3,2:nx+1,2) + vcu(3,3:nx+2,3)!a
      cu(3,nx:nx+1,2) = cu(3,nx:nx+1,2) + vcu(3,1:2,3)!e
      cu(3,1:2,2) = cu(3,1:2,2) + vcu(3,nx+3:nx+4,3)!b
      else
!For regular 2D case 
!For j_x
      cu(1,2:nx+1,2:ny+1) = cu(1,2:nx+1,2:ny+1) + vcu(1,3:nx+2,3:ny+2)!a
      cu(1,nx:nx+1,2:ny+1) = cu(1,nx:nx+1,2:ny+1) + vcu(1,1:2,3:ny+2)!b
      cu(1,2:nx+1,ny:ny+1) = cu(1,2:nx+1,ny:ny+1) + vcu(1,3:nx+2,1:2)!e
      cu(1,2:nx+1,2:3) = cu(1,2:nx+1,2:3) + vcu(1,3:nx+2,ny+3:ny+4)!c
      cu(1,nx:nx+1,ny:ny+1) = cu(1,nx:nx+1,ny:ny+1)  + vcu(1,1:2,1:2)!f
      cu(1,nx:nx+1,2:3) = cu(1,nx:nx+1,2:3) + vcu(1,1:2,ny+3:ny+4)!d
      cu(1,2,2:ny+1) = cu(1,2,2:ny+1) + vcu(1,nx+3,3:ny+2)!i
      cu(1,2,ny:ny+1) = cu(1,2,ny:ny+1) + vcu(1,nx+3,1:2)!l
      cu(1,2,2:3) = cu(1,2,2:3) + vcu(1,nx+3,ny+3:ny+4)!h
!for j_y
      cu(2,2:nx+1,2:ny+1) = cu(2,2:nx+1,2:ny+1) + vcu(2,3:nx+2,3:ny+2)!a
      cu(2,2:nx+1,ny:ny+1) = cu(2,2:nx+1,ny:ny+1) + vcu(2,3:nx+2,1:2)!b
      cu(2,nx:nx+1,2:ny+1) = cu(2,nx:nx+1,2:ny+1) + vcu(2,1:2,3:ny+2)!e
      cu(2,2:3,2:ny+1) = cu(2,2:3,2:ny+1) + vcu(2,nx+3:nx+4,3:ny+2)!c
      cu(2,nx:nx+1,ny:ny+1) = cu(2,nx:nx+1,ny:ny+1)  + vcu(2,1:2,1:2)!f
      cu(2,2:3,ny:ny+1) = cu(2,2:3,ny:ny+1) + vcu(2,nx+3:nx+4,1:2)!d
      cu(2,2:nx+1,2) = cu(2,2:nx+1,2) + vcu(2,3:nx+2,ny+3)!i
      cu(2,nx:nx+1,2) = cu(2,nx:nx+1,2) + vcu(2,1:2,ny+3)!l
      cu(2,2:3,2) = cu(2,2:3,2) + vcu(2,nx+3:nx+4,ny+3)!h
!for j_z
      cu(3,2:nx+1,2:ny+1) = cu(3,2:nx+1,2:ny+1) + vcu(3,3:nx+2,3:ny+2)!a
      cu(3,2:3,2:ny+1) = cu(3,2:3,2:ny+1) + vcu(3,nx+3:nx+4,3:ny+2)!c
      cu(3,nx:nx+1,ny:ny+1) = cu(3,nx:nx+1,ny:ny+1) + vcu(3,1:2,1:2)!f
      cu(3,nx:nx+1,2:ny+1) = cu(3,nx:nx+1,2:ny+1) + vcu(3,1:2,3:ny+2)!e
      cu(3,2:nx+1,2:3) = cu(3,2:nx+1,2:3) + vcu(3,3:nx+2,ny+3:ny+4)!i
      cu(3,2:nx+1,ny:ny+1) = cu(3,2:nx+1,ny:ny+1) + vcu(3,3:nx+2,1:2)!b
      cu(3,2:3,2:3) = cu(3,2:3,2:3) + vcu(3,nx+3:nx+4,ny+3:ny+4)!h
      cu(3,nx:nx+1,2:3) = cu(3,nx:nx+1,2:3) + vcu(3,1:2,ny+3:ny+4)!l
      cu(3,2:3,ny:ny+1) = cu(3,2:3,ny:ny+1) + vcu(3,nx+3:nx+4,1:2)!d
      endif
      return
      end
      subroutine AVCGUARD2C(cu,vcu,nx,ny,nxe,nye)
c accumulate extended periodic current
c cubic interpolation
      implicit none
      real cu, vcu
      integer nx, ny, nxe, nye, j
      dimension cu(3,nxe,nye), vcu(3,nx+5,ny+5)
      if(ny==1) then
!consider 1D case.
        do 10 j = 1, nx+5
           vcu(1,j,3) = sum(vcu(1,j,:))
           vcu(2,j,3) = sum(vcu(2,j,:))
           vcu(3,j,3) = sum(vcu(3,j,:))
   10 continue
!j_x
      cu(1,3:nx+2,3) = cu(1,3:nx+2,3) + vcu(1,3:nx+2,3)!a
      cu(1,nx+1:nx+2,3) = cu(1,nx+1:nx+2,3) + vcu(1,1:2,3)!b
      cu(1,3:4,3) = cu(1,3:4,3) + vcu(1,nx+3:nx+4,3)!i
!j_y
      cu(2,3:nx+2,3) = cu(2,3:nx+2,3) + vcu(2,3:nx+2,3)!a
      cu(2,nx+1:nx+2,3) = cu(2,nx+1:nx+2,3) + vcu(2,1:2,3)!e
      cu(2,3:5,3) = cu(2,3:5,3) + vcu(2,nx+3:nx+5,3)!c
!j_z
      cu(3,3:nx+2,3) = cu(3,3:nx+2,3) + vcu(3,3:nx+2,3)!a
      cu(3,nx+1:nx+2,3) = cu(3,nx+1:nx+2,3) + vcu(3,1:2,3)!e
      cu(3,3:5,3) = cu(3,3:5,3) + vcu(3,nx+3:nx+5,3)!c
      else 
!regular 2D case
!For j_x
      cu(1,3:nx+2,3:ny+2) = cu(1,3:nx+2,3:ny+2) + vcu(1,3:nx+2,
     1 3:ny+2)!a
      cu(1,nx+1:nx+2,3:ny+2) = cu(1,nx+1:nx+2,3:ny+2) + vcu(1,1:2,
     1 3:ny+2)!b
      cu(1,3:nx+2,ny+1:ny+2) = cu(1,3:nx+2,ny+1:ny+2) + vcu(1,3:nx+2,
     1 1:2)!e
      cu(1,3:nx+2,3:5) = cu(1,3:nx+2,3:5) + vcu(1,3:nx+2,ny+3:ny+5)!c
      cu(1,nx+1:nx+2,ny+1:ny+2) = cu(1,nx+1:nx+2,ny+1:ny+2) 
     1 + vcu(1,1:2,1:2)!f
      cu(1,3:4,3:ny+2) = cu(1,3:4,3:ny+2) + vcu(1,nx+3:nx+4,3:ny+2)!i
      cu(1,3:4,3:5) = cu(1,3:4,3:5) + vcu(1,nx+3:nx+4,ny+3:ny+5)!k
      cu(1,3:4,ny+1:ny+2) = cu(1,3:4,ny+1:ny+2) 
     1 + vcu(1,nx+3:nx+4,1:2)!q
      cu(1,nx+1:nx+2,3:5) = cu(1,nx+1:nx+2,3:5) + vcu(1,1:2,
     1 ny+3:ny+5)!j
!for j_y
      cu(2,3:nx+2,3:ny+2) = cu(2,3:nx+2,3:ny+2) + vcu(2,3:nx+2,
     1 3:ny+2)!a
      cu(2,3:nx+2,ny+1:ny+2) = cu(2,3:nx+2,ny+1:ny+2) + vcu(2,3:nx+2,
     1 1:2)!b
      cu(2,nx+1:nx+2,3:ny+2) = cu(2,nx+1:nx+2,3:ny+2) + vcu(2,1:2,
     1 3:ny+2)!e
      cu(2,3:5,3:ny+2) = cu(2,3:5,3:ny+2) + vcu(2,nx+3:nx+5,3:ny+2)!c
      cu(2,nx+1:nx+2,ny+1:ny+2) = cu(2,nx+1:nx+2,ny+1:ny+2) 
     1 + vcu(2,1:2,1:2)!f
      cu(2,3:nx+2,3:4) = cu(2,3:nx+2,3:4) + vcu(2,3:nx+2,ny+3:ny+4)!i
      cu(2,3:5,3:4) = cu(2,3:5,3:4) + vcu(2,nx+3:nx+5,ny+3:ny+4)!k
      cu(2,nx+1:nx+2,3:4) = cu(2,nx+1:nx+2,3:4) 
     1 + vcu(2,1:2,ny+3:ny+4)!q
      cu(2,3:5,ny+1:ny+2) = cu(2,3:5,ny+1:ny+2) + vcu(2,nx+3:nx+5,
     1 1:2)!j
!for j_z
      cu(3,3:nx+2,3:ny+2) = cu(3,3:nx+2,3:ny+2) + vcu(3,3:nx+2,3:ny+2)!a
      cu(3,3:nx+2,ny+1:ny+2) = cu(3,3:nx+2,ny+1:ny+2) + vcu(3,3:nx+2,
     1 1:2)!b
      cu(3,nx+1:nx+2,3:ny+2) = cu(3,nx+1:nx+2,3:ny+2) + vcu(3,1:2,
     1 3:ny+2)!e
      cu(3,3:5,3:ny+2) = cu(3,3:5,3:ny+2) + vcu(3,nx+3:nx+5,3:ny+2)!c
      cu(3,3:nx+2,3:5) = cu(3,3:nx+2,3:5) + vcu(3,3:nx+2,ny+3:ny+5)!i
      cu(3,nx+1:nx+2,3:5) = cu(3,nx+1:nx+2,3:5) + vcu(3,1:2,ny+3:ny+5)!u
      cu(3,3:5,3:5) = cu(3,3:5,3:5) + vcu(3,nx+3:nx+5,ny+3:ny+5)!k
      cu(3,nx+1:nx+2,ny+1:ny+2) = cu(3,nx+1:nx+2,ny+1:ny+2) + 
     1 vcu(3,1:2,1:2)!f
      cu(3,3:5,ny+1:ny+2) = cu(3,3:5,ny+1:ny+2) + 
     1 vcu(3,nx+3:nx+5,1:2)!d
! For the special case of nx or ny <=2, we need to double wrap
      if(ny.le.2) then
        cu(:,:,3) = cu(:,:,3) + cu(:,:,5)
        cu(:,:,5) = 0.0
      endif
      if(nx.le.2) then
        cu(:,3,:) = cu(:,3,:) + cu(:,5,:)
        cu(:,5,:) = 0.0
      endif
      endif
      return
      end
      subroutine AGUARD2L(q,nx,ny,nxe,nye)
c accumulate extended periodic scalar field
c linear interpolation
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c accumulate edges of extended field
      do 10 k = 1, ny
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1) = q(j,1) + q(j,ny+1)
      q(j,ny+1) = 0.0
   20 continue
      q(1,1) = q(1,1) + q(nx+1,ny+1)
      q(nx+1,ny+1) = 0.0
      return
      end
      subroutine CGUARD2C(fxy,nx,ny,nxe,nye,ndim)
c replicate extended periodic field
c cubic interpolation
      implicit none
      real fxy
      integer nx, ny, nxe, nye, ndim
      dimension fxy(ndim,nxe,nye)
c local data
      integer i, j, k
      do 20 k = 1, ny
      do 10 i = 1, ndim
      fxy(i,1,k+2) = fxy(i,nx+1,k+2)
      fxy(i,2,k+2) = fxy(i,nx+2,k+2)
      fxy(i,nx+3,k+2) = fxy(i,3,k+2)
      fxy(i,nx+4,k+2) = fxy(i,4,k+2)
      fxy(i,nx+5,k+2) = fxy(i,5,k+2)
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      fxy(i,j+2,1) = fxy(i,j+2,ny+1)
      fxy(i,j+2,2) = fxy(i,j+2,ny+2)
      fxy(i,j+2,ny+3) = fxy(i,j+2,3)
      fxy(i,j+2,ny+4) = fxy(i,j+2,4)
      fxy(i,j+2,ny+5) = fxy(i,j+2,5)
   30 continue
   40 continue
      do 50 i = 1, ndim
      fxy(i,1,1) = fxy(i,nx+1,ny+1)
      fxy(i,2,1) = fxy(i,nx+2,ny+1)
      fxy(i,nx+3,1) = fxy(i,3,ny+1)
      fxy(i,nx+4,1) = fxy(i,4,ny+1)
      fxy(i,nx+5,1) = fxy(i,5,ny+1)
      fxy(i,1,2) = fxy(i,nx+1,ny+2)
      fxy(i,2,2) = fxy(i,nx+2,ny+2)
      fxy(i,nx+3,2) = fxy(i,3,ny+2)
      fxy(i,nx+4,2) = fxy(i,4,ny+2)
      fxy(i,nx+5,2) = fxy(i,5,ny+2)
      fxy(i,1,ny+3) = fxy(i,nx+1,3)
      fxy(i,2,ny+3) = fxy(i,nx+2,3)
      fxy(i,nx+3,ny+3) = fxy(i,3,3)
      fxy(i,nx+4,ny+3) = fxy(i,4,3)
      fxy(i,nx+5,ny+3) = fxy(i,5,3)
      fxy(i,1,ny+4) = fxy(i,nx+1,4)
      fxy(i,2,ny+4) = fxy(i,nx+2,4)
      fxy(i,nx+3,ny+4) = fxy(i,3,4)
      fxy(i,nx+4,ny+4) = fxy(i,4,4)
      fxy(i,nx+5,ny+4) = fxy(i,5,4)
      fxy(i,1,ny+5) = fxy(i,nx+1,5)
      fxy(i,2,ny+5) = fxy(i,nx+2,5)
      fxy(i,nx+3,ny+5) = fxy(i,3,5)
      fxy(i,nx+4,ny+5) = fxy(i,4,5)
      fxy(i,nx+5,ny+5) = fxy(i,5,5)
   50 continue
      return
      end
      subroutine DGUARD2C(q,nx,ny,nxe,nye)
c replicate extended periodic field
c cubic interpolation
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
      do 10 k = 1, ny
      q(1,k+2) = q(nx+1,k+2)
      q(2,k+2) = q(nx+2,k+2)
      q(nx+3,k+2) = q(3,k+2)
      q(nx+4,k+2) = q(4,k+2)
      q(nx+5,k+2) = q(5,k+2)
   10 continue
      do 20 j = 1, nx
      q(j+2,1) = q(j+2,ny+1)
      q(j+2,2) = q(j+2,ny+2)
      q(j+2,ny+3) = q(j+2,3)
      q(j+2,ny+4) = q(j+2,4)
      q(j+2,ny+5) = q(j+2,5)
   20 continue
      q(1,1) = q(nx+1,ny+1)
      q(2,1) = q(nx+2,ny+1)
      q(nx+3,1) = q(3,ny+1)
      q(nx+4,1) = q(4,ny+1)
      q(nx+5,1) = q(5,ny+1)
      q(1,2) = q(nx+1,ny+2)
      q(2,2) = q(nx+2,ny+2)
      q(nx+3,2) = q(3,ny+2)
      q(nx+4,2) = q(4,ny+2)
      q(nx+5,2) = q(5,ny+2)
      q(1,ny+3) = q(nx+1,3)
      q(2,ny+3) = q(nx+2,3)
      q(nx+3,ny+3) = q(3,3)
      q(nx+4,ny+3) = q(4,3)
      q(nx+5,ny+3) = q(5,3)
      q(1,ny+4) = q(nx+1,4)
      q(2,ny+4) = q(nx+2,4)
      q(nx+3,ny+4) = q(3,4)
      q(nx+4,ny+4) = q(4,4)
      q(nx+5,ny+4) = q(5,4)
      q(1,ny+5) = q(nx+1,5)
      q(2,ny+5) = q(nx+2,5)
      q(nx+3,ny+5) = q(3,5)
      q(nx+4,ny+5) = q(4,5)
      q(nx+5,ny+5) = q(5,5)
      return
      end
      subroutine SCGUARD2C(cu,xj0,yj0,zj0,nx,ny,nxe,nye)
c initialize extended periodic vector field
c cubic interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      cu(1,j+2,k+2) = xj0
      cu(2,j+2,k+2) = yj0
      cu(3,j+2,k+2) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+2) = 0.
      cu(i,2,k+2) = 0.
      cu(i,nx+3,k+2) = 0.
      cu(i,nx+4,k+2) = 0.
      cu(i,nx+5,k+2) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 3
      cu(i,j+2,1) = 0.
      cu(i,j+2,2) = 0.
      cu(i,j+2,ny+3) = 0.
      cu(i,j+2,ny+4) = 0.
      cu(i,j+2,ny+5) = 0.
   40 continue
   50 continue
      do 60 i = 1, 3
      cu(i,1,1) = 0.
      cu(i,2,1) = 0.
      cu(i,nx+3,1) = 0.
      cu(i,nx+4,1) = 0.
      cu(i,nx+5,1) = 0.
      cu(i,1,2) = 0.
      cu(i,2,2) = 0.
      cu(i,nx+3,2) = 0.
      cu(i,nx+4,2) = 0.
      cu(i,nx+5,2) = 0.
      cu(i,1,ny+3) = 0.
      cu(i,2,ny+3) = 0.
      cu(i,nx+3,ny+3) = 0.
      cu(i,nx+4,ny+3) = 0.
      cu(i,nx+5,ny+3) = 0.
      cu(i,1,ny+4) = 0.
      cu(i,2,ny+4) = 0.
      cu(i,nx+3,ny+4) = 0.
      cu(i,nx+4,ny+4) = 0.
      cu(i,nx+5,ny+4) = 0.
      cu(i,1,ny+5) = 0.
      cu(i,2,ny+5) = 0.
      cu(i,nx+3,ny+5) = 0.
      cu(i,nx+4,ny+5) = 0.
      cu(i,nx+5,ny+5) = 0.
   60 continue
      return
      end
      subroutine SCGUARD22C(cu,xj0,yj0,nx,ny,nxe,nye)
c initialize extended periodic vector field
c cubic interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      cu(1,j+2,k+2) = xj0
      cu(2,j+2,k+2) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+2) = 0.
      cu(i,2,k+2) = 0.
      cu(i,nx+3,k+2) = 0.
      cu(i,nx+4,k+2) = 0.
      cu(i,nx+5,k+2) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 2
      cu(i,j+2,1) = 0.
      cu(i,j+2,2) = 0.
      cu(i,j+2,ny+3) = 0.
      cu(i,j+2,ny+4) = 0.
      cu(i,j+2,ny+5) = 0.
   40 continue
   50 continue
      do 60 i = 1, 2
      cu(i,1,1) = 0.
      cu(i,2,1) = 0.
      cu(i,nx+3,1) = 0.
      cu(i,nx+4,1) = 0.
      cu(i,nx+5,1) = 0.
      cu(i,1,2) = 0.
      cu(i,2,2) = 0.
      cu(i,nx+3,2) = 0.
      cu(i,nx+4,2) = 0.
      cu(i,nx+5,2) = 0.
      cu(i,1,ny+3) = 0.
      cu(i,2,ny+3) = 0.
      cu(i,nx+3,ny+3) = 0.
      cu(i,nx+4,ny+3) = 0.
      cu(i,nx+5,ny+3) = 0.
      cu(i,1,ny+4) = 0.
      cu(i,2,ny+4) = 0.
      cu(i,nx+3,ny+4) = 0.
      cu(i,nx+4,ny+4) = 0.
      cu(i,nx+5,ny+4) = 0.
      cu(i,1,ny+5) = 0.
      cu(i,2,ny+5) = 0.
      cu(i,nx+3,ny+5) = 0.
      cu(i,nx+4,ny+5) = 0.
      cu(i,nx+5,ny+5) = 0.
   60 continue
      return
      end
      subroutine SGUARD2C(q,qi0,nx,ny,nxe,nye)
c initialize extended periodic scalar field
c cubic interpolation
      implicit none
      real q, qi0
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c initialize extended field, with zero in the edges
      do 20 k = 1, ny
      do 10 j = 1, nx
      q(j+2,k+2) = qi0
   10 continue
      q(1,k+2) = 0.
      q(2,k+2) = 0.
      q(nx+3,k+2) = 0.
      q(nx+4,k+2) = 0.
      q(nx+5,k+2) = 0.
   20 continue
      do 30 j = 1, nx
      q(j+2,1) = 0.
      q(j+2,2) = 0.
      q(j+2,ny+3) = 0.
      q(j+2,ny+4) = 0.
      q(j+2,ny+5) = 0.
   30 continue
      q(1,1) = 0.
      q(2,1) = 0.
      q(nx+3,1) = 0.
      q(nx+4,1) = 0.
      q(nx+5,1) = 0.
      q(1,2) = 0.
      q(2,2) = 0.
      q(nx+3,2) = 0.
      q(nx+4,2) = 0.
      q(nx+5,2) = 0.
      q(1,ny+3) = 0.
      q(2,ny+3) = 0.
      q(nx+3,ny+3) = 0.
      q(nx+4,ny+3) = 0.
      q(nx+5,ny+3) = 0.
      q(1,ny+4) = 0.
      q(2,ny+4) = 0.
      q(nx+3,ny+4) = 0.
      q(nx+4,ny+4) = 0.
      q(nx+5,ny+4) = 0.
      q(1,ny+5) = 0.
      q(2,ny+5) = 0.
      q(nx+3,ny+5) = 0.
      q(nx+4,ny+5) = 0.
      q(nx+5,ny+5) = 0.
      return
      end
      subroutine ACGUARD2C(cu,nx,ny,nxe,nye,ndim)
c accumulate extended periodic vector field
c cubic interpolation
      implicit none
      real cu
      integer nx, ny, nxe, nye, ndim
      dimension cu(ndim,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, ndim
      cu(i,3,k+2) = cu(i,3,k+2) + cu(i,nx+3,k+2)
      cu(i,4,k+2) = cu(i,4,k+2) + cu(i,nx+4,k+2)
      cu(i,5,k+2) = cu(i,5,k+2) + cu(i,nx+5,k+2)
      cu(i,nx+1,k+2) = cu(i,nx+1,k+2) + cu(i,1,k+2)
      cu(i,nx+2,k+2) = cu(i,nx+2,k+2) + cu(i,2,k+2)
      cu(i,nx+3,k+2) = 0.0
      cu(i,nx+4,k+2) = 0.0
      cu(i,nx+5,k+2) = 0.0
      cu(i,1,k+2) = 0.0
      cu(i,2,k+2) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      cu(i,j+2,3) = cu(i,j+2,3) + cu(i,j+2,ny+3)
      cu(i,j+2,4) = cu(i,j+2,4) + cu(i,j+2,ny+4)
      cu(i,j+2,5) = cu(i,j+2,5) + cu(i,j+2,ny+5)
      cu(i,j+2,ny+1) = cu(i,j+2,ny+1) + cu(i,j+2,1)
      cu(i,j+2,ny+2) = cu(i,j+2,ny+2) + cu(i,j+2,2)
      cu(i,j+2,ny+3) = 0.0
      cu(i,j+2,ny+4) = 0.0
      cu(i,j+2,ny+5) = 0.0
      cu(i,j+2,1) = 0.0
      cu(i,j+2,2) = 0.0
   30 continue
   40 continue
      do 50 i = 1, ndim
      cu(i,3,3) = cu(i,3,3) + cu(i,nx+3,ny+3)
      cu(i,4,3) = cu(i,4,3) + cu(i,nx+4,ny+3)
      cu(i,5,3) = cu(i,5,3) + cu(i,nx+5,ny+3)
      cu(i,nx+1,3) = cu(i,nx+1,3) + cu(i,1,ny+3)
      cu(i,nx+2,3) = cu(i,nx+2,3) + cu(i,2,ny+3)
      cu(i,nx+3,ny+3) = 0.0
      cu(i,nx+4,ny+3) = 0.0
      cu(i,nx+5,ny+3) = 0.0
      cu(i,1,ny+3) = 0.0
      cu(i,2,ny+3) = 0.0
      cu(i,3,4) = cu(i,3,4) + cu(i,nx+3,ny+4)
      cu(i,4,4) = cu(i,4,4) + cu(i,nx+4,ny+4)
      cu(i,5,4) = cu(i,5,4) + cu(i,nx+5,ny+4)
      cu(i,nx+1,4) = cu(i,nx+1,4) + cu(i,1,ny+4)
      cu(i,nx+2,4) = cu(i,nx+2,4) + cu(i,2,ny+4)
      cu(i,nx+3,ny+4) = 0.0
      cu(i,nx+4,ny+4) = 0.0
      cu(i,nx+5,ny+4) = 0.0
      cu(i,1,ny+4) = 0.0
      cu(i,2,ny+4) = 0.0
      cu(i,3,5) = cu(i,3,5) + cu(i,nx+3,ny+5)
      cu(i,4,5) = cu(i,4,5) + cu(i,nx+4,ny+5)
      cu(i,5,5) = cu(i,5,5) + cu(i,nx+5,ny+5)
      cu(i,nx+1,5) = cu(i,nx+1,5) + cu(i,1,ny+5)
      cu(i,nx+2,5) = cu(i,nx+2,5) + cu(i,2,ny+5)
      cu(i,nx+3,ny+5) = 0.0
      cu(i,nx+4,ny+5) = 0.0
      cu(i,nx+5,ny+5) = 0.0
      cu(i,1,ny+5) = 0.0
      cu(i,2,ny+5) = 0.0
      cu(i,3,ny+1) = cu(i,3,ny+1) + cu(i,nx+3,1)
      cu(i,4,ny+1) = cu(i,4,ny+1) + cu(i,nx+4,1)
      cu(i,5,ny+1) = cu(i,5,ny+1) + cu(i,nx+5,1)
      cu(i,nx+1,ny+1) = cu(i,nx+1,ny+1) + cu(i,1,1)
      cu(i,nx+2,ny+1) = cu(i,nx+2,ny+1) + cu(i,2,1)
      cu(i,nx+3,1) = 0.0
      cu(i,nx+4,1) = 0.0
      cu(i,nx+5,1) = 0.0
      cu(i,1,1) = 0.0
      cu(i,2,1) = 0.0
      cu(i,3,ny+2) = cu(i,3,ny+2) + cu(i,nx+3,2)
      cu(i,4,ny+2) = cu(i,4,ny+2) + cu(i,nx+4,2)
      cu(i,5,ny+2) = cu(i,5,ny+2) + cu(i,nx+5,2)
      cu(i,nx+1,ny+2) = cu(i,nx+1,ny+2) + cu(i,1,2)
      cu(i,nx+2,ny+2) = cu(i,nx+2,ny+2) + cu(i,2,2)
      cu(i,nx+3,2) = 0.0
      cu(i,nx+4,2) = 0.0
      cu(i,nx+5,2) = 0.0
      cu(i,1,2) = 0.0
      cu(i,2,2) = 0.0
   50 continue
      return
      end
      subroutine AGUARD2C(q,nx,ny,nxe,nye)
c accumulate extended periodic scalar field
c cubic interpolation
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k
c accumulate edges of extended field
      do 10 k = 1, ny
      q(3,k+2) = q(3,k+2) + q(nx+3,k+2)
      q(4,k+2) = q(4,k+2) + q(nx+4,k+2)
      q(5,k+2) = q(5,k+2) + q(nx+5,k+2)
      q(nx+1,k+2) = q(nx+1,k+2) + q(1,k+2)
      q(nx+2,k+2) = q(nx+2,k+2) + q(2,k+2)
      q(nx+3,k+2) = 0.0
      q(nx+4,k+2) = 0.0
      q(nx+5,k+2) = 0.0
      q(1,k+2) = 0.0
      q(2,k+2) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j+2,3) = q(j+2,3) + q(j+2,ny+3)
      q(j+2,4) = q(j+2,4) + q(j+2,ny+4)
      q(j+2,5) = q(j+2,5) + q(j+2,ny+5)
      q(j+2,ny+1) = q(j+2,ny+1) + q(j+2,1)
      q(j+2,ny+2) = q(j+2,ny+2) + q(j+2,2)
      q(j+2,ny+3) = 0.0
      q(j+2,ny+4) = 0.0
      q(j+2,ny+5) = 0.0
      q(j+2,1) = 0.0
      q(j+2,2) = 0.0
   20 continue
      q(3,3) = q(3,3) + q(nx+3,ny+3)
      q(4,3) = q(4,3) + q(nx+4,ny+3)
      q(5,3) = q(5,3) + q(nx+5,ny+3)
      q(nx+1,3) = q(nx+1,3) + q(1,ny+3)
      q(nx+2,3) = q(nx+2,3) + q(2,ny+3)
      q(nx+3,ny+3) = 0.0
      q(nx+4,ny+3) = 0.0
      q(nx+5,ny+3) = 0.0
      q(1,ny+3) = 0.0
      q(2,ny+3) = 0.0
      q(3,4) = q(3,4) + q(nx+3,ny+4)
      q(4,4) = q(4,4) + q(nx+4,ny+4)
      q(5,4) = q(5,4) + q(nx+5,ny+4)
      q(nx+1,4) = q(nx+1,4) + q(1,ny+4)
      q(nx+2,4) = q(nx+2,4) + q(2,ny+4)
      q(nx+3,ny+4) = 0.0
      q(nx+4,ny+4) = 0.0
      q(nx+5,ny+4) = 0.0
      q(1,ny+4) = 0.0
      q(2,ny+4) = 0.0
      q(3,5) = q(3,5) + q(nx+3,ny+5)
      q(4,5) = q(4,5) + q(nx+4,ny+5)
      q(5,5) = q(5,5) + q(nx+5,ny+5)
      q(nx+1,5) = q(nx+1,5) + q(1,ny+5)
      q(nx+2,5) = q(nx+2,5) + q(2,ny+5)
      q(nx+3,ny+5) = 0.0
      q(nx+4,ny+5) = 0.0
      q(nx+5,ny+5) = 0.0
      q(1,ny+5) = 0.0
      q(2,ny+5) = 0.0
      q(3,ny+1) = q(3,ny+1) + q(nx+3,1)
      q(4,ny+1) = q(4,ny+1) + q(nx+4,1)
      q(5,ny+1) = q(5,ny+1) + q(nx+5,1)
      q(nx+1,ny+1) = q(nx+1,ny+1) + q(1,1)
      q(nx+2,ny+1) = q(nx+2,ny+1) + q(2,1)
      q(nx+3,1) = 0.0
      q(nx+4,1) = 0.0
      q(nx+5,1) = 0.0
      q(1,1) = 0.0
      q(2,1) = 0.0
      q(3,ny+2) = q(3,ny+2) + q(nx+3,2)
      q(4,ny+2) = q(4,ny+2) + q(nx+4,2)
      q(5,ny+2) = q(5,ny+2) + q(nx+5,2)
      q(nx+1,ny+2) = q(nx+1,ny+2) + q(1,2)
      q(nx+2,ny+2) = q(nx+2,ny+2) + q(2,2)
      q(nx+3,2) = 0.0
      q(nx+4,2) = 0.0
      q(nx+5,2) = 0.0
      q(1,2) = 0.0
      q(2,2) = 0.0
!Add all y indicies together for each x index for the special case
!of one cell wide in y.
      if(ny==1) then
        do 25 j = 1, nxe
        q(j,3) = sum(q(j,:))
        q(j,1:2) = 0.0
        q(j,4:nye) = 0.0
   25 continue
      endif 
      return
      end
      subroutine POISP2(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,nxv,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffc
c for isign = -1, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fx,fy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c for isign = 1, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fx,we
c approximate flop count is: 14*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fy
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c cmplx(q(2*j-1,k),q(2*j,k)) = complex charge density
c for fourier mode (j-1,k-1)
c cmplx(fx(2*j-1,k),fx(2*j,k)) = x component of complex force/charge,
c cmplx(fy(2*j-1,k),fy(2*j,k)) = y component of complex force/charge,
c for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c ffc(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffc(2*j-1,k) = potential green's function g for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      dimension q(nxv,ny+1), fx(nxv,ny+1), fy(nxv,ny+1)
      dimension ffc(nxv,nyhd)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      if (at3.eq.0.) then
         ffc(2*j,k) = 1.0
         ffc(2*j-1,k) = affp
      else
         ffc(2*j,k) = exp(-.5*((dkx*ax)**2 + at2))
         ffc(2*j-1,k) = affp*ffc(2*j,k)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = ffc(2*j-1,k)*ffc(2*j,k)
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      fx(2*j-1,k) = at2*q(2*j,k)
      fx(2*j,k) = -at2*q(2*j-1,k)
      fx(2*j-1,k1) = at2*q(2*j,k1)
      fx(2*j,k1) = -at2*q(2*j-1,k1)
      fy(2*j-1,k) = at3*q(2*j,k)
      fy(2*j,k) = -at3*q(2*j-1,k)
      fy(2*j-1,k1) = -at3*q(2*j,k1)
      fy(2*j,k1) = at3*q(2*j-1,k1)
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2 + q(2*j-1,k1)**2 + q(2*
     1j,k1)**2)
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ffc(1,k)*ffc(2,k)
      at3 = dny*float(k - 1)*at1
      fx(1,k) = 0.0
      fx(2,k) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      fy(1,k) = at3*q(2,k)
      fy(2,k) = -at3*q(1,k)
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = ffc(2*j-1,1)*ffc(2*j,1)
      at2 = dnx*float(j - 1)*at1
      fx(2*j-1,1) = at2*q(2*j,1)
      fx(2*j,1) = -at2*q(2*j-1,1)
      fx(2*j-1,k1) = 0.0
      fx(2*j,k1) = 0.0
      fy(2*j-1,1) = 0.0
      fy(2*j,1) = 0.0
      fy(2*j-1,k1) = 0.0
      fy(2*j,k1) = 0.0
      wp = wp + at1*(q(2*j-1,1)**2 + q(2*j,1)**2)
   70 continue
      fx(1,1) = 0.0
      fx(2,1) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      fy(1,1) = 0.0
      fy(2,1) = 0.0
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
      we = float(nx*ny)*wp
      return
c calculate potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at2 = ffc(2*j-1,k)
      at1 = at2*ffc(2*j,k)
      fx(2*j-1,k) = at2*q(2*j-1,k)
      fx(2*j,k) = at2*q(2*j,k)
      fx(2*j-1,k1) = at2*q(2*j-1,k1)
      fx(2*j,k1) = at2*q(2*j,k1)
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2 + q(2*j-1,k1)**2 + q(2*
     1j,k1)**2)
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = ffc(1,k)
      at1 = at2*ffc(2,k)
      fx(1,k) = at2*q(1,k)
      fx(2,k) = at2*q(2,k)
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at2 = ffc(2*j-1,1)
      at1 = at2*ffc(2*j,1)
      fx(2*j-1,1) = at2*q(2*j-1,1)
      fx(2*j,1) = at2*q(2*j,1)
      fx(2*j-1,k1) = 0.0
      fx(2*j,k1) = 0.0
      wp = wp + at1*(q(2*j-1,1)**2 + q(2*j,1)**2)
  120 continue
      fx(1,1) = 0.0
      fx(2,1) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      we = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      k1 = ny2 - k
      do 140 j = 2, nxh
      at1 = ffc(2*j,k)
      fy(2*j-1,k) = at1*q(2*j-1,k)
      fy(2*j,k) = at1*q(2*j,k)
      fy(2*j-1,k1) = at1*q(2*j-1,k1)
      fy(2*j,k1) = at1*q(2*j,k1)
  140 continue
  150 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = ffc(2,k)
      fy(1,k) = at1*q(1,k)
      fy(2,k) = at1*q(2,k)
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nxh
      at1 = ffc(2*j,1)
      fy(2*j-1,1) = at1*q(2*j-1,1)
      fy(2*j,1) = at1*q(2*j,1)
      fy(2*j-1,k1) = 0.0
      fy(2*j,k1) = 0.0
  170 continue
      fy(1,1) = ffc(2,1)*q(1,1)
      fy(2,1) = 0.0
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
      return
      end
      subroutine POISP2X(q,fx,fy,isign,ffc,ax,ay,affp,we,nx,ny,nxv,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with periodic boundary conditions.
c explicitly handles the nx/2 + 1 mode
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffc
c for isign = -1, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fx,fy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c for isign = 1, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fx,we
c approximate flop count is: 14*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fy
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c cmplx(q(2*j-1,k),q(2*j,k)) = complex charge density
c for fourier mode (j-1,k-1)
c cmplx(fx(2*j-1,k),fx(2*j,k)) = x component of complex force/charge,
c cmplx(fy(2*j-1,k),fy(2*j,k)) = y component of complex force/charge,
c for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c ffc(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffc(2*j-1,k) = potential green's function g for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+2
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxv, nyhd
      real ax, ay, affp, we
      real q, fx, fy, ffc
      dimension q(nxv,ny+1), fx(nxv,ny+1), fy(nxv,ny+1)
      dimension ffc(nxv,nyhd)
c local data
      integer j, k, nxh, nyh, ny2, k1
      real dnx, dny, dkx, dky, at1, at2, at3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      if (at3.eq.0.) then
         ffc(2*j,k) = 1.0
         ffc(2*j-1,k) = affp
      else
         ffc(2*j,k) = exp(-.5*((dkx*ax)**2 + at2))
         ffc(2*j-1,k) = affp*ffc(2*j,k)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = ffc(2*j-1,k)*ffc(2*j,k)
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      fx(2*j-1,k) = at2*q(2*j,k)
      fx(2*j,k) = -at2*q(2*j-1,k)
      fx(2*j-1,k1) = at2*q(2*j,k1)
      fx(2*j,k1) = -at2*q(2*j-1,k1)
      fy(2*j-1,k) = at3*q(2*j,k)
      fy(2*j,k) = -at3*q(2*j-1,k)
      fy(2*j-1,k1) = -at3*q(2*j,k1)
      fy(2*j,k1) = at3*q(2*j-1,k1)
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2 + q(2*j-1,k1)**2 + q(2*
     1j,k1)**2)
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ffc(1,k)*ffc(2,k)
      at3 = dny*float(k - 1)*at1
      at2 = at3*q(2,k)
      at3 = at3*q(1,k)
      fx(1,k) = 0.0
      fx(2,k) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      fy(1,k) = at2
      fy(2,k) = -at3
      fy(1,k1) = at2
      fy(2,k1) = at3
      fx(nx+1,k) = 0.0
      fx(nx+2,k) = 0.0
      fx(nx+1,k1) = 0.0
      fx(nx+2,k1) = 0.0
      fy(nx+1,k) = 0.0
      fy(nx+2,k) = 0.0
      fy(nx+1,k1) = 0.0
      fy(nx+2,k1) = 0.0
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = ffc(2*j-1,1)*ffc(2*j,1)
      at2 = dnx*float(j - 1)*at1
      fx(2*j-1,1) = at2*q(2*j,1)
      fx(2*j,1) = -at2*q(2*j-1,1)
      fx(2*j-1,k1) = 0.0
      fx(2*j,k1) = 0.0
      fy(2*j-1,1) = 0.0
      fy(2*j,1) = 0.0
      fy(2*j-1,k1) = 0.0
      fy(2*j,k1) = 0.0
      wp = wp + at1*(q(2*j-1,1)**2 + q(2*j,1)**2)
   70 continue
      fx(1,1) = 0.0
      fx(2,1) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      fy(1,1) = 0.0
      fy(2,1) = 0.0
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
      fx(nx+1,1) = 0.0
      fx(nx+2,1) = 0.0
      fx(nx+1,k1) = 0.0
      fx(nx+2,k1) = 0.0
      fy(nx+1,1) = 0.0
      fy(nx+2,1) = 0.0
      fy(nx+1,k1) = 0.0
      fy(nx+2,k1) = 0.0
      we = float(nx*ny)*wp
      return
c calculate potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at2 = ffc(2*j-1,k)
      at1 = at2*ffc(2*j,k)
      fx(2*j-1,k) = at2*q(2*j-1,k)
      fx(2*j,k) = at2*q(2*j,k)
      fx(2*j-1,k1) = at2*q(2*j-1,k1)
      fx(2*j,k1) = at2*q(2*j,k1)
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2 + q(2*j-1,k1)**2 + q(2*
     1j,k1)**2)
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = ffc(1,k)
      at1 = at2*ffc(2,k)
      at3 = at2*q(1,k)
      at2 = at2*q(2,k)
      fx(1,k) = at3
      fx(2,k) = at2
      fx(1,k1) = at3
      fx(2,k1) = -at2
      fx(nx+1,k) = 0.0
      fx(nx+2,k) = 0.0
      fx(nx+1,k1) = 0.0
      fx(nx+2,k1) = 0.0
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at2 = ffc(2*j-1,1)
      at1 = at2*ffc(2*j,1)
      fx(2*j-1,1) = at2*q(2*j-1,1)
      fx(2*j,1) = at2*q(2*j,1)
      fx(2*j-1,k1) = 0.0
      fx(2*j,k1) = 0.0
      wp = wp + at1*(q(2*j-1,1)**2 + q(2*j,1)**2)
  120 continue
      fx(1,1) = 0.0
      fx(2,1) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      fx(nx+1,1) = 0.0
      fx(nx+2,1) = 0.0
      fx(nx+1,k1) = 0.0
      fx(nx+2,k1) = 0.0
      we = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      k1 = ny2 - k
      do 140 j = 2, nxh
      at1 = ffc(2*j,k)
      fy(2*j-1,k) = at1*q(2*j-1,k)
      fy(2*j,k) = at1*q(2*j,k)
      fy(2*j-1,k1) = at1*q(2*j-1,k1)
      fy(2*j,k1) = at1*q(2*j,k1)
  140 continue
  150 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = ffc(2,k)
      at2 = at1*q(1,k)
      at3 = at1*q(2,k)
      fy(1,k) = at2
      fy(2,k) = at3
      fy(1,k1) = at2
      fy(2,k1) = -at3
      fy(nx+1,k) = 0.0
      fy(nx+2,k) = 0.0
      fy(nx+1,k1) = 0.0
      fy(nx+2,k1) = 0.0
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nxh
      at1 = ffc(2*j,1)
      fy(2*j-1,1) = at1*q(2*j-1,1)
      fy(2*j,1) = at1*q(2*j,1)
      fy(2*j-1,k1) = 0.0
      fy(2*j,k1) = 0.0
  170 continue
      fy(1,1) = ffc(2,1)*q(1,1)
      fy(2,1) = 0.0
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
      fy(nx+1,1) = 0.0
      fy(nx+2,1) = 0.0
      fy(nx+1,k1) = 0.0
      fy(nx+2,k1) = 0.0
      return
      end
      subroutine POIS22cf(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,n
     1xhd,nyhd,smoothfactor,smoothmask)
c this subroutine initialize ffc in a "cut-off" way
      double precision wp,kke,jje
      complex q, fxy, ffc, zero, zt1, zt2
      dimension q(nxvh,nyv), fxy(2,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.0,0.0)
      print *, 'ffc cut-off'
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      kke = real(k**2)/real( (smoothmask*nyh)**2 )
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      jje = real(j**2)/real( (smoothmask*nxh)**2 )
      !if ( kke + jje .le. 1 ) then
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         if ((kke + jje) .le. 1.0) then
            ffc(j,k) = cmplx(affp*at4/at3,at4)
         else
            ffc(j,k) = cmplx(affp*at4/at3*smoothfactor,at4*smoothfactor)
         endif
      endif
   10 continue
   20 continue
   30 return
      end
      subroutine POIS22(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nxh
     1d,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffc, zero, zt1, zt2
      integer yee
      real kxop, kyop, thetax, thetay, dkx, dky
      dimension q(nxvh,nyv), fxy(2,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
!      if(yee == 1) then
!        thetay = dky/2.0
!        kyop = 2.0*sin(thetay)
!        at1 = kyop*kyop
!      else
        at1 = dky*dky
!      endif
      !print*, 'at1 =', at1
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
!      if(yee == 1) then
!        thetax = dkx/2.0
!        kxop = 2.0*sin(thetax)
!        at3 = kxop*kxop + at1
!      else
        at3 = dkx*dkx + at1
!      endif
      !print*, 'at3 =', at3
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
      !print*, 'affp =', affp
!      print*, 'k =', k, 'thetay = ', thetay, 'j =', j, 
!     1 'thetax =', thetax, 'at3 =', at3, 'ffc =', ffc(j,k)
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine FLDSHIFT2(exy,bxy,isign,nx,ny,nxvh,nyv)
!This subroutine shifts the fourier space fields between the 
!integer grid points where the charge is defined, and the 
!locations on the yee cube where they are defined. isign = -1 
!shifts them to the integer points, and isign = 1 shifts to 
!the yee cell positions.
      double precision wp
      complex zero, exy, bxy
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      real kxop, kyop, thetax, thetay, dkx, dky
      complex phix, phiy
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.0,0.0)

c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = float(isign)*dky/2.0
      kyop = 2.0*sin(thetay)
      phiy = cmplx(cos(thetay),sin(thetay))
      do 40 j = 2, nxh
      dkx = dnx*(j - 1)
      thetax = float(isign)*dkx/2.0
      kxop = 2.0*sin(thetax)
      phix = cmplx(cos(thetax),sin(thetax))
      exy(1,j,k) = exy(1,j,k)*phix
      exy(2,j,k) = exy(2,j,k)*phiy
      exy(3,j,k) = exy(3,j,k)
      exy(1,j,k1) = exy(1,j,k1)*phix
      exy(2,j,k1) = exy(2,j,k1)*conjg(phiy)
      exy(3,j,k1) = exy(3,j,k1)
      bxy(1,j,k) = bxy(1,j,k)*phiy
      bxy(2,j,k) = bxy(2,j,k)*phix
      bxy(3,j,k) = bxy(3,j,k)*phix*phiy
      bxy(1,j,k1) = bxy(1,j,k1)*conjg(phiy)
      bxy(2,j,k1) = bxy(2,j,k1)*phix
      bxy(3,j,k1) = bxy(3,j,k1)*phix*conjg(phiy)
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      dkx = dnx*float(nx/2)
      thetax = float(isign)*dkx/2.0
      kxop = 2.0*sin(thetax)
      phix = cmplx(cos(thetax),sin(thetax))     
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = float(isign)*dky/2.0
      kyop = 2.0*sin(thetay)
      phiy = cmplx(cos(thetay),sin(thetay))
      exy(1,1,k) = exy(1,1,k)
      exy(2,1,k) = exy(2,1,k)*phiy
      exy(3,1,k) = exy(3,1,k)
      exy(1,1,k1) = exy(1,1,k1)*conjg(phix)
      exy(2,1,k1) = exy(2,1,k1)*conjg(phiy)
      exy(3,1,k1) = exy(3,1,k1)
      bxy(1,1,k) = bxy(1,1,k)*phiy
      bxy(2,1,k) = bxy(2,1,k)
      bxy(3,1,k) = bxy(3,1,k)*phiy
      bxy(1,1,k1) = bxy(1,1,k1)*conjg(phiy)
      bxy(2,1,k1) = bxy(2,1,k1)*conjg(phix)
      bxy(3,1,k1) = bxy(3,1,k1)*conjg(phix)*conjg(phiy)
   60 continue
c mode numbers ky = 0, ny/2 (blocks e and k)
      dky = dny*float(ny/2)
      thetay = float(isign)*dky/2.0
      kyop = 2.0*sin(thetay)
      phiy = cmplx(cos(thetay),sin(thetay))
      k1 = nyh + 1
      do 70 j = 2, nxh
      dkx = dnx*(j - 1)
      thetax = float(isign)*dkx/2.0
      kxop = 2.0*sin(thetax)
      phix = cmplx(cos(thetax),sin(thetax))
      exy(1,j,1) = exy(1,j,1)*phix
      exy(2,j,1) = exy(2,j,1)
      exy(3,j,1) = exy(3,j,1)
      exy(1,j,k1) = exy(1,j,k1)*phix
      exy(2,j,k1) = exy(2,j,k1)*conjg(phiy)
      exy(3,j,k1) = exy(3,j,k1)
      bxy(1,j,1) = bxy(1,j,1)
      bxy(2,j,1) = bxy(2,j,1)*phix
      bxy(3,j,1) = bxy(3,j,1)*phix
      bxy(1,j,k1) = bxy(1,j,k1)*conjg(phiy)
      bxy(2,j,k1) = bxy(2,j,k1)*phix
      bxy(3,j,k1) = bxy(3,j,k1)*phix*conjg(phiy)
   70 continue
!      exy(1,1,1) = exy(1,1,1) ! (kx=0,ky=0;kx=nx/2,ky=0) => (no shift, zero)
      exy(1,1,1) = cmplx(real(exy(1,1,1)),0.0)      
      exy(2,1,1) = exy(2,1,1) ! (0,0;nx/2,0) => (no shift, no shift)
      exy(3,1,1) = exy(3,1,1) ! (0,0;nx/2,0) => (no shift, no shift)
!      exy(1,1,k1) = exy(1,1,k1) ! (0,ny/2;nx/2,ny/2) => (no shift, zero)
      exy(1,1,k1) = cmplx(real(exy(1,1,k1)),0.0)  
!      exy(2,1,k1) = exy(2,1,k1) ! (0,ny/2;nx/2,ny/2) => (zero, zero)
      exy(2,1,k1) = zero
      exy(3,1,k1) = exy(3,1,k1) ! (0,ny/2;nx/2,ny/2) => (no shift, no shift)
      bxy(1,1,1) = bxy(1,1,1) ! (0,0;nx/2,0) => (no shift, no shift)
!      bxy(2,1,1) = bxy(2,1,1) ! (0,0;nx/2,0) => (no shift, zero)
      bxy(2,1,1) = cmplx(real(bxy(2,1,1)),0.0) 
!      bxy(3,1,1) = bxy(3,1,1) ! (0,0;nx/2,0) => (no shift, zero)
      bxy(3,1,1) = cmplx(real(bxy(3,1,1)),0.0) 
!      bxy(1,1,k1) = bxy(1,1,k1) ! (0,ny/2;nx/2,ny/2) => (zero, zero)
      bxy(1,1,k1) = zero
!      bxy(2,1,k1) = bxy(2,1,k1) ! (0,ny/2;nx/2,ny/2) => (no shift, zero)
      bxy(2,1,k1) = cmplx(real(bxy(2,1,k1)),0.0) 
!      bxy(3,1,k1) = bxy(3,1,k1) ! (0,ny/2;nx/2,ny/2) => (zero, zero)
      bxy(3,1,k1) = zero
      return
      end
      subroutine FLDMODINIT2(fxy,exy,bxy,nx,ny,nxvh,nyv)
!This subroutine takes the force fxyz from the poisson solver at
!the initial time step and multiplies it by k/[k] so it satisfies
!the finite difference consinuity equation in Fourier space. It
!then copies fxyz to exyz for later use in the maxwel2 solver
!and shifts exyz to its yee cell position. It also shifts bxyz
!from the integer grid points to the yee cell positions
      double precision wp
      complex q, fxy, zero, exy, bxy
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      real kxop, kyop, thetax, thetay, dkx, dky
      complex phix, phiy, f
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.0,0.0)
   30 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      kyop = 2.0*sin(thetay)
      phiy = cmplx(cos(thetay),sin(thetay))
      do 40 j = 2, nxh
      dkx = dnx*(j - 1)
      thetax = dkx/2.0
      kxop = 2.0*sin(thetax)
      phix = cmplx(cos(thetax),sin(thetax))
!      fxy(1,j,k) = fxy(1,j,k)*dkx/kxop*phix
!      fxy(2,j,k) = fxy(2,j,k)*dky/kyop*phiy
!      fxy(3,j,k) = zero
!      fxy(1,j,k1) = fxy(1,j,k1)*dkx/kxop*phix
!      fxy(2,j,k1) = fxy(2,j,k1)*dky/kyop*conjg(phiy)
!      fxy(3,j,k1) = zero
      fxy(1,j,k) = fxy(1,j,k)*dkx/kxop
      fxy(2,j,k) = fxy(2,j,k)*dky/kyop
      fxy(3,j,k) = fxy(3,j,k)
      fxy(1,j,k1) = fxy(1,j,k1)*dkx/kxop
      fxy(2,j,k1) = fxy(2,j,k1)*dky/kyop
      fxy(3,j,k1) =  fxy(3,j,k1)
      exy(1,j,k) = fxy(1,j,k)*phix
      exy(2,j,k) = fxy(2,j,k)*phiy
      exy(3,j,k) = fxy(3,j,k1)
      exy(1,j,k1) = fxy(1,j,k1)*phix
      exy(2,j,k1) = fxy(2,j,k1)*conjg(phiy)
      exy(3,j,k1) = fxy(3,j,k1)
      bxy(1,j,k) = bxy(1,j,k)*phiy
      bxy(2,j,k) = bxy(2,j,k)*phix
      bxy(3,j,k) = bxy(3,j,k)*phix*phiy
      bxy(1,j,k1) = bxy(1,j,k1)*conjg(phiy)
      bxy(2,j,k1) = bxy(2,j,k1)*phix
      bxy(3,j,k1) = bxy(3,j,k1)*phix*conjg(phiy)
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      dkx = dnx*float(nx/2)
      thetax = dkx/2.0
      kxop = 2.0*sin(thetax)
      phix = cmplx(cos(thetax),sin(thetax)) 
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      kyop = 2.0*sin(thetay)
      phiy = cmplx(cos(thetay),sin(thetay))
!      fxy(1,1,k) = zero
!      fxy(2,1,k) = fxy(2,1,k)*dky/kyop!*phiy
!      fxy(3,1,k) = zero
!      fxy(1,1,k1) = fxy(1,1,k1)*dkx/kxop!*conjg(phix)
!      !print*, 'ex(k) =', fxy(1,1,k)
!      !print*, 'ex(k1) =', fxy(1,1,k1)
!      !print*, 'kxEx =', -kxop*fxy(1,1,k)
!      !print*, 'k1xEx =', -kxop*fxy(1,1,k1)
!      fxy(2,1,k1) = fxy(2,1,k1)*dky/kyop!*conjg(phiy)
!      fxy(3,1,k1) = zero
      fxy(1,1,k) = fxy(1,1,k)!zero
      fxy(2,1,k) = fxy(2,1,k)*dky/kyop
      fxy(3,1,k) = fxy(3,1,k)!zero
      fxy(1,1,k1) = fxy(1,1,k1)*dkx/kxop
      fxy(2,1,k1) = fxy(2,1,k1)*dky/kyop
      fxy(3,1,k1) = fxy(3,1,k1)!zero
      exy(1,1,k) = fxy(1,1,k)
      exy(2,1,k) = fxy(2,1,k)*phiy
      exy(3,1,k) = fxy(3,1,k)
      exy(1,1,k1) = fxy(1,1,k1)*conjg(phix)
      exy(2,1,k1) = fxy(2,1,k1)*conjg(phiy)
      exy(3,1,k1) = fxy(3,1,k1)
      bxy(1,1,k) = bxy(1,1,k)*phiy
      bxy(2,1,k) = bxy(2,1,k)
      bxy(3,1,k) = bxy(3,1,k)*phiy
      bxy(1,1,k1) = bxy(1,1,k1)*conjg(phiy)
      bxy(2,1,k1) = bxy(2,1,k1)*conjg(phix)
      bxy(3,1,k1) = bxy(3,1,k1)*conjg(phix)*conjg(phiy)
   60 continue
c mode numbers ky = 0, ny/2 (blocks e and k)
      dky = dny*float(ny/2)
      thetay = dky/2.0
      kyop = 2.0*sin(thetay)
      phiy = cmplx(cos(thetay),sin(thetay))
      k1 = nyh + 1
      do 70 j = 2, nxh
      dkx = dnx*(j - 1)
      thetax = dkx/2.0
      kxop = 2.0*sin(thetax)
      phix = cmplx(cos(thetax),sin(thetax))
!      fxy(1,j,1) = fxy(1,j,1)*dkx/kxop*phix
!      fxy(2,j,1) = zero
!      fxy(3,j,1) = zero
!      fxy(1,j,k1) = fxy(1,j,k1)*dkx/kxop*phix
!      fxy(2,j,k1) = fxy(2,j,k1)*dky/kyop*conjg(phiy)
!      fxy(3,j,k1) = zero
      fxy(1,j,1) = fxy(1,j,1)*dkx/kxop
      fxy(2,j,1) = fxy(2,j,1)!zero
      fxy(3,j,1) = fxy(3,j,1)!zero
      fxy(1,j,k1) = fxy(1,j,k1)*dkx/kxop
      fxy(2,j,k1) = fxy(2,j,k1)*dky/kyop
      fxy(3,j,k1) = fxy(3,j,k1)!zero
      exy(1,j,1) = fxy(1,j,1)*phix
      exy(2,j,1) = fxy(2,j,1)
      exy(3,j,1) = fxy(3,j,1)
      exy(1,j,k1) = fxy(1,j,k1)*phix
      exy(2,j,k1) = fxy(2,j,k1)*conjg(phiy)
      exy(3,j,k1) = fxy(3,j,k1)
      bxy(1,j,1) = bxy(1,j,1)
      bxy(2,j,1) = bxy(2,j,1)*phix
      bxy(3,j,1) = bxy(3,j,1)*phix
      bxy(1,j,k1) = bxy(1,j,k1)*conjg(phiy)
      bxy(2,j,k1) = bxy(2,j,k1)*phix
      bxy(3,j,k1) = bxy(3,j,k1)*phix*conjg(phiy)
   70 continue!Can't calculate these correctly, so zero them all out
      dkx = dnx*float(nx/2)
      kxop = 2.0*sin(dkx/2.0)
      phix = cmplx(cos(thetax),sin(thetax))
      !f = fxy(1,1,1)
      fxy(1,1,1) = zero! cmplx(0.0,aimag(f)*dkx/kxop)!*conjg(phix)
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      !f = fxy(1,1,k1)
      fxy(1,1,k1) = zero! cmplx(0.0,aimag(f)*dkx/kxop)!*conjg(phix)
      !f = fxy(2,1,k1)
      fxy(2,1,k1) = zero! cmplx(real(f)*dky/kyop,aimag(f)*dky/kyop)
      !1 *conjg(phiy)
      fxy(3,1,k1) = zero
      exy(1,1,1) = fxy(1,1,1)
      exy(2,1,1) = fxy(2,1,1)
      exy(3,1,1) = fxy(3,1,1)
      exy(1,1,k1) = fxy(1,1,k1)
      exy(2,1,k1) = fxy(2,1,k1)
      exy(3,1,k1) = fxy(3,1,k1)
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      return
      end
      subroutine POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nxh
     1d,nyhd,yee)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.  Zeros out z component.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh

!*****
!Now tests for if yee = 1. If it does, the solver executes a block
!of code identical to what it did before (i.e., with no [k]s), except
!now it also calculates the nx/2, ny/2 and 0 modes. The sub blocks are
!labeled according to fft2Dpacking.ods
!*****
      double precision wp
      complex q, fxy, ffc, zero, zt1, zt2
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      integer yee
      real kxop, kyop, thetax, thetay, dkx, dky, ffcT, affp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0

      if(yee==1) then
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 51 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 41 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   41 continue
   51 continue
c mode numbers kx = 0, nx/2 (blocks a and i)
cdir$ ivdep
      dkx = dnx*float(nx/2)
      do 61 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      ffcT = at1*dky*dky/(dky*dky + dkx*dkx)
      !at3 = dny*float(k - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      zt2 = cmplx(aimag(q(1,k1)),-real(q(1,k1)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
!We can't calculate these modes correctly because of how fft data is stored
      fxy(1,1,k1) = zero
      !fxy(1,1,k1) = -dkx*ffcT*zt2
      fxy(2,1,k1) = zero
      !fxy(2,1,k1) = -dky*ffcT*zt2
      fxy(3,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k))+ q(1,k1)*conjg(q(1,k1)))
   61 continue
c mode numbers ky = 0, ny/2 (blocks e and k)
      k1 = nyh + 1
      dky = dnx*float(ny/2)
      do 71 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      ffcT = at1*dkx*dkx/(dkx*dkx + dky*dky)
      at2 = dkx*at1
      at3 = dkx*fccT
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero!dkx*ffcT*zt2
      fxy(2,j,k1) = zero!-dky*ffcT*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   71 continue
!Unfortionatly, we can't calculate these correctly when we have to shift,
!so just zero them out for now. If we don't have to shift, then it works
!out, as explaned below.
      affp = at1*dkx*dkx
      dkx = dnx*float(nx/2)
      !Note that the aimag part is actually equal to
      ! -i*affp*(-1)/dkx*real(q(nx/2,0))*phix, where phix = -i
      fxy(1,1,1) = zero!cmplx(0.0, affp/dkx*aimag(q(1,1)))!
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      ffcT = affp/(dkx*dkx + dky*dky)
      !Note that the aimag part is actually equal to
      ! -i*affp*(-dkx)/dk^2*real(q(nx/2,ny/2))*phix, where phix = -i
      fxy(1,1,k1) = zero !cmplx(0.0, dkx*ffcT*aimag(q(1,k1)))
      !Note that the real part is actually equal to
      ! -i*affp*(-1)/dky*real(q(0,ny/2))*phiy, where phiy = -i
      !The aimag part is equal to
      !  -i*affp*(-dky)/dk^2*real(q(nx/2,ny/2))*phiy, where phiy = -i
      fxy(2,1,k1) = zero!cmplx(affp/dky*real(q(1,k1)),
!     1 dky*ffcT*aimag(q(1,k1)))
      fxy(3,1,k1) = zero
      we = float(nx*ny)*wp

      else

c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      !print*, fxy(3,j,1)
      !print*, fxy(3,j,k1)
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      !print*, fxy(3,1,k1)
      we = float(nx*ny)*wp

      !print*, fxy(3,:,:)
      endif
      !do j = 1, nxh+1
      !do k = 1, ny+1
      !if(fxy(3,j,k).ne.zero) print*, j,k,fxy(3,j,k) 
      !enddo
      !enddo
      !print*, 'pois', fxy(3,140,513)

      return
      end
      subroutine POIS23X(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,nx
     1hd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.  Zeros out z component.
c explicitly handles the nx/2 + 1 mode
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh+1
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer j, k, nxh, nyh, ny2, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = at3*conjg(zt1)
      fxy(3,1,k1) = zero
      fxy(1,nxh+1,k) = zero
      fxy(2,nxh+1,k) = zero
      fxy(3,nxh+1,k) = zero
      fxy(1,nxh+1,k1) = zero
      fxy(2,nxh+1,k1) = zero
      fxy(3,nxh+1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      fxy(1,nxh+1,1) = zero
      fxy(2,nxh+1,1) = zero
      fxy(3,nxh+1,1) = zero
      fxy(1,nxh+1,k1) = zero
      fxy(2,nxh+1,k1) = zero
      fxy(3,nxh+1,k1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine POISP2T(qt,fxt,fyt,isign,ffct,ax,ay,affp,we,nx,ny,nxvh,
     1nyv,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyv,nxhd,nyhd,
c output: ffct
c for isign = -1, input: qt,ffct,isign,nx,ny,nxvh,nyv,nxhd,nyhd,
c output: fxt,fyt,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c for isign = 1, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fx,we
c approximate flop count is: 14*nxc*nyc + 8*(nxc + nyc)
c for isign = 2, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fy
c approximate flop count is: 4*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c cmplx(qt(2*k-1,j),q(2*k,j)) = complex charge density
c for fourier mode (k-1,j-1)
c cmplx(fxt(2*k-1,j),fxt(2*k,j)) = x component of complex force/charge,
c cmplx(fyt(2*k-1,j),fyt(2*k,j)) = y component of complex force/charge,
c for fourier mode (k-1,j-1)
c if isign = 0, form factor array is prepared
c ffct(2*k,j) = finite-size particle shape factor s
c for fourier mode (k-1,j-1)
c ffct(2*k-1,j) = potential green's function g for fourier mode
c (k-1,j-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh+1
c nyv = first dimension of field arrays, must be >= ny
c nxhd = second dimension of form factor array, must be >= nxh
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      real qt, fxt, fyt, ffct
      dimension qt(2*nyv,nxvh), fxt(2*nyv,nxvh), fyt(2*nyv,nxvh)
      dimension ffct(2*nyhd,nxhd)
c local data
      integer j, k, nxh, nyh, ny2, j1, k1
      real dnx, dny, dkx, dky, at1, at2, at3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 j = 1, nxh
      dkx = dnx*real(j - 1)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      if (at3.eq.0.) then
         ffct(2*k,j) = 1.0
         ffct(2*k-1,j) = affp
      else
         ffct(2*k,j) = exp(-.5*((dky*ay)**2 + at2))
         ffct(2*k-1,j) = affp*ffct(2*k,j)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 j = 2, nxh
      dkx = dnx*real(j - 1)
      do 40 k = 2, nyh
      k1 = ny2 - k
      at1 = ffct(2*k-1,j)*ffct(2*k,j)
      at2 = dkx*at1
      at3 = dny*real(k - 1)*at1
      fxt(2*k-1,j) = at2*qt(2*k,j)
      fxt(2*k,j) = -at2*qt(2*k-1,j)
      fxt(2*k1-1,j) = at2*qt(2*k1,j)
      fxt(2*k1,j) = -at2*qt(2*k1-1,j)
      fyt(2*k-1,j) = at3*qt(2*k,j)
      fyt(2*k,j) = -at3*qt(2*k-1,j)
      fyt(2*k1-1,j) = -at3*qt(2*k1,j)
      fyt(2*k1,j) = at3*qt(2*k1-1,j)
      wp = wp + at1*(qt(2*k-1,j)**2 + qt(2*k,j)**2 + qt(2*k1-1,j)**2 + q
     1t(2*k1,j)**2)
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      at1 = ffct(1,j)*ffct(2,j)
      at2 = dnx*real(j - 1)*at1
      fxt(1,j) = at2*qt(2,j)
      fxt(2,j) = -at2*qt(1,j)
      fxt(2*k1-1,j) = 0.0
      fxt(2*k1,j) = 0.0
      fyt(1,j) = 0.0
      fyt(2,j) = 0.0
      fyt(2*k1-1,j) = 0.0
      fyt(2*k1,j) = 0.0
      wp = wp + at1*(qt(1,j)**2 + qt(2,j)**2)
   60 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 70 k = 2, nyh
      k1 = ny2 - k
      at1 = ffct(2*k-1,1)*ffct(2*k,1)
      at3 = dny*real(k - 1)*at1
      at2 = at3*qt(2*k,1)
      at3 = at3*qt(2*k-1,1)
      fxt(2*k-1,1) = 0.0
      fxt(2*k,1) = 0.0
      fxt(2*k1-1,1) = 0.0
      fxt(2*k1,1) = 0.0
      fyt(2*k-1,1) = at2
      fyt(2*k,1) = -at3
      fyt(2*k1-1,1) = at2
      fyt(2*k1,1) = at3
      fxt(2*k-1,j1) = 0.0
      fxt(2*k,j1) = 0.0
      fxt(2*k1-1,j1) = 0.0
      fxt(2*k1,j1) = 0.0
      fyt(2*k-1,j1) = 0.0
      fyt(2*k,j1) = 0.0
      fyt(2*k1-1,j1) = 0.0
      fyt(2*k1,j1) = 0.0
      wp = wp + at1*(qt(2*k-1,1)**2 + qt(2*k,1)**2)
   70 continue
      fxt(1,1) = 0.0
      fxt(2,1) = 0.0
      fxt(ny2-1,1) = 0.0
      fxt(ny2,1) = 0.0
      fyt(1,1) = 0.0
      fyt(2,1) = 0.0
      fyt(ny2-1,1) = 0.0
      fyt(ny2,1) = 0.0
      fxt(1,j1) = 0.0
      fxt(2,j1) = 0.0
      fxt(ny2-1,j1) = 0.0
      fxt(ny2,j1) = 0.0
      fyt(1,j1) = 0.0
      fyt(2,j1) = 0.0
      fyt(ny2-1,j1) = 0.0
      fyt(ny2,j1) = 0.0
      we = real(nx*ny)*wp
      return
c calculate potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 j = 2, nxh
      do 90 k = 2, nyh
      k1 = ny2 - k
      at2 = ffct(2*k-1,j)
      at1 = at2*ffct(2*k,j)
      fxt(2*k-1,j) = at2*qt(2*k-1,j)
      fxt(2*k,j) = at2*qt(2*k,j)
      fxt(2*k1-1,j) = at2*qt(2*k1-1,j)
      fxt(2*k1,j) = at2*qt(2*k1,j)
      wp = wp + at1*(qt(2*k-1,j)**2 + qt(2*k,j)**2 + qt(2*k1-1,j)**2 + q
     1t(2*k1,j)**2)
   90 continue
  100 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 110 j = 2, nxh
      at2 = ffct(1,j)
      at1 = at2*ffct(2,j)
      fxt(1,j) = at2*qt(1,j)
      fxt(2,j) = at2*qt(2,j)
      fxt(2*k1-1,j) = 0.0
      fxt(2*k1,j) = 0.0
      wp = wp + at1*(qt(1,j)**2 + qt(2,j)**2)
  110 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 120 k = 2, nyh
      k1 = ny2 - k
      at2 = ffct(2*k-1,1)
      at1 = at2*ffct(2*k,1)
      at3 = at2*qt(2*k-1,1)
      at2 = at2*qt(2*k,1)
      fxt(2*k-1,1) = at3
      fxt(2*k,1) = at2
      fxt(2*k1-1,1) = at3
      fxt(2*k1,1) = -at2
      fxt(2*k-1,j1) = 0.0
      fxt(2*k,j1) = 0.0
      fxt(2*k1-1,j1) = 0.0
      fxt(2*k1,j1) = 0.0
      wp = wp + at1*(qt(2*k-1,1)**2 + qt(2*k,1)**2)
  120 continue
      fxt(1,1) = 0.0
      fxt(2,1) = 0.0
      fxt(ny2-1,1) = 0.0
      fxt(ny2,1) = 0.0
      fxt(1,j1) = 0.0
      fxt(2,j1) = 0.0
      fxt(ny2-1,j1) = 0.0
      fxt(ny2,j1) = 0.0
      we = real(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  130 do 150 j = 2, nxh
      do 140 k = 2, nyh
      k1 = ny2 - k
      at1 = ffct(2*k,j)
      fyt(2*k-1,j) = at1*qt(2*k-1,j)
      fyt(2*k,j)= at1*qt(2*k,j)
      fyt(2*k1-1,j) = at1*qt(2*k1-1,j)
      fyt(2*k1,j) = at1*qt(2*k1,j)
  140 continue
  150 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 160 j = 2, nxh
      at1 = ffct(2,j)
      fyt(1,j) = at1*qt(1,j)
      fyt(2,j) = at1*qt(2,j)
      fyt(2*k1-1,j) = 0.0
      fyt(2*k1,j) = 0.0
  160 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 170 k = 2, nyh
      k1 = ny2 - k
      at1 = ffct(2*k,1)
      at2 = at1*qt(2*k-1,1)
      at3 = at1*qt(2*k,1)
      fyt(2*k-1,1) = at2
      fyt(2*k,1) = at3
      fyt(2*k1-1,1) = at2
      fyt(2*k1,1) = -at3
      fyt(2*k-1,j1) = 0.0
      fyt(2*k,j1) = 0.0
      fyt(2*k1-1,j1) = 0.0
      fyt(2*k1,j1) = 0.0
  170 continue
      fyt(1,1) = ffct(2,1)*qt(1,1)
      fyt(2,1) = 0.0
      fyt(ny2-1,1) = 0.0
      fyt(ny2,1) = 0.0
      fyt(1,j1) = 0.0
      fyt(2,j1) = 0.0
      fyt(ny2-1,j1) = 0.0
      fyt(ny2,j1) = 0.0
      return
      end
      subroutine POIS23T(qt,fxyt,isign,ffct,ax,ay,affp,we,nx,ny,nxvh,nyv
     1,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions.  Zeros out z component.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffct
c for isign /= 0, input: qt,ffct,isign,nx,ny,nxvh,nyhd, output: fxyt,we
c approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c qt(k,j) = complex charge density for fourier mode (k-1,j-1)
c fxyt(k,1,j) = x component of complex force/charge,
c fxyt(k,2,j) = y component of complex force/charge,
c fxyt(k,3,j) = zero,
c all for fourier mode (k-1,j-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffct(k,j)) = finite-size particle shape factor s
c for fourier mode (k-1,j-1)
c real(ffct(k,j)) = potential green's function g
c for fourier mode (k-1,j-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh+1
c nyv = first dimension of field arrays, must be >= ny
c nxhd = second dimension of form factor array, must be >= nxh
c nyhd = first dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex qt, fxyt, ffct
      dimension qt(nyv,nxvh), fxyt(nyv,3,nxvh)
      dimension ffct(nyhd,nxhd)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 j = 1, nxh
      dkx = dnx*real(j - 1)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, nyh
      dky = dny*real(k - 1)
      at3 = dky*dky + at1
      at4 = exp(-.5*((dky*ay)**2 + at2))
      if (at3.eq.0.) then
         ffct(k,j) = cmplx(affp,1.0)
      else
         ffct(k,j) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 j = 2, nxh
      dkx = dnx*real(j - 1)
      do 40 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffct(k,j))*aimag(ffct(k,j))
      at2 = dkx*at1
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(aimag(qt(k,j)),-real(qt(k,j)))
      zt2 = cmplx(aimag(qt(k1,j)),-real(qt(k1,j)))
      fxyt(k,1,j) = at2*zt1
      fxyt(k,2,j) = at3*zt1
      fxyt(k,3,j) = zero
      fxyt(k1,1,j) = at2*zt2
      fxyt(k1,2,j) = -at3*zt2
      fxyt(k1,3,j) = zero
      wp = wp + at1*(qt(k,j)*conjg(qt(k,j)) + qt(k1,j)*conjg(qt(k1,j)))
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      at1 = real(ffct(1,j))*aimag(ffct(1,j))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(qt(1,j)),-real(qt(1,j)))
      fxyt(1,1,j) = at2*zt1
      fxyt(1,2,j) = zero
      fxyt(1,3,j) = zero
      fxyt(k1,1,j) = zero
      fxyt(k1,2,j) = zero
      fxyt(k1,3,j) = zero
      wp = wp + at1*(qt(1,j)*conjg(qt(1,j)))
   60 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 70 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffct(k,1))*aimag(ffct(k,1))
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(aimag(qt(k,1)),-real(qt(k,1)))
      fxyt(k,1,1) = zero
      fxyt(k,2,1) = at3*zt1
      fxyt(k,3,1) = zero
      fxyt(k1,1,1) = zero
      fxyt(k1,2,1) = at3*conjg(zt1)
      fxyt(k1,3,1) = zero
      fxyt(k,1,j1) = zero
      fxyt(k,2,j1) = zero
      fxyt(k,3,j1) = zero
      fxyt(k1,1,j1) = zero
      fxyt(k1,2,j1) = zero
      fxyt(k1,3,j1) = zero
      wp = wp + at1*(qt(k,1)*conjg(qt(k,1)))
   70 continue
      k1 = nyh + 1
      fxyt(1,1,1) = zero
      fxyt(1,2,1) = zero
      fxyt(1,3,1) = zero
      fxyt(k1,1,1) = zero
      fxyt(k1,2,1) = zero
      fxyt(k1,3,1) = zero
      fxyt(1,1,j1) = zero
      fxyt(1,2,j1) = zero
      fxyt(1,3,j1) = zero
      fxyt(k1,1,j1) = zero
      fxyt(k1,2,j1) = zero
      fxyt(k1,3,j1) = zero
      we = real(nx*ny)*wp
      return
      end
      subroutine DIVF2(f,df,nx,ny,ndim,nxvh,nyv)
c this subroutine calculates the divergence in fourier space
c input: all except df, output: df
c approximate flop count is: 16*nxc*nyc + 5*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex f, df, zero, zt1
      dimension f(ndim,nxvh,nyv), df(nxvh,nyv)
      if (ndim.lt.2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the divergence
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt1 = dkx*f(1,j,k) + dky*f(2,j,k)
      df(j,k) = cmplx(-aimag(zt1),real(zt1))
      zt1 = dkx*f(1,j,k1) - dky*f(2,j,k1)
      df(j,k1) = cmplx(-aimag(zt1),real(zt1))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      df(1,k) = dky*cmplx(-aimag(f(2,1,k)),real(f(2,1,k)))
      df(1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      df(j,1) = dkx*cmplx(-aimag(f(1,j,1)),real(f(1,j,1)))
      df(j,k1) = zero
   40 continue
      df(1,1) = zero
      df(1,k1) = zero
      return
      end
      subroutine GRADF2(df,f,nx,ny,ndim,nxvh,nyv)
c this subroutine calculates the gradient in fourier space
c input: all except f, output: f
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex df, f, zero, zt1
      dimension df(nxvh,nyv), f(ndim,nxvh,nyv)
      if (ndim.lt.2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the gradient
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt1 = cmplx(-aimag(df(j,k)),real(df(j,k)))
      f(1,j,k) = dkx*zt1
      f(2,j,k) = dky*zt1
      zt1 = cmplx(-aimag(df(j,k1)),real(df(j,k1)))
      f(1,j,k1) = dkx*zt1
      f(2,j,k1) = -dky*zt1
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      f(1,1,k) = zero
      f(2,1,k) = dky*cmplx(-aimag(df(1,k)),real(df(1,k)))
      f(1,1,k1) = zero
      f(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      f(1,j,1) = dkx*cmplx(-aimag(df(j,1)),real(df(j,1)))
      f(2,j,1) = zero
      f(1,j,k1) = zero
      f(2,j,k1) = zero
   40 continue
      f(1,1,1) = zero
      f(2,1,1) = zero
      f(1,1,k1) = zero
      f(2,1,k1) = zero
      if (ndim.eq.2) return
c handle case of ndim = 3
      do 60 k = 2, nyh
      k1 = ny2 - k
      do 50 j = 2, nxh
      f(3,j,k) = zero
      f(3,j,k1) = zero
   50 continue
   60 continue
      k1 = nyh + 1
      do 70 j = 2, nxh
      f(3,j,1) = zero
      f(3,j,k1) = zero
   70 continue
      f(3,1,1) = zero
      f(3,1,k1) = zero
      return
      end
      subroutine CURLF2(f,g,nx,ny,nxvh,nyv)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 32*nxc*nyc + 10*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex f, g, zero, zt1, zt2, zt3
      dimension f(3,nxvh,nyv), g(3,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the curl
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt1 = cmplx(-aimag(f(3,j,k)),real(f(3,j,k)))
      zt2 = cmplx(-aimag(f(2,j,k)),real(f(2,j,k)))
      zt3 = cmplx(-aimag(f(1,j,k)),real(f(1,j,k)))
      g(1,j,k) = dky*zt1
      g(2,j,k) = -dkx*zt1
      g(3,j,k) = dkx*zt2 - dky*zt3
      zt1 = cmplx(-aimag(f(3,j,k1)),real(f(3,j,k1)))
      zt2 = cmplx(-aimag(f(2,j,k1)),real(f(2,j,k1)))
      zt3 = cmplx(-aimag(f(1,j,k1)),real(f(1,j,k1)))
      g(1,j,k1) = -dky*zt1
      g(2,j,k1) = -dkx*zt1
      g(3,j,k1) = dkx*zt2 + dky*zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      zt1 = cmplx(-aimag(f(3,1,k)),real(f(3,1,k)))
      zt3 = cmplx(-aimag(f(1,1,k)),real(f(1,1,k)))
      g(1,1,k) = dky*zt1
      g(2,1,k) = zero
      g(3,1,k) = -dky*zt3
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      g(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt1 = cmplx(-aimag(f(3,j,1)),real(f(3,j,1)))
      zt2 = cmplx(-aimag(f(2,j,1)),real(f(2,j,1)))
      g(1,j,1) = zero
      g(2,j,1) = -dkx*zt1
      g(3,j,1) = dkx*zt2
      g(1,j,k1) = zero
      g(2,j,k1) = zero
      g(3,j,k1) = zero
   40 continue
      g(1,1,1) = zero
      g(2,1,1) = zero
      g(3,1,1) = zero
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      g(3,1,k1) = zero
      return
      end
      subroutine CURLF22(f,g,nx,ny,nxvh,nyv)
c this subroutine calculates the curl in fourier space
c input: all except g, output: g
c approximate flop count is: 18*nxc*nyc + 6*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c g(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex f, g, zero, zt2, zt3
      dimension f(2,nxvh,nyv), g(nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the curl
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt2 = cmplx(-aimag(f(2,j,k)),real(f(2,j,k)))
      zt3 = cmplx(-aimag(f(1,j,k)),real(f(1,j,k)))
      g(j,k) = dkx*zt2 - dky*zt3
      zt2 = cmplx(-aimag(f(2,j,k1)),real(f(2,j,k1)))
      zt3 = cmplx(-aimag(f(1,j,k1)),real(f(1,j,k1)))
      g(j,k1) = dkx*zt2 + dky*zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      zt3 = cmplx(-aimag(f(1,1,k)),real(f(1,1,k)))
      g(1,k) = -dky*zt3
      g(1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt2 = cmplx(-aimag(f(2,j,1)),real(f(2,j,1)))
      g(j,1) = dkx*zt2
      g(j,k1) = zero
   40 continue
      g(1,1) = zero
      g(1,k1) = zero
      return
      end
      subroutine LAPLACE23(f,g,nx,ny,nxvh,nyv)
c this subroutine calculates the vector laplacian in fourier space
c input: all except g, output: g
c approximate flop count is: 16*nxc*nyc + 9*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the laplacian is calculated using the equations:
c gx(kx,ky) = -(kx*kx+ky*ky)*fx(kx,ky)
c gy(kx,ky) = -(kx*kx+ky*ky)*fy(kx,ky)
c gz(kx,ky) = -(kx*kx+ky*ky)*fz(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex f, g, zero
      dimension f(3,nxvh,nyv), g(3,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the laplacian
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      dk2 = dkx*dkx + dky2
      g(1,j,k) = -dk2*f(1,j,k)
      g(2,j,k) = -dk2*f(2,j,k)
      g(3,j,k) = -dk2*f(3,j,k)
      g(1,j,k1) = -dk2*f(1,j,k1)
      g(2,j,k1) = -dk2*f(2,j,k1)
      g(3,j,k1) = -dk2*f(3,j,k1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dk2 = dky*dky
      g(1,1,k) = -dk2*f(1,1,k) 
      g(2,1,k) = -dk2*f(2,1,k)
      g(3,1,k) = -dk2*f(3,1,k)
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      g(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      dk2 = dkx*dkx
      g(1,j,1) = -dk2*f(1,j,1)
      g(2,j,1) = -dk2*f(2,j,1)
      g(3,j,1) = -dk2*f(3,j,1)
      g(1,j,k1) = zero
      g(2,j,k1) = zero
      g(3,j,k1) = zero
   40 continue
      g(1,1,1) = zero
      g(2,1,1) = zero
      g(3,1,1) = zero
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      g(3,1,k1) = zero
      return
      end
      subroutine LAPLACE22(f,g,nx,ny,nxvh,nyv)
c this subroutine calculates the vector laplacian in fourier space
c input: all except g, output: g
c approximate flop count is: 12*nxc*nyc + 7*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the laplacian is calculated using the equations:
c gx(kx,ky) = -(kx*kx+ky*ky)*fx(kx,ky)
c gy(kx,ky) = -(kx*kx+ky*ky)*fy(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0.
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex f, g, zero
      dimension f(2,nxvh,nyv), g(2,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the laplacian
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      dk2 = dkx*dkx + dky2
      g(1,j,k) = -dk2*f(1,j,k)
      g(2,j,k) = -dk2*f(2,j,k)
      g(1,j,k1) = -dk2*f(1,j,k1)
      g(2,j,k1) = -dk2*f(2,j,k1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dk2 = dky*dky
      g(1,1,k) = -dk2*f(1,1,k) 
      g(2,1,k) = -dk2*f(2,1,k)
      g(1,1,k1) = zero
      g(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      dk2 = dkx*dkx
      g(1,j,1) = -dk2*f(1,j,1)
      g(2,j,1) = -dk2*f(2,j,1)
      g(1,j,k1) = zero
      g(2,j,k1) = zero
   40 continue
      g(1,1,1) = zero
      g(2,1,1) = zero
      g(1,1,k1) = zero
      g(2,1,k1) = zero
      return
      end
      subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nxh
c nyv = second dimension of current array, must be >= ny
      complex cu, zero, zt1
      dimension cu(3,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
      zt1 = at1*(dkx*cu(1,j,k1) - dky*cu(2,j,k1))
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt1
      cu(2,j,k1) = cu(2,j,k1) + dky*zt1
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   40 continue
      !cu(1,1,1) = zero
      !cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
      subroutine MODECHOP2(fxy,kxcutoff,kycutoff,nx,ny,nxvh,nyv)
!*****
!This subroutine cuts off kx and ky above kxcutoff and kycutoff
!for debugging by setting the modes to zero. The blocks of code
!are labled according to fft2Dpacking.ods
!*****
      complex fxy, zero, zt1
      dimension fxy(3,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      jcutoff = kxcutoff
      kcutoff = kycutoff
c mode numbers kxcutoff < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = jcutoff, nxh
      fxy(1,j,k) = zero
      fxy(2,j,k) = zero
      fxy(3,j,k) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
   10 continue
   20 continue
c mode numbers 0 < kx < nx/2 and kycutoff < ky < ny/2
      do 21 k = kcutoff, nyh
      k1 = ny2 - k
      do 11 j = 2, nxh
      fxy(1,j,k) = zero
      fxy(2,j,k) = zero
      fxy(3,j,k) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
   11 continue
   21 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = kcutoff, nyh
      k1 = ny2 - k
!block a
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(3,1,k) = zero
!block i
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = jcutoff, nxh
!block e
      fxy(1,j,1) = zero
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
!block k
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
   40 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      return
      end
      subroutine PARFORCE2(fxy,fp,nx,ny,nxvh,nyv)
!*****
!This subroutine calculates the parallel component of fxy with
!the smae method as CUPERP2, except it keeps the parallel part
!instead of subtracting it from the total field. The blocks are
!labeled according to fft2Dpacking.ods.
!*****
      complex fxy, fp, zero, zt1
      dimension fxy(3,nxvh,nyv), fp(2,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate parallel part of force (blocks b and h) 
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*fxy(1,j,k) + dky*fxy(2,j,k))
      fp(1,j,k) = dkx*zt1
      fp(2,j,k) = dky*zt1
      zt1 = at1*(dkx*fxy(1,j,k1) - dky*fxy(2,j,k1))
      fp(1,j,k1) = dkx*zt1
      fp(2,j,k1) = - dky*zt1
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
!block a
      at1 = 1./dky2
      zt1 = at1*dky*fxy(2,1,k)
      fp(1,1,k) = zero
      fp(2,1,k) = dky*zt1
!block i
      dkx = dnx*float(nx/2)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*fxy(1,1,k1) - dky*fxy(2,1,k1))
      fp(1,1,k1) = dkx*zt1
      fp(2,1,k1) = - dky*zt1
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
!block e
      at1 = 1./(dkx*dkx)
      zt1 = at1*dkx*fxy(1,j,1)
      fp(1,j,1) = dkx*zt1
      fp(2,j,1) = zero
!block k
      dky = dny*float(ny/2)
      dky2 = dky*dky
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*fxy(1,j,k1) - dky*fxy(2,j,k1))
      fp(1,j,k1) = dkx*zt1
      fp(2,j,k1) = - dky*zt1
   40 continue
      fp(1,1,1) = zero
      fp(2,1,1) = zero
      
      fp(1,1,k1) = zero
      fp(2,1,k1) = zero

!      print*, 'fxF =', fxy(1,1:nxh,1)
!      print*, 'fxF =', fxy(1,1:nxh,2)
!      print*, 'fxF =', fxy(1,1:nxh,3)
!      print*, 'fxF =', fxy(1,1:nxh,4)

      print*, 'fpxF =', fp(1,1:nxh,1)
      print*, 'fpxF =', fp(1,1:nxh,2)
      print*, 'fpxF =', fp(1,1:nxh,3)
      print*, 'fpxF =', fp(1,1:nxh,4)

      print*, 'fpyF =', fp(2,1:nxh,1)
      print*, 'fpyF =', fp(2,1:nxh,2)
      print*, 'fpyF =', fp(2,1:nxh,3)
      print*, 'fpyF =', fp(2,1:nxh,4)

      return
      end
      subroutine CUSHIFT2(cu,nx,ny,nxvh,nyv)
!*****
!This subroutine shifts each foruier component of the charge
!conserving current from it's position on the yee cell to the
!integer grid points, where the charge is. It also multiplies
!the ccc by [k]/k so that is satisfies the finite difference
!continuity equation.
!*****
      complex cu, zero, zt1, phix, phiy
      dimension cu(3,nxvh,nyv)
!      real thetax, thetay
!      integer yee
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)

! mode numbers 0 < kx < nx/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      phiy = cmplx(cos(thetay),sin(thetay))
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      phix = cmplx(cos(thetax),sin(thetax))
! 0 < ky < ny/2 (block b)
      cu(1,j,k) = cu(1,j,k)*conjg(phix)
      cu(1,j,k) = cu(1,j,k)*sin(thetax)/thetax
      cu(2,j,k) = cu(2,j,k)*conjg(phiy)
      cu(2,j,k) = cu(2,j,k)*sin(thetay)/thetay
! -ny/2 <= ky < 0 (block h)
      cu(1,j,k1) = cu(1,j,k1)*conjg(phix)
      cu(1,j,k1) = cu(1,j,k1)*sin(thetax)/thetax
      cu(2,j,k1) = cu(2,j,k1)*phiy
      cu(2,j,k1) = cu(2,j,k1)*sin(thetay)/thetay
   10 continue
   20 continue
c mode numbers kx = 0
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      phiy = cmplx(cos(thetay),sin(thetay))
! 0 < ky < ny/2 (block a)
      cu(1,1,k) = cu(1,1,k)
      cu(2,1,k) = cu(2,1,k)*conjg(phiy)
      cu(2,1,k) = cu(2,1,k)*sin(thetay)/thetay
! -ny/2 <= ky < 0 (block i)
      cu(1,1,k1) = cu(1,1,k1)
      cu(2,1,k1) = cu(2,1,k1)*phiy
      cu(2,1,k1) = cu(2,1,k1)*sin(thetay)/thetay
   30 continue
c mode numbers ! 0 < x < nx/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      dky = dny*float(k1 - 1)
      phix = cmplx(cos(thetax),sin(thetax))
      thetay = dky/2.0
      phiy = cmplx(cos(thetay),sin(thetay))
! ky = 0 (block e)
      cu(1,j,1) = cu(1,j,1)*conjg(phix)
      cu(1,j,1) = cu(1,j,1)*sin(thetax)/thetax
      cu(2,j,1) = cu(2,j,1)
      cu(2,j,1) = cu(2,j,1)
! ky = -ny/2 (block k)
      cu(1,j,k1) = cu(1,j,k1)*conjg(phix)
      cu(1,j,k1) = cu(1,j,k1)*sin(thetax)/thetax
      cu(2,1,k1) = cu(2,1,k1)*phiy
      cu(2,1,k1) = cu(2,1,k1)*sin(thetay)/thetay
   40 continue
!ky = 0; kx = 0, nx/2
      cu(1,1,1) = cmplx(real(cu(1,1,1)),0.0)! d,f
      cu(2,1,1) = cu(2,1,1)!d,f
!ky = -ny/2; kx = 0, nx/2
      cu(1,1,k1) = cmplx(real(cu(1,1,k1)),0.0)! j,l
      cu(2,1,k1) = zero!j,l
      return
      end
      subroutine SAVCUZ2(cu,cu0,nxvh,nyv)
c save initial net current values for k=0 electric field calculation
c cu(i,j) = complex current density for component i, fourier mode (j-1)
c cu0(1:3) = net initial current
c cux0 = net current in x direction
c nxvh = first dimension of field arrays, must be >= nx/2
      integer nxvh, nyv
      real cu0
      dimension cu0(3)
      complex cu
      dimension cu(3,nxvh,nyv)
c save initial net current
      cu0(1) = real(cu(1,1,1))
      cu0(2) = real(cu(2,1,1))
      cu0(3) = real(cu(3,1,1))
      return
      end
      subroutine MAXZ2(exy,cu,cu0,ffc,dt,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine adds k=0 mode to the solution of maxwell's equation in
c fourier space for electric fields with periodic boundary conditions.
c input: all, output: wf, fx, eyz
c the electric field is updated a whole step using the equations:
c fx(0) = fx(0) - affp*dt*(cux0-cu0(1))*s(0)
c ey(0) = ey(0) - affp*dt*(cuy(0)-cu0(2))*s(0)
c ez(0) = ez(0) - affp*dt*(cuz(0)-cu0(3))*s(0)
c and s(0) = 1.0
c fx(j) = complex force/charge for fourier mode j-1
c eyz(i,j) = complex transverse electric field
c cu(i,j) = complex current density
c for component i, all for fourier mode (j-1)
c cu0(1:3) = net initial current
c cux0 = net current in x direction
c real(ffc(1)) = affp = normalization constant = nx/np,
c where np=number of particles
c dt = time interval between successive calculations
c transverse electric field energy is incremented, using
c wf = wf + 0.5*nx*sum((1/affp)*|eyz(0)|**2)
c nx = system length in x direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxhd = first dimension of form factor array, must be >= nx/2
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real cu0, dt, wf
      dimension cu0(3)
      complex exy, cu, ffc
      dimension exy(3,nxvh,nyv), cu(3,nxvh,nyv), ffc(nxhd,nyhd)
c local data
      real affp, anorm, afdt
      complex zt1, zt2, zt3
      double precision ws
      affp = real(ffc(1,1))
      anorm = 1.0/affp
      afdt = affp*dt*aimag(ffc(1,1))
      !print *,afdt,affp,dt
      
c update electromagnetic field and sum field energies
      ws = 0.0d0
c calculate the electromagnetic fields for mode number kx = 0

      zt1 = exy(1,1,1) - afdt*(real(cu(1,1,1)) - cu0(1))
      zt2 = exy(2,1,1) - afdt*(real(cu(2,1,1)) - cu0(2))
      zt3 = exy(3,1,1) - afdt*(real(cu(3,1,1)) - cu0(3))
      ws = ws + anorm*(zt1*conjg(zt1) + zt2*conjg(zt2) +zt3*conjg(zt3))
      exy(1,1,1) = zt1
      exy(2,1,1) = zt2
      exy(3,1,1) = zt3
      
      wf = wf + 0.5*real(nx*ny)*ws
      return
      end
      subroutine CUPERP22(cu,nx,ny,nxvh,nyv)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cu
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nxh
c nyv = second dimension of current array, must be >= ny
      complex cu, zero, zt1
      dimension cu(2,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
      zt1 = at1*(dkx*cu(1,j,k1) - dky*cu(2,j,k1))
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt1
      cu(2,j,k1) = cu(2,j,k1) + dky*zt1
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   40 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
      subroutine CUPERP2T(cut,nx,ny,nxvh,nyv)
c this subroutine calculates the transverse current in fourier space
c input: all, output: cut
c approximate flop count is: 36*nxc*nyc
c and nxc*nyc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
c and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
c cut(k,i,j) = complex current density for fourier mode (k-1,j-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nxh
c nyv = second dimension of current array, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex cut
      dimension cut(nyv,3,nxvh)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      real dnx, dny, dkx, dkx2, dky, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
c calculate transverse part of current
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 2, nxh
      dkx = dnx*real(j - 1)
      dkx2 = dkx*dkx
      do 10 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      at1 = 1.0/(dkx2 + dky*dky)
      zt1 = at1*(dkx*cut(k,1,j) + dky*cut(k,2,j))
      cut(k,1,j) = cut(k,1,j) - dkx*zt1
      cut(k,2,j) = cut(k,2,j) - dky*zt1
      zt1 = at1*(dkx*cut(k1,1,j) - dky*cut(k1,2,j))
      cut(k1,1,j) = cut(k1,1,j) - dkx*zt1
      cut(k1,2,j) = cut(k1,2,j) + dky*zt1
   10 continue
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      cut(1,1,j) = zero
      cut(k1,1,j) = zero
      cut(k1,2,j) = zero
   30 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 40 k = 2, nyh
      k1 = ny2 - k
      cut(k,2,1) = zero
      cut(k1,1,1) = conjg(cut(k,1,1))
      cut(k1,2,1) = zero
      cut(k,1,j1) = zero
      cut(k,2,j1) = zero
      cut(k1,1,j1) = zero
      cut(k1,2,j1) = zero
   40 continue
      k1 = nyh + 1
      cut(1,1,1) = zero
      cut(1,2,1) = zero
      cut(k1,1,1) = zero
      cut(k1,2,1) = zero
      cut(1,1,j1) = zero
      cut(1,2,j1) = zero
      cut(k1,1,j1) = zero
      cut(k1,2,j1) = zero
      return
      end
      subroutine BPOIS23(cu,bxy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,nxvh,ny
     1v,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd, output: bxy,wm
c approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
c for isign = 1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd, output: bxy,wm
c approximate flop count is: 68*nxc*nyc + 21*(nxc + nyc)
c for isign = 2, input: cu,ffc,isign,nx,ny,nxvh,nyhd, output: bxy
c approximate flop count is: 12*nxc*nyc + 6*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
c = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(1,j,k) = x component of complex magnetic field
c bxy(2,j,k) = y component of complex magnetic field
c bxy(3,j,k) = z component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffc, zero, zt1, zt2, zt3
      dimension cu(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      bxy(1,j,k1) = -at3*zt1
      bxy(2,j,k1) = -at2*zt1
      bxy(3,j,k1) = at2*zt2 + at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
     2 cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bxy(1,j,1) = zero
      bxy(2,j,1) = -at2*zt1
      bxy(3,j,1) = at2*zt2
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
   70 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at2 = ci2*real(ffc(j,k))
      at1 = at2*aimag(ffc(j,k))
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(2,j,k) = at2*cu(2,j,k)
      bxy(3,j,k) = at2*cu(3,j,k)
      bxy(1,j,k1) = at2*cu(1,j,k1)
      bxy(2,j,k1) = at2*cu(2,j,k1)
      bxy(3,j,k1) = at2*cu(3,j,k1)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
     2 cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = ci2*real(ffc(1,k))
      at1 = at2*aimag(ffc(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = at2*cu(2,1,k)
      bxy(3,1,k) = at2*cu(3,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at2 = ci2*real(ffc(j,1))
      at1 = at2*aimag(ffc(j,1))
      bxy(1,j,1) = at2*cu(1,j,1)
      bxy(2,j,1) = at2*cu(2,j,1)
      bxy(3,j,1) = at2*cu(3,j,1)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
  120 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      k1 = ny2 - k
      do 140 j = 2, nxh
      at1 = aimag(ffc(j,k))
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(2,j,k) = at1*cu(2,j,k)
      bxy(3,j,k) = at1*cu(3,j,k)
      bxy(1,j,k1) = at1*cu(1,j,k1)
      bxy(2,j,k1) = at1*cu(2,j,k1)
      bxy(3,j,k1) = at1*cu(3,j,k1)
  140 continue
  150 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffc(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = at1*cu(2,1,k)
      bxy(3,1,k) = at1*cu(3,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nxh
      at1 = aimag(ffc(j,1))
      bxy(1,j,1) = at1*cu(1,j,1)
      bxy(2,j,1) = at1*cu(2,j,1)
      bxy(3,j,1) = at1*cu(3,j,1)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
  170 continue
      at1 = aimag(ffc(1,1))
      bxy(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
      bxy(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
      bxy(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      return
      end
      subroutine BPOIS22(cu,bxy,bz,isign,ffc,ax,ay,affp,ci,wm,nx,ny,nxvh
     1,nyv,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd, output: bz,wm
c approximate flop count is: 60*nxc*nyc + 30*(nxc + nyc)
c for isign = 1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd, output: bxy,wm
c approximate flop count is: 50*nxc*nyc + 24*(nxc + nyc)
c for isign = 2, input: cu,ffc,isign,nx,ny,nxvh,nyhd, output: bxy
c approximate flop count is: 8*nxc*nyc + 4*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bz(kx=pi) = bz(ky=pi) = 0, and bz(kx=0,ky=0) = 0.
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(1,j,k) = x component of complex vector potential
c bxy(2,j,k) = y component of complex vector potential
c bz(j,k) = z component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, bz, ffc, zero, zt2, zt3
      dimension cu(2,nxvh,nyv), bxy(2,nxvh,nyv), bz(nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bz(j,k) = at2*zt2 - at3*zt3
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      bz(j,k1) = at2*zt2 + at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))
     2)
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bz(1,k) = -at3*zt3
      bz(1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bz(j,1) = at2*zt2
      bz(j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)))
   70 continue
      bz(1,1) = zero
      bz(1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at2 = ci2*real(ffc(j,k))
      at1 = at2*aimag(ffc(j,k))
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(2,j,k) = at2*cu(2,j,k)
      bxy(1,j,k1) = at2*cu(1,j,k1)
      bxy(2,j,k1) = at2*cu(2,j,k1)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))
     2)
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = ci2*real(ffc(1,k))
      at1 = at2*aimag(ffc(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = at2*cu(2,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at2 = ci2*real(ffc(j,1))
      at1 = at2*aimag(ffc(j,1))
      bxy(1,j,1) = at2*cu(1,j,1)
      bxy(2,j,1) = at2*cu(2,j,1)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)))
  120 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      k1 = ny2 - k
      do 140 j = 2, nxh
      at1 = aimag(ffc(j,k))
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(2,j,k) = at1*cu(2,j,k)
      bxy(1,j,k1) = at1*cu(1,j,k1)
      bxy(2,j,k1) = at1*cu(2,j,k1)
  140 continue
  150 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffc(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = at1*cu(2,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nxh
      at1 = aimag(ffc(j,1))
      bxy(1,j,1) = at1*cu(1,j,1)
      bxy(2,j,1) = at1*cu(2,j,1)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
  170 continue
      at1 = aimag(ffc(1,1))
      bxy(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
      bxy(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      return
      end
      subroutine IBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd,yee)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions.
c input: cu,ffc,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
c = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
!*****
!Checks to see if yee = 1. If it does, then a separate block
!of code is executed that accounts for the phase shift between
!the charge conserving current and the components of B on a yee cube.
!*****
      double precision wp
      complex cu, bxy, ffc, zero, zt1, zt2, zt3
      integer yee
      real dkx, dky, thetax, thetay, kxop, kyop
      complex phix, phiy
      dimension cu(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0

!      if(yee==1) then
!!The Darwin B is calculated from Ampere's law with no E.
!!We need to modify what viktor has so that B = -affp*[k]/[k]^2 \cross J
!!And J comes from OSIRIS, so it's at the grid points of E
!!So we need to shift B to its grid points on the Yee mesh too

!c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!      do 21 k = 2, nyh
!      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      thetay = dky/2.0
!      kyop = 2.0*sin(thetay)
!      phiy = cmplx(cos(thetay),sin(thetay))
!      do 11 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      thetax = dkx/2.0
!      kxop = 2.0*sin(thetax)
!      phix = cmplx(cos(thetax),sin(thetax))
!      at1 = ci2*real(ffc(j,k))
!      !need at1 = c2*affp/[k]^2; [k] = [kx] + [ky]
!      at1 = at1*(dkx*dkx + dky*dky)/(kxop*kxop + kyop*kyop)
!      at2 = kxop*at1
!      at3 = kyop*at1
!      at1 = at1*aimag(ffc(j,k))
!      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
!      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
!      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
!      bxy(1,j,k) = at3*zt1*phiy
!      bxy(2,j,k) = -at2*zt1*phix
!      bxy(3,j,k) = at2*zt2*phix - at3*zt3*phiy
!      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
!      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
!      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
!      bxy(1,j,k1) = -at3*zt1*conjg(phiy)
!      bxy(2,j,k1) = -at2*zt1*phix
!      bxy(3,j,k1) = at2*zt2*phix + at3*zt3*conjg(phiy)
!      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
!     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
!     2 cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
!   11 continue
!   21 continue
!c mode numbers kx = 0, nx/2
!cdir$ ivdep
!      do 31 k = 2, nyh
!      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      thetay = dky/2.0
!      kyop = 2.0*sin(thetay)
!      phiy = cmplx(cos(thetay),sin(thetay))
!      at1 = ci2*real(ffc(1,k))
!      at1 = at1*dky*dky/(kyop*kyop)!dkx = kyop = 0 only
!      at3 = kyop*at1
!      at1 = at1*aimag(ffc(1,k))
!      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
!      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
!      bxy(1,1,k) = at3*zt1*phiy
!      bxy(2,1,k) = zero
!      bxy(3,1,k) = -at3*zt3*phiy
!      bxy(1,1,k1) = zero
!      bxy(2,1,k1) = zero
!      bxy(3,1,k1) = zero
!      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
!     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
!   31 continue
!c mode numbers ky = 0, ny/2
!      k1 = nyh + 1
!      do 41 j = 2, nxh
!      dkx = dnx*float(nx/2)
!      thetax = dkx/2.0
!      kxop = 2.0*sin(thetax)
!      phix = cmplx(cos(thetax),sin(thetax))
!      at1 = ci2*real(ffc(j,1))
!      at1 = at1*dkx*dkx/(kxop*kxop)!dky = kxop = 0 only
!      at2 = kxop*at1
!      at1 = at1*aimag(ffc(j,1))
!      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
!      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
!      bxy(1,j,1) = zero
!      bxy(2,j,1) = -at2*zt1*phix
!      bxy(3,j,1) = at2*zt2*phix
!      bxy(1,j,k1) = zero
!      bxy(2,j,k1) = zero
!      bxy(3,j,k1) = zero
!      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
!     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
!   41 continue
!      bxy(1,1,1) = zero
!      bxy(2,1,1) = zero
!      bxy(3,1,1) = zero
!      bxy(1,1,k1) = zero
!      bxy(2,1,k1) = zero
!      bxy(3,1,k1) = zero
!      wm = float(nx*ny)*wp

!      else

c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffc(j,k))
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      bxy(1,j,k1) = -at3*zt1
      bxy(2,j,k1) = -at2*zt1
      bxy(3,j,k1) = at2*zt2 + at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
     2 cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      at1 = at1*aimag(ffc(1,k))
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffc(j,1))
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bxy(1,j,1) = zero
      bxy(2,j,1) = -at2*zt1
      bxy(3,j,1) = at2*zt2
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = float(nx*ny)*wp

!      endif

      return
      end
      subroutine IBPOIS23T(cut,bxyt,ffct,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with periodic boundary conditions.
c input: cut,ffct,ci,nx,ny,nxv,nyhd, output: bxyt,wm
c approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
c = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cut(k,i,j) = complex current density for fourier mode (k-1,j-1)
c bxyt(k,i,j) = i component of complex magnetic field
c all for fourier mode (k-1,j-1)
c aimag(ffct(k,j)) = finite-size particle shape factor s
c for fourier mode (k-1,j-1)
c real(ffct(k,j)) = potential green's function g
c for fourier mode (k-1,j-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, wm
      complex cut, bxyt, ffct
      dimension cut(nyv,3,nxvh), bxyt(nyv,3,nxvh)
      dimension ffct(nyhd,nxhd)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      real dnx, dny, dkx, ci2, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 2, nxh
      dkx = dnx*real(j - 1)
      do 10 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffct(k,j))
      at2 = dkx*at1
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffct(k,j))
      zt1 = cmplx(-aimag(cut(k,3,j)),real(cut(k,3,j)))
      zt2 = cmplx(-aimag(cut(k,2,j)),real(cut(k,2,j)))
      zt3 = cmplx(-aimag(cut(k,1,j)),real(cut(k,1,j)))
      bxyt(k,1,j) = at3*zt1
      bxyt(k,2,j) = -at2*zt1
      bxyt(k,3,j) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cut(k1,3,j)),real(cut(k1,3,j)))
      zt2 = cmplx(-aimag(cut(k1,2,j)),real(cut(k1,2,j)))
      zt3 = cmplx(-aimag(cut(k1,1,j)),real(cut(k1,1,j)))
      bxyt(k1,1,j) = -at3*zt1
      bxyt(k1,2,j) = -at2*zt1
      bxyt(k1,3,j) = at2*zt2 + at3*zt3
      wp = wp + at1*(cut(k,1,j)*conjg(cut(k,1,j)) + cut(k,2,j)*conjg(cut
     1(k,2,j)) + cut(k,3,j)*conjg(cut(k,3,j)) + cut(k1,1,j)*conjg(cut(k1
     2,1,j)) + cut(k1,2,j)*conjg(cut(k1,2,j)) + cut(k1,3,j)*conjg(cut(k1
     3,3,j)))
   10 continue
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      at1 = ci2*real(ffct(1,j))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffct(1,j))
      zt1 = cmplx(-aimag(cut(1,3,j)),real(cut(1,3,j)))
      zt2 = cmplx(-aimag(cut(1,2,j)),real(cut(1,2,j)))
      bxyt(1,1,j) = zero
      bxyt(1,2,j) = -at2*zt1
      bxyt(1,3,j) = at2*zt2
      bxyt(k1,1,j) = zero
      bxyt(k1,2,j) = zero
      bxyt(k1,3,j) = zero
      wp = wp + at1*(cut(1,1,j)*conjg(cut(1,1,j)) + cut(1,2,j)*conjg(cut
     1(1,2,j)) + cut(1,3,j)*conjg(cut(1,3,j)))
   30 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 40 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffct(k,1))
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffct(k,1))
      zt1 = cmplx(-aimag(cut(k,3,1)),real(cut(k,3,1)))
      zt3 = cmplx(-aimag(cut(k,1,1)),real(cut(k,1,1)))
      bxyt(k,1,1) = at3*zt1
      bxyt(k,2,1) = zero
      bxyt(k,3,1) = -at3*zt3
      bxyt(k1,1,1) = at3*conjg(zt1)
      bxyt(k1,2,1) = zero
      bxyt(k1,3,1) = -at3*conjg(zt3)
      bxyt(k,1,j1) = zero
      bxyt(k,2,j1) = zero
      bxyt(k,3,j1) = zero
      bxyt(k1,1,j1) = zero
      bxyt(k1,2,j1) = zero
      bxyt(k1,3,j1) = zero
      wp = wp + at1*(cut(k,1,1)*conjg(cut(k,1,1)) + cut(k,2,1)*conjg(cut
     1(k,2,1)) + cut(k,3,1)*conjg(cut(k,3,1)))
   40 continue
      k1 = nyh + 1
      bxyt(1,1,1) = zero
      bxyt(1,2,1) = zero
      bxyt(1,3,1) = zero
      bxyt(k1,1,1) = zero
      bxyt(k1,2,1) = zero
      bxyt(k1,3,1) = zero
      bxyt(1,1,j1) = zero
      bxyt(1,2,j1) = zero
      bxyt(1,3,j1) = zero
      bxyt(k1,1,j1) = zero
      bxyt(k1,2,j1) = zero
      bxyt(k1,3,j1) = zero
      wm = real(nx*ny)*wp
      return
      end
      subroutine MAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nxhd,
     1nyhd,parmax,yee)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh

!*****
!Now tests for if yee = 1. If it does, the solver executes a block
!of code that is based off of what it did before. This new block has 
!the relative phase shifts between the fields on the yee mesh in the 
!update parts of the equations being used to calculate the fields
!below. These equations also have [k] (krop) replacing k (dkr) at the
!needed places. The nx/2, ny/2 and 0 modes are also claclulated. The 
!sub blocks are labeled according to fft2Dpacking.ods
!If yee = 0, but we still want to use the ccc, then we calculate
!kx = 0, ky = [-ny/2+1, ny/2-1], and kx = [-nx/2+1, nx/2-1], ky = 0.
!*****
      double precision wp, ws
      complex exy, bxy, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      complex phix, phiy
      complex e1T, e2T, e3T, b1T, b2T, b3T, cu1T, cu2T, cu3T
      integer parmax, yee
      real kxop, kyop, thetax, thetay
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0

      !print*, exy(3,:,1)
      if(yee==1) then

      do 25 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 15 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      thetay = dky/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      phiy = cmplx(cos(thetay),sin(thetay))
      kxop = 2.0*sin(thetax)
      kyop = 2.0*sin(thetay)      
      afdt = adt*aimag(ffc(j,k))
c update magnetic field half time step, ky > 0 (block b)
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(kyop*zt1*phiy)
      zt5 = bxy(2,j,k) + dth*(kxop*zt1*phix)
      !print*, 'before', bxy(3,j,k)
      !print*, 'delta1 = ', kxop*zt2*phix, 'delta2 =', kyop*zt3*phiy
      zt6 = bxy(3,j,k) - dth*(kxop*zt2*phix - kyop*zt3*phiy)
      !print*, 'after', zt6
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(kyop*zt1*conjg(phiy))
     1 - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(kxop*zt2*conjg(phix) - 
     1 kyop*zt3*conjg(phiy)) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(kyop*zt1*phiy)
      zt5 = zt5 + dth*(kxop*zt1*phix)
      zt6 = zt6 - dth*(kxop*zt2*phix - kyop*zt3*phiy)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c update magnetic field half time step, ky < 0 (block h)
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(kyop*zt1*conjg(phiy))
      zt5 = bxy(2,j,k1) + dth*(kxop*zt1*phix)
      zt6 = bxy(3,j,k1) - dth*(kxop*zt2*phix + kyop*zt3*conjg(phiy))
      !print*, 'B1', j, k, zt4, zt5, zt6
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(kyop*zt1*phiy) - afdt*cu(1,j,k1)!*phix
      zt8 = exy(2,j,k1) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(kxop*zt2*conjg(phix) + kyop*zt3*phiy)
     1 - afdt*cu(3,j,k1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      !print*, 'E', j, k, zt7, zt8, zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 + dth*(kyop*zt1*conjg(phiy))
      zt5 = zt5 + dth*(kxop*zt1*phix)
      zt6 = zt6 - dth*(kxop*zt2*phix + kyop*zt3*conjg(phiy))
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      !print*, j, k, bxy(3,j,k), bxy(3,j,k1)
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   15 continue
   25 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 35 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      phiy = cmplx(cos(thetay),sin(thetay))
      kyop = 2.0*sin(thetay) 
      afdt = adt*aimag(ffc(1,k))
c kx = 0 
c update magnetic field half time step, ky > 0 (block a)
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(kyop*zt1*phiy)
      zt5 = bxy(2,1,k)
      zt6 = bxy(3,1,k) + dth*(kyop*zt3*phiy)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(kyop*zt1*conjg(phiy)) - afdt*cu(1,1,k)
      if(parmax==1) then
      zt8 = exy(2,1,k) - afdt*cu(2,1,k)
      else
      zt8 = zero
      endif
      zt9 = exy(3,1,k) - cdt*(kyop*zt3*conjg(phiy)) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,1,k) = zt7
      exy(2,1,k) = zt8
      exy(3,1,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
     1 zt9*conjg(zt9))
      zt4 = zt4 - dth*(kyop*zt1*phiy)
      zt6 = zt6 + dth*(kyop*zt3*phiy)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zt5
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + 
     1 zt6*conjg(zt6))
!*****
!We cannot calculate these correctly either, because of how the fft data for
!the nx/2 mode is stored 
!*****
!c kx = nx/2
!c update magnetic field half time step, ky < 0 (block i)
!      dkx = dnx*float(nx/2)
!      thetax = dkx/2.0
!      phix = cmplx(cos(thetax),sin(thetax))
!      kxop = 2.0*sin(thetax)       
!      zt1 = cmplx(-aimag(exy(3,1,k1)),real(exy(3,1,k1)))
!      zt2 = cmplx(-aimag(exy(2,1,k1)),real(exy(2,1,k1)))
!      zt3 = cmplx(-aimag(exy(1,1,k1)),real(exy(1,1,k1)))
!      zt4 = bxy(1,1,k1) - dth*(-kyop*zt1*conjg(phiy))
!      zt5 = bxy(2,1,k1) - dth*kxop*zt1*conjg(phix)!
!      zt6 = bxy(3,1,k1) + dth*(-kyop*zt3*conjg(phiy) + 
!     1 kxop*zt5*conjg(phix))!
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      zt7 = exy(1,1,k1) - cdt*(kyop*zt1*phiy) - afdt*cu(1,1,k1)
!      if(parmax==1) then
!      zt8 = exy(2,1,k1) + cdt*kxop*zt1*phix - afdt*cu(2,1,k1)!
!      else
!      zt8 = zero
!      endif
!      zt9 = exy(3,1,k1) + cdt*(-kxop*zt2*phix
!     1 + kyop*zt3*phiy) - afdt*cu(3,1,k1)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      exy(1,1,k1) = zt7
!      exy(2,1,k1) = zt8
!      exy(3,1,k1) = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
!     1 zt9*conjg(zt9))
!      zt4 = bxy(1,1,k1) - dth*(-kyop*zt1*conjg(phiy))
!      zt5 = bxy(2,1,k1) - dth*kxop*zt1*conjg(phix)!
!      zt6 = bxy(3,1,k1) + dth*(-kyop*zt3*conjg(phiy) + 
!     1 kxop*zt2*conjg(phix))!
!      bxy(1,1,k1) = zt4
!      bxy(2,1,k1) = zt5
!      bxy(3,1,k1) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + 
!     1 zt6*conjg(zt6))
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
   35 continue
c mode numbers 0 < kx < nx/2
      k1 = nyh + 1
      do 45 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      kxop = 2.0*sin(thetax)
      afdt = adt*aimag(ffc(j,1))
c update magnetic field half time step, ky = 0 (block e)
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
!      zt3 = cmplx(-aimag(exy(1,j,1)),real(exy(1,j,1)))
      zt4 = bxy(1,j,1)
      zt5 = bxy(2,j,1) + dth*(kxop*zt1*phix)
      zt6 = bxy(3,j,1) - dth*(kxop*zt2*phix)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      if(parmax==1) then
      zt7 = exy(1,j,1) - afdt*cu(1,j,1)!*phix
      else
      zt7 = zero
      endif
      zt8 = exy(2,j,1) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(kxop*zt2*conjg(phix)) - afdt*cu(3,j,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exy(1,j,1) = zt7
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + 
     1 zt9*conjg(zt9))
      zt5 = zt5 + dth*(kxop*zt1*phix)
      zt6 = zt6 - dth*(kxop*zt2*phix)
      bxy(1,j,1) = zt4
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) +
     1 zt6*conjg(zt6))
!*****
!We cannot calculate these correctly either, because of how the fft data for
!the ny/2 mode is stored 
!*****
!c update magnetic field half time step, ky = ny/2 (block k)
!      dky = dnx*float(ny/2)
!      thetay = dky/2.0
!      phiy = cmplx(cos(thetay),sin(thetay))
!      kyop = 2.0*sin(thetay)
!      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
!      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
!      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
!      zt4 = bxy(1,j,k1) + dth*kyop*zt1*conjg(phiy)
!      zt5 = bxy(2,j,k1) + dth*(kxop*zt1*phix)
!      zt6 = bxy(3,j,k1) - dth*(kxop*zt2*phix + kyop*zt1*conjg(phiy))
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      if(parmax==1) then
!      zt7 = exy(1,j,k1) - cdt*kyop*zt1*phiy - afdt*cu(1,j,k1)!*phix
!      else
!      zt7 = zero
!      endif
!      zt8 = exy(2,j,k1) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,k1)
!      zt9 = exy(3,j,k1) + cdt*(kxop*zt2*conjg(phix) + kyop*zt3*phiy)
!     1 - afdt*cu(3,j,k1)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      exy(1,j,k1) = zt7
!      exy(2,j,k1) = zt8
!      exy(3,j,k1) = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + 
!     1 zt9*conjg(zt9))
!      zt4 = bxy(1,j,k1) + dth*kyop*zt1*conjg(phiy)
!      zt5 = bxy(2,j,k1) + dth*(kxop*zt1*phix)
!      zt6 = bxy(3,j,k1) - dth*(kxop*zt2*phix + kyop*zt1*conjg(phiy))
!      bxy(1,j,k1) = zt4
!      bxy(2,j,k1) = zt5
!      bxy(3,j,k1) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) +
!     1 zt6*conjg(zt6))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   45 continue
!*****
!Unfortionatly, when we have to shift for the pusher, later, these modes
!become completely imaginary, we we can't calculate them correctly with
!way the fft data is stored. We comment these lines out, but they can be
!used in a future version that modefies the pusher
!*****

!kx = -nx/2, ky = 0 (block f)*
!      dkx = dnx*float(nx/2)
!      thetax = dkx/2.0
!      phix = cmplx(cos(thetax),sin(thetax))
!      kxop = 2.0*sin(thetax)

!      e3T = cmplx(aimag(exy(3,1,1)),0.0)
!      e2T = cmplx(aimag(exy(2,1,1)),0.0)
!      e1T = cmplx(aimag(exy(1,1,1)),0.0)
!      b3T = cmplx(aimag(bxy(3,1,1)),0.0)
!      b2T = cmplx(aimag(bxy(2,1,1)),0.0)
!      b1T = cmplx(aimag(bxy(1,1,1)),0.0)
!      cu3T = cmplx(aimag(cu(3,1,1)),0.0)
!      cu2T = cmplx(aimag(cu(2,1,1)),0.0)
!      cu1T = cmplx(aimag(cu(1,1,1)),0.0)

!      zt1 = cmplx(-aimag(e3T),real(e3T))
!      zt2 = cmplx(-aimag(e2T),real(e2T))
!      zt3 = cmplx(-aimag(e1T),real(e1T))

!      zt4 = b1T
!      zt5 = b2T - dth*(kxop*zt1*conjg(phix))!phix = i*sin(-pi/2) = -i
!      zt6 = b3T + dth*(kxop*zt2*conjg(phix))

!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      if(parmax==1) then
!      zt7 = e1T - afdt*cu1T
!      else
!      zt7 = zero
!      endif
!      zt8 = e2T + cdt*(kxop*zt1*phix) - afdt*cu2T
!      zt9 = e3T - cdt*(kxop*zt2*phix) - afdt*cu3T
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      e1T = zt7
!      e2T = zt8
!      e3T = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + 
!     1 zt9*conjg(zt9))
!      zt4 = zt4
!      zt5 = zt5 - dth*(kxop*zt1*conjg(phix))
!      zt6 = zt6 + dth*(kxop*zt2*conjg(phix))
!      b1T = zt4
!      b2T = zt5
!      b3T = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) +
!     1 zt6*conjg(zt6))

!      exy(1,1,1) = cmplx(0.0,real(e1T))
!      exy(2,1,1) = cmplx(0.0,real(e2T))
!      exy(3,1,1) = cmplx(0.0,real(e3T))
!      bxy(1,1,1) = cmplx(0.0,real(b1T))
!      bxy(2,1,1) = cmplx(0.0,real(b2T))
!      bxy(3,1,1) = cmplx(0.0,real(b3T))

!!kx = -nx/2, ky = -ny/2 (block l)
!      dky = dny*float(ny/2)
!      thetay = dky/2.0
!      phiy = cmplx(cos(thetay),sin(thetay))
!      kyop = 2.0*sin(thetay)

!      e3T = cmplx(aimag(exy(3,1,k1)),0.0)
!      e2T = cmplx(aimag(exy(2,1,k1)),0.0)
!      e1T = cmplx(aimag(exy(1,1,k1)),0.0)
!      b3T = cmplx(aimag(bxy(3,1,k1)),0.0)
!      b2T = cmplx(aimag(bxy(2,1,k1)),0.0)
!      b1T = cmplx(aimag(bxy(1,1,k1)),0.0)
!      cu3T = cmplx(aimag(cu(3,1,k1)),0.0)
!      cu2T = cmplx(aimag(cu(2,1,k1)),0.0)
!      cu1T = cmplx(aimag(cu(1,1,k1)),0.0)

!      zt1 = cmplx(-aimag(e3T),real(e3T))
!      zt2 = cmplx(-aimag(e2T),real(e2T))
!      zt3 = cmplx(-aimag(e1T),real(e1T))

!      zt4 = b1T + dth*(kyop*zt1*conjg(phiy))
!      zt5 = b2T - dth*(kxop*zt1*conjg(phix))!phix = i*sin(pi/2) = i
!      zt6 = b3T + dth*(kxop*zt2*conjg(phix) - kyop*zt3*conjg(phiy))

!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      if(parmax==1) then
!      zt7 = e1T - cdt*(kyop*zt1*phiy) - afdt*cu1T
!      else
!      zt7 = zero
!      endif
!      zt8 = e2T + cdt*(kxop*zt1*phix) - afdt*cu2T
!      zt9 = e3T - cdt*(kxop*zt2*phix - kyop*zt3*phiy)
!     1 - afdt*cu3T
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      e1T = zt7
!      e2T = zt8
!      e3T = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + 
!     1 zt9*conjg(zt9))
!      zt4 = zt4 + dth*(kyop*zt1*conjg(phiy))
!      zt5 = zt5 - dth*(kxop*zt1*conjg(phix))
!      zt6 = zt6 + dth*(kxop*zt2*conjg(phix) - 
!     1 kyop*zt3*conjg(phiy))
!      b1T = zt4
!      b2T = zt5
!      b3T = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) +
!     1 zt6*conjg(zt6))

!      exy(1,1,k1) = cmplx(real(exy(1,1,k1)),real(e1T))
!      exy(2,1,k1) = cmplx(real(exy(2,1,k1)),real(e2T))
!      exy(3,1,k1) = cmplx(real(exy(3,1,k1)),real(e3T))
!      bxy(1,1,k1) = cmplx(real(bxy(1,1,k1)),real(b1T))
!      bxy(2,1,k1) = cmplx(real(bxy(2,1,k1)),real(b2T))
!      bxy(3,1,k1) = cmplx(real(bxy(3,1,k1)),real(b3T))

!!kx = 0, ky = -ny/2 (block j)
!      dkx = 0.0
!      thetax = dkx/2.0
!      phix = cmplx(1.0,0.0)! cmplx(cos(thetax),sin(thetax))
!      kxop = 2.0*sin(thetax)

!      e3T = cmplx(real(exy(3,1,k1)),0.0)
!      e2T = cmplx(real(exy(2,1,k1)),0.0)
!      e1T = cmplx(real(exy(1,1,k1)),0.0)
!      b3T = cmplx(real(bxy(3,1,k1)),0.0)
!      b2T = cmplx(real(bxy(2,1,k1)),0.0)
!      b1T = cmplx(real(bxy(1,1,k1)),0.0)
!      cu3T = cmplx(real(cu(3,1,k1)),0.0)
!      cu2T = cmplx(real(cu(2,1,k1)),0.0)
!      cu1T = cmplx(real(cu(1,1,k1)),0.0)    

!      zt1 = cmplx(-aimag(e3T),real(e3T))
!      zt2 = cmplx(-aimag(e2T),real(e2T))
!      zt3 = cmplx(-aimag(e1T),real(e1T))

!      zt4 = b1T + dth*(kyop*zt1*conjg(phiy))
!      zt5 = b2T
!      zt6 = b3T + dth*(-kyop*zt3*conjg(phiy))

!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      if(parmax==1) then
!      zt7 = e1T - cdt*(kyop*zt1*phiy) - afdt*cu1T
!      else
!      zt7 = zero
!      endif
!      zt8 = e2T - afdt*cu2T
!      zt9 = e3T - cdt*(-kyop*zt3*phiy) - afdt*cu3T
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      e1T = zt7
!      e2T = zt8
!      e3T = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + 
!     1 zt9*conjg(zt9))
!      zt4 = zt4 + dth*(kyop*zt1*conjg(phiy))
!      zt5 = zt5
!      zt6 = zt6 + dth*(-kyop*zt3*conjg(phiy))
!      b1T = zt4
!      b2T = zt5
!      b3T = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) +
!     1 zt6*conjg(zt6))

!      exy(1,1,k1) = cmplx(real(e1T),aimag(exy(1,1,k1)))
!      exy(2,1,k1) = cmplx(real(e2T),aimag(exy(2,1,k1)))
!      exy(3,1,k1) = cmplx(real(e3T),aimag(exy(3,1,k1)))
!      bxy(1,1,k1) = cmplx(real(b1T),aimag(bxy(1,1,k1)))
!      bxy(2,1,k1) = cmplx(real(b2T),aimag(bxy(2,1,k1)))
!      bxy(3,1,k1) = cmplx(real(b3T),aimag(bxy(3,1,k1)))

      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero

      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp

      else

c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)     
      afdt = adt*aimag(ffc(j,k))
c update magnetic field half time step, ky > 0 (block b)
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
      !print*, 'delta1 = ', dkx*zt2, 'delta2 =', dky*zt3
      zt6 = bxy(3,j,k) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) - cdt*(dkx*zt1) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c update magnetic field half time step, ky < 0 (block h)
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 + dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      !print*, j, k, bxy(3,j,k), bxy(3,j,k1)
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffc(1,k))
!kx = 0
c update magnetic field half time step, ky > 0 (block a)
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
      zt5 = bxy(2,1,k)
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
      if(parmax==1) then
      zt8 = exy(2,1,k) - afdt*cu(2,1,k)
      else
      zt8 = zero
      endif
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,1,k) = zt7
      exy(2,1,k) = zt8
      exy(3,1,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
     1 zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt6 = zt6 + dth*(dky*zt3)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zt5
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + 
     1 zt6*conjg(zt6))
!kx = nx/2
c update magnetic field half time step, ky > 0 (block i)
!      zt1 = cmplx(-aimag(exy(3,1,k1)),real(exy(3,1,k1)))
!      zt3 = cmplx(-aimag(exy(1,1,k1)),real(exy(1,1,k1)))
!      zt4 = bxy(1,1,k1) + dth*(dky*zt1)
!      zt5 = bxy(2,1,k1)
!      zt6 = bxy(3,1,k1) - dth*(dky*zt3)
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      zt7 = exy(1,1,k1) - cdt*(dky*zt1) - afdt*cu(1,1,k1)
!      if(parmax==1) then
!      zt8 = exy(2,1,k1) - afdt*cu(2,1,k1)
!      else
!      zt8 = zero
!      endif
!      zt9 = exy(3,1,k1) + cdt*(dky*zt3) - afdt*cu(3,1,k1)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      exy(1,1,k1) = zt7
!      exy(2,1,k1) = zt8
!      exy(3,1,k1) = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
!     1 zt9*conjg(zt9))
!      zt4 = zt4 + dth*(dky*zt1)
!      zt6 = zt6 - dth*(dky*zt3)
!      bxy(1,1,k1) = zt4
!      bxy(2,1,k1) = zt5
!      bxy(3,1,k1) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + 
!     1 zt6*conjg(zt6))
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
   30 continue
c mode numbers 0 < kx < nx/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffc(j,1))
c update magnetic field half time step, ky = 0 (block e)
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
!      zt3 = cmplx(-aimag(exy(1,j,1)),real(exy(1,j,1)))
      zt4 = bxy(1,j,1)
      zt5 = bxy(2,j,1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,1) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      if(parmax==1) then
      zt7 = exy(1,j,1) - afdt*cu(1,j,1)
      else
      zt7 = zero
      endif
      zt8 = exy(2,j,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exy(1,j,1) = zt7
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + 
     1 zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxy(1,j,1) = zt4
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) +
     1 zt6*conjg(zt6))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero


   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp

!c calculate the electromagnetic fields
!c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!      do 20 k = 2, nyh
!      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      do 10 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      afdt = adt*aimag(ffc(j,k))
!c update magnetic field half time step, ky > 0
!      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
!      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
!      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
!      zt4 = bxy(1,j,k) - dth*(dky*zt1)
!      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
!      zt6 = bxy(3,j,k) - dth*(dkx*zt2 - dky*zt3)
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
!      zt8 = exy(2,j,k) - cdt*(dkx*zt1) - afdt*cu(2,j,k)
!      zt9 = exy(3,j,k) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      exy(1,j,k) = zt7
!      exy(2,j,k) = zt8
!      exy(3,j,k) = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
!      zt4 = zt4 - dth*(dky*zt1)
!      zt5 = zt5 + dth*(dkx*zt1)
!      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
!      bxy(1,j,k) = zt4
!      bxy(2,j,k) = zt5
!      bxy(3,j,k) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
!c update magnetic field half time step, ky < 0
!      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
!      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
!      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
!      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
!      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
!      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
!      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
!      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      exy(1,j,k1) = zt7
!      exy(2,j,k1) = zt8
!      exy(3,j,k1) = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
!      zt4 = zt4 + dth*(dky*zt1)
!      zt5 = zt5 + dth*(dkx*zt1)
!      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
!      bxy(1,j,k1) = zt4
!      bxy(2,j,k1) = zt5
!      bxy(3,j,k1) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
!   10 continue
!   20 continue
!c mode numbers kx = 0, nx/2
!cdir$ ivdep
!      do 30 k = 2, nyh
!      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      afdt = adt*aimag(ffc(1,k))
!c update magnetic field half time step
!      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
!      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
!      zt4 = bxy(1,1,k) - dth*(dky*zt1)
!      zt6 = bxy(3,1,k) + dth*(dky*zt3)
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
!      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt3 = cmplx(-aimag(zt7),real(zt7))
!      exy(1,1,k) = zt7
!      exy(2,1,k) = zero
!      exy(3,1,k) = zt9
!      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
!      zt4 = zt4 - dth*(dky*zt1)
!      zt6 = zt6 + dth*(dky*zt3)
!      bxy(1,1,k) = zt4
!      bxy(2,1,k) = zero
!      bxy(3,1,k) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
!      bxy(1,1,k1) = zero
!      bxy(2,1,k1) = zero
!      bxy(3,1,k1) = zero
!      exy(1,1,k1) = zero
!      exy(2,1,k1) = zero
!      exy(3,1,k1) = zero
!   30 continue
!c mode numbers ky = 0, ny/2
!      k1 = nyh + 1
!      do 40 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      afdt = adt*aimag(ffc(j,1))
!c update magnetic field half time step
!      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
!      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
!      zt5 = bxy(2,j,1) + dth*(dkx*zt1)
!      zt6 = bxy(3,j,1) - dth*(dkx*zt2)
!c update electric field whole time step
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt8 = exy(2,j,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1)
!      zt9 = exy(3,j,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1)
!c update magnetic field half time step and store electric field
!      zt1 = cmplx(-aimag(zt9),real(zt9))
!      zt2 = cmplx(-aimag(zt8),real(zt8))
!      exy(1,j,1) = zero
!      exy(2,j,1) = zt8
!      exy(3,j,1) = zt9
!      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
!      zt5 = zt5 + dth*(dkx*zt1)
!      zt6 = zt6 - dth*(dkx*zt2)
!      bxy(1,j,1) = zero
!      bxy(2,j,1) = zt5
!      bxy(3,j,1) = zt6
!      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
!      bxy(1,j,k1) = zero
!      bxy(2,j,k1) = zero
!      bxy(3,j,k1) = zero
!      exy(1,j,k1) = zero
!      exy(2,j,k1) = zero
!      exy(3,j,k1) = zero
!   40 continue
!      bxy(1,1,1) = zero
!      bxy(2,1,1) = zero
!      bxy(3,1,1) = zero
!      exy(1,1,1) = zero
!      exy(2,1,1) = zero
!      exy(3,1,1) = zero
!      bxy(1,1,k1) = zero
!      bxy(2,1,k1) = zero
!      bxy(3,1,k1) = zero
!      exy(1,1,k1) = zero
!      exy(2,1,k1) = zero
!      exy(3,1,k1) = zero
!      wf = float(nx*ny)*ws
!      wm = float(nx*ny)*c2*wp

      endif

      return
      end
!*****
      subroutine MAXWELBYEE2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,
     1nxhd,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse magnetic field for a half time step
c input: all, output: wf, wm, bxy
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c j,k = fourier mode numbers, except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
c bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      complex phix, phiy
      real kxop, kyop, thetax, thetay
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      thetay = dky/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      phiy = cmplx(cos(thetay),sin(thetay))
      kxop = 2.0*sin(thetax)
      kyop = 2.0*sin(thetay)   
      afdt = adt*aimag(ffc(j,k))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(kyop*zt1*phiy)
      zt5 = bxy(2,j,k) + dth*(kxop*zt1*phix)
      zt6 = bxy(3,j,k) - dth*(kxop*zt2*phix - kyop*zt3*phiy)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      !print*, j, k, zt5
c update magnetic field half time step, ky < 0
      !print*, 'E', j, k, exy(1,j,k1), exy(2,j,k1), exy(3,j,k1)
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(kyop*zt1*conjg(phiy))
      zt5 = bxy(2,j,k1) + dth*(kxop*zt1*phix)
      zt6 = bxy(3,j,k1) - dth*(kxop*zt2*phix + kyop*zt3*conjg(phiy))
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      !print*, 'B', j, k, zt4, zt5, zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6)) 
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      dkx = float(nx/2)
      thetax = dkx/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      kxop = 2.0*sin(thetax)
!kx = 0
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      phiy = cmplx(cos(thetay),sin(thetay))
      kyop = 2.0*sin(thetay) 
      afdt = adt*aimag(ffc(1,k))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(kyop*zt1*phiy)
      zt5 = bxy(2,1,k)
      zt6 = bxy(3,1,k) + dth*(kyop*zt3*phiy)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zt5
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
!kx = nx/2
!      zt1 = cmplx(-aimag(exy(3,1,k1)),real(exy(3,1,k1)))
!      zt2 = cmplx(-aimag(exy(2,1,k1)),real(exy(2,1,k1)))
!      zt3 = cmplx(-aimag(exy(1,1,k1)),real(exy(1,1,k1)))
!      zt4 = bxy(1,1,k1) + dth*(kyop*zt1*conjg(phiy))
!      zt5 = bxy(2,1,k1) - dth*(kxop*zt1*conjg(phix))
!      zt6 = bxy(3,1,k1) + dth*(kxop*zt2*conjg(phix) - 
!     1 kyop*zt3*conjg(phiy))
!      bxy(1,1,k1) = zt4
!      bxy(2,1,k1) = zt5
!      bxy(3,1,k1) = zt6
!      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
!      print*, 'dBx =', dth*(kyop*zt1*conjg(phiy))
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      kxop = 2.0*sin(thetax)
      afdt = adt*aimag(ffc(j,1))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
      zt4 = bxy(1,j,1)
      zt5 = bxy(2,j,1) + dth*(kxop*zt1*phix)
      zt6 = bxy(3,j,1) - dth*(kxop*zt2*phix)
      bxy(1,j,1) = zt4
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      !print*, j, k, k1, bxy(2,j,1)
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   40 continue
!      bxy(1,1,1) = zero
!      bxy(2,1,1) = zero
!      bxy(3,1,1) = zero
!      exy(1,1,1) = zero
!      exy(2,1,1) = zero
!      exy(3,1,1) = zero
      bxy(1,1,1) = cmplx(real(bxy(1,1,1)),0.0)
      bxy(2,1,1) = cmplx(real(bxy(2,1,1)),0.0)
      bxy(3,1,1) = cmplx(real(bxy(3,1,1)),0.0)
      exy(1,1,1) = cmplx(real(exy(1,1,1)),0.0)
      exy(2,1,1) = cmplx(real(exy(2,1,1)),0.0)
      exy(3,1,1) = cmplx(real(exy(3,1,1)),0.0) 
      
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
!
      subroutine MAXWELEYEE2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,
     1nxhd,nyhd,parmax)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric field for a full time step with periodic boundary
c conditions.
c input: all, output: wf, wm, exy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      complex phix, phiy
      integer parmax
      real kxop, kyop, thetax, thetay
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      thetay = dky/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      phiy = cmplx(cos(thetay),sin(thetay))
      kxop = 2.0*sin(thetax)
      kyop = 2.0*sin(thetay)
      afdt = adt*aimag(ffc(j,k))
c update electric field whole time step, ky > 0
      zt4 = bxy(1,j,k)
      zt5 = bxy(2,j,k) 
      zt6 = bxy(3,j,k)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(kyop*zt1*conjg(phiy))
     1 - afdt*cu(1,j,k)!*phix
      zt8 = exy(2,j,k) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(kxop*zt2*conjg(phix) - 
     1 kyop*zt3*conjg(phiy)) - afdt*cu(3,j,k)
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
c update electric field whole time step, ky < 0
      zt4 = bxy(1,j,k1)
      zt5 = bxy(2,j,k1)
      zt6 = bxy(3,j,k1)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(kyop*zt1*phiy) - afdt*cu(1,j,k1)!*phix
      zt8 = exy(2,j,k1) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(kxop*zt2*conjg(phix) + kyop*zt3*phiy)
     1 - afdt*cu(3,j,k1)
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      !print*, j, k, exy(1,j,k1), exy(2,j,k1), exy(3,j,k1)
   10 continue
   20 continue
c mode numbers kx = 0
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      thetay = dky/2.0
      phiy = cmplx(cos(thetay),sin(thetay))
      kyop = 2.0*sin(thetay) 
      afdt = adt*aimag(ffc(1,k))
      afdt = adt*aimag(ffc(1,k))
c update electric field whole time step
      zt4 = bxy(1,1,k)
      zt6 = bxy(3,1,k)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(kyop*zt1*conjg(phiy)) - afdt*cu(1,1,k)
      if(parmax==1) then
      zt8 = exy(2,1,k) - afdt*cu(2,1,k)
      else
      zt8 = zero
      endif
      zt9 = exy(3,1,k) - cdt*(kyop*zt3*conjg(phiy)) - afdt*cu(3,1,k)
      exy(1,1,k) = zt7
      exy(2,1,k) = zt8
      exy(3,1,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
     1 zt9*conjg(zt9))
! mode numbers kx = nx/2
!      dkx = float(nx/2)
!      thetax = dkx/2.0
!      phix = cmplx(cos(thetax),sin(thetax))
!      kxop = 2.0*sin(thetax)
!      zt4 = bxy(1,1,k1)
!      zt5 = bxy(2,1,k1)
!      zt6 = bxy(3,1,k1)
!      zt1 = cmplx(-aimag(zt6),real(zt6))
!      zt2 = cmplx(-aimag(zt5),real(zt5))
!      zt3 = cmplx(-aimag(zt4),real(zt4))
!      zt7 = exy(1,1,k1) - cdt*(kyop*zt1*phiy) - afdt*cu(1,1,k1)
!      zt8 = exy(2,1,k1) + cdt*(kxop*zt1*phix) - afdt*cu(2,1,k1)
!      zt9 = exy(3,1,k1) - cdt*(kxop*zt2*phix - kyop*zt3*phiy) - 
!     1 afdt*cu(3,1,k1)
!      exy(1,1,k1) = zt7
!      exy(2,1,k1) = zt8
!      exy(3,1,k1) = zt9
!      !print*, exy(3,1,k1)
!      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      thetax = dkx/2.0
      phix = cmplx(cos(thetax),sin(thetax))
      kxop = 2.0*sin(thetax)
      afdt = adt*aimag(ffc(j,1))
c update electric field whole time step
      zt5 = bxy(2,j,1)
      zt6 = bxy(3,j,1)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      if(parmax==1) then
      zt7 = exy(1,j,1) - afdt*cu(1,j,1)!*phix
      else
      zt7 = zero
      endif
      zt8 = exy(2,j,1) - cdt*(kxop*zt1*conjg(phix)) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(kxop*zt2*conjg(phix)) - afdt*cu(3,j,1)
      exy(1,j,1) = zt7
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      !print*, j, k, zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
   40 continue
!      exy(1,1,1) = zero
!      exy(2,1,1) = zero
!      exy(3,1,1) = zero
!      bxy(1,1,1) = zero
!      bxy(2,1,1) = zero
!      bxy(3,1,1) = zero
      bxy(1,1,1) = cmplx(real(bxy(1,1,1)),0.0)
      bxy(2,1,1) = cmplx(real(bxy(2,1,1)),0.0)
      bxy(3,1,1) = cmplx(real(bxy(3,1,1)),0.0)
      exy(1,1,1) = cmplx(real(exy(1,1,1)),0.0)
      exy(2,1,1) = cmplx(real(exy(2,1,1)),0.0)
      exy(3,1,1) = cmplx(real(exy(3,1,1)),0.0) 
    
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
      subroutine MAXWELB2(exy,bxy,ci,dt,wm,nx,ny,nxvh,nyv,
     1nxhd,nyhd)
!      subroutine MAXWELB2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,
!     1nxhd,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse magnetic field for a half time step
c input: all, output: wf, wm, bxy
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c j,k = fourier mode numbers, except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
c bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex exy, bxy!, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv)!, cu(3,nxvh,nyv)
      !dimension ffc(nxhd,nyhd)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)  
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k) - dth*(dkx*zt2 - dky*zt3)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6)) 
c update magnetic field half time step, ky < 0
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6)) 
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
!kx = 0
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
      zt5 = bxy(2,1,k)
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zt5
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
!kx = nx/2
C      dkx = dnx*float(nx/2)
C      zt1 = cmplx(-aimag(exy(3,1,k1)),real(exy(3,1,k1)))
C      zt2 = cmplx(-aimag(exy(2,1,k1)),real(exy(2,1,k1)))
C      zt3 = cmplx(-aimag(exy(1,1,k1)),real(exy(1,1,k1)))
C      zt4 = bxy(1,1,k1) + dth*(dky*zt1)
C      zt5 = bxy(2,1,k1) - dth*(dkx*zt1)
C      zt6 = bxy(3,1,k1) + dth*(dkx*zt2 - dky*zt3)
C      bxy(1,1,k1) = zt4
C      bxy(2,1,k1) = zt5
C      bxy(3,1,k1) = zt6
      bxy(1,1,k1) = zero ! (0,ny/2;nx/2,ny/2) => (zero, zero)
      bxy(2,1,k1) = cmplx(real(bxy(2,1,k1)),0.0) ! (0,ny/2;nx/2,ny/2) => (no shift, zero)
      bxy(3,1,k1) = zero ! (0,ny/2;nx/2,ny/2) => (zero, zero)
      wp = wp + anorm*(bxy(2,1,k1)*conjg(bxy(2,1,k1)))
!      bxy(1,1,k1) = zero
!      bxy(2,1,k1) = zero
!      bxy(3,1,k1) = zero

   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
      zt4 = bxy(1,j,1)
      zt5 = bxy(2,j,1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,1) - dth*(dkx*zt2)
      bxy(1,j,1) = zt4
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
! ky = ny/2
C      dky = dny*float(ny/2)
C      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
C      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
C      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
C      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
C      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
C      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
C      bxy(1,j,k1) = zt4
C      bxy(2,j,k1) = zt5
C      bxy(3,j,k1) = zt6
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
   40 continue
!      bxy(1,1,1) = zero
!      bxy(2,1,1) = zero
!      bxy(3,1,1) = zero
!      bxy(1,1,1) = cmplx(real(bxy(1,1,1)),0.0)
      bxy(1,1,1) = bxy(1,1,1) ! (0,0;nx/2,0) => (no shift, no shift)
      bxy(2,1,1) = cmplx(real(bxy(2,1,1)),0.0) ! (0,0;nx/2,0) => (no shift, zero)
      bxy(3,1,1) = cmplx(real(bxy(3,1,1)),0.0) ! (0,0;nx/2,0) => (no shift, zero)
      zt4 = bxy(1,j,1)
      zt5 = bxy(2,j,1)
      zt6 = bxy(3,j,1)
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6)) 
      wm = float(nx*ny)*c2*wp
      return
      end
!
      subroutine MAXWELE2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,
     1nxhd,nyhd,parmax)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric field for a full time step with periodic boundary
c conditions.
c input: all, output: wf, wm, exy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      integer parmax
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffc(j,k))
c update electric field whole time step, ky > 0
      zt4 = bxy(1,j,k)
      zt5 = bxy(2,j,k) 
      zt6 = bxy(3,j,k)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) - cdt*(dkx*zt1) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k)
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
c update electric field whole time step, ky < 0
      zt4 = bxy(1,j,k1)
      zt5 = bxy(2,j,k1)
      zt6 = bxy(3,j,k1)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
   10 continue
   20 continue
c mode numbers kx = 0
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffc(1,k))
c update electric field whole time step
      zt4 = bxy(1,1,k)
      zt6 = bxy(3,1,k)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
      if(parmax==1) then
      zt8 = exy(2,1,k) - afdt*cu(2,1,k)
      else
      zt8 = zero
      endif
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
      exy(1,1,k) = zt7
      exy(2,1,k) = zt8
      exy(3,1,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
     1 zt9*conjg(zt9))
! mode numbers kx = nx/2
C      dkx = dnx*float(nx/2)
C      zt4 = bxy(1,1,k1)
C      zt5 = bxy(2,1,k1)
C      zt6 = bxy(3,1,k1)
C      zt1 = cmplx(-aimag(zt6),real(zt6))
C      zt2 = cmplx(-aimag(zt5),real(zt5))
C      zt3 = cmplx(-aimag(zt4),real(zt4))
C      zt7 = exy(1,1,k1) + cdt*(dky*zt1) - afdt*cu(1,1,k1)
C      zt8 = exy(2,1,k1) + cdt*(dkx*zt1) - afdt*cu(2,1,k1)
C      zt9 = exy(3,1,k1) - cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,1,k1)
C      exy(1,1,k1) = zt7
C      exy(2,1,k1) = zt8
C      exy(3,1,k1) = zt9
!      !print*, exy(3,1,k1)
!      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
!      exy(1,1,k1) = zero
!      exy(2,1,k1) = zero
!      exy(3,1,k1) = zero
      exy(1,1,k1) = cmplx(real(exy(1,1,k1)),0.0)! (0,ny/2;nx/2,ny/2) => (no shift, zero)
      exy(2,1,k1) = zero ! (0,ny/2;nx/2,ny/2) => (zero, zero)
      exy(3,1,k1) = exy(3,1,k1) ! (0,ny/2;nx/2,ny/2) => (no shift, no shift)
      zt7 = exy(1,j,k1)
      zt9 = exy(3,j,k1)
      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffc(j,1))
c update electric field whole time step
! ky = 0
      zt5 = bxy(2,j,1)
      zt6 = bxy(3,j,1)
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      if(parmax==1) then
      zt7 = exy(1,j,1) - afdt*cu(1,j,1)
      else
      zt7 = zero
      endif
      zt8 = exy(2,j,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1)
      exy(1,j,1) = zt7
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      !print*, j, k, zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
! ky = ny/2
C      dky = dny*float(ny/2)
C      zt4 = bxy(1,j,k1)
C      zt5 = bxy(2,j,k1)
C      zt6 = bxy(3,j,k1)
C      zt1 = cmplx(-aimag(zt6),real(zt6))
C      zt2 = cmplx(-aimag(zt5),real(zt5))
C      zt3 = cmplx(-aimag(zt4),real(zt4))
C      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
C      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
C      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
C      exy(1,j,k1) = zt7
C      exy(2,j,k1) = zt8
C      exy(3,j,k1) = zt9
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   40 continue
!      exy(1,1,1) = cmplx(real(exy(1,1,1)),0.0)
!      exy(2,1,1) = cmplx(real(exy(2,1,1)),0.0)
!      exy(3,1,1) = cmplx(real(exy(3,1,1)),0.0) 
      exy(1,1,1) = cmplx(real(exy(1,1,1)),0.0)      
      exy(2,1,1) = exy(2,1,1) ! (0,0;nx/2,0) => (no shift, no shift)
      exy(3,1,1) = exy(3,1,1) ! (0,0;nx/2,0) => (no shift, no shift)
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) +
     1 zt9*conjg(zt9))
      wf = float(nx*ny)*ws
      return
      end
!*****
      subroutine MAXWEL2PS(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nxh
     1d,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      real knorm
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c1 = 1./ci
      c1dt = c1*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      knorm = sqrt(dkx*dkx+dky*dky)
      ss = sin(knorm*c1dt)
      cc = cos(knorm*c1dt)
      p1 = ss / (knorm*c1)
      p2 = (1 - cc) / (knorm*knorm*c2)
      afdt = adt*aimag(ffc(j,k))
      affpd = affp*aimag(ffc(j,k))
c update electric field whole time step, ky > 0
      zt3 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt1 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      zt4 = cc*exy(1,j,k) + p1*c2*(dky*zt3) - p1*affpd*cu(1,j,k)
      zt5 = cc*exy(2,j,k) - p1*c2*(dkx*zt3) - p1*affpd*cu(2,j,k)
      zt6 = cc*exy(3,j,k) + p1*c2*(dkx*zt2-dky*zt1) - p1*affpd*cu(3,j,k)
c update magnetic field whole time step, ky > 0
      zt1 = p1*exy(1,j,k) - p2*affpd*cu(1,j,k)
      zt2 = p1*exy(2,j,k) - p2*affpd*cu(2,j,k)
      zt3 = p1*exy(3,j,k) - p2*affpd*cu(3,j,k)
      exy(1,j,k) = zt4
      exy(2,j,k) = zt5
      exy(3,j,k) = zt6  
      ws = ws + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      zt7 = cmplx(-aimag(zt1),real(zt1))
      zt8 = cmplx(-aimag(zt2),real(zt2))
      zt9 = cmplx(-aimag(zt3),real(zt3))
      zt4 = cc*bxy(1,j,k) - dky*zt9
      zt5 = cc*bxy(2,j,k) + dkx*zt9
      zt6 = cc*bxy(3,j,k) - (dkx*zt8-dky*zt7)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c update electric field whole time step, ky < 0
      zt3 = cmplx(-aimag(bxy(3,j,k1)),real(bxy(3,j,k1)))
      zt2 = cmplx(-aimag(bxy(2,j,k1)),real(bxy(2,j,k1)))
      zt1 = cmplx(-aimag(bxy(1,j,k1)),real(bxy(1,j,k1)))
      zt4 = cc*exy(1,j,k1) + p1*c2*(dky*zt3) - p1*affpd*cu(1,j,k1)
      zt5 = cc*exy(2,j,k1) - p1*c2*(dkx*zt3) - p1*affpd*cu(2,j,k1)
      zt6 = cc*exy(3,j,k1)+p1*c2*(dkx*zt2-dky*zt1) - p1*affpd*cu(3,j,k1)
c update magnetic field whole time step, ky > 0
      zt1 = p1*exy(1,j,k1) - p2*affpd*cu(1,j,k1)
      zt2 = p1*exy(2,j,k1) - p2*affpd*cu(2,j,k1)
      zt3 = p1*exy(3,j,k1) - p2*affpd*cu(3,j,k1)
      exy(1,j,k1) = zt4
      exy(2,j,k1) = zt5
      exy(3,j,k1) = zt6  
      ws = ws + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      zt7 = cmplx(-aimag(zt1),real(zt1))
      zt8 = cmplx(-aimag(zt2),real(zt2))
      zt9 = cmplx(-aimag(zt3),real(zt3))
      zt4 = cc*bxy(1,j,k1) - dky*zt9
      zt5 = cc*bxy(2,j,k1) + dkx*zt9
      zt6 = cc*bxy(3,j,k1) - (dkx*zt8-dky*zt7)
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dkx = 0.0
      afdt = adt*aimag(ffc(1,k))
      affpd = affp*aimag(ffc(1,k))
      knorm = sqrt(dkx*dkx+dky*dky)
      ss = sin(knorm*c1dt)
      cc = cos(knorm*c1dt)
      p1 = ss / (knorm*c1)
      p2 = (1 - cc) / (knorm*knorm*c2)
      zt3 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      zt1 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      zt4 = cc*exy(1,1,k) + p1*c2*(dky*zt3) - p1*affpd*cu(1,1,k)
      zt6 = cc*exy(3,1,k) - p1*c2*(dky*zt1) - p1*affpd*cu(3,1,k)
c update magnetic field whole time step, ky > 0
      zt1 = p1*exy(1,1,k) - p2*affpd*cu(1,1,k)
      zt3 = p1*exy(3,1,k) - p2*affpd*cu(3,1,k)
      exy(1,1,k) = zt4
      exy(2,1,k) = zero
      exy(3,1,k) = zt6  
      ws = ws + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
      zt7 = cmplx(-aimag(zt1),real(zt1))
      zt9 = cmplx(-aimag(zt3),real(zt3))
      zt4 = cc*bxy(1,1,k) - dky*zt9
      zt6 = cc*bxy(3,1,k) + dky*zt7
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      dky = 0.0
      afdt = adt*aimag(ffc(j,1))
      affpd = affp*aimag(ffc(j,1))
      knorm = sqrt(dkx*dkx+dky*dky)
      ss = sin(knorm*c1dt)
      cc = cos(knorm*c1dt)
      p1 = ss / (knorm*c1)
      p2 = (1 - cc) / (knorm*knorm*c2)
      zt3 = cmplx(-aimag(bxy(3,j,1)),real(bxy(3,j,1)))
      zt2 = cmplx(-aimag(bxy(2,j,1)),real(bxy(2,j,1)))
      zt5 = cc*exy(2,j,1) - p1*c2*(dkx*zt3) - p1*affpd*cu(2,j,1)
      zt6 = cc*exy(3,j,1) + p1*c2*(dkx*zt2) - p1*affpd*cu(3,j,1)
c update magnetic field whole time step, ky > 0
      zt2 = p1*exy(2,j,1) - p2*affpd*cu(2,j,1)
      zt3 = p1*exy(3,j,1) - p2*affpd*cu(3,j,1)
      exy(1,j,1) = zero
      exy(2,j,1) = zt5
      exy(3,j,1) = zt6  
      ws = ws + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      zt8 = cmplx(-aimag(zt2),real(zt2))
      zt9 = cmplx(-aimag(zt3),real(zt3))
      zt5 = cc*bxy(2,j,1) + dkx*zt9
      zt6 = cc*bxy(3,j,1) - dkx*zt8
      bxy(1,j,1) = zero
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   40 continue
      !bxy(1,1,1) = zero
      !bxy(2,1,1) = zero
      !bxy(3,1,1) = zero
      !exy(1,1,1) = zero
      !exy(2,1,1) = zero
      !exy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
      subroutine MAXWEL2T(exyt,bxyt,cut,ffct,ci,dt,wf,wm,nx,ny,nxvh,nyv,
     1nxhd,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
c by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
c bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
c                       - affp*dt*cux(kx,ky)*s(kx,ky)
c ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
c                       - affp*dt*cuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*dt*cuz(kx,ky)*s(kx,ky)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c and similarly for bx, by, bz.
c cut(k,i,j) = complex current density
c exyt(k,i,j) = complex transverse electric field
c bxyt(k,i,j) = complex magnetic field
c for component i, all for fourier mode (k-1,j-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffct(k,j)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (k-1,j-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, dt, wf, wm
      complex exyt, bxyt, cut, ffct
      dimension exyt(nyv,3,nxvh), bxyt(nyv,3,nxvh), cut(nyv,3,nxvh)
      dimension ffct(nyhd,nxhd)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      real dnx, dny, dth, c2, cdt, affp, anorm, dkx, dky, afdt, adt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffct(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 2, nxh
      dkx = dnx*real(j - 1)
      do 10 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      afdt = adt*aimag(ffct(k,j))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exyt(k,3,j)),real(exyt(k,3,j)))
      zt2 = cmplx(-aimag(exyt(k,2,j)),real(exyt(k,2,j)))
      zt3 = cmplx(-aimag(exyt(k,1,j)),real(exyt(k,1,j)))
      zt4 = bxyt(k,1,j) - dth*(dky*zt1)
      zt5 = bxyt(k,2,j) + dth*(dkx*zt1)
      zt6 = bxyt(k,3,j) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyt(k,1,j) + cdt*(dky*zt1) - afdt*cut(k,1,j)
      zt8 = exyt(k,2,j) - cdt*(dkx*zt1) - afdt*cut(k,2,j)
      zt9 = exyt(k,3,j) + cdt*(dkx*zt2 - dky*zt3) - afdt*cut(k,3,j)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyt(k,1,j) = zt7
      exyt(k,2,j) = zt8
      exyt(k,3,j) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxyt(k,1,j) = zt4
      bxyt(k,2,j) = zt5
      bxyt(k,3,j) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
c update magnetic field half time step, ky < 0
      zt1 = cmplx(-aimag(exyt(k1,3,j)),real(exyt(k1,3,j)))
      zt2 = cmplx(-aimag(exyt(k1,2,j)),real(exyt(k1,2,j)))
      zt3 = cmplx(-aimag(exyt(k1,1,j)),real(exyt(k1,1,j)))
      zt4 = bxyt(k1,1,j) + dth*(dky*zt1)
      zt5 = bxyt(k1,2,j) + dth*(dkx*zt1)
      zt6 = bxyt(k1,3,j) - dth*(dkx*zt2 + dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyt(k1,1,j) - cdt*(dky*zt1) - afdt*cut(k1,1,j)
      zt8 = exyt(k1,2,j) - cdt*(dkx*zt1) - afdt*cut(k1,2,j)
      zt9 = exyt(k1,3,j) + cdt*(dkx*zt2 + dky*zt3) - afdt*cut(k1,3,j)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyt(k1,1,j) = zt7
      exyt(k1,2,j) = zt8
      exyt(k1,3,j) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 + dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxyt(k1,1,j) = zt4
      bxyt(k1,2,j) = zt5
      bxyt(k1,3,j) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffct(1,j))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exyt(1,3,j)),real(exyt(1,3,j)))
      zt2 = cmplx(-aimag(exyt(1,2,j)),real(exyt(1,2,j)))
      zt5 = bxyt(1,2,j) + dth*(dkx*zt1)
      zt6 = bxyt(1,3,j) - dth*(dkx*zt2)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = exyt(1,2,j) - cdt*(dkx*zt1) - afdt*cut(1,2,j)
      zt9 = exyt(1,3,j) + cdt*(dkx*zt2) - afdt*cut(1,3,j)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exyt(1,1,j) = zero
      exyt(1,2,j) = zt8
      exyt(1,3,j) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxyt(1,1,j) = zero
      bxyt(1,2,j) = zt5
      bxyt(1,3,j) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      bxyt(k1,1,j) = zero
      bxyt(k1,2,j) = zero
      bxyt(k1,3,j) = zero
      exyt(k1,1,j) = zero
      exyt(k1,2,j) = zero
      exyt(k1,3,j) = zero
   30 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 40 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      afdt = adt*aimag(ffct(k,1))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exyt(k,3,1)),real(exyt(k,3,1)))
      zt3 = cmplx(-aimag(exyt(k,1,1)),real(exyt(k,1,1)))
      zt4 = bxyt(k,1,1) - dth*(dky*zt1)
      zt6 = bxyt(k,3,1) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exyt(k,1,1) + cdt*(dky*zt1) - afdt*cut(k,1,1)
      zt9 = exyt(k,3,1) - cdt*(dky*zt3) - afdt*cut(k,3,1)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exyt(k,1,1) = zt7
      exyt(k,2,1) = zero
      exyt(k,3,1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt6 = zt6 + dth*(dky*zt3)
      bxyt(k,1,1) = zt4
      bxyt(k,2,1) = zero
      bxyt(k,3,1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
      bxyt(k1,1,1) = conjg(zt4)
      bxyt(k1,2,1) = zero
      bxyt(k1,3,1) = conjg(zt6)
      exyt(k1,1,1) = conjg(zt7)
      exyt(k1,2,1) = zero
      exyt(k1,3,1) = conjg(zt9)
      bxyt(k,1,j1) = zero
      bxyt(k,2,j1) = zero
      bxyt(k,3,j1) = zero
      exyt(k,1,j1) = zero
      exyt(k,2,j1) = zero
      exyt(k,3,j1) = zero
      bxyt(k1,1,j1) = zero
      bxyt(k1,2,j1) = zero
      bxyt(k1,3,j1) = zero
      exyt(k1,1,j1) = zero
      exyt(k1,2,j1) = zero
      exyt(k1,3,j1) = zero
   40 continue
      k1 = nyh + 1
      bxyt(1,1,1) = zero
      bxyt(1,2,1) = zero
      bxyt(1,3,1) = zero
      exyt(1,1,1) = zero
      exyt(1,2,1) = zero
      exyt(1,3,1) = zero
      bxyt(k1,1,1) = zero
      bxyt(k1,2,1) = zero
      bxyt(k1,3,1) = zero
      exyt(k1,1,1) = zero
      exyt(k1,2,1) = zero
      exyt(k1,3,1) = zero
      bxyt(1,1,j1) = zero
      bxyt(1,2,j1) = zero
      bxyt(1,3,j1) = zero
      exyt(1,1,j1) = zero
      exyt(1,2,j1) = zero
      exyt(1,3,j1) = zero
      bxyt(k1,1,j1) = zero
      bxyt(k1,2,j1) = zero
      bxyt(k1,3,j1) = zero
      exyt(k1,1,j1) = zero
      exyt(k1,2,j1) = zero
      exyt(k1,3,j1) = zero
      wf = real(nx*ny)*ws
      wm = real(nx*ny)*c2*wp
      return
      end
      subroutine EMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
!*****
!Now includes option to add fxy to exy if isign = 0.
!*****
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      integer i, j, k, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
c add the fields
         !print*, 'bef', isign, fxy(3,140,513)
      if (isign.gt.0) then
         !print*, 'bef', fxy(3,140,513)
         do 30 k = 2, nyh
         k1 = ny2 - k
         do 20 j = 1, nxh
         at1 = aimag(ffc(j,k))
         do 10 i = 1, 3
         fxy(i,j,k) = fxy(i,j,k) + exy(i,j,k)*at1
         fxy(i,j,k1) = fxy(i,j,k1) + exy(i,j,k1)*at1
         !if(k1.gt.512) print*, k1, fxy(i,j,k1)
   10    continue
   20    continue
   30    continue
         !print*, 'aft', fxy(3,140,513)
         k1 = nyh + 1
         do 50 j = 1, nxh
         at1 = aimag(ffc(j,1))
         do 40 i = 1, 3
         fxy(i,j,1) = fxy(i,j,1) + exy(i,j,1)*at1
         fxy(i,j,k1) = fxy(i,j,k1) + exy(i,j,k1)*at1
   40    continue
   50    continue
         !print*, 'aft', fxy(3,140,513)
c copy the fields
      else if (isign.lt.0) then
         do 80 k = 2, nyh
         k1 = ny2 - k
         do 70 j = 1, nxh
         at1 = aimag(ffc(j,k))
         do 60 i = 1, 3
         fxy(i,j,k) = exy(i,j,k)*at1
         fxy(i,j,k1) = exy(i,j,k1)*at1
   60    continue
   70    continue
   80    continue
         k1 = nyh + 1
         do 100 j = 1, nxh
         at1 = aimag(ffc(j,1))
         do 90 i = 1, 3
         fxy(i,j,1) = exy(i,j,1)*at1
         fxy(i,j,k1) = exy(i,j,k1)*at1
   90    continue
  100    continue
!add fxy to exy for later use in maxwel2 solver
      else if (isign==0) then
         do 130 k = 2, nyh
         k1 = ny2 - k
         do 120 j = 1, nxh
         do 110 i = 1, 3
         exy(i,j,k) = exy(i,j,k) + fxy(i,j,k)
         exy(i,j,k1) = exy(i,j,k1) + fxy(i,j,k1)
  110    continue
  120    continue
  130    continue
         k1 = nyh + 1
         do 150 j = 1, nxh
         do 140 i = 1, 3
         exy(i,j,1) = exy(i,j,1) + fxy(i,j,1)
         exy(i,j,k1) = exy(i,j,k1) + fxy(i,j,k1)
  140    continue
  150    continue
      endif
      !print*, 'aft', isign, fxy(3,140,513)
      return
      end
      subroutine LSREMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      integer i, j, k, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
c add the fields
      if (isign.gt.0) then
!      print*, '---'
!      print*, fxy(3,:,1)
!      print*, fxy(3,:,2)
!      print*, fxy(3,:,3)
!      print*, fxy(3,:,4)
!      print*, fxy(3,:,5)
!      print*, fxy(3,:,6)
!      print*, fxy(3,:,7)
!      print*, fxy(3,:,8)
!      print*, fxy(3,:,9)
!      print*, fxy(3,:,10)
!      print*, fxy(3,:,11)
!      print*, fxy(3,:,12)
!      print*, fxy(3,:,13)
!      print*, fxy(3,:,14)
!      print*, fxy(3,:,15)
!      print*, fxy(3,:,16)
!      print*, fxy(3,:,17)
!      print*, fxy(3,:,18)
!      print*, fxy(3,:,19)
!      print*, fxy(3,:,20)
!      print*, fxy(3,:,21)
!      print*, fxy(3,:,22)
!      print*, fxy(3,:,23)
!      print*, fxy(3,:,24)
!      print*, fxy(3,:,25)
!      print*, fxy(3,:,26)
!      print*, fxy(3,:,27)
!      print*, fxy(3,:,28)
!      print*, fxy(3,:,29)
!      print*, fxy(3,:,30)
!      print*, fxy(3,:,31)
!      print*, fxy(3,:,32)
!      print*, fxy(3,:,33)
         do 30 k = 2, nyh
         k1 = ny2 - k
         do 20 j = 1, nxh
         !at1 = aimag(ffc(j,k))
         do 10 i = 1, 3
         exy(i,j,k) =  exy(i,j,k) + fxy(i,j,k)
         exy(i,j,k1) =  exy(i,j,k1) + fxy(i,j,k1)
   10    continue
   20    continue
   30    continue
         k1 = nyh + 1
         do 50 j = 1, nxh
         !at1 = aimag(ffc(j,1))
         do 40 i = 1, 3
         exy(i,j,1) =  exy(i,j,1) + fxy(i,j,1)
         exy(i,j,k1) = exy(i,j,k1) + fxy(i,j,k1)
   40    continue
   50    continue
!      print*, '---'
!      print*, exy(3,:,1)
!      print*, exy(3,:,2)
!      print*, exy(3,:,3)
!      print*, exy(3,:,4)
!      print*, exy(3,:,5)
!      print*, exy(3,:,6)
!      print*, exy(3,:,7)
!      print*, exy(3,:,8)
!      print*, exy(3,:,9)
!      print*, exy(3,:,10)
!      print*, exy(3,:,11)
!      print*, exy(3,:,12)
!      print*, exy(3,:,13)
!      print*, exy(3,:,14)
!      print*, exy(3,:,15)
!      print*, exy(3,:,16)
!      print*, exy(3,:,17)
!      print*, exy(3,:,18)
!      print*, exy(3,:,19)
!      print*, exy(3,:,20)
!      print*, exy(3,:,21)
!      print*, exy(3,:,22)
!      print*, exy(3,:,23)
!      print*, exy(3,:,24)
!      print*, exy(3,:,25)
!      print*, exy(3,:,26)
!      print*, exy(3,:,27)
!      print*, exy(3,:,28)
!      print*, exy(3,:,29)
!      print*, exy(3,:,30)
!      print*, exy(3,:,31)
!      print*, exy(3,:,32)
!      print*, exy(3,:,33)
c copy the fields
      else if (isign.lt.0) then
         print *, 'isign shoud be greater than 1'
      endif
      return
      end
      subroutine LSREMFIELD4(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
c exy(half,Fourier) to fxy(full,Fourier) [substitution]
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      integer i, j, k, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
c add the fields
      if (isign.gt.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         do 20 j = 1, nxh
         !at1 = aimag(ffc(j,k))
         do 10 i = 1, 3
         fxy(i,j,k) = exy(i,j,k)
         fxy(i,j,k1) = exy(i,j,k1)
   10    continue
   20    continue
   30    continue
         k1 = nyh + 1
         do 50 j = 1, nxh
         !at1 = aimag(ffc(j,1))
         do 40 i = 1, 3
         fxy(i,j,1) = exy(i,j,1)
         fxy(i,j,k1) = exy(i,j,k1)
   40    continue
   50    continue
c copy the fields
      else if (isign.lt.0) then
         print *, 'isign shoud be greater than 1'
      endif
      return
      end
      subroutine LSREMFIELD3(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
      integer i, j, k, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
c add the fields
      if (isign.gt.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         do 20 j = 1, nxh
         !at1 = aimag(ffc(j,k))
         do 10 i = 1, 3
         exy(i,j,k) =  fxy(i,j,k)
         exy(i,j,k1) =  fxy(i,j,k1)
   10    continue
   20    continue
   30    continue
         k1 = nyh + 1
         do 50 j = 1, nxh
         !at1 = aimag(ffc(j,1))
         do 40 i = 1, 3
         exy(i,j,1) =  fxy(i,j,1)
         exy(i,j,k1) = fxy(i,j,k1)
   40    continue
   50    continue
c copy the fields
      else if (isign.lt.0) then
         print *, 'isign shoud be greater than 1'
      endif
      return
      end
      subroutine EMFIELDR2(fxy,exy,ffd,isign,nx,ny,nxv,nxe,nye)
c this subroutine either adds real vector fields if isign > 0
c or copies real vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, nxv, nxe, nye
      complex ffd
      real fxy, exy
      dimension fxy(3,nxe,nye), exy(3,nxe,nye)
      dimension ffd(nxv,ny)
c local data
      integer j, k, nx1
      real at1
      nx1 = nx + 1
c smooth the em field and add
      if (isign.gt.0) then
         do 20 k = 1, ny
         do 10 j = 1, nx
         at1 = aimag(ffd(j,k))
         fxy(1,j,k) = fxy(1,j,k) + exy(1,j,k)*at1
         fxy(2,j,k) = fxy(2,j,k) + exy(2,j,k)*at1
         fxy(3,j,k) = fxy(3,j,k) + exy(3,j,k)*at1
   10    continue
         fxy(1,nx+1,k) = fxy(1,nx+1,k) + exy(1,nx+1,k)
         fxy(2,nx+1,k) = fxy(2,nx+1,k) + exy(2,nx+1,k)
         fxy(3,nx+1,k) = fxy(3,nx+1,k) + exy(3,nx+1,k)
   20    continue
         do 30 j = 1, nx1
         fxy(1,j,ny+1) = fxy(1,j,ny+1) + exy(1,j,ny+1)
         fxy(2,j,ny+1) = fxy(2,j,ny+1) + exy(2,j,ny+1)
         fxy(3,j,ny+1) = fxy(3,j,ny+1) + exy(3,j,ny+1)
   30    continue
c copy and smooth the magnetic fields
      else if (isign.lt.0) then
         do 50 k = 1, ny
         do 40 j = 1, nx
         at1 = aimag(ffd(j,k))
         fxy(1,j,k) = exy(1,j,k)*at1
         fxy(2,j,k) = exy(2,j,k)*at1
         fxy(3,j,k) = exy(3,j,k)*at1
   40    continue
         fxy(1,nx+1,k) = exy(1,nx+1,k)
         fxy(2,nx+1,k) = exy(2,nx+1,k)
         fxy(3,nx+1,k) = exy(3,nx+1,k)
   50    continue
         do 60 j = 1, nx1
         fxy(1,j,ny+1) = exy(1,j,ny+1)
         fxy(2,j,ny+1) = exy(2,j,ny+1)
         fxy(3,j,ny+1) = exy(3,j,ny+1)
   60    continue
c copy the electric fields
      else
         do 80 k = 1, ny
         do 70 j = 1, nx
         fxy(1,j,k) = exy(1,j,k)
         fxy(2,j,k) = exy(2,j,k)
         fxy(3,j,k) = exy(3,j,k)
   70    continue
         fxy(1,nx+1,k) = exy(1,nx+1,k)
         fxy(2,nx+1,k) = exy(2,nx+1,k)
         fxy(3,nx+1,k) = exy(3,nx+1,k)
   80    continue
         do 90 j = 1, nx1
         fxy(1,j,ny+1) = exy(1,j,ny+1)
         fxy(2,j,ny+1) = exy(2,j,ny+1)
         fxy(3,j,ny+1) = exy(3,j,ny+1)
   90    continue
      endif
      return
      end
      subroutine EMFIELDC2(fxy,exy,ffb,isign,nx,ny,nxv,nxe,nyeh)
c this subroutine either adds real complex fields if isign > 0
c or copies real vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, nxv, nxe, nyeh
      complex fxy, exy, ffb
      dimension fxy(3,nxe,nyeh), exy(3,nxe,nyeh)
      dimension ffb(nxv,ny/2)
c local data
      integer j, k, nx1, nyh
      real at1
      nyh = max(1,ny/2)
      nx1 = nx + 1
c smooth the em field and add
      if (isign.gt.0) then
         do 20 k = 2, nyh
         do 10 j = 1, nx
         at1 = aimag(ffb(j,k))
         fxy(1,j,k) = fxy(1,j,k) + exy(1,j,k)*at1
         fxy(2,j,k) = fxy(2,j,k) + exy(2,j,k)*at1
         fxy(3,j,k) = fxy(3,j,k) + exy(3,j,k)*at1
   10    continue
         fxy(1,nx+1,k) = fxy(1,nx+1,k) + exy(1,nx+1,k)
         fxy(2,nx+1,k) = fxy(2,nx+1,k) + exy(2,nx+1,k)
         fxy(3,nx+1,k) = fxy(3,nx+1,k) + exy(3,nx+1,k)
   20    continue
         do 30 j = 1, nx
         at1 = aimag(ffb(j,1))
         fxy(1,j,1) = fxy(1,j,1) + cmplx(real(exy(1,j,1))*at1,aimag(exy(
     11,j,1)))
         fxy(2,j,1) = fxy(2,j,1) + cmplx(real(exy(2,j,1))*at1,aimag(exy(
     12,j,1)))
         fxy(3,j,1) = fxy(3,j,1) + cmplx(real(exy(3,j,1))*at1,aimag(exy(
     13,j,1)))
   30    continue
         fxy(1,nx+1,1) = fxy(1,nx+1,1) + exy(1,nx+1,1)
         fxy(2,nx+1,1) = fxy(2,nx+1,1) + exy(2,nx+1,1)
         fxy(3,nx+1,1) = fxy(3,nx+1,1) + exy(3,nx+1,1)
c copy and smooth the magnetic fields
      else if (isign.lt.0) then
         do 50 k = 2, nyh
         do 40 j = 1, nx
         at1 = aimag(ffb(j,k))
         fxy(1,j,k) = exy(1,j,k)*at1
         fxy(2,j,k) = exy(2,j,k)*at1
         fxy(3,j,k) = exy(3,j,k)*at1
   40    continue
         fxy(1,nx+1,k) = exy(1,nx+1,k)
         fxy(2,nx+1,k) = exy(2,nx+1,k)
         fxy(3,nx+1,k) = exy(3,nx+1,k)
   50    continue
         do 60 j = 1, nx
         at1 = aimag(ffb(j,1))
         fxy(1,j,1) = cmplx(real(exy(1,j,1))*at1,aimag(exy(1,j,1)))
         fxy(2,j,1) = cmplx(real(exy(2,j,1))*at1,aimag(exy(2,j,1)))
         fxy(3,j,1) = cmplx(real(exy(3,j,1))*at1,aimag(exy(3,j,1)))
   60    continue
         fxy(1,nx+1,1) = exy(1,nx+1,1)
         fxy(2,nx+1,1) = exy(2,nx+1,1)
         fxy(3,nx+1,1) = exy(3,nx+1,1)
c copy the electric fields
      else
         do 80 k = 1, nyh
         do 70 j = 1, nx1
         fxy(1,j,k) = exy(1,j,k)
         fxy(2,j,k) = exy(2,j,k)
         fxy(3,j,k) = exy(3,j,k)
   70    continue
   80    continue
      endif
      return
      end
      subroutine EMFIELD2T(fxyt,exyt,ffct,isign,nx,ny,nxvh,nyv,nxhd,nyhd
     1)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxyt, exyt, ffct
      dimension fxyt(nyv,3,nxvh), exyt(nyv,3,nxvh)
      dimension ffct(nyhd,nxhd)
c local data
      integer i, j, k, nxh, nyh, ny2, j1, k1
      real at1
      complex zero
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
      zero = cmplx(0.0,0.0)
c add the fields
      if (isign.gt.0) then
         do 30 j = 1, nxh
         do 20 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffct(k,j))
         do 10 i = 1, 3
         fxyt(k,i,j) = fxyt(k,i,j) + exyt(k,i,j)*at1
         fxyt(k1,i,j) = fxyt(k1,i,j) + exyt(k1,i,j)*at1
   10    continue
   20    continue
   30    continue
         k1 = nyh + 1
         do 50 j = 1, nxh
         at1 = aimag(ffct(1,j))
         do 40 i = 1, 3
         fxyt(1,i,j) = fxyt(1,i,j) + exyt(1,i,j)*at1
         fxyt(k1,i,j) = fxyt(k1,i,j) + exyt(k1,i,j)*at1
   40    continue
   50    continue
         do 70 k = 2, nyh
         k1 = ny2 - k
         do 60 i = 1, 3
         fxyt(k,i,j1) = zero
         fxyt(k1,i,j1) = zero
   60    continue
   70    continue
         k1 = nyh + 1
         do 80 i = 1, 3
         fxyt(1,i,j1) = zero
         fxyt(k1,i,j1) = zero
   80    continue
c copy the fields
      else if (isign.lt.0) then
         do 110 j = 1, nxh
         do 100 k = 2, nyh
         k1 = ny2 - k
         at1 = aimag(ffct(k,j))
         do 90 i = 1, 3
         fxyt(k,i,j) = exyt(k,i,j)*at1
         fxyt(k1,i,j) = exyt(k1,i,j)*at1
   90    continue
  100    continue
  110    continue
         k1 = nyh + 1
         do 130 j = 1, nxh
         at1 = aimag(ffct(1,j))
         do 120 i = 1, 3
         fxyt(1,i,j) = exyt(1,i,j)*at1
         fxyt(k1,i,j) = exyt(k1,i,j)*at1
  120    continue
  130    continue
         do 150 k = 2, nyh
         k1 = ny2 - k
         do 140 i = 1, 3
         fxyt(k,i,j1) = zero
         fxyt(k1,i,j1) = zero
  140    continue
  150    continue
         k1 = nyh + 1
         do 160 i = 1, 3
         fxyt(1,i,j1) = zero
         fxyt(k1,i,j1) = zero
  160    continue
      endif
      return
      end
      subroutine AVPOT23(bxy,axy,nx,ny,nxvh,nyv)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with periodic boundary conditions.
c input: bxy, nx, ny, nxvh, nyv, output: axy
c approximate flop count is: 38*nxc*nyc + 10*(nxc + nyc)
c and nxc*nyc divides, where nxc = nx/2 - 1, nyc = ny/2 - 1
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = ax(ky=pi) = ay(ky=pi) = az(ky=pi) 
c = 0, and ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c bxy(i,j,k) = i component of complex magnetic field
c axy(i,j,k) = i component of complex vector potential
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      complex bxy, axy, zero, zt1, zt2, zt3
      dimension bxy(3,nxvh,nyv), axy(3,nxvh,nyv)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate vector potential
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      axy(1,j,k) = at3*zt1
      axy(2,j,k) = -at2*zt1
      axy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(bxy(3,j,k1)),real(bxy(3,j,k1)))
      zt2 = cmplx(-aimag(bxy(2,j,k1)),real(bxy(2,j,k1)))
      zt3 = cmplx(-aimag(bxy(1,j,k1)),real(bxy(1,j,k1)))
      axy(1,j,k1) = -at3*zt1
      axy(2,j,k1) = -at2*zt1
      axy(3,j,k1) = at2*zt2 + at3*zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      zt1 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      axy(1,1,k) = at3*zt1
      axy(2,1,k) = zero
      axy(3,1,k) = -at3*zt3
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      at2 = 1.0/dkx
      zt1 = cmplx(-aimag(bxy(3,j,1)),real(bxy(3,j,1)))
      zt2 = cmplx(-aimag(bxy(2,j,1)),real(bxy(2,j,1)))
      axy(1,j,1) = zero
      axy(2,j,1) = -at2*zt1
      axy(3,j,1) = at2*zt2
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      axy(3,j,k1) = zero
   40 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      return
      end
      subroutine AVRPOT23(axy,bxy,ffc,ci,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c the radiative part of the vector potential
c with periodic boundary conditions.
c input: all, output: axy
c approximate flop count is: 68*nxc*nyc + 20*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the radiative vector potential is updated using the equations:
c ax(kx,ky) = (sqrt(-1)*ky*bz(kx,ky)
c                       - affp*ci2*cux(kx,ky)*s(kx,ky)/(kx*kx+ky*ky)
c ay(kx,ky) = -(sqrt(-1)*kx*bz(kx,ky)
c                       + affp*ci2*cuy(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = (sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
c                       - affp*ci2*cuz(kx,ky)*s(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, ci2 = ci*ci
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c j,k = fourier mode numbers, except for
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
c ax(kx=0,ky=0) = ay(kx=0,ky=0) = az(kx=0,ky=0) = 0.
c axy(i,j,k) = on entry, complex current density cu
c axy(i,j,k) = on exit, complex radiative vector potential
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffc(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci
      complex axy, bxy, ffc
      dimension axy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, afc2, dkx, dky, dky2, at1, at2
      complex zero, zt1, zt2, zt3
      if (ci.le.0.) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      afc2 = real(ffc(1,1))*ci*ci
      zero = cmplx(0.,0.)
c calculate the radiative vector potential
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = 1.0/(dkx*dkx + dky2)
      at2 = afc2*aimag(ffc(j,k))
c update radiative vector potential, ky > 0
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      axy(1,j,k) = at1*(dky*zt1 - at2*axy(1,j,k))
      axy(2,j,k) = -at1*(dkx*zt1 + at2*axy(2,j,k))
      axy(3,j,k) = at1*((dkx*zt2 - dky*zt3) - at2*axy(3,j,k))
c update radiative vector potential, ky < 0
      zt1 = cmplx(-aimag(bxy(3,j,k1)),real(bxy(3,j,k1)))
      zt2 = cmplx(-aimag(bxy(2,j,k1)),real(bxy(2,j,k1)))
      zt3 = cmplx(-aimag(bxy(1,j,k1)),real(bxy(1,j,k1)))
      axy(1,j,k1) = -at1*(dky*zt1 + at2*axy(1,j,k1))
      axy(2,j,k1) = -at1*(dkx*zt1 + at2*axy(2,j,k1))
      axy(3,j,k1) = at1*((dkx*zt2 + dky*zt3) - at2*axy(3,j,k1))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at1 = 1.0/(dky*dky)
      at2 = afc2*aimag(ffc(1,k))
c update radiative vector potential
      zt1 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      axy(1,1,k) = at1*(dky*zt1 - at2*axy(1,1,k))
      axy(2,1,k) = zero
      axy(3,1,k) = -at1*(dky*zt3 + at2*axy(3,1,k))
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      at1 = 1.0/(dkx*dkx)
      at2 = afc2*aimag(ffc(j,1))
c update radiative vector potential
      zt1 = cmplx(-aimag(bxy(3,j,1)),real(bxy(3,j,1)))
      zt2 = cmplx(-aimag(bxy(2,j,1)),real(bxy(2,j,1)))
      axy(1,j,1) = zero
      axy(2,j,1) = -at1*(dkx*zt1 + at2*axy(2,j,1))
      axy(3,j,1) = at1*(dkx*zt2 - at2*axy(3,j,1))
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      axy(3,j,k1) = zero
   40 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      return
      end
      subroutine GTMODES2(pot,pott,nx,ny,it,modesx,modesy,nxe,nye,nt2,mo
     1desxd,modesyd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c nxe = first dimension of input array pot, nxe >= nx
c nye = second dimension of input array pot, nye >= ny
c nt2 = first dimension of output array pott, nt2 >= 2*it
c modesxd = second dimension of output array pott, modesxd >= modesx
c modesyd = third dimension of output array pott, modesyd  = 2*modesy
      implicit none
      integer nx, ny, it, modesx, modesy, nxe, nye, nt2
      integer modesxd, modesyd
      real pot, pott
      dimension pot(nxe,nye), pott(nt2,modesxd,modesyd)
c local data
      integer i1, i2, nxh, nyh, kmax, jmax, ny2, j, k, j1, k1
      i2 = it + it
      i1 = i2 - 1
      if (i2.gt.nt2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pott(i1,j,2*k-2) = pot(2*j-1,k)
      pott(i2,j,2*k-2) = pot(2*j,k)
      pott(i1,j,2*k-1) = pot(2*j-1,k1)
      pott(i2,j,2*k-1) = pot(2*j,k1)
   10 continue
c mode numbers kx = 0, nx/2
      pott(i1,1,2*k-2) = pot(1,k)
      pott(i2,1,2*k-2) = pot(2,k)
      pott(i1,1,2*k-1) = pot(1,k)
      pott(i2,1,2*k-1) = -pot(2,k)
      if (modesx.gt.nxh) then
         pott(i1,j1,2*k-2) = pot(1,k1)
         pott(i2,j1,2*k-2) = -pot(2,k1)
         pott(i1,j1,2*k-1) = pot(1,k1)
         pott(i2,j1,2*k-1) = pot(2,k1)
      endif
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, jmax
      pott(i1,j,1) = pot(2*j-1,1)
      pott(i2,j,1) = pot(2*j,1)
   30 continue
      pott(i1,1,1) = pot(1,1)
      pott(i2,1,1) = 0.
      if (modesx.gt.nxh) then
         pott(i1,j1,1) = pot(2,1)
         pott(i2,j1,1) = 0.
      endif
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 40 j = 2, jmax
         pott(i1,j,ny) = pot(2*j-1,k1)
         pott(i2,j,ny) = pot(2*j,k1)
   40    continue
         pott(i1,1,ny) = pot(1,k1)
         pott(i2,1,ny) = 0.
         if (modesx.gt.nxh) then
            pott(i1,j1,ny) = pot(2,k1)
            pott(i2,j1,ny) = 0.
         endif
      endif
      return
      end
      subroutine PTMODES2(pot,pott,nx,ny,it,modesx,modesy,nxe,nye,nt2,mo
     1desxd,modesyd)
c this subroutine extracts lowest order modes from a location in a time
c history array pott and stores them into complex array pot
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c nxe = first dimension of input array pot, nxe >= nx
c nye = second dimension of input array pot, nye >= ny
c nt2 = first dimension of output array pott, nt2 >= 2*it
c modesxd = second dimension of output array pott, modesxd >= modesx
c modesyd = third dimension of output array pott, modesyd  = 2*modesy
      implicit none
      integer nx, ny, it, modesx, modesy, nxe, nye, nt2
      integer modesxd, modesyd
      real pot, pott
      dimension pot(nxe,nye), pott(nt2,modesxd,modesyd)
c local data
      integer i1, i2, nxh, nyh, kmax, jmax, ny2, j, k, j1, k1
      i2 = it + it
      i1 = i2 - 1
      if (i2.gt.nt2) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, kmax
      k1 = ny2 - k
      do 10 j = 2, jmax
      pot(2*j-1,k) = pott(i1,j,2*k-2)
      pot(2*j,k) = pott(i2,j,2*k-2)
      pot(2*j-1,k1) = pott(i1,j,2*k-1)
      pot(2*j,k1) = pott(i2,j,2*k-1)
   10 continue
      do 20 j = jmax+1, nxh
      pot(2*j-1,k) = 0.
      pot(2*j,k) = 0.
      pot(2*j-1,k1) = 0.
      pot(2*j,k1) = 0.
   20 continue
c mode numbers kx = 0, nx/2
      pot(1,k) = pott(i1,1,2*k-2)
      pot(2,k) = pott(i2,1,2*k-2)
      pot(1,k1) = 0.
      pot(2,k1) = 0.
      if (modesx.gt.nxh) then
         pot(1,k1) = pott(i1,j1,2*k-2)
         pot(2,k1) = -pott(i2,j1,2*k-2)
      endif
   30 continue
      do 50 k = kmax+1, nyh
      k1 = ny2 - k
      do 40 j = 1, nxh
      pot(2*j-1,k) = 0.
      pot(2*j,k) = 0.
      pot(2*j-1,k1) = 0.
      pot(2*j,k1) = 0.
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, jmax
      pot(2*j-1,1) = pott(i1,j,1)
      pot(2*j,1) = pott(i2,j,1)
      pot(2*j-1,k1) = 0.
      pot(2*j,k1) = 0.
   60 continue
      do 70 j = jmax+1, nxh
      pot(2*j-1,1) = 0.
      pot(2*j,1) = 0.
      pot(2*j-1,k1) = 0.
      pot(2*j,k1) = 0.
   70 continue
      pot(1,1) = pott(i1,1,1)
      pot(2,1) = 0.
      pot(1,k1) = 0.
      pot(2,k1) = 0.
      if (modesx.gt.nxh) then
         pot(2,1) = pott(i1,j1,1)
      endif
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 80 j = 2, jmax
         pot(2*j-1,k1) = pott(i1,j,ny)
         pot(2*j,k1) = pott(i2,j,ny)
   80    continue
         pot(1,k1) = pott(i1,1,ny)
         pot(2,k1) = 0.
         if (modesx.gt.nxh) then
            pot(2,k1) = pott(i1,j1,ny)
         endif
      endif
      return
      end
      subroutine GTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,nxvh,n
     1yv,nt,modesxd,modesyd)
c this subroutine extracts lowest order modes from complex vector array
c vpot and stores them into a location in a time history array vpott
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c nyv = third dimension of input array vpot, nyv >= ny
c nt = first dimension of output array vpott, nt >= it
c modesxd = second dimension of output array vpott, modesxd >= modesx
c modesyd = third dimension of output array vpott, modesyd  = 2*modesy
      implicit none
      integer nx, ny, it, modesx, modesy, ndim, nxvh, nyv, nt
      integer modesxd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nxvh,nyv), vpott(nt,ndim,modesxd,modesyd)
c local data
      integer nxh, nyh, kmax, jmax, ny2, i, j, k, j1, k1
      if (it.gt.nt) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 40 k = 2, kmax
      k1 = ny2 - k
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpott(it,i,j,2*k-2) = vpot(i,j,k)
      vpott(it,i,j,2*k-1) = vpot(i,j,k1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 i = 1, ndim
      vpott(it,i,1,2*k-2) = vpot(i,1,k)
      vpott(it,i,1,2*k-1) = conjg(vpot(i,1,k))
      if (modesx.gt.nxh) then
         vpott(it,i,j1,2*k-2) = conjg(vpot(i,1,k1))
         vpott(it,i,j1,2*k-1) = vpot(i,1,k1)
      endif
   30 continue
   40 continue
c mode numbers ky = 0, ny/2
      do 60 j = 2, jmax
      do 50 i = 1, ndim
      vpott(it,i,j,1) = vpot(i,j,1)
   50 continue
   60 continue
      do 70 i = 1, ndim
      vpott(it,i,1,1) = cmplx(real(vpot(i,1,1)),0.0)
      if (modesx.gt.nxh) then
         vpott(it,i,j1,1) = cmplx(aimag(vpot(i,1,1)),0.0)
      endif
   70 continue
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 90 j = 2, jmax
         do 80 i = 1, ndim
         vpott(it,i,j,ny) = vpot(i,j,k1)
   80    continue
   90    continue
         do 100 i = 1, ndim
         vpott(it,i,1,ny) = cmplx(real(vpot(i,1,k1)),0.0)
         if (modesx.gt.nxh) then
            vpott(it,i,j1,ny) = cmplx(aimag(vpot(i,1,k1)),0.0)
         endif
  100    continue
      endif
      return
      end
      subroutine PTVMODES2(vpot,vpott,nx,ny,it,modesx,modesy,ndim,nxvh,n
     1yv,nt,modesxd,modesyd)
c this subroutine extracts lowest order modes from a location in a time
c history array vpott and stores them into complex vector array vpot
c modes stored: kx=(0,1,...,NX/2), ky=(0,+-1,+-2,...,+-(NY/2-1),NY/2)
c nx/ny = system length in x/y direction
c it = current time
c modesx/modesy = number of modes to store in x/y direction,
c where modesx <= nx/2+1, modesy <= ny/2+1
c ndim = number of field arrays, must be >= 1
c nxvh = second dimension of input array vpot, nxvh >= nx/2
c nyv = third dimension of input array vpot, nyv >= ny
c nt = first dimension of output array vpott, nt >= it
c modesxd = second dimension of output array vpott, modesxd >= modesx
c modesyd = third dimension of output array vpott, modesyd  = 2*modesy
      implicit none
      integer nx, ny, it, modesx, modesy, ndim, nxvh, nyv, nt
      integer modesxd, modesyd
      complex vpot, vpott
      dimension vpot(ndim,nxvh,nyv), vpott(nt,ndim,modesxd,modesyd)
c local data
      integer nxh, nyh, kmax, jmax, ny2, i, j, k, j1, k1
      complex zero
      if (it.gt.nt) return
      nxh = nx/2
      nyh = max(1,ny/2)
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      ny2 = ny + 2
      jmax = min0(modesx,nxh)
      kmax = min0(modesy,nyh)
      j1 = nxh + 1
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 60 k = 2, kmax
      k1 = ny2 - k
      do 20 j = 2, jmax
      do 10 i = 1, ndim
      vpot(i,j,k) = vpott(it,i,j,2*k-2)
      vpot(i,j,k1) = vpott(it,i,j,2*k-1)
   10 continue
   20 continue
      do 40 j = jmax+1, nxh
      do 30 i = 1, ndim
      vpot(i,j,k) = zero
      vpot(i,j,k1) = zero
   30 continue
   40 continue
c mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      vpot(i,1,k) = vpott(it,i,1,2*k-2)
      vpot(i,1,k1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,k1) = conjg(vpott(it,i,j1,2*k-2))
      endif
   50 continue
   60 continue
      do 90 k = kmax+1, nyh
      k1 = ny2 - k
      do 80 j = 1, nxh
      do 70 i = 1, ndim
      vpot(i,j,k) = zero
      vpot(i,j,k1) = zero
   70 continue
   80 continue
   90 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 110 j = 2, jmax
      do 100 i = 1, ndim
      vpot(i,j,1) = vpott(it,i,j,1)
      vpot(i,j,k1) = zero
  100 continue
  110 continue
      do 130 j = jmax+1, nxh
      do 120 i = 1, ndim
      vpot(i,j,1) = zero
      vpot(i,j,k1) = zero
  120 continue
  130 continue
      do 140 i = 1, ndim
      vpot(i,1,1) = cmplx(real(vpott(it,i,1,1)),0.0)
      vpot(i,1,k1) = zero
      if (modesx.gt.nxh) then
         vpot(i,1,1) = cmplx(real(vpot(i,1,1)),real(vpott(it,i,j1,1)))
      endif
  140 continue
      if (modesy.gt.nyh) then
         k1 = nyh + 1
         do 160 j = 2, jmax
         do 150 i = 1, ndim
         vpot(i,j,k1) = vpott(it,i,j,ny)
  150    continue
  160    continue
         do 170 i = 1, ndim
         vpot(i,1,k1) = cmplx(real(vpott(it,i,1,ny)),0.0)
         if (modesx.gt.nxh) then
            vpot(i,1,k1) = cmplx(real(vpot(i,1,k1)),real(vpott(it,i,j1,n
     1y)))
         endif
  170    continue
      endif
      return
      end
      subroutine POYNT2(q,exy,bxy,ffc,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine calculates the momentum in the electromagnetic field
c given by the poynting flux.  inputs are the charge density, transverse
c electric field, and magnetic field.  outputs are sx, sy, sz
c equation used is:
c sx = sum((fy(j,k)+exy(2,j,k))*conjg(bxy(3,j,k))-exy(3,j,k)*
c conjg(bxy(2,j,k)))
c sy = sum(exy(3,j,k)*conjg(bxy(1,j,k))-(fx(j,k)+exy(1,j,k))*
c conjg(bxy(3,j,k)))
c sz = sum((fx(j,k)+exy(1,j,k))*conjg(bxy(2,j,k))-(fy(j,k)+exy(2,j,k))*
c conjg(bxy(1,j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c exy(1,j,k) = x component of transverse electric field,
c exy(2,j,k) = y component of transverse electric field,
c exy(3,j,k) = z component of transverse electric field,
c bxy(1,j,k) = x component of magnetic field,
c bxy(2,j,k) = y component of magnetic field,
c bxy(3,j,k) = z component of magnetic field,
c all for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c sx/sy/sz = x/y/z components of field momentum
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real sx, sy, sz
      complex q, exy, bxy, ffc
      dimension q(nxvh,nyv), exy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, affp, dky, at1, at2, at3
      complex zt1, zt2
      double precision wx, wy, wz
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      affp = real(ffc(1,1))
c calculate force/charge and sum field energy
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      at1 = real(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = at3*zt1 + exy(2,j,k)
      zt1 = at2*zt1 + exy(1,j,k)
      wx = wx + (zt2*conjg(bxy(3,j,k))-exy(3,j,k)*conjg(bxy(2,j,k)))
      wy = wy + (exy(3,j,k)*conjg(bxy(1,j,k))-zt1*conjg(bxy(3,j,k)))
      wz = wz + (zt1*conjg(bxy(2,j,k))-zt2*conjg(bxy(1,j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      zt1 = at2*zt2 + exy(1,j,k1)
      zt2 = -at3*zt2 + exy(2,j,k1)
      wx = wx + (zt2*conjg(bxy(3,j,k1))-exy(3,j,k1)*conjg(bxy(2,j,k1)))
      wy = wy + (exy(3,j,k1)*conjg(bxy(1,j,k1))-zt1*conjg(bxy(3,j,k1)))
      wz = wz + (zt1*conjg(bxy(2,j,k1))-zt2*conjg(bxy(1,j,k1)))
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      at1 = real(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      zt2 = at3*zt1 + exy(2,1,k)
      zt1 = exy(1,1,k)
      wx = wx + (zt2*conjg(bxy(3,1,k))-exy(3,1,k)*conjg(bxy(2,1,k)))
      wy = wy + (exy(3,1,k)*conjg(bxy(1,1,k))-zt1*conjg(bxy(3,1,k)))
      wz = wz + (zt1*conjg(bxy(2,1,k))-zt2*conjg(bxy(1,1,k)))
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nxh
      at1 = real(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      zt2 = exy(2,j,1)
      zt1 = at2*zt1 + exy(1,j,1)
      wx = wx + (zt2*conjg(bxy(3,j,1))-exy(3,j,1)*conjg(bxy(2,j,1)))
      wy = wy + (exy(3,j,1)*conjg(bxy(1,j,1))-zt1*conjg(bxy(3,j,1)))
      wz = wz + (zt1*conjg(bxy(2,j,1))-zt2*conjg(bxy(1,j,1)))
   40 continue
      at1 = 2.0*float(nx*ny)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
      subroutine DPOYNT2(q,cu,ffc,ci,sx,sy,sz,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine calculates the momentum in the darwin field given by
c the poynting flux in 2-1/2d.  inputs are the charge density and
c current density.  outputs are sx, sy, sz
c equation used is:
c sx = sum(fy(j,k)*conjg(bz(j,k))
c sy = sum(-fx(j,k)*conjg(bz(j,k)))
c sz = sum(fx(j,k)*conjg(by(j,k))-fy(j,k)*conjg(bx(j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky)
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c all for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c sx/sy/sz = x/y/z components of field momentum
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, sx, sy, sz
      complex q, cu, ffc
      dimension q(nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, affp, dky, at1, at2, at3, at4, at5
      complex zt1, zt2, zt3, zt4, zt5
      double precision wx, wy, wz
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      ci2 = ci*ci
      affp = real(ffc(1,1))
c sum darwin field momentum
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      at1 = real(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt4 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt5 = at3*zt4
      zt4 = at2*zt4
      at4 = ci2*at2
      at5 = ci2*at3
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      zt3 = conjg(at4*zt2 - at5*zt3)
      zt2 = conjg(at4*zt1)
      zt1 = conjg(at5*zt1)
      wx = wx + zt5*zt3
      wy = wy - zt4*zt3
      wz = wz - (zt4*zt2+zt5*zt1)
      zt5 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      zt4 = at2*zt5
      zt5 = at3*zt5
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      zt3 = conjg(at4*zt2 + at5*zt3)
      zt2 = conjg(at4*zt1)
      zt1 = conjg(at5*zt1)
      wx = wx - zt5*zt3
      wy = wy - zt4*zt3
      wz = wz - (zt4*zt2+zt5*zt1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      at1 = real(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt4 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      zt5 = at3*zt4
      at3 = ci2*at3
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      zt1 = conjg(at3*zt1)
      zt3 = conjg(at3*zt3)
      wx = wx - zt5*zt3
      wz = wz - zt5*zt1
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nxh
      at1 = real(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt4 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      zt4 = at2*zt4
      at2 = ci2*at2
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      zt3 = conjg(at2*zt2)
      zt2 = conjg(at2*zt1)
      wy = wy - zt4*zt3
      wz = wz - zt4*zt2
   40 continue
      at1 = 2.0*float(nx*ny)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
      subroutine DPOYNT22(q,cu,ffc,ci,sx,sy,nx,ny,nxvh,nyv,nxhd,nyhd)
c this subroutine calculates the momentum in the darwin field given by
c the poynting flux in 2d.  inputs are the charge density and
c current density.  outputs are sx, sy
c equation used is:
c sx = sum(fy(j,k)*conjg(bz(j,k))
c sy = sum(-fx(j,k)*conjg(bz(j,k))), where
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))
c kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky)
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c all for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c sx/sy = x/y components of field momentum
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, sx, sy
      complex q, cu, ffc
      dimension q(nxvh,nyv), cu(2,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, affp, dky, at1, at2, at3, at4, at5
      complex zt2, zt3, zt4, zt5
      double precision wx, wy
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      ci2 = ci*ci
      affp = real(ffc(1,1))
c sum darwin field momentum
      wx = 0.0d0
      wy = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nxh
      at1 = real(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt4 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt5 = at3*zt4
      zt4 = at2*zt4
      at4 = ci2*at2
      at5 = ci2*at3
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      zt3 = conjg(at4*zt2 - at5*zt3)
      wx = wx + zt5*zt3
      wy = wy - zt4*zt3
      zt5 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      zt4 = at2*zt5
      zt5 = at3*zt5
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      zt3 = conjg(at4*zt2 + at5*zt3)
      wx = wx - zt5*zt3
      wy = wy - zt4*zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      at1 = real(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt4 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      zt5 = at3*zt4
      at3 = ci2*at3
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      zt3 = conjg(at3*zt3)
      wx = wx - zt5*zt3
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nxh
      at1 = real(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt4 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      zt4 = at2*zt4
      at2 = ci2*at2
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      zt3 = conjg(at2*zt2)
      wy = wy - zt4*zt3
   40 continue
      at1 = 2.0*float(nx*ny)/affp
      sx = at1*wx
      sy = at1*wy
      return
      end
      subroutine SCFGUARD2(cus,cu,q2m0,nx,ny,nxe,nye)
c initialize extended periodic field with scaled field
c quadratic interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(3,nxe,nye), cu(3,nxe,nye)
      integer i, j, k
c initialize extended field, with zero in the edges
      do 40 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 3
      cus(i,j+1,k+1) = -q2m0*cu(i,j+1,k+1)
   10 continue
   20 continue
      do 30 i = 1, 3
      cus(i,1,k+1) = 0.
      cus(i,nx+2,k+1) = 0.
      cus(i,nx+3,k+1) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx
      do 50 i = 1, 3
      cus(i,j+1,1) = 0.
      cus(i,j+1,ny+2) = 0.
      cus(i,j+1,ny+3) = 0.
   50 continue
   60 continue
      do 70 i = 1, 3
      cus(i,1,1) = 0.
      cus(i,nx+2,1) = 0.
      cus(i,nx+3,1) = 0.
      cus(i,1,ny+2) = 0.
      cus(i,nx+2,ny+2) = 0.
      cus(i,nx+3,ny+2) = 0.
      cus(i,1,ny+3) = 0.
      cus(i,nx+2,ny+3) = 0.
      cus(i,nx+3,ny+3) = 0.
   70 continue
      return
      end
      subroutine SCFGUARD22(cus,cu,q2m0,nx,ny,nxe,nye)
c initialize extended periodic field with scaled field
c quadratic interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(2,nxe,nye), cu(2,nxe,nye)
      integer i, j, k
c initialize extended field, with zero in the edges
      do 40 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 2
      cus(i,j+1,k+1) = -q2m0*cu(i,j+1,k+1)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,1,k+1) = 0.
      cus(i,nx+2,k+1) = 0.
      cus(i,nx+3,k+1) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx
      do 50 i = 1, 2
      cus(i,j+1,1) = 0.
      cus(i,j+1,ny+2) = 0.
      cus(i,j+1,ny+3) = 0.
   50 continue
   60 continue
      do 70 i = 1, 2
      cus(i,1,1) = 0.
      cus(i,nx+2,1) = 0.
      cus(i,nx+3,1) = 0.
      cus(i,1,ny+2) = 0.
      cus(i,nx+2,ny+2) = 0.
      cus(i,nx+3,ny+2) = 0.
      cus(i,1,ny+3) = 0.
      cus(i,nx+2,ny+3) = 0.
      cus(i,nx+3,ny+3) = 0.
   70 continue
      return
      end
      subroutine SMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
c initialize extended periodic field
c quadratic interpolation
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nx, ny, nxe, nye
      dimension amu(4,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      amu(1,j+1,k+1) = x2y2m0
      amu(2,j+1,k+1) = xym0
      amu(3,j+1,k+1) = zxm0
      amu(4,j+1,k+1) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,1,k+1) = 0.
      amu(i,nx+2,k+1) = 0.
      amu(i,nx+3,k+1) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 4
      amu(i,j+1,1) = 0.
      amu(i,j+1,ny+2) = 0.
      amu(i,j+1,ny+3) = 0.
   40 continue
   50 continue
      do 60 i = 1, 4
      amu(i,1,1) = 0.
      amu(i,nx+2,1) = 0.
      amu(i,nx+3,1) = 0.
      amu(i,1,ny+2) = 0.
      amu(i,nx+2,ny+2) = 0.
      amu(i,nx+3,ny+2) = 0.
      amu(i,1,ny+3) = 0.
      amu(i,nx+2,ny+3) = 0.
      amu(i,nx+3,ny+3) = 0.
   60 continue
      return
      end
      subroutine SMCGUARD22(amu,x2y2m0,xym0,nx,ny,nxe,nye)
c initialize extended periodic field
c quadratic interpolation
      implicit none
      real amu, x2y2m0, xym0
      integer nx, ny, nxe, nye
      dimension amu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      amu(1,j+1,k+1) = x2y2m0
      amu(2,j+1,k+1) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,1,k+1) = 0.
      amu(i,nx+2,k+1) = 0.
      amu(i,nx+3,k+1) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 2
      amu(i,j+1,1) = 0.
      amu(i,j+1,ny+2) = 0.
      amu(i,j+1,ny+3) = 0.
   40 continue
   50 continue
      do 60 i = 1, 2
      amu(i,1,1) = 0.
      amu(i,nx+2,1) = 0.
      amu(i,nx+3,1) = 0.
      amu(i,1,ny+2) = 0.
      amu(i,nx+2,ny+2) = 0.
      amu(i,nx+3,ny+2) = 0.
      amu(i,1,ny+3) = 0.
      amu(i,nx+2,ny+3) = 0.
      amu(i,nx+3,ny+3) = 0.
   60 continue
      return
      end
      subroutine AMCGUARD2(amu,nx,ny,nxe,nye,ndim)
c accumulate extended periodic tensor field
c quadratic interpolation
c handles special case of ny=1 properly
      implicit none
      real amu
      integer nx, ny, nxe, nye, ndim
      dimension amu(ndim,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, ndim
      amu(i,2,k+1) = amu(i,2,k+1) + amu(i,nx+2,k+1)
      amu(i,3,k+1) = amu(i,3,k+1) + amu(i,nx+3,k+1)
      amu(i,nx+1,k+1) = amu(i,nx+1,k+1) + amu(i,1,k+1)
      amu(i,nx+2,k+1) = 0.0
      amu(i,nx+3,k+1) = 0.0
      amu(i,1,k+1) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx+3
      do 30 i = 1, ndim
      amu(i,j,3) = amu(i,j,3) + amu(i,j,ny+3)
      amu(i,j,2) = amu(i,j,2) + amu(i,j,ny+2)
      amu(i,j,ny+1) = amu(i,j,ny+1) + amu(i,j,1)
      amu(i,j,ny+3) = 0.0
      amu(i,j,ny+2) = 0.0
      amu(i,j,1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, ndim
      amu(i,2,3) = amu(i,2,3) + amu(i,nx+2,3)
      amu(i,3,3) = amu(i,3,3) + amu(i,nx+3,3)
      amu(i,nx+1,3) = amu(i,nx+1,3) + amu(i,1,3)
      amu(i,nx+2,3) = 0.0
      amu(i,nx+3,3) = 0.0
      amu(i,1,3) = 0.0
      amu(i,2,2) = amu(i,2,2) + amu(i,nx+2,2)
      amu(i,3,2) = amu(i,3,2) + amu(i,nx+3,2)
      amu(i,nx+1,2) = amu(i,nx+1,2) + amu(i,1,2)
      amu(i,nx+2,2) = 0.0
      amu(i,nx+3,2) = 0.0
      amu(i,1,2) = 0.0
      amu(i,2,ny+1) = amu(i,2,ny+1) + amu(i,nx+2,ny+1)
      amu(i,3,ny+1) = amu(i,3,ny+1) + amu(i,nx+3,ny+1)
      amu(i,nx+1,ny+1) = amu(i,nx+1,ny+1) + amu(i,1,ny+1)
      amu(i,nx+2,ny+1) = 0.0
      amu(i,nx+3,ny+1) = 0.0
      amu(i,1,ny+1) = 0.0
   50 continue
      return
      end
      subroutine SCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
c initialize extended periodic field with scaled field
c linear interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(3,nxe,nye), cu(3,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 40 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 3
      cus(i,j,k) = -q2m0*cu(i,j,k)
   10 continue
   20 continue
      do 30 i = 1, 3
      cus(i,nx+1,k) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx
      do 50 i = 1, 3
      cus(i,j,ny+1) = 0.
   50 continue
   60 continue
      do 70 i = 1, 3
      cus(i,nx+1,ny+1) = 0.
   70 continue
      return
      end
      subroutine SCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
c initialize extended periodic field with scaled field
c linear interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(2,nxe,nye), cu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 40 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 2
      cus(i,j,k) = -q2m0*cu(i,j,k)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,nx+1,k) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx
      do 50 i = 1, 2
      cus(i,j,ny+1) = 0.
   50 continue
   60 continue
      do 70 i = 1, 2
      cus(i,nx+1,ny+1) = 0.
   70 continue
      return
      end
      subroutine SMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
c initialize extended periodic field
c linear interpolation
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nx, ny, nxe, nye
      dimension amu(4,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      amu(1,j,k) = x2y2m0
      amu(2,j,k) = xym0
      amu(3,j,k) = zxm0
      amu(4,j,k) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,nx+1,k) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 4
      amu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 i = 1, 4
      amu(i,nx+1,ny+1) = 0.
   60 continue
      return
      end
      subroutine SMCGUARD22L(amu,x2y2m0,xym0,nx,ny,nxe,nye)
c initialize extended periodic field
c linear interpolation
      implicit none
      real amu, x2y2m0, xym0
      integer nx, ny, nxe, nye
      dimension amu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      amu(1,j,k) = x2y2m0
      amu(2,j,k) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,nx+1,k) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 2
      amu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 i = 1, 2
      amu(i,nx+1,ny+1) = 0.
   60 continue
      return
      end
      subroutine AMCGUARD2L(amu,nx,ny,nxe,nye,ndim)
c accumulate extended periodic tensor field
c linear interpolation
      implicit none
      real amu
      integer nx, ny, nxe, nye, ndim
      dimension amu(ndim,nxe,nye)
c local data
      integer i, j, k
c accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, ndim
      amu(i,1,k) = amu(i,1,k) + amu(i,nx+1,k)
      amu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      amu(i,j,1) = amu(i,j,1) + amu(i,j,ny+1)
      amu(i,j,ny+1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, ndim
      amu(i,1,1) = amu(i,1,1) + amu(i,nx+1,ny+1)
      amu(i,nx+1,ny+1) = 0.0
   50 continue
      return
      end
      subroutine SMCGUARD2C(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,nxe,nye)
c initialize extended periodic tensor field
c cubic interpolation
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nx, ny, nxe, nye
      dimension amu(4,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      amu(1,j+2,k+2) = x2y2m0
      amu(2,j+2,k+2) = xym0
      amu(3,j+2,k+2) = zxm0
      amu(4,j+2,k+2) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,1,k+2) = 0.
      amu(i,2,k+2) = 0.
      amu(i,nx+3,k+2) = 0.
      amu(i,nx+4,k+2) = 0.
      amu(i,nx+5,k+2) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 4
      amu(i,j+2,1) = 0.
      amu(i,j+2,2) = 0.
      amu(i,j+2,ny+3) = 0.
      amu(i,j+2,ny+4) = 0.
      amu(i,j+2,ny+5) = 0.
   40 continue
   50 continue
      do 60 i = 1, 4
      amu(i,1,1) = 0.
      amu(i,2,1) = 0.
      amu(i,nx+3,1) = 0.
      amu(i,nx+4,1) = 0.
      amu(i,nx+5,1) = 0.
      amu(i,1,2) = 0.
      amu(i,2,2) = 0.
      amu(i,nx+3,2) = 0.
      amu(i,nx+4,2) = 0.
      amu(i,nx+5,2) = 0.
      amu(i,1,ny+3) = 0.
      amu(i,2,ny+3) = 0.
      amu(i,nx+3,ny+3) = 0.
      amu(i,nx+4,ny+3) = 0.
      amu(i,nx+5,ny+3) = 0.
      amu(i,1,ny+4) = 0.
      amu(i,2,ny+4) = 0.
      amu(i,nx+3,ny+4) = 0.
      amu(i,nx+4,ny+4) = 0.
      amu(i,nx+5,ny+4) = 0.
      amu(i,1,ny+5) = 0.
      amu(i,2,ny+5) = 0.
      amu(i,nx+3,ny+5) = 0.
      amu(i,nx+4,ny+5) = 0.
      amu(i,nx+5,ny+5) = 0.
   60 continue
      return
      end
      subroutine SMCGUARD22C(amu,x2y2m0,xym0,nx,ny,nxe,nye)
c initialize extended periodic tensor field
c cubic interpolation
      implicit none
      real amu, x2y2m0, xym0
      integer nx, ny, nxe, nye
      dimension amu(2,nxe,nye)
c local data
      integer i, j, k
c initialize extended field, with zero in the edges
      do 30 k = 1, ny
      do 10 j = 1, nx
      amu(1,j+2,k+2) = x2y2m0
      amu(2,j+2,k+2) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,1,k+2) = 0.
      amu(i,2,k+2) = 0.
      amu(i,nx+3,k+2) = 0.
      amu(i,nx+4,k+2) = 0.
      amu(i,nx+5,k+2) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 2
      amu(i,j+2,1) = 0.
      amu(i,j+2,2) = 0.
      amu(i,j+2,ny+3) = 0.
      amu(i,j+2,ny+4) = 0.
      amu(i,j+2,ny+5) = 0.
   40 continue
   50 continue
      do 60 i = 1, 2
      amu(i,1,1) = 0.
      amu(i,2,1) = 0.
      amu(i,nx+3,1) = 0.
      amu(i,nx+4,1) = 0.
      amu(i,nx+5,1) = 0.
      amu(i,1,2) = 0.
      amu(i,2,2) = 0.
      amu(i,nx+3,2) = 0.
      amu(i,nx+4,2) = 0.
      amu(i,nx+5,2) = 0.
      amu(i,1,ny+3) = 0.
      amu(i,2,ny+3) = 0.
      amu(i,nx+3,ny+3) = 0.
      amu(i,nx+4,ny+3) = 0.
      amu(i,nx+5,ny+3) = 0.
      amu(i,1,ny+4) = 0.
      amu(i,2,ny+4) = 0.
      amu(i,nx+3,ny+4) = 0.
      amu(i,nx+4,ny+4) = 0.
      amu(i,nx+5,ny+4) = 0.
      amu(i,1,ny+5) = 0.
      amu(i,2,ny+5) = 0.
      amu(i,nx+3,ny+5) = 0.
      amu(i,nx+4,ny+5) = 0.
      amu(i,nx+5,ny+5) = 0.
   60 continue
      return
      end
      subroutine DCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2-1/2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c amu(3,j,k) = zx component of complex momentum flux
c amu(4,j,k) = zy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = second dimension of field arrays, must be >= nxh
c nyv = third dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k)),-real(amu(3,j,k)))
      zt2 = cmplx(aimag(amu(4,j,k)),-real(amu(4,j,k)))
      dcu(3,j,k) = dkx*zt1 + dky*zt2
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k1)),-real(amu(3,j,k1)))
      zt2 = cmplx(aimag(amu(4,j,k1)),-real(amu(4,j,k1)))
      dcu(3,j,k1) = dkx*zt1 - dky*zt2
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dky*zt2
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dkx*zt1
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   40 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(3,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
      subroutine DCUPERP22(dcu,amu,nx,ny,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex dcu, amu
      dimension dcu(2,nxvh,nyv), amu(2,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dkxy*zt1 - dkxy2*zt2)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dky*zt2
      dcu(2,1,k) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dkx*zt2
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
   40 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      return
      end
      subroutine ADCUPERP23(dcu,amu,nx,ny,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2-1/2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = dcu(3,kx,ky)-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on input:
c dcu(i,j,k) = complex acceleration density for fourier mode (j-1,k-1)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c amu(3,j,k) = zx component of complex momentum flux
c amu(4,j,k) = zy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex dcu, amu
      dimension dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dky*dcu(1,j,k) - dkx*dcu(2,j,k) + dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k)),-real(amu(3,j,k)))
      zt2 = cmplx(aimag(amu(4,j,k)),-real(amu(4,j,k)))
      dcu(3,j,k) = dcu(3,j,k) + dkx*zt1 + dky*zt2
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dky*dcu(1,j,k1) + dkx*dcu(2,j,k1) + dkxy*zt1 - dkxy2*zt
     12)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
      zt1 = cmplx(aimag(amu(3,j,k1)),-real(amu(3,j,k1)))
      zt2 = cmplx(aimag(amu(4,j,k1)),-real(amu(4,j,k1)))
      dcu(3,j,k1) = dcu(3,j,k1) + dkx*zt1 - dky*zt2
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dcu(1,1,k) + dky*zt2
      dcu(2,1,k) = zero
      zt2 = cmplx(aimag(amu(4,1,k)),-real(amu(4,1,k)))
      dcu(3,1,k) = dcu(3,1,k) + dky*zt2
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dcu(2,j,1) + dkx*zt2
      zt1 = cmplx(aimag(amu(3,j,1)),-real(amu(3,j,1)))
      dcu(3,j,1) = dcu(3,j,1) + dkx*zt1
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
      dcu(3,j,k1) = zero
   40 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
      return
      end
      subroutine ADCUPERP22(dcu,amu,nx,ny,nxvh,nyv)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2d with periodic boundary conditions.
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for dcu(i,kx=pi) = dcu(i,ky=pi) =  dcu(i,kx=0,ky=0) = 0.
c the transverse part is calculated using the equation:
c dcu(1,kx,ky) = dcu(1,kx,ky)-kx*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c dcu(2,kx,ky) = dcu(2,kx,ky)-ky*(kx*dcu(1,kx,ky)+ky*dcu(2,kx,ky))/
c               (kx*kx+ky*ky)
c on input:
c dcu(i,j,k) = complex acceleration density for fourier mode (j-1,k-1)
c on output:
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c amu(1,j,k) = xx component of complex momentum flux
c amu(2,j,k) = xy component of complex momentum flux
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex dcu, amu
      dimension dcu(2,nxvh,nyv), amu(2,nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2, at1
      complex zero, zt1, zt2, zt3
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*float(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      zt1 = cmplx(aimag(amu(1,j,k)),-real(amu(1,j,k)))
      zt2 = cmplx(aimag(amu(2,j,k)),-real(amu(2,j,k)))
      zt3 = at1*(dky*dcu(1,j,k) - dkx*dcu(2,j,k) + dkxy*zt1 + dkxy2*zt2)
      dcu(1,j,k) = dky*zt3
      dcu(2,j,k) = -dkx*zt3
      zt1 = cmplx(aimag(amu(1,j,k1)),-real(amu(1,j,k1)))
      zt2 = cmplx(aimag(amu(2,j,k1)),-real(amu(2,j,k1)))
      zt3 = at1*(dky*dcu(1,j,k1) + dkx*dcu(2,j,k1) + dkxy*zt1 - dkxy2*zt
     12)
      dcu(1,j,k1) = dky*zt3
      dcu(2,j,k1) = dkx*zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      zt2 = cmplx(aimag(amu(2,1,k)),-real(amu(2,1,k)))
      dcu(1,1,k) = dcu(1,1,k) + dky*zt2
      dcu(2,1,k) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*float(j - 1)
      zt2 = cmplx(aimag(amu(2,j,1)),-real(amu(2,j,1)))
      dcu(1,j,1) = zero
      dcu(2,j,1) = dcu(2,j,1) + dkx*zt2
      dcu(1,j,k1) = zero
      dcu(2,j,k1) = zero
   40 continue
      dcu(1,1,1) = zero
      dcu(2,1,1) = zero
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      return
      end
      subroutine EPOIS23(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,nx
     1vh,nyv,nxhd,nyhd,yee)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,nxvh,nyhd, output:ffe
c for isign /= 0, input: dcu,ffe,isign,ci,nx,ny,nxvh,nyv,nxhd,nyhd,
c output: exy,wf
c approximate flop count is: 68*nxc*nyc + 33*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2+wp0*ci2*s(kx,ky)**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = ex(ky=pi) = ey(ky=pi) = ez(ky=pi) 
c = 0, and ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c exy(1,j,k) = x component of complex transverse electric field
c exy(2,j,k) = y component of complex transverse electric field
c exy(3,j,k) = z component of complex transverse electric field
c all for fourier mode (j-1,k-1)
c aimag(ffe(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffe(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprical of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
c    |dcu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
!*****
!This has been modefied to partially accomodate the yee cell, but
!the changes have been commented out. Thre are some corrections
!that need to be made to at2 if we want to do this correctly.
!*****
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, wp0, ci, wf, yee
      complex dcu, exy, ffe
      dimension dcu(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffe(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp
!      real thetax, thetay, kxop, kyop
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      wpc = wp0*ci2
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffe(j,k) = cmplx(affp,1.)
      else
         ffe(j,k) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
c calculate smoothed transverse electric field and sum field energy
   30 if (isign.gt.0) go to 80
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      thetay = dky/2.0
!      kyop = 2.0*sin(thetay)
      do 40 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      thetax = dkx/2.0
!      kxop = 2.0*sin(thetax)
      at2 = -ci2*real(ffe(j,k))
!      at2 = at2*(dkx*dkx + dky*dky)/(kxop*kxop + kyop*kyop)
      at1 = at2*aimag(ffe(j,k))
      at2 = at2*at2
      exy(1,j,k) = at1*dcu(1,j,k)
      exy(2,j,k) = at1*dcu(2,j,k)
      exy(3,j,k) = at1*dcu(3,j,k)
      exy(1,j,k1) = at1*dcu(1,j,k1)
      exy(2,j,k1) = at1*dcu(2,j,k1)
      exy(3,j,k1) = at1*dcu(3,j,k1)
      wp = wp + at2*(dcu(1,j,k)*conjg(dcu(1,j,k)) + dcu(2,j,k)*conjg(dcu
     1(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k)) + dcu(1,j,k1)*conjg(dcu(1,
     2j,k1)) + dcu(2,j,k1)*conjg(dcu(2,j,k1)) + dcu(3,j,k1)*conjg(dcu(3,
     3j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      thetay = dky/2.0
!      kyop = 2.0*sin(thetay)
      at2 = -ci2*real(ffe(1,k))
!      at2 = at2*(dky*dky)/(kyop*kyop)!kxop = dkx = 0
      at1 = at2*aimag(ffe(1,k))
      at2 = at2*at2
      exy(1,1,k) = at1*dcu(1,1,k)
      exy(2,1,k) = at1*dcu(2,1,k)
      exy(3,1,k) = at1*dcu(3,1,k)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wp = wp + at2*(dcu(1,1,k)*conjg(dcu(1,1,k)) + dcu(2,1,k)*conjg(dcu
     1(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      thetax = dkx/2.0
!      kxop = 2.0*sin(thetax)
      at2 = -ci2*real(ffe(j,1))
!      at2 = at2*(dkx*dkx)/(kxop*kxop)!kyop = dky = 0
      at1 = at2*aimag(ffe(j,1))
      at2 = at2*at2
      exy(1,j,1) = at1*dcu(1,j,1)
      exy(2,j,1) = at1*dcu(2,j,1)
      exy(3,j,1) = at1*dcu(3,j,1)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at2*(dcu(1,j,1)*conjg(dcu(1,j,1)) + dcu(2,j,1)*conjg(dcu
     1(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
   70 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*wp/real(ffe(1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   80 wp = 0.0d0

!      if(yee==1) then
!c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
!      do 101 k = 2, nyh
!      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      thetay = dky/2.0
!      kyop = 2.0*sin(thetay)
!      do 91 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      thetax = dkx/2.0
!      kxop = 2.0*sin(thetax)
!      at2 = -ci2*real(ffe(j,k))
!      at2 = at2*(dkx*dkx + dky*dky)/(kxop*kxop + kyop*kyop)
!      at1 = at2*at2
!      exy(1,j,k) = at2*dcu(1,j,k)
!      exy(2,j,k) = at2*dcu(2,j,k)
!      exy(3,j,k) = at2*dcu(3,j,k)
!      exy(1,j,k1) = at2*dcu(1,j,k1)
!      exy(2,j,k1) = at2*dcu(2,j,k1)
!      exy(3,j,k1) = at2*dcu(3,j,k1)
!      wp = wp + at1*(dcu(1,j,k)*conjg(dcu(1,j,k)) + dcu(2,j,k)*conjg(dcu
!     1(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k)) + dcu(1,j,k1)*conjg(dcu(1,
!     2j,k1)) + dcu(2,j,k1)*conjg(dcu(2,j,k1)) + dcu(3,j,k1)*conjg(dcu(3,
!     3j,k1)))
!   91 continue
!  101 continue
!c mode numbers kx = 0, nx/2
!cdir$ ivdep
!      do 111 k = 2, nyh
!      k1 = ny2 - k
!      dky = dny*float(k - 1)
!      thetay = dky/2.0
!      kyop = 2.0*sin(thetay)
!      at2 = -ci2*real(ffe(1,k))
!      at2 = at2*(dky*dky)/(kyop*kyop)!kxop = dkx = 0
!      at1 = at2*at2
!      exy(1,1,k) = at2*dcu(1,1,k)
!      exy(2,1,k) = at2*dcu(2,1,k)
!      exy(3,1,k) = at2*dcu(3,1,k)
!      exy(1,1,k1) = zero
!      exy(2,1,k1) = zero
!      exy(3,1,k1) = zero
!      wp = wp + at1*(dcu(1,1,k)*conjg(dcu(1,1,k)) + dcu(2,1,k)*conjg(dcu
!     1(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
!  111 continue
!c mode numbers ky = 0, ny/2
!      k1 = nyh + 1
!      do 121 j = 2, nxh
!      dkx = dnx*float(j - 1)
!      thetax = dkx/2.0
!      kxop = 2.0*sin(thetax)
!      at2 = -ci2*real(ffe(j,1))
!      at2 = at2*(dkx*dkx)/(kxop*kxop)!kyop = dky = 0
!      at1 = at2*at2
!      exy(1,j,1) = at2*dcu(1,j,1)
!      exy(2,j,1) = at2*dcu(2,j,1)
!      exy(3,j,1) = at2*dcu(3,j,1)
!      exy(1,j,k1) = zero
!      exy(2,j,k1) = zero
!      exy(3,j,k1) = zero
!      wp = wp + at1*(dcu(1,j,1)*conjg(dcu(1,j,1)) + dcu(2,j,1)*conjg(dcu
!     1(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
!  121 continue
!      exy(1,1,1) = zero
!      exy(2,1,1) = zero
!      exy(3,1,1) = zero
!      exy(1,1,k1) = zero
!      exy(2,1,k1) = zero
!      exy(3,1,k1) = zero
!      wf = float(nx*ny)*wp/real(ffe(1,1))

!      else

c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*at2
      exy(1,j,k) = at2*dcu(1,j,k)
      exy(2,j,k) = at2*dcu(2,j,k)
      exy(3,j,k) = at2*dcu(3,j,k)
      exy(1,j,k1) = at2*dcu(1,j,k1)
      exy(2,j,k1) = at2*dcu(2,j,k1)
      exy(3,j,k1) = at2*dcu(3,j,k1)
      wp = wp + at1*(dcu(1,j,k)*conjg(dcu(1,j,k)) + dcu(2,j,k)*conjg(dcu
     1(2,j,k)) + dcu(3,j,k)*conjg(dcu(3,j,k)) + dcu(1,j,k1)*conjg(dcu(1,
     2j,k1)) + dcu(2,j,k1)*conjg(dcu(2,j,k1)) + dcu(3,j,k1)*conjg(dcu(3,
     3j,k1)))
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*at2
      exy(1,1,k) = at2*dcu(1,1,k)
      exy(2,1,k) = at2*dcu(2,1,k)
      exy(3,1,k) = at2*dcu(3,1,k)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wp = wp + at1*(dcu(1,1,k)*conjg(dcu(1,1,k)) + dcu(2,1,k)*conjg(dcu
     1(2,1,k)) + dcu(3,1,k)*conjg(dcu(3,1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*at2
      exy(1,j,1) = at2*dcu(1,j,1)
      exy(2,j,1) = at2*dcu(2,j,1)
      exy(3,j,1) = at2*dcu(3,j,1)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at1*(dcu(1,j,1)*conjg(dcu(1,j,1)) + dcu(2,j,1)*conjg(dcu
     1(2,j,1)) + dcu(3,j,1)*conjg(dcu(3,j,1)))
  120 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*wp/real(ffe(1,1))

!      endif
!      print*, maxval(abs(exy(1,:,:)))
!      print*, maxval(abs(exy(2,:,:)))
!      print*, maxval(abs(exy(2,:,:)))
      return
      end
      subroutine EPOIS22(dcu,exy,isign,ffe,ax,ay,affp,wp0,ci,wf,nx,ny,nx
     1vh,nyv,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with periodic boundary conditions.
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,nxvh,nyhd, output:ffe
c for isign =/ 0, input: dcu,ffe,isign,ci,nx,ny,nxvh,nyv,nxhd,nyhd,
c output: exy,wf
c approximate flop count is: 48*nxc*nyc + 24*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2+wp0*ci2*s(kx,ky)**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c ex(kx=pi) = ey(kx=pi) = ex(ky=pi) = ey(ky=pi)  = 0.0, and
c ex(kx=0,ky=0) = ey(kx=0,ky=0) = 0.
c if isign = 1, unsmoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)
c dcu(i,j,k) = transverse part of complex derivative of current for
c fourier mode (j-1,k-1)
c exy(1,j,k) = x component of complex transverse electric field
c exy(2,j,k) = y component of complex transverse electric field
c all for fourier mode (j-1,k-1)
c aimag(ffe(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffe(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c wp0 = normalized total plasma frequency squared
c ci = reciprical of velocity of light
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((affp/((kx**2+ky**2)*ci*ci)**2)
c    |dcu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the derivative of current is
c divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, ffe
      dimension dcu(2,nxvh,nyv), exy(2,nxvh,nyv)
      dimension ffe(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      wpc = wp0*ci2
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffe(j,k) = cmplx(affp,1.)
      else
         ffe(j,k) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
c calculate smoothed transverse electric field and sum field energy
   30 if (isign.gt.0) go to 80
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      do 40 j = 2, nxh
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*aimag(ffe(j,k))
      at2 = at2*at2
      exy(1,j,k) = at1*dcu(1,j,k)
      exy(2,j,k) = at1*dcu(2,j,k)
      exy(1,j,k1) = at1*dcu(1,j,k1)
      exy(2,j,k1) = at1*dcu(2,j,k1)
      wp = wp + at2*(dcu(1,j,k)*conjg(dcu(1,j,k)) + dcu(2,j,k)*conjg(dcu
     1(2,j,k)) + dcu(1,j,k1)*conjg(dcu(1,j,k1)) + dcu(2,j,k1)*conjg(dcu(
     22,j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*aimag(ffe(1,k))
      at2 = at2*at2
      exy(1,1,k) = at1*dcu(1,1,k)
      exy(2,1,k) = at1*dcu(2,1,k)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      wp = wp + at2*(dcu(1,1,k)*conjg(dcu(1,1,k)) + dcu(2,1,k)*conjg(dcu
     1(2,1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*aimag(ffe(j,1))
      at2 = at2*at2
      exy(1,j,1) = at1*dcu(1,j,1)
      exy(2,j,1) = at1*dcu(2,j,1)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      wp = wp + at2*(dcu(1,j,1)*conjg(dcu(1,j,1)) + dcu(2,j,1)*conjg(dcu
     1(2,j,1)))
   70 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      wf = float(nx*ny)*wp/real(ffe(1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   80 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at2 = -ci2*real(ffe(j,k))
      at1 = at2*at2
      exy(1,j,k) = at2*dcu(1,j,k)
      exy(2,j,k) = at2*dcu(2,j,k)
      exy(1,j,k1) = at2*dcu(1,j,k1)
      exy(2,j,k1) = at2*dcu(2,j,k1)
      wp = wp + at2*(dcu(1,j,k)*conjg(dcu(1,j,k)) + dcu(2,j,k)*conjg(dcu
     1(2,j,k)) + dcu(1,j,k1)*conjg(dcu(1,j,k1)) + dcu(2,j,k1)*conjg(dcu(
     22,j,k1)))
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = -ci2*real(ffe(1,k))
      at1 = at2*at2
      exy(1,1,k) = at2*dcu(1,1,k)
      exy(2,1,k) = at2*dcu(2,1,k)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      wp = wp + at1*(dcu(1,1,k)*conjg(dcu(1,1,k)) + dcu(2,1,k)*conjg(dcu
     1(2,1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at2 = -ci2*real(ffe(j,1))
      at1 = at2*at2
      exy(1,j,1) = at2*dcu(1,j,1)
      exy(2,j,1) = at2*dcu(2,j,1)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      wp = wp + at1*(dcu(1,j,1)*conjg(dcu(1,j,1)) + dcu(2,j,1)*conjg(dcu
     1(2,j,1)))
  120 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      wf = float(nx*ny)*wp/real(ffe(1,1))
      return
      end
      subroutine DCGUARD2(daxy,nx,ny,nxe,nye,ndim)
c replicate extended periodic field
c quadratic interpolation
      implicit none
      real daxy
      integer nx, ny, nxe, nye, ndim
      dimension daxy(ndim,nxe,nye)
      integer i, j, k
      do 20 k = 1, ny
      do 10 i = 1, ndim
      daxy(i,1,k+1) = daxy(i,nx+1,k+1)
      daxy(i,nx+2,k+1) = daxy(i,2,k+1)
      daxy(i,nx+3,k+1) = daxy(i,3,k+1)
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      daxy(i,j+1,1) = daxy(i,j+1,ny+1)
      daxy(i,j+1,ny+2) = daxy(i,j+1,2)
      daxy(i,j+1,ny+3) = daxy(i,j+1,3)
   30 continue
   40 continue
      do 50 i = 1, ndim
      daxy(i,1,1) = daxy(i,nx+1,ny+1)
      daxy(i,nx+2,1) = daxy(i,2,ny+1)
      daxy(i,nx+3,1) = daxy(i,3,ny+1)
      daxy(i,1,ny+2) = daxy(i,nx+1,2)
      daxy(i,nx+2,ny+2) = daxy(i,2,2)
      daxy(i,nx+3,ny+2) = daxy(i,3,2)
      daxy(i,1,ny+3) = daxy(i,nx+1,3)
      daxy(i,nx+2,ny+3) = daxy(i,2,3)
      daxy(i,nx+3,ny+3) = daxy(i,3,3)
   50 continue
      return
      end
      subroutine DCGUARD2L(daxy,nx,ny,nxe,nye,ndim)
c replicate extended periodic field
c linear interpolation
      implicit none
      real daxy
      integer nx, ny, nxe, nye, ndim
      dimension daxy(ndim,nxe,nye)
c local data
      integer i, j, k
      do 20 k = 1, ny
      do 10 i = 1, ndim
      daxy(i,nx+1,k) = daxy(i,1,k)
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, ndim
      daxy(i,j,ny+1) = daxy(i,j,1)
   30 continue
   40 continue
      do 50 i = 1, ndim
      daxy(i,nx+1,ny+1) = daxy(i,1,1)
   50 continue
      return
      end
      subroutine APOIS23(cu,daxy,axy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,nx
     1vh,nyv,nxhd,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c smoothed derivative of vector potential, or smoothed vector potential,
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd,
c output: daxy,wm
c approximate flop count is: 90*nxc*nyc + 36*(nxc + nyc)
c for isign = 1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd, output: axy,wm
c approximate flop count is: 68*nxc*nyc + 21*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c daxy(1,kx,ky)=dax/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c daxy(2,kx,ky)=day/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cuy(kx,ky)*s(kx,ky)
c daxy(3,kx,ky)=dax/dy = ci*ci*sqrt(-1)*ky*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c daxy(4,kx,ky)=daz/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c daxy(5,kx,ky)=daz/dy = ci*ci*sqrt(-1)*ky*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for daxy(kx=pi) = daxy(ky=pi) = daxy(kx=0,ky=0) = 0.
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c daxy(kx=pi) = daxy(ky=pi) = daxy(kx=0,ky=0) = 0.
c if isign = 1, smoothed vector potential is calculated using:
c axy(1,kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c axy(2,kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)*s(kx,ky)
c axy(3,kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c except for axy(kx=pi) = axy(ky=pi) = axy(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c daxy(i,j,k) = i component of smoothed derivative of vector potential
c axy(i,j,k) = i component of smoothed vector potential
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, ci, wm
      complex cu, daxy, axy, ffc
      dimension cu(3,nxvh,nyv), daxy(5,nxvh,nyv), axy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2, zt3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate derivative of vector potential and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      daxy(1,j,k) = at2*zt3
      daxy(2,j,k) = at2*zt2
      daxy(3,j,k) = at3*zt3
      daxy(4,j,k) = at2*zt1
      daxy(5,j,k) = at3*zt1
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      daxy(1,j,k1) = at2*zt3
      daxy(2,j,k1) = at2*zt2
      daxy(3,j,k1) = -at3*zt3
      daxy(4,j,k1) = at2*zt1
      daxy(5,j,k1) = -at3*zt1
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
     2 cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      daxy(1,1,k) = zero
      daxy(2,1,k) = zero
      daxy(3,1,k) = at3*zt3
      daxy(4,1,k) = zero
      daxy(5,1,k) = at3*zt1
      daxy(1,1,k1) = zero
      daxy(2,1,k1) = zero
      daxy(3,1,k1) = zero
      daxy(4,1,k1) = zero
      daxy(5,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      zt3 = cmplx(-aimag(cu(1,j,1)),real(cu(1,j,1)))
      daxy(1,j,1) = at2*zt3
      daxy(2,j,1) = at2*zt2
      daxy(3,j,1) = zero
      daxy(4,j,1) = at2*zt1
      daxy(5,j,1) = zero
      daxy(1,j,k1) = zero
      daxy(2,j,k1) = zero
      daxy(3,j,k1) = zero
      daxy(4,j,k1) = zero
      daxy(5,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
   70 continue
      daxy(1,1,1) = zero
      daxy(2,1,1) = zero
      daxy(3,1,1) = zero
      daxy(4,1,1) = zero
      daxy(5,1,1) = zero
      daxy(1,1,k1) = zero
      daxy(2,1,k1) = zero
      daxy(3,1,k1) = zero
      daxy(4,1,k1) = zero
      daxy(5,1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothed vector potential and sum field energy
   80 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      axy(1,j,k) = at1*cu(1,j,k)
      axy(2,j,k) = at1*cu(2,j,k)
      axy(3,j,k) = at1*cu(3,j,k)
      axy(1,j,k1) = at1*cu(1,j,k1)
      axy(2,j,k1) = at1*cu(2,j,k1)
      axy(3,j,k1) = at1*cu(3,j,k1)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
     2 cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      axy(1,1,k) = at1*cu(1,1,k)
      axy(2,1,k) = at1*cu(2,1,k)
      axy(3,1,k) = at1*cu(3,1,k)
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      axy(1,j,1) = at1*cu(1,j,1)
      axy(2,j,1) = at1*cu(2,j,1)
      axy(3,j,1) = at1*cu(3,j,1)
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      axy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
  120 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      wm = float(nx*ny)*wp
      return
      end
      subroutine APOIS22(cu,daxy,axy,isign,ffc,ax,ay,affp,ci,wm,nx,ny,nx
     1vh,nyv,nxhd,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c smoothed derivative of vector potential, or smoothed vector potential,
c with periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nyhd, output: ffc
c for isign = -1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd,
c output: daxy,wm
c approximate flop count is: 63*nxc*nyc + 25*(nxc + nyc)
c for isign = 1, input: cu,ffc,isign,ci,nx,ny,nxvh,nyhd, output: axy,wm
c approximate flop count is: 48*nxc*nyc + 11*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c daxy(1,kx,ky)=dax/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c daxy(2,kx,ky)=day/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cuy(kx,ky)*s(kx,ky)
c daxy(3,kx,ky)=dax/dy = ci*ci*sqrt(-1)*ky*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c daxy(kx=pi) = daxy(ky=pi) = daxy(kx=0,ky=0) = 0.
c if isign = 1, smoothed vector potential is calculated using:
c axy(1,kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c axy(2,kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)*s(kx,ky)
c axy(3,kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c except for axy(kx=pi) = axy(ky=pi) = axy(kx=0,ky=0) = 0.
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c daxy(i,j,k) = i component of smoothed derivative of vector potential
c axy(i,j,k) = i component of smoothed vector potential
c all for fourier mode (j-1,k-1)
c aimag(ffc(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffc(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
c nxhd = first dimension of form factor array, must be >= nxh
c nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, ci, wm
      complex cu, daxy, axy, ffc
      dimension cu(2,nxvh,nyv), daxy(3,nxvh,nyv), axy(2,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, ci2, dkx, dky, at1, at2, at3, at4
      complex zero, zt2, zt3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffc(j,k) = cmplx(affp,1.)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate derivative of vector potential and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      daxy(1,j,k) = at2*zt3
      daxy(2,j,k) = at2*zt2
      daxy(3,j,k) = at3*zt3
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      daxy(1,j,k1) = at2*zt3
      daxy(2,j,k1) = at2*zt2
      daxy(3,j,k1) = -at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))
     2)
   40 continue
   50 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*float(k - 1)*at1
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      daxy(1,1,k) = zero
      daxy(2,1,k) = zero
      daxy(3,1,k) = at3*zt3
      daxy(1,1,k1) = zero
      daxy(2,1,k1) = zero
      daxy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*float(j - 1)*at1
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      zt3 = cmplx(-aimag(cu(1,j,1)),real(cu(1,j,1)))
      daxy(1,j,1) = at2*zt3
      daxy(2,j,1) = at2*zt2
      daxy(3,j,1) = zero
      daxy(1,j,k1) = zero
      daxy(2,j,k1) = zero
      daxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)))
   70 continue
      daxy(1,1,1) = zero
      daxy(2,1,1) = zero
      daxy(3,1,1) = zero
      daxy(1,1,k1) = zero
      daxy(2,1,k1) = zero
      daxy(3,1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothed vector potential and sum field energy
   80 wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nxh
      at1 = ci2*real(ffc(j,k))*aimag(ffc(j,k))
      axy(1,j,k) = at1*cu(1,j,k)
      axy(2,j,k) = at1*cu(2,j,k)
      axy(1,j,k1) = at1*cu(1,j,k1)
      axy(2,j,k1) = at1*cu(2,j,k1)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))
     2)
   90 continue
  100 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))*aimag(ffc(1,k))
      axy(1,1,k) = at1*cu(1,1,k)
      axy(2,1,k) = at1*cu(2,1,k)
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nxh
      at1 = ci2*real(ffc(j,1))*aimag(ffc(j,1))
      axy(1,j,1) = at1*cu(1,j,1)
      axy(2,j,1) = at1*cu(2,j,1)
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)))
  120 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      wm = float(nx*ny)*wp
      return
      end
      subroutine ADDQEI2(qe,qi,nx,ny,nxe,nye)
c adds electron and ion densities
c assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c nx/ny = system length in x/y direction
c nxe = first dimension of charge arrays, nxe must be >= nx
c nye = second dimension of charge arrays, nye must be >= ny
      implicit none
      real qe, qi
      integer nx, ny, nxe, nye
      dimension qe(nxe,nye), qi(nxe,nye)
c local data
      integer j, k
      do 20 k = 1, ny
      do 10 j = 1, nx
      qe(j,k) = qe(j,k) + qi(j,k)
   10 continue
   20 continue
      return
      end
      subroutine ADDQEI2X(qe,qi,qbme,qbmi,wpmax,wpmin,nx,ny,nxe,nye)
c adds electron and ion densities, and calculates maximum and minimum
c plasma frequency.  assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c qbme/qbmi = charge/mass ratio for electrons/ions
c wpmax/wpmin = maximum/minimum plasma frequency
c nx/ny = system length in x/y direction
c nxe = first dimension of charge arrays, nxe must be >= nx
c nye = second dimension of charge arrays, nye must be >= ny
      implicit none
      real qe, qi, qbme, qbmi, wpmax, wpmin
      integer nx, ny, nxe, nye
      dimension qe(nxe,nye), qi(nxe,nye)
c local data
      integer j, k
      real at1
      wpmax = qbme*qe(1,1) + qbmi*qi(1,1)
      wpmin = wpmax
      do 20 k = 1, ny
      do 10 j = 1, nx
      at1 = qbme*qe(j,k) + qbmi*qi(j,k)
      qe(j,k) = qe(j,k) + qi(j,k)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
      return
      end
      subroutine BADDEXT2(bxy,omx,omy,omz,nx,ny,nxe,nye)
c adds constant to magnetic field for 2-1/2d code
c bxy = magnetic field
c omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
c nx/ny = system length in x/y direction
c nxe = second dimension of magnetic field array, nxe must be >= nx
c nye = third dimension of magnetic field array, nye must be >= ny
      implicit none
      real bxy, omx, omy, omz
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer j, k
      do 20 k = 1, ny
      do 10 j = 1, nx
      bxy(1,j,k) = bxy(1,j,k) + omx
      bxy(2,j,k) = bxy(2,j,k) + omy
      bxy(3,j,k) = bxy(3,j,k) + omz
   10 continue
   20 continue
      return
      end
      subroutine BADDEXT22(bz,omz,nx,ny,nxe,nye)
c adds constant to magnetic field for 2d code
c bz = magnetic field
c omz = magnetic field electron cyclotron frequency in z 
c nx/ny = system length in x/y direction
c nxe = first dimension of magnetic field array, nxe must be >= nx
c nye = second dimension of magnetic field array, nye must be >= ny
      implicit none
      real bz, omz
      integer nx, ny, nxe, nye
      dimension bz(nxe,nye)
c local data
      integer j, k
      do 20 k = 1, ny
      do 10 j = 1, nx
      bz(j,k) = bz(j,k) + omz
   10 continue
   20 continue
      return
      end
      subroutine IMOMENT2(qi,fxy,pxi,pyi,pzi,dt,nx,ny,nxe,nye)
c calculate ion momentum for 2-1/2d code from integral of qi*fxy
c assumes guard cells have already been added
      implicit none
      real qi, fxy, pxi, pyi, pzi, dt
      integer nx, ny, nxe, nye
      dimension qi(nxe,nye), fxy(3,nxe,nye)
c local data
      integer j, k
      double precision dt1, sum1, sum2, sum3
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 20 k = 1, ny
      do 10 j = 1, nx
      dt1 = dble(qi(j,k))
      sum1 = sum1 + dt1*fxy(1,j,k)
      sum2 = sum2 + dt1*fxy(2,j,k)
      sum3 = sum3 + dt1*fxy(3,j,k)
   10 continue
   20 continue
      pxi = pxi + sum1*dt
      pyi = pyi + sum2*dt
      pzi = pzi + sum3*dt
      return
      end
      subroutine IMOMENT22(qi,fxy,pxi,pyi,dt,nx,ny,nxe,nye)
c calculate ion momentum for 2d code from integral of qi*fxy
c assumes guard cells have already been added
      implicit none
      real qi, fxy, pxi, pyi, dt
      integer nx, ny, nxe, nye
      dimension qi(nxe,nye), fxy(2,nxe,nye)
c local data
      integer j, k
      double precision dt1, sum1, sum2
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 20 k = 1, ny
      do 10 j = 1, nx
      dt1 = dble(qi(j,k))
      sum1 = sum1 + dt1*fxy(1,j,k)
      sum2 = sum2 + dt1*fxy(2,j,k)
   10 continue
   20 continue
      pxi = pxi + sum1*dt
      pyi = pyi + sum2*dt
      return
      end
      subroutine POIS23GL(q,fxy,affp,nx,ny,nxvh,nyv)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge with periodic boundary conditions.  Zeros out z component
c input: q,nx,ny,nxv, output: fxy
c approximate flop count is: 20*nxc*nyc + 16*(nxc + nyc)
c and nxc*nyc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2)), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
c fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c affp = normalization constant for field
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nxh
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex q, fxy, affp
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
c local data
      integer j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, dky2, at1, at2, at3
      complex zt1, zt2, zero
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      zero = cmplx(0.,0.)
c calculate force/charge
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j-1)
      at1 = affp/(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
   10 continue
c mode numbers kx = 0, nx/2
      dkx = dnx*real(nxh)
      at1 = affp/(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      at1 = affp/dky
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      zt2 = cmplx(aimag(q(1,k1)),-real(q(1,k1)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at1*zt1
      fxy(3,1,k) = zero
      fxy(1,1,k1) = at2*zt2
      fxy(2,1,k1) = -at3*zt2
      fxy(3,1,k1) = zero
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)
      dky2 = dky*dky
      do 30 j = 2, nxh
      dkx = dnx*real(j-1)
      at1 = affp/(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      at1 = affp/dkx
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,1) = at1*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = at3*zt2
      fxy(3,j,k1) = zero
   30 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      return
      end
      subroutine VCCOPY2(f,g,nx,ny,ndim,nxv,nyv)
c this subroutine copies (nx,ny) complex vector array elements
c to another complex array, nxv >= nx; nyv >= ny
      implicit none
      integer nx, ny, ndim, nxv, nyv
      complex f, g
      dimension f(ndim,nxv,nyv), g(ndim,nxv,nyv)
c local data
      integer i, j, k
      do 30 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, ndim
      g(i,j,k) = f(i,j,k)
   10 continue
   20 continue
   30 continue
      return
      end
      subroutine RTOC2FILL(fc,nx,ny,nxvh,nyv)
c this subroutine fills in place data from a real to complex packed
c format to a half-size complex array format
c nx/ny = system length in x/y direction
c nxv = first dimension of complex array fc, nxvh >= nx/2+1
c nyv = second dimension of complex array fc, nyv >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex fc
      dimension fc(nxvh,nyv)
c local data
      integer nxh, nyh, ny2, k, j1, k1
      complex zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers kx = nx/2
      do 10 k = 2, nyh
      k1 = ny2 - k
      fc(j1,k1) = fc(1,k1)
      fc(j1,k) = conjg(fc(1,k1))
      fc(1,k1) = conjg(fc(1,k))
   10 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      zt1 = cmplx(aimag(fc(1,1)),0.0)
      fc(1,1) = cmplx(real(fc(1,1)),0.0)
      fc(j1,1) = zt1
      zt1 = cmplx(aimag(fc(1,k1)),0.0)
      fc(1,k1) = cmplx(real(fc(1,k1)),0.0)
      fc(j1,k1) = zt1
      return
      end
      subroutine CTOR2FILL(fc,nx,ny,nxvh,nyv)
c this subroutine fills in place data from a half-size complex array
c format to a real to complex packed format
c nx/ny = system length in x/y direction
c nxvh = first dimension of complex array fc, nxvh >= nx/2+1
c nyv = second dimension of complex array fc, nyv >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex fc
      dimension fc(nxvh,nyv)
c local data
      integer nxh, nyh, ny2, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers kx = nx/2
      do 10 k = 2, nyh
      k1 = ny2 - k
      fc(1,k1) = fc(j1,k1)
   10 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      fc(1,1) = cmplx(real(fc(1,1)),real(fc(j1,1)))
      fc(1,k1) = cmplx(real(fc(1,k1)),real(fc(j1,k1)))
      return
      end
      subroutine RTOC2PACK(f,fc,nx,ny,nxe,nye,nxvh,nyv)
c this subroutine copies data stored in a real to complex packed format
c to a half-size complex array
c nx/ny = system length in x/y direction
c nxe = first dimension of input real array f, nxe >= nx
c nye = second dimension of input real array f, nye >= ny
c nxvh = first dimension of output complex array fc, nxvh >= nx/2+1
c nyv = second dimension of output complex array fc, nyv >= ny
      implicit none
      integer nx, ny, nxe, nye, nxvh, nyv
      real f
      complex fc
      dimension f(nxe,nye), fc(nxvh,nyv)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      fc(j,k) = cmplx(f(2*j-1,k),f(2*j,k))
      fc(j,k1) = cmplx(f(2*j-1,k1),f(2*j,k1))
   10 continue
c mode numbers kx = 0, nx/2
      fc(1,k) = cmplx(f(1,k),f(2,k))
      fc(1,k1) = cmplx(f(1,k),-f(2,k))
      fc(j1,k) = cmplx(f(1,k1),-f(2,k1))
      fc(j1,k1) = cmplx(f(1,k1),f(2,k1))
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      fc(j,1) = cmplx(f(2*j-1,1),f(2*j,1))
      fc(j,k1) = cmplx(f(2*j-1,k1),f(2*j,k1))
   30 continue
      fc(1,1) = cmplx(f(1,1),0.0)
      fc(j1,1) = cmplx(f(2,1),0.0)
      fc(1,k1) = cmplx(f(1,k1),0.0)
      fc(j1,k1) = cmplx(f(2,k1),0.0)
      return
      end
      subroutine CTOR2PACK(f,fc,nx,ny,nxe,nye,nxv,nyv)
c this subroutine copies data stored in a complex array to a
c real to complex packed format
c nx/ny = system length in x/y direction
c nxe = first dimension of output real array f, nxe >= nx
c nye = second dimension of output real array f, nye >= ny
c nxv = first dimension of input complex array fc, nxv >= nx/2+1
c nyv = second dimension of input complex array fc, nyv >= ny
      implicit none
      integer nx, ny, nxe, nye, nxv, nyv
      real f
      complex fc
      dimension f(nxe,nye), fc(nxv,nyv)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      f(2*j-1,k) = real(fc(j,k))
      f(2*j,k) = aimag(fc(j,k))
      f(2*j-1,k1) = real(fc(j,k1))
      f(2*j,k1) = aimag(fc(j,k1))
   10 continue
c mode numbers kx = 0, nx/2
      f(1,k) = real(fc(1,k))
      f(2,k) = aimag(fc(1,k))
      f(1,k1) = real(fc(j1,k1))
      f(2,k1) = aimag(fc(j1,k1))
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      f(2*j-1,1) = real(fc(j,1))
      f(2*j,1) = aimag(fc(j,1))
      f(2*j-1,k1) = real(fc(j,k1))
      f(2*j,k1) = aimag(fc(j,k1))
   30 continue
      f(1,1) = real(fc(1,1))
      f(2,1) = real(fc(j1,1))
      f(1,k1) = real(fc(1,k1))
      f(2,k1) = real(fc(j1,k1))
      return
      end
      subroutine RTOC2PACKT(f,gc,nx,ny,nxe,nye,nxvh,nyv)
c this subroutine copies data stored in a real to complex packed format
c to a transposed complex array
c nx/ny = system length in x/y direction
c nxe = first dimension of input real array f, nxe >= nx
c nye = second dimension of input real array f, nye >= ny
c nxvh = second dimension of output complex array gc, nxvh >= nx/2+1
c nyv = first dimension of output complex array gc, nyv >= ny
      implicit none
      integer nx, ny, nxe, nye, nxvh, nyv
      real f
      complex gc
      dimension f(nxe,nye), gc(nyv,nxvh)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      gc(k,j) = cmplx(f(2*j-1,k),f(2*j,k))
      gc(k1,j) = cmplx(f(2*j-1,k1),f(2*j,k1))
   10 continue
c mode numbers kx = 0, nx/2
      gc(k,1) = cmplx(f(1,k),f(2,k))
      gc(k1,1) = cmplx(f(1,k),-f(2,k))
      gc(k,j1) = cmplx(f(1,k1),-f(2,k1))
      gc(k1,j1) = cmplx(f(1,k1),f(2,k1))
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      gc(1,j) = cmplx(f(2*j-1,1),f(2*j,1))
      gc(k1,j) = cmplx(f(2*j-1,k1),f(2*j,k1))
   30 continue
      gc(1,1) = cmplx(f(1,1),0.0)
      gc(1,j1) = cmplx(f(2,1),0.0)
      gc(k1,1) = cmplx(f(1,k1),0.0)
      gc(k1,j1) = cmplx(f(2,k1),0.0)
      return
      end
      subroutine CTOR2PACKT(f,gc,nx,ny,nxe,nye,nxvh,nyv)
c this subroutine copies data stored in a transposed complex array to a
c real to complex packed format
c nx/ny = system length in x/y direction
c nxe = first dimension of output real array f, nxe >= nx
c nye = second dimension of output real array f, nye >= ny
c nxvh = second dimension of input complex array gc, nxvh >= nx/2+1
c nyv = first dimension of input complex array gc, nyv >= ny
      implicit none
      integer nx, ny, nxe, nye, nxvh, nyv
      real f
      complex gc
      dimension f(nxe,nye), gc(nyv,nxvh)
c local data
      integer nxh, nyh, ny2, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 j = 2, nxh
      do 10 k = 2, nyh
      k1 = ny2 - k
      f(2*j-1,k) = real(gc(k,j))
      f(2*j,k) = aimag(gc(k,j))
      f(2*j-1,k1) = real(gc(k1,j))
      f(2*j,k1) = aimag(gc(k1,j))
   10 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      f(2*j-1,1) = real(gc(1,j))
      f(2*j,1) = aimag(gc(1,j))
      f(2*j-1,k1) = real(gc(k1,j))
      f(2*j,k1) = aimag(gc(k1,j))
   20 continue
c mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      f(1,k) = real(gc(k,1))
      f(2,k) = aimag(gc(k,1))
      f(1,k1) = real(gc(k1,j1))
      f(2,k1) = aimag(gc(k1,j1))
   30 continue
      k1 = nyh + 1
      f(1,1) = real(gc(1,1))
      f(2,1) = real(gc(1,j1))
      f(1,k1) = real(gc(k1,1))
      f(2,k1) = real(gc(k1,j1))
      return
      end
      subroutine RTOCN2PACKT(fc,gc,nx,ny,ndim,nxeh,nye,nxvh,nyv)
c this subroutine copies mutiple data stored in a real to complex packed
c format to a multiple transposed complex array
c nx/ny = system length in x/y direction
c ndim = leading dimension of array fc
c nxeh = second dimension of output real array fc, nxeh >= nx/2
c nye = second dimension of input real array fc, nye >= ny
c nxvh = second dimension of output complex array gc, nxvh >= nx/2+1
c nyv = first dimension of output complex array gc, nyv >= ny
      implicit none
      integer nx, ny, ndim, nxeh, nye, nxvh, nyv
      complex fc, gc
      dimension fc(ndim,nxeh,nye), gc(nyv,ndim,nxvh)
c local data
      integer nxh, nyh, ny2, i, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 40 k = 2, nyh
      k1 = ny2 - k
      do 20 j = 2, nxh
      do 10 i = 1, ndim
      gc(k,i,j) = fc(i,j,k)
      gc(k1,i,j) = fc(i,j,k1)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
      do 30 i = 1, ndim
      gc(k,i,1) = fc(i,1,k)
      gc(k1,i,1) = conjg(fc(i,1,k))
      gc(k,i,j1) = conjg(fc(i,1,k1))
      gc(k1,i,j1) = fc(i,1,k1)
   30 continue
   40 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nxh
      do 50 i = 1, ndim
      gc(1,i,j) = fc(i,j,1)
      gc(k1,i,j) = fc(i,j,k1)
   50 continue
   60 continue
      do 70 i = 1, ndim
      gc(1,i,1) = cmplx(real(fc(i,1,1)),0.0)
      gc(1,i,j1) = cmplx(aimag(fc(i,1,1)),0.0)
      gc(k1,i,1) = cmplx(real(fc(i,1,k1)),0.0)
      gc(k1,i,j1) = cmplx(aimag(fc(i,1,k1)),0.0)
   70 continue
      return
      end
      subroutine CTORN2PACKT(fc,gc,nx,ny,ndim,nxeh,nye,nxvh,nyv)
c this subroutine copies multiple data stored in a transposed complex
c array to multiple real to complex packed formats
c nx/ny = system length in x/y direction
c ndim = leading dimension of array fc
c nxeh = second dimension of output real array fc, nxeh >= nx/2
c nye = third dimension of output real array fc, nye >= ny
c nxvh = third dimension of input complex array gc, nxvh >= nx/2+1
c nyv = second dimension of input complex array gc, nyv >= ny
      implicit none
      integer nx, ny, ndim, nxeh, nye, nxvh, nyv
      complex fc, gc
      dimension fc(ndim,nxeh,nye), gc(nyv,ndim,nxvh)
c local data
      integer nxh, nyh, ny2, i, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      j1 = nxh + 1
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 j = 2, nxh
      do 20 i = 1, ndim
      do 10 k = 2, nyh
      k1 = ny2 - k
      fc(i,j,k) = gc(k,i,j)
      fc(i,j,k1) = gc(k1,i,j)
   10 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      fc(i,j,1) = gc(1,i,j)
      fc(i,j,k1) = gc(k1,i,j)
   20 continue
   30 continue
c mode numbers kx = 0, nx/2
      do 50 i = 1, ndim
      do 40 k = 2, nyh
      k1 = ny2 - k
      fc(i,1,k) = gc(k,i,1)
      fc(i,1,k1) = gc(k1,i,j1)
   40 continue
      k1 = nyh + 1
      fc(i,1,1) = cmplx(real(gc(1,i,1)),real(gc(1,i,j1)))
      fc(i,1,k1) = cmplx(real(gc(k1,i,1)),real(gc(k1,i,j1)))
   50 continue
      return
      end
      subroutine RTOC2PACKD(f,fc,nx,ny,nxe,nye,nxv,nyv)
c this subroutine copies data stored in a real to complex packed format
c to a complex array
c nx/ny = system length in x/y direction
c nxe = first dimension of input real array f, nxe >= nx
c nye = second dimension of input real array f, nye >= ny
c nxv = first dimension of output complex array fc, nxv >= nx
c nyv = second dimension of output complex array fc, nyv >= ny
      implicit none
      integer nx, ny, nxe, nye, nxv, nyv
      real f
      complex fc
      dimension f(nxe,nye), fc(nxv,nyv)
c local data
      integer nxh, nyh, nx2, ny2, j, k, j1, k1
      nxh = nx/2
      nyh = max(1,ny/2)
      nx2 = nx + 2
      ny2 = ny + 2
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nxh
      j1 = nx2 - j
      fc(j,k) = cmplx(f(2*j-1,k),f(2*j,k))
      fc(j,k1) = cmplx(f(2*j-1,k1),f(2*j,k1))
      fc(j1,k) = conjg(fc(j,k1))
      fc(j1,k1) = conjg(fc(j,k))
   10 continue
c mode numbers kx = 0, nx/2
      j1 = nxh + 1
      fc(1,k) = cmplx(f(1,k),f(2,k))
      fc(1,k1) = cmplx(f(1,k),-f(2,k))
      fc(j1,k) = cmplx(f(1,k1),-f(2,k1))
      fc(j1,k1) = cmplx(f(1,k1),f(2,k1))
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nxh
      j1 = nx2 - j
      fc(j,1) = cmplx(f(2*j-1,1),f(2*j,1))
      fc(j1,1) = conjg(fc(j,1))
      fc(j,k1) = cmplx(f(2*j-1,k1),f(2*j,k1))
      fc(j1,k1) = conjg(fc(j,k1))
   30 continue
      j1 = nxh + 1
      fc(1,1) = cmplx(f(1,1),0.0)
      fc(j1,1) = cmplx(f(2,1),0.0)
      fc(1,k1) = cmplx(f(1,k1),0.0)
      fc(j1,k1) = cmplx(f(2,k1),0.0)
      return
      end
      subroutine ADDVRFIELD2(a,b,c,ndim,nxe,nye)
c this subroutine calculates a = b + c for real vector fields
      implicit none
      integer ndim, nxe, nye
      real a, b, c
      dimension a(ndim,nxe,nye), b(ndim,nxe,nye), c(ndim,nxe,nye)
c local data
      integer i, j, k
      do 30 k = 1, nye
      do 20 j = 1, nxe
      do 10 i = 1, ndim
      a(i,j,k) = b(i,j,k) + c(i,j,k)
   10 continue
   20 continue
   30 continue
      return
      end


