c-----------------------------------------------------------------------
c 2d PIC library for solving field equations with dirichlet boundary
c conditions
c dfield2lib.f contains procedures to manage guard cells and solve
c              fields equations in fourier space for dirichlet
c              boundary conditions:
c LCGUARD2 replicates fields for 2 component vector array to replace
c          quadratic with linear interpolation at the edges.
c LDGUARD2 replicates fields for scalar array to replace quadratic with
c          linear interpolation at the edges.
c LBGUARD2 replicates fields for 3 component vector array to replace
c          quadratic with linear interpolation at the edges.
c LSCGUARD2 initialize field for 3 component non-periodic vector array,
c           quadratic interpolation.
c LSCGUARD22 initialize field for 2 component non-periodic vector array,
c            quadratic interpolation.
c LSGUARD2 initialize field for non-periodic scalar array, quadratic
c          interpolation.
c LACGUARD2 add guard cells for 3 component non-periodic vector array,
c           to replace quadratic with linear interpolation at the edges.
c LACGUARD22 add guard cells for 2 component non-periodic vector array,
c            to replace quadratic with linear interpolation at the edges.
c LAGUARD2 add guard cells for non-periodic scalar array, to replace
c          quadratic with linear interpolation at the edges.
c LSCGUARD2L initialize field for 3 component non-periodic vector array,
c            quadratic interpolation.
c LSCGUARD22L initialize field for 2 component non-periodic vector
c             array, quadratic interpolation.
c LSGUARD2L initialize field for non-periodic scalar array, quadratic
c          interpolation.
c DBLSIN2C creates doubled array for 2 component vector data to enable
c          various sine/cosine transforms to be perfomed with ffts.
c DBLSIN2D creates doubled array for scalar vector data to enable
c          various sine/cosine transforms to be perfomed with ffts.
c DBLSIN2B creates doubled array for 3 component vector data to enable
c          various sine/cosine transforms to be perfomed with ffts.
c HAFDBL2C extracts data from doubled array for 2 component vector data.
c HAFDBL2D extracts data from doubled array for scalar data.
c HAFDBL2B extracts data from doubled array for 3 component vector data.
c POISDX2 solve 2d poisson equation for electric force, potential, or
c         smoothing with dirichlet boundary conditions, using doubled
c         ffts.
c POISD2 solve 2d poisson equation for electric force, potential, or
c        smoothing with dirichlet boundary conditions, using sine/cosine
c        transforms.
c POISDX22 solve 2d poisson equation for electric force with dirichlet
c          boundary conditions, using doubled ffts.
c POISD22 solve 2d poisson equation for electric force with dirichlet
c         boundary conditions, using sine/cosine transforms.
c POISDX23 solve 2-1/2d poisson equation for electric force with
c          dirichlet boundary conditions, using doubled ffts.
c POISD23 solve 2-1/2d poisson equation for electric force with
c         dirichlet boundary conditions, using sine/cosine transforms.
c POISDX2T solve 2d poisson equation for electric force, potential, or
c          smoothing with dirichlet boundary conditions for transposed
c          data, using doubled ffts.
c POISD2T solve 2d poisson equation for electric force, potential, or
c         smoothing with dirichlet boundary conditions for transposed
c         data, using sine/cosine transforms.
c DIVFD2 calculates 2d divergence of n component vector in fourier space
c        with dirichlet boundary conditions using sine/cosine 
c        transforms.
c GRADFD2 calculates 2d gradient of scalar field in fourier space,
c         with dirichlet boundary conditions using sine/cosine
c         transforms.
c CURLFD2 calculates 2d divergence of 3 component vector in fourier
c         space with dirichlet boundary conditions using sine/cosine
c         transforms.
c CURLFD22 calculates 2d divergence of 2 component vector in fourier
c          space with dirichlet boundary conditions using sine/cosine
c          transforms.
c CUPERPDX2 calculates 2d transverse current of 3 component vector in
c           fourier space with dirichlet boundary conditions, using
c           doubled ffts.
c CUPERPD2 calculates 2d transverse current of 3 component vector in
c          fourier space with dirichlet boundary conditions, using
c          sine/cosine transforms.
c CUPERPDX22 calculates 2d transverse current of 2 component vector in
c           fourier space with dirichlet boundary conditions, using
c           doubled ffts.
c CUPERPD22 calculates 2d transverse current of 2 component vector in
c          fourier space with dirichlet boundary conditions, using
c          sine/cosine transforms.
c BPOISDX23 solve 2-1/2d vector poisson equation for magnetic force,
c           vector potential, or smoothing with dirichlet boundary
c           conditions, using doubled ffts.
c BPOISD23 solve 2-1/2d vector poisson equation for magnetic force,
c          vector potential, or smoothing with dirichlet boundary
c          conditions, using sine/cosine transforms.
c BPOISDX22 solve 2d vector poisson equation for magnetic force,
c           vector potential, or smoothing with dirichlet boundary
c           conditions, using doubled ffts.
c BPOISD22 solve 2d vector poisson equation for magnetic force,
c          vector potential, or smoothing with dirichlet boundary
c          conditions, using sine/cosine transforms.
c IBPOISDX23 solve 2-1/2d vector poisson equation for magnetic field
c            with dirichlet boundary conditions, using doubled ffts.
c IBPOISD23 solve 2-1/2d vector poisson equation for magnetic field
c           with dirichlet boundary conditions, using sine/cosine
c           transforms.
c MAXWELDX2 solve 2d maxwell equation for electric and magnetic fields
c           with dirichlet boundary conditions, using doubled ffts.
c MAXWELD2 solve 2d maxwell equation for electric and magnetic fields
c          with dirichlet boundary conditions, using sine/cosine
c          transforms.
c DMFIELDD2 copies scalar data from doubled fft format to sine format.
c FMFIELDD2 copies scalar data from doubled fft format to cosine format.
c CMFIELDD2 copies 3 component vector data from doubled fft format to
c           sine/cosine format.
c EMFIELDD2 combines and smooths 2d electric or magnetic forces in
c           sine/cosine or cosine/sine format to doubled fft format.
c PMFIELDD2 copies scalar data from sine format to doubled fft format.
c BMFIELDD2 copies 3 component vector data from cosine/sine format to
c           doubled fft format.
c CPFIELDD2 copies 2d electric field in sine/cosine to doubled fft
c           format.
c AVPOTDX23 calculate 2-1/2d vector potential from magnetic field
c           with dirichlet boundary conditions, using doubled ffts.
c AVPOTD23 calculate 2-1/2d vector potential from magnetic field
c          with dirichlet boundary conditions, using sine/cosine
c          transforms.
c DPOYNTDX2 calculate electromagnetic flux in 2-1/2d Darwin field
c           with dirichlet boundary conditions, using doubled ffts.
c LACFGUARD2 increment 3 component non-periodic field with scaled vector
c            array, quadratic interpolation.
c LACFGUARD22 increment 2 component non-periodic field with scaled
c             vector array, quadratic interpolation.
c LSCFGUARD2L initialize 3 component non-periodic field with scaled
c             vector array,linear interpolation.
c LSCFGUARD22L initialize 2 component non-periodic field with scaled
c              vector array, linear interpolation.
c LSMCGUARD2 initialize non-periodic field for 4 component tensor array,
c            quadratic interpolation.
c LSMCGUARD22 initialize non-periodic field for 2 component tensor
c             array, quadratic interpolation.
c LAMCGUARD2 add guard cells for 4 component non-periodic tensor array,
c            quadratic interpolation.
c LSMCGUARD2L initialize non-periodic field for 4 component tensor
c             array, linear interpolation.
c LSMCGUARD22L initialize non-periodic field for 2 component tensor
c              array, linear interpolation.
c DBLSIN2M creates doubled array for 4 component tensor data to enable
c          various sine/cosine transforms to be perfomed with ffts.
c DBLSIN22M creates doubled array for 2 component tensor data to enable
c           various sine/cosine transforms to be perfomed with ffts.
c DCUPERPDX23 calculate 2-1/d transverse derivative of current density
c             from momentum flux with dirichlet boundary conditions,
c             using doubled ffts.
c ADCUPERPDX23 calculate 2-1/2d transverse derivative of current density
c             from momentum flux and acceleration density with dirichlet
c             boundary conditions, using doubled ffts.
c EPOISDX23 solve vector poisson equation for 2-1/2d transverse electric
c           field or force with dirichlet boundary conditions, using
c           doubled ffts.
c HAFDBL2N extracts data from doubled array for n component tensor data.
c LDCGUARD2 replicates fields for n component vector array to replace
c           quadratic with linear interpolation at the edges.
c APOISDX23 solves 2-1/2d poisson equation for smoothed derivative of
c           vector potential or smoothed vector potential with dirichlet
c           boundary conditions, using doubled ffts.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: july 21, 2010
c-----------------------------------------------------------------------
      subroutine LCGUARD2(fxy,nx,ny,nxe,nye)
c this subroutine replicates field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = first dimension of input array fxy, must be >= nx+3
c nye = second dimension of input array fxy, must be >= ny+3
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer j, k, ny1, nx3
      ny1 = ny + 1
      nx3 = nx + 3
      do 10 k = 1, ny1
      fxy(1,1,k+1) = 2.*fxy(1,2,k+1) - fxy(1,3,k+1)
      fxy(2,1,k+1) = 2.*fxy(2,2,k+1) - fxy(2,3,k+1)
      fxy(1,nx3,k+1) = 2.*fxy(1,nx+2,k+1) - fxy(1,nx+1,k+1)
      fxy(2,nx3,k+1) = 2.*fxy(2,nx+2,k+1) - fxy(2,nx+1,k+1)
   10 continue
      do 20 j = 1, nx3
      fxy(1,j,1) = 2.*fxy(1,j,2) - fxy(1,j,3)
      fxy(2,j,1) = 2.*fxy(2,j,2) - fxy(2,j,3)
      fxy(1,j,ny+3) = 2.*fxy(1,j,ny+2) - fxy(1,j,ny+1)
      fxy(2,j,ny+3) = 2.*fxy(2,j,ny+2) - fxy(2,j,ny+1)
   20 continue
      return
      end
      subroutine LDGUARD2(q,nx,ny,nxe,nye)
c this subroutine replicates scalar field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = first dimension of input array q, must be >= nx+3
c nye = second dimension of input array q, must be >= ny+3
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k, ny1, nx3
      ny1 = ny + 1
      nx3 = nx + 3
      do 10 k = 1, ny1
      q(1,k+1) = 2.*q(2,k+1) - q(3,k+1)
      q(nx3,k+1) = 2.*q(nx+2,k+1) - q(nx+1,k+1)
   10 continue
      do 20 j = 1, nx3
      q(j,1) = 2.*q(j,2) - q(j,3)
      q(j,ny+3) = 2.*q(j,ny+2) - q(j,ny+1)
   20 continue
      return
      end
      subroutine LBGUARD2(bxy,nx,ny,nxe,nye)
c this subroutine replicates vector field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = first dimension of input array bxy, must be >= nx+3
c nye = second dimension of input array bxy, must be >= ny+3
      implicit none
      real bxy
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer j, k, ny1, nx3
      ny1 = ny + 1
      nx3 = nx + 3
      do 10 k = 1, ny1
      bxy(1,1,k+1) = 2.*bxy(1,2,k+1) - bxy(1,3,k+1)
      bxy(2,1,k+1) = 2.*bxy(2,2,k+1) - bxy(2,3,k+1)
      bxy(3,1,k+1) = 2.*bxy(3,2,k+1) - bxy(3,3,k+1)
      bxy(1,nx3,k+1) = 2.*bxy(1,nx+2,k+1) - bxy(1,nx+1,k+1)
      bxy(2,nx3,k+1) = 2.*bxy(2,nx+2,k+1) - bxy(2,nx+1,k+1)
      bxy(3,nx3,k+1) = 2.*bxy(3,nx+2,k+1) - bxy(3,nx+1,k+1)
   10 continue
      do 20 j = 1, nx3
      bxy(1,j,1) = 2.*bxy(1,j,2) - bxy(1,j,3)
      bxy(2,j,1) = 2.*bxy(2,j,2) - bxy(2,j,3)
      bxy(3,j,1) = 2.*bxy(3,j,2) - bxy(3,j,3)
      bxy(1,j,ny+3) = 2.*bxy(1,j,ny+2) - bxy(1,j,ny+1)
      bxy(2,j,ny+3) = 2.*bxy(2,j,ny+2) - bxy(2,j,ny+1)
      bxy(3,j,ny+3) = 2.*bxy(3,j,ny+2) - bxy(3,j,ny+1)
   20 continue
      return
      end
      subroutine LSCGUARD2(cu,xj0,yj0,zj0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension cu(3,nxe,nye)
      integer i, j, k, nxg, nyg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx3 = nx + 3
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      cu(1,j+ngx+1,k+ngy+1) = xj0
      cu(2,j+ngx+1,k+ngy+1) = yj0
      cu(3,j+ngx+1,k+ngy+1) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+ngy+1) = 0.
      cu(i,2,k+ngy+1) = 0.
      cu(i,nx+2,k+ngy+1) = 0.
      cu(i,nx+3,k+ngy+1) = 0.
   20 continue
      cu(1,ngx+2,k+ngy+1) = .5*xj0
      cu(2,ngx+2,k+ngy+1) = .5*yj0
      cu(3,ngx+2,k+ngy+1) = .5*zj0
      cu(1,nx-ngx+2,k+ngy+1) = .5*xj0
      cu(2,nx-ngx+2,k+ngy+1) = .5*yj0
      cu(3,nx-ngx+2,k+ngy+1) = .5*zj0
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 3
      cu(i,j,1) = 0.
      cu(i,j,2) = 0.
      cu(i,j,ny+2) = 0.
      cu(i,j,ny+3) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      cu(1,j+ngx+1,ngy+2) = .5*xj0
      cu(2,j+ngx+1,ngy+2) = .5*yj0
      cu(3,j+ngx+1,ngy+2) = .5*zj0
      cu(1,j+ngx+1,ny-ngy+2) = .5*xj0
      cu(2,j+ngx+1,ny-ngy+2) = .5*yj0
      cu(3,j+ngx+1,ny-ngy+2) = .5*zj0
   60 continue
      do 70 i = 1, 3
      cu(i,1,ngy+2) = 0.
      cu(i,2,ngy+2) = 0.
      cu(i,nx+2,ngy+2) = 0.
      cu(i,nx+3,ngy+2) = 0.
      cu(i,1,ny-ngy+2) = 0.
      cu(i,2,ny-ngy+2) = 0.
      cu(i,nx+2,ny-ngy+2) = 0.
      cu(i,nx+3,ny-ngy+2) = 0.
   70 continue
      cu(1,ngx+2,ngy+2) = .25*xj0
      cu(2,ngx+2,ngy+2) = .25*yj0
      cu(3,ngx+2,ngy+2) = .25*zj0
      cu(1,nx-ngx+2,ngy+2) = .25*xj0
      cu(2,nx-ngx+2,ngy+2) = .25*yj0
      cu(3,nx-ngx+2,ngy+2) = .25*zj0
      cu(1,ngx+2,ny-ngy+2) = .25*xj0
      cu(2,ngx+2,ny-ngy+2) = .25*yj0
      cu(3,ngx+2,ny-ngy+2) = .25*zj0
      cu(1,nx-ngx+2,ny-ngy+2) = .25*xj0
      cu(2,nx-ngx+2,ny-ngy+2) = .25*yj0
      cu(3,nx-ngx+2,ny-ngy+2) = .25*zj0
      return
      end
      subroutine LSCGUARD22(cu,xj0,yj0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension cu(2,nxe,nye)
      integer i, j, k, nxg, nyg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx3 = nx + 3
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      cu(1,j+ngx+1,k+ngy+1) = xj0
      cu(2,j+ngx+1,k+ngy+1) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+ngy+1) = 0.
      cu(i,2,k+ngy+1) = 0.
      cu(i,nx+2,k+ngy+1) = 0.
      cu(i,nx+3,k+ngy+1) = 0.
   20 continue
      cu(1,ngx+2,k+ngy+1) = .5*xj0
      cu(2,ngx+2,k+ngy+1) = .5*yj0
      cu(1,nx-ngx+2,k+ngy+1) = .5*xj0
      cu(2,nx-ngx+2,k+ngy+1) = .5*yj0
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 2
      cu(i,j,1) = 0.
      cu(i,j,2) = 0.
      cu(i,j,ny+2) = 0.
      cu(i,j,ny+3) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      cu(1,j+ngx+1,ngy+2) = .5*xj0
      cu(2,j+ngx+1,ngy+2) = .5*yj0
      cu(1,j+ngx+1,ny-ngy+2) = .5*xj0
      cu(2,j+ngx+1,ny-ngy+2) = .5*yj0
   60 continue
      do 70 i = 1, 2
      cu(i,1,ngy+2) = 0.
      cu(i,2,ngy+2) = 0.
      cu(i,nx+2,ngy+2) = 0.
      cu(i,nx+3,ngy+2) = 0.
      cu(i,1,ny-ngy+2) = 0.
      cu(i,2,ny-ngy+2) = 0.
      cu(i,nx+2,ny-ngy+2) = 0.
      cu(i,nx+3,ny-ngy+2) = 0.
   70 continue
      cu(1,ngx+2,ngy+2) = .25*xj0
      cu(2,ngx+2,ngy+2) = .25*yj0
      cu(1,nx-ngx+2,ngy+2) = .25*xj0
      cu(2,nx-ngx+2,ngy+2) = .25*yj0
      cu(1,ngx+2,ny-ngy+2) = .25*xj0
      cu(2,ngx+2,ny-ngy+2) = .25*yj0
      cu(1,nx-ngx+2,ny-ngy+2) = .25*xj0
      cu(2,nx-ngx+2,ny-ngy+2) = .25*yj0
      return
      end
      subroutine LSGUARD2(q,qi0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic scalar field
c ngx/ngy = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real q, qi0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension q(nxe,nye)
      integer j, k, nxg, nyg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx3 = nx + 3
      do 20 k = 2, nyg
      do 10 j = 2, nxg
      q(j+ngx+1,k+ngy+1) = qi0
   10 continue
      q(1,k+ngy+1) = 0.
      q(2,k+ngy+1) = 0.
      q(nx+2,k+ngy+1) = 0.
      q(nx+3,k+ngy+1) = 0.
      q(ngx+2,k+ngy+1) = .5*qi0
      q(nx-ngx+2,k+ngy+1) = .5*qi0
   20 continue
      do 30 j = 1, nx3
      q(j,1) = 0.
      q(j,2) = 0.
      q(j,ny+2) = 0.
      q(j,ny+3) = 0.
   30 continue
      do 40 j = 2, nxg
      q(j+ngx+1,ngy+2) = .5*qi0
      q(j+ngx+1,ny-ngy+2) = .5*qi0
   40 continue
      q(1,ngy+2) = 0.
      q(2,ngy+2) = 0.
      q(nx+2,ngy+2) = 0.
      q(nx+3,ngy+2) = 0.
      q(1,ny-ngy+2) = 0.
      q(2,ny-ngy+2) = 0.
      q(nx+2,ny-ngy+2) = 0.
      q(nx+3,ny-ngy+2) = 0.
      q(ngx+2,ngy+2) = .25*qi0
      q(nx-ngx+2,ngy+2) = .25*qi0
      q(ngx+2,ny-ngy+2) = .25*qi0
      q(nx-ngx+2,ny-ngy+2) = .25*qi0
      return
      end
      subroutine LACGUARD2(cu,nx,ny,nxe,nye)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = second dimension of input array cu, must be >= nx+3
c nye = third dimension of input array cu, must be >= ny+3
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
c local data
      integer i, j, k, nx1, ny3
      nx1 = nx + 1
      ny3 = ny + 3
c add up guard cells
      do 20 k = 1, ny3
      do 10 i = 1, 3
      cu(i,2,k) = cu(i,2,k) + 2.*cu(i,1,k)
      cu(i,3,k) = cu(i,3,k) - cu(i,1,k)
      cu(i,nx+1,k) = cu(i,nx+1,k) - cu(i,nx+3,k)
      cu(i,nx+2,k) = cu(i,nx+2,k) + 2.*cu(i,nx+3,k)
      cu(i,1,k) = 0.
      cu(i,nx+3,k) = 0.
   10 continue
   20 continue
      do 40 j = 1, nx1
      do 30 i = 1, 3
      cu(i,j+1,2) = cu(i,j+1,2) + 2.*cu(i,j+1,1)
      cu(i,j+1,3) = cu(i,j+1,3) - cu(i,j+1,1)
      cu(i,j+1,ny+1) = cu(i,j+1,ny+1) - cu(i,j+1,ny+3)
      cu(i,j+1,ny+2) = cu(i,j+1,ny+2) + 2.*cu(i,j+1,ny+3)
      cu(i,j+1,1) = 0.
      cu(i,j+1,ny+3) = 0.
   30 continue
   40 continue
      return
      end
      subroutine LACGUARD22(cu,nx,ny,nxe,nye)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = second dimension of input array cu, must be >= nx+3
c nye = third dimension of input array cu, must be >= ny+3
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
c local data
      integer i, j, k, nx1, ny3
      nx1 = nx + 1
      ny3 = ny + 3
c add up guard cells
      do 20 k = 1, ny3
      do 10 i = 1, 2
      cu(i,2,k) = cu(i,2,k) + 2.*cu(i,1,k)
      cu(i,3,k) = cu(i,3,k) - cu(i,1,k)
      cu(i,nx+1,k) = cu(i,nx+1,k) - cu(i,nx+3,k)
      cu(i,nx+2,k) = cu(i,nx+2,k) + 2.*cu(i,nx+3,k)
      cu(i,1,k) = 0.
      cu(i,nx+3,k) = 0.
   10 continue
   20 continue
      do 40 j = 1, nx1
      do 30 i = 1, 2
      cu(i,j+1,2) = cu(i,j+1,2) + 2.*cu(i,j+1,1)
      cu(i,j+1,3) = cu(i,j+1,3) - cu(i,j+1,1)
      cu(i,j+1,ny+1) = cu(i,j+1,ny+1) - cu(i,j+1,ny+3)
      cu(i,j+1,ny+2) = cu(i,j+1,ny+2) + 2.*cu(i,j+1,ny+3)
      cu(i,j+1,1) = 0.
      cu(i,j+1,ny+3) = 0.
   30 continue
   40 continue
      return
      end
      subroutine LAGUARD2(q,nx,ny,nxe,nye)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c for scalar field
c nx/ny = system length in x/y direction
c nxe = first dimension of input array q, must be >= nx+3
c nye = second dimension of input array q, must be >= ny+3
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k, nx1, ny3
      nx1 = nx + 1
      ny3 = ny + 3
c add up guard cells
      do 10 k = 1, ny3
      q(2,k) = q(2,k) + 2.*q(1,k)
      q(3,k) = q(3,k) - q(1,k)
      q(nx+1,k) = q(nx+1,k) - q(nx+3,k)
      q(nx+2,k) = q(nx+2,k) + 2.*q(nx+3,k)
      q(1,k) = 0.
      q(nx+3,k) = 0.
   10 continue
      do 20 j = 1, nx1
      q(j+1,2) = q(j+1,2) + 2.*q(j+1,1)
      q(j+1,3) = q(j+1,3) - q(j+1,1)
      q(j+1,ny+1) = q(j+1,ny+1) - q(j+1,ny+3)
      q(j+1,ny+2) = q(j+1,ny+2) + 2.*q(j+1,ny+3)
      q(j+1,1) = 0.
      q(j+1,ny+3) = 0.
   20 continue
      return
      end
      subroutine LSCGUARD2L(cu,xj0,yj0,zj0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension cu(3,nxe,nye)
      integer i, j, k, nxg, nyg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx1 = nx + 1
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      cu(1,j+ngx,k+ngy) = xj0
      cu(2,j+ngx,k+ngy) = yj0
      cu(3,j+ngx,k+ngy) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+ngy) = 0.
      cu(i,nx+1,k+ngy) = 0.
   20 continue
      cu(1,ngx+1,k+ngy) = .5*xj0
      cu(2,ngx+1,k+ngy) = .5*yj0
      cu(3,ngx+1,k+ngy) = .5*zj0
      cu(1,nx-ngx+1,k+ngy) = .5*xj0
      cu(2,nx-ngx+1,k+ngy) = .5*yj0
      cu(3,nx-ngx+1,k+ngy) = .5*zj0
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 3
      cu(i,j,1) = 0.
      cu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      cu(1,j+ngx,ngy+1) = .5*xj0
      cu(2,j+ngx,ngy+1) = .5*yj0
      cu(3,j+ngx,ngy+1) = .5*zj0
      cu(1,j+ngx,ny-ngy+1) = .5*xj0
      cu(2,j+ngx,ny-ngy+1) = .5*yj0
      cu(3,j+ngx,ny-ngy+1) = .5*zj0
   60 continue
      do 70 i = 1, 3
      cu(i,1,ngy+1) = 0.
      cu(i,nx+1,ngy+1) = 0.
      cu(i,1,ny-ngy+1) = 0.
      cu(i,nx+1,ny-ngy+1) = 0.
   70 continue
      cu(1,ngx+1,ngy+1) = .25*xj0
      cu(2,ngx+1,ngy+1) = .25*yj0
      cu(3,ngx+1,ngy+1) = .25*zj0
      cu(1,nx-ngx+1,ngy+1) = .25*xj0
      cu(2,nx-ngx+1,ngy+1) = .25*yj0
      cu(3,nx-ngx+1,ngy+1) = .25*zj0
      cu(1,ngx+1,ny-ngy+1) = .25*xj0
      cu(2,ngx+1,ny-ngy+1) = .25*yj0
      cu(3,ngx+1,ny-ngy+1) = .25*zj0
      cu(1,nx-ngx+1,ny-ngy+1) = .25*xj0
      cu(2,nx-ngx+1,ny-ngy+1) = .25*yj0
      cu(3,nx-ngx+1,ny-ngy+1) = .25*zj0
      return
      end
      subroutine LSCGUARD22L(cu,xj0,yj0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension cu(2,nxe,nye)
      integer i, j, k, nxg, nyg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx1 = nx + 1
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      cu(1,j+ngx,k+ngy) = xj0
      cu(2,j+ngx,k+ngy) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+ngy) = 0.
      cu(i,nx+1,k+ngy) = 0.
   20 continue
      cu(1,ngx+1,k+ngy) = .5*xj0
      cu(2,ngx+1,k+ngy) = .5*yj0
      cu(1,nx-ngx+1,k+ngy) = .5*xj0
      cu(2,nx-ngx+1,k+ngy) = .5*yj0
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 2
      cu(i,j,1) = 0.
      cu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      cu(1,j+ngx,ngy+1) = .5*xj0
      cu(2,j+ngx,ngy+1) = .5*yj0
      cu(1,j+ngx,ny-ngy+1) = .5*xj0
      cu(2,j+ngx,ny-ngy+1) = .5*yj0
   60 continue
      do 70 i = 1, 2
      cu(i,1,ngy+1) = 0.
      cu(i,nx+1,ngy+1) = 0.
      cu(i,1,ny-ngy+1) = 0.
      cu(i,nx+1,ny-ngy+1) = 0.
   70 continue
      cu(1,ngx+1,ngy+1) = .25*xj0
      cu(2,ngx+1,ngy+1) = .25*yj0
      cu(1,nx-ngx+1,ngy+1) = .25*xj0
      cu(2,nx-ngx+1,ngy+1) = .25*yj0
      cu(1,ngx+1,ny-ngy+1) = .25*xj0
      cu(2,ngx+1,ny-ngy+1) = .25*yj0
      cu(1,nx-ngx+1,ny-ngy+1) = .25*xj0
      cu(2,nx-ngx+1,ny-ngy+1) = .25*yj0
      return
      end
      subroutine LSGUARD2L(q,qi0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic scalar field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real q, qi0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension q(nxe,nye)
      integer j, k, nxg, nyg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx1 = nx + 1
      do 20 k = 2, nyg
      do 10 j = 2, nxg
      q(j+ngx,k+ngy) = qi0
   10 continue
      q(1,k+ngy) = 0.
      q(nx+1,k+ngy) = 0.
      q(ngx+1,k+ngy) = .5*qi0
      q(nx-ngx+1,k+ngy) = .5*qi0
   20 continue
      do 30 j = 1, nx1
      q(j,1) = 0.
      q(j,ny+1) = 0.
   30 continue
      do 40 j = 2, nxg
      q(j+ngx,ngy+1) = .5*qi0
      q(j+ngx,ny-ngy+1) = .5*qi0
   40 continue
      q(1,ngy+1) = 0.
      q(nx+1,ngy+1) = 0.
      q(1,ny-ngy+1) = 0.
      q(nx+1,ny-ngy+1) = 0.
      q(ngx+1,ngy+1) = .25*qi0
      q(nx-ngx+1,ngy+1) = .25*qi0
      q(ngx+1,ny-ngy+1) = .25*qi0
      q(nx-ngx+1,ny-ngy+1) = .25*qi0
      return
      end
      subroutine DBLSIN2C(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx+1
c nyv = third dimension of input array cu, must be >= ny+1
c nx2v = second dimension of output array cu2, must be >= 2*nx
c ny2 = third dimension of output array cu2, must be >= 2*ny
      implicit none
      real cu, cu2
      integer nx, ny, nxv, nyv, nx2v, ny2
      dimension cu(2,nxv,nyv), cu2(2,nx2v,ny2)
c local data
      integer j, k, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 20 k = 1, nys
      do 10 j = 1, nxs
      cu2(1,j+1,k+1) = cu(1,j+1,k+1)
      cu2(2,j+1,k+1) = cu(2,j+1,k+1)
      cu2(1,nx+j+1,k+1) = cu(1,nx-j+1,k+1)
      cu2(2,nx+j+1,k+1) = -cu(2,nx-j+1,k+1)
      cu2(1,j+1,ny+k+1) = -cu(1,j+1,ny-k+1)
      cu2(2,j+1,ny+k+1) = cu(2,j+1,ny-k+1)
      cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
      cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
   10 continue
      cu2(1,1,k+1) = cu(1,1,k+1)
      cu2(2,1,k+1) = 0.
      cu2(1,nx+1,k+1) = cu(1,nx+1,k+1)
      cu2(2,nx+1,k+1) = 0.
      cu2(1,1,k+ny+1) = -cu(1,1,ny-k+1)
      cu2(2,1,k+ny+1) = 0.
      cu2(1,nx+1,k+ny+1) = -cu(1,nx+1,ny-k+1)
      cu2(2,nx+1,k+ny+1) = 0.
   20 continue
      do 30 j = 1, nxs
      cu2(1,j+1,1) = 0.
      cu2(2,j+1,1) = cu(2,j+1,1)
      cu2(1,j+nx+1,1) = 0.
      cu2(2,j+nx+1,1) = -cu(2,nx-j+1,1)
      cu2(1,j+1,ny+1) = 0.
      cu2(2,j+1,ny+1) = cu(2,j+1,ny+1)
      cu2(1,j+nx+1,ny+1) = 0.
      cu2(2,j+nx+1,ny+1) = -cu(2,nx-j+1,ny+1)
   30 continue
      cu2(1,1,1) = 0.
      cu2(2,1,1) = 0.
      cu2(1,nx+1,1) = 0.
      cu2(2,nx+1,1) = 0.
      cu2(1,1,ny+1) = 0.
      cu2(2,1,ny+1) = 0.
      cu2(1,nx+1,ny+1) = 0.
      cu2(2,nx+1,ny+1) = 0.
      return
      end
      subroutine DBLSIN2D(q,q2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates an odd array q2 from an array q, so that
c a 2d sine transform can be performed with a 2d real to complex fft.
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = first dimension of input array q, must be >= nx
c nyv = second dimension of input array q, must be >= ny
c nx2v = first dimension of output array q2, must be >= 2*nx
c ny2 = second dimension of output array q2, must be >= 2*ny
      implicit none
      real q, q2
      integer nx, ny, nxv, nyv, nx2v, ny2
      dimension q(nxv,nyv), q2(nx2v,ny2)
c local data
      integer j, k, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 20 k = 1, nys
      do 10 j = 1, nxs
      q2(j+1,k+1) = q(j+1,k+1)
      q2(nx+j+1,k+1) = -q(nx-j+1,k+1)
      q2(j+1,ny+k+1) = -q(j+1,ny-k+1)
      q2(nx+j+1,ny+k+1) = q(nx-j+1,ny-k+1)
   10 continue
      q2(1,k+1) = 0.
      q2(nx+1,k+1) = 0.
      q2(1,k+ny+1) = 0.
      q2(nx+1,k+ny+1) = 0.
   20 continue
      do 30 j = 1, nx
      q2(j,1) = 0.
      q2(j+nx,1) = 0.
      q2(j,ny+1) = 0.
      q2(j+nx,ny+1) = 0.
   30 continue
      return
      end
      subroutine DBLSIN2B(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in y,
c y component is an odd function in x, and the z component is an odd
c function in both x and y.  Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx+1
c nyv = third dimension of input array cu, must be >= ny+1
c nx2v = second dimension of output array cu2, must be >= 2*nx
c ny2 = third dimension of output array cu2, must be >= 2*ny
      implicit none
      real cu, cu2
      integer nx, ny, nxv, nyv, nx2v, ny2
      dimension cu(3,nxv,nyv), cu2(3,nx2v,ny2)
c local data
      integer j, k, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 20 k = 1, nys
      do 10 j = 1, nxs
      cu2(1,j+1,k+1) = cu(1,j+1,k+1)
      cu2(2,j+1,k+1) = cu(2,j+1,k+1)
      cu2(3,j+1,k+1) = cu(3,j+1,k+1)
      cu2(1,nx+j+1,k+1) = cu(1,nx-j+1,k+1)
      cu2(2,nx+j+1,k+1) = -cu(2,nx-j+1,k+1)
      cu2(3,nx+j+1,k+1) = -cu(3,nx-j+1,k+1)
      cu2(1,j+1,ny+k+1) = -cu(1,j+1,ny-k+1)
      cu2(2,j+1,ny+k+1) = cu(2,j+1,ny-k+1)
      cu2(3,j+1,ny+k+1) = -cu(3,j+1,ny-k+1)
      cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
      cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
      cu2(3,nx+j+1,ny+k+1) = cu(3,nx-j+1,ny-k+1)
   10 continue
      cu2(1,1,k+1) = cu(1,1,k+1)
      cu2(2,1,k+1) = 0.
      cu2(3,1,k+1) = 0.
      cu2(1,nx+1,k+1) = cu(1,nx+1,k+1)
      cu2(2,nx+1,k+1) = 0.
      cu2(3,nx+1,k+1) = 0.
      cu2(1,1,k+ny+1) = -cu(1,1,ny-k+1)
      cu2(2,1,k+ny+1) = 0.
      cu2(3,1,k+ny+1) = 0.
      cu2(1,nx+1,k+ny+1) = -cu(1,nx+1,ny-k+1)
      cu2(2,nx+1,k+ny+1) = 0.
      cu2(3,nx+1,k+ny+1) = 0.
   20 continue
      do 30 j = 1, nxs
      cu2(1,j+1,1) = 0.
      cu2(2,j+1,1) = cu(2,j+1,1)
      cu2(3,j+1,1) = 0.
      cu2(1,j+nx+1,1) = 0.
      cu2(2,j+nx+1,1) = -cu(2,nx-j+1,1)
      cu2(3,j+nx+1,1) = 0.
      cu2(1,j+1,ny+1) = 0.
      cu2(2,j+1,ny+1) = cu(2,j+1,ny+1)
      cu2(3,j+1,ny+1) = 0.
      cu2(1,j+nx+1,ny+1) = 0.
      cu2(2,j+nx+1,ny+1) = -cu(2,nx-j+1,ny+1)
      cu2(3,j+nx+1,ny+1) = 0.
   30 continue
      cu2(1,1,1) = 0.
      cu2(2,1,1) = 0.
      cu2(3,1,1) = 0.
      cu2(1,nx+1,1) = 0.
      cu2(2,nx+1,1) = 0.
      cu2(3,nx+1,1) = 0.
      cu2(1,1,ny+1) = 0.
      cu2(2,1,ny+1) = 0.
      cu2(3,1,ny+1) = 0.
      cu2(1,nx+1,ny+1) = 0.
      cu2(2,nx+1,ny+1) = 0.
      cu2(3,nx+1,ny+1) = 0.
      return
      end
      subroutine HAFDBL2C(fxy,fxy2,nx,ny,nxe,nye,nx2v,ny2)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx/ny = system length in x/y direction
c nxe = second dimension of output array fxy, must be >= nx+1
c nye = third dimension of ouput array fxy, must be >= ny+1
c nx2v = second dimension of input array fxy2, must be >= 2*nx
c ny2 = third dimension of input array fxy2, must be >= 2*ny
      implicit none
      real fxy, fxy2
      integer nx, ny, nxe, nye, nx2v, ny2
      dimension fxy(2,nxe,nye), fxy2(2,nx2v,ny2)
c local data
      integer j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      fxy(1,j,k) = fxy2(1,j,k)
      fxy(2,j,k) = fxy2(2,j,k)
   10 continue
   20 continue
      return
      end
      subroutine HAFDBL2D(q,q2,nx,ny,nxe,nye,nx2v,ny2)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array q, must be >= nx+1
c nye = second dimension of ouput array q, must be >= ny+1
c nx2v = first dimension of input array q2, must be >= 2*nx
c ny2 = second dimension of input array q2, must be >= 2*ny
      implicit none
      real q, q2
      integer nx, ny, nxe, nye, nx2v, ny2
      dimension q(nxe,nye), q2(nx2v,ny2)
c local data
      integer j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      q(j,k) = q2(j,k)
   10 continue
   20 continue
      return
      end
      subroutine HAFDBL2B(bxy,bxy2,nx,ny,nxe,nye,nx2v,ny2)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx/ny = system length in x/y direction
c nxe = second dimension of output array bxy, must be >= nx+1
c nye = third dimension of ouput array bxy, must be >= ny+1
c nx2v = second dimension of input array bxy2, must be >= 2*nx
c ny2 = third dimension of input array bxy2, must be >= 2*ny
      implicit none
      real bxy, bxy2
      integer nx, ny, nxe, nye, nx2v, ny2
      dimension bxy(3,nxe,nye), bxy2(3,nx2v,ny2)
c local data
      integer j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      bxy(1,j,k) = bxy2(1,j,k)
      bxy(2,j,k) = bxy2(2,j,k)
      bxy(3,j,k) = bxy2(3,j,k)
   10 continue
   20 continue
      return
      end
      subroutine POISDX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,ny2d
     1)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nx2v,ny2d, output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fx,fy,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fx,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fy
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
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
c ffd(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffd(2*j-1,k) = potential green's function g for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nx2v = first dimension of field arrays, must be >= 2*nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp
      dimension q(nx2v,ny2d), fx(nx2v,ny2d), fy(nx2v,ny2d)
      dimension ffd(nx2v,ny)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      if (at3.eq.0.) then
         ffd(2*j,k) = 1.0
         ffd(2*j-1,k) = affp
      else
         ffd(2*j,k) = exp(-.5*((dkx*ax)**2 + at2))
         ffd(2*j-1,k) = affp*ffd(2*j,k)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 50 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ffd(2*j-1,k)*ffd(2*j,k)
      at3 = -at1*q(2*j-1,k)
      at2 = dnx*float(j - 1)*at3
      at3 = dky*at3
      fx(2*j-1,k) = 0.0
      fx(2*j,k) = at2
      fx(2*j-1,k1) = 0.0
      fx(2*j,k1) = -at2
      fy(2*j-1,k) = 0.0
      fy(2*j,k) = at3
      fy(2*j-1,k1) = 0.0
      fy(2*j,k1) = at3
      wp = wp + at1*q(2*j-1,k)**2
   40 continue
      fx(1,k) = 0.0
      fx(2,k) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
      fy(1,k) = 0.0
      fy(2,k) = 0.0
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
   50 continue
      do 60 j = 1, nx
      fx(2*j-1,1) = 0.0
      fx(2*j,1) = 0.0
      fx(2*j-1,ny+1) = 0.0
      fx(2*j,ny+1) = 0.0
      fy(2*j-1,1) = 0.0
      fy(2*j,1) = 0.0
      fy(2*j-1,ny+1) = 0.0
      fy(2*j,ny+1) = 0.0
   60 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
      do 90 k = 2, ny
      k1 = ny2 - k
      do 80 j = 2, nx
      at2 = ffd(2*j-1,k)
      at1 = at2*ffd(2*j,k)
      at3 = at2*q(2*j-1,k)
      fx(2*j-1,k) = at3
      fx(2*j,k) = 0.0
      fx(2*j-1,k1) = -at3
      fx(2*j,k1) = 0.0
      wp = wp + at1*q(2*j-1,k)**2
   80 continue
      fx(1,k) = 0.0
      fx(2,k) = 0.0
      fx(1,k1) = 0.0
      fx(2,k1) = 0.0
   90 continue
      do 100 j = 1, nx
      fx(2*j-1,1) = 0.0
      fx(2*j,1) = 0.0
      fx(2*j-1,ny+1) = 0.0
      fx(2*j,ny+1) = 0.0
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  110 do 130 k = 2, ny
      k1 = ny2 - k
      do 120 j = 2, nx
      at1 = ffd(2*j,k)
      at2 = at1*q(2*j-1,k)
      fy(2*j-1,k) = at2
      fy(2*j,k) = 0.0
      fy(2*j-1,k1) = -at2
      fy(2*j,k1) = 0.0
  120 continue
      fy(1,k) = 0.0
      fy(2,k) = 0.0
      fy(1,k1) = 0.0
      fy(2,k1) = 0.0
  130 continue
      do 140 j = 1, nx
      fy(2*j-1,1) = 0.0
      fy(2*j,1) = 0.0
      fy(2*j-1,ny+1) = 0.0
      fy(2*j,ny+1) = 0.0
  140 continue
      return
      end
      subroutine POISD2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye,nx
     12v)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nx2v,ny2d, output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,nx2v, output: fx,fy,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffd,isign,nx,ny,nxe,ny,nx2v, output: fx,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffd,isign,nx,ny,nxe,ny,nx2v, output: fy
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(j,k) = charge density for fourier mode (j-1,k-1)
c fx(j,k) = x component of force/charge,
c fy(j,k) = y component of force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c ffd(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffd(2*j-1,k) = potential green's function g for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nx2v = first dimension of form factor array, must be >= 2*nx
      double precision wp
      dimension q(nxe,nye), fx(nxe,nye), fy(nxe,nye)
      dimension ffd(nx2v,ny)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      if (at3.eq.0.) then
         ffd(2*j,k) = 1.0
         ffd(2*j-1,k) = affp
      else
         ffd(2*j,k) = exp(-.5*((dkx*ax)**2 + at2))
         ffd(2*j-1,k) = affp*ffd(2*j,k)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 50 k = 2, ny
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ffd(2*j-1,k)*ffd(2*j,k)
      at3 = -at1*q(j,k)
      at2 = dnx*float(j - 1)*at3
      at3 = dky*at3
      fx(j,k) = at2
      fy(j,k) = at3
      wp = wp + at1*q(j,k)**2
   40 continue
      fx(1,k) = 0.0
      fy(1,k) = 0.0
      fx(nx+1,k) = 0.0
      fy(nx+1,k) = 0.0
   50 continue
      do 60 j = 1, nx1
      fx(j,1) = 0.0
      fx(j,ny+1) = 0.0
      fy(j,1) = 0.0
      fy(j,ny+1) = 0.0
   60 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
      do 90 k = 2, ny
      do 80 j = 2, nx
      at2 = ffd(2*j-1,k)
      at1 = at2*ffd(2*j,k)
      at3 = at2*q(j,k)
      fx(j,k) = at3
      wp = wp + at1*q(j,k)**2
   80 continue
      fx(1,k) = 0.0
      fx(nx+1,k) = 0.0
   90 continue
      do 100 j = 1, nx1
      fx(j,1) = 0.0
      fx(j,ny+1) = 0.0
  100 continue
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  110 do 130 k = 2, ny
      do 120 j = 2, nx
      at1 = ffd(2*j,k)
      at2 = at1*q(j,k)
      fy(j,k) = at2
  120 continue
      fy(1,k) = 0.0
      fy(nx+1,k) = 0.0
  130 continue
      do 140 j = 1, nx1
      fy(j,1) = 0.0
      fy(j,ny+1) = 0.0
  140 continue
      return
      end
      subroutine POISDX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,nxv,ny2d, output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(nxv,ny2d), fxy(2,nxv,ny2d)
      dimension ffd(nxv,ny)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.0)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
      do 50 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffd(j,k))*aimag(ffd(j,k))
      at3 = -at1*real(q(j,k))
      at2 = dnx*float(j - 1)*at3
      at3 = dky*at3
      fxy(1,j,k) = cmplx(0.0,at2)
      fxy(2,j,k) = cmplx(0.0,at3)
      fxy(1,j,k1) = cmplx(0.0,-at2)
      fxy(2,j,k1) = cmplx(0.0,at3)
      wp = wp + at1*real(q(j,k))**2
   40 continue
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
   50 continue
      do 60 j = 1, nx
      fxy(1,j,1) = zero
      fxy(2,j,1) = zero
      fxy(1,j,ny+1) = zero
      fxy(2,j,ny+1) = zero
   60 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
      subroutine POISD22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye,nxv
     1)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nx2v,ny2d, output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c q(j,k) = charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of force/charge,
c fxy(2,j,k) = y component of force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nxv = first dimension of form factor array, must be >= nx
      double precision wp
      complex ffd
      dimension q(nxe,nye), fxy(2,nxe,nye)
      dimension ffd(nxv,ny)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.0)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
      do 50 k = 2, ny
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffd(j,k))*aimag(ffd(j,k))
      at3 = -at1*q(j,k)
      at2 = dnx*float(j - 1)*at3
      at3 = dky*at3
      fxy(1,j,k) = at2
      fxy(2,j,k) = at3
      wp = wp + at1*q(j,k)**2
   40 continue
      fxy(1,1,k) = 0.0
      fxy(2,1,k) = 0.0
      fxy(1,nx+1,k) = 0.0
      fxy(2,nx+1,k) = 0.0
   50 continue
      do 60 j = 1, nx1
      fxy(1,j,1) = 0.0
      fxy(2,j,1) = 0.0
      fxy(1,j,ny+1) = 0.0
      fxy(2,j,ny+1) = 0.0
   60 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
      subroutine POISDX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential).
c Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,nxv,ny2d, output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp
      complex q, fxy, ffd, zero
      dimension q(nxv,ny2d), fxy(3,nxv,ny2d)
      dimension ffd(nxv,ny)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.0)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
      do 50 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffd(j,k))*aimag(ffd(j,k))
      at3 = -at1*real(q(j,k))
      at2 = dnx*float(j - 1)*at3
      at3 = dky*at3
      fxy(1,j,k) = cmplx(0.0,at2)
      fxy(2,j,k) = cmplx(0.0,at3)
      fxy(3,j,k) = zero
      fxy(1,j,k1) = cmplx(0.0,-at2)
      fxy(2,j,k1) = cmplx(0.0,at3)
      fxy(3,j,k1) = zero
      wp = wp + at1*real(q(j,k))**2
   40 continue
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
   50 continue
      do 60 j = 1, nx
      fxy(1,j,1) = zero
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,ny+1) = zero
      fxy(2,j,ny+1) = zero
      fxy(3,j,ny+1) = zero
   60 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
      subroutine POISD23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxe,nye,nxv
     1)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential).
c Zeros out z component
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nx2v,ny2d, output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fxy,we
c approximate flop count is: 11*nx*ny
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c q(j,k) = charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of force/charge,
c fxy(2,j,k) = y component of force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
c nxv = first dimension of form factor array, must be >= nx
      double precision wp
      complex ffd
      dimension q(nxe,nye), fxy(3,nxe,nye)
      dimension ffd(nxv,ny)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.0)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
      do 50 k = 2, ny
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffd(j,k))*aimag(ffd(j,k))
      at3 = -at1*q(j,k)
      at2 = dnx*float(j - 1)*at3
      at3 = dky*at3
      fxy(1,j,k) = at2
      fxy(2,j,k) = at3
      fxy(3,j,k) = 0.0
      wp = wp + at1*q(j,k)**2
   40 continue
      fxy(1,1,k) = 0.0
      fxy(2,1,k) = 0.0
      fxy(3,1,k) = 0.0
      fxy(1,nx+1,k) = 0.0
      fxy(2,nx+1,k) = 0.0
      fxy(3,nx+1,k) = 0.0
   50 continue
      do 60 j = 1, nx1
      fxy(1,j,1) = 0.0
      fxy(2,j,1) = 0.0
      fxy(3,j,1) = 0.0
      fxy(1,j,ny+1) = 0.0
      fxy(2,j,ny+1) = 0.0
      fxy(3,j,ny+1) = 0.0
   60 continue
      we = 2.0*float(nx*ny)*wp
      return
      end
      subroutine POISDX2T(qt,fxt,fyt,isign,ffdt,ax,ay,affp,we,nx,ny,nxv,
     1ny2d,nyd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-sin, sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffdt
c for isign = -1, input: qt,ffdt,isign,nx,ny,nxv,ny2d,
c output: fxt,fyt,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fx,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffc,isign,nx,ny,nxv,nyhd, output: fy
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c cmplx(qt(2*k-1,j),q(2*k,j)) = complex charge density
c for fourier mode (k-1,j-1)
c cmplx(fxt(2*k-1,j),fxt(2*k,j))) = x component of complex force/charge,
c cmplx(fyt(2*k-1,j),fyt(2*k,j))) = y component of complex force/charge,
c for fourier mode (k-1,j-1)
c if isign = 0, form factor array is prepared
c ffdt(2*k,j) = finite-size particle shape factor s
c for fourier mode (k-1,j-1)
c ffdt(2*k-1,j) = potential green's function g for fourier mode
c (k-1,j-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = second dimension of field arrays, must be >= nx+1
c ny2d = first dimension of field arrays, must be >= 2*ny
c nyd = first dimension of form factor array, must be >= ny
      implicit none
      integer isign, nx, ny, nxv, ny2d, nyd
      real ax, ay, affp, we
      real qt, fxt, fyt, ffdt
      dimension qt(2*ny2d,nxv), fxt(2*ny2d,nxv), fyt(2*ny2d,nxv)
      dimension ffdt(2*nyd,nx)
c local data
      integer j, k, ny2, j1, k1
      real dnx, dny, dkx, dky, at1, at2, at3
      double precision wp
      ny2 = 2*ny + 2
      j1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 j = 1, nx
      dkx = dnx*float(j - 1)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      if (at3.eq.0.) then
         ffdt(2*k,j) = 1.0
         ffdt(2*k-1,j) = affp
      else
         ffdt(2*k,j) = exp(-.5*((dky*ay)**2 + at2))
         ffdt(2*k-1,j) = affp*ffdt(2*k,j)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 50 j = 2, nx
      dkx = dnx*float(j - 1)
      do 40 k = 2, ny
      k1 = ny2 - k
      at1 = ffdt(2*k-1,j)*ffdt(2*k,j)
      at3 = -at1*qt(2*k-1,j)
      at2 = dkx*at3
      at3 = dny*float(k - 1)*at3
      fxt(2*k-1,j) = 0.0
      fxt(2*k,j) = at2
      fxt(2*k1-1,j) = 0.0
      fxt(2*k1,j) = -at2
      fyt(2*k-1,j) = 0.0
      fyt(2*k,j) = at3
      fyt(2*k1-1,j) = 0.0
      fyt(2*k1,j) = at3
      wp = wp + at1*qt(2*k-1,j)**2
   40 continue
      fxt(1,j) = 0.0
      fxt(2,j) = 0.0
      fxt(ny2-1,j) = 0.0
      fxt(ny2,j) = 0.0
      fyt(1,j) = 0.0
      fyt(2,j) = 0.0
      fyt(ny2-1,j) = 0.0
      fyt(ny2,j) = 0.0
   50 continue
      do 60 k = 2, ny
      k1 = ny2 - k
      fxt(2*k-1,1) = 0.0
      fxt(2*k,1) = 0.0
      fxt(2*k1-1,1) = 0.0
      fxt(2*k1,1) = 0.0
      fyt(2*k-1,1) = 0.0
      fyt(2*k,1) = 0.0
      fyt(2*k1-1,1) = 0.0
      fyt(2*k1,1) = 0.0
      fxt(2*k-1,j1) = 0.0
      fxt(2*k,j1) = 0.0
      fxt(2*k1-1,j1) = 0.0
      fxt(2*k1,j1) = 0.0
      fyt(2*k-1,j1) = 0.0
      fyt(2*k,j1) = 0.0
      fyt(2*k1-1,j1) = 0.0
      fyt(2*k1,j1) = 0.0
   60 continue
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
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
      do 90 j = 2, nx
      do 80 k = 2, ny
      k1 = ny2 - k
      at2 = ffdt(2*k-1,j)
      at1 = at2*ffdt(2*k,j)
      at3 = at2*qt(2*k-1,j)
      fxt(2*k-1,j) = at3
      fxt(2*k,j) = 0.0
      fxt(2*k1-1,j) = -at3
      fxt(2*k1,j) = 0.0
      wp = wp + at1*qt(2*k-1,j)**2
   80 continue
      fxt(1,j) = 0.0
      fxt(2,j) = 0.0
      fxt(ny2-1,j) = 0.0
      fxt(ny2,j) = 0.0
   90 continue
      do 100 k = 2, ny
      k1 = ny2 - k
      fxt(2*k-1,1)= 0.0
      fxt(2*k,1) = 0.0
      fxt(2*k1-1,1) = 0.0
      fxt(2*k1,1) = 0.0
      fxt(2*k-1,j1)= 0.0
      fxt(2*k,j1) = 0.0
      fxt(2*k1-1,j1) = 0.0
      fxt(2*k1,j1) = 0.0
  100 continue
      fxt(1,1) = 0.0
      fxt(2,1) = 0.0
      fxt(ny2-1,1) = 0.0
      fxt(ny2,1) = 0.0
      fxt(1,j1) = 0.0
      fxt(2,j1) = 0.0
      fxt(ny2-1,j1) = 0.0
      fxt(ny2,j1) = 0.0
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  110 do 130 j = 2, nx
      do 120 k = 2, ny
      k1 = ny2 - k
      at1 = ffdt(2*k,j)
      fyt(2*k-1,j) = at1*qt(2*k-1,j)
      fyt(2*k,j)= at1*qt(2*k,j)
      fyt(2*k1-1,j) = at1*qt(2*k1-1,j)
      fyt(2*k1,j) = at1*qt(2*k1,j)
  120 continue
      fyt(1,j) = 0.0
      fyt(2,j) = 0.0
      fyt(ny2-1,j) = 0.0
      fyt(ny2,j) = 0.0
  130 continue
      do 140 k = 2, ny
      k1 = ny2 - k
      fyt(2*k-1,1) = 0.0
      fyt(2*k,1) = 0.0
      fyt(2*k1-1,1) = 0.0
      fyt(2*k1,1) = 0.0
      fyt(2*k-1,j1) = 0.0
      fyt(2*k,j1) = 0.0
      fyt(2*k1-1,j1) = 0.0
      fyt(2*k1,j1) = 0.0
  140 continue
      fyt(1,1) = 0.0
      fyt(2,1) = 0.0
      fyt(ny2-1,1) = 0.0
      fyt(ny2,1) = 0.0
      fyt(1,j1) = 0.0
      fyt(2,j1) = 0.0
      fyt(ny2-1,j1) = 0.0
      fyt(ny2,j1) = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine POISD2T(qt,fxt,fyt,isign,ffdt,ax,ay,affp,we,nx,ny,nxe,n
     1ye,ny2v)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c mode ny+1 is zero and is not written
c for isign = 0, input: isign,ax,ay,affp,nx,ny,ny2v, output: ffdt
c for isign = -1, input: q,ffdt,isign,nx,ny,nxe,nye,ny2v,
c output: fxt,fyt,we
c approximate flop count is: 11*nx*ny
c for isign = 1, input: q,ffdt,isign,nx,ny,nxe,nye,ny2v, output: fxt,we
c approximate flop count is: 6*nx*ny
c for isign = 2, input: q,ffdt,isign,nx,ny,nxe,nye,ny2v, output: fyt
c approximate flop count is: 2*nx*ny
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)*s(kx,ky)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c q(k,j) = charge density for fourier mode (j-1,k-1)
c fx(k,j) = x component of force/charge,
c fy(k,j) = y component of force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c ffdt(2*k,j) = finite-size particle shape factor s
c ffdt(2*k-1,j) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = 2*nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = second dimension of field arrays, must be >= nx+1
c nye = first dimension of field arrays, must be = ny
c ny2v = first dimension of form factor array, must be >= 2*ny
      implicit none
      integer isign, nx, ny, nxe, nye, ny2v
      real ax, ay, affp, we
      real qt, fxt, fyt, ffdt
      dimension qt(nye,nxe), fxt(nye,nxe), fyt(nye,nxe)
      dimension ffdt(ny2v,nx)
c local data
      integer j, k, nx1
      real dnx, dny, dkx, dky, at1, at2, at3
      double precision wp
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 j = 1, nx
      dkx = dnx*float(j - 1)
      at1 = dkx*dkx
      at2 = (dkx*ax)**2
      do 10 k = 1, ny
      dky = dny*float(k - 1)
      at3 = dky*dky + at1
      if (at3.eq.0.) then
         ffdt(2*k,j) = 1.0
         ffdt(2*k-1,j) = affp
      else
         ffdt(2*k,j) = exp(-.5*((dky*ay)**2 + at2))
         ffdt(2*k-1,j) = affp*ffdt(2*k,j)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
      do 50 j = 2, nx
      dkx = dnx*float(j - 1)
      do 40 k = 2, ny
      at1 = ffdt(2*k-1,j)*ffdt(2*k,j)
      at3 = -at1*qt(k,j)
      at2 = dkx*at3
      at3 = dny*float(k - 1)*at3
      fxt(k,j) = at2
      fyt(k,j) = at3
      wp = wp + at1*qt(k,j)**2
   40 continue
      fxt(1,j) = 0.0
c     fxt(ny+1,j) = 0.0
      fyt(1,j) = 0.0
c     fyt(ny+1,j) = 0.0
   50 continue
      do 60 k = 1, ny
      fxt(k,1) = 0.0
      fyt(k,1) = 0.0
      fxt(k,nx+1) = 0.0
      fyt(k,nx+1) = 0.0
   60 continue
c     fxt(ny+1,1) = 0.0
c     fyt(ny+1,1) = 0.0
c     fxt(ny+1,nx+1) = 0.0
c     fyt(ny+1,nx+1) = 0.0
      we = 2.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
      do 90 j = 2, nx
      do 80 k = 2, ny
      at2 = ffdt(2*k-1,j)
      at1 = at2*ffdt(2*k,j)
      at3 = at2*qt(k,j)
      fxt(k,j) = at3
      wp = wp + at1*qt(k,j)**2
   80 continue
      fxt(1,j) = 0.0
c     fxt(ny+1,j) = 0.0
   90 continue
      do 100 k = 1, ny
      fxt(k,1) = 0.0
      fxt(k,nx+1) = 0.0
  100 continue
c     fxt(ny+1,1) = 0.0
c     fxt(ny+1,nx+1) = 0.0
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  110 do 130 k = 2, ny
      do 120 j = 2, nx
      at1 = ffdt(2*k,j)
      at2 = at1*qt(k,j)
      fyt(k,j) = at2
  120 continue
      fyt(1,j) = 0.0
c     fyt(ny+1,j) = 0.0
  130 continue
      do 140 k = 1, ny
      fyt(k,1) = 0.0
      fyt(k,nx+1) = 0.0
  140 continue
c     fyt(ny+1,1) = 0.0
c     fyt(ny+1,nx+1) = 0.0
      return
      end
      subroutine DIVFD2(f,df,nx,ny,ndim,nxv,nyv)
c this subroutine calculates the divergence in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the charge density from the electric field
c input: all except df, output: df
c approximate flop count is: 6*nxc*nyc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      real f, df
      dimension f(ndim,nxv,nyv), df(nxv,nyv)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the gradient
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      df(j,k) = -(dkx*f(1,j,k) + dky*f(2,j,k))
   10 continue
c mode numbers kx = 0, nx
      df(1,k) = 0.
      df(nx+1,k) = 0.
   20 continue
c mode numbers ky = 0, ny
      do 30 j = 1, nx1
      df(j,1) = 0.
      df(j,ny+1) = 0.
   30 continue
      return
      end
      subroutine GRADFD2(df,f,nx,ny,ndim,nxv,nyv)
c this subroutine calculates the gradient in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the electric field from the potential
c input: all except f, output: f
c approximate flop count is: 4*nxc*nyc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      real df, f
      dimension df(nxv,nyv), f(ndim,nxv,nyv)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the gradient
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      f(1,j,k) = dkx*df(j,k)
      f(2,j,k) = dky*df(j,k)
   10 continue
c mode numbers kx = 0, nx
      f(1,1,k) = 0.
      f(2,1,k) = 0.
      f(1,nx+1,k) = 0.
      f(2,nx+1,k) = 0.
   20 continue
c mode numbers ky = 0, ny
      do 30 j = 1, nx1
      f(1,j,1) = 0.
      f(2,j,1) = 0.
      f(1,j,ny+1) = 0.
      f(2,j,ny+1) = 0.
   30 continue
      if (ndim.eq.2) return
c handle case of ndim = 3
      do 50 k = 2, ny
      do 40 j = 2, nx
      f(3,j,k) = 0.0
   40 continue
      f(3,1,k) = 0.0
      f(3,nx+1,k) = 0.0
   50 continue
      do 60 j = 1, nx1
      f(3,j,1) = 0.0
      f(3,j,ny+1) = 0.0
   60 continue
      return
      end
      subroutine CURLFD2(f,g,nx,ny,nxv,nyv)
c this subroutine calculates the gradient in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 8*nxc*nyc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      real f, g
      dimension f(3,nxv,nyv), g(3,nxv,nyv)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the curl
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      g(1,j,k) = dky*f(3,j,k)
      g(2,j,k) = -dkx*f(3,j,k)
      g(3,j,k) = dkx*f(2,j,k) - dky*f(1,j,k)
   10 continue
   20 continue
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      dky = dny*float(k - 1)
      g(1,1,k) = 0.
      g(2,1,k) = 0.
      g(3,1,k) = -dky*f(1,1,k)
      g(1,nx+1,k) = 0.
      g(2,nx+1,k) = 0.
      g(3,nx+1,k) = 0.
   30 continue
c mode numbers ky = 0, ny
      do 40 j = 2, nx
      g(1,j,1) = 0.
      g(2,j,1) = 0.
      g(3,j,1) = dkx*f(2,j,1)
      g(1,j,ny+1) = 0.
      g(2,j,ny+1) = 0.
      g(3,j,ny+1) = 0.
   40 continue
      g(1,1,1) = 0.
      g(2,1,1) = 0.
      g(3,1,1) = 0.
      g(1,nx+1,1) = 0.
      g(2,nx+1,1) = 0.
      g(3,nx+1,1) = 0.
      g(1,1,ny+1) = 0.
      g(2,1,ny+1) = 0.
      g(3,1,ny+1) = 0.
      g(1,nx+1,ny+1) = 0.
      g(2,nx+1,ny+1) = 0.
      g(3,nx+1,ny+1) = 0.
      return
      end
      subroutine CURLFD22(f,g,nx,ny,nxv,nyv)
c this subroutine calculates the gradient in fourier space
c with dirichlet boundary conditions (zero potential)
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 8*nxc*nyc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      real f, g
      dimension f(2,nxv,nyv), g(nxv,nyv)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate the curl
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      g(j,k) = dkx*f(2,j,k) - dky*f(1,j,k)
   10 continue
   20 continue
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      dky = dny*float(k - 1)
      g(1,k) = -dky*f(1,1,k)
      g(nx+1,k) = 0.
   30 continue
c mode numbers ky = 0, ny
      do 40 j = 2, nx
      g(j,1) = dkx*f(2,j,1)
      g(j,ny+1) = 0.
   40 continue
      g(1,1) = 0.
      g(nx+1,1) = 0.
      g(1,ny+1) = 0.
      g(nx+1,ny+1) = 0.
      return
      end
      subroutine CUPERPDX2(cu,nx,ny,nxv,ny2d)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      complex cu, zero
      dimension cu(3,nxv,ny2d)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = at1*(dkx*aimag(cu(1,j,k)) + dky*aimag(cu(2,j,k)))
      at3 = aimag(cu(1,j,k)) - dkx*at2
      at4 = aimag(cu(2,j,k)) - dky*at2
      cu(1,j,k) = cmplx(0.,at3)
      cu(2,j,k) = cmplx(0.,at4)
      cu(1,j,k1) = cmplx(0.,-at3)
      cu(2,j,k1) = cmplx(0.,at4)
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   20 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx
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
      subroutine CUPERPD2(cu,nx,ny,nxe,nye)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxe = first dimension of current array, must be >= nx+1
c nye = second dimension of current array, must be >= ny+1
      dimension cu(3,nxe,nye)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      at3 = cu(1,j,k) - dkx*at2
      at4 = cu(2,j,k) - dky*at2
      cu(1,j,k) = at3
      cu(2,j,k) = at4
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = 0.
   20 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx1
      cu(1,j,1) = 0.
   40 continue
      cu(2,1,1) = 0.
      cu(2,1,ny+1) = 0.
      return
      end
      subroutine CUPERPDX22(cu,nx,ny,nxv,ny2d)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      complex cu, zero
      dimension cu(2,nxv,ny2d)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = at1*(dkx*aimag(cu(1,j,k)) + dky*aimag(cu(2,j,k)))
      at3 = aimag(cu(1,j,k)) - dkx*at2
      at4 = aimag(cu(2,j,k)) - dky*at2
      cu(1,j,k) = cmplx(0.,at3)
      cu(2,j,k) = cmplx(0.,at4)
      cu(1,j,k1) = cmplx(0.,-at3)
      cu(2,j,k1) = cmplx(0.,at4)
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   20 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx
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
      subroutine CUPERPD22(cu,nx,ny,nxe,nye)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c input: all, output: cu
c approximate flop count is: 13*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxe = first dimension of current array, must be >= nx+1
c nye = second dimension of current array, must be >= ny+1
      dimension cu(2,nxe,nye)
      nx1 = nx + 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      at3 = cu(1,j,k) - dkx*at2
      at4 = cu(2,j,k) - dky*at2
      cu(1,j,k) = at3
      cu(2,j,k) = at4
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = 0.
   20 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx1
      cu(1,j,1) = 0.
   40 continue
      cu(2,1,1) = 0.
      cu(2,1,ny+1) = 0.
      return
      end
      subroutine BPOISDX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxv,n
     1y2d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 20*nx*ny + 11*(nx + ny)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 14*nx*ny + 8*(nx + ny)
c for isign = 2, input: cu,ffd,isign,nx,ny,nxv,ny2d, output: bxy
c approximate flop count is: 8*nx*ny
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
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
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp
      complex cu, bxy, ffd, zero
      dimension cu(3,nxv,ny2d), bxy(3,nxv,ny2d)
      dimension ffd(nxv,ny)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 50 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffd(j,k))*aimag(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      bxy(1,j,k) = cmplx(0.,at3*real(cu(3,j,k)))
      bxy(2,j,k) = cmplx(0.,-at2*real(cu(3,j,k)))
      bxy(3,j,k) = cmplx(at3*aimag(cu(1,j,k))-at2*aimag(cu(2,j,k)),0.)
      bxy(1,j,k1) = bxy(1,j,k)
      bxy(2,j,k1) = -bxy(2,j,k)
      bxy(3,j,k1) = bxy(3,j,k)
      wp = wp + at1*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2 + real(cu
     1(3,j,k))**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      at3 = dky*at1
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = cmplx(at3*aimag(cu(1,1,k)),0.)      
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(aimag(cu(1,1,k))**2 + real(cu(3,1,k))**2)
   60 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 70 j = 2, nx
      at1 = ci2*real(ffd(j,1))*aimag(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      bxy(1,j,1) = zero
      bxy(2,j,1) = zero
      bxy(3,j,1) = cmplx(-at2*aimag(cu(2,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(2,j,1))**2 + real(cu(3,j,1))**2)
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 100 k = 2, ny
      k1 = ny2 - k
      do 90 j = 2, nx
      at2 = ci2*real(ffd(j,k))
      at1 = at2*aimag(ffd(j,k))
      bxy(1,j,k) = cmplx(0.,at2*aimag(cu(1,j,k)))
      bxy(2,j,k) = cmplx(0.,at2*aimag(cu(2,j,k)))
      bxy(3,j,k) = cmplx(at2*real(cu(3,j,k)),0.)
      bxy(1,j,k1) = -bxy(1,j,k)
      bxy(2,j,k1) = bxy(2,j,k)
      bxy(3,j,k1) = -bxy(3,j,k)
      wp = wp + at1*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2 + real(cu
     1(3,j,k))**2)
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, ny
      k1 = ny2 - k
      at2 = ci2*real(ffd(1,k))
      at1 = at2*aimag(ffd(1,k))
      bxy(1,1,k) = cmplx(0.,at2*aimag(cu(1,1,k)))
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(aimag(cu(1,1,k))**2 + real(cu(3,1,k))**2)
  110 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 120 j = 2, nx
      at2 = ci2*real(ffd(j,1))
      at1 = at2*aimag(ffd(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,at2*aimag(cu(2,j,1)))
      bxy(3,j,1) = zero
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(2,j,1))**2 + real(cu(3,j,1))**2)
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
c mode numbers 0 < kx < nx and 0 < ky < ny
  130 do 150 k = 2, ny
      k1 = ny2 - k
      do 140 j = 2, nx
      at1 = aimag(ffd(j,k))
      bxy(1,j,k) = cmplx(0.,at1*aimag(cu(1,j,k)))
      bxy(2,j,k) = cmplx(0.,at1*aimag(cu(2,j,k)))
      bxy(3,j,k) = cmplx(at1*real(cu(3,j,k)),0.)
      bxy(1,j,k1) = -bxy(1,j,k)
      bxy(2,j,k1) = bxy(2,j,k)
      bxy(3,j,k1) = -bxy(3,j,k)
  140 continue
c mode numbers kx = 0, nx
      at1 = aimag(ffd(1,k))
      bxy(1,1,k) = cmplx(0.,at1*aimag(cu(1,1,k)))
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
  150 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 160 j = 2, nx
      at1 = aimag(ffd(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,at1*aimag(cu(2,j,1)))
      bxy(3,j,1) = zero
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
  160 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      return
      end
      subroutine BPOISD23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxe,ny
     1e,nxv)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 20*nx*ny + 11*(nx + ny)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 14*nx*ny + 8*(nx + ny)
c for isign = 2, input: cu,ffd,isign,nx,ny,nxv,ny2d, output: bxy
c approximate flop count is: 8*nx*ny
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
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
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
      double precision wp
      complex ffd
      dimension cu(3,nxe,nye), bxy(3,nxe,nye)
      dimension ffd(nxv,ny)
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 50 k = 2, ny
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffd(j,k))*aimag(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      bxy(1,j,k) = at3*cu(3,j,k)
      bxy(2,j,k) = -at2*cu(3,j,k)
      bxy(3,j,k) = at2*cu(2,j,k) - at3*cu(1,j,k)
      wp = wp + at1*(cu(1,j,k)**2 + cu(2,j,k)**2 + cu(3,j,k)**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, ny
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      at3 = dky*at1
      bxy(1,1,k) = 0.
      bxy(2,1,k) = 0.
      bxy(3,1,k) = -at3*cu(1,1,k)
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      bxy(3,nx+1,k) = 0.
      wp = wp + at1*(cu(1,1,k)**2 + cu(3,1,k)**2)
   60 continue
c mode numbers ky = 0, ny
      do 70 j = 2, nx
      at1 = ci2*real(ffd(j,1))*aimag(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      bxy(1,j,1) = 0.
      bxy(2,j,1) = 0.
      bxy(3,j,1) = at2*cu(2,j,1)
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      bxy(3,j,ny+1) = 0.
      wp = wp + at1*(cu(2,j,1)**2 + cu(3,j,1)**2)
   70 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(3,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(3,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(3,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      bxy(3,nx+1,ny+1) = 0.
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 100 k = 2, ny
      do 90 j = 2, nx
      at2 = ci2*real(ffd(j,k))
      at1 = at2*aimag(ffd(j,k))
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(2,j,k) = at2*cu(2,j,k)
      bxy(3,j,k) = at2*cu(3,j,k)
      wp = wp + at1*(cu(1,j,k)**2 + cu(2,j,k)**2 + cu(3,j,k)**2)
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, ny
      at2 = ci2*real(ffd(1,k))
      at1 = at2*aimag(ffd(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = 0.
      bxy(3,1,k) = 0.
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      bxy(3,nx+1,k) = 0.
      wp = wp + at1*(cu(1,1,k)**2 + cu(3,1,k))**2
  110 continue
c mode numbers ky = 0, ny
      do 120 j = 2, nx
      at2 = ci2*real(ffd(j,1))
      at1 = at2*aimag(ffd(j,1))
      bxy(1,j,1) = 0.
      bxy(2,j,1) = at2*cu(2,j,1)
      bxy(3,j,1) = 0.
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      bxy(3,j,ny+1) = 0.
      wp = wp + at1*(cu(2,j,1)**2 + cu(3,j,1)**2)
  120 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(3,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(3,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(3,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      bxy(3,nx+1,ny+1) = 0.
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny
  130 do 150 k = 2, ny
      do 140 j = 2, nx
      at1 = aimag(ffd(j,k))
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(2,j,k) = at1*cu(2,j,k)
      bxy(3,j,k) = at1*cu(3,j,k)
  140 continue
c mode numbers kx = 0, nx
      at1 = aimag(ffd(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = 0.
      bxy(3,1,k) = 0.
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      bxy(3,nx+1,k) = 0.
  150 continue
c mode numbers ky = 0, ny
      do 160 j = 2, nx
      at1 = aimag(ffd(j,1))
      bxy(1,j,1) = 0.
      bxy(2,j,1) = at1*cu(2,j,1)
      bxy(3,j,1) = 0.
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      bxy(3,j,ny+1) = 0.
  160 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(3,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(3,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(3,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      bxy(3,nx+1,ny+1) = 0.
      return
      end
      subroutine BPOISDX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nx
     1v,ny2d)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bz,wm
c approximate flop count is: 15*nx*ny + 9*(nx + ny)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 10*nx*ny + 6*(nx + ny)
c for isign = 2, input: cu,ffd,isign,nx,ny,nxv,ny2d, output: bxy
c approximate flop count is: 6*nx*ny
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
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
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp
      complex cu, bxy, bz, ffd, zero
      dimension cu(2,nxv,ny2d), bxy(2,nxv,ny2d), bz(nxv,ny2d)
      dimension ffd(nxv,ny)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 50 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffd(j,k))*aimag(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      bz(j,k) = cmplx(at3*aimag(cu(1,j,k))-at2*aimag(cu(2,j,k)),0.)
      bz(j,k1) = bz(j,k)
      wp = wp + at1*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      at3 = dky*at1
      bz(1,k) = cmplx(at3*aimag(cu(1,1,k)),0.)      
      bz(1,k1) = zero
      wp = wp + at1*aimag(cu(1,1,k))**2
   60 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 70 j = 2, nx
      at1 = ci2*real(ffd(j,1))*aimag(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      bz(j,1) = cmplx(-at2*aimag(cu(2,j,1)),0.)
      bz(j,k1) = zero
      wp = wp + at1*aimag(cu(2,j,1))**2
   70 continue
      bz(1,1) = zero
      bz(1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 100 k = 2, ny
      k1 = ny2 - k
      do 90 j = 2, nx
      at2 = ci2*real(ffd(j,k))
      at1 = at2*aimag(ffd(j,k))
      bxy(1,j,k) = cmplx(0.,at2*aimag(cu(1,j,k)))
      bxy(2,j,k) = cmplx(0.,at2*aimag(cu(2,j,k)))
      bxy(1,j,k1) = -bxy(1,j,k)
      bxy(2,j,k1) = bxy(2,j,k)
      wp = wp + at1*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2)
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, ny
      k1 = ny2 - k
      at2 = ci2*real(ffd(1,k))
      at1 = at2*aimag(ffd(1,k))
      bxy(1,1,k) = cmplx(0.,at2*aimag(cu(1,1,k)))
      bxy(2,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      wp = wp + at1*aimag(cu(1,1,k))**2
  110 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 120 j = 2, nx
      at2 = ci2*real(ffd(j,1))
      at1 = at2*aimag(ffd(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,at2*aimag(cu(2,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      wp = wp + at1*aimag(cu(2,j,1))**2
  120 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny
  130 do 150 k = 2, ny
      k1 = ny2 - k
      do 140 j = 2, nx
      at1 = aimag(ffd(j,k))
      bxy(1,j,k) = cmplx(0.,at1*aimag(cu(1,j,k)))
      bxy(2,j,k) = cmplx(0.,at1*aimag(cu(2,j,k)))
      bxy(1,j,k1) = -bxy(1,j,k)
      bxy(2,j,k1) = bxy(2,j,k)
  140 continue
c mode numbers kx = 0, nx
      at1 = aimag(ffd(1,k))
      bxy(1,1,k) = cmplx(0.,at1*aimag(cu(1,1,k)))
      bxy(2,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
  150 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 160 j = 2, nx
      at1 = aimag(ffd(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,at1*aimag(cu(2,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
  160 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      return
      end
      subroutine BPOISD22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxe
     1,nye,nxv)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 20*nx*ny + 11*(nx + ny)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 14*nx*ny + 8*(nx + ny)
c for isign = 2, input: cu,ffd,isign,nx,ny,nxv,ny2d, output: bxy
c approximate flop count is: 8*nx*ny
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
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
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
      double precision wp
      complex ffd
      dimension cu(2,nxe,nye), bxy(2,nxe,nye), bz(nxe,nye)
      dimension ffd(nxv,ny)
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 50 k = 2, ny
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffd(j,k))*aimag(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      bz(j,k) = at2*cu(2,j,k) - at3*cu(1,j,k)
      wp = wp + at1*(cu(1,j,k)**2 + cu(2,j,k)**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, ny
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      at3 = dky*at1
      bz(1,k) = -at3*cu(1,1,k)
      bz(nx+1,k) = 0.
      wp = wp + at1*cu(1,1,k)**2
   60 continue
c mode numbers ky = 0, ny
      do 70 j = 2, nx
      at1 = ci2*real(ffd(j,1))*aimag(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      bz(j,1) = at2*cu(2,j,1)
      bz(j,ny+1) = 0.
      wp = wp + at1*cu(2,j,1)**2
   70 continue
      bz(1,1) = 0.
      bz(nx+1,1) = 0.
      bz(1,ny+1) = 0.
      bz(nx+1,ny+1) = 0.
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 100 k = 2, ny
      do 90 j = 2, nx
      at2 = ci2*real(ffd(j,k))
      at1 = at2*aimag(ffd(j,k))
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(2,j,k) = at2*cu(2,j,k)
      wp = wp + at1*(cu(1,j,k)**2 + cu(2,j,k)**2)
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, ny
      at2 = ci2*real(ffd(1,k))
      at1 = at2*aimag(ffd(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = 0.
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      wp = wp + at1*cu(1,1,k)**2
  110 continue
c mode numbers ky = 0, ny
      do 120 j = 2, nx
      at2 = ci2*real(ffd(j,1))
      at1 = at2*aimag(ffd(j,1))
      bxy(1,j,1) = 0.
      bxy(2,j,1) = at2*cu(2,j,1)
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      wp = wp + at1*cu(2,j,1)**2
  120 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny
  130 do 150 k = 2, ny
      do 140 j = 2, nx
      at1 = aimag(ffd(j,k))
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(2,j,k) = at1*cu(2,j,k)
  140 continue
c mode numbers kx = 0, nx
      at1 = aimag(ffd(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = 0.
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
  150 continue
c mode numbers ky = 0, ny
      do 160 j = 2, nx
      at1 = aimag(ffd(j,1))
      bxy(1,j,1) = 0.
      bxy(2,j,1) = at1*cu(2,j,1)
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
  160 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      return
      end
      subroutine IBPOISDX23(cu,bxy,ffd,ci,wm,nx,ny,nxv,ny2d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c input: cu,ffd,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 20*nx*ny + 11*(nx + ny)
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*|cu(kx,ky)*s(kx,ky)|**2),
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp
      complex cu, bxy, ffd, zero
      dimension cu(3,nxv,ny2d), bxy(3,nxv,ny2d)
      dimension ffd(nxv,ny)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      zero = cmplx(0.,0.)
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      at1 = ci2*real(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffd(j,k))
      bxy(1,j,k) = cmplx(0.,at3*real(cu(3,j,k)))
      bxy(2,j,k) = cmplx(0.,-at2*real(cu(3,j,k)))
      bxy(3,j,k) = cmplx(at3*aimag(cu(1,j,k))-at2*aimag(cu(2,j,k)),0.)
      bxy(1,j,k1) = bxy(1,j,k)
      bxy(2,j,k1) = -bxy(2,j,k)
      bxy(3,j,k1) = bxy(3,j,k)
      wp = wp + at1*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2 + real(cu
     1(3,j,k))**2)
   10 continue
   20 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))
      at3 = dky*at1
      at1 = at1*aimag(ffd(1,k))
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = cmplx(at3*aimag(cu(1,1,k)),0.)      
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(aimag(cu(1,1,k))**2 + real(cu(3,1,k))**2)
   30 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx
      at1 = ci2*real(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffd(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = zero
      bxy(3,j,1) = cmplx(-at2*aimag(cu(2,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(2,j,1))**2 + real(cu(3,j,1))**2)
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = float(nx*ny)*wp
      return
      end
      subroutine IBPOISD23(cu,bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with dirichlet boundary conditions (zero potential).
c fourier coefficients are constructed to perform the appropriate
c sin-sin, sin-cos, or cos-sin transforms
c input: cu,ffd,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 20*nx*ny + 11*(nx + ny)
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
      double precision wp
      complex ffd
      dimension cu(3,nxe,nye), bxy(3,nxe,nye)
      dimension ffd(nxv,ny)
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      at1 = ci2*real(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffd(j,k))
      bxy(1,j,k) = at3*cu(3,j,k)
      bxy(2,j,k) = -at2*cu(3,j,k)
      bxy(3,j,k) = at2*cu(2,j,k) - at3*cu(1,j,k)
      wp = wp + at1*(cu(1,j,k)**2 + cu(2,j,k)**2 + cu(3,j,k)**2)
   10 continue
   20 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))
      at3 = dky*at1
      at1 = at1*aimag(ffd(1,k))
      bxy(1,1,k) = 0.
      bxy(2,1,k) = 0.
      bxy(3,1,k) = -at3*cu(1,1,k)
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      bxy(3,nx+1,k) = 0.
      wp = wp + at1*(cu(1,1,k)**2 + cu(3,1,k)**2)
   30 continue
c mode numbers ky = 0, ny
      do 40 j = 2, nx
      at1 = ci2*real(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffd(j,1))
      bxy(1,j,1) = 0.
      bxy(2,j,1) = 0.
      bxy(3,j,1) = at2*cu(2,j,1)
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      bxy(3,j,ny+1) = 0.
      wp = wp + at1*(cu(2,j,1)**2 + cu(3,j,1)**2)
   40 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(3,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(3,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(3,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      bxy(3,nx+1,ny+1) = 0.
      wm = float(nx*ny)*wp
      return
      end
      subroutine MAXWELDX2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxv,ny2d)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with dirichlet boundary
c conditions (zero potential).
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 60*nxc*nyc + 36*(nxc + nyc)
c where nxc = nx - 1, nyc = ny - 1
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
c where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffd(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffd(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      double precision wp, ws
      complex exy, bxy, cu, ffd
      complex zero
      dimension exy(3,nxv,ny2d), bxy(3,nxv,ny2d), cu(3,nxv,ny2d)
      dimension ffd(nxv,ny)
      if (ci.le.0.) return
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffd(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffd(j,k))
      at7 = aimag(exy(1,j,k))
      at8 = aimag(exy(2,j,k))
      at9 = real(exy(3,j,k))
c update magnetic field half time step, ky > 0
      at4 = aimag(bxy(1,j,k)) - dth*(dky*at9)
      at5 = aimag(bxy(2,j,k)) + dth*(dkx*at9)
      at6 = real(bxy(3,j,k)) + dth*(dkx*at8 - dky*at7)
c update electric field whole time step
      at7 = at7 + cdt*(dky*at6) - afdt*aimag(cu(1,j,k))
      at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(2,j,k))
      at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*real(cu(3,j,k))
c update magnetic field half time step and store electric field
      at4 = at4 - dth*(dky*at9)
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 + dth*(dkx*at8 - dky*at7)
      ws = ws + anorm*(at7*at7 + at8*at8 + at9*at9)
      wp = wp + anorm*(at4*at4 + at5*at5 + at6*at6)
      exy(1,j,k) = cmplx(0.,at7)
      exy(2,j,k) = cmplx(0.,at8)
      exy(3,j,k) = cmplx(at9,0.)
      bxy(1,j,k) = cmplx(0.,at4)
      bxy(2,j,k) = cmplx(0.,at5)
      bxy(3,j,k) = cmplx(at6,0.)
c update electric and magnetic fields, ky < 0
      exy(1,j,k1) = cmplx(0.,-at7)
      exy(2,j,k1) = cmplx(0.,at8)
      exy(3,j,k1) = cmplx(-at9,0.)
      bxy(1,j,k1) = cmplx(0.,at4)
      bxy(2,j,k1) = cmplx(0.,-at5)
      bxy(3,j,k1) = cmplx(at6,0.)
   10 continue
   20 continue
      ws = ws + ws
      wp = wp + wp
c mode numbers kx = 0, nx/2
      do 30 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffd(1,k))
      at7 = aimag(exy(1,1,k))
      at9 = real(exy(3,1,k))
c update magnetic field half time step, ky > 0
      at4 = aimag(bxy(1,1,k)) - dth*(dky*at9)
      at6 = real(bxy(3,1,k)) - dth*(dky*at7)
c update electric field whole time step
      at7 = at7 + cdt*(dky*at6) - afdt*aimag(cu(1,1,k))
      at9 = at9 + cdt*(dky*at4) - afdt*real(cu(3,1,k))
c update magnetic field half time step and store electric field
      at4 = at4 - dth*(dky*at9)
      at6 = at6 - dth*(dky*at7)
      ws = ws + anorm*(at7*at7 + at9*at9)
      wp = wp + anorm*(at4*at4 + at6*at6)
      exy(1,1,k) = cmplx(0.,at7)
      exy(2,1,k) = zero
      exy(3,1,k) = cmplx(at9,0.)
      bxy(1,1,k) = cmplx(0.,at4)
      bxy(2,1,k) = zero
      bxy(3,1,k) = cmplx(at6,0.)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffd(j,1))
      at8 = aimag(exy(2,j,1))
      at9 = real(exy(3,j,1))
c update magnetic field half time step, ky > 0
      at5 = aimag(bxy(2,j,1)) + dth*(dkx*at9)
      at6 = real(bxy(3,j,1)) + dth*(dkx*at8)
c update electric field whole time step
      at8 = at8 - cdt*(dkx*at6) - afdt*aimag(cu(2,j,1))
      at9 = at9 - cdt*(dkx*at5) - afdt*real(cu(3,j,1))
c update magnetic field half time step and store electric field
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 + dth*(dkx*at8)
      ws = ws + anorm*(at8*at8 + at9*at9)
      wp = wp + anorm*(at5*at5 + at6*at6)
      exy(1,j,1) = zero
      exy(2,j,1) = cmplx(0.,at8)
      exy(3,j,1) = cmplx(at9,0.)
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,at5)
      bxy(3,j,1) = cmplx(at6,0.)
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
      return
      end
      subroutine MAXWELD2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxe,nye,nxv)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with  dirichlet boundary
c conditions (zero potential).
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 57*nxc*nyc + 36*(nxc + nyc)
c where nxc = nx - 1, nyc = ny - 1
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
c where kx = pi*j/nx, ky = pi*k/ny, c2 = 1./(ci*ci)
c and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffd(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffd(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
      double precision wp, ws
      complex ffd
      dimension exy(3,nxe,nye), bxy(3,nxe,nye), cu(3,nxe,nye)
      dimension ffd(nxv,ny)
      if (ci.le.0.) return
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffd(1,1))
      adt = affp*dt
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffd(j,k))
      at7 = exy(1,j,k)
      at8 = exy(2,j,k)
      at9 = exy(3,j,k)
c update magnetic field half time step, ky > 0
      at4 = bxy(1,j,k) - dth*(dky*at9)
      at5 = bxy(2,j,k) + dth*(dkx*at9)
      at6 = bxy(3,j,k) + dth*(dky*at7 - dkx*at8)
c update electric field whole time step
      at7 = at7 - cdt*(dky*at6) - afdt*cu(1,j,k)
      at8 = at8 + cdt*(dkx*at6) - afdt*cu(2,j,k)
      at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      at4 = at4 - dth*(dky*at9)
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 + dth*(dky*at7 - dkx*at8)
      ws = ws + anorm*(at7*at7 + at8*at8 + at9*at9)
      wp = wp + anorm*(at4*at4 + at5*at5 + at6*at6)
      bxy(1,j,k) = at4
      bxy(2,j,k) = at5
      bxy(3,j,k) = at6
      exy(1,j,k) = at7
      exy(2,j,k) = at8
      exy(3,j,k) = at9
   10 continue
   20 continue
      ws = ws + ws
      wp = wp + wp
c mode numbers kx = 0, nx/2
      do 30 k = 2, ny
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffd(1,k))
      at7 = exy(1,1,k)
      at9 = exy(3,1,k)
c update magnetic field half time step, ky > 0
      at4 = bxy(1,1,k) - dth*(dky*at9)
      at6 = bxy(3,1,k) + dth*(dky*at7)
c update electric field whole time step
      at7 = at7 - cdt*(dky*at6) - afdt*cu(1,1,k)
      at9 = at9 + cdt*(dky*at4) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      at4 = at4 - dth*(dky*at9)
      at6 = at6 + dth*(dky*at7)
      ws = ws + anorm*(at7*at7 + at9*at9)
      wp = wp + anorm*(at4*at4 + at6*at6)
      exy(1,1,k) = at7
      exy(2,1,k) = 0.
      exy(3,1,k) = at9
      exy(1,nx+1,k) = 0.
      exy(2,nx+1,k) = 0.
      exy(3,nx+1,k) = 0.
      bxy(1,1,k) = at4
      bxy(2,1,k) = 0.
      bxy(3,1,k) = at6
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      bxy(3,nx+1,k) = 0.
   30 continue
c mode numbers ky = 0, ny
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffd(j,1))
      at8 = exy(2,j,1)
      at9 = exy(3,j,1)
c update magnetic field half time step, ky > 0
      at5 = bxy(2,j,1) + dth*(dkx*at9)
      at6 = bxy(3,j,1) - dth*(dkx*at8)
c update electric field whole time step
      at8 = at8 + cdt*(dkx*at6) - afdt*cu(2,j,1)
      at9 = at9 - cdt*(dkx*at5) - afdt*cu(3,j,1)
c update magnetic field half time step and store electric field
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 - dth*(dkx*at8)
      ws = ws + anorm*(at8*at8 + at9*at9)
      wp = wp + anorm*(at5*at5 + at6*at6)
      bxy(1,j,1) = 0.
      bxy(2,j,1) = at5
      bxy(3,j,1) = at6
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      bxy(3,j,ny+1) = 0.
      exy(1,j,1) = 0.
      exy(2,j,1) = at8
      exy(3,j,1) = at9
      exy(1,j,ny+1) = 0.
      exy(2,j,ny+1) = 0.
      exy(3,j,ny+1) = 0.
   40 continue
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(3,1,1) = 0.
      bxy(1,nx+1,1) = 0.
      bxy(2,nx+1,1) = 0.
      bxy(3,nx+1,1) = 0.
      bxy(1,1,ny+1) = 0.
      bxy(2,1,ny+1) = 0.
      bxy(3,1,ny+1) = 0.
      bxy(1,nx+1,ny+1) = 0.
      bxy(2,nx+1,ny+1) = 0.
      bxy(3,nx+1,ny+1) = 0.
      exy(1,1,1) = 0.
      exy(2,1,1) = 0.
      exy(3,1,1) = 0.
      exy(1,nx+1,1) = 0.
      exy(2,nx+1,1) = 0.
      exy(3,nx+1,1) = 0.
      exy(1,1,ny+1) = 0.
      exy(2,1,ny+1) = 0.
      exy(3,1,ny+1) = 0.
      exy(1,nx+1,ny+1) = 0.
      exy(2,nx+1,ny+1) = 0.
      exy(3,nx+1,ny+1) = 0.
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
      subroutine DMFIELDD2(q2,q,nx,ny,nxv,ny2d,nxe,nye)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast sine transform in x and y
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex q2
      real q
      dimension q2(nxv,ny2d), q(nxe,nye)
      integer j, k, k1, ny2
      ny2 = 2*ny + 2
      do 20 k = 2, ny
      k1 = ny2 - k
      do 10 j = 1, nx
      q(j,k) = -real(q2(j,k))
   10 continue
      q(nx+1,k) = -real(q2(1,k1))
   20 continue
      do 30 j = 1, nx
      q(j,1) = -real(q2(j,1))
      q(j,ny+1) = -real(q2(j,ny+1))
   30 continue
      q(nx+1,1) = aimag(q2(1,1))
      q(nx+1,ny+1) = aimag(q2(1,ny+1))
      return
      end
      subroutine FMFIELDD2(q2,q,nx,ny,nxv,ny2d,nxe,nye)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast cosine transform in x and y
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex q2
      real q
      dimension q2(nxv,ny2d), q(nxe,nye)
      integer j, k, k1, ny2
      ny2 = 2*ny + 2
      do 20 k = 2, ny
      k1 = ny2 - k
      do 10 j = 1, nx
      q(j,k) = real(q2(j,k))
   10 continue
      q(nx+1,k) = real(q2(1,k1))
   20 continue
      do 30 j = 1, nx
      q(j,1) = real(q2(j,1))
      q(j,ny+1) = real(q2(j,ny+1))
   30 continue
      q(nx+1,1) = aimag(q2(1,1))
      q(nx+1,ny+1) = aimag(q2(1,ny+1))
      return
      end
      subroutine CMFIELDD2(cu2,cu,nx,ny,nxv,ny2d,nxe,nye)
c this subroutine copies the current into a smaller array
c which would have been created by fast sine/cosine transforms in x
c and y
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex cu2
      real cu
      dimension cu2(3,nxv,ny2d), cu(3,nxe,nye)
      integer j, k, k1, ny2
      ny2 = 2*ny + 2
      do 20 k = 2, ny
      k1 = ny2 - k
      do 10 j = 1, nx
      cu(1,j,k) = -aimag(cu2(1,j,k))
      cu(2,j,k) = -aimag(cu2(2,j,k))
      cu(3,j,k) = -real(cu2(3,j,k))
   10 continue
      cu(1,nx+1,k) = aimag(cu2(1,1,k1))
      cu(2,nx+1,k) = aimag(cu2(2,1,k1))
      cu(3,nx+1,k) = -real(cu2(3,1,k1))
   20 continue
      do 30 j = 2, nx
      cu(1,j,1) = -aimag(cu2(1,j,1))
      cu(2,j,1) = -aimag(cu2(2,j,1))
      cu(3,j,1) = -real(cu2(3,j,1))
      cu(1,j,ny+1) = -aimag(cu2(1,j,ny+1))
      cu(2,j,ny+1) = -aimag(cu2(2,j,ny+1))
      cu(3,j,ny+1) = -real(cu2(3,j,ny+1))
   30 continue
      cu(1,1,1) = -real(cu2(1,1,1))
      cu(2,1,1) = -real(cu2(2,1,1))
      cu(3,1,1) = -real(cu2(3,1,1))
      cu(1,nx+1,1) = aimag(cu2(1,1,1))
      cu(2,nx+1,1) = aimag(cu2(2,1,1))
      cu(3,nx+1,1) = aimag(cu2(3,1,1))
      cu(1,1,ny+1) = -real(cu2(1,1,ny+1))
      cu(2,1,ny+1) = -real(cu2(2,1,ny+1))
      cu(3,1,ny+1) = -real(cu2(3,1,ny+1))
      cu(1,nx+1,ny+1) = aimag(cu2(1,1,ny+1))
      cu(2,nx+1,ny+1) = aimag(cu2(2,1,ny+1))
      cu(3,nx+1,ny+1) = aimag(cu2(3,1,ny+1))
      return
      end
      subroutine EMFIELDD2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, nxv, ny2d, nxe, nye
      complex fxy, ffd
      real exy
      dimension fxy(3,nxv,ny2d), exy(3,nxe,nye)
      dimension ffd(nxv,ny)
      complex zero
      integer j, k, k1, ny2
      real at1
      ny2 = 2*ny + 2
      zero = cmplx(0.,0.)
c add the fields
      if (isign.gt.0) then
         do 20 k = 2, ny
         k1 = ny2 - k
         do 10 j = 2, nx
         at1 = aimag(ffd(j,k))
         fxy(1,j,k) = fxy(1,j,k) + cmplx(0.,-exy(1,j,k)*at1)
         fxy(2,j,k) = fxy(2,j,k) + cmplx(0.,-exy(2,j,k)*at1)
         fxy(3,j,k) = fxy(3,j,k) + cmplx(-exy(3,j,k)*at1,0.)
         fxy(1,j,k1) = fxy(1,j,k1) + cmplx(0.,exy(1,j,k)*at1)
         fxy(2,j,k1) = fxy(2,j,k1) + cmplx(0.,-exy(2,j,k)*at1)
         fxy(3,j,k1) = fxy(3,j,k1) + cmplx(exy(3,j,k)*at1,0.)
   10    continue
         at1 = aimag(ffd(1,k))
         fxy(1,1,k) = fxy(1,1,k) + cmplx(0.,-exy(1,1,k)*at1)
         fxy(2,1,k) = fxy(2,1,k) + cmplx(0.,-exy(2,1,k)*at1)
         fxy(3,1,k) = fxy(3,1,k) + cmplx(-exy(3,1,k)*at1,0.)
   20    continue
         do 30 j = 1, nx
         at1 = aimag(ffd(j,1))
         fxy(1,j,1) = fxy(1,j,1) + cmplx(0.,-exy(1,j,1)*at1)
         fxy(2,j,1) = fxy(2,j,1) + cmplx(0.,-exy(2,j,1)*at1)
         fxy(3,j,1) = fxy(3,j,1) + cmplx(-exy(3,j,1)*at1,0.)
   30    continue
c copy the magnetic fields
      else if (isign.lt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         do 40 j = 2, nx
         at1 = aimag(ffd(j,k))
         fxy(1,j,k) = cmplx(0.,-exy(1,j,k)*at1)
         fxy(2,j,k) = cmplx(0.,-exy(2,j,k)*at1)
         fxy(3,j,k) = cmplx(exy(3,j,k)*at1,0.)
         fxy(1,j,k1) = cmplx(0.,-exy(1,j,k)*at1)
         fxy(2,j,k1) = cmplx(0.,exy(2,j,k)*at1)
         fxy(3,j,k1) = cmplx(exy(3,j,k)*at1,0.)
   40    continue
         at1 = aimag(ffd(1,k))
         fxy(1,1,k) = cmplx(0.,-exy(1,1,k)*at1)
         fxy(2,1,k) = cmplx(0.,-exy(2,1,k)*at1)
         fxy(3,1,k) = cmplx(exy(3,1,k)*at1,0.)
         fxy(1,1,k1) = zero
         fxy(2,1,k1) = zero
         fxy(3,1,k1) = zero
   50    continue
         k1 = ny + 1
         do 60 j = 1, nx
         at1 = aimag(ffd(j,1))
         fxy(1,j,1) = cmplx(0.,-exy(1,j,1)*at1)
         fxy(2,j,1) = cmplx(0.,-exy(2,j,1)*at1)
         fxy(3,j,1) = cmplx(exy(3,j,1)*at1,0.)
         fxy(1,j,k1) = zero
         fxy(2,j,k1) = zero
         fxy(3,j,k1) = zero
   60    continue
c copy the electric fields
      else
         do 80 k = 2, ny
         k1 = ny2 - k
         do 70 j = 2, nx
         fxy(1,j,k) = cmplx(0.,-exy(1,j,k))
         fxy(2,j,k) = cmplx(0.,-exy(2,j,k))
         fxy(3,j,k) = cmplx(-exy(3,j,k),0.)
         fxy(1,j,k1) = cmplx(0.,exy(1,j,k))
         fxy(2,j,k1) = cmplx(0.,-exy(2,j,k))
         fxy(3,j,k1) = cmplx(exy(3,j,k),0.)
   70    continue
         fxy(1,1,k) = cmplx(0.,-exy(1,1,k))
         fxy(2,1,k) = cmplx(0.,-exy(2,1,k))
         fxy(3,1,k) = cmplx(-exy(3,1,k),0.)
         fxy(1,1,k1) = cmplx(0.,exy(1,nx+1,k))
         fxy(2,1,k1) = cmplx(0.,-exy(2,nx+1,k))
         fxy(3,1,k1) = cmplx(exy(3,nx+1,k),0.)
   80    continue
         k1 = ny + 1
         do 90 j = 2, nx
         fxy(1,j,1) = cmplx(0.,-exy(1,j,1))
         fxy(2,j,1) = cmplx(0.,-exy(2,j,1))
         fxy(3,j,1) = cmplx(-exy(3,j,1),0.)
         fxy(1,j,k1) = cmplx(0.,-exy(1,j,k1))
         fxy(2,j,k1) = cmplx(0.,-exy(2,j,k1))
         fxy(3,j,k1) = cmplx(-exy(3,j,k1),0.)
   90    continue
         fxy(1,1,1) = cmplx(-exy(1,1,1),exy(1,nx+1,1))
         fxy(2,1,1) = cmplx(-exy(2,1,1),exy(2,nx+1,1))
         fxy(3,1,1) = cmplx(-exy(3,1,1),exy(3,nx+1,1))
         fxy(1,1,k1) = cmplx(-exy(1,1,k1),exy(1,nx+1,k1))
         fxy(2,1,k1) = cmplx(-exy(2,1,k1),exy(2,nx+1,k1))
         fxy(3,1,k1) = cmplx(-exy(3,1,k1),exy(3,nx+1,k1))
      endif
      return
      end
      subroutine PMFIELDD2(pot2,pot,nx,ny,nxv,ny2d,nxe,nye)
c copies image charges appropriate for potential
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex pot2
      real pot
      dimension pot2(nxv,ny2d), pot(nxe,nye)
      complex zero
      integer j, k, k1, ny2
      ny2 = 2*ny + 2
      zero = cmplx(0.,0.)
      do 20 k = 2, ny
      k1 = ny2 - k
      do 10 j = 2, nx
      pot2(j,k) = cmplx(-pot(j,k),0.)
      pot2(j,k1) = cmplx(pot(j,k),0.)
   10 continue
      pot2(1,k) = cmplx(-pot(1,k),0.)
      pot2(1,k1) = cmplx(-pot(nx+1,k),0.)
   20 continue
      k1 = ny + 1
      do 30 j = 2, nx
      pot2(j,1) = cmplx(-pot(j,1),0.)
      pot2(j,k1) = cmplx(-pot(j,k1),0.)
   30 continue
      pot2(1,1) = cmplx(-pot(1,1),pot(nx+1,1))
      pot2(1,k1) = cmplx(-pot(1,k1),pot(nx+1,k1))
      return
      end
      subroutine BMFIELDD2(fxy,bxy,nx,ny,nxv,ny2d,nxe,nye)
c copies image charges appropriate for magnetic field
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex fxy
      real bxy
      dimension fxy(3,nxv,ny2d), bxy(3,nxe,nye)
      integer j, k, k1, ny2
      ny2 = 2*ny + 2
c copy the magnetic fields
      do 20 k = 2, ny
      k1 = ny2 - k
      do 10 j = 2, nx
      fxy(1,j,k) = cmplx(0.,-bxy(1,j,k))
      fxy(2,j,k) = cmplx(0.,-bxy(2,j,k))
      fxy(3,j,k) = cmplx(bxy(3,j,k),0.)
      fxy(1,j,k1) = cmplx(0.,-bxy(1,j,k))
      fxy(2,j,k1) = cmplx(0.,bxy(2,j,k))
      fxy(3,j,k1) = cmplx(bxy(3,j,k),0.)
   10 continue
      fxy(1,1,k) = cmplx(0.,-bxy(1,1,k))
      fxy(2,1,k) = cmplx(0.,-bxy(2,1,k))
      fxy(3,1,k) = cmplx(bxy(3,1,k),0.)
      fxy(1,1,k1) = cmplx(0.,-bxy(1,nx+1,k))
      fxy(2,1,k1) = cmplx(0.,bxy(2,nx+1,k))
      fxy(3,1,k1) = cmplx(bxy(3,nx+1,k),0.)
   20 continue
      k1 = ny + 1
      do 30 j = 2, nx
      fxy(1,j,1) = cmplx(0.,-bxy(1,j,1))
      fxy(2,j,1) = cmplx(0.,-bxy(2,j,1))
      fxy(3,j,1) = cmplx(bxy(3,j,1),0.)
      fxy(1,j,k1) = cmplx(0.,-bxy(1,j,k1))
      fxy(2,j,k1) = cmplx(0.,-bxy(2,j,k1))
      fxy(3,j,k1) = cmplx(bxy(3,j,k1),0.)
   30 continue
      fxy(1,1,1) = cmplx(-bxy(1,1,1),bxy(1,nx+1,1))
      fxy(2,1,1) = cmplx(-bxy(2,1,1),bxy(2,nx+1,1))
      fxy(3,1,1) = cmplx(bxy(3,1,1),bxy(3,nx+1,1))
      fxy(1,1,k1) = cmplx(-bxy(1,1,k1),bxy(1,nx+1,k1))
      fxy(2,1,k1) = cmplx(-bxy(2,1,k1),bxy(2,nx+1,k1))
      fxy(3,1,k1) = cmplx(bxy(3,1,k1),bxy(3,nx+1,k1))
      return
      end
      subroutine CPFIELDD2(fxy,exy,nx,ny,nxv,ny2d,nxe,nye)
c this subroutine copies complex vector fields and adds image charges
c appropriate for electric field
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex fxy
      real exy
      dimension fxy(3,nxv,ny2d), exy(3,nxe,nye)
c local data
      integer isign
      complex ffd
      dimension ffd(1,1)
      isign = 0
      call EMFIELDD2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
      return
      end
      subroutine AVPOTDX23(bxy,axy,nx,ny,nxv,ny2d)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential).
c input: bxy, nx, ny, nxv, ny2d, output: axy
c approximate flop count is: 14*nxc*nyc + 4*(nxc + nyc)
c and nxc*nyc divides, where nxc = nx - 1, nyc = ny - 1
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
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      complex bxy, axy, zero
      dimension bxy(3,nxv,ny2d), axy(3,nxv,ny2d)
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c calculate vector potential
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      at4 = real(bxy(3,j,k))
      at5 = aimag(bxy(2,j,k))
      at6 = aimag(bxy(1,j,k))
      at6 = at3*at6 - at2*at5
      at5 = -at2*at4
      at4 = at3*at4
      axy(1,j,k) = cmplx(0.,at4)
      axy(2,j,k) = cmplx(0.,at5)
      axy(3,j,k) = cmplx(at6,0.)
      axy(1,j,k1) = cmplx(0.,-at4)
      axy(2,j,k1) = cmplx(0.,at5)
      axy(3,j,k1) = cmplx(-at6,0.)
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      at4 = real(bxy(3,1,k))
      at6 = aimag(bxy(1,1,k))
      axy(1,1,k) = cmplx(0.,at3*at4)
      axy(2,1,k) = zero
      axy(3,1,k) = cmplx(at3*at6,0.)
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = ny + 1
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      at2 = 1.0/dkx
      at4 = real(bxy(3,j,1))
      at5 = aimag(bxy(2,j,1))
      axy(1,j,1) = zero
      axy(2,j,1) = cmplx(0.,-at2*at4)
      axy(3,j,1) = cmplx(-at2*at5,0.)
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
      subroutine AVPOTD23(bxy,axy,nx,ny,nxe,nye)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential).
c input: bxy, nx, ny, nxe, nye, output: axy
c approximate flop count is: 12*nxc*nyc + 4*(nxc + nyc)
c and nxc*nyc divides, where nxc = nx - 1, nyc = ny - 1
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
c nxe = first dimension of field arrays, must be >= nx+1
c nye = second dimension of field arrays, must be >= ny+1
      dimension bxy(3,nxe,nye), axy(3,nxe,nye)
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
c calculate vector potential
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      at4 = bxy(3,j,k)
      at5 = bxy(2,j,k)
      at6 = bxy(1,j,k)
      axy(1,j,k) = -at3*at4
      axy(2,j,k) = at2*at4
      axy(3,j,k) = at3*at6 - at2*at5
   10 continue
   20 continue
c mode numbers kx = 0, nx/2
cdir$ ivdep
      do 30 k = 2, ny
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      at4 = bxy(3,1,k)
      at6 = bxy(1,1,k)
      axy(1,1,k) = -at3*at4
      axy(2,1,k) = 0.
      axy(3,1,k) = at3*at6
      axy(1,nx+1,k) = 0.
      axy(2,nx+1,k) = 0.
      axy(3,nx+1,k) = 0.
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      at2 = 1.0/dkx
      at4 = bxy(3,j,1)
      at5 = bxy(2,j,1)
      axy(1,j,1) = 0.
      axy(2,j,1) = at2*at4
      axy(3,j,1) = -at2*at5
      axy(1,j,ny+1) = 0.
      axy(2,j,ny+1) = 0.
      axy(3,j,ny+1) = 0.
   40 continue
      axy(1,1,1) = 0.
      axy(2,1,1) = 0.
      axy(3,1,1) = 0.
      axy(1,nx+1,1) = 0.
      axy(2,nx+1,1) = 0.
      axy(3,nx+1,1) = 0.
      axy(1,1,ny+1) = 0.
      axy(2,1,ny+1) = 0.
      axy(3,1,ny+1) = 0.
      axy(1,nx+1,ny+1) = 0.
      axy(2,nx+1,ny+1) = 0.
      axy(3,nx+1,ny+1) = 0.
      return
      end
      subroutine DPOYNTDX2(q,cu,ffd,ci,sx,sy,sz,nx,ny,nxv,ny2d)
c this subroutine calculates the momentum in the darwin field given by
c the poynting flux in 2-1/2d.
c with dirichlet boundary conditions (zero potential).
c inputs are the charge density and current density.
c outputs are sx, sy, sz
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
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c sx/sy/sz = x/y/z components of field momentum
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      implicit none
      integer nx, ny, nxv, ny2d
      real ci, sx, sy, sz
      complex q, cu, ffd
      dimension q(nxv,ny2d), cu(3,nxv,ny2d)
      dimension ffd(nxv,ny)
      integer ny2, j, k
      real dnx, dny, ci2, affp, dky
      real at1, at2, at3, at4, at5
      double precision wx, wy, wz
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      ci2 = ci*ci
      affp = real(ffd(1,1))
c calculate force/charge and sum field energy
      wx = 0.0d0
      wy = 0.0d0
      wz = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      at1 = real(ffd(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      at1 = real(q(j,k))
      at4 = at2*at1
      at1 = at3*at1
      at5 = ci2*real(cu(3,j,k))
      at2 = at2*at5
      at3 = at3*at5
      wz = wz + (at4*at2 + at1*at3)
   10 continue
   20 continue
      wz = wz + wz
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      at1 = real(ffd(1,k))
      at3 = dny*float(k - 1)*at1
      at1 = real(q(j,k))
      at4 = 0.0
      at1 = at3*at1
      at3 = ci2*at3*real(cu(3,j,k))
      wz = wz + at1*at3
   30 continue
c mode numbers ky = 0, ny
      do 40 j = 2, nx
      at1 = real(ffd(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = real(q(j,k))
      at4 = at2*at1
      at2 = ci2*at2*real(cu(3,j,k))
      wz = wz + at4*at2
   40 continue
      at1 = 2.0*float(nx*ny)/affp
      sx = at1*wx
      sy = at1*wy
      sz = at1*wz
      return
      end
      subroutine LACFGUARD2(cus,cu,q2m0,nx,ny,nxe,nye)
c increment extended field with scaled field
c quadratic interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(3,nxe,nye), cu(3,nxe,nye)
      integer i, j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
c initialize extended field, with zero in the edges
      do 40 k = 1, ny1
      do 20 j = 1, nx1
      do 10 i = 1, 3
      cus(i,j+1,k+1) = cus(i,j+1,k+1) - q2m0*cu(i,j+1,k+1)
   10 continue
   20 continue
      do 30 i = 1, 3
      cus(i,1,k+1) = 0.
      cus(i,nx+3,k+1) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx1
      do 50 i = 1, 3
      cus(i,j+1,1) = 0.
      cus(i,j+1,ny+3) = 0.
   50 continue
   60 continue
      do 70 i = 1, 3
      cus(i,1,1) = 0.
      cus(i,nx+3,1) = 0.
      cus(i,1,ny+3) = 0.
      cus(i,nx+3,ny+3) = 0.
   70 continue
      return
      end
      subroutine LACFGUARD22(cus,cu,q2m0,nx,ny,nxe,nye)
c increment extended field with scaled field
c quadratic interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(2,nxe,nye), cu(2,nxe,nye)
      integer i, j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
c initialize extended field, with zero in the edges
      do 40 k = 1, ny1
      do 20 j = 1, nx1
      do 10 i = 1, 2
      cus(i,j+1,k+1) = cus(i,j+1,k+1) - q2m0*cu(i,j+1,k+1)
   10 continue
   20 continue
      do 30 i = 1, 2
      cus(i,1,k+1) = 0.
      cus(i,nx+3,k+1) = 0.
   30 continue
   40 continue
      do 60 j = 1, nx1
      do 50 i = 1, 2
      cus(i,j+1,1) = 0.
      cus(i,j+1,ny+3) = 0.
   50 continue
   60 continue
      do 70 i = 1, 2
      cus(i,1,1) = 0.
      cus(i,nx+3,1) = 0.
      cus(i,1,ny+3) = 0.
      cus(i,nx+3,ny+3) = 0.
   70 continue
      return
      end
      subroutine LSCFGUARD2L(cus,cu,q2m0,nx,ny,nxe,nye)
c initialize extended field with scaled field
c linear interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(3,nxe,nye), cu(3,nxe,nye)
      integer i, j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
c initialize extended field
      do 30 k = 1, ny1
      do 20 j = 1, nx1
      do 10 i = 1, 3
      cus(i,j,k) = -q2m0*cu(i,j,k)
   10 continue
   20 continue
   30 continue
      return
      end
      subroutine LSCFGUARD22L(cus,cu,q2m0,nx,ny,nxe,nye)
c initialize extended field with scaled field
c linear interpolation
      implicit none
      real cus, cu, q2m0
      integer nx, ny, nxe, nye
      dimension cus(2,nxe,nye), cu(2,nxe,nye)
      integer i, j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
c initialize extended field
      do 30 k = 1, ny1
      do 20 j = 1, nx1
      do 10 i = 1, 2
      cus(i,j,k) = -q2m0*cu(i,j,k)
   10 continue
   20 continue
   30 continue
      return
      end
      subroutine LSMCGUARD2(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ngx,ngy,nxe,
     1nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension amu(4,nxe,nye)
      integer i, j, k, nxg, nyg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx3 = nx + 3
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      amu(1,j+ngx+1,k+ngy+1) = x2y2m0
      amu(2,j+ngx+1,k+ngy+1) = xym0
      amu(3,j+ngx+1,k+ngy+1) = zxm0
      amu(4,j+ngx+1,k+ngy+1) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,1,k+ngy+1) = 0.
      amu(i,2,k+ngy+1) = 0.
      amu(i,nx+2,k+ngy+1) = 0.
      amu(i,nx+3,k+ngy+1) = 0.
   20 continue
      amu(1,ngx+2,k+ngy+1) = .5*x2y2m0
      amu(2,ngx+2,k+ngy+1) = .5*xym0
      amu(3,ngx+2,k+ngy+1) = .5*zxm0
      amu(4,ngx+2,k+ngy+1) = .5*zym0
      amu(1,nx-ngx+2,k+ngy+1) = .5*x2y2m0
      amu(2,nx-ngx+2,k+ngy+1) = .5*xym0
      amu(3,nx-ngx+2,k+ngy+1) = .5*zxm0
      amu(4,nx-ngx+2,k+ngy+1) = .5*zym0
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 4
      amu(i,j,1) = 0.
      amu(i,j,2) = 0.
      amu(i,j,ny+2) = 0.
      amu(i,j,ny+3) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      amu(1,j+ngx+1,ngy+2) = .5*x2y2m0
      amu(2,j+ngx+1,ngy+2) = .5*xym0
      amu(3,j+ngx+1,ngy+2) = .5*zxm0
      amu(4,j+ngx+1,ngy+2) = .5*zym0
      amu(1,j+ngx+1,ny-ngy+2) = .5*x2y2m0
      amu(2,j+ngx+1,ny-ngy+2) = .5*xym0
      amu(3,j+ngx+1,ny-ngy+2) = .5*zxm0
      amu(4,j+ngx+1,ny-ngy+2) = .5*zym0
   60 continue
      do 70 i = 1, 4
      amu(i,1,ngy+2) = 0.
      amu(i,2,ngy+2) = 0.
      amu(i,nx+2,ngy+2) = 0.
      amu(i,nx+3,ngy+2) = 0.
      amu(i,1,ny-ngy+2) = 0.
      amu(i,2,ny-ngy+2) = 0.
      amu(i,nx+2,ny-ngy+2) = 0.
      amu(i,nx+3,ny-ngy+2) = 0.
   70 continue
      amu(1,ngx+2,ngy+2) = .25*x2y2m0
      amu(2,ngx+2,ngy+2) = .25*xym0
      amu(3,ngx+2,ngy+2) = .25*zxm0
      amu(4,ngx+2,ngy+2) = .25*zym0
      amu(1,nx-ngx+2,ngy+2) = .25*x2y2m0
      amu(2,nx-ngx+2,ngy+2) = .25*xym0
      amu(3,nx-ngx+2,ngy+2) = .25*zxm0
      amu(4,nx-ngx+2,ngy+2) = .25*zym0
      amu(1,ngx+2,ny-ngy+2) = .25*x2y2m0
      amu(2,ngx+2,ny-ngy+2) = .25*xym0
      amu(3,ngx+2,ny-ngy+2) = .25*zxm0
      amu(4,ngx+2,ny-ngy+2) = .25*zym0
      amu(1,nx-ngx+2,ny-ngy+2) = .25*x2y2m0
      amu(2,nx-ngx+2,ny-ngy+2) = .25*xym0
      amu(3,nx-ngx+2,ny-ngy+2) = .25*zxm0
      amu(4,nx-ngx+2,ny-ngy+2) = .25*zym0
      return
      end
      subroutine LSMCGUARD22(amu,x2y2m0,xym0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real amu, x2y2m0, xym0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension amu(2,nxe,nye)
      integer i, j, k, nxg, nyg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx3 = nx + 3
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      amu(1,j+ngx+1,k+ngy+1) = x2y2m0
      amu(2,j+ngx+1,k+ngy+1) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,1,k+ngy+1) = 0.
      amu(i,2,k+ngy+1) = 0.
      amu(i,nx+2,k+ngy+1) = 0.
      amu(i,nx+3,k+ngy+1) = 0.
   20 continue
      amu(1,ngx+2,k+ngy+1) = .5*x2y2m0
      amu(2,ngx+2,k+ngy+1) = .5*xym0
      amu(1,nx-ngx+2,k+ngy+1) = .5*x2y2m0
      amu(2,nx-ngx+2,k+ngy+1) = .5*xym0
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 2
      amu(i,j,1) = 0.
      amu(i,j,2) = 0.
      amu(i,j,ny+2) = 0.
      amu(i,j,ny+3) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      amu(1,j+ngx+1,ngy+2) = .5*x2y2m0
      amu(2,j+ngx+1,ngy+2) = .5*xym0
      amu(1,j+ngx+1,ny-ngy+2) = .5*x2y2m0
      amu(2,j+ngx+1,ny-ngy+2) = .5*xym0
   60 continue
      do 70 i = 1, 2
      amu(i,1,ngy+2) = 0.
      amu(i,2,ngy+2) = 0.
      amu(i,nx+2,ngy+2) = 0.
      amu(i,nx+3,ngy+2) = 0.
      amu(i,1,ny-ngy+2) = 0.
      amu(i,2,ny-ngy+2) = 0.
      amu(i,nx+2,ny-ngy+2) = 0.
      amu(i,nx+3,ny-ngy+2) = 0.
   70 continue
      amu(1,ngx+2,ngy+2) = .25*x2y2m0
      amu(2,ngx+2,ngy+2) = .25*xym0
      amu(1,nx-ngx+2,ngy+2) = .25*x2y2m0
      amu(2,nx-ngx+2,ngy+2) = .25*xym0
      amu(1,ngx+2,ny-ngy+2) = .25*x2y2m0
      amu(2,ngx+2,ny-ngy+2) = .25*xym0
      amu(1,nx-ngx+2,ny-ngy+2) = .25*x2y2m0
      amu(2,nx-ngx+2,ny-ngy+2) = .25*xym0
      return
      end
      subroutine LAMCGUARD2(amu,nx,ny,nxe,nye,ndim)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = second dimension of input array amu, must be >= nx+3
c nxe = third dimension of input array amu, must be >= ny+3
      implicit none
      real amu
      integer nx, ny, nxe, nye, ndim
      dimension amu(ndim,nxe,nye)
c local data
      integer i, j, k, nx1, ny3
      nx1 = nx + 1
      ny3 = ny + 3
c add up guard cells
      do 20 k = 1, ny3
      do 10 i = 1, ndim
      amu(i,2,k) = amu(i,2,k) + 2.*amu(i,1,k)
      amu(i,3,k) = amu(i,3,k) - amu(i,1,k)
      amu(i,nx+1,k) = amu(i,nx+1,k) - amu(i,nx+3,k)
      amu(i,nx+2,k) = amu(i,nx+2,k) + 2.*amu(i,nx+3,k)
      amu(i,1,k) = 0.
      amu(i,nx+3,k) = 0.
   10 continue
   20 continue
      do 40 j = 1, nx1
      do 30 i = 1, ndim
      amu(i,j+1,2) = amu(i,j+1,2) + 2.*amu(i,j+1,1)
      amu(i,j+1,3) = amu(i,j+1,3) - amu(i,j+1,1)
      amu(i,j+1,ny+1) = amu(i,j+1,ny+1) - amu(i,j+1,ny+3)
      amu(i,j+1,ny+2) = amu(i,j+1,ny+2) + 2.*amu(i,j+1,ny+3)
      amu(i,j+1,1) = 0.
      amu(i,j+1,ny+3) = 0.
   30 continue
   40 continue
      return
      end
      subroutine LSMCGUARD2L(amu,x2y2m0,xym0,zxm0,zym0,nx,ny,ngx,ngy,nxe
     1,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real amu, x2y2m0, xym0, zxm0, zym0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension amu(4,nxe,nye)
      integer i, j, k, nxg, nyg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx1 = nx + 1
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      amu(1,j+ngx,k+ngy) = x2y2m0
      amu(2,j+ngx,k+ngy) = xym0
      amu(3,j+ngx,k+ngy) = zxm0
      amu(4,j+ngx,k+ngy) = zym0
   10 continue
      do 20 i = 1, 4
      amu(i,1,k+ngy) = 0.
      amu(i,nx+1,k+ngy) = 0.
   20 continue
      amu(1,ngx+1,k+ngy) = .5*x2y2m0
      amu(2,ngx+1,k+ngy) = .5*xym0
      amu(3,ngx+1,k+ngy) = .5*zxm0
      amu(4,ngx+1,k+ngy) = .5*zym0
      amu(1,nx-ngx+1,k+ngy) = .5*x2y2m0
      amu(2,nx-ngx+1,k+ngy) = .5*xym0
      amu(3,nx-ngx+1,k+ngy) = .5*zxm0
      amu(4,nx-ngx+1,k+ngy) = .5*zym0
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 4
      amu(i,j,1) = 0.
      amu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      amu(1,j+ngx,ngy+1) = .5*x2y2m0
      amu(2,j+ngx,ngy+1) = .5*xym0
      amu(3,j+ngx,ngy+1) = .5*zxm0
      amu(4,j+ngx,ngy+1) = .5*zym0
      amu(1,j+ngx,ny-ngy+1) = .5*x2y2m0
      amu(2,j+ngx,ny-ngy+1) = .5*xym0
      amu(3,j+ngx,ny-ngy+1) = .5*zxm0
      amu(4,j+ngx,ny-ngy+1) = .5*zym0
   60 continue
      do 70 i = 1, 4
      amu(i,1,ngy+1) = 0.
      amu(i,nx+1,ngy+1) = 0.
      amu(i,1,ny-ngy+1) = 0.
      amu(i,nx+1,ny-ngy+1) = 0.
   70 continue
      amu(1,ngx+1,ngy+1) = .25*x2y2m0
      amu(2,ngx+1,ngy+1) = .25*xym0
      amu(3,ngx+1,ngy+1) = .25*zxm0
      amu(4,ngx+1,ngy+1) = .25*zym0
      amu(1,nx-ngx+1,ngy+1) = .25*x2y2m0
      amu(2,nx-ngx+1,ngy+1) = .25*xym0
      amu(3,nx-ngx+1,ngy+1) = .25*zxm0
      amu(4,nx-ngx+1,ngy+1) = .25*zym0
      amu(1,ngx+1,ny-ngy+1) = .25*x2y2m0
      amu(2,ngx+1,ny-ngy+1) = .25*xym0
      amu(3,ngx+1,ny-ngy+1) = .25*zxm0
      amu(4,ngx+1,ny-ngy+1) = .25*zym0
      amu(1,nx-ngx+1,ny-ngy+1) = .25*x2y2m0
      amu(2,nx-ngx+1,ny-ngy+1) = .25*xym0
      amu(3,nx-ngx+1,ny-ngy+1) = .25*zxm0
      amu(4,nx-ngx+1,ny-ngy+1) = .25*zym0
      return
      end
      subroutine LSMCGUARD22L(amu,x2y2m0,xym0,nx,ny,ngx,ngy,nxe,nye)
c initialize extended non-periodic field
c ngx/ngy = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real amu, x2y2m0, xym0
      integer nx, ny, ngx, ngy, nxe, nye
      dimension amu(2,nxe,nye)
      integer i, j, k, nxg, nyg, nx1
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nyg = ny - 2*ngy
      nx1 = nx + 1
      do 30 k = 2, nyg
      do 10 j = 2, nxg
      amu(1,j+ngx,k+ngy) = x2y2m0
      amu(2,j+ngx,k+ngy) = xym0
   10 continue
      do 20 i = 1, 2
      amu(i,1,k+ngy) = 0.
      amu(i,nx+1,k+ngy) = 0.
   20 continue
      amu(1,ngx+1,k+ngy) = .5*x2y2m0
      amu(2,ngx+1,k+ngy) = .5*xym0
      amu(1,nx-ngx+1,k+ngy) = .5*x2y2m0
      amu(2,nx-ngx+1,k+ngy) = .5*xym0
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 2
      amu(i,j,1) = 0.
      amu(i,j,ny+1) = 0.
   40 continue
   50 continue
      do 60 j = 2, nxg
      amu(1,j+ngx,ngy+1) = .5*x2y2m0
      amu(2,j+ngx,ngy+1) = .5*xym0
      amu(1,j+ngx,ny-ngy+1) = .5*x2y2m0
      amu(2,j+ngx,ny-ngy+1) = .5*xym0
   60 continue
      do 70 i = 1, 2
      amu(i,1,ngy+1) = 0.
      amu(i,nx+1,ngy+1) = 0.
      amu(i,1,ny-ngy+1) = 0.
      amu(i,nx+1,ny-ngy+1) = 0.
   70 continue
      amu(1,ngx+1,ngy+1) = .25*x2y2m0
      amu(2,ngx+1,ngy+1) = .25*xym0
      amu(1,nx-ngx+1,ngy+1) = .25*x2y2m0
      amu(2,nx-ngx+1,ngy+1) = .25*xym0
      amu(1,ngx+1,ny-ngy+1) = .25*x2y2m0
      amu(2,ngx+1,ny-ngy+1) = .25*xym0
      amu(1,nx-ngx+1,ny-ngy+1) = .25*x2y2m0
      amu(2,nx-ngx+1,ny-ngy+1) = .25*xym0
      return
      end
      subroutine DBLSIN2M(amu,amu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a doubled tensor array amu2 from a tensor array
c amu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the xx-yy component is an odd function in
c both x and y, the xy component is an even function in both x and y, the
c zx component is an even function in x and an odd function in y, the zy
c component is an odd function in x and and even function in y.
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array amu, must be >= nx+1
c nyv = third dimension of input array amu, must be >= ny+1
c nx2v = second dimension of output array amu2, must be >= 2*nx
c ny2 = third dimension of output array amu2, must be >= 2*ny
      implicit none
      real amu, amu2
      integer nx, ny, nxv, nyv, nx2v, ny2
      dimension amu(4,nxv,nyv), amu2(4,nx2v,ny2)
c local data
      integer j, k, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 20 k = 1, nys
      do 10 j = 1, nxs
      amu2(1,j+1,k+1) = amu(1,j+1,k+1)
      amu2(2,j+1,k+1) = amu(2,j+1,k+1)
      amu2(3,j+1,k+1) = amu(3,j+1,k+1)
      amu2(4,j+1,k+1) = amu(4,j+1,k+1)
      amu2(1,nx+j+1,k+1) = -amu(1,nx-j+1,k+1)
      amu2(2,nx+j+1,k+1) = amu(2,nx-j+1,k+1)
      amu2(3,nx+j+1,k+1) = amu(3,nx-j+1,k+1)
      amu2(4,nx+j+1,k+1) = -amu(4,nx-j+1,k+1)
      amu2(1,j+1,ny+k+1) = -amu(1,j+1,ny-k+1)
      amu2(2,j+1,ny+k+1) = amu(2,j+1,ny-k+1)
      amu2(3,j+1,ny+k+1) = -amu(3,j+1,ny-k+1)
      amu2(4,j+1,ny+k+1) = amu(4,j+1,ny-k+1)
      amu2(1,nx+j+1,ny+k+1) = amu(1,nx-j+1,ny-k+1)
      amu2(2,nx+j+1,ny+k+1) = amu(2,nx-j+1,ny-k+1)
      amu2(3,nx+j+1,ny+k+1) = -amu(3,nx-j+1,ny-k+1)
      amu2(4,nx+j+1,ny+k+1) = -amu(4,nx-j+1,ny-k+1)
   10 continue
      amu2(1,1,k+1) = 0.
      amu2(2,1,k+1) = amu(2,1,k+1)
      amu2(3,1,k+1) = amu(3,1,k+1)
      amu2(4,1,k+1) = 0.
      amu2(1,nx+1,k+1) = 0.
      amu2(2,nx+1,k+1) = amu(2,nx+1,k+1)
      amu2(3,nx+1,k+1) = amu(3,nx+1,k+1)
      amu2(4,nx+1,k+1) = 0.
      amu2(1,1,k+ny+1) = 0.
      amu2(2,1,k+ny+1) = amu(2,1,ny-k+1)
      amu2(3,1,k+ny+1) = -amu(3,1,ny-k+1)
      amu2(4,1,k+ny+1) = 0.
      amu2(1,nx+1,k+ny+1) = 0.
      amu2(2,nx+1,k+ny+1) = amu(2,nx+1,ny-k+1)
      amu2(3,nx+1,k+ny+1) = -amu(3,nx+1,ny-k+1)
      amu2(4,nx+1,k+ny+1) = 0.
   20 continue
      do 30 j = 1, nxs
      amu2(1,j+1,1) = 0.
      amu2(2,j+1,1) = amu(2,j+1,1)
      amu2(3,j+1,1) = 0.
      amu2(4,j+1,1) = amu(4,j+1,1)
      amu2(1,j+nx+1,1) = 0.
      amu2(2,nx+j+1,1) = amu(2,nx-j+1,1)
      amu2(3,j+nx+1,1) = 0.
      amu2(4,j+nx+1,1) = -amu(4,nx-j+1,1)
      amu2(1,j+1,ny+1) = 0.
      amu2(2,j+1,ny+1) = amu(2,j+1,ny+1)
      amu2(3,j+1,ny+1) = 0.
      amu2(4,j+1,ny+1) = amu(4,j+1,ny+1)
      amu2(1,j+nx+1,ny+1) = 0.
      amu2(2,nx+j+1,ny+1) = amu(2,nx-j+1,ny+1)
      amu2(3,j+nx+1,ny+1) = 0.
      amu2(4,j+nx+1,ny+1) = -amu(4,nx-j+1,ny+1)
   30 continue
      amu2(1,1,1) = 0.
      amu2(2,1,1) = amu(2,1,1)
      amu2(3,1,1) = 0.
      amu2(4,1,1) = 0.
      amu2(1,nx+1,1) = 0.
      amu2(2,nx+1,1) = amu(2,nx+1,1)
      amu2(3,nx+1,1) = 0.
      amu2(4,nx+1,1) = 0.
      amu2(1,1,ny+1) = 0.
      amu2(2,1,ny+1) = amu(2,1,ny+1)
      amu2(3,1,ny+1) = 0.
      amu2(4,1,ny+1) = 0.
      amu2(1,nx+1,ny+1) = 0.
      amu2(2,nx+1,ny+1) = amu(2,nx+1,ny+1)
      amu2(3,nx+1,ny+1) = 0.
      amu2(4,nx+1,ny+1) = 0.
      return
      end
      subroutine DBLSIN22M(amu,amu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a doubled tensor array amu2 from a tensor array
c amu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the xx-yy component is an odd function in
c both x and y, the xy component is an even function in both x and y
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array amu, must be >= nx+1
c nyv = third dimension of input array amu, must be >= ny+1
c nx2v = second dimension of output array amu2, must be >= 2*nx
c ny2 = third dimension of output array amu2, must be >= 2*ny
      implicit none
      real amu, amu2
      integer nx, ny, nxv, nyv, nx2v, ny2
      dimension amu(2,nxv,nyv), amu2(2,nx2v,ny2)
c local data
      integer j, k, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 20 k = 1, nys
      do 10 j = 1, nxs
      amu2(1,j+1,k+1) = amu(1,j+1,k+1)
      amu2(2,j+1,k+1) = amu(2,j+1,k+1)
      amu2(1,nx+j+1,k+1) = -amu(1,nx-j+1,k+1)
      amu2(2,nx+j+1,k+1) = amu(2,nx-j+1,k+1)
      amu2(1,j+1,ny+k+1) = -amu(1,j+1,ny-k+1)
      amu2(2,j+1,ny+k+1) = amu(2,j+1,ny-k+1)
      amu2(1,nx+j+1,ny+k+1) = amu(1,nx-j+1,ny-k+1)
      amu2(2,nx+j+1,ny+k+1) = amu(2,nx-j+1,ny-k+1)
   10 continue
      amu2(1,1,k+1) = 0.
      amu2(2,1,k+1) = amu(2,1,k+1)
      amu2(1,nx+1,k+1) = 0.
      amu2(2,nx+1,k+1) = amu(2,nx+1,k+1)
      amu2(1,1,k+ny+1) = 0.
      amu2(2,1,k+ny+1) = amu(2,1,ny-k+1)
      amu2(1,nx+1,k+ny+1) = 0.
      amu2(2,nx+1,k+ny+1) = amu(2,nx+1,ny-k+1)
   20 continue
      do 30 j = 1, nxs
      amu2(1,j+1,1) = 0.
      amu2(2,j+1,1) = amu(2,j+1,1)
      amu2(1,j+nx+1,1) = 0.
      amu2(2,nx+j+1,1) = amu(2,nx-j+1,1)
      amu2(1,j+1,ny+1) = 0.
      amu2(2,j+1,ny+1) = amu(2,j+1,ny+1)
      amu2(1,j+nx+1,ny+1) = 0.
      amu2(2,nx+j+1,ny+1) = amu(2,nx-j+1,ny+1)
   30 continue
      amu2(1,1,1) = 0.
      amu2(2,1,1) = amu(2,1,1)
      amu2(1,nx+1,1) = 0.
      amu2(2,nx+1,1) = amu(2,nx+1,1)
      amu2(1,1,ny+1) = 0.
      amu2(2,1,ny+1) = amu(2,1,ny+1)
      amu2(1,nx+1,ny+1) = 0.
      amu2(2,nx+1,ny+1) = amu(2,nx+1,ny+1)
      return
      end
      subroutine DCUPERPDX23(dcu,amu,nx,ny,nxv,ny2d)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux
c in 2-1/2d with dirichlet boundary conditions (zero potential)
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = -sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = -sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = -sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
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
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      implicit none
      integer nx, ny, nxv, ny2d
      complex dcu, amu
      dimension dcu(3,nxv,ny2d), amu(4,nxv,ny2d)
      integer ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2
      real at1, at2, at3, at4
      complex zero
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      at3 = at1*(dkxy*real(amu(1,j,k)) + dkxy2*real(amu(2,j,k)))
      at2 = dky*at3
      at3 = dkx*at3
      dcu(1,j,k) = cmplx(0.0,-at2)
      dcu(2,j,k) = cmplx(0.0,at3)
      at4 = dkx*aimag(amu(3,j,k)) + dky*aimag(amu(4,j,k))
      dcu(3,j,k) = cmplx(at4,0.0)
      dcu(1,j,k1) = cmplx(0.0,at2)
      dcu(2,j,k1) = cmplx(0.0,at3)
      dcu(3,j,k1) = cmplx(-at4,0.0)
   10 continue
   20 continue
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dcu(1,1,k) = cmplx(0.0,-dky*real(amu(2,1,k)))
      dcu(2,1,k) = zero
      dcu(3,1,k) = cmplx(dky*aimag(amu(4,1,k)),0.0)
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      dcu(1,j,1) = zero
      dcu(2,j,1) = cmplx(0.0,-dkx*real(amu(2,j,1)))
      dcu(3,j,1) = cmplx(dkx*aimag(amu(3,j,1)),0.0)
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
      subroutine ADCUPERPDX23(dcu,amu,nx,ny,nxv,ny2d)
c this subroutine calculates transverse part of the derivative of
c the current density from the momentum flux and acceleration density
c in 2-1/2d with dirichlet boundary conditions (zero potential)
c the derivative of the current is calculated using the equations:
c dcu(1,kx,ky) = dcu(1,kx,ky)-sqrt(-1)*(kx*vx*vx+ky*vx*vy)
c dcu(2,kx,ky) = dcu(2,kx,ky)-sqrt(-1)*(kx*vx*vy+ky*vy*vy)
c dcu(3,kx,ky) = dcu(3,kx,ky)-sqrt(-1)*(kx*vx*vz+ky*vy*vz)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
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
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      implicit none
      integer nx, ny, nxv, ny2d
      complex dcu, amu
      dimension dcu(3,nxv,ny2d), amu(4,nxv,ny2d)
      integer ny2, j, k, k1
      real dnx, dny, dky, dky2, dkx, dkx2, dkxy, dkxy2
      real at1, at2, at3, at4
      complex zero
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
c mode numbers 0 < kx < nx/ and 0 < ky < ny/
      do 20 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      dkx2 = dkx*dkx
      dkxy = dkx*dky
      dkxy2 = dky2 - dkx2
      at1 = 1.0/(dkx2 + dky2)
      at3 = at1*(dky*aimag(dcu(1,j,k)) - dkx*aimag(dcu(2,j,k)) - dkxy*re
     1al(amu(1,j,k)) - dkxy2*real(amu(2,j,k)))
      at2 = dky*at3
      at3 = -dkx*at3
      dcu(1,j,k) = cmplx(0.0,at2)
      dcu(2,j,k) = cmplx(0.0,at3)
      at4 = real(dcu(3,j,k)) + dkx*aimag(amu(3,j,k)) + dky*aimag(amu(4,j
     1,k))
      dcu(3,j,k) = cmplx(at4,0.0)
      dcu(1,j,k1) = cmplx(0.0,-at2)
      dcu(2,j,k1) = cmplx(0.0,at3)
      dcu(3,j,k1) = cmplx(-at4,0.0)
   10 continue
   20 continue
c mode numbers kx = 0, nx
      do 30 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dcu(1,1,k) = cmplx(0.0,aimag(dcu(1,1,k))-dky*real(amu(2,1,k)))
      dcu(2,1,k) = zero
      dcu(3,1,k) = cmplx(real(dcu(3,1,k))+dky*aimag(amu(4,1,k)),0.0)
      dcu(1,1,k1) = zero
      dcu(2,1,k1) = zero
      dcu(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      dcu(1,j,1) = zero
      dcu(2,j,1) = cmplx(0.0,aimag(dcu(2,j,1))-dkx*real(amu(2,j,1)))
      dcu(3,j,1) = cmplx(real(dcu(3,j,1))+dkx*aimag(amu(3,j,1)),0.0)
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
      subroutine EPOISDX23(dcu,exy,isign,fff,ax,ay,affp,wp0,ci,wf,nx,ny,
     1nxv,ny2d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c transverse electric field (or convolution of transverse electric field
c over particle shape), with dirichlet boundary conditions
c (zero potential).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c using algorithm described in J. Busnardo-Neto, P. L. Pritchett,
c A. T. Lin, and J. M. Dawson, J. Computational Phys. 23, 300 (1977).
c for isign = 0, input: isign,ax,ay,affp,wp0,nx,ny,nxv,ny2d, output:fff
c for isign =/ 0, input: dcu,fff,isign,ci,nx,ny,nxv,ny2d, output: exy,wf
c approximate flop count is: 16*nxc*nyc + 14*(nxc + nyc)
c where nxc = nx - 1, nyc = ny - 1
c if isign = 0, form factor array is prepared
c if isign = -1, smoothed transverse electric field is calculated
c using the equation:
c ex(kx,ky) = -ci*ci*g(kx,ky)*dcux(kx,ky)*s(kx,ky)
c ey(kx,ky) = -ci*ci*g(kx,ky)*dcuy(kx,ky)*s(kx,ky)
c ez(kx,ky) = -ci*ci*g(kx,ky)*dcuz(kx,ky)*s(kx,ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2+wp0*ci2*s(kx,ky)**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)), except for
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
c aimag(fff(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(fff(j,k)) = potential green's function g
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
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      implicit none
      integer isign, nx, ny, nxv, ny2d
      real ax, ay, affp, wp0, ci, wf
      complex dcu, exy, fff
      dimension dcu(3,nxv,ny2d), exy(3,nxv,ny2d)
      dimension fff(nxv,ny)
      integer ny2, j, k, k1
      real dnx, dny, ci2, wpc, dkx, dky, at1, at2, at3, at4
      complex zero
      double precision wp
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
      wpc = wp0*ci2
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         fff(j,k) = cmplx(affp,1.)
      else
         fff(j,k) = cmplx(affp*at4/(at3 + wpc*at4*at4),at4)
      endif
   10 continue
   20 continue
      return
c calculate smoothed transverse electric field and sum field energy
   30 if (isign.gt.0) go to 80
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 50 k = 2, ny
      k1 = ny2 - k
      do 40 j = 2, nx
      at2 = -ci2*real(fff(j,k))
      at1 = at2*aimag(fff(j,k))
      at2 = at2*at2
      exy(1,j,k) = cmplx(0.0,at1*aimag(dcu(1,j,k)))
      exy(2,j,k) = cmplx(0.0,at1*aimag(dcu(2,j,k)))
      exy(3,j,k) = cmplx(at1*real(dcu(3,j,k)),0.0)
      exy(1,j,k1) = -exy(1,j,k)
      exy(2,j,k1) = exy(2,j,k)
      exy(3,j,k1) = -exy(3,j,k)
      wp = wp + at2*(aimag(dcu(1,j,k))**2 + aimag(dcu(2,j,k))**2 + real(
     1dcu(3,j,k))**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 60 k = 2, ny
      k1 = ny2 - k
      at2 = -ci2*real(fff(1,k))
      at1 = at2*aimag(fff(1,k))
      at2 = at2*at2
      exy(1,1,k) = cmplx(0.0,at1*aimag(dcu(1,1,k)))
      exy(2,1,k) = cmplx(0.0,at1*aimag(dcu(2,1,k)))
      exy(3,1,k) = cmplx(at1*real(dcu(3,1,k)),0.0)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wp = wp + at2*(aimag(dcu(1,1,k))**2 + aimag(dcu(2,1,k))**2 + real(
     1dcu(3,1,k))**2)
   60 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 70 j = 2, nx
      at2 = -ci2*real(fff(j,1))
      at1 = at2*aimag(fff(j,1))
      at2 = at2*at2
      exy(1,j,1) = cmplx(0.0,at1*aimag(dcu(1,j,1)))
      exy(2,j,1) = cmplx(0.0,at1*aimag(dcu(2,j,1)))
      exy(3,j,1) = cmplx(at1*real(dcu(3,j,1)),0.0)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at2*(aimag(dcu(1,j,1))**2 + aimag(dcu(2,j,1))**2 + real(
     1dcu(3,j,1))**2)
   70 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*wp/real(fff(1,1))
      return
c calculate unsmoothed transverse electric field and sum field energy
   80 wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 100 k = 2, ny
      k1 = ny2 - k
      do 90 j = 2, nx
      at2 = -ci2*real(fff(j,k))
      at1 = at2*at2
      exy(1,j,k) = cmplx(0.0,at2*aimag(dcu(1,j,k)))
      exy(2,j,k) = cmplx(0.0,at2*aimag(dcu(2,j,k)))
      exy(3,j,k) = cmplx(at2*real(dcu(3,j,k)),0.0)
      exy(1,j,k1) = -exy(1,j,k)
      exy(2,j,k1) = exy(2,j,k)
      exy(3,j,k1) = -exy(3,j,k)
      wp = wp + at1*(aimag(dcu(1,j,k))**2 + aimag(dcu(2,j,k))**2 + real(
     1dcu(3,j,k))**2)
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 110 k = 2, ny
      k1 = ny2 - k
      at2 = -ci2*real(fff(1,k))
      at1 = at2*at2
      exy(1,1,k) = cmplx(0.0,at2*aimag(dcu(1,1,k)))
      exy(2,1,k) = cmplx(0.0,at2*aimag(dcu(2,1,k)))
      exy(3,1,k) = cmplx(at2*real(dcu(3,1,k)),0.0)
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wp = wp + at1*(aimag(dcu(1,1,k))**2 + aimag(dcu(2,1,k))**2 + real(
     1dcu(3,1,k))**2)
  110 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 120 j = 2, nx
      at2 = -ci2*real(fff(j,1))
      at1 = at2*at2
      exy(1,j,1) = cmplx(0.0,at2*aimag(dcu(1,j,1)))
      exy(2,j,1) = cmplx(0.0,at2*aimag(dcu(2,j,1)))
      exy(3,j,1) = cmplx(at2*real(dcu(3,j,1)),0.0)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      wp = wp + at1*(aimag(dcu(1,j,1))**2 + aimag(dcu(2,j,1))**2 + real(
     1dcu(3,j,1))**2)
  120 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = float(nx*ny)*wp/real(fff(1,1))
      return
      end
      subroutine HAFDBL2N(daxy,daxy2,nx,ny,nxe,nye,nx2v,ny2,ndim)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c with variable first dimension
c nx/ny = system length in x/y direction
c nxe = second dimension of output array daxy, must be >= nx+1
c nye = third dimension of ouput array daxy, must be >= ny+1
c nx2v = second dimension of input array daxy2, must be >= 2*nx
c ny2 = third dimension of input array daxy2, must be >= 2*ny
c ndim = first dimension of arrays daxy, daxy2
      implicit none
      real daxy, daxy2
      integer nx, ny, nxe, nye, nx2v, ny2, ndim
      dimension daxy(ndim,nxe,nye), daxy2(ndim,nx2v,ny2)
c local data
      integer i, j, k, nx1, ny1
      nx1 = nx + 1
      ny1 = ny + 1
      do 30 k = 1, ny1
      do 20 j = 1, nx1
      do 10 i = 1, ndim
      daxy(i,j,k) = daxy2(i,j,k)
   10 continue
   20 continue
   30 continue
      return
      end
      subroutine LDCGUARD2(daxy,nx,ny,nxe,nye,ndim)
c this subroutine replicates vector field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation on both boundaries
c nx/ny = system length in x/y direction
c nxe = first dimension of input array bxy, must be >= nx+3
c nxe = second dimension of input array bxy, must be >= ny+3
      implicit none
      real daxy
      integer nx, ny, nxe, nye, ndim
      dimension daxy(ndim,nxe,nye)
c local data
      integer i, j, k, ny1, nx3
      ny1 = ny + 1
      nx3 = nx + 3
      do 20 k = 1, ny1
      do 10 i = 1, ndim
      daxy(i,1,k+1) = 2.*daxy(i,2,k+1) - daxy(i,3,k+1)
      daxy(i,nx3,k+1) = 2.*daxy(i,nx+2,k+1) - daxy(i,nx+1,k+1)
   10 continue
   20 continue
      do 40 j = 1, nx3
      do 30 i = 1, ndim
      daxy(i,j,1) = 2.*daxy(i,j,2) - daxy(i,j,3)
      daxy(i,j,ny+3) = 2.*daxy(i,j,ny+2) - daxy(i,j,ny+1)
   30 continue
   40 continue
      return
      end
      subroutine APOISDX23(cu,daxy,axy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,
     1nxv,ny2d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c smoothed derivative of vector potential, or smoothed vector potential,
c with dirichlet boundary conditions (zero potential).
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,ny2d, output: ffd
c for isign = -1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d,output: daxy,wm
c approximate flop count is: 22*nxc*nyc + 15*(nxc + nyc)
c for isign = 1, input: cu,ffd,isign,ci,nx,ny,nxv,ny2d, output: axy,wm
c approximate flop count is: 15*nxc*nyc + 12*(nxc + nyc)
c where nxc = nx - 1, nyc = ny - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c daxy(1,kx,ky)=dax/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c daxy(2,kx,ky)=day/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cuy(kx,ky)*s(kx,ky)
c daxy(3,kx,ky)=dax/dy = ci*ci*sqrt(-1)*ky*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c daxy(4,kx,ky)=daz/dx = ci*ci*sqrt(-1)*kx*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c daxy(5,kx,ky)=daz/dy = ci*ci*sqrt(-1)*ky*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), 
c if isign = 1, smoothed vector potential is calculated using:
c axy(1,kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)*s(kx,ky)
c axy(2,kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)*s(kx,ky)
c axy(3,kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)*s(kx,ky)
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c daxy(i,j,k) = i component of smoothed derivative of vector potential
c axy(i,j,k) = i component of smoothed vector potential
c all for fourier mode (j-1,k-1)
c aimag(ffd(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffd(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
      implicit none
      integer isign, nx, ny, nxv, ny2d
      real ax, ay, affp, ci, wm
      complex cu, daxy, axy, ffd
      dimension cu(3,nxv,ny2d), daxy(5,nxv,ny2d), axy(3,nxv,ny2d)
      dimension ffd(nxv,ny)
      integer ny2, j, k, k1
      real dnx, dny, ci2, dkx, dky, at0, at1, at2, at3, at4, at5
      complex zero
      double precision wp
      ny2 = 2*ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny + ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, ny
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffd(j,k) = cmplx(affp,1.)
      else
         ffd(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate derivative of vector potential and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 50 k = 2, ny
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at0 = ci2*real(ffd(j,k))*aimag(ffd(j,k))
      at2 = dnx*float(j - 1)*at0
      at3 = dky*at0
      at5 = at3*real(cu(3,j,k))
      at4 = at2*real(cu(3,j,k))
      at3 = -at3*aimag(cu(1,j,k))
      at1 = at2*aimag(cu(1,j,k))
      at2 = -at2*aimag(cu(2,j,k))
      daxy(1,j,k) = cmplx(-at1,0.0)
      daxy(2,j,k) = cmplx(at2,0.0)
      daxy(3,j,k) = cmplx(at3,0.0)
      daxy(4,j,k) = cmplx(0.0,at4)
      daxy(5,j,k) = cmplx(0.0,at5)
      daxy(1,j,k1) = cmplx(at1,0.0)
      daxy(2,j,k1) = cmplx(at2,0.0)
      daxy(3,j,k1) = cmplx(at3,0.0)
      daxy(4,j,k1) = cmplx(0.0,-at4)
      daxy(5,j,k1) = cmplx(0.0,at5)
      wp = wp + at0*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2 + real(cu
     1(3,j,k))**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 60 k = 2, ny
      k1 = ny2 - k
      at0 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      at3 = dny*float(k - 1)*at0
      daxy(1,1,k) = zero
      daxy(2,1,k) = zero
      daxy(3,1,k) = cmplx(-at3*aimag(cu(1,1,k)),0.0)
      daxy(4,1,k) = zero
      daxy(5,1,k) = cmplx(0.0,at3*real(cu(3,1,k)))
      daxy(1,1,k1) = zero
      daxy(2,1,k1) = zero
      daxy(3,1,k1) = zero
      daxy(4,1,k1) = zero
      daxy(5,1,k1) = zero
      wp = wp + at0*(aimag(cu(1,1,k))**2 + aimag(cu(2,1,k))**2 + real(cu
     1(3,1,k))**2)
   60 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 70 j = 2, nx
      at0 = ci2*real(ffd(j,1))*aimag(ffd(j,1))
      at2 = dnx*float(j - 1)*at0
      daxy(1,j,1) = cmplx(-at2*aimag(cu(1,j,1)),0.0)
      daxy(2,j,1) = cmplx(-at2*aimag(cu(2,j,1)),0.0)
      daxy(3,j,1) = zero
      daxy(4,j,1) = cmplx(0.0,at2*real(cu(3,j,1)))
      daxy(5,j,1) = zero
      daxy(1,j,k1) = zero
      daxy(2,j,k1) = zero
      daxy(3,j,k1) = zero
      daxy(4,j,k1) = zero
      daxy(5,j,k1) = zero
      wp = wp + at0*(aimag(cu(1,j,1))**2 + aimag(cu(2,j,1))**2 + real(cu
     1(3,j,1))**2)
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
c mode numbers 0 < kx < nx and 0 < ky < ny
      do 100 k = 2, ny
      k1 = ny2 - k
      do 90 j = 2, nx
      at0 = ci2*real(ffd(j,k))*aimag(ffd(j,k))
      at3 = at0*real(cu(3,j,k))
      at2 = at0*aimag(cu(2,j,k))
      at1 = at0*aimag(cu(1,j,k))
      axy(1,j,k) = cmplx(0.0,at1)
      axy(2,j,k) = cmplx(0.0,at2)
      axy(3,j,k) = cmplx(at3,0.0)
      axy(1,j,k1) = cmplx(0.0,-at1)
      axy(2,j,k1) = cmplx(0.0,at2)
      axy(3,j,k1) = cmplx(-at3,0.0)
      wp = wp + at0*(aimag(cu(1,j,k))**2 + aimag(cu(2,j,k))**2 + real(cu
     1(3,j,k))**2)
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 110 k = 2, ny
      k1 = ny2 - k
      at1 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      axy(1,1,k) = cmplx(0.0,at1*aimag(cu(1,1,k)))
      axy(2,1,k) = cmplx(0.0,at1*aimag(cu(2,1,k)))
      axy(3,1,k) = cmplx(at1*real(cu(3,1,k)),0.0)
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
      wp = wp + at1*(aimag(cu(1,1,k))**2 + aimag(cu(2,1,k))**2 + real(cu
     1(3,1,k))**2)
  110 continue
c mode numbers ky = 0, ny
      k1 = ny + 1
      do 120 j = 2, nx
      at1 = ci2*real(ffd(j,1))*aimag(ffd(j,1))
      axy(1,j,1) = cmplx(0.0,at1*aimag(cu(1,j,1)))
      axy(2,j,1) = cmplx(0.0,at1*aimag(cu(2,j,1)))
      axy(3,j,1) = cmplx(at1*real(cu(3,j,1)),0.0)
      axy(1,j,k1) = zero
      axy(2,j,k1) = zero
      axy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(1,j,1))**2 + aimag(cu(2,j,1))**2 + real(cu
     1(3,j,1))**2)
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
