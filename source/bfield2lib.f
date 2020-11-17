c 2d PIC library for solving field equations with mixed
c dirichlet/periodic boundary conditions
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: june 29, 2007
      subroutine MCGUARD2(fxy,nx,ny,nxe,nye)
c this subroutine replicates field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c nx/ny = system length in x/y direction
c nxe = first dimension of input array fxy, must be >= nx+3
c nxe = second dimension of input array fxy, must be >= ny+3
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer i, j, k, nx3
      nx3 = nx + 3
      do 10 k = 1, ny
      fxy(1,1,k+1) = 2.*fxy(1,2,k+1) - fxy(1,3,k+1)
      fxy(2,1,k+1) = 2.*fxy(2,2,k+1) - fxy(2,3,k+1)
      fxy(1,nx3,k+1) = 2.*fxy(1,nx+2,k+1) - fxy(1,nx+1,k+1)
      fxy(2,nx3,k+1) = 2.*fxy(2,nx+2,k+1) - fxy(2,nx+1,k+1)
   10 continue
c copy periodic edge
      do 30 j = 1, nx3
      do 20 i = 1, 2
      fxy(i,j,1) = fxy(i,j,ny+1)
      fxy(i,j,ny+2) = fxy(i,j,2)
      fxy(i,j,ny+3) = fxy(i,j,3)
   20 continue
   30 continue
      return
      end
      subroutine MDGUARD2(q,nx,ny,nxe,nye)
c this subroutine replicates scalar field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c nx/ny = system length in x/y direction
c nxe = first dimension of input array q, must be >= nx+3
c nxe = second dimension of input array q, must be >= ny+3
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, k, nx3
      nx3 = nx + 3
      do 10 k = 1, ny
      q(1,k+1) = 2.*q(2,k+1) - q(3,k+1)
      q(nx3,k+1) = 2.*q(nx+2,k+1) - q(nx+1,k+1)
   10 continue
c copy periodic edge
      do 20 j = 1, nx3
      q(j,1) = q(j,ny+1)
      q(j,ny+2) = q(j,2)
      q(j,ny+3) = q(j,3)
   20 continue
      return
      end
      subroutine MBGUARD2(bxy,nx,ny,nxe,nye)
c this subroutine replicates vector field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c nx/ny = system length in x/y direction
c nxe = first dimension of input array bxy, must be >= nx+3
c nxe = second dimension of input array bxy, must be >= ny+3
      implicit none
      real bxy
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
c local data
      integer i, j, k, nx3
      nx3 = nx + 3
      do 10 k = 1, ny
      bxy(1,1,k+1) = 2.*bxy(1,2,k+1) - bxy(1,3,k+1)
      bxy(2,1,k+1) = 2.*bxy(2,2,k+1) - bxy(2,3,k+1)
      bxy(3,1,k+1) = 2.*bxy(3,2,k+1) - bxy(3,3,k+1)
      bxy(1,nx3,k+1) = 2.*bxy(1,nx+2,k+1) - bxy(1,nx+1,k+1)
      bxy(2,nx3,k+1) = 2.*bxy(2,nx+2,k+1) - bxy(2,nx+1,k+1)
      bxy(3,nx3,k+1) = 2.*bxy(3,nx+2,k+1) - bxy(3,nx+1,k+1)
   10 continue
c copy periodic edge
      do 50 j = 1, nx3
      do 40 i = 1, 3
      bxy(i,j,1) = bxy(i,j,ny+1)
      bxy(i,j,ny+2) = bxy(i,j,2)
      bxy(i,j,ny+3) = bxy(i,j,3)
   40 continue
   50 continue
      return
      end
      subroutine MSCGUARD2(cu,xj0,yj0,zj0,nx,ny,ngx,nxe,nye)
c initialize extended non-periodic field in x and periodic in y
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, ngx, nxe, nye
      dimension cu(3,nxe,nye)
      integer i, j, k, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 30 k = 1, ny
      do 10 j = 2, nxg
      cu(1,j+ngx+1,k+1) = xj0
      cu(2,j+ngx+1,k+1) = yj0
      cu(3,j+ngx+1,k+1) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+1) = 0.
      cu(i,2,k+1) = 0.
      cu(i,nx+2,k+1) = 0.
      cu(i,nx+3,k+1) = 0.
   20 continue
      cu(1,ngx+2,k+1) = .5*xj0
      cu(2,ngx+2,k+1) = .5*yj0
      cu(3,ngx+2,k+1) = .5*zj0
      cu(1,nx-ngx+2,k+1) = .5*xj0
      cu(2,nx-ngx+2,k+1) = .5*yj0
      cu(3,nx-ngx+2,k+1) = .5*zj0
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 3
      cu(i,j,1) = 0.
      cu(i,j,ny+2) = 0.
      cu(i,j,ny+3) = 0.
   40 continue
   50 continue
      return
      end
      subroutine MSCGUARD22(cu,xj0,yj0,nx,ny,ngx,nxe,nye)
c initialize extended non-periodic field in x and periodic in y
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, ngx, nxe, nye
      dimension cu(2,nxe,nye)
      integer i, j, k, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 30 k = 1, ny
      do 10 j = 2, nxg
      cu(1,j+ngx+1,k+1) = xj0
      cu(2,j+ngx+1,k+1) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k+1) = 0.
      cu(i,2,k+1) = 0.
      cu(i,nx+2,k+1) = 0.
      cu(i,nx+3,k+1) = 0.
   20 continue
      cu(1,ngx+2,k+1) = .5*xj0
      cu(2,ngx+2,k+1) = .5*yj0
      cu(1,nx-ngx+2,k+1) = .5*xj0
      cu(2,nx-ngx+2,k+1) = .5*yj0
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 2
      cu(i,j,1) = 0.
      cu(i,j,ny+2) = 0.
      cu(i,j,ny+3) = 0.
   40 continue
   50 continue
      return
      end
      subroutine MSGUARD2(q,qi0,nx,ny,ngx,nxe,nye)
c initialize extended non-periodic scalar field in x and periodic in y
c ngx = (0,1) = number of grid cells away from edge
c quadratic interpolation
      implicit none
      real q, qi0
      integer nx, ny, ngx, nxe, nye
      dimension q(nxe,nye)
      integer j, k, nxg, nx3
      if ((ngx.lt.0).or.(ngx.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      do 20 k = 1, ny
      do 10 j = 2, nxg
      q(j+ngx+1,k+1) = qi0
   10 continue
      q(1,k+1) = 0.
      q(2,k+1) = 0.
      q(nx+2,k+1) = 0.
      q(nx+3,k+1) = 0.
      q(ngx+2,k+1) = .5*qi0
      q(nx-ngx+2,k+1) = .5*qi0
   20 continue
      do 30 j = 1, nx3
      q(j,1) = 0.
      q(j,ny+2) = 0.
      q(j,ny+3) = 0.
   30 continue
      return
      end
      subroutine MACGUARD2(cu,nx,ny,nxe,nye)
c this subroutine adds adds up guard cells, bounded in x, periodic in y
c with quadratic interpolation, disabling quadratic interpolation
c in the x direction within half a cell of the the edges and reducing
c it to linear interpolation. 
c in the y direction, periodicity is assumed.
c nx/ny = system length in x/y direction
c nxe = first dimension of input array q, must be >= nx+3
c nye = second dimension of input array q, must be >= ny+3
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
      cu(i,j+1,2) = cu(i,j+1,2) + cu(i,j+1,ny+2)
      cu(i,j+1,3) = cu(i,j+1,3) + cu(i,j+1,ny+3)
      cu(i,j+1,ny+1) = cu(i,j+1,ny+1) + cu(i,j+1,1)
      cu(i,j+1,ny+2) = 0.0
      cu(i,j+1,ny+3) = 0.0
      cu(i,j+1,1) = 0.0
   30 continue
   40 continue
      return
      end
      subroutine MACGUARD22(cu,nx,ny,nxe,nye)
c this subroutine adds adds up guard cells, bounded in x, periodic in y
c with quadratic interpolation, disabling quadratic interpolation
c in the x direction within half a cell of the the edges and reducing
c it to linear interpolation. 
c in the y direction, periodicity is assumed.
c nx/ny = system length in x/y direction
c nxe = first dimension of input array q, must be >= nx+3
c nye = second dimension of input array q, must be >= ny+3
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
      cu(i,j+1,2) = cu(i,j+1,2) + cu(i,j+1,ny+2)
      cu(i,j+1,3) = cu(i,j+1,3) + cu(i,j+1,ny+3)
      cu(i,j+1,ny+1) = cu(i,j+1,ny+1) + cu(i,j+1,1)
      cu(i,j+1,ny+2) = 0.0
      cu(i,j+1,ny+3) = 0.0
      cu(i,j+1,1) = 0.0
   30 continue
   40 continue
      return
      end
      subroutine MAGUARD2(q,nx,ny,nxe,nye)
c this subroutine adds adds up guard cells, bounded in x, periodic in y
c with quadratic interpolation, disabling quadratic interpolation
c in the x direction within half a cell of the the edges and reducing
c it to linear interpolation, for scalar field
c in the y direction, periodicity is assumed.
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
      q(j+1,2) = q(j+1,2) + q(j+1,ny+2)
      q(j+1,3) = q(j+1,3) + q(j+1,ny+3)
      q(j+1,ny+1) = q(j+1,ny+1) + q(j+1,1)
      q(j+1,ny+2) = 0.0
      q(j+1,ny+3) = 0.0
      q(j+1,1) = 0.0
   20 continue
      return
      end
      subroutine MCGUARD2L(fxy,nx,ny,nxe,nye)
c this subroutine replicates extended field, periodic in y,
c with linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array fxy, must be >= nx+1
c nye = first dimension of output array fxy, must be >= ny+1
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye)
c local data
      integer j, nx1
      nx1 = nx + 1
      do 10 j = 1, nx1
      fxy(1,j,ny+1) = fxy(1,j,1)
      fxy(2,j,ny+1) = fxy(2,j,1)
   10 continue
      return
      end
      subroutine MDGUARD2L(q,nx,ny,nxe,nye)
c this subroutine replicates extended scalar field, periodic in y,
c with linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array q, must be >= nx+1
c nye = first dimension of output array q, must be >= ny+1
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
c local data
      integer j, nx1
      nx1 = nx + 1
      do 10 j = 1, nx1
      q(j,ny+1) = q(j,1)
   10 continue
      return
      end
      subroutine MBGUARD2L(fxy,nx,ny,nxe,nye)
c this subroutine replicates extended field, periodic in y,
c with linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array fxy, must be >= nx+1
c nye = first dimension of output array fxy, must be >= ny+1
      implicit none
      real fxy
      integer nx, ny, nxe, nye
      dimension fxy(3,nxe,nye)
c local data
      integer j, nx1
      nx1 = nx + 1
      do 10 j = 1, nx1
      fxy(1,j,ny+1) = fxy(1,j,1)
      fxy(2,j,ny+1) = fxy(2,j,1)
      fxy(3,j,ny+1) = fxy(3,j,1)
   10 continue
      return
      end
      subroutine MSCGUARD2L(cu,xj0,yj0,zj0,nx,ny,ngx,nxe,nye)
c initialize extended semi-periodic field
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real cu, xj0, yj0, zj0
      integer nx, ny, ngx, nxe, nye
      dimension cu(3,nxe,nye)
      integer i, j, k, nxg, nx1
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 30 k = 1, ny
      do 10 j = 2, nxg
      cu(1,j+ngx,k) = xj0
      cu(2,j+ngx,k) = yj0
      cu(3,j+ngx,k) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k) = 0.
      cu(i,nx+1,k) = 0.
   20 continue
      cu(1,ngx+1,k) = .5*xj0
      cu(2,ngx+1,k) = .5*yj0
      cu(3,ngx+1,k) = .5*zj0
      cu(1,nx-ngx+1,k) = .5*xj0
      cu(2,nx-ngx+1,k) = .5*yj0
      cu(3,nx-ngx+1,k) = .5*zj0
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 3
      cu(i,j,ny+1) = 0.
   40 continue
   50 continue
      return
      end
      subroutine MSCGUARD22L(cu,xj0,yj0,nx,ny,ngx,nxe,nye)
c initialize extended semi-periodic field
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real cu, xj0, yj0
      integer nx, ny, ngx, nxe, nye
      dimension cu(2,nxe,nye)
      integer i, j, k, nxg, nx1
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 30 k = 1, ny
      do 10 j = 2, nxg
      cu(1,j+ngx,k) = xj0
      cu(2,j+ngx,k) = yj0
   10 continue
      do 20 i = 1, 2
      cu(i,1,k) = 0.
      cu(i,nx+1,k) = 0.
   20 continue
      cu(1,ngx+1,k) = .5*xj0
      cu(2,ngx+1,k) = .5*yj0
      cu(1,nx-ngx+1,k) = .5*xj0
      cu(2,nx-ngx+1,k) = .5*yj0
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 2
      cu(i,j,ny+1) = 0.
   40 continue
   50 continue
      return
      end
      subroutine MSGUARD2L(q,qi0,nx,ny,ngx,nxe,nye)
c initialize extended semi-periodic field
c ngx = (0,1) = number of grid cells away from edge
c linear interpolation
      implicit none
      real q, qi0
      integer nx, ny, ngx, nxe, nye
      dimension q(nxe,nye)
      integer j, k, nxg, nx1
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      do 20 k = 1, ny
      do 10 j = 2, nxg
      q(j+ngx,k) = qi0
   10 continue
      q(1,k) = 0.
      q(nx+1,k) = 0.
      q(ngx+1,k) = .5*qi0
      q(nx-ngx+1,k) = .5*qi0
   20 continue
      do 30 j = 1, nx1
      q(j,ny+1) = 0.
   30 continue
      return
      end
      subroutine MACGUARD2L(cu,nx,ny,nxe,nye)
c accumulate extended periodic field in y direction
c linear interpolation
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
      integer i, j
c accumulate edges of extended field
      do 20 j = 1, nx
      do 10 i = 1, 3
      cu(i,j,1) = cu(i,j,1) + cu(i,j,ny+1)
      cu(i,j,ny+1) = 0.0
   10 continue
   20 continue
      return
      end
      subroutine MACGUARD22L(cu,nx,ny,nxe,nye)
c accumulate extended periodic field in y direction
c linear interpolation
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(2,nxe,nye)
      integer i, j
c accumulate edges of extended field
      do 20 j = 1, nx
      do 10 i = 1, 2
      cu(i,j,1) = cu(i,j,1) + cu(i,j,ny+1)
      cu(i,j,ny+1) = 0.0
   10 continue
   20 continue
      return
      end
      subroutine MAGUARD2L(q,nx,ny,nxe,nye)
c accumulate extended periodic field in y direction
c linear interpolation, for scalar field
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
      integer j
c accumulate edges of extended field
      do 10 j = 1, nx
      q(j,1) = q(j,1) + q(j,ny+1)
      q(j,ny+1) = 0.0
   10 continue
      return
      end
      subroutine SGLSIN2C(cu,cu1,nx,ny,nxv,nyv,nx2v)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an even function in x,
c and y component is an odd function in x.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx+1
c nyv = third dimension of input array cu, must be >= ny
c nx2v = second dimension of output array cu1, must be >= 2*nx
      implicit none
      real cu, cu1
      integer nx, ny, nxv, nyv, nx2v
      dimension cu(2,nxv,nyv), cu1(2,nx2v,nyv)
c local data
      integer j, k, nxs
c copy to double array
      nxs = nx - 1
      do 20 k = 1, ny
      do 10 j = 1, nxs
      cu1(1,j+1,k) = cu(1,j+1,k)
      cu1(2,j+1,k) = cu(2,j+1,k)
      cu1(1,nx+j+1,k) = cu(1,nx-j+1,k)
      cu1(2,nx+j+1,k) = -cu(2,nx-j+1,k)
   10 continue
      cu1(1,1,k) = cu(1,1,k)
      cu1(2,1,k) = 0.
      cu1(1,nx+1,k) = cu(1,nx+1,k)
      cu1(2,nx+1,k) = 0.
   20 continue
      return
      end
      subroutine SGLSIN2D(q,q1,nx,ny,nxv,nyv,nx2v)
c this subroutine creates an odd array q1 from an array q, so that
c a 1d sine transform can be performed with a 2d real to complex fft.
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = first dimension of input array q, must be >= nx+1
c nyv = second dimension of input array q, must be >= ny
c nx2v = first dimension of output array q1, must be >= 2*nx
      implicit none
      real q, q1
      integer nx, ny, nxv, nyv, nx2v
      dimension q(nxv,nyv), q1(nx2v,nyv)
c local data
      integer j, k, nxs
c copy to double array
      nxs = nx - 1
      do 20 k = 1, ny
      do 10 j = 1, nxs
      q1(j+1,k) = q(j+1,k)
      q1(nx+j+1,k) = -q(nx-j+1,k)
   10 continue
      q1(1,k) = 0.
      q1(nx+1,k) = 0.
   20 continue
      return
      end
      subroutine SGLSIN2B(cu,cu1,nx,ny,nxv,nyv,nx2v)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an even function in x,
c y and z components are an odd function in x.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx+1
c nyv = third dimension of input array cu, must be >= ny
c nx2v = second dimension of output array cu1, must be >= 2*nx
      implicit none
      real cu, cu1
      integer nx, ny, nxv, nyv, nx2v
      dimension cu(3,nxv,nyv), cu1(3,nx2v,nyv)
c local data
      integer j, k, nxs
c copy to double array
      nxs = nx - 1
      do 20 k = 1, ny
      do 10 j = 1, nxs
      cu1(1,j+1,k) = cu(1,j+1,k)
      cu1(2,j+1,k) = cu(2,j+1,k)
      cu1(3,j+1,k) = cu(3,j+1,k)
      cu1(1,nx+j+1,k) = cu(1,nx-j+1,k)
      cu1(2,nx+j+1,k) = -cu(2,nx-j+1,k)
      cu1(3,nx+j+1,k) = -cu(3,nx-j+1,k)
   10 continue
      cu1(1,1,k) = cu(1,1,k)
      cu1(2,1,k) = 0.
      cu1(3,1,k) = 0.
      cu1(1,nx+1,k) = cu(1,nx+1,k)
      cu1(2,nx+1,k) = 0.
      cu1(3,nx+1,k) = 0.
   20 continue
      return
      end
      subroutine HAFSGL2C(fxy,fxy1,nx,ny,nxe,nye,nx2v)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array fxy, must be >= nx+1
c nye = second dimension of ouput array fxy, must be >= ny+1
c nx2v = first dimension of input array fxy1, must be >= 2*nx
      implicit none
      real fxy, fxy1
      integer nx, ny, nxe, nye, nx2v
      dimension fxy(2,nxe,nye), fxy1(2,nx2v,nye)
c local data
      integer j, k, nx1
      nx1 = nx + 1
      do 20 k = 1, ny
      do 10 j = 1, nx1
      fxy(1,j,k) = fxy1(1,j,k)
      fxy(2,j,k) = fxy1(2,j,k)
   10 continue
   20 continue
      return
      end
      subroutine HAFSGL2D(q,q1,nx,ny,nxe,nye,nx2v)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array q, must be >= nx+1
c nye = second dimension of ouput array q, must be >= ny+1
c nx2v = first dimension of input array q1, must be >= 2*nx
      implicit none
      real q, q1
      integer nx, ny, nxe, nye, nx2v
      dimension q(nxe,nye), q1(nx2v,nye)
c local data
      integer j, k, nx1
      nx1 = nx + 1
      do 20 k = 1, ny
      do 10 j = 1, nx1
      q(j,k) = q1(j,k)
   10 continue
   20 continue
      return
      end
      subroutine HAFSGL2B(bxy,bxy1,nx,ny,nxe,nye,nx2v)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx/ny = system length in x/y direction
c nxe = first dimension of output array bxy, must be >= nx+1
c nye = second dimension of ouput array bxy, must be >= ny+1
c nx2v = first dimension of input array bxy1, must be >= 2*nx
      implicit none
      real bxy, bxy1
      integer nx, ny, nxe, nye, nx2v
      dimension bxy(3,nxe,nye), bxy1(3,nx2v,nye)
c local data
      integer j, k, nx1
      nx1 = nx + 1
      do 20 k = 1, ny
      do 10 j = 1, nx1
      bxy(1,j,k) = bxy1(1,j,k)
      bxy(2,j,k) = bxy1(2,j,k)
      bxy(3,j,k) = bxy1(3,j,k)
   10 continue
   20 continue
      return
      end
      subroutine POISMX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,nyv,
     1nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign = -1, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd,
c output: fx,fy,we
c approximate flop count is: 18*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fx,we
c approximate flop count is: 9*nxc*nyc + 5*(nxc + nyc)
c for isign = 2, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fy
c approximate flop count is: 3*nxc*nyc + 1*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
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
c ffb(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffb(2*j-1,k) = potential green's function g for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nx2v = first dimension of field arrays, must be >= 2*nx
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      dimension q(nx2v,nyv), fx(nx2v,nyv), fy(nx2v,nyv)
      dimension ffb(nx2v,nyhd)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      if (at3.eq.0.) then
         ffb(2*j,k) = 1.
         ffb(2*j-1,k) = affp
      else
         ffb(2*j,k) = exp(-.5*((dkx*ax)**2 + at2))
         ffb(2*j-1,k) = affp*ffb(2*j,k)/at3
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ffb(2*j-1,k)*ffb(2*j,k)
      at2 = dnx*float(j - 1)*at1
      at3 = at2*q(2*j,k)
      at2 = at2*q(2*j-1,k)
      fx(2*j-1,k) = at3
      fx(2*j,k) = -at2
      fx(2*j-1,k1) = at3
      fx(2*j,k1) = at2
      at3 = dky*at1
      at2 = at3*q(2*j,k)
      at3 = at3*q(2*j-1,k)
      fy(2*j-1,k) = at2
      fy(2*j,k) = -at3
      fy(2*j-1,k1) = -at2
      fy(2*j,k1) = -at3
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2)
   40 continue
c mode numbers kx = 0, nx
      fx(1,k) = 0.
      fx(2,k) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      fy(1,k) = 0.
      fy(2,k) = 0.
      fy(1,k1) = 0.
      fy(2,k1) = 0.
   50 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nx
      at1 = ffb(2*j-1,1)*ffb(2*j,1)
      at3 = at1*q(2*j,1)
      at2 = dnx*float(j - 1)*at3
      fx(2*j-1,1) = at2
      fx(2*j,1) = 0.
      fx(2*j-1,k1) = 0.
      fx(2*j,k1) = 0.
      fy(2*j-1,1) = 0.
      fy(2*j,1) = 0.
      fy(2*j-1,k1) = 0.
      fy(2*j,k1) = 0.
      wp = wp + at1*q(2*j,1)**2
   60 continue
      fx(1,1) = 0.
      fx(2,1) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      fy(1,1) = 0.
      fy(2,1) = 0.
      fy(1,k1) = 0.
      fy(2,k1) = 0.
      we = float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 90 k = 2, nyh
      k1 = ny2 - k
      do 80 j = 2, nx
      at2 = ffb(2*j-1,k)
      at1 = at2*ffb(2*j,k)
      at3 = at2*q(2*j-1,k)
      at2 = at2*q(2*j,k)
      fx(2*j-1,k) = at3
      fx(2*j,k) = at2
      fx(2*j-1,k1) = -at3
      fx(2*j,k1) = at2
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2)
   80 continue
c mode numbers kx = 0, nx
      fx(1,k) = 0.
      fx(2,k) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
   90 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 100 j = 2, nx
      at2 = ffb(2*j-1,1)
      at1 = at2*ffb(2*j,1)
      fx(2*j-1,1) = 0.
      fx(2*j,1) = at2*q(2*j,1)
      fx(2*j-1,k1) = 0.
      fx(2*j,k1) = 0.
      wp = wp + at1*q(2*j,1)**2
  100 continue
      fx(1,1) = 0.
      fx(2,1) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      we = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  110 do 130 k = 2, nyh
      k1 = ny2 - k
      do 120 j = 2, nx
      at1 = ffb(2*j,k)
      at3 = at1*q(2*j-1,k)
      at2 = at1*q(2*j,k)
      fy(2*j-1,k) = at3
      fy(2*j,k) = at2
      fy(2*j-1,k1) = -at3
      fy(2*j,k1) = at2
  120 continue
c mode numbers kx = 0, nx
      fy(1,k) = 0.
      fy(2,k) = 0.
      fy(1,k1) = 0.
      fy(2,k1) = 0.
  130 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 140 j = 2, nx
      at1 = ffb(2*j,1)
      fy(2*j-1,1) = 0.
      fy(2*j,1) = at1*q(2*j,1)
      fy(2*j-1,k1) = 0.
      fy(2*j,k1) = 0.
  140 continue
      fy(1,1) = 0.
      fy(2,1) = 0.
      fy(1,k1) = 0.
      fy(2,k1) = 0.
      return
      end
      subroutine POISM2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nyeh,n
     1xv,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign = -1, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd,
c output: fx,fy,we
c approximate flop count is: 19*nxc*nyc + 8*(nxc + nyc)
c for isign = 1, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fx,we
c approximate flop count is: 9*nxc*nyc + 5*(nxc + nyc)
c for isign = 2, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fy
c approximate flop count is: 3*nxc*nyc + 1*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
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
c ffb(2*j,k) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c ffb(2*j-1,k) = potential green's function g for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= ny/2
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fx, fy, ffb, zero
      dimension q(nxe,nyeh), fx(nxe,nyeh), fy(nxe,nyeh)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      fx(j,k) = -at2*q(j,k)
      fy(j,k) = at3*cmplx(aimag(q(j,k)),-real(q(j,k)))
      wp = wp + at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      fx(1,k) = zero
      fy(1,k) = zero
      fx(nx+1,k) = zero
      fy(nx+1,k) = zero
   50 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fx(j,1) = cmplx(-at2,0.)
      fy(j,1) = zero
      wp = wp + at1*real(q(j,1))**2
   60 continue
      fx(1,1) = zero
      fy(1,1) = zero
      fx(nx+1,1) = zero
      fy(nx+1,1) = zero
      we = float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 90 k = 2, nyh
      do 80 j = 2, nx
      at2 = real(ffb(j,k))
      at1 = at2*aimag(ffb(j,k))
      fx(j,k) = at2*q(j,k)
      wp = wp + at1*q(j,k)*conjg(q(j,k))
   80 continue
c mode numbers kx = 0, nx
      fx(1,k) = zero
      fx(nx+1,k) = zero
   90 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      do 100 j = 2, nx
      at2 = real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      at2 = at2*real(q(j,1))
      fx(j,1) = cmplx(at2,0)
      wp = wp + at1*real(q(j,1))**2
  100 continue
      fx(1,1) = zero
      fx(nx+1,1) = zero
      we = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  110 do 130 k = 2, nyh
      do 120 j = 2, nx
      at1 = aimag(ffb(j,k))
      fy(j,k) = at1*q(j,k)
  120 continue
c mode numbers kx = 0, nx
      fy(1,k) = zero
      fy(nx+1,k) = zero
  130 continue
c mode numbers ky = 0, ny/2
      do 140 j = 2, nx
      at1 = aimag(ffb(j,1))
      at2 = at1*real(q(j,1))
      fy(j,1) = cmplx(at2,0.)
  140 continue
      fy(1,1) = zero
      fy(nx+1,1) = zero
      return
      end
      subroutine POISMX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,ny
     1hd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 8*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero, zt1, zt2
      dimension q(nxv,nyv), fxy(2,nxv,nyv)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = conjg(zt1)
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      wp = wp + at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
   50 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*aimag(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(at2,0.)
      fxy(2,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      wp = wp + at1*aimag(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine POISM22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nyeh,nx
     1v,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fxy,we
c approximate flop count is: 19*nxc*nyc + 8*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= ny/2
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero
      dimension q(nxe,nyeh), fxy(2,nxe,nyeh)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      fxy(1,j,k) = -at2*q(j,k)
      fxy(2,j,k) = at3*cmplx(aimag(q(j,k)),-real(q(j,k)))
      wp = wp + at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(1,nx+1,k) = zero
      fxy(2,nx+1,k) = zero
   50 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(-at2,0.)
      fxy(2,j,1) = zero
      wp = wp + at1*real(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,nx+1,1) = zero
      fxy(2,nx+1,1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine POISMX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,ny
     1hd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions.
c Zeros out z component
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fxy,we
c approximate flop count is: 26*nxc*nyc + 8*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero, zt1, zt2
      dimension q(nxv,nyv), fxy(3,nxv,nyv)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = conjg(zt1)
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
   50 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*aimag(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(at2,0.)
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*aimag(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine POISM23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxe,nyeh,nx
     1v,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed dirichlet/periodic boundary conditions.
c Zeros out z component
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign /= 0, input: q,ffb,isign,nx,ny,nxv,nyv,nyhd, output: fxy,we
c approximate flop count is: 19*nxc*nyc + 8*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c equation used is:
c fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0.
c q(j,k) = complex charge density for fourier mode (j-1,k-1)
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c all for fourier mode (j-1,k-1)
c if isign = 0, form factor array is prepared
c if isign is not equal to 0, force/charge is calculated
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= ny/2
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex q, fxy, ffb, zero
      dimension q(nxe,nyeh), fxy(3,nxe,nyeh)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
c calculate force/charge and sum field energy
   30 wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      fxy(1,j,k) = -at2*q(j,k)
      fxy(2,j,k) = at3*cmplx(aimag(q(j,k)),-real(q(j,k)))
      fxy(3,j,k) = zero
      wp = wp + at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      fxy(1,1,k) = zero
      fxy(2,1,k) = zero
      fxy(3,1,k) = zero
      fxy(1,nx+1,k) = zero
      fxy(2,nx+1,k) = zero
      fxy(3,nx+1,k) = zero
   50 continue
      wp = wp + wp
c mode numbers ky = 0, ny/2
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(-at2,0.)
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      wp = wp + at1*real(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,nx+1,1) = zero
      fxy(2,nx+1,1) = zero
      fxy(3,nx+1,1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine DIVFM2(f,df,nx,ny,ndim,nxe,nyeh)
c this subroutine calculates the divergence in fourier space
c with mixed dirichlet/periodic boundary conditions
c intended for calculating the charge density from the electric field
c input: all except df, output: df
c approximate flop count is: 9*nxc*nyc + 4*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the divergence is calculated using the equation:
c df(kx,ky) = sqrt(-1)*(kx*fx(kx,ky)+ky*fy(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for df(kx=pi) = df(ky=pi) = df(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= ny/2
      complex f, df, zero, zt1
      dimension f(ndim,nxe,nyeh), df(nxe,nyeh)
      if (ndim.lt.2) return
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the divergence
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      zt1 = cmplx(-aimag(f(2,j,k)),real(f(2,j,k)))
      df(j,k) = dky*zt1 - dkx*f(1,j,k)
   10 continue
c mode numbers kx = 0, nx
      df(1,k) = dky*cmplx(-aimag(f(2,1,k)),real(f(2,1,k)))
      df(nx+1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      dkx = dnx*float(j - 1)
      df(j,1) = cmplx(-dkx*real(f(1,j,1)),0.)
   30 continue
      df(1,1) = zero
      df(nx+1,1) = zero
      return
      end
      subroutine GRADFM2(df,f,nx,ny,ndim,nxe,nyeh)
c this subroutine calculates the gradient in fourier space
c with mixed dirichlet/periodic boundary conditions
c intended for calculating the electric field from the potential
c input: all except f, output: f
c approximate flop count is: 7*nxc*nyc + 3*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the gradient is calculated using the equations:
c fx(kx,ky) = sqrt(-1)*kx*df(kx,ky)
c fy(kx,ky) = sqrt(-1)*ky*df(kx,ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for fx(kx=pi) = fy(kx=pi) = 0, fx(ky=pi) = fy(ky=pi) = 0,
c and fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c ndim = number of field arrays, must be >= 2
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= ny/2
      complex df, f, zero
      dimension df(nxe,nyeh), f(ndim,nxe,nyeh)
      if (ndim.lt.2) return
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the gradient
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      f(1,j,k) = dkx*df(j,k)
      f(2,j,k) = dky*cmplx(-aimag(df(j,k)),real(df(j,k)))
   10 continue
c mode numbers kx = 0, nx
      f(1,1,k) = zero
      f(2,1,k) = zero
      f(1,nx+1,k) = zero
      f(2,nx+1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      dkx = dnx*float(j - 1)
      f(1,j,1) = cmplx(dkx*real(df(j,1)),0.)
      f(2,j,1) = zero
   30 continue
      f(1,1,1) = zero
      f(2,1,1) = zero
      f(1,nx+1,1) = zero
      f(2,nx+1,1) = zero
      if (ndim.eq.2) return
c handle case of ndim = 3
      do 50 k = 1, nyh
      do 40 j = 1, nx
      f(3,j,k) = zero
   40 continue
      f(3,nx+1,k) = zero
   50 continue
      return
      end
      subroutine CURLFM2(f,g,nx,ny,nxe,nyeh)
c this subroutine calculates the gradient in fourier space
c with mixed dirichlet/periodic boundary conditions
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 16*nxc*nyc + 5*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
      complex f, g, zero, zt1, zt3
      dimension f(3,nxe,nyeh), g(3,nxe,nyeh)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the curl
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      zt1 = cmplx(-aimag(f(3,j,k)),real(f(3,j,k)))
      zt3 = cmplx(-aimag(f(1,j,k)),real(f(1,j,k)))
      g(1,j,k) = dky*zt1
      g(2,j,k) = -dkx*f(3,j,k)
      g(3,j,k) = dkx*f(2,j,k) - dky*zt3
   10 continue
      zt3 = cmplx(-aimag(f(1,1,k)),real(f(1,1,k)))
      g(1,1,k) = zero
      g(2,1,k) = zero
      g(3,1,k) = -dky*zt3
      g(1,nx+1,k) = zero
      g(2,nx+1,k) = zero
      g(3,nx+1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      dkx = dnx*float(j - 1)
      g(1,j,1) = zero
      g(2,j,1) = cmplx(-dkx*real(f(3,j,1)),0.)
      g(3,j,1) = cmplx(dkx*real(f(2,j,1)),0.)
   30 continue
      g(1,1,1) = zero
      g(2,1,1) = zero
      g(3,1,1) = zero
      g(1,nx+1,1) = zero
      g(2,nx+1,1) = zero
      g(3,nx+1,1) = zero
      return
      end
      subroutine CURLFM22(f,g,nx,ny,nxe,nyeh)
c this subroutine calculates the gradient in fourier space
c with mixed dirichlet/periodic boundary conditions
c intended for calculating the magnetic field from the vector potential
c input: all except g, output: g
c approximate flop count is: 9*nxc*nyc + 3*(nxc + nyc)
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c the curl is calculated using the equations:
c gx(kx,ky) = sqrt(-1)*ky*fz(kx,ky)
c gy(kx,ky) = -sqrt(-1)*kx*fz(kx,ky)
c gz(kx,ky) = sqrt(-1)*(kx*fy(kx,ky)-ky*fx(kx,ky))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k,l = fourier mode numbers,
c except for gx(kx=pi) = gy(kx=pi) = 0, gx(ky=pi) = gy(ky=pi) = 0,
c and gx(kx=0,ky=0) = gy(kx=0,ky=0) = 0.
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
      complex f, g, zero, zt3
      dimension f(2,nxe,nyeh), g(nxe,nyeh)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate the curl
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      zt3 = cmplx(-aimag(f(1,j,k)),real(f(1,j,k)))
      g(j,k) = dkx*f(2,j,k) - dky*zt3
   10 continue
      zt3 = cmplx(-aimag(f(1,1,k)),real(f(1,1,k)))
      g(1,k) = -dky*zt3
      g(nx+1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      dkx = dnx*float(j - 1)
      g(j,1) = cmplx(dkx*real(f(2,j,1)),0.)
   30 continue
      g(1,1) = zero
      g(nx+1,1) = zero
      return
      end
      subroutine CUPERPMX2(cu,nx,ny,nxv,nyv)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions.
c input: all, output: cu
c approximate flop count is: 29*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(ky=pi) = cuy(ky=pi) = 0,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx
c nyv = second dimension of current array, must be >= ny
      complex cu, zero, zt1, zt2
      dimension cu(3,nxv,nyv)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
      zt2 = conjg(zt1)
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt2
      cu(2,j,k1) = cu(2,j,k1) + dky*zt2
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nx
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   30 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
      subroutine CUPERPM2(cu,nx,ny,nxe,nyeh)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions.
c input: all, output: cu
c approximate flop count is: 20*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(ky=pi) = cuy(ky=pi) = 0,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxe = first dimension of current array, must be >= nx+1
c nyeh = second dimension of current array, must be >= nyh
      complex cu, zero, zt1
      dimension cu(3,nxe,nyeh)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = cmplx(aimag(cu(2,j,k)),-real(cu(2,j,k)))
      zt1 = at1*(dkx*cu(1,j,k) + dky*zt1)
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*cmplx(-aimag(zt1),real(zt1))
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
      cu(1,nx+1,k) = zero
      cu(2,nx+1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      cu(1,j,1) = zero
      cu(2,j,1) = cmplx(real(cu(2,j,1)),0.0)
   30 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,nx+1,1) = zero
      cu(2,nx+1,1) = zero
      return
      end
      subroutine CUPERPMX22(cu,nx,ny,nxv,nyv)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions.
c input: all, output: cu
c approximate flop count is: 29*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(ky=pi) = cuy(ky=pi) = 0,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx
c nyv = second dimension of current array, must be >= ny
      complex cu, zero, zt1, zt2
      dimension cu(2,nxv,nyv)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
      zt2 = conjg(zt1)
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt2
      cu(2,j,k1) = cu(2,j,k1) + dky*zt2
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   20 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 30 j = 2, nx
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   30 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
      subroutine CUPERPM22(cu,nx,ny,nxe,nyeh)
c this subroutine calculates the transverse current in fourier space
c with mixed dirichlet/periodic boundary conditions.
c input: all, output: cu
c approximate flop count is: 20*nxc*nyc
c and nxc*nyc divides
c where nxc = nx - 1, nyc = ny/2 - 1
c the transverse current is calculated using the equation:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c except for cux(ky=pi) = cuy(ky=pi) = 0,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxe = first dimension of current array, must be >= nx+1
c nyeh = second dimension of current array, must be >= nyh
      complex cu, zero, zt1
      dimension cu(2,nxe,nyeh)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate transverse part of current
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = cmplx(aimag(cu(2,j,k)),-real(cu(2,j,k)))
      zt1 = at1*(dkx*cu(1,j,k) + dky*zt1)
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*cmplx(-aimag(zt1),real(zt1))
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
      cu(1,nx+1,k) = zero
      cu(2,nx+1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      cu(1,j,1) = zero
      cu(2,j,1) = cmplx(real(cu(2,j,1)),0.0)
   30 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,nx+1,1) = zero
      cu(2,nx+1,1) = zero
      return
      end
      subroutine BPOISMX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxv,n
     1yv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 52*nxc*nyc + 19*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 38*nxc*nyc + 13*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,nxv,nyhd, output: bxy
c approximate flop count is: 11*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(ky=pi) = by(ky=pi) = bz(ky=pi).
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
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
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
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
      dimension cu(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      zt3 = at2*zt2 - at3*zt3
      zt2 = -at2*zt1
      zt1 = at3*zt1
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(3,j,k) = zt3
      bxy(1,j,k1) = -conjg(zt1)
      bxy(2,j,k1) = conjg(zt2)
      bxy(3,j,k1) = conjg(zt3)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(at2*aimag(cu(3,j,1)),0.)
      bxy(3,j,1) = cmplx(-at2*aimag(cu(2,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(2,j,1))**2 + aimag(cu(3,j,1))**2)
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
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nx
      at2 = ci2*real(ffb(j,k))
      at1 = at2*aimag(ffb(j,k))
      zt1 = at2*cu(1,j,k)
      zt2 = at2*cu(2,j,k)
      zt3 = at2*cu(3,j,k)
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(3,j,k) = zt3
      bxy(1,j,k1) = conjg(zt1)
      bxy(2,j,k1) = -conjg(zt2)
      bxy(3,j,k1) = -conjg(zt3)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = ci2*real(ffb(1,k))
      at1 = at2*aimag(ffb(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(at2*real(cu(1,j,1)),0.)
      bxy(2,j,1) = cmplx(0.,at2*aimag(cu(2,j,1)))
      bxy(3,j,1) = cmplx(0.,at2*aimag(cu(3,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(2,j,1))**2 + aimag(cu(3,j,1))**2)
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
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      k1 = ny2 - k
      do 140 j = 2, nx
      at1 = aimag(ffb(j,k))
      zt1 = at1*cu(1,j,k)
      zt2 = at1*cu(2,j,k)
      zt3 = at1*cu(3,j,k)
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(3,j,k) = zt3
      bxy(1,j,k1) = conjg(zt1)
      bxy(2,j,k1) = -conjg(zt2)
      bxy(3,j,k1) = -conjg(zt3)
  140 continue
  150 continue
c mode numbers kx = 0, nx
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(at1*real(cu(1,j,1)),0.)
      bxy(2,j,1) = cmplx(0.,at1*aimag(cu(2,j,1)))
      bxy(3,j,1) = cmplx(0.,at1*aimag(cu(3,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
  170 continue
      at1 = aimag(ffb(1,1))
      bxy(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      return
      end
      subroutine BPOISM23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxe,ny
     1eh,nxv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 49*nxc*nyc + 19*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 37*nxc*nyc + 13*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,nxv,nyhd, output: bxy
c approximate flop count is: 6*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(ky=pi) = by(ky=pi) = bz(ky=pi).
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
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt3
      dimension cu(3,nxe,nyeh), bxy(3,nxe,nyeh)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*cu(3,j,k)
      bxy(3,j,k) = at2*cu(2,j,k) - at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, nyh
      at1 = ci2*real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
      bxy(3,nx+1,k) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
   60 continue
c mode numbers ky = 0, ny/2
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(-at2*real(cu(3,j,1)),0.)
      bxy(3,j,1) = cmplx(at2*real(cu(2,j,1)),0.)
      wp = wp + at1*(real(cu(2,j,1))**2 + real(cu(3,j,1))**2)
   70 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      bxy(3,nx+1,1) = zero
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 100 k = 2, nyh
      do 90 j = 2, nx
      at2 = ci2*real(ffb(j,k))
      at1 = at2*aimag(ffb(j,k))
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(2,j,k) = at2*cu(2,j,k)
      bxy(3,j,k) = at2*cu(3,j,k)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, nyh
      at2 = ci2*real(ffb(1,k))
      at1 = at2*aimag(ffb(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
      bxy(3,nx+1,k) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
  110 continue
c mode numbers ky = 0, ny/2
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = at2*cmplx(real(cu(2,j,1)),0.)
      bxy(3,j,1) = at2*cmplx(real(cu(3,j,1)),0.)
      wp = wp + at1*(real(cu(2,j,1))**2 + real(cu(3,j,1))**2)
  120 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      bxy(3,nx+1,1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      do 140 j = 2, nx
      at1 = aimag(ffb(j,k))
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(2,j,k) = at1*cu(2,j,k)
      bxy(3,j,k) = at1*cu(3,j,k)
  140 continue
  150 continue
c mode numbers kx = 0, nx
      do 160 k = 2, nyh
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
      bxy(3,nx+1,k) = zero
  160 continue
c mode numbers ky = 0, ny/2
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = at1*cmplx(real(cu(1,j,1)),0.)
      bxy(2,j,1) = at1*cmplx(real(cu(2,j,1)),0.)
      bxy(3,j,1) = at1*cmplx(real(cu(3,j,1)),0.)
  170 continue
      at1 = aimag(ffb(1,1))
      bxy(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      bxy(3,nx+1,1) = zero
      return
      end
      subroutine BPOISMX22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nx
     1v,nyv,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bz,wm
c approximate flop count is: 63*nxc*nyc + 19*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 38*nxc*nyc + 13*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,nxv,nyhd, output: bxy
c approximate flop count is: 11*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bz(kx=pi) = 0, bz(ky=pi).
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
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
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
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, bz, ffb, zero, zt1, zt2
      dimension cu(2,nxv,nyv), bxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt2 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      zt1 = at2*zt1 - at3*zt2
      bz(j,k) = zt1
      bz(j,k1) = conjg(zt1)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)))
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      zt2 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bz(1,k) = -at3*zt2
      bz(1,k1) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
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
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 2, nx
      at2 = ci2*real(ffb(j,k))
      at1 = at2*aimag(ffb(j,k))
      zt1 = at2*cu(1,j,k)
      zt2 = at2*cu(2,j,k)
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(1,j,k1) = conjg(zt1)
      bxy(2,j,k1) = -conjg(zt2)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)))
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 110 k = 2, nyh
      k1 = ny2 - k
      at2 = ci2*real(ffb(1,k))
      at1 = at2*aimag(ffb(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(at2*real(cu(1,j,1)),0.)
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
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      k1 = ny2 - k
      do 140 j = 2, nx
      at1 = aimag(ffb(j,k))
      zt1 = at1*cu(1,j,k)
      zt2 = at1*cu(2,j,k)
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(1,j,k1) = conjg(zt1)
      bxy(2,j,k1) = -conjg(zt2)
  140 continue
  150 continue
c mode numbers kx = 0, nx
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(at1*real(cu(1,j,1)),0.)
      bxy(2,j,1) = cmplx(0.,at1*aimag(cu(2,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
  170 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      return
      end
      subroutine BPOISM22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxe
     1,nyeh,nxv,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed dirichlet/periodic boundary conditions.
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nxv,nyhd, output: ffb
c for isign = -1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 49*nxc*nyc + 19*(nxc + nyc)
c for isign = 1, input: cu,ffb,isign,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 37*nxc*nyc + 13*(nxc + nyc)
c for isign = 2, input: cu,ffb,isign,nx,ny,nxv,nyhd, output: bxy
c approximate flop count is: 6*nxc*nyc + 2*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c if isign = 0, form factor array is prepared
c if isign < 0, the magnetic field is calculated using the equations:
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky))*
c            s(kx,ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bz(ky=pi).
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
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ax/ay = half-width of particle in x/y direction
c affp = normalization constant = nx*ny/np, where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, bz, ffb, zero, zt3
      dimension cu(2,nxe,nyeh), bxy(2,nxe,nyeh), bz(nxe,nyeh)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 30
c prepare form factor array
      do 20 k = 1, nyh
      dky = dny*float(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nx
      dkx = dnx*float(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.) then
         ffb(j,k) = cmplx(affp,1.)
      else
         ffb(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
   30 if (isign.gt.0) go to 80
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 50 k = 2, nyh
      dky = dny*float(k - 1)
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,k))*aimag(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bz(j,k) = at2*cu(2,j,k) - at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)))
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, nyh
      at1 = ci2*real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bz(1,k) = -at3*zt3
      bz(nx+1,k) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
   60 continue
c mode numbers ky = 0, ny/2
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      bz(j,1) = cmplx(at2*real(cu(2,j,1)),0.)
      wp = wp + at1*real(cu(2,j,1))**2
   70 continue
      bz(1,1) = zero
      bz(nx+1,1) = zero
      wm = float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   80 if (isign.gt.1) go to 130
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 100 k = 2, nyh
      do 90 j = 2, nx
      at2 = ci2*real(ffb(j,k))
      at1 = at2*aimag(ffb(j,k))
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(2,j,k) = at2*cu(2,j,k)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)))
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, nyh
      at2 = ci2*real(ffb(1,k))
      at1 = at2*aimag(ffb(1,k))
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
  110 continue
c mode numbers ky = 0, ny/2
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = at2*cmplx(real(cu(2,j,1)),0.)
      wp = wp + at1*real(cu(2,j,1))**2
  120 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      do 140 j = 2, nx
      at1 = aimag(ffb(j,k))
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(2,j,k) = at1*cu(2,j,k)
  140 continue
  150 continue
c mode numbers kx = 0, nx
      do 160 k = 2, nyh
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = zero
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
  160 continue
c mode numbers ky = 0, ny/2
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = at1*cmplx(real(cu(1,j,1)),0.)
      bxy(2,j,1) = at1*cmplx(real(cu(2,j,1)),0.)
  170 continue
      at1 = aimag(ffb(1,1))
      bxy(1,1,1) = cmplx(at1*real(cu(1,1,1)),0.)
      bxy(2,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      return
      end
      subroutine IBPOISMX23(cu,bxy,ffb,ci,wm,nx,ny,nxv,nyv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with mixed dirichlet/periodic boundary conditions.
c input: cu,ffb,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 63*nxc*nyc + 19*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(ky=pi) = by(ky=pi) = bz(ky=pi).
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
      dimension cu(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      at1 = ci2*real(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffb(j,k))
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = conjg(zt1)
      zt2 = conjg(zt2)
      zt3 = conjg(zt3)
      bxy(1,j,k1) = -at3*zt1
      bxy(2,j,k1) = -at2*zt1
      bxy(3,j,k1) = at2*zt2 - at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   10 continue
   20 continue
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      at1 = at1*aimag(ffb(1,k))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffb(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(at2*aimag(cu(3,j,1)),0.)
      bxy(3,j,1) = cmplx(-at2*aimag(cu(2,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(2,j,1))**2 + aimag(cu(3,j,1))**2)
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
      subroutine IBPOISM23(cu,bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with mixed dirichlet/periodic boundary conditions.
c input: cu,ffb,ci,nx,ny,nxv,nyhd, output: bxy,wm
c approximate flop count is: 49*nxc*nyc + 19*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
c bx(ky=pi) = by(ky=pi) = bz(ky=pi).
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c bxy(i,j,k) = i component of complex magnetic field
c all for fourier mode (j-1,k-1)
c aimag(ffb(j,k)) = finite-size particle shape factor s
c for fourier mode (j-1,k-1)
c real(ffb(j,k)) = potential green's function g
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
c    |cu(kx,ky)*s(kx,ky)|**2), where
c affp = normalization constant = nx*ny/np, where np=number of particles
c this expression is valid only if the current is divergence-free
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt3
      dimension cu(3,nxe,nyeh), bxy(3,nxe,nyeh)
      dimension ffb(nxv,nyhd)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      at1 = ci2*real(ffb(j,k))
      at2 = dnx*float(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffb(j,k))
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*cu(3,j,k)
      bxy(3,j,k) = at2*cu(2,j,k) - at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   10 continue
   20 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 30 k = 2, nyh
      at1 = ci2*real(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      at1 = at1*aimag(ffb(1,k))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
      bxy(3,nx+1,k) = zero
      wp = wp + at1*cu(1,1,k)*conjg(cu(1,1,k))
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffb(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(-at2*real(cu(3,j,1)),0.)
      bxy(3,j,1) = cmplx(at2*real(cu(2,j,1)),0.)
      wp = wp + at1*(real(cu(2,j,1))**2 + real(cu(3,j,1))**2)
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      bxy(3,nx+1,1) = zero
      wm = float(nx*ny)*wp
      return
      end
      subroutine MAXWELMX2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxv,nyv,nyhd
     1)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields
c with mixed dirichlet/periodic boundary conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 146*nxc*nyc + 45*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
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
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffb(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffb(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffb
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nxv,nyv), bxy(3,nxv,nyv), cu(3,nxv,nyv)
      dimension ffb(nxv,nyhd)
      if (ci.le.0.) return
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffb(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffb(j,k))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
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
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
c update electric and  magnetic fields, ky < 0
      exy(1,j,k1) = conjg(zt7)
      exy(2,j,k1) = -conjg(zt8)
      exy(3,j,k1) = -conjg(zt9)
      bxy(1,j,k1) = -conjg(zt4)
      bxy(2,j,k1) = conjg(zt5)
      bxy(3,j,k1) = conjg(zt6)
   10 continue
   20 continue
      ws = ws + ws
      wp = wp + wp
c mode numbers kx = 0, nx
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffb(1,k))
c update magnetic field half time step
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
c update magnetic field half time step and store electric field
      zt3 = cmplx(-aimag(zt7),real(zt7))
      zt6 = zt6 + dth*(dky*zt3)
      ws = ws + anorm*zt7*conjg(zt7)
      wp = wp + anorm*zt6*conjg(zt6)
      exy(1,1,k) = zt7
      exy(2,1,k) = zero
      exy(3,1,k) = zero
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffb(j,1))
c update magnetic field half time step
      at1 = aimag(exy(3,j,1))
      at2 = aimag(exy(2,j,1))
      at5 = real(bxy(2,j,1)) - dth*(dkx*at1)
      at6 = real(bxy(3,j,1)) + dth*(dkx*at2)
c update electric field whole time step
      at8 = at2 - cdt*(dkx*at6) - afdt*aimag(cu(2,j,1))
      at9 = at1 + cdt*(dkx*at5) - afdt*aimag(cu(3,j,1))
c update magnetic field half time step and store electric field
      at5 = at5 - dth*(dkx*at9)
      at6 = at6 + dth*(dkx*at8)
      ws = ws + anorm*(at8*at8 + at9*at9)
      wp = wp + anorm*(at5*at5 + at6*at6)
      exy(1,j,1) = zero
      exy(2,j,1) = cmplx(0.,at8)
      exy(3,j,1) = cmplx(0.,at9)
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(at5,0.)
      bxy(3,j,1) = cmplx(at6,0.)
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
   40 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
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
      subroutine MAXWELM2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxe,nyeh,nxv,
     1nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields
c with mixed dirichlet/periodic boundary conditions.
c input: all, output: wf, wm, exy, bxy
c approximate flop count is: 134*nxc*nyc + 45*(nxc + nyc)
c where nxc = nx - 1, nyc = ny/2 - 1
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
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0.
c and similarly for bx, by, bz.
c cu(i,j,k) = complex current density
c exy(i,j,k) = complex transverse electric field
c bxy(i,j,k) = complex magnetic field
c for component i, all for fourier mode (j-1,k-1)
c real(ffb(1,1)) = affp = normalization constant = nx*ny/np,
c where np=number of particles
c aimag(ffb(j,k)) = finite-size particle shape factor s,
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
c for fourier mode (j-1,k-1)
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffb
      complex zero, zt1, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exy(3,nxe,nyeh), bxy(3,nxe,nyeh), cu(3,nxe,nyeh)
      dimension ffb(nxv,nyhd)
      if (ci.le.0.) return
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      affp = real(ffb(1,1))
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffb(j,k))
c update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*exy(3,j,k))
      zt6 = bxy(3,j,k) - dth*(dkx*exy(2,j,k) - dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) + cdt*(dkx*zt6) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) - cdt*(dkx*zt5 + dky*zt3) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt9)
      zt6 = zt6 - dth*(dkx*zt8 - dky*zt3)
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
   10 continue
   20 continue
      ws = ws + ws
      wp = wp + wp
c mode numbers kx = 0, nx
      do 30 k = 2, nyh
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffb(1,k))
c update magnetic field half time step
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
c update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
c update magnetic field half time step and store electric field
      zt3 = cmplx(-aimag(zt7),real(zt7))
      zt6 = zt6 + dth*(dky*zt3)
      ws = ws + anorm*zt7*conjg(zt7)
      wp = wp + anorm*zt6*conjg(zt6)
      exy(1,1,k) = zt7
      exy(2,1,k) = zero
      exy(3,1,k) = zero
      exy(1,nx+1,k) = zero
      exy(2,nx+1,k) = zero
      exy(3,nx+1,k) = zero
      bxy(1,1,k) = zero
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      bxy(1,nx+1,k) = zero
      bxy(2,nx+1,k) = zero
      bxy(3,nx+1,k) = zero
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffb(j,1))
c update magnetic field half time step
      at1 = real(exy(3,j,1))
      at2 = real(exy(2,j,1))
      at5 = real(bxy(2,j,1)) + dth*(dkx*at1)
      at6 = real(bxy(3,j,1)) - dth*(dkx*at2)
c update electric field whole time step
      at8 = at2 + cdt*(dkx*at6) - afdt*real(cu(2,j,1))
      at9 = at1 - cdt*(dkx*at5) - afdt*real(cu(3,j,1))
c update magnetic field half time step and store electric field
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 - dth*(dkx*at8)
      ws = ws + anorm*(at8*at8 + at9*at9)
      wp = wp + anorm*(at5*at5 + at6*at6)
      exy(1,j,1) = zero
      exy(2,j,1) = cmplx(at8,0.)
      exy(3,j,1) = cmplx(at9,0.)
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(at5,0.)
      bxy(3,j,1) = cmplx(at6,0.)
   40 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      exy(1,nx+1,1) = zero
      exy(2,nx+1,1) = zero
      exy(3,nx+1,1) = zero
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,nx+1,1) = zero
      bxy(2,nx+1,1) = zero
      bxy(3,nx+1,1) = zero
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
      subroutine DMFIELDM2(q1,q,nx,ny,nxv,nyv,nxe,nyeh)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast sine transform in x
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex q1, q
      dimension q1(nxv,nyv), q(nxe,nyeh)
      integer j, k, nx1, nyh
      nx1 = nx + 1
      nyh = ny/2
      do 20 k = 2, nyh
      do 10 j = 1, nx
      q(j,k) = cmplx(-aimag(q1(j,k)),real(q1(j,k)))
   10 continue
   20 continue
      do 30 j = 2, nx
      q(j,1) = cmplx(-aimag(q1(j,1)),-aimag(q1(j,nyh+1)))
   30 continue
      do 40 k = 2, nyh
      q(nx1,k) = cmplx(real(q1(1,ny+2-k)),-aimag(q1(1,ny+2-k)))
   40 continue
      q(1,1) = cmplx(real(q1(1,1)),real(q1(1,nyh+1)))
      q(nx1,1) = cmplx(aimag(q1(1,1)),aimag(q1(1,nyh+1)))
      return
      end
      subroutine FMFIELDM2(q1,q,nx,ny,nxv,nyv,nxe,nyeh)
c this subroutine copies the charge density into a smaller array
c which would have been created by a fast cosine transform in x
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex q1, q
      dimension q1(nxv,nyv), q(nxe,nyeh)
      integer j, k, nx1, nyh
      nx1 = nx + 1
      nyh = ny/2
      do 20 k = 2, nyh
      do 10 j = 1, nx
      q(j,k) = q1(j,k)
   10 continue
   20 continue
      do 30 j = 2, nx
      q(j,1) = cmplx(real(q1(j,1)),real(q1(j,nyh+1)))
   30 continue
      do 40 k = 2, nyh
      q(nx1,k) = cmplx(real(q1(1,ny+2-k)),-aimag(q1(1,ny+2-k)))
   40 continue
      q(1,1) = cmplx(real(q1(1,1)),real(q1(1,nyh+1)))
      q(nx1,1) = cmplx(aimag(q1(1,1)),aimag(q1(1,nyh+1)))
      return
      end
      subroutine CMFIELDM2(cu1,cu,nx,ny,nxv,nyv,nxe,nyeh)
c this subroutine copies the current into a smaller array
c which would have been created by a fast sine/cosine transforms in x
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex cu1, cu
      dimension cu1(3,nxv,nyv), cu(3,nxe,nyeh)
      integer j, k, nx1, nyh
      nx1 = nx + 1
      nyh = ny/2
      do 20 k = 2, nyh
      do 10 j = 1, nx
      cu(1,j,k) = cu1(1,j,k)
      cu(2,j,k) = cmplx(-aimag(cu1(2,j,k)),real(cu1(2,j,k)))
      cu(3,j,k) = cmplx(-aimag(cu1(3,j,k)),real(cu1(3,j,k)))
   10 continue
   20 continue
      do 30 j = 2, nx
      cu(1,j,1) = cmplx(real(cu1(1,j,1)),real(cu1(1,j,nyh+1)))
      cu(2,j,1) = cmplx(-aimag(cu1(2,j,1)),-aimag(cu1(2,j,nyh+1)))
      cu(3,j,1) = cmplx(-aimag(cu1(3,j,1)),-aimag(cu1(3,j,nyh+1)))
   30 continue
      do 40 k = 2, nyh
      cu(1,nx1,k) = cmplx(real(cu1(1,1,ny+2-k)),-aimag(cu1(1,1,ny+2-k)))
      cu(2,nx1,k) = cmplx(real(cu1(2,1,ny+2-k)),-aimag(cu1(2,1,ny+2-k)))
      cu(3,nx1,k) = cmplx(real(cu1(3,1,ny+2-k)),-aimag(cu1(3,1,ny+2-k)))
   40 continue
      cu(1,1,1) = cmplx(real(cu1(1,1,1)),real(cu1(1,1,nyh+1)))
      cu(2,1,1) = cmplx(real(cu1(2,1,1)),real(cu1(2,1,nyh+1)))
      cu(3,1,1) = cmplx(real(cu1(3,1,1)),real(cu1(3,1,nyh+1)))
      cu(1,nx1,1) = cmplx(aimag(cu1(1,1,1)),aimag(cu1(1,1,nyh+1)))
      cu(2,nx1,1) = cmplx(aimag(cu1(2,1,1)),aimag(cu1(2,1,nyh+1)))
      cu(3,nx1,1) = cmplx(aimag(cu1(3,1,1)),aimag(cu1(3,1,nyh+1)))
      return
      end
      subroutine EMFIELDM2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,nyhd
     1)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, nxv, nyv, nxe, nyeh, nyhd
      complex fxy, exy, ffb
      dimension fxy(3,nxv,nyv), exy(3,nxe,nyeh)
      dimension ffb(nxv,nyhd)
      complex zero, zt2, zt3, zt4
      integer j, k, nyh, ny2, k1
      real at1
      nyh = ny/2
      ny2 = ny + 2
      zero = cmplx(0.,0.)
c add the fields
      if (isign.gt.0) then
         do 20 k = 2, nyh
         k1 = ny2 - k
         do 10 j = 2, nx
         at1 = aimag(ffb(j,k))
         zt2 = exy(1,j,k)*at1
         zt3 = exy(2,j,k)*at1
         zt4 = exy(3,j,k)*at1
         zt3 = cmplx(aimag(zt3),-real(zt3))
         zt4 = cmplx(aimag(zt4),-real(zt4))
         fxy(1,j,k) = fxy(1,j,k) + zt2
         fxy(2,j,k) = fxy(2,j,k) + zt3
         fxy(3,j,k) = fxy(3,j,k) + zt4
         fxy(1,j,k1) = fxy(1,j,k1) + conjg(zt2)
         fxy(2,j,k1) = fxy(2,j,k1) - conjg(zt3)
         fxy(3,j,k1) = fxy(3,j,k1) - conjg(zt4)
   10    continue
         at1 = aimag(ffb(1,k))
         zt3 = exy(2,1,k)*at1
         zt4 = exy(3,1,k)*at1
         fxy(1,1,k) = fxy(1,1,k) + exy(1,1,k)*at1
         fxy(2,1,k) = fxy(2,1,k) + cmplx(aimag(zt3),-real(zt3))
         fxy(3,1,k) = fxy(3,1,k) + cmplx(aimag(zt4),-real(zt4))
   20    continue
         k1 = nyh + 1
         do 30 j = 1, nx
         at1 = aimag(ffb(j,1))
         fxy(1,j,1) = fxy(1,j,1) + cmplx(real(exy(1,j,1))*at1,0.0)
         fxy(2,j,1) = fxy(2,j,1) + cmplx(0.0,-real(exy(2,j,1))*at1)
         fxy(3,j,1) = fxy(3,j,1) + cmplx(0.0,-real(exy(3,j,1))*at1)
   30    continue
c copy the magnetic fields
      else if (isign.lt.0) then
         do 50 k = 2, nyh
         k1 = ny2 - k
         do 40 j = 2, nx
         at1 = aimag(ffb(j,k))
         zt2 = exy(1,j,k)*at1
         zt3 = exy(2,j,k)*at1
         zt4 = exy(3,j,k)*at1
         zt2 = cmplx(aimag(zt2),-real(zt2))
         fxy(1,j,k) = zt2
         fxy(2,j,k) = zt3
         fxy(3,j,k) = zt4
         fxy(1,j,k1) = -conjg(zt2)
         fxy(2,j,k1) = conjg(zt3)
         fxy(3,j,k1) = conjg(zt4)
   40    continue
         at1 = aimag(ffb(1,k))
         zt2 = exy(1,1,k)*at1
         fxy(1,1,k) = cmplx(aimag(zt2),-real(zt2))
         fxy(2,1,k) = exy(2,1,k)*at1
         fxy(3,1,k) = exy(3,1,k)*at1
         fxy(1,1,k1) = zero
         fxy(2,1,k1) = zero
         fxy(3,1,k1) = zero
   50    continue
         k1 = nyh + 1
         do 60 j = 1, nx
         at1 = aimag(ffb(j,1))
         fxy(1,j,1) = cmplx(0.0,-real(exy(1,j,1))*at1)
         fxy(2,j,1) = cmplx(real(exy(2,j,1))*at1,0.0)
         fxy(3,j,1) = cmplx(real(exy(3,j,1))*at1,0.0)
         fxy(1,j,k1) = zero
         fxy(2,j,k1) = zero
         fxy(3,j,k1) = zero
   60    continue
c copy the electric fields
      else
         do 80 k = 2, nyh
         k1 = ny2 - k
         do 70 j = 2, nx
         zt2 = exy(1,j,k)
         zt3 = cmplx(aimag(exy(2,j,k)),-real(exy(2,j,k)))
         zt4 = cmplx(aimag(exy(3,j,k)),-real(exy(3,j,k)))
         fxy(1,j,k) = zt2
         fxy(2,j,k) = zt3
         fxy(3,j,k) = zt4
         fxy(1,j,k1) = conjg(zt2)
         fxy(2,j,k1) = -conjg(zt3)
         fxy(3,j,k1) = -conjg(zt4)
   70    continue
         fxy(1,1,k) = exy(1,1,k)
         fxy(2,1,k) = cmplx(aimag(exy(2,1,k)),-real(exy(2,1,k)))
         fxy(3,1,k) = cmplx(aimag(exy(3,1,k)),-real(exy(3,1,k)))
         fxy(1,1,k1) = conjg(exy(1,nx+1,k))
         fxy(2,1,k1) = cmplx(real(exy(2,nx+1,k)),-aimag(exy(2,nx+1,k)))
         fxy(3,1,k1) = cmplx(real(exy(3,nx+1,k)),-aimag(exy(3,nx+1,k)))
   80    continue
         k1 = nyh + 1
         do 90 j = 2, nx
         fxy(1,j,1) = cmplx(real(exy(1,j,1)),0.0)
         fxy(2,j,1) = cmplx(0.0,-real(exy(2,j,1)))
         fxy(3,j,1) = cmplx(0.0,-real(exy(3,j,1)))
         fxy(1,j,k1) = cmplx(aimag(exy(1,j,1)),0.0)
         fxy(2,j,k1) = cmplx(0.0,-aimag(exy(2,j,1)))
         fxy(3,j,k1) = cmplx(0.0,-aimag(exy(3,j,1)))
   90    continue
         fxy(1,1,1) = cmplx(real(exy(1,1,1)),real(exy(1,nx+1,1)))
         fxy(2,1,1) = cmplx(real(exy(2,1,1)),real(exy(2,nx+1,1)))
         fxy(3,1,1) = cmplx(real(exy(3,1,1)),real(exy(3,nx+1,1)))
         fxy(1,1,k1) = cmplx(aimag(exy(1,1,1)),aimag(exy(1,nx+1,1)))
         fxy(2,1,k1) = cmplx(aimag(exy(2,1,1)),aimag(exy(2,nx+1,1)))
         fxy(3,1,k1) = cmplx(aimag(exy(3,1,1)),aimag(exy(3,nx+1,1)))
      endif
      return
      end
      subroutine PMFIELDM2(pot1,pot,nx,ny,nxv,nyv,nxe,nyeh)
c copies image charges appropriate for potential
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex pot1, pot
      dimension pot1(nxv,nyv), pot(nxe,nyeh)
      complex zero, zt1
      integer j, k, nyh, ny2, k1
      nyh = ny/2
      ny2 = ny + 2
      zero = cmplx(0.,0.)
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nx
      zt1 = cmplx(aimag(pot(j,k)),-real(pot(j,k)))
      pot1(j,k) = zt1
      pot1(j,k1) = -conjg(zt1)
   10 continue
      pot1(1,k) = cmplx(aimag(pot(1,k)),-real(pot(1,k)))
      pot1(1,k1) = cmplx(real(pot(nx+1,k)),-aimag(pot(nx+1,k)))
   20 continue
      k1 = nyh + 1
      do 30 j = 2, nx
      pot1(j,1) = cmplx(0.0,-real(pot(j,1)))
      pot1(j,k1) = cmplx(0.0,-aimag(pot(j,1)))
   30 continue
      pot1(1,1) = cmplx(real(pot(1,1)),real(pot(nx+1,1)))
      pot1(1,k1) = cmplx(aimag(pot(1,1)),aimag(pot(nx+1,1)))
      return
      end
      subroutine BMFIELDM2(fxy,bxy,nx,ny,nxv,nyv,nxe,nyeh)
c copies image charges appropriate for magnetic field
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex fxy, bxy
      dimension fxy(3,nxv,nyv), bxy(3,nxe,nyeh)
      complex zt2, zt3, zt4
      integer j, k, nyh, ny2, k1
      nyh = ny/2
      ny2 = ny + 2
c copy the magnetic fields
      do 20 k = 2, nyh
      k1 = ny2 - k
      do 10 j = 2, nx
      zt2 = cmplx(aimag(bxy(1,j,k)),-real(bxy(1,j,k)))
      zt3 = bxy(2,j,k)
      zt4 = bxy(3,j,k)
      fxy(1,j,k) = zt2
      fxy(2,j,k) = zt3
      fxy(3,j,k) = zt4
      fxy(1,j,k1) = -conjg(zt2)
      fxy(2,j,k1) = conjg(zt3)
      fxy(3,j,k1) = conjg(zt4)
   10 continue
      fxy(1,1,k) = cmplx(aimag(bxy(1,1,k)),-real(bxy(1,1,k)))
      fxy(2,1,k) = bxy(2,1,k)
      fxy(3,1,k) = bxy(3,1,k)
      fxy(1,1,k1) = cmplx(real(bxy(1,nx+1,k)),-aimag(bxy(1,nx+1,k)))
      fxy(2,1,k1) = conjg(bxy(2,nx+1,k))
      fxy(3,1,k1) = conjg(bxy(3,nx+1,k))
   20 continue
      k1 = nyh + 1
      do 30 j = 2, nx
      fxy(1,j,1) = cmplx(0.0,-real(bxy(1,j,1)))
      fxy(2,j,1) = cmplx(real(bxy(2,j,1)),0.0)
      fxy(3,j,1) = cmplx(real(bxy(3,j,1)),0.0)
      fxy(1,j,k1) = cmplx(0.0,-aimag(bxy(1,j,1)))
      fxy(2,j,k1) = cmplx(aimag(bxy(2,j,1)),0.0)
      fxy(3,j,k1) = cmplx(aimag(bxy(3,j,1)),0.0)
   30 continue
      fxy(1,1,1) = cmplx(real(bxy(1,1,1)),real(bxy(1,nx+1,1)))
      fxy(2,1,1) = cmplx(real(bxy(2,1,1)),real(bxy(2,nx+1,1)))
      fxy(3,1,1) = cmplx(real(bxy(3,1,1)),real(bxy(3,nx+1,1)))
      fxy(1,1,k1) = cmplx(aimag(bxy(1,1,1)),aimag(bxy(1,nx+1,1)))
      fxy(2,1,k1) = cmplx(aimag(bxy(2,1,1)),aimag(bxy(2,nx+1,1)))
      fxy(3,1,k1) = cmplx(aimag(bxy(3,1,1)),aimag(bxy(3,nx+1,1)))
      return
      end
      subroutine CPFIELDM2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
c this subroutine copies complex vector fields and adds image charges
c appropriate for electric field
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex fxy, exy
      dimension fxy(3,nxv,nyv), exy(3,nxe,nyeh)
c local data
      integer isign, nyhd
      complex ffb
      dimension ffb(1,1)
      isign = 0
      nyhd = 1
      call EMFIELDM2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,nyhd)
      return
      end
      subroutine AVPOTMX23(bxy,axy,nx,ny,nxv,nyv)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with mixed dirichlet/periodic boundary conditions.
c input: bxy, nx, ny, nxv, nyv, output: axy
c approximate flop count is: 26*nxc*nyc + 5*(nxc + nyc)
c and nxc*nyc divides, where nxc = nx - 1, nyc = ny/2 - 1
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) .
c bxy(i,j,k) = i component of complex magnetic field
c axy(i,j,k) = i component of complex vector potential
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
      complex bxy, axy, zero, zt1, zt2, zt3
      dimension bxy(3,nxv,nyv), axy(3,nxv,nyv)
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate vector potential
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      zt3 = at2*zt2 - at3*zt3
      zt2 = -at2*zt1
      zt1 = at3*zt1
      axy(1,j,k) = zt1
      axy(2,j,k) = zt2
      axy(3,j,k) = zt3
      axy(1,j,k1) = conjg(zt1)
      axy(2,j,k1) = -conjg(zt2)
      axy(3,j,k1) = -conjg(zt3)
   10 continue
   20 continue
c mode numbers kx = 0, nx
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      zt1 = cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      axy(1,1,k) = at3*zt1
      axy(2,1,k) = zero
      axy(3,1,k) = zero
      axy(1,1,k1) = zero
      axy(2,1,k1) = zero
      axy(3,1,k1) = zero
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      at2 = 1.0/dkx
      axy(1,j,1) = zero
      axy(2,j,1) = cmplx(0.,-at2*real(bxy(3,j,1)))
      axy(3,j,1) = cmplx(0.,at2*real(bxy(2,j,1)))
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
      subroutine AVPOTM23(bxy,axy,nx,ny,nxe,nyeh)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with mixed dirichlet/periodic boundary conditions.
c input: bxy, nx, ny, nxv, nyv, output: axy
c approximate flop count is: 21*nxc*nyc + 5*(nxc + nyc)
c and nxc*nyc divides, where nxc = nx - 1, nyc = ny/2 - 1
c the vector potential is calculated using the equations:
c ax(kx,ky) = sqrt(-1)*(ky*bz(kx,ky))/(kx*kx+ky*ky)
c ay(kx,ky) = -sqrt(-1)*(kx*bz(kx,ky))/(kx*kx+ky*ky)
c az(kx,ky) = sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))/(kx*kx+ky*ky),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
c ax(ky=pi) = ay(ky=pi) = az(ky=pi) .
c bxy(i,j,k) = i component of complex magnetic field
c axy(i,j,k) = i component of complex vector potential
c all for fourier mode (j-1,k-1)
c nx/ny = system length in x/y direction
c nxe = first dimension of field arrays, must be >= nx+1
c nyeh = second dimension of field arrays, must be >= nyh
      complex bxy, axy, zero
      dimension bxy(3,nxe,nyeh), axy(3,nxe,nyeh)
      nyh = ny/2
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(ny)
      zero = cmplx(0.,0.)
c calculate vector potential
c mode numbers 0 < kx < nx and 0 < ky < ny/2
      do 20 k = 2, nyh
      dky = dny*float(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nx
      dkx = dnx*float(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      at2 = dkx*at1
      at3 = dky*at1
      axy(1,j,k) = at3*cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      axy(2,j,k) = at2*bxy(3,j,k)
      axy(3,j,k) = at3*cmplx(aimag(bxy(1,j,k)),-real(bxy(1,j,k))) - at2*
     1bxy(2,j,k)
   10 continue
   20 continue
c mode numbers kx = 0, nx
      do 30 k = 2, nyh
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      axy(1,1,k) = at3*cmplx(-aimag(bxy(3,1,k)),real(bxy(3,1,k)))
      axy(2,1,k) = zero
      axy(3,1,k) = zero
      axy(1,nx+1,k) = zero
      axy(2,nx+1,k) = zero
      axy(3,nx+1,k) = zero
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      at2 = 1.0/dkx
      axy(1,j,1) = zero
      axy(2,j,1) = cmplx(at2*real(bxy(3,j,1)),0.0)
      axy(3,j,1) = cmplx(-at2*real(bxy(2,j,1)),0.0)
   40 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      axy(1,nx+1,1) = zero
      axy(2,nx+1,1) = zero
      axy(3,nx+1,1) = zero
      return
      end
