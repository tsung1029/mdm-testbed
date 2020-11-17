c 2d PIC library for solving field equations with neumann
c boundary conditions
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: april 10, 2006
      subroutine DBLCOS2C(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in x,
c and y component is an odd function in y.
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx
c nyv = third dimension of input array cu, must be >= ny
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
      cu2(1,nx+j+1,k+1) = -cu(1,nx-j+1,k+1)
      cu2(2,nx+j+1,k+1) = cu(2,nx-j+1,k+1)
      cu2(1,j+1,ny+k+1) = cu(1,j+1,ny-k+1)
      cu2(2,j+1,ny+k+1) = -cu(2,j+1,ny-k+1)
      cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
      cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
   10 continue
      cu2(1,1,k+1) = 0.
      cu2(2,1,k+1) = cu(2,1,k+1)
      cu2(1,nx+1,k+1) = 0.
      cu2(2,nx+1,k+1) = cu(2,nx+1,k+1)
      cu2(1,1,k+ny+1) = 0.
      cu2(2,1,k+ny+1) = -cu(2,1,ny-k+1)
      cu2(1,nx+1,k+ny+1) = 0.
      cu2(2,nx+1,k+ny+1) = -cu(2,nx+1,ny-k+1)
   20 continue
      do 30 j = 1, nxs
      cu2(1,j+1,1) = cu(1,j+1,1)
      cu2(2,j+1,1) = 0.
      cu2(1,j+nx+1,1) = -cu(1,nx-j+1,1)
      cu2(2,j+nx+1,1) = 0.
      cu2(1,j+1,ny+1) = cu(1,j+1,ny+1)
      cu2(2,j+1,ny+1) = 0.
      cu2(1,j+nx+1,ny+1) = -cu(1,nx-j+1,ny+1)
      cu2(2,j+nx+1,ny+1) = 0.
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
      subroutine DBLCOS2D(q,q2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates an even array q2 from an array q, so that
c a 2d cosine transform can be performed with a 2d real to complex fft.
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
      q2(nx+j+1,k+1) = q(nx-j+1,k+1)
      q2(j+1,ny+k+1) = q(j+1,ny-k+1)
      q2(nx+j+1,ny+k+1) = q(nx-j+1,ny-k+1)
   10 continue
      q2(1,k+1) = q(1,k+1)
      q2(nx+1,k+1) = q(nx+1,k+1)
      q2(1,k+ny+1) = q(1,ny-k+1)
      q2(nx+1,k+ny+1) = q(nx+1,ny-k+1)
   20 continue
      do 30 j = 1, nxs
      q2(j+1,1) = q(j+1,1)
      q2(nx+j+1,1) = q(nx-j+1,1)
      q2(j+1,ny+1) = q(j+1,ny+1)
      q2(nx+j+1,ny+1) = q(nx-j+1,ny+1)
   30 continue
      q2(1,1) = q(1,1)
      q2(nx+1,1) = q(nx+1,1)
      q2(1,ny+1) = q(1,ny+1)
      q2(nx+1,ny+1) = q(nx+1,ny+1)
      return
      end
      subroutine DBLCOS2B(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in x,
c y component is an odd function in y, and the z component is an even
c function in both x and y.
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx
c nyv = third dimension of input array cu, must be >= ny
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
      cu2(1,nx+j+1,k+1) = -cu(1,nx-j+1,k+1)
      cu2(2,nx+j+1,k+1) = cu(2,nx-j+1,k+1)
      cu2(3,nx+j+1,k+1) = cu(3,nx-j+1,k+1)
      cu2(1,j+1,ny+k+1) = cu(1,j+1,ny-k+1)
      cu2(2,j+1,ny+k+1) = -cu(2,j+1,ny-k+1)
      cu2(3,j+1,ny+k+1) = cu(3,j+1,ny-k+1)
      cu2(1,nx+j+1,ny+k+1) = -cu(1,nx-j+1,ny-k+1)
      cu2(2,nx+j+1,ny+k+1) = -cu(2,nx-j+1,ny-k+1)
      cu2(3,nx+j+1,ny+k+1) = cu(3,nx-j+1,ny-k+1)
   10 continue
      cu2(1,1,k+1) = 0.
      cu2(2,1,k+1) = cu(2,1,k+1)
      cu2(3,1,k+1) = cu(3,1,k+1)
      cu2(1,nx+1,k+1) = 0.
      cu2(2,nx+1,k+1) = cu(2,nx+1,k+1)
      cu2(3,nx+1,k+1) = cu(3,nx+1,k+1)
      cu2(1,1,k+ny+1) = 0.
      cu2(2,1,k+ny+1) = -cu(2,1,ny-k+1)
      cu2(3,1,k+ny+1) = cu(3,1,ny-k+1)
      cu2(1,nx+1,k+ny+1) = 0.
      cu2(2,nx+1,k+ny+1) = -cu(2,nx+1,ny-k+1)
      cu2(3,nx+1,k+ny+1) = cu(3,nx+1,ny-k+1)
   20 continue
      do 30 j = 1, nxs
      cu2(1,j+1,1) = cu(1,j+1,1)
      cu2(2,j+1,1) = 0.
      cu2(3,j+1,1) = cu(3,j+1,1)
      cu2(1,j+nx+1,1) = -cu(1,nx-j+1,1)
      cu2(2,j+nx+1,1) = 0.
      cu2(3,j+nx+1,1) = cu(3,nx-j+1,1)
      cu2(1,j+1,ny+1) = cu(1,j+1,ny+1)
      cu2(2,j+1,ny+1) = 0.
      cu2(3,j+1,ny+1) = cu(3,j+1,ny+1)
      cu2(1,j+nx+1,ny+1) = -cu(1,nx-j+1,ny+1)
      cu2(2,j+nx+1,ny+1) = 0.
      cu2(3,j+nx+1,ny+1) = cu(3,nx-j+1,ny+1)
   30 continue
      cu2(1,1,1) = 0.
      cu2(2,1,1) = 0.
      cu2(3,1,1) = cu(3,1,1)
      cu2(1,nx+1,1) = 0.
      cu2(2,nx+1,1) = 0.
      cu2(3,nx+1,1) = cu(3,nx+1,1)
      cu2(1,1,ny+1) = 0.
      cu2(2,1,ny+1) = 0.
      cu2(3,1,ny+1) = cu(3,1,ny+1)
      cu2(1,nx+1,ny+1) = 0.
      cu2(2,nx+1,ny+1) = 0.
      cu2(3,nx+1,ny+1) = cu(3,nx+1,ny+1)
      return
      end
      subroutine SGLVCS2C(cu,cu1,nx,ny,ncsx,ncsy,nxv,nyv,nx2v)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c ncsx/ncsy = (0,1) = (sine,cosine) transform in x/y direction
c nxv = second dimension of input array cu, must be >= nx
c nyv = third dimension of input array cu, must be >= ny
c nx2v = second dimension of output array cu1, must be >= 2*nx
      implicit none
      real cu, cu1
      integer nx, ny, ncsx, ncsy, nxv, nyv, nx2v
      dimension cu(2,nxv,nyv), cu1(2,nx2v,nyv)
c local data
      integer j, k, nxs
      real csx, csy
c copy to double array
      nxs = nx - 1
      csx = 2*ncsx - 1
      csy = 2*ncsy - 1
      do 20 k = 1, ny
      do 10 j = 1, nxs
      cu1(1,j+1,k) = cu(1,j+1,k)
      cu1(2,j+1,k) = cu(2,j+1,k)
      cu1(1,nx+j+1,k) = csx*cu(1,nx-j+1,k)
      cu1(2,nx+j+1,k) = csy*cu(2,nx-j+1,k)
   10 continue
      cu1(1,1,k) = 0.
      cu1(2,1,k) = 0.
      cu1(1,nx+1,k) = 0.
      cu1(2,nx+1,k) = 0.
   20 continue
      return
      end
      subroutine SGLVCS2D(q,q1,nx,ny,ncs,nxv,nyv,nx2v)
c this subroutine creates an even array q1 from an array q, so that
c a 1d sine transform can be performed with a 2d real to complex fft.
c linear interpolation
c nx/ny = system length in x/y direction
c ncs = (0,1) = (sine,cosine) transform
c nxv = first dimension of input array q, must be >= nx
c nyv = second dimension of input array q, must be >= ny
c nx2v = first dimension of output array q1, must be >= 2*nx
      implicit none
      real q, q1
      integer nx, ny, ncs, nxv, nyv, nx2v
      dimension q(nxv,nyv), q1(nx2v,nyv)
c local data
      integer j, k, nxs
      real cs
c copy to double array
      nxs = nx - 1
      cs = 2*ncs - 1
      do 20 k = 1, ny
      do 10 j = 1, nxs
      q1(j+1,k) = q(j+1,k)
      q1(nx+j+1,k) = cs*q(nx-j+1,k)
   10 continue
      q1(1,k) = 0.
      q1(nx+1,k) = 0.
   20 continue
      return
      end
      subroutine SGLVCS2B(cu,cu1,nx,ny,ncsx,ncsy,ncsz,nxv,nyv,nx2v)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c ncsx/ncsy/ncsz = (0,1) = (sine,cosine) transform in x/y/z direction
c nxv = second dimension of input array cu, must be >= nx
c nyv = third dimension of input array cu, must be >= ny
c nx2v = second dimension of output array cu1, must be >= 2*nx
      implicit none
      real cu, cu1
      integer nx, ny, ncsx, ncsy, ncsz, nxv, nyv, nx2v
      dimension cu(3,nxv,nyv), cu1(3,nx2v,nyv)
c local data
      integer j, k, nxs
      real csx, csy, csz
c copy to double array
      nxs = nx - 1
      csx = 2*ncsx - 1
      csy = 2*ncsy - 1
      csz = 2*ncsz - 1
      do 20 k = 1, ny
      do 10 j = 1, nxs
      cu1(1,j+1,k) = cu(1,j+1,k)
      cu1(2,j+1,k) = cu(2,j+1,k)
      cu1(3,j+1,k) = cu(3,j+1,k)
      cu1(1,nx+j+1,k) = csx*cu(1,nx-j+1,k)
      cu1(2,nx+j+1,k) = csy*cu(2,nx-j+1,k)
      cu1(3,nx+j+1,k) = csz*cu(3,nx-j+1,k)
   10 continue
      cu1(1,1,k) = 0.
      cu1(2,1,k) = 0.
      cu1(3,1,k) = 0.
      cu1(1,nx+1,k) = 0.
      cu1(2,nx+1,k) = 0.
      cu1(3,nx+1,k) = 0.
   20 continue
      return
      end
      subroutine POISNX2(q,fx,fy,isign,ffd,ax,ay,affp,we,nx,ny,nx2v,ny2d
     1)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with neumann boundary conditions (zero normal efield).
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
c s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
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
         ffd(2*j,k) = 1.
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
      fx(2*j-1,k) = 0.
      fx(2*j,k) = at2
      fx(2*j-1,k1) = 0.
      fx(2*j,k1) = at2
      fy(2*j-1,k) = 0.
      fy(2*j,k) = at3
      fy(2*j-1,k1) = 0.
      fy(2*j,k1) = -at3
      wp = wp + at1*q(2*j-1,k)**2
   40 continue
      at1 = ffd(1,k)*ffd(2,k)
      at3 = -dky*at1*q(1,k)
      fx(1,k) = 0.
      fx(2,k) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      fy(1,k) = 0.
      fy(2,k) = at3
      fy(1,k1) = 0.
      fy(2,k1) = 0.
      wp = wp + .5*at1*q(1,k)**2
   50 continue
      do 60 j = 2, nx
      at1 = ffd(2*j-1,1)*ffd(2*j,1)
      at3 = -at1*q(2*j-1,1)
      at2 = dnx*float(j - 1)*at3
      fx(2*j-1,1) = 0.
      fx(2*j,1) = at2
      fx(2*j-1,ny+1) = 0.
      fx(2*j,ny+1) = 0.
      fy(2*j-1,1) = 0.
      fy(2*j,1) = 0.
      fy(2*j-1,ny+1) = 0.
      fy(2*j,ny+1) = 0.
      wp = wp + .5*at1*q(2*j-1,1)**2
   60 continue
      fx(1,1) = 0.
      fx(2,1) = 0.
      fx(1,ny+1) = 0.
      fx(2,ny+1) = 0.
      fy(1,1) = 0.
      fy(2,1) = 0.
      fy(1,ny+1) = 0.
      fy(2,ny+1) = 0.
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
      fx(2*j,k) = 0.
      fx(2*j-1,k1) = at3
      fx(2*j,k1) = 0.
      wp = wp + at1*q(2*j-1,k)**2
   80 continue
      at2 = ffd(1,k)
      at1 = at2*ffd(2,k)
      at3 = at2*q(1,k)
      fx(1,k) = at3
      fx(2,k) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      wp = wp + .5*at1*q(1,k)**2
   90 continue
      do 100 j = 2, nx
      at2 = ffd(2*j-1,1)
      at1 = at2*ffd(2*j,1)
      at3 = at2*q(2*j-1,1)
      fx(2*j-1,1) = at3
      fx(2*j,1) = 0.
      fx(2*j-1,ny+1) = 0.
      fx(2*j,ny+1) = 0.
      wp = wp + .5*at1*q(2*j-1,1)**2
  100 continue
      fx(1,1) = 0.
      fx(2,1) = 0.
      fx(1,ny+1) = 0.
      fx(2,ny+1) = 0.
      we = 2.0*float(nx*ny)*wp
      return
c calculate smoothing
  110 do 130 k = 2, ny
      k1 = ny2 - k
      do 120 j = 2, nx
      at1 = ffd(2*j,k)
      at2 = at1*q(2*j-1,k)
      fy(2*j-1,k) = at2
      fy(2*j,k) = 0.
      fy(2*j-1,k1) = at2
      fy(2*j,k1) = 0
  120 continue
      at1 = ffd(2,k)
      at2 = at1*q(1,k)
      fy(1,k) = at2
      fy(2,k) = 0.
      fy(1,k1) = 0.
      fy(2,k1) = 0
  130 continue
      do 140 j = 2, nx
      at1 = ffd(2*j,1)
      at2 = at1*q(2*j-1,1)
      fy(2*j-1,1) = at2
      fy(2*j,1) = 0.
      fy(2*j-1,ny+1) = 0.
      fy(2*j,ny+1) = 0.
  140 continue
      at1 = ffd(2,1)
      at2 = at1*q(1,1)
      fy(1,1) = at2
      fy(2,1) = 0.
      fy(1,ny+1) = 0.
      fy(2,ny+1) = 0.
      return
      end
      subroutine POISNX22(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with neumann boundary conditions (zero normal efield).
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nx2v,ny2d, output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fxy,we
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
      fxy(1,j,k) = cmplx(0.,at2)
      fxy(2,j,k) = cmplx(0.,at3)
      fxy(1,j,k1) = cmplx(0.,at2)
      fxy(2,j,k1) = cmplx(0.,-at3)
      wp = wp + at1*real(q(j,k))**2
   40 continue
      at1 = real(ffd(1,k))*aimag(ffd(1,k))
      at3 = -dky*at1*real(q(1,k))
      fxy(1,1,k) = zero
      fxy(2,1,k) = cmplx(0.,at3)
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      wp = wp + .5*at1*real(q(1,k))**2
   50 continue
      do 60 j = 2, nx
      at1 = real(ffd(j,1))*aimag(ffd(j,1))
      at3 = -at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(0.,at2)
      fxy(2,j,1) = zero
      fxy(1,j,ny+1) = zero
      fxy(2,j,ny+1) = zero
      wp = wp + .5*at1*real(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,1,ny+1) = zero
      fxy(2,1,ny+1) = zero
      we = 2.0*float(nx*ny)*wp
      return
      end
      subroutine POISNX23(q,fxy,isign,ffd,ax,ay,affp,we,nx,ny,nxv,ny2d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with neumann boundary conditions (zero normal efield).
c Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,affp,nx,ny,nx2v,ny2d, output: ffd
c for isign /= 0, input: q,ffd,isign,nx,ny,nx2v,ny2d, output: fxy,we
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
      fxy(1,j,k) = cmplx(0.,at2)
      fxy(2,j,k) = cmplx(0.,at3)
      fxy(3,j,k) = zero
      fxy(1,j,k1) = cmplx(0.,at2)
      fxy(2,j,k1) = cmplx(0.,-at3)
      fxy(3,j,k1) = zero
      wp = wp + at1*real(q(j,k))**2
   40 continue
      at1 = real(ffd(1,k))*aimag(ffd(1,k))
      at3 = -dky*at1*real(q(1,k))
      fxy(1,1,k) = zero
      fxy(2,1,k) = cmplx(0.,at3)
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      wp = wp + .5*at1*real(q(1,k))**2
   50 continue
      do 60 j = 2, nx
      at1 = real(ffd(j,1))*aimag(ffd(j,1))
      at3 = -at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(0.,at2)
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,ny+1) = zero
      fxy(2,j,ny+1) = zero
      fxy(3,j,ny+1) = zero
      wp = wp + .5*at1*real(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,ny+1) = zero
      fxy(2,1,ny+1) = zero
      fxy(3,1,ny+1) = zero
      we = 2.0*float(nx*ny)*wp
      return
      end
      subroutine CUPERPNX2(cu,nx,ny,nxv,ny2d)
c this subroutine calculates the transverse current in fourier space
c with neumann boundary conditions (zero normal efield).
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
      cu(1,j,k1) = cmplx(0.,at3)
      cu(2,j,k1) = cmplx(0.,-at4)
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
      subroutine CUPERPN2(cu,nx,ny,nxe,nye)
c this subroutine calculates the transverse current in fourier space
c with neumann boundary conditions (zero normal efield).
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
      subroutine CUPERPNX22(cu,nx,ny,nxv,ny2d)
c this subroutine calculates the transverse current in fourier space
c with neumann boundary conditions (zero normal efield).
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
      cu(1,j,k1) = cmplx(0.,at3)
      cu(2,j,k1) = cmplx(0.,-at4)
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
      subroutine BPOISNX23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxv,n
     1y2d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with neumann boundary conditions (zero normal efield).
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
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
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
      bxy(1,j,k1) = -bxy(1,j,k)
      bxy(2,j,k1) = bxy(2,j,k)
      bxy(3,j,k1) = -bxy(3,j,k)
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
      bxy(1,1,k) = cmplx(0.,at3*real(cu(3,1,k)))
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
      bxy(2,j,1) = cmplx(0.,-at2*real(cu(3,j,1)))
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
      bxy(1,j,k1) = bxy(1,j,k)
      bxy(2,j,k1) = -bxy(2,j,k)
      bxy(3,j,k1) = bxy(3,j,k)
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
      bxy(3,1,k) = cmplx(at2*real(cu(3,1,k)),0.)
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
      bxy(3,j,1) = cmplx(at2*real(cu(3,j,1)),0.)
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
      bxy(1,j,k1) = bxy(1,j,k)
      bxy(2,j,k1) = -bxy(2,j,k)
      bxy(3,j,k1) = bxy(3,j,k)
  140 continue
c mode numbers kx = 0, nx
      at1 = aimag(ffd(1,k))
      bxy(1,1,k) = cmplx(0.,at1*aimag(cu(1,1,k)))
      bxy(2,1,k) = zero
      bxy(3,1,k) = cmplx(at1*real(cu(3,1,k)),0.)
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
      bxy(3,j,1) = cmplx(at1*real(cu(3,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
  160 continue
      at1 = aimag(ffd(1,1))
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      return
      end
      subroutine BPOISN23(cu,bxy,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nxe,ny
     1e,nxv)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with neumann boundary conditions (zero normal efield).
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
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky)*s(kx,ky),
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
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
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
      bxy(3,j,k) = at3*cu(1,j,k) - at2*cu(2,j,k)
      wp = wp + at1*(cu(1,j,k)**2 + cu(2,j,k)**2 + cu(3,j,k)**2)
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, ny
      dky = dny*float(k - 1)
      at1 = ci2*real(ffd(1,k))*aimag(ffd(1,k))
      at3 = dky*at1
      bxy(1,1,k) = at3*cu(3,1,k)
      bxy(2,1,k) = 0.
      bxy(3,1,k) = at3*cu(1,1,k)
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
      bxy(2,j,1) = -at2*cu(3,j,1)
      bxy(3,j,1) = -at2*cu(2,j,1)
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
      bxy(3,1,k) = at2*cu(3,1,k)
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
      bxy(3,j,1) = at2*cu(3,j,1)
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
      bxy(3,1,k) = at1*cu(3,1,k)
      bxy(1,nx+1,k) = 0.
      bxy(2,nx+1,k) = 0.
      bxy(3,nx+1,k) = 0.
  150 continue
c mode numbers ky = 0, ny
      do 160 j = 2, nx
      at1 = aimag(ffd(j,1))
      bxy(1,j,1) = 0.
      bxy(2,j,1) = at1*cu(2,j,1)
      bxy(3,j,1) = at1*cu(3,j,1)
      bxy(1,j,ny+1) = 0.
      bxy(2,j,ny+1) = 0.
      bxy(3,j,ny+1) = 0.
  160 continue
      at1 = aimag(ffd(1,1))
      bxy(1,1,1) = 0.
      bxy(2,1,1) = 0.
      bxy(3,1,1) = at1*cu(3,1,1)
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
      subroutine BPOISNX22(cu,bxy,bz,isign,ffd,ax,ay,affp,ci,wm,nx,ny,nx
     1v,ny2d)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with neumann boundary conditions (zero normal efield).
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
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
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
      bz(j,k1) = -bz(j,k)
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
      bxy(1,j,k1) = bxy(1,j,k)
      bxy(2,j,k1) = -bxy(2,j,k)
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
      bxy(1,j,k1) = bxy(1,j,k)
      bxy(2,j,k1) = -bxy(2,j,k)
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
      subroutine IBPOISNX23(cu,bxy,ffd,ci,wm,nx,ny,nxv,ny2d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with neumann boundary conditions (zero normal efield)
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c input: cu,ffd,ci,nx,ny,nxv,ny2d, output: bxy,wm
c approximate flop count is: 20*nx*ny + 11*(nx + ny)
c the magnetic field is calculated using the equations:
c bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
c by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*-kx*cuz(kx,ky),
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
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci, where
c    |cu(kx,ky)*s(kx,ky)|**2), where
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
      bxy(1,j,k1) = -bxy(1,j,k)
      bxy(2,j,k1) = bxy(2,j,k)
      bxy(3,j,k1) = -bxy(3,j,k)
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
      bxy(1,1,k) = cmplx(0.,at3*real(cu(3,1,k)))
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
      bxy(2,j,1) = cmplx(0.,-at2*real(cu(3,j,1)))
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
      subroutine IBPOISN23(cu,bxy,ffd,ci,wm,nx,ny,nxe,nye,nxv)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with neumann boundary conditions (zero normal efield)
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
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
c wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci
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
      bxy(3,j,k) = at3*cu(1,j,k) - at2*cu(2,j,k)
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
      bxy(1,1,k) = at3*cu(3,1,k)
      bxy(2,1,k) = 0.
      bxy(3,1,k) = at3*cu(1,1,k)
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
      bxy(2,j,1) = -at2*cu(3,j,1)
      bxy(3,j,1) = -at2*cu(2,j,1)
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
      subroutine MAXWELNX2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxv,ny2d)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with neumann boundary
c conditions (zero normal efield).
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
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
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
      exy(1,j,k1) = cmplx(0.,at7)
      exy(2,j,k1) = cmplx(0.,-at8)
      exy(3,j,k1) = cmplx(at9,0.)
      bxy(1,j,k1) = cmplx(0.,-at4)
      bxy(2,j,k1) = cmplx(0.,at5)
      bxy(3,j,k1) = cmplx(-at6,0.)
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
      subroutine MAXWELN2(exy,bxy,cu,ffd,ci,dt,wf,wm,nx,ny,nxe,nye,nxv)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields with neumann boundary
c conditions (zero normal efield).
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
c wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
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
      at6 = bxy(3,j,k) + dth*(dkx*at8 - dky*at7)
c update electric field whole time step
      at7 = at7 + cdt*(dky*at6) - afdt*cu(1,j,k)
      at8 = at8 - cdt*(dkx*at6) - afdt*cu(2,j,k)
      at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,j,k)
c update magnetic field half time step and store electric field
      at4 = at4 - dth*(dky*at9)
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 + dth*(dkx*at8 - dky*at7)
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
      at6 = bxy(3,1,k) - dth*(dky*at7)
c update electric field whole time step
      at7 = at7 + cdt*(dky*at6) - afdt*cu(1,1,k)
      at9 = at9 + cdt*(dky*at4) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      at4 = at4 - dth*(dky*at9)
      at6 = at6 - dth*(dky*at7)
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
      at6 = bxy(3,j,1) + dth*(dkx*at8)
c update electric field whole time step
      at8 = at8 - cdt*(dkx*at6) - afdt*cu(2,j,1)
      at9 = at9 - cdt*(dkx*at5) - afdt*cu(3,j,1)
c update magnetic field half time step and store electric field
      at5 = at5 + dth*(dkx*at9)
      at6 = at6 + dth*(dkx*at8)
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
      subroutine CMFIELDN2(cu2,cu,nx,ny,nxv,ny2d,nxe,nye)
c this subroutine copies the current into a smaller array
      implicit none
      integer nx, ny, nxv, ny2d, nxe, nye
      complex cu2
      real cu
      dimension cu2(3,nxv,ny2d), cu(3,nxe,nye)
      integer j, k
      do 20 k = 1, ny
      do 10 j = 1, nx
      cu(1,j,k) = aimag(cu2(1,j,k))
      cu(2,j,k) = aimag(cu2(2,j,k))
      cu(3,j,k) = real(cu2(3,j,k))
   10 continue
   20 continue
      return
      end
      subroutine EMFIELDN2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
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
         fxy(1,j,k) = fxy(1,j,k) + cmplx(0.,exy(1,j,k)*at1)
         fxy(2,j,k) = fxy(2,j,k) + cmplx(0.,exy(2,j,k)*at1)
         fxy(3,j,k) = fxy(3,j,k) + cmplx(exy(3,j,k)*at1,0.)
         fxy(1,j,k1) = fxy(1,j,k1) + cmplx(0.,exy(1,j,k)*at1)
         fxy(2,j,k1) = fxy(2,j,k1) + cmplx(0.,-exy(2,j,k)*at1)
         fxy(3,j,k1) = fxy(3,j,k1) + cmplx(exy(3,j,k)*at1,0.)
   10    continue
         at1 = aimag(ffd(1,k))
         fxy(1,1,k) = fxy(1,1,k) + cmplx(0.,exy(1,1,k)*at1)
         fxy(2,1,k) = fxy(2,1,k) + cmplx(0.,exy(2,1,k)*at1)
         fxy(3,1,k) = fxy(3,1,k) + cmplx(exy(3,1,k)*at1,0.)
   20    continue
         do 30 j = 1, nx
         at1 = aimag(ffd(j,1))
         fxy(1,j,1) = fxy(1,j,1) + cmplx(0.,exy(1,j,1)*at1)
         fxy(2,j,1) = fxy(2,j,1) + cmplx(0.,exy(2,j,1)*at1)
         fxy(3,j,1) = fxy(3,j,1) + cmplx(exy(3,j,1)*at1,0.)
   30    continue
c copy the magnetic fields
      else if (isign.lt.0) then
         do 50 k = 2, ny
         k1 = ny2 - k
         do 40 j = 2, nx
         at1 = aimag(ffd(j,k))
         fxy(1,j,k) = cmplx(0.,exy(1,j,k)*at1)
         fxy(2,j,k) = cmplx(0.,exy(2,j,k)*at1)
         fxy(3,j,k) = cmplx(exy(3,j,k)*at1,0.)
         fxy(1,j,k1) = cmplx(0.,-exy(1,j,k)*at1)
         fxy(2,j,k1) = cmplx(0.,exy(2,j,k)*at1)
         fxy(3,j,k1) = cmplx(-exy(3,j,k)*at1,0.)
   40    continue
         at1 = aimag(ffd(1,k))
         fxy(1,1,k) = cmplx(0.,exy(1,1,k)*at1)
         fxy(2,1,k) = cmplx(0.,exy(2,1,k)*at1)
         fxy(3,1,k) = cmplx(exy(3,1,k)*at1,0.)
         fxy(1,1,k1) = zero
         fxy(2,1,k1) = zero
         fxy(3,1,k1) = zero
   50    continue
         k1 = ny + 1
         do 60 j = 1, nx
         at1 = aimag(ffd(j,1))
         fxy(1,j,1) = cmplx(0.,exy(1,j,1)*at1)
         fxy(2,j,1) = cmplx(0.,exy(2,j,1)*at1)
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
         fxy(1,j,k) = cmplx(0.,exy(1,j,k))
         fxy(2,j,k) = cmplx(0.,exy(2,j,k))
         fxy(3,j,k) = cmplx(exy(3,j,k),0.)
         fxy(1,j,k1) = cmplx(0.,exy(1,j,k))
         fxy(2,j,k1) = cmplx(0.,-exy(2,j,k))
         fxy(3,j,k1) = cmplx(exy(3,j,k),0.)
   70    continue
         fxy(1,1,k) = cmplx(0.,exy(1,1,k))
         fxy(2,1,k) = cmplx(0.,exy(2,1,k))
         fxy(3,1,k) = cmplx(exy(3,1,k),0.)
         fxy(1,1,k1) = zero
         fxy(2,1,k1) = zero
         fxy(3,1,k1) = zero
   80    continue
         k1 = ny + 1
         do 90 j = 1, nx
         fxy(1,j,1) = cmplx(0.,exy(1,j,1))
         fxy(2,j,1) = cmplx(0.,exy(2,j,1))
         fxy(3,j,1) = cmplx(exy(3,j,1),0.)
         fxy(1,j,k1) = zero
         fxy(2,j,k1) = zero
         fxy(3,j,k1) = zero
   90    continue
      endif
      return
      end
      subroutine CPFIELDN2(fxy,exy,nx,ny,nxv,ny2d,nxe,nye)
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
      call EMFIELDN2(fxy,exy,ffd,isign,nx,ny,nxv,ny2d,nxe,nye)
      return
      end
      subroutine AVPOTNX23(bxy,axy,nx,ny,nxv,ny2d)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with neumann boundary conditions (zero normal efield)
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
      axy(1,j,k1) = cmplx(0.,at4)
      axy(2,j,k1) = cmplx(0.,-at5)
      axy(3,j,k1) = cmplx(at6,0.)
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
      subroutine AVPOTN23(bxy,axy,nx,ny,nxe,nye)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with neumann boundary conditions (zero normal efield)
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
      axy(1,j,k) = at3*at4
      axy(2,j,k) = -at2*at4
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
      axy(1,1,k) = at3*at4
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
      axy(2,j,1) = -at2*at4
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
