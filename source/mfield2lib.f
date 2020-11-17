c 2d PIC library for solving field equations with mixed
c neumann/periodic boundary conditions
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: april 21, 2006
      subroutine SGLCOS2C(cu,cu1,nx,ny,nxv,nyv,nx2v)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in x,
c and y component is an even function in x.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx
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
      cu1(1,nx+j+1,k) = -cu(1,nx-j+1,k)
      cu1(2,nx+j+1,k) = cu(2,nx-j+1,k)
   10 continue
      cu1(1,1,k) = 0.
      cu1(2,1,k) = cu(2,1,k)
      cu1(1,nx+1,k) = 0.
      cu1(2,nx+1,k) = cu(2,nx+1,k)
   20 continue
      return
      end
      subroutine SGLCOS2D(q,q1,nx,ny,nxv,nyv,nx2v)
c this subroutine creates an even array q1 from an array q, so that
c a 1d sine transform can be performed with a 2d real to complex fft.
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = first dimension of input array q, must be >= nx
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
      q1(nx+j+1,k) = q(nx-j+1,k)
   10 continue
      q1(1,k) = q(1,k)
      q1(nx+1,k) = q(nx+1,k)
   20 continue
      return
      end
      subroutine SGLCOS2B(cu,cu1,nx,ny,nxv,nyv,nx2v)
c this subroutine creates a doubled vector array cu1 from a vector array
c cu, so that various 1d sine/cosine transforms can be performed with a
c 2d real to complex fft.  the x component is an odd function in x,
c y and z components are an even function in x.
c Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny = system length in x/y direction
c nxv = second dimension of input array cu, must be >= nx
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
      cu1(1,nx+j+1,k) = -cu(1,nx-j+1,k)
      cu1(2,nx+j+1,k) = cu(2,nx-j+1,k)
      cu1(3,nx+j+1,k) = cu(3,nx-j+1,k)
   10 continue
      cu1(1,1,k) = 0.
      cu1(2,1,k) = cu(2,1,k)
      cu1(3,1,k) = cu(3,1,k)
      cu1(1,nx+1,k) = 0.
      cu1(2,nx+1,k) = cu(2,nx+1,k)
      cu1(3,nx+1,k) = cu(2,nx+1,k)
   20 continue
      return
      end
      subroutine POISMNX2(q,fx,fy,isign,ffb,ax,ay,affp,we,nx,ny,nx2v,nyv
     1,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function.
c with mixed neumann/periodic boundary conditions.
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
      fx(2*j-1,k1) = -at3
      fx(2*j,k1) = -at2
      at3 = dky*at1
      at2 = at3*q(2*j,k)
      at3 = at3*q(2*j-1,k)
      fy(2*j-1,k) = at2
      fy(2*j,k) = -at3
      fy(2*j-1,k1) = at2
      fy(2*j,k1) = at3
      wp = wp + 2.0*at1*(q(2*j-1,k)**2 + q(2*j,k)**2)
   40 continue
c mode numbers kx = 0, nx
      at1 = ffb(1,k)*ffb(2,k)
      at3 = dky*at1
      at2 = at3*q(2,k)
      at3 = at3*q(1,k)
      fx(1,k) = 0.
      fx(2,k) = 0.
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      fy(1,k) = at2
      fy(2,k) = -at3
      fy(1,k1) = 0.
      fy(2,k1) = 0.
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nx
      at1 = ffb(2*j-1,1)*ffb(2*j,1)
      at3 = at1*q(2*j-1,1)
      at2 = dnx*float(j - 1)*at3
      fx(2*j-1,1) = 0.
      fx(2*j,1) = -at2
      fx(2*j-1,k1) = 0.
      fx(2*j,k1) = 0.
      fy(2*j-1,1) = 0.
      fy(2*j,1) = 0.
      fy(2*j-1,k1) = 0.
      fy(2*j,k1) = 0.
      wp = wp + at1*q(2*j-1,1)**2
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
      fx(2*j-1,k1) = at3
      fx(2*j,k1) = -at2
      wp = wp + 2.0*at1*(q(2*j-1,k)**2 + q(2*j,k)**2)
   80 continue
c mode numbers kx = 0, nx
      at2 = ffb(1,k)
      at1 = at2*ffb(2,k)
      at3 = at2*q(1,k)
      at2 = at2*q(2,k)
      fx(1,k) = at3
      fx(2,k) = at2
      fx(1,k1) = 0.
      fx(2,k1) = 0.
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
   90 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 100 j = 2, nx
      at2 = ffb(2*j-1,1)
      at1 = at2*ffb(2*j,1)
      fx(2*j-1,1) = at2*q(2*j-1,1)
      fx(2*j,1) = 0.
      fx(2*j-1,k1) = 0.
      fx(2*j,k1) = 0.
      wp = wp + at1*q(2*j-1,1)**2
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
      fy(2*j-1,k1) = at3
      fy(2*j,k1) = -at2
  120 continue
c mode numbers kx = 0, nx
      at1 = ffb(2,k)
      at3 = at1*q(1,k)
      at2 = at1*q(2,k)
      fy(1,k) = at3
      fy(2,k) = at2
      fy(1,k1) = 0.
      fy(2,k1) = 0.
  130 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 140 j = 2, nx
      at1 = ffb(2*j,1)
      fy(2*j-1,1) = at1*q(2*j-1,1)
      fy(2*j,1) = 0.
      fy(2*j-1,k1) = 0.
      fy(2*j,k1) = 0.
  140 continue
      fy(1,1) = ffb(2,1)*q(1,1)
      fy(2,1) = 0.
      fy(1,k1) = 0.
      fy(2,k1) = 0.
      return
      end
      subroutine POISMNX22(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,n
     1yhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed neumann/periodic boundary conditions.
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
      fxy(1,j,k1) = -at2*zt2
      fxy(2,j,k1) = at3*zt2
      wp = wp + 2.0*at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      at1 = real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dky*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(0.,-at2)
      fxy(2,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      wp = wp + at1*real(q(j,1))**2
   60 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      we = float(nx*ny)*wp
      return
      end
      subroutine POISMNX23(q,fxy,isign,ffb,ax,ay,affp,we,nx,ny,nxv,nyv,n
     1yhd)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with mixed neumann/periodic boundary conditions.
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
      fxy(1,j,k1) = -at2*zt2
      fxy(2,j,k1) = at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + 2.0*at1*q(j,k)*conjg(q(j,k))
   40 continue
c mode numbers kx = 0, nx
      at1 = real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dky*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, nx
      at1 = real(ffb(j,1))*aimag(ffb(j,1))
      at3 = at1*real(q(j,1))
      at2 = dnx*float(j - 1)*at3
      fxy(1,j,1) = cmplx(0.,-at2)
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*real(q(j,1))**2
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
      subroutine CUPERPMNX2(cu,nx,ny,nxv,nyv)
c this subroutine calculates the transverse current in fourier space
c with mixed neumann/periodic boundary conditions.
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
      cu(1,j,k1) = cu(1,j,k1) + dkx*zt2
      cu(2,j,k1) = cu(2,j,k1) - dky*zt2
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
      subroutine CUPERPMN2(cu,nx,ny,nxe,nyeh)
c this subroutine calculates the transverse current in fourier space
c with mixed neumann/periodic boundary conditions.
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
c nxe = first dimension of current array, must be >= nx
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
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
   10 continue
c mode numbers kx = 0, nx
      cu(2,1,k) = zero
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, nx
      cu(1,j,1) = zero
   30 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      return
      end
      subroutine CUPERPMNX22(cu,nx,ny,nxv,nyv)
c this subroutine calculates the transverse current in fourier space
c with mixed neumann/periodic boundary conditions.
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
      cu(1,j,k1) = cu(1,j,k1) + dkx*zt2
      cu(2,j,k1) = cu(2,j,k1) - dky*zt2
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
      subroutine BPOISMNX23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxv,
     1nyv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed neumann/periodic boundary conditions.
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
      bxy(1,j,k1) = conjg(zt1)
      bxy(2,j,k1) = -conjg(zt2)
      bxy(3,j,k1) = -conjg(zt3)
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
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1
     1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,-at2*real(cu(3,j,1)))
      bxy(3,j,1) = cmplx(0.,at2*real(cu(2,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2 + real(cu(
     13,j,1))**2)
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
      bxy(1,j,k1) = -conjg(zt1)
      bxy(2,j,k1) = conjg(zt2)
      bxy(3,j,k1) = conjg(zt3)
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
      bxy(1,1,k) = zero
      bxy(2,1,k) = at2*cu(2,1,k)
      bxy(3,1,k) = at2*cu(3,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1
     1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(0.,at2*aimag(cu(1,j,1)))
      bxy(2,j,1) = cmplx(at2*real(cu(2,j,1)),0.)
      bxy(3,j,1) = cmplx(at2*real(cu(3,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2 + real(cu(
     13,j,1))**2)
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
      bxy(1,j,k1) = -conjg(zt1)
      bxy(2,j,k1) = conjg(zt2)
      bxy(3,j,k1) = conjg(zt3)
  140 continue
  150 continue
c mode numbers kx = 0, nx
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = zero
      bxy(2,1,k) = at1*cu(2,1,k)
      bxy(3,1,k) = at1*cu(3,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(0.,at1*aimag(cu(1,j,1)))
      bxy(2,j,1) = cmplx(at1*real(cu(2,j,1)),0.)
      bxy(3,j,1) = cmplx(at1*real(cu(3,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
  170 continue
      at1 = aimag(ffb(1,1))
      bxy(1,1,1) = zero
      bxy(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
      bxy(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      return
      end
      subroutine BPOISMN23(cu,bxy,isign,ffb,ax,ay,affp,ci,wm,nx,ny,nxe,n
     1yeh,nxv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed neumann/periodic boundary conditions.
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
c nxe = first dimension of field arrays, must be >= nx
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
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
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   40 continue
   50 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 60 k = 2, nyh
      at1 = ci2*real(ffb(1,k))*aimag(ffb(1,k))
      at3 = dny*float(k - 1)*at1
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1
     1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,-at2*real(cu(3,j,1)))
      bxy(3,j,1) = cmplx(0.,at2*real(cu(2,j,1)))
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2 + real(cu(
     13,j,1))**2)
   70 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
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
      zt1 = at2*cu(1,j,k)
      zt2 = at2*cu(2,j,k)
      zt3 = at2*cu(3,j,k)
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(3,j,k) = zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)))
   90 continue
  100 continue
      wp = wp + wp
c mode numbers kx = 0, nx
      do 110 k = 2, nyh
      at2 = ci2*real(ffb(1,k))
      at1 = at2*aimag(ffb(1,k))
      bxy(1,1,k) = zero
      bxy(2,1,k) = at2*cu(2,1,k)
      bxy(3,1,k) = at2*cu(3,1,k)
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1
     1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(0.,at2*aimag(cu(1,j,1)))
      bxy(2,j,1) = cmplx(at2*real(cu(2,j,1)),0.)
      bxy(3,j,1) = cmplx(at2*real(cu(3,j,1)),0.)
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2 + real(cu(
     13,j,1))**2)
  120 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      wm = float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx and 0 < ky < ny/2
  130 do 150 k = 2, nyh
      do 140 j = 2, nx
      at1 = aimag(ffb(j,k))
      zt1 = at1*cu(1,j,k)
      zt2 = at1*cu(2,j,k)
      zt3 = at1*cu(3,j,k)
      bxy(1,j,k) = zt1
      bxy(2,j,k) = zt2
      bxy(3,j,k) = zt3
  140 continue
  150 continue
c mode numbers kx = 0, nx
      do 160 k = 2, nyh
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = zero
      bxy(2,1,k) = at1*cu(2,1,k)
      bxy(3,1,k) = at1*cu(3,1,k)
  160 continue
c mode numbers ky = 0, ny/2
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(0.,at1*aimag(cu(1,j,1)))
      bxy(2,j,1) = cmplx(at1*real(cu(2,j,1)),0.)
      bxy(3,j,1) = cmplx(at1*real(cu(3,j,1)),0.)
  170 continue
      at1 = aimag(ffb(1,1))
      bxy(1,1,1) = zero
      bxy(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
      bxy(3,1,1) = cmplx(at1*real(cu(3,1,1)),0.)
      return
      end
      subroutine BPOISMNX22(cu,bxy,bz,isign,ffb,ax,ay,affp,ci,wm,nx,ny,n
     1xv,nyv,nyhd)
c this subroutine solves 2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with mixed neumann/periodic boundary conditions.
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
      bz(j,k1) = -conjg(zt1)
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
      bz(1,k) = zero
      bz(1,k1) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)))
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nx
      at1 = ci2*real(ffb(j,1))*aimag(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      bz(j,1) = cmplx(0.,at2*real(cu(2,j,1)))
      bz(j,k1) = zero
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2)
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
      bxy(1,j,k1) = -conjg(zt1)
      bxy(2,j,k1) = conjg(zt2)
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
      bxy(1,1,k) = zero
      bxy(2,1,k) = at2*cu(2,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)))
  110 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 120 j = 2, nx
      at2 = ci2*real(ffb(j,1))
      at1 = at2*aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(0.,at2*aimag(cu(1,j,1)))
      bxy(2,j,1) = cmplx(at2*real(cu(2,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2)
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
      bxy(1,j,k1) = -conjg(zt1)
      bxy(2,j,k1) = conjg(zt2)
  140 continue
  150 continue
c mode numbers kx = 0, nx
cdir$ ivdep
      do 160 k = 2, nyh
      k1 = ny2 - k
      at1 = aimag(ffb(1,k))
      bxy(1,1,k) = zero
      bxy(2,1,k) = at1*cu(2,1,k)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
  160 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 170 j = 2, nx
      at1 = aimag(ffb(j,1))
      bxy(1,j,1) = cmplx(0.,at1*aimag(cu(1,j,1)))
      bxy(2,j,1) = cmplx(at1*real(cu(2,j,1)),0.)
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
  170 continue
      at1 = aimag(ffb(1,1))
      bxy(1,1,1) = zero
      bxy(2,1,1) = cmplx(at1*real(cu(2,1,1)),0.)
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      return
      end
      subroutine IBPOISMNX23(cu,bxy,ffb,ci,wm,nx,ny,nxv,nyv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with mixed neumann/periodic boundary conditions.
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
      bxy(1,j,k1) = at3*zt1
      bxy(2,j,k1) = at2*zt1
      bxy(3,j,k1) = at3*zt3 - at2*zt2
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
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1
     1,k)))
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffb(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,-at2*real(cu(3,j,1)))
      bxy(3,j,1) = cmplx(0.,at2*real(cu(2,j,1)))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2 + real(cu(
     13,j,1))**2)
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
      subroutine IBPOISMN23(cu,bxy,ffb,ci,wm,nx,ny,nxe,nyeh,nxv,nyhd)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field, with mixed neumann/periodic boundary conditions.
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
c nxe = first dimension of field arrays, must be >= nx
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp
      complex cu, bxy, ffb, zero, zt1, zt2, zt3
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
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
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
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
      wp = wp + at1*(cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1
     1,k)))
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      at1 = ci2*real(ffb(j,1))
      at2 = dnx*float(j - 1)*at1
      at1 = at1*aimag(ffb(j,1))
      bxy(1,j,1) = zero
      bxy(2,j,1) = cmplx(0.,-at2*real(cu(3,j,1)))
      bxy(3,j,1) = cmplx(0.,at2*real(cu(2,j,1)))
      wp = wp + at1*(aimag(cu(1,j,1))**2 + real(cu(2,j,1))**2 + real(cu(
     13,j,1))**2)
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      wm = float(nx*ny)*wp
      return
      end
      subroutine MAXWELMNX2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxv,nyv,nyh
     1d)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields
c with mixed neumann/periodic boundary conditions.
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
      exy(1,j,k1) = -conjg(zt7)
      exy(2,j,k1) = conjg(zt8)
      exy(3,j,k1) = conjg(zt9)
      bxy(1,j,k1) = conjg(zt4)
      bxy(2,j,k1) = -conjg(zt5)
      bxy(3,j,k1) = -conjg(zt6)
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
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
c update electric field whole time step
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt8 = exy(2,1,k) - afdt*cu(2,1,k)
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      wp = wp + anorm*zt4*conjg(zt4)
      exy(1,1,k) = zero
      exy(2,1,k) = zt8
      exy(3,1,k) = zt9
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
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
      at1 = real(exy(3,j,1))
      at2 = real(exy(2,j,1))
      at5 = aimag(bxy(2,j,1)) + dth*(dkx*at1)
      at6 = aimag(bxy(3,j,1)) - dth*(dkx*at2)
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
      bxy(2,j,1) = cmplx(0.,at5)
      bxy(3,j,1) = cmplx(0.,at6)
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
      subroutine MAXWELMN2(exy,bxy,cu,ffb,ci,dt,wf,wm,nx,ny,nxe,nyeh,nxv
     1,nyhd)
c this subroutine solves 2d maxwell's equation in fourier space for
c transverse electric and magnetic fields
c with mixed neumann/periodic boundary conditions.
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
c nxe = first dimension of field arrays, must be >= nx
c nyeh = second dimension of field arrays, must be >= nyh
c nxv = first dimension of form factor array, must be >= nx
c nyhd = second dimension of form factor array, must be >= nyh
      double precision wp, ws
      complex exy, bxy, cu, ffb
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
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
   10 continue
   20 continue
      ws = ws + ws
      wp = wp + wp
c mode numbers kx = 0, nx
      do 30 k = 2, nyh
      dky = dny*float(k - 1)
      afdt = adt*aimag(ffb(1,k))
c update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
c update electric field whole time step
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt8 = exy(2,1,k) - afdt*cu(2,1,k)
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
c update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      wp = wp + anorm*zt4*conjg(zt4)
      exy(1,1,k) = zero
      exy(2,1,k) = zt8
      exy(3,1,k) = zt9
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zero
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      afdt = adt*aimag(ffb(j,1))
c update magnetic field half time step
      at1 = real(exy(3,j,1))
      at2 = real(exy(2,j,1))
      at5 = aimag(bxy(2,j,1)) + dth*(dkx*at1)
      at6 = aimag(bxy(3,j,1)) - dth*(dkx*at2)
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
      bxy(2,j,1) = cmplx(0.,at5)
      bxy(3,j,1) = cmplx(0.,at6)
   40 continue
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      wf = float(nx*ny)*ws
      wm = float(nx*ny)*c2*wp
      return
      end
      subroutine CMFIELDMN2(cu1,cu,nx,ny,nxv,nyv,nxe,nyeh)
c this subroutine copies the current into a smaller array
      implicit none
      integer nx, ny, nxv, nyv, nxe, nyeh
      complex cu1, cu
      dimension cu1(3,nxv,nyv), cu(3,nxe,nyeh)
      integer j, k, nyh
      nyh = ny/2
      do 20 k = 1, nyh
      do 10 j = 1, nx
      cu(1,j,k) = cu1(1,j,k)
      cu(2,j,k) = cu1(2,j,k)
      cu(3,j,k) = cu1(3,j,k)
   10 continue
   20 continue
      return
      end
      subroutine EMFIELDMN2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,nyh
     1d)
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
         fxy(1,j,k) = fxy(1,j,k) + zt2
         fxy(2,j,k) = fxy(2,j,k) + zt3
         fxy(3,j,k) = fxy(3,j,k) + zt4
         fxy(1,j,k1) = fxy(1,j,k1) - conjg(zt2)
         fxy(2,j,k1) = fxy(2,j,k1) + conjg(zt3)
         fxy(3,j,k1) = fxy(3,j,k1) + conjg(zt4)
   10    continue
         at1 = aimag(ffb(1,k))
         fxy(1,1,k) = fxy(1,1,k) + exy(1,1,k)*at1
         fxy(2,1,k) = fxy(2,1,k) + exy(2,1,k)*at1
         fxy(3,1,k) = fxy(3,1,k) + exy(3,1,k)*at1
   20    continue
         k1 = nyh + 1
         do 30 j = 1, nx
         at1 = aimag(ffb(j,1))
         fxy(1,j,1) = fxy(1,j,1) + exy(1,j,1)*at1
         fxy(2,j,1) = fxy(2,j,1) + exy(2,j,1)*at1
         fxy(3,j,1) = fxy(3,j,1) + exy(3,j,1)*at1
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
         fxy(1,j,k) = zt2
         fxy(2,j,k) = zt3
         fxy(3,j,k) = zt4
         fxy(1,j,k1) = conjg(zt2)
         fxy(2,j,k1) = -conjg(zt3)
         fxy(3,j,k1) = -conjg(zt4)
   40    continue
         at1 = aimag(ffb(1,k))
         fxy(1,1,k) = exy(1,1,k)*at1
         fxy(2,1,k) = exy(2,1,k)*at1
         fxy(3,1,k) = exy(3,1,k)*at1
         fxy(1,1,k1) = zero
         fxy(2,1,k1) = zero
         fxy(3,1,k1) = zero
   50    continue
         k1 = nyh + 1
         do 60 j = 1, nx
         at1 = aimag(ffb(j,1))
         fxy(1,j,1) = exy(1,j,1)*at1
         fxy(2,j,1) = exy(2,j,1)*at1
         fxy(3,j,1) = exy(3,j,1)*at1
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
         zt3 = exy(2,j,k)
         zt4 = exy(3,j,k)
         fxy(1,j,k) = zt2
         fxy(2,j,k) = zt3
         fxy(3,j,k) = zt4
         fxy(1,j,k1) = -conjg(zt2)
         fxy(2,j,k1) = conjg(zt3)
         fxy(3,j,k1) = conjg(zt4)
   70    continue
         fxy(1,1,k) = exy(1,1,k)
         fxy(2,1,k) = exy(2,1,k)
         fxy(3,1,k) = exy(3,1,k)
         fxy(1,1,k1) = zero
         fxy(2,1,k1) = zero
         fxy(3,1,k1) = zero
   80    continue
         k1 = nyh + 1
         do 90 j = 1, nx
         fxy(1,j,1) = exy(1,j,1)
         fxy(2,j,1) = exy(2,j,1)
         fxy(3,j,1) = exy(3,j,1)
         fxy(1,j,k1) = zero
         fxy(2,j,k1) = zero
         fxy(3,j,k1) = zero
   90    continue
      endif
      return
      end
      subroutine CPFIELDMN2(fxy,exy,nx,ny,nxv,nyv,nxe,nyeh)
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
      call EMFIELDMN2(fxy,exy,ffb,isign,nx,ny,nxv,nyv,nxe,nyeh,nyhd)
      return
      end
      subroutine AVPOTMNX23(bxy,axy,nx,ny,nxv,nyv)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with mixed neumann/periodic boundary conditions.
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
      axy(1,j,k1) = -conjg(zt1)
      axy(2,j,k1) = conjg(zt2)
      axy(3,j,k1) = conjg(zt3)
   10 continue
   20 continue
c mode numbers kx = 0, nx
cdir$ ivdep
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      axy(1,1,k) = zero
      axy(2,1,k) = zero
      axy(3,1,k) = -at3*zt3
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
      axy(2,j,1) = cmplx(at2*aimag(bxy(3,j,1)),0.)
      axy(3,j,1) = cmplx(-at2*aimag(bxy(2,j,1)),0.)
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
      subroutine AVPOTMN23(bxy,axy,nx,ny,nxe,nyeh)
c this subroutine calculates 2-1/2d vector potential from magnetic field
c in fourier space with mixed neumann/periodic boundary conditions.
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
c nxe = first dimension of field arrays, must be >= nx
c nyeh = second dimension of field arrays, must be >= nyh
      complex bxy, axy, zero, zt1, zt2, zt3
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
      zt1 = cmplx(-aimag(bxy(3,j,k)),real(bxy(3,j,k)))
      zt2 = cmplx(-aimag(bxy(2,j,k)),real(bxy(2,j,k)))
      zt3 = cmplx(-aimag(bxy(1,j,k)),real(bxy(1,j,k)))
      zt3 = at2*zt2 - at3*zt3
      zt2 = -at2*zt1
      zt1 = at3*zt1
      axy(1,j,k) = zt1
      axy(2,j,k) = zt2
      axy(3,j,k) = zt3
   10 continue
   20 continue
c mode numbers kx = 0, nx
      do 30 k = 2, nyh
      dky = dny*float(k - 1)
      at3 = 1.0/dky
      zt3 = cmplx(-aimag(bxy(1,1,k)),real(bxy(1,1,k)))
      axy(1,1,k) = zero
      axy(2,1,k) = zero
      axy(3,1,k) = -at3*zt3
   30 continue
c mode numbers ky = 0, ny/2
      do 40 j = 2, nx
      dkx = dnx*float(j - 1)
      at2 = 1.0/dkx
      axy(1,j,1) = zero
      axy(2,j,1) = cmplx(at2*aimag(bxy(3,j,1)),0.)
      axy(3,j,1) = cmplx(-at2*aimag(bxy(2,j,1)),0.)
   40 continue
      axy(1,1,1) = zero
      axy(2,1,1) = zero
      axy(3,1,1) = zero
      return
      end
