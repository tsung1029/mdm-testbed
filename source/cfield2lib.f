c 2d PIC library for solving field equations with open (vacuum) boundary
c conditions
c written by viktor k. decyk, ucla
c copyright 1991, regents of the university of california
c update: october 17, 2007
      subroutine ZDBL2C(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a zeroed array cu2 from an array cu, so that
c a 2d real to complex transform will perform a correct convolution.
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
      integer i, j, k
c copy to double array
      do 30 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 2
      cu2(i,j,k) = cu(i,j,k)
      cu2(i,nx+j,k) = 0.
      cu2(i,j,ny+k) = 0.
      cu2(i,nx+j,ny+k) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
      subroutine ZDBL2D(q,q2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a zeroed array q2 from an array q, so that
c a 2d real to complex transform will perform a correct convolution.
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
      integer j, k
c copy to double array
      do 20 k = 1, ny
      do 10 j = 1, nx
      q2(j,k) = q(j,k)
      q2(nx+j,k) = 0.
      q2(j,ny+k) = 0.
      q2(nx+j,ny+k) = 0.
   10 continue
   20 continue
      return
      end
      subroutine ZDBL2B(cu,cu2,nx,ny,nxv,nyv,nx2v,ny2)
c this subroutine creates a zeroed array cu2 from an array cu, so that
c a 2d real to complex transform will perform a correct convolution.
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
      integer i, j, k
c copy to double array
      do 30 k = 1, ny
      do 20 j = 1, nx
      do 10 i = 1, 3
      cu2(i,j,k) = cu(i,j,k)
      cu2(i,nx+j,k) = 0.
      cu2(i,j,ny+k) = 0.
      cu2(i,nx+j,ny+k) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
      subroutine FORMC2(ffg,f,fpotc,mixup2,sct2,affp,ar,indx1,indy1,nx1d
     1,ny1d,nxv,ny2v,nxhy2,nxyh2)
c this subroutine calculates the form factor array ffg needed by field
c solvers with open (vacuum) boundary conditions using hockney's method.
c the four green's functions calculated are:
c g(kx,ky) = affp*inverse FFT of potr
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c gx(kx,ky) = affp*s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = affp*s(kx,ky)*inverse FFT of (y/r)*Er
c where the fields due to the finite-sized particles are given by fpotc
c input: mixup2,sct2,affp,ar,indx1,indy1,nx1d,ny1d,nx2v,ny2v,nxhy2,nxyh2
c output: ffg, f
c ffg(1,j,k) = potential green's function g
c ffg(2,j,k) = finite-size particle shape factor s
c ffg(3,j,k) = x component of electric field green's function gx
c ffg(4,j,k) = y component of electric field green's function gy
c f = scratch array used by FFT
c mixup2/sct2 = bit-reverse and sine-cosine table used by FFT
c affp = normalization constant = nx*ny/np, where np=number of particles
c ar = half-width of particle in r direction
c indx1/indy1 = exponent which determines FFT length in x/y direction,
c where 2*nx=2**indx1, 2*ny=2**indy1
c nx1d = second dimension of ffg arrays, must be >= nx+1
c ny1d = third dimension of field arrays, must be >= ny+1
c nxv = half of first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nxhy2 = maximum of (nx,2*ny)
c nxyh2 = maximum of (nx,ny)
      implicit none
      real ffg, f
      integer mixup2
      complex sct2
      real affp, ar
      integer indx1, indy1, nx1d, ny1d, nxv, ny2v, nxhy2, nxyh2
      dimension ffg(4,nx1d,ny1d), f(2*nxv,ny2v)
      dimension mixup2(nxhy2), sct2(nxyh2)
      real fpotc
      external fpotc
c local data
      integer nx, ny, nx1, ny1, nx2, ny2, isign, j, k, j1, k1, ifun
      real an, ari, at1, x, y, r
      real POTC2, erfn
      external POTC2, erfn
      nx2 = 2**(indx1)
      ny2 = 2**(indy1)
      nx = nx2/2
      ny = ny2/2
      nx1 = nx + 1
      ny1 = ny + 1
      ari = 0.0
      if (ar.gt.0.) ari = 1.0/ar
      an = float(nx2*ny2)
c calculate potential green's function
      ifun = 1
      do 20 k = 1, ny2
      k1 = k - 1
      if (k1.gt.ny) k1 = k1 - ny2
      at1 = float(k1)**2
      do 10 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at1 + float(j1)**2)
      f(j,k) = fpotc(r,affp,ari,1)
   10 continue
   20 continue
      isign = -1
      call WFFT2RX(f,isign,mixup2,sct2,indx1,indy1,nxv,ny2v,nxhy2,nxyh2)
      do 40 k = 1, ny1
      do 30 j = 1, nx
      ffg(ifun,j,k) = an*f(2*j-1,k)
   30 continue
   40 continue
      do 50 k = 2, ny
      k1 = ny2 + 2 - k
      ffg(ifun,nx1,k) = an*f(1,k1)
   50 continue
      ffg(ifun,nx1,1) = an*f(2,1)
      ffg(ifun,nx1,ny1) = an*f(2,ny1)
c calculate particle smoothing function
      ifun = ifun + 1
      do 70 k = 1, ny2
      k1 = k - 1
      if (k1.gt.ny) k1 = k1 - ny2
      at1 = float(k1)**2
      do 60 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at1 + float(j1)**2)
      f(j,k) = POTC2(r,affp,ari,2)
   60 continue
   70 continue
      isign = -1
      call WFFT2RX(f,isign,mixup2,sct2,indx1,indy1,nxv,ny2v,nxhy2,nxyh2)
      do 90 k = 1, ny1
      do 80 j = 1, nx
      ffg(ifun,j,k) = an*f(2*j-1,k)
   80 continue
   90 continue
      do 100 k = 2, ny
      k1 = ny2 + 2 - k
      ffg(ifun,nx1,k) = an*f(1,k1)
  100 continue
      ffg(ifun,nx1,1) = an*f(2,1)
      ffg(ifun,nx1,ny1) = an*f(2,ny1)
c calculate green's function for x component of electric field
      ifun = ifun + 1
      do 120 k = 1, ny2
      k1 = k - 1
      if (k1.gt.ny) k1 = k1 - ny2
      at1 = float(k1)**2
      do 110 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      x = float(j1)
      r = sqrt(at1 + x*x)
      f(j,k) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k) = f(j,k)*(x/r)
  110 continue
  120 continue
      isign = -1
      call WFFT2RX(f,isign,mixup2,sct2,indx1,indy1,nxv,ny2v,nxhy2,nxyh2)
      do 140 k = 1, ny1
      do 130 j = 2, nx
      ffg(ifun,j,k) = an*f(2*j,k)
  130 continue
      ffg(ifun,1,k) = an*f(1,k)
  140 continue
      do 150 k = 2, ny
      k1 = ny2 + 2 - k
      ffg(ifun,nx1,k) = an*f(1,k1)
  150 continue
      ffg(ifun,nx1,1) = an*f(2,1)
      ffg(ifun,nx1,ny1) = an*f(2,ny1)
c calculate green's function for y component of electric field
      ifun = ifun + 1
      do 170 k = 1, ny2
      k1 = k - 1
      if (k1.gt.ny) k1 = k1 - ny2
      y = float(k1)
      at1 = y*y
      do 160 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at1 + float(j1)**2)
      f(j,k) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k) = f(j,k)*(y/r)
  160 continue
  170 continue
      isign = -1
      call WFFT2RX(f,isign,mixup2,sct2,indx1,indy1,nxv,ny2v,nxhy2,nxyh2)
      do 190 k = 2, ny
      k1 = ny2 + 2 - k
      do 180 j = 1, nx
      ffg(ifun,j,k) = an*f(2*j,k)
  180 continue
      ffg(ifun,nx1,k) = an*f(2,k1)
  190 continue
      do 200 j = 1, nx
      ffg(ifun,j,1) = an*f(2*j-1,1)
      ffg(ifun,j,ny1) = an*f(2*j-1,ny1)
  200 continue
      ffg(ifun,nx1,1) = an*f(2,1)
      ffg(ifun,nx1,ny1) = an*f(2,ny1)
      return
      end
      subroutine POISC2(q,fx,fy,isign,ffg,we,nx,ny,nx2v,ny2d,nx1d,ny1d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides a smoothing function
c with open (vacuum) boundary conditions using hockney's method.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c for isign = -1, input: q,ffg,isign,nx,ny,nx1d,ny1d,nx2v,ny2d
c output: fx,fy,we
c approximate flop count is: 35*nx*ny + 28*(nx + ny)
c for isign = 1, input: q,ffg,isign,nx,ny,,nx1d,ny1d,nx2v,ny2d
c output: fx,we
c approximate flop count is: 14*nx*ny + 16*(nx + ny)
c for isign = 2, input: q,ffg,isign,nx,ny,nx1d,ny1dnx2v,ny2d
c output: fy
c approximate flop count is: 4*nx*ny + 4*(nx + ny)
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky) = gx(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = gy(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c gx(kx,ky) = s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = s(kx,ky)*inverse FFT of (y/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky) = g(kx,ky)*q(kx,ky)
c where g(kx,ky) = affp*inverse FFT of potr
c where potr is the potential of a single finite-sized particle
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky) = q(kx,ky)*s(kx,ky)
c cmplx(q(2*j-1,k),q(2*j,k)) = complex charge density
c cmplx(fx(2*j-1,k),fx(2*j,k)) = x component of complex force/charge,
c cmplx(fy(2*j-1,k),fy(2*j,k)) = y component of complex force/charge,
c ffg(1,j,k) = potential green's function g
c ffg(2,j,k) = finite-size particle shape factor s
c ffg(3,j,k) = x component of electric field green's function gx
c ffg(4,j,k) = y component of electric field green's function gy
c all for fourier mode (j-1,k-1)
c the ffg array is calculated by the subroutine FORMC2
c electric field energy is also calculated and returned in we
c nx/ny = system length in x/y direction
c nx2v = first dimension of field arrays, must be >= 2*nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nx1d = second dimension of ffg arrays, must be >= nx+1
c ny1d = third dimension of ffg arrays, must be >= ny+1
      implicit none
      real ffg, q, fx, fy
      integer isign, nx, ny, nx2v, ny2d, nx1d, ny1d
      real we
      dimension q(nx2v,ny2d), fx(nx2v,ny2d), fy(nx2v,ny2d)
      dimension ffg(4,nx1d,ny1d)
c local data
      double precision wp
      integer nx1, ny1, ny22, j, k, k1
      real at1, at2, at3, at4, at5
      if (isign.eq.0) return
      nx1 = nx + 1
      ny1 = ny + 1
      ny22 = ny + ny + 2
      if (isign.gt.0) go to 70
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      at5 = 1.0
      do 20 k = 2, ny
      k1 = ny22 - k
      at1 = ffg(3,1,k)
      at5 = -at5
      do 10 j = 2, nx
      at1 = -at1
      at2 = ffg(3,j,k)
      at3 = at5*ffg(4,j,1)
      at4 = ffg(4,j,k)
      fx(2*j-1,k) = at1*q(2*j-1,k) - at2*q(2*j,k)
      fx(2*j,k) = at2*q(2*j-1,k) + at1*q(2*j,k)
      fx(2*j-1,k1) = at1*q(2*j-1,k1) - at2*q(2*j,k1)
      fx(2*j,k1) = at2*q(2*j-1,k1) + at1*q(2*j,k1)
      fy(2*j-1,k) = at3*q(2*j-1,k) - at4*q(2*j,k)
      fy(2*j,k) = at4*q(2*j-1,k) + at3*q(2*j,k)
      fy(2*j-1,k1) = at3*q(2*j-1,k1) + at4*q(2*j,k1)
      fy(2*j,k1) = at3*q(2*j,k1) - at4*q(2*j-1,k1)
      wp = wp + ffg(1,j,k)*(q(2*j-1,k)**2 + q(2*j,k)**2 + q(2*j-1,k1)**2
     1 + q(2*j,k1)**2)
   10 continue
   20 continue
c mode number kx = 0
      at3 = ffg(4,1,1)
      do 30 k = 2, ny
      at1 = ffg(3,1,k)
      at3 = -at3
      at4 = ffg(4,1,k)
      fx(1,k) = at1*q(1,k)
      fx(2,k) = at1*q(2,k)
      fy(1,k) = at3*q(1,k) - at4*q(2,k)
      fy(2,k) = at4*q(1,k) + at3*q(2,k)
      wp = wp + ffg(1,1,k)*(q(1,k)**2 + q(2,k)**2)
   30 continue
c mode number kx = nx/2
      at3 = ffg(4,nx1,1)
      do 40 k = 2, ny
      k1 = ny22 - k
      at1 = ffg(3,nx1,k)
      at3 = -at3
      at4 = ffg(4,nx1,k)
      fx(1,k1) = at1*q(1,k1)
      fx(2,k1) = at1*q(2,k1)
      fy(1,k1) = at3*q(1,k1) - at4*q(2,k1)
      fy(2,k1) = at4*q(1,k1) + at3*q(2,k1)
      wp = wp + ffg(1,nx1,k)*(q(1,k1)**2 + q(2,k1)**2)
   40 continue
c mode number ky = 0
      at1 = ffg(3,1,1)
      do 50 j = 2, nx
      at1 = -at1
      at2 = ffg(3,j,1)
      at4 = ffg(4,j,1)
      fx(2*j-1,1) = at1*q(2*j-1,1) - at2*q(2*j,1)
      fx(2*j,1) = at2*q(2*j-1,1) + at1*q(2*j,1)
      fy(2*j-1,1) = at4*q(2*j-1,1)
      fy(2*j,1) = at4*q(2*j,1)
      wp = wp + ffg(1,j,1)*(q(2*j-1,1)**2 + q(2*j,1)**2)
   50 continue
c mode number ky = ny/2
      at1 = ffg(3,1,ny1)
      do 60 j = 2, nx
      at1 = -at1
      at2 = ffg(3,j,ny1)
      at4 = ffg(4,j,ny1)
      fx(2*j-1,ny1) = at1*q(2*j-1,ny1) - at2*q(2*j,ny1)
      fx(2*j,ny1) = at2*q(2*j-1,ny1) + at1*q(2*j,ny1)
      fy(2*j-1,ny1) = at4*q(2*j-1,ny1)
      fy(2*j,ny1) = at4*q(2*j,ny1)
      wp = wp + ffg(1,j,ny1)*(q(2*j-1,ny1)**2 + q(2*j,ny1)**2)
   60 continue
c mode numbers ky = 0, kx = 0, nx/2
      fx(1,1) = ffg(3,1,1)*q(1,1)
      fx(2,1) = ffg(3,nx1,1)*q(2,1)
      fy(1,1) = ffg(4,1,1)*q(1,1)
      fy(2,1) = ffg(4,nx1,1)*q(2,1)
      wp = wp + .5*(ffg(1,1,1)*q(1,1)**2 + ffg(1,nx1,1)*q(2,1)**2)
c mode numbers ky = ny/2, kx = 0, nx/2
      fx(1,ny1) = ffg(3,1,ny1)*q(1,ny1)
      fx(2,ny1) = ffg(3,nx1,ny1)*q(2,ny1)
      fy(1,ny1) = ffg(4,1,ny1)*q(1,ny1)
      fy(2,ny1) = ffg(4,nx1,ny1)*q(2,ny1)
      wp = wp + .5*(ffg(1,1,ny1)*q(1,ny1)**2 + ffg(1,nx1,ny1)*q(2,ny1)**
     12)
      we = 4.0*float(nx*ny)*wp
      return
c calculate potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 90 k = 2, ny
      k1 = ny22 - k
      do 80 j = 2, nx
      at2 = ffg(1,j,k)
c     at1 = at2*ffg(2,j,k)
      at1 = at2
      fx(2*j-1,k) = at2*q(2*j-1,k)
      fx(2*j,k) = at2*q(2*j,k)
      fx(2*j-1,k1) = at2*q(2*j-1,k1)
      fx(2*j,k1) = at2*q(2*j,k1)
      wp = wp + at1*(q(2*j-1,k)**2 + q(2*j,k)**2 + q(2*j-1,k1)**2 + q(2*
     1j,k1)**2)
   80 continue
c mode number kx = 0
      at2 = ffg(1,1,k)
c     at1 = at2*ffg(2,1,k)
      at1 = at2
      fx(1,k) = at2*q(1,k)
      fx(2,k) = at2*q(2,k)
      wp = wp + at1*(q(1,k)**2 + q(2,k)**2)
c mode number kx = nx/2
      at2 = ffg(1,nx1,k)
c     at1 = at2*ffg(2,nx1,k)
      at1 = at2
      fx(1,k1) = at2*q(1,k1)
      fx(2,k1) = at2*q(2,k1)
      wp = wp + at1*(q(1,k1)**2 + q(2,k1)**2)
   90 continue
c mode numbers ky = 0, ny/2
      do 100 j = 2, nx
      at2 = ffg(1,j,1)
c     at1 = at2*ffg(2,j,1)
      at1 = at2
      fx(2*j-1,1) = at2*q(2*j-1,1)
      fx(2*j,1) = at2*q(2*j,1)
      wp = wp + at1*(q(2*j-1,1)**2 + q(2*j,1)**2)
      at2 = ffg(1,j,ny1)
c     at1 = at2*ffg(2,j,ny1)
      at1 = at2
      fx(2*j-1,ny1) = at2*q(2*j-1,ny1)
      fx(2*j,ny1) = at2*q(2*j,ny1)
      wp = wp + at1*(q(2*j-1,ny1)**2 + q(2*j,ny1)**2)
  100 continue
c mode numbers ky = 0, kx = 0
      at2 = ffg(1,1,1)
c     at1 = at2*ffg(2,1,1)
      at1 = at2
      fx(1,1) = at2*q(1,1)
      wp = wp + .5*at1*q(1,1)**2
c mode numbers ky = 0, kx = nx/2
      at2 = ffg(1,nx1,1)
c     at1 = at2*ffg(2,nx1,1)
      at1 = at2
      fx(2,1) = at2*q(2,1)
      wp = wp + .5*at1*q(2,1)**2
c mode numbers ky = ny/2, kx = 0
      at2 = ffg(1,1,ny1)
c     at1 = at2*ffg(2,1,ny1)
      at1 = at2
      fx(1,ny1) = at2*q(1,ny1)
      wp = wp + .5*at1*q(1,ny1)**2
c mode numbers ky = ny/2, kx = nx/2
      at2 = ffg(1,nx1,ny1)
c     at1 = at2*ffg(2,nx1,ny1)
      at1 = at2
      fx(2,ny1) = at2*q(2,ny1)
      wp = wp + .5*at1*q(2,ny1)**2
      we = 4.0*float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  110 do 130 k = 2, ny
      k1 = ny22 - k
      do 120 j = 2, nx
      at1 = ffg(2,j,k)
      fy(2*j-1,k) = at1*q(2*j-1,k)
      fy(2*j,k) = at1*q(2*j,k)
      fy(2*j-1,k1) = at1*q(2*j-1,k1)
      fy(2*j,k1) = at1*q(2*j,k1)
  120 continue
c mode number kx = 0
      at1 = ffg(2,1,k)
      fy(1,k) = at1*q(1,k)
      fy(2,k) = at1*q(2,k)
c mode number kx = nx/2
      at1 = ffg(2,nx1,k)
      fy(1,k1) = at1*q(1,k1)
      fy(2,k1) = at1*q(2,k1)
  130 continue
c mode numbers ky = 0, ny/2
      do 140 j = 2, nx
      at1 = ffg(2,j,1)
      fy(2*j-1,1) = at1*q(2*j-1,1)
      fy(2*j,1) = at1*q(2*j,1)
      at1 = ffg(2,j,ny1)
      fy(2*j-1,ny1) = at1*q(2*j-1,ny1)
      fy(2*j,ny1) = at1*q(2*j,ny1)
  140 continue
c mode numbers ky = 0, kx = 0, nx/2
      fy(1,1) = ffg(2,1,1)*q(1,1)
      fy(2,1) = ffg(2,nx1,1)*q(2,1)
c mode numbers ky = ny/2, kx = 0, nx/2
      fy(1,ny1) = ffg(2,1,ny1)*q(1,ny1)
      fy(2,ny1) = ffg(2,nx1,ny1)*q(2,ny1)
      return
      end
      subroutine POISC22(q,fxy,ffg,we,nx,ny,nxv,ny2d,nx1d,ny1d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with open (vacuum) boundary conditions using hockney's method.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c input: q,ffg,isign,nx,ny,nx1d,ny1d,nx2v,ny2d, output: fxy,we
c approximate flop count is: 44*nx*ny + 36*(nx + ny)
c force/charge is calculated using the equations:
c fx(kx,ky) = gx(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = gy(kx,ky)*s(kx,ky)*q(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c gx(kx,ky) = s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = s(kx,ky)*inverse FFT of (y/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c q(j,k) = complex charge density
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c ffg(1,j,k) = potential green's function g
c ffg(2,j,k) = finite-size particle shape factor s
c ffg(3,j,k) = x component of electric field green's function gx
c ffg(4,j,k) = y component of electric field green's function gy
c all for fourier mode (j-1,k-1)
c the ffg array is calculated by the subroutine FORMC2
c electric field energy is also calculated and returned in we
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nx1d = second dimension of ffg arrays, must be >= nx+1
c ny1d = third dimension of ffg arrays, must be >= ny+1
      implicit none
      complex q, fxy
      real ffg
      integer nx, ny, nxv, ny2d, nx1d, ny1d
      real we
      dimension q(nxv,ny2d), fxy(2,nxv,ny2d)
      dimension ffg(4,nx1d,ny1d)
c local data
      double precision wp
      integer nx1, ny1, ny22, j, k, k1
      real at1, at2, at3
      complex zt1, zt2
      nx1 = nx + 1
      ny1 = ny + 1
      ny22 = ny + ny + 2
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      at3 = 1.0
      do 20 k = 2, ny
      k1 = ny22 - k
      at1 = ffg(3,1,k)
      at3 = -at3
      do 10 j = 2, nx
      at1 = -at1
      at2 = at3*ffg(4,j,1)
      zt1 = cmplx(at1,ffg(3,j,k))
      zt2 = cmplx(at2,ffg(4,j,k))
      fxy(1,j,k) = zt1*q(j,k)
      fxy(2,j,k) = zt2*q(j,k)
      fxy(1,j,k1) = zt1*q(j,k1)
      fxy(2,j,k1) = conjg(zt2)*q(j,k1)
      wp = wp + ffg(1,j,k)*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)
     1))
   10 continue
   20 continue
c mode number kx = 0
      at3 = ffg(4,1,1)
      do 30 k = 2, ny
      at3 = -at3
      zt1 = cmplx(at3,ffg(4,1,k))
      fxy(1,1,k) = ffg(3,1,k)*q(1,k)
      fxy(2,1,k) = zt1*q(1,k)
      wp = wp + ffg(1,1,k)*(q(1,k)*conjg(q(1,k)))
   30 continue
c mode number kx = nx/2
      at3 = ffg(4,nx1,1)
      do 40 k = 2, ny
      k1 = ny22 - k
      at3 = -at3
      zt1 = cmplx(at3,ffg(4,nx1,k))
      fxy(1,1,k1) = ffg(3,nx1,k)*q(1,k1)
      fxy(2,1,k1) = zt1*q(1,k1)
      wp = wp + ffg(1,nx1,k)*(q(1,k1)*conjg(q(1,k1)))
   40 continue
c mode number ky = 0
      at1 = ffg(3,1,1)
      do 50 j = 2, nx
      at1 = -at1
      zt1 = cmplx(at1,ffg(3,j,1))
      fxy(1,j,1) = zt1*q(j,1)
      fxy(2,j,1) = ffg(4,j,1)*q(j,1)
      wp = wp + ffg(1,j,1)*(q(j,1)*conjg(q(j,1)))
   50 continue
c mode number ky = ny/2
      at1 = ffg(3,1,ny1)
      do 60 j = 2, nx
      at1 = -at1
      zt1 = cmplx(at1,ffg(3,j,ny1))
      fxy(1,j,ny1) = zt1*q(j,ny1)
      fxy(2,j,ny1) = ffg(4,j,ny1)*q(j,ny1)
      wp = wp + ffg(1,j,ny1)*(q(j,ny1)*conjg(q(j,ny1)))
   60 continue
c mode numbers ky = 0, kx = 0, nx/2
      fxy(1,1,1) = cmplx(ffg(3,1,1)*real(q(1,1)),ffg(3,nx1,1)*aimag(q(1,
     11)))
      fxy(2,1,1) = cmplx(ffg(4,1,1)*real(q(1,1)),ffg(4,nx1,1)*aimag(q(1,
     11)))
      wp = wp + .5*(ffg(1,1,1)*real(q(1,1))**2 + ffg(1,nx1,1)*aimag(q(1,
     11))**2)
c mode numbers ky = ny/2, kx = 0, nx/2
      fxy(1,1,ny1) = cmplx(ffg(3,1,ny1)*real(q(1,ny1)),ffg(3,nx1,ny1)*ai
     1mag(q(1,ny1)))
      fxy(2,1,ny1) = cmplx(ffg(4,1,ny1)*real(q(1,ny1)),ffg(4,nx1,ny1)*ai
     1mag(q(1,ny1)))
      wp = wp + .5*(ffg(1,1,ny1)*real(q(1,ny1))**2 + ffg(1,nx1,ny1)*aima
     1g(q(1,ny1))**2)
      we = 4.0*float(nx*ny)*wp
      return
      end
      subroutine POISC23(q,fxy,ffg,we,nx,ny,nxv,ny2d,nx1d,ny1d)
c this subroutine solves 2d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with open (vacuum) boundary conditions using hockney's method.
c Zeros out z component
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c input: q,ffg,isign,nx,ny,nx1d,ny1d,nx2v,ny2d, output: fxy,we
c approximate flop count is: 44*nx*ny + 36*(nx + ny)
c force/charge is calculated using the equations:
c fx(kx,ky) = gx(kx,ky)*s(kx,ky)*q(kx,ky),
c fy(kx,ky) = gy(kx,ky)*s(kx,ky)*q(kx,ky),
c fz(kx,ky) = zero,
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c gx(kx,ky) = s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = s(kx,ky)*inverse FFT of (y/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c q(j,k) = complex charge density
c fxy(1,j,k) = x component of complex force/charge,
c fxy(2,j,k) = y component of complex force/charge,
c fxy(3,j,k) = zero,
c ffg(1,j,k) = potential green's function g
c ffg(2,j,k) = finite-size particle shape factor s
c ffg(3,j,k) = x component of electric field green's function gx
c ffg(4,j,k) = y component of electric field green's function gy
c all for fourier mode (j-1,k-1)
c the ffg array is calculated by the subroutine FORMC2
c electric field energy is also calculated and returned in we
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nx1d = second dimension of ffg arrays, must be >= nx+1
c ny1d = third dimension of ffg arrays, must be >= ny+1
      implicit none
      complex q, fxy
      real ffg
      integer nx, ny, nxv, ny2d, nx1d, ny1d
      real we
      dimension q(nxv,ny2d), fxy(3,nxv,ny2d)
      dimension ffg(4,nx1d,ny1d)
c local data
      double precision wp
      integer nx1, ny1, ny22, j, k, k1
      real at1, at2, at3
      complex zt1, zt2, zero
      nx1 = nx + 1
      ny1 = ny + 1
      ny22 = ny + ny + 2
      zero = cmplx(0.,0.)
c calculate force/charge and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      at3 = 1.0
      do 20 k = 2, ny
      k1 = ny22 - k
      at1 = ffg(3,1,k)
      at3 = -at3
      do 10 j = 2, nx
      at1 = -at1
      at2 = at3*ffg(4,j,1)
      zt1 = cmplx(at1,ffg(3,j,k))
      zt2 = cmplx(at2,ffg(4,j,k))
      fxy(1,j,k) = zt1*q(j,k)
      fxy(2,j,k) = zt2*q(j,k)
      fxy(3,j,k) = zero
      fxy(1,j,k1) = zt1*q(j,k1)
      fxy(2,j,k1) = conjg(zt2)*q(j,k1)
      fxy(3,j,k1) = zero
      wp = wp + ffg(1,j,k)*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)
     1))
   10 continue
   20 continue
c mode number kx = 0
      at3 = ffg(4,1,1)
      do 30 k = 2, ny
      at3 = -at3
      zt1 = cmplx(at3,ffg(4,1,k))
      fxy(1,1,k) = ffg(3,1,k)*q(1,k)
      fxy(2,1,k) = zt1*q(1,k)
      fxy(3,1,k) = zero
      wp = wp + ffg(1,1,k)*(q(1,k)*conjg(q(1,k)))
   30 continue
c mode number kx = nx/2
      at3 = ffg(4,nx1,1)
      do 40 k = 2, ny
      k1 = ny22 - k
      at3 = -at3
      zt1 = cmplx(at3,ffg(4,nx1,k))
      fxy(1,1,k1) = ffg(3,nx1,k)*q(1,k1)
      fxy(2,1,k1) = zt1*q(1,k1)
      fxy(3,1,k1) = zero
      wp = wp + ffg(1,nx1,k)*(q(1,k1)*conjg(q(1,k1)))
   40 continue
c mode number ky = 0
      at1 = ffg(3,1,1)
      do 50 j = 2, nx
      at1 = -at1
      zt1 = cmplx(at1,ffg(3,j,1))
      fxy(1,j,1) = zt1*q(j,1)
      fxy(2,j,1) = ffg(4,j,1)*q(j,1)
      fxy(3,j,1) = zero
      wp = wp + ffg(1,j,1)*(q(j,1)*conjg(q(j,1)))
   50 continue
c mode number ky = ny/2
      at1 = ffg(3,1,ny1)
      do 60 j = 2, nx
      at1 = -at1
      zt1 = cmplx(at1,ffg(3,j,ny1))
      fxy(1,j,ny1) = zt1*q(j,ny1)
      fxy(2,j,ny1) = ffg(4,j,ny1)*q(j,ny1)
      fxy(3,j,ny1) = zero
      wp = wp + ffg(1,j,ny1)*(q(j,ny1)*conjg(q(j,ny1)))
   60 continue
c mode numbers ky = 0, kx = 0, nx/2
      fxy(1,1,1) = cmplx(ffg(3,1,1)*real(q(1,1)),ffg(3,nx1,1)*aimag(q(1,
     11)))
      fxy(2,1,1) = cmplx(ffg(4,1,1)*real(q(1,1)),ffg(4,nx1,1)*aimag(q(1,
     11)))
      fxy(3,1,1) = zero
      wp = wp + .5*(ffg(1,1,1)*real(q(1,1))**2 + ffg(1,nx1,1)*aimag(q(1,
     11))**2)
c mode numbers ky = ny/2, kx = 0, nx/2
      fxy(1,1,ny1) = cmplx(ffg(3,1,ny1)*real(q(1,ny1)),ffg(3,nx1,ny1)*ai
     1mag(q(1,ny1)))
      fxy(2,1,ny1) = cmplx(ffg(4,1,ny1)*real(q(1,ny1)),ffg(4,nx1,ny1)*ai
     1mag(q(1,ny1)))
      fxy(3,1,ny1) = zero
      wp = wp + .5*(ffg(1,1,ny1)*real(q(1,ny1))**2 + ffg(1,nx1,ny1)*aima
     1g(q(1,ny1))**2)
      we = 4.0*float(nx*ny)*wp
      return
      end
      subroutine CUPERPC2(cu,ffg,nx,ny,nxv,ny2d,nx1d,ny1d)
c this subroutine calculates the transverse current in fourier space
c with open (vacuum) boundary conditions using hockney's method.
c input: all, output: cu
c approximate flop count is: 99*nxc*nyc + 72*(nxc + nyc)
c and 2*nx*ny divides
c where nxc = nx - 1, nyc = ny - 1
c the transverse current is calculated using an equation analogous to:
c cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
c ffg(3,j,k) = x component of electric field green's function gx
c ffg(4,j,k) = y component of electric field green's function gy
c all for fourier mode (j-1,k-1)
c the ffg array is calculated by the subroutine FORMC2
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nx1d = second dimension of ffg arrays, must be >= nx+1
c ny1d = third dimension of ffg arrays, must be >= ny+1
      implicit none
      complex cu
      real ffg
      integer nx, ny, nxv, ny2d, nx1d, ny1d
      dimension cu(3,nxv,ny2d), ffg(4,nx1d,ny1d)
c local data
      integer nx1, ny1, ny22, j, k, k1
      real at1, at2, at3, at4, at5
      complex zt1, zt2, zt3
      nx1 = nx + 1
      ny1 = ny + 1
      ny22 = ny + ny + 2
c calculate transverse part of current
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      at3 = 1.0
      do 20 k = 2, ny
      k1 = ny22 - k
      at1 = ffg(3,1,k)
      at3 = -at3
      do 10 j = 2, nx
      at1 = -at1
      at2 = at3*ffg(4,j,1)
      zt1 = cmplx(at1,ffg(3,j,k))
      zt2 = cmplx(at2,ffg(4,j,k))
      zt3 = (zt1*cu(1,j,k) + zt2*cu(2,j,k))/(zt1*zt1 + zt2*zt2)
      cu(1,j,k) = cu(1,j,k) - zt1*zt3
      cu(2,j,k) = cu(2,j,k) - zt2*zt3
      zt2 = conjg(zt2)
      zt3 = (zt1*cu(1,j,k1) + zt2*cu(2,j,k1))/(zt1*zt1 + zt2*zt2)
      cu(1,j,k1) = cu(1,j,k1) - zt1*zt3
      cu(2,j,k1) = cu(2,j,k1) - zt2*zt3
   10 continue
   20 continue
c mode number kx = 0
      at3 = ffg(4,1,1)
      do 30 k = 2, ny
      at3 = -at3
      at1 = ffg(3,1,k)
      zt2 = cmplx(at3,ffg(4,1,k))
      zt3 = (at1*cu(1,1,k) + zt2*cu(2,1,k))/(at1*at1 + zt2*zt2)
      cu(1,1,k) = cu(1,1,k) - at1*zt3
      cu(2,1,k) = cu(2,1,k) - zt2*zt3
   30 continue
c mode number kx = nx/2
      at3 = ffg(4,nx1,1)
      do 40 k = 2, ny
      k1 = ny22 - k
      at3 = -at3
      at1 = ffg(3,nx1,k)
      zt2 = cmplx(at3,ffg(4,nx1,k))
      zt3 = (at1*cu(1,1,k1) + zt2*cu(2,1,k1))/(at1*at1 + zt2*zt2)
      cu(1,1,k1) = cu(1,1,k1) - at1*zt3
      cu(2,1,k1) = cu(2,1,k1) - zt2*zt3
   40 continue
c mode number ky = 0
      at1 = ffg(3,1,1)
      do 50 j = 2, nx
      at1 = -at1
      at2 = ffg(4,j,1)
      zt1 = cmplx(at1,ffg(3,j,1))
      zt3 = (zt1*cu(1,j,1) + at2*cu(2,j,1))/(zt1*zt1 + at2*at2)
      cu(1,j,1) = cu(1,j,1) - zt1*zt3
      cu(2,j,1) = cu(2,j,1) - at2*zt3
   50 continue
c mode number ky = ny/2
      at1 = ffg(3,1,ny1)
      do 60 j = 2, nx
      at1 = -at1
      at2 = ffg(4,j,ny1)
      zt1 = cmplx(at1,ffg(3,j,ny1))
      zt3 = (zt1*cu(1,j,ny1) + at2*cu(2,j,ny1))/(zt1*zt1 + at2*at2)
      cu(1,j,ny1) = cu(1,j,ny1) - zt1*zt3
      cu(2,j,ny1) = cu(2,j,ny1) - at2*zt3
   60 continue
c mode numbers ky = 0, kx = 0, nx/2
      at1 = ffg(3,1,1)
      at2 = ffg(4,1,1)
      at3 = (at1*real(cu(1,1,1)) + at2*real(cu(2,1,1)))/(at1*at1 + at2*a
     1t2)
      at1 = at1*at3
      at2 = at2*at3
      at3 = ffg(3,nx1,1)
      at4 = ffg(4,nx1,1)
      at5 = (at3*aimag(cu(1,1,1)) + at4*aimag(cu(2,1,1)))/(at3*at3 + at4
     1*at4)
      cu(1,1,1) = cu(1,1,1) - cmplx(at1,at3*at5)
      cu(2,1,1) = cu(2,1,1) - cmplx(at2,at4*at5)
c mode numbers ky = ny/2, kx = 0, nx/2
      at1 = ffg(3,1,ny1)
      at2 = ffg(4,1,ny1)
      at3 = (at1*real(cu(1,1,ny1)) + at2*real(cu(2,1,ny1)))/(at1*at1 + a
     1t2*at2)
      at1 = at1*at3
      at2 = at2*at3
      at3 = ffg(3,nx1,ny1)
      at4 = ffg(4,nx1,ny1)
      at5 = (at3*aimag(cu(1,1,ny1)) + at4*aimag(cu(2,1,ny1)))/(at3*at3 +
     1at4*at4)
      cu(1,1,ny1) = cu(1,1,ny1) - cmplx(at1,at3*at5)
      cu(2,1,ny1) = cu(2,1,ny1) - cmplx(at2,at4*at5)
      return
      end
      subroutine BPOISC23(cu,bxy,isign,ffg,ci,wm,nx,ny,nxv,ny2d,nx1d,ny1
     1d)
c this subroutine solves 2-1/2d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with open (vacuum) boundary conditions using hockney's method.
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c for isign = -1
c input: cu,ffg,isign,ci,nx,ny,nx1d,ny1d,nx2v,ny2d, output: bxy,wm
c approximate flop count is: 144*nx*ny + 104*(nx + ny)
c for isign = 1
c input: cu,ffg,isign,ci,nx,ny,nx1d,ny1d,nx2v,ny2d, output: bxy,wm
c approximate flop count is: 63*nx*ny + 66*(nx + ny)
c for isign = 2
c input: cu,ffg,isign,ci,nx,ny,nx1d,ny1d,nx2v,ny2d, output: bxy
c approximate flop count is: 12*nx*ny + 12*(nx + ny)
c if isign < 0, the magnetic field is calculated using the equations:
c bx(kx,ky) = -ci*ci*sqrt(-1)*gy(kx,ky)*cuz(kx,ky)*s(kx,ky),
c by(kx,ky) = ci*ci*sqrt(-1)*gx(kx,ky)*cuz(kx,ky)*s(kx,ky),
c bz(kx,ky) = ci*ci*sqrt(-1)*(gy(kx,ky)*cux(kx,ky)-gx(kx,ky)*cuy(kx,ky))
c            *s(kx,ky),
c where kx = pi*j/nx, ky = pi*k/ny, and j,k = fourier mode numbers,
c gx(kx,ky) = s(kx,ky)*inverse FFT of (x/r)*Er
c gy(kx,ky) = s(kx,ky)*inverse FFT of (y/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky) = inverse FFT of the density of a finite-sized particle
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky) = ci*ci*g(kx,ky)*cux(kx,ky)
c by(kx,ky) = ci*ci*g(kx,ky)*cuy(kx,ky)
c bz(kx,ky) = ci*ci*g(kx,ky)*cuz(kx,ky)
c where g(kx,ky) = affp*inverse FFT of potr
c where potr is the potential of a single finite-sized particle
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky) = cux(kx,ky)*s(kx,ky)
c by(kx,ky) = cuy(kx,ky)*s(kx,ky)
c bz(kx,ky) = cuz(kx,ky)*s(kx,ky)
c cu(i,j,k) = complex current density
c bxy(1,j,k) = x component of complex magnetic field
c bxy(2,j,k) = y component of complex magnetic field
c bxy(3,j,k) = z component of complex magnetic field
c all for fourier mode (j-1,k-1)
c ffg(1,j,k) = potential green's function g
c ffg(2,j,k) = finite-size particle shape factor s
c ffg(3,j,k) = x component of electric field green's function gx
c ffg(4,j,k) = y component of electric field green's function gy
c all for fourier mode (j-1,k-1)
c the ffg array is calculated by the subroutine FORMC2
c ci = reciprical of velocity of light
c magnetic field energy is also calculated and returned in wm
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nx1d = second dimension of ffg arrays, must be >= nx+1
c ny1d = third dimension of ffg arrays, must be >= ny+1
      implicit none
      complex cu, bxy
      real ffg
      integer isign, nx, ny, nxv, ny2d, nx1d, ny1d
      real ci, wm
      dimension cu(3,nxv,ny2d), bxy(3,nxv,ny2d)
      dimension ffg(4,nx1d,ny1d)
c local data
      double precision wp
      integer nx1, ny1, ny22, j, k, k1
      real ci2, at1, at2, at3, at4
      complex zt1, zt2
      if (isign.eq.0) return
      nx1 = nx + 1
      ny1 = ny + 1
      ny22 = ny + ny + 2
      ci2 = ci*ci
      if (isign.gt.0) go to 70
c calculate magnetic field and sum field energy
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      at3 = 1.0
      do 20 k = 2, ny
      k1 = ny22 - k
      at1 = ffg(3,1,k)
      at3 = -at3
      do 10 j = 2, nx
      at1 = -at1
      at2 = at3*ffg(4,j,1)
      zt1 = cmplx(at1,ffg(3,j,k))*ci2
      zt2 = cmplx(at2,ffg(4,j,k))*ci2
      bxy(1,j,k) = -zt2*cu(3,j,k)
      bxy(2,j,k) = zt1*cu(3,j,k)
      bxy(3,j,k) = zt2*cu(1,j,k) - zt1*cu(2,j,k)
      zt2 = conjg(zt2)
      bxy(1,j,k1) = -zt2*cu(3,j,k1)
      bxy(2,j,k1) = zt1*cu(3,j,k1)
      bxy(3,j,k1) = zt2*cu(1,j,k1) - zt1*cu(2,j,k1)
      wp = wp + ffg(1,j,k)*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg
     1(cu(2,j,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j
     2,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)
     3))*ci2
   10 continue
   20 continue
c mode number kx = 0
      at3 = ffg(4,1,1)
      do 30 k = 2, ny
      at3 = -at3
      at1 = ffg(3,1,k)*ci2
      zt1 = cmplx(at3,ffg(4,1,k))*ci2
      bxy(1,1,k) = -zt1*cu(3,1,k)
      bxy(2,1,k) = at1*cu(3,1,k)
      bxy(3,1,k) = zt1*cu(1,1,k) - at1*cu(2,1,k)
      wp = wp + ffg(1,1,k)*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg
     1(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))*ci2
   30 continue
c mode number kx = nx/2
      at3 = ffg(4,nx1,1)
      do 40 k = 2, ny
      k1 = ny22 - k
      at3 = -at3
      at1 = ffg(3,nx1,k)*ci2
      zt1 = cmplx(at3,ffg(4,nx1,k))*ci2
      bxy(1,1,k1) = -zt1*cu(3,1,k1)
      bxy(2,1,k1) = at1*cu(3,1,k1)
      bxy(3,1,k1) = zt1*cu(1,1,k1) - at1*cu(2,1,k1)
      wp = wp + ffg(1,nx1,k)*(cu(1,1,k1)*conjg(cu(1,1,k1)) + cu(2,1,k1)*
     1conjg(cu(2,1,k1)) + cu(3,1,k1)*conjg(cu(3,1,k1)))*ci2
   40 continue
c mode number ky = 0
      at1 = ffg(3,1,1)
      do 50 j = 2, nx
      at1 = -at1
      at3 = ffg(4,j,1)*ci2
      zt1 = cmplx(at1,ffg(3,j,1))*ci2
      bxy(1,j,1) = -at3*cu(3,j,1)
      bxy(2,j,1) = zt1*cu(3,j,1)
      bxy(3,j,1) = at3*cu(1,j,1) - zt1*cu(2,j,1)
      wp = wp + ffg(1,j,1)*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg
     1(cu(2,j,1)) + cu(3,j,1)*conjg(cu(3,j,1)))*ci2
   50 continue
c mode number ky = ny/2
      at1 = ffg(3,1,ny1)
      do 60 j = 2, nx
      at1 = -at1
      at3 = ffg(4,j,ny1)*ci2
      zt1 = cmplx(at1,ffg(3,j,ny1))*ci2
      bxy(1,j,ny1) = -at3*cu(3,j,ny1)
      bxy(2,j,ny1) = zt1*cu(3,j,ny1)
      bxy(3,j,ny1) = at3*cu(1,j,ny1) - zt1*cu(2,j,ny1)
      wp = wp + ffg(1,j,ny1)*(cu(1,j,ny1)*conjg(cu(1,j,ny1)) + cu(2,j,ny
     11)*conjg(cu(2,j,ny1)) + cu(3,j,ny1)*conjg(cu(3,j,ny1)))*ci2
   60 continue
c mode numbers ky = 0, kx = 0, nx/2
      at1 = ffg(3,1,1)*ci2
      at2 = ffg(3,nx1,1)*ci2
      at3 = ffg(4,1,1)*ci2
      at4 = ffg(4,nx1,1)*ci2
      bxy(1,1,1) = -cmplx(at3*real(cu(3,1,1)),at4*aimag(cu(3,1,1)))
      bxy(2,1,1) = cmplx(at1*real(cu(3,1,1)),at2*aimag(cu(3,1,1)))
      bxy(3,1,1) = cmplx(at3*real(cu(1,1,1))-at1*real(cu(2,1,1)),at4*aim
     1ag(cu(1,1,1))-at2*aimag(cu(2,1,1)))
      wp = wp + .5*ffg(1,1,1)*(real(cu(1,1,1))**2 + real(cu(2,1,1))**2 +
     1 real(cu(3,1,1))**2)*ci2 + .5*ffg(1,nx1,1)*(aimag(cu(1,1,1))**2 + 
     2aimag(cu(2,1,1))**2 + aimag(cu(3,1,1))**2)*ci2
c mode numbers ky = ny/2, kx = 0, nx/2
      at1 = ffg(3,1,ny1)*ci2
      at2 = ffg(3,nx1,ny1)*ci2
      at3 = ffg(4,1,ny1)*ci2
      at4 = ffg(4,nx1,ny1)*ci2
      bxy(1,1,ny1) = -cmplx(at3*real(cu(3,1,ny1)),at4*aimag(cu(3,1,ny1))
     1)
      bxy(2,1,ny1) = cmplx(at1*real(cu(3,1,ny1)),at2*aimag(cu(3,1,ny1)))
      bxy(3,1,ny1) = cmplx(at3*real(cu(1,1,ny1))-at1*real(cu(2,1,ny1)),a
     1t4*aimag(cu(1,1,ny1))-at2*aimag(cu(2,1,ny1)))
      wp = wp + .5*ffg(1,1,ny1)*(real(cu(1,1,ny1))**2 + real(cu(2,1,ny1)
     1)**2 + real(cu(3,1,ny1))**2)*ci2 + .5*ffg(1,nx1,ny1)*(aimag(cu(1,1
     2,ny1))**2 + aimag(cu(2,1,ny1))**2 + aimag(cu(3,1,ny1))**2)*ci2
      wm = 4.0*float(nx*ny)*wp
      return
c calculate vector potential and sum field energy
   70 if (isign.gt.1) go to 110
      wp = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 90 k = 2, ny
      k1 = ny22 - k
      do 80 j = 2, nx
      at2 = ffg(1,j,k)*ci2
c     at1 = at2*ffg(2,j,k)
      at1 = at2
      bxy(1,j,k) = at2*cu(1,j,k)
      bxy(1,j,k1) = at2*cu(1,j,k1)
      bxy(2,j,k) = at2*cu(2,j,k)
      bxy(2,j,k1) = at2*cu(2,j,k1)
      bxy(3,j,k) = at2*cu(3,j,k)
      bxy(3,j,k1) = at2*cu(3,j,k1)
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k)) + cu(2,j,k)*conjg(cu(2,j
     1,k)) + cu(3,j,k)*conjg(cu(3,j,k)) + cu(1,j,k1)*conjg(cu(1,j,k1)) +
     2cu(2,j,k1)*conjg(cu(2,j,k1)) + cu(3,j,k1)*conjg(cu(3,j,k1)))
   80 continue
c mode number kx = 0
      at2 = ffg(1,1,k)*ci2
c     at1 = at2*ffg(2,1,k)
      at1 = at2
      bxy(1,1,k) = at2*cu(1,1,k)
      bxy(2,1,k) = at2*cu(2,1,k)
      bxy(3,1,k) = at2*cu(3,1,k)
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k)) + cu(2,1,k)*conjg(cu(2,1
     1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
c mode number kx = nx/2
      at2 = ffg(1,nx1,k)*ci2
c     at1 = at2*ffg(2,nx1,k)
      at1 = at2
      bxy(1,1,k1) = at2*cu(1,1,k1)
      bxy(2,1,k1) = at2*cu(2,1,k1)
      bxy(3,1,k1) = at2*cu(3,1,k1)
      wp = wp + at1*(cu(1,1,k1)*conjg(cu(1,1,k1)) + cu(2,1,k1)*conjg(cu(
     12,1,k1)) + cu(3,1,k1)*conjg(cu(3,1,k1)))
   90 continue
c mode numbers ky = 0, ny/2
      do 100 j = 2, nx
      at2 = ffg(1,j,1)*ci2
c     at1 = at2*ffg(2,j,1)
      at1 = at2
      bxy(1,j,1) = at2*cu(1,j,1)
      bxy(2,j,1) = at2*cu(2,j,1)
      bxy(3,j,1) = at2*cu(3,j,1)
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1)) + cu(2,j,1)*conjg(cu(2,j
     1,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
      at2 = ffg(1,j,ny1)*ci2
c     at1 = at2*ffg(2,j,ny1)
      at1 = at2
      bxy(1,j,ny1) = at2*cu(1,j,ny1)
      bxy(2,j,ny1) = at2*cu(2,j,ny1)
      bxy(3,j,ny1) = at2*cu(3,j,ny1)
      wp = wp + at1*(cu(1,j,ny1)*conjg(cu(1,j,ny1)) + cu(2,j,ny1)*conjg(
     1cu(2,j,ny1)) + cu(3,j,ny1)*conjg(cu(3,j,ny1)))
  100 continue
c mode numbers ky = 0, kx = 0, and kx = nx/2
      at2 = ffg(1,1,1)*ci2
c     at1 = at2*ffg(2,1,1)
      at1 = at2
      at4 = ffg(1,nx1,1)*ci2
c     at3 = at4*ffg(2,nx1,1)
      at3 = at4
      bxy(1,1,1) = cmplx(at2*real(cu(1,1,1)),at4*aimag(cu(1,1,1)))
      bxy(2,1,1) = cmplx(at2*real(cu(2,1,1)),at4*aimag(cu(2,1,1)))
      bxy(3,1,1) = cmplx(at2*real(cu(3,1,1)),at4*aimag(cu(3,1,1)))
      wp = wp + .5*at1*(real(cu(1,1,1))**2 + real(cu(2,1,1))**2 + real(c
     1u(3,1,1))**2) + .5*at3*(aimag(cu(1,1,1))**2 + aimag(cu(2,1,1))**2 
     2+ aimag(cu(3,1,1))**2)
c mode numbers ky = ny/2, kx = 0, and kx = nx/2
      at2 = ffg(1,1,ny1)*ci2
c     at1 = at2*ffg(2,1,ny1)
      at1 = at2
      at4 = ffg(1,nx1,ny1)*ci2
c     at3 = at4*ffg(2,nx1,ny1)
      at3 = at4
      bxy(1,1,ny1) = cmplx(at2*real(cu(1,1,ny1)),at4*aimag(cu(1,1,ny1)))
      bxy(2,1,ny1) = cmplx(at2*real(cu(2,1,ny1)),at4*aimag(cu(2,1,ny1)))
      bxy(3,1,ny1) = cmplx(at2*real(cu(3,1,ny1)),at4*aimag(cu(3,1,ny1)))
      wp = wp + .5*at1*(real(cu(1,1,ny1))**2 + real(cu(2,1,ny1))**2 + re
     1al(cu(3,1,ny1))**2) + .5*at3*(aimag(cu(1,1,ny1))**2 + aimag(cu(2,1
     2,ny1))**2 + aimag(cu(3,1,ny1))**2)
      wm = 4.0*float(nx*ny)*wp
      return
c calculate smoothing
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
  110 do 130 k = 2, ny
      k1 = ny22 - k
      do 120 j = 2, nx
      at1 = ffg(2,j,k)
      bxy(1,j,k) = at1*cu(1,j,k)
      bxy(1,j,k1) = at1*cu(1,j,k1)
      bxy(2,j,k) = at1*cu(2,j,k)
      bxy(2,j,k1) = at1*cu(2,j,k1)
      bxy(3,j,k) = at1*cu(3,j,k)
      bxy(3,j,k1) = at1*cu(3,j,k1)
  120 continue
c mode number kx = 0
      at1 = ffg(2,1,k)
      bxy(1,1,k) = at1*cu(1,1,k)
      bxy(2,1,k) = at1*cu(2,1,k)
      bxy(3,1,k) = at1*cu(3,1,k)
c mode number kx = nx/2
      at1 = ffg(2,nx1,k)
      bxy(1,1,k1) = at1*cu(1,1,k1)
      bxy(2,1,k1) = at1*cu(2,1,k1)
      bxy(3,1,k1) = at1*cu(3,1,k1)
  130 continue
c mode numbers ky = 0, ny/2
      do 140 j = 2, nx
      at1 = ffg(2,j,1)
      bxy(1,j,1) = at1*cu(1,j,1)
      bxy(2,j,1) = at1*cu(2,j,1)
      bxy(3,j,1) = at1*cu(3,j,1)
      at1 = ffg(2,j,ny1)
      bxy(1,j,ny1) = at1*cu(1,j,ny1)
      bxy(2,j,ny1) = at1*cu(2,j,ny1)
      bxy(3,j,ny1) = at1*cu(3,j,ny1)
  140 continue
c mode numbers ky = 0, kx = 0, and kx = nx/2
      at1 = ffg(2,1,1)
      at3 = ffg(2,nx1,1)
      bxy(1,1,1) = cmplx(at1*real(cu(1,1,1)),at3*aimag(cu(1,1,1)))
      bxy(2,1,1) = cmplx(at1*real(cu(2,1,1)),at3*aimag(cu(2,1,1)))
      bxy(3,1,1) = cmplx(at1*real(cu(3,1,1)),at3*aimag(cu(3,1,1)))
c mode numbers ky = ny/2, kx = 0, and kx = nx/2
      at1 = ffg(2,1,ny1)
      at3 = ffg(2,nx1,ny1)
      bxy(1,1,ny1) = cmplx(at1*real(cu(1,1,ny1)),at3*aimag(cu(1,1,ny1)))
      bxy(2,1,ny1) = cmplx(at1*real(cu(2,1,ny1)),at3*aimag(cu(2,1,ny1)))
      bxy(3,1,ny1) = cmplx(at1*real(cu(3,1,ny1)),at3*aimag(cu(3,1,ny1)))
      return
      end
      function POTC3(r,affp,ari,ifun)
c this function calculates the fields for finite-size gaussian particles
c in 3D:
c if ifun = 1, calculate potential function
c POTC3 = (affp/(4*pi))*erfn(r/(ar*sqrt(2.)))/r, for r > 0.
c POTC3 = (affp/(4*pi))*sqrt(2./3.14159265358979)/ar, for r = 0.
c if ifun = 2, calculate particle shape function
c POTC3 = exp(-(r/(sqrt(2.)*ar))**2)/(sqrt(2.*pi)*ar)**3, for r > 0.
c POTC3 = 1./(sqrt(2.*pi)*ar)**3, for r = 0.
c if ifun = 3, calculate radial electric field
c POTC3 = (affp/(4*pi))*(1/r)*(erf(r/(sqrt(2.)*ar))/r -
c exp(-(r/(sqrt(2.)*ar))**2)*sqrt(2./3.14159265358979)/ar, for r > 0.
c POTC3 = 0.0, for r = 0.
c where erfn is the error function
c and where the finite-size particle density is given by:
c rho(r) = exp(-(r/sqrt(2)*ar)**2)/(sqrt(2*pi)*ar)**3
c affp = 4*pi*e**2/(me*(omega0**2)*delta**3) = 1/(n0*delta**3)
c where n0*delta**3 = number density per grid
c r = radial coordinate
c affp = normalization constant
c ari = 1/ar = inverse of particle size function
c (ari = 0., means use point particle result)
c ifun = (1,2,3) = calculate (potential,shape,electric field)
      implicit none
      real r, affp, ari
      integer ifun
c local data
c pi4i = 1/4*pi, sqt2i = 1./sqrt(2.), sqt2pi = sqrt(2./pi)
      real pi4i, sqt2i, sqt2pi
      parameter(pi4i=0.5/6.28318530717959)
      parameter(sqt2i=0.707106781186548,sqt2pi=0.797884560802865)
      real POTC3, erfn
      external erfn
      real anorm, at1, ri
      anorm = affp*pi4i
c calculate potential function
      if (ifun.eq.1) then
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = anorm*sqt2pi*ari
            else
               POTC3 = anorm*erfn(r*sqt2i*ari)/r
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               POTC3 = anorm/r
            endif
         endif
c calculate particle shape function
      else if (ifun.eq.2) then
         anorm = affp*(.5*sqt2pi*ari)**3
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = anorm
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC3 = anorm*exp(-(at1*at1))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = affp
            else
               POTC3 = 0.0
            endif
         endif
c calculate radial electric field
      else if (ifun.eq.3) then
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               ri = 1.0/r
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC3 = anorm*ri*(erfn(at1)*ri - sqt2pi*ari*exp(-(at1*at1
     1)))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               POTC3 = anorm/(r*r)
            endif
         endif
      endif
      return
      end
      function POTC2(r,affp,ari,ifun)
c this function calculates the fields for finite-size gaussian particles
c in 2D:
c if ifun = 1, calculate potential function
c POTC2 = -(affp/(4*pi))*(e1(r**2/(2*ar**2)) + ln(r**2)), for r > 0.
c POTC2 = -(affp/(4*pi))*(ln(2) - gamma + 2*ln(ar), for r = 0.
c if ifun = 2, calculate particle shape function
c POTC2 = exp(-(r/(sqrt(2.)*ar))**2)/(sqrt(2.*pi)*ar)**2, for r > 0.
c POTC2 = 1./(sqrt(2.*pi)*ar)**2, for r = 0.
c if ifun = 3, calculate radial electric field
c POTC2 = 2*(1 - exp(-(r/(sqrt(2.)*ar))**2)/r, for r > 0.
c POTC2 = 0.0, for r = 0.
c where e1 is the exponential integral
c and where the finite-size particle density is given by:
c rho(r) = exp(-(r/sqrt(2)*ar)**2)/(2*pi*ar**2), qm = q/e
c affp = 4*pi*e**2/(me*(omega0**2)*delta**2) = 1/(n0*delta**2)
c where n0*delta**2 = number density per grid
c r = radial coordinate
c affp = normalization constant
c ari = 1/ar = inverse of particle size function
c (ari = 0., means use point particle result)
c ifun = (1,2,3) = calculate (potential,shape,electric field)
      implicit none
      real r, affp, ari
      integer ifun
c local data
c pi4i = 1/4*pi, sqt2i = 1./sqrt(2.), sqt2pi = sqrt(2./pi)
      real pi4i, sqt2i, sqt2pi
      parameter(pi4i=0.5/6.28318530717959)
      parameter(sqt2i=0.707106781186548,sqt2pi=0.797884560802865)
      real POTC2, e1ln
      external e1ln
      real anorm, at1
c calculate potential function
      if (ifun.eq.1) then
         anorm = -affp*pi4i
c finite-size particles
         if (ari.gt.0.) then
            POTC2 = anorm*(e1ln((r*sqt2i*ari)**2) - 2.0*alog(sqt2i*ari))
c point particles
         else
            if (r.eq.0.) then
               POTC2 = 0.0
            else
               POTC2 = 2.0*anorm*alog(r)
            endif
         endif
c calculate particle shape function
      else if (ifun.eq.2) then
         anorm = affp*(.5*sqt2pi*ari)**2
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC2 = anorm
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC2 = anorm*exp(-(at1*at1))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC2 = affp
            else
               POTC2 = 0.0
            endif
         endif
c calculate radial electric field
      else if (ifun.eq.3) then
         anorm = 2.*affp*pi4i
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC2 = 0.0
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC2 = anorm*(1.0 - exp(-(at1*at1)))/r
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC2 = 0.0
            else
               POTC2 = anorm/r
            endif
         endif
      endif
      return
      end
