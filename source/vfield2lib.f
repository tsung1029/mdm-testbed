c 2d PIC library for solving field equations with vacuum boundary
c conditions
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: may 5, 2004
      subroutine NCGUARD2(fxy,bv,nx,ny,nxe,nye)
c this subroutine replicates extended field, bounded in x,
c periodic in y, with quadratic interpolation
c first the boundary values of the fields which are stored in bv are
c copied to the coordinate j=nx+2.  then guard cells are added in the
c x direction at j=1,nx+3 to disable quadratic interpolation within half
c a cell of the the edges and reduce it to linear interpolation.
c finally, in the y direction, periodicity is assumed.
c bv(k,5) = ex(x=Lx,y), bv(k,6) = ey(x=Lx,y)
c nx/ny = system length in x/y direction
c nxe = first dimension of output array fxy, must be >= nx+3
c nye = first dimension of output array fxy, must be >= ny+3
      implicit none
      real fxy, bv
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye), bv(ny,6)
c local data
      integer i, j, k, nx3
      nx3 = nx + 3
      do 20 k = 1, ny
      do 10 i = 1, 2
      fxy(i,1,k+1) = 2.*fxy(i,2,k+1) - fxy(i,3,k+1)
      fxy(i,nx+2,k+1) = bv(k,i+4)
      fxy(i,nx+3,k+1) = 2.*bv(k,i+4) - fxy(i,nx+1,k+1)
   10 continue
   20 continue
      do 40 j = 1, nx3
      do 30 i = 1, 2
      fxy(i,j,1) = fxy(i,j,ny+1)
      fxy(i,j,ny+2) = fxy(i,j,2)
      fxy(i,j,ny+3) = fxy(i,j,3)
   30 continue
   40 continue
      return
      end
      subroutine NDGUARD2(q,bv,nx,ny,nxe,nye)
c this subroutine replicates extended scalar field, bounded in x,
c periodic in y, with quadratic interpolation
c first the boundary values of the field which is stored in bv is
c copied to the coordinate j=nx+2.  then guard cells are added in the
c x direction at j=1,nx+3 to disable quadratic interpolation within half
c a cell of the the edges and reduce it to linear interpolation.
c finally, in the y direction, periodicity is assumed.
c bv(k,6) = phi(x=Lx,y)
c nx/ny = system length in x/y direction
c nxe = first dimension of output array q, must be >= nx+3
c nye = first dimension of output array q, must be >= ny+3
      implicit none
      real q, bv
      integer nx, ny, nxe, nye
      dimension q(nxe,nye), bv(ny,6)
c local data
      integer j, k, nx3
      nx3 = nx + 3
      do 10 k = 1, ny
      q(nx+2,k+1) = bv(k,6)
      q(1,k+1) = 2.*q(2,k+1) - q(3,k+1)
      q(nx+3,k+1) = 2.*bv(k,6) - q(nx+1,k+1)
   10 continue
      do 20 j = 1, nx3
      q(j,1) = q(j,ny+1)
      q(j,ny+2) = q(j,2)
      q(j,ny+3) = q(j,3)
   20 continue
      return
      end
      subroutine NCGUARD2L(fxy,bv,nx,ny,nxe,nye)
c this subroutine replicates extended field, bounded in x,
c periodic in y, with linear interpolation
c first the boundary values of the fields which are stored in bv are
c copied to the coordinate j=nx+1.
c in the y direction, periodicity is assumed.
c bv(k,5) = ex(x=Lx,y), bv(k,6) = ey(x=Lx,y)
c nx/ny = system length in x/y direction
c nxe = first dimension of output array fxy, must be >= nx+1
c nye = first dimension of output array fxy, must be >= ny+1
      implicit none
      real fxy, bv
      integer nx, ny, nxe, nye
      dimension fxy(2,nxe,nye), bv(ny,6)
c local data
      integer j, k, nx1
      nx1 = nx + 1
      do 10 k = 1, ny
      fxy(1,nx+1,k) = bv(k,5)
      fxy(2,nx+1,k) = bv(k,6)
   10 continue
      do 20 j = 1, nx1
      fxy(1,j,ny+1) = fxy(1,j,1)
      fxy(2,j,ny+1) = fxy(2,j,1)
   20 continue
      return
      end
      subroutine NDGUARD2L(q,bv,nx,ny,nxe,nye)
c this subroutine replicates extended scalar field, bounded in x,
c periodic in y, with linear interpolation
c first the boundary values of the field which is stored in bv is
c copied to the coordinate j=nx+1.
c in the y direction, periodicity is assumed.
c bv(k,5) = ex(x=Lx,y), bv(k,6) = ey(x=Lx,y)
c nx/ny = system length in x/y direction
c nxe = first dimension of output array q, must be >= nx+1
c nye = first dimension of output array q, must be >= ny+1
      implicit none
      real q, bv
      integer nx, ny, nxe, nye
      dimension q(nxe,nye), bv(ny,6)
c local data
      integer j, k, nx1
      nx1 = nx + 1
      do 10 k = 1, ny
      q(nx+1,k) = bv(k,6)
   10 continue
      do 20 j = 1, nx1
      q(j,ny+1) = q(j,1)
   20 continue
      return
      end
      subroutine BNDRYV2(q,ffc,bv,nx,ny,nxv,nxd,nyd,nyhd)
c this subroutine calculates the boundary values of electric field of
c the periodic solution of poisson's equation in Fourier space from the
c charge density.  The results are used in calculating the solution of a
c laplacian in order to satisfy non-periodic boundary conditions.
c algorithm used in described in V. K. Decyk and J. M. Dawson,
c Journal of Computational Physics 30, 407 (1979).
c input: q, ffp, nx, ny, nxvh, nxhd, nyhd, output: bv
c approximate flop count = 16*nxc*nyc + 8*nyc + 7*nxc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c q(j,k) = input complex charge density,
c real(ffc(j,k)) = finite-size particle shape factor s,
c aimag(ffc(j,k)) = potential green's function g,
c all for for fourier mode (j-1,k-1)
c bv = boundary fields, bv(k,3) = KmPm and bv(k,4) = PIm, except
c imag(PI0) = net charge density rho
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nxd/nyd = must be >= nx/ny
c nyhd = must be >= ny/2
      implicit none
      integer nx, ny, nxv, nxd, nyd, nyhd
      real q, ffc, bv
      dimension q(nxv,ny), ffc(nxd,nyhd), bv(nyd,4)
c local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, sum1, sum2, sum3, sum4, at1, at2
      nxh = nx/2
      nyh = ny/2
      ny2 = ny + 2
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
c calculate KmPm and PIm
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
      do 10 j = 2, nxh
      at2 = ffc(2*j-1,k)
      at1 = dky*at2
      at2 = dnx*float(j - 1)*at2
      sum1 = sum1 + at1*(q(2*j-1,k) + q(2*j-1,k1))
      sum2 = sum2 + at1*(q(2*j,k) - q(2*j,k1))
      sum3 = sum3 + at2*(q(2*j,k) + q(2*j,k1))
      sum4 = sum4 + at2*(q(2*j-1,k1) - q(2*j-1,k))
   10 continue
      bv(2*k-1,3) = sum1 + dky*ffc(1,k)*q(1,k)
      bv(2*k,3) = sum2 + dky*ffc(1,k)*q(2,k)
      bv(2*k-1,4) = sum3
      bv(2*k,4) = sum4
   20 continue
c calculate P0 and PI0
      sum1 = 0.
      sum2 = 0.
      do 30 j = 2, nxh
      at1 = ffc(2*j-1,1)
      sum1 = sum1 + at1*q(2*j-1,1)
      sum2 = sum2 + (dnx*float(j - 1)*at1)*q(2*j,1)
   30 continue
      bv(1,3) = sum1 + sum1
      bv(2,3) = 0.
c imaginary part of bv(1,4) contains net charge rho00
      bv(1,4) = sum2 + sum2
      bv(2,4) = ffc(1,1)*q(1,1)
      return
      end
      subroutine POISB2(fx,fy,isign,ffb,bv,bcd,mixup,sct,t,we,affp,indx,
     1ny,nxv,nxd,nxhd,nyd,nyhd)
c this subroutine finds corrections to 2d poisson's equation for 
c force/charge or potential with vacuum boundary conditions and 
c external surface charge in one direction, and periodic in the other.
c average potential across system is zero.
c a periodic solution is assumed to have been found first with poisp2,
c and boundary values with bndryv2
c algorithm used in described in V. K. Decyk and J. M. Dawson,
c Journal of Computational Physics 30, 407 (1979).
c for isign = 0, input: isign,indx,ny,nxv,nxd,nxhd,nyd,nyhd
c                output: ffb,bcd
c                scratch: mixup,sct,t
c for isign = -1, input:  isign,fx,fy,ffb,bv,bcd,affp,indx,nxh,ny
c                         nxv,nxd,nxhd,nyd,nyhd
c                 output: fx, fy, bv, we
c approximate flop count is: 48*nxc*nyc + 88*nyc + 3*nxc + nyc divides
c for isign = 1, input:  isign,fx,ffb,bv,bcd,affp,indx,ny
c                        nxv,nxd,nxhd,nyd,nyhd
c                output: fx, bv, we
c approximate flop count is: 24*nxc*nyc + + 74*nyc + 6*nxc + nyc divides
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c if isign < 0, the force/charge correction is calculated:
c exc(km,x) = -Am*exp(-km*(Lx-x)) + Bm*exp(-km*x)
c eyc(km,x) = -sqrt(-1)*(Am*exp(-km*(Lx-x)) + Bm*exp(-km*x))
c exc(k0,x) = -4*pi*rho00*(Lx/2-x) - A0
c eyc(k0,x) = 0.
c where Am = .5*(4*pi*sigma(x=Lx,k) + PIm - km*Pm),
c       Bm = .5*(4*pi*sigma(x=0,k) - PIm - km*Pm),
c and   A0 = 2*pi*sigma(x=Lx) - 2*pi*sigma(x=0) + PI0
c where PIm and Pm are the periodic ex and phi at the boundaries
c the calculations are done in fourier space and are added to the
c periodic forces already in fx, fy
c on output, bv = value of electric fields on right boundary:
c bv(k,5) = ex(x=Lx), bv(k,6) = ey(x=Lx)
c if isign = 1, potential correction is calculated:
c potc(km,x) = (Am*exp(-km*(Lx-x)) + Bm*exp(-km*x)/km
c potc(k0,x) = 2*pi*rho00*x*(Lx-x) - A0*(Lx/2-x) - P0
c the calculation is done in fourier space and is added to the
c periodic potential already in fx.
c on output, bv = value of potential on right boundary:
c bv(k,6) = phi(x=Lx)
c if isign = 0, form factor arrays ffb and bcd are prepared
c on input, fx and/or fy contain periodic part of solution
c on output, fx and/or fy contain total solution
c ffb(j,k) = (1/nx)*inverse fft(exp(-dky*float(nx + 1 - j))))
c real(ffb(j,1)) = (1/nx)*inverse fft((j - 1)*(nx + 1 - j))))
c aimag(ffb(j,1)) = (1/nx)*inverse fft((nx/2 + 1 - j)))
c on input, bv = input surface charge and boundary values
c for fourier mode k-1:
c bv(k,1) = 4*pi*sigma(x=0), bv(k,2) = 4*pi*sigma(x=Lx)
c bv(k,3) = KmPm, bv(k,4) = PIm
c both are normalized in the same way as the electric field.
c bcd(k) = exp(-ky*Lx)
c mixup = array of bit reversed addresses for fft
c sct = sine/cosine table for fft
c t = complex scratch array, used during initialiation of fft tables
c we = bounded corrections to periodic electric field energy
c affp = normalization constant for poisson's equation
c indx = exponent which determines length in x direction, where nx=2**indx
c ny = system length in y direction
c nxv = first dimension of field arrays, must be >= nx
c nxd = must be >= nx
c nxhd = must be >= nx/2
c nyd = must be >= ny
c nyhd = must be >= ny/2
      implicit none
      complex sct
      integer isign, mixup, indx, ny, nxv, nxd, nxhd, nyd, nyhd
      real fx, fy, ffb, bv, bcd, we, affp, t
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension ffb(nxd,nyhd), bv(nyd,6), bcd(nyhd)
      dimension mixup(nxhd), sct(nxhd), t(nxd)
c local data
      double precision wp, wb
      real cr, ci, dr, di, t1r, t1i, t2r, t2i, t3r, t3i, t4r, t4i
      integer nx, nxh, nx1, nyh, ny2
      integer is, j, j1, k, k1
      real dny, anx, anxi, dky, at1, at2, at3, at4, rho, rholx, dkyi
      real sum1, sum2, sum3, sum4
      nx = 2**indx
      nxh = nx/2
      nx1 = nx + 1
      nyh = ny/2
      ny2 = ny + 2
      dny = 6.28318530717959/float(ny)
      anx = float(nx)
c initialization
      if (isign.ne.0) go to 50
c prepare fft tables
      is = 0
      call FFT1RX(ffb,t,is,mixup,sct,indx,nx,nxh)
      is = -1
c prepare form factor array
      do 10 j = 1, nx
      j1 = j - 1
      ffb(j,1) = float(j1*(nx - j1))
      ffb(j,2) = float(nxh - j1)
   10 continue
      call FFT1RX(ffb(1,1),t,is,mixup,sct,indx,nx,nxh)
      call FFT1RX(ffb(1,2),t,is,mixup,sct,indx,nx,nxh)
      do 20 j = 1, nxh
      ffb(2*j,1) = ffb(2*j,2)
   20 continue
      do 40 k = 2, nyh
      dky = dny*float(k - 1)
      do 30 j = 1, nx
      ffb(j,k) = exp(-amin1(50.,dky*float(nx1 - j)))
   30 continue
      bcd(k) = exp(-amin1(50.,dky*anx))
      call FFT1RX(ffb(1,k),t,is,mixup,sct,indx,nx,nxh)
   40 continue
      return
   50 if (isign.gt.0) go to 90
c calculate force/charge and sum field energy
      anxi = 1./anx
      wp = 0.0d0
      wb = 0.0d0
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
c find constants for solution of homogeneous equation
      cr = .5*(bv(2*k-1,1) - bv(2*k-1,3) - bv(2*k-1,4))
      ci = .5*(bv(2*k,1) - bv(2*k,3) - bv(2*k,4))
      dr = .5*(bv(2*k-1,2) - bv(2*k-1,3) + bv(2*k-1,4))
      di = .5*(bv(2*k,2) - bv(2*k,3) + bv(2*k,4))
      t1r = dr - cr
      t1i = di - ci
      t2r = dr + cr
      t2i = di + ci
c boundary fields
      t3r = cr + dr*bcd(k)
      t3i = ci + di*bcd(k)
      t4r = cr*bcd(k) + dr
      t4i = ci*bcd(k) + di
c calculate internal and boundary energy corrections
      at2 = anxi/dky
      wp = wp + (t2r*bv(2*k-1,3) + t2i*bv(2*k,3) + t1r*bv(2*k-1,4) + t1i
     1*bv(2*k,4))*at2*(1. - bcd(k))
      wb = wb + (bv(2*k-1,1)*(bv(2*k-1,3) + t3r) + bv(2*k,1)*(bv(2*k,3) 
     1+ t3i) + bv(2*k-1,2)*(bv(2*k-1,3) + t4r) + bv(2*k,2)*(bv(2*k,3) + 
     1t4i))*at2
c homogenous electric field in x direction at x = Lx
      bv(2*k-1,5) = cr*bcd(k) - dr
      bv(2*k,5) = ci*bcd(k) - di
c homogenous electric field in y direction at x = Lx
      bv(2*k-1,6) = t4i
      bv(2*k,6) = -t4r
c homogenous electric field in x direction at x = 0
c     bv(2*k-1,7) = cr - dr*bcd(k)
c     bv(2*k,7) = ci - di*bcd(k)
c homogenous electric field in y direction at x = 0
c     bv(2*k-1,8) = t3i
c     bv(2*k,8) = -t3r
c calculate extra term in homogeneous solution
      cr = cr*(1. - bcd(k))*anxi
      ci = ci*(1. - bcd(k))*anxi
      dr = ci
      di = -cr
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
c add solutions of homogeneous equation to periodic solution
      do 60 j = 2, nxh
      sum1 = sum1 + (fx(2*j-1,k) + fx(2*j-1,k1))
      sum2 = sum2 + (fx(2*j,k) - fx(2*j,k1))
      sum3 = sum3 + (fy(2*j-1,k) + fy(2*j-1,k1))
      sum4 = sum4 + (fy(2*j,k) - fy(2*j,k1))
      fx(2*j-1,k) = fx(2*j-1,k) - (t1r*ffb(2*j-1,k) - t2i*ffb(2*j,k)) + 
     1cr
      fx(2*j,k) = fx(2*j,k) - (t1i*ffb(2*j-1,k) + t2r*ffb(2*j,k)) + ci
      fx(2*j-1,k1) = fx(2*j-1,k1) - (t1r*ffb(2*j-1,k) + t2i*ffb(2*j,k)) 
     1+ cr
      fx(2*j,k1) = fx(2*j,k1) + (t1i*ffb(2*j-1,k) - t2r*ffb(2*j,k)) - ci
      fy(2*j-1,k) = fy(2*j-1,k) + (t2i*ffb(2*j-1,k) + t1r*ffb(2*j,k)) + 
     1dr
      fy(2*j,k) = fy(2*j,k) - (t2r*ffb(2*j-1,k) - t1i*ffb(2*j,k)) + di
      fy(2*j-1,k1) = fy(2*j-1,k1) + (t2i*real(ffb(2*j-1,k)) - t1r*ffb(2*
     1j,k)) + dr
      fy(2*j,k1) = fy(2*j,k1) + (t2r*real(ffb(2*j-1,k)) + t1i*ffb(2*j,k)
     1) - di
   60 continue
c modes with n = 0, nx/2 are special
      sum1 = sum1 + fx(1,k)
      sum2 = sum2 + fx(2,k)
      sum3 = sum3 + fy(1,k)
      sum4 = sum4 + fy(2,k)
      fx(1,k) = -t1r*ffb(1,k) + cr
      fx(2,k) = -t1i*ffb(1,k) + ci
      fx(1,k1) = -t1r*ffb(2,k) + cr
      fx(2,k1) = t1i*ffb(2,k) - ci
      fy(1,k) = fy(1,k) + t2i*ffb(1,k) + dr
      fy(2,k) = fy(2,k) - t2r*ffb(1,k) + di
      fy(1,k1) = t2i*ffb(2,k) + dr
      fy(2,k1) = t2r*ffb(2,k) - di
c electric field in x direction at x = Lx
      bv(2*k-1,5) = bv(2*k-1,5) + sum1
      bv(2*k,5) = bv(2*k,5) + sum2
c electric field in y direction at x = Lx
      bv(2*k-1,6) = bv(2*k-1,6) + sum3
      bv(2*k,6) = bv(2*k,6) + sum4
c electric field in x direction at x = 0
c     bv(2*k-1,7) = bv(2*k-1,7) + sum1
c     bv(2*k,7) = bv(2*k,7) + sum2
c electric field in y direction at x = 0
c     bv(2*k-1,8) = bv(2*k-1,8) + sum3
c     bv(2*k,8) = bv(2*k,8) + sum4
   70 continue
c find constants for solution of homogeneous equation
      rho = bv(2,4)
      rholx = .5*rho*anx
c find constants for solution of homogeneous equation
      at1 = rho*ffb(2,1)
      at2 = -(.5*(bv(1,2) - bv(1,1)) + bv(1,4))
      at3 = .5*at2*anx
c calculate energies
      wp = wp - .5*(at2*bv(1,4) - rho*(rholx*anx/6. - 2.*bv(1,3)))
      wb = wb - .5*(bv(1,2) - bv(1,1))*at2
      we = anx*float(ny)*(wp + wb)/affp
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
c add solution of homogeneous equation to periodic solution
      do 80 j = 2, nxh
      sum1 = sum1 + fx(2*j-1,1)
      sum2 = sum2 + fy(2*j-1,1)
      if (rho.ne.0.) then
         fx(2*j-1,1) = fx(2*j-1,1) - at1
         fx(2*j,1) = fx(2*j,1) - rho*ffb(2*j,1)
      endif
   80 continue
      sum1 = sum1 + sum1
      sum2 = sum2 + sum2
      fx(1,1) = at2 - at1
      fx(2,1) = -at1
c electric field in x direction at x = Lx
      bv(1,5) = (at2 + rholx) + sum1
      bv(2,5) = 0.
c electric field in y direction at x = Lx
      bv(1,6) = sum2
      bv(2,6) = 0.
c electric field in x direction at x = 0
c     bv(1,7) = (at2 - rholx) + sum1
c     bv(2,7) = 0.
c electric field in y direction at x = 0
c     bv(1,8) = sum2
c     bv(2,8) = 0.
      return
c calculate potential and sum field energy
   90 anxi = 1./anx
      wp = 0.0d0
      wb = 0.0d0
      do 110 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
      dkyi = .5/dky
c find constants for solution of homogeneous equation
      cr = dkyi*(bv(2*k-1,1) - bv(2*k-1,3) - bv(2*k-1,4))
      ci = dkyi*(bv(2*k,1) - bv(2*k,3) - bv(2*k,4))
      dr = dkyi*(bv(2*k-1,2) - bv(2*k-1,3) + bv(2*k-1,4))
      di = dkyi*(bv(2*k,2) - bv(2*k,3) + bv(2*k,4))
      t1r = dr + cr
      t1i = di + ci
      t2r = dr - cr
      t2i = di - ci
c boundary potentials
      t3r = cr + dr*bcd(k)
      t3i = ci + di*bcd(k)
      t4r = cr*bcd(k) + dr
      t4i = ci*bcd(k) + di
c calculate internal and boundary energy corrections
      at2 = anxi/dky
      wp = wp + (t1r*bv(2*k-1,3) + t1i*bv(2*k,3) + t2r*bv(2*k-1,4) + t2i
     1*bv(2*k,4))*anxi*(1. - bcd(k))
      wb = wb + (bv(2*k-1,1)*(bv(2*k-1,3) + t3r*dky) + bv(2*k,1)*(bv(2*k
     1,3) + t3i*dky) + bv(2*k-1,2)*(bv(2*k-1,3) + t4r*dky) + bv(2*k,2)*(
     1bv(2*k,3) + t4i*dky))*at2
c calculate extra term in homogeneous solution
      cr = cr*(1. - bcd(k))*anxi
      ci = ci*(1. - bcd(k))*anxi
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
c add solutions of homogeneous equation to periodic solution
      do 100 j = 2, nxh
      sum1 = sum1 + (fx(2*j-1,k) + fx(2*j-1,k1))
      sum2 = sum2 + (fx(2*j,k) - fx(2*j,k1))
      fx(2*j-1,k) = fx(2*j-1,k) + t1r*ffb(2*j-1,k) - t2i*ffb(2*j,k) + cr
      fx(2*j,k) = fx(2*j,k) + t1i*ffb(2*j-1,k) + t2r*ffb(2*j,k) + ci
      fx(2*j-1,k1) = fx(2*j-1,k1) + t1r*ffb(2*j-1,k) + t2i*ffb(2*j,k) + 
     1cr
      fx(2*j,k1) = fx(2*j,k1) - t1i*ffb(2*j-1,k) + t2r*ffb(2*j,k) - ci
  100 continue
c modes with n = 0, nx/2 are special
      sum1 = sum1 + fx(1,k)
      sum2 = sum2 + fx(2,k)
      fx(1,k) = fx(1,k) + t1r*ffb(1,k) + cr
      fx(2,k) = fx(2,k) + t1i*ffb(1,k) + ci
      fx(1,k1) = t1r*ffb(2,k) + cr
      fx(2,k1) = -t1i*ffb(2,k) - ci
c potential at x = Lx
      bv(2*k-1,6) = t4r + sum1
      bv(2*k,6) = t4i + sum2
c potential at x = 0
c     bv(2*k-1,5) = t3r + sum1
c     bv(2*k,5) = t3i + sum2
  110 continue
c find constants for solution of homogeneous equation
      rho = bv(2,4)
      rholx = .5*rho*anx
c find constants for solution of homogeneous equation
      at2 = -(.5*(bv(1,2) - bv(1,1)) + bv(1,4))
      at3 = .5*at2*anx
      at1 = at2*ffb(2,1)
      at4 = -bv(1,3)
c calculate energies
      wp = wp - .5*(at2*bv(1,4) - rho*(rholx*anx/6. + 2.*at4))
      wb = wb - .5*(bv(1,2) - bv(1,1))*at2
      we = anx*float(ny)*(wp + wb)/affp
c homogenous potential at x = Lx
      bv(1,6) = -at3
      bv(2,6) = 0.
c homogenous potential at x = 0
c     bv(1,5) = at3
c     bv(2,5) = 0.
      at3 = .5*rho
c find boundary values of periodic solution
      sum1 = 0.
c add solution of homogeneous equation to periodic solution
      do 120 j = 2, nxh
      sum1 = sum1 + fx(2*j-1,1)
      fx(2*j-1,1) = fx(2*j-1,1) + at1 + at3*ffb(2*j-1,1)
      fx(2*j,1) = fx(2*j,1) + at2*ffb(2*j,1)
  120 continue
      sum1 = sum1 + sum1
      fx(1,1) = at3*ffb(1,1) + at1 + at4
      fx(2,1) = -at3*ffb(2,1) + at1
c potential at x = Lx
      bv(1,6) = bv(1,6) + sum1 + at4
c potential at x = 0
c     bv(1,5) = bv(1,5) + sum1 + at4
      return
      end
      subroutine POISB22(fxy,isign,ffb,bv,bcd,mixup,sct,t,we,affp,indx,n
     1y,nxvh,nxhd,nyhd)
c this subroutine finds corrections to 2d poisson's equation for 
c force/charge or potential with vacuum boundary conditions and 
c external surface charge in one direction, and periodic in the other.
c average potential across system is zero.
c a periodic solution is assumed to have been found first with poispp2,
c and boundary values with bndryv2
c algorithm used in described in V. K. Decyk and J. M. Dawson,
c Journal of Computational Physics 30, 407 (1979).
c for isign = 0, input: isign,indx,nxh,ny,nyh,nxvh
c                output: ffb,bcd
c                scratch: mixup,sct,t
c for isign /= 0, input:  isign,fxy,ffb,bv,bcd,affp,indx,ny,
c                         nxvh,nxhd,nyhd
c                 output: fxy, bv, we
c approximate flop count is: 54*nxc*nyc + + 115*nyc + 3*nxc
c where nxc = nx/2 - 1, nyc = ny/2 - 1
c equations used are:
c exc(km,x) = -Am*exp(-km*(Lx-x)) + Bm*exp(-km*x)
c eyc(km,x) = -sqrt(-1)*(Am*exp(-km*(Lx-x)) + Bm*exp(-km*x))
c exc(k0,x) = -4*pi*rho00*(Lx/2-x) - A0
c eyc(k0,x) = 0.
c where Am = .5*(4*pi*sigma(x=Lx,k) + PIm - km*Pm),
c       Bm = .5*(4*pi*sigma(x=0,k) - PIm - km*Pm),
c and   A0 = 2*pi*sigma(x=Lx) - 2*pi*sigma(x=0) + PI0
c where PIm and Pm are the periodic ex and phi at the boundaries
c the calculations are done in fourier space and are added to the
c periodic forces already in fxy
c on output, bv = value of electric fields on right boundary:
c bv(k,5) = ex(x=Lx), bv(k,6) = ey(x=Lx)
c if isign = 0, form factor arrays ffb and bcd are prepared
c if isign is not equal to 0, force/charge is calculated
c on input, fxy contain periodic part of solution
c on output, fxy contain total solution
c ffb(j,k) = (1/nx)*inverse fft(exp(-dky*float(nx + 1 - j))))
c real(ffb(j,1)) = (1/nx)*inverse fft((j - 1)*(nx + 1 - j))))
c aimag(ffb(j,1)) = (1/nx)*inverse fft((nx/2 + 1 - j)))
c on input, bv = input surface charge and boundary values
c for fourier mode k-1:
c bv(k,1) = 4*pi*sigma(x=0), bv(k,2) = 4*pi*sigma(x=Lx)
c bv(k,3) = KmPm, bv(k,4) = PIm
c both are normalized in the same way as the electric field.
c bcd(k) = exp(-ky*Lx)
c mixup = array of bit reversed addresses for fft
c sct = sine/cosine table for fft
c t = complex scratch array, used during initialiation of fft tables
c we = bounded corrections to periodic electric field energy
c affp = normalization constant for poisson's equation
c indx = exponent which determines length in x direction, where nx=2**indx
c ny = system length in y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxhd = must be >= nx/2
c nyhd = must be >= ny/2
      implicit none
      complex fxy, ffb, bv, sct, t
      integer isign, mixup, indx, ny, nxvh, nxhd, nyhd
      real bcd, we, affp
      dimension fxy(2,nxvh,ny)
      dimension ffb(nxhd,nyhd), bv(nyhd,6), bcd(nyhd)
      dimension mixup(nxhd), sct(nxhd), t(nxhd)
c local data
      double precision wp, wb
      complex zc, zd, zt1, zt2, zt3, zt4
      integer nx, nxh, nx1, nyh, ny2
      integer is, j, j1, j2, j3, k, k1
      real dny, anx, anxi, dky, at1, at2, at3, rho, rholx
      real sum1, sum2, sum3, sum4
      nx = 2**indx
      nxh = nx/2
      nx1 = nx + 1
      nyh = ny/2
      ny2 = ny + 2
      dny = 6.28318530717959/float(ny)
      anx = float(nx)
c initialization
      if (isign.ne.0) go to 50
c prepare fft tables
      is = 0
      call FFT1RX(ffb,t,is,mixup,sct,indx,nx,nxh)
      is = -1
c prepare form factor array
      do 10 j = 1, nxh
      j1 = j - 1
      j2 = j1 + j1
      j3 = j2 + 1
      ffb(j,1) = cmplx(float(j2*(nx - j2)),float(j3*(nx - j3)))
      ffb(j,2) = cmplx(float(nxh - j2),float(nxh - j3))
   10 continue
      call FFT1RX(ffb(1,1),t,is,mixup,sct,indx,nx,nxh)
      call FFT1RX(ffb(1,2),t,is,mixup,sct,indx,nx,nxh)
      do 20 j = 1, nxh
      ffb(j,1) = cmplx(real(ffb(j,1)),aimag(ffb(j,2)))
   20 continue
      do 40 k = 2, nyh
      dky = dny*float(k - 1)
      do 30 j = 1, nxh
      j2 = j + j
      j1 = j2 - 1
      ffb(j,k) = cmplx(exp(-amin1(50.,dky*float(nx1 - j1))),exp(-amin1(5
     10.,dky*float(nx1 - j2))))
   30 continue
      bcd(k) = exp(-amin1(50.,dky*anx))
      call FFT1RX(ffb(1,k),t,is,mixup,sct,indx,nx,nxh)
   40 continue
      return
c calculate force/charge and sum field energy
   50 wp = 0.0d0
      wb = 0.0d0
      anxi = 1./anx
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*float(k - 1)
c find constants for solution of homogeneous equation
      zc = .5*(bv(k,1) - bv(k,3) - bv(k,4))
      zd = .5*(bv(k,2) - bv(k,3) + bv(k,4))
      zt1 = zc - zd
      zt2 = cmplx(0.,1.)*conjg(zc + zd)
c boundary fields
      zt3 = zc + zd*bcd(k)
      zt4 = zc*bcd(k) + zd
c calculate internal and boundary energy corrections
      at2 = anxi/dky
      wp = wp + (aimag(zt2*bv(k,3)) - real(conjg(zt1)*bv(k,4)))*at2*(1. 
     1- bcd(k))
      wb = wb + (conjg(bv(k,1))*(bv(k,3) + zt3) + conjg(bv(k,2))*(bv(k,3
     1) + zt4))*at2
c homogenous electric field in x direction at x = Lx
      bv(k,5) = zc*bcd(k) - zd
c homogenous electric field in y direction at x = Lx
      bv(k,6) = -cmplx(0.,1.)*zt4
c homogenous electric field in x direction at x = 0
c     bv(k,7) = zc - zd*bcd(k)
c homogenous electric field in y direction at x = 0
c     bv(k,8) = -cmplx(0.,1.)*zt3
c calculate extra term in homogeneous solution
      zc = zc*(1. - bcd(k))*anxi
      zd = -cmplx(0.,1.)*zc
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
      sum3 = 0.
      sum4 = 0.
c add solutions of homogeneous equation to periodic solution
      do 60 j = 2, nxh
      sum1 = sum1 + real(fxy(1,j,k) + fxy(1,j,k1))
      sum2 = sum2 + aimag(fxy(1,j,k) - fxy(1,j,k1))
      sum3 = sum3 + real(fxy(2,j,k) + fxy(2,j,k1))
      sum4 = sum4 + aimag(fxy(2,j,k) - fxy(2,j,k1))
      fxy(1,j,k) = fxy(1,j,k) + zt1*real(ffb(j,k)) + conjg(zt2)*aimag(ff
     1b(j,k)) + zc
      fxy(1,j,k1) = fxy(1,j,k1) + conjg(zt1)*real(ffb(j,k)) - zt2*aimag(
     1ffb(j,k)) + conjg(zc)
      fxy(2,j,k) = fxy(2,j,k) + conjg(zt2)*real(ffb(j,k)) - zt1*aimag(ff
     1b(j,k)) + zd
      fxy(2,j,k1) = fxy(2,j,k1) + zt2*real(ffb(j,k)) + conjg(zt1)*aimag(
     1ffb(j,k)) + conjg(zd)
   60 continue
c modes with n = 0, nx/2 are special
      sum1 = sum1 + real(fxy(1,1,k))
      sum2 = sum2 + aimag(fxy(1,1,k))
      sum3 = sum3 + real(fxy(2,1,k))
      sum4 = sum4 + aimag(fxy(2,1,k))
      fxy(1,1,k) = zt1*real(ffb(1,k)) + zc
      fxy(1,1,k1) = conjg(zt1)*aimag(ffb(1,k)) + conjg(zc)
      fxy(2,1,k) = fxy(2,1,k) + conjg(zt2)*real(ffb(1,k)) + zd
      fxy(2,1,k1) = zt2*aimag(ffb(1,k)) + conjg(zd)
c electric field in x direction at x = Lx
      bv(k,5) = bv(k,5) + cmplx(sum1,sum2)
c electric field in y direction at x = Lx
      bv(k,6) = bv(k,6) + cmplx(sum3,sum4)
c electric field in x direction at x = 0
c     bv(k,7) = bv(k,7) + cmplx(sum1,sum2)
c electric field in y direction at x = 0
c     bv(k,8) = bv(k,8) + cmplx(sum3,sum4)
   70 continue
c find constants for solution of homogeneous equation
      rho = aimag(bv(1,4))
      rholx = .5*rho*anx
c find constants for solution of homogeneous equation
      at1 = rho*aimag(ffb(1,1))
      at2 = -(.5*(bv(1,2) - bv(1,1)) + bv(1,4))
      at3 = .5*at2*anx
c calculate energies
      wp = wp - .5*(at2*bv(1,4) - rho*(rholx*anx/6. - 2.*real(bv(1,3))))
      wb = wb - .5*(bv(1,2) - bv(1,1))*at2
      we = anx*float(ny)*(wp + wb)/affp
c find boundary values of periodic solution
      sum1 = 0.
      sum2 = 0.
c add solution of homogeneous equation to periodic solution
      do 80 j = 2, nxh
      sum1 = sum1 + real(fxy(1,j,1))
      sum2 = sum2 + real(fxy(2,j,1))
      if (rho.ne.0.) fxy(1,j,1) = fxy(1,j,1) - cmplx(at1,rho*aimag(ffb(j
     1,1)))
   80 continue
      sum1 = sum1 + sum1
      sum2 = sum2 + sum2
      fxy(1,1,1) = cmplx(at2 - at1,-at1)
c electric field in x direction at x = Lx
      bv(1,5) = cmplx(at2+rholx,0.) + cmplx(sum1,0.)
c electric field in y direction at x = Lx
      bv(1,6) = cmplx(sum2,0.)
c electric field in x direction at x = 0
c     bv(1,7) = cmplx(at2-rholx,0.) + cmplx(sum1,0.)
c electric field in y direction at x = 0
c     bv(1,8) = cmplx(sum2,0.)
      return
      end
