c-----------------------------------------------------------------------
c 2d PIC library for initialization
c init2lib.f contains procedures to initialize particle
c            co-ordinates:
c DISTR2 initializes x, y and vx, vy co-ordinates for 2d code, with
c        uniform distribution in space.
c DISTR2H initializes x, y and vx, vy, vz co-ordinates for magnetized
c         2-1/2d codes, with uniform distribution in space.
c LDISTR2 initializes x, y  co-ordinates for 2d code, with bi-linear
c         distribution in space.
c FDISTR2 initializes x, y co-ordinates for 2d code, with general
c         distribution in space.
c RDISTR2 randomizes x, y co-ordinates for 2d code, with general
c         distribution in space previously calculated by FDISTR2.
c VDISTR2 initializes vx, vy co-ordinates for 2d code, with maxwellian
c         velocity distribution with drift.
c VDISTR2H initializes vx, vy, vz co-ordinates for 2-1/2d code, with
c          maxwellian velocity distribution with drift.
c GBDISTR2L calculates guiding centers for magnetized 2-1/2d codes.
c GBZDISTR2L calculates guiding centers for magnetized 2d codes.
c GRBDISTR2L calculates guiding centers for relativistic, magnetized
c            2-1/2d codes.
c GRBZDISTR2L calculates guiding centers for relativistic, magnetized
c             2d codes.
c REVDISTR2 recalculates 2d particle velocities with maxwellian velocity
c           with drift for fraction f of particles, (0.<f<1.)
c REVDISTR2H recalculates 2-1/2d particle velocities with maxwellian
c            velocity with drift for fraction f of particles, (0.<f<1.)
c ranorm = generates gaussian random numbers.
c randum = generates uniform random numbers.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: january 8, 2010
c-----------------------------------------------------------------------
      subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,ipb
     1c)
c for 2d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      double precision ranorm
      double precision dsum1, dsum2
      dimension part(idimp,nop)
      npxy = npx*npy
c set boundary values
      edgelx = 0.
      edgely = 0.
      at1 = float(nx)/float(npx)
      at2 = float(ny)/float(npy)
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny)/float(npy)
      endif
c uniform density profile
      do 20 k = 1, npy
      k1 = npx*(k - 1)
      at3 = edgely + at2*(float(k) - .5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(float(j) - .5)
      part(2,j+k1) = at3
   10 continue
   20 continue
c maxwellian velocity distribution
      do 30 j = 1, npxy
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      at1 = 1./float(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
   50 continue
      return
      end
      subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,nop,
     1nx,ny,ipbc)
c for 2-1/2d code, this subroutine calculates initial particle
c co-ordinates and velocities with uniform density and maxwellian
c velocity with drift
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ranorm = gaussian random number with zero mean and unit variance
      double precision ranorm
      double precision dsum1, dsum2, dsum3
      dimension part(idimp,nop)
      npxy = npx*npy
c set boundary values
      edgelx = 0.
      edgely = 0.
      at1 = float(nx)/float(npx)
      at2 = float(ny)/float(npy)
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny)/float(npy)
      endif
c uniform density profile
      do 20 k = 1, npy
      k1 = npx*(k - 1)
      at3 = edgely + at2*(float(k) - .5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(float(j) - .5)
      part(2,j+k1) = at3
   10 continue
   20 continue
c maxwellian velocity distribution
      do 30 j = 1, npxy
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
      part(5,j) = vtz*ranorm()
   30 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
      dsum3 = dsum3 + part(5,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./float(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
      part(5,j) = part(5,j) - sum3
   50 continue
      return
      end
      subroutine LDISTR2(part,anlx,anly,npx,npy,idimp,nop,nx,ny,ipbc)
c for 2d code, this subroutine calculates initial particle
c co-ordinates with the following bi-linear density profile:
c n(x,y) = n(x)*n(y), where n(x) = n0x*(1. + anlx*(x/nx - .5)) and 
c n(y) = n0y*(1. + anly*(y/ny - .5)) and where
c n0x = npx/(nx - 2*edgelx) and n0y = npy/(ny - 2*edgely)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c anlx/anly = initial linear density weight in x/y direction
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4 or 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real anlx, anly, part
      dimension part(idimp,nop)
c local data
      integer j, k, k1
      real edgelx, edgely, at1, at2, at3, bt1, bt2, antx, anty
c set boundary values
      edgelx = 0.
      edgely = 0.
      at1 = float(nx)/float(npx)
      at2 = float(ny)/float(npy)
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny)/float(npy)
      endif
      if (anly.ne.0.) then
         anty = anly/float(ny)
         at2 = 2.*anty*at2
         bt2 = 1. - .5*anty*(float(ny) - 2.*edgely)
      endif
      if (anlx.ne.0.) then
         antx = anlx/float(nx)
         at1 = 2.*antx*at1
         bt1 = 1. - .5*antx*(float(nx) - 2.*edgelx)
      endif
c linear or uniform density profile
      do 30 k = 1, npy
      k1 = npx*(k - 1)
c linear density in y
      if (anly.ne.0.) then
         at3 = edgely + (sqrt(bt2*bt2 + at2*(float(k) - .5)) - bt2)/anty
c uniform density in y
      else
         at3 = edgely + at2*(float(k) - .5)
      endif
c linear density in x
      if (anlx.ne.0.) then
         do 10 j = 1, npx
         part(1,j+k1) = edgelx + (sqrt(bt1*bt1 + at1*(float(j) - .5)) - 
     1bt1)/antx
         part(2,j+k1) = at3
   10    continue
c uniform density in x
      else
         do 20 j = 1, npx
         part(1,j+k1) = edgelx + at1*(float(j) - .5)
         part(2,j+k1) = at3
   20    continue
      endif
   30 continue
      return
      end
      subroutine FDISTR2(part,fnx,argx1,argx2,argx3,fny,argy1,argy2,argy
     13,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
c for 2d code, this subroutine calculates initial particle co-ordinates
c with general density profile n(x,y) = n(x)*n(y), 
c where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
c and integral of the density is given by fnx(x,argx1,argx2,argx3,1)
c and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
c and integral of the density is given by fny(y,argy1,argy2,argy3,1)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c fnx/fny = density and density integral function in x/y direction
c argx1,argx2,argx3 = arguments to fnx
c argy1,argy2,argy3 = arguments to fny
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4 or 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc, ierr
      real argx1, argx2, argx3, argy1, argy2, argy3, part
      dimension part(idimp,nop)
      real fnx, fny
      external fnx, fny
c local data
      integer imax, i, j, k, k1
      real edgelx, edgely, anx, any, bnx, bny, xt, yt, xt0, yt0, x0, y0
      real xn, yn, eps, big, f, fp
      ierr = 0
c eps = convergence criterion
      imax = 20
      eps = 0.0001
      big = 0.5
c set boundary values
      edgelx = 0.
      edgely = 0.
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
      endif
c find normalization for function
      anx = float(nx) - edgelx
      any = float(ny) - edgely
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      y0 = fny(edgely,argy1,argy2,argy3,1)
      bnx = float(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      bny = float(npy)/(fny(any,argy1,argy2,argy3,1) - y0)
      x0 = bnx*x0 - .5
      y0 = bny*y0 - .5
c density profile in x
      xt0 = edgelx
      xt = xt0 + 0.5/(bnx*fnx(xt0,argx1,argx2,argx3,0))
      do 20 j = 1, npx
      xn = float(j) + x0
c guess next value for xt
      if (j.gt.1) xt = xt + 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
      xt = max(edgelx,min(xt,anx))
      i = 0
   10 f = bnx*fnx(xt,argx1,argx2,argx3,1) - xn
c find improved value for xt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bnx*fnx(xt,argx1,argx2,argx3,0)
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(xt - xt0)
            xt = xt0 + fp
         else
            fp = xt - xt0
c           xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 10
         write (2,*) 'newton iteration max exceeded, xt = ', xt
         ierr = ierr + 1
      endif
      part(1,j) = xt
      xt0 = xt
   20 continue
c quit if error
      if (ierr.ne.0) return
c density profile in y
      yt0 = edgely
      yt = yt0 + 0.5/(bny*fny(yt0,argy1,argy2,argy3,0))
      do 50 k = 1, npy
      k1 = npx*(k - 1)
      yn = float(k) + y0
c guess next value for yt
      if (k.gt.1) yt = yt + 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
      yt = max(edgely,min(yt,any))
      i = 0
   30 f = bny*fny(yt,argy1,argy2,argy3,1) - yn
c find improved value for yt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bny*fny(yt,argy1,argy2,argy3,0)
            yt0 = yt
            yt = yt - f/fp
            yt = max(edgely,min(yt,any))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(yt - yt0)
            yt = yt0 + fp
         else
            fp = yt - yt0
c           yt0 = yt
            yt = yt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 30
         write (2,*) 'newton iteration max exceeded, yt = ', yt
         ierr = ierr + 1
      endif
c store co-ordinates
      do 40 j = 1, npx
      part(1,j+k1) = part(1,j)
      part(2,j+k1) = yt
   40 continue
      yt0 = yt
   50 continue
      return
      end
      subroutine FDISTR2B(part,fnx,argx1,argx2,argx3,fny,argy1,argy2,arg
     1y3,npx,npy,idimp,nop,nx,ny,ipbc,ierr,plasmabstart,plasmabend,plasm   
     1abystart,plasmabyend)
c for 2d code, this subroutine calculates initial particle co-ordinates
c with general density profile n(x,y) = n(x)*n(y), 
c where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
c and integral of the density is given by fnx(x,argx1,argx2,argx3,1)
c and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
c and integral of the density is given by fny(y,argy1,argy2,argy3,1)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c fnx/fny = density and density integral function in x/y direction
c argx1,argx2,argx3 = arguments to fnx
c argy1,argy2,argy3 = arguments to fny
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4 or 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc, ierr
      integer plasmabstart, plasmabend
      integer plasmabystart, plasmabyend
      real argx1, argx2, argx3, argy1, argy2, argy3, part
      dimension part(idimp,nop)
      real fnx, fny
      external fnx, fny
c local data
      integer imax, i, j, k, k1
      real edgelx, edgely, anx, any, bnx, bny, xt, yt, xt0, yt0, x0, y0
      real xn, yn, eps, big, f, fp
      ierr = 0
c eps = convergence criterion
      imax = 20
      eps = 0.0001
      big = 0.5
c set boundary values
      edgelx = float(plasmabstart)
      edgely = float(plasmabystart)
      
c find normalization for function
      anx = float(nx) - float(plasmabend)
      any = float(ny) - float(plasmabyend)
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      y0 = fny(edgely,argy1,argy2,argy3,1)
      bnx = float(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      bny = float(npy)/(fny(any,argy1,argy2,argy3,1) - y0)
      x0 = bnx*x0 - .5
      y0 = bny*y0 - .5
c density profile in x
      xt0 = edgelx
      xt = xt0 + 0.5/(bnx*fnx(xt0,argx1,argx2,argx3,0))
      do 20 j = 1, npx
      xn = float(j) + x0
c guess next value for xt
      if (j.gt.1) xt = xt + 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
      xt = max(edgelx,min(xt,anx))
      i = 0
   10 f = bnx*fnx(xt,argx1,argx2,argx3,1) - xn
c find improved value for xt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bnx*fnx(xt,argx1,argx2,argx3,0)
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(xt - xt0)
            xt = xt0 + fp
         else
            fp = xt - xt0
c           xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 10
         write (2,*) 'newton iteration max exceeded, xt = ', xt
         ierr = ierr + 1
      endif
      part(1,j) = xt
      xt0 = xt
   20 continue
c quit if error
      if (ierr.ne.0) return
c density profile in y
      yt0 = edgely
      yt = yt0 + 0.5/(bny*fny(yt0,argy1,argy2,argy3,0))
      do 50 k = 1, npy
      k1 = npx*(k - 1)
      yn = float(k) + y0
c guess next value for yt
      if (k.gt.1) yt = yt + 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
      yt = max(edgely,min(yt,any))
      i = 0
   30 f = bny*fny(yt,argy1,argy2,argy3,1) - yn
c find improved value for yt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bny*fny(yt,argy1,argy2,argy3,0)
            yt0 = yt
            yt = yt - f/fp
            yt = max(edgely,min(yt,any))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(yt - yt0)
            yt = yt0 + fp
         else
            fp = yt - yt0
c           yt0 = yt
            yt = yt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 30
         write (2,*) 'newton iteration max exceeded, yt = ', yt
         ierr = ierr + 1
      endif
c store co-ordinates
      do 40 j = 1, npx
      part(1,j+k1) = part(1,j)
      part(2,j+k1) = yt
   40 continue
      yt0 = yt
   50 continue
      return
      end
      subroutine RDISTR2(part,fnx,argx1,argx2,argx3,fny,argy1,argy2,argy
     13,npx,npy,idimp,nop,nx,ny,ipbc)
c for 2d code, this subroutine randomizes initial particle co-ordinates
c which have been previously calculated using a general density profile
c n(x,y) = n(x)*n(y), 
c where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
c and integral of the density is given by fnx(x,argx1,argx2,argx3,1)
c and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
c and integral of the density is given by fny(y,argy1,argy2,argy3,1)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c fnx/fny = density and density integral function in x/y direction
c argx1,argx2,argx3 = arguments to fnx
c argy1,argy2,argy3 = arguments to fny
c npx/npy = initial number of particles distributed in x/y direction
c idimp = size of phase space = 4 or 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real argx1, argx2, argx3, argy1, argy2, argy3, part
      dimension part(idimp,nop)
      real fnx, fny
      external fnx, fny
c local data
      integer j, k, k1
      real edgelx, edgely, anx, any, bnx, bny, xt, yt, xt0, yt0, x0, y0
      real xn, yn
      double precision randum
c set boundary values
      edgelx = 0.
      edgely = 0.
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
      endif
c find normalization for function
      anx = float(nx) - edgelx
      any = float(ny) - edgely
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      y0 = fny(edgely,argy1,argy2,argy3,1)
      bnx = float(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      bny = float(npy)/(fny(any,argy1,argy2,argy3,1) - y0)
c randomize co-ordinates
      do 70 k = 1, npy
      k1 = npx*(k - 1)
      yt = part(2,k1+1)
      yt0 = 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
      do 60 j = 1, npx
      xt = part(1,j+k1)
      xt0 = 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
      xn = xt + xt0*(randum() - 0.5)
      xn = max(edgelx,min(xn,anx))
      part(1,j+k1) = xn
      yn = yt + yt0*(randum() - 0.5)
      yn = max(edgely,min(yn,any))
      part(2,j+k1) = yn
   60 continue
   70 continue
      return
      end
      function FLDISTR1(x,anlx,anxi,shift,intg)
c this function calculates either a density function or its integral
c for a linear density profile.  Used in initializing particle
c coordinates.  The three parameters are redundant, and one can set one
c of them arbitrarily.  A convenient choice is to set  anxi = 1/Lx,
c anlx = NH - NL, shift = (1 - NL)/(NH - NL), where NL is the density
c at the left, and NH at the right compared to the average density
c if intg = 0, n(x) = 1. + anlx*(x*anxi - shift)
c if intg = 1, n(x) = x + .5*anlx*x*(x*anxi - 2.*shift)
      implicit none
      integer intg
      real x, anlx, anxi, shift
c local data
      real FLDISTR1, f
      if (intg.eq.0) then
         f = 1.0 + anlx*(x*anxi - shift)
      else if (intg.eq.1) then
         if (anxi.eq.0.) then
            f = x
         else
            f = x + .5*anlx*x*(x*anxi - 2.*shift)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FLDISTR1 Error: f = ', f
      FLDISTR1 = f
      return
      end
      function FLDISTR1B(x,anlx,anxi,shift,intg)
c this function calculates either a density function or its integral
c for a linear density profile.  Used in initializing particle
c coordinates.  The three parameters are redundant, and one can set one
c of them arbitrarily.  A convenient choice is to set  anxi = 1/Lx,
c anlx = NH - NL, shift = (1 - NL)/(NH - NL), where NL is the density
c at the left, and NH at the right compared to the average density
c if intg = 0, n(x) = 1. + anlx*(x*anxi - shift)
c if intg = 1, n(x) = x + .5*anlx*x*(x*anxi - 2.*shift)
      implicit none
      integer intg
      real x, anlx, anxi, shift
c local data
      real FLDISTR1B, f
      if (intg.eq.0) then
         f = 1.0 + anlx*(x*anxi - shift)
      else if (intg.eq.1) then
         if (anxi.eq.0.) then
            f = x
         else
            f = x + .5*anlx*x*(x*anxi - 2.*shift)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FLDISTR1 Error: f = ', f
      FLDISTR1B = f
      return
      end
      function FSDISTR1(x,ans,dkx,phase,intg)
c this function calculates either a density function or its integral
c for a sinusoidal density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ans*sin(dkx*x - phase)
c if intg = 1, n(x) = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
      implicit none
      integer intg
      real x, ans, dkx, phase
c local data
      real FSDISTR1, f
      if (intg.eq.0) then
         f = 1.0 + ans*sin(dkx*x - phase)
      else if (intg.eq.1) then
         if (dkx.eq.0.) then
            f = x - ans*sin(phase)*x
         else
            f = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FSDISTR1 Error: f = ', f
      FSDISTR1 = f
      return
      end
      function FGDISTR1(x,ang,wi,x0,intg)
c this function calculates either a density function or its integral
c for a gaussian density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ang*exp(-((x-x0)*wi)**2/2.)
c if intg = 1, n(x) = x + (ang*sqrt(pi/2)/wi)*
c                         (erf((x-x0)*wi/sqrt(2)) + erf(x0*wi/sqrt(2)))
      implicit none
      integer intg
      real x, ang, x0, wi
c local data
      real FGDISTR1, f, sqrt2i, sqtpih, aw, t, erfn
      external erfn
      data sqrt2i, sqtpih /0.7071067811865476,1.253314137397325/
      save sqrt2i, sqtpih
      aw = wi*sqrt2i
      t = (x - x0)*aw
      if (intg.eq.0) then
         if (abs(t).lt.8.) then
            f = 1.0 + ang*exp(-t**2)
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + ang)*x
         else
            f = x + (ang*sqtpih/wi)*(erfn(t) + erfn(x0*aw))
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FGDISTR1 Error: f = ', f
      FGDISTR1 = f
      return
      end
      function FHDISTR1(x,anh,wi,x0,intg)
c this function calculates either a density function or its integral
c for a hyperbolic secant squared density profile.  Used in initializing
c particle coordinates.
c if intg = 0, n(x) = 1.0 + anh*sech((x-x0)*wi)**2
c if intg = 1, n(x) = x + (anh/wi)*(tanh((x-x0)*wi) + tanh(x0*wi))
      implicit none
      integer intg
      real x, anh, x0, wi
c local data
      real FHDISTR1, f, g, t, u
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (abs(t).lt.32.) then
            u = exp(-abs(t))
            f = 1.0 + anh*(2.*u/(1.0 + u*u))**2
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + anh)*x
         else
            if (abs(t).lt.32.) then
               u = exp(-abs(t))**2
               f = (1.0 - u)/(1.0 + u)
            else
               f = 1.0
            endif
            if (t.lt.0.) f = -f
            t = x0*wi
            if (abs(t).lt.32.) then
               u = exp(-abs(t))**2
               g = (1.0 - u)/(1.0 + u)
            else
               g = 1.0
            endif
            if (t.lt.0.) g = -g
            f = x + (anh/wi)*(f + g)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FHDISTR1 Error: f = ', f
      FHDISTR1 = f
      return
      end
      function FEDISTR1(x,ane,wi,x0,intg)
c this function calculates either a density function or its integral
c for an exponential density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ane*exp((x-x0)*wi)
c if intg = 1, n(x) = x + (ane/wi)*(exp((x-x0)*wi) - exp(-x0*wi)
      implicit none
      integer intg
      real x, ane, x0, wi
c local data
      real FEDISTR1, f, t
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (t.gt.(-64.)) then
            f = 1.0 + ane*exp(t)
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + ane)*x
         else
            f = x0*wi
            if (f.gt.64) then
               f = 0.
            else
               f = exp(-f)
            endif
            f = x + (ane/wi)*(exp(t) - f)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FEDISTR1 Error: f = ', f
      FEDISTR1 = f
      return
      end
      subroutine VDISTR2(part,vtx,vty,vdx,vdy,idimp,nop)
c for 2d code, this subroutine calculates initial particle
c velocities with maxwellian velocity with drift
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c idimp = size of phase space = 4
c nop = number of particles
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer idimp, nop
      real vtx, vty, vdx, vdy, part
      dimension part(idimp,nop)
c local data
      integer j
      real sum1, sum2, at1
      double precision ranorm
      double precision dsum1, dsum2
      if (nop.eq.0) return
c maxwellian velocity distribution
      do 10 j = 1, nop
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
   10 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      do 20 j = 1, nop
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
   20 continue
      sum1 = dsum1
      sum2 = dsum2
      at1 = 1./float(nop)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      do 30 j = 1, nop
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
   30 continue
      return
      end
      subroutine VDISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,idimp,nop)
c for 2-1/2d code, this subroutine calculates initial particle
c velocities with maxwellian velocity with drift
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c idimp = size of phase space = 5
c nop = number of particles
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer idimp, nop
      real vtx, vty, vtz, vdx, vdy, vdz, part
      dimension part(idimp,nop)
c local data
      integer j
      real sum1, sum2, sum3, at1
      double precision ranorm
      double precision dsum1, dsum2, dsum3
      if (nop.eq.0) return
c maxwellian velocity distribution
      do 10 j = 1, nop
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
      part(5,j) = vtz*ranorm()
   10 continue
c add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 20 j = 1, nop
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
      dsum3 = dsum3 + part(5,j)
   20 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./float(nop)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 30 j = 1, nop
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
      part(5,j) = part(5,j) - sum3
   30 continue
      return
      end
      subroutine BDISTR2L(part,bx,by,bz,qbm,idimp,nop,nx,ny,nxv)
c for 2-1/2d code, this subroutine reinterprets current particle
c positions as positions of guiding centers, and calculates the actual
c particle positions, with periodic boundary conditions
c in converting from guiding center to actual co-ordinates,
c the following equations are used:
c       x(t) = xg(t) - (vy(t)*omz - vz(t)*omy)/om**2
c       y(t) = yg(t) - (vz(t)*omx - vx(t)*omz)/om**2
c where omx = (q/m)*bx(xg(t),yg(t)),
c       omy = (q/m)*by(xg(t),yg(t)),
c and   omz = (q/m)*bz(xg(t),yg(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c bx(x,y) = (1-dy)*((1-dx)*bx(n,m)+dx*bx(n+1,m)) + dy*((1-dx)*bx(n,m+1)
c    + dx*bx(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c bx(j,k) = x component of magnetic field at grid (j,k)
c by(j,k) = y component of magnetic field at grid (j,k)
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      dimension part(idimp,nop)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
c calculate actual position from guiding center
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qbm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      amx = qbm - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      mm = mm + 1
      amy = 1. - dyp
      mp = mm + 1
      if (mp.gt.ny) mp = mp - ny
c find magnetic field
      omx = amy*(amx*bx(nn,mm) + dxp*bx(np,mm)) + dyp*(amx*bx(nn,mp) + d
     1xp*bx(np,mp))
      omy = amy*(amx*by(nn,mm) + dxp*by(np,mm)) + dyp*(amx*by(nn,mp) + d
     1xp*by(np,mp))
      omz = amy*(amx*bz(nn,mm) + dxp*bz(np,mm)) + dyp*(amx*bz(nn,mp) + d
     1xp*bz(np,mp))
      at3 = sqrt(omx*omx + omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3
      omxt = omx*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j) - (part(4,j)*omzt - part(5,j)*omyt)
      dy = part(2,j) - (part(5,j)*omxt - part(3,j)*omzt)
c periodic boundary conditions
      n = abs(dx)/anx
      if (dx.lt.zero) dx = dx + float(n + 1)*anx
      if (dx.ge.anx) dx = dx - float(n)*anx
      part(1,j) = dx
      m = abs(dy)/any
      if (dy.lt.zero) dy = dy + float(m + 1)*any
      if (dy.ge.any) dy = dy - float(m)*any
      part(2,j) = dy
   10 continue
      return
      end
      subroutine GBDISTR2L(part,bxy,qbm,idimp,nop,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine reinterprets current particle
c positions as positions of guiding centers, and calculates the actual
c particle positions
c in converting from guiding center to actual co-ordinates,
c the following equations are used:
c       x(t) = xg(t) - (vy(t)*omz - vz(t)*omy)/om**2
c       y(t) = yg(t) - (vz(t)*omx - vx(t)*omz)/om**2
c where omx = (q/m)*bxy(1,xg(t),yg(t)),
c       omy = (q/m)*bxy(2,xg(t),yg(t)),
c and   omz = (q/m)*bxy(3,xg(t),yg(t)),
c and the magnetic field components bxyz(i,x(t),y(t)) are approximated
c by interpolation from the nearest grid points:
c bxy(i,x,y) = (1-dy)*((1-dx)*bxy(i,n,m)+dx*bxy(i,n+1,m)) + 
c               dy*((1-dx)*bxy(i,n,m+1) + dx*bxy(i,n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c bxyz(i,j,k) = i component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
      dimension bxy(3,nxv,nyv)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgerx = float(nx)
         edgery = float(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         edgerx = float(nx-1)
         edgery = float(ny)
      endif
c calculate actual position from guiding center
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qbm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qbm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find magnetic field
      omx = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp)) + amy*(dxp*bxy(1,n
     1p,mm) + amx*bxy(1,nn,mm))
      omy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp)) + amy*(dxp*bxy(2,n
     1p,mm) + amx*bxy(2,nn,mm))
      omz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp)) + amy*(dxp*bxy(3,n
     1p,mm) + amx*bxy(3,nn,mm))
      at3 = sqrt(omx*omx + omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3
      omxt = omx*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j) - (part(4,j)*omzt - part(5,j)*omyt)
      dy = part(2,j) - (part(5,j)*omxt - part(3,j)*omzt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + float(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - float(n)*edgerx
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery)) then
            if ((dy.ge.edgely).and.(dy.lt.edgery)) then
c if x co-ordinate only is out of bounds, try switching vy
               dx = part(1,j) + (part(4,j)*omzt + part(5,j)*omyt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(4,j) = -part(4,j)
               else
c otherwise, try switching both vy and vz
                  dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
                  dy = part(2,j) + (part(5,j)*omxt + part(3,j)*omzt)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(4,j) = -part(4,j)
                     part(5,j) = -part(5,j)
                  endif
               endif
            else if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
c if y co-ordinate only is out of bounds, try switching vx
               dy = part(2,j) - (part(5,j)*omxt + part(3,j)*omzt)
               if ((dy.ge.edgely).and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
               else
c otherwise, try switching both vx and vz
                  dx = part(1,j) - (part(4,j)*omzt + part(5,j)*omyt)
                  dy = part(2,j) + (part(5,j)*omxt - part(3,j)*omzt)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(3,j) = -part(3,j)
                     part(5,j) = -part(5,j)
                  endif
               endif
            endif
c if both co-ordinates are out of bounds, try switching vx, vy, vz
            if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(d
     1y.ge.edgery)) then
               dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
               dy = part(2,j) + (part(5,j)*omxt - part(3,j)*omzt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).
     1and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
                  part(4,j) = -part(4,j)
                  part(5,j) = -part(5,j)
               else
c give up if larmor radius is too large
                  dx = part(1,j)
                  dy = part(2,j)
               endif
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
c rotate particle position by reversing velocity in y and z
            dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
            dy = part(2,j) + (part(5,j)*omxt + part(3,j)*omzt)
c give up if larmor radius is too large
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j)
               dy = part(2,j)
            else
               part(4,j) = -part(4,j)
               part(5,j) = -part(5,j)
            endif
         endif
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
      subroutine GBZDISTR2L(part,bz,qbm,idimp,nop,nx,ny,nxv,nyv,ipbc)
c for 2d code, this subroutine reinterprets current particle
c positions as positions of guiding centers, and calculates the actual
c particle positions
c in converting from guiding center to actual co-ordinates,
c the following equations are used:
c       x(t) = xg(t) - (vy(t)*omz)/om**2
c       y(t) = yg(t) + (vx(t)*omz)/om**2
c where omz = (q/m)*bz(xg(t),yg(t)),
c and the magnetic field component bz(x(t),y(t)) is approximated
c by interpolation from the nearest grid points:
c bz(x,y) = (1-dy)*((1-dx)*bz(n,m)+dx*bz(n+1,m)) + 
c               dy*((1-dx)*bz(n,m+1) + dx*bz(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
      dimension bz(nxv,nyv)
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgerx = float(nx)
         edgery = float(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         edgerx = float(nx-1)
         edgery = float(ny)
      endif
c calculate actual position from guiding center
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qbm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qbm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find magnetic field
      omz = dyp*(dxp*bz(np,mp) + amx*bz(nn,mp)) + amy*(dxp*bz(np,mm) + a
     1mx*bz(nn,mm))
      at3 = abs(omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3
      omzt = omz*at3
c correct position
      dx = part(1,j) - part(4,j)*omzt
      dy = part(2,j) + part(3,j)*omzt
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + float(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - float(n)*edgerx
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery)) then
            if ((dy.ge.edgely).and.(dy.lt.edgery)) then
c if x co-ordinate only is out of bounds, try switching vy
               dx = part(1,j) + part(4,j)*omzt
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(4,j) = -part(4,j)
               endif
            else if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
c if y co-ordinate only is out of bounds, try switching vx
               dy = part(2,j) - part(3,j)*omzt
               if ((dy.ge.edgely).and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
               endif
            endif
c if both co-ordinates are out of bounds, try switching vx, vy
            if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(d
     1y.ge.edgery)) then
               dx = part(1,j) + part(4,j)*omzt
               dy = part(2,j) - part(3,j)*omzt
               if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).
     1and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
                  part(4,j) = -part(4,j)
               else
c give up if larmor radius is too large
                  dx = part(1,j)
                  dy = part(2,j)
               endif
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
c rotate particle position by reversing velocity in y
            dx = part(1,j) + part(4,j)*omzt
            dy = part(2,j) + part(3,j)*omzt
c give up if larmor radius is too large
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j)
               dy = part(2,j)
            else
               part(4,j) = -part(4,j)
            endif
         endif
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
      subroutine GRBDISTR2L(part,bxy,qbm,ci,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2-1/2d code, this subroutine reinterprets current particle
c positions as positions of guiding centers, and calculates the actual
c particle positions for relativistic particles
c in converting from guiding center to actual co-ordinates,
c the following equations are used:
c       x(t) = xg(t) - gami*(py(t)*omz - pz(t)*omy)/om**2
c       y(t) = yg(t) - gami*(pz(t)*omx - px(t)*omz)/om**2
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c and omx = (q/m)*bxy(1,xg(t),yg(t)),
c     omy = (q/m)*bxy(2,xg(t),yg(t)),
c     omz = (q/m)*bxy(2,xg(t),yg(t)),
c and the magnetic field components bxyz(i,x(t),y(t)) are approximated
c by interpolation from the nearest grid points:
c bxy(i,x,y) = (1-dy)*((1-dx)*bxy(i,n,m)+dx*bxy(i,n+1,m)) + 
c               dy*((1-dx)*bxy(i,n,m+1) + dx*bxy(i,n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c bxyz(i,j,k) = i component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c ci = reciprical of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
      dimension bxy(3,nxv,nyv)
      ci2 = ci*ci
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgerx = float(nx)
         edgery = float(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         edgerx = float(nx-1)
         edgery = float(ny)
      endif
c calculate actual position from guiding center
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qbm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qbm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      omx = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp)) + amy*(dxp*bxy(1,n
     1p,mm) + amx*bxy(1,nn,mm))
      omy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp)) + amy*(dxp*bxy(2,n
     1p,mm) + amx*bxy(2,nn,mm))
      omz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp)) + amy*(dxp*bxy(3,n
     1p,mm) + amx*bxy(3,nn,mm))
      at3 = sqrt(omx*omx + omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3*gami
      omxt = omx*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j) - (part(4,j)*omzt - part(5,j)*omyt)
      dy = part(2,j) - (part(5,j)*omxt - part(3,j)*omzt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + float(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - float(n)*edgerx
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery)) then
            if ((dy.ge.edgely).and.(dy.lt.edgery)) then
c if x co-ordinate only is out of bounds, try switching vy
               dx = part(1,j) + (part(4,j)*omzt + part(5,j)*omyt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(4,j) = -part(4,j)
               else
c otherwise, try switching both vy and vz
                  dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
                  dy = part(2,j) + (part(5,j)*omxt + part(3,j)*omzt)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(4,j) = -part(4,j)
                     part(5,j) = -part(5,j)
                  endif
               endif
            else if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
c if y co-ordinate only is out of bounds, try switching vx
               dy = part(2,j) - (part(5,j)*omxt + part(3,j)*omzt)
               if ((dy.ge.edgely).and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
               else
c otherwise, try switching both vx and vz
                  dx = part(1,j) - (part(4,j)*omzt + part(5,j)*omyt)
                  dy = part(2,j) + (part(5,j)*omxt - part(3,j)*omzt)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(3,j) = -part(3,j)
                     part(5,j) = -part(5,j)
                  endif
               endif
            endif
c if both co-ordinates are out of bounds, try switching vx, vy, vz
            if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(d
     1y.ge.edgery)) then
               dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
               dy = part(2,j) + (part(5,j)*omxt - part(3,j)*omzt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).
     1and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
                  part(4,j) = -part(4,j)
                  part(5,j) = -part(5,j)
               else
c give up if larmor radius is too large
                  dx = part(1,j)
                  dy = part(2,j)
               endif
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
c rotate particle position by reversing velocity in y and z
            dx = part(1,j) + (part(4,j)*omzt - part(5,j)*omyt)
            dy = part(2,j) + (part(5,j)*omxt + part(3,j)*omzt)
c give up if larmor radius is too large
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j)
               dy = part(2,j)
            else
               part(4,j) = -part(4,j)
               part(5,j) = -part(5,j)
            endif
         endif
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
      subroutine GRBZDISTR2L(part,bz,qbm,ci,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine reinterprets current particle
c positions as positions of guiding centers, and calculates the actual
c particle positions for relativistic particles
c in converting from guiding center to actual co-ordinates,
c the following equations are used:
c       x(t) = xg(t) - gami*(py(t)*omz)/om**2
c       y(t) = yg(t) + gami*(px(t)*omz)/om**2
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c and omz = (q/m)*bz(xg(t),yg(t)),
c and the magnetic field component bz(x(t),y(t)) is approximated
c by interpolation from the nearest grid points:
c bz(x,y) = (1-dy)*((1-dx)*bz(n,m)+dx*bz(n+1,m)) + 
c               dy*((1-dx)*bz(n,m+1) + dx*bz(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c ci = reciprical of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
      dimension bz(nxv,nyv)
      ci2 = ci*ci
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgerx = float(nx)
         edgery = float(ny)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 0.
         edgerx = float(nx-1)
         edgery = float(ny)
      endif
c calculate actual position from guiding center
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qbm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qbm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      omz = dyp*(dxp*bz(np,mp) + amx*bz(nn,mp)) + amy*(dxp*bz(np,mm) + a
     1mx*bz(nn,mm))
      at3 = abs(omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3*gami
      omzt = omz*at3
c correct position
      dx = part(1,j) - part(4,j)*omzt
      dy = part(2,j) + part(3,j)*omzt
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + float(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - float(n)*edgerx
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery)) then
            if ((dy.ge.edgely).and.(dy.lt.edgery)) then
c if x co-ordinate only is out of bounds, try switching vy
               dx = part(1,j) + part(4,j)*omzt
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(4,j) = -part(4,j)
               endif
            else if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
c if y co-ordinate only is out of bounds, try switching vx
               dy = part(2,j) - part(3,j)*omzt
               if ((dy.ge.edgely).and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
               endif
            endif
c if both co-ordinates are out of bounds, try switching vx, vy
            if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(d
     1y.ge.edgery)) then
               dx = part(1,j) + part(4,j)*omzt
               dy = part(2,j) - part(3,j)*omzt
               if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).
     1and.(dy.lt.edgery)) then
                  part(3,j) = -part(3,j)
                  part(4,j) = -part(4,j)
               else
c give up if larmor radius is too large
                  dx = part(1,j)
                  dy = part(2,j)
               endif
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
c rotate particle position by reversing velocity in y
            dx = part(1,j) + part(4,j)*omzt
            dy = part(2,j) + part(3,j)*omzt
c give up if larmor radius is too large
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j)
               dy = part(2,j)
            else
               part(4,j) = -part(4,j)
            endif
         endif
         m = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(m + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(m)*edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
      subroutine REVDISTR2(part,vtx,vty,vdx,vdy,f,nstart,npxy,idimp,nop)
c for 2d code, this subroutine recalculates particle velocities with
c maxwellian velocity with drift for fraction f of particles, (0.<f<1.)
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c vtx/vty = thermal velocity of electrons in x/y direction
c vdx/vdy = drift velocity of beam electrons in x/y direction
c f = fraction of particles to be redistributed
c nstart = initial redistributed particle
c npxy = number of particles redistributed
c idimp = size of phase space = 4
c nop = number of particles
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer nstart, npxy, idimp, nop
      real vtx, vty, vdx, vdy, f
      real part
      dimension part(idimp,nop)
c local data
      integer j, ns, np
      real fc, r
      double precision randum, ranorm
      fc = 1.0 - f
      ns = nstart - 1
      np = min0(npxy+ns,nop)
c maxwellian velocity distribution with drift
      do 10 j = 1, np
      r = randum()
      if (r.ge.fc) then
         part(3,j+ns) = vtx*ranorm() + vdx
         part(4,j+ns) = vty*ranorm() + vdy
      endif
   10 continue
      return
      end
      subroutine REVDISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,f,nstart,npxy,i
     1dimp,nop)
c for 2-1/2d code, this subroutine recalculates particle velocities with
c maxwellian velocity with drift for fraction f of particles, (0.<f<1.)
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c f = fraction of particles to be redistributed
c nstart = initial redistributed particle
c npxy = number of particles redistributed
c idimp = size of phase space = 5
c nop = number of particles
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer nstart, npxy, idimp, nop
      real vtx, vty, vtz, vdx, vdy, vdz, f
      real part
      dimension part(idimp,nop)
c local data
      integer j, ns, np
      real fc, r
      double precision randum, ranorm
      fc = 1.0 - f
      ns = nstart - 1
      np = min0(npxy+ns,nop)
c maxwellian velocity distribution with drift
      do 10 j = 1, np
      r = randum()
      if (r.ge.fc) then
         part(3,j+ns) = vtx*ranorm() + vdx
         part(4,j+ns) = vty*ranorm() + vdy
         part(5,j+ns) = vtz*ranorm() + vdz
      endif
   10 continue
      return
      end
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      integer r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      integer r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
      function erfn(x)
c this function calculates the real error function, according to the
c formulae given in Abramowitz and Stegun, Handbook of Mathematical
c Functions, p. 299.  Error is < 1.5 x 10-7.
      implicit none
      real x
c local data
      real erfn, p, a1, a2, a3, a4, a5, t, f
      data p, a1, a2 /0.3275911,0.254829592,-0.284496736/
      data a3, a4, a5 /1.421413741,-1.453152027,1.061405429/
      save p, a1, a2, a3, a4, a5
      f = abs(x)
      t = 1.0/(1.0 + p*f)
      if (f.le.8.) then
         erfn = 1.0 - t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))*exp(-x*x)
      else
         erfn = 1.0
      endif
      if (x.lt.0.) erfn = -erfn
      return
      end
      function e1ln(x)
c this function calculates the sum of the exponential integral and the
c natural logarithm, according to the formulae given in Abramowitz and
c Stegun, Handbook of Mathematical Functions, p. 231.
c Error is < 2.0 x 10-7.
      implicit none
      real x
c local data
      real e1ln, a0, a1, a2, a3, a4, a5, b1, b2, b3, b4, c1, c2, c3, c4
      data a0, a1, a2 /-0.57721566,0.99999193,-0.24991055/
      data a3, a4, a5 /0.05519968,-0.00976004,0.00107857/
      data b1, b2, b3 /8.5733287401,18.0590169730,8.6347608925/
      data c1, c2, c3 /9.5733223454,25.6329561486,21.0996530827/
      data b4, c4 /0.2677737343,3.9584969228/
      save 
      if (x.le.1.0) then
         e1ln = a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5))))
      else if (x.lt.50.0) then
         e1ln = alog(x) + (exp(-x)/x)*((b4 + x*(b3 + x*(b2 + x*(b1 + x))
     1))/(c4 + x*(c3 + x*(c2 + x*(c1 + x)))))
      else
         e1ln = alog(x)
      endif
      return
      end

