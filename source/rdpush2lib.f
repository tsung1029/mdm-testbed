c-----------------------------------------------------------------------
c 2d PIC library for pushing relativistic particles with darwin electric
c and magnetic fields and depositing current and derivative of current
c rdpush2lib.f contains procedures to process relativistic particles
c              with darwin electric and magnetic fields:
c GRMJPOST2 deposits momentum flux for 2-1/2d code, quadratic
c           interpolation, STANDARD optimization, for relativistic
c           particles.
c GSRMJPOST2 deposits momentum flux for 2-1/2d code, quadratic
c            interpolation, LOOKAHEAD optimization, for relativistic
c            particles.
c GRDCJPOST2 deposits momentum flux, acceleration density and current
c            density for 2-1/2d code, quadratic interpolation, STANDARD
c            optimization, for relativistic particles.
c GSRDCJPOST2 deposits momentum flux, acceleration density and current
c             density for 2-1/2d code, quadratic interpolation,
c             LOOKAHEAD optimization, for relativistic particles.
c GRMJPOST22 deposits momentum flux for 2d code, quadratic
c            interpolation, STANDARD optimization, for relativistic
c            particles.
c GSRMJPOST22 deposits momentum flux for 2d code, quadratic
c             interpolation, LOOKAHEAD optimization, for relativistic
c             particles.
c GRDCJPOST22 deposits momentum flux, acceleration density and current
c             density for 2d code, quadratic interpolation, STANDARD
c             optimization, for relativistic particles.
c GSRDCJPOST22 deposits momentum flux, acceleration density and current
c              density for 2d code, quadratic interpolation, LOOKAHEAD
c              optimization, for relativistic particles.
c GRMJPOST2L deposits momentum flux for 2-1/2d code, linear
c            interpolation, STANDARD optimization, for relativistic
c            particles.
c GSRMJPOST2L deposits momentum flux for 2-1/2d code, linear
c             interpolation, LOOKAHEAD optimization, for relativistic
c             particles.
c GRDCJPOST2L deposits momentum flux, acceleration density and current
c             density for 2-1/2d code, linear interpolation, STANDARD
c             optimization, for relativistic particles.
c GSRDCJPOST2L deposits momentum flux, acceleration density and current
c              density for 2-1/2d code, linear interpolation, LOOKAHEAD
c              optimization, for relativistic particles.
c GRMJPOST22L deposits momentum flux for 2d code, linear interpolation,
c             STANDARD optimization, for relativistic particles.
c GSRMJPOST22L deposits momentum flux for 2d code, linear interpolation,
c              LOOKAHEAD optimization, for relativistic particles.
c GRDCJPOST22L deposits momentum flux, acceleration density and current
c              density for 2d code, linear interpolation, STANDARD
c              optimization, for relativistic particles.
c GSRDCJPOST22L deposits momentum flux, acceleration density and current
c               density for 2d code, linear interpolation, LOOKAHEAD
c               optimization, for relativistic particles.
c GRMJPOST2C deposits momentum flux for 2-1/2d code, cubic
c            interpolation, STANDARD optimization, for relativistic
c            particles.
c GRDCJPOST2C deposits momentum flux, acceleration density and current
c             density for 2-1/2d code, cubic interpolation, STANDARD
c             optimization, for relativistic particles.
c GRMJPOST22C deposits momentum flux for 2 code, cubic interpolation,
c             STANDARD optimization, for relativistic particles.
c GRDCJPOST22C deposits momentum flux, acceleration density and current
c              density for 2d code, cubic interpolation, LOOKAHEAD
c              optimization, for relativistic particles.
c GRMJPOST2GL deposits momentum flux, using gridless method, for
c             relativistic particles.
c GRDCJPOST2GL deposits momentum flux, acceleration density and current
c              density, using gridless method, for relativistic
c              particles.
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: september 16, 2011
c-----------------------------------------------------------------------
      subroutine GRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation for relativistic particles
c scalar version using guard cells
c 123 flops/particle, 1 divide, 41 loads, 36 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,j+1,k+1) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nxv = first dimension of flux array, must be >= nx+3
c nyv = second dimension of flux array, must be >= ny+3
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(4,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qmh, ci2, gami2, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, vz, p2, v1, v2, v3, v4
      qmh = .5*qm
      ci2 = ci*ci
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = 1.0/(1.0 + p2*ci2)
c calculate weights
      nl = nn + 1
      amx = qm*(.75 - dxp*dxp)
      mm = mm + 2
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = qmh*(.5 - dxp)**2
      np = nl + 2
      dxp = qmh*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vz*vx)*gami2
      v4 = (vz*vy)*gami2
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      amu(3,nl,mm) = amu(3,nl,mm) + v3*dx
      amu(4,nl,mm) = amu(4,nl,mm) + v4*dx
      dx = dxl*dyl
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dy = amx*dyl
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      amu(3,np,mm) = amu(3,np,mm) + v3*dz
      amu(4,np,mm) = amu(4,np,mm) + v4*dz
      dz = dxp*dyl
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      amu(3,nl,ml) = amu(3,nl,ml) + v3*dx
      amu(4,nl,ml) = amu(4,nl,ml) + v4*dx
      dx = dxl*dyp
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      amu(3,nn,ml) = amu(3,nn,ml) + v3*dy
      amu(4,nn,ml) = amu(4,nn,ml) + v4*dy
      dy = amx*dyp
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      amu(3,np,ml) = amu(3,np,ml) + v3*dz
      amu(4,np,ml) = amu(4,np,ml) + v4*dz
      dz = dxp*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      amu(3,nl,mp) = amu(3,nl,mp) + v3*dx
      amu(4,nl,mp) = amu(4,nl,mp) + v4*dx
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      amu(3,np,mp) = amu(3,np,mp) + v3*dz
      amu(4,np,mp) = amu(4,np,mp) + v4*dz
   10 continue
      return
      end
      subroutine GSRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation for relativistic particles
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 123 flops/particle, 1 divide, 41 loads, 36 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,n) = ith component of momentum flux at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyv = dimension of current array, must be >= nxv*(ny+3)
      implicit none
      integer nop, idimp, nxv, nxyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(4,nxyv)
      integer nnn, mmn, j, nn, mm, ml, mn, mp
      real dxn, dyn, qmh, ci2, gami2, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vxn, vyn, vzn, vx, vy, vz, p2, v1, v2, v3, v4
      real dx1, dy1, dx2, dy2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      qmh = .5*qm
      ci2 = ci*ci
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      vzn = part(5,1)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami2 = 1.0/(1.0 + p2*ci2)
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j) + .5
      mmn = part(2,j) + .5
      dxp = dxn
      dyp = dyn
      dxn = part(1,j) - float(nnn)
      dyn = part(2,j) - float(mmn)
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = mm + nn
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv
c calculate weights
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = (vxn*vxn - vyn*vyn)*gami2
      v2 = (vxn*vyn)*gami2
      v3 = (vzn*vxn)*gami2
      v4 = (vzn*vyn)*gami2
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      vzn = part(5,j)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit momentum flux
      dx1 = amu(1,mn) + v1*dx
      dy1 = amu(2,mn) + v2*dx
      amy = amu(3,mn) + v3*dx
      vx = amu(4,mn) + v4*dx
      dx2 = amu(1,mn+1) + v1*dy
      dy2 = amu(2,mn+1) + v2*dy
      dx3 = amu(3,mn+1) + v3*dy
      vy = amu(4,mn+1) + v4*dy
      dx = amu(1,mn+2) + v1*dz
      dy = amu(2,mn+2) + v2*dz
      vz = amu(3,mn+2) + v3*dz
      dz = amu(4,mn+2) + v4*dz
      amu(1,mn) = dx1
      amu(2,mn) = dy1
      amu(3,mn) = amy
      amu(4,mn) = vx
      amu(1,mn+1) = dx2
      amu(2,mn+1) = dy2
      amu(3,mn+1) = dx3
      amu(4,mn+1) = vy
      amu(1,mn+2) = dx
      amu(2,mn+2) = dy
      amu(3,mn+2) = vz
      amu(4,mn+2) = dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml) + v1*dx
      dy1 = amu(2,ml) + v2*dx
      amy = amu(3,ml) + v3*dx
      vx = amu(4,ml) + v4*dx
      dx2 = amu(1,ml+1) + v1*dy
      dy2 = amu(2,ml+1) + v2*dy
      dyl = amu(3,ml+1) + v3*dy
      vy = amu(4,ml+1) + v4*dy
      dx = amu(1,ml+2) + v1*dz
      dy = amu(2,ml+2) + v2*dz
      vz = amu(3,ml+2) + v3*dz
      dz = amu(4,ml+2) + v4*dz
      amu(1,ml) = dx1
      amu(2,ml) = dy1
      amu(3,ml) = amy
      amu(4,ml) = vx
      amu(1,ml+1) = dx2
      amu(2,ml+1) = dy2
      amu(3,ml+1) = dyl
      amu(4,ml+1) = vy
      amu(1,ml+2) = dx
      amu(2,ml+2) = dy
      amu(3,ml+2) = vz
      amu(4,ml+2) = dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp) + v1*dx
      dy1 = amu(2,mp) + v2*dx
      amy = amu(3,mp) + v3*dx
      vx = amu(4,mp) + v4*dx
      dxl = amu(1,mp+1) + v1*dy
      amx = amu(2,mp+1) + v2*dy
      dxp = amu(3,mp+1) + v3*dy
      vy = amu(4,mp+1) + v4*dy
      dx = amu(1,mp+2) + v1*dz
      dy = amu(2,mp+2) + v2*dz
      vz = amu(3,mp+2) + v3*dz
      dz = amu(4,mp+2) + v4*dz
      amu(1,mp) = dx1
      amu(2,mp) = dy1
      amu(3,mp) = amy
      amu(4,mp) = vx
      amu(1,mp+1) = dxl
      amu(2,mp+1) = amx
      amu(3,mp+1) = dxp
      amu(4,mp+1) = vy
      amu(1,mp+2) = dx
      amu(2,mp+2) = dy
      amu(3,mp+2) = vz
      amu(4,mp+2) = dz
c find inverse gamma for next particle
      gami2 = 1.0/(1.0 + p2*ci2)
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      amx = qm*(.75 - dxn*dxn)
      amy = .75 - dyn*dyn
      ml = mm + nn
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
      mp = mn + nxv
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = (vxn*vxn - vyn*vyn)*gami2
      v2 = (vxn*vyn)*gami2
      v3 = (vzn*vxn)*gami2
      v4 = (vzn*vyn)*gami2
      amu(1,mn) = amu(1,mn) + v1*dx
      amu(2,mn) = amu(2,mn) + v2*dx
      amu(3,mn) = amu(3,mn) + v3*dx
      amu(4,mn) = amu(4,mn) + v4*dx
      amu(1,mn+1) = amu(1,mn+1) + v1*dy
      amu(2,mn+1) = amu(2,mn+1) + v2*dy
      amu(3,mn+1) = amu(3,mn+1) + v3*dy
      amu(4,mn+1) = amu(4,mn+1) + v4*dy
      amu(1,mn+2) = amu(1,mn+2) + v1*dz
      amu(2,mn+2) = amu(2,mn+2) + v2*dz
      amu(3,mn+2) = amu(3,mn+2) + v3*dz
      amu(4,mn+2) = amu(4,mn+2) + v4*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml) = amu(1,ml) + v1*dx
      amu(2,ml) = amu(2,ml) + v2*dx
      amu(3,ml) = amu(3,ml) + v3*dx
      amu(4,ml) = amu(4,ml) + v4*dx
      amu(1,ml+1) = amu(1,ml+1) + v1*dy
      amu(2,ml+1) = amu(2,ml+1) + v2*dy
      amu(3,ml+1) = amu(3,ml+1) + v3*dy
      amu(4,ml+1) = amu(4,ml+1) + v4*dy
      amu(1,ml+2) = amu(1,ml+2) + v1*dz
      amu(2,ml+2) = amu(2,ml+2) + v2*dz
      amu(3,ml+2) = amu(3,ml+2) + v3*dz
      amu(4,ml+2) = amu(4,ml+2) + v4*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp) = amu(1,mp) + v1*dx
      amu(2,mp) = amu(2,mp) + v2*dx
      amu(3,mp) = amu(3,mp) + v3*dx
      amu(4,mp) = amu(4,mp) + v4*dx
      amu(1,mp+1) = amu(1,mp+1) + v1*dy
      amu(2,mp+1) = amu(2,mp+1) + v2*dy
      amu(3,mp+1) = amu(3,mp+1) + v3*dy
      amu(4,mp+1) = amu(4,mp+1) + v4*dy
      amu(1,mp+2) = amu(1,mp+2) + v1*dz
      amu(2,mp+2) = amu(2,mp+2) + v2*dz
      amu(3,mp+2) = amu(3,mp+2) + v3*dz
      amu(4,mp+2) = amu(4,mp+2) + v4*dz
      return
      end
      subroutine GRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp,n
     1op,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation for relativistic particles.
c scalar version using guard cells
c 430 flops/particle, 2 divide, 1 sqrt, 150 loads, 80 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami.
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c part(5,n) = momentum pz of particle n at t - dt/2
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c fxy(3,j+1,k+1) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j+1,k+1) = ith component of current density
c at grid point j,k for i = 1, 3
c dcu(i,j+1,k+1) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j+1,k+1) = ith component of momentum flux
c at grid point j,k for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of flux array, must be >= nx+3
c nyv = second dimension of flux array, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, p2, v1, v2, v3, v4
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nl = nn + 1
      amx = .75 - dxp*dxp
      mm = mm + 2
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))
      dz = amy*(dxl*fxy(3,nl,mm) + amx*fxy(3,nn,mm) + dxp*fxy(3,np,mm))
     1+ dyl*(dxl*fxy(3,nl,ml) + amx*fxy(3,nn,ml) + dxp*fxy(3,np,ml)) + d
     2yp*(dxl*fxy(3,nl,mp) + amx*fxy(3,nn,mp) + dxp*fxy(3,np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = amy*(dxl*bxy(1,nl,mm) + amx*bxy(1,nn,mm) + dxp*bxy(1,np,mm))
     1+ dyl*(dxl*bxy(1,nl,ml) + amx*bxy(1,nn,ml) + dxp*bxy(1,np,ml)) + d
     2yp*(dxl*bxy(1,nl,mp) + amx*bxy(1,nn,mp) + dxp*bxy(1,np,mp))              
      oy = amy*(dxl*bxy(2,nl,mm) + amx*bxy(2,nn,mm) + dxp*bxy(2,np,mm))
     1+ dyl*(dxl*bxy(2,nl,ml) + amx*bxy(2,nn,ml) + dxp*bxy(2,np,ml)) + d
     2yp*(dxl*bxy(2,nl,mp) + amx*bxy(2,nn,mp) + dxp*bxy(2,np,mp))
      oz = amy*(dxl*bxy(3,nl,mm) + amx*bxy(3,nn,mm) + dxp*bxy(3,np,mm))
     1+ dyl*(dxl*bxy(3,nl,ml) + amx*bxy(3,nn,ml) + dxp*bxy(3,np,ml)) + d
     2yp*(dxl*bxy(3,nl,mp) + amx*bxy(3,nn,mp) + dxp*bxy(3,np,mp))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      amu(3,nl,mm) = amu(3,nl,mm) + v3*dx
      amu(4,nl,mm) = amu(4,nl,mm) + v4*dx
      dcu(1,nl,mm) = dcu(1,nl,mm) + vx*dx
      dcu(2,nl,mm) = dcu(2,nl,mm) + vy*dx
      dcu(3,nl,mm) = dcu(3,nl,mm) + vz*dx
      cu(1,nl,mm) = cu(1,nl,mm) + ox*dx
      cu(2,nl,mm) = cu(2,nl,mm) + oy*dx
      cu(3,nl,mm) = cu(3,nl,mm) + oz*dx
      dx = dxl*dyl
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      dcu(3,nn,mm) = dcu(3,nn,mm) + vz*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + oz*dy
      dy = amx*dyl
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      amu(3,np,mm) = amu(3,np,mm) + v3*dz
      amu(4,np,mm) = amu(4,np,mm) + v4*dz
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dz
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dz
      dcu(3,np,mm) = dcu(3,np,mm) + vz*dz
      cu(1,np,mm) = cu(1,np,mm) + ox*dz
      cu(2,np,mm) = cu(2,np,mm) + oy*dz
      cu(3,np,mm) = cu(3,np,mm) + oz*dz
      dz = dxp*dyl
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      amu(3,nl,ml) = amu(3,nl,ml) + v3*dx
      amu(4,nl,ml) = amu(4,nl,ml) + v4*dx
      dcu(1,nl,ml) = dcu(1,nl,ml) + vx*dx
      dcu(2,nl,ml) = dcu(2,nl,ml) + vy*dx
      dcu(3,nl,ml) = dcu(3,nl,ml) + vz*dx
      cu(1,nl,ml) = cu(1,nl,ml) + ox*dx
      cu(2,nl,ml) = cu(2,nl,ml) + oy*dx
      cu(3,nl,ml) = cu(3,nl,ml) + oz*dx
      dx = dxl*dyp
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      amu(3,nn,ml) = amu(3,nn,ml) + v3*dy
      amu(4,nn,ml) = amu(4,nn,ml) + v4*dy
      dcu(1,nn,ml) = dcu(1,nn,ml) + vx*dy
      dcu(2,nn,ml) = dcu(2,nn,ml) + vy*dy
      dcu(3,nn,ml) = dcu(3,nn,ml) + vz*dy
      cu(1,nn,ml) = cu(1,nn,ml) + ox*dy
      cu(2,nn,ml) = cu(2,nn,ml) + oy*dy
      cu(3,nn,ml) = cu(3,nn,ml) + oz*dy
      dy = amx*dyp
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      amu(3,np,ml) = amu(3,np,ml) + v3*dz
      amu(4,np,ml) = amu(4,np,ml) + v4*dz
      dcu(1,np,ml) = dcu(1,np,ml) + vx*dz
      dcu(2,np,ml) = dcu(2,np,ml) + vy*dz
      dcu(3,np,ml) = dcu(3,np,ml) + vz*dz
      cu(1,np,ml) = cu(1,np,ml) + ox*dz
      cu(2,np,ml) = cu(2,np,ml) + oy*dz
      cu(3,np,ml) = cu(3,np,ml) + oz*dz
      dz = dxp*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      amu(3,nl,mp) = amu(3,nl,mp) + v3*dx
      amu(4,nl,mp) = amu(4,nl,mp) + v4*dx
      dcu(1,nl,mp) = dcu(1,nl,mp) + vx*dx
      dcu(2,nl,mp) = dcu(2,nl,mp) + vy*dx
      dcu(3,nl,mp) = dcu(3,nl,mp) + vz*dx
      cu(1,nl,mp) = cu(1,nl,mp) + ox*dx
      cu(2,nl,mp) = cu(2,nl,mp) + oy*dx
      cu(3,nl,mp) = cu(3,nl,mp) + oz*dx
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      dcu(3,nn,mp) = dcu(3,nn,mp) + vz*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + oz*dy
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      amu(3,np,mp) = amu(3,np,mp) + v3*dz
      amu(4,np,mp) = amu(4,np,mp) + v4*dz
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dz
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dz
      dcu(3,np,mp) = dcu(3,np,mp) + vz*dz
      cu(1,np,mp) = cu(1,np,mp) + ox*dz
      cu(2,np,mp) = cu(2,np,mp) + oy*dz
      cu(3,np,mp) = cu(3,np,mp) + oz*dz
   10 continue
      return
      end
      subroutine GSRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation for relativistic particles.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 430 flops/particle, 2 divide, 1 sqrt, 150 loads, 80 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami.
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c part(5,n) = momentum pz of particle n at t - dt/2
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c fxy(3,j+1,k+1) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k, for i = 1, 3
c dcu(i,n) = ith component of acceleration density at grid point j,k
c where n = j + nxv*k, for i = 1, 3
c amu(i,n) = ith component of momentum flux at grid point j,k
c where n = j + nxv*k, for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxn, dyn
      real dxp, dyp, amx, amy, dxl, dyl, dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, p2, v1, v2, v3, v4
      real dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4, dx5, dy5, dx6
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1) + .5
      mmn = part(2,j+1) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1) - float(nnn)
      dyn = part(2,j+1) - float(mmn)
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      ml = mm + nn
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mn = ml + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mp = mn + nxv
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))
      dz = amy*(dxl*fxy(3,mn) + amx*fxy(3,mn+1) + dxp*fxy(3,mn+2)) + dyl
     1*(dxl*fxy(3,ml) + amx*fxy(3,ml+1) + dxp*fxy(3,ml+2)) + dyp*(dxl*fx
     2y(3,mp) + amx*fxy(3,mp+1) + dxp*fxy(3,mp+2)) 
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = amy*(dxl*bxy(1,mn) + amx*bxy(1,mn+1) + dxp*bxy(1,mn+2)) + dyl
     1*(dxl*bxy(1,ml) + amx*bxy(1,ml+1) + dxp*bxy(1,ml+2)) + dyp*(dxl*bx
     2y(1,mp) + amx*bxy(1,mp+1) + dxp*bxy(1,mp+2))              
      oy = amy*(dxl*bxy(2,mn) + amx*bxy(2,mn+1) + dxp*bxy(2,mn+2)) + dyl
     1*(dxl*bxy(2,ml) + amx*bxy(2,ml+1) + dxp*bxy(2,ml+2)) + dyp*(dxl*bx
     2y(2,mp) + amx*bxy(2,mp+1) + dxp*bxy(2,mp+2)) 
      oz = amy*(dxl*bxy(3,mn) + amx*bxy(3,mn+1) + dxp*bxy(3,mn+2)) + dyl
     1*(dxl*bxy(3,ml) + amx*bxy(3,ml+1) + dxp*bxy(3,ml+2)) + dyp*(dxl*bx
     2y(3,mp) + amx*bxy(3,mp+1) + dxp*bxy(3,mp+2))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      dx1 = amu(1,mn) + v1*dx
      dy1 = amu(2,mn) + v2*dx
      amy = amu(3,mn) + v3*dx
      dx2 = amu(4,mn) + v4*dx
      dy2 = amu(1,mn+1) + v1*dy
      dx3 = amu(2,mn+1) + v2*dy
      dy3 = amu(3,mn+1) + v3*dy
      dx4 = amu(4,mn+1) + v4*dy
      dy4 = amu(1,mn+2) + v1*dz
      dx5 = amu(2,mn+2) + v2*dz
      dy5 = amu(3,mn+2) + v3*dz
      dx6 = amu(4,mn+2) + v4*dz
      amu(1,mn) = dx1
      amu(2,mn) = dy1
      amu(3,mn) = amy
      amu(4,mn) = dx2
      amu(1,mn+1) = dy2
      amu(2,mn+1) = dx3
      amu(3,mn+1) = dy3
      amu(4,mn+1) = dx4
      amu(1,mn+2) = dy4
      amu(2,mn+2) = dx5
      amu(3,mn+2) = dy5
      amu(4,mn+2) = dx6
      dx1 = dcu(1,mn) + vx*dx
      dy1 = dcu(2,mn) + vy*dx
      amy = dcu(3,mn) + vz*dx
      dx2 = dcu(1,mn+1) + vx*dy
      dy2 = dcu(2,mn+1) + vy*dy
      dx3 = dcu(3,mn+1) + vz*dy
      dy3 = dcu(1,mn+2) + vx*dz
      dx4 = dcu(2,mn+2) + vy*dz
      dy4 = dcu(3,mn+2) + vz*dz
      dcu(1,mn) = dx1
      dcu(2,mn) = dy1
      dcu(3,mn) = amy
      dcu(1,mn+1) = dx2
      dcu(2,mn+1) = dy2
      dcu(3,mn+1) = dx3
      dcu(1,mn+2) = dy3
      dcu(2,mn+2) = dx4
      dcu(3,mn+2) = dy4
      dx1 = cu(1,mn) + ox*dx
      dy1 = cu(2,mn) + oy*dx
      amy = cu(3,mn) + oz*dx
      dx2 = cu(1,mn+1) + ox*dy
      dy2 = cu(2,mn+1) + oy*dy
      dx3 = cu(3,mn+1) + oz*dy
      dy3 = cu(1,mn+2) + ox*dz
      dx4 = cu(2,mn+2) + oy*dz
      dy4 = cu(3,mn+2) + oz*dz
      cu(1,mn) = dx1
      cu(2,mn) = dy1
      cu(3,mn) = amy
      cu(1,mn+1) = dx2
      cu(2,mn+1) = dy2
      cu(3,mn+1) = dx3
      cu(1,mn+2) = dy3
      cu(2,mn+2) = dx4
      cu(3,mn+2) = dy4
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml) + v1*dx
      dy1 = amu(2,ml) + v2*dx
      amy = amu(3,ml) + v3*dx
      dx2 = amu(4,ml) + v4*dx
      dy2 = amu(1,ml+1) + v1*dy
      dx3 = amu(2,ml+1) + v2*dy
      dy3 = amu(3,ml+1) + v3*dy
      dx4 = amu(4,ml+1) + v4*dy
      dy4 = amu(1,ml+2) + v1*dz
      dx5 = amu(2,ml+2) + v2*dz
      dy5 = amu(3,ml+2) + v3*dz
      dx6 = amu(4,ml+2) + v4*dz
      amu(1,ml) = dx1
      amu(2,ml) = dy1
      amu(3,ml) = amy
      amu(4,ml) = dx2
      amu(1,ml+1) = dy2
      amu(2,ml+1) = dx3
      amu(3,ml+1) = dy3
      amu(4,ml+1) = dx4
      amu(1,ml+2) = dy4
      amu(2,ml+2) = dx5
      amu(3,ml+2) = dy5
      amu(4,ml+2) = dx6
      dx1 = dcu(1,ml) + vx*dx
      dy1 = dcu(2,ml) + vy*dx
      amy = dcu(3,ml) + vz*dx
      dx2 = dcu(1,ml+1) + vx*dy
      dy2 = dcu(2,ml+1) + vy*dy
      dx3 = dcu(3,ml+1) + vz*dy
      dy3 = dcu(1,ml+2) + vx*dz
      dx4 = dcu(2,ml+2) + vy*dz
      dy4 = dcu(3,ml+2) + vz*dz
      dcu(1,ml) = dx1
      dcu(2,ml) = dy1
      dcu(3,ml) = amy
      dcu(1,ml+1) = dx2
      dcu(2,ml+1) = dy2
      dcu(3,ml+1) = dx3
      dcu(1,ml+2) = dy3
      dcu(2,ml+2) = dx4
      dcu(3,ml+2) = dy4
      dx1 = cu(1,ml) + ox*dx
      dy1 = cu(2,ml) + oy*dx
      amy = cu(3,ml) + oz*dx
      dx2 = cu(1,ml+1) + ox*dy
      dy2 = cu(2,ml+1) + oy*dy
      dx3 = cu(3,ml+1) + oz*dy
      dy3 = cu(1,ml+2) + ox*dz
      dx4 = cu(2,ml+2) + oy*dz
      dy4 = cu(3,ml+2) + oz*dz
      cu(1,ml) = dx1
      cu(2,ml) = dy1
      cu(3,ml) = amy
      cu(1,ml+1) = dx2
      cu(2,ml+1) = dy2
      cu(3,ml+1) = dx3
      cu(1,ml+2) = dy3
      cu(2,ml+2) = dx4
      cu(3,ml+2) = dy4
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp) + v1*dx
      dy1 = amu(2,mp) + v2*dx
      amy = amu(3,mp) + v3*dx
      dx2 = amu(4,mp) + v4*dx
      dy2 = amu(1,mp+1) + v1*dy
      dx3 = amu(2,mp+1) + v2*dy
      dy3 = amu(3,mp+1) + v3*dy
      dx4 = amu(4,mp+1) + v4*dy
      dy4 = amu(1,mp+2) + v1*dz
      dx5 = amu(2,mp+2) + v2*dz
      dy5 = amu(3,mp+2) + v3*dz
      dx6 = amu(4,mp+2) + v4*dz
      amu(1,mp) = dx1
      amu(2,mp) = dy1
      amu(3,mp) = amy
      amu(4,mp) = dx2
      amu(1,mp+1) = dy2
      amu(2,mp+1) = dx3
      amu(3,mp+1) = dy3
      amu(4,mp+1) = dx4
      amu(1,mp+2) = dy4
      amu(2,mp+2) = dx5
      amu(3,mp+2) = dy5
      amu(4,mp+2) = dx6
      dx1 = dcu(1,mp) + vx*dx
      dy1 = dcu(2,mp) + vy*dx
      amy = dcu(3,mp) + vz*dx
      dx2 = dcu(1,mp+1) + vx*dy
      dy2 = dcu(2,mp+1) + vy*dy
      dx3 = dcu(3,mp+1) + vz*dy
      dy3 = dcu(1,mp+2) + vx*dz
      dx4 = dcu(2,mp+2) + vy*dz
      dy4 = dcu(3,mp+2) + vz*dz
      dcu(1,mp) = dx1
      dcu(2,mp) = dy1
      dcu(3,mp) = amy
      dcu(1,mp+1) = dx2
      dcu(2,mp+1) = dy2
      dcu(3,mp+1) = dx3
      dcu(1,mp+2) = dy3
      dcu(2,mp+2) = dx4
      dcu(3,mp+2) = dy4
      dx1 = cu(1,mp) + ox*dx
      dy1 = cu(2,mp) + oy*dx
      amy = cu(3,mp) + oz*dx
      dx2 = cu(1,mp+1) + ox*dy
      dy2 = cu(2,mp+1) + oy*dy
      dx3 = cu(3,mp+1) + oz*dy
      dy3 = cu(1,mp+2) + ox*dz
      dx4 = cu(2,mp+2) + oy*dz
      dy4 = cu(3,mp+2) + oz*dz
      cu(1,mp) = dx1
      cu(2,mp) = dy1
      cu(3,mp) = amy
      cu(1,mp+1) = dx2
      cu(2,mp+1) = dy2
      cu(3,mp+1) = dx3
      cu(1,mp+2) = dy3
      cu(2,mp+2) = dx4
      cu(3,mp+2) = dy4
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      ml = mm + nn
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
      mp = mn + nxv
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
      dz = amy*(dxl*fxy(3,mn) + amx*fxy(3,mn+1) + dxp*fxy(3,mn+2)) + dyl
     1*(dxl*fxy(3,ml) + amx*fxy(3,ml+1) + dxp*fxy(3,ml+2)) + dyp*(dxl*fx
     2y(3,mp) + amx*fxy(3,mp+1) + dxp*fxy(3,mp+2))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = amy*(dxl*bxy(1,mn) + amx*bxy(1,mn+1) + dxp*bxy(1,mn+2)) + dyl
     1*(dxl*bxy(1,ml) + amx*bxy(1,ml+1) + dxp*bxy(1,ml+2)) + dyp*(dxl*bx
     2y(1,mp) + amx*bxy(1,mp+1) + dxp*bxy(1,mp+2))              
      oy = amy*(dxl*bxy(2,mn) + amx*bxy(2,mn+1) + dxp*bxy(2,mn+2)) + dyl
     1*(dxl*bxy(2,ml) + amx*bxy(2,ml+1) + dxp*bxy(2,ml+2)) + dyp*(dxl*bx
     2y(2,mp) + amx*bxy(2,mp+1) + dxp*bxy(2,mp+2))  
      oz = amy*(dxl*bxy(3,mn) + amx*bxy(3,mn+1) + dxp*bxy(3,mn+2)) + dyl
     1*(dxl*bxy(3,ml) + amx*bxy(3,ml+1) + dxp*bxy(3,ml+2)) + dyp*(dxl*bx
     2y(3,mp)+ amx*bxy(3,mp+1) + dxp*bxy(3,mp+2)) 
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      amu(1,mn) = amu(1,mn) + v1*dx
      amu(2,mn) = amu(2,mn) + v2*dx
      amu(3,mn) = amu(3,mn) + v3*dx
      amu(4,mn) = amu(4,mn) + v4*dx
      amu(1,mn+1) = amu(1,mn+1) + v1*dy
      amu(2,mn+1) = amu(2,mn+1) + v2*dy
      amu(3,mn+1) = amu(3,mn+1) + v3*dy
      amu(4,mn+1) = amu(4,mn+1) + v4*dy
      amu(1,mn+2) = amu(1,mn+2) + v1*dz
      amu(2,mn+2) = amu(2,mn+2) + v2*dz
      amu(3,mn+2) = amu(3,mn+2) + v3*dz
      amu(4,mn+2) = amu(4,mn+2) + v4*dz
      dcu(1,mn) = dcu(1,mn) + vx*dx
      dcu(2,mn) = dcu(2,mn) + vy*dx
      dcu(3,mn) = dcu(3,mn) + vz*dx
      dcu(1,mn+1) = dcu(1,mn+1) + vx*dy
      dcu(2,mn+1) = dcu(2,mn+1) + vy*dy
      dcu(3,mn+1) = dcu(3,mn+1) + vz*dy
      dcu(1,mn+2) = dcu(1,mn+2) + vx*dz
      dcu(2,mn+2) = dcu(2,mn+2) + vy*dz
      dcu(3,mn+2) = dcu(3,mn+2) + vz*dz
      cu(1,mn) = cu(1,mn) + ox*dx
      cu(2,mn) = cu(2,mn) + oy*dx
      cu(3,mn) = cu(3,mn) + oz*dx
      cu(1,mn+1) = cu(1,mn+1) + ox*dy
      cu(2,mn+1) = cu(2,mn+1) + oy*dy
      cu(3,mn+1) = cu(3,mn+1) + oz*dy
      cu(1,mn+2) = cu(1,mn+2) + ox*dz
      cu(2,mn+2) = cu(2,mn+2) + oy*dz
      cu(3,mn+2) = cu(3,mn+2) + oz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml) = amu(1,ml) + v1*dx
      amu(2,ml) = amu(2,ml) + v2*dx
      amu(3,ml) = amu(3,ml) + v3*dx
      amu(4,ml) = amu(4,ml) + v4*dx
      amu(1,ml+1) = amu(1,ml+1) + v1*dy
      amu(2,ml+1) = amu(2,ml+1) + v2*dy
      amu(3,ml+1) = amu(3,ml+1) + v3*dy
      amu(4,ml+1) = amu(4,ml+1) + v4*dy
      amu(1,ml+2) = amu(1,ml+2) + v1*dz
      amu(2,ml+2) = amu(2,ml+2) + v2*dz
      amu(3,ml+2) = amu(3,ml+2) + v3*dz
      amu(4,ml+2) = amu(4,ml+2) + v4*dz
      dcu(1,ml) = dcu(1,ml) + vx*dx
      dcu(2,ml) = dcu(2,ml) + vy*dx
      dcu(3,ml) = dcu(3,ml) + vz*dx
      dcu(1,ml+1) = dcu(1,ml+1) + vx*dy
      dcu(2,ml+1) = dcu(2,ml+1) + vy*dy
      dcu(3,ml+1) = dcu(3,ml+1) + vz*dy
      dcu(1,ml+2) = dcu(1,ml+2) + vx*dz
      dcu(2,ml+2) = dcu(2,ml+2) + vy*dz
      dcu(3,ml+2) = dcu(3,ml+2) + vz*dz
      cu(1,ml) = cu(1,ml) + ox*dx
      cu(2,ml) = cu(2,ml) + oy*dx
      cu(3,ml) = cu(3,ml) + oz*dx
      cu(1,ml+1) = cu(1,ml+1) + ox*dy
      cu(2,ml+1) = cu(2,ml+1) + oy*dy
      cu(3,ml+1) = cu(3,ml+1) + oz*dy
      cu(1,ml+2) = cu(1,ml+2) + ox*dz
      cu(2,ml+2) = cu(2,ml+2) + oy*dz
      cu(3,ml+2) = cu(3,ml+2) + oz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp) = amu(1,mp) + v1*dx
      amu(2,mp) = amu(2,mp) + v2*dx
      amu(3,mp) = amu(3,mp) + v3*dx
      amu(4,mp) = amu(4,mp) + v4*dx
      amu(1,mp+1) = amu(1,mp+1) + v1*dy
      amu(2,mp+1) = amu(2,mp+1) + v2*dy
      amu(3,mp+1) = amu(3,mp+1) + v3*dy
      amu(4,mp+1) = amu(4,mp+1) + v4*dy
      amu(1,mp+2) = amu(1,mp+2) + v1*dz
      amu(2,mp+2) = amu(2,mp+2) + v2*dz
      amu(3,mp+2) = amu(3,mp+2) + v3*dz
      amu(4,mp+2) = amu(4,mp+2) + v4*dz
      dcu(1,mp) = dcu(1,mp) + vx*dx
      dcu(2,mp) = dcu(2,mp) + vy*dx
      dcu(3,mp) = dcu(3,mp) + vz*dx
      dcu(1,mp+1) = dcu(1,mp+1) + vx*dy
      dcu(2,mp+1) = dcu(2,mp+1) + vy*dy
      dcu(3,mp+1) = dcu(3,mp+1) + vz*dy
      dcu(1,mp+2) = dcu(1,mp+2) + vx*dz
      dcu(2,mp+2) = dcu(2,mp+2) + vy*dz
      dcu(3,mp+2) = dcu(3,mp+2) + vz*dz
      cu(1,mp) = cu(1,mp) + ox*dx
      cu(2,mp) = cu(2,mp) + oy*dx
      cu(3,mp) = cu(3,mp) + oz*dx
      cu(1,mp+1) = cu(1,mp+1) + ox*dy
      cu(2,mp+1) = cu(2,mp+1) + oy*dy
      cu(3,mp+1) = cu(3,mp+1) + oz*dy
      cu(1,mp+2) = cu(1,mp+2) + ox*dz
      cu(2,mp+2) = cu(2,mp+2) + oy*dz
      cu(3,mp+2) = cu(3,mp+2) + oz*dz
      return
      end
      subroutine GRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation for relativistic particles
c scalar version using guard cells
c 81 flops/particle, 1 divide, 22 loads, 18 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy, for i = 1, 2
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c amu(i,j+1,k+1) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of flux array, must be >= nx+3
c nyv = second dimension of flux array, must be >= ny+3
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qmh, ci2, gami2, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, p2, v1, v2
      qmh = .5*qm
      ci2 = ci*ci
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami2 = 1.0/(1.0 + p2*ci2)
c calculate weights
      nl = nn + 1
      amx = qm*(.75 - dxp*dxp)
      mm = mm + 2
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = qmh*(.5 - dxp)**2
      np = nl + 2
      dxp = qmh*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      dx = dxl*dyl
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      dy = amx*dyl
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      dz = dxp*dyl
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      dx = dxl*dyp
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      dy = amx*dyp
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      dz = dxp*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
   10 continue
      return
      end
      subroutine GSRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle momentum flux
c using second-order spline interpolation for relativistic particles
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 81 flops/particle, 1 divide, 22 loads, 18 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy, for i = 1, 2
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c amu(i,n) = ith component of momentum flux at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyv = dimension of current array, must be >= nxv*(ny+3)
      implicit none
      integer nop, idimp, nxv, nxyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(2,nxyv)
      integer nnn, mmn, j, nn, mm, ml, mn, mp
      real dxn, dyn, qmh, ci2, gami2, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, p2, v1, v2
      real dx1, dy1, dx2, dy2
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      qmh = .5*qm
      ci2 = ci*ci
c find inverse gamma
      vx = part(3,1)
      vy = part(4,1)
      p2 = vx*vx + vy*vy
      gami2 = 1.0/(1.0 + p2*ci2)
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j) + .5
      mmn = part(2,j) + .5
      dxp = dxn
      dyp = dyn
      dxn = part(1,j) - float(nnn)
      dyn = part(2,j) - float(mmn)
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = mm + nn
      dxl = qmh*(.5 - dxp)**2
      dxp = qmh*(.5 + dxp)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv
c calculate weights
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
c get momentum for next particle
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
c deposit momentum flux
      dx1 = amu(1,mn) + v1*dx
      dy1 = amu(2,mn) + v2*dx
      dx2 = amu(1,mn+1) + v1*dy
      dy2 = amu(2,mn+1) + v2*dy
      dx = amu(1,mn+2) + v1*dz
      dy = amu(2,mn+2) + v2*dz
      amu(1,mn) = dx1
      amu(2,mn) = dy1
      amu(1,mn+1) = dx2
      amu(2,mn+1) = dy2
      amu(1,mn+2) = dx
      amu(2,mn+2) = dy
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml) + v1*dx
      dy1 = amu(2,ml) + v2*dx
      dx2 = amu(1,ml+1) + v1*dy
      dy2 = amu(2,ml+1) + v2*dy
      dx = amu(1,ml+2) + v1*dz
      dy = amu(2,ml+2) + v2*dz
      amu(1,ml) = dx1
      amu(2,ml) = dy1
      amu(1,ml+1) = dx2
      amu(2,ml+1) = dy2
      amu(1,ml+2) = dx
      amu(2,ml+2) = dy
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp) + v1*dx
      dy1 = amu(2,mp) + v2*dx
      dxl = amu(1,mp+1) + v1*dy
      amx = amu(2,mp+1) + v2*dy
      dx = amu(1,mp+2) + v1*dz
      dy = amu(2,mp+2) + v2*dz
      amu(1,mp) = dx1
      amu(2,mp) = dy1
      amu(1,mp+1) = dxl
      amu(2,mp+1) = amx
      amu(1,mp+2) = dx
      amu(2,mp+2) = dy
c find inverse gamma for next particle
      gami2 = 1.0/(1.0 + p2*ci2)
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      amx = qm*(.75 - dxn*dxn)
      amy = .75 - dyn*dyn
      ml = mm + nn
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
      mp = mn + nxv
c deposit momentum flux
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      amu(1,mn) = amu(1,mn) + v1*dx
      amu(2,mn) = amu(2,mn) + v2*dx
      amu(1,mn+1) = amu(1,mn+1) + v1*dy
      amu(2,mn+1) = amu(2,mn+1) + v2*dy
      amu(1,mn+2) = amu(1,mn+2) + v1*dz
      amu(2,mn+2) = amu(2,mn+2) + v2*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml) = amu(1,ml) + v1*dx
      amu(2,ml) = amu(2,ml) + v2*dx
      amu(1,ml+1) = amu(1,ml+1) + v1*dy
      amu(2,ml+1) = amu(2,ml+1) + v2*dy
      amu(1,ml+2) = amu(1,ml+2) + v1*dz
      amu(2,ml+2) = amu(2,ml+2) + v2*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp) = amu(1,mp) + v1*dx
      amu(2,mp) = amu(2,mp) + v2*dx
      amu(1,mp+1) = amu(1,mp+1) + v1*dy
      amu(2,mp+1) = amu(2,mp+1) + v2*dy
      amu(1,mp+2) = amu(1,mp+2) + v1*dz
      amu(2,mp+2) = amu(2,mp+2) + v2*dz
      return
      end
      subroutine GRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp,n
     1op,nxv,nyv)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation for relativistic particles.
c scalar version using guard cells
c 252 flops/particle, 2 divide, 1 sqrt, 58 loads, 54 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj*gami/dt, where j = x,y, for i = 1, 2
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy, for i = 1, 2
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = -rot(2)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t))*gami, and
c gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j+1,k+1) = ith component of current density
c at grid point j,k for i = 1, 2
c dcu(i,j+1,k+1) = ith component of acceleration density
c at grid point j,k for i = 1, 2
c amu(i,j+1,k+1) = ith component of momentum flux
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, p2, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nl = nn + 1
      amx = .75 - dxp*dxp
      mm = mm + 2
      amy = .75 - dyp*dyp
      nn = nl + 1
      dxl = .5*(.5 - dxp)**2
      np = nl + 2
      dxp = .5*(.5 + dxp)**2
      ml = mm - 1
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      dyp = .5*(.5 + dyp)**2
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = amy*(dxl*bz(nl,mm) + amx*bz(nn,mm) + dxp*bz(np,mm)) + dyl*(dx
     1l*bz(nl,ml) + amx*bz(nn,ml) + dxp*bz(np,ml)) + dyp*(dxl*bz(nl,mp) 
     2+ amx*bz(nn,mp) + dxp*bz(np,mp))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      dcu(1,nl,mm) = dcu(1,nl,mm) + vx*dx
      dcu(2,nl,mm) = dcu(2,nl,mm) + vy*dx
      cu(1,nl,mm) = cu(1,nl,mm) + ox*dx
      cu(2,nl,mm) = cu(2,nl,mm) + oy*dx
      dx = dxl*dyl
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
      dy = amx*dyl
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dz
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dz
      cu(1,np,mm) = cu(1,np,mm) + ox*dz
      cu(2,np,mm) = cu(2,np,mm) + oy*dz
      dz = dxp*dyl
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      dcu(1,nl,ml) = dcu(1,nl,ml) + vx*dx
      dcu(2,nl,ml) = dcu(2,nl,ml) + vy*dx
      cu(1,nl,ml) = cu(1,nl,ml) + ox*dx
      cu(2,nl,ml) = cu(2,nl,ml) + oy*dx
      dx = dxl*dyp
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      dcu(1,nn,ml) = dcu(1,nn,ml) + vx*dy
      dcu(2,nn,ml) = dcu(2,nn,ml) + vy*dy
      cu(1,nn,ml) = cu(1,nn,ml) + ox*dy
      cu(2,nn,ml) = cu(2,nn,ml) + oy*dy
      dy = amx*dyp
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      dcu(1,np,ml) = dcu(1,np,ml) + vx*dz
      dcu(2,np,ml) = dcu(2,np,ml) + vy*dz
      cu(1,np,ml) = cu(1,np,ml) + ox*dz
      cu(2,np,ml) = cu(2,np,ml) + oy*dz
      dz = dxp*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      dcu(1,nl,mp) = dcu(1,nl,mp) + vx*dx
      dcu(2,nl,mp) = dcu(2,nl,mp) + vy*dx
      cu(1,nl,mp) = cu(1,nl,mp) + ox*dx
      cu(2,nl,mp) = cu(2,nl,mp) + oy*dx
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dz
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dz
      cu(1,np,mp) = cu(1,np,mp) + ox*dz
      cu(2,np,mp) = cu(2,np,mp) + oy*dz
   10 continue
      return
      end
      subroutine GSRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nxyv)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using second-order spline
c interpolation for relativistic particles.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 252 flops/particle, 2 divide, 1 sqrt, 58 loads, 54 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c dcu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c dcu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c dcu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c dcu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c dcu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c dcu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c dcu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c dcu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*dvj*gami/dt, where j = x,y, for i = 1, 2
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(.75-dx**2)*(.75-dy**2)
c amu(i,n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c amu(i,n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c amu(i,n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c amu(i,n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c amu(i,n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c amu(i,n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c amu(i,n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c amu(i,n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy, for i = 1, 2
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = -rot(2)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t))*gami, and
c gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j+1,k+1) = ith component of current density
c at grid point j,k for i = 1, 2
c dcu(i,j+1,k+1) = ith component of acceleration density
c at grid point j,k for i = 1, 2
c amu(i,j+1,k+1) = ith component of momentum flux
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      dimension cu(2,nxyv), dcu(2,nxyv), amu(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxn, dyn
      real dxp, dyp, amx, amy, dxl, dyl, dx, dy, dz, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, p2, v1, v2
      real dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4, dx5
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1) + .5
      mmn = part(2,j+1) + .5
      dx = dxn
      dy = dyn
      dxn = part(1,j+1) - float(nnn)
      dyn = part(2,j+1) - float(mmn)
      amx = .75 - dx*dx
      amy = .75 - dy*dy
      ml = mm + nn
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      mn = ml + nxv
      dyl = .5*(.5 - dy)**2
      dyp = .5*(.5 + dy)**2
      mp = mn + nxv
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = amy*(dxl*bz(mn) + amx*bz(mn+1) + dxp*bz(mn+2)) + dyl*(dxl*bz(
     1ml) + amx*bz(ml+1) + dxp*bz(ml+2)) + dyp*(dxl*bz(mp) + amx*bz(mp+1
     2) + dxp*bz(mp+2))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      dx1 = amu(1,mn) + v1*dx
      dy1 = amu(2,mn) + v2*dx
      dy2 = amu(1,mn+1) + v1*dy
      dx3 = amu(2,mn+1) + v2*dy
      dy4 = amu(1,mn+2) + v1*dz
      dx5 = amu(2,mn+2) + v2*dz
      amu(1,mn) = dx1
      amu(2,mn) = dy1
      amu(1,mn+1) = dy2
      amu(2,mn+1) = dx3
      amu(1,mn+2) = dy4
      amu(2,mn+2) = dx5
      dx1 = dcu(1,mn) + vx*dx
      dy1 = dcu(2,mn) + vy*dx
      dx2 = dcu(1,mn+1) + vx*dy
      dy2 = dcu(2,mn+1) + vy*dy
      dy3 = dcu(1,mn+2) + vx*dz
      dx4 = dcu(2,mn+2) + vy*dz
      dcu(1,mn) = dx1
      dcu(2,mn) = dy1
      dcu(1,mn+1) = dx2
      dcu(2,mn+1) = dy2
      dcu(1,mn+2) = dy3
      dcu(2,mn+2) = dx4
      dx1 = cu(1,mn) + ox*dx
      dy1 = cu(2,mn) + oy*dx
      dx2 = cu(1,mn+1) + ox*dy
      dy2 = cu(2,mn+1) + oy*dy
      dy3 = cu(1,mn+2) + ox*dz
      dx4 = cu(2,mn+2) + oy*dz
      cu(1,mn) = dx1
      cu(2,mn) = dy1
      cu(1,mn+1) = dx2
      cu(2,mn+1) = dy2
      cu(1,mn+2) = dy3
      cu(2,mn+2) = dx4
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = amu(1,ml) + v1*dx
      dy1 = amu(2,ml) + v2*dx
      dy2 = amu(1,ml+1) + v1*dy
      dx3 = amu(2,ml+1) + v2*dy
      dy4 = amu(1,ml+2) + v1*dz
      dx5 = amu(2,ml+2) + v2*dz
      amu(1,ml) = dx1
      amu(2,ml) = dy1
      amu(1,ml+1) = dy2
      amu(2,ml+1) = dx3
      amu(1,ml+2) = dy4
      amu(2,ml+2) = dx5
      dx1 = dcu(1,ml) + vx*dx
      dy1 = dcu(2,ml) + vy*dx
      dx2 = dcu(1,ml+1) + vx*dy
      dy2 = dcu(2,ml+1) + vy*dy
      dy3 = dcu(1,ml+2) + vx*dz
      dx4 = dcu(2,ml+2) + vy*dz
      dcu(1,ml) = dx1
      dcu(2,ml) = dy1
      dcu(1,ml+1) = dx2
      dcu(2,ml+1) = dy2
      dcu(1,ml+2) = dy3
      dcu(2,ml+2) = dx4
      dx1 = cu(1,ml) + ox*dx
      dy1 = cu(2,ml) + oy*dx
      dx2 = cu(1,ml+1) + ox*dy
      dy2 = cu(2,ml+1) + oy*dy
      dy3 = cu(1,ml+2) + ox*dz
      dx4 = cu(2,ml+2) + oy*dz
      cu(1,ml) = dx1
      cu(2,ml) = dy1
      cu(1,ml+1) = dx2
      cu(2,ml+1) = dy2
      cu(1,ml+2) = dy3
      cu(2,ml+2) = dx4
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = amu(1,mp) + v1*dx
      dy1 = amu(2,mp) + v2*dx
      dy2 = amu(1,mp+1) + v1*dy
      dx3 = amu(2,mp+1) + v2*dy
      dy4 = amu(1,mp+2) + v1*dz
      dx5 = amu(2,mp+2) + v2*dz
      amu(1,mp) = dx1
      amu(2,mp) = dy1
      amu(1,mp+1) = dy2
      amu(2,mp+1) = dx3
      amu(1,mp+2) = dy4
      amu(2,mp+2) = dx5
      dx1 = dcu(1,mp) + vx*dx
      dy1 = dcu(2,mp) + vy*dx
      dx2 = dcu(1,mp+1) + vx*dy
      dy2 = dcu(2,mp+1) + vy*dy
      dy3 = dcu(1,mp+2) + vx*dz
      dx4 = dcu(2,mp+2) + vy*dz
      dcu(1,mp) = dx1
      dcu(2,mp) = dy1
      dcu(1,mp+1) = dx2
      dcu(2,mp+1) = dy2
      dcu(1,mp+2) = dy3
      dcu(2,mp+2) = dx4
      dx1 = cu(1,mp) + ox*dx
      dy1 = cu(2,mp) + oy*dx
      dx2 = cu(1,mp+1) + ox*dy
      dy2 = cu(2,mp+1) + oy*dy
      dy3 = cu(1,mp+2) + ox*dz
      dx4 = cu(2,mp+2) + oy*dz
      cu(1,mp) = dx1
      cu(2,mp) = dy1
      cu(1,mp+1) = dx2
      cu(2,mp+1) = dy2
      cu(1,mp+2) = dy3
      cu(2,mp+2) = dx4
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      amx = .75 - dxn*dxn
      amy = .75 - dyn*dyn
      ml = mm + nn
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyn)**2
      dyp = .5*(.5 + dyn)**2
      mp = mn + nxv
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,nop)
      vy = part(4,nop)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = amy*(dxl*bz(mn) + amx*bz(mn+1) + dxp*bz(mn+2)) + dyl*(dxl*bz(
     1ml) + amx*bz(ml+1) + dxp*bz(ml+2)) + dyp*(dxl*bz(mp)+ amx*bz(mp+1)
     2 + dxp*bz(mp+2)) 
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      amu(1,mn) = amu(1,mn) + v1*dx
      amu(2,mn) = amu(2,mn) + v2*dx
      amu(1,mn+1) = amu(1,mn+1) + v1*dy
      amu(2,mn+1) = amu(2,mn+1) + v2*dy
      amu(1,mn+2) = amu(1,mn+2) + v1*dz
      amu(2,mn+2) = amu(2,mn+2) + v2*dz
      dcu(1,mn) = dcu(1,mn) + vx*dx
      dcu(2,mn) = dcu(2,mn) + vy*dx
      dcu(1,mn+1) = dcu(1,mn+1) + vx*dy
      dcu(2,mn+1) = dcu(2,mn+1) + vy*dy
      dcu(1,mn+2) = dcu(1,mn+2) + vx*dz
      dcu(2,mn+2) = dcu(2,mn+2) + vy*dz
      cu(1,mn) = cu(1,mn) + ox*dx
      cu(2,mn) = cu(2,mn) + oy*dx
      cu(1,mn+1) = cu(1,mn+1) + ox*dy
      cu(2,mn+1) = cu(2,mn+1) + oy*dy
      cu(1,mn+2) = cu(1,mn+2) + ox*dz
      cu(2,mn+2) = cu(2,mn+2) + oy*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      amu(1,ml) = amu(1,ml) + v1*dx
      amu(2,ml) = amu(2,ml) + v2*dx
      amu(1,ml+1) = amu(1,ml+1) + v1*dy
      amu(2,ml+1) = amu(2,ml+1) + v2*dy
      amu(1,ml+2) = amu(1,ml+2) + v1*dz
      amu(2,ml+2) = amu(2,ml+2) + v2*dz
      dcu(1,ml) = dcu(1,ml) + vx*dx
      dcu(2,ml) = dcu(2,ml) + vy*dx
      dcu(1,ml+1) = dcu(1,ml+1) + vx*dy
      dcu(2,ml+1) = dcu(2,ml+1) + vy*dy
      dcu(1,ml+2) = dcu(1,ml+2) + vx*dz
      dcu(2,ml+2) = dcu(2,ml+2) + vy*dz
      cu(1,ml) = cu(1,ml) + ox*dx
      cu(2,ml) = cu(2,ml) + oy*dx
      cu(1,ml+1) = cu(1,ml+1) + ox*dy
      cu(2,ml+1) = cu(2,ml+1) + oy*dy
      cu(1,ml+2) = cu(1,ml+2) + ox*dz
      cu(2,ml+2) = cu(2,ml+2) + oy*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      amu(1,mp) = amu(1,mp) + v1*dx
      amu(2,mp) = amu(2,mp) + v2*dx
      amu(1,mp+1) = amu(1,mp+1) + v1*dy
      amu(2,mp+1) = amu(2,mp+1) + v2*dy
      amu(1,mp+2) = amu(1,mp+2) + v1*dz
      amu(2,mp+2) = amu(2,mp+2) + v2*dz
      dcu(1,mp) = dcu(1,mp) + vx*dx
      dcu(2,mp) = dcu(2,mp) + vy*dx
      dcu(1,mp+1) = dcu(1,mp+1) + vx*dy
      dcu(2,mp+1) = dcu(2,mp+1) + vy*dy
      dcu(1,mp+2) = dcu(1,mp+2) + vx*dz
      dcu(2,mp+2) = dcu(2,mp+2) + vy*dz
      cu(1,mp) = cu(1,mp) + ox*dx
      cu(2,mp) = cu(2,mp) + oy*dx
      cu(1,mp+1) = cu(1,mp+1) + ox*dy
      cu(2,mp+1) = cu(2,mp+1) + oy*dy
      cu(1,mp+2) = cu(1,mp+2) + ox*dz
      cu(2,mp+2) = cu(2,mp+2) + oy*dz
      return
      end
      subroutine GRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation for relativistic particles
c scalar version using guard cells
c 62 flops/particle, 1 divide, 21 loads, 16 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,j,k) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nxv = first dimension of flux array, must be >= nx+1
c nyv = second dimension of flux array, must be >= ny+1
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(4,nxv,nyv)
      integer j, nn, mm, np, mp
      real ci2, gami2, dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz, p2, v1, v2, v3, v4
      ci2 = ci*ci
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = 1.0/(1.0 + p2*ci2)
c calculate weights
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit momentum flux
      dx = dxp*dyp
      dy = amx*dyp
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vz*vx)*gami2
      v4 = (vz*vy)*gami2
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      amu(3,np,mp) = amu(3,np,mp) + v3*dx
      amu(4,np,mp) = amu(4,np,mp) + v4*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      amu(3,np,mm) = amu(3,np,mm) + v3*dx
      amu(4,np,mm) = amu(4,np,mm) + v4*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
   10 continue
      return
      end
      subroutine GSRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using first-order linear interpolation for relativistic particles
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 62 flops/particle, 1 divide, 21 loads, 16 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,n) = ith component of momentum flux at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
      implicit none
      integer nop, idimp, nxv, nxyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(4,nxyv)
      integer nnn, mmn, nn, mm, mp, j
      real dxn, dyn, ci2, gami2, dxp, dyp, amx, amy, dx, dy, dz
      real vxn, vyn, vzn, vx, vy, p2, v1, v2, v3, v4, dx1, dy1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      ci2 = ci*ci
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      vzn = part(5,1)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami2 = 1.0/(1.0 + p2*ci2)
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - float(nnn)
      dyn = part(2,j) - float(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
c calculate weights
      dx = dxp*dyp
      dz = amx*dyp
      v1 = (vxn*vxn - vyn*vyn)*gami2
      v2 = (vxn*vyn)*gami2
      v3 = (vzn*vxn)*gami2
      v4 = (vzn*vyn)*gami2
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      vzn = part(5,j)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit momentum flux
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dyp = amu(3,mp+1) + v3*dx
      vx = amu(4,mp+1) + v4*dx
      dx = amu(1,mp) + v1*dz
      dy = amu(2,mp) + v2*dz
      vy = amu(3,mp) + v3*dz
      dz = amu(4,mp) + v4*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(3,mp+1) = dyp
      amu(4,mp+1) = vx
      amu(1,mp) = dx
      amu(2,mp) = dy
      amu(3,mp) = vy
      amu(4,mp) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dyp = amu(3,mm+1) + v3*dx
      vx = amu(4,mm+1) + v4*dx
      dx = amu(1,mm) + v1*dz
      dy = amu(2,mm) + v2*dz
      vy = amu(3,mm) + v3*dz
      dz = amu(4,mm) + v4*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(3,mm+1) = dyp
      amu(4,mm+1) = vx
      amu(1,mm) = dx
      amu(2,mm) = dy
      amu(3,mm) = vy
      amu(4,mm) = dz
c find inverse gamma for next particle
      gami2 = 1.0/(1.0 + p2*ci2)
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit momentum flux
      dx = dxp*dyn
      dy = amx*dyn
      v1 = (vxn*vxn - vyn*vyn)*gami2
      v2 = (vxn*vyn)*gami2
      v3 = (vzn*vxn)*gami2
      v4 = (vzn*vyn)*gami2
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(3,mp+1) = amu(3,mp+1) + v3*dx
      amu(4,mp+1) = amu(4,mp+1) + v4*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      amu(3,mp) = amu(3,mp) + v3*dy
      amu(4,mp) = amu(4,mp) + v4*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(3,mm+1) = amu(3,mm+1) + v3*dx
      amu(4,mm+1) = amu(4,mm+1) + v4*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      amu(3,mm) = amu(3,mm) + v3*dy
      amu(4,mm) = amu(4,mm) + v4*dy
      return
      end
      subroutine GRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation for relativistic particles.
c scalar version using guard cells
c 241 flops/particle, 2 divide, 1 sqrt, 69 loads, 40 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami.
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c part(5,n) = momentum pz of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j,k) = ith component of current density
c at grid point j,k for i = 1, 3
c dcu(i,j,k) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j,k) = ith component of momentum flux
c at grid point j,k for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, p2, v1, v2, v3, v4
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp)) + amy*(dxp*fxy(3,np
     1,mm) + amx*fxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp)) + amy*(dxp*bxy(1,np
     1,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp)) + amy*(dxp*bxy(2,np
     1,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp)) + amy*(dxp*bxy(3,np
     1,mm) + amx*bxy(3,nn,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxp*dyp
      dy = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      amu(3,np,mp) = amu(3,np,mp) + v3*dx
      amu(4,np,mp) = amu(4,np,mp) + v4*dx
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dx
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dx
      dcu(3,np,mp) = dcu(3,np,mp) + vz*dx
      cu(1,np,mp) = cu(1,np,mp) + ox*dx
      cu(2,np,mp) = cu(2,np,mp) + oy*dx
      cu(3,np,mp) = cu(3,np,mp) + oz*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      dcu(3,nn,mp) = dcu(3,nn,mp) + vz*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + oz*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      amu(3,np,mm) = amu(3,np,mm) + v3*dx
      amu(4,np,mm) = amu(4,np,mm) + v4*dx
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dx
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dx
      dcu(3,np,mm) = dcu(3,np,mm) + vz*dx
      cu(1,np,mm) = cu(1,np,mm) + ox*dx
      cu(2,np,mm) = cu(2,np,mm) + oy*dx
      cu(3,np,mm) = cu(3,np,mm) + oz*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      dcu(3,nn,mm) = dcu(3,nn,mm) + vz*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + oz*dy
   10 continue
      return
      end
      subroutine GSRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nxyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using first-order spline
c interpolation for relativistic particles.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 241 flops/particle, 2 divide, 1 sqrt, 69 loads, 40 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c ddcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami.
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c part(5,n) = momentum pz of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c dcu(i,n) = ith component of acceleration density at grid point j,k
c where n = j + nxv*(k-1)
c amu(i,n) = ith component of momentum at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      integer nnn, mmn, nop1, j, nn, mm, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, p2, v1, v2, v3, v4, dx1, dy1, dx2, dy2, dx3, dy3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - float(nnn)
      dyn = part(2,j+1) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxp*dyp
      dz = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dyp = amu(3,mp+1) + v3*dx
      dx2 = amu(4,mp+1) + v4*dx
      dy2 = amu(1,mp) + v1*dz
      dx3 = amu(2,mp) + v2*dz
      dy3 = amu(3,mp) + v3*dz
      dy = amu(4,mp) + v4*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(3,mp+1) = dyp
      amu(4,mp+1) = dx2
      amu(1,mp) = dy2
      amu(2,mp) = dx3
      amu(3,mp) = dy3
      amu(4,mp) = dy
      dx1 = dcu(1,mp+1) + vx*dx
      dy1 = dcu(2,mp+1) + vy*dx
      dyp = dcu(3,mp+1) + vz*dx
      dx2 = dcu(1,mp) + vx*dz
      dy2 = dcu(2,mp) + vy*dz
      dy = dcu(3,mp) + vz*dz
      dcu(1,mp+1) = dx1
      dcu(2,mp+1) = dy1
      dcu(3,mp+1) = dyp
      dcu(1,mp) = dx2
      dcu(2,mp) = dy2
      dcu(3,mp) = dy
      dx1 = cu(1,mp+1) + ox*dx
      dy1 = cu(2,mp+1) + oy*dx
      dyp = cu(3,mp+1) + oz*dx
      dx2 = cu(1,mp) + ox*dz
      dy2 = cu(2,mp) + oy*dz
      dy = cu(3,mp) + oz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx2
      cu(2,mp) = dy2
      cu(3,mp) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dyp = amu(3,mm+1) + v3*dx
      dx1 = amu(4,mm+1) + v4*dx
      dy1 = amu(1,mm) + v1*dz
      dx2 = amu(2,mm) + v2*dz
      dy2 = amu(3,mm) + v3*dz
      dy = amu(4,mm) + v4*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(3,mm+1) = dyp
      amu(4,mm+1) = dx1
      amu(1,mm) = dy1
      amu(2,mm) = dx2
      amu(3,mm) = dy2
      amu(4,mm) = dy
      dxp = dcu(1,mm+1) + vx*dx
      amx = dcu(2,mm+1) + vy*dx
      dyp = dcu(3,mm+1) + vz*dx
      dx1 = dcu(1,mm) + vx*dz
      dy1 = dcu(2,mm) + vy*dz
      dy = dcu(3,mm) + vz*dz
      dcu(1,mm+1) = dxp
      dcu(2,mm+1) = amx
      dcu(3,mm+1) = dyp
      dcu(1,mm) = dx1
      dcu(2,mm) = dy1
      dcu(3,mm) = dy
      dxp = cu(1,mm+1) + ox*dx
      amx = cu(2,mm+1) + oy*dx
      dyp = cu(3,mm+1) + oz*dx
      dx1 = cu(1,mm) + ox*dz
      dy1 = cu(2,mm) + oy*dz
      dy = cu(3,mm) + oz*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(3,mm+1) = dyp
      cu(1,mm) = dx1
      cu(2,mm) = dy1
      cu(3,mm) = dy
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxn
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxp*dyn
      dy = amx*dyn
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(3,mp+1) = amu(3,mp+1) + v3*dx
      amu(4,mp+1) = amu(4,mp+1) + v4*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      amu(3,mp) = amu(3,mp) + v3*dy
      amu(4,mp) = amu(4,mp) + v4*dy
      dcu(1,mp+1) = dcu(1,mp+1) + vx*dx
      dcu(2,mp+1) = dcu(2,mp+1) + vy*dx
      dcu(3,mp+1) = dcu(3,mp+1) + vz*dx
      dcu(1,mp) = dcu(1,mp) + vx*dy
      dcu(2,mp) = dcu(2,mp) + vy*dy
      dcu(3,mp) = dcu(3,mp) + vz*dy
      cu(1,mp+1) = cu(1,mp+1) + ox*dx
      cu(2,mp+1) = cu(2,mp+1) + oy*dx
      cu(3,mp+1) = cu(3,mp+1) + oz*dx
      cu(1,mp) = cu(1,mp) + ox*dy
      cu(2,mp) = cu(2,mp) + oy*dy
      cu(3,mp) = cu(3,mp) + oz*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(3,mm+1) = amu(3,mm+1) + v3*dx
      amu(4,mm+1) = amu(4,mm+1) + v4*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      amu(3,mm) = amu(3,mm) + v3*dy
      amu(4,mm) = amu(4,mm) + v4*dy
      dcu(1,mm+1) = dcu(1,mm+1) + vx*dx
      dcu(2,mm+1) = dcu(2,mm+1) + vy*dx
      dcu(3,mm+1) = dcu(3,mm+1) + vz*dx
      dcu(1,mm) = dcu(1,mm) + vx*dy
      dcu(2,mm) = dcu(2,mm) + vy*dy
      dcu(3,mm) = dcu(3,mm) + vz*dy
      cu(1,mm+1) = cu(1,mm+1) + ox*dx
      cu(2,mm+1) = cu(2,mm+1) + oy*dx
      cu(3,mm+1) = cu(3,mm+1) + oz*dx
      cu(1,mm) = cu(1,mm) + ox*dy
      cu(2,mm) = cu(2,mm) + oy*dy
      cu(3,mm) = cu(3,mm) + oz*dy
      return
      end
      subroutine GRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle momentum flux
c using first-order spline interpolation for relativistic particles
c scalar version using guard cells
c 40 flops/particle, 1 divide, 12 loads, 8 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy, for i = 1, 2
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c amu(i,j,k) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of flux array, must be >= nx+1
c nyv = second dimension of flux array, must be >= ny+1
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real ci2, gami2, dxp, dyp, amx, amy
      real dx, dy, vx, vy, p2, v1, v2
      ci2 = ci*ci
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami2 = 1.0/(1.0 + p2*ci2)
c calculate weights
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit momentum flux
      dx = dxp*dyp
      dy = amx*dyp
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
   10 continue
      return
      end
      subroutine GSRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle momentum flux
c using first-order linear interpolation for relativistic particles
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 40 flops/particle, 1 divide, 12 loads, 8 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy, for i = 1, 2
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c amu(i,n) = ith component of momentum flux at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
      implicit none
      integer nop, idimp, nxv, nxyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(2,nxyv)
      integer nnn, mmn, nn, mm, mp, j
      real dxn, dyn, ci2, gami2, dxp, dyp, amx, amy, dx, dy, dz, vx, vy
      real p2, v1, v2, dx1, dy1
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      ci2 = ci*ci
c find inverse gamma
      vx = part(3,1)
      vy = part(4,1)
      p2 = vx*vx + vy*vy
      gami2 = 1.0/(1.0 + p2*ci2)
c find interpolation weights
      do 10 j = 2, nop
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j)
      mmn = part(2,j)
      dxp = qm*dxn
      dyp = dyn
      dxn = part(1,j) - float(nnn)
      dyn = part(2,j) - float(mmn)
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyp
c calculate weights
      dx = dxp*dyp
      dz = amx*dyp
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
c get momentum for next particle
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
c deposit momentum flux
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dx = amu(1,mp) + v1*dz
      dy = amu(2,mp) + v2*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(1,mp) = dx
      amu(2,mp) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dx = amu(1,mm) + v1*dz
      dy = amu(2,mm) + v2*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(1,mm) = dx
      amu(2,mm) = dy
c find inverse gamma for next particle
      gami2 = 1.0/(1.0 + p2*ci2)
   10 continue
c deposit momentum flux for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit momentum flux
      dx = dxp*dyn
      dy = amx*dyn
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      return
      end
      subroutine GRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nyv)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density and current density using first-order spline
c interpolation for relativistic particles.
c scalar version using guard cells
c 137 flops/particle, 2 divide, 1 sqrt, 28 loads, 24 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj*gami/dt, where j = x,y, for i = 1, 2
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy, for i = 1, 2
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = -rot(2)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t))*gami, and
c gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j,k) = ith component of current density
c at grid point j,k for i = 1, 2
c dcu(i,j,k) = ith component of acceleration density
c at grid point j,k for i = 1, 2
c amu(i,j,k) = ith component of momentum flux
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxp, dyp, amx, amy, dx, dy
      real ox, oy, oz, acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, p2, v1, v2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyp*(dxp*bz(np,mp) + amx*bz(nn,mp)) + amy*(dxp*bz(np,mm) + am
     1x*bz(nn,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxp*dyp
      dy = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      amu(1,np,mp) = amu(1,np,mp) + v1*dx
      amu(2,np,mp) = amu(2,np,mp) + v2*dx
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dx
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dx
      cu(1,np,mp) = cu(1,np,mp) + ox*dx
      cu(2,np,mp) = cu(2,np,mp) + oy*dx
      dx = dxp*amy
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      dy = amx*amy
      amu(1,np,mm) = amu(1,np,mm) + v1*dx
      amu(2,np,mm) = amu(2,np,mm) + v2*dx
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dx
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dx
      cu(1,np,mm) = cu(1,np,mm) + ox*dx
      cu(2,np,mm) = cu(2,np,mm) + oy*dx
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
   10 continue
      return
      end
      subroutine GSRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nxyv)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using first-order spline
c interpolation for relativistic particles.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 137 flops/particle, 2 divide, 1 sqrt, 28 loads, 24 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c dcu(i,n+1,m)=qci*dx*(1.-dy)
c dcu(i,n,m+1)=qci*(1.-dx)*dy
c dcu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*dvj*gami/dt, where j = x,y, for i = 1, 2
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the grid points
c amu(i,n,m)=qci*(1.-dx)*(1.-dy)
c amu(i,n+1,m)=qci*dx*(1.-dy)
c amu(i,n,m+1)=qci*(1.-dx)*dy
c amu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy, for i = 1, 2
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = -rot(2)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t))*gami, and
c gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c dcu(i,n) = ith component of acceleration density at grid point j,k
c where n = j + nxv*(k-1)
c amu(i,n) = ith component of momentum at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      dimension cu(2,nxyv), dcu(2,nxyv), amu(2,nxyv)
      integer nnn, mmn, nop1, j, nn, mm, mp
      real qtmh, dti, ci2, gami, qtmg, gh, dxn, dyn, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, p2, v1, v2, dx1, dy1, dx2, dy2, dx3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmn
      nnn = part(1,j+1)
      mmn = part(2,j+1)
      dxp = dxn
      dyp = dyn
      dxn = part(1,j+1) - float(nnn)
      dyn = part(2,j+1) - float(mmn)
      mm = mm + nn
      amx = 1. - dxp
      mp = mm + nxv
      amy = 1. - dyp
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyp*(dxp*bz(mp+1) + amx*bz(mp)) + amy*(dxp*bz(mm+1) + amx*bz
     1(mm))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxp
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxp*dyp
      dz = amx*dyp
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      dx1 = amu(1,mp+1) + v1*dx
      dy1 = amu(2,mp+1) + v2*dx
      dy2 = amu(1,mp) + v1*dz
      dx3 = amu(2,mp) + v2*dz
      amu(1,mp+1) = dx1
      amu(2,mp+1) = dy1
      amu(1,mp) = dy2
      amu(2,mp) = dx3
      dx1 = dcu(1,mp+1) + vx*dx
      dy1 = dcu(2,mp+1) + vy*dx
      dx2 = dcu(1,mp) + vx*dz
      dy2 = dcu(2,mp) + vy*dz
      dcu(1,mp+1) = dx1
      dcu(2,mp+1) = dy1
      dcu(1,mp) = dx2
      dcu(2,mp) = dy2
      dx1 = cu(1,mp+1) + ox*dx
      dy1 = cu(2,mp+1) + oy*dx
      dx2 = cu(1,mp) + ox*dz
      dy2 = cu(2,mp) + oy*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(1,mp) = dx2
      cu(2,mp) = dy2
      dx = dxp*amy
      dz = amx*amy
      dxp = amu(1,mm+1) + v1*dx
      amx = amu(2,mm+1) + v2*dx
      dy1 = amu(1,mm) + v1*dz
      dx2 = amu(2,mm) + v2*dz
      amu(1,mm+1) = dxp
      amu(2,mm+1) = amx
      amu(1,mm) = dy1
      amu(2,mm) = dx2
      dxp = dcu(1,mm+1) + vx*dx
      amx = dcu(2,mm+1) + vy*dx
      dx1 = dcu(1,mm) + vx*dz
      dy1 = dcu(2,mm) + vy*dz
      dcu(1,mm+1) = dxp
      dcu(2,mm+1) = amx
      dcu(1,mm) = dx1
      dcu(2,mm) = dy1
      dxp = cu(1,mm+1) + ox*dx
      amx = cu(2,mm+1) + oy*dx
      dx1 = cu(1,mm) + ox*dz
      dy1 = cu(2,mm) + oy*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(1,mm) = dx1
      cu(2,mm) = dy1
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,nop)
      vy = part(4,nop)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyn*(dxn*bz(mp+1) + amx*bz(mp)) + amy*(dxn*bz(mm+1) + amx*bz
     1(mm))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      amx = qm*amx
      dxp = qm*dxn
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxp*dyn
      dy = amx*dyn
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      amu(1,mp+1) = amu(1,mp+1) + v1*dx
      amu(2,mp+1) = amu(2,mp+1) + v2*dx
      amu(1,mp) = amu(1,mp) + v1*dy
      amu(2,mp) = amu(2,mp) + v2*dy
      dcu(1,mp+1) = dcu(1,mp+1) + vx*dx
      dcu(2,mp+1) = dcu(2,mp+1) + vy*dx
      dcu(1,mp) = dcu(1,mp) + vx*dy
      dcu(2,mp) = dcu(2,mp) + vy*dy
      cu(1,mp+1) = cu(1,mp+1) + ox*dx
      cu(2,mp+1) = cu(2,mp+1) + oy*dx
      cu(1,mp) = cu(1,mp) + ox*dy
      cu(2,mp) = cu(2,mp) + oy*dy
      dx = dxp*amy
      dy = amx*amy
      amu(1,mm+1) = amu(1,mm+1) + v1*dx
      amu(2,mm+1) = amu(2,mm+1) + v2*dx
      amu(1,mm) = amu(1,mm) + v1*dy
      amu(2,mm) = amu(2,mm) + v2*dy
      dcu(1,mm+1) = dcu(1,mm+1) + vx*dx
      dcu(2,mm+1) = dcu(2,mm+1) + vy*dx
      dcu(1,mm) = dcu(1,mm) + vx*dy
      dcu(2,mm) = dcu(2,mm) + vy*dy
      cu(1,mm+1) = cu(1,mm+1) + ox*dx
      cu(2,mm+1) = cu(2,mm+1) + oy*dx
      cu(1,mm) = cu(1,mm) + ox*dy
      cu(2,mm) = cu(2,mm) + oy*dy
      return
      end
      subroutine GRMJPOST2C(part,amu,qm,ci,nop,idimp,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using third-order spline interpolation for relativistic particles
c scalar version using guard cells
c 185 flops/particle, 1 divide, 69 loads, 64 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c amu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c amu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c amu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c amu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c amu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c amu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c amu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,j+2,k+2) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nxv = first dimension of flux array, must be >= nx+5
c nyv = second dimension of flux array, must be >= ny+5
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(4,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2, dx, dy, dz, dw, vx, vy, vz, v1, v2, v3, v4
      real ci2, gami2, p2
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
      ci2 = ci*ci
c find interpolation weights
      do 10 j = 1, nop
      nl = part(1,j)
      ml = part(2,j)
      dxp = part(1,j) - float(nl)
      dyp = part(2,j) - float(ml)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = 1.0/(1.0 + p2*ci2)
c calculate weights
      nl = nl + 1
      ml = ml + 1
      dxn = 1.0 - dxp
      dyn = 1.0 - dyp
      dxq = dxp*dxp
      dyq = dyp*dyp
      dxl = dxn*dxn
      dyl = dyn*dyn
      nn = nl + 1
      mm = ml + 1
      at1 = dxl*(1.0 + dxp)
      at2 = dxq*(1.0 + dxn)
      dxq = qms*dxp*dxq
      dxl = qms*dxn*dxl
      dxp = qm23 - qmh*at1
      dxn = qm23 - qmh*at2
      mp = mm + 1
      np = nn + 1
      at1 = dyl*(1.0 + dyp)
      at2 = dyq*(1.0 + dyn)
      dyq = sixth*dyp*dyq
      dyl = sixth*dyn*dyl
      mq = mp + 1
      nq = np + 1
      dyp = two3rds - .5*at1
      dyn = two3rds - .5*at2
c deposit momentum flux
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vz*vx)*gami2
      v4 = (vz*vy)*gami2
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      amu(3,nl,ml) = amu(3,nl,ml) + v3*dx
      amu(4,nl,ml) = amu(4,nl,ml) + v4*dx
      dx = dxl*dyn
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      amu(3,nn,ml) = amu(3,nn,ml) + v3*dy
      amu(4,nn,ml) = amu(4,nn,ml) + v4*dy
      dy = dxn*dyn
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      amu(3,np,ml) = amu(3,np,ml) + v3*dz
      amu(4,np,ml) = amu(4,np,ml) + v4*dz
      dz = dxp*dyn
      amu(1,nq,ml) = amu(1,nq,ml) + v1*dw
      amu(2,nq,ml) = amu(2,nq,ml) + v2*dw
      amu(3,nq,ml) = amu(3,nq,ml) + v3*dw
      amu(4,nq,ml) = amu(4,nq,ml) + v4*dw
      dw = dxq*dyn
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      amu(3,nl,mm) = amu(3,nl,mm) + v3*dx
      amu(4,nl,mm) = amu(4,nl,mm) + v4*dx
      dx = dxl*dyp
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dy = dxn*dyp
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      amu(3,np,mm) = amu(3,np,mm) + v3*dz
      amu(4,np,mm) = amu(4,np,mm) + v4*dz
      dz = dxp*dyp
      amu(1,nq,mm) = amu(1,nq,mm) + v1*dw
      amu(2,nq,mm) = amu(2,nq,mm) + v2*dw
      amu(3,nq,mm) = amu(3,nq,mm) + v3*dw
      amu(4,nq,mm) = amu(4,nq,mm) + v4*dw
      dw = dxq*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      amu(3,nl,mp) = amu(3,nl,mp) + v3*dx
      amu(4,nl,mp) = amu(4,nl,mp) + v4*dx
      dx = dxl*dyq
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dy = dxn*dyq
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      amu(3,np,mp) = amu(3,np,mp) + v3*dz
      amu(4,np,mp) = amu(4,np,mp) + v4*dz
      dz = dxp*dyq
      amu(1,nq,mp) = amu(1,nq,mp) + v1*dw
      amu(2,nq,mp) = amu(2,nq,mp) + v2*dw
      amu(3,nq,mp) = amu(3,nq,mp) + v3*dw
      amu(4,nq,mp) = amu(4,nq,mp) + v4*dw
      dw = dxq*dyq
      amu(1,nl,mq) = amu(1,nl,mq) + v1*dx
      amu(2,nl,mq) = amu(2,nl,mq) + v2*dx
      amu(3,nl,mq) = amu(3,nl,mq) + v3*dx
      amu(4,nl,mq) = amu(4,nl,mq) + v4*dx
      amu(1,nn,mq) = amu(1,nn,mq) + v1*dy
      amu(2,nn,mq) = amu(2,nn,mq) + v2*dy
      amu(3,nn,mq) = amu(3,nn,mq) + v3*dy
      amu(4,nn,mq) = amu(4,nn,mq) + v4*dy
      amu(1,np,mq) = amu(1,np,mq) + v1*dz
      amu(2,np,mq) = amu(2,np,mq) + v2*dz
      amu(3,np,mq) = amu(3,np,mq) + v3*dz
      amu(4,np,mq) = amu(4,np,mq) + v4*dz
      amu(1,nq,mq) = amu(1,nq,mq) + v1*dw
      amu(2,nq,mq) = amu(2,nq,mq) + v2*dw
      amu(3,nq,mq) = amu(3,nq,mq) + v3*dw
      amu(4,nq,mq) = amu(4,nq,mq) + v4*dw
   10 continue
      return
      end
      subroutine GRDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nyv)
c for 2-1/2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using third-order spline
c interpolation for relativistic particles.
c scalar version using guard cells
c 680 flops/particle, 2 divides, 1 sqrt. 261 loads, 160 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c cu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c cu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c cu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c cu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c cu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c cu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c cu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c and qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c dcu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c dcu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c dcu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c dcu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c dcu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c dcu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c dcu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c and qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c amu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c amu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c amu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c amu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c amu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c amu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c amu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami.
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1./6.)((1-dy)**3)*((1./6.)*((1-dx)**3)*fx(n-1,m-1) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m-1) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m-1) + (1./6.)*(dx**3)*fx(n+2,m-1)) + 
c (2./3.-.5*dy*dy*(2-dy))*((1./6.)*((1-dx)**3)*fx(n-1,m) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m) + (1./6.)*(dx**3)*fx(n+2,m)) + 
c (2./3.-.5*((1-dx)**2)*((1./6.)*((1-dx)**3)*fx(n-1,m+1) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m+1) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m+1) + (1./6.)*(dx**3)*fx(n+2,m+1)) + 
c (1./6.)*(dy**3)*((1./6.)*((1-dx)**3)*fx(n-1,m+1) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m+1) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m+1) + (1./6.)*(dx**3)*fx(n+2,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c part(5,n) = momentum pz of particle n at t - dt/2
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c fxy(3,j+2,k+2) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+2,k+2) = x component of magnetic field at grid (j,k)
c bxy(2,j+2,k+2) = y component of magnetic field at grid (j,k)
c bxy(3,j+2,k+2) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j+2,k+2) = ith component of current density
c at grid point j,k for i = 1, 3
c dcu(i,j+2,k+2) = ith component of acceleration density
c at grid point j,k for i = 1, 3
c amu(i,j+2,k+2) = ith component of momentum flux
c at grid point j,k for i = 1, 4
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nxv = first dimension of flux array, must be >= nx+5
c nyv = second dimension of flux array, must be >= ny+5
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, dti, dxl, dyl, dxn, dyn, dxp, dyp, dxq, dyq, at1, at2
      real dx, dy, dz, dw, ox, oy, oz
      real acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real vx, vy, vz, v1, v2, v3, v4, ci2, gami, qtmg, gh, p2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop
c find interpolation weights
      nl = part(1,j)
      ml = part(2,j)
      dxp = part(1,j) - float(nl)
      dyp = part(2,j) - float(ml)
      nl = nl + 1
      ml = ml + 1
      dxn = 1.0 - dxp
      dyn = 1.0 - dyp
      dxq = dxp*dxp
      dyq = dyp*dyp
      dxl = dxn*dxn
      dyl = dyn*dyn
      nn = nl + 1
      mm = ml + 1
      at1 = dxl*(1.0 + dxp)
      at2 = dxq*(1.0 + dxn)
      dxq = sixth*dxp*dxq
      dxl = sixth*dxn*dxl
      dxp = two3rds - .5*at1
      dxn = two3rds - .5*at2
      mp = mm + 1
      np = nn + 1
      at1 = dyl*(1.0 + dyp)
      at2 = dyq*(1.0 + dyn)
      dyq = sixth*dyp*dyq
      dyl = sixth*dyn*dyl
      mq = mp + 1
      nq = np + 1
      dyp = two3rds - .5*at1
      dyn = two3rds - .5*at2
c find electric field
      dx = dyl*(dxl*fxy(1,nl,ml) + dxn*fxy(1,nn,ml) +  dxp*fxy(1,np,ml) 
     1+  dxq*fxy(1,nq,ml)) + dyn*(dxl*fxy(1,nl,mm) + dxn*fxy(1,nn,mm) + 
     2dxp*fxy(1,np,mm) +  dxq*fxy(1,nq,mm)) + dyp*(dxl*fxy(1,nl,mp) + dx
     3n*fxy(1,nn,mp) +  dxp*fxy(1,np,mp) +  dxq*fxy(1,nq,mp)) + dyq*(dxl
     4*fxy(1,nl,mq) + dxn*fxy(1,nn,mq) +  dxp*fxy(1,np,mq) +  dxq*fxy(1,
     5nq,mq))
      dy = dyl*(dxl*fxy(2,nl,ml) + dxn*fxy(2,nn,ml) +  dxp*fxy(2,np,ml) 
     1+  dxq*fxy(2,nq,ml)) + dyn*(dxl*fxy(2,nl,mm) + dxn*fxy(2,nn,mm) + 
     2dxp*fxy(2,np,mm) +  dxq*fxy(2,nq,mm)) + dyp*(dxl*fxy(2,nl,mp) + dx
     3n*fxy(2,nn,mp) +  dxp*fxy(2,np,mp) +  dxq*fxy(2,nq,mp)) + dyq*(dxl
     4*fxy(2,nl,mq) + dxn*fxy(2,nn,mq) +  dxp*fxy(2,np,mq) +  dxq*fxy(2,
     5nq,mq))
      dz = dyl*(dxl*fxy(3,nl,ml) + dxn*fxy(3,nn,ml) +  dxp*fxy(3,np,ml) 
     1+  dxq*fxy(3,nq,ml)) + dyn*(dxl*fxy(3,nl,mm) + dxn*fxy(3,nn,mm) + 
     2dxp*fxy(3,np,mm) +  dxq*fxy(3,nq,mm)) + dyp*(dxl*fxy(3,nl,mp) + dx
     3n*fxy(3,nn,mp) +  dxp*fxy(3,np,mp) +  dxq*fxy(3,nq,mp)) + dyq*(dxl
     4*fxy(3,nl,mq) + dxn*fxy(3,nn,mq) +  dxp*fxy(3,np,mq) +  dxq*fxy(3,
     5nq,mq))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      ox = dyl*(dxl*bxy(1,nl,ml) + dxn*bxy(1,nn,ml) +  dxp*bxy(1,np,ml) 
     1+  dxq*bxy(1,nq,ml)) + dyn*(dxl*bxy(1,nl,mm) + dxn*bxy(1,nn,mm) + 
     2dxp*bxy(1,np,mm) +  dxq*bxy(1,nq,mm)) + dyp*(dxl*bxy(1,nl,mp) + dx
     3n*bxy(1,nn,mp) +  dxp*bxy(1,np,mp) +  dxq*bxy(1,nq,mp)) + dyq*(dxl
     4*bxy(1,nl,mq) + dxn*bxy(1,nn,mq) +  dxp*bxy(1,np,mq) +  dxq*bxy(1,
     5nq,mq))
      oy = dyl*(dxl*bxy(2,nl,ml) + dxn*bxy(2,nn,ml) +  dxp*bxy(2,np,ml) 
     1+  dxq*bxy(2,nq,ml)) + dyn*(dxl*bxy(2,nl,mm) + dxn*bxy(2,nn,mm) + 
     2dxp*bxy(2,np,mm) +  dxq*bxy(2,nq,mm)) + dyp*(dxl*bxy(2,nl,mp) + dx
     3n*bxy(2,nn,mp) +  dxp*bxy(2,np,mp) +  dxq*bxy(2,nq,mp)) + dyq*(dxl
     4*bxy(2,nl,mq) + dxn*bxy(2,nn,mq) +  dxp*bxy(2,np,mq) +  dxq*bxy(2,
     5nq,mq))
      oz = dyl*(dxl*bxy(3,nl,ml) + dxn*bxy(3,nn,ml) +  dxp*bxy(3,np,ml) 
     1+  dxq*bxy(3,nq,ml)) + dyn*(dxl*bxy(3,nl,mm) + dxn*bxy(3,nn,mm) + 
     2dxp*bxy(3,np,mm) +  dxq*bxy(3,nq,mm)) + dyp*(dxl*bxy(3,nl,mp) + dx
     3n*bxy(3,nn,mp) +  dxp*bxy(3,np,mp) +  dxq*bxy(3,nq,mp)) + dyq*(dxl
     4*bxy(3,nl,mq) + dxn*bxy(3,nn,mq) +  dxp*bxy(3,np,mq) +  dxq*bxy(3,
     5nq,mq))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux, acceleration density, and current density
      dxl = qm*dxl
      dxn = qm*dxn
      dxp = qm*dxp
      dxq = qm*dxq
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      v3 = oz*ox
      v4 = oz*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      amu(3,nl,ml) = amu(3,nl,ml) + v3*dx
      amu(4,nl,ml) = amu(4,nl,ml) + v4*dx
      dcu(1,nl,ml) = dcu(1,nl,ml) + vx*dx
      dcu(2,nl,ml) = dcu(2,nl,ml) + vy*dx
      dcu(3,nl,ml) = dcu(3,nl,ml) + vz*dx
      cu(1,nl,ml) = cu(1,nl,ml) + ox*dx
      cu(2,nl,ml) = cu(2,nl,ml) + oy*dx
      cu(3,nl,ml) = cu(3,nl,ml) + oz*dx
      dx = dxl*dyn
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      amu(3,nn,ml) = amu(3,nn,ml) + v3*dy
      amu(4,nn,ml) = amu(4,nn,ml) + v4*dy
      dcu(1,nn,ml) = dcu(1,nn,ml) + vx*dy
      dcu(2,nn,ml) = dcu(2,nn,ml) + vy*dy
      dcu(3,nn,ml) = dcu(3,nn,ml) + vz*dy
      cu(1,nn,ml) = cu(1,nn,ml) + ox*dy
      cu(2,nn,ml) = cu(2,nn,ml) + oy*dy
      cu(3,nn,ml) = cu(3,nn,ml) + oz*dy
      dy = dxn*dyn
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      amu(3,np,ml) = amu(3,np,ml) + v3*dz
      amu(4,np,ml) = amu(4,np,ml) + v4*dz
      dcu(1,np,ml) = dcu(1,np,ml) + vx*dz
      dcu(2,np,ml) = dcu(2,np,ml) + vy*dz
      dcu(3,np,ml) = dcu(3,np,ml) + vz*dz
      cu(1,np,ml) = cu(1,np,ml) + ox*dz
      cu(2,np,ml) = cu(2,np,ml) + oy*dz
      cu(3,np,ml) = cu(3,np,ml) + oz*dz
      dz = dxp*dyn
      amu(1,nq,ml) = amu(1,nq,ml) + v1*dw
      amu(2,nq,ml) = amu(2,nq,ml) + v2*dw
      amu(3,nq,ml) = amu(3,nq,ml) + v3*dw
      amu(4,nq,ml) = amu(4,nq,ml) + v4*dw
      dcu(1,nq,ml) = dcu(1,nq,ml) + vx*dw
      dcu(2,nq,ml) = dcu(2,nq,ml) + vy*dw
      dcu(3,nq,ml) = dcu(3,nq,ml) + vz*dw
      cu(1,nq,ml) = cu(1,nq,ml) + ox*dw
      cu(2,nq,ml) = cu(2,nq,ml) + oy*dw
      cu(3,nq,ml) = cu(3,nq,ml) + oz*dw
      dw = dxq*dyn
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      amu(3,nl,mm) = amu(3,nl,mm) + v3*dx
      amu(4,nl,mm) = amu(4,nl,mm) + v4*dx
      dcu(1,nl,mm) = dcu(1,nl,mm) + vx*dx
      dcu(2,nl,mm) = dcu(2,nl,mm) + vy*dx
      dcu(3,nl,mm) = dcu(3,nl,mm) + vz*dx
      cu(1,nl,mm) = cu(1,nl,mm) + ox*dx
      cu(2,nl,mm) = cu(2,nl,mm) + oy*dx
      cu(3,nl,mm) = cu(3,nl,mm) + oz*dx
      dx = dxl*dyp
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      amu(3,nn,mm) = amu(3,nn,mm) + v3*dy
      amu(4,nn,mm) = amu(4,nn,mm) + v4*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      dcu(3,nn,mm) = dcu(3,nn,mm) + vz*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + oz*dy
      dy = dxn*dyp
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      amu(3,np,mm) = amu(3,np,mm) + v3*dz
      amu(4,np,mm) = amu(4,np,mm) + v4*dz
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dz
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dz
      dcu(3,np,mm) = dcu(3,np,mm) + vz*dz
      cu(1,np,mm) = cu(1,np,mm) + ox*dz
      cu(2,np,mm) = cu(2,np,mm) + oy*dz
      cu(3,np,mm) = cu(3,np,mm) + oz*dz
      dz = dxp*dyp
      amu(1,nq,mm) = amu(1,nq,mm) + v1*dw
      amu(2,nq,mm) = amu(2,nq,mm) + v2*dw
      amu(3,nq,mm) = amu(3,nq,mm) + v3*dw
      amu(4,nq,mm) = amu(4,nq,mm) + v4*dw
      dcu(1,nq,mm) = dcu(1,nq,mm) + vx*dw
      dcu(2,nq,mm) = dcu(2,nq,mm) + vy*dw
      dcu(3,nq,mm) = dcu(3,nq,mm) + vz*dw
      cu(1,nq,mm) = cu(1,nq,mm) + ox*dw
      cu(2,nq,mm) = cu(2,nq,mm) + oy*dw
      cu(3,nq,mm) = cu(3,nq,mm) + oz*dw
      dw = dxq*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      amu(3,nl,mp) = amu(3,nl,mp) + v3*dx
      amu(4,nl,mp) = amu(4,nl,mp) + v4*dx
      dcu(1,nl,mp) = dcu(1,nl,mp) + vx*dx
      dcu(2,nl,mp) = dcu(2,nl,mp) + vy*dx
      dcu(3,nl,mp) = dcu(3,nl,mp) + vz*dx
      cu(1,nl,mp) = cu(1,nl,mp) + ox*dx
      cu(2,nl,mp) = cu(2,nl,mp) + oy*dx
      cu(3,nl,mp) = cu(3,nl,mp) + oz*dx
      dx = dxl*dyq
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      amu(3,nn,mp) = amu(3,nn,mp) + v3*dy
      amu(4,nn,mp) = amu(4,nn,mp) + v4*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      dcu(3,nn,mp) = dcu(3,nn,mp) + vz*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + oz*dy
      dy = dxn*dyq
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      amu(3,np,mp) = amu(3,np,mp) + v3*dz
      amu(4,np,mp) = amu(4,np,mp) + v4*dz
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dz
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dz
      dcu(3,np,mp) = dcu(3,np,mp) + vz*dz
      cu(1,np,mp) = cu(1,np,mp) + ox*dz
      cu(2,np,mp) = cu(2,np,mp) + oy*dz
      cu(3,np,mp) = cu(3,np,mp) + oz*dz
      dz = dxp*dyq
      amu(1,nq,mp) = amu(1,nq,mp) + v1*dw
      amu(2,nq,mp) = amu(2,nq,mp) + v2*dw
      amu(3,nq,mp) = amu(3,nq,mp) + v3*dw
      amu(4,nq,mp) = amu(4,nq,mp) + v4*dw
      dcu(1,nq,mp) = dcu(1,nq,mp) + vx*dw
      dcu(2,nq,mp) = dcu(2,nq,mp) + vy*dw
      dcu(3,nq,mp) = dcu(3,nq,mp) + vz*dw
      cu(1,nq,mp) = cu(1,nq,mp) + ox*dw
      cu(2,nq,mp) = cu(2,nq,mp) + oy*dw
      cu(3,nq,mp) = cu(3,nq,mp) + oz*dw
      dw = dxq*dyq
      amu(1,nl,mq) = amu(1,nl,mq) + v1*dx
      amu(2,nl,mq) = amu(2,nl,mq) + v2*dx
      amu(3,nl,mq) = amu(3,nl,mq) + v3*dx
      amu(4,nl,mq) = amu(4,nl,mq) + v4*dx
      dcu(1,nl,mq) = dcu(1,nl,mq) + vx*dx
      dcu(2,nl,mq) = dcu(2,nl,mq) + vy*dx
      dcu(3,nl,mq) = dcu(3,nl,mq) + vz*dx
      cu(1,nl,mq) = cu(1,nl,mq) + ox*dx
      cu(2,nl,mq) = cu(2,nl,mq) + oy*dx
      cu(3,nl,mq) = cu(3,nl,mq) + oz*dx
      amu(1,nn,mq) = amu(1,nn,mq) + v1*dy
      amu(2,nn,mq) = amu(2,nn,mq) + v2*dy
      amu(3,nn,mq) = amu(3,nn,mq) + v3*dy
      amu(4,nn,mq) = amu(4,nn,mq) + v4*dy
      dcu(1,nn,mq) = dcu(1,nn,mq) + vx*dy
      dcu(2,nn,mq) = dcu(2,nn,mq) + vy*dy
      dcu(3,nn,mq) = dcu(3,nn,mq) + vz*dy
      cu(1,nn,mq) = cu(1,nn,mq) + ox*dy
      cu(2,nn,mq) = cu(2,nn,mq) + oy*dy
      cu(3,nn,mq) = cu(3,nn,mq) + oz*dy
      amu(1,np,mq) = amu(1,np,mq) + v1*dz
      amu(2,np,mq) = amu(2,np,mq) + v2*dz
      amu(3,np,mq) = amu(3,np,mq) + v3*dz
      amu(4,np,mq) = amu(4,np,mq) + v4*dz
      dcu(1,np,mq) = dcu(1,np,mq) + vx*dz
      dcu(2,np,mq) = dcu(2,np,mq) + vy*dz
      dcu(3,np,mq) = dcu(3,np,mq) + vz*dz
      cu(1,np,mq) = cu(1,np,mq) + ox*dz
      cu(2,np,mq) = cu(2,np,mq) + oy*dz
      cu(3,np,mq) = cu(3,np,mq) + oz*dz
      amu(1,nq,mq) = amu(1,nq,mq) + v1*dw
      amu(2,nq,mq) = amu(2,nq,mq) + v2*dw
      amu(3,nq,mq) = amu(3,nq,mq) + v3*dw
      amu(4,nq,mq) = amu(4,nq,mq) + v4*dw
      dcu(1,nq,mq) = dcu(1,nq,mq) + vx*dw
      dcu(2,nq,mq) = dcu(2,nq,mq) + vy*dw
      dcu(3,nq,mq) = dcu(3,nq,mq) + vz*dw
      cu(1,nq,mq) = cu(1,nq,mq) + ox*dw
      cu(2,nq,mq) = cu(2,nq,mq) + oy*dw
      cu(3,nq,mq) = cu(3,nq,mq) + oz*dw
   10 continue
      return
      end
      subroutine GRMJPOST22C(part,amu,qm,ci,nop,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle momentum flux
c using third-order spline interpolation for relativistic particles
c scalar version using guard cells
c 115 flops/particle, 1 divide, 36 loads, 32 stores
c input: all, output: amu
c momentum flux is approximated by values at the nearest grid points
c amu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c amu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c amu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c amu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c amu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c amu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c amu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c amu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy, for i = 1, 2
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c amu(i,j+2,k+2) = ith component of momentum flux at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of flux array, must be >= nx+5
c nyv = second dimension of flux array, must be >= ny+5
      implicit none
      integer nop, idimp, nxv, nyv
      real part, amu, qm, ci
      dimension part(idimp,nop), amu(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2, dx, dy, dz, dw, vx, vy, v1, v2, ci2, gami2, p2
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
      ci2 = ci*ci
c find interpolation weights
      do 10 j = 1, nop
      nl = part(1,j)
      ml = part(2,j)
      dxp = part(1,j) - float(nl)
      dyp = part(2,j) - float(ml)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami2 = 1.0/(1.0 + p2*ci2)
c calculate weights
      nl = nl + 1
      ml = ml + 1
      dxn = 1.0 - dxp
      dyn = 1.0 - dyp
      dxq = dxp*dxp
      dyq = dyp*dyp
      dxl = dxn*dxn
      dyl = dyn*dyn
      nn = nl + 1
      mm = ml + 1
      at1 = dxl*(1.0 + dxp)
      at2 = dxq*(1.0 + dxn)
      dxq = qms*dxp*dxq
      dxl = qms*dxn*dxl
      dxp = qm23 - qmh*at1
      dxn = qm23 - qmh*at2
      mp = mm + 1
      np = nn + 1
      at1 = dyl*(1.0 + dyp)
      at2 = dyq*(1.0 + dyn)
      dyq = sixth*dyp*dyq
      dyl = sixth*dyn*dyl
      mq = mp + 1
      nq = np + 1
      dyp = two3rds - .5*at1
      dyn = two3rds - .5*at2
c deposit momentum flux
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      dx = dxl*dyn
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      dy = dxn*dyn
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      dz = dxp*dyn
      amu(1,nq,ml) = amu(1,nq,ml) + v1*dw
      amu(2,nq,ml) = amu(2,nq,ml) + v2*dw
      dw = dxq*dyn
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      dx = dxl*dyp
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      dy = dxn*dyp
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      dz = dxp*dyp
      amu(1,nq,mm) = amu(1,nq,mm) + v1*dw
      amu(2,nq,mm) = amu(2,nq,mm) + v2*dw
      dw = dxq*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      dx = dxl*dyq
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      dy = dxn*dyq
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      dz = dxp*dyq
      amu(1,nq,mp) = amu(1,nq,mp) + v1*dw
      amu(2,nq,mp) = amu(2,nq,mp) + v2*dw
      dw = dxq*dyq
      amu(1,nl,mq) = amu(1,nl,mq) + v1*dx
      amu(2,nl,mq) = amu(2,nl,mq) + v2*dx
      amu(1,nn,mq) = amu(1,nn,mq) + v1*dy
      amu(2,nn,mq) = amu(2,nn,mq) + v2*dy
      amu(1,np,mq) = amu(1,np,mq) + v1*dz
      amu(2,np,mq) = amu(2,np,mq) + v2*dz
      amu(1,nq,mq) = amu(1,nq,mq) + v1*dw
      amu(2,nq,mq) = amu(2,nq,mq) + v2*dw
   10 continue
      return
      end
      subroutine GRDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nyv)
c for 2d code, this subroutine calculates particle momentum flux,
c acceleration density, and current density using third-order spline
c interpolation for relativistic particles.
c scalar version using guard cells
c 401 flops/particle, 2 divides, 1sqrt, 148 loads, 96 stores
c input: all, output: cu, dcu, amu
c current density is approximated by values at the nearest grid points
c cu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c cu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c cu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c cu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c cu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c cu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c cu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c cu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c cu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c cu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is approximated by values at the nearest grid
c points
c dcu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c dcu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c dcu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c dcu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c dcu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c dcu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c dcu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c dcu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c dcu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c dcu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c and qci = qm*dvj*gami/dt, where j = x,y, for i = 1, 2
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is approximated by values at the nearest grid points
c amu(i,n-1,m-1)=(qci/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c amu(i,n,m-1)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c amu(i,n+1,m-1)=qci*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c amu(i,n+2,m-1)=(qci/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c amu(i,n-1,m)=(qci/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n,m)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+1,m)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n+2,m)=(qci/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c amu(i,n-1,m+1)=(qci/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n,m+1)=qci*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+1,m+1)=qci*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n+2,m+1)=(qci/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c amu(i,n-1,m+2)=(qci/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c amu(i,n,m+2)=qci*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c amu(i,n+1,m+2)=qci*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c amu(i,n+2,m+2)=(qci/6.)*(dx**3)*(1./6.)*(dy**3)
c and qci = qm*pj*pk*gami**2, where jk = xx-yy,xy, for i = 1, 2
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = -rot(2)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t))*gami, and
c gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1./6.)((1-dy)**3)*((1./6.)*((1-dx)**3)*fx(n-1,m-1) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m-1) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m-1) + (1./6.)*(dx**3)*fx(n+2,m-1)) + 
c (2./3.-.5*dy*dy*(2-dy))*((1./6.)*((1-dx)**3)*fx(n-1,m) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m) + (1./6.)*(dx**3)*fx(n+2,m)) + 
c (2./3.-.5*((1-dx)**2)*((1./6.)*((1-dx)**3)*fx(n-1,m+1) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m+1) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m+1) + (1./6.)*(dx**3)*fx(n+2,m+1)) + 
c (1./6.)*(dy**3)*((1./6.)*((1-dx)**3)*fx(n-1,m+1) +
c (2./3.-.5*dx*dx*(2-dx))*fx(n,m+1) + (2./3.-.5*((1-dx)**2)*(1+dx))*
c fx(n+1,m+1) + (1./6.)*(dx**3)*fx(n+2,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+2,k+2) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c cu(i,j+2,k+2) = ith component of current density
c at grid point j,k for i = 1, 2
c dcu(i,j+2,k+2) = ith component of acceleration density
c at grid point j,k for i = 1, 2
c amu(i,j+2,k+2) = ith component of momentum flux
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of flux array, must be >= nx+5
c nyv = second dimension of flux array, must be >= ny+5
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, dti, dxl, dyl, dxn, dyn, dxp, dyp, dxq, dyq, at1, at2
      real dx, dy, dz, dw, ox, oy, oz
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real vx, vy, v1, v2, ci2, gami, qtmg, gh, p2
      qtmh = .5*qbm*dt
      dti = 1.0/dt
      ci2 = ci*ci
      do 10 j = 1, nop
c find interpolation weights
      nl = part(1,j)
      ml = part(2,j)
      dxp = part(1,j) - float(nl)
      dyp = part(2,j) - float(ml)
      nl = nl + 1
      ml = ml + 1
      dxn = 1.0 - dxp
      dyn = 1.0 - dyp
      dxq = dxp*dxp
      dyq = dyp*dyp
      dxl = dxn*dxn
      dyl = dyn*dyn
      nn = nl + 1
      mm = ml + 1
      at1 = dxl*(1.0 + dxp)
      at2 = dxq*(1.0 + dxn)
      dxq = sixth*dxp*dxq
      dxl = sixth*dxn*dxl
      dxp = two3rds - .5*at1
      dxn = two3rds - .5*at2
      mp = mm + 1
      np = nn + 1
      at1 = dyl*(1.0 + dyp)
      at2 = dyq*(1.0 + dyn)
      dyq = sixth*dyp*dyq
      dyl = sixth*dyn*dyl
      mq = mp + 1
      nq = np + 1
      dyp = two3rds - .5*at1
      dyn = two3rds - .5*at2
c find electric field
      dx = dyl*(dxl*fxy(1,nl,ml) + dxn*fxy(1,nn,ml) +  dxp*fxy(1,np,ml) 
     1+  dxq*fxy(1,nq,ml)) + dyn*(dxl*fxy(1,nl,mm) + dxn*fxy(1,nn,mm) + 
     2dxp*fxy(1,np,mm) +  dxq*fxy(1,nq,mm)) + dyp*(dxl*fxy(1,nl,mp) + dx
     3n*fxy(1,nn,mp) +  dxp*fxy(1,np,mp) +  dxq*fxy(1,nq,mp)) + dyq*(dxl
     4*fxy(1,nl,mq) + dxn*fxy(1,nn,mq) +  dxp*fxy(1,np,mq) +  dxq*fxy(1,
     5nq,mq))
      dy = dyl*(dxl*fxy(2,nl,ml) + dxn*fxy(2,nn,ml) +  dxp*fxy(2,np,ml) 
     1+  dxq*fxy(2,nq,ml)) + dyn*(dxl*fxy(2,nl,mm) + dxn*fxy(2,nn,mm) + 
     2dxp*fxy(2,np,mm) +  dxq*fxy(2,nq,mm)) + dyp*(dxl*fxy(2,nl,mp) + dx
     3n*fxy(2,nn,mp) +  dxp*fxy(2,np,mp) +  dxq*fxy(2,nq,mp)) + dyq*(dxl
     4*fxy(2,nl,mq) + dxn*fxy(2,nn,mq) +  dxp*fxy(2,np,mq) +  dxq*fxy(2,
     5nq,mq))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      vx = part(3,j)
      vy = part(4,j)
      acx = vx + dx
      acy = vy + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyl*(dxl*bz(nl,ml) + dxn*bz(nn,ml) +  dxp*bz(np,ml) + dxq*bz(
     1nq,ml)) + dyn*(dxl*bz(nl,mm) + dxn*bz(nn,mm) + dxp*bz(np,mm) +  dx
     2q*bz(nq,mm)) + dyp*(dxl*bz(nl,mp) + dxn*bz(nn,mp) +  dxp*bz(np,mp)
     3+  dxq*bz(nq,mp)) + dyq*(dxl*bz(nl,mq) + dxn*bz(nn,mq) +  dxp*bz(n
     4p,mq) +  dxq*bz(nq,mq))
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      v1 = (rot1*acx + rot2*acy)*anorm + dx
      v2 = (rot1*acy - rot2*acx)*anorm + dy
c deposit momentum flux, acceleration density, and current density
      dxl = qm*dxl
      dxn = qm*dxn
      dxp = qm*dxp
      dxq = qm*dxq
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      vx = v1 - vx
      vy = v2 - vy
      gh = 2.0*ci2*(ox*dx + oy*dy)
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      v1 = ox*ox - oy*oy
      v2 = ox*oy
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      amu(1,nl,ml) = amu(1,nl,ml) + v1*dx
      amu(2,nl,ml) = amu(2,nl,ml) + v2*dx
      dcu(1,nl,ml) = dcu(1,nl,ml) + vx*dx
      dcu(2,nl,ml) = dcu(2,nl,ml) + vy*dx
      cu(1,nl,ml) = cu(1,nl,ml) + ox*dx
      cu(2,nl,ml) = cu(2,nl,ml) + oy*dx
      dx = dxl*dyn
      amu(1,nn,ml) = amu(1,nn,ml) + v1*dy
      amu(2,nn,ml) = amu(2,nn,ml) + v2*dy
      dcu(1,nn,ml) = dcu(1,nn,ml) + vx*dy
      dcu(2,nn,ml) = dcu(2,nn,ml) + vy*dy
      cu(1,nn,ml) = cu(1,nn,ml) + ox*dy
      cu(2,nn,ml) = cu(2,nn,ml) + oy*dy
      dy = dxn*dyn
      amu(1,np,ml) = amu(1,np,ml) + v1*dz
      amu(2,np,ml) = amu(2,np,ml) + v2*dz
      dcu(1,np,ml) = dcu(1,np,ml) + vx*dz
      dcu(2,np,ml) = dcu(2,np,ml) + vy*dz
      cu(1,np,ml) = cu(1,np,ml) + ox*dz
      cu(2,np,ml) = cu(2,np,ml) + oy*dz
      dz = dxp*dyn
      amu(1,nq,ml) = amu(1,nq,ml) + v1*dw
      amu(2,nq,ml) = amu(2,nq,ml) + v2*dw
      dcu(1,nq,ml) = dcu(1,nq,ml) + vx*dw
      dcu(2,nq,ml) = dcu(2,nq,ml) + vy*dw
      cu(1,nq,ml) = cu(1,nq,ml) + ox*dw
      cu(2,nq,ml) = cu(2,nq,ml) + oy*dw
      dw = dxq*dyn
      amu(1,nl,mm) = amu(1,nl,mm) + v1*dx
      amu(2,nl,mm) = amu(2,nl,mm) + v2*dx
      dcu(1,nl,mm) = dcu(1,nl,mm) + vx*dx
      dcu(2,nl,mm) = dcu(2,nl,mm) + vy*dx
      cu(1,nl,mm) = cu(1,nl,mm) + ox*dx
      cu(2,nl,mm) = cu(2,nl,mm) + oy*dx
      dx = dxl*dyp
      amu(1,nn,mm) = amu(1,nn,mm) + v1*dy
      amu(2,nn,mm) = amu(2,nn,mm) + v2*dy
      dcu(1,nn,mm) = dcu(1,nn,mm) + vx*dy
      dcu(2,nn,mm) = dcu(2,nn,mm) + vy*dy
      cu(1,nn,mm) = cu(1,nn,mm) + ox*dy
      cu(2,nn,mm) = cu(2,nn,mm) + oy*dy
      dy = dxn*dyp
      amu(1,np,mm) = amu(1,np,mm) + v1*dz
      amu(2,np,mm) = amu(2,np,mm) + v2*dz
      dcu(1,np,mm) = dcu(1,np,mm) + vx*dz
      dcu(2,np,mm) = dcu(2,np,mm) + vy*dz
      cu(1,np,mm) = cu(1,np,mm) + ox*dz
      cu(2,np,mm) = cu(2,np,mm) + oy*dz
      dz = dxp*dyp
      amu(1,nq,mm) = amu(1,nq,mm) + v1*dw
      amu(2,nq,mm) = amu(2,nq,mm) + v2*dw
      dcu(1,nq,mm) = dcu(1,nq,mm) + vx*dw
      dcu(2,nq,mm) = dcu(2,nq,mm) + vy*dw
      cu(1,nq,mm) = cu(1,nq,mm) + ox*dw
      cu(2,nq,mm) = cu(2,nq,mm) + oy*dw
      dw = dxq*dyp
      amu(1,nl,mp) = amu(1,nl,mp) + v1*dx
      amu(2,nl,mp) = amu(2,nl,mp) + v2*dx
      dcu(1,nl,mp) = dcu(1,nl,mp) + vx*dx
      dcu(2,nl,mp) = dcu(2,nl,mp) + vy*dx
      cu(1,nl,mp) = cu(1,nl,mp) + ox*dx
      cu(2,nl,mp) = cu(2,nl,mp) + oy*dx
      dx = dxl*dyq
      amu(1,nn,mp) = amu(1,nn,mp) + v1*dy
      amu(2,nn,mp) = amu(2,nn,mp) + v2*dy
      dcu(1,nn,mp) = dcu(1,nn,mp) + vx*dy
      dcu(2,nn,mp) = dcu(2,nn,mp) + vy*dy
      cu(1,nn,mp) = cu(1,nn,mp) + ox*dy
      cu(2,nn,mp) = cu(2,nn,mp) + oy*dy
      dy = dxn*dyq
      amu(1,np,mp) = amu(1,np,mp) + v1*dz
      amu(2,np,mp) = amu(2,np,mp) + v2*dz
      dcu(1,np,mp) = dcu(1,np,mp) + vx*dz
      dcu(2,np,mp) = dcu(2,np,mp) + vy*dz
      cu(1,np,mp) = cu(1,np,mp) + ox*dz
      cu(2,np,mp) = cu(2,np,mp) + oy*dz
      dz = dxp*dyq
      amu(1,nq,mp) = amu(1,nq,mp) + v1*dw
      amu(2,nq,mp) = amu(2,nq,mp) + v2*dw
      dcu(1,nq,mp) = dcu(1,nq,mp) + vx*dw
      dcu(2,nq,mp) = dcu(2,nq,mp) + vy*dw
      cu(1,nq,mp) = cu(1,nq,mp) + ox*dw
      cu(2,nq,mp) = cu(2,nq,mp) + oy*dw
      dw = dxq*dyq
      amu(1,nl,mq) = amu(1,nl,mq) + v1*dx
      amu(2,nl,mq) = amu(2,nl,mq) + v2*dx
      dcu(1,nl,mq) = dcu(1,nl,mq) + vx*dx
      dcu(2,nl,mq) = dcu(2,nl,mq) + vy*dx
      cu(1,nl,mq) = cu(1,nl,mq) + ox*dx
      cu(2,nl,mq) = cu(2,nl,mq) + oy*dx
      amu(1,nn,mq) = amu(1,nn,mq) + v1*dy
      amu(2,nn,mq) = amu(2,nn,mq) + v2*dy
      dcu(1,nn,mq) = dcu(1,nn,mq) + vx*dy
      dcu(2,nn,mq) = dcu(2,nn,mq) + vy*dy
      cu(1,nn,mq) = cu(1,nn,mq) + ox*dy
      cu(2,nn,mq) = cu(2,nn,mq) + oy*dy
      amu(1,np,mq) = amu(1,np,mq) + v1*dz
      amu(2,np,mq) = amu(2,np,mq) + v2*dz
      dcu(1,np,mq) = dcu(1,np,mq) + vx*dz
      dcu(2,np,mq) = dcu(2,np,mq) + vy*dz
      cu(1,np,mq) = cu(1,np,mq) + ox*dz
      cu(2,np,mq) = cu(2,np,mq) + oy*dz
      amu(1,nq,mq) = amu(1,nq,mq) + v1*dw
      amu(2,nq,mq) = amu(2,nq,mq) + v2*dw
      dcu(1,nq,mq) = dcu(1,nq,mq) + vx*dw
      dcu(2,nq,mq) = dcu(2,nq,mq) + vy*dw
      cu(1,nq,mq) = cu(1,nq,mq) + ox*dw
      cu(2,nq,mq) = cu(2,nq,mq) + oy*dw
   10 continue
      return
      end
      subroutine GRMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxvh,ny
     1v)
c for 2-1/2d code, this subroutine calculates particle momentum flux
c using gridless spectral version for relativistic particles,
c and periodic boundaries
c baseline scalar version
c 44*(NX/2)*(NY/2) + 10 flops/particle, 2*NX*NY loads and stores
c plus (NX/2 + NY/2) sines and cosines/particle
c input: all, output: amu
c momentum flux is calculated from the expression:
c amu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c and qci = qm*pj*pk*gami2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = pj(t-dt/2) and pk = pk(t-dt/2)
c where gami2 = 1./(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x momentum of particle n at t - dt/2
c part(4,n) = y momentum of particle n at t - dt/2
c part(5,n) = z momentum of particle n at t - dt/2
c amu(i,j,k) = ith component of momentum flux
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nx/2
c nyv = second dimension of current array, must be >= ny
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv
      real part, qm, ci
      complex amu, sctx
      dimension part(idimp,nop), amu(4,nxvh,nyv), sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, ci2, gami2, dnx, dny, dkx, dky, at1, at3
      real vx, vy, vz, p2, v1, v2, v3, v4
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qmn = qm/real(nx*ny)
      ci2 = ci*ci
c find fourier components
      do 50 i = 1, nop
      vx = part(3,i)
      vy = part(4,i)
      vz = part(5,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami2 = qmn/(1.0 + p2*ci2)
      v1 = (vx*vx - vy*vy)*gami2
      v2 = (vx*vy)*gami2
      v3 = (vz*vx)*gami2
      v4 = (vz*vy)*gami2
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      amu(1,j,k) = amu(1,j,k) + v1*zt3
      amu(2,j,k) = amu(2,j,k) + v2*zt3
      amu(3,j,k) = amu(3,j,k) + v3*zt3
      amu(4,j,k) = amu(4,j,k) + v4*zt3
      amu(1,j,k1) = amu(1,j,k1) + v1*zt4
      amu(2,j,k1) = amu(2,j,k1) + v2*zt4
      amu(3,j,k1) = amu(3,j,k1) + v3*zt4
      amu(4,j,k1) = amu(4,j,k1) + v4*zt4
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = at3*zt2
      amu(1,1,k) = amu(1,1,k) + v1*zt1
      amu(2,1,k) = amu(2,1,k) + v2*zt1
      amu(3,1,k) = amu(3,1,k) + v3*zt1
      amu(4,1,k) = amu(4,1,k) + v4*zt1
      amu(1,1,k1) = amu(1,1,k1) + v1*zt2
      amu(2,1,k1) = amu(2,1,k1) + v2*zt2
      amu(3,1,k1) = amu(3,1,k1) + v3*zt2
      amu(4,1,k1) = amu(4,1,k1) + v4*zt2
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      amu(1,j,1) = amu(1,j,1) + v1*zt3
      amu(2,j,1) = amu(2,j,1) + v2*zt3
      amu(3,j,1) = amu(3,j,1) + v3*zt3
      amu(4,j,1) = amu(4,j,1) + v4*zt3
      amu(1,j,k1) = amu(1,j,k1) + v1*zt1
      amu(2,j,k1) = amu(2,j,k1) + v2*zt1
      amu(3,j,k1) = amu(3,j,k1) + v3*zt1
      amu(4,j,k1) = amu(4,j,k1) + v4*zt1
   40 continue
      at3 = real(sctx(nxh))
      amu(1,1,1) = cmplx(real(amu(1,1,1))+v1,aimag(amu(1,1,1))+v1*at3)
      amu(2,1,1) = cmplx(real(amu(2,1,1))+v2,aimag(amu(2,1,1))+v2*at3)
      amu(3,1,1) = cmplx(real(amu(3,1,1))+v3,aimag(amu(3,1,1))+v3*at3)
      amu(4,1,1) = cmplx(real(amu(4,1,1))+v4,aimag(amu(4,1,1))+v4*at3)
      at3 = at1*at3
      amu(1,1,k1) = cmplx(real(amu(1,1,k1))+v1*at1,aimag(amu(1,1,k1))+v1
     1*at3)
      amu(2,1,k1) = cmplx(real(amu(2,1,k1))+v2*at1,aimag(amu(2,1,k1))+v2
     1*at3)
      amu(3,1,k1) = cmplx(real(amu(3,1,k1))+v3*at1,aimag(amu(3,1,k1))+v3
     1*at3)
      amu(4,1,k1) = cmplx(real(amu(4,1,k1))+v4*at1,aimag(amu(4,1,k1))+v4
     1*at3)
   50 continue
      return
      end
      subroutine GRDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,ci,
     1idimp,nop,nx,ny,nxvh,nyv)
c for 2-1/2d code, this subroutine updates particle momentum flux,
c acceleration density and current density
c using gridless spectral version for relativistic particles,
c with periodic boundaries,
c baseline scalar version
c 132*(NX/2)*(NY/2) + 105 flops/particle,
c 19*NX*NY/2 loads and 13*NX*NY/2 stores
c input: all, output: cu, dcu, amu
c current density is calculated from the expression:
c cu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c where qci = qm*pj*gami, where j = x,y,z, for i = 1, 3
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c acceleration density is calculated from the expression:
c dcu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c where qci = qm*dvj*gami/dt, where j = x,y,z, for i = 1, 3
c where dvj = dpj - pj*gami*dgamma, dpj = (pj(t+dt/2)-pj(t-dt/2)), 
c pj = .5*(pj(t+dt/2)+pj(t-dt/2)),
c dgamma = (q/m)*ci*ci*gami*(sum(pj*Ej))*dt,
c and Ej = jth component of electric field
c momentum flux is calculated from the expression:
c amu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c where qci = qm*vj*vk, where jk = xx-y*y,xy,zx,zy, for i = 1, 4
c where qci = qm*pj*pk*gami**2, where jk = xx-yy,xy,zx,zy, for i = 1, 4
c where pj = 0.5*(pj(t+dt/2)+pj(t-dt/2),
c and pk = 0.5*(pk(t+dt/2)+pk(t-dt/2))
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
c omz = (q/m)*bz(x(t),y(t))*gami.
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,j,k) = i component of magnetic field at fourier grid point n,m
c cu(i,n,m) = current density at fourier grid point n,m
c dcu(i,n,m) = acceleration density at fourier grid point n,m
c amu(i,n,m) = momentum flux at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv
      real part, qm, qbm, dt, ci
      complex fxy, bxy, cu, dcu, amu, sctx
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension cu(3,nxvh,nyv), dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3, qtmh, qmn, dti
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9, anorm
      real ci2, gami, qtmg, gh, vx, vy, vz, p2, v1, v2, v3, v4
      double precision exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qtmh = 0.5*qbm*dt
      qmn = qm/real(nx*ny)
      dti = qmn/dt
      ci2 = ci*ci
      do 90 i = 1, nop
c find electric and magnetic field
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      bx = 0.0d0
      by = 0.0d0
      bz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      bxl = bxl + (real(bxy(1,j,k)*zt3) + real(bxy(1,j,k1)*zt4))
      byl = byl + (real(bxy(2,j,k)*zt3) + real(bxy(2,j,k1)*zt4))
      bzl = bzl + (real(bxy(3,j,k)*zt3) + real(bxy(3,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
      eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
      bxl = bxl + (real(bxy(1,1,k)*zt1) + real(bxy(1,1,k1)*zt2))
      byl = byl + (real(bxy(2,1,k)*zt1) + real(bxy(2,1,k1)*zt2))
      bzl = bzl + (real(bxy(3,1,k)*zt1) + real(bxy(3,1,k1)*zt2))
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      bx = bx + bxl
      by = by + byl
      bz = bz + bzl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      bxl = 0.0d0
      byl = 0.0d0
      bzl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
      eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
      ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
      bxl = bxl + (real(bxy(1,j,1)*zt3) + real(bxy(1,j,k1)*zt1))
      byl = byl + (real(bxy(2,j,1)*zt3) + real(bxy(2,j,k1)*zt1))
      bzl = bzl + (real(bxy(3,j,1)*zt3) + real(bxy(3,j,k1)*zt1))
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      bx = 2.0d0*(bx + bxl)
      by = 2.0d0*(by + byl)
      bz = 2.0d0*(bz + bzl)
      at3 = real(sctx(nxh))
      ex = ex + (real(fxy(1,1,1)) + real(fxy(1,1,k1))*at1)
      ey = ey + (real(fxy(2,1,1)) + real(fxy(2,1,k1))*at1)
      ez = ez + (real(fxy(3,1,1)) + real(fxy(3,1,k1))*at1)
      bx = bx + (real(bxy(1,1,1)) + real(bxy(1,1,k1))*at1)
      by = by + (real(bxy(2,1,1)) + real(bxy(2,1,k1))*at1)
      bz = bz + (real(bxy(3,1,1)) + real(bxy(3,1,k1))*at1)
      at1 = at1*at3
      dx = ex + (aimag(fxy(1,1,1))*at3 + aimag(fxy(1,1,k1))*at1)
      dy = ey + (aimag(fxy(2,1,1))*at3 + aimag(fxy(2,1,k1))*at1)
      dz = ez + (aimag(fxy(3,1,1))*at3 + aimag(fxy(3,1,k1))*at1)
      ox = bx + (aimag(bxy(1,1,1))*at3 + aimag(bxy(1,1,k1))*at1)
      oy = by + (aimag(bxy(2,1,1))*at3 + aimag(bxy(2,1,k1))*at1)
      oz = bz + (aimag(bxy(3,1,1))*at3 + aimag(bxy(3,1,k1))*at1)
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      vx = part(3,i)
      vy = part(4,i)
      vz = part(5,i)
      acx = vx + dx
      acy = vy + dy
      acz = vz + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
      gh = 0.5*gami
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
      qtmg = dti*gami
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = rot4 + omzt
      rot4 = rot4 - omzt
      rot3 = rot7 - omyt
      rot7 = rot7 + omyt
      rot6 = rot8 + omxt
      rot8 = rot8 - omxt
c new velocity
      v1 = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      v2 = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      v3 = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
c deposit momentum flux and acceleration density
      ox = gh*(v1 + vx)
      oy = gh*(v2 + vy)
      oz = gh*(v3 + vz)
      gh = 2.0*ci2*(ox*dx + oy*dy + oz*dz)
      vx = v1 - vx
      vy = v2 - vy
      vz = v3 - vz
      vx = qtmg*(vx - ox*gh)
      vy = qtmg*(vy - oy*gh)
      vz = qtmg*(vz - oz*gh)
      v1 = qmn*(ox*ox - oy*oy)
      ox = qmn*ox
      v2 = ox*oy
      oy = qmn*oy
      v3 = oz*ox
      v4 = oz*oy
      oz = qmn*oz
      do 50 j = 1, nxh
      sctx(j) = conjg(sctx(j))
   50 continue
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 70 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 60 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      amu(1,j,k) = amu(1,j,k) + v1*zt3
      amu(2,j,k) = amu(2,j,k) + v2*zt3
      amu(3,j,k) = amu(3,j,k) + v3*zt3
      amu(4,j,k) = amu(4,j,k) + v4*zt3
      dcu(1,j,k) = dcu(1,j,k) + vx*zt3
      dcu(2,j,k) = dcu(2,j,k) + vy*zt3
      dcu(3,j,k) = dcu(3,j,k) + vz*zt3
      cu(1,j,k) = cu(1,j,k) + ox*zt3
      cu(2,j,k) = cu(2,j,k) + oy*zt3
      cu(3,j,k) = cu(3,j,k) + oz*zt3
      amu(1,j,k1) = amu(1,j,k1) + v1*zt4
      amu(2,j,k1) = amu(2,j,k1) + v2*zt4
      amu(3,j,k1) = amu(3,j,k1) + v3*zt4
      amu(4,j,k1) = amu(4,j,k1) + v4*zt4
      dcu(1,j,k1) = dcu(1,j,k1) + vx*zt4
      dcu(2,j,k1) = dcu(2,j,k1) + vy*zt4
      dcu(3,j,k1) = dcu(3,j,k1) + vz*zt4
      cu(1,j,k1) = cu(1,j,k1) + ox*zt4
      cu(2,j,k1) = cu(2,j,k1) + oy*zt4
      cu(3,j,k1) = cu(3,j,k1) + oz*zt4
   60 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = at3*zt2
      amu(1,1,k) = amu(1,1,k) + v1*zt1
      amu(2,1,k) = amu(2,1,k) + v2*zt1
      amu(3,1,k) = amu(3,1,k) + v3*zt1
      amu(4,1,k) = amu(4,1,k) + v4*zt1
      dcu(1,1,k) = dcu(1,1,k) + vx*zt1
      dcu(2,1,k) = dcu(2,1,k) + vy*zt1
      dcu(3,1,k) = dcu(3,1,k) + vz*zt1
      cu(1,1,k) = cu(1,1,k) + ox*zt1
      cu(2,1,k) = cu(2,1,k) + oy*zt1
      cu(3,1,k) = cu(3,1,k) + oz*zt1
      amu(1,1,k1) = amu(1,1,k1) + v1*zt2
      amu(2,1,k1) = amu(2,1,k1) + v2*zt2
      amu(3,1,k1) = amu(3,1,k1) + v3*zt2
      amu(4,1,k1) = amu(4,1,k1) + v4*zt2
      dcu(1,1,k1) = dcu(1,1,k1) + vx*zt2
      dcu(2,1,k1) = dcu(2,1,k1) + vy*zt2
      dcu(3,1,k1) = dcu(3,1,k1) + vz*zt2
      cu(1,1,k1) = cu(1,1,k1) + ox*zt2
      cu(2,1,k1) = cu(2,1,k1) + oy*zt2
      cu(3,1,k1) = cu(3,1,k1) + oz*zt2
   70 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      do 80 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      amu(1,j,1) = amu(1,j,1) + v1*zt3
      amu(2,j,1) = amu(2,j,1) + v2*zt3
      amu(3,j,1) = amu(3,j,1) + v3*zt3
      amu(4,j,1) = amu(4,j,1) + v4*zt3
      dcu(1,j,1) = dcu(1,j,1) + vx*zt3
      dcu(2,j,1) = dcu(2,j,1) + vy*zt3
      dcu(3,j,1) = dcu(3,j,1) + vz*zt3
      cu(1,j,1) = cu(1,j,1) + ox*zt3
      cu(2,j,1) = cu(2,j,1) + oy*zt3
      cu(3,j,1) = cu(3,j,1) + oz*zt3
      amu(1,j,k1) = amu(1,j,k1) + v1*zt1
      amu(2,j,k1) = amu(2,j,k1) + v2*zt1
      amu(3,j,k1) = amu(3,j,k1) + v3*zt1
      amu(4,j,k1) = amu(4,j,k1) + v4*zt1
      dcu(1,j,k1) = dcu(1,j,k1) + vx*zt1
      dcu(2,j,k1) = dcu(2,j,k1) + vy*zt1
      dcu(3,j,k1) = dcu(3,j,k1) + vz*zt1
      cu(1,j,k1) = cu(1,j,k1) + ox*zt1
      cu(2,j,k1) = cu(2,j,k1) + oy*zt1
      cu(3,j,k1) = cu(3,j,k1) + oz*zt1
   80 continue
      at3 = real(sctx(nxh))
      amu(1,1,1) = cmplx(real(amu(1,1,1))+v1,aimag(amu(1,1,1))+v1*at3)
      amu(2,1,1) = cmplx(real(amu(2,1,1))+v2,aimag(amu(2,1,1))+v2*at3)
      amu(3,1,1) = cmplx(real(amu(3,1,1))+v3,aimag(amu(3,1,1))+v3*at3)
      amu(4,1,1) = cmplx(real(amu(4,1,1))+v4,aimag(amu(4,1,1))+v4*at3)
      dcu(1,1,1) = cmplx(real(dcu(1,1,1))+vx,aimag(dcu(1,1,1))+vx*at3)
      dcu(2,1,1) = cmplx(real(dcu(2,1,1))+vy,aimag(dcu(2,1,1))+vy*at3)
      dcu(3,1,1) = cmplx(real(dcu(3,1,1))+vz,aimag(dcu(3,1,1))+vz*at3)
      cu(1,1,1) = cmplx(real(cu(1,1,1))+ox,aimag(cu(1,1,1))+ox*at3)
      cu(2,1,1) = cmplx(real(cu(2,1,1))+oy,aimag(cu(2,1,1))+oy*at3)
      cu(3,1,1) = cmplx(real(cu(3,1,1))+oz,aimag(cu(3,1,1))+oz*at3)
      at3 = at1*at3
      amu(1,1,k1) = cmplx(real(amu(1,1,k1))+v1*at1,aimag(amu(1,1,k1))+v1
     1*at3)
      amu(2,1,k1) = cmplx(real(amu(2,1,k1))+v2*at1,aimag(amu(2,1,k1))+v2
     1*at3)
      amu(3,1,k1) = cmplx(real(amu(3,1,k1))+v3*at1,aimag(amu(3,1,k1))+v3
     1*at3)
      amu(4,1,k1) = cmplx(real(amu(4,1,k1))+v4*at1,aimag(amu(4,1,k1))+v4
     1*at3)
      dcu(1,1,k1) = cmplx(real(dcu(1,1,k1))+vx*at1,aimag(dcu(1,1,k1))+vx
     1*at3)
      dcu(2,1,k1) = cmplx(real(dcu(2,1,k1))+vy*at1,aimag(dcu(2,1,k1))+vy
     1*at3)
      dcu(3,1,k1) = cmplx(real(dcu(3,1,k1))+vz*at1,aimag(dcu(3,1,k1))+vz
     1*at3)
      cu(1,1,k1) = cmplx(real(cu(1,1,k1))+ox*at1,aimag(cu(1,1,k1))+ox*at
     13)
      cu(2,1,k1) = cmplx(real(cu(2,1,k1))+oy*at1,aimag(cu(2,1,k1))+oy*at
     13)
      cu(3,1,k1) = cmplx(real(cu(3,1,k1))+oz*at1,aimag(cu(3,1,k1))+oz*at
     13)
   90 continue
      return
      end
