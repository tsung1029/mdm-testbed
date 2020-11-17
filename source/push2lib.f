c-----------------------------------------------------------------------
c 2d PIC library for pushing particles and depositing charge
c push2lib.f contains procedures to process particles:
c GPOST2 deposits charge density, quadratic interpolation, STANDARD
c        optimization.
c GSPOST2 deposits charge density, quadratic interpolation, LOOKAHEAD
c         optimization.
c GSPOST2X deposits charge density, quadratic interpolation, VECTOR
c          optimization.
c GPOST2L deposits charge density, linear interpolation, STANDARD
c         optimization.
c GSPOST2L deposits charge density, linear interpolation, LOOKAHEAD
c          optimization.
c GSPOST2XL deposits charge density, linear interpolation, VECTOR
c           optimization.
c GPOST2C deposits charge density, cubic interpolation, STANDARD
c         optimization.
c GPUSH2 push particles, quadratic interpolation, STANDARD optimization.
c GSPUSH2 push particles, quadratic interpolation, LOOKAHEAD
c         optimization.
c GPUSH2L push particles, linear interpolation, STANDARD optimization.
c GSPUSH2L push particles, linear interpolation, LOOKAHEAD optimization.
c GPUSH2C push particles, cubic interpolation, STANDARD optimization.
c SORTP2Y sort particles by y grid, quadratic interpolation, memory
c         conserving algorithm.
c SORTP2YL sort particles by y grid, linear interpolation, memory
c          conserving algorithm.
c DSORTP2Y sort particles by y grid, quadratic interpolation, high
c          performance algorithm.
c DSORTP2YL sort particles by y grid, linear interpolation, high
c           performance algorithm.
c SORTP2 sort particles by grid, quadratic interpolation, memory
c        conserving algorithm.
c SORTP2L sort particles by grid, linear interpolation, memory
c         conserving algorithm.
c RMOVE2 remove particles instead of reflecting at boundary.
c PUSH2ZF update particle co-ordinates for particles with fixed
c         velocities.
c DPOST2GL deposits charge density, using gridless method.
c DPOST2GLX deposits charge density, using optimized gridless method.
c PUSH2GL push particles, using gridless method.
c PUSH2GLX push particles, using optimized gridless method.
c GCJPOST2 deposits time-centered particle current density, quadratic
c          interpolation.
c GCJPOST2L deposits time-centered particle current density, linear
c           interpolation.
c GCJPOST2C deposits time-centered particle current density, cubic
c           interpolation.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: august 28, 2010
c-----------------------------------------------------------------------
      subroutine DPOST2(part,q,qm,nop,idimp,nx,ny,nxv)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c baseline scalar version
c 43 flops/particle, 11 loads, 9 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of charge array, must be >= nx
      dimension part(idimp,nop), q(nxv,ny)
      qmh = .5*qm
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      mm = mm + 1
      if (mm.gt.ny) mm = mm - ny
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = qmh*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = qmh*(.5 + dxp)**2
      ml = mm - 1
      if (ml.lt.1) ml = ml + ny
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      if (mp.gt.ny) mp = mp - ny
      dyp = .5*(.5 + dyp)**2
c deposit charge
      q(nl,mm) = q(nl,mm) + dxl*amy
      q(nn,mm) = q(nn,mm) + amx*amy
      q(np,mm) = q(np,mm) + dxp*amy
      q(nl,ml) = q(nl,ml) + dxl*dyl
      q(nn,ml) = q(nn,ml) + amx*dyl
      q(np,ml) = q(np,ml) + dxp*dyl
      q(nl,mp) = q(nl,mp) + dxl*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mp) = q(np,mp) + dxp*dyp
   10 continue
      return
      end
      subroutine GPOST2(part,q,qm,nop,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c scalar version using guard cells
c 43 flops/particle, 11 loads, 9 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j+1,k+1) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of charge array, must be >= nx+3
c nyv = second dimension of charge array, must be >= ny+3
      dimension part(idimp,nop), q(nxv,nyv)
      qmh = .5*qm
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
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
c deposit charge
      q(nl,mm) = q(nl,mm) + dxl*amy
      q(nn,mm) = q(nn,mm) + amx*amy
      q(np,mm) = q(np,mm) + dxp*amy
      q(nl,ml) = q(nl,ml) + dxl*dyl
      q(nn,ml) = q(nn,ml) + amx*dyl
      q(np,ml) = q(np,ml) + dxp*dyl
      q(nl,mp) = q(nl,mp) + dxl*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mp) = q(np,mp) + dxp*dyp
   10 continue
      return
      end
      subroutine GSPOST2(part,q,qm,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 43 flops/particle, 11 loads, 9 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j+1,k+1) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of charge array, must be >= nx+3
c nxyv = dimension of charge array, must be >= nxv*(ny+3)
      dimension part(idimp,nop), q(nxyv)
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      qmh = .5*qm
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
c deposit charge
      dx = q(mn) + dxl*amy
      dy = q(mn+1) + amx*amy
      amy = q(mn+2) + dxp*amy
      dx1 = q(ml) + dxl*dyl
      dy1 = q(ml+1) + amx*dyl
      dyl = q(ml+2) + dxp*dyl
      dxl = q(mp) + dxl*dyp
      amx = q(mp+1) + amx*dyp
      dyp = q(mp+2) + dxp*dyp
      q(mn) = dx
      q(mn+1) = dy
      q(mn+2) = amy
      q(ml) = dx1
      q(ml+1) = dy1
      q(ml+2) = dyl
      q(mp) = dxl
      q(mp+1) = amx
      q(mp+2) = dyp
   10 continue
c deposit charge for last particle
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
c deposit charge
      q(mn) = q(mn) + dxl*amy
      q(mn+1) = q(mn+1) + amx*amy
      q(mn+2) = q(mn+2) + dxp*amy
      q(ml) = q(ml) + dxl*dyl
      q(ml+1) = q(ml+1) + amx*dyl
      q(ml+2) = q(ml+2) + dxp*dyl
      q(mp) = q(mp) + dxl*dyp
      q(mp+1) = q(mp+1) + amx*dyp
      q(mp+2) = q(mp+2) + dxp*dyp
      return
      end
      subroutine SPOST2X(part,q,qm,nop,idimp,nx,ny,nxv,nxvy)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with 1d addressing
c 43 flops/particle, 29 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n) = charge density at grid point j,k where n = j + nxv*(k - 1)
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of charge array, must be >= nx
c nxvy = nxv*ny, dimension of charge array
c vectorized version with 1d addressing of case 2 in
c v.k.decyk et al, computers in physics 10, 290  (1996).
      parameter(npp=512)
      dimension part(idimp,nop), q(nxvy)
      dimension nn(9,npp), amxy(9,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      qmh = .5*qm
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb) + .5
      dxl = part(1,j+jb) - float(n)
      n = n + 1
      if (n.gt.nx) n = n - nx
      amx = qm*(.75 - dxl*dxl)
      np = n + 1
      if (np.gt.nx) np = np - nx
      dxp = qmh*(.5 + dxl)**2
      nl = n - 1
      if (nl.lt.1) nl = nl + nx
      dxl = qmh*(.5 - dxl)**2
      m = part(2,j+jb) + .5
      dyl = part(2,j+jb) - float(m)
      if (m.ge.ny) m = m - ny
      amy = .75 - dyl*dyl
      mp = m + 1
      if (mp.ge.ny) mp = mp - ny
      dyp = .5*(.5 + dyl)**2
      ml = m - 1
      if (ml.lt.0) ml = ml + ny
      dyl = .5*(.5 - dyl)**2
      m = nxv*m
      mp = nxv*mp
      ml = nxv*ml
      nn(1,j) = n + m
      nn(2,j) = np + m
      nn(3,j) = nl + m
      nn(4,j) = n + mp
      nn(5,j) = np + mp
      nn(6,j) = nl + mp
      nn(7,j) = n + ml
      nn(8,j) = np + ml
      nn(9,j) = nl + ml
      amxy(1,j) = amx*amy
      amxy(2,j) = dxp*amy
      amxy(3,j) = dxl*amy
      amxy(4,j) = amx*dyp
      amxy(5,j) = dxp*dyp
      amxy(6,j) = dxl*dyp
      amxy(7,j) = amx*dyl
      amxy(8,j) = dxp*dyl
      amxy(9,j) = dxl*dyl
   10 continue
c deposit charge
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 9
      q(nn(i,j)) = q(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GSPOST2X(part,q,qm,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with guard cells and 1d addressing
c 43 flops/particle, 29 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(.75-dx**2)*(.75-dy**2)
c q(n+1,m)=.5*qm*((.5+dx)**2)*(.75-dy**2)
c q(n-1,m)=.5*qm*((.5-dx)**2)*(.75-dy**2)
c q(n,m+1)=.5*qm*(.75-dx**2)*(.5+dy)**2
c q(n+1,m+1)=.25*qm*((.5+dx)**2)*(.5+dy)**2
c q(n-1,m+1)=.25*qm*((.5-dx)**2)*(.5+dy)**2
c q(n,m-1)=.5*qm*(.75-dx**2)*(.5-dy)**2
c q(n+1,m-1)=.25*qm*((.5+dx)**2)*(.5-dy)**2
c q(n-1,m-1)=.25*qm*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n) = charge density at grid point j,k where n = j + nxv*k + 1
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of charge array, must be >= nx+3
c nxyv = dimension of charge array, must be >= nxv*(ny+3)
c vectorized version case 6 in
c v.k.decyk et al, computers in physics 10, 290  (1996).
      parameter(npp=512)
      dimension part(idimp,nop), q(nxyv)
      dimension nn(9,npp), amxy(9,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      qmh = .5*qm
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb) + .5
      m = part(2,j+jb) + .5
      dxp = part(1,j+jb) - float(n)
      dyp = part(2,j+jb) - float(m)
      n = n + 1
      m = nxv*m
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = m + n
      dxl = .5*qm*(.5 - dxp)**2
      dxp = .5*qm*(.5 + dxp)**2
      mn = ml + nxv
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv
      nn(1,j) = mn
      nn(2,j) = mn + 1
      nn(3,j) = mn + 2
      nn(4,j) = ml
      nn(5,j) = ml + 1
      nn(6,j) = ml + 2
      nn(7,j) = mp
      nn(8,j) = mp + 1
      nn(9,j) = mp + 2
      amxy(1,j) = dxl*amy
      amxy(2,j) = amx*amy
      amxy(3,j) = dxp*amy
      amxy(4,j) = dxl*dyl
      amxy(5,j) = amx*dyl
      amxy(6,j) = dxp*dyl
      amxy(7,j) = dxl*dyp
      amxy(8,j) = amx*dyp
      amxy(9,j) = dxp*dyp
   10 continue
c deposit charge
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 9
      q(nn(i,j)) = q(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine DPOST2L(part,q,qm,nop,idimp,nx,ny,nxv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c baseline scalar version
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of charge array, must be >= nx
      dimension part(idimp,nop), q(nxv,ny)
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
      if (mp.gt.ny) mp = mp - ny
c deposit charge
      q(nn,mm) = q(nn,mm) + amx*amy
      q(np,mm) = q(np,mm) + dxp*amy
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mp) = q(np,mp) + dxp*dyp
   10 continue
      return
      end
      subroutine GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of charge array, must be >= nx+1
c nxy = second dimension of charge array, must be >= ny+1
      dimension part(idimp,nop), q(nxv,nyv)
c find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit charge
      !print*, 'mp =', mp
      q(np,mp) = q(np,mp) + dxp*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mm) = q(np,mm) + dxp*amy
      q(nn,mm) = q(nn,mm) + amx*amy
   10 continue
      return
      end
      subroutine GSPOST2L(part,q,qm,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 17 flops/particle, 6 loads, 4 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of charge array, must be >= nx+1
c nxyv = dimension of charge array, must be >= nxv*(ny+1)
      dimension part(idimp,nop), q(nxyv)
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
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
c deposit charge
      dx1 = q(mp+1) + dxp*dyp
      dyp = q(mp) + amx*dyp
      dxp = q(mm+1) + dxp*amy
      amy = q(mm) + amx*amy
      q(mp+1) = dx1
      q(mp) = dyp
      q(mm+1) = dxp
      q(mm) = amy
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit charge
      q(mp+1) = q(mp+1) + dxp*dyn
      q(mp) = q(mp) + amx*dyn
      q(mm+1) = q(mm+1) + dxp*amy
      q(mm) = q(mm) + amx*amy
      return
      end
      subroutine SPOST2XL(part,q,qm,nop,idimp,nx,ny,nxv,nxvy)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with 1d addressing
c 17 flops/particle, 14 loads, 12 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n) = charge density at grid point j,k where n = j + nxv*(k - 1)
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of charge array, must be >= nx
c nxvy = nxv*ny, dimension of charge array
      parameter(npp=128)
      dimension part(idimp,nop), q(nxvy)
      dimension nn(4,npp), amxy(4,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb)
      dxp = qm*(part(1,j+jb) - float(n))
      n = n + 1
      amx = qm - dxp
      np = n + 1
      if (np.gt.nx) np = np - nx
      m = part(2,j+jb)
      dyp = part(2,j+jb) - float(m)
      amy = 1. - dyp
      mp = m + 1
      if (mp.ge.ny) mp = mp - ny
      m = nxv*m
      mp = nxv*mp
      nn(1,j) = n + m
      nn(2,j) = np + m
      nn(3,j) = n + mp
      nn(4,j) = np + mp
      amxy(1,j) = amx*amy
      amxy(2,j) = dxp*amy
      amxy(3,j) = amx*dyp
      amxy(4,j) = dxp*dyp
   10 continue
c deposit charge
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 4
      q(nn(i,j)) = q(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GSPOST2XL(part,q,qm,nop,idimp,nxv,nxyv)
c for 2d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with guard cells and 1d addressing
c 17 flops/particle, 14 loads, 12 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m)=qm*(1.-dx)*(1.-dy)
c q(n+1,m)=qm*dx*(1.-dy)
c q(n,m+1)=qm*(1.-dx)*dy
c q(n+1,m+1)=qm*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n) = charge density at grid point j,k where n = j + nxv*k + 1
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first virtual dimension of charge array, must be >= nx+1
c nxyv = dimension of charge array, must be >= nxv*(ny+1)
      parameter(npp=128)
      dimension part(idimp,nop), q(nxyv)
      dimension nn(4,npp), amxy(4,npp)
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
c find interpolation weights
      do 10 j = 1, npb
      n = part(1,j+jb)
      m = part(2,j+jb)
      dxp = qm*(part(1,j+jb) - float(n))
      dyp = part(2,j+jb) - float(m)
      n = n + 1
      m = nxv*m
      amx = qm - dxp
      m = m + n
      amy = 1. - dyp
      mp = m + nxv
      nn(4,j) = m
      nn(3,j) = m + 1
      nn(2,j) = mp
      nn(1,j) = mp + 1
      amxy(1,j) = dxp*dyp
      amxy(2,j) = amx*dyp
      amxy(3,j) = dxp*amy
      amxy(4,j) = amx*amy
   10 continue
c deposit charge
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 4
      q(nn(i,j)) = q(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GPOST2C(part,q,qm,nop,idimp,nxv,nyv)
c for 2d code, this subroutine calculates particle charge density
c using third-order spline interpolation, periodic boundaries
c scalar version using guard cells
c 68 flops/particle, 18 loads, 16 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n-1,m-1)=(qm/6.)*((1-dx)**3)*(1./6.)*((1-dy)**3)
c q(n,m-1)=(qm)*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*((1-dy)**3)
c q(n+1,m-1)=(qm)*(2./3.-.5*((1-dx)**2)*(1+dx))*(1./6.)*((1-dy)**3)
c q(n+2,m-1)=(qm/6.)*(dx**3)*(1./6.)*((1-dy)**3)
c q(n-1,m)=(qm/6.)*((1-dx)**3)*(2./3.-.5*dy*dy*(2-dy))
c q(n,m)=(qm)*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*dy*dy*(2-dy))
c q(n+1,m)=(qm)*(2./3.-.5*((1-dx)**2)*(2./3.-.5*dy*dy*(2-dy))
c q(n+2,m)=(qm/6.)*(dx**3)*(2./3.-.5*dy*dy*(2-dy))
c q(n-1,m+1)=(qm/6.)*((1-dx)**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c q(n,m+1)=(qm)*(2./3.-.5*dx*dx*(2-dx))*(2./3.-.5*((1-dy)**2)*(1+dy))
c q(n+1,m+1)=(qm)*(2./3.-.5*((1-dx)**2)*(2./3.-.5*((1-dy)**2)*(1+dy))
c q(n+2,m+1)=(qm/6.)*(dx**3)*(2./3.-.5*((1-dy)**2)*(1+dy))
c q(n-1,m+2)=(qm/6.)*((1-dx)**3)*(1./6.)*(dy**3)
c q(n,m+2)=(qm)*(2./3.-.5*dx*dx*(2-dx))*(1./6.)*(dy**3)
c q(n+1,m+2)=(qm)*(2./3.-.5*((1-dx)**2)*(1./6.)*(dy**3)
c q(n+2,m+2)=(qm/6.)*(dx**3)*(1./6.)*(dy**3)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(j,k) = charge density at grid point j,k
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nxv = first dimension of charge array, must be >= nx+5
c nxy = second dimension of charge array, must be >= ny+5
      implicit none
      integer nop, idimp, nxv, nyv
      real part, q, qm
      dimension part(idimp,nop), q(nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
c find interpolation weights
      do 10 j = 1, nop
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
c deposit charge
      q(nl,ml) = q(nl,ml) + dxl*dyl
      q(nn,ml) = q(nn,ml) + dxn*dyl
      q(np,ml) = q(np,ml) + dxp*dyl
      q(nq,ml) = q(nq,ml) + dxq*dyl
      q(nl,mm) = q(nl,mm) + dxl*dyn
      q(nn,mm) = q(nn,mm) + dxn*dyn
      q(np,mm) = q(np,mm) + dxp*dyn
      q(nq,mm) = q(nq,mm) + dxq*dyn
      q(nl,mp) = q(nl,mp) + dxl*dyp
      q(nn,mp) = q(nn,mp) + dxn*dyp
      q(np,mp) = q(np,mp) + dxp*dyp
      q(nq,mp) = q(nq,mp) + dxq*dyp
      q(nl,mq) = q(nl,mq) + dxl*dyq
      q(nn,mq) = q(nn,mq) + dxn*dyq
      q(np,mq) = q(np,mq) + dxp*dyq
      q(nq,mq) = q(nq,mq) + dxq*dyq
   10 continue
      return
      end
      subroutine PUSH2(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with periodic boundary conditions.
c baseline scalar version
c 82 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      double precision sum1
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      qtm = qbm*dt
      sum1 = 0.0d0
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      mm = mm + 1
      if (mm.gt.ny) mm = mm - ny
      amx = .75 - dxp*dxp
      amy = .75 - dyp*dyp
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dxl = .5*(.5 - dxp)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*(.5 + dxp)**2
      ml = mm - 1
      if (ml.lt.1) ml = ml + ny
      dyl = .5*(.5 - dyp)**2
      mp = mm + 1
      if (mp.gt.ny) mp = mp - ny
      dyp = .5*(.5 + dyp)**2
c find acceleration
      dx = amy*(dxl*fx(nl,mm) + amx*fx(nn,mm) + dxp*fx(np,mm)) + dyl*(dx
     1l*fx(nl,ml) + amx*fx(nn,ml) + dxp*fx(np,ml)) + dyp*(dxl*fx(nl,mp)
     2+ amx*fx(nn,mp) + dxp*fx(np,mp))
      dy = amy*(dxl*fy(nl,mm) + amx*fy(nn,mm) + dxp*fy(np,mm)) + dyl*(dx
     1l*fy(nl,ml) + amx*fy(nn,ml) + dxp*fy(np,ml)) + dyp*(dxl*fy(nl,mp)
     2+ amx*fy(nn,mp) + dxp*fy(np,mp))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions.
c scalar version using guard cells
c 82 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = second dimension of field array, must be >= nx+3
c nyv = third dimension of field array, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      qtm = qbm*dt
      sum1 = 0.0d0
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
c find acceleration
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))              
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GSPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,ipb
     1c)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with various boundary conditions.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 82 flops/particle, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = second virtual dimension of field array, must be >= nx+3
c nxyv = dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtm = qbm*dt
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
c find acceleration
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))              
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
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
c find acceleration
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp)+ amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))              
c new velocity
      dx = part(3,nop) + qtm*dx
      dy = part(4,nop) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop))**2 + (dy + part(4,nop))**2
      part(3,nop) = dx
      part(4,nop) = dy
c new position
      dx = part(1,nop) + dx*dt
      dy = part(2,nop) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c normalize kinetic energy
   20 ek = ek + .125*sum1
      return
      end
      subroutine PUSH2L(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions.
c baseline scalar version
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      double precision sum1
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      qtm = qbm*dt
      sum1 = 0.0d0
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = 1. - dyp
      mp = mm + 1
      if (mp.gt.ny) mp = mp - ny
c find acceleration
      dx = amy*(amx*fx(nn,mm) + dxp*fx(np,mm)) + dyp*(amx*fx(nn,mp) + dx
     1p*fx(np,mp))
      dy = amy*(amx*fy(nn,mm) + dxp*fy(np,mm)) + dyp*(amx*fy(nn,mp) + dx
     1p*fy(np,mp))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c scalar version using guard cells
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxy = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      qtm = qbm*dt
      sum1 = 0.0d0
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
c find acceleration
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with various boundary conditions.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 44 flops/particle, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = second virtual dimension of field arrays, must be >= nx+1
c nxyv = dimension of field arrays, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtm = qbm*dt
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
c find acceleration
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c push last particle
      nn = nnn + 1
      mm = nxv*mmn
      mm = mm + nn
      amx = 1. - dxn
      mp = mm + nxv
      amy = 1. - dyn
c find acceleration
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c new velocity
      dx = part(3,nop) + qtm*dx
      dy = part(4,nop) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,nop))**2 + (dy + part(4,nop))**2
      part(3,nop) = dx
      part(4,nop) = dy
c new position
      dx = part(1,nop) + dx*dt
      dy = part(2,nop) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(3,nop) = -part(3,nop)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
c normalize kinetic energy
   20 ek = ek + .125*sum1
      return
      end
      subroutine GPUSH2C(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order spline
c interpolation in space, with various boundary conditions.
c scalar version using guard cells
c 120 flops/particle, 36 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
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
c and similarly for fy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nxy = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, qbm, dt, ek
      dimension part(idimp,nop), fxy(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtm, edgelx, edgely, edgerx, edgery
      real dx, dy, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq, at1, at2
      double precision sum1
      qtm = qbm*dt
      sum1 = 0.0d0
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
c find acceleration
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
c new velocity
      dx = part(3,j) + qtm*dx
      dy = part(4,j) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,j))**2 + (dy + part(4,j))**2
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine SORTP2Y(part,pt,ip,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c quadratic interpolation
c part = particle array
c part(2,n) = position y of particle n
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(ny1)
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = part(2,j) + 0.5
      m = m + 1
      npic(m) = npic(m) + 1
      ip(j) = m
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip(j)
      npic(m) = npic(m) + 1
      ip(j) = npic(m)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, nop
      pt(ip(j)) = part(i,j)
   50 continue
      do 60 j = 1, nop
      part(i,j) = pt(j)
   60 continue
   70 continue
      return
      end
      subroutine SORTP2YL(part,pt,ip,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c linear interpolation
c part = particle array
c part(2,n) = position y of particle n
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(ny1)
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = part(2,j)
      m = m + 1
      npic(m) = npic(m) + 1
      ip(j) = m
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip(j)
      npic(m) = npic(m) + 1
      ip(j) = npic(m)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, nop
      pt(ip(j)) = part(i,j)
   50 continue
      do 60 j = 1, nop
      part(i,j) = pt(j)
   60 continue
   70 continue
      return
      end
      subroutine DSORTP2Y(parta,partb,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c quadratic interpolation
c parta/partb = input/output particle arrays
c parta(2,n) = position y of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      implicit none
      integer npic, idimp, nop, ny1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
c local data
      integer i, j, k, m, isum, ist, ip
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(2,j) + 0.5
      m = m + 1
      npic(m) = npic(m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      m = parta(2,j) + 0.5
      m = m + 1
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(m) = ip
   50 continue
      return
      end
      subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c linear interpolation
c parta/partb = input/output particle arrays
c parta(2,n) = position y of particle n
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      implicit none
      integer npic, idimp, nop, ny1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
c local data
      integer i, j, k, m, isum, ist, ip
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      m = parta(2,j)
      m = m + 1
      npic(m) = npic(m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      m = parta(2,j)
      m = m + 1
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(m) = ip
   50 continue
      return
      end
      subroutine IPSORTP2Y(part,pt2,ip2,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c quadratic interpolation
c part = particle array
c part(2,n) = position y of particle n
c pt2 = scratch array for reordering particles
c ip2 = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      implicit none
      integer idimp, nop, ny1
      real part, pt2
      integer ip2, npic
      dimension part(idimp,nop), pt2(idimp,2)
      dimension ip2(2,nop), npic(ny1)
c local data
      integer i, j, k, n, m, isum, ist
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid, and mark unprocessed locations
      do 20 j = 1, nop
      m = part(2,j) + 0.5
      m = m + 1
      npic(m) = npic(m) + 1
      ip2(1,j) = m
      ip2(2,j) = 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip2(1,j)
      npic(m) = npic(m) + 1
      ip2(1,j) = npic(m)
   40 continue
c reorder particles following trail of reordered particles
      n = 0
c find first particle which is moving to another grid
   50 n = n + 1
c check to see if we are done
      if (n.gt.nop) return
c check to see if next one points to itself
      if (n.eq.ip2(1,n)) ip2(2,n) = 0
      if (ip2(2,n).eq.0) go to 50
      m = n
      ip2(2,m) = 0
      do 60 i = 1, idimp
      pt2(i,1) = part(i,m)
   60 continue
      do 100 j = 1, nop
c swap components
      m = ip2(1,m)
      do 70 i = 1, idimp
      pt2(i,2) = part(i,m)
      part(i,m) = pt2(i,1)
      pt2(i,1) = pt2(i,2)
   70 continue
c if next location is not processed, mark it as (about to be) processed
      if (ip2(2,m).eq.1) then
         ip2(2,m) = 0
c start new cycle
      else
c find new starting location for cycle
         if (ip2(2,n).eq.0) then
   80       n = n + 1
c check to see if we are done
            if (n.gt.nop) return
c check to see if next one points to itself
            if (n.eq.ip2(1,n)) ip2(2,n) = 0
            if (ip2(2,n).eq.0) go to 80
         endif
         m = n
         ip2(2,m) = 0
         do 90 i = 1, idimp
         pt2(i,1) = part(i,m)
   90    continue
      endif
  100 continue
      return
      end
      subroutine IPSORTP2YL(part,pt2,ip2,npic,idimp,nop,ny1)
c this subroutine sorts particles by y grid
c linear interpolation
c part = particle array
c part(2,n) = position y of particle n
c pt2 = scratch array for reordering particles
c ip2 = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c ny1 = system length in y direction + 1
      implicit none
      integer idimp, nop, ny1
      real part, pt2
      integer ip2, npic
      dimension part(idimp,nop), pt2(idimp,2)
      dimension ip2(2,nop), npic(ny1)
c local data
      integer i, j, k, n, m, isum, ist
c clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
c find how many particles in each grid, and mark unprocessed locations
      do 20 j = 1, nop
      m = part(2,j)
      m = m + 1
      npic(m) = npic(m) + 1
      ip2(1,j) = m
      ip2(2,j) = 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip2(1,j)
      npic(m) = npic(m) + 1
      ip2(1,j) = npic(m)
   40 continue
c reorder particles following trail of reordered particles
      n = 0
c find first particle which is moving to another grid
   50 n = n + 1
c check to see if we are done
      if (n.gt.nop) return
c check to see if next one points to itself
      if (n.eq.ip2(1,n)) ip2(2,n) = 0
      if (ip2(2,n).eq.0) go to 50
      m = n
      ip2(2,m) = 0
      do 60 i = 1, idimp
      pt2(i,1) = part(i,m)
   60 continue
      do 100 j = 1, nop
c swap components
      m = ip2(1,m)
      do 70 i = 1, idimp
      pt2(i,2) = part(i,m)
      part(i,m) = pt2(i,1)
      pt2(i,1) = pt2(i,2)
   70 continue
c if next location is not processed, mark it as (about to be) processed
      if (ip2(2,m).eq.1) then
         ip2(2,m) = 0
c start new cycle
      else
c find new starting location for cycle
         if (ip2(2,n).eq.0) then
   80       n = n + 1
c check to see if we are done
            if (n.gt.nop) return
c check to see if next one points to itself
            if (n.eq.ip2(1,n)) ip2(2,n) = 0
            if (ip2(2,n).eq.0) go to 80
         endif
         m = n
         ip2(2,m) = 0
         do 90 i = 1, idimp
         pt2(i,1) = part(i,m)
   90    continue
      endif
  100 continue
      return
      end
      subroutine SORTP2(part,pt,ip,npic,idimp,nop,nx1,nxy1)
c this subroutine sorts particles by x,y grid
c quadratic interpolation
c part = particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
c nxy1 = nx1*ny1, where ny1 = system length in y direction + 1
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(nxy1)
c clear counter array
      do 10 k = 1, nxy1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = part(1,j) + 0.5
      n = n + 1
      m = part(2,j) + 0.5
      m = n + nx1*m
      npic(m) = npic(m) + 1
      ip(j) = m
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nxy1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip(j)
      npic(m) = npic(m) + 1
      ip(j) = npic(m)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, nop
      pt(ip(j)) = part(i,j)
   50 continue
      do 60 j = 1, nop
      part(i,j) = pt(j)
   60 continue
   70 continue
      return
      end
      subroutine SORTP2L(part,pt,ip,npic,idimp,nop,nx1,nxy1)
c this subroutine sorts particles by x,y grid
c linear interpolation
c part = particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
c nxy1 = nx1*ny1, where ny1 = system length in y direction + 1
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(nxy1)
c clear counter array
      do 10 k = 1, nxy1
      npic(k) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, nop
      n = part(1,j)
      n = n + 1
      m = part(2,j)
      m = n + nx1*m
      npic(m) = npic(m) + 1
      ip(j) = m
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nxy1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip(j)
      npic(m) = npic(m) + 1
      ip(j) = npic(m)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, nop
      pt(ip(j)) = part(i,j)
   50 continue
      do 60 j = 1, nop
      part(i,j) = pt(j)
   60 continue
   70 continue
      return
      end
      subroutine IPSORTP2(part,pt2,ip2,npic,idimp,nop,nx1,nxy1)
c this subroutine sorts particles by x,y grid
c quadratic interpolation
c part = particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c pt2 = scratch array for reordering particles
c ip2 = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
c nxy1 = nx1*ny1, where ny1 = system length in y direction + 1
      implicit none
      integer idimp, nop, nx1, nxy1
      real part, pt2
      integer ip2, npic
      dimension part(idimp,nop), pt2(idimp,2)
      dimension ip2(2,nop), npic(nxy1)
c local data
      integer i, j, k, n, m, isum, ist
c clear counter array
      do 10 k = 1, nxy1
      npic(k) = 0
   10 continue
c find how many particles in each grid, and mark unprocessed locations
      do 20 j = 1, nop
      n = part(1,j) + 0.5
      n = n + 1
      m = part(2,j) + 0.5
      m = n + nx1*m
      npic(m) = npic(m) + 1
      ip2(1,j) = m
      ip2(2,j) = 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nxy1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip2(1,j)
      npic(m) = npic(m) + 1
      ip2(1,j) = npic(m)
   40 continue
c reorder particles following trail of reordered particles
      n = 0
c find first particle which is moving to another grid
   50 n = n + 1
c check to see if we are done
      if (n.gt.nop) return
c check to see if next one points to itself
      if (n.eq.ip2(1,n)) ip2(2,n) = 0
      if (ip2(2,n).eq.0) go to 50
      m = n
      ip2(2,m) = 0
      do 60 i = 1, idimp
      pt2(i,1) = part(i,m)
   60 continue
      do 100 j = 1, nop
c swap components
      m = ip2(1,m)
      do 70 i = 1, idimp
      pt2(i,2) = part(i,m)
      part(i,m) = pt2(i,1)
      pt2(i,1) = pt2(i,2)
   70 continue
c if next location is not processed, mark it as (about to be) processed
      if (ip2(2,m).eq.1) then
         ip2(2,m) = 0
c start new cycle
      else
c find new starting location for cycle
         if (ip2(2,n).eq.0) then
   80       n = n + 1
c check to see if we are done
            if (n.gt.nop) return
c check to see if next one points to itself
            if (n.eq.ip2(1,n)) ip2(2,n) = 0
            if (ip2(2,n).eq.0) go to 80
         endif
         m = n
         ip2(2,m) = 0
         do 90 i = 1, idimp
         pt2(i,1) = part(i,m)
   90    continue
      endif
  100 continue
      return
      end
      subroutine IPSORTP2L(part,pt2,ip2,npic,idimp,nop,nx1,nxy1)
c this subroutine sorts particles by x,y grid
c linear interpolation
c part = particle array
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c pt2 = scratch array for reordering particles
c ip2 = index array for reordering particles
c npic = address offset for reordering particles
c idimp = size of phase space = 4
c nop = number of particles
c nx1 = system length in x direction + 1
c nxy1 = nx1*ny1, where ny1 = system length in y direction + 1
      implicit none
      integer idimp, nop, nx1, nxy1
      real part, pt2
      integer ip2, npic
      dimension part(idimp,nop), pt2(idimp,2)
      dimension ip2(2,nop), npic(nxy1)
c local data
      integer i, j, k, n, m, isum, ist
c clear counter array
      do 10 k = 1, nxy1
      npic(k) = 0
   10 continue
c find how many particles in each grid, and mark unprocessed locations
      do 20 j = 1, nop
      n = part(1,j)
      n = n + 1
      m = part(2,j)
      m = n + nx1*m
      npic(m) = npic(m) + 1
      ip2(1,j) = m
      ip2(2,j) = 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nxy1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, nop
      m = ip2(1,j)
      npic(m) = npic(m) + 1
      ip2(1,j) = npic(m)
   40 continue
c reorder particles following trail of reordered particles
      n = 0
c find first particle which is moving to another grid
   50 n = n + 1
c check to see if we are done
      if (n.gt.nop) return
c check to see if next one points to itself
      if (n.eq.ip2(1,n)) ip2(2,n) = 0
      if (ip2(2,n).eq.0) go to 50
      m = n
      ip2(2,m) = 0
      do 60 i = 1, idimp
      pt2(i,1) = part(i,m)
   60 continue
      do 100 j = 1, nop
c swap components
      m = ip2(1,m)
      do 70 i = 1, idimp
      pt2(i,2) = part(i,m)
      part(i,m) = pt2(i,1)
      pt2(i,1) = pt2(i,2)
   70 continue
c if next location is not processed, mark it as (about to be) processed
      if (ip2(2,m).eq.1) then
         ip2(2,m) = 0
c start new cycle
      else
c find new starting location for cycle
         if (ip2(2,n).eq.0) then
   80       n = n + 1
c check to see if we are done
            if (n.gt.nop) return
c check to see if next one points to itself
            if (n.eq.ip2(1,n)) ip2(2,n) = 0
            if (ip2(2,n).eq.0) go to 80
         endif
         m = n
         ip2(2,m) = 0
         do 90 i = 1, idimp
         pt2(i,1) = part(i,m)
   90    continue
      endif
  100 continue
      return
      end
      subroutine RMOVE2(part,ihole,nx,ny,idimp,nop,ntmax,ipbc)
c this subroutine removes particles which would normally be reflected
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c ihole = location of holes left in particle arrays
c nx/ny = system length in x/y direction
c idimp = size of phase space = 4
c nop = number of particles
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,-1,-2,-3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      real part
      integer ihole
      integer nx, ny, idimp, nop, ntmax, ipbc
      dimension part(idimp,nop)
      dimension ihole(ntmax)
c local data
      integer idps
      parameter(idps=2)
      integer i, j, j1, j2, nter, np, jss
      dimension jss(idps)
      real edgelx, edgely, edgerx, edgery, dx, dy
c set boundary values
      if (ipbc.eq.(-1)) then
         edgelx = 0.
         edgely = 0.
         edgerx = float(nx)
         edgery = float(ny)
      else if (ipbc.eq.(-2)) then
         edgelx = 1.
         edgely = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
      else if (ipbc.eq.(-3)) then
         edgelx = 1.
         edgely = 0.
         edgerx = float(nx-1)
         edgery = float(ny)
      endif
      nter = 0
      np = nop
c buffer outgoing particles
   10 jss(1) = 0
      jss(2) = 0
      do 20 j = 1, np
      dx = part(1,j)
      dy = part(2,j)
c periodic boundary conditions
      if (ipbc.eq.(-1)) then
         if (dx.lt.edgelx) part(1,j) = dx + edgerx
         if (dx.ge.edgerx) part(1,j) = part(1,j) - edgerx
         if (dy.lt.edgely) part(2,j) = dy + edgery
         if (dy.ge.edgery) part(2,j) = part(2,j) - edgery
c reflecting boundary conditions
      else if (ipbc.eq.(-2)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1).lt.ntmax) then
               jss(1) = jss(1) + 1
               ihole(jss(1)) = j
            else
               jss(2) = 1
               go to 30
            endif
         else if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            if (jss(1).lt.ntmax) then
               jss(1) = jss(1) + 1
               ihole(jss(1)) = j
            else
               jss(2) = 1
               go to 30
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.(-3)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1).lt.ntmax) then
               jss(1) = jss(1) + 1
               ihole(jss(1)) = j
            else
               jss(2) = 1
               go to 30
            endif
         endif
         if (dy.lt.edgely) part(2,j) = dy + edgery
         if (dy.ge.edgery) part(2,j) = part(2,j) - edgery
      endif
   20 continue
   30 continue
c fill up holes in particle array with particles from bottom
      do 50 j = 1, jss(1)
      j1 = np - j + 1
      j2 = jss(1) - j + 1
      if (j1.gt.ihole(j2)) then
c move particle only if it is below current hole
         do 40 i = 1, idimp
         part(i,ihole(j2)) = part(i,j1)
   40    continue
      endif
   50 continue
      np = np - jss(1)
c check if buffer overflowed and more particles remain to be checked
      if (jss(2).gt.0) then
         nter = nter + 1
         go to 10
      endif
c information
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, ntmax=', ntmax
      endif
      nop = np
      return
      end
      subroutine PUSH2ZF(part,dt,ek,idimp,nop,nx,ny,ipbc)
c for 2d code, this subroutine updates particle co-ordinates for
c particles with fixed velocities, with various boundary conditions.
c 7 flops/particle, 4 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t+dt/2)**2+vy(t+dt/2)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, ipbc
      real part, dt, ek
      dimension part(idimp,nop)
      integer j
      real edgelx, edgely, edgerx, edgery, dx, dy
      double precision sum1
      sum1 = 0.0d0
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
      do 10 j = 1, nop
c velocity
      dx = part(3,j)
      dy = part(4,j)
c average kinetic energy
      sum1 = sum1 + dx**2 + dy**2
c new position
      dx = part(1,j) + dx*dt
      dy = part(2,j) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
      subroutine DPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxvh,nyv)
c for 2d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 16*(NX/2)*(NY/2) flops/particle
c plus (NX/2 + NY/2) sines and cosines/particle
c input: all, output: q
c charge density is calculated from the expression:
c q(n,m) = sum(qm*exp(-sqrt(-1)*2*n*pi*x/nx)*exp(-sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxvh = first dimension of charge array, must be >= nx/2
c nyv = second dimension of charge array, must be >= ny
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv
      real part, qm
      complex q, sctx
      dimension part(idimp,nop), q(nxvh,nyv), sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3
      complex zt1, zt2, zt3
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qmn = qm/real(nx*ny)
c find fourier components
      do 50 i = 1, nop
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),-sin(dkx))
   10 continue
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = qmn*cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      q(j,k) = q(j,k) + zt1*zt3
      q(j,k1) = q(j,k1) + zt2*zt3
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      q(1,k) = q(1,k) + zt1
      q(1,k1) = q(1,k1) + at3*zt2
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = qmn*cos(dky)
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      q(j,1) = q(j,1) + qmn*zt3
      q(j,k1) = q(j,k1) + at1*zt3
   40 continue
      at3 = real(sctx(nxh))
      q(1,1) = cmplx(real(q(1,1))+qmn,aimag(q(1,1))+qmn*at3)
      q(1,k1) = cmplx(real(q(1,k1))+at1,aimag(q(1,k1))+at1*at3)
   50 continue
      return
      end
      subroutine DPOST2GLX(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,npp)
c for 2d code, this subroutine calculates particle charge density
c using gridless spectral version, periodic boundaries
c baseline scalar version
c 16*(NX/2)*(NY/2) flops/particle
c plus (NX/2 + NY/2) sines and cosines/particle
c input: all, output: q
c charge density is calculated from the expression:
c q(n,m) = sum(qm*exp(-sqrt(-1)*2*n*pi*x/nx)*exp(-sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c q(n,m) = charge density at fourier grid point n,m
c sctxp = scratch array for sines and cosines
c qm = charge on particle, in units of e
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxvh = first dimension of charge array, must be >= nx/2
c nyv = second dimension of charge array, must be >= ny
c npp = number of particles to be processed together
c optimized version
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv, npp
      real part, qm
      complex q, sctxp
      dimension part(idimp,nop), q(nxvh,nyv), sctxp(nxvh,npp)
      integer i, j, k, l, k1, nxh, nyh, ny2, npb, ipp, ie, ib
      real qmn, dnx, dny, dkx, dky, at1, at3
      complex zt1, zt2, zt3
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      qmn = qm/real(nx*ny)
c find fourier components
c outer loop over blocks of particles
      do 80 l = 1, ipp
      ie = l*npb
      ib = ie - npb
      if (ie.gt.nop) npb = nop - ib
      do 20 i = 1, npb
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i+ib)
      sctxp(j,i) = cmplx(cos(dkx),-sin(dkx))
   10 continue
   20 continue
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      do 40 i = 1, npb
      dky = dny*real(k-1)*part(2,i+ib)
      zt1 = qmn*cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 30 j = 2, nxh
      zt3 = sctxp(j-1,i)
      q(j,k) = q(j,k) + zt1*zt3
      q(j,k1) = q(j,k1) + zt2*zt3
   30 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctxp(nxh,i))
      q(1,k) = q(1,k) + zt1
      q(1,k1) = q(1,k1) + at3*zt2
   40 continue
   50 continue
      do 70 i = 1, npb
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i+ib)
      at1 = qmn*cos(dky)
      do 60 j = 2, nxh
      zt3 = sctxp(j-1,i)
      q(j,1) = q(j,1) + qmn*zt3
      q(j,k1) = q(j,k1) + at1*zt3
   60 continue
      at3 = real(sctxp(nxh,i))
      q(1,1) = cmplx(real(q(1,1))+qmn,aimag(q(1,1))+qmn*at3)
      q(1,k1) = cmplx(real(q(1,k1))+at1,aimag(q(1,k1))+at1*at3)
   70 continue
   80 continue
      return
      end
      subroutine PUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxvh,ny
     1v,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, periodic boundaries and various particle boundary conditions.
c scalar version using guard cells
c 28*(NX/2)*(NY/2) + 14 flops/particle, NX*NY loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c fy(x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxy = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt, ek
      complex fxy, sctx
      dimension part(idimp,nop), fxy(2,nxvh,nyv), sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtm, dx, dy
      double precision sum1, exl, eyl, ex, ey
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qtm = qbm*dt
      sum1 = 0.0d0
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
      do 50 i = 1, nop
c find electric field
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ex = 0.0d0
      ey = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
      eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      ex = ex + exl
      ey = ey + eyl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
      eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      at3 = real(sctx(nxh))
      ex = ex + (real(fxy(1,1,1)) + real(fxy(1,1,k1))*at1)
      ey = ey + (real(fxy(2,1,1)) + real(fxy(2,1,k1))*at1)
      at1 = at1*at3
      dx = ex + (aimag(fxy(1,1,1))*at3 + aimag(fxy(1,1,k1))*at1)
      dy = ey + (aimag(fxy(2,1,1))*at3 + aimag(fxy(2,1,k1))*at1)
c new velocity
      dx = part(3,i) + qtm*dx
      dy = part(4,i) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,i))**2 + (dy + part(4,i))**2
      part(3,i) = dx
      part(4,i) = dy
c new position
      dx = part(1,i) + dx*dt
      dy = part(2,i) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(3,i) = -part(3,i)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   50 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine PUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ek,idimp,nop,nx,ny,
     1nxvh,nyv,npp,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, periodic boundaries and various particle boundary conditions.
c scalar version using guard cells
c 28*(NX/2)*(NY/2) + 14 flops/particle, NX*NY loads, 4 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c fy(x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c sctxp = scratch array for sines and cosines
c exyp = scratch array for particle forces
c qbm = particle charge/mass
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxy = second dimension of field arrays, must be >= ny
c npp = number of particles to be processed together
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
c optimized version
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv, npp, ipbc
      real part, qbm, dt, ek
      double precision exyp
      complex fxy, sctxp
      dimension part(idimp,nop), fxy(2,nxvh,nyv), sctxp(nxvh,npp)
      dimension exyp(2,npp)
      integer i, j, k, l, k1, nxh, nyh, ny2, npb, ipp, ie, ib
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtm, dx, dy
      double precision sum1, exl, eyl, ex, ey
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      qtm = qbm*dt
      sum1 = 0.0d0
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
c outer loop over blocks of particles
      do 80 l = 1, ipp
      ie = l*npb
      ib = ie - npb
      if (ie.gt.nop) npb = nop - ib
      do 20 i = 1, npb
c find electric field
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i+ib)
      sctxp(j,i) = cmplx(cos(dkx),sin(dkx))
   10 continue
      exyp(1,i) = 0.0d0
      exyp(2,i) = 0.0d0
   20 continue
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      do 40 i = 1, npb
      dky = dny*real(k-1)*part(2,i+ib)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      do 30 j = 2, nxh
      zt3 = sctxp(j-1,i)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
   30 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctxp(nxh,i))
      zt2 = zt2*at3
      exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
      eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      exyp(1,i) = exyp(1,i) + exl
      exyp(2,i) = exyp(2,i) + eyl
   40 continue
   50 continue
      do 70 i = 1, npb
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i+ib)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      do 60 j = 2, nxh
      zt3 = sctxp(j-1,i)
      zt1 = at1*zt3
      exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
      eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
   60 continue
      ex = 2.0d0*(exyp(1,i) + exl)
      ey = 2.0d0*(exyp(2,i) + eyl)
      at3 = real(sctxp(nxh,i))
      ex = ex + (real(fxy(1,1,1)) + real(fxy(1,1,k1))*at1)
      ey = ey + (real(fxy(2,1,1)) + real(fxy(2,1,k1))*at1)
      at1 = at1*at3
      dx = ex + (aimag(fxy(1,1,1))*at3 + aimag(fxy(1,1,k1))*at1)
      dy = ey + (aimag(fxy(2,1,1))*at3 + aimag(fxy(2,1,k1))*at1)
c new velocity
      dx = part(3,i+ib) + qtm*dx
      dy = part(4,i+ib) + qtm*dy
c average kinetic energy
      sum1 = sum1 + (dx + part(3,i+ib))**2 + (dy + part(4,i+ib))**2
      part(3,i+ib) = dx
      part(4,i+ib) = dy
c new position
      dx = part(1,i+ib) + dx*dt
      dy = part(2,i+ib) + dy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+ib)
            part(3,i+ib) = -part(3,i+ib)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+ib)
            part(4,i+ib) = -part(4,i+ib)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+ib)
            part(3,i+ib) = -part(3,i+ib)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i+ib) = dx
      part(2,i+ib) = dy
   70 continue
   80 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
      subroutine GCJPOST2(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation.
c scalar version using guard cells
c 96 flops/particle, 40 loads, 18 stores
c input: all, output: cu
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c fy(x,y) = (.75-dy**2)*((.75-dx**2)*fy(n,m)+(.5*(.5+dx)**2)*fy(n+1,m)+
c (.5*(.5-dx)**2)*fy(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fy(n,m+1)+
c (.5*(.5+dx)**2)*fy(n+1,m+1)+(.5*(.5-dx)**2)*fy(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fy(n,m-1)+(.5*(.5+dx)**2)*fy(n+1,m-1)+
c (.5*(.5-dx)**2)*fy(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x velocity of particle n at t - dt/2
c part(4,n) = y velocity of particle n at t - dt/2
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c cu(i,j+1,k+1) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of current density array, must be >= nx+3
c nyv = second dimension of current density array, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), cu(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy
      qtmh = 0.5*qbm*dt
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
c new velocity
      vx = part(3,j) + qtmh*dx
      vy = part(4,j) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      cu(1,nl,mm) = cu(1,nl,mm) + vx*dx
      cu(2,nl,mm) = cu(2,nl,mm) + vy*dx
      dx = dxl*dyl
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      dy = amx*dyl
      cu(1,np,mm) = cu(1,np,mm) + vx*dz
      cu(2,np,mm) = cu(2,np,mm) + vy*dz
      dz = dxp*dyl
      cu(1,nl,ml) = cu(1,nl,ml) + vx*dx
      cu(2,nl,ml) = cu(2,nl,ml) + vy*dx
      dx = dxl*dyp
      cu(1,nn,ml) = cu(1,nn,ml) + vx*dy
      cu(2,nn,ml) = cu(2,nn,ml) + vy*dy
      dy = amx*dyp
      cu(1,np,ml) = cu(1,np,ml) + vx*dz
      cu(2,np,ml) = cu(2,np,ml) + vy*dz
      dz = dxp*dyp
      cu(1,nl,mp) = cu(1,nl,mp) + vx*dx
      cu(2,nl,mp) = cu(2,nl,mp) + vy*dx
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(1,np,mp) = cu(1,np,mp) + vx*dz
      cu(2,np,mp) = cu(2,np,mp) + vy*dz
   10 continue
      return
      end
      subroutine GCJPOST2L(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates particle current density
c using first-order spline interpolation.
c scalar version using guard cells
c 52 flops/particle, 20 loads, 0 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1.m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c cu(i,j,k) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of current density array, must be >= nx+1
c nyv = second dimension of current density array, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), cu(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, dxp, dyp, amx, amy, dx, dy, vx, vy
      qtmh = 0.5*qbm*dt
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
c new velocity
      vx = part(3,j) + qtmh*dx
      vy = part(4,j) + qtmh*dy
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
      dx = dxp*dyp
      dy = amx*dyp
      cu(1,np,mp) = cu(1,np,mp) + vx*dx
      cu(2,np,mp) = cu(2,np,mp) + vy*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + vx*dx
      cu(2,np,mm) = cu(2,np,mm) + vy*dx
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
   10 continue
      return
      end
      subroutine GCJPOST2C(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates particle current density
c using third-order spline interpolation.
c scalar version using guard cells
c 194 flops/particle, 68 loads, 32 stores
c input: all, output: cu
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c velocity equations at t=t+dt/2 are calculated from:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
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
c similarly for fy(x,y).
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = x velocity of particle n at t - dt/2
c part(4,n) = y velocity of particle n at t - dt/2
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c cu(i,j+2,k+2) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of current density array, must be >= nx+5
c nyv = second dimension of current density array, must be >= ny+5
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), cu(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, dxl, dyl, dxn, dyn, dxp, dyp, dxq, dyq, at1, at2
      real dx, dy, dz, dw, vx, vy
      qtmh = 0.5*qbm*dt
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
c new velocity
      vx = part(3,j) + qtmh*dx
      vy = part(4,j) + qtmh*dy
c deposit current density
      dxl = qm*dxl
      dxn = qm*dxn
      dxp = qm*dxp
      dxq = qm*dxq
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      cu(1,nl,ml) = cu(1,nl,ml) + vx*dx
      cu(2,nl,ml) = cu(2,nl,ml) + vy*dx
      dx = dxl*dyn
      cu(1,nn,ml) = cu(1,nn,ml) + vx*dy
      cu(2,nn,ml) = cu(2,nn,ml) + vy*dy
      dy = dxn*dyn
      cu(1,np,ml) = cu(1,np,ml) + vx*dz
      cu(2,np,ml) = cu(2,np,ml) + vy*dz
      dz = dxp*dyn
      cu(1,nq,ml) = cu(1,nq,ml) + vx*dw
      cu(2,nq,ml) = cu(2,nq,ml) + vy*dw
      dw = dxq*dyn
      cu(1,nl,mm) = cu(1,nl,mm) + vx*dx
      cu(2,nl,mm) = cu(2,nl,mm) + vy*dx
      dx = dxl*dyp
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      dy = dxn*dyp
      cu(1,np,mm) = cu(1,np,mm) + vx*dz
      cu(2,np,mm) = cu(2,np,mm) + vy*dz
      dz = dxp*dyp
      cu(1,nq,mm) = cu(1,nq,mm) + vx*dw
      cu(2,nq,mm) = cu(2,nq,mm) + vy*dw
      dw = dxq*dyp
      cu(1,nl,mp) = cu(1,nl,mp) + vx*dx
      cu(2,nl,mp) = cu(2,nl,mp) + vy*dx
      dx = dxl*dyq
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      dy = dxn*dyq
      cu(1,np,mp) = cu(1,np,mp) + vx*dz
      cu(2,np,mp) = cu(2,np,mp) + vy*dz
      dz = dxp*dyq
      cu(1,nq,mp) = cu(1,nq,mp) + vx*dw
      cu(2,nq,mp) = cu(2,nq,mp) + vy*dw
      dw = dxq*dyq
      cu(1,nl,mq) = cu(1,nl,mq) + vx*dx
      cu(2,nl,mq) = cu(2,nl,mq) + vy*dx
      cu(1,nn,mq) = cu(1,nn,mq) + vx*dy
      cu(2,nn,mq) = cu(2,nn,mq) + vy*dy
      cu(1,np,mq) = cu(1,np,mq) + vx*dz
      cu(2,np,mq) = cu(2,np,mq) + vy*dz
      cu(1,nq,mq) = cu(1,nq,mq) + vx*dw
      cu(2,nq,mq) = cu(2,nq,mq) + vy*dw
   10 continue
      return
      end
