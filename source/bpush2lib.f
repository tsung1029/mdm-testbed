c-----------------------------------------------------------------------
c 2d PIC library for pushing particles with magnetic field
c and depositing current
c bpush2lib.f contains procedures to process particles with magnetic
c             fields:
c GJPOST2 deposits 3 component current density, quadratic interpolation,
c         STANDARD optimization.
c GSJPOST2 deposits 3 component current density, quadratic
c          interpolation, LOOKAHEAD optimization.
c GSJPOST2X deposits 3 component current density, quadratic
c           interpolation, VECTOR optimization.
c GJPOST2L deposits 3 component current density, linear interpolation,
c          STANDARD optimization.
c GSJPOST2L deposits 3 component current density, linear interpolation,
c           LOOKAHEAD optimization.
c GSJPOST2XL deposits 3 component current density, linear interpolation,
c            VECTOR optimization.
c GJPOST22 deposits 2 component current density, quadratic
c          interpolation, STANDARD optimization.
c GSJPOST22 deposits 2 component current density, quadratic
c           interpolation, LOOKAHEAD optimization.
c GSJPOST22X deposits 2 component current density, quadratic
c            interpolation, VECTOR optimization.
c GJPOST22L deposits 2 component current density, linear interpolation,
c           STANDARD optimization.
c GSJPOST22L deposits 2 component current density, linear interpolation,
c            LOOKAHEAD optimization.
c GSJPOST22XL deposits 2 component current density, linear interpolation,
c             VECTOR optimization.
c GJPOST2C deposits 3 component current density, cubic interpolation,
c          STANDARD optimization.
c GJPOST22C deposits 2 component current density, cubic interpolation,
c           STANDARD optimization.
c GBPUSH2 push particles with 3 component magnetic field, 2 component
c         electric field, quadratic interpolation, STANDARD
c         optimization.
c GSBPUSH2 push particles with 3 component magnetic field, 2 component
c          electric field, quadratic interpolation, LOOKAHEAD
c          optimization.
c GBPUSH2L push particles with 3 component magnetic field, 2 component
c          electric field, linear interpolation, STANDARD optimization.
c GSBPUSH2L push particles with 3 component magnetic field, 2 component
c           electric field, linear interpolation, LOOKAHEAD optimization
c GBPUSH23 push particles with 3 component magnetic field, 3 component
c          electric field, quadratic interpolation, STANDARD
c          optimization.
c GSBPUSH23 push particles with 3 component magnetic field, 3 component
c           electric field, quadratic interpolation, LOOKAHEAD
c           optimization.
c GBPUSH23L push particles with 3 component magnetic field, 3 component
c           electric field, linear interpolation, STANDARD optimization.
c GSBPUSH23L push particles with 3 component magnetic field, 3 component
c            electric field, linear interpolation, LOOKAHEAD
c            optimization.
c GBPUSH22 push particles with 1 component magnetic field, 2 component
c          electric field, quadratic interpolation, STANDARD
c          optimization.
c GSBPUSH22 push particles with 1 component magnetic field, 2 component
c           electric field, quadratic interpolation, LOOKAHEAD
c           optimization.
c GBPUSH22L push particles with 1 component magnetic field, 2 component
c           electric field, linear interpolation, STANDARD optimization.
c GSBPUSH22L push particles with 1 component magnetic field, 2 component
c            electric field, linear interpolation, LOOKAHEAD
c            optimization.
c GBPUSH2C push particles with 3 component magnetic field, 2 component
c          electric field, cubic interpolation, STANDARD
c          optimization.
c GBPUSH23C push particles with 3 component magnetic field, 3 component
c           electric field, cubic interpolation, STANDARD
c           optimization.
c GBPUSH22C push particles with 1 component magnetic field, 2 component
c           electric field, cubic interpolation, STANDARD optimization.
c RETARD2 retard particle position a half time-step for 2-1/2d code.
c DJPOST2GL deposits 3 component current density, using gridless method.
c BPUSH23GL push particles with 3 component magnetic field, 3 component
c           electric field, using gridless method.
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: august 11, 2009
c-----------------------------------------------------------------------
      subroutine DJPOST2(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 78 flops/particle, 32 loads, 29 stores
c input: all, output: part, cux, cuy, cuz
c current density is approximated by values at the nearest grid points
c cui(n,m)=qci*(.75-dx**2)*(.75-dy**2)
c cui(n+1,m)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cui(n-1,m)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cui(n,m+1)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cui(n+1,m+1)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cui(n-1,m+1)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cui(n,m-1)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cui(n+1,m-1)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cui(n-1,m-1)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cui(j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx
      dimension part(idimp,nop), cux(nxv,ny), cuy(nxv,ny), cuz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cux(nl,mm) = cux(nl,mm) + vx*dx
      cux(nn,mm) = cux(nn,mm) + vx*dy
      cux(np,mm) = cux(np,mm) + vx*dz
      cuy(nl,mm) = cuy(nl,mm) + vy*dx
      cuy(nn,mm) = cuy(nn,mm) + vy*dy
      cuy(np,mm) = cuy(np,mm) + vy*dz
      cuz(nl,mm) = cuz(nl,mm) + vz*dx
      cuz(nn,mm) = cuz(nn,mm) + vz*dy
      cuz(np,mm) = cuz(np,mm) + vz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cux(nl,ml) = cux(nl,ml) + vx*dx
      cux(nn,ml) = cux(nn,ml) + vx*dy
      cux(np,ml) = cux(np,ml) + vx*dz
      cuy(nl,ml) = cuy(nl,ml) + vy*dx
      cuy(nn,ml) = cuy(nn,ml) + vy*dy
      cuy(np,ml) = cuy(np,ml) + vy*dz
      cuz(nl,ml) = cuz(nl,ml) + vz*dx
      cuz(nn,ml) = cuz(nn,ml) + vz*dy
      cuz(np,ml) = cuz(np,ml) + vz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cux(nl,mp) = cux(nl,mp) + vx*dx
      cux(nn,mp) = cux(nn,mp) + vx*dy
      cux(np,mp) = cux(np,mp) + vx*dz
      cuy(nl,mp) = cuy(nl,mp) + vy*dx
      cuy(nn,mp) = cuy(nn,mp) + vy*dy
      cuy(np,mp) = cuy(np,mp) + vy*dz
      cuz(nl,mp) = cuz(nl,mp) + vz*dx
      cuz(nn,mp) = cuz(nn,mp) + vz*dy
      cuz(np,mp) = cuz(np,mp) + vz*dz
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j) = dy
   10 continue
      return
      end
      subroutine GJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 78 flops/particle, 32 loads, 29 stores
c input: all, output: part, cu
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
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+3
c nyv = second dimension of current array, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(3,nxv,nyv)
      qmh = .5*qm
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cu(1,nl,mm) = cu(1,nl,mm) + vx*dx
      cu(2,nl,mm) = cu(2,nl,mm) + vy*dx
      cu(3,nl,mm) = cu(3,nl,mm) + vz*dx
      dx = dxl*dyl
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + vz*dy
      dy = amx*dyl
      cu(1,np,mm) = cu(1,np,mm) + vx*dz
      cu(2,np,mm) = cu(2,np,mm) + vy*dz
      cu(3,np,mm) = cu(3,np,mm) + vz*dz
      dz = dxp*dyl
      cu(1,nl,ml) = cu(1,nl,ml) + vx*dx
      cu(2,nl,ml) = cu(2,nl,ml) + vy*dx
      cu(3,nl,ml) = cu(3,nl,ml) + vz*dx
      dx = dxl*dyp
      cu(1,nn,ml) = cu(1,nn,ml) + vx*dy
      cu(2,nn,ml) = cu(2,nn,ml) + vy*dy
      cu(3,nn,ml) = cu(3,nn,ml) + vz*dy
      dy = amx*dyp
      cu(1,np,ml) = cu(1,np,ml) + vx*dz
      cu(2,np,ml) = cu(2,np,ml) + vy*dz
      cu(3,np,ml) = cu(3,np,ml) + vz*dz
      dz = dxp*dyp
      cu(1,nl,mp) = cu(1,nl,mp) + vx*dx
      cu(2,nl,mp) = cu(2,nl,mp) + vy*dx
      cu(3,nl,mp) = cu(3,nl,mp) + vz*dx
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + vz*dy
      cu(1,np,mp) = cu(1,np,mp) + vx*dz
      cu(2,np,mp) = cu(2,np,mp) + vy*dz
      cu(3,np,mp) = cu(3,np,mp) + vz*dz
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
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
      return
      end
      subroutine GSJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 78 flops/particle, 32 loads, 29 stores
c input: all, output: part, cu
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
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyv = dimension of current array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(3,nxyv)
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      qmh = .5*qm
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j-1)
      vy = part(4,j-1)
      vz = part(5,j-1)
      dx1 = cu(1,mn) + vx*dx
      dy1 = cu(2,mn) + vy*dx
      amy = cu(3,mn) + vz*dx
      dx2 = cu(1,mn+1) + vx*dy
      dy2 = cu(2,mn+1) + vy*dy
      dx3 = cu(3,mn+1) + vz*dy
      dx = cu(1,mn+2) + vx*dz
      dy = cu(2,mn+2) + vy*dz
      dz = cu(3,mn+2) + vz*dz
      cu(1,mn) = dx1
      cu(2,mn) = dy1
      cu(3,mn) = amy
      cu(1,mn+1) = dx2
      cu(2,mn+1) = dy2
      cu(3,mn+1) = dx3
      cu(1,mn+2) = dx
      cu(2,mn+2) = dy
      cu(3,mn+2) = dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = cu(1,ml) + vx*dx
      dy1 = cu(2,ml) + vy*dx
      amy = cu(3,ml) + vz*dx
      dx2 = cu(1,ml+1) + vx*dy
      dy2 = cu(2,ml+1) + vy*dy
      dyl = cu(3,ml+1) + vz*dy
      dx = cu(1,ml+2) + vx*dz
      dy = cu(2,ml+2) + vy*dz
      dz = cu(3,ml+2) + vz*dz
      cu(1,ml) = dx1
      cu(2,ml) = dy1
      cu(3,ml) = amy
      cu(1,ml+1) = dx2
      cu(2,ml+1) = dy2
      cu(3,ml+1) = dyl
      cu(1,ml+2) = dx
      cu(2,ml+2) = dy
      cu(3,ml+2) = dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = cu(1,mp) + vx*dx
      dy1 = cu(2,mp) + vy*dx
      amy = cu(3,mp) + vz*dx
      dxl = cu(1,mp+1) + vx*dy
      amx = cu(2,mp+1) + vy*dy
      dxp = cu(3,mp+1) + vz*dy
      dx = cu(1,mp+2) + vx*dz
      dy = cu(2,mp+2) + vy*dz
      dz = cu(3,mp+2) + vz*dz
      cu(1,mp) = dx1
      cu(2,mp) = dy1
      cu(3,mp) = amy
      cu(1,mp+1) = dxl
      cu(2,mp+1) = amx
      cu(3,mp+1) = dxp
      cu(1,mp+2) = dx
      cu(2,mp+2) = dy
      cu(3,mp+2) = dz
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
c deposit current for last particle
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      cu(1,mn) = cu(1,mn) + vx*dx
      cu(2,mn) = cu(2,mn) + vy*dx
      cu(3,mn) = cu(3,mn) + vz*dx
      cu(1,mn+1) = cu(1,mn+1) + vx*dy
      cu(2,mn+1) = cu(2,mn+1) + vy*dy
      cu(3,mn+1) = cu(3,mn+1) + vz*dy
      cu(1,mn+2) = cu(1,mn+2) + vx*dz
      cu(2,mn+2) = cu(2,mn+2) + vy*dz
      cu(3,mn+2) = cu(3,mn+2) + vz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cu(1,ml) = cu(1,ml) + vx*dx
      cu(2,ml) = cu(2,ml) + vy*dx
      cu(3,ml) = cu(3,ml) + vz*dx
      cu(1,ml+1) = cu(1,ml+1) + vx*dy
      cu(2,ml+1) = cu(2,ml+1) + vy*dy
      cu(3,ml+1) = cu(3,ml+1) + vz*dy
      cu(1,ml+2) = cu(1,ml+2) + vx*dz
      cu(2,ml+2) = cu(2,ml+2) + vy*dz
      cu(3,ml+2) = cu(3,ml+2) + vz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cu(1,mp) = cu(1,mp) + vx*dx
      cu(2,mp) = cu(2,mp) + vy*dx
      cu(3,mp) = cu(3,mp) + vz*dx
      cu(1,mp+1) = cu(1,mp+1) + vx*dy
      cu(2,mp+1) = cu(2,mp+1) + vy*dy
      cu(3,mp+1) = cu(3,mp+1) + vz*dy
      cu(1,mp+2) = cu(1,mp+2) + vx*dz
      cu(2,mp+2) = cu(2,mp+2) + vy*dz
      cu(3,mp+2) = cu(3,mp+2) + vz*dz
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
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
      return
      end
      subroutine SJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with 1d addressing
c 78 flops/particle, 87 loads, 83 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(n,m,i)=qci*(.75-dx**2)*(.75-dy**2)
c cu(n+1,m,i)=.5*qci*((.5+dx)**2)*(.75-dy**2)
c cu(n-1,m,i)=.5*qci*((.5-dx)**2)*(.75-dy**2)
c cu(n,m+1,i)=.5*qci*(.75-dx**2)*(.5+dy)**2
c cu(n+1,m+1,i)=.25*qci*((.5+dx)**2)*(.5+dy)**2
c cu(n-1,m+1,i)=.25*qci*((.5-dx)**2)*(.5+dy)**2
c cu(n,m-1,i)=.5*qci*(.75-dx**2)*(.5-dy)**2
c cu(n+1,m-1,i)=.25*qci*((.5+dx)**2)*(.5-dy)**2
c cu,n-1,m-1,i)=.25*qci*((.5-dx)**2)*(.5-dy)**2
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(n,i) = ith component of current density at grid point j,k
c where n = j + nxv*(k - 1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx
c nxvy = nxv*ny, dimension of current array
      parameter(npp=512)
      dimension part(idimp,nop), cu(3*nxvy)
      dimension nn(27,npp), amxy(27,npp)
      nxvy2 = nxvy + nxvy
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      zero = 0.
      anx = float(nx)
      any = float(ny)
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
      in = n + m
      nn(1,j) = in
      nn(1+9,j) = in + nxvy
      nn(1+18,j) = in + nxvy2
      in = np + m
      nn(2,j) = in
      nn(2+9,j) = in + nxvy
      nn(2+18,j) = in + nxvy2
      in = nl + m
      nn(3,j) = in
      nn(3+9,j) = in + nxvy
      nn(3+18,j) = in + nxvy2
      in = n + mp
      nn(4,j) = in
      nn(4+9,j) = in + nxvy
      nn(4+18,j) = in + nxvy2
      in = np + mp
      nn(5,j) = in
      nn(5+9,j) = in + nxvy
      nn(5+18,j) = in + nxvy2
      in = nl + mp
      nn(6,j) = in
      nn(6+9,j) = in + nxvy
      nn(6+18,j) = in + nxvy2
      in = n + ml
      nn(7,j) = in
      nn(7+9,j) = in + nxvy
      nn(7+18,j) = in + nxvy2
      in = np + ml
      nn(8,j) = in
      nn(8+9,j) = in + nxvy
      nn(8+18,j) = in + nxvy2
      in = nl + ml
      nn(9,j) = in
      nn(9+9,j) = in + nxvy
      nn(9+18,j) = in + nxvy2
      dx = amx*amy
      dy = dxp*amy
      dz = dxl*amy
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      vz = part(5,j+jb)
      amxy(1,j) = vx*dx
      amxy(2,j) = vx*dy
      amxy(3,j) = vx*dz
      amxy(1+9,j) = vy*dx
      amxy(2+9,j) = vy*dy
      amxy(3+9,j) = vy*dz
      amxy(1+18,j) = vz*dx
      amxy(2+18,j) = vz*dy
      amxy(3+18,j) = vz*dz
      dx = amx*dyp
      dy = dxp*dyp
      dz = dxl*dyp
      amxy(4,j) = vx*dx
      amxy(5,j) = vx*dy
      amxy(6,j) = vx*dz
      amxy(4+9,j) = vy*dx
      amxy(5+9,j) = vy*dy
      amxy(6+9,j) = vy*dz
      amxy(4+18,j) = vz*dx
      amxy(5+18,j) = vz*dy
      amxy(6+18,j) = vz*dz
      dx = amx*dyl
      dy = dxp*dyl
      dz = dxl*dyl
      amxy(7,j) = vx*dx
      amxy(8,j) = vx*dy
      amxy(9,j) = vx*dz
      amxy(7+9,j) = vy*dx
      amxy(8+9,j) = vy*dy
      amxy(9+9,j) = vy*dz
      amxy(7+18,j) = vz*dx
      amxy(8+18,j) = vz*dy
      amxy(9+18,j) = vz*dz
c advance position half a time-step
      dx = part(1,j+jb) + vx*dt
      dy = part(2,j+jb) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j+jb) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j+jb) = dy
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 27
      cu(nn(i,j)) = cu(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 78 flops/particle, 87 loads, 83 stores
c input: all, output: part, cu
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
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyv = dimension of current array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      parameter(npp=512)
      dimension part(idimp,nop), cu(3*nxyv)
      dimension nn(27,npp), amxy(27,npp)
      nxv3 = 3*nxv
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      qmh = .5*qm
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
      n3 = 3*n + 1
      m = nxv3*m
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = m + n3
      dxl = .5*qm*(.5 - dxp)**2
      dxp = .5*qm*(.5 + dxp)**2
      mn = ml + nxv3
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv3
      nn(1,j) = mn
      nn(2,j) = mn + 1
      nn(3,j) = mn + 2
      nn(4,j) = mn + 3
      nn(5,j) = mn + 4
      nn(6,j) = mn + 5
      nn(7,j) = mn + 6
      nn(8,j) = mn + 7
      nn(9,j) = mn + 8
      nn(10,j) = ml
      nn(11,j) = ml + 1
      nn(12,j) = ml + 2
      nn(13,j) = ml + 3
      nn(14,j) = ml + 4
      nn(15,j) = ml + 5
      nn(16,j) = ml + 6
      nn(17,j) = ml + 7
      nn(18,j) = ml + 8
      nn(19,j) = mp
      nn(20,j) = mp + 1
      nn(21,j) = mp + 2
      nn(22,j) = mp + 3
      nn(23,j) = mp + 4
      nn(24,j) = mp + 5
      nn(25,j) = mp + 6
      nn(26,j) = mp + 7
      nn(27,j) = mp + 8
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      vz = part(5,j+jb)
      amxy(1,j) = vx*dx
      amxy(2,j) = vy*dx
      amxy(3,j) = vz*dx
      dx = dxl*dyl
      amxy(4,j) = vx*dy
      amxy(5,j) = vy*dy
      amxy(6,j) = vz*dy
      dy = amx*dyl
      amxy(7,j) = vx*dz
      amxy(8,j) = vy*dz
      amxy(9,j) = vz*dz
      dz = dxp*dyl
      amxy(10,j) = vx*dx
      amxy(11,j) = vy*dx
      amxy(12,j) = vz*dx
      dx = dxl*dyp
      amxy(13,j) = vx*dy
      amxy(14,j) = vy*dy
      amxy(15,j) = vz*dy
      dy = amx*dyp
      amxy(16,j) = vx*dz
      amxy(17,j) = vy*dz
      amxy(18,j) = vz*dz
      dz = dxp*dyp
      amxy(19,j) = vx*dx
      amxy(20,j) = vy*dx
      amxy(21,j) = vz*dx
      amxy(22,j) = vx*dy
      amxy(23,j) = vy*dy
      amxy(24,j) = vz*dy
      amxy(25,j) = vx*dz
      amxy(26,j) = vy*dz
      amxy(27,j) = vz*dz
c advance position half a time-step
      dx = part(1,j+jb) + vx*dt
      dy = part(2,j+jb) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j+jb)
            part(4,j+jb) = -part(4,j+jb)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j+jb) = dx
      part(2,j+jb) = dy
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 27
      cu(nn(i,j)) = cu(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine DJPOST2L(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 37 flops/particle, 17 loads, 14 stores
c input: all, output: part, cux, cuy, cuz
c current density is approximated by values at the nearest grid points
c cui(n,m)=qci*(1.-dx)*(1.-dy)
c cu(n+1,m)=qci*dx*(1.-dy)
c cu(n,m+1)=qci*(1.-dx)*dy
c cu(n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cui(j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx
      dimension part(idimp,nop), cux(nxv,ny), cuy(nxv,ny), cuz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
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
c deposit current
      dx = amx*amy
      dy = dxp*amy
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cux(nn,mm) = cux(nn,mm) + vx*dx
      cux(np,mm) = cux(np,mm) + vx*dy
      cuy(nn,mm) = cuy(nn,mm) + vy*dx
      cuy(np,mm) = cuy(np,mm) + vy*dy
      cuz(nn,mm) = cuz(nn,mm) + vz*dx
      cuz(np,mm) = cuz(np,mm) + vz*dy
      dx = amx*dyp
      dy = dxp*dyp
      cux(nn,mp) = cux(nn,mp) + vx*dx
      cux(np,mp) = cux(np,mp) + vx*dy
      cuy(nn,mp) = cuy(nn,mp) + vy*dx
      cuy(np,mp) = cuy(np,mp) + vy*dy
      cuz(nn,mp) = cuz(nn,mp) + vz*dx
      cuz(np,mp) = cuz(np,mp) + vz*dy
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j) = dy
   10 continue
      return
      end
      subroutine GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 37 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+1
c nyv = second dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(3,nxv,nyv)
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
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cu(1,np,mp) = cu(1,np,mp) + vx*dx
      cu(2,np,mp) = cu(2,np,mp) + vy*dx
      cu(3,np,mp) = cu(3,np,mp) + vz*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + vz*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + vx*dx
      cu(2,np,mm) = cu(2,np,mm) + vy*dx
      cu(3,np,mm) = cu(3,np,mm) + vz*dx
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + vz*dy
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
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
      return
      end
      subroutine GSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 37 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(3,nxyv)
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
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
c deposit current
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1)
      vy = part(4,j-1)
      vz = part(5,j-1)
      dx1 = cu(1,mp+1) + vx*dx
      dy1 = cu(2,mp+1) + vy*dx
      dyp = cu(3,mp+1) + vz*dx
      dx = cu(1,mp) + vx*dz
      dy = cu(2,mp) + vy*dz
      dz = cu(3,mp) + vz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx
      cu(2,mp) = dy
      cu(3,mp) = dz
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + vx*dx
      amx = cu(2,mm+1) + vy*dx
      dyp = cu(3,mm+1) + vz*dx
      dx = cu(1,mm) + vx*dz
      dy = cu(2,mm) + vy*dz
      dz = cu(3,mm) + vz*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(3,mm+1) = dyp
      cu(1,mm) = dx
      cu(2,mm) = dy
      cu(3,mm) = dz
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop)
      vy = part(4,nop)
      vz = part(5,nop)
      cu(1,mp+1) = cu(1,mp+1) + vx*dx
      cu(2,mp+1) = cu(2,mp+1) + vy*dx
      cu(3,mp+1) = cu(3,mp+1) + vz*dx
      cu(1,mp) = cu(1,mp) + vx*dy
      cu(2,mp) = cu(2,mp) + vy*dy
      cu(3,mp) = cu(3,mp) + vz*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + vx*dx
      cu(2,mm+1) = cu(2,mm+1) + vy*dx
      cu(3,mm+1) = cu(3,mm+1) + vz*dx
      cu(1,mm) = cu(1,mm) + vx*dy
      cu(2,mm) = cu(2,mm) + vy*dy
      cu(3,mm) = cu(3,mm) + vz*dy
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
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
      return
      end
      subroutine SJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with 1d addressing
c 37 flops/particle, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(n,m,i)=qci*(1.-dx)*(1.-dy)
c cu(n+1,m,i)=qci*dx*(1.-dy)
c cu(n,m+1,i)=qci*(1.-dx)*dy
c cu(n+1,m+1,i)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(n,i) = ith component of current density at grid point j,k
c where n = j + nxv*(k - 1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx
c nxvy = nxv*ny, dimension of current array
      parameter(npp=128)
      dimension part(idimp,nop), cu(3*nxvy)
      dimension nn(12,npp), amxy(12,npp)
      nxvy2 = nxvy + nxvy
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      zero = 0.
      anx = float(nx)
      any = float(ny)
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
      in = n + m
      nn(1,j) = in
      nn(1+4,j) = in + nxvy
      nn(1+8,j) = in + nxvy2
      in = np + m
      nn(2,j) = in
      nn(2+4,j) = in + nxvy
      nn(2+8,j) = in + nxvy2
      in = n + mp
      nn(3,j) = in
      nn(3+4,j) = in + nxvy
      nn(3+8,j) = in + nxvy2
      in = np + mp
      nn(4,j) = in
      nn(4+4,j) = in + nxvy
      nn(4+8,j) = in + nxvy2
      dx = amx*amy
      dy = dxp*amy
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      vz = part(5,j+jb)
      amxy(1,j) = vx*dx
      amxy(2,j) = vx*dy
      amxy(1+4,j) = vy*dx
      amxy(2+4,j) = vy*dy
      amxy(1+8,j) = vz*dx
      amxy(2+8,j) = vz*dy
      dx = amx*dyp
      dy = dxp*dyp
      amxy(3,j) = vx*dx
      amxy(4,j) = vx*dy
      amxy(3+4,j) = vy*dx
      amxy(4+4,j) = vy*dy
      amxy(3+8,j) = vz*dx
      amxy(4+8,j) = vz*dy
c advance position half a time-step
      dx = part(1,j+jb) + vx*dt
      dy = part(2,j+jb) + vy*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j+jb) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j+jb) = dy
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 12
      cu(nn(i,j)) = cu(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 37 flops/particle, 41 loads, 33 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      parameter(npp=128)
      dimension part(idimp,nop), cu(3*nxyv)
      dimension nn(12,npp), amxy(12,npp)
      nxv3 = 3*nxv
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
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
      n3 = 3*n + 1
      m = nxv3*m
      amx = qm - dxp
      m = m + n3
      amy = 1. - dyp
      mp = m + nxv3
      nn(10,j) = m
      nn(11,j) = m + 1
      nn(12,j) = m + 2
      nn(7,j) = m + 3
      nn(8,j) = m + 4
      nn(9,j) = m + 5
      nn(4,j) = mp
      nn(5,j) = mp + 1
      nn(6,j) = mp + 2
      nn(1,j) = mp + 3
      nn(2,j) = mp + 4
      nn(3,j) = mp + 5
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      vz = part(5,j+jb)
      amxy(1,j) = vx*dx
      amxy(2,j) = vy*dx
      amxy(3,j) = vz*dx
      amxy(4,j) = vx*dy
      amxy(5,j) = vy*dy
      amxy(6,j) = vz*dy
      dx = dxp*amy
      dy = amx*amy
      amxy(7,j) = vx*dx
      amxy(8,j) = vy*dx
      amxy(9,j) = vz*dx
      amxy(10,j) = vx*dy
      amxy(11,j) = vy*dy
      amxy(12,j) = vz*dy
c advance position half a time-step
      dx = part(1,j+jb) + vx*dt
      dy = part(2,j+jb) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j+jb)
            part(4,j+jb) = -part(4,j+jb)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j+jb) = dx
      part(2,j+jb) = dy
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 12
      cu(nn(i,j)) = cu(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 60 flops/particle, 22 loads, 20 stores
c input: all, output: part, cu
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
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+3
c nyv = second dimension of current array, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(2,nxv,nyv)
      qmh = .5*qm
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j)
      vy = part(4,j)
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
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
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
      return
      end
      subroutine GSJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 60 flops/particle, 22 loads, 20 stores
c input: all, output: part, cu
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
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyv = dimension of current array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(2,nxyv)
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      qmh = .5*qm
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j-1)
      vy = part(4,j-1)
      dx1 = cu(1,mn) + vx*dx
      dy1 = cu(2,mn) + vy*dx
      dx2 = cu(1,mn+1) + vx*dy
      dy2 = cu(2,mn+1) + vy*dy
      dx = cu(1,mn+2) + vx*dz
      dy = cu(2,mn+2) + vy*dz
      cu(1,mn) = dx1
      cu(2,mn) = dy1
      cu(1,mn+1) = dx2
      cu(2,mn+1) = dy2
      cu(1,mn+2) = dx
      cu(2,mn+2) = dy
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = cu(1,ml) + vx*dx
      dy1 = cu(2,ml) + vy*dx
      dx2 = cu(1,ml+1) + vx*dy
      dy2 = cu(2,ml+1) + vy*dy
      dx = cu(1,ml+2) + vx*dz
      dy = cu(2,ml+2) + vy*dz
      cu(1,ml) = dx1
      cu(2,ml) = dy1
      cu(1,ml+1) = dx2
      cu(2,ml+1) = dy2
      cu(1,ml+2) = dx
      cu(2,ml+2) = dy
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = cu(1,mp) + vx*dx
      dy1 = cu(2,mp) + vy*dx
      dxl = cu(1,mp+1) + vx*dy
      amx = cu(2,mp+1) + vy*dy
      dx = cu(1,mp+2) + vx*dz
      dy = cu(2,mp+2) + vy*dz
      cu(1,mp) = dx1
      cu(2,mp) = dy1
      cu(1,mp+1) = dxl
      cu(2,mp+1) = amx
      cu(1,mp+2) = dx
      cu(2,mp+2) = dy
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
c deposit current for last particle
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,nop)
      vy = part(4,nop)
      cu(1,mn) = cu(1,mn) + vx*dx
      cu(2,mn) = cu(2,mn) + vy*dx
      cu(1,mn+1) = cu(1,mn+1) + vx*dy
      cu(2,mn+1) = cu(2,mn+1) + vy*dy
      cu(1,mn+2) = cu(1,mn+2) + vx*dz
      cu(2,mn+2) = cu(2,mn+2) + vy*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cu(1,ml) = cu(1,ml) + vx*dx
      cu(2,ml) = cu(2,ml) + vy*dx
      cu(1,ml+1) = cu(1,ml+1) + vx*dy
      cu(2,ml+1) = cu(2,ml+1) + vy*dy
      cu(1,ml+2) = cu(1,ml+2) + vx*dz
      cu(2,ml+2) = cu(2,ml+2) + vy*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cu(1,mp) = cu(1,mp) + vx*dx
      cu(2,mp) = cu(2,mp) + vy*dx
      cu(1,mp+1) = cu(1,mp+1) + vx*dy
      cu(2,mp+1) = cu(2,mp+1) + vy*dy
      cu(1,mp+2) = cu(1,mp+2) + vx*dz
      cu(2,mp+2) = cu(2,mp+2) + vy*dz
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
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
      return
      end
      subroutine GSJPOST22X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 60 flops/particle, 58 loads, 56 stores
c input: all, output: part, cu
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
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+3
c nxyv = dimension of current array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      parameter(npp=512)
      dimension part(idimp,nop), cu(2*nxyv)
      dimension nn(18,npp), amxy(18,npp)
      nxv2 = 2*nxv
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
      qmh = .5*qm
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
      n2 = 2*n + 1
      m = nxv2*m
      amx = qm*(.75 - dxp*dxp)
      amy = .75 - dyp*dyp
      ml = m + n2
      dxl = .5*qm*(.5 - dxp)**2
      dxp = .5*qm*(.5 + dxp)**2
      mn = ml + nxv2
      dyl = .5*(.5 - dyp)**2
      dyp = .5*(.5 + dyp)**2
      mp = mn + nxv2
      nn(1,j) = mn
      nn(2,j) = mn + 1
      nn(3,j) = mn + 2
      nn(4,j) = mn + 3
      nn(5,j) = mn + 4
      nn(6,j) = mn + 5
      nn(7,j) = ml
      nn(8,j) = ml + 1
      nn(9,j) = ml + 2
      nn(10,j) = ml + 3
      nn(11,j) = ml + 4
      nn(12,j) = ml + 5
      nn(13,j) = mp
      nn(14,j) = mp + 1
      nn(15,j) = mp + 2
      nn(16,j) = mp + 3
      nn(17,j) = mp + 4
      nn(18,j) = mp + 5
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      amxy(1,j) = vx*dx
      amxy(2,j) = vy*dx
      dx = dxl*dyl
      amxy(3,j) = vx*dy
      amxy(4,j) = vy*dy
      dy = amx*dyl
      amxy(5,j) = vx*dz
      amxy(6,j) = vy*dz
      dz = dxp*dyl
      amxy(7,j) = vx*dx
      amxy(8,j) = vy*dx
      dx = dxl*dyp
      amxy(9,j) = vx*dy
      amxy(10,j) = vy*dy
      dy = amx*dyp
      amxy(11,j) = vx*dz
      amxy(12,j) = vy*dz
      dz = dxp*dyp
      amxy(13,j) = vx*dx
      amxy(14,j) = vy*dx
      amxy(15,j) = vx*dy
      amxy(16,j) = vy*dy
      amxy(17,j) = vx*dz
      amxy(18,j) = vy*dz
c advance position half a time-step
      dx = part(1,j+jb) + vx*dt
      dy = part(2,j+jb) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j+jb)
            part(4,j+jb) = -part(4,j+jb)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j+jb) = dx
      part(2,j+jb) = dy
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 18
      cu(nn(i,j)) = cu(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 29 flops/particle, 12 loads, 10 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+1
c nyv = second dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(2,nxv,nyv)
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
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j)
      vy = part(4,j)
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
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
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
      return
      end
      subroutine GSJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 29 flops/particle, 12 loads, 10 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(2,nxyv)
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
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
c deposit current
      dx = dxp*dyp
      dz = amx*dyp
      vx = part(3,j-1)
      vy = part(4,j-1)
      dx1 = cu(1,mp+1) + vx*dx
      dy1 = cu(2,mp+1) + vy*dx
      dx = cu(1,mp) + vx*dz
      dy = cu(2,mp) + vy*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(1,mp) = dx
      cu(2,mp) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + vx*dx
      amx = cu(2,mm+1) + vy*dx
      dx = cu(1,mm) + vx*dz
      dy = cu(2,mm) + vy*dz
      cu(1,mm+1) = dxp
      cu(2,mm+1) = amx
      cu(1,mm) = dx
      cu(2,mm) = dy
c advance position half a time-step
      dx = part(1,j-1) + vx*dt
      dy = part(2,j-1) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1)
            part(4,j-1) = -part(4,j-1)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1)
            part(3,j-1) = -part(3,j-1)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j-1) = dx
      part(2,j-1) = dy
   10 continue
c deposit current for last particle
      nn = nnn + 1
      mm = nxv*mmn
      dxp = qm*dxn
      mm = mm + nn
      amx = qm - dxp
      mp = mm + nxv
      amy = 1. - dyn
c deposit current
      dx = dxp*dyn
      dy = amx*dyn
      vx = part(3,nop)
      vy = part(4,nop)
      cu(1,mp+1) = cu(1,mp+1) + vx*dx
      cu(2,mp+1) = cu(2,mp+1) + vy*dx
      cu(1,mp) = cu(1,mp) + vx*dy
      cu(2,mp) = cu(2,mp) + vy*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + vx*dx
      cu(2,mm+1) = cu(2,mm+1) + vy*dx
      cu(1,mm) = cu(1,mm) + vx*dy
      cu(2,mm) = cu(2,mm) + vy*dy
c advance position half a time-step
      dx = part(1,nop) + vx*dt
      dy = part(2,nop) + vy*dt
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
      return
      end
      subroutine GSJPOST22XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc
     1)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 29 flops/particle, 29 loads, 26 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first virtual dimension of current array, must be >= nx+1
c nxyv = dimension of current array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      parameter(npp=128)
      dimension part(idimp,nop), cu(2*nxyv)
      dimension nn(8,npp), amxy(8,npp)
      nxv2 = 2*nxv
      npb = npp
      if (npb.gt.nop) npb = nop
      ipp = float(nop - 1)/float(npb) + 1.
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
      n2 = 2*n + 1
      m = nxv2*m
      amx = qm - dxp
      m = m + n2
      amy = 1. - dyp
      mp = m + nxv2
      nn(7,j) = m
      nn(8,j) = m + 1
      nn(5,j) = m + 2
      nn(6,j) = m + 3
      nn(3,j) = mp
      nn(4,j) = mp + 1
      nn(1,j) = mp + 2
      nn(2,j) = mp + 3
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      amxy(1,j) = vx*dx
      amxy(2,j) = vy*dx
      amxy(3,j) = vx*dy
      amxy(4,j) = vy*dy
      dx = dxp*amy
      dy = amx*amy
      amxy(5,j) = vx*dx
      amxy(6,j) = vy*dx
      amxy(7,j) = vx*dy
      amxy(8,j) = vy*dy
c advance position half a time-step
      dx = part(1,j+jb) + vx*dt
      dy = part(2,j+jb) + vy*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j+jb)
            part(4,j+jb) = -part(4,j+jb)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j+jb)
            part(3,j+jb) = -part(3,j+jb)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,j+jb) = dx
      part(2,j+jb) = dy
   10 continue
c deposit current
      do 30 j = 1, npb
cdir$ ivdep
      do 20 i = 1, 8
      cu(nn(i,j)) = cu(nn(i,j)) + amxy(i,j)
   20 continue
   30 continue
   40 continue
      return
      end
      subroutine GJPOST2C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using third-order spline interpolation, periodic boundaries
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 152 flops/particle, 53 loads, 50 stores
c input: all, output: part, cu
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
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of charge array, must be >= nx+5
c nxy = second dimension of charge array, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(3,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2
      real edgelx, edgely, edgerx, edgery, dx, dy, dz, dw, vx, vy, vz
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
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
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cu(1,nl,ml) = cu(1,nl,ml) + vx*dx
      cu(2,nl,ml) = cu(2,nl,ml) + vy*dx
      cu(3,nl,ml) = cu(3,nl,ml) + vz*dx
      dx = dxl*dyn
      cu(1,nn,ml) = cu(1,nn,ml) + vx*dy
      cu(2,nn,ml) = cu(2,nn,ml) + vy*dy
      cu(3,nn,ml) = cu(3,nn,ml) + vz*dy
      dy = dxn*dyn
      cu(1,np,ml) = cu(1,np,ml) + vx*dz
      cu(2,np,ml) = cu(2,np,ml) + vy*dz
      cu(3,np,ml) = cu(3,np,ml) + vz*dz
      dz = dxp*dyn
      cu(1,nq,ml) = cu(1,nq,ml) + vx*dw
      cu(2,nq,ml) = cu(2,nq,ml) + vy*dw
      cu(3,nq,ml) = cu(3,nq,ml) + vz*dw
      dw = dxq*dyn
      cu(1,nl,mm) = cu(1,nl,mm) + vx*dx
      cu(2,nl,mm) = cu(2,nl,mm) + vy*dx
      cu(3,nl,mm) = cu(3,nl,mm) + vz*dx
      dx = dxl*dyp
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + vz*dy
      dy = dxn*dyp
      cu(1,np,mm) = cu(1,np,mm) + vx*dz
      cu(2,np,mm) = cu(2,np,mm) + vy*dz
      cu(3,np,mm) = cu(3,np,mm) + vz*dz
      dz = dxp*dyp
      cu(1,nq,mm) = cu(1,nq,mm) + vx*dw
      cu(2,nq,mm) = cu(2,nq,mm) + vy*dw
      cu(3,nq,mm) = cu(3,nq,mm) + vz*dw
      dw = dxq*dyp
      cu(1,nl,mp) = cu(1,nl,mp) + vx*dx
      cu(2,nl,mp) = cu(2,nl,mp) + vy*dx
      cu(3,nl,mp) = cu(3,nl,mp) + vz*dx
      dx = dxl*dyq
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + vz*dy
      dy = dxn*dyq
      cu(1,np,mp) = cu(1,np,mp) + vx*dz
      cu(2,np,mp) = cu(2,np,mp) + vy*dz
      cu(3,np,mp) = cu(3,np,mp) + vz*dz
      dz = dxp*dyq
      cu(1,nq,mp) = cu(1,nq,mp) + vx*dw
      cu(2,nq,mp) = cu(2,nq,mp) + vy*dw
      cu(3,nq,mp) = cu(3,nq,mp) + vz*dw
      dw = dxq*dyq
      cu(1,nl,mq) = cu(1,nl,mq) + vx*dx
      cu(2,nl,mq) = cu(2,nl,mq) + vy*dx
      cu(3,nl,mq) = cu(3,nl,mq) + vz*dx
      cu(1,nn,mq) = cu(1,nn,mq) + vx*dy
      cu(2,nn,mq) = cu(2,nn,mq) + vy*dy
      cu(3,nn,mq) = cu(3,nn,mq) + vz*dy
      cu(1,np,mq) = cu(1,np,mq) + vx*dz
      cu(2,np,mq) = cu(2,np,mq) + vy*dz
      cu(3,np,mq) = cu(3,np,mq) + vz*dz
      cu(1,nq,mq) = cu(1,nq,mq) + vx*dw
      cu(2,nq,mq) = cu(2,nq,mq) + vy*dw
      cu(3,nq,mq) = cu(3,nq,mq) + vz*dw
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
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
      return
      end
      subroutine GJPOST22C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2d code, this subroutine calculates particle current density
c using third-order spline interpolation, periodic boundaries
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 120 flops/particle, 36 loads, 34 stores
c input: all, output: part, cu
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
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*vi, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of charge array, must be >= nx+5
c nxy = second dimension of charge array, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real part, cu, qm, dt
      dimension part(idimp,nop), cu(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2
      real edgelx, edgely, edgerx, edgery, dx, dy, dz, dw, vx, vy
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
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
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      vx = part(3,j)
      vy = part(4,j)
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
c advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
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
      return
      end
      subroutine BPUSH2(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny,nx
     1v)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover.
c baseline scalar version
c 188 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bx(j,k) = x component of magnetic field at grid (j,k)
c by(j,k) = y component of magnetic field at grid (j,k)
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      double precision sum1
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fx(nl,mm) + amx*fx(nn,mm) + dxp*fx(np,mm)) + dyl*(dx
     1l*fx(nl,ml) + amx*fx(nn,ml) + dxp*fx(np,ml)) + dyp*(dxl*fx(nl,mp)
     2+ amx*fx(nn,mp) + dxp*fx(np,mp))
      dy = amy*(dxl*fy(nl,mm) + amx*fy(nn,mm) + dxp*fy(np,mm)) + dyl*(dx
     1l*fy(nl,ml) + amx*fy(nn,ml) + dxp*fy(np,ml)) + dyp*(dxl*fy(nl,mp)
     2+ amx*fy(nn,mp) + dxp*fy(np,mp))
c find magnetic field
      ox = amy*(dxl*bx(nl,mm) + amx*bx(nn,mm) + dxp*bx(np,mm)) + dyl*(dx
     1l*bx(nl,ml) + amx*bx(nn,ml) + dxp*bx(np,ml)) + dyp*(dxl*bx(nl,mp)
     2+ amx*bx(nn,mp) + dxp*bx(np,mp))
      oy = amy*(dxl*by(nl,mm) + amx*by(nn,mm) + dxp*by(np,mm)) + dyl*(dx
     1l*by(nl,ml) + amx*by(nn,ml) + dxp*by(np,ml)) + dyp*(dxl*by(nl,mp)
     2+ amx*by(nn,mp) + dxp*by(np,mp))
      oz = amy*(dxl*bz(nl,mm) + amx*bz(nn,mm) + dxp*bz(np,mm)) + dyl*(dx
     1l*bz(nl,ml) + amx*bz(nn,ml) + dxp*bz(np,ml)) + dyp*(dxl*bz(nl,mp)
     2+ amx*bz(nn,mp) + dxp*bz(np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
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
      ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv,
     1nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 188 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp)) 
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GSBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 188 flops/particle, 1 divide, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bxy(3,nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
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
   20 ek = ek + .5*sum1
      return
      end
      subroutine BPUSH2L(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover.
c baseline scalar version
c 107 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bx(j,k) = x component of magnetic field at grid (j,k)
c by(j,k) = y component of magnetic field at grid (j,k)
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      double precision sum1
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(amx*fx(nn,mm) + dxp*fx(np,mm)) + dyp*(amx*fx(nn,mp) + dx
     1p*fx(np,mp))
      dy = amy*(amx*fy(nn,mm) + dxp*fy(np,mm)) + dyp*(amx*fy(nn,mp) + dx
     1p*fy(np,mp))
c find magnetic field
      ox = amy*(amx*bx(nn,mm) + dxp*bx(np,mm)) + dyp*(amx*bx(nn,mp) + dx
     1p*bx(np,mp))
      oy = amy*(amx*by(nn,mm) + dxp*by(np,mm)) + dyp*(amx*by(nn,mp) + dx
     1p*by(np,mp))
      oz = amy*(amx*bz(nn,mm) + dxp*bz(np,mm)) + dyp*(amx*bz(nn,mp) + dx
     1p*bz(np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
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
      ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 107 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp)) + amy*(dxp*bxy(1,np
     1,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp)) + amy*(dxp*bxy(2,np
     1,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp)) + amy*(dxp*bxy(3,np
     1,mm) + amx*bxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GSBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 107 flops/particle, 1 divide, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bxy(3,nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
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
   20 ek = ek + .5*sum1
      return
      end
      subroutine BPUSH2CQ(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field and periodic boundary
c conditions.  Using the Boris Mover with correction.
c baseline scalar version
c 204 flops/particle, 2 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bx(j,k) = x component of magnetic field at grid (j,k)
c by(j,k) = y component of magnetic field at grid (j,k)
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      double precision sum1
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fx(nl,mm) + amx*fx(nn,mm) + dxp*fx(np,mm)) + dyl*(dx
     1l*fx(nl,ml) + amx*fx(nn,ml) + dxp*fx(np,ml)) + dyp*(dxl*fx(nl,mp)
     2+ amx*fx(nn,mp) + dxp*fx(np,mp))
      dy = amy*(dxl*fy(nl,mm) + amx*fy(nn,mm) + dxp*fy(np,mm)) + dyl*(dx
     1l*fy(nl,ml) + amx*fy(nn,ml) + dxp*fy(np,ml)) + dyp*(dxl*fy(nl,mp)
     2+ amx*fy(nn,mp) + dxp*fy(np,mp))
c find magnetic field
      ox = amy*(dxl*bx(nl,mm) + amx*bx(nn,mm) + dxp*bx(np,mm)) + dyl*(dx
     1l*bx(nl,ml) + amx*bx(nn,ml) + dxp*bx(np,ml)) + dyp*(dxl*bx(nl,mp)
     2+ amx*bx(nn,mp) + dxp*bx(np,mp))
      oy = amy*(dxl*by(nl,mm) + amx*by(nn,mm) + dxp*by(np,mm)) + dyl*(dx
     1l*by(nl,ml) + amx*by(nn,ml) + dxp*by(np,ml)) + dyp*(dxl*by(nl,mp)
     2+ amx*by(nn,mp) + dxp*by(np,mp))
      oz = amy*(dxl*bz(nl,mm) + amx*bz(nn,mm) + dxp*bz(np,mm)) + dyl*(dx
     1l*bz(nl,ml) + amx*bz(nn,ml) + dxp*bz(np,ml)) + dyp*(dxl*bz(nl,mp)
     2+ amx*bz(nn,mp) + dxp*bz(np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j) + (dx*dtt + omxt*omt)
      dy = part(2,j) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH2CQ(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ny
     1v,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field.  Using the Boris Mover
c with correction.
c scalar version using guard cells
c 204 flops/particle, 2 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp)) 
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j) + (dx*dtt + omxt*omt)
      dy = part(2,j) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            omt = .5*omzt*dt
            dx = part(1,j)
            dy = dy + omt*part(3,j)
            if ((dy.lt.edgely).or.(dy.ge.edgery)) then
               dx = dx - omt*part(4,j)
               dy = part(2,j) + omt*part(3,j)
               part(4,j) = -part(4,j)
            endif
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            omt = .5*omzt*dt
            dy = part(2,j)
            dx = dx - omt*part(4,j)
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j) - omt*part(4,j)
               dy = dy + omt*part(3,j)
               part(3,j) = -part(3,j)
            endif
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            dy = dy + .5*omzt*dt*part(3,j)
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
      subroutine BPUSH2CL(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover with correction.
c baseline scalar version
c 123 flops/particle, 2 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bx(j,k) = x component of magnetic field at grid (j,k)
c by(j,k) = y component of magnetic field at grid (j,k)
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx
      double precision sum1
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(amx*fx(nn,mm) + dxp*fx(np,mm)) + dyp*(amx*fx(nn,mp) + dx
     1p*fx(np,mp))
      dy = amy*(amx*fy(nn,mm) + dxp*fy(np,mm)) + dyp*(amx*fy(nn,mp) + dx
     1p*fy(np,mp))
c find magnetic field
      ox = amy*(amx*bx(nn,mm) + dxp*bx(np,mm)) + dyp*(amx*bx(nn,mp) + dx
     1p*bx(np,mp))
      oy = amy*(amx*by(nn,mm) + dxp*by(np,mm)) + dyp*(amx*by(nn,mp) + dx
     1p*by(np,mp))
      oz = amy*(amx*bz(nn,mm) + dxp*bz(np,mm)) + dyp*(amx*bz(nn,mp) + dx
     1p*bz(np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j) + (dx*dtt + omxt*omt)
      dy = part(2,j) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j) = dx
      if (dy.lt.zero) dy = dy + any
      if (dy.ge.any) dy = dy - any
      part(2,j) = dy
   10 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH2CL(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ny
     1v,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover
c with correction.
c scalar version using guard cells
c 123 flops/particle, 2 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c where vplx = omx*vpl/om and vply = omy*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp)) + amy*(dxp*bxy(1,np
     1,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp)) + amy*(dxp*bxy(2,np
     1,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp)) + amy*(dxp*bxy(3,np
     1,mm) + amx*bxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j) + (dx*dtt + omxt*omt)
      dy = part(2,j) + (dy*dtt + omyt*omt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            omt = .5*omzt*dt
            dx = part(1,j)
            dy = dy + omt*part(3,j)
            if ((dy.lt.edgely).or.(dy.ge.edgery)) then
               dx = dx - omt*part(4,j)
               dy = part(2,j) + omt*part(3,j)
               part(4,j) = -part(4,j)
            endif
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            omt = .5*omzt*dt
            dy = part(2,j)
            dx = dx - omt*part(4,j)
            if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
               dx = part(1,j) - omt*part(4,j)
               dy = dy + omt*part(3,j)
               part(3,j) = -part(3,j)
            endif
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            dy = dy + .5*omzt*dt*part(3,j)
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
      subroutine GBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 211 flops/particle, 1 divide, 59 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c fxy(3,j+1,k+1) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      qtmh = .5*qbm*dt
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GSBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 211 flops/particle, 1 divide, 59 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c fxy(3,j+1,k+1) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+1,k+1) = x component of magnetic field at grid (j,k)
c bxy(2,j+1,k+1) = y component of magnetic field at grid (j,k)
c bxy(3,j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
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
   20 ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp)) + amy*(dxp*fxy(3,np
     1,mm) + amx*fxy(3,nn,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp)) + amy*(dxp*bxy(1,np
     1,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp)) + amy*(dxp*bxy(2,np
     1,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp)) + amy*(dxp*bxy(3,np
     1,mm) + amx*bxy(3,nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,n
     1xv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 119 flops/particle, 1 divide, 29 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c fxy(3,j,k) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j,k) = x component of magnetic field at grid (j,k)
c bxy(2,j,k) = y component of magnetic field at grid (j,k)
c bxy(3,j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyp*(dxp*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxp*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyp*(dxp*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxp*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyp*(dxp*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxp*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find magnetic field
      ox = dyn*(dxn*bxy(1,mp+1) + amx*bxy(1,mp)) + amy*(dxn*bxy(1,mm+1)
     1+ amx*bxy(1,mm))
      oy = dyn*(dxn*bxy(2,mp+1) + amx*bxy(2,mp)) + amy*(dxn*bxy(2,mm+1)
     1+ amx*bxy(2,mm))
      oz = dyn*(dxn*bxy(3,mp+1) + amx*bxy(3,mp)) + amy*(dxn*bxy(3,mm+1)
     1+ amx*bxy(3,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,nop) = dx
      part(4,nop) = dy
      part(5,nop) = dz
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
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
   20 ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv,
     1nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field in z direction. 
c Using the Boris Mover.
c scalar version using guard cells
c 143 flops/particle, 1 divide, 31 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp)) 
c find magnetic field
      oz = amy*(dxl*bz(nl,mm) + amx*bz(nn,mm) + dxp*bz(np,mm)) + dyl*(dx
     1l*bz(nl,ml) + amx*bz(nn,ml) + dxp*bz(np,ml)) + dyp*(dxl*bz(nl,mp) 
     2+ amx*bz(nn,mp) + dxp*bz(np,mp))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GSBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nxyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field in z direction. 
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 143 flops/particle, 1 divide, 31 loads, 4 stores
c input: all, output: part, ek
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (.75-dy**2)*((.75-dx**2)*fx(n,m)+(.5*(.5+dx)**2)*fx(n+1,m)+
c (.5*(.5-dx)**2)*fx(n-1,m)) + (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1)+(.5*(.5-dx)**2)*fx(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1)+(.5*(.5+dx)**2)*fx(n+1,m-1)+
c (.5*(.5-dx)**2)*fx(n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
c find magnetic field
      oz = amy*(dxl*bz(mn) + amx*bz(mn+1) + dxp*bz(mn+2)) + dyl*(dxl*bz(
     1ml) + amx*bz(ml+1) + dxp*bz(ml+2)) + dyp*(dxl*bz(mp) + amx*bz(mp+1
     2) + dxp*bz(mp+2))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
c find electric field
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
c find magnetic field
      oz = amy*(dxl*bz(mn) + amx*bz(mn+1) + dxp*bz(mn+2)) + dyl*(dxl*bz(
     1ml) + amx*bz(ml+1) + dxp*bz(ml+2)) + dyp*(dxl*bz(mp)+ amx*bz(mp+1)
     2 + dxp*bz(mp+2)) 
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,nop) = dx
      part(4,nop) = dy
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
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
   20 ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order spline
c interpolation in space, with magnetic field in z direction. 
c Using the Boris Mover.
c scalar version using guard cells
c 64 flops/particle, 1 divide, 16 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c find magnetic field
      oz = dyp*(dxp*bz(np,mp) + amx*bz(nn,mp)) + amy*(dxp*bz(np,mm) + am
     1x*bz(nn,mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GSBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nxyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order spline
c interpolation in space, with magnetic field in z direction. 
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 64 flops/particle, 1 divide, 16 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      double precision sum1
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find electric field
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find magnetic field
      oz = dyp*(dxp*bz(mp+1) + amx*bz(mp)) + amy*(dxp*bz(mm+1) + amx*bz
     1(mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
c find electric field
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find magnetic field
      oz = dyn*(dxn*bz(mp+1) + amx*bz(mp)) + amy*(dxn*bz(mm+1) + amx*bz
     1(mm))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,nop) = dx
      part(4,nop) = dy
c new position
      dx = part(1,nop) + dx*dtc
      dy = part(2,nop) + dy*dtc
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
   20 ek = ek + .5*sum1
      return
      end
      subroutine GBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 273 flops/particle, 1 divide, 85 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2))
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)),
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
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+2,k+2) = x component of magnetic field at grid (j,k)
c bxy(2,j+2,k+2) = y component of magnetic field at grid (j,k)
c bxy(3,j+2,k+2) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nyv = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, bxy, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, edgelx, edgely, edgerx, edgery
      real dx, dy, dz, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq, at1, at2
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order linear
c interpolation in space, with magnetic field. Using the Boris Mover.
c scalar version using guard cells
c 311 flops/particle, 1 divide, 101 loads, 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c fxy(3,j+2,k+2) = z component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bxy(1,j+2,k+2) = x component of magnetic field at grid (j,k)
c bxy(2,j+2,k+2) = y component of magnetic field at grid (j,k)
c bxy(3,j+2,k+2) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nyv = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, bxy, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, edgelx, edgely, edgerx, edgery
      real dx, dy, dz, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq, at1, at2
      real ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = .5*qbm*dt
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine GBPUSH22C(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order spline
c interpolation in space, with magnetic field in z direction. 
c Using the Boris Mover.
c scalar version using guard cells
c 168 flops/particle, 1 divide, 52 loads, 4 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = -rot(2)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(1)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (omz*dt/2)**2/(1 + (omz*dt/2)**2)
c    rot(2) = 2*(omz*dt/2)/(1 + (omz*dt/2)**2)
c where omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
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
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(+2j,k+2) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nyv = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, bz, qbm, dt, dtc, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, edgelx, edgely, edgerx, edgery
      real dx, dy, oz, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq, at1, at2
      real acx, acy, omzt, omt, anorm, rot1, rot2
      double precision sum1
      qtmh = .5*qbm*dt
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
c find magnetic field
      oz = dyl*(dxl*bz(nl,ml) + dxn*bz(nn,ml) +  dxp*bz(np,ml) + dxq*bz(
     1nq,ml)) + dyn*(dxl*bz(nl,mm) + dxn*bz(nn,mm) + dxp*bz(np,mm) +  dx
     2q*bz(nq,mm)) + dyp*(dxl*bz(nl,mp) + dxn*bz(nn,mp) +  dxp*bz(np,mp)
     3+  dxq*bz(nq,mp)) + dyq*(dxl*bz(nl,mq) + dxn*bz(nn,mq) +  dxp*bz(n
     4p,mq) +  dxq*bz(nq,mq))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy)
c calculate cyclotron frequency
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omzt*omzt
      anorm = 2./(1. + omt)
      rot1 = .5*(1. - omt)
      rot2 = omzt
c new velocity
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot1*acy - rot2*acx)*anorm + dy
      part(3,j) = dx
      part(4,j) = dy
c new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
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
      subroutine RETARD2(part,dtc,idimp,nop,nx,ny,ipbc)
c for 2-1/2d code, particle positions are retarded a half time-step
c input: all, output: part
c equations used are:
c x(t+dt) = x(t) - vx(t+dt/2)*dtc, y(t+dt) = y(t) - vy(t+dt/2)*dtc,
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c dtc = time interval between successive co-ordinate calculations
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
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
c retard position half a time-step for current deposit
      dx = part(1,j) - part(3,j)*dtc
      dy = part(2,j) - part(4,j)*dtc
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
      return
      end
      subroutine DJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxvh,nyv,i
     1pbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using gridless spectral version, periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 36*(NX/2)*(NY/2) + 7 flops/particle, 3*NX*NY/2 loads and stores
c plus (NX/2 + NY/2) sines and cosines/particle
c input: all, output: part, cu
c charge density is calculated from the expression:
c cu(i,n,m) = sum(qmi*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                exp(-sqrt(-1)*2*m*pi*y/ny))
c where qmi = qm*vi, where i = x,y,z
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x velocity of particle n
c part(4,n) = y velocity of particle n
c part(5,n) = z velocity of particle n
c cu(i,n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nx/2
c nyv = second dimension of current array, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc
      real part, qm, dt
      complex cu, sctx
      dimension part(idimp,nop), cu(3,nxvh,nyv), sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3
      real vx, vy, vz, edgelx, edgely, edgerx, edgery, dx, dy
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qmn = qm/real(nx*ny)
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
c find fourier components
      do 50 i = 1, nop
      vx = qmn*part(3,i)
      vy = qmn*part(4,i)
      vz = qmn*part(5,i)
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
      cu(1,j,k) = cu(1,j,k) + vx*zt3
      cu(2,j,k) = cu(2,j,k) + vy*zt3
      cu(3,j,k) = cu(3,j,k) + vz*zt3
      cu(1,j,k1) = cu(1,j,k1) + vx*zt4
      cu(2,j,k1) = cu(2,j,k1) + vy*zt4
      cu(3,j,k1) = cu(3,j,k1) + vz*zt4
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = at3*zt2
      cu(1,1,k) = cu(1,1,k) + vx*zt1
      cu(2,1,k) = cu(2,1,k) + vy*zt1
      cu(3,1,k) = cu(3,1,k) + vz*zt1
      cu(1,1,k1) = cu(1,1,k1) + vx*zt2
      cu(2,1,k1) = cu(2,1,k1) + vy*zt2
      cu(3,1,k1) = cu(3,1,k1) + vz*zt2
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      cu(1,j,1) = cu(1,j,1) + vx*zt3
      cu(2,j,1) = cu(2,j,1) + vy*zt3
      cu(3,j,1) = cu(3,j,1) + vz*zt3
      cu(1,j,k1) = cu(1,j,k1) + vx*zt1
      cu(2,j,k1) = cu(2,j,k1) + vy*zt1
      cu(3,j,k1) = cu(3,j,k1) + vz*zt1
   40 continue
      at3 = real(sctx(nxh))
      cu(1,1,1) = cmplx(real(cu(1,1,1))+vx,aimag(cu(1,1,1))+vx*at3)
      cu(2,1,1) = cmplx(real(cu(2,1,1))+vy,aimag(cu(2,1,1))+vy*at3)
      cu(3,1,1) = cmplx(real(cu(3,1,1))+vz,aimag(cu(3,1,1))+vz*at3)
      at3 = at1*at3
      cu(1,1,k1) = cmplx(real(cu(1,1,k1))+vx*at1,aimag(cu(1,1,k1))+vx*at
     13)
      cu(2,1,k1) = cmplx(real(cu(2,1,k1))+vy*at1,aimag(cu(2,1,k1))+vy*at
     13)
      cu(3,1,k1) = cmplx(real(cu(3,1,k1))+vz*at1,aimag(cu(3,1,k1))+vz*at
     13)
c advance position half a time-step
      dx = part(1,i) + part(3,i)*dt
      dy = part(2,i) + part(4,i)*dt
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
      return
      end
      subroutine BPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,nx,
     1ny,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, with magnetic field, periodic boundaries and various particle
c boundary conditions.  Using the Boris Mover.
c scalar version using guard cells
c 60*(NX/2)*(NY/2) + 57 flops/particle, 1 divide, 3*NX*NY + 5 loads,
c 5 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
c omz = (q/m)*bz(x(t),y(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c bxy(i,j,k) = i component of magnetic field at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qtmh = .5*qbm*dt
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
      acx = part(3,i) + dx
      acy = part(4,i) + dy
      acz = part(5,i) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
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
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,i) = dx
      part(4,i) = dy
      part(5,i) = dz
c new position
      dx = part(1,i) + dx*dtc
      dy = part(2,i) + dy*dtc
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
      ek = ek + .5*sum1
      return
      end
