c 2d PIC library for pushing particles with darwin vector potential
c using hamiltonian formulation
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: november 18, 2009
      subroutine CR8PCM23(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv)
c for 2-1/2d code, this subroutine calculates the canonical momentum of
c a particle from the vector potential at t - dt/2, and adds this
c information to the particle array, using quadratic interpolation
c canonical momentum p(t-dt/2) = v(t-dt/2) + qi*A(t-dt/2)/(c*mi)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c on entry
c part(6,n) = position x of particle n at t - dt/2
c part(7,n) = position y of particle n at t - dt/2
c on exit
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c axy(1,j+1,k+1) = x component of vector potential at grid (j,k)
c axy(2,j+1,k+1) = y component of vector potential at grid (j,k)
c axy(3,j+1,k+1) = z component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c sx/sy/sz = x/y/z components of field momentum
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field array, must be >= nx+3
c nyv = second dimension of field array, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, axy, qbm, sx, sy, sz
      dimension part(idimp,nop)
      dimension axy(3,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dxp, dyp, amx, amy, dxl, dyl, ax, ay, az
      double precision sum1, sum2, sum3
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 10 j = 1, nop
c find interpolation weights
      nn = part(6,j) + .5
      mm = part(7,j) + .5
      dxp = part(6,j) - float(nn)
      dyp = part(7,j) - float(mm)
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
c find vector potential at t - dt/2
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
      az = amy*(dxl*axy(3,nl,mm) + amx*axy(3,nn,mm) + dxp*axy(3,np,mm))
     1+ dyl*(dxl*axy(3,nl,ml) + amx*axy(3,nn,ml) + dxp*axy(3,np,ml)) + d
     2yp*(dxl*axy(3,nl,mp) + amx*axy(3,nn,mp) + dxp*axy(3,np,mp))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      sum1 = sum1 + ax
      sum2 = sum2 + ay
      sum3 = sum3 + az
      part(6,j) = part(3,j) + ax
      part(7,j) = part(4,j) + ay
      part(8,j) = part(5,j) + az
   10 continue
      sx = sum1
      sy = sum2
      sz = sum3
      return
      end
      subroutine GHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv,n
     1yv)
c for 2-1/2d code, this subroutine calculates current density, using
c second-order spline interpolation.
c scalar version using guard cells
c 380 flops/particle, 1 divide, 131 loads, 27 stores
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j+1,k+1) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension cu(3,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl, dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,nl,mm) + amx*daxy(1,nn,mm) + dxp*daxy(1,np,m
     1m)) + dyl*(dxl*daxy(1,nl,ml) + amx*daxy(1,nn,ml) + dxp*daxy(1,np,m
     2l)) + dyp*(dxl*daxy(1,nl,mp) + amx*daxy(1,nn,mp) + dxp*daxy(1,np,m
     3p))
      oyx = amy*(dxl*daxy(2,nl,mm) + amx*daxy(2,nn,mm) + dxp*daxy(2,np,m
     1m)) + dyl*(dxl*daxy(2,nl,ml) + amx*daxy(2,nn,ml) + dxp*daxy(2,np,m
     2l)) + dyp*(dxl*daxy(2,nl,mp) + amx*daxy(2,nn,mp) + dxp*daxy(2,np,m
     3p))  
      oxy = amy*(dxl*daxy(3,nl,mm) + amx*daxy(3,nn,mm) + dxp*daxy(3,np,m
     1m)) + dyl*(dxl*daxy(3,nl,ml) + amx*daxy(3,nn,ml) + dxp*daxy(3,np,m
     2l)) + dyp*(dxl*daxy(3,nl,mp) + amx*daxy(3,nn,mp) + dxp*daxy(3,np,m
     3p))
      ozx = amy*(dxl*daxy(4,nl,mm) + amx*daxy(4,nn,mm) + dxp*daxy(4,np,m
     1m)) + dyl*(dxl*daxy(4,nl,ml) + amx*daxy(4,nn,ml) + dxp*daxy(4,np,m
     2l)) + dyp*(dxl*daxy(4,nl,mp) + amx*daxy(4,nn,mp) + dxp*daxy(4,np,m
     3p))  
      ozy = amy*(dxl*daxy(5,nl,mm) + amx*daxy(5,nn,mm) + dxp*daxy(5,np,m
     1m)) + dyl*(dxl*daxy(5,nl,ml) + amx*daxy(5,nn,ml) + dxp*daxy(5,np,m
     2l)) + dyp*(dxl*daxy(5,nl,mp) + amx*daxy(5,nn,mp) + dxp*daxy(5,np,m
     3p))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))
      dz = amy*(dxl*fxy(3,nl,mm) + amx*fxy(3,nn,mm) + dxp*fxy(3,np,mm))
     1+ dyl*(dxl*fxy(3,nl,ml) + amx*fxy(3,nn,ml) + dxp*fxy(3,np,ml)) + d
     2yp*(dxl*fxy(3,nl,mp) + amx*fxy(3,nn,mp) + dxp*fxy(3,np,mp)) 
c find vector potential at t
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
      az = amy*(dxl*axy(3,nl,mm) + amx*axy(3,nn,mm) + dxp*axy(3,np,mm))
     1+ dyl*(dxl*axy(3,nl,ml) + amx*axy(3,nn,ml) + dxp*axy(3,np,ml)) + d
     2yp*(dxl*axy(3,nl,mp) + amx*axy(3,nn,mp) + dxp*axy(3,np,mp))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      cu(1,nl,mm) = cu(1,nl,mm) + px*dx
      cu(2,nl,mm) = cu(2,nl,mm) + py*dx
      cu(3,nl,mm) = cu(3,nl,mm) + pz*dx
      dx = dxl*dyl
      cu(1,nn,mm) = cu(1,nn,mm) + px*dy
      cu(2,nn,mm) = cu(2,nn,mm) + py*dy
      cu(3,nn,mm) = cu(3,nn,mm) + pz*dy
      dy = amx*dyl
      cu(1,np,mm) = cu(1,np,mm) + px*dz
      cu(2,np,mm) = cu(2,np,mm) + py*dz
      cu(3,np,mm) = cu(3,np,mm) + pz*dz
      dz = dxp*dyl
      cu(1,nl,ml) = cu(1,nl,ml) + px*dx
      cu(2,nl,ml) = cu(2,nl,ml) + py*dx
      cu(3,nl,ml) = cu(3,nl,ml) + pz*dx
      dx = dxl*dyp
      cu(1,nn,ml) = cu(1,nn,ml) + px*dy
      cu(2,nn,ml) = cu(2,nn,ml) + py*dy
      cu(3,nn,ml) = cu(3,nn,ml) + pz*dy
      dy = amx*dyp
      cu(1,np,ml) = cu(1,np,ml) + px*dz
      cu(2,np,ml) = cu(2,np,ml) + py*dz
      cu(3,np,ml) = cu(3,np,ml) + pz*dz
      dz = dxp*dyp
      cu(1,nl,mp) = cu(1,nl,mp) + px*dx
      cu(2,nl,mp) = cu(2,nl,mp) + py*dx
      cu(3,nl,mp) = cu(3,nl,mp) + pz*dx
      cu(1,nn,mp) = cu(1,nn,mp) + px*dy
      cu(2,nn,mp) = cu(2,nn,mp) + py*dy
      cu(3,nn,mp) = cu(3,nn,mp) + pz*dy
      cu(1,np,mp) = cu(1,np,mp) + px*dz
      cu(2,np,mp) = cu(2,np,mp) + py*dz
      cu(3,np,mp) = cu(3,np,mp) + pz*dz
   10 continue
      return
      end
      subroutine GSHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv,
     1nxyv)
c for 2-1/2d code, this subroutine calculates current density, using
c second-order spline interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 380 flops/particle, 1 divide, 131 loads, 27 stores
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
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j+1,k+1) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension cu(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm, dx1, dy1, dx2, dy2, dx3, dy3, dx4, dy4
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2)) 
      ozx = amy*(dxl*daxy(4,mn) + amx*daxy(4,mn+1) + dxp*daxy(4,mn+2)) +
     1 dyl*(dxl*daxy(4,ml) + amx*daxy(4,ml+1) + dxp*daxy(4,ml+2)) + dyp*
     2(dxl*daxy(4,mp) + amx*daxy(4,mp+1) + dxp*daxy(4,mp+2)) 
      ozy = amy*(dxl*daxy(5,mn) + amx*daxy(5,mn+1) + dxp*daxy(5,mn+2)) +
     1 dyl*(dxl*daxy(5,ml) + amx*daxy(5,ml+1) + dxp*daxy(5,ml+2)) + dyp*
     2(dxl*daxy(5,mp) + amx*daxy(5,mp+1) + dxp*daxy(5,mp+2)) 
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
      dz = amy*(dxl*fxy(3,mn) + amx*fxy(3,mn+1) + dxp*fxy(3,mn+2)) + dyl
     1*(dxl*fxy(3,ml) + amx*fxy(3,ml+1) + dxp*fxy(3,ml+2)) + dyp*(dxl*fx
     2y(3,mp) + amx*fxy(3,mp+1) + dxp*fxy(3,mp+2))
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      dx1 = cu(1,mn) + px*dx
      dy1 = cu(2,mn) + py*dx
      amy = cu(3,mn) + pz*dx
      dx2 = cu(1,mn+1) + px*dy
      dy2 = cu(2,mn+1) + py*dy
      dx3 = cu(3,mn+1) + pz*dy
      dy3 = cu(1,mn+2) + px*dz
      dx4 = cu(2,mn+2) + py*dz
      dy4 = cu(3,mn+2) + pz*dz
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
      dx1 = cu(1,ml) + px*dx
      dy1 = cu(2,ml) + py*dx
      amy = cu(3,ml) + pz*dx
      dx2 = cu(1,ml+1) + px*dy
      dy2 = cu(2,ml+1) + py*dy
      dx3 = cu(3,ml+1) + pz*dy
      dy3 = cu(1,ml+2) + px*dz
      dx4 = cu(2,ml+2) + py*dz
      dy4 = cu(3,ml+2) + pz*dz
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
      dx1 = cu(1,mp) + px*dx
      dy1 = cu(2,mp) + py*dx
      amy = cu(3,mp) + pz*dx
      dx2 = cu(1,mp+1) + px*dy
      dy2 = cu(2,mp+1) + py*dy
      dx3 = cu(3,mp+1) + pz*dy
      dy3 = cu(1,mp+2) + px*dz
      dx4 = cu(2,mp+2) + py*dz
      dy4 = cu(3,mp+2) + pz*dz
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2)) 
      ozx = amy*(dxl*daxy(4,mn) + amx*daxy(4,mn+1) + dxp*daxy(4,mn+2)) +
     1 dyl*(dxl*daxy(4,ml) + amx*daxy(4,ml+1) + dxp*daxy(4,ml+2)) + dyp*
     2(dxl*daxy(4,mp) + amx*daxy(4,mp+1) + dxp*daxy(4,mp+2)) 
      ozy = amy*(dxl*daxy(5,mn) + amx*daxy(5,mn+1) + dxp*daxy(5,mn+2)) +
     1 dyl*(dxl*daxy(5,ml) + amx*daxy(5,ml+1) + dxp*daxy(5,ml+2)) + dyp*
     2(dxl*daxy(5,mp) + amx*daxy(5,mp+1) + dxp*daxy(5,mp+2)) 
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
      dz = amy*(dxl*fxy(3,mn) + amx*fxy(3,mn+1) + dxp*fxy(3,mn+2)) + dyl
     1*(dxl*fxy(3,ml) + amx*fxy(3,ml+1) + dxp*fxy(3,ml+2)) + dyp*(dxl*fx
     2y(3,mp) + amx*fxy(3,mp+1) + dxp*fxy(3,mp+2))
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,nop)
      py = part(7,nop)
      pz = part(8,nop)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      cu(1,mn) = cu(1,mn) + px*dx
      cu(2,mn) = cu(2,mn) + py*dx
      cu(3,mn) = cu(3,mn) + pz*dx
      cu(1,mn+1) = cu(1,mn+1) + px*dy
      cu(2,mn+1) = cu(2,mn+1) + py*dy
      cu(3,mn+1) = cu(3,mn+1) + pz*dy
      cu(1,mn+2) = cu(1,mn+2) + px*dz
      cu(2,mn+2) = cu(2,mn+2) + py*dz
      cu(3,mn+2) = cu(3,mn+2) + pz*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cu(1,ml) = cu(1,ml) + px*dx
      cu(2,ml) = cu(2,ml) + py*dx
      cu(3,ml) = cu(3,ml) + pz*dx
      cu(1,ml+1) = cu(1,ml+1) + px*dy
      cu(2,ml+1) = cu(2,ml+1) + py*dy
      cu(3,ml+1) = cu(3,ml+1) + pz*dy
      cu(1,ml+2) = cu(1,ml+2) + px*dz
      cu(2,ml+2) = cu(2,ml+2) + py*dz
      cu(3,ml+2) = cu(3,ml+2) + pz*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cu(1,mp) = cu(1,mp) + px*dx
      cu(2,mp) = cu(2,mp) + py*dx
      cu(3,mp) = cu(3,mp) + pz*dx
      cu(1,mp+1) = cu(1,mp+1) + px*dy
      cu(2,mp+1) = cu(2,mp+1) + py*dy
      cu(3,mp+1) = cu(3,mp+1) + pz*dy
      cu(1,mp+2) = cu(1,mp+2) + px*dz
      cu(2,mp+2) = cu(2,mp+2) + py*dz
      cu(3,mp+2) = cu(3,mp+2) + pz*dz
      return
      end
      subroutine GHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c second-order spline interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 314 flops/particle, 1 divide, 104 loads, 8 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fx(1,n,m)+(.5*(.5+dx)**2)*fx(1,n+1,m)+
c (.5*(.5-dx)**2)*fx(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fx(1,n,m+1)+
c (.5*(.5+dx)**2)*fx(1,n+1,m+1)+(.5*(.5-dx)**2)*fx1,(n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(1,n,m-1)+
c (.5*(.5+dx)**2)*fx(1,n+1,m-1)+(.5*(.5-dx)**2)*fx(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2 + vz(t)**2)
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl, dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm, dth, edgelx, edgely, edgerx, edgery
      double precision sum1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,nl,mm) + amx*daxy(1,nn,mm) + dxp*daxy(1,np,m
     1m)) + dyl*(dxl*daxy(1,nl,ml) + amx*daxy(1,nn,ml) + dxp*daxy(1,np,m
     2l)) + dyp*(dxl*daxy(1,nl,mp) + amx*daxy(1,nn,mp) + dxp*daxy(1,np,m
     3p))
      oyx = amy*(dxl*daxy(2,nl,mm) + amx*daxy(2,nn,mm) + dxp*daxy(2,np,m
     1m)) + dyl*(dxl*daxy(2,nl,ml) + amx*daxy(2,nn,ml) + dxp*daxy(2,np,m
     2l)) + dyp*(dxl*daxy(2,nl,mp) + amx*daxy(2,nn,mp) + dxp*daxy(2,np,m
     3p))  
      oxy = amy*(dxl*daxy(3,nl,mm) + amx*daxy(3,nn,mm) + dxp*daxy(3,np,m
     1m)) + dyl*(dxl*daxy(3,nl,ml) + amx*daxy(3,nn,ml) + dxp*daxy(3,np,m
     2l)) + dyp*(dxl*daxy(3,nl,mp) + amx*daxy(3,nn,mp) + dxp*daxy(3,np,m
     3p))
      ozx = amy*(dxl*daxy(4,nl,mm) + amx*daxy(4,nn,mm) + dxp*daxy(4,np,m
     1m)) + dyl*(dxl*daxy(4,nl,ml) + amx*daxy(4,nn,ml) + dxp*daxy(4,np,m
     2l)) + dyp*(dxl*daxy(4,nl,mp) + amx*daxy(4,nn,mp) + dxp*daxy(4,np,m
     3p))  
      ozy = amy*(dxl*daxy(5,nl,mm) + amx*daxy(5,nn,mm) + dxp*daxy(5,np,m
     1m)) + dyl*(dxl*daxy(5,nl,ml) + amx*daxy(5,nn,ml) + dxp*daxy(5,np,m
     2l)) + dyp*(dxl*daxy(5,nl,mp) + amx*daxy(5,nn,mp) + dxp*daxy(5,np,m
     3p))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))
      dz = amy*(dxl*fxy(3,nl,mm) + amx*fxy(3,nn,mm) + dxp*fxy(3,np,mm))
     1+ dyl*(dxl*fxy(3,nl,ml) + amx*fxy(3,nn,ml) + dxp*fxy(3,np,ml)) + d
     2yp*(dxl*fxy(3,nl,mp) + amx*fxy(3,nn,mp) + dxp*fxy(3,np,mp)) 
c find vector potential at t
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
      az = amy*(dxl*axy(3,nl,mm) + amx*axy(3,nn,mm) + dxp*axy(3,np,mm))
     1+ dyl*(dxl*axy(3,nl,ml) + amx*axy(3,nn,ml) + dxp*axy(3,np,ml)) + d
     2yp*(dxl*axy(3,nl,mp) + amx*axy(3,nn,mp) + dxp*axy(3,np,mp))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,j) = dx
      part(7,j) = dy
      part(8,j) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c second-order spline interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 314 flops/particle, 1 divide, 104 loads, 8 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2 + vz(t)**2)
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm, dth, edgelx, edgely, edgerx, edgery
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2)) 
      ozx = amy*(dxl*daxy(4,mn) + amx*daxy(4,mn+1) + dxp*daxy(4,mn+2)) +
     1 dyl*(dxl*daxy(4,ml) + amx*daxy(4,ml+1) + dxp*daxy(4,ml+2)) + dyp*
     2(dxl*daxy(4,mp) + amx*daxy(4,mp+1) + dxp*daxy(4,mp+2)) 
      ozy = amy*(dxl*daxy(5,mn) + amx*daxy(5,mn+1) + dxp*daxy(5,mn+2)) +
     1 dyl*(dxl*daxy(5,ml) + amx*daxy(5,ml+1) + dxp*daxy(5,ml+2)) + dyp*
     2(dxl*daxy(5,mp) + amx*daxy(5,mp+1) + dxp*daxy(5,mp+2)) 
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
      dz = amy*(dxl*fxy(3,mn) + amx*fxy(3,mn+1) + dxp*fxy(3,mn+2)) + dyl
     1*(dxl*fxy(3,ml) + amx*fxy(3,ml+1) + dxp*fxy(3,ml+2)) + dyp*(dxl*fx
     2y(3,mp) + amx*fxy(3,mp+1) + dxp*fxy(3,mp+2))
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,j) = dx
      part(7,j) = dy
      part(8,j) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2)) 
      ozx = amy*(dxl*daxy(4,mn) + amx*daxy(4,mn+1) + dxp*daxy(4,mn+2)) +
     1 dyl*(dxl*daxy(4,ml) + amx*daxy(4,ml+1) + dxp*daxy(4,ml+2)) + dyp*
     2(dxl*daxy(4,mp) + amx*daxy(4,mp+1) + dxp*daxy(4,mp+2)) 
      ozy = amy*(dxl*daxy(5,mn) + amx*daxy(5,mn+1) + dxp*daxy(5,mn+2)) +
     1 dyl*(dxl*daxy(5,ml) + amx*daxy(5,ml+1) + dxp*daxy(5,ml+2)) + dyp*
     2(dxl*daxy(5,mp) + amx*daxy(5,mp+1) + dxp*daxy(5,mp+2)) 
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
      dz = amy*(dxl*fxy(3,mn) + amx*fxy(3,mn+1) + dxp*fxy(3,mn+2)) + dyl
     1*(dxl*fxy(3,ml) + amx*fxy(3,ml+1) + dxp*fxy(3,ml+2)) + dyp*(dxl*fx
     2y(3,mp) + amx*fxy(3,mp+1) + dxp*fxy(3,mp+2))
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,nop)
      py = part(7,nop)
      pz = part(8,nop)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,nop) = dx
      part(7,nop) = dy
      part(8,nop) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,nop) = px
      part(4,nop) = py
      part(5,nop) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(7,nop) = -part(7,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
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
      subroutine GHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 94 flops/particle, 32 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dth, dxp, dyp, amx, amy, dxl, dyl, dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
      az = amy*(dxl*axy(3,nl,mm) + amx*axy(3,nn,mm) + dxp*axy(3,np,mm))
     1+ dyl*(dxl*axy(3,nl,ml) + amx*axy(3,nn,ml) + dxp*axy(3,np,ml)) + d
     2yp*(dxl*axy(3,nl,mp) + amx*axy(3,nn,mp) + dxp*axy(3,np,mp))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 94 flops/particle, 32 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real dth, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,nop) - ax
      py = part(7,nop) - ay
      pz = part(8,nop) - az
c new position at time t + dt/2
      dx = part(1,nop) + (px - part(3,nop))*dth
      dy = part(2,nop) + (py - part(4,nop))*dth
      part(3,nop) = px
      part(4,nop) = py
      part(5,nop) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,nop) = -part(7,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      return
      end
      subroutine GHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 94 flops/particle, 32 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dxp, dyp, amx, amy, dxl, dyl, dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
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
c find vector potential at t
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
      az = amy*(dxl*axy(3,nl,mm) + amx*axy(3,nn,mm) + dxp*axy(3,np,mm))
     1+ dyl*(dxl*axy(3,nl,ml) + amx*axy(3,nn,ml) + dxp*axy(3,np,ml)) + d
     2yp*(dxl*axy(3,nl,mp) + amx*axy(3,nn,mp) + dxp*axy(3,np,mp))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
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
      subroutine GSHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ipb
     1c)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 94 flops/particle, 32 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
      az = amy*(dxl*axy(3,mn) + amx*axy(3,mn+1) + dxp*axy(3,mn+2)) + dyl
     1*(dxl*axy(3,ml) + amx*axy(3,ml+1) + dxp*axy(3,ml+2)) + dyp*(dxl*ax
     2y(3,mp) + amx*axy(3,mp+1) + dxp*axy(3,mp+2))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,nop) - ax
      py = part(7,nop) - ay
      pz = part(8,nop) - az
c new position at time t + dt
      dx = part(1,nop) + (px - 0.5*part(3,j))*dt
      dy = part(2,nop) + (py - 0.5*part(4,j))*dt
      part(3,nop) = px
      part(4,nop) = py
      part(5,nop) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,nop) = -part(7,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
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
      subroutine CR8PCM22(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates the canonical momentum of
c a particle from the vector potential at t - dt/2, and adds this
c information to the particle array, using quadratic interpolation
c canonical momentum p(t-dt/2) = v(t-dt/2) + qi*A(t-dt/2)/(c*mi)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c on entry
c part(5,n) = position x of particle n at t - dt/2
c part(6,n) = position y of particle n at t - dt/2
c on exit
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c axy(1,j+1,k+1) = x component of vector potential at grid (j,k)
c axy(2,j+1,k+1) = y component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c sx/sy = x/ycomponents of field momentum
c idimp = size of phase space = 6
c nop = number of particles
c nxv = first dimension of field array, must be >= nx+3
c nyv = second dimension of field array, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, axy, qbm, sx, sy
      dimension part(idimp,nop)
      dimension axy(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dxp, dyp, amx, amy, dxl, dyl, ax, ay
      double precision sum1, sum2
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 10 j = 1, nop
c find interpolation weights
      nn = part(5,j) + .5
      mm = part(6,j) + .5
      dxp = part(5,j) - float(nn)
      dyp = part(6,j) - float(mm)
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
c find vector potential at t - dt/2
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      sum1 = sum1 + ax
      sum2 = sum2 + ay
      part(5,j) = part(3,j) + ax
      part(6,j) = part(4,j) + ay
   10 continue
      sx = sum1
      sy = sum2
      return
      end
      subroutine GHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv,
     1nyv)
c for 2d code, this subroutine calculates current density, using
c second-order spline interpolation.
c scalar version using guard cells
c 257 flops/particle, 1 divide, 85 loads, 18 stores
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
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j+1,k+1) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension cu(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, dxp, dyp, amx, amy, dxl, dyl, dx, dy, dz, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,nl,mm) + amx*daxy(1,nn,mm) + dxp*daxy(1,np,m
     1m)) + dyl*(dxl*daxy(1,nl,ml) + amx*daxy(1,nn,ml) + dxp*daxy(1,np,m
     2l)) + dyp*(dxl*daxy(1,nl,mp) + amx*daxy(1,nn,mp) + dxp*daxy(1,np,m
     3p))
      oyx = amy*(dxl*daxy(2,nl,mm) + amx*daxy(2,nn,mm) + dxp*daxy(2,np,m
     1m)) + dyl*(dxl*daxy(2,nl,ml) + amx*daxy(2,nn,ml) + dxp*daxy(2,np,m
     2l)) + dyp*(dxl*daxy(2,nl,mp) + amx*daxy(2,nn,mp) + dxp*daxy(2,np,m
     3p))  
      oxy = amy*(dxl*daxy(3,nl,mm) + amx*daxy(3,nn,mm) + dxp*daxy(3,np,m
     1m)) + dyl*(dxl*daxy(3,nl,ml) + amx*daxy(3,nn,ml) + dxp*daxy(3,np,m
     2l)) + dyp*(dxl*daxy(3,nl,mp) + amx*daxy(3,nn,mp) + dxp*daxy(3,np,m
     3p))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))
c find vector potential at t
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      cu(1,nl,mm) = cu(1,nl,mm) + px*dx
      cu(2,nl,mm) = cu(2,nl,mm) + py*dx
      dx = dxl*dyl
      cu(1,nn,mm) = cu(1,nn,mm) + px*dy
      cu(2,nn,mm) = cu(2,nn,mm) + py*dy
      dy = amx*dyl
      cu(1,np,mm) = cu(1,np,mm) + px*dz
      cu(2,np,mm) = cu(2,np,mm) + py*dz
      dz = dxp*dyl
      cu(1,nl,ml) = cu(1,nl,ml) + px*dx
      cu(2,nl,ml) = cu(2,nl,ml) + py*dx
      dx = dxl*dyp
      cu(1,nn,ml) = cu(1,nn,ml) + px*dy
      cu(2,nn,ml) = cu(2,nn,ml) + py*dy
      dy = amx*dyp
      cu(1,np,ml) = cu(1,np,ml) + px*dz
      cu(2,np,ml) = cu(2,np,ml) + py*dz
      dz = dxp*dyp
      cu(1,nl,mp) = cu(1,nl,mp) + px*dx
      cu(2,nl,mp) = cu(2,nl,mp) + py*dx
      cu(1,nn,mp) = cu(1,nn,mp) + px*dy
      cu(2,nn,mp) = cu(2,nn,mp) + py*dy
      cu(1,np,mp) = cu(1,np,mp) + px*dz
      cu(2,np,mp) = cu(2,np,mp) + py*dz
   10 continue
      return
      end
      subroutine GSHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv
     1,nxyv)
c for 2d code, this subroutine calculates current density, using
c second-order spline interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 257 flops/particle, 1 divide, 85 loads, 18 stores
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
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j+1,k+1) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension cu(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real qtmh, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      real dx1, dy1, dx2, dy2, dx3, dy3
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2))  
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      dx1 = cu(1,mn) + px*dx
      dy1 = cu(2,mn) + py*dx
      dx2 = cu(1,mn+1) + px*dy
      dy2 = cu(2,mn+1) + py*dy
      dy3 = cu(1,mn+2) + px*dz
      dx3 = cu(2,mn+2) + py*dz
      cu(1,mn) = dx1
      cu(2,mn) = dy1
      cu(1,mn+1) = dx2
      cu(2,mn+1) = dy2
      cu(1,mn+2) = dy3
      cu(2,mn+2) = dx3
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      dx1 = cu(1,ml) + px*dx
      dy1 = cu(2,ml) + py*dx
      dx2 = cu(1,ml+1) + px*dy
      dy2 = cu(2,ml+1) + py*dy
      dy3 = cu(1,ml+2) + px*dz
      dx3 = cu(2,ml+2) + py*dz
      cu(1,ml) = dx1
      cu(2,ml) = dy1
      cu(1,ml+1) = dx2
      cu(2,ml+1) = dy2
      cu(1,ml+2) = dy3
      cu(2,ml+2) = dx3
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      dx1 = cu(1,mp) + px*dx
      dy1 = cu(2,mp) + py*dx
      dx2 = cu(1,mp+1) + px*dy
      dy2 = cu(2,mp+1) + py*dy
      dy3 = cu(1,mp+2) + px*dz
      dx3 = cu(2,mp+2) + py*dz
      cu(1,mp) = dx1
      cu(2,mp) = dy1
      cu(1,mp+1) = dx2
      cu(2,mp+1) = dy2
      cu(1,mp+2) = dy3
      cu(2,mp+2) = dx3
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2))  
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,nop)
      py = part(6,nop)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      cu(1,mn) = cu(1,mn) + px*dx
      cu(2,mn) = cu(2,mn) + py*dx
      cu(1,mn+1) = cu(1,mn+1) + px*dy
      cu(2,mn+1) = cu(2,mn+1) + py*dy
      cu(1,mn+2) = cu(1,mn+2) + px*dz
      cu(2,mn+2) = cu(2,mn+2) + py*dz
      dx = dxl*dyl
      dy = amx*dyl
      dz = dxp*dyl
      cu(1,ml) = cu(1,ml) + px*dx
      cu(2,ml) = cu(2,ml) + py*dx
      cu(1,ml+1) = cu(1,ml+1) + px*dy
      cu(2,ml+1) = cu(2,ml+1) + py*dy
      cu(1,ml+2) = cu(1,ml+2) + px*dz
      cu(2,ml+2) = cu(2,ml+2) + py*dz
      dx = dxl*dyp
      dy = amx*dyp
      dz = dxp*dyp
      cu(1,mp) = cu(1,mp) + px*dx
      cu(2,mp) = cu(2,mp) + py*dx
      cu(1,mp+1) = cu(1,mp+1) + px*dy
      cu(2,mp+1) = cu(2,mp+1) + py*dy
      cu(1,mp+2) = cu(1,mp+2) + px*dz
      cu(2,mp+2) = cu(2,mp+2) + py*dz
      return
      end
      subroutine GHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c second-order spline interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 217 flops/particle, 1 divide, 67 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2)
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dth, qtmh, dxp, dyp, amx, amy, dxl, dyl, dx, dy, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      real edgelx, edgely, edgerx, edgery
      double precision sum1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,nl,mm) + amx*daxy(1,nn,mm) + dxp*daxy(1,np,m
     1m)) + dyl*(dxl*daxy(1,nl,ml) + amx*daxy(1,nn,ml) + dxp*daxy(1,np,m
     2l)) + dyp*(dxl*daxy(1,nl,mp) + amx*daxy(1,nn,mp) + dxp*daxy(1,np,m
     3p))
      oyx = amy*(dxl*daxy(2,nl,mm) + amx*daxy(2,nn,mm) + dxp*daxy(2,np,m
     1m)) + dyl*(dxl*daxy(2,nl,ml) + amx*daxy(2,nn,ml) + dxp*daxy(2,np,m
     2l)) + dyp*(dxl*daxy(2,nl,mp) + amx*daxy(2,nn,mp) + dxp*daxy(2,np,m
     3p))  
      oxy = amy*(dxl*daxy(3,nl,mm) + amx*daxy(3,nn,mm) + dxp*daxy(3,np,m
     1m)) + dyl*(dxl*daxy(3,nl,ml) + amx*daxy(3,nn,ml) + dxp*daxy(3,np,m
     2l)) + dyp*(dxl*daxy(3,nl,mp) + amx*daxy(3,nn,mp) + dxp*daxy(3,np,m
     3p))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = amy*(dxl*fxy(1,nl,mm) + amx*fxy(1,nn,mm) + dxp*fxy(1,np,mm))
     1+ dyl*(dxl*fxy(1,nl,ml) + amx*fxy(1,nn,ml) + dxp*fxy(1,np,ml)) + d
     2yp*(dxl*fxy(1,nl,mp) + amx*fxy(1,nn,mp) + dxp*fxy(1,np,mp))              
      dy = amy*(dxl*fxy(2,nl,mm) + amx*fxy(2,nn,mm) + dxp*fxy(2,np,mm))
     1+ dyl*(dxl*fxy(2,nl,ml) + amx*fxy(2,nn,ml) + dxp*fxy(2,np,ml)) + d
     2yp*(dxl*fxy(2,nl,mp) + amx*fxy(2,nn,mp) + dxp*fxy(2,np,mp))
c find vector potential at t
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
      part(5,j) = dx
      part(6,j) = dy
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      part(3,j) = px
      part(4,j) = py
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,nxyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c second-order spline interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 217 flops/particle, 1 divide, 67 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*fxy(1,n,m)+(.5*(.5+dx)**2)*fxy(1,n+1,m)+
c (.5*(.5-dx)**2)*fxy(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*fxy(1,n,m+1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m+1)+(.5*(.5-dx)**2)*fxy(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fxy(1,n,m-1)+
c (.5*(.5+dx)**2)*fxy(1,n+1,m-1)+(.5*(.5-dx)**2)*fxy(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j+1,k+1) = ith component of force/charge at grid (j,k)
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c daxy(i,j+1,k+1) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2)
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real dth, qtmh, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      real edgelx, edgely, edgerx, edgery
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2))  
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp) + amx*fxy(2,mp+1) + dxp*fxy(2,mp+2)) 
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
      part(5,j) = dx
      part(6,j) = dy
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      part(3,j) = px
      part(4,j) = py
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
c find vector potential gradients at t
      oxx = amy*(dxl*daxy(1,mn) + amx*daxy(1,mn+1) + dxp*daxy(1,mn+2)) +
     1 dyl*(dxl*daxy(1,ml) + amx*daxy(1,ml+1) + dxp*daxy(1,ml+2)) + dyp*
     2(dxl*daxy(1,mp) + amx*daxy(1,mp+1) + dxp*daxy(1,mp+2))
      oyx = amy*(dxl*daxy(2,mn) + amx*daxy(2,mn+1) + dxp*daxy(2,mn+2)) +
     1 dyl*(dxl*daxy(2,ml) + amx*daxy(2,ml+1) + dxp*daxy(2,ml+2)) + dyp*
     2(dxl*daxy(2,mp) + amx*daxy(2,mp+1) + dxp*daxy(2,mp+2)) 
      oxy = amy*(dxl*daxy(3,mn) + amx*daxy(3,mn+1) + dxp*daxy(3,mn+2)) +
     1 dyl*(dxl*daxy(3,ml) + amx*daxy(3,ml+1) + dxp*daxy(3,ml+2)) + dyp*
     2(dxl*daxy(3,mp) + amx*daxy(3,mp+1) + dxp*daxy(3,mp+2))  
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = amy*(dxl*fxy(1,mn) + amx*fxy(1,mn+1) + dxp*fxy(1,mn+2)) + dyl
     1*(dxl*fxy(1,ml) + amx*fxy(1,ml+1) + dxp*fxy(1,ml+2)) + dyp*(dxl*fx
     2y(1,mp) + amx*fxy(1,mp+1) + dxp*fxy(1,mp+2))              
      dy = amy*(dxl*fxy(2,mn) + amx*fxy(2,mn+1) + dxp*fxy(2,mn+2)) + dyl
     1*(dxl*fxy(2,ml) + amx*fxy(2,ml+1) + dxp*fxy(2,ml+2)) + dyp*(dxl*fx
     2y(2,mp)+ amx*fxy(2,mp+1) + dxp*fxy(2,mp+2))  
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,nop)
      py = part(6,nop)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
      part(5,nop) = dx
      part(6,nop) = dy
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      part(3,nop) = px
      part(4,nop) = py
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py)
c new position at time t + dt/2
      dx = part(1,nop) + px*dth
      dy = part(2,nop) + py*dth
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
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(6,nop) = -part(6,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
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
      subroutine GHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 74 flops/particle, 22 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dth, dxp, dyp, amx, amy, dxl, dyl, dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 74 flops/particle, 22 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real dth, dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
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
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,nop) - ax
      py = part(6,nop) - ay
c new position at time t + dt/2
      dx = part(1,nop) + (px - part(3,nop))*dth
      dy = part(2,nop) + (py - part(4,nop))*dth
      part(3,nop) = px
      part(4,nop) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,nop) = -part(6,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      return
      end
      subroutine GHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 74 flops/particle, 22 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nyv = second dimension of field arrays, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real dxp, dyp, amx, amy, dxl, dyl, dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
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
c find vector potential at t+dt/2
      ax = amy*(dxl*axy(1,nl,mm) + amx*axy(1,nn,mm) + dxp*axy(1,np,mm))
     1+ dyl*(dxl*axy(1,nl,ml) + amx*axy(1,nn,ml) + dxp*axy(1,np,ml)) + d
     2yp*(dxl*axy(1,nl,mp) + amx*axy(1,nn,mp) + dxp*axy(1,np,mp))              
      ay = amy*(dxl*axy(2,nl,mm) + amx*axy(2,nn,mm) + dxp*axy(2,np,mm))
     1+ dyl*(dxl*axy(2,nl,ml) + amx*axy(2,nn,ml) + dxp*axy(2,np,ml)) + d
     2yp*(dxl*axy(2,nl,mp) + amx*axy(2,nn,mp) + dxp*axy(2,np,mp))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
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
      subroutine GSHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ipb
     1c)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 74 flops/particle, 22 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = 
c (.75-dy**2)*((.75-dx**2)*ax(1,n,m)+(.5*(.5+dx)**2)*ax(1,n+1,m)+
c (.5*(.5-dx)**2)*ax(1,n-1,m))
c + (.5*(.5+dy)**2)*((.75-dx**2)*ax(1,n,m+1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m+1)+(.5*(.5-dx)**2)*ax(1,n-1,m+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*ax(1,n,m-1)+
c (.5*(.5+dx)**2)*ax(1,n+1,m-1)+(.5*(.5-dx)**2)*ax(1,n-1,m-1))
c where n,m = nearest grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j+1,k+1) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+3
c nxyv = actual dimension of field array, must be >= nxv*(ny+3)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, ml, mn, mp
      real dxn, dyn, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1) + .5
      mmn = part(2,1) + .5
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
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
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
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
c find vector potential at t
      ax = amy*(dxl*axy(1,mn) + amx*axy(1,mn+1) + dxp*axy(1,mn+2)) + dyl
     1*(dxl*axy(1,ml) + amx*axy(1,ml+1) + dxp*axy(1,ml+2)) + dyp*(dxl*ax
     2y(1,mp) + amx*axy(1,mp+1) + dxp*axy(1,mp+2))              
      ay = amy*(dxl*axy(2,mn) + amx*axy(2,mn+1) + dxp*axy(2,mn+2)) + dyl
     1*(dxl*axy(2,ml) + amx*axy(2,ml+1) + dxp*axy(2,ml+2)) + dyp*(dxl*ax
     2y(2,mp) + amx*axy(2,mp+1) + dxp*axy(2,mp+2)) 
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,nop) - ax
      py = part(6,nop) - ay
c new position at time t + dt
      dx = part(1,nop) + (px - 0.5*part(3,nop))*dt
      dy = part(2,nop) + (py - 0.5*part(4,nop))*dt
      part(3,nop) = px
      part(4,nop) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,nop) = -part(6,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
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
      subroutine CR8PCM23L(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv)
c for 2-1/2d code, this subroutine calculates the canonical momentum of
c a particle from the vector potential at t - dt/2, and adds this
c information to the particle array, using linear interpolation
c canonical momentum p(t-dt/2) = v(t-dt/2) + qi*A(t-dt/2)/(c*mi)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c on entry
c part(6,n) = position x of particle n at t - dt/2
c part(7,n) = position y of particle n at t - dt/2
c on exit
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c axy(1,j,k) = x component of vector potential at grid (j,k)
c axy(2,j,k) = y component of vector potential at grid (j,k)
c axy(3,j,k) = z component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c sx/sy/sz = x/y/z components of field momentum
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field array, must be >= nx+1
c nyv = second dimension of field array, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, axy, qbm, sx, sy, sz
      dimension part(idimp,nop)
      dimension axy(3,nxv,nyv)
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy, ax, ay, az
      double precision sum1, sum2, sum3
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 10 j = 1, nop
c find interpolation weights
      nn = part(6,j)
      mm = part(7,j)
      dxp = part(6,j) - float(nn)
      dyp = part(7,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find vector potential at t - dt/2
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
      az = dyp*(dxp*axy(3,np,mp) + amx*axy(3,nn,mp)) + amy*(dxp*axy(3,np
     1,mm) + amx*axy(3,nn,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      sum1 = sum1 + ax
      sum2 = sum2 + ay
      sum3 = sum3 + az
      part(6,j) = part(3,j) + ax
      part(7,j) = part(4,j) + ay
      part(8,j) = part(5,j) + az
   10 continue
      sx = sum1
      sy = sum2
      sz = sum3
      return
      end
      subroutine GHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv,
     1nyv)
c for 2-1/2d code, this subroutine calculates current density, using
c first-order linear interpolation.
c scalar version using guard cells
c 208 flops/particle, 1 divide, 61 loads, 12 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j,k) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension cu(3,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, dxp, dyp, amx, amy, dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,np,mp) + amx*daxy(1,nn,mp)) + amy*(dxp*daxy(
     11,np,mm) + amx*daxy(1,nn,mm))
      oyx = dyp*(dxp*daxy(2,np,mp) + amx*daxy(2,nn,mp)) + amy*(dxp*daxy(
     12,np,mm) + amx*daxy(2,nn,mm))
      oxy = dyp*(dxp*daxy(3,np,mp) + amx*daxy(3,nn,mp)) + amy*(dxp*daxy(
     13,np,mm) + amx*daxy(3,nn,mm))
      ozx = dyp*(dxp*daxy(4,np,mp) + amx*daxy(4,nn,mp)) + amy*(dxp*daxy(
     14,np,mm) + amx*daxy(4,nn,mm))
      ozy = dyp*(dxp*daxy(5,np,mp) + amx*daxy(5,nn,mp)) + amy*(dxp*daxy(
     15,np,mm) + amx*daxy(5,nn,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp)) + amy*(dxp*fxy(3,np
     1,mm) + amx*fxy(3,nn,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
      az = dyp*(dxp*axy(3,np,mp) + amx*axy(3,nn,mp)) + amy*(dxp*axy(3,np
     1,mm) + amx*axy(3,nn,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      dx = dxp*dyp
      dy = amx*dyp
      cu(1,np,mp) = cu(1,np,mp) + px*dx
      cu(2,np,mp) = cu(2,np,mp) + py*dx
      cu(3,np,mp) = cu(3,np,mp) + pz*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + px*dy
      cu(2,nn,mp) = cu(2,nn,mp) + py*dy
      cu(3,nn,mp) = cu(3,nn,mp) + pz*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + px*dx
      cu(2,np,mm) = cu(2,np,mm) + py*dx
      cu(3,np,mm) = cu(3,np,mm) + pz*dx
      cu(1,nn,mm) = cu(1,nn,mm) + px*dy
      cu(2,nn,mm) = cu(2,nn,mm) + py*dy
      cu(3,nn,mm) = cu(3,nn,mm) + pz*dy
   10 continue
      return
      end
      subroutine GSHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv
     1,nxyv)
c for 2-1/2d code, this subroutine calculates current density, using
c first-order linear interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 208 flops/particle, 1 divide, 61 loads, 12 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 3
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j,k) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension cu(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real qtmh, dxn, dyn, dxp, dyp, amx, amy, dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm, dx1, dy1, dx2, dy2
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxp*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyp*(dxp*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxp*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyp*(dxp*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxp*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
      ozx = dyp*(dxp*daxy(4,mp+1) + amx*daxy(4,mp)) + amy*(dxp*daxy(4,mm
     1+1)+ amx*daxy(4,mm))
      ozy = dyp*(dxp*daxy(5,mp+1) + amx*daxy(5,mp)) + amy*(dxp*daxy(5,mm
     1+1)+ amx*daxy(5,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyp*(dxp*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxp*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      dx = dxp*dyp
      dz = amx*dyp
      dx1 = cu(1,mp+1) + px*dx
      dy1 = cu(2,mp+1) + py*dx
      dyp = cu(3,mp+1) + pz*dx
      dx2 = cu(1,mp) + px*dz
      dy2 = cu(2,mp) + py*dz
      dy = cu(3,mp) + pz*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(3,mp+1) = dyp
      cu(1,mp) = dx2
      cu(2,mp) = dy2
      cu(3,mp) = dy
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + px*dx
      amx = cu(2,mm+1) + py*dx
      dyp = cu(3,mm+1) + pz*dx
      dx1 = cu(1,mm) + px*dz
      dy1 = cu(2,mm) + py*dz
      dy = cu(3,mm) + pz*dz
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
c find vector potential gradients at t
      oxx = dyn*(dxn*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxn*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyn*(dxn*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxn*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyn*(dxn*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxn*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
      ozx = dyn*(dxn*daxy(4,mp+1) + amx*daxy(4,mp)) + amy*(dxn*daxy(4,mm
     1+1)+ amx*daxy(4,mm))
      ozy = dyn*(dxn*daxy(5,mp+1) + amx*daxy(5,mp)) + amy*(dxn*daxy(5,mm
     1+1)+ amx*daxy(5,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find vector potential at t
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyn*(dxn*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxn*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,nop)
      py = part(7,nop)
      pz = part(8,nop)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c deposit current density
      amx = qm*amx
      dxp = qm*dxn
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      dx = dxp*dyn
      dy = amx*dyn
      cu(1,mp+1) = cu(1,mp+1) + px*dx
      cu(2,mp+1) = cu(2,mp+1) + py*dx
      cu(3,mp+1) = cu(3,mp+1) + pz*dx
      cu(1,mp) = cu(1,mp) + px*dy
      cu(2,mp) = cu(2,mp) + py*dy
      cu(3,mp) = cu(3,mp) + pz*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + px*dx
      cu(2,mm+1) = cu(2,mm+1) + py*dx
      cu(3,mm+1) = cu(3,mm+1) + pz*dx
      cu(1,mm) = cu(1,mm) + px*dy
      cu(2,mm) = cu(2,mm) + py*dy
      cu(3,mm) = cu(3,mm) + pz*dy
      return
      end
      subroutine GHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c first-order linear interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 188 flops/particle, 1 divide, 49 loads, 8 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2 + vz(t)**2)
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, dxp, dyp, amx, amy, dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm, dth, edgelx, edgely, edgerx, edgery
      double precision sum1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,np,mp) + amx*daxy(1,nn,mp)) + amy*(dxp*daxy(
     11,np,mm) + amx*daxy(1,nn,mm))
      oyx = dyp*(dxp*daxy(2,np,mp) + amx*daxy(2,nn,mp)) + amy*(dxp*daxy(
     12,np,mm) + amx*daxy(2,nn,mm))
      oxy = dyp*(dxp*daxy(3,np,mp) + amx*daxy(3,nn,mp)) + amy*(dxp*daxy(
     13,np,mm) + amx*daxy(3,nn,mm))
      ozx = dyp*(dxp*daxy(4,np,mp) + amx*daxy(4,nn,mp)) + amy*(dxp*daxy(
     14,np,mm) + amx*daxy(4,nn,mm))
      ozy = dyp*(dxp*daxy(5,np,mp) + amx*daxy(5,nn,mp)) + amy*(dxp*daxy(
     15,np,mm) + amx*daxy(5,nn,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp)) + amy*(dxp*fxy(3,np
     1,mm) + amx*fxy(3,nn,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
      az = dyp*(dxp*axy(3,np,mp) + amx*axy(3,nn,mp)) + amy*(dxp*axy(3,np
     1,mm) + amx*axy(3,nn,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,j) = dx
      part(7,j) = dy
      part(8,j) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c first-order linear interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 188 flops/particle, 1 divide, 49 loads, 8 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2 + vz(t)**2)
c idimp = size of phase space = 8
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real qtmh, dxn, dyn, dxp, dyp, amx, amy, dx, dy, dz, ax, ay, az
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      real anorm, dth, edgelx, edgely, edgerx, edgery
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxp*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyp*(dxp*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxp*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyp*(dxp*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxp*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
      ozx = dyp*(dxp*daxy(4,mp+1) + amx*daxy(4,mp)) + amy*(dxp*daxy(4,mm
     1+1)+ amx*daxy(4,mm))
      ozy = dyp*(dxp*daxy(5,mp+1) + amx*daxy(5,mp)) + amy*(dxp*daxy(5,mm
     1+1)+ amx*daxy(5,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyp*(dxp*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxp*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyp*(dxp*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxp*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,j)
      py = part(7,j)
      pz = part(8,j)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,j) = dx
      part(7,j) = dy
      part(8,j) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
c find vector potential gradients at t
      oxx = dyn*(dxn*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxn*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyn*(dxn*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxn*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyn*(dxn*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxn*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
      ozx = dyn*(dxn*daxy(4,mp+1) + amx*daxy(4,mp)) + amy*(dxn*daxy(4,mm
     1+1)+ amx*daxy(4,mm))
      ozy = dyn*(dxn*daxy(5,mp+1) + amx*daxy(5,mp)) + amy*(dxn*daxy(5,mm
     1+1)+ amx*daxy(5,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field at t
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
      dz = dyn*(dxn*fxy(3,mp+1) + amx*fxy(3,mp)) + amy*(dxn*fxy(3,mm+1)
     1+ amx*fxy(3,mm))
c find vector potential at t
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyn*(dxn*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxn*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,nop)
      py = part(7,nop)
      pz = part(8,nop)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,nop) = dx
      part(7,nop) = dy
      part(8,nop) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,nop) = px
      part(4,nop) = py
      part(5,nop) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,nop) + px*dth
      dy = part(2,nop) + py*dth
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
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(7,nop) = -part(7,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
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
      subroutine GHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipb
     1c)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 47 flops/particle, 19 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxv,nyv)
      integer j, nn, mm, np, mp
      real dth, dxp, dyp, amx, amy, dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
      az = dyp*(dxp*axy(3,np,mp) + amx*axy(3,nn,mp)) + amy*(dxp*axy(3,np
     1,mm) + amx*axy(3,nn,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 47 flops/particle, 19 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real dth, dxn, dyn, dxp, dyp, amx, amy, dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyp*(dxp*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxp*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = 0.0
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
c find vector potential at t+dt/2
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyn*(dxn*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxn*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,nop) - ax
      py = part(7,nop) - ay
      pz = part(8,nop) - az
c new position at time t + dt/2
      dx = part(1,nop) + (px - part(3,nop))*dth
      dy = part(2,nop) + (py - part(4,nop))*dth
      part(3,nop) = px
      part(4,nop) = py
      part(5,nop) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,nop) = -part(7,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
            part(3,nop) = 0.0
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      return
      end
      subroutine GHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 47 flops/particle, 19 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxv,nyv)
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy, dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
      az = dyp*(dxp*axy(3,np,mp) + amx*axy(3,nn,mp)) + amy*(dxp*axy(3,np
     1,mm) + amx*axy(3,nn,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
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
      subroutine GSHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 47 flops/particle, 19 loads, 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real dxn, dyn, dxp, dyp, amx, amy, dx, dy, ax, ay, az
      real px, py, pz, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyp*(dxp*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxp*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,j) - ax
      py = part(7,j) - ay
      pz = part(8,j) - az
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,j) = -part(7,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,j) = -part(6,j)
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
c find vector potential at t+dt/2
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
      az = dyn*(dxn*axy(3,mp+1) + amx*axy(3,mp)) + amy*(dxn*axy(3,mm+1)
     1+ amx*axy(3,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,nop) - ax
      py = part(7,nop) - ay
      pz = part(8,nop) - az
c new position at time t + dt
      dx = part(1,nop) + (px - 0.5*part(3,nop))*dt
      dy = part(2,nop) + (py - 0.5*part(4,nop))*dt
      part(3,nop) = px
      part(4,nop) = py
      part(5,nop) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,nop) = -part(7,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,nop) = -part(6,nop)
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
      subroutine CR8PCM22L(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates the canonical momentum of
c a particle from the vector potential at t - dt/2, and adds this
c information to the particle array, using linear interpolation
c canonical momentum p(t-dt/2) = v(t-dt/2) + qi*A(t-dt/2)/(c*mi)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c on entry
c part(5,n) = position x of particle n at t - dt/2
c part(6,n) = position y of particle n at t - dt/2
c on exit
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c axy(1,j,k) = x component of vector potential at grid (j,k)
c axy(2,j,k) = y component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c sx/sy = x/y components of field momentum
c idimp = size of phase space = 6
c nop = number of particles
c nxv = first dimension of field array, must be >= nx+1
c nyv = second dimension of field array, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, axy, qbm, sx, sy
      dimension part(idimp,nop)
      dimension axy(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy, ax, ay
      double precision sum1, sum2
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 10 j = 1, nop
c find interpolation weights
      nn = part(5,j)
      mm = part(6,j)
      dxp = part(5,j) - float(nn)
      dyp = part(6,j) - float(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1. - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c find vector potential at t - dt/2
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      sum1 = sum1 + ax
      sum2 = sum2 + ay
      part(5,j) = part(3,j) + ax
      part(6,j) = part(4,j) + ay
   10 continue
      sx = sum1
      sy = sum2
      return
      end
      subroutine GHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv
     1,nyv)
c for 2d code, this subroutine calculates current density, using
c first-order linear interpolation.
c scalar version using guard cells
c 136 flops/particle, 1 divide, 38 loads, 8 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j,k) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension cu(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, dxp, dyp, amx, amy, dx, dy, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,np,mp) + amx*daxy(1,nn,mp)) + amy*(dxp*daxy(
     11,np,mm) + amx*daxy(1,nn,mm))
      oyx = dyp*(dxp*daxy(2,np,mp) + amx*daxy(2,nn,mp)) + amy*(dxp*daxy(
     12,np,mm) + amx*daxy(2,nn,mm))
      oxy = dyp*(dxp*daxy(3,np,mp) + amx*daxy(3,nn,mp)) + amy*(dxp*daxy(
     13,np,mm) + amx*daxy(3,nn,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      dx = dxp*dyp
      dy = amx*dyp
      cu(1,np,mp) = cu(1,np,mp) + px*dx
      cu(2,np,mp) = cu(2,np,mp) + py*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + px*dy
      cu(2,nn,mp) = cu(2,nn,mp) + py*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + px*dx
      cu(2,np,mm) = cu(2,np,mm) + py*dx
      cu(1,nn,mm) = cu(1,nn,mm) + px*dy
      cu(2,nn,mm) = cu(2,nn,mm) + py*dy
   10 continue
      return
      end
      subroutine GSHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nx
     1v,nxyv)
c for 2d code, this subroutine calculates current density, using
c first-order linear interpolation.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 136 flops/particle, 1 divide, 38 loads, 8 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*vj, where j = x,y,z, for i = 1, 2
c where vj = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c cu(i,j,k) = ith component of current density at grid (j,k)
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
      implicit none
      integer idimp, nop, nxv, nxyv
      real part, fxy, axy, daxy, cu, qm, qbm, dt
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension cu(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real qtmh, dxn, dyn, dxp, dyp, amx, amy, dx, dy, dz, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      real dx1, dy1, dx2, dy2
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      qtmh = .5*qbm*dt
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxp*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyp*(dxp*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxp*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyp*(dxp*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxp*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      dx = dxp*dyp
      dz = amx*dyp
      dx1 = cu(1,mp+1) + px*dx
      dy1 = cu(2,mp+1) + py*dx
      dx2 = cu(1,mp) + px*dz
      dy2 = cu(2,mp) + py*dz
      cu(1,mp+1) = dx1
      cu(2,mp+1) = dy1
      cu(1,mp) = dx2
      cu(2,mp) = dy2
      dx = dxp*amy
      dz = amx*amy
      dxp = cu(1,mm+1) + px*dx
      amx = cu(2,mm+1) + py*dx
      dx1 = cu(1,mm) + px*dz
      dy1 = cu(2,mm) + py*dz
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
c find vector potential gradients at t
      oxx = dyn*(dxn*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxn*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyn*(dxn*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxn*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyn*(dxn*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxn*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find vector potential at t
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,nop)
      py = part(6,nop)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
c deposit current density
      amx = qm*amx
      dxp = qm*dxn
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      dx = dxp*dyn
      dy = amx*dyn
      cu(1,mp+1) = cu(1,mp+1) + px*dx
      cu(2,mp+1) = cu(2,mp+1) + py*dx
      cu(1,mp) = cu(1,mp) + px*dy
      cu(2,mp) = cu(2,mp) + py*dy
      dx = dxp*amy
      dy = amx*amy
      cu(1,mm+1) = cu(1,mm+1) + px*dx
      cu(2,mm+1) = cu(2,mm+1) + py*dx
      cu(1,mm) = cu(1,mm) + px*dy
      cu(2,mm) = cu(2,mm) + py*dy
      return
      end
      subroutine GHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c first-order linear interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 124 flops/particle, 1 divide, 32 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2)
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      integer j, nn, mm, np, mp
      real dth, qtmh, dxp, dyp, amx, amy, dx, dy, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      real edgelx, edgely, edgerx, edgery
      double precision sum1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,np,mp) + amx*daxy(1,nn,mp)) + amy*(dxp*daxy(
     11,np,mm) + amx*daxy(1,nn,mm))
      oyx = dyp*(dxp*daxy(2,np,mp) + amx*daxy(2,nn,mp)) + amy*(dxp*daxy(
     12,np,mm) + amx*daxy(2,nn,mm))
      oxy = dyp*(dxp*daxy(3,np,mp) + amx*daxy(3,nn,mp)) + amy*(dxp*daxy(
     13,np,mm) + amx*daxy(3,nn,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp)) + amy*(dxp*fxy(1,np
     1,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp)) + amy*(dxp*fxy(2,np
     1,mm) + amx*fxy(2,nn,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
      part(5,j) = dx
      part(6,j) = dy
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      part(3,j) = px
      part(4,j) = py
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,nxyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time and
c first-order linear interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 124 flops/particle, 1 divide, 32 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(3)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(4)*(py(t-dt/2) + .5*(q/m)*Fy*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(4) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c,
c Fx = Ex - 2*(Ax*omxx + Ay*omyx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx)/c
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)),
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fxy(1,x,y) = (1-dy)*((1-dx)*fxy(1,n,m)+dx*fxy(1,n+1,m)) +
c dy*((1-dx)*fxy(1,n,m+1)+ dx*fxy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = canonical momentum px of particle n at t - dt/2
c part(6,n) = canonical momentum py of particle n at t - dt/2
c fxy(i,j,k) = ith component of force/charge at grid (j,k)
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c daxy(i,j,k) = ith component of derivative of vector potential
c at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2)
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, fxy, axy, daxy, qbm, dt, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real dth, qtmh, dxn, dyn, dxp, dyp, amx, amy, dx, dy, ax, ay
      real oxx, oyx, oxy, omxxt, omyxt, omxyt
      real px, py, acx, acy, rot0, rot1, rot2, rot3, rot4, anorm
      real edgelx, edgely, edgerx, edgery
      double precision sum1
      sum1 = 0.0d0
      if (nop.lt.1) go to 20
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
      qtmh = qbm*dth
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
c find vector potential gradients at t
      oxx = dyp*(dxp*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxp*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyp*(dxp*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxp*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyp*(dxp*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxp*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = dyp*(dxp*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxp*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyp*(dxp*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxp*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find vector potential at t
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,j)
      py = part(6,j)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
      part(5,j) = dx
      part(6,j) = dy
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      part(3,j) = px
      part(4,j) = py
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py)
c new position at time t + dt/2
      dx = part(1,j) + px*dth
      dy = part(2,j) + py*dth
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
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
c find vector potential gradients at t
      oxx = dyn*(dxn*daxy(1,mp+1) + amx*daxy(1,mp)) + amy*(dxn*daxy(1,mm
     1+1)+ amx*daxy(1,mm))
      oyx = dyn*(dxn*daxy(2,mp+1) + amx*daxy(2,mp)) + amy*(dxn*daxy(2,mm
     1+1)+ amx*daxy(2,mm))
      oxy = dyn*(dxn*daxy(3,mp+1) + amx*daxy(3,mp)) + amy*(dxn*daxy(3,mm
     1+1)+ amx*daxy(3,mm))
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
c find electric field at t
      dx = dyn*(dxn*fxy(1,mp+1) + amx*fxy(1,mp)) + amy*(dxn*fxy(1,mm+1)
     1+ amx*fxy(1,mm))
      dy = dyn*(dxn*fxy(2,mp+1) + amx*fxy(2,mp)) + amy*(dxn*fxy(2,mm+1)
     1+ amx*fxy(2,mm))
c find vector potential at t
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      dx = qtmh*dx - (ax*omxxt + ay*omyxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt)
c half acceleration
      px = part(5,nop)
      py = part(6,nop)
      acx = px + dx
      acy = py + dy
c calculate "rotation" matrix
      rot0 = omyxt*omxyt
      rot1 = 1.0 + omxxt
      rot4 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot3 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot0)
      rot1 = rot1*rot1 + rot0
      rot4 = rot4*rot4 + rot0
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy)*anorm + dx
      dy = (rot3*acx + rot4*acy)*anorm + dy
      part(5,nop) = dx
      part(6,nop) = dy
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      part(3,nop) = px
      part(4,nop) = py
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py)
c new position at time t + dt/2
      dx = part(1,nop) + px*dth
      dy = part(2,nop) + py*dth
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
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop)
            part(6,nop) = -part(6,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop)
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
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
      subroutine GHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipb
     1c)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 36 flops/particle, 14 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real dth, dxp, dyp, amx, amy, dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
      subroutine GSHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 36 flops/particle, 14 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real dth, dxn, dyn, dxp, dyp, amx, amy, dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt/2
      dx = part(1,j) + (px - part(3,j))*dth
      dy = part(2,j) + (py - part(4,j))*dth
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = 0.0
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
c find vector potential at t+dt/2
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,nop) - ax
      py = part(6,nop) - ay
c new position at time t + dt/2
      dx = part(1,nop) + (px - part(3,nop))*dth
      dy = part(2,nop) + (py - part(4,nop))*dth
      part(3,nop) = px
      part(4,nop) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,nop) = -part(6,nop)
            part(4,nop) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
            part(3,nop) = 0.0
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,nop) = dx
      part(2,nop) = dy
      return
      end
      subroutine GHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells
c 36 flops/particle, 14 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nyv = second dimension of field arrays, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy, dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,np,mp) + amx*axy(1,nn,mp)) + amy*(dxp*axy(1,np
     1,mm) + amx*axy(1,nn,mm))
      ay = dyp*(dxp*axy(2,np,mp) + amx*axy(2,nn,mp)) + amy*(dxp*axy(2,np
     1,mm) + amx*axy(2,nn,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
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
      subroutine GSHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, periodic boundaries
c other boundary conditions are only approximate.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 36 flops/particle, 14 loads, 4 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) is approximated by interpolation from the nearest
c grid points:
c axy(1,x,y) = (1-dy)*((1-dx)*axy(1,n,m)+dx*axy(1,n+1,m)) +
c dy*((1-dx)*axy(1,n,m+1)+ dx*axy(1,n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for the other components
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = canonical momentum px of particle n at t + dt/2
c part(6,n) = canonical momentum py of particle n at t + dt/2
c axy(i,j,k) = ith component of vector potential at grid (j,k)
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 6
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+1
c nxyv = actual dimension of field array, must be >= nxv*(ny+1)
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc
      real part, axy, qbm, dt
      dimension part(idimp,nop), axy(2,nxyv)
      integer nop1, j, nnn, mmn, nn, mm, mp
      real dxn, dyn, dxp, dyp, amx, amy, dx, dy, ax, ay
      real px, py, edgelx, edgely, edgerx, edgery
      if (nop.lt.1) return
c begin first particle
      nnn = part(1,1)
      mmn = part(2,1)
      dxn = part(1,1) - float(nnn)
      dyn = part(2,1) - float(mmn)
      nop1 = nop - 1
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
c find vector potential at t+dt/2
      ax = dyp*(dxp*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxp*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyp*(dxp*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxp*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,j) - ax
      py = part(6,j) - ay
c new position at time t + dt
      dx = part(1,j) + (px - 0.5*part(3,j))*dt
      dy = part(2,j) + (py - 0.5*part(4,j))*dt
      part(3,j) = px
      part(4,j) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,j) = -part(6,j)
            part(4,j) = -part(4,j)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,j) = -part(5,j)
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
c find vector potential at t+dt/2
      ax = dyn*(dxn*axy(1,mp+1) + amx*axy(1,mp)) + amy*(dxn*axy(1,mm+1)
     1+ amx*axy(1,mm))
      ay = dyn*(dxn*axy(2,mp+1) + amx*axy(2,mp)) + amy*(dxn*axy(2,mm+1)
     1+ amx*axy(2,mm))
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
c velocity at time t + dt/2
      px = part(5,nop) - ax
      py = part(6,nop) - ay
c new position at time t + dt
      dx = part(1,nop) + (px - 0.5*part(3,nop))*dt
      dy = part(2,nop) + (py - 0.5*part(4,nop))*dt
      part(3,nop) = px
      part(4,nop) = py
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
            part(3,nop) = -part(3,nop)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(6,nop) = -part(6,nop)
            part(4,nop) = -part(4,nop)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(5,nop) = -part(5,nop)
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
      subroutine CR8PCM23GL(part,axy,sctx,qbm,sx,sy,sz,idimp,nop,nx,ny,n
     1xvh,nyv)
c for 2-1/2d code, this subroutine calculates the canonical momentum of
c a particle from the vector potential at t - dt/2, and adds this
c information to the particle array, 
c using gridless spectral version with periodic boundaries,
c canonical momentum p(t-dt/2) = v(t-dt/2) + qi*A(t-dt/2)/(c*mi)
c axy(i,x(t),y(t)) are calculated from the expression
c axy(i,x,y) = sum(axy(i,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c on entry
c part(6,n) = position x of particle n at t - dt/2
c part(7,n) = position y of particle n at t - dt/2
c on exit
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c axy(i,j,k) = i component of vector potential at fourier grid point n,m
c qbm = particle charge/mass ratio
c sx/sy/sz = x/y/z components of field momentum
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv
      real part, qbm, sx, sy, sz
      complex axy, sctx
      dimension part(idimp,nop), axy(3,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3
      real ax, ay, az
      double precision axl, ayl, azl, ox, oy, oz
      double precision sum1, sum2, sum3
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 50 i = 1, nop
c find vector potential at t - dt/2
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(6,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ox = 0.0d0
      oy = 0.0d0
      oz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(7,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      axl = axl + (real(axy(1,j,k)*zt3) + real(axy(1,j,k1)*zt4))
      ayl = ayl + (real(axy(2,j,k)*zt3) + real(axy(2,j,k1)*zt4))
      azl = azl + (real(axy(3,j,k)*zt3) + real(axy(3,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      axl = axl + (real(axy(1,1,k)*zt1) + real(axy(1,1,k1)*zt2))
      ayl = ayl + (real(axy(2,1,k)*zt1) + real(axy(2,1,k1)*zt2))
      azl = azl + (real(axy(3,1,k)*zt1) + real(axy(3,1,k1)*zt2))
      ox = ox + axl
      oy = oy + ayl
      oz = oz + azl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(7,i)
      at1 = cos(dky)
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      axl = axl + (real(axy(1,j,1)*zt3) + real(axy(1,j,k1)*zt1))
      ayl = ayl + (real(axy(2,j,1)*zt3) + real(axy(2,j,k1)*zt1))
      azl = azl + (real(axy(3,j,1)*zt3) + real(axy(3,j,k1)*zt1))
   40 continue
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      oz = 2.0d0*(oz + azl)
      at3 = real(sctx(nxh))
      ox = ox + (real(axy(1,1,1)) + real(axy(1,1,k1))*at1)
      oy = oy + (real(axy(2,1,1)) + real(axy(2,1,k1))*at1)
      oz = oz + (real(axy(3,1,1)) + real(axy(3,1,k1))*at1)
      at1 = at1*at3
      ax = ox + (aimag(axy(1,1,1))*at3 + aimag(axy(1,1,k1))*at1)
      ay = oy + (aimag(axy(2,1,1))*at3 + aimag(axy(2,1,k1))*at1)
      az = oz + (aimag(axy(3,1,1))*at3 + aimag(axy(3,1,k1))*at1)
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      sum1 = sum1 + ax
      sum2 = sum2 + ay
      sum3 = sum3 + az
      part(6,i) = part(3,i) + ax
      part(7,i) = part(4,i) + ay
      part(8,i) = part(5,i) + az
   50 continue
      sx = sum1
      sy = sum2
      sz = sum3
      return
      end
      subroutine GHJPOST2GL(part,fxy,axy,daxy,cu,sctx,qm,qbm,dt,idimp,no
     1p,nx,ny,nxvh,nyv)
c for 2-1/2d code, this subroutine calculates current density,
c using gridless spectral version with periodic boundaries,
c baseline scalar version
c 124*(NX/2)*(NY/2) + 68 flops/particle,
c 28*NX*NY/2 loads and 6*NX*NY/2 stores
c input: all, output: cu
c current density is calculated from the expression:
c cu(i,n,m) = sum(qci*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                  exp(-sqrt(-1)*2*m*pi*y/ny))
c where qci = qm*vj, where j = x,y,z, for i = 1, 3
c and where vj = .5*(vj(t+dt/2)+vj(t-dt/2))
c momentum equations at t=t+dt/2 are calculated from:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are calculated from the expression
c fxy(i,x,y) = sum(fxy(i,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for axy(i,x,y), and daxy(i,x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c axy(i,j,k) = i component of vector potential at fourier grid point n,m
c daxy(i,j,k) = i component of derivative of vector potential at fourier
c grid point n,m
c cu(i,n,m) = current density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv
      real part, qm, qbm, dt
      complex fxy, axy, daxy, cu, sctx
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), axy(3,nxvh,nyv), daxy(5,nxvh,nyv)
      dimension cu(3,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3, qtmh, qmn
      real dx, dy, dz, ax, ay, az, anorm
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      double precision exl, eyl, ezl, ex, ey, ez
      double precision axl, ayl, azl, ox, oy, oz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qtmh = 0.5*qbm*dt
      qmn = qm/real(nx*ny)
      do 120 i = 1, nop
c find fields at t
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
c find vector potential gradients
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      ox = 0.0d0
      oy = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(daxy(1,j,k)*zt3) + real(daxy(1,j,k1)*zt4))
      eyl = eyl + (real(daxy(2,j,k)*zt3) + real(daxy(2,j,k1)*zt4))
      ezl = ezl + (real(daxy(3,j,k)*zt3) + real(daxy(3,j,k1)*zt4))
      axl = axl + (real(daxy(4,j,k)*zt3) + real(daxy(4,j,k1)*zt4))
      ayl = ayl + (real(daxy(5,j,k)*zt3) + real(daxy(5,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      exl = exl + (real(daxy(1,1,k)*zt1) + real(daxy(1,1,k1)*zt2))
      eyl = eyl + (real(daxy(2,1,k)*zt1) + real(daxy(2,1,k1)*zt2))
      ezl = ezl + (real(daxy(3,1,k)*zt1) + real(daxy(3,1,k1)*zt2))
      axl = axl + (real(daxy(4,1,k)*zt1) + real(daxy(4,1,k1)*zt2))
      ayl = ayl + (real(daxy(5,1,k)*zt1) + real(daxy(5,1,k1)*zt2))
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      ox = ox + axl
      oy = oy + ayl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      exl = exl + (real(daxy(1,j,1)*zt3) + real(daxy(1,j,k1)*zt1))
      eyl = eyl + (real(daxy(2,j,1)*zt3) + real(daxy(2,j,k1)*zt1))
      ezl = ezl + (real(daxy(3,j,1)*zt3) + real(daxy(3,j,k1)*zt1))
      axl = axl + (real(daxy(4,j,1)*zt3) + real(daxy(4,j,k1)*zt1))
      ayl = ayl + (real(daxy(5,j,1)*zt3) + real(daxy(5,j,k1)*zt1))
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      at3 = real(sctx(nxh))
      ex = ex + (real(daxy(1,1,1)) + real(daxy(1,1,k1))*at1)
      ey = ey + (real(daxy(2,1,1)) + real(daxy(2,1,k1))*at1)
      ez = ez + (real(daxy(3,1,1)) + real(daxy(3,1,k1))*at1)
      ox = ox + (real(daxy(4,1,1)) + real(daxy(4,1,k1))*at1)
      oy = oy + (real(daxy(5,1,1)) + real(daxy(5,1,k1))*at1)
      at1 = at1*at3
      oxx = ex + (aimag(daxy(1,1,1))*at3 + aimag(daxy(1,1,k1))*at1)
      oyx = ey + (aimag(daxy(2,1,1))*at3 + aimag(daxy(2,1,k1))*at1)
      oxy = ez + (aimag(daxy(3,1,1))*at3 + aimag(daxy(3,1,k1))*at1)
      ozx = ox + (aimag(daxy(4,1,1))*at3 + aimag(daxy(4,1,k1))*at1)
      ozy = oy + (aimag(daxy(5,1,1))*at3 + aimag(daxy(5,1,k1))*at1)
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field and vector potential
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      ox = 0.0d0
      oy = 0.0d0
      oz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 50 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      axl = axl + (real(axy(1,j,k)*zt3) + real(axy(1,j,k1)*zt4))
      ayl = ayl + (real(axy(2,j,k)*zt3) + real(axy(2,j,k1)*zt4))
      azl = azl + (real(axy(3,j,k)*zt3) + real(axy(3,j,k1)*zt4))
   50 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
      eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
      axl = axl + (real(axy(1,1,k)*zt1) + real(axy(1,1,k1)*zt2))
      ayl = ayl + (real(axy(2,1,k)*zt1) + real(axy(2,1,k1)*zt2))
      azl = azl + (real(axy(3,1,k)*zt1) + real(axy(3,1,k1)*zt2))
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      ox = ox + axl
      oy = oy + ayl
      oz = oz + azl
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 70 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
      eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
      ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
      axl = axl + (real(axy(1,j,1)*zt3) + real(axy(1,j,k1)*zt1))
      ayl = ayl + (real(axy(2,j,1)*zt3) + real(axy(2,j,k1)*zt1))
      azl = azl + (real(axy(3,j,1)*zt3) + real(axy(3,j,k1)*zt1))
   70 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      oz = 2.0d0*(oz + azl)
      at3 = real(sctx(nxh))
      ex = ex + (real(fxy(1,1,1)) + real(fxy(1,1,k1))*at1)
      ey = ey + (real(fxy(2,1,1)) + real(fxy(2,1,k1))*at1)
      ez = ez + (real(fxy(3,1,1)) + real(fxy(3,1,k1))*at1)
      ox = ox + (real(axy(1,1,1)) + real(axy(1,1,k1))*at1)
      oy = oy + (real(axy(2,1,1)) + real(axy(2,1,k1))*at1)
      oz = oz + (real(axy(3,1,1)) + real(axy(3,1,k1))*at1)
      at1 = at1*at3
      dx = ex + (aimag(fxy(1,1,1))*at3 + aimag(fxy(1,1,k1))*at1)
      dy = ey + (aimag(fxy(2,1,1))*at3 + aimag(fxy(2,1,k1))*at1)
      dz = ez + (aimag(fxy(3,1,1))*at3 + aimag(fxy(3,1,k1))*at1)
      ax = ox + (aimag(axy(1,1,1))*at3 + aimag(axy(1,1,k1))*at1)
      ay = oy + (aimag(axy(2,1,1))*at3 + aimag(axy(2,1,k1))*at1)
      az = oz + (aimag(axy(3,1,1))*at3 + aimag(axy(3,1,k1))*at1)
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,i)
      py = part(7,i)
      pz = part(8,i)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
c deposit current density
      px = qmn*px
      py = qmn*py
      pz = qmn*pz
      do 80 j = 1, nxh
      sctx(j) = conjg(sctx(j))
   80 continue
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 100 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),-sin(dky))
      zt2 = conjg(zt1)
      do 90 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      cu(1,j,k) = cu(1,j,k) + px*zt3
      cu(2,j,k) = cu(2,j,k) + py*zt3
      cu(3,j,k) = cu(3,j,k) + pz*zt3
      cu(1,j,k1) = cu(1,j,k1) + px*zt4
      cu(2,j,k1) = cu(2,j,k1) + py*zt4
      cu(3,j,k1) = cu(3,j,k1) + pz*zt4
   90 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = at3*zt2
      cu(1,1,k) = cu(1,1,k) + px*zt1
      cu(2,1,k) = cu(2,1,k) + py*zt1
      cu(3,1,k) = cu(3,1,k) + pz*zt1
      cu(1,1,k1) = cu(1,1,k1) + px*zt2
      cu(2,1,k1) = cu(2,1,k1) + py*zt2
      cu(3,1,k1) = cu(3,1,k1) + pz*zt2
  100 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      do 110 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      cu(1,j,1) = cu(1,j,1) + px*zt3
      cu(2,j,1) = cu(2,j,1) + py*zt3
      cu(3,j,1) = cu(3,j,1) + pz*zt3
      cu(1,j,k1) = cu(1,j,k1) + px*zt1
      cu(2,j,k1) = cu(2,j,k1) + py*zt1
      cu(3,j,k1) = cu(3,j,k1) + pz*zt1
  110 continue
      at3 = real(sctx(nxh))
      cu(1,1,1) = cmplx(real(cu(1,1,1))+px,aimag(cu(1,1,1))+px*at3)
      cu(2,1,1) = cmplx(real(cu(2,1,1))+py,aimag(cu(2,1,1))+py*at3)
      cu(3,1,1) = cmplx(real(cu(3,1,1))+pz,aimag(cu(3,1,1))+pz*at3)
      at3 = at1*at3
      cu(1,1,k1) = cmplx(real(cu(1,1,k1))+px*at1,aimag(cu(1,1,k1))+px*at
     13)
      cu(2,1,k1) = cmplx(real(cu(2,1,k1))+py*at1,aimag(cu(2,1,k1))+py*at
     13)
      cu(3,1,k1) = cmplx(real(cu(3,1,k1))+pz*at1,aimag(cu(3,1,k1))+pz*at
     13)
  120 continue
      return
      end
      subroutine GHPUSH23GL(part,fxy,axy,daxy,sctx,qbm,dt,ek,idimp,nop,n
     1x,ny,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates,
c canonical momenta, and velocities using leap-frog scheme in time with
c gridless spectral version with periodic boundaries,
c other boundary conditions are only approximate.
c baseline scalar version
c 100*(NX/2)*(NY/2) + 1 divide + 75 flops/particle,
c 22*NX*NY/2 + 5 loads and 8 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fx*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*Fx*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*Fy*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*Fz*dt) + .5*(q/m)*Fy(x(t),y(t))*dt)
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*Fz*dt
c where q/m is charge/mass, and the "rotation matrix" is given by:
c    rot(1) = ((1 + omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(2) = 2*omyx/(1 - omxx**2 - omyx*omxy)
c    rot(3) = 2*(omzx*(1 + omxx) + omyx*omzy)/(1 - omxx**2 - omyx*omxy)
c    rot(4) = 2*omxy/(1 - omxx**2 - omyx*omxy)
c    rot(5) = ((1 - omxx)**2 + omyx*omxy)/(1 - omxx**2 - omyx*omxy)
c    rot(6) = 2*(omzy*(1 - omxx) + omxy*omzx)/(1 - omxx**2 - omyx*omxy)
c omxx = .5*(q/m)*(dAx/dx)/c, omyx = .5*(q/m)*(dAy/dx)/c,
c omxy = .5*(q/m)*(dAx/dy)/c, omzx = .5*(q/m)*(dAz/dx)/c,
c omzy = .5*(q/m)*(dAz/dy)/c
c Fx = Ex - 2*(Ax*omxx + Ay*omyx + Az*omzx)/c
c Fy = Ey - 2*(Ax*omxy - Ay*omxx + Az*omzy)/c
c Fz = Ez
c Ex = fxy(1,x(t),y(t)), Ey = fxy(2,x(t),y(t)), Ez = fxy(3,x(t),y(t)),
c Ax = axy(1,x(t),y(t)), Ay = axy(2,x(t),y(t)), Az = axy(3,x(t),y(t)),
c dAx/x = daxy(1,x(t),y(t)), dAy/dx = daxy(2,x(t),y(t)),
c dAx/y = daxy(3,x(t),y(t)), dAz/dx = daxy(4,x(t),y(t)),
c dAz/y = daxy(5,x(t),y(t))
c position equations used are:
c x(t+dt/2) = x(t) + .5*vx(t)*dt
c y(t+dt/2) = y(t) + .5*vy(t)*dt
c where vj(t) = .5*(pj(t+dt/2)+pj(t-dt/2)) - (q/m)*A(t)/c
c fxy(i,x(t),y(t)), axy(i,x(t),y(t)), and daxy(i,x(t),y(t))
c are calculated from the expression
c fxy(i,x,y) = sum(fxy(i,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c similarly for axy(i,x,y), and daxy(i,x,y)
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = velocity vx of particle n at t - dt/2
c part(4,n) = velocity vy of particle n at t - dt/2
c part(5,n) = velocity vz of particle n at t - dt/2
c part(6,n) = canonical momentum px of particle n at t - dt/2
c part(7,n) = canonical momentum py of particle n at t - dt/2
c part(8,n) = canonical momentum pz of particle n at t - dt/2
c fxy(i,n,m) = i component of force/charge at fourier grid point n,m
c axy(i,j,k) = i component of vector potential at fourier grid point n,m
c daxy(i,j,k) = i component of derivative of vector potential at fourier
c grid point n,m
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum(vx(t)**2 + vy(t)**2 + vz(t)**2)
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt, ek
      complex fxy, axy, daxy, sctx
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), axy(3,nxvh,nyv), daxy(5,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3, dth, qtmh
      real dx, dy, dz, ax, ay, az, anorm, edgelx, edgely, edgerx, edgery
      real oxx, oyx, ozx, oxy, ozy, omxxt, omyxt, omzxt, omxyt, omzyt
      real px, py, pz, acx, acy, acz, rot1, rot2, rot3, rot4, rot5, rot6
      double precision exl, eyl, ezl, ex, ey, ez
      double precision axl, ayl, azl, ox, oy, oz
      double precision sum1
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dth = 0.5*dt
      qtmh = qbm*dth
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
      do 80 i = 1, nop
c find fields at t
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
c find vector potential gradients
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      ox = 0.0d0
      oy = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(daxy(1,j,k)*zt3) + real(daxy(1,j,k1)*zt4))
      eyl = eyl + (real(daxy(2,j,k)*zt3) + real(daxy(2,j,k1)*zt4))
      ezl = ezl + (real(daxy(3,j,k)*zt3) + real(daxy(3,j,k1)*zt4))
      axl = axl + (real(daxy(4,j,k)*zt3) + real(daxy(4,j,k1)*zt4))
      ayl = ayl + (real(daxy(5,j,k)*zt3) + real(daxy(5,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      exl = exl + (real(daxy(1,1,k)*zt1) + real(daxy(1,1,k1)*zt2))
      eyl = eyl + (real(daxy(2,1,k)*zt1) + real(daxy(2,1,k1)*zt2))
      ezl = ezl + (real(daxy(3,1,k)*zt1) + real(daxy(3,1,k1)*zt2))
      axl = axl + (real(daxy(4,1,k)*zt1) + real(daxy(4,1,k1)*zt2))
      ayl = ayl + (real(daxy(5,1,k)*zt1) + real(daxy(5,1,k1)*zt2))
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      ox = ox + axl
      oy = oy + ayl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      exl = exl + (real(daxy(1,j,1)*zt3) + real(daxy(1,j,k1)*zt1))
      eyl = eyl + (real(daxy(2,j,1)*zt3) + real(daxy(2,j,k1)*zt1))
      ezl = ezl + (real(daxy(3,j,1)*zt3) + real(daxy(3,j,k1)*zt1))
      axl = axl + (real(daxy(4,j,1)*zt3) + real(daxy(4,j,k1)*zt1))
      ayl = ayl + (real(daxy(5,j,1)*zt3) + real(daxy(5,j,k1)*zt1))
   40 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      at3 = real(sctx(nxh))
      ex = ex + (real(daxy(1,1,1)) + real(daxy(1,1,k1))*at1)
      ey = ey + (real(daxy(2,1,1)) + real(daxy(2,1,k1))*at1)
      ez = ez + (real(daxy(3,1,1)) + real(daxy(3,1,k1))*at1)
      ox = ox + (real(daxy(4,1,1)) + real(daxy(4,1,k1))*at1)
      oy = oy + (real(daxy(5,1,1)) + real(daxy(5,1,k1))*at1)
      at1 = at1*at3
      oxx = ex + (aimag(daxy(1,1,1))*at3 + aimag(daxy(1,1,k1))*at1)
      oyx = ey + (aimag(daxy(2,1,1))*at3 + aimag(daxy(2,1,k1))*at1)
      oxy = ez + (aimag(daxy(3,1,1))*at3 + aimag(daxy(3,1,k1))*at1)
      ozx = ox + (aimag(daxy(4,1,1))*at3 + aimag(daxy(4,1,k1))*at1)
      ozy = oy + (aimag(daxy(5,1,1))*at3 + aimag(daxy(5,1,k1))*at1)
c calculate "cyclotron" frequency
      omxxt = qtmh*oxx
      omyxt = qtmh*oyx
      omxyt = qtmh*oxy
      omzxt = qtmh*ozx
      omzyt = qtmh*ozy
c find electric field and vector potential
      ex = 0.0d0
      ey = 0.0d0
      ez = 0.0d0
      ox = 0.0d0
      oy = 0.0d0
      oz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 60 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 50 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      exl = exl + (real(fxy(1,j,k)*zt3) + real(fxy(1,j,k1)*zt4))
      eyl = eyl + (real(fxy(2,j,k)*zt3) + real(fxy(2,j,k1)*zt4))
      ezl = ezl + (real(fxy(3,j,k)*zt3) + real(fxy(3,j,k1)*zt4))
      axl = axl + (real(axy(1,j,k)*zt3) + real(axy(1,j,k1)*zt4))
      ayl = ayl + (real(axy(2,j,k)*zt3) + real(axy(2,j,k1)*zt4))
      azl = azl + (real(axy(3,j,k)*zt3) + real(axy(3,j,k1)*zt4))
   50 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      exl = exl + (real(fxy(1,1,k)*zt1) + real(fxy(1,1,k1)*zt2))
      eyl = eyl + (real(fxy(2,1,k)*zt1) + real(fxy(2,1,k1)*zt2))
      ezl = ezl + (real(fxy(3,1,k)*zt1) + real(fxy(3,1,k1)*zt2))
      axl = axl + (real(axy(1,1,k)*zt1) + real(axy(1,1,k1)*zt2))
      ayl = ayl + (real(axy(2,1,k)*zt1) + real(axy(2,1,k1)*zt2))
      azl = azl + (real(axy(3,1,k)*zt1) + real(axy(3,1,k1)*zt2))
      ex = ex + exl
      ey = ey + eyl
      ez = ez + ezl
      ox = ox + axl
      oy = oy + ayl
      oz = oz + azl
   60 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      exl = 0.0d0
      eyl = 0.0d0
      ezl = 0.0d0
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 70 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      exl = exl + (real(fxy(1,j,1)*zt3) + real(fxy(1,j,k1)*zt1))
      eyl = eyl + (real(fxy(2,j,1)*zt3) + real(fxy(2,j,k1)*zt1))
      ezl = ezl + (real(fxy(3,j,1)*zt3) + real(fxy(3,j,k1)*zt1))
      axl = axl + (real(axy(1,j,1)*zt3) + real(axy(1,j,k1)*zt1))
      ayl = ayl + (real(axy(2,j,1)*zt3) + real(axy(2,j,k1)*zt1))
      azl = azl + (real(axy(3,j,1)*zt3) + real(axy(3,j,k1)*zt1))
   70 continue
      ex = 2.0d0*(ex + exl)
      ey = 2.0d0*(ey + eyl)
      ez = 2.0d0*(ez + ezl)
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      oz = 2.0d0*(oz + azl)
      at3 = real(sctx(nxh))
      ex = ex + (real(fxy(1,1,1)) + real(fxy(1,1,k1))*at1)
      ey = ey + (real(fxy(2,1,1)) + real(fxy(2,1,k1))*at1)
      ez = ez + (real(fxy(3,1,1)) + real(fxy(3,1,k1))*at1)
      ox = ox + (real(axy(1,1,1)) + real(axy(1,1,k1))*at1)
      oy = oy + (real(axy(2,1,1)) + real(axy(2,1,k1))*at1)
      oz = oz + (real(axy(3,1,1)) + real(axy(3,1,k1))*at1)
      at1 = at1*at3
      dx = ex + (aimag(fxy(1,1,1))*at3 + aimag(fxy(1,1,k1))*at1)
      dy = ey + (aimag(fxy(2,1,1))*at3 + aimag(fxy(2,1,k1))*at1)
      dz = ez + (aimag(fxy(3,1,1))*at3 + aimag(fxy(3,1,k1))*at1)
      ax = ox + (aimag(axy(1,1,1))*at3 + aimag(axy(1,1,k1))*at1)
      ay = oy + (aimag(axy(2,1,1))*at3 + aimag(axy(2,1,k1))*at1)
      az = oz + (aimag(axy(3,1,1))*at3 + aimag(axy(3,1,k1))*at1)
c calculate half impulse
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
      dx = qtmh*dx - (ax*omxxt + ay*omyxt + az*omzxt)
      dy = qtmh*dy - (ax*omxyt - ay*omxxt + az*omzyt)
      dz = qtmh*dz
c half acceleration
      px = part(6,i)
      py = part(7,i)
      pz = part(8,i)
      acx = px + dx
      acy = py + dy
      acz = pz + dz
c calculate "rotation" matrix
      rot5 = omyxt*omxyt
      rot3 = 1.0 + omxxt
      rot6 = 1.0 - omxxt
      rot2 = omyxt + omyxt
      rot4 = omxyt + omxyt
      anorm = 1.0/(1.0 - omxxt*omxxt - rot5)
      rot1 = rot3*rot3 + rot5
      rot5 = rot6*rot6 + rot5
      rot3 = 2.0*(omzxt*rot3 + omyxt*omzyt)
      rot6 = 2.0*(omzyt*rot6 + omxyt*omzxt)
c new canonical momentum at t + dt/2
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = acz + dz
      part(6,i) = dx
      part(7,i) = dy
      part(8,i) = dz
c velocity at time t
      px = 0.5*(px + dx) - ax
      py = 0.5*(py + dy) - ay
      pz = 0.5*(pz + dz) - az
      part(3,i) = px
      part(4,i) = py
      part(5,i) = pz
c time-centered kinetic energy
      sum1 = sum1 + (px*px + py*py + pz*pz)
c new position at time t + dt/2
      dx = part(1,i) + px*dth
      dy = part(2,i) + py*dth
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
            part(6,i) = -part(6,i)
            part(3,i) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i)
            part(7,i) = -part(7,i)
            part(4,i) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i)
            part(6,i) = -part(6,i)
            part(3,i) = 0.0
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
c set new position
      part(1,i) = dx
      part(2,i) = dy
   80 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
      subroutine GHPUSHXH23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxvh,
     1nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time with
c gridless spectral version with periodic boundaries,
c other boundary conditions are only approximate.
c baseline scalar version
c 27*(NX/2)*(NY/2) + 12 flops/particle,
c 3*NX*NY/2 + 8 loads and 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + .5*(vxnew(t+dt/2) - vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + .5*(vynew(t+dt/2) - vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) are calculated from the expression
c axy(i,x,y) = sum(axy(i,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j,k) = i component of vector potential at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt
      complex axy, sctx
      dimension part(idimp,nop), axy(3,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3, dth
      real dx, dy, ax, ay, az, edgelx, edgely, edgerx, edgery
      real px, py, pz
      double precision axl, ayl, azl, ox, oy, oz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dth = 0.5*dt
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
c find vector potential at t+dt/2
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ox = 0.0d0
      oy = 0.0d0
      oz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      axl = axl + (real(axy(1,j,k)*zt3) + real(axy(1,j,k1)*zt4))
      ayl = ayl + (real(axy(2,j,k)*zt3) + real(axy(2,j,k1)*zt4))
      azl = azl + (real(axy(3,j,k)*zt3) + real(axy(3,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      axl = axl + (real(axy(1,1,k)*zt1) + real(axy(1,1,k1)*zt2))
      ayl = ayl + (real(axy(2,1,k)*zt1) + real(axy(2,1,k1)*zt2))
      azl = azl + (real(axy(3,1,k)*zt1) + real(axy(3,1,k1)*zt2))
      ox = ox + axl
      oy = oy + ayl
      oz = oz + azl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      axl = axl + (real(axy(1,j,1)*zt3) + real(axy(1,j,k1)*zt1))
      ayl = ayl + (real(axy(2,j,1)*zt3) + real(axy(2,j,k1)*zt1))
      azl = azl + (real(axy(3,j,1)*zt3) + real(axy(3,j,k1)*zt1))
   40 continue
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      oz = 2.0d0*(oz + azl)
      at3 = real(sctx(nxh))
      ox = ox + (real(axy(1,1,1)) + real(axy(1,1,k1))*at1)
      oy = oy + (real(axy(2,1,1)) + real(axy(2,1,k1))*at1)
      oz = oz + (real(axy(3,1,1)) + real(axy(3,1,k1))*at1)
      at1 = at1*at3
      ax = ox + (aimag(axy(1,1,1))*at3 + aimag(axy(1,1,k1))*at1)
      ay = oy + (aimag(axy(2,1,1))*at3 + aimag(axy(2,1,k1))*at1)
      az = oz + (aimag(axy(3,1,1))*at3 + aimag(axy(3,1,k1))*at1)
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,i) - ax
      py = part(7,i) - ay
      pz = part(8,i) - az
c new position at time t + dt/2
      dx = part(1,i) + (px - part(3,i))*dth
      dy = part(2,i) + (py - part(4,i))*dth
      part(3,i) = px
      part(4,i) = py
      part(5,i) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,i) = -part(6,i)
            part(3,i) = 0.0
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dth
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,i) = -part(7,i)
            part(4,i) = 0.0
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dth
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,i) = -part(6,i)
            part(3,i) = 0.0
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
      subroutine GHPUSHX23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxvh,n
     1yv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time with
c gridless spectral version with periodic boundaries,
c other boundary conditions are only approximate.
c baseline scalar version
c 27*(NX/2)*(NY/2) + 12 flops/particle,
c 3*NX*NY/2 + 8 loads and 5 stores
c input: all, output: part
c position equations used are:
c x(t+dt/2) = x(t+dt/2) + (vxnew(t+dt/2) - .5*vxold(t+dt/2))*dt
c y(t+dt/2) = y(t+dt/2) + (vynew(t+dt/2) - .5*vyold(t+dt/2))*dt
c where vj(t+dt/2) = pj(t+dt/2) - (q/m)*A(t+dt/2)/c
c axy(i,x(t),y(t)) are calculated from the expression
c axy(i,x,y) = sum(axy(i,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n at t + dt/2
c part(2,n) = position y of particle n at t + dt/2
c part(3,n) = velocity vx of particle n at t + dt/2
c part(4,n) = velocity vy of particle n at t + dt/2
c part(5,n) = velocity vz of particle n at t + dt/2
c part(6,n) = canonical momentum px of particle n at t + dt/2
c part(7,n) = canonical momentum py of particle n at t + dt/2
c part(8,n) = canonical momentum pz of particle n at t + dt/2
c axy(i,j,k) = i component of vector potential at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c idimp = size of phase space = 8
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt
      complex axy, sctx
      dimension part(idimp,nop), axy(3,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real dnx, dny, dkx, dky, at1, at3
      real dx, dy, ax, ay, az, edgelx, edgely, edgerx, edgery
      real px, py, pz
      double precision axl, ayl, azl, ox, oy, oz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
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
c find vector potential at t+dt/2
      do 10 j = 1, nxh
      dkx = dnx*real(j)*part(1,i)
      sctx(j) = cmplx(cos(dkx),sin(dkx))
   10 continue
      ox = 0.0d0
      oy = 0.0d0
      oz = 0.0d0
c mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k-1)*part(2,i)
      zt1 = cmplx(cos(dky),sin(dky))
      zt2 = conjg(zt1)
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 20 j = 2, nxh
      zt3 = sctx(j-1)
      zt4 = zt2*zt3
      zt3 = zt1*zt3
      axl = axl + (real(axy(1,j,k)*zt3) + real(axy(1,j,k1)*zt4))
      ayl = ayl + (real(axy(2,j,k)*zt3) + real(axy(2,j,k1)*zt4))
      azl = azl + (real(axy(3,j,k)*zt3) + real(axy(3,j,k1)*zt4))
   20 continue
c mode numbers kx = 0, nx/2
      at3 = real(sctx(nxh))
      zt2 = zt2*at3
      axl = axl + (real(axy(1,1,k)*zt1) + real(axy(1,1,k1)*zt2))
      ayl = ayl + (real(axy(2,1,k)*zt1) + real(axy(2,1,k1)*zt2))
      azl = azl + (real(axy(3,1,k)*zt1) + real(axy(3,1,k1)*zt2))
      ox = ox + axl
      oy = oy + ayl
      oz = oz + azl
   30 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      dky = dny*real(nyh)*part(2,i)
      at1 = cos(dky)
      axl = 0.0d0
      ayl = 0.0d0
      azl = 0.0d0
      do 40 j = 2, nxh
      zt3 = sctx(j-1)
      zt1 = at1*zt3
      axl = axl + (real(axy(1,j,1)*zt3) + real(axy(1,j,k1)*zt1))
      ayl = ayl + (real(axy(2,j,1)*zt3) + real(axy(2,j,k1)*zt1))
      azl = azl + (real(axy(3,j,1)*zt3) + real(axy(3,j,k1)*zt1))
   40 continue
      ox = 2.0d0*(ox + axl)
      oy = 2.0d0*(oy + ayl)
      oz = 2.0d0*(oz + azl)
      at3 = real(sctx(nxh))
      ox = ox + (real(axy(1,1,1)) + real(axy(1,1,k1))*at1)
      oy = oy + (real(axy(2,1,1)) + real(axy(2,1,k1))*at1)
      oz = oz + (real(axy(3,1,1)) + real(axy(3,1,k1))*at1)
      at1 = at1*at3
      ax = ox + (aimag(axy(1,1,1))*at3 + aimag(axy(1,1,k1))*at1)
      ay = oy + (aimag(axy(2,1,1))*at3 + aimag(axy(2,1,k1))*at1)
      az = oz + (aimag(axy(3,1,1))*at3 + aimag(axy(3,1,k1))*at1)
c calculate total canonical momentum
      ax = qbm*ax
      ay = qbm*ay
      az = qbm*az
c velocity at time t + dt/2
      px = part(6,i) - ax
      py = part(7,i) - ay
      pz = part(8,i) - az
c new position at time t + dt
      dx = part(1,i) + (px - 0.5*part(3,i))*dt
      dy = part(2,i) + (py - 0.5*part(4,i))*dt
      part(3,i) = px
      part(4,i) = py
      part(5,i) = pz
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,i) = -part(6,i)
            part(3,i) = -part(3,i)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = dy - py*dt
            if (dy.lt.edgely) then
               dy = edgely
            elseif (dy.ge.edgery) then
               dy = 0.999999*edgery
            endif
            part(7,i) = -part(7,i)
            part(4,i) = -part(4,i)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = dx - px*dt
            if (dx.lt.edgelx) then
               dx = edgelx
            elseif (dx.ge.edgerx) then
               dx = 0.999999*edgerx
            endif
            part(6,i) = -part(6,i)
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
