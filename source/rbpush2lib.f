c-----------------------------------------------------------------------
c 2d PIC library for pushing relativistic particles with magnetic field
c and depositing current
c rbpush2lib.f contains procedures to process relativistic particles
c              with and without magnetic fields:
c GRJPOST2 deposits 3 component current density, quadratic
c          interpolation, STANDARD optimization, for relativistic
c          particles.
c GSRJPOST2 deposits 3 component current density, quadratic
c           interpolation, LOOKAHEAD optimization, for relativistic
c           particles.
c GSRJPOST2X deposits 3 component current density, quadratic
c            interpolation, VECTOR optimization, for relativistic
c            particles.
c GRJPOST2L deposits 3 component current density, linear interpolation,
c           STANDARD optimization, for relativistic particles.
c GSRJPOST2L deposits 3 component current density, linear interpolation,
c            LOOKAHEAD optimization, for relativistic particles.
c GSRJPOST2XL deposits 3 component current density, linear
c             interpolation, VECTOR optimization, for relativistic
c             particles.
c GRJPOST22 deposits 2 component current density, quadratic
c           interpolation, STANDARD optimization, for relativistic
c           particles.
c GSRJPOST22 deposits 2 component current density, quadratic
c            interpolation, LOOKAHEAD optimization, for relativistic
c            particles.
c GSRJPOST22X deposits 2 component current density, quadratic
c             interpolation, VECTOR optimization, for relativistic
c             particles.
c GRJPOST22L deposits 2 component current density, linear interpolation,
c            STANDARD optimization, for relativistic particles.
c GSRJPOST22L deposits 2 component current density, linear
c             interpolation, LOOKAHEAD optimization, for relativistic
c             particles.
c GSRJPOST22XL deposits 2 component current density, linear
c              interpolation, VECTOR optimization, for relativistic
c              particles.
c GRJPOST2C deposits 3 component current density, cubic interpolation,
c           STANDARD optimization, for relativistic particles.
c GRJPOST22C deposits 2 component current density, cubic interpolation,
c            STANDARD optimization, for relativistic particles.
c GRPUSH2 push particles with 2 component electric field, quadratic
c         interpolation, STANDARD optimization, for relativistic
c         particles.
c GSRPUSH2 push particles with 2 component electric field, quadratic
c          interpolation, LOOKAHEAD optimization, for relativistic
c          particles.
c GRPUSH2L push particles with 2 component electric field, linear
c          interpolation, STANDARD optimization, for relativistic
c          particles.
c GSRPUSH2L push particles with 2 component electric field, linear
c           interpolation, LOOKAHEAD optimization, for relativistic
c           particles.
c GRBPUSH2 push particles with 3 component magnetic field, 2 component
c          electric field, quadratic interpolation, STANDARD
c          optimization, for relativistic particles.
c GSRBPUSH2 push particles with 3 component magnetic field, 2 component
c           electric field, quadratic interpolation, LOOKAHEAD
c           optimization, for relativistic particles.
c GRBPUSH2L push particles with 3 component magnetic field, 2 component
c           electric field, linear interpolation, STANDARD optimization,
c           for relativistic particles.
c GSRBPUSH2L push particles with 3 component magnetic field, 2 component
c            electric field, linear interpolation, LOOKAHEAD
c            optimization, for relativistic particles.
c GRBPUSH23 push particles with 3 component magnetic field, 3 component
c           electric field, quadratic interpolation, STANDARD
c           optimization, for relativistic particles.
c GSRBPUSH23 push particles with 3 component magnetic field, 3 component
c            electric field, quadratic interpolation, LOOKAHEAD
c            optimization, for relativistic particles.
c GRBPUSH23L push particles with 3 component magnetic field, 3 component
c            electric field, linear interpolation, STANDARD
c            optimization, for relativistic particles.
c GSRBPUSH23L push particles with 3 component magnetic field,
c             3 component electric field, linear interpolation,
c             LOOKAHEAD optimization, for relativistic particles.
c GRBPUSH22 push particles with 1 component magnetic field, 2 component
c           electric field, quadratic interpolation, STANDARD
c           optimization, for relativistic particles.
c GSRBPUSH22 push particles with 1 component magnetic field, 2 component
c            electric field, quadratic interpolation, LOOKAHEAD
c            optimization, for relativistic particles.
c GRBPUSH22L push particles with 1 component magnetic field, 2 component
c            electric field, linear interpolation, STANDARD
c            optimization, for relativistic particles.
c GSRBPUSH22L push particles with 1 component magnetic field,
c             2 component electric field, linear interpolation,
c             LOOKAHEAD optimization, for relativistic particles.
c GRPUSH2C push particles with 2 component electric field, cubic
c          interpolation, STANDARD optimization, for relativistic
c          particles.
c GRBPUSH2C push particles with 3 component magnetic field, 2 component
c          electric field, cubic interpolation, STANDARD
c          optimization, for relativistic particles.
c GRBPUSH23C push particles with 3 component magnetic field, 3 component
c            electric field, cubic interpolation, STANDARD
c            optimization, for relativistic particles.
c GRBPUSH22C push particles with 1 component magnetic field, 2 component
c            electric field, cubic interpolation, STANDARD optimization,
c            for relativistic particles.
c RRETARD2 retard particle position a half time-step for 2-1/2d code,
c          and relativistic particles.
c RRETARD22 retard particle position a half time-step for 2d code,
c           and relativistic particles.
c CPTOV2 converts momentum to velocity for relativistic particles for
c        2-1/2d code.
c CPTOV22 converts momentum to velocity for relativistic particles for
c         2d code.
c RPUSH2ZF update particle co-ordinates for particles with fixed
c          velocities, for 2d code, and relativistic particles.
c RPUSH23ZF update particle co-ordinates for particles with fixed
c           velocities, for 2-1/2d code, and relativistic particles.
c RDJPOST2GL deposits 3 component current density, using gridless
c            method, for relativistic particles.
c RPUSH2GL push particles with 2 component electric field, using
c          gridless method, for relativistic particles.
c RPUSH2GLX push particles with 2 component electric field, using
c           optimized gridless method, for relativistic particles.
c RBPUSH23GL push particles with 3 component magnetic field, 3 component
c            electric field, using gridless method, for relativistic
c            particles.
c GRCJPOST2 deposits time-centered particle current density, quadratic
c           interpolation, for relativistic particles.
c GRCJPOST2L deposits time-centered particle current density, linear
c            interpolation, for relativistic particles.
c GRCJPOST2C deposits time-centered particle current density, cubic
c            interpolation, for relativistic particles.
!*****
!  PSPLIT2 splits the real particles to virtual particles that do not
!          cross cell boundaries
!  OGRJPOST2, 2L, 2C deposits the virtual particle current using
!                    quadratic, linear, or cubic interpolation.
!*****
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: september 16, 2011
c-----------------------------------------------------------------------
      subroutine PSPLIT2(part,vpart,nop,dt,nsplit,ci,idimp,nx,
     1ny,ipbc,order)
!*****
!This subroutine splits the real particles into virtual particles
!that do not cross cell boundaries. It is based on the subroutine
!split_2d from OSIRIS.
!*****
      implicit none
      real :: part, vpart
      real :: dt, ci
      integer :: nop, nsplit, idimp, nx, ny, ipbc, order
      dimension part(idimp,nop), vpart(idimp+4,4*nop)
      integer :: ntime
! local vars
      integer :: j, k, dxi, dyi, n0, n1, m0, m1
      real :: x0, x1, vx, vy, vz, px, py, pz, p2, delta, xint, gami
      real :: y0, y1, yint
      real :: vxint, vyint, vzint, xint2, yint2

      k = 0
      do 10 j = 1, nop
                px = part(3,j)
                py = part(4,j)
                pz = part(5,j)
                p2 = px*px + py*py + pz*pz
                gami = 1.0/sqrt(1.0 + p2*ci*ci)
                vx = px*gami
                vy = py*gami
                vz = pz*gami

!x parameters
                if(order==1.or.order==3) then
                  x0 = part(1,j) + 0.5*(3.0 - dt*vx)
                  x1 = part(1,j) + 0.5*(3.0 + dt*vx)
                elseif(order==2) then
                  x0 = part(1,j) + 0.5*(2.0 - dt*vx)
                  x1 = part(1,j) + 0.5*(2.0 + dt*vx)
                endif
! find the nearest grid point for the initial and final positions                
                n0 = nint(x0)
                n1 = nint(x1)
! convert x1, x0 to OSIRIS coordinates so we can directly follow the deposition algorithm
                x0 = x0 - real(n0)
                x1 = x1 - real(n0)
                n0 = n0
                n1 = n1
! follow OSIRIS splitter
                dxi = n1 - n0
!y parameters
                if(order==1.or.order==3) then
                  y0 = part(2,j) + 0.5*(3.0 - dt*vy)
                  y1 = part(2,j) + 0.5*(3.0 + dt*vy)
                elseif(order==2) then
                  y0 = part(2,j) + 0.5*(2.0 - dt*vy)
                  y1 = part(2,j) + 0.5*(2.0 + dt*vy)
                endif

                m0 = nint(y0)
                m1 = nint(y1)
                y0 = y0 - real(m0)
                y1 = y1 - real(m0)
                m0 = m0! + 1
                m1 = m1! + 1
                dyi = m1 - m0

!vpart(:,k) = start x pos, end x pos, start cell x index,
             !start y pos, end y pos, start cell y index, vx, vy, vz
!no splitting
                if(dxi==0.and.dyi==0) then
                  k = k + 1
                  vpart(1,k) = x0
                  vpart(2,k) = x1
                  vpart(3,k) = n0
                  vpart(4,k) = y0
                  vpart(5,k) = y1
                  vpart(6,k) = m0
! OSIRIS uses v (not p) in the virtual particle array.
                  vpart(7,k) = vx
                  vpart(8,k) = vy
                  vpart(9,k) = vz

! x cross only
                elseif((dxi.ne.0).and.(dyi==0)) then
                  xint = 0.5*dxi
		  delta = (xint - x0) / (x1 - x0)
                  yint = y0 + (y1 - y0)*delta            
! parameters for the virtual particle that stays in the first x cell
                  k = k + 1		
	          vpart(1,k) = x0  
	          vpart(2,k) = xint
                  vpart(3,k) = n0
                  vpart(4,k) = y0
                  vpart(5,k) = yint
                  vpart(6,k) = m0
		  vpart(7,k) = vx * delta
		  vpart(8,k) = vy * delta
		  vpart(9,k) = vz * delta
		      
! parameters for the virtual particle that moves to a new x cell
		  k = k + 1
		  vpart(1,k) = - xint      
		  vpart(2,k) = x1 - dxi
                  vpart(3,k) = n1
                  vpart(4,k) = yint
                  vpart(5,k) = y1
                  vpart(6,k) = m0
                  vpart(7,k) = vx * (1-delta)
		  vpart(8,k) = vy * (1-delta)
		  vpart(9,k) = vz * (1-delta)

! y cross only
                elseif((dxi==0).and.(dyi.ne.0)) then
                  yint = 0.5*dyi
		  delta = (yint - y0) / (y1 - y0)
                  xint = x0 + (x1 - x0)*delta         
! parameters for the virtual particle that stays in the first y cell
                  k = k + 1		
		  vpart(1,k) = x0  
		  vpart(2,k) = xint
                  vpart(3,k) = n0
                  vpart(4,k) = y0
                  vpart(5,k) = yint
                  vpart(6,k) = m0
		  vpart(7,k) = vx * delta
		  vpart(8,k) = vy * delta
		  vpart(9,k) = vz * delta
		      
! parameters for the virtual particle that moves to a new y cell
		  k = k + 1
		  vpart(1,k) = xint      
		  vpart(2,k) = x1
                vpart(3,k) = n0
                vpart(4,k) = -yint
                vpart(5,k) = y1 - dyi
                vpart(6,k) = m1
                vpart(7,k) = vx * (1-delta)
		  vpart(8,k) = vy * (1-delta)
		  vpart(9,k) = vz * (1-delta)

!x and y cross
                  elseif((dxi.ne.0).and.(dyi.ne.0)) then
! first consider x direction
                    xint = 0.5*dxi
		    delta = (xint - x0) / (x1 - x0)
                    yint = y0 + (y1 - y0)*delta
                   
                    if((yint >= -0.5).and.(yint < 0.5)) then
! no y cross on 1st vp
			k = k + 1  
			vpart(1,k) = x0
                     vpart(2,k) = xint 
                     vpart(3,k) = n0      
			vpart(4,k) = y0          
			vpart(5,k) = yint
			vpart(6,k) = m0
			vpart(7,k) = vx*delta
                     vpart(8,k) = vy*delta
			vpart(9,k) = vz*delta
      
                     vxint = vx*(1-delta)
                     vyint = vy*(1-delta)
                     vzint = vz*(1-delta)
			
! y split 2nd vp
! first vp from original vp
			k = k + 1

			yint2 = 0.5*dyi 
			delta = ( yint2 - yint ) / ( y1 - yint )
			xint2 = -xint + (x1 - xint)*delta  
			
			vpart(1,k) = -xint
                        vpart(2,k) = xint2
                        vpart(3,k) = n1
			vpart(4,k) = yint	
			vpart(5,k) = yint2 
                        vpart(6,k) = m0 
			vpart(7,k) = vxint*delta
                        vpart(8,k) = vyint*delta
			vpart(9,k) = vzint*delta

!second vp from original vp			 
			k = k + 1
              
                        vpart(1,k) = xint2
                        vpart(2,k) = x1 - dxi
                        vpart(3,k) = n1
			vpart(4,k) = -yint2	
			vpart(5,k) = y1 - dyi 
                        vpart(6,k) = m1 
			vpart(7,k) = vxint*(1-delta)
                        vpart(8,k) = vyint*(1-delta)
			vpart(9,k) = vzint*(1-delta)
			
		 else
!first vp does cross y (and x)		
			vzint = vz * delta
			
			! y split 1st vp
			yint2 = 0.5*dyi
			delta = ( yint2 - y0 ) / ( yint - y0)
			xint2 = x0 + (xint - x0) * delta
			
			k = k + 1
!vp from original vp that doesn't cross y
                        vpart(1,k) = x0
                        vpart(2,k) = xint2
                        vpart(3,k) = n0
			vpart(4,k) = y0	
			vpart(5,k) = yint2 
                        vpart(6,k) = m0 
			vpart(7,k) = vxint*delta
                        vpart(8,k) = vyint*delta
			vpart(9,k) = vzint*delta
			
			k = k + 1
!second vp from original vp that does cross y
                        vpart(1,k) = xint2
                        vpart(2,k) = xint
                        vpart(3,k) = n0
			vpart(4,k) = -yint2
			vpart(5,k) = yint - dyi 
                        vpart(6,k) = m1
			vpart(7,k) = vxint*(1-delta)
                        vpart(8,k) = vyint*(1-delta)
			vpart(9,k) = vzint*(1-delta)
			 
			! no y cross on second vp
                        k = k + 1
                        vpart(1,k) = -xint
                        vpart(2,k) = x1 - dxi
                        vpart(3,k) = n1
			vpart(4,k) = yint - dyi
			vpart(5,k) = y1 - dyi 
                        vpart(6,k) = m1
			vpart(7,k) = vx - vxint
                        vpart(8,k) = vy - vyint
			vpart(9,k) = vz - vzint		
		  endif
	        endif
   10 continue

      nsplit = k

      return
      end
      subroutine GRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 88 flops/particle, 1 divide, 1 sqrt, 32 loads, 29 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+3
c nyv = second dimension of current array, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(3,nxv,nyv)
      qmh = .5*qm
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
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine OGRJPOST2(part,vpart,vcu,qm,dth,ci,nsplit,nop, !!!
     1idimp,nx,ny,nxv,nyv,ipbc,dt)
!*****
!This subroutine deposits the charge conserving current for quadratic
!particle shape into the virtual current array. It is based on the
!OSIRIS subroutine getjr_2d_s2. It also does a half position update, 
!like GRJPOST2.
!*****
      dimension part(idimp,nop), vcu(3,nx+4,ny+4) 
      dimension vpart(idimp+4,4*nop)
      dimension S0x(3), S1x(3), S0y(3), S1y(3)
      dimension wp1(3), wp2(3), wl1(2), wl2(2)
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
      !print*, 'nsplit =', nsplit

      if(ny==1) then
        do 5 j = 1, nsplit
c find interpolation weights   
        x0 = vpart(1,j)
        x1 = vpart(2,j)
        ix = vpart(3,j) + 1
        y0 = vpart(4,j)
        y1 = vpart(5,j)
        iy = vpart(6,j) + 1

        qnx = 1.0/(4.0*dt)*qm
        qvy = vpart(idimp+3,j)*qm/3.0
        qvz = vpart(idimp+4,j)*qm/3.0

	! get spline weitghts for x and y
	  S0x(1) = 0.5*(0.5 - x0)**2
	  S0x(2) = 0.75 - x0**2
	  S0x(3) = 0.5*(0.5 + x0)**2
	
	  S1x(1) = 0.5*(0.5 - x1)**2
	  S1x(2) = 0.75 - x1**2
	  S1x(3) = 0.5*(0.5 + x1)**2
	
	  S0y(1) = 0.5*(0.5 - y0)**2
	  S0y(2) = 0.75 - y0**2
	  S0y(3) = 0.5*(0.5 + y0)**2
	
	  S1y(1) =  0.5*(0.5 - y1)**2
	  S1y(2) = 0.75 - y1**2
	  S1y(3) = 0.5*(0.5 + y1)**2
	
	
	! get longitudinal motion weights
	  wl1(1) = (qnx*(x0 - x1)*(-1 + x0 + x1))
	  wl1(2) = (qnx*(-(x0*(1 + x0)) + x1*(1 + x1)))

	! get perpendicular motion weights
	  wp1(1) = (S0y(1)+S1y(1))
	  wp1(2) = (S0y(2)+S1y(2))
	  wp1(3) = (S0y(3)+S1y(3))
	
	! accumulate j1
	  vcu(1,ix-1,iy-1) = vcu(1,ix-1,iy-1) + wl1(1) * wp1(1)
	  vcu(1,ix  ,iy-1) = vcu(1,ix  ,iy-1) + wl1(2) * wp1(1)
	  vcu(1,ix-1,iy  ) = vcu(1,ix-1,iy  ) + wl1(1) * wp1(2)
	  vcu(1,ix  ,iy  ) = vcu(1,ix  ,iy  ) + wl1(2) * wp1(2)
	  vcu(1,ix-1,iy+1) = vcu(1,ix-1,iy+1) + wl1(1) * wp1(3)
	  vcu(1,ix  ,iy+1) = vcu(1,ix  ,iy+1) + wl1(2) * wp1(3)
	
	! accumulate j2
	  vcu(2,ix-1,iy-1) = vcu(2,ix-1,iy-1) + qvy * (S0x(1)*S0y(1)+
     1  S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	  vcu(2,ix  ,iy-1) = vcu(2,ix  ,iy-1) + qvy * (S0x(2)*S0y(1)+
     1  S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	  vcu(2,ix+1,iy-1) = vcu(2,ix+1,iy-1) + qvy * (S0x(3)*S0y(1)+
     1  S1x(3)*S1y(1)+(S0x(3)*S1y(1)+S1x(3)*S0y(1))/2.)
	  vcu(2,ix-1,iy  ) = vcu(2,ix-1,iy  ) + qvy * (S0x(1)*S0y(2)+
     1  S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	  vcu(2,ix  ,iy  ) = vcu(2,ix  ,iy  ) + qvy * (S0x(2)*S0y(2)+
     1  S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	  vcu(2,ix+1,iy  ) = vcu(2,ix+1,iy  ) + qvy * (S0x(3)*S0y(2)+
     1  S1x(3)*S1y(2)+(S0x(3)*S1y(2)+S1x(3)*S0y(2))/2.)
	  vcu(2,ix-1,iy+1) = vcu(2,ix-1,iy+1) + qvy * (S0x(1)*S0y(3)+
     1  S1x(1)*S1y(3)+(S0x(1)*S1y(3)+S1x(1)*S0y(3))/2.)
	  vcu(2,ix  ,iy+1) = vcu(2,ix  ,iy+1) + qvy * (S0x(2)*S0y(3)+
     1  S1x(2)*S1y(3)+(S0x(2)*S1y(3)+S1x(2)*S0y(3))/2.)
	  vcu(2,ix+1,iy+1) = vcu(2,ix+1,iy+1) + qvy * (S0x(3)*S0y(3)+
     1  S1x(3)*S1y(3)+(S0x(3)*S1y(3)+S1x(3)*S0y(3))/2.)
	
	! accumulate j3
	  vcu(3,ix-1,iy-1) = vcu(3,ix-1,iy-1) + qvz * (S0x(1)*S0y(1)+
     1  S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	  vcu(3,ix  ,iy-1) = vcu(3,ix  ,iy-1) + qvz * (S0x(2)*S0y(1)+
     1  S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	  vcu(3,ix+1,iy-1) = vcu(3,ix+1,iy-1) + qvz * (S0x(3)*S0y(1)+
     1  S1x(3)*S1y(1)+(S0x(3)*S1y(1)+S1x(3)*S0y(1))/2.)
	  vcu(3,ix-1,iy  ) = vcu(3,ix-1,iy  ) + qvz * (S0x(1)*S0y(2)+
     1  S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	  vcu(3,ix  ,iy  ) = vcu(3,ix  ,iy  ) + qvz * (S0x(2)*S0y(2)+
     1  S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	  vcu(3,ix+1,iy  ) = vcu(3,ix+1,iy  ) + qvz * (S0x(3)*S0y(2)+
     1  S1x(3)*S1y(2)+(S0x(3)*S1y(2)+S1x(3)*S0y(2))/2.)
	  vcu(3,ix-1,iy+1) = vcu(3,ix-1,iy+1) + qvz * (S0x(1)*S0y(3)+
     1  S1x(1)*S1y(3)+(S0x(1)*S1y(3)+S1x(1)*S0y(3))/2.)
	  vcu(3,ix  ,iy+1) = vcu(3,ix  ,iy+1) + qvz * (S0x(2)*S0y(3)+
     1  S1x(2)*S1y(3)+(S0x(2)*S1y(3)+S1x(2)*S0y(3))/2.)
	  vcu(3,ix+1,iy+1) = vcu(3,ix+1,iy+1) + qvz * (S0x(3)*S0y(3)+
     1  S1x(3)*S1y(3)+(S0x(3)*S1y(3)+S1x(3)*S0y(3))/2.)
    5 continue    
      else
        do 10 j = 1, nsplit
c find interpolation weights
      
        x0 = vpart(1,j)
        x1 = vpart(2,j)
        ix = vpart(3,j) + 1
        y0 = vpart(4,j)
        y1 = vpart(5,j)
        iy = vpart(6,j) + 1

        qnx = 1.0/(4.0*dt)*qm
        qny = 1.0/(4.0*dt)*qm
        qvz = vpart(idimp+4,j)*qm/3.0
      !print*, 'qvz =', qvz

	! get spline weitghts for x and y
	  S0x(1) = 0.5*(0.5 - x0)**2
	  S0x(2) = 0.75 - x0**2
	  S0x(3) = 0.5*(0.5 + x0)**2
	
	  S1x(1) = 0.5*(0.5 - x1)**2
	  S1x(2) = 0.75 - x1**2
	  S1x(3) = 0.5*(0.5 + x1)**2
	
	  S0y(1) = 0.5*(0.5 - y0)**2
	  S0y(2) = 0.75 - y0**2
	  S0y(3) = 0.5*(0.5 + y0)**2
	
	  S1y(1) =  0.5*(0.5 - y1)**2
	  S1y(2) = 0.75 - y1**2
	  S1y(3) = 0.5*(0.5 + y1)**2
	
	! get longitudinal motion weights
	  wl1(1) = (qnx*(x0 - x1)*(-1 + x0 + x1))
	  wl1(2) = (qnx*(-(x0*(1 + x0)) + x1*(1 + x1)))
	
	  wl2(1) = (qny*(y0 - y1)*(-1 + y0 + y1))
	  wl2(2) = (qny*(-(y0*(1 + y0)) + y1*(1 + y1)))
	
	! get perpendicular motion weights
	  wp1(1) = (S0y(1)+S1y(1))
	  wp1(2) = (S0y(2)+S1y(2))
	  wp1(3) = (S0y(3)+S1y(3))
	
	  wp2(1) = (S0x(1)+S1x(1))
	  wp2(2) = (S0x(2)+S1x(2))
	  wp2(3) = (S0x(3)+S1x(3))
	
	! accumulate j1
	  vcu(1,ix-1,iy-1) = vcu(1,ix-1,iy-1) + wl1(1) * wp1(1)
	  vcu(1,ix  ,iy-1) = vcu(1,ix  ,iy-1) + wl1(2) * wp1(1)
	  vcu(1,ix-1,iy  ) = vcu(1,ix-1,iy  ) + wl1(1) * wp1(2)
	  vcu(1,ix  ,iy  ) = vcu(1,ix  ,iy  ) + wl1(2) * wp1(2)
	  vcu(1,ix-1,iy+1) = vcu(1,ix-1,iy+1) + wl1(1) * wp1(3)
	  vcu(1,ix  ,iy+1) = vcu(1,ix  ,iy+1) + wl1(2) * wp1(3)
	
	! accumulate j2
	  vcu(2,ix-1,iy-1) = vcu(2,ix-1,iy-1) + wl2(1) * wp2(1)
	  vcu(2,ix  ,iy-1) = vcu(2,ix  ,iy-1) + wl2(1) * wp2(2)
	  vcu(2,ix+1,iy-1) = vcu(2,ix+1,iy-1) + wl2(1) * wp2(3)
	  vcu(2,ix-1,iy  ) = vcu(2,ix-1,iy  ) + wl2(2) * wp2(1)
	  vcu(2,ix  ,iy  ) = vcu(2,ix  ,iy  ) + wl2(2) * wp2(2)
	  vcu(2,ix+1,iy  ) = vcu(2,ix+1,iy  ) + wl2(2) * wp2(3)
	
	! accumulate j3
	  vcu(3,ix-1,iy-1) = vcu(3,ix-1,iy-1) + qvz * (S0x(1)*S0y(1)+
     1  S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	  vcu(3,ix  ,iy-1) = vcu(3,ix  ,iy-1) + qvz * (S0x(2)*S0y(1)+
     1  S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	  vcu(3,ix+1,iy-1) = vcu(3,ix+1,iy-1) + qvz * (S0x(3)*S0y(1)+
     1  S1x(3)*S1y(1)+(S0x(3)*S1y(1)+S1x(3)*S0y(1))/2.)
	  vcu(3,ix-1,iy  ) = vcu(3,ix-1,iy  ) + qvz * (S0x(1)*S0y(2)+
     1  S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	  vcu(3,ix  ,iy  ) = vcu(3,ix  ,iy  ) + qvz * (S0x(2)*S0y(2)+
     1  S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	  vcu(3,ix+1,iy  ) = vcu(3,ix+1,iy  ) + qvz * (S0x(3)*S0y(2)+
     1  S1x(3)*S1y(2)+(S0x(3)*S1y(2)+S1x(3)*S0y(2))/2.)
	  vcu(3,ix-1,iy+1) = vcu(3,ix-1,iy+1) + qvz * (S0x(1)*S0y(3)+
     1  S1x(1)*S1y(3)+(S0x(1)*S1y(3)+S1x(1)*S0y(3))/2.)
	  vcu(3,ix  ,iy+1) = vcu(3,ix  ,iy+1) + qvz * (S0x(2)*S0y(3)+
     1  S1x(2)*S1y(3)+(S0x(2)*S1y(3)+S1x(2)*S0y(3))/2.)
	  vcu(3,ix+1,iy+1) = vcu(3,ix+1,iy+1) + qvz * (S0x(3)*S0y(3)+
     1  S1x(3)*S1y(3)+(S0x(3)*S1y(3)+S1x(3)*S0y(3))/2.)
   10 continue
      endif
      do 20 j = 1, nop
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
c advance position half a time-step
      dx = part(1,j) + vx*dth
      dy = part(2,j) + vy*dth
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
   20 continue
      return
      end !!!
      subroutine GSRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,ipb
     1c)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 88 flops/particle, 1 divide, 1 sqrt, 32 loads, 29 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
      ci2 = ci*ci
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      vzn = part(5,1)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      do 10 j = 2, nop
c find interpolation weights
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      vzn = part(5,j)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
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
      subroutine GSRJPOST2X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,ip
     1bc)
c for 2-1/2d code, this subroutine calculates particle current density
c using second-order spline interpolation, for relativistic particles,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 88 flops/particle, 1 divide, 1 sqrt, 87 loads, 83 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
      do 10 j = 1, npb
c find interpolation weights
      n = part(1,j+jb) + .5
      m = part(2,j+jb) + .5
      dxp = part(1,j+jb) - float(n)
      dyp = part(2,j+jb) - float(m)
c find inverse gamma
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      vz = part(5,j+jb)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc
     1)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+1
c nyv = second dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(3,nxv,nyv)
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
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine OGRJPOST2L(part,vpart,vcu,qm,dth,ci,nsplit,nop, !!!
     1idimp,nx,ny,nxv,nyv,ipbc,dt)
!*****
!This subroutine deposits the charge conserving current for linear
!particle shape into the virtual current array. It is based on the
!OSIRIS subroutine getjr_2d_s1. It also does a half position update, 
!like GRJPOST2L.
!*****
      dimension part(idimp,nop), vcu(3,nx+3,ny+3) 
      dimension vpart(idimp+4,4*nop)
      dimension S0x(2), S1x(2), S0y(2), S1y(2)
      dimension wp1(2), wp2(2)
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

      if(ny==1) then
        do 5 j = 1, nsplit
c find interpolation weights
     
        x0 = vpart(1,j)
        x1 = vpart(2,j)
        ix = vpart(3,j)
        y0 = vpart(4,j)
        y1 = vpart(5,j)
        iy = vpart(6,j)

        qnx = 1.0/(2.0*dt)*qm
        qvy = vpart(idimp+3,j)*qm/3.0
        qvz = vpart(idimp+4,j)*qm/3.0

        S0x(1) = 0.5 - x0
        S0x(2) = 0.5 + x0 
        S1x(1) = 0.5 - x1
        S1x(2) = 0.5 + x1    
 
        S0y(1) = 0.5 - y0
        S0y(2) = 0.5 + y0 
        S1y(1) = 0.5 - y1
        S1y(2) = 0.5 + y1

        wl1 = qnx*(-x0 + x1) 
        wl2 = qny*(-y0 + y1)

        wp1(1) = S0y(1) + S1y(1)
        wp1(2) = S0y(2) + S1y(2)

! x current
        vcu(1,ix,iy) = vcu(1,ix,iy) + wl1*wp1(1)
        vcu(1,ix,iy+1) = vcu(1,ix,iy+1) + wl1*wp1(2)

! y current
        vcu(2,ix  ,iy  ) = vcu(2,ix  ,iy  ) + qvy * (S0x(1)*S0y(1)
     1  +S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)

        vcu(2,ix+1,iy  ) = vcu(2,ix+1,iy  ) + qvy * (S0x(2)*S0y(1)
     1  +S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)

        vcu(2,ix  ,iy+1) = vcu(2,ix  ,iy+1) + qvy * (S0x(1)*S0y(2)
     1  +S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)

        vcu(2,ix+1,iy+1) = vcu(2,ix+1,iy+1) + qvy * (S0x(2)*S0y(2)
     1  +S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)

! z current
        vcu(3,ix  ,iy  ) = vcu(3,ix  ,iy  ) + qvz * (S0x(1)*S0y(1)
     1  +S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)

        vcu(3,ix+1,iy  ) = vcu(3,ix+1,iy  ) + qvz * (S0x(2)*S0y(1)
     1  +S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)

        vcu(3,ix  ,iy+1) = vcu(3,ix  ,iy+1) + qvz * (S0x(1)*S0y(2)
     1  +S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)

        vcu(3,ix+1,iy+1) = vcu(3,ix+1,iy+1) + qvz * (S0x(2)*S0y(2)
     1  +S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
    5 continue
      else
        do 10 j = 1, nsplit
c find interpolation weights
     
        x0 = vpart(1,j)
        x1 = vpart(2,j)
        ix = vpart(3,j)
        y0 = vpart(4,j)
        y1 = vpart(5,j)
        iy = vpart(6,j)

        qnx = 1.0/(2.0*dt)*qm
        qny = 1.0/(2.0*dt)*qm
        qvz = vpart(idimp+4,j)*qm/3.0

        S0x(1) = 0.5 - x0
        S0x(2) = 0.5 + x0 
        S1x(1) = 0.5 - x1
        S1x(2) = 0.5 + x1    
 
        S0y(1) = 0.5 - y0
        S0y(2) = 0.5 + y0 
        S1y(1) = 0.5 - y1
        S1y(2) = 0.5 + y1

        wl1 = qnx*(-x0 + x1) 
        wl2 = qny*(-y0 + y1)

        wp1(1) = S0y(1) + S1y(1)
        wp1(2) = S0y(2) + S1y(2)
        wp2(1) = S0x(1) + S1x(1)
        wp2(2) = S0x(2) + S1x(2)

! x current
        vcu(1,ix,iy) = vcu(1,ix,iy) + wl1*wp1(1)
        vcu(1,ix,iy+1) = vcu(1,ix,iy+1) + wl1*wp1(2)

! y current
        vcu(2,ix,iy) = vcu(2,ix,iy) + wl2*wp2(1)
        vcu(2,ix+1,iy) = vcu(2,ix+1,iy) + wl2*wp2(2)

! z current
        vcu(3,ix  ,iy  ) = vcu(3,ix  ,iy  ) + qvz * (S0x(1)*S0y(1)
     1  +S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)

        vcu(3,ix+1,iy  ) = vcu(3,ix+1,iy  ) + qvz * (S0x(2)*S0y(1)
     1  +S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)

        vcu(3,ix  ,iy+1) = vcu(3,ix  ,iy+1) + qvz * (S0x(1)*S0y(2)
     1  +S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)

        vcu(3,ix+1,iy+1) = vcu(3,ix+1,iy+1) + qvz * (S0x(2)*S0y(2)
     1  +S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
   10 continue
      endif     
      do 20 j = 1, nop
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
c advance position half a time-step
      dx = part(1,j) + vx*dth
      dy = part(2,j) + vy*dth
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
   20 continue
      return
      end !!!
      subroutine GSRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,ip
     1bc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
      ci2 = ci*ci
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      vzn = part(5,1)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      do 10 j = 2, nop
c find interpolation weights
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      vzn = part(5,j)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
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
      subroutine GSRJPOST2XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,i
     1pbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using first-order linear interpolation, for relativistic particles,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 47 flops/particle, 1 divide, 1 sqrt,, 41 loads, 33 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
      do 10 j = 1, npb
c find interpolation weights
      n = part(1,j+jb)
      m = part(2,j+jb)
      dxp = qm*(part(1,j+jb) - float(n))
      dyp = part(2,j+jb) - float(m)
c find inverse gamma
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      vz = part(5,j+jb)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine GRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc
     1)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 67 flops/particle, 1 divide, 1 sqrt, 22 loads, 20 stores
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
c and qci = qm*pi*gami, where i = x,y
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+3
c nyv = second dimension of current array, must be >= ny+3
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(2,nxv,nyv)
      qmh = .5*qm
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
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j) + .5
      mm = part(2,j) + .5
      dxp = part(1,j) - float(nn)
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c deposit current
      dx = dxl*amy
      dy = amx*amy
      dz = dxp*amy
      vx = vx*gami
      vy = vy*gami
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
      subroutine GSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,ip
     1bc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 67 flops/particle, 1 divide, 1 sqrt, 22 loads, 20 stores
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
c and qci = qm*pi*gami, where i = x,y
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
      ci2 = ci*ci
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      p2 = vxn*vxn + vyn*vyn
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      do 10 j = 2, nop
c find interpolation weights
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
      vx = vxn*gami
      vy = vyn*gami
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      p2 = vxn*vxn + vyn*vyn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      vx = vxn*gami
      vy = vyn*gami
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
      subroutine GSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,i
     1pbc)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation, for relativistic particles,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 67 flops/particle, 1 divide, 1 sqrt, 58 loads, 56 stores
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
c and qci = qm*pi*gami, where i = x,y
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
      do 10 j = 1, npb
c find interpolation weights
      n = part(1,j+jb) + .5
      m = part(2,j+jb) + .5
      dxp = part(1,j+jb) - float(n)
      dyp = part(2,j+jb) - float(m)
c find inverse gamma
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
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
      subroutine GRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipb
     1c)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 36 flops/particle, 1 divide, 1 sqrt, 12 loads, 10 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,j,k) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of current array, must be >= nx+1
c nyv = second dimension of current array, must be >= ny+1
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop), cu(2,nxv,nyv)
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
      do 10 j = 1, nop
c find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - float(nn))
      dyp = part(2,j) - float(mm)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
c deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = vx*gami
      vy = vy*gami
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
      subroutine GSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,i
     1pbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 36 flops/particle, 1 divide, 1 sqrt, 12 loads, 10 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
      ci2 = ci*ci
c find inverse gamma
      vxn = part(3,1)
      vyn = part(4,1)
      p2 = vxn*vxn + vyn*vyn
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      do 10 j = 2, nop
c find interpolation weights
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
      vx = vxn*gami
      vy = vyn*gami
c get momentum for next particle
      vxn = part(3,j)
      vyn = part(4,j)
      p2 = vxn*vxn + vyn*vyn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      vx = vxn*gami
      vy = vyn*gami
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
      subroutine GSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,
     1ipbc)
c for 2d code, this subroutine calculates particle current density
c using first-order linear interpolation, for relativistic particles,
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 36 flops/particle, 1 divide, 1 sqrt, 29 loads, 26 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c where n,m = nearest grid points and dx = x-n, dy = y-m
c and qci = qm*pi*gami, where i = x,y
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,n) = ith component of current density at grid point j,k
c where n = j + nxv*(k-1)
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
c outer loop over blocks of particles
      do 40 k = 1, ipp
      je = k*npb
      jb = je - npb
      if (je.gt.nop) npb = nop - jb
      do 10 j = 1, npb
c find interpolation weights
      n = part(1,j+jb)
      m = part(2,j+jb)
      dxp = qm*(part(1,j+jb) - float(n))
      dyp = part(2,j+jb) - float(m)
c find inverse gamma
      vx = part(3,j+jb)
      vy = part(4,j+jb)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
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
      subroutine GRJPOST2C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc
     1)
c for 2-1/2d code, this subroutine calculates particle current density
c using third-order spline interpolation or relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 162 flops/particle, 1 divide, 1 sqrt, 53 loads, 50 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxv = first dimension of charge array, must be >= nx+5
c nxy = second dimension of charge array, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real part, cu, qm, dt, ci
      dimension part(idimp,nop), cu(3,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2, ci2, p2, gami
      real edgelx, edgely, edgerx, edgery, dx, dy, dz, dw, vx, vy, vz
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
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
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c deposit charge
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine OGRJPOST2C(part,vpart,vcu,qm,dth,ci,nsplit,nop, !!!
     1idimp,nx,ny,nxv,nyv,ipbc,dt)
!*****
!This subroutine deposits the charge conserving current for cubic
!particle shape into the virtual current array. It is based on the
!OSIRIS subroutine getjr_2d_s3. It also does a half position update, 
!like GRJPOST2C.
!*****
      dimension part(idimp,nop), vcu(3,nx+5,ny+5) 
      dimension vpart(idimp+4,4*nop)
      dimension S0x(4), S1x(4), S0y(4), S1y(4)
      dimension wp1(4), wp2(4), wl1(3), wl2(3)
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

      if(ny==1) then
        do 5 j = 1, nsplit
c find interpolation weights
      
        x0 = vpart(1,j)
        x1 = vpart(2,j)
        ix = vpart(3,j) + 1
        y0 = vpart(4,j)
        y1 = vpart(5,j)
        iy = vpart(6,j) + 1

        qnx = 1.0/dt*qm
        qvy = vpart(idimp+3,j)*qm/3.0
        qvz = vpart(idimp+4,j)*qm/3.0

	! get spline weitghts for x and y
	  S0x(1) = -(-0.5 + x0)**3/6.
	  S0x(2) = (4 - 6*(0.5 + x0)**2 + 3*(0.5 + x0)**3)/6.
	  S0x(3) = (23 + 30*x0 - 12*x0**2 - 24*x0**3)/48.
	  S0x(4) = (0.5 + x0)**3/6.
	
	  S1x(1) = -(-0.5 + x1)**3/6.
	  S1x(2) = (4 - 6*(0.5 + x1)**2 + 3*(0.5 + x1)**3)/6.
	  S1x(3) = (23 + 30*x1 - 12*x1**2 - 24*x1**3)/48.
	  S1x(4) = (0.5 + x1)**3/6.
	
	  S0y(1) = -(-0.5 + y0)**3/6.
	  S0y(2) = (4 - 6*(0.5 + y0)**2 + 3*(0.5 + y0)**3)/6.
	  S0y(3) = (23 + 30*y0 - 12*y0**2 - 24*y0**3)/48.
	  S0y(4) = (0.5 + y0)**3/6.
	
	  S1y(1) = -(-0.5 + y1)**3/6.
	  S1y(2) = (4 - 6*(0.5 + y1)**2 + 3*(0.5 + y1)**3)/6.
	  S1y(3) = (23 + 30*y1 - 12*y1**2 - 24*y1**3)/48.
	  S1y(4) = (0.5 + y1)**3/6.
	
	
	! get longitudinal motion weights
	  wl1(1) = (qnx*(-(-0.5 + x0)**3 + (-0.5 + x1)**3))/6.
	  wl1(2) = (qnx*(-9*x0 + 4*x0**3 + 9*x1 - 4*x1**3))/12.
	  wl1(3) = (qnx*(-(x0*(3 + 6*x0 + 4*x0**2)) + x1*(3 + 6*x1
     1   + 4*x1**2)))/24.

	
	! get perpendicular motion weights
	  wp1(1) = 0.5*(S0y(1)+S1y(1))
	  wp1(2) = 0.5*(S0y(2)+S1y(2))
	  wp1(3) = 0.5*(S0y(3)+S1y(3))
	  wp1(4) = 0.5*(S0y(4)+S1y(4))

	
	! accumulate j1
	  vcu(1,ix-1,iy-1) = vcu(1,ix-1,iy-1) + wl1(1) * wp1(1)
	  vcu(1,ix  ,iy-1) = vcu(1,ix  ,iy-1) + wl1(2) * wp1(1)
	  vcu(1,ix+1,iy-1) = vcu(1,ix+1,iy-1) + wl1(3) * wp1(1)
	  vcu(1,ix-1,iy  ) = vcu(1,ix-1,iy  ) + wl1(1) * wp1(2)
	  vcu(1,ix  ,iy  ) = vcu(1,ix  ,iy  ) + wl1(2) * wp1(2)
	  vcu(1,ix+1,iy  ) = vcu(1,ix+1,iy  ) + wl1(3) * wp1(2)
	  vcu(1,ix-1,iy+1) = vcu(1,ix-1,iy+1) + wl1(1) * wp1(3)
	  vcu(1,ix  ,iy+1) = vcu(1,ix  ,iy+1) + wl1(2) * wp1(3)
	  vcu(1,ix+1,iy+1) = vcu(1,ix+1,iy+1) + wl1(3) * wp1(3)
	  vcu(1,ix-1,iy+2) = vcu(1,ix-1,iy+2) + wl1(1) * wp1(4)
	  vcu(1,ix  ,iy+2) = vcu(1,ix  ,iy+2) + wl1(2) * wp1(4)
	  vcu(1,ix+1,iy+2) = vcu(1,ix+1,iy+2) + wl1(3) * wp1(4)
	
	! accumulate j2
	  vcu(2,ix-1,iy-1) = vcu(2,ix-1,iy-1) + qvy * (S0x(1)*S0y(1)
     1  +S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	  vcu(2,ix  ,iy-1) = vcu(2,ix  ,iy-1) + qvy * (S0x(2)*S0y(1)
     1  +S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	  vcu(2,ix+1,iy-1) = vcu(2,ix+1,iy-1) + qvy * (S0x(3)*S0y(1)
     1  +S1x(3)*S1y(1)+(S0x(3)*S1y(1)+S1x(3)*S0y(1))/2.)
	  vcu(2,ix+2,iy-1) = vcu(2,ix+2,iy-1) + qvy * (S0x(4)*S0y(1)
     1  +S1x(4)*S1y(1)+(S0x(4)*S1y(1)+S1x(4)*S0y(1))/2.)
	  vcu(2,ix-1,iy  ) = vcu(2,ix-1,iy  ) + qvy * (S0x(1)*S0y(2)
     1  +S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	  vcu(2,ix  ,iy  ) = vcu(2,ix  ,iy  ) + qvy * (S0x(2)*S0y(2)
     1  +S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	  vcu(2,ix+1,iy  ) = vcu(2,ix+1,iy  ) + qvy * (S0x(3)*S0y(2)
     1  +S1x(3)*S1y(2)+(S0x(3)*S1y(2)+S1x(3)*S0y(2))/2.)
	  vcu(2,ix+2,iy  ) = vcu(2,ix+2,iy  ) + qvy * (S0x(4)*S0y(2)
     1  +S1x(4)*S1y(2)+(S0x(4)*S1y(2)+S1x(4)*S0y(2))/2.)
	  vcu(2,ix-1,iy+1) = vcu(2,ix-1,iy+1) + qvy * (S0x(1)*S0y(3)
     1  +S1x(1)*S1y(3)+(S0x(1)*S1y(3)+S1x(1)*S0y(3))/2.)
	  vcu(2,ix  ,iy+1) = vcu(2,ix  ,iy+1) + qvy * (S0x(2)*S0y(3)
     1  +S1x(2)*S1y(3)+(S0x(2)*S1y(3)+S1x(2)*S0y(3))/2.)
	  vcu(2,ix+1,iy+1) = vcu(2,ix+1,iy+1) + qvy * (S0x(3)*S0y(3)
     1  +S1x(3)*S1y(3)+(S0x(3)*S1y(3)+S1x(3)*S0y(3))/2.)
	  vcu(2,ix+2,iy+1) = vcu(2,ix+2,iy+1) + qvy * (S0x(4)*S0y(3)
     1  +S1x(4)*S1y(3)+(S0x(4)*S1y(3)+S1x(4)*S0y(3))/2.)
	  vcu(2,ix-1,iy+2) = vcu(2,ix-1,iy+2) + qvy * (S0x(1)*S0y(4)
     1  +S1x(1)*S1y(4)+(S0x(1)*S1y(4)+S1x(1)*S0y(4))/2.)
	  vcu(2,ix  ,iy+2) = vcu(2,ix  ,iy+2) + qvy * (S0x(2)*S0y(4)
     1  +S1x(2)*S1y(4)+(S0x(2)*S1y(4)+S1x(2)*S0y(4))/2.)
	  vcu(2,ix+1,iy+2) = vcu(2,ix+1,iy+2) + qvy * (S0x(3)*S0y(4)
     1  +S1x(3)*S1y(4)+(S0x(3)*S1y(4)+S1x(3)*S0y(4))/2.)
	  vcu(2,ix+2,iy+2) = vcu(2,ix+2,iy+2) + qvy * (S0x(4)*S0y(4)
     1  +S1x(4)*S1y(4)+(S0x(4)*S1y(4)+S1x(4)*S0y(4))/2.)
	
	! accumulate j3
	  vcu(3,ix-1,iy-1) = vcu(3,ix-1,iy-1) + qvz * (S0x(1)*S0y(1)
     1  +S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	  vcu(3,ix  ,iy-1) = vcu(3,ix  ,iy-1) + qvz * (S0x(2)*S0y(1)
     1  +S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	  vcu(3,ix+1,iy-1) = vcu(3,ix+1,iy-1) + qvz * (S0x(3)*S0y(1)
     1  +S1x(3)*S1y(1)+(S0x(3)*S1y(1)+S1x(3)*S0y(1))/2.)
	  vcu(3,ix+2,iy-1) = vcu(3,ix+2,iy-1) + qvz * (S0x(4)*S0y(1)
     1  +S1x(4)*S1y(1)+(S0x(4)*S1y(1)+S1x(4)*S0y(1))/2.)
	  vcu(3,ix-1,iy  ) = vcu(3,ix-1,iy  ) + qvz * (S0x(1)*S0y(2)
     1  +S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	  vcu(3,ix  ,iy  ) = vcu(3,ix  ,iy  ) + qvz * (S0x(2)*S0y(2)
     1  +S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	  vcu(3,ix+1,iy  ) = vcu(3,ix+1,iy  ) + qvz * (S0x(3)*S0y(2)
     1  +S1x(3)*S1y(2)+(S0x(3)*S1y(2)+S1x(3)*S0y(2))/2.)
	  vcu(3,ix+2,iy  ) = vcu(3,ix+2,iy  ) + qvz * (S0x(4)*S0y(2)
     1  +S1x(4)*S1y(2)+(S0x(4)*S1y(2)+S1x(4)*S0y(2))/2.)
	  vcu(3,ix-1,iy+1) = vcu(3,ix-1,iy+1) + qvz * (S0x(1)*S0y(3)
     1  +S1x(1)*S1y(3)+(S0x(1)*S1y(3)+S1x(1)*S0y(3))/2.)
	  vcu(3,ix  ,iy+1) = vcu(3,ix  ,iy+1) + qvz * (S0x(2)*S0y(3)
     1  +S1x(2)*S1y(3)+(S0x(2)*S1y(3)+S1x(2)*S0y(3))/2.)
	  vcu(3,ix+1,iy+1) = vcu(3,ix+1,iy+1) + qvz * (S0x(3)*S0y(3)
     1  +S1x(3)*S1y(3)+(S0x(3)*S1y(3)+S1x(3)*S0y(3))/2.)
	  vcu(3,ix+2,iy+1) = vcu(3,ix+2,iy+1) + qvz * (S0x(4)*S0y(3)
     1  +S1x(4)*S1y(3)+(S0x(4)*S1y(3)+S1x(4)*S0y(3))/2.)
	  vcu(3,ix-1,iy+2) = vcu(3,ix-1,iy+2) + qvz * (S0x(1)*S0y(4)
     1  +S1x(1)*S1y(4)+(S0x(1)*S1y(4)+S1x(1)*S0y(4))/2.)
	  vcu(3,ix  ,iy+2) = vcu(3,ix  ,iy+2) + qvz * (S0x(2)*S0y(4)
     1  +S1x(2)*S1y(4)+(S0x(2)*S1y(4)+S1x(2)*S0y(4))/2.)
	  vcu(3,ix+1,iy+2) = vcu(3,ix+1,iy+2) + qvz * (S0x(3)*S0y(4)
     1  +S1x(3)*S1y(4)+(S0x(3)*S1y(4)+S1x(3)*S0y(4))/2.)
	  vcu(3,ix+2,iy+2) = vcu(3,ix+2,iy+2) + qvz * (S0x(4)*S0y(4)
     1  +S1x(4)*S1y(4)+(S0x(4)*S1y(4)+S1x(4)*S0y(4))/2.)
    5 continue 	
      else
        do 10 j = 1, nsplit
c find interpolation weights
      
        x0 = vpart(1,j)
        x1 = vpart(2,j)
        ix = vpart(3,j) + 1
        y0 = vpart(4,j)
        y1 = vpart(5,j)
        iy = vpart(6,j) + 1

        qnx = 1.0/dt*qm
        qny = 1.0/dt*qm
        qvz = vpart(idimp+4,j)*qm/3.0

	! get spline weitghts for x and y
	  S0x(1) = -(-0.5 + x0)**3/6.
	  S0x(2) = (4 - 6*(0.5 + x0)**2 + 3*(0.5 + x0)**3)/6.
	  S0x(3) = (23 + 30*x0 - 12*x0**2 - 24*x0**3)/48.
	  S0x(4) = (0.5 + x0)**3/6.
	
	  S1x(1) = -(-0.5 + x1)**3/6.
	  S1x(2) = (4 - 6*(0.5 + x1)**2 + 3*(0.5 + x1)**3)/6.
	  S1x(3) = (23 + 30*x1 - 12*x1**2 - 24*x1**3)/48.
	  S1x(4) = (0.5 + x1)**3/6.
	
	  S0y(1) = -(-0.5 + y0)**3/6.
	  S0y(2) = (4 - 6*(0.5 + y0)**2 + 3*(0.5 + y0)**3)/6.
	  S0y(3) = (23 + 30*y0 - 12*y0**2 - 24*y0**3)/48.
	  S0y(4) = (0.5 + y0)**3/6.
	
	  S1y(1) = -(-0.5 + y1)**3/6.
	  S1y(2) = (4 - 6*(0.5 + y1)**2 + 3*(0.5 + y1)**3)/6.
	  S1y(3) = (23 + 30*y1 - 12*y1**2 - 24*y1**3)/48.
	  S1y(4) = (0.5 + y1)**3/6.
	
	
	! get longitudinal motion weights
	  wl1(1) = (qnx*(-(-0.5 + x0)**3 + (-0.5 + x1)**3))/6.
	  wl1(2) = (qnx*(-9*x0 + 4*x0**3 + 9*x1 - 4*x1**3))/12.
	  wl1(3) = (qnx*(-(x0*(3 + 6*x0 + 4*x0**2)) + x1*(3 + 6*x1
     1   + 4*x1**2)))/24.
	
	  wl2(1) = (qny*(-(-0.5 + y0)**3 + (-0.5 + y1)**3))/6.
	  wl2(2) = (qny*(-9*y0 + 4*y0**3 + 9*y1 - 4*y1**3))/12.
	  wl2(3) = (qny*(-(y0*(3 + 6*y0 + 4*y0**2)) + y1*(3 + 6*y1
     1   + 4*y1**2)))/24.
	
	! get perpendicular motion weights
	  wp1(1) = 0.5*(S0y(1)+S1y(1))
	  wp1(2) = 0.5*(S0y(2)+S1y(2))
	  wp1(3) = 0.5*(S0y(3)+S1y(3))
	  wp1(4) = 0.5*(S0y(4)+S1y(4))
	
	  wp2(1) = 0.5*(S0x(1)+S1x(1))
	  wp2(2) = 0.5*(S0x(2)+S1x(2))
	  wp2(3) = 0.5*(S0x(3)+S1x(3))
	  wp2(4) = 0.5*(S0x(4)+S1x(4))
	
	! accumulate j1
	  vcu(1,ix-1,iy-1) = vcu(1,ix-1,iy-1) + wl1(1) * wp1(1)
	  vcu(1,ix  ,iy-1) = vcu(1,ix  ,iy-1) + wl1(2) * wp1(1)
	  vcu(1,ix+1,iy-1) = vcu(1,ix+1,iy-1) + wl1(3) * wp1(1)
	  vcu(1,ix-1,iy  ) = vcu(1,ix-1,iy  ) + wl1(1) * wp1(2)
	  vcu(1,ix  ,iy  ) = vcu(1,ix  ,iy  ) + wl1(2) * wp1(2)
	  vcu(1,ix+1,iy  ) = vcu(1,ix+1,iy  ) + wl1(3) * wp1(2)
	  vcu(1,ix-1,iy+1) = vcu(1,ix-1,iy+1) + wl1(1) * wp1(3)
	  vcu(1,ix  ,iy+1) = vcu(1,ix  ,iy+1) + wl1(2) * wp1(3)
	  vcu(1,ix+1,iy+1) = vcu(1,ix+1,iy+1) + wl1(3) * wp1(3)
	  vcu(1,ix-1,iy+2) = vcu(1,ix-1,iy+2) + wl1(1) * wp1(4)
	  vcu(1,ix  ,iy+2) = vcu(1,ix  ,iy+2) + wl1(2) * wp1(4)
	  vcu(1,ix+1,iy+2) = vcu(1,ix+1,iy+2) + wl1(3) * wp1(4)
	
	! accumulate j2
	  vcu(2,ix-1,iy-1) = vcu(2,ix-1,iy-1) + wl2(1) * wp2(1)
	  vcu(2,ix  ,iy-1) = vcu(2,ix  ,iy-1) + wl2(1) * wp2(2)
	  vcu(2,ix+1,iy-1) = vcu(2,ix+1,iy-1) + wl2(1) * wp2(3)
	  vcu(2,ix+2,iy-1) = vcu(2,ix+2,iy-1) + wl2(1) * wp2(4)
	  vcu(2,ix-1,iy  ) = vcu(2,ix-1,iy  ) + wl2(2) * wp2(1)
	  vcu(2,ix  ,iy  ) = vcu(2,ix  ,iy  ) + wl2(2) * wp2(2)
	  vcu(2,ix+1,iy  ) = vcu(2,ix+1,iy  ) + wl2(2) * wp2(3)
	  vcu(2,ix+2,iy  ) = vcu(2,ix+2,iy  ) + wl2(2) * wp2(4)
	  vcu(2,ix-1,iy+1) = vcu(2,ix-1,iy+1) + wl2(3) * wp2(1)
	  vcu(2,ix  ,iy+1) = vcu(2,ix  ,iy+1) + wl2(3) * wp2(2)
	  vcu(2,ix+1,iy+1) = vcu(2,ix+1,iy+1) + wl2(3) * wp2(3)
	  vcu(2,ix+2,iy+1) = vcu(2,ix+2,iy+1) + wl2(3) * wp2(4)
	
	! accumulate j3
	  vcu(3,ix-1,iy-1) = vcu(3,ix-1,iy-1) + qvz * (S0x(1)*S0y(1)
     1  +S1x(1)*S1y(1)+(S0x(1)*S1y(1)+S1x(1)*S0y(1))/2.)
	  vcu(3,ix  ,iy-1) = vcu(3,ix  ,iy-1) + qvz * (S0x(2)*S0y(1)
     1  +S1x(2)*S1y(1)+(S0x(2)*S1y(1)+S1x(2)*S0y(1))/2.)
	  vcu(3,ix+1,iy-1) = vcu(3,ix+1,iy-1) + qvz * (S0x(3)*S0y(1)
     1  +S1x(3)*S1y(1)+(S0x(3)*S1y(1)+S1x(3)*S0y(1))/2.)
	  vcu(3,ix+2,iy-1) = vcu(3,ix+2,iy-1) + qvz * (S0x(4)*S0y(1)
     1  +S1x(4)*S1y(1)+(S0x(4)*S1y(1)+S1x(4)*S0y(1))/2.)
	  vcu(3,ix-1,iy  ) = vcu(3,ix-1,iy  ) + qvz * (S0x(1)*S0y(2)
     1  +S1x(1)*S1y(2)+(S0x(1)*S1y(2)+S1x(1)*S0y(2))/2.)
	  vcu(3,ix  ,iy  ) = vcu(3,ix  ,iy  ) + qvz * (S0x(2)*S0y(2)
     1  +S1x(2)*S1y(2)+(S0x(2)*S1y(2)+S1x(2)*S0y(2))/2.)
	  vcu(3,ix+1,iy  ) = vcu(3,ix+1,iy  ) + qvz * (S0x(3)*S0y(2)
     1  +S1x(3)*S1y(2)+(S0x(3)*S1y(2)+S1x(3)*S0y(2))/2.)
	  vcu(3,ix+2,iy  ) = vcu(3,ix+2,iy  ) + qvz * (S0x(4)*S0y(2)
     1  +S1x(4)*S1y(2)+(S0x(4)*S1y(2)+S1x(4)*S0y(2))/2.)
	  vcu(3,ix-1,iy+1) = vcu(3,ix-1,iy+1) + qvz * (S0x(1)*S0y(3)
     1  +S1x(1)*S1y(3)+(S0x(1)*S1y(3)+S1x(1)*S0y(3))/2.)
	  vcu(3,ix  ,iy+1) = vcu(3,ix  ,iy+1) + qvz * (S0x(2)*S0y(3)
     1  +S1x(2)*S1y(3)+(S0x(2)*S1y(3)+S1x(2)*S0y(3))/2.)
	  vcu(3,ix+1,iy+1) = vcu(3,ix+1,iy+1) + qvz * (S0x(3)*S0y(3)
     1  +S1x(3)*S1y(3)+(S0x(3)*S1y(3)+S1x(3)*S0y(3))/2.)
	  vcu(3,ix+2,iy+1) = vcu(3,ix+2,iy+1) + qvz * (S0x(4)*S0y(3)
     1  +S1x(4)*S1y(3)+(S0x(4)*S1y(3)+S1x(4)*S0y(3))/2.)
	  vcu(3,ix-1,iy+2) = vcu(3,ix-1,iy+2) + qvz * (S0x(1)*S0y(4)
     1  +S1x(1)*S1y(4)+(S0x(1)*S1y(4)+S1x(1)*S0y(4))/2.)
	  vcu(3,ix  ,iy+2) = vcu(3,ix  ,iy+2) + qvz * (S0x(2)*S0y(4)
     1  +S1x(2)*S1y(4)+(S0x(2)*S1y(4)+S1x(2)*S0y(4))/2.)
	  vcu(3,ix+1,iy+2) = vcu(3,ix+1,iy+2) + qvz * (S0x(3)*S0y(4)
     1  +S1x(3)*S1y(4)+(S0x(3)*S1y(4)+S1x(3)*S0y(4))/2.)
	  vcu(3,ix+2,iy+2) = vcu(3,ix+2,iy+2) + qvz * (S0x(4)*S0y(4)
     1  +S1x(4)*S1y(4)+(S0x(4)*S1y(4)+S1x(4)*S0y(4))/2.)
	
   10 continue
      endif
     
      do 20 j = 1, nop
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
c advance position half a time-step
      dx = part(1,j) + vx*dth
      dy = part(2,j) + vy*dth
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
   20 continue
      return
      end !!!
      subroutine GRJPOST22C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipb
     1c)
c for 2d code, this subroutine calculates particle current density
c using third-order spline interpolation for relativistic particles
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 127 flops/particle, 1 divide, 1 sqrt, 36 loads, 34 stores
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
c and qci = qm*pi*gami, where i = x,y
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c cu(i,j+1,k+1) = ith component of current density at grid point j,k
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
c nx/ny = system length in x/y direction
c nxv = first dimension of charge array, must be >= nx+5
c nxy = second dimension of charge array, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real part, cu, qm, dt, ci
      dimension part(idimp,nop), cu(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qms, qmh, qm23, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq
      real at1, at2, ci2, p2, gami
      real edgelx, edgely, edgerx, edgery, dx, dy, dz, dw, vx, vy
      qms = sixth*qm
      qmh = 0.5*qm
      qm23 = 4.0*qms
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
      do 10 j = 1, nop
c find interpolation weights
      nl = part(1,j)
      ml = part(2,j)
      dxp = part(1,j) - float(nl)
      dyp = part(2,j) - float(ml)
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c deposit charge
      dx = dxl*dyl
      dy = dxn*dyl
      dz = dxp*dyl
      dw = dxq*dyl
      vx = vx*gami
      vy = vy*gami
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
      subroutine GRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv,i
     1pbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles
c with various boundary conditions.
c scalar version using guard cells
c 90 flops/particle, 2 divides, 2 sqrts, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,j) = dx
      part(4,j) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nxyv
     1,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles
c with various boundary conditions.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290  (1996).
c 90 flops/particle, 2 divides, 2 sqrts, 22 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
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
      qtmh = .5*qbm*dt
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,j) = dx
      part(4,j) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,nop) = dx
      part(4,nop) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv,
     1ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles
c with various boundary conditions.
c scalar version using guard cells
c 52 flops/particle, 2 divides, 2 sqrts, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,j) = dx
      part(4,j) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nxy
     1v,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles
c with various boundary conditions.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 52 flops/particle, 2 divides, 2 sqrts, 12 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
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
      qtmh = .5*qbm*dt
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,j) = dx
      part(4,j) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,nop) = dx
      part(4,nop) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny,
     1nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 200 flops/particle, 4 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2))
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + (pz(t-dt/2)**2)/
c      (1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 200 flops/particle, 4 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2))
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + (pz(t-dt/2)**2)/
c      (1. + gami)
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
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
      subroutine GRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 117 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2))
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + (pz(t-dt/2)**2)/
c      (1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 117 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2))
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)),
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + (pz(t-dt/2)**2)/
c      (1. + gami)
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
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 223 flops/particle, 4 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 223 flops/particle, 4 divides, 2 sqrts, 50 loads, 5 stores
c input: all, output: part, ek
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,
     1ny,nxv,nxyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
c input: all, output: part, ek
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
c bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
      acz = part(5,nop) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny,
     1nxv,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c in z direction, using the Boris Mover.
c scalar version using guard cells
c 153 flops/particle, 4 divides, 2 sqrts, 31 loads, 4 stores
c input: all, output: part, ek
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
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = amy*(dxl*bz(nl,mm) + amx*bz(nn,mm) + dxp*bz(np,mm)) + dyl*(dx
     1l*bz(nl,ml) + amx*bz(nn,ml) + dxp*bz(np,ml)) + dyp*(dxl*bz(nl,mp) 
     2+ amx*bz(nn,mp) + dxp*bz(np,mp))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nxyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c in z direction, using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 153 flops/particle, 4 divides, 2 sqrts, 31 loads, 4 stores
c input: all, output: part, ek
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
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j+1,k+1) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gami)
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
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = amy*(dxl*bz(mn) + amx*bz(mn+1) + dxp*bz(mn+2)) + dyl*(dxl*bz(
     1ml) + amx*bz(ml+1) + dxp*bz(ml+2)) + dyp*(dxl*bz(mp) + amx*bz(mp+1
     2) + dxp*bz(mp+2))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = amy*(dxl*bz(mn) + amx*bz(mn+1) + dxp*bz(mn+2)) + dyl*(dxl*bz(
     1ml) + amx*bz(ml+1) + dxp*bz(ml+2)) + dyp*(dxl*bz(mp)+ amx*bz(mp+1)
     2 + dxp*bz(mp+2)) 
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c in z direction, using the Boris Mover.
c scalar version using guard cells
c 74 flops/particle, 4 divides, 2 sqrts, 16 loads, 4 stores
c input: all, output: part, ek
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
c where omz = (q/m)*bz(x(t),y(t))*gami,
c and gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyp*(dxp*bz(np,mp) + amx*bz(nn,mp)) + amy*(dxp*bz(np,mm) + am
     1x*bz(nn,mm))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GSRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nxyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles with magnetic field
c in z direction, using the Boris Mover.
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 74 flops/particle, 4 divides, 2 sqrts, 16 loads, 4 stores
c input: all, output: part, ek
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
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c fx(x(t),y(t)), fy(x(t),y(t)), bz(x(t),y(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1,m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c similarly for fy(x,y), bx(x,y), by(x,y), bz(x,y)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(j,k) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gami)
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
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyp*(dxp*bz(mp+1) + amx*bz(mp)) + amy*(dxp*bz(mm+1) + amx*bz
     1(mm))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,nop) + dx
      acy = part(4,nop) + dy
c find inverse gamma
      p2 = acx*acx + acy*acy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      oz = dyn*(dxn*bz(mp+1) + amx*bz(mp)) + amy*(dxn*bz(mm+1) + amx*bz
     1(mm))
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop) + dx*dtg
      dy = part(2,nop) + dy*dtg
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
   20 ek = ek + sum1
      return
      end
      subroutine GRPUSH2C(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv,
     1ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order spline
c interpolation in space, for relativistic particles
c with various boundary conditions.
c scalar version using guard cells
c 128 flops/particle, 2 divides, 2 sqrts, 36 loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fx(j,k) = x component of force/charge at grid (j,k)
c fy(j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c qbm = particle charge/mass
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nxy = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, qbm, dt, ci, ek
      dimension part(idimp,nop), fxy(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, edgelx, edgely, edgerx, edgery
      real dx, dy, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq, at1, at2
      real ci2, acx, acy, p2, dtg
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,j) = dx
      part(4,j) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GRBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 285 flops/particle, 4 divides, 2 sqrts, 85 loads, 5 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(3)*(pz(t-dt/2)) +.5*(q/m)*fx(x(t),y(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(6)*(pz(t-dt/2)) + .5*(q/m)*fy(x(t),y(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
c    rot(9)*(pz(t-dt/2))
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + (pz(t-dt/2)**2)/
c      (1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nyv = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, bxy, qbm, dt, dtc, ci, ek
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
      real ci2, p2, gami, qtmg, dtg
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j)
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GRBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order linear
c interpolation in space, for relativistic particles with magnetic field
c Using the Boris Mover.
c scalar version using guard cells
c 323 flops/particle, 4 divides, 2 sqrts. 101 loads, 5 stores
c input: all, output: part, ek
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
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
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nyv = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, bxy, qbm, dt, dtc, ci, ek
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
      real ci2, p2, gami, qtmg, dtg
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GRBPUSH22C(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and third-order spline
c interpolation in space, for relativistic particles with magnetic field
c in z direction, using the Boris Mover.
c scalar version using guard cells
c 178 flops/particle, 4 divides, 2 sqrts, 52 loads, 4 stores
c input: all, output: part, ek
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
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
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
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c bz(+2j,k+2) = z component of magnetic field at grid (j,k)
c that is, the convolution of magnetic field over particle shape
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxv = first dimension of field arrays, must be >= nx+5
c nyv = second dimension of field arrays, must be >= ny+5
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real part, fxy, bz, qbm, dt, dtc, ci, ek
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, edgelx, edgely, edgerx, edgery
      real dx, dy, oz, dxl, dxn, dxp, dxq, dyl, dyn, dyp, dyq, at1, at2
      real acx, acy, omzt, omt, anorm, rot1, rot2
      real ci2, p2, gami, qtmg, dtg
      double precision sum1
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
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
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine RRETARD2(part,dtc,ci,idimp,nop,nx,ny,ipbc)
c for 2-1/2d code, particle positions are retarded a half time-step,
c for relativistic particles
c input: all, output: part
c equations used are:
c x(t+dt) = x(t) - px(t+dt/2)*dtg
c y(t+dt) = y(t) - py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
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
      do 10 j = 1, nop
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c retard position half a time-step for current deposit
      dx = part(1,j) - vx*dtg
      dy = part(2,j) - vy*dtg
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
      subroutine RRETARD22(part,dtc,ci,idimp,nop,nx,ny,ipbc)
c for 2d code, particle positions are retarded a half time-step
c input: all, output: part
c equations used are:
c x(t+dt) = x(t) - px(t+dt/2)*dtg
c y(t+dt) = y(t) - py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      dimension part(idimp,nop)
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
      do 10 j = 1, nop
c find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      dtg = dtc/sqrt(1.0 + p2*ci2)
c retard position half a time-step for current deposit
      dx = part(1,j) - vx*dtg
      dy = part(2,j) - vy*dtg
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
      subroutine CPTOV2(part,ci,nop,idimp)
c for 2-1/2d code, this subroutine converts momentum to velocity for
c relativistic particles
c part(i,n) = part(i,n)/sqrt(1.+sum(part(i,n)**2)*ci*ci)
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
      implicit none
      integer nop, idimp
      real part, ci
      dimension part(idimp,nop)
      integer j
      real ci2, vx, vy, vz, p2, gami
      ci2 = ci*ci
c convert from momentum to velocity
      do 10 j = 1, nop
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      part(3,j) = vx*gami
      part(4,j) = vy*gami
      part(5,j) = vz*gami
   10 continue
      return
      end
      subroutine CPTOV22(part,ci,nop,idimp)
c for 2d code, this subroutine converts momentum to velocity for
c relativistic particles
c part(i,n) = part(i,n)/sqrt(1.+sum(part(i,n)**2)*ci*ci)
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 4
      implicit none
      integer nop, idimp
      real part, ci
      dimension part(idimp,nop)
      integer j
      real ci2, vx, vy, p2, gami
      ci2 = ci*ci
c convert from momentum to velocity
      do 10 j = 1, nop
      vx = part(3,j)
      vy = part(4,j)
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
      part(3,j) = vx*gami
      part(4,j) = vy*gami
   10 continue
      return
      end
      subroutine RPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
c for 2d code, this subroutine updates particle co-ordinates for
c relativistic particles with fixed velocities,
c with various boundary conditions.
c 11 flops/particle, 2 divides, 1 sqrts, 4 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg, y(t+dt) = y(t) + py(t+dt/2)*dtg,
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2))**2 + (py(t-dt/2))**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t+dt/2)**2+py(t+dt/2)**2)*ci*ci)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, ipbc
      real part, dt, ci, ek
      dimension part(idimp,nop)
      integer j
      real ci2, edgelx, edgely, edgerx, edgery, dx, dy, p2, gam, dtg
      double precision sum1
      ci2 = ci*ci
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
c kinetic energy
      p2 = dx*dx + dy*dy
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
      dtg = dt/gam
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine RPUSH23ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates for
c relativistic particles with fixed velocities,
c with various boundary conditions.
c 13 flops/particle, 2 divides, 1 sqrts, 5 loads, 2 stores
c input: all, output: part, ek
c equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg, y(t+dt) = y(t) + py(t+dt/2)*dtg,
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c part(5,n) = momentum pz of particle n
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2))**2+(py(t-dt/2))**2+(pz(t-dt/2))**2)/(1. + gamma)
c where gamma=sqrt(1.+(px(t+dt/2)**2+py(t+dt/2)**2+pz(t+dt/2)**2)*ci*ci)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, ipbc
      real part, dt, ci, ek
      dimension part(idimp,nop)
      integer j
      real ci2, edgelx, edgely, edgerx, edgery, dx, dy, dz, p2, gam, dtg
      double precision sum1
      ci2 = ci*ci
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
      dz = part(5,j)
c kinetic energy
      p2 = dx*dx + dy*dy + dz*dz
      gam = sqrt(1.0 + p2*ci2)
      sum1 = sum1 + p2/(1.0 + gam)
      dtg = dt/gam
c new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine RDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nxvh,n
     1yv,ipbc)
c for 2-1/2d code, this subroutine calculates particle current density
c using gridless spectral version for relativistic particles,
c periodic boundaries
c in addition, particle positions are advanced a half time-step
c baseline scalar version
c 36*(NX/2)*(NY/2) + 17 flops/particle, 1 divide, 1, sqrt,
c 3*NX*NY/2 loads and stores
c plus (NX/2 + NY/2) sines and cosines/particle
c input: all, output: part, cu
c charge density is calculated from the expression:
c cu(i,n,m) = sum(qmi*exp(-sqrt(-1)*2*n*pi*x/nx)*
c                exp(-sqrt(-1)*2*m*pi*y/ny))
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = x momentum of particle n
c part(4,n) = y momentum of particle n
c part(5,n) = z momentum of particle n
c cu(i,n,m) = charge density at fourier grid point n,m
c sctx = scratch array for sines and cosines
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c nop = number of particles
c idimp = size of phase space = 5
c nx/ny = system length in x/y direction
c nxvh = first dimension of current array, must be >= nx/2
c nyv = second dimension of current array, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc
      real part, qm, dt, ci
      complex cu, sctx
      dimension part(idimp,nop), cu(3,nxvh,nyv), sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real qmn, dnx, dny, dkx, dky, at1, at3
      real vx, vy, vz, edgelx, edgely, edgerx, edgery, dx, dy
      real ci2, p2, gami, px, py, pz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qmn = qm/real(nx*ny)
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
c find fourier components
      do 50 i = 1, nop
c find inverse gamma
      vx = part(3,i)
      vy = part(4,i)
      vz = part(5,i)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
      px = vx*gami
      py = vy*gami
      pz = vz*gami
      vx = qmn*px
      vy = qmn*py
      vz = qmn*pz
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
      dx = part(1,i) + px*dt
      dy = part(2,i) + py*dt
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
      subroutine RPUSH2GL(part,fxy,sctx,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv
     1h,nyv,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, periodic boundaries and various particle boundary conditions,
c for relativistic particles.
c scalar version using guard cells
c 28*(NX/2)*(NY/2) + 22 flops/particle, 2 divides, 2 sqrts
c NX*NY loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c fy(x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
c idimp = size of phase space = 4
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nxy = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt, ci, ek
      complex fxy, sctx
      dimension part(idimp,nop), fxy(2,nxvh,nyv), sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, acx, acy, p2, ci2, dtg
      double precision sum1, exl, eyl, ex, ey
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,i) + dx
      acy = part(4,i) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,i) = dx
      part(4,i) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,i) + dx*dtg
      dy = part(2,i) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine RPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,nop,nx
     1,ny,nxvh,nyv,npp,ipbc)
c for 2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, periodic boundaries and various particle boundary conditions,
c for relativistic particles.
c scalar version using guard cells
c 28*(NX/2)*(NY/2) + 22 flops/particle, 2 divides, 2 sqrts
c NX*NY loads, 4 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2))*ci*ci)
c fx(x(t),y(t)) and fy(x(t),y(t)) are calculated from the expression
c fx(x,y) = sum(fxy(1,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c fy(x,y) = sum(fxy(2,n,m)*exp(sqrt(-1)*2*n*pi*x/nx)*
c               exp(sqrt(-1)*2*m*pi*y/ny))
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = momentum px of particle n
c part(4,n) = momentum py of particle n
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c sctxp = scratch array for sines and cosines
c exyp = scratch array for particle forces
c qbm = particle charge/mass
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t))*ci*ci)
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
      real part, qbm, dt, ci, ek
      double precision exyp
      complex fxy, sctxp
      dimension part(idimp,nop), fxy(2,nxvh,nyv), sctxp(nxvh,npp)
      dimension exyp(2,npp)
      integer i, j, k, l, k1, nxh, nyh, ny2, npb, ipp, ie, ib
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, acx, acy, p2, ci2, dtg
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
c half acceleration
      acx = part(3,i+ib) + dx
      acy = part(4,i+ib) + dy
c time-centered kinetic energy
      p2 = acx*acx + acy*acy
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      part(3,i+ib) = dx
      part(4,i+ib) = dy
c update inverse gamma
      p2 = dx*dx + dy*dy
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,i+ib) + dx*dtg
      dy = part(2,i+ib) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine RBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,nop
     1,nx,ny,nxvh,nyv,ipbc)
c for 2-1/2d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time using gridless spectral
c version, for relativistic particles with magnetic field, periodic
c boundaries and various particle boundary conditions.
c Using the Boris Mover.
c scalar version using guard cells
c 60*(NX/2)*(NY/2) + 69 flops/particle, 4 divides, 2 sqrts,
c 3*NX*NY + 5 loads, 5 stores
c input: all, output: part, ek
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
c omz = (q/m)*bz(x(t),y(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c sctx = scratch array for sines and cosines
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c ci = reciprocal of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
c idimp = size of phase space = 5
c nop = number of particles
c nx/ny = system length in x/y direction
c nxvh = first dimension of field arrays, must be >= nx/2
c nyv = second dimension of field arrays, must be >= ny
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc
      real part, qbm, dt, dtc, ci, ek
      complex fxy, bxy, sctx
      dimension part(idimp,nop), fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension sctx(nxvh)
      integer i, j, k, k1, nxh, nyh, ny2
      real edgelx, edgely, edgerx, edgery, dnx, dny, dkx, dky, at1, at3
      real qtmh, dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt
      real omt, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real anorm, ci2, p2, gami, qtmg, dtg
      double precision sum1, exl, eyl, ezl, ex, ey, ez
      double precision bxl, byl, bzl, bx, by, bz
      complex zt1, zt2, zt3, zt4
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,i) + dx*dtg
      dy = part(2,i) + dy*dtg
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
      ek = ek + sum1
      return
      end
      subroutine GRCJPOST2(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates particle current density
c using second-order spline interpolation for relativistic particles.
c scalar version using guard cells
c 103 flops/particle, 1 sqrt, 40 loads, 18 stores
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
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
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
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j+1,k+1) = x component of force/charge at grid (j,k)
c fxy(2,j+1,k+1) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c cu(i,j+1,k+1) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of current density array, must be >= nx+3
c nyv = second dimension of current density array, must be >= ny+3
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, cu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), cu(2,nxv,nyv)
      integer j, nn, mm, nl, np, ml, mp
      real qtmh, ci2, gami, dxp, dyp, amx, amy, dxl, dyl
      real dx, dy, dz, vx, vy, p2
      qtmh = 0.5*qbm*dt
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
c new momentum
      vx = part(3,j) + qtmh*dx
      vy = part(4,j) + qtmh*dy
c find inverse gamma
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c deposit current density
      amx = qm*amx
      dxl = qm*dxl
      dxp = qm*dxp
      vx = gami*vx
      vy = gami*vy
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
      subroutine GRCJPOST2L(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates particle current density
c using first-order spline interpolation for relativistic particles.
c scalar version using guard cells
c 59 flops/particle, 1 sqrt, 20 loads, 0 stores
c input: all, output: cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m)=qci*(1.-dx)*(1.-dy)
c cu(i,n+1,m)=qci*dx*(1.-dy)
c cu(i,n,m+1)=qci*(1.-dx)*dy
c cu(i,n+1,m+1)=qci*dx*dy
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c where n,m = nearest grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
c fx(x(t),y(t)) and fy(x(t),y(t)) are approximated by interpolation from
c the nearest grid points:
c fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
c    + dx*fx(n+1.m+1))
c fy(x,y) = (1-dy)*((1-dx)*fy(n,m)+dx*fy(n+1,m)) + dy*((1-dx)*fy(n,m+1)
c    + dx*fy(n+1.m+1))
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j,k) = x component of force/charge at grid (j,k)
c fxy(2,j,k) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c cu(i,j,k) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of current density array, must be >= nx+1
c nyv = second dimension of current density array, must be >= ny+1
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, cu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), cu(2,nxv,nyv)
      integer j, nn, mm, np, mp
      real qtmh, ci2, gami, dxp, dyp, amx, amy, dx, dy, vx, vy, p2
      qtmh = 0.5*qbm*dt
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
c new momentum
      vx = part(3,j) + qtmh*dx
      vy = part(4,j) + qtmh*dy
c find inverse gamma
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c deposit current density
      amx = qm*amx
      dxp = qm*dxp
      vx = gami*vx
      vy = gami*vy
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
      subroutine GRCJPOST2C(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv)
c for 2d code, this subroutine calculates particle current density
c using third-order spline interpolation for relativistic particles.
c scalar version using guard cells
c 201 flops/particle, 1 sqrt, 68 loads, 32 stores
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
c and qci = qm*pj*gami, where j = x,y, for i = 1, 2
c where pj = .5*(pj(t+dt/2)+pj(t-dt/2))
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c momentum equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t))*dt,
c where q/m is charge/mass, and
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
c and similarly for fy(x,y).
c where n,m = leftmost grid points and dx = x-n, dy = y-m
c part(1,n) = position x of particle n at t
c part(2,n) = position y of particle n at t
c part(3,n) = momentum px of particle n at t - dt/2
c part(4,n) = momentum py of particle n at t - dt/2
c fxy(1,j+2,k+2) = x component of force/charge at grid (j,k)
c fxy(2,j+2,k+2) = y component of force/charge at grid (j,k)
c that is, convolution of electric field over particle shape
c cu(i,j+2,k+2) = ith component of current density
c at grid point j,k for i = 1, 2
c qm = charge on particle, in units of e
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c ci = reciprocal of velocity of light
c idimp = size of phase space = 4
c nop = number of particles
c nxv = first dimension of current density array, must be >= nx+5
c nyv = second dimension of current density array, must be >= ny+5
      implicit none
      integer idimp, nop, nxv, nyv
      real part, fxy, cu, qm, qbm, dt, ci
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), cu(2,nxv,nyv)
c local data
      real sixth, two3rds
      parameter(sixth=0.166666666666667,two3rds=4.0*sixth)
      integer j, nl, nn, np, nq, ml, mm, mp, mq
      real qtmh, dxl, dyl, dxn, dyn, dxp, dyp, dxq, dyq, at1, at2
      real dx, dy, dz, dw, vx, vy, ci2, gami, p2
      qtmh = 0.5*qbm*dt
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
c new momentum
      vx = part(3,j) + qtmh*dx
      vy = part(4,j) + qtmh*dy
c find inverse gamma
      p2 = vx*vx + vy*vy
      gami = 1.0/sqrt(1.0 + p2*ci2)
c deposit current density
      dxl = qm*dxl
      dxn = qm*dxn
      dxp = qm*dxp
      dxq = qm*dxq
      vx = gami*vx
      vy = gami*vy
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
