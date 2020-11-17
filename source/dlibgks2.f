c 2-1/2d PIC library for diagnostics graphics
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: march 1, 2012
      subroutine GRASP23(part,label,itime,isc,nx,ny,iyp,ixp,idimp,npxy,n
     1p,irc)
c for 2-1/2d code, this subroutine displays (iyp-ixp) phase space
c plots background particles in blue and beam particles in red
c part(1,n) = position x of particle n
c part(2,n) = position y of particle n
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c label = species label
c itime = current time step
c isc = power of 2 scale of range of values of velocity
c nx/ny = system length in x/y direction
c iyp/ixp = phase space coordinates to be displayed
c idimp = size of phase space = 5
c npxy = number of background particles
c np = number of particles
c irc = return code (0 = normal return)
c nxbs = length of scratch variable for plotting
      parameter(nxbs=65)
      dimension part(idimp,np)
      character*(*) label
c idwk = workstation identifier
c ncols = number of foreground colors available for line plotting
c rx, ry = ndc coordinates of upper-right corner of workstation window
c iplot = plot location on page, 0 <= iplot < nplot
c nplot = number of plots per page
c iclr = (-1,0,1) = (no,default,yes) erase plot (default=when iplot=0)
c iupd = (-1,0,1) = (no,default,yes) end plot
c (default=when iplot=nplot-1)
c kprime = table of color indices for prime colors
      common /plotcm/ idwk,ncols,rx,ry,iplot,nplot,iclr,iupd,idstr,idloc
     1,nclsp,ifrg,isx,isy,kprime(8)
c x,y = scratch arrays for plotting
      dimension x(nxbs), y(nxbs)
      character*26 chr
      character*2 lblsp(5)
      character*10 chrs(2)
      save lblsp
   91 format (1x,a2,' VERSUS ',a2,', T = ',i7)
      data lblsp /' X',' Y','VX','VY','VZ'/
c dv = scale will be set in powers of this parameter
      data dv /2.0/
c smin/smax = range of x values of plotting window
c tmin/tmax = range of y values of plotting window
      data smin,tmin,smax,tmax /.25,.14,.975,.975/
c ntx/nty = number of ticks in grid in x/y direction
      data ntx,nty /9,9/
c csize = vertical size of characters
      data zero,csize /0.,0.034/
c sample labels
      data chrs /'BACKGROUND','   BEAM   '/
c set return code to normal
      irc = 0
c exit if plots are suppressed (nplot = 0)
      if (nplot.lt.1) return
      algdvi = 1./alog(dv)
c find y scale for plot
      if (iyp.le.1) then
         if (iyp.lt.1) iyp = 1
         ymin = zero
         ymax = float(nx)
      elseif (iyp.eq.2) then
         ymin = zero
         ymax = float(ny)
      elseif (iyp.gt.2) then
         if (iyp.gt.idimp) iyp = idimp
         is = isc
         if (abs(is).gt.116) then
            ymax = abs(part(iyp,1))
            do 10 j = 1, np
            ymax = amax1(ymax,abs(part(iyp,j)))
   10       continue
            if (ymax.eq.0.) ymax = 1.0e-35
            is = alog(ymax)*algdvi
            if (ymax.ge.1.) is = is + 1
            if (ymax.le.dv**(is-1)) is = is - 1
         endif
         ymax = dv**is
         ymin = -ymax
      endif
c find x scale for plot
      if (ixp.le.1) then
         if (ixp.lt.1) ixp = 1
         xmin = zero
         xmax = float(nx)
      elseif (ixp.eq.2) then
         xmin = zero
         xmax = float(ny)
      elseif (ixp.gt.2) then
         if (ixp.gt.idimp) ixp = idimp
         is = isc
         if (abs(is).gt.116) then
            xmax = abs(part(ixp,1))
            do 20 j = 1, np
            xmax = amax1(xmax,abs(part(ixp,j)))
   20       continue
            if (xmax.eq.0.) xmax = 1.0e-35
            is = alog(xmax)*algdvi
            if (xmax.ge.1.) is = is + 1
            if (xmax.le.dv**(is-1)) is = is - 1
         endif
         xmax = dv**is
         xmin = -xmax
      endif
c find location for plot
      npl1 = sqrt(float(nplot-1)) + 0.0001
      npl = npl1 + 1
      apl = 1./float(npl)
      iy = iplot/npl
      ix = iplot - iy*npl
      aplx = apl*rx
      aply = apl*ry
      orx = aplx*float(ix)
      ory = aply*float(npl1 - iy)
      smn = orx + aplx*smin
      smx = orx + aplx*smax
      tmn = ory + aply*tmin
      tmx = ory + aply*tmax
      chh = aply*csize
c initiate plot
      if (((iplot.eq.0).and.(iclr.eq.0)).or.(iclr.eq.1)) then
c clear workstation, always
         call gclrwk(idwk,1)
      endif
c write labels
      write (chr,91) lblsp(iyp), lblsp(ixp), itime
c select point as marker
      mks = 1
c draw grid and labels, call identity transformation
      call tickd(xmin,xmax,ymin,ymax,orx,ory,smn,smx,tmn,tmx,ntx,nty,lab
     1el,chr,chh)
c display sample markers
      ngs = 2
      if ((npxy.eq.0).or.(np.eq.npxy)) ngs = 1
      call dsmpln(orx,smn,tmn,tmx,ngs,1,1,mks,chrs,chh)
c define transformation number 2
      nrt = 1
c set window
      call gswn(nrt,xmin,xmax,ymin,ymax)
c set viewport
      call gsvp(nrt,smn,smx,tmn,tmx)
c select normalization transformation
      call gselnt(nrt)
c set marker size scale factor, 1.0 = nominal
      call gsmksc(1.0)
c set marker type, 1 = dot, 2 = plus, 3 = star, 4 = circle, 5 = cross
      call gsmk(mks)
c set clipping indicator, 1 = on
      call gsclip(1)
c plot particles
      do 50 k = 1, 2
c determine how many particles to plot
      if (k.eq.1) then
         nd = npxy
      elseif (k.eq.2) then
         nd = np - npxy
      endif
      if (nd.eq.0) go to 50
      koff = npxy*(k - 1)
      icol = k + 1 - ncols*(k/ncols)
      icol = kprime(icol+1)
c set polymarker color index
c 1=foreground, 2=blue, 3=red, 4=yellow, 5=cyan, 6=magenta, 7=green
      call gspmci(icol)
c parameters for plots
      nxs = nxbs - 1
c blocking parameters for plots
      nxb = (nd - 1)/nxs + 1
      npts = nxs
c length of last block
      nptl = nd - nxs*(nxb - 1)
c plot polymarkers
      npt = npts
c loop over number of blocks
      do 40 j = 1, nxb
      js = nxs*(j - 1) + koff
      if (j.eq.nxb) npt = nptl
c calculate x,y axes for block
      do 30 i = 1, npt
      x(i) = part(ixp,i+js)
      y(i) = part(iyp,i+js)
   30 continue
c treat dots by drawing a line to itself
c     call spdots(x,y,npt,icol,nxbs)
c treat dots by drawing a line to itself with non-zero width in x
      dd = (xmax - xmin)*float(npl)*2.0e-3
      call spddots(x,y,dd,npt,icol,nxbs)
   40 continue
   50 continue
c update workstation, perform
      call guwk(idwk,1)
c update plot number
      iplot = iplot + 1
      if (iplot.eq.nplot) iplot = 0
      if (((iplot.eq.0).and.(iupd.eq.0)).or.(iupd.eq.1)) then
c update workstation, perform
         call guwk(idwk,1)
c read code from input device, if present
         call readrc(irc)
      endif
c reset defaults
      iclr = 0
      iupd = 0
      return
      end
      subroutine spddots(x,y,dd,npt,icol,nxbs)
c this subroutine draws dot markers by drawing a line to itself
c with non-zero width in x
c x, y = arrays to be plotted
c dd = smallest visible size, such as (xmax - xmin)*4.0e-3
c npt = number of points to be plotted
c icol = color index
c nxbs = dimension of x, y arrays
      dimension x(nxbs), y(nxbs)
c xs, ys = scratch arrays for plotting
      dimension xs(2), ys(2)
c set polyline color index
      call gsplci(icol)
      do 10 j = 1, npt
      xs(1) = x(j)
      ys(1) = y(j)
      xs(2) = xs(1) + dd
      ys(2) = ys(1)
c draw polyline
      call gpl(2,xs,ys)
   10 continue
      return
      end
