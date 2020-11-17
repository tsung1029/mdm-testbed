!-----------------------------------------------------------------------
! * * * periodic 2d darwin particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with periodic
! electromagnetic forces obtained by solving non-radiative form of
! maxwell's equation with fast fourier transforms
! using hamiltonian algorithm described in C.W. Nielsen and H.R. Lewis,
! Methods in Computational Physics, vol 16, p. 367 (1976).
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 2006, regents of the university of california
! update: march 1, 2012
      program hdbeps2
      use init2d
      use empush2d
      use ehdpush2d
      use field2d
      use diag2d
      use hdsimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: nmv = 40, nma = 40, ipbc = 1
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, iua = 15, ium = 21, iuer = 2
      integer :: npxy, npxyb, np, npxyi, npxybi, npi, npxy1
      integer :: nx, ny, nxh, nyh, nxe, nye
      integer :: nxye, nxeh, nxyh, nxhy, nx1, ny1, ntmax, idimp, ndc2
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: j, k, it, jt, itw, modesy2d, modesy2p, modesy2a
      integer :: iur1, iur2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, wf = 0.0, time = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0 
      real :: tfft = 0.0, tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: tdhjpost = 0.0, tdhjposti = 0.0, ts = 0.0
      real :: qbme, qbmi, affp, dth, qi0, q2m0, wp0, wpmax, wpmin
      real :: etx, we, wm, wke
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: sx, sy, sz, pxe, pye, pze, wx, wy, wz
      real :: wef, wtot, vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      real, dimension(:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:), pointer :: qe, qi
      real, dimension(:,:,:), pointer :: fxyze, axyze, daxyze
      real, dimension(:,:,:), pointer :: cu, cus
      complex, dimension(:,:,:), pointer :: q2m
      complex, dimension(:,:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:,:), pointer :: qt, sfield
      complex, dimension(:,:), pointer :: dent, pott
      real, dimension(:,:,:), pointer :: vfield
      complex, dimension(:,:,:), pointer :: vpott
      real, dimension(:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: fa, fam, fai, fami
      real, dimension(:,:), pointer :: wt
      character(len=10) :: cdrun
      character(len=32) :: fname
      integer, external :: NDIAN, NDPREC, IDPREC
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
  993 format (' electric(l,t), magnetic energies = ',3e14.7)
  996 format (' total momentum = ',3e14.7)
! get unit number for error file
      iuer = get_funit(iuer)
! read namelist
      iuin = get_funit(iuin)
      open(unit=iuin,file='input2',form='formatted',status='old')
      read (iuin,input2)
      if (movion==1) then
         rewind iuin
         read (iuin,ions2)
      endif
! override input data
      psolve = 1
! ndc2 = number of corrections in secondary darwin iteration
      ndc2 = ndc
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      iuot = get_funit(iuot)
      fname = 'output2.'//cdrun
      open(unit=iuot,file=trim(fname),form='formatted',status='replace')
! open initial diagnostic metafile
      iudm = get_funit(iudm)
      fname = 'diag2.init.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! np = total number of electrons in simulation
      npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb
! npi = total number of ions in simulation
      npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 4; nye = ny + 3
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871
         nxe = nx + 2; nye = ny + 1
      endif
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
      call MP_SETSTACK(262144)
      nxye = nxe*nye; nxeh = nxe/2;
! nxyh = maximum(nx,ny)/2
      nxyh = max(nx,ny)/2
! nxhy = maximum(nx/2,ny)
      nxhy = max(nxh,ny)
! dimensions for index and sorting arrays
      nx1 = nx + 1; ny1 = ny + 1
! ntmax = size of buffer for removing particles
      ntmax = 1 + (4*(npxy*vty + npxyb*vtdy) + 2.8*npxyb*abs(vdy))*dt/ny
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! idimp = dimension of phase space = 6 or 8
      idimp = 2 + 2*ndim
      if (relativity==1) idimp = idimp + 1
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
! part(6,n) = canonical momentum px of particle n
! part(7,n) = canonical momentum py of particle n
! part(8,n) = canonical momentum pz of particle n
      allocate(part(idimp,np))
! in real space, qe(j,k) = charge density at grid point (j,k)
! in real space, qi(j,k) = ion charge density at grid point (j,k)
      allocate(qe(nxe,nye),qi(nxe,nye))
! fxyze(i,j,k) = i component of longitudinal electric force/charge
! axyze(i,j,k) = i component of smoothed vector potential
      allocate(fxyze(ndim,nxe,nye),axyze(ndim,nxe,nye))
! daxyze(i,j,k) = i component of derivative of smoothed vector potential
      allocate(daxyze(2*ndim-1,nxe,nye))
! cu(i,j,k) = i component of current density
! cus(i,j,k) = i component of modified current density
      allocate(cu(ndim,nxe,nye),cus(ndim,nxe,nye))
! store shift constants
      allocate(q2m(1,1,1))
! ffc, ffe = form factor array for poisson solvers
      allocate(ffc(nxeh,nyh),ffe(nxeh,nyh))
! mixup, sct = arrays for fft
      allocate(mixup(nxhy),sct(nxyh))
! sorting arrays
      allocate(pt(max(np,npi)),ip(max(np,npi)),npic(ny1))
      if (sortime > 0) then
         allocate(part2(idimp,np))
      else
         allocate(part2(0,0))
      endif
!
! open graphics device
      call GROPEN
      call SETNPLT(nplot,irc)
      call STPALIT(idpal)
! initialize timer
      call wtimer(time,ltime,-1)
! initialize constants
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
      qbme = qme
      affp = real(nx*ny)/real(np)
      qbmi = zero
      dth = .5*dt
      qbmi = zero
      q2m0 = qbme*qme*real(np)/real(nx*ny)
      if (movion==1) then
         qbmi = qmi/rmass
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
         vtdzi = vtdz/sqrt(rmass*rtempdzi)
         q2m0 = q2m0 + qbmi*qmi*real(npi)/real(nx*ny)
      endif
      wp0 = q2m0*affp
      q2m = cmplx(q2m0,wp0)
!
! diagnostics
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! ion density or potential diagnostics
      if ((ntp > 0) .or. (ndp > 0) .or. (ntd > 0) .or. (ndd > 0)) then
         allocate(sfield(nxe,nye))
      endif
! ion density diagnostics
      if (ntd > 0) then
         if (modesxd > nxh) modesxd = nxh
         if (modesyd > nyh) modesyd = nyh
         allocate(qt(nxe,nye))
         modesy2d = 2*modesyd - 1
         allocate(dent(modesxd,modesy2d))
         ceng = 0.0
         fdname = 'denk2.'//cdrun
         write (iudm,den2d,iostat=irc)
      endif
! velocity diagnostics
      if ((ntv > 0) .or. (ndv > 0)) then
         allocate(fv(2*nmv+2,ndim),fvm(2,ndim))
         allocate(fa(2*nma+2,ndim),fam(2,ndim))
! fix velocity range
         if (ndim==2) then
            fv(1,:) = 8.*max(vtx,vty)
            fa(1,:) = 0.125*max(vtx,vty)
         elseif (ndim==3) then
            fv(1,:) = 8.*max(vtx,vty,vtz)
            fa(1,:) = 0.125*max(vtx,vty,vtz)
         endif
         if (ntv > 0) then
            iuv = get_funit(iuv)
            fname = 'fv2.'//cdrun
            open(unit=iuv,file=trim(fname),form='formatted',status='unkn&
     &own')
! write captions
            if (ndim==2) then
               write (iuv,*) 'it vdx vdy vtx vty'
            else if (ndim==3) then
               write (iuv,*) 'it vdx vdy vdz vtx vty vtz'
            endif
         endif
         if (movion==1) then
            allocate(fvi(2*nmv+2,ndim),fvmi(2,ndim))
            allocate(fai(2*nma+2,ndim),fami(2,ndim))
! fix velocity range 
            if (ndim==2) then
               fvi(1,:) = 8.*max(vtxi,vtyi)
               fai(1,:) = 0.125*max(vtxi,vtyi)
            else if (ndim==3) then
               fvi(1,:) = 8.*max(vtxi,vtyi,vtzi)
               fai(1,:) = 0.125*max(vtxi,vtyi,vtzi)
            endif
            if (ntv > 0) then
               iuvi = get_funit(iuvi)
               fname = 'fvi2.'//cdrun
               open(unit=iuvi,file=trim(fname),form='formatted',status='&
     &unknown')
! write captions
               if (ndim==2) then
                  write (iuvi,*) 'it vdxi vdyi vtxi vtyi'
               else if (ndim==3) then
                  write (iuvi,*) 'it vdxi vdyi vdzi vtxi vtyi vtzi'
               endif
            endif
         endif
      endif
! potential diagnostics
      if (ntp > 0) then
         if (modesxp > nxh) modesxp = nxh
         if (modesyp > nyh) modesyp = nyh
         modesy2p = 2*modesyp - 1
         allocate(pott(modesxp,modesy2p))
         ceng = affp
         fpname = 'potk2.'//cdrun
         write (iudm,pot2d,iostat=irc)
      endif
! vector potential diagnostics
      if ((nta > 0) .or. (nda > 0)) then
         allocate(vfield(ndim,nxe,nye))
         if (nta > 0) then
            if (modesxa > nxh) modesxa = nxh
            if (modesya > nyh) modesya = nyh
            modesy2a = 2*modesya - 1
            allocate(vpott(ndim,modesxa,modesy2a))
            ceng = affp
            faname = 'vpotk2.'//cdrun
            write (iudm,vpot2d,iostat=irc)
         endif
      endif
! energy diagnostics
      if (ndw > 0) then
         allocate(wt((nloop-1)/ndw-(itime0/ndw)+1,7))
         itw = 0
      endif
! momentum diagnostics
      if (ntm > 0) then
         ium = get_funit(ium)
         fname = 'momentum2.'//cdrun
         open(unit=ium,file=trim(fname),form='formatted',status='unknown&
     &')
      endif
! write out input file
      write (iudm,input2,iostat=irc)
! close diagnostic metafile
      close(unit=iudm)
! open restart files
      if (nustrt==0) then
         call restart_open(nustrt,ntr,idrun0,iur1,iur2,iuer)
      else
         call restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
      endif
      if ((iur1 < 0) .and. (iur2 < 0)) go to 3000
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny)
      call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny)
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npi))
         allocate(parti2(0,0))
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         npxy1 = npxy + 1
! background electrons
!        if (npxy > 0) call distr(part,1,npxy,vtx,vty,vx0,vy0,npx,npy,nx&
!    &,ny,ipbc)
!        if (npxy > 0) call distr(part,1,npxy,vtx,vty,vtz,vx0,vy0,vz0,np&
!    &x,npy,nx,ny,ipbc)
         if (npxy > 0) then
            call fdistr(part,1,npxy,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,npx,npy,nx,ny,ipbc,ndprof,nsrand)
            call vdistr(part,1,npxy,vtx,vty,vtz,vx0,vy0,vz0,ndim)
         endif
! beam electrons
!        if (npxyb > 0) call distr(part,npxy1,npxyb,vtdx,vtdy,vdx,vdy,np&
!    &xb,npyb,nx,ny,ipbc)
!        if (npxyb > 0) call distr(part,npxy1,npxyb,vtdx,vtdy,vtdz,vdx,v&
!    &dy,vdz,npxb,npyb,nx,ny,ipbc)
         if (npxyb > 0) then
            call fdistr(part,npxy1,npxyb,ampdx,scaledx,shiftdx,ampdy,sca&
     &ledy,shiftdy,npxb,npyb,nx,ny,ipbc,ndprof,nsrand)
            call vdistr(part,npxy1,npxyb,vtdx,vtdy,vtdz,vdx,vdy,vdz,ndim&
     &)
         endif
! fix guiding centers
!        if (relativity==1) then
!           call distr(part,bxyze,np,qbme,ci,nx,ny,ipbc,inorder)
!        else
!           call distr(part,bxyze,np,qbme,nx,ny,ipbc,inorder)
!        endif
! calculate initial electron momentum
         if (ntm > 0) call initmomt2(part,np,pxe,pye,pze,ndim)
! initialize ions
         if (movion==1) then
            npxy1 = npxyi + 1
! background ions
!           if (npxyi > 0) call distr(parti,1,npxyi,vtxi,vtyi,vxi0,vyi0,&
!    &npxi,npyi,nx,ny,ipbc)
!           if (npxyi > 0) call distr(parti,1,npxyi,vtxi,vtyi,vtzi,vxi0,&
!    &vyi0,vzi0,npxi,npyi,nx,ny,ipbc)
            if (npxyi > 0) then
               call fdistr(parti,1,npxyi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,npxi,npyi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,1,npxyi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,n&
     &dim)
            endif
! beam ions
!           if (npxybi > 0) call distr(parti,npxy1,npxybi,vtdxi,vtdyi,vd&
!    &xi,vdyi,npxbi,npybi,nx,ny,ipbc)
!           if (npxybi > 0) call distr(parti,npxy1,npxybi,vtdxi,vtdyi,vt&
!    &dzi,vdxi,vdyi,vdzi,npxbi,npybi,nx,ny,ipbc)
            if (npxybi > 0) then
               call fdistr(parti,npxy1,npxybi,ampdxi,scaledxi,shiftdxi,a&
     &mpdyi,scaledyi,shiftdyi,npxbi,npybi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,npxy1,npxybi,vtdxi,vtdyi,vtdzi,vdxi,vdy&
     &i,vdzi,ndim)
            endif
! fix guiding centers
!           if (relativity==1) then
!              call distr(parti,bxyze,npi,qbmi,ci,nx,ny,ipbc,inorder)
!           else
!              call distr(parti,bxyze,npi,qbmi,nx,ny,ipbc,inorder)
!           endif
! calculate initial ion momentum
            if (ntm > 0) call initmomt2(parti,npi,pxi,pyi,pzi,ndim)
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
            call dpostg(part,qi,np,nx,ny,-qme,tdpost,inorder,dopt)
! debug
!           call sguard(qi,qi0,nx,ny,inorder)
! freeze the ions now
         else if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density
            call dpostg(parti,qi,npi,nx,ny,qmi,tdposti,inorder,dopt)
! delete ions
            deallocate(parti)
            movion = 0
         endif
! retard electron positionss to deposit current
         call retardg(part,np,dth,ci,nx,ny,ipbc,relativity)
! save retarded co-ordinates
         do j = 1, np
         part(ndim+3,j) = part(1,j)
         part(ndim+4,j) = part(2,j)
         enddo
! initialize current density to background
         call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit current at t-dt/2
         call djpostg(part,cu,np,qme,dth,ci,tdjpost,nx,ny,ipbc,relativit&
     &y,inorder,djopt)
! retard ion positions to deposit current
         if (movion==1) then
            call retardg(parti,npi,dth,ci,nx,ny,ipbc,relativity)
! save retarded ion co-ordinates
            do j = 1, npi
            parti(ndim+3,j) = parti(1,j)
            parti(ndim+4,j) = parti(2,j)
            enddo
! deposit ion current at t-dt/2
            call djpostg(parti,cu,npi,qmi,dth,ci,tdjposti,nx,ny,ipbc,rel&
     &ativity,inorder,djopt)
         endif
! add guard cells for current
         call acguard(cu,nx,ny,inorder)
! transform current to fourier space
         isign = -1
         call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
!        call fftn(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
         call cuperp(cu,nx,ny,inorder)
! calculate smoothed vector potential in fourier space at time t - dt/2
         call sapois(cu,axyze,ffc,ci,tfield,nx,ny,inorder)
! transform vector potential to real space
         isign = 1
         call fft(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
!        call fftn(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
         call cguard(axyze,nx,ny,inorder)
! create electron array based on canonical momentum
         call cr8pcmg(part,axyze,np,qbme,ci,relativity,inorder)
! create ion array based on canonical momentum
         if (movion==1) then
            call cr8pcmg(parti,axyze,npi,qbmi,ci,relativity,inorder)
         endif
! restart
      else
! read primary and, if necessary, secondary restart files
         do j = 1, 2
            if (j==1) then
               it = iur1
            else
               write (iuer,*) 'Initial Restart Error, irc = ', irc
               rewind it
               it = iur2
            endif
            call restart_bread(it,itime,itime0,np,part,movion,npi,parti,&
     &qi,irc,iuer,q2m)
            if (irc /= 0) cycle
            wpmin = wp0
            q2m0 = real(q2m(1,1,1)); wp0 = aimag(q2m(1,1,1))
            if (wp0 /= wpmin) then
               call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny)
            endif
! initiate momentum diagnostic
            if (ntm > 0) then
               call initmomt2(part,np,pxe,pye,pze)
               if (movion==1) call initmomt2(parti,npi,pxi,pyi,pzi)
            endif
! extend run
            if (nustrt==0) then
               itime0 = itime + itime0
               t0 = dt*real(itime0)
               itime = 0
               ntime = itime + itime0
               if (iur1 >= 0) close (unit=iur1)
               if (iur2 >= 0) close (unit=iur2)
               call restart_open(1,ntr,idrun,iur1,iur2,iuer)
               exit
            endif
! read diagnostics
            call restart_dread(it,itime,itw,wt,iud,ndrec,fdname,iup,npre&
     &c,fpname,irc,iuer,iua,narec,faname)
            if (irc /= 0) cycle
            ntime = itime + itime0
            t0 = dt*real(itime0)
            rewind it
            exit
         enddo
         if (irc /= 0) then
            write (iuer,*) 'Restart Error, irc = ', irc
            go to 3000
         endif
      endif
! record time
      call wtimer(time,ltime)
      write (iuot,*) 'initialization wall clock time = ', time, 'sec'
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
      write (iuot,991) ntime
! initialize current density to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit current for electrons
      call djpostg(part,cu,np,qme,zero,ci,tdjpost,nx,ny,ipbc,relativity,&
     &inorder,djopt)
! deposit electron charge density
      call dpostg(part,qe,np,nx,ny,qme,tdpost,inorder,dopt)
! deposit current and charge density for ions
      if (movion==1) then
         call djpostg(parti,cu,npi,qmi,zero,ci,tdjposti,nx,ny,ipbc,relat&
     &ivity,inorder,djopt)
         call dpostg(parti,qi,npi,nx,ny,qmi,tdposti,inorder,dopt)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti)
            movion = 0
         endif
      endif
! add guard cells for current
      call acguard(cu,nx,ny,inorder)
! add electron and ion densities
      call addqei(qe,qi,qbme,qbmi,wpmax,wpmin,nx,ny,inorder)
      wp0 = 0.5*(wpmax + wpmin)
! recalculate form factors
      if ((wp0 > 1.15*q2m0) .or. (wp0 < 0.85*q2m0)) then
         q2m0 = wp0
         wp0 = affp*wp0
         q2m = cmplx(q2m0,wp0)
         call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny)
         write (iuer,*) ntime, 'new shift constants,q2m0,wp0=', q2m
      endif
! ion density diagnostic
      if ((ntd > 0) .or. (ndd > 0)) then
         it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
         jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
         if ((it==0) .or. (jt==0)) then
            qt = qi
! transform ion density to fourier space
            isign = -1
            call fft(qt,isign,mixup,sct,tfft,indx,indy,inorder)
! calculate smoothing in fourier space
            call spois(qt,sfield,ffc,nx,ny,inorder)
! store selected fourier modes
            if (it==0) then
               call gtmodes(sfield,dent,nx,ny,modesxd,modesyd,inorder)
! write diagnostic output
               if (ndrec==0) then
                  ndrec = -1; iud = get_funit(iud)
                  call bfopen(dent,modesxd,modesy2d,iud,ndrec,trim(fdnam&
     &e))
               endif
               call writebf(dent,modesxd,modesy2d,iud,ndrec,order=LINEAR&
     &)
            endif
! transform ion density to real space
            if (jt==0) then
               isign = 1
               call fft(sfield,isign,mixup,sct,tfft,indx,indy,inorder)
               call cguard(sfield,nx,ny,inorder)
! display ion density
               call displays(sfield,' I DENSITY',ntime,999,2,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! velocity diagnostic
      if ((ntv > 0) .or. (ndv > 0)) then
         it = -1; if (ntv > 0) it = ntime - ntv*(ntime/ntv)
         jt = -1; if (ndv > 0) jt = ntime - ndv*(ntime/ndv)
         if ((it==0) .or. (jt==0)) then
! calculate electron distribution function and moments
            call vdist(part,fv,fvm,np,nmv)
! calculate electron field momentum distribution function and moments
            call adist(part,fa,fam,np,nma)
! print out velocity and field momentum moments
            if (it==0) then
               write (iuv,*) it, fvm(1,:), fvm(2,:)
               write (10,*) it, fam(1,:), fam(2,:)
            endif
! display velocity and field momentum distributions
            if (jt==0) then
               call displayfv(fv,fvm,' ELECTRON',ntime,nmv,2,irc)
               if (irc==1) go to 2000
               call displayfa(fa,fam,' ELECTRON',ntime,nma,2,irc)
               if (irc==1) go to 2000
            endif
            if (movion==1) then
! calculate ion distribution function and moments
               call vdist(parti,fvi,fvmi,npi,nmv)
! print out velocity moments
               if (it==0) then
                  write (iuvi,*) it, fvm(1,:), fvm(2,:)
               endif
! display velocity distribution
               if (jt==0) then
                  call displayfv(fvi,fvmi,' ION',ntime,nmv,2,irc)
                  if (irc==1) go to 2000
               endif
            endif
         endif
      endif
! phase space diagnostic
      fname = ' ELECTRON PHASE SPACE'
      call phasediag(part,nts,nds,nx,ny,np,npxy,ntime,fname,irc)
      if (movion==1) then
         fname = ' ION PHASE SPACE'
         call phasediag(parti,nts,nds,nx,ny,npi,npxyi,ntime,fname,irc)
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,isign,mixup,sct,tfft,indx,indy,inorder)
! potential diagnostic
      if ((ntp > 0) .or. (ndp > 0)) then
         it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
         jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
         if ((it==0) .or. (jt==0)) then
! calculate potential in fourier space
            call pois(qe,sfield,ffc,we,nx,ny,inorder)
! store selected fourier modes
            if (it==0) then
               call gtmodes(sfield,pott,nx,ny,modesxp,modesyp,inorder)
! write diagnostic output
               if (nprec==0) then
                  nprec = -1; iup = get_funit(iup)
                  call bfopen(pott,modesxp,modesy2p,iup,nprec,trim(fpnam&
     &e))
               endif
               call writebf(pott,modesxp,modesy2p,iup,nprec,order=LINEAR&
     &)
            endif
! transform potential to real space
            if (jt==0) then
               isign = 1
               call fft(sfield,isign,mixup,sct,tfft,indx,indy,inorder)
               call cguard(sfield,nx,ny,inorder)
! display potential
               call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! calculate longitudinal electric field in fourier space
      call pois3(qe,fxyze,ffc,we,tfield,nx,ny,inorder)
! transform longitudinal electric field to real space
      isign = 1
      call fft(fxyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(fxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(fxyze,nx,ny,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! calculate smoothed derivatives of vector potential in fourier space
      call dapois(cu,daxyze,ffc,ci,wm,tfield,nx,ny,inorder)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poynt(qe,cu,ffc,ci,sx,sy,sz,nx,ny,inorder)
         endif
      endif
! transform derivative of vector potential to real space
      isign = 1
      call fftn(daxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(daxyze,nx,ny,inorder)
! calculate smoothed vector potential in fourier space
      call sapois(cu,axyze,ffc,ci,tfield,nx,ny,inorder)
      if (ndc > 0) cus = q2m0*axyze
! transform vector potential to real space
      isign = 1
      call fft(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(axyze,nx,ny,inorder)
!
! main iteration loop
      do k = 1, ndc
! initialize current density to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit current
      call dhjpostg(part,fxyze,axyze,daxyze,cu,np,qme,qbme,dt,ci,tdhjpos&
     &t,relativity,inorder,djopt)
      if (movion==1) then
         call dhjpostg(parti,fxyze,axyze,daxyze,cu,npi,qmi,qbmi,dt,ci,td&
     &hjposti,relativity,inorder,djopt)
      endif
! add guard cells for current
      call acguard(cu,nx,ny,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! modify current for stability
      cus = cus + cu
! calculate smoothed derivatives of vector potential in fourier space
      call dapois(cu,daxyze,ffc,ci,wm,tfield,nx,ny,inorder)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poynt(qe,cu,ffc,ci,sx,sy,sz,nx,ny,inorder)
         endif
      endif
! transform derivative of vector potential to real space
      isign = 1
      call fftn(daxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(daxyze,nx,ny,inorder)
! calculate smoothed vector potential in fourier space
      call sapois(cus,axyze,ffe,ci,tfield,nx,ny,inorder)
      if (k.lt.ndc) cus = q2m0*axyze
! transform vector potential to real space
      isign = 1
      call fft(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(axyze,nx,ny,inorder)
      enddo
!
! vector potential diagnostic
      if ((nta > 0) .or. (nda > 0)) then
         it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
         jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
         if ((it==0) .or. (jt==0)) then
! calculate vector potential in fourier space
            call apois(cu,vfield,ffc,ci,wm,nx,ny,inorder)
! store selected fourier modes
            if (it==0) then
               call gtmodes(vfield,vpott,nx,ny,modesxa,modesya,inorder)
! write diagnostic output
               if (narec==0) then
                  narec = -1; iua = get_funit(iua)
                  call bfopen(vpott,modesxa,modesy2a,iua,narec,trim(fana&
     &me))
               endif
               call writebf(vpott,modesxa,modesy2a,iua,narec,order=LINEA&
     &R)
            endif
! transform vector potential to real space
            if (jt==0) then
               isign = 1
               call fft(vfield,isign,mixup,sct,tfft,indx,indy,inorder)
               call cguard(vfield,nx,ny,inorder)
! display absolute value of vector potential
               call displayv(vfield,' VECTOR POTENTIAL',ntime,999,1,ndst&
     &yle,nx,ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxyze(1,:,:) = fxyze(1,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! advance position and momentum at time t + dt/2
      call hpushg(part,fxyze,axyze,daxyze,np,qbme,dt,ci,wke,tpush,nx,ny,&
     &ipbc,relativity,inorder,popt)
      if (movion==1) then
         wki = 0.
         call hpushg(parti,fxyze,axyze,daxyze,npi,qbmi,dt,ci,wki,tpushi,&
     &nx,ny,ipbc,relativity,inorder,popt)
         wki = wki*rmass
      endif
! initialize current density to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit current for electrons with non-relativistic algorithm
      call djpost(part,cu,np,qme,zero,tdjpost,nx,ny,ipbc,inorder,djopt)
! deposit current for ions with non-relativistic algorithm
      if (movion==1) then
         call djpost(parti,cu,npi,qmi,zero,tdjposti,nx,ny,ipbc,inorder,d&
     &jopt)
      endif
! add guard cells for current
      call acguard(cu,nx,ny,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! calculate smoothed vector potential in fourier space
      call sapois(cu,axyze,ffc,ci,tfield,nx,ny,inorder)
      cus = q2m0*axyze
! transform vector potential to real space
      isign = 1
      call fft(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(axyze,nx,ny,inorder)
!
! secondary iteration loop
      do k = 1, ndc2
! advance position and velocity at time t + dt/2
      call hpushxhg(part,axyze,np,qbme,dt,ci,tpush,nx,ny,ipbc,relativity&
     &,inorder,popt)
      if (movion==1) then
         call hpushxhg(parti,axyze,npi,qbmi,dt,ci,tpushi,nx,ny,ipbc,rela&
     &tivity,inorder,popt)
      endif
! initialize current density to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit current for electrons with non-relativistic algorithm
      call djpost(part,cu,np,qme,zero,tdjpost,nx,ny,ipbc,inorder,djopt)
! deposit current for ions with non-relativistic algorithm
      if (movion==1) then
         call djpost(parti,cu,npi,qmi,zero,tdjposti,nx,ny,ipbc,inorder,d&
     &jopt)
      endif
! add guard cells for current
      call acguard(cu,nx,ny,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! modify current for stability
      cus = cus + cu
! calculate smoothed vector potential in fourier space
      call sapois(cus,axyze,ffe,ci,tfield,nx,ny,inorder)
      if (k.lt.ndc) cus = q2m0*axyze
! transform vector potential to real space
      isign = 1
      call fft(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(axyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(axyze,nx,ny,inorder)
      enddo
!
! advance position and velocity at time t + dt
      call hpushxg(part,axyze,np,qbme,dt,ci,tpush,nx,ny,ipbc,relativity,&
     &inorder,popt)
      if (movion==1) then
         call hpushxg(parti,axyze,npi,qbmi,dt,ci,tpushi,nx,ny,ipbc,relat&
     &ivity,inorder,popt)
      endif
! calculate electron and field momentum
      if (ntm > 0) then
         it = ntime/ntm
         it = ntime - ntm*it + 1
         if (it > 1) it = it - ntm
         if (it >= 0) then
            call premoment2(part,ntime,np,ium,pxe,pye,pze,sx,sy,sz,wx,wy&
     &,wz,ndim,nprint=it)
!           call precmoment2(part,ntime,np,ium,wx,wy,wz,ndim)
! calculate ion momentum
            if (movion==0) then
               if (it==1) then
! exyze not defined
!                 call imoment(qi,exyze,ium,pxi,pyi,pzi,dt,wx,wy,wz,nx,n&
!    &y,inorder)
                  continue
               endif
            else if (movion==1) then
               call primoment2(parti,npi,ium,rmass,pxi,pyi,pzi,wx,wy,wz,&
     &ndim,nprint=it)
!              call pricmoment2(parti,npi,ium,rmass,wx,wy,wz,ndim)
            endif
            if (it==1) then
! print total momentum
               write (ium,996) wx, wy, wz
            endif
         endif
      endif
! sort electrons
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
!           call sortp(part,pt,ip,np,npic,tsort,inorder)
            call sortp(part,part2,np,npic,tsort,inorder)
         endif
      endif
! sort ions
      if ((movion==1).and.(sortimi > 0)) then
         if (mod(ntime,sortimi)==0) then
            call sortp(parti,pt,ip,npi,npic,tsorti,inorder)
!           call sortp(parti,parti2,npi,npic,tsorti,inorder)
         endif
      endif
! energy diagnostic
      if ((ntw > 0) .or. (ndw > 0)) then
         it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
         jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
         if ((it==0) .or. (jt==0)) then
            wef = we + wm
            wtot = wef + wke + wki
            if (it==0) then
               write (iuot,992) wef, wke, wtot
               write (iuot,993) we, wf, wm
            endif
            if (jt==0) then
               itw = itw + 1
               wt(itw,:) = (/wef,wke,wki,wtot,we,wf,wm/)
            endif
         endif
      endif
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
! write primary file
            call restart_bwrite(iur1,itime,itime0,np,part,movion,npi,par&
     &ti,qi,q2m)
            call restart_dwrite(iur1,itime,itw,wt,ndrec,fdname,nprec,fpn&
     &ame,narec,faname)
! write secondary file
            call restart_bwrite(iur2,itime,itime0,np,part,movion,npi,par&
     &ti,qi,q2m)
            call restart_dwrite(iur2,itime,itw,wt,ndrec,fdname,nprec,fpn&
     &ame,narec,faname)
         endif
      endif
      call wtimer(tloop,ltime)
      time = time + tloop
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! energy diagnostic
      if (ndw > 0) then
         ts = t0 + dt*real(ntw)*(itime0-(itime0/ntw)*ntw)
         call displayw(wt,ts,dt*real(ntw),itw,irc)
! check error return code
         if (irc==1) go to 3000
      endif
! accumulate timings
      write (iuot,*) 'hamiltonian darwin code hdbeps2'
      write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost + tdjpost  + tdhjpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'electron current deposit time = ', tdjpost, 'sec'
      write (iuot,*) 'electron current iteration deposit time = ', tdhjp&
     &ost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      if (movion==1) then
         totpushi = tpushi + tdposti + tdjposti  + tdhjposti
         write (iuot,*) 'ion push time = ', tpushi, 'sec'
         write (iuot,*) 'ion charge deposit time = ', tdposti, 'sec'
         write (iuot,*) 'ion current deposit time = ', tdjposti, 'sec'
         write (iuot,*) 'ion current iteration deposit time = ', tdhjpos&
     &ti, 'sec'
         write (iuot,*) 'total ion push time = ', totpushi, 'sec'
         write (iuot,*) 'ion sort time = ', tsorti
         totpushi = totpushi + tsorti
         write (iuot,*) 'total ion time = ', totpushi, 'sec'
      endif
      write (iuot,*) 'total fft time=', tfft, 'sec'
      write (iuot,*) 'field time=', tfield
      time = time - (totpush + totpushi + tfft + tfield)
      write (iuot,*) 'other time=', time, 'sec'
! write final diagnostic metafile
      close(unit=iudm)
      fname = 'diag2.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! ion density diagnostics
      if (ntd > 0) then
         ndrec = ndrec - 1
         ceng = zero
         write (iudm,den2d,iostat=irc)
         if (irc /= 0) write (iuer,*) 'den2d namelist not written'
      endif
! potential diagnostics
      if (ntp > 0) then
         nprec = nprec - 1
         ceng = affp
         write (iudm,pot2d,iostat=irc)
         if (irc /= 0) write (iuer,*) 'pot2d namelist not written'
      endif
! vector potential diagnostics
      if (nta > 0) then
         narec = narec - 1
         ceng = affp
         write (iudm,vpot2d,iostat=irc)
         if (irc /= 0) write (iuer,*) 'vpot2d namelist not written'
      endif
! write out input file
      write (iudm,input2,iostat=irc)
      if (irc /= 0) write (iuer,*) 'input2 namelist not written'
! done
      write (iuot,*) '* * * q.e.d. * * *'
      close(unit=iudm)
      close(unit=iuot)
      close(unit=iuer)
! close graphics device
 3000 call GRCLOSE
      call MP_END
      stop
      end program hdbeps2
