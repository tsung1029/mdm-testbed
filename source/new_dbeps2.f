!-----------------------------------------------------------------------
! * * * periodic 2d darwin particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with periodic
! electromagnetic forces obtained by solving non-radiative form of
! maxwell's equation with fast fourier transforms
! using algorithm similar to that described in
! J. Busnardo-Neto, P. L. Pritchett, A. T. Lin, and J. M. Dawson,
! J. Computational Phys. 23, 300 (1977).
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 2006, regents of the university of california
! update: march 2, 2012
      program dbeps2
      use init2d
      use empush2d
      use field2d
      use diag2d
      use emsimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: nmv = 40, ipbc = 1
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, iua = 15, iuj = 25, ium = 21
      integer :: iuer = 2
      integer :: npxy, npxyb, np, npxyi, npxybi, npi, npxy1
      integer :: nx, ny, nxh, nyh, nxe, nye
      integer :: nxye, nxeh, nxyh, nxhy, nx1, ny1, ntmax, idimp
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: j, k, it, itw, iur1, iur2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, time = 0.0, tloop = 0.0, ts = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0 
      real :: tfft = 0.0, tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: tdcjpost = 0.0, tdcjposti = 0.0
      real :: qbme, qbmi, affp, qi0, omt, q2m0, wp0, wpmax, wpmin
      real :: etx, we, wf, wm, wke
      real :: pxe = 0.0, pye = 0.0, pze = 0.0
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: sx, sy, sz, wx, wy, wz
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      real, dimension(:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:), pointer :: qe, qi
      real, dimension(:,:,:), pointer :: fxyze, exyze, bxyze
      real, dimension(:,:,:), pointer :: cu, cus, amu
      complex, dimension(:,:,:), pointer :: q2m
      complex, dimension(:,:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:,:), pointer :: qt, sfield
      complex, dimension(:,:), pointer :: dent, pott
      real, dimension(:,:,:), pointer :: cut, vfield
      complex, dimension(:,:,:), pointer :: vpott, vcurt
      real, dimension(:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
      character(len=10) :: cdrun
      character(len=32) :: fname
      integer, external :: NDIAN, NDPREC, IDPREC
  991 format (' T = ',i7)
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
      idcode = 3
      psolve = 1
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      iuot = get_funit(iuot)
      fname = 'output2.'//cdrun
      open(unit=iuot,file=trim(fname),form='formatted',status='replace')
! np = total number of electrons in simulation
      npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb
! npi = total number of ions in simulation
      npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi
      nx = 2**indx; ny = 2**indy; nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 4; nye = ny + 3
!     ax = .866025; ay = .866025
      if (inorder==LINEAR) then
!        ax = .912871; ay = .912871
         nxe = nx + 2; nye = ny + 1
      else if (inorder==CUBIC) then
         nxe = nx + 6; nye = ny + 5
!        ax = .816497; ay = .816497
      endif
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
      call MP_SETSTACK(262144)
      nxye = nxe*nye; nxeh = nxe/2
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
! idimp = dimension of phase space = 4 or 5
      idimp = 2 + ndim
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
      allocate(part(idimp,np))
! in real space, qe(j,k) = charge density at grid point (j,k)
! in real space, qi(j,k) = ion charge density at grid point (j,k)
      allocate(qe(nxe,nye),qi(nxe,nye))
! fxyze(i,j,k) = i component of longitudinal electric force/charge
! exyze(i,j,k) = i component of transverse electric force/charge
      allocate(fxyze(ndim,nxe,nye),exyze(ndim,nxe,nye))
! bxyze(i,j,k) = i component of convolution of magnetic force/charge
      allocate(bxyze(2*ndim-3,nxe,nye))
! cu(i,j,k) = i component of current density
! cus(i,j,k) = i component of acceleration density
      allocate(cu(ndim,nxe,nye),cus(ndim,nxe,nye))
! amu(i,j,k) = i component of momentum flux
      allocate(amu(2*ndim-2,nxe,nye))
! store shift constants
      allocate(q2m(1,1,1))
! ffc, ffe = form factor arrays for poisson solvers
      allocate(ffc(nxeh,nyh),ffe(nxeh,nyh))
! mixup, sct = arrays for fft
      allocate(mixup(nxhy),sct(nxyh))
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
      omt = sqrt(omx*omx + omy*omy + omz*omz)
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
! set initial time
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! set default diagnostic file names
      if (ntd > 0) fdname = 'denk2.'//cdrun
      if (ntp > 0) fpname = 'potk2.'//cdrun
      if (nta > 0) faname = 'vpotk2.'//cdrun
      if (ntj > 0) fjname = 'vcurk2.'//cdrun
! energy diagnostics
      if (ndw > 0) then
         allocate(wt((nloop-1)/ndw-(itime0/ndw)+1,7))
         itw = 0
      endif
! open restart files
      if (nustrt==0) then
         call restart_open(nustrt,ntr,idrun0,iur1,iur2,iuer)
      else
         call restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
      endif
      if ((iur1 < 0) .and. (iur2 < 0)) go to 3000
! initialize external magnetic field
      bxyze = 0.0
      if (omt > 0) then
         call baddext(bxyze,omx,omy,omz,nx,ny,inorder)
         call cguard(bxyze,nx,ny,inorder)
      endif
      fxyze = 0.0; cu = 0.0
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
         if (npxy > 0) then
            call fdistr(part,1,npxy,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,npx,npy,nx,ny,ipbc,ndprof,nsrand)
            call vdistr(part,1,npxy,vtx,vty,vtz,vx0,vy0,vz0,ndim)
         endif
! beam electrons
         if (npxyb > 0) then
            call fdistr(part,npxy1,npxyb,ampdx,scaledx,shiftdx,ampdy,   &
     &scaledy,shiftdy,npxb,npyb,nx,ny,ipbc,ndprof,nsrand)
            call vdistr(part,npxy1,npxyb,vtdx,vtdy,vtdz,vdx,vdy,vdz,ndim&
     &)
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
            call dpostg(part,qi,np,nx,ny,-qme,tdpost,inorder,dopt)
! debug
!           call sguard(qi,qi0,nx,ny,inorder)
         endif
! fix guiding centers for electrons
         if (omt > 0) then
            if (relativity==1) then
               call distr(part,bxyze,np,qbme,ci,nx,ny,ipbc,inorder)
            else
               call distr(part,bxyze,np,qbme,nx,ny,ipbc,inorder)
            endif
         endif
! calculate initial electron momentum
         if (ntm > 0) call initmomt2(part,np,pxe,pye,pze,ndim)
! initialize ions
         if (movion==1) then
            npxy1 = npxyi + 1
! background ions
            if (npxyi > 0) then
               call fdistr(parti,1,npxyi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,npxi,npyi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,1,npxyi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0, &
     &ndim)
            endif
! beam ions
            if (npxybi > 0) then
               call fdistr(parti,npxy1,npxybi,ampdxi,scaledxi,shiftdxi, &
     &ampdyi,scaledyi,shiftdyi,npxbi,npybi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,npxy1,npxybi,vtdxi,vtdyi,vtdzi,vdxi,   &
     &vdyi,vdzi,ndim)
            endif
! fix guiding centers for ions
            if (omt > 0) then
               if (relativity==1) then
                  call distr(parti,bxyze,npi,qbmi,ci,nx,ny,ipbc,inorder)
               else
                  call distr(parti,bxyze,npi,qbmi,nx,ny,ipbc,inorder)
               endif
            endif
! calculate initial ion momentum
            if (ntm > 0) call initmomt2(parti,npi,pxi,pyi,pzi,ndim) 
         endif
! freeze the ions now
         if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density
            call dpostg(parti,qi,npi,nx,ny,qmi,tdposti,inorder,dopt)
! delete ions
            deallocate(parti)
            movion = 0
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
               call initmomt2(part,np,pxe,pye,pze,ndim)
               if (movion==1) call initmomt2(parti,npi,pxi,pyi,pzi,ndim)
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
            call restart_dread(it,itime,itw,wt,iud,ndrec,fdname,iup,    &
     &nprec,fpname,irc,iuer,iua,narec,faname,iuj,njrec,fjname)
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
!
! sorting arrays
      allocate(pt(max(np,npi)),ip(max(np,npi)),npic(ny1))
      if (sortime > 0) then
         allocate(part2(idimp,np))
      else
         allocate(part2(0,0))
      endif
! initialize diagnostics
! open initial diagnostic metafile
      iudm = get_funit(iudm)
      fname = 'diag2.init.'//cdrun
      open(unit=iudm,file=trim(fname),form='formatted',status='replace')
! ion density or potential diagnostics
      if ((ntp > 0) .or. (ndp > 0) .or. (ntd > 0) .or. (ndd > 0)) then
         allocate(sfield(nxe,nye))
      endif
! ion density diagnostics
      call initmodediag(dent,ntd,nxh,nyh,modesxd,modesyd,iud,ndrec,     &
     &fdname)
      if (ntd > 0) then
         allocate(qt(nxe,nye))
         ceng = zero
         write (iudm,den2d,iostat=irc)
      endif
! velocity diagnostics
      fname = 'fv2.'//cdrun
      call initveldiag(fv,fvm,vtx,vty,zero,ntv,ndv,nmv,ndim,iuv,fname)
      if (movion==1) then
         fname = 'fvi2.'//cdrun
         call initveldiag(fvi,fvmi,vtxi,vtyi,zero,ntv,ndv,nmv,2,iuvi,   &
     &fname)
      endif
! potential diagnostics
      call initmodediag(pott,ntp,nxh,nyh,modesxp,modesyp,iup,nprec,     &
     &fpname)
      if (ntp > 0) then
         ceng = affp
         write (iudm,pot2d,iostat=irc)
      endif
! vector potential or ion current diagnostics
      if ((nta > 0) .or. (nda > 0) .or. (ntj > 0) .or. (ndj > 0)) then
         allocate(vfield(ndim,nxe,nye))
         vfield = 0.0
      endif
! vector potential diagnostics
      call initvmodediag(vpott,nta,nxh,nyh,ndim,modesxa,modesya,iua,    &
     &narec,faname)
      if (nta > 0) then
         ceng = affp
         write (iudm,vpot2d,iostat=irc)
      endif
! ion current diagnostics
      call initvmodediag(vcurt,ntj,nxh,nyh,ndim,modesxj,modesyj,iuj,    &
     &njrec,fjname)
      if (ntj > 0) then
         allocate(cut(ndim,nxe,nye))
         ceng = zero
         write (iudm,vcur2d,iostat=irc)
      endif
! momentum diagnostics
      fname = 'momentum2.'//cdrun
      if (ntm > 0) then
         ium = get_funit(ium )
         open(unit=ium,file=trim(fname),form='formatted',status='unknown&
     &')
      endif
! write out input file
      write (iudm,input2,iostat=irc)
! close diagnostic metafile
      close(unit=iudm)
! record time
      call wtimer(time,ltime)
      write (iuot,*) 'initialization wall clock time = ', time, 'sec'
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
      write (iuot,991) ntime
! initialize current density and momentum flux to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
      call sguard(amu,zero,zero,zero,zero,nx,ny,inorder)
! deposit current and momentum flux for electrons
      call djpostg(part,cu,np,qme,zero,ci,tdjpost,nx,ny,ipbc,relativity,&
     &inorder,djopt)
      call dmjpostg(part,amu,np,qme,ci,tdcjpost,relativity,inorder,djopt&
     &)
! deposit electron charge density
      call dpostg(part,qe,np,nx,ny,qme,tdpost,inorder,dopt)
! save electron current for ion current diagnostic
      if (ndc==0) call vpotdiagprep(cu,vfield,ntj,ndj,ntime)
! deposit current, momentum flux, and charge density for ions
      if (movion==1) then
         call djpostg(parti,cu,npi,qmi,zero,ci,tdjposti,nx,ny,ipbc,     &
     &relativity,inorder,djopt)
         call dmjpostg(parti,amu,npi,qmi,ci,tdcjposti,relativity,inorder&
     &,djopt)
         call dpostg(parti,qi,npi,nx,ny,qmi,tdposti,inorder,dopt)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti)
            movion = 0
         endif
      endif
! ion current diagnostic
      if ((ndc==0)) then
         call vcurdiag(cut,cu,vfield,vcurt,ffc,mixup,sct,tfft,ntj,ndj,nx&
     &,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,ndstyle,1,irc,      &
     &inorder)
         if (irc==1) go to 2000
      endif
! add guard cells for current and momentum flux
      call acguard(cu,nx,ny,inorder)
      call amcguard(amu,nx,ny,inorder)
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
      call dendiag(qt,qi,sfield,dent,ffc,mixup,sct,tfft,ntd,ndd,nx,ny,  &
     &modesxd,modesyd,iud,ndrec,indx,indy,ntime,ndstyle,1,irc,inorder)
      if (irc==1) go to 2000
! velocity diagnostic
      call veldiag(part,fv,fvm,ntv,ndv,np,nmv,iuv,ntime,' ELECTRON',irc)
      if (irc==1) go to 2000
      if (movion==1) then
         call veldiag(parti,fvi,fvmi,ntv,ndv,npi,nmv,iuvi,ntime,' ION', &
     &irc)
         if (irc==1) go to 2000
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
      call potdiag(qe,sfield,pott,ffc,mixup,sct,tfft,ntp,ndp,nx,ny,     &
     &modesxp,modesyp,iup,nprec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! calculate longitudinal electric force in fourier space
      call pois3(qe,fxyze,ffc,we,tfield,nx,ny,inorder)
! transform longitudinal electric force to real space
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
! calculate magnetic field in fourier space
      call bpois(cu,bxyze,ffc,ci,wm,tfield,nx,ny,inorder)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poynt(qe,cu,ffc,ci,sx,sy,sz,nx,ny,inorder)
         endif
      endif
! vector potential diagnostic
      if (ndc==0) call vpotdiagprep(cu,vfield,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(bxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      if (omt > 0) call baddext(bxyze,omx,omy,omz,nx,ny,inorder)
      call cguard(bxyze,nx,ny,inorder)
! transform momentum flux to fourier space
      isign = -1
      call fftn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! debug
!     call iwfft2rn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of time derivative of current
      call dcuperp(cus,amu,nx,ny,inorder)
! calculate convective part of transverse electric field
      call epois(cus,cu,ffe,ci,wf,tfield,nx,ny,inorder)
! transform transverse electric field to real space
      isign = 1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(cu,nx,ny,inorder)
! add longitudinal and transverse fields
!     exyze = cu + fxyze
      call addfields(exyze,cu,fxyze)
!
! inner iteration loop
      do k = 1, ndc
! initialize current, acceleration density and momentum flux
      call sguard(cus,cu,q2m0,nx,ny,inorder)
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
      call sguard(amu,zero,zero,zero,zero,nx,ny,inorder)
! deposit electron current and acceleration density and momentum flux
      call dcjpostg(part,exyze,bxyze,cu,cus,amu,np,qme,qbme,dt,ci,      &
     &tdcjpost,relativity,inorder,popt)
! save electron current for ion current diagnostic
      if (k==ndc) call vpotdiagprep(cu,vfield,ntj,ndj,ntime)
! deposit ion current and acceleration density and momentum flux
      if (movion==1) then
         call dcjpostg(parti,exyze,bxyze,cu,cus,amu,npi,qmi,qbmi,dt,ci, &
     &tdcjposti,relativity,inorder,popt)
      endif
! ion current diagnostic
      if (k==ndc) then
         call vcurdiag(cut,cu,vfield,vcurt,ffc,mixup,sct,tfft,ntj,ndj,nx&
     &,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,ndstyle,1,irc,      &
     &inorder)
         if (irc==1) go to 2000
      endif
! add guard cells for current, acceleration density, and momentum flux
      call acguard(cu,nx,ny,inorder)
      call acguard(cus,nx,ny,inorder)
      call amcguard(amu,nx,ny,inorder)
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! calculate magnetic field in fourier space
      call bpois(cu,bxyze,ffc,ci,wm,tfield,nx,ny,inorder)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poynt(qe,cu,ffc,ci,sx,sy,sz,nx,ny,inorder)
         endif
      endif
! save current for vector potential diagnostic
      if (k==ndc) call vpotdiagprep(cu,vfield,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      if (omt > 0) call baddext(bxyze,omx,omy,omz,nx,ny,inorder)
      call cguard(bxyze,nx,ny,inorder)
! transform acceleration density and momentum flux to fourier space
      isign = -1
      call fft(cus,isign,mixup,sct,tfft,indx,indy,inorder)
      call fftn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! debug
!     call iwfft2rn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of time derivative of current
      call adcuperp(cus,amu,nx,ny,inorder)
! calculate transverse electric field
      call epois(cus,cu,ffe,ci,wf,tfield,nx,ny,inorder)
! transform transverse electric field to real space
      isign = 1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(cu,nx,ny,inorder)
! add longitudinal and transverse fields
!     exyze = cu + fxyze
      call addfields(exyze,cu,fxyze)
      enddo
!
! vector potential diagnostic
      call vpotdiag(vfield,cus,vpott,ffc,mixup,sct,ci,tfft,nta,nda,nx,ny&
     &,modesxa,modesya,iua,narec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        exyze(1,:,:) = exyze(1,:,:) + etx
!     endif
! push particles
      wke = 0.
      call push3g(part,exyze,bxyze,np,qbme,dt,dt,ci,wke,tpush,nx,ny,ipbc&
     &,relativity,inorder,popt)
! debug: zero force
!     call push3zfg(part,np,dt,ci,wke,tpush,nx,ny,ipbc,ndim,relativity)
      if (movion==1) then
         wki = 0.
         call push3g(parti,exyze,bxyze,npi,qbmi,dt,dt,ci,wki,tpushi,nx, &
     &ny,ipbc,relativity,inorder,popt)
! debug: zero force
!        call push3zfg(parti,npi,dt,ci,wki,tpushi,nx,ny,ipbc,ndim,      &
!    &relativity)
         wki = wki*rmass
      endif
! calculate electron momentum
      call emomtdiag(part,pxe,pye,pze,sx,sy,sz,wx,wy,wz,ntm,np,ium,ntime&
     &,ndim)
! calculate ion and total momentum
      call imomtdiag(parti,qi,exyze,dt,rmass,pxi,pyi,pzi,wx,wy,wz,ntm,  &
     &movion,npi,ium,nx,ny,ntime,ndim,inorder)
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
      call dmenergy(wt,we,wf,wm,wke,wki,ntw,ndw,itw,iuot,ntime)
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
! write primary file
            call restart_bwrite(iur1,itime,itime0,np,part,movion,npi,   &
     &parti,qi,q2m)
            call restart_dwrite(iur1,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname,narec,faname,njrec,fjname)
! write secondary file
            call restart_bwrite(iur2,itime,itime0,np,part,movion,npi,   &
     &parti,qi,q2m)
            call restart_dwrite(iur2,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname,narec,faname,njrec,fjname)
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
      write (iuot,*) 'darwin code dbeps2'
      write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost + tdjpost + tdcjpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'electron current deposit time = ', tdjpost, 'sec'
      write (iuot,*) 'electron current derivative deposit time = ',     &
     &tdcjpost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      if (movion==1) then
         totpushi = tpushi + tdposti + tdjposti + tdcjposti
         write (iuot,*) 'ion push time = ', tpushi, 'sec'
         write (iuot,*) 'ion charge deposit time = ', tdposti, 'sec'
         write (iuot,*) 'ion current deposit time = ', tdjposti, 'sec'
         write (iuot,*) 'ion current derivative deposit time = ', tdcjpo&
     &sti, 'sec'
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
! ion current diagnostics
      if (ntj > 0) then
         njrec = njrec - 1
         ceng = zero
         write (iudm,vcur2d,iostat=irc)
         if (irc /= 0) write (iuer,*) 'vcur2d namelist not written'
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
      end program dbeps2
