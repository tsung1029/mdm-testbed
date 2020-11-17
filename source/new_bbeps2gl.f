!-----------------------------------------------------------------------
! * * * periodic 2d electromagnetic particle simulation kernel code * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with periodic
! electromagnetic forces obtained by solving maxwells's equation with
! gridless implementation.
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G3
! copyright 1999, regents of the university of california
! update: march 1, 2012
      program bbeps2gl
      use init2d
      use empush2d
      use field2d
      use diag2d
      use emsimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! idimp = dimension of phase space = 5
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      integer :: idimp = 5, nmv = 40, ipbc = 1
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, iua = 15, iue = 25, ium = 21
      integer :: iuer = 2
      integer :: npxy, npxyb, np, npxyi, npxybi, npi, npxy1
      integer :: nx, ny, nxh, nyh, nxe, nye
      integer :: nxye, nxeh, nxyh, nxhy, nx1, ny1, ntmax
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: j, it, itw, iur1, iur2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, time = 0.0, tloop = 0.0, ts = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0 
      real :: tfft = 0.0, tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, dth, qi0, q2m0, wp0
      real :: etx, we, wf, wm, wke
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: sx, sy, sz, pxe, pye, pze, wx, wy, wz
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      real, dimension(:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:), pointer :: qe, qi
      real, dimension(:,:,:), pointer :: cu, amu, fxyze, bxyze
      complex, dimension(:,:,:), pointer :: exyz, bxyz
      complex, dimension(:,:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:,:), pointer :: qt, sfield
      complex, dimension(:,:), pointer :: dent, pott
      real, dimension(:,:,:), pointer :: vfield
      complex, dimension(:,:,:), pointer :: vpott, vpotr
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
      idcode = 12
      psolve = 1
      ndim = 3
      inorder = 1
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
      nxe = nx + 2; nye = ny + 1
      ax = 1.0; ay = 1.0
!     ax = 0.0; ay = 0.0
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
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
      allocate(part(idimp,np))
! in real space, qe(j,k) = charge density at grid point (j,k)
! in real space, qi(j,k) = ion charge density at grid point (j,k)
      allocate(qe(nxe,nye),qi(nxe,nye))
! in real space, fxyze(i,j,k) = i component of force/charge at grid (j,k)
! in other words, fxyze are the convolutions of the electric field
! over the particle shape
      allocate(fxyze(3,nxe,nye))
! cu(i,j,k) = i component of current at grid (j,k).
! bxyze(i,j,k) = i component of magnetic field at grid (j,k).  
! bxyze is the convolution of the magnetic field over the particle shape
      allocate(cu(3,nxe,nye),bxyze(3,nxe,nye))
! in fourier space, exyz = transverse electric field
! in fourier space, bxyz = magnetic field
      allocate(exyz(3,nxeh,nye),bxyz(3,nxeh,nye))
! ffc = form factor array for poisson solver
      allocate(ffc(nxeh,nyh))
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
      affp = float(nx*ny)/float(np)
      q2m0 = qbme*qme*real(np)/real(nx*ny)
      dth = 0.0
! debug
!     dth = .5*dt
! end debug
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
      if (nte > 0) fename = 'vpotrk2.'//cdrun
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
! initialize electromagnetic fields
      bxyze = 0.0
      bxyz = cmplx(0.0,0.0)
      exyz = cmplx(0.0,0.0)
      cu = 0.0
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny)
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
            call vdistr(part,1,npxy,vtx,vty,vtz,vx0,vy0,vz0)
         endif
! beam electrons
         if (npxyb > 0) then
            call fdistr(part,npxy1,npxyb,ampdx,scaledx,shiftdx,ampdy,   &
     &scaledy,shiftdy,npxb,npyb,nx,ny,ipbc,ndprof,nsrand)
            call vdistr(part,npxy1,npxyb,vtdx,vtdy,vtdz,vdx,vdy,vdz)
         endif
! fix guiding centers
!        if (relativity==1) then
!           call distr(part,bxyze,np,qbme,ci,nx,ny,ipbc,inorder)
!        else
!           call distr(part,bxyze,np,qbme,nx,ny,ipbc,inorder)
!        endif
! calculate initial electron momentum
         if (ntm > 0) call initmomt2(part,np,pxe,pye,pze)
! initialize ions
         if (movion==1) then
            npxy1 = npxyi + 1
! background ions
            if (npxyi > 0) then
               call fdistr(parti,1,npxyi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,npxi,npyi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,1,npxyi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0)
            endif
! beam ions
            if (npxybi > 0) then
               call fdistr(parti,npxy1,npxybi,ampdxi,scaledxi,shiftdxi, &
     &ampdyi,scaledyi,shiftdyi,npxbi,npybi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,npxy1,npxybi,vtdxi,vtdyi,vtdzi,vdxi,   &
     &vdyi,vdzi)
            endif
! fix guiding centers
!           if (relativity==1) then
!              call distr(parti,bxyze,npi,qbmi,ci,nx,ny,ipbc,inorder)
!           else
!              call distr(parti,bxyze,npi,qbmi,nx,ny,ipbc,inorder)
!           endif
! calculate initial ion momentum
            if (ntm > 0) call initmomt2(parti,npi,pxi,pyi,pzi)  
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
            call sguard(qi,zero,nx,ny,inorder)
            call dpostgl(part,qi,np,-qme,nx,ny,tdpost)
! debug
!           call sguard(qi,qi0,nx,ny,inorder)
! freeze the ions now
         else if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density to background
            call sguard(qi,zero,nx,ny,inorder)
! initialize ion charge density
            call dpostgl(parti,qi,npi,qmi,nx,ny,tdposti)
! delete ions
            deallocate(parti)
            movion = 0
         endif
! retard electron positions to deposit current
!        call retardg(part,np,dth,ci,nx,ny,ipbc,relativity)
! retard ion positions to deposit current
!        if (movion==1) then
!           call retardg(parti,npi,dth,ci,nx,ny,ipbc,relativity)
!        endif
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
     &qi,irc,iuer,exyz,bxyz)
            if (irc /= 0) cycle
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
     &c,fpname,irc,iuer,iua,narec,faname,iue,nerec,fename)
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
      if (ntime > 0) dth = .5*dt
!
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
      call initveldiag(fv,fvm,vtx,vty,zero,ntv,ndv,nmv,2,iuv,fname)
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
! vector potential or electromagnetic diagnostics
      if ((nta > 0) .or. (nda > 0) .or. (nte > 0) .or. (nde > 0)) then
         allocate(vfield(3,nxe,nye))
         vfield = 0.0
      endif
! vector potential diagnostics
      call initvmodediag(vpott,nta,nxh,nyh,ndim,modesxa,modesya,iua,    &
     &narec,faname)
      if (nta > 0) then
         ceng = affp
         write (iudm,vpot2d,iostat=irc)
      endif
! electromagnetic diagnostics
      call initvmodediag(vpotr,nte,nxh,nyh,ndim,modesxe,modesye,iue,    &
     &nerec,fename)
      if (nte > 0) then
         ceng = affp
         write (iudm,em2d,iostat=irc)
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
!
! electromagnetic diagnostic
      if (ntime > 0) call vpotdiagprep(cu,vfield,nte,nde,ntime)
! initialize current density to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit electron current
      call djpostglg(part,cu,np,qme,dth,ci,nx,ny,ipbc,tdjpost,relativity&
     &)
! initialize electron charge density to background
      call sguard(qe,zero,nx,ny,inorder)
! deposit electron charge density
      call dpostgl(part,qe,np,qme,nx,ny,tdpost)
! add ion current and deposit ion charge density
      if (movion==1) then
         call djpostglg(parti,cu,npi,qmi,dth,ci,nx,ny,ipbc,tdjposti,    &
     &relativity)
! initialize ion charge density to background
         call sguard(qi,zero,nx,ny,inorder)
! deposit ion charge density
         call dpostgl(parti,qi,npi,qmi,nx,ny,tdposti)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti)
            movion = 0
         endif
      endif
! add electron and ion densities
      call addqei(qe,qi,nx,ny,inorder)
! ion density diagnostic
      call dendiag(qt,qi,sfield,dent,ffc,mixup,sct,tfft,ntd,ndd,nx,ny,  &
     &modesxd,modesyd,iud,ndrec,indx,indy,ntime,ndstyle,0,irc,inorder)
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
! potential diagnostic
      call potdiag(qe,sfield,pott,ffc,mixup,sct,tfft,ntp,ndp,nx,ny,     &
     &modesxp,modesyp,iup,nprec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! calculate longitudinal electric field in fourier space
      call pois3(qe,fxyze,ffc,we,tfield,nx,ny,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! electromagnetic diagnostic
      call vpotrdiag(bxyz,cu,vfield,vpotr,ffc,mixup,sct,tfft,ci,nte,nde,&
     &nx,ny,modesxe,modesye,iue,nerec,indx,indy,ntime,ndstyle,irc,      &
     &inorder)
      if (irc==1) go to 2000
! calculate electromagnetic fields in fourier space
      if (ntime==0) then
! calculate initial darwin magnetic field
         call ibpois(cu,bxyz,ffc,ci,wm,nx,ny,inorder)
         wf = 0.
! calculate initial darwin electric field
         allocate(amu(4,nxe,nye),ffe(nxeh,nyh))
! deposit momentum flux
         call sguard(amu,zero,zero,zero,zero,nx,ny,inorder)
         call dmjpostglg(part,amu,np,qme,ci,nx,ny,ts,relativity)
         if (movion==1) then
            call dmjpostglg(parti,amu,npi,qmi,ci,nx,ny,ts,relativity)
         endif
! solve for darwin electric field
         bxyze = cu
         call dcuperp(bxyze,amu,nx,ny,inorder)
         call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny)
         call iepois(bxyze,exyz,ffe,ci,wf,nx,ny,inorder)
         deallocate(amu,ffe)
         dth = .5*dt
! calculate electromagnetic fields
      else
         call maxwel(exyz,bxyz,cu,ffc,ci,dt,wf,wm,tfield,nx,ny,inorder)
      endif
! vector potential diagnostic
      call avpotdiag(bxyz,vfield,vpott,mixup,sct,tfft,nta,nda,nx,ny,    &
     &modesxa,modesya,iua,narec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! add longitudinal and transverse electric fields
      isign = 1
      call emfield(fxyze,exyz,ffc,isign,nx,ny,inorder)
! copy magnetic field
      isign = -1
      call emfield(bxyze,bxyz,ffc,isign,nx,ny,inorder)
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxyze(1,:,:) = fxyze(1,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! push particles
      call push3glg(part,fxyze,bxyze,np,qbme,dt,dth,ci,wke,nx,ny,ipbc,  &
     &tpush,relativity)
! debug: zero force
!     call push3zfg(part,np,dt,ci,wke,tpush,nx,ny,ipbc,ndim,relativity)
      if (movion==1) then
         wki = 0.
         call push3glg(parti,fxyze,bxyze,npi,qbmi,dt,dth,ci,wki,nx,ny,  &
     &ipbc,tpushi,relativity)
         wki = wki*rmass
      endif
! calculate electron and field momentum
      call fmomtdiag(part,qe,ffc,exyz,bxyz,pxe,pye,pze,sx,sy,sz,wx,wy,wz&
     &,ntm,np,ium,nx,ny,ntime,inorder)
! calculate ion and total momentum
      call imomtdiag(parti,qi,fxyze,dt,rmass,pxi,pyi,pzi,wx,wy,wz,ntm,  &
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
      call emenergy(wt,we,wf,wm,wke,wki,ntw,ndw,itw,iuot,ntime)
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
! write primary file
            call restart_bwrite(iur1,itime,itime0,np,part,movion,npi,   &
     &parti,qi,exyz,bxyz)
            call restart_dwrite(iur1,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname,narec,faname,nerec,fename)
! write secondary file
            call restart_bwrite(iur2,itime,itime0,np,part,movion,npi,   &
     &parti,qi,exyz,bxyz)
            call restart_dwrite(iur2,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname,narec,faname,nerec,fename)
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
      write (iuot,*) 'gridless electromagnetic code bbeps2gl'
      write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost + tdjpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'electron current deposit time = ', tdjpost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      if (movion==1) then
         totpushi = tpushi + tdposti + tdjposti
         write (iuot,*) 'ion push time = ', tpushi, 'sec'
         write (iuot,*) 'ion charge deposit time = ', tdposti, 'sec'
         write (iuot,*) 'ion current deposit time = ', tdjposti, 'sec'
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
! electromagnetic diagnostics
      if (nte > 0) then
         nerec = nerec - 1
         ceng = affp
         write (iudm,em2d,iostat=irc)
         if (irc /= 0) write (iuer,*) 'em2d namelist not written'
      endif
! write out input file
      write (iudm,input2,iostat=irc)
      if (irc /= 0) write (iuer,*) 'input2 namelist not written'
! done
      write (iuot,*) '* * * q.e.d. * * *'
! close graphics device
 3000 call GRCLOSE
      call MP_END
      stop
      end program bbeps2gl
