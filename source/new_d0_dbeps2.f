!-----------------------------------------------------------------------
! * * * periodic 2d darwin particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with 
! electromagnetic forces obtained by solving non-radiative form of
! maxwell's equation with fast fourier transforms for various boundary
! conditions, using algorithm similar to that described in
! J. Busnardo-Neto, P. L. Pritchett, A. T. Lin, and J. M. Dawson,
! J. Computational Phys. 23, 300 (1977).
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 2006, regents of the university of california
! update: march 1, 2012
      program d0_dbeps2
      use init2d
      use empush2d
      use npfield2d
! debug
      use dfield2d, only: icuperpdx2
! end debug
      use diag2d
      use emsimul2d
      use mp0d, only: mpinit, ncpus
      implicit none
! nmv = number of segments in v for velocity distribution
      integer :: nmv = 40
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, iua = 15, ium = 21, iuer = 2
      integer :: npxy, npxyb, np, npxyi, npxybi, npi, npxy1
      integer :: nx, ny, nxh, nyh, nxe, nye
      integer :: nxye, nxeh, nxyh, nxhy, nx1, ny1, ntmax, idimp, ipbc
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: j, k, it, jt, itw, iur1, iur2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, time = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0, ts = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0 
      real :: tfft = 0.0, tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: tdcjpost = 0.0, tdcjposti = 0.0
      real :: qbme, qbmi, affp, qi0, q2m0, wp0, wpmax, wpmin
      real :: etx, we, wf, wm, wke
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: sx, sy, sz, pxe, pye, pze, wx, wy, wz
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      real, dimension(:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:), pointer :: qe, qi
      real, dimension(:,:,:), pointer :: fxyze, exyze, bxyze
      real, dimension(:,:,:), pointer :: cu, cus, amu, ecu
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
      real, dimension(:,:), pointer :: wt
! dirichlet or neumann boundary conditions
      integer :: indx1, nx2, nxv, nx2v
      integer :: indy1, ny2
      real, dimension(:,:), pointer :: q2, sfield2
      real, dimension(:,:,:), pointer :: cu2, cus2, amu2
      real, dimension(:,:,:), pointer :: fxy2, exy2, bxy2
      complex, dimension(:,:), pointer :: ffd, fff
      integer, dimension(:), pointer :: mixup2
      complex, dimension(:), pointer :: sct2
!
      character(len=10) :: cdrun
      character(len=32) :: fname
      integer, external :: NDIAN, NDPREC, IDPREC
! debug
      integer :: i, jj, kk, kk2
      real :: wg, eps, epsmax, sum1, at1
      complex zt1
      real, dimension(:,:,:), pointer :: gu, gus, gmu
      real, dimension(:,:,:), pointer :: cuk, cuk2, gus2
      real, dimension(:,:,:), pointer :: gxyze
      complex, dimension(:,:,:), pointer :: cuc, guc, cuc2, guc2
! end debug
  991 format (' T = ',i7)
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
      idcode = 6
!     ndim = 3
      if (psolve > 2) psolve = 2
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
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871
         nxe = nx + 2; nye = ny + 1
      endif
! boundary conditions
      ipbc = psolve
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
! idimp = dimension of phase space = 4 or 5
      idimp = 2 + ndim
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
      allocate(part(idimp,np))
! in real space, qe(j,k) = electron charge density at grid point (j,k)
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
!
      allocate(ecu(ndim,nxe,nye))
!
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
      fxyze = 0.0; cu = 0.0
! debug
      allocate(gu(ndim,nxe,nye),gus(ndim,nxe,nye))
      allocate(gmu(2*ndim-2,nxe,nye),gxyze(ndim,nxe,nye))
      allocate(cuk(ndim,nxe,nye),cuk2(ndim,2*(nx+2),2*ny))
      allocate(gus2(ndim,2*(nx+2),2*ny))
      allocate(cuc(ndim,nxeh,ny),guc(ndim,nxeh,ny))
      allocate(cuc2(ndim,nx+2,2*ny),guc2(ndim,nx+2,2*ny))
! end debug
! non-periodic boundary conditions
      indx1 = indx + 1; indy1 = indy + 1
      nxv = nx + 2; nx2v = 2*nxv; ny2 = 2*ny; nx2 = 2*nx
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
         allocate(q2(nx2v,ny2),sfield2(nx2v,ny2))
         allocate(cu2(ndim,nx2v,ny2),cus2(ndim,nx2v,ny2))
         allocate(amu2(2*ndim-2,nx2v,ny2))
         allocate(fxy2(ndim,nx2v,ny2),exy2(ndim,nx2v,ny2))
         allocate(bxy2(2*ndim-3,nx2v,ny2))
         allocate(ffd(nxv,ny),fff(nxv,ny))
         allocate(mixup2(2*nxhy))
         allocate(sct2(2*nxyh))
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
      if (ipbc==1) then
         affp = float(nx*ny)/float(np)
      else if (ipbc==2) then
         affp = float((nx-2)*(ny-2))/float(np)
      else if (ipbc==3) then
         affp = float((nx-2)*ny)/float(np)
      endif
      qbmi = zero
      q2m0 = qbme*qme/affp
      if (movion==1) then
         qbmi = qmi/rmass
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
         vtdzi = vtdz/sqrt(rmass*rtempdzi)
         if (ipbc==1) then
            q2m0 = q2m0 + qbmi*qmi*real(npi)/real(nx*ny)
         else if (ipbc==2) then
            q2m0 = q2m0 + qbmi*qmi*real(npi)/float((nx-2)*(ny-2))
         else if (ipbc==3) then
            q2m0 = q2m0 + qbmi*qmi*real(npi)/float((nx-2)*ny)
         endif
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
! energy diagnostics
      if (ndw > 0) then
         allocate(wt((nloop-1)/ntw-(itime0/ntw)+1,7))
         itw = 0
      endif
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
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
         call fst_init(mixup,sct2,indx,indy)
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisd_init(ffd,ax,ay,affp,nx,ny)
         call epoisd_init(fff,ax,ay,affp,wp0,ci,nx,ny)
      endif
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npi))
         allocate(parti2(0,0))
      endif
! debug
      if (ipbc==2) vdy = 0.
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
            call vdistr(part,npxy1,npxyb,vtdy,vtdy,vtdz,vdx,vdy,vdz,ndim&
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
            call sguardp(qi,zero,nx,ny,ipbc,inorder)
            call dpost(part,qi,np,-qme,tdpost,inorder,dopt)
            call aguardp(qi,nx,ny,ipbc,inorder)
! debug
!           call sguardp(qi,qi0,nx,ny,ipbc,inorder)
! freeze the ions now
         else if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density
            call sguardp(qi,zero,nx,ny,ipbc,inorder)
            call dpost(parti,qi,npi,qmi,tdposti,inorder,dopt)
            call aguardp(qi,nx,ny,ipbc,inorder)
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
            call restart_dread(it,itime,itw,wt,iud,ndrec,fdname,iup,    &
     &nprec,fpname,irc,iuer,iua,narec,faname)
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
! vector potential diagnostics
      if ((nta > 0) .or. (nda > 0)) then
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
! debug
      write (71,991) ntime
! end debug
!
! initialize current density and momentum flux to background
      call sguardp(cu,zero,zero,zero,nx,ny,ipbc,inorder)
      call sguardp(amu,zero,zero,zero,zero,nx,ny,ipbc,inorder)
! deposit current and momentum flux for electrons
      call djpostg(part,cu,np,qme,zero,ci,tdjpost,nx,ny,ipbc,relativity,&
     &inorder,djopt)
      call dmjpostg(part,amu,np,qme,ci,tdcjpost,relativity,inorder,djopt&
     &)
! deposit electron charge density
      call sguardp(qe,zero,nx,ny,ipbc,inorder)
      call dpost(part,qe,np,qme,tdpost,inorder,dopt)
      call aguardp(qe,nx,ny,ipbc,inorder)
! deposit current, momentum flux, and charge density for ions
      if (movion==1) then
         call djpostg(parti,cu,npi,qmi,zero,ci,tdjposti,nx,ny,ipbc,     &
     &relativity,inorder,djopt)
         call dmjpostg(parti,amu,npi,qmi,ci,tdcjposti,relativity,inorder&
     &,djopt)
         call sguardp(qi,zero,nx,ny,ipbc,inorder)
         call dpost(parti,qi,npi,qmi,tdposti,inorder,dopt)
         call aguardp(qi,nx,ny,ipbc,inorder)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti)
            movion = 0
         endif
      endif
! add guard cells for current and momentum flux
      call aguardp(cu,nx,ny,ipbc,inorder)
      call amguardp(amu,nx,ny,ipbc,inorder)
! add electron and ion densities
      call addqei(qe,qi,qbme,qbmi,wpmax,wpmin,nx,ny,inorder)
!     wp0 = 0.5*(wpmax + wpmin)
! recalculate form factors
!     if ((wp0 > 1.15*q2m0) .or. (wp0 < 0.85*q2m0)) then
!        q2m0 = wp0
!        wp0 = affp*wp0
!        q2m = cmplx(q2m0,wp0)
!        if (psolve==DIRICHLET_2D) then
!           call epoisd_init(fff,ax,ay,affp,wp0,ci,nx,ny)
!        else if (psolve==PERIODIC_2D) then
!           call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny)
!        endif
!        write (iuer,*) ntime, 'new shift constants,q2m0,wp0=', q2m
!     endif
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
!
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_2D) then
! ion density diagnostic
      if ((ntd > 0) .or. (ndd > 0)) then
         it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
         jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
         if ((it==0) .or. (jt==0)) then
            qt = qi
! transform ion density to fourier space
!           call dblsin(qt,q2,nx,ny,inorder)
            isign = -1
!           call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
            call fsst(qt,isign,mixup,sct2,tfft,indx,indy,inorder)
! calculate smoothing in fourier space
!           call poisdx(q2,sfield2,ffd,nx,ny)
            call poisd(qt,sfield,ffd,nx,ny,inorder)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(sfield,dent,nx,ny,modesxd,modesyd,inorder)
! write diagnostic output
!              if (ndrec==0) then
!                 ndrec = -1; iud = get_funit(12)
!                 call bfopen(dent,modesxd,modesy2d,iud,ndrec,trim(fdnam&
!    &e))
!              endif
!              call writebf(dent,modesxd,modesy2d,iud,ndrec,order=LINEAR&
!    &)
!           endif
! transform ion density to real space
            if (jt==0) then
               isign = 1
!              call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,LINEA&
!    &R)
!              call hafdbl(sfield,sfield2,nx,ny,inorder)
               call fsst(sfield,isign,mixup,sct2,tfft,indx,indy,inorder)
               call cguardp(sfield,nx,ny,ipbc,inorder)
! display ion density
               call displays(sfield,' I DENSITY',ntime,999,2,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! debug
      call dblsin(qe,q2,nx,ny,inorder)
      isign = -1
      call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! end debug
! transform charge to fourier space
!     call dblsin(qe,q2,nx,ny,inorder)
      isign = -1
!     call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call fsst(qe,isign,mixup,sct2,tfft,indx,indy,inorder)
! potential diagnostic
      if ((ntp > 0) .or. (ndp > 0)) then
         it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
         jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
         if ((it==0) .or. (jt==0)) then
! solve for potential
!           call poisdx(q2,sfield2,ffd,we,nx,ny)
            call poisd(qe,sfield,ffd,we,nx,ny,inorder)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(sfield,pott,nx,ny,modesxp,modesyp,inorder)
! write diagnostic output
!              if (nprec==0) then
!                 nprec = -1; iup = get_funit(11)
!                 call bfopen(pott,modesxp,modesy2p,iup,nprec,trim(fpnam&
!    &e))
!              endif
!              call writebf(pott,modesxp,modesy2p,iup,nprec,order=LINEAR&
!    &)
!           endif
! transform potential to real space
            if (jt==0) then
               isign = 1
!              call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,LINEA&
!    &R)
!              call hafdbl(sfield,sfield2,nx,ny,inorder)
               call fsst(sfield,isign,mixup,sct2,tfft,indx,indy,inorder)
               call cguardp(sfield,nx,ny,ipbc,inorder)
! display potential
               call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! calculate longitudinal electric field in fourier space
!     call poisdx(q2,fxy2,ffd,we,nx,ny)
      call poisd(qe,fxyze,ffd,we,nx,ny,inorder)
! transform longitudinal electric field to real space
      isign = 1
      call fcst(fxyze,isign,mixup,sct2,tfft,indx,indy,inorder)
!     call fft(fxy2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
!     call hafdbl(fxyze,fxy2,nx,ny,inorder)
      call cguardp(fxyze,nx,ny,ipbc,inorder)
! debug
      call dblsin(cu,cu2,nx,ny,inorder)
      isign = -1
      call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call icuperpdx2(cu2,nx,ny)
! end debug
! transform current to fourier space
!     call dblsin(cu,cu2,nx,ny,inorder)
      isign = -1
!     call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call fcst(cu,isign,mixup,sct2,tfft,indx,indy,inorder)
! take transverse part of current
!     call cmfieldd(cu2,cu,nx,ny)
      call cuperpd(cu,nx,ny,inorder)
! calculate magnetic field in fourier space
!     call bpoisdx(cu2,bxy2,ffd,ci,wm,nx,ny)
      call bpoisd(cu,bxyze,ffd,ci,wm,nx,ny,inorder)
! calculate the momentum in the darwin field
      if (ntm > 0) then
         it = ntime/ntm
         if (ntime==ntm*it) then
            call poyntdx(q2,cu2,ffd,ci,sx,sy,sz,nx,ny)
         endif
      endif
! debug
      write (*,*) 'wm,sx,sy,sz=',wm,sx,sy,sz
! end debug
! vector potential diagnostic
      if (ndc==0) call vpotdiagprep(cu,vfield,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fsct(bxyze,isign,mixup,sct2,tfft,indx,indy,inorder)
!     call fft(bxy2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
!     call hafdbl(bxyze,bxy2,nx,ny,inorder)
      call cguardp(bxyze,nx,ny,ipbc,inorder)
! transform momentum flux to fourier space
      call dblsinm(amu,amu2,nx,ny,inorder)
      isign = -1
      call fftn(amu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! take transverse part of time derivative of current
      call dcuperpdx(cus2,amu2,nx,ny)
! calculate convective part of transverse electric field
      call epoisdx(cus2,cu2,fff,ci,wf,nx,ny)
! debug
      write (*,*) 'wf=',wf
      cuk2 = cu2
      cu = 0.0
! end debug
! transform transverse electric field to real space
      isign = 1
!     call fcst(cu,isign,mixup,sct2,tfft,indx,indy,inorder)
      call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call hafdbl(cu,cu2,nx,ny,inorder)
      call cguardp(cu,nx,ny,ipbc,inorder)
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_2D) then
! ion density diagnostic
      call dendiag(qt,qi,sfield,dent,ffc,mixup,sct,tfft,ntd,ndd,nx,ny,  &
     &modesxd,modesyd,iud,ndrec,indx,indy,ntime,ndstyle,1,irc,inorder)
      if (irc==1) go to 2000
! transform charge to fourier space
      isign = -1
      call fft(qe,isign,mixup,sct,tfft,indx,indy,inorder)
! potential diagnostic
      call potdiag(qe,sfield,pott,ffc,mixup,sct,tfft,ntp,ndp,nx,ny,     &
     &modesxp,modesyp,iup,nprec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
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
      call cguard(bxyze,nx,ny,inorder)
! transform momentum flux to fourier space
      isign = -1
      call fftn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of time derivative of current
      call dcuperp(cus,amu,nx,ny,inorder)
! calculate convective part of transverse electric field
      call epois(cus,cu,ffe,ci,wf,tfield,nx,ny,inorder)
! debug
      cuk = cu
! end debug
! transform transverse electric field to real space
      isign = 1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(cu,nx,ny,inorder)
!
      endif
!
! add longitudinal and transverse fields
      exyze = cu + fxyze
!
! main iteration loop
      do k = 1, ndc
! debug
      call displayv(cu,' ET FIELD',ntime,999,1,ndstyle,nx,ny,irc,inorder&
     &)
      if (irc==1) go to 2000
      call sguardp(gus,zero,zero,zero,nx,ny,ipbc,inorder)
      call sguardp(gu,zero,zero,zero,nx,ny,ipbc,inorder)
      call sguardp(gmu,zero,zero,zero,zero,nx,ny,ipbc,inorder)
      call dcjpostg(part,exyze,bxyze,gu,gus,gmu,np,qme,qbme,dt,ci,tdcjpo&
     &st,relativity,inorder,popt)
      call aguardp(gu,nx,ny,ipbc,inorder)
      call aguardp(gus,nx,ny,ipbc,inorder)
      call amguardp(gmu,nx,ny,ipbc,inorder)
! end debug
! initialize current, acceleration density and momentum flux
      call sfguardp(cus,cu,ecu,q2m0,nx,ny,ipbc,inorder)
      call sguardp(cu,zero,zero,zero,nx,ny,ipbc,inorder)
      call sguardp(amu,zero,zero,zero,zero,nx,ny,ipbc,inorder)
! deposit current and acceleration density and momentum flux
      call dcjpostg(part,exyze,bxyze,cu,cus,amu,np,qme,qbme,dt,ci,      &
     &tdcjpost,relativity,inorder,popt)
      if (movion==1) then
         call dcjpostg(parti,exyze,bxyze,cu,cus,amu,npi,qmi,qbmi,dt,ci, &
     &tdcjposti,relativity,inorder,popt)
      endif
! add guard cells for current, acceleration density, and momentum flux
      call aguardp(cu,nx,ny,ipbc,inorder)
      call aguardp(cus,nx,ny,ipbc,inorder)
      call afguardp(cus,ecu,q2m0,nx,ny,ipbc,inorder)
      call amguardp(amu,nx,ny,ipbc,inorder)
!
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_2D) then
! debug
      call dblsin(cu,cu2,nx,ny,inorder)
      isign = -1
      call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call icuperpdx2(cu2,nx,ny)
! end debug
! transform current to fourier space
!     call dblsin(cu,cu2,nx,ny,inorder)
      isign = -1
!     call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call fcst(cu,isign,mixup,sct2,tfft,indx,indy,inorder)
! take transverse part of current
!     call cmfieldd(cu2,cu,nx,ny)
      call cuperpd(cu,nx,ny,inorder)
! calculate magnetic field in fourier space
!     call bpoisdx(cu2,bxy2,ffd,ci,wm,nx,ny)
      call bpoisd(cu,bxyze,ffd,ci,wm,nx,ny,inorder)
! calculate the momentum in the darwin field
       call poyntdx(q2,cu2,ffd,ci,sx,sy,sz,nx,ny)
! debug
      write (*,*) 'wm,sx,sy,sz=',k,wm,sx,sy,sz
! end debug
! vector potential diagnostic
      if (k==ndc) call vpotdiagprep(cu,vfield,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fsct(bxyze,isign,mixup,sct2,tfft,indx,indy,inorder)
!     call fft(bxy2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
!     call hafdbl(bxyze,bxy2,nx,ny,inorder)
      call cguardp(bxyze,nx,ny,ipbc,inorder)
! transform acceleration density and momentum flux to fourier space
      call dblsin(cus,cus2,nx,ny,inorder)
      isign = -1
      call fft(cus2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call dblsinm(amu,amu2,nx,ny,inorder)
      isign = -1
      call fftn(amu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! debug
!     call ivccopy(cus2,cuc2,nx,ny2,LINEAR)
!     call dblsin(gus,gus2,nx,ny,inorder)
!     isign = -1
!     call fft(gus2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
!     gus2 = gus2 - q2m0*cuk2
!     call ivccopy(gus2,guc2,nx,ny2,LINEAR)
!     epsmax = 0.0
!     do kk = 1, ny2
!     do jj = 1, nx
!     do i = 1, 3
!     eps = abs(cuc2(i,jj,kk)-guc2(i,jj,kk))
!     if (eps > epsmax) then
!        write (71,*) i,jj,kk,cuc2(i,jj,kk),guc2(i,jj,kk),eps
!        epsmax = eps
!     endif
!     enddo
!     enddo
!     enddo
!     write (71,*) 'cus2 epsmax=',k,epsmax
! end debug
! take transverse part of time derivative of current
      call adcuperpdx(cus2,amu2,nx,ny)
! calculate convective part of transverse electric field
      call epoisdx(cus2,cu2,fff,ci,wf,nx,ny)
! debug
      write (*,*) 'wf=',k,wf
      cuk2 = cu2
! end debug
! transform transverse electric field to real space
      isign = 1
!     call fcst(cu,isign,mixup,sct2,tfft,indx,indy,inorder)
      call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call hafdbl(cu,cu2,nx,ny,inorder)
      call cguardp(cu,nx,ny,ipbc,inorder)
! debug
      call dblsin(gus,cus2,nx,ny,inorder)
      isign = -1
      call fft(cus2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call dblsinm(gmu,amu2,nx,ny,inorder)
      isign = -1
      call fftn(amu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call adcuperpdx(cus2,amu2,nx,ny)
      call epoisdx(cus2,cu2,ffd,ci,wg,nx,ny)
      write (*,*) 'wg=',k,wg   
      isign = 1
      call fft(cu2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call hafdbl(gu,cu2,nx,ny,inorder)
      call cguardp(gu,nx,ny,ipbc,inorder)
! end debug
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_2D) then
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! calculate magnetic field in fourier space
      call bpois(cu,bxyze,ffc,ci,wm,tfield,nx,ny,inorder)
! calculate the momentum in the darwin field
      call poynt(qe,cu,ffc,ci,sx,sy,sz,nx,ny,inorder)
! vector potential diagnostic
      if (k==ndc) call vpotdiagprep(cu,vfield,nta,nda,ntime)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(bxyze,nx,ny,inorder)
! transform acceleration density and momentum flux to fourier space
      isign = -1
      call fft(cus,isign,mixup,sct,tfft,indx,indy,inorder)
      call fftn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! debug
!     call ivccopy(cus,cuc,nxh,ny,inorder)
!     isign = -1
!     call fft(gus,isign,mixup,sct,tfft,indx,indy,inorder)
!     cuk = gus - q2m0*cuk
!     call ivccopy(cuk,guc,nxh,ny,inorder)
!     epsmax = 0.0
!     do kk = 1, ny
!     do jj = 1, nxh
!     do i = 1, 3
!     eps = abs(cuc(i,jj,kk)-guc(i,jj,kk))
!     if (eps > epsmax) then
!        write (71,*) i,jj,kk,cuc(i,jj,kk),guc(i,jj,kk),eps
!        epsmax = eps
!     endif
!     enddo
!     enddo
!     enddo
!     write (71,*) 'cus epsmax=',k,epsmax
! end debug
! take transverse part of time derivative of current
      call adcuperp(cus,amu,nx,ny,inorder)
! calculate convective part of transverse electric field
      call epois(cus,cu,ffe,ci,wf,tfield,nx,ny,inorder)
! debug
      cuk = cu
! end debug
! transform transverse electric field to real space
      isign = 1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(cu,nx,ny,inorder)
! debug
      isign = -1
      call fft(gus,isign,mixup,sct,tfft,indx,indy,inorder)
      call fftn(gmu,isign,mixup,sct,tfft,indx,indy,inorder)
      call adcuperp(gus,gmu,nx,ny,inorder)
      call epois(gus,gu,ffc,ci,wf,tfield,nx,ny,inorder)
      isign = 1
      call fft(gu,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(gu,nx,ny,inorder)
! end debug
!
      endif
!
! add longitudinal and transverse fields
      exyze = cu + fxyze
! debug
      call displayv(cu,' ET FIELD',ntime,999,1,ndstyle,nx,ny,irc,inorder&
     &)
      if (irc==1) go to 2000
      epsmax = 0.0
      sum1 = 0.0
      do kk = 1, nye
      do jj = 1, nxe
      do i = 1, 3
      eps = abs(cu(i,jj,kk)-gu(i,jj,kk))
      sum1 = sum1 + eps
      if (eps > epsmax) then
!        write (71,*) i,jj,kk,cu(i,jj,kk),gu(i,jj,kk),eps
         epsmax = eps
      endif
      enddo
      enddo
      enddo
      write (71,*) 'exyze epsmax,sum1=',k,epsmax,sum1
! end debug
      enddo
!
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_2D) then
! vector potential diagnostic
      if ((nta > 0) .or. (nda > 0)) then
         it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
         jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
         if ((it==0) .or. (jt==0)) then
! calculate vector potential in fourier space
            call apoisd(vfield,cus,ffd,ci,wm,nx,ny,inorder)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(cus,vpott,nx,ny,modesxa,modesya,inorder)
! write diagnostic output
!              if (narec==0) then
!                 narec = -1; iua = get_funit(15)
!                 call bfopen(vpott,modesxa,modesy2a,iua,narec,trim(fana&
!    &me))
!              endif
!              call writebf(vpott,modesxa,modesy2a,iua,narec,order=LINEA&
!    &R)
!           endif
! transform vector potential to real space
            if (jt==0) then
               isign = 1
               call fcst(cus,isign,mixup,sct2,tfft,indx,indy,inorder)
               call cguardp(cus,nx,ny,ipbc,inorder)
! display absolute value of vector potential
               call displayv(cus,' VECTOR POTENTIAL',ntime,999,1,ndstyle&
     &,nx,ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_2D) then
! vector potential diagnostic
      call vpotdiag(vfield,cus,vpott,ffc,mixup,sct,ci,tfft,nta,nda,nx,ny&
     &,modesxa,modesya,iua,narec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
!
      endif
!
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        exyze(1,:,:) = exyze(1,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! push particles
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
! calculate electron and field momentum
      if (ntm > 0) then
         it = ntime/ntm
         it = ntime - ntm*it + 1
         if (it > 1) it = it - ntm
         if (it >= 0) then
            call premoment2(part,ntime,np,ium,pxe,pye,pze,sx,sy,sz,wx,wy&
     &,wz,ndim,nprint=it)
! calculate ion momentum
            if (movion==0) then
               if (it==1) then
                  call imoment(qi,exyze,ium,pxi,pyi,pzi,dt,wx,wy,wz,nx, &
     &ny,inorder)
               endif
            else if (movion==1) then
               call primoment2(parti,npi,ium,rmass,pxi,pyi,pzi,wx,wy,wz,&
     &ndim,nprint=it)
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
     &fpname,narec,faname)
! write secondary file
            call restart_bwrite(iur2,itime,itime0,np,part,movion,npi,   &
     &parti,qi,q2m)
            call restart_dwrite(iur2,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname,narec,faname)
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
      write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost + tdjpost + tdcjpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'electron current deposit time = ', tdjpost, 'sec'
      write (iuot,*) 'electron current derivative deposit time = ', tdcj&
     &post, 'sec'
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
      time = time - (totpush + totpushi + tfft)
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
! write out input file
      write (iudm,input2,iostat=irc)
      if (irc /= 0) write (iuer,*) 'input2 namelist not written'
! done
      write (iuot,*) '* * * q.e.d. * * *'
! close graphics device
 3000 call GRCLOSE
      call MP_END
      stop
      end program d0_dbeps2
