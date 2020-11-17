    !-----------------------------------------------------------------------
! * * * periodic 2d electromagnetic particle simulation kernel code * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves electrons and ions, with periodic
! electromagnetic forces obtained by solving maxwell's equation with
! fast fourier transforms.
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G3
! copyright 1999, regents of the university of california
! update: march 2, 2012
      program bbeps2
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
      integer :: temp_pei
      integer :: iuvi = 20, iup = 11, iua = 15, iuj = 25, iue = 26
      integer :: ium = 21, iuer = 2
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
      real :: qbme, qbmi, affp, dth, qi0, omt, q2m0, wp0
      real :: etx, we, wf, wm, wke
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: sx, sy, sz, pxe, pye, pze, wx, wy, wz
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      real, dimension(:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:), pointer :: qe, qi
      real, dimension(:,:,:), pointer :: cu, amu, fxyze, bxyze
      real, dimension(:,:,:), pointer :: lsrfxyze,lsrbxyze
      real, dimension(:,:,:), pointer :: dampfxyze,dampbxyze
      real, dimension(:,:), pointer :: dampfactor
      integer :: dampinit_flag
      integer :: lsrii,lsrjj
      complex, dimension(:,:,:), pointer :: exyz, bxyz
      real, dimension(:,:,:), pointer :: exyzr
      complex, dimension(:,:), pointer :: ffc, ffe
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      real, dimension(:), pointer :: cu0
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:,:), pointer :: qt, sfield
      complex, dimension(:,:), pointer :: dent, pott
      real, dimension(:,:,:), pointer :: cui, cut, vfield
      complex, dimension(:,:,:), pointer :: vpott, vcurt, vpotr
      real, dimension(:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=32) :: fe1name,fe2name,fe3name
      character(len=32) :: fb1name,fb2name,fb3name
      character(len=32) :: fecdname,fline1name,fline2name,fline3name,ficdname
      character(len=32) :: flinb1name,flinb2name,flinb3name
      character(len=32) :: fecu1name,fecu2name,fecu3name
      integer :: plasmalength, plasmaylength,beamlength, beamylength
      
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
      idcode = 2
      psolve = 1
      ndim = 3
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file

      temp_pei = system("mkdir MS");
      temp_pei = system("mkdir MS/FLD");
      temp_pei = system("mkdir MS/FLD/e1");
      temp_pei = system("mkdir MS/FLD/e2");
      temp_pei = system("mkdir MS/FLD/e3");
      temp_pei = system("mkdir MS/FLD/b1");
      temp_pei = system("mkdir MS/FLD/b2");
      temp_pei = system("mkdir MS/FLD/b3");
      temp_pei = system("mkdir MS/FLD/e1r");
      temp_pei = system("mkdir MS/FLD/e2r");
      temp_pei = system("mkdir MS/FLD/e3r");
      temp_pei = system("mkdir MS/FLD/b1r");
      temp_pei = system("mkdir MS/FLD/b2r");
      temp_pei = system("mkdir MS/FLD/b3r");
      temp_pei = system("mkdir MS/FLD/e1line");
      temp_pei = system("mkdir MS/FLD/e2line");
      temp_pei = system("mkdir MS/FLD/e3line");
      temp_pei = system("mkdir MS/DENSITY");
      temp_pei = system("mkdir MS/DENSITY/electron");
      temp_pei = system("mkdir MS/DENSITY/ion");
      temp_pei = system("mkdir MS/PART");
      temp_pei = system("mkdir MS/PART/electron");
      
      iuot = get_funit(iuot)
      fname = 'MS//output2.'//cdrun
      open(unit=iuot,file=trim(fname),form='formatted',status='replace')
! np = total number of electrons in simulation
      nx = 2**indx; ny = 2**indy; 
      if (ndprof == 5) then
        plasmalength = nx - plasmabstart - plasmabend
        plasmaylength = ny - plasmabystart - plasmabyend
        beamlength = nx - beamstart - beamend
        beamylength = ny - beamystart - beamyend
        npx = npxg * plasmalength; npy = npyg * plasmaylength
        npxb = npxbg * beamlength; npyb = npybg * beamylength
      endif  
      npxy = npx*npy; npxyb = npxb*npyb; np = npxy + npxyb
      print *,' number of electron = ', np
! npi = total number of ions in simulation
      if (ndprof == 5) then
        plasmalength = nx - plasmabstart - plasmabend
        plasmaylength = ny - plasmabystart - plasmabyend
        npxi = npxg * plasmalength; npyi = npyg * plasmaylength
      endif
      npxyi = npxi*npyi; npxybi = npxbi*npybi; npi = npxyi + npxybi
      if (movion==1) print *,' number of ion = ', npi
      nxh = nx/2; nyh = max(1,ny/2)
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
      !ntmax = 100
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
      allocate(part(idimp,np))
      part = 0.0
      allocate(cu0(3))
! in real space, qe(j,k) = charge density at grid point (j,k)
! in real space, qi(j,k) = ion charge density at grid point (j,k)
      allocate(qe(nxe,nye),qi(nxe,nye))
! in real space, fxyze(i,j,k) = i component of force/charge at grid (j,k)
! in other words, fxyze are the convolutions of the electric field
! over the particle shape
      allocate(fxyze(3,nxe,nye))
      if (launchlas == 1) then
        allocate(lsrfxyze(3,nxe,nye)); allocate(lsrbxyze(3,nxe,nye))
      end if
      if (dampfield == 1) then
        allocate(dampfxyze(3,nxe,nye)); allocate(dampbxyze(3,nxe,nye))
        allocate(dampfactor(nxe,nye))
        dampinit_flag = 0
      endif
        
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
!
! open graphics device
      call GROPEN
      call SETNPLT(nplot,irc)
      call STPALIT(idpal)
! initialize timer
      call wtimer(time,ltime,-1)
! initialize constants
      if (ndprof == 5) then
        qme = -1.0 * edensity/real(npxg*npyg)
      endif
      print*, 'electron charge = ', qme
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
      if (ndprof == 5) then
        qbme = -1.0   ! electron q/m ratio
        !affp = real(plasmalength*ny)/real(np)*r2wpw0*(real(npxg*npyg)/edensity)
        affp = affpin
      else
        qbme = qme
        affp = real(nx*ny)/real(np)
      endif
      print*, 'electron charge/mass = ', qbme
      if (laseronly == 1) then
        qbme = zero
      endif
      !print *, "affp = ",affp
      qbmi = zero
      omt = sqrt(omx*omx + omy*omy + omz*omz)
      !q2m0 = qbme*qme*real(np)/real(nx*ny)
      q2m0 = qbme*qme/affp
      dth = 0.0
      !print *, qme, affp,np,npx,npy
! debug
!     dth = .5*dt
! end debug
      if (movion==1) then
         if (ndprof == 5) then
            ndprofi = ndprof
            qmi = - qme
            qbmi = 1.0/rmass
         else
            qbmi = qmi/rmass
         endif
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
         vtdzi = vtdz/sqrt(rmass*rtempdzi)
         !q2m0 = q2m0 + qbmi*qmi*real(npi)/real(nx*ny)
         q2m0 = q2m0 + qbmi*qmi/affp
      endif
      print*, 'ion charge = ', qmi
      print*, 'ion charge/mass = ', qbmi
      wp0 = q2m0*affp
! set initial time
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! set default diagnostic file names
      if (ntd > 0) fdname = 'MS//denk2.'//cdrun
      if (ntp > 0) fpname = 'MS//potk2.'//cdrun
      if (nta > 0) faname = 'MS//vpotk2.'//cdrun
      if (ntj > 0) fjname = 'MS//vcurk2.'//cdrun
      if (nte > 0) fename = 'MS//vpotrk2.'//cdrun
! EM diagnostics
      !fe1name = 'MS//e1.'//cdrun
      !fe2name = 'MS//e2.'//cdrun
      !fe3name = 'MS//e3.'//cdrun
      !fb1name = 'MS//b1.'//cdrun
      !fb2name = 'MS//b2.'//cdrun
      !fb3name = 'MS//b3.'//cdrun
      !fecdname = 'MS//elecCharDen.'//cdrun
      !ficdname = 'MS//ionCharDen.'//cdrun
      !fecu1name = 'MS//elecCurr1.'//cdrun
      !fecu2name = 'MS//elecCurr2.'//cdrun
      !fecu3name = 'MS//elecCurr3.'//cdrun
      !fline1name = 'MS//e1line.'//cdrun
      !fline2name = 'MS//e2line.'//cdrun
      !fline3name = 'MS//e3line.'//cdrun
      !flinb1name = 'MS//b1line.'//cdrun
      !flinb2name = 'MS//b2line.'//cdrun
      !flinb3name = 'MS//b3line.'//cdrun
      
      !open (unit = 12341, file = fe1name)
      !open (unit = 12342, file = fe2name)
      !open (unit = 12343, file = fe3name)
      !open (unit = 12344, file = fb1name)
      !open (unit = 12345, file = fb2name)
      !open (unit = 12346, file = fb3name)
      !open (unit = 12351, file = fecdname)
      !open (unit = 12352, file = ficdname)
      !open (unit = 12361, file = fecu1name)
      !open (unit = 12362, file = fecu2name)
      !open (unit = 12363, file = fecu3name)
      !open (unit = 12381, file = fline1name)
      !open (unit = 12382, file = fline2name)
      !open (unit = 12383, file = fline3name)
      !open (unit = 12384, file = flinb1name)
      !open (unit = 12385, file = flinb2name)
      !open (unit = 12386, file = flinb3name)
      open (unit = 12371, file = 'MS//lsrfxyze.20')
      !open (unit = 12372, file = 'MS//lsrbxyze.20')
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
      if (omt > 0) then
         call baddext(bxyze,omx,omy,omz,nx,ny,inorder)
         call cguard(bxyze,nx,ny,inorder)
      endif
      bxyz = cmplx(0.0,0.0)
      exyz = cmplx(0.0,0.0)
      cu = 0.0
! prepare fft tables
      call fft_init(mixup,sct,indx,indy)
! calculate form factors
      call pois_init(ffc,ax,ay,affp,nx,ny,smoothtype,smoothfactor,      &
     &smoothmask)
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npi))
         parti = 0.0
         allocate(parti2(0,0))
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         npxy1 = npxy + 1
! background electrons
         if (npxy > 0) then
            call fdistr(part,1,npxy,ampdx,scaledx,shiftdx,ampdy,scaledy,&
     &shiftdy,npx,npy,nx,ny,ipbc,ndprof,nsrand,1)
            call vdistr(part,1,npxy,vtx,vty,vtz,vx0,vy0,vz0)
         endif

! beam electrons       
         if (npxyb > 0) then
            call fdistr(part,npxy1,npxyb,ampdx,scaledx,shiftdx,ampdy,   &
     &scaledy,shiftdy,npxb,npyb,nx,ny,ipbc,ndprof,nsrand,2)
            call vdistr(part,npxy1,npxyb,vtdx,vtdy,vtdz,vdx,vdy,vdz)
         endif

! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
            call dpostg(part,qi,npxy,nx,ny,-qme,tdpost,inorder,dopt)
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
         if (ntm > 0) call initmomt2(part,np,pxe,pye,pze)
! initialize ions
         if (movion==1) then
            npxy1 = npxyi + 1
! background ions
            if (npxyi > 0) then
               call fdistr(parti,1,npxyi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,npxi,npyi,nx,ny,ipbc,ndprofi,nsrandi,1)
               call vdistr(parti,1,npxyi,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0)
            endif
! beam ions
            if (npxybi > 0) then
               call fdistr(parti,npxy1,npxybi,ampdxi,scaledxi,shiftdxi, &
     &ampdyi,scaledyi,shiftdyi,npxbi,npybi,nx,ny,ipbc,ndprofi,nsrandi,2)
               call vdistr(parti,npxy1,npxybi,vtdxi,vtdyi,vtdzi,vdxi,   &
     &vdyi,vdzi)
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
            if (ntm > 0) call initmomt2(parti,npi,pxi,pyi,pzi)  
         endif
! freeze the ions now
         if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density
            call dpostg(parti,qi,npi,nx,ny,qmi,tdposti,inorder,dopt)
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
            call restart_dread(it,itime,itw,wt,iud,ndrec,fdname,iup,    &
     &nprec,fpname,irc,iuer,iua,narec,faname,iuj,njrec,fjname,iue,nerec,&
     &fename)
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
      fname = 'MS//diag2.init.'//cdrun
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
         call initveldiag(fvi,fvmi,vtxi,vtyi,zero,ntv,ndv,nmv,ndim,iuvi,&
     &fname)
      endif
! potential diagnostics
      call initmodediag(pott,ntp,nxh,nyh,modesxp,modesyp,iup,nprec,     &
     &fpname)
      if (ntp > 0) then
         ceng = affp
         write (iudm,pot2d,iostat=irc)
      endif
! vector potential, ion current or electromagnetic diagnostics
      if ((nta > 0) .or. (nda > 0) .or. (ntj > 0) .or. (ndj > 0) .or.   &
     &(nte > 0) .or. (nde > 0)) then
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
! ion current diagnostics
      call initvmodediag(vcurt,ntj,nxh,nyh,ndim,modesxj,modesyj,iuj,    &
     &njrec,fjname)
      if (ntj > 0) then
         allocate(cui(ndim,nxe,nye),cut(ndim,nxe,nye))
         cui = 0.0; cut = 0.0
         ceng = zero
         write (iudm,vcur2d,iostat=irc)
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
!!!
!open(unit = 40, file = 'qe.bin', access = 'stream', status = 'new') !!!
!open(unit = 50, file = 'qi.bin', access = 'stream', status = 'new') !!!
!!!
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
      write (iuot,991) ntime
      print *, 'now at ntime/nloop = ', ntime,'/',nloop
!
! electromagnetic diagnostic
      call vpotdiagprep(cu,vfield,nte,nde,ntime-1)
! initialize current density to background
      call sguard(cu,zero,zero,zero,nx,ny,inorder)
! deposit electron current
      call djpostg(part,cu,np,qme,dth,ci,tdjpost,nx,ny,ipbc,relativity, &
     &inorder,djopt)
     
! current diagnostics
      !call ECurrdiag(cu,ntime,ntcu,nxe,nye,inorder)
! deposit electron charge density
      call dpostg(part,qe,np,nx,ny,qme,tdpost,inorder,dopt)
! electron charge density diagnostics
     ! write(40) qe(1:nx,1) !!!
      call ECharDendiag(qe, ntime, ntqe, nxe, nye, inorder)
! save electron current for ion current diagnostic
      call vpotrdiagprep(cu,cut,ntj,ndj,ntime)
! add ion current and deposit ion charge density
      if (movion==1) then
         call djpostg(parti,cu,npi,qmi,dth,ci,tdjposti,nx,ny,ipbc,      &
     &relativity,inorder,djopt)
! deposit ion charge density
         call dpostg(parti,qi,npi,nx,ny,qmi,tdposti,inorder,dopt)
        ! write(50) qi(1:nx,1) !!!
         call ICharDendiag(qi, ntime, ntqi, nxe, nye, inorder)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti)
            movion = 0
         endif
      endif
! save ion current
      call vpotrdiagprep(cu,cut,ntj,ndj,ntime,diff=.true.)
! add guard cells for current
      call acguard(cu,nx,ny,inorder)
! add electron and ion densities
      call addqei(qe,qi,nx,ny,inorder)
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
! transform current to fourier space
      isign = -1
      call fft(cu,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(cu,isign,mixup,sct,tfft,indx,indy,inorder)
! take transverse part of current
      call cuperp(cu,nx,ny,inorder)
! electromagnetic diagnostic
      call vpotrdiag(bxyz,cu,vfield,vpotr,ffc,mixup,sct,tfft,ci,nte,nde,&
     &nx,ny,modesxe,modesye,iue,nerec,indx,indy,ntime,ndstyle,irc,      &
     &inorder)
      if (irc==1) go to 2000
! ion current diagnostic
      call avcurdiag(cut,cui,vfield,vcurt,ffc,mixup,sct,tfft,ntj,ndj,nx,&
     &ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,ndstyle,1,irc,inorder&
     &)
      if (irc==1) go to 2000
      
! calculate electromagnetic fields in fourier space
      if (ntime==0) then
! calculate initial darwin magnetic field
         call ibpois(cu,bxyz,ffc,ci,wm,nx,ny,inorder)
         wf = 0.
         ! save current soub
! calculate initial darwin electric field
         allocate(amu(4,nxe,nye),ffe(nxeh,nyh))
! deposit momentum flux
         call sguard(amu,zero,zero,zero,zero,nx,ny,inorder)
         call dmjpostg(part,amu,np,qme,ci,ts,relativity,inorder,djopt)
         if (movion==1) then
            call dmjpostg(parti,amu,npi,qmi,ci,ts,relativity,inorder,   &
     &djopt)
         endif
         call amcguard(amu,nx,ny,inorder)
! solve for darwin electric field
         isign = -1
         call fftn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
! debug
!        call iwfft2rn(amu,isign,mixup,sct,tfft,indx,indy,inorder)
         bxyze = cu
         call dcuperp(bxyze,amu,nx,ny,inorder)
         call epois_init(ffe,ax,ay,affp,wp0,ci,nx,ny)
         call iepois(bxyze,exyz,ffe,ci,wf,nx,ny,inorder)
         !    SAVCUZ2(cu,cu0,nxvh,nyv)
         call isavcuz2(cu,cu0,inorder)

         deallocate(amu,ffe)
         dth = .5*dt
! calculate electromagnetic fields
      else
         call maxwel(exyz,bxyz,cu,ffc,ci,dt,wf,wm,tfield,nx,ny,inorder,solvertype)
      endif
      !    MAXZ2(exy,cu,cu0,ffc,dt,wf,nx,ny,nxvh,nyv,nxhd,nyhd)
      call imaxz2(exyz,cu,cu0,ffc,dt,wf,nx,ny,inorder)
      
      call EMdiagF(exyz,bxyz,ntime,nteme1,nteme2,nteme3, ntemb1,ntemb2, &
     &ntemb3, ntle1,ntle2,ntle3, ntlb1,ntlb2, ntlb3,nxe,nye,inorder)

!     launch a laser in real space in lsrfxyze and lsrbxyze
      if (launchlas == 1) then
        if (launchtype == 1) then
         !lsrfxyze = 0.0; lsrbxyze = 0.0
         call launchlaser(lsrfxyze,lsrbxyze,nxe,nye,nx,ny,ci,launchlas, &
     &lase0, laspol, lasxstart,lasxend,launchtype,lask, ntime, inorder, &
     &dt,ntlas,launchshape,lasw0,lasfocus,rayl,gamma,antlength)
         
!        transfer it to fourier space
         isign = -1
         call fft(lsrfxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         call fft(lsrbxyze,isign,mixup,sct,tfft,indx,indy,inorder) 
         !write (12371,'(ES24.14E2)') ((lsrfxyze(3,lsrii,lsrjj),lsrii=1,nx),lsrjj=1,ny)
         !write (12372,'(ES24.14E2)') ((lsrbxyze(2,lsrii,lsrjj),lsrii=1,nx),lsrjj=1,ny)
        
!        add them to exyz and bxyz
         call ilsremfield2(lsrfxyze,exyz,ffc,1,nx,ny,inorder)
         call ilsremfield2(lsrbxyze,bxyz,ffc,1,nx,ny,inorder)
        else if (launchtype == 2) then
         !print *, 'antenna'
         ! Fourier space (exyz and bxyz) to Fourier space (lsrfxyze and lsrbxyze)
         call ilsremfield4(lsrfxyze,exyz,ffc,1,nx,ny,inorder)
         call ilsremfield4(lsrbxyze,bxyz,ffc,1,nx,ny,inorder)
         ! Fourier space (lsrfxyze and lsrbxyze) to real space (lsrfxyze and lsrbxyze)
         isign = 1
         call fft(lsrfxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         call fft(lsrbxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         ! launch laser at specific section
         call launchlaser(lsrfxyze,lsrbxyze,nxe,nye,nx,ny,ci,launchlas, &
     &lase0, laspol, lasxstart,lasxend,launchtype,lask, ntime, inorder, &
     &dt,ntlas,launchshape,lasw0,lasfocus,rayl,gamma,antlength)
         ! real space (lsrfxyze and lsrbxyze) to Fourier space (exyz and bxyz)
         isign = -1
         call fft(lsrfxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         call fft(lsrbxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         ! lsrfxyze(full,Fourier) to exyz(half,Fourier)
         call ilsremfield3(lsrfxyze,exyz,ffc,1,nx,ny,inorder)
         call ilsremfield3(lsrbxyze,bxyz,ffc,1,nx,ny,inorder)
        endif
      else if (launchlas == 2) then
         deallocate(lsrfxyze); deallocate(lsrbxyze)
         launchlas = 0
      endif
!     finish launching a laser, EM fields are added to exyz and bxyz 

!   damp the electric field
         ! exyz(half,Fourier) to exyz(full,Fourier)
      if (dampfield == 1) then
         print *,'damping the E-field'
         ! initialize damping form factor
         call dampfactorinit(dampfactor,nxe,nye,nx,ny,inorder,&
      &dampystart,dampyend,dampinit_flag, dampdrop)
         ! damp the E-field
         dampfxyze = 0.0
         call ilsremfield4(dampfxyze,exyz,ffc,1,nx,ny,inorder)
         !call ilsremfield4(dampbxyze,bxyz,ffc,1,nx,ny,inorder)
         ! exyz(full,Fourier) to exyz(full,real)
         isign = 1
         call fft(dampfxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         !call fft(dampbxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         ! damp it
         !call dampEM
         call dampEM(dampfxyze,dampfactor,nxe,nye,nx,ny,inorder,dampe1,dampe2,dampe3)
         ! exyz(full,real) to exyz(full, Fourier)
         isign = -1
         call fft(dampfxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         !call fft(dampbxyze,isign,mixup,sct,tfft,indx,indy,inorder)
         ! exyz(full,Fourier) to exyz(half,Fourier)
         call ilsremfield3(dampfxyze,exyz,ffc,1,nx,ny,inorder)
         !call ilsremfield3(dampbxyze,bxyz,ffc,1,nx,ny,inorder)
      endif
!   finish damping the electric field

! vector potential diagnostic
      call avpotdiag(bxyz,vfield,vpott,mixup,sct,tfft,nta,nda,nx,ny,    &
     &modesxa,modesya,iua,narec,indx,indy,ntime,ndstyle,irc,inorder)
      if (irc==1) go to 2000
! calculate longitudinal electric field in fourier space
      call pois3(qe,fxyze,ffc,we,tfield,nx,ny,inorder)
      ! correct the EM field for k = 0
      
! add longitudinal and transverse electric fields
      isign = 1
      call emfield(fxyze,exyz,ffc,isign,nx,ny,inorder)
! copy magnetic field
      isign = -1
      call emfield(bxyze,bxyz,ffc,isign,nx,ny,inorder)
! transform force/charge to real space
      isign = 1
      call fft(fxyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(fxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      call cguard(fxyze,nx,ny,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,isign,mixup,sct,tfft,indx,indy,inorder)
!     call fftn(bxyze,isign,mixup,sct,tfft,indx,indy,inorder)
      if (omt > 0) call baddext(bxyze,omx,omy,omz,nx,ny,inorder)
      call cguard(bxyze,nx,ny,inorder)
         
! EM diagnostics in real space
      call EMdiag(fxyze,bxyze,ntime,nteme1,nteme2,nteme3, ntemb1,ntemb2, &
     &ntemb3, ntle1,ntle2,ntle3, ntlb1,ntlb2, ntlb3,nxe,nye,inorder)
      
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxyze(1,:,:) = fxyze(1,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! push particles
      call push3g(part,fxyze,bxyze,np,qbme,dt,dth,ci,wke,tpush,nx,ny,   &
     &ipbc,relativity,inorder,popt)
! debug: zero force
!     call push3zfg(part,np,dth,ci,wke,tpush,nx,ny,ipbc,ndim,relativity)
      if (movion==1) then
         wki = 0.
         call push3g(parti,fxyze,bxyze,npi,qbmi,dt,dth,ci,wki,tpushi,nx,&
     &ny,ipbc,relativity,inorder,popt)
! debug: zero force
!        call push3zfg(parti,npi,dth,ci,wki,tpushi,nx,ny,ipbc,ndim,     &
!    &relativity)
         wki = wki*rmass
      endif

      call PartDiag(part,ntime,ntpart,idimp, np)

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
     &fpname,narec,faname,njrec,fjname,nerec,fename)
! write secondary file
            call restart_bwrite(iur2,itime,itime0,np,part,movion,npi,   &
     &parti,qi,exyz,bxyz)
            call restart_dwrite(iur2,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname,narec,faname,njrec,fjname,nerec,fename)
         endif
      endif
      call wtimer(tloop,ltime)
      time = time + tloop
!!!      
     ! isign = 1
     ! call fft(qe,isign,mixup,sct,tfft,indx,indy,inorder)
     ! call fft(qi,isign,mixup,sct,tfft,indx,indy,inorder)
     ! write(40) qe(1:nx,1) !!!
     ! write(50) qi(1:nx,1) !!!
     ! isign = -1
     ! call fft(qe,isign,mixup,sct,tfft,indx,indy,inorder)
     ! call fft(qi,isign,mixup,sct,tfft,indx,indy,inorder)
!!!      
      
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
      write (iuot,*) 'electromagnetic code bbeps2'
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
      fname = 'MS//diag2.'//cdrun
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
! electromagnetic diagnostics
      if (nte > 0) then
         nerec = nerec - 1
         ceng = affp
         write (iudm,em2d,iostat=irc)
         if (irc /= 0) write (iuer,*) 'em2d namelist not written'
      endif
      
!!!
     ! close(40) !!!
    !  close(50)
!!!      
! write out input file
      write (iudm,input2,iostat=irc)
      if (irc /= 0) write (iuer,*) 'input2 namelist not written'
! done
      write (iuot,*) '* * * q.e.d. * * *'
      close(unit=iudm)
      close(unit=iuot)
      close(unit=iuer)
      
!     cleanup
      if (dampfield == 1) then
        deallocate(dampfxyze,dampbxyze,dampfactor)
      endif
      
! close graphics device
 3000 call GRCLOSE
      call MP_END
      stop
      end program bbeps2
