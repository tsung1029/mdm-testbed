!-----------------------------------------------------------------------
! * * * periodic 2d electrostatic particle simulation kernel code * * *
! this is a simple 2d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves electrons and ions, with electrostatic
! forces obtained by solving poisson's equation with fast fourier
! transforms for various boundary conditions
! written by viktor k. decyk, ucla
! Fortran 90 for Macintosh G5
! copyright 1999, regents of the university of california
! update: march 1, 2012
      program d0_beps2
      use init2d
      use espush2d
      use npfield2d
      use vfield2d
      use diag2d
      use simul2d
      use fft12d
      use mp0d, only: mpinit, ncpus
      implicit none
! idimp = dimension of phase space = 4
! nmv = number of segments in v for velocity distribution
      integer :: idimp = 4, nmv = 40
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, iuer = 2
      integer :: npxy, npxyb, np, npxyi, npxybi, npi, npxy1
      integer :: nx, ny, nxh, nyh, nxe, nye
      integer :: nxye, nxeh, nxyh, nxhy, nx1, ny1, ntmax, ipbc
      integer :: nloop, itime, itime0, ntime, ltime, isign, irc
      integer :: j, it, jt, itw, iur1, iur2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, time = 0.0, tloop = 0.0, ts = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0 
      real :: tfft = 0.0, tfield = 0.0, totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, qi0, etx, we, wke
      real :: vtxi, vtyi, vtdxi, vtdyi
      real, dimension(:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:), pointer :: qe, qi
      real, dimension(:,:,:), pointer :: fxye
      complex, dimension(:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:), pointer :: pt
      integer, dimension(:), pointer :: ip, npic
      real, dimension(:,:), pointer :: qt, sfield
      complex, dimension(:,:), pointer :: dent, pott
      real, dimension(:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! semi-periodic boundary conditions
      integer :: indx1, nxv, nx2v
      real, dimension(:,:), pointer :: q1, sfield1
      real, dimension(:,:,:), pointer :: fxy1
      integer, dimension(:), pointer :: mixup1
      complex, dimension(:), pointer :: sct1
      complex, dimension(:,:), pointer :: ffb
! dirichlet or neumann boundary conditions
      integer :: indy1, ny2
      real, dimension(:,:), pointer :: q2, sfield2
      real, dimension(:,:,:), pointer :: fxy2
      integer, dimension(:), pointer :: mixup2
      complex, dimension(:), pointer :: sct2
      complex, dimension(:,:), pointer :: ffd
! vacuum boundary conditions
      real, dimension(:,:,:), pointer :: ffg
! semi-periodic mixed boundary conditions
      integer :: indx2, nx2, nx4v
      real, dimension(:,:), pointer :: q3, sfield3
      real, dimension(:,:,:), pointer :: fxy3
      integer, dimension(:), pointer :: mixup3
      complex, dimension(:), pointer :: sct3
      complex, dimension(:), pointer :: sct2x
      complex, dimension(:,:), pointer :: ffh
! semi-periodic vacuum boundary conditions
      real :: wb
      integer, dimension(:), pointer :: mixupx, mixupy
      complex, dimension(:), pointer :: sctx, scty
      complex, dimension(:,:), pointer :: ffv
      real, dimension(:,:), pointer :: bv
      real, dimension(:), pointer :: bcd, tx, ty
!
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
      idcode = 4
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
         ax = .912871; ay = .912871
         nxe = nx + 2; nye = ny + 1
      endif
      nxye = nxe*nye; nxeh = nxe/2
! boundary conditions
      ipbc = psolve
      if ((psolve==VACUUM_2D).or.(psolve==VACUUM_3D)) ipbc = 2
      if (psolve==NEUMANN_2D) ipbc = 2
      if ((psolve==NEUMANN_PERIODIC_2D).or.                             &
     &(psolve==DIRICHLET_NEUMANN_PERIODIC_2D)) ipbc = 3
      if (psolve==9) ipbc = 3
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
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
      allocate(part(idimp,np))
! in real space, qe(j,k) = charge density at grid point (j,k)
! in real space, qi(j,k) = ion charge density at grid point (j,k)
      allocate(qe(nxe,nye),qi(nxe,nye))
! in real space, fxye(i,j,k) = i component of force/charge at grid (j,k)
! in other words, fxye are the convolutions of the electric field
! over the particle shape
      allocate(fxye(2,nxe,nye))
! ffc = form factor array for poisson solver
      allocate(ffc(nxeh,nyh))
! mixup, sct = arrays for fft
      allocate(mixup(nxhy),sct(nxyh))
!
! non-periodic boundary conditions
      indx1 = indx + 1; indy1 = indy + 1; indx2 = indx + 2
      nxv = nx + 2; nx2v = 2*nxv; ny2 = 2*ny; nx2 = 2*nx; nx4v = 2*nx2v
! dirichlet conditions
      if (psolve==DIRICHLET_2D) then
         allocate(q2(nx2v,ny2),fxy2(2,nx2v,ny2))
         allocate(sfield2(nx2v,ny2))
         allocate(ffd(nxv,ny))
         allocate(mixup2(2*nxhy))
         allocate(sct2(2*nxyh))
! semi-periodic dirichlet boundary conditions
      else if (psolve==DIRICHLET_PERIODIC_2D) then
         allocate(q1(nx2v,ny),fxy1(2,nx2v,ny))
         allocate(sfield1(nx2v,ny))
         allocate(ffb(nxv,nyh))
         allocate(mixup1(max(nx,ny)),sct1(max(nx,nyh)))
         allocate(sct2(2*nxyh))
! vacuum boundary conditions
      else if ((psolve==VACUUM_2D).or.(psolve==VACUUM_3D)) then
         allocate(q2(nx2v,ny2),fxy2(2,nx2v,ny2))
         allocate(sfield2(nx2v,ny2))
         allocate(ffg(4,nx1,ny1))
         allocate(mixup2(2*nxhy),sct2(2*nxyh))
! neumann boundary conditions
      else if (psolve==NEUMANN_2D) then
         allocate(q2(nx2v,ny2),fxy2(2,nx2v,ny2))
         allocate(sfield2(nx2v,ny2))
         allocate(ffd(nxv,ny))
         allocate(mixup2(2*nxhy),sct2(2*nxyh))
! semi-periodic neumann boundary conditions
      else if (psolve==NEUMANN_PERIODIC_2D) then
         allocate(q1(nx2v,ny),fxy1(2,nx2v,ny))
         allocate(sfield1(nx2v,ny))
         allocate(ffb(nxv,nyh))
         allocate(mixup1(max(nx,ny)),sct1(max(nx,nyh)))
         allocate(sct2(2*nxyh))
! semi-periodic dirichlet-neumann boundary conditions
      else if (psolve==DIRICHLET_NEUMANN_PERIODIC_2D) then
         allocate(q3(nx4v,ny),fxy1(2,nx2v,ny),fxy3(2,nx4v,ny))
         allocate(sfield1(nx2v,ny),sfield3(nx4v,ny))
         allocate(ffh(nxv,nyh))
         allocate(mixup3(max(2*nx,ny)),sct3(max(2*nx,nyh)))
         allocate(sct2(2*nxyh),sct2x(nx))
! semi-periodic vacuum boundary conditions
      else if (psolve==9) then
         allocate(mixupx(nxh),mixupy(nyh),sctx(nxh),scty(nyh))
         allocate(ffv(nxh,nyh),bv(ny,6),bcd(nyh),tx(nx),ty(ny))
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
      if (movion==1) then
         qbmi = qmi/rmass
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtdxi = vtdx/sqrt(rmass*rtempdxi)
         vtdyi = vtdy/sqrt(rmass*rtempdyi)
      endif
! set initial time
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! set default diagnostic file names
      if (ntd > 0) fdname = 'denk2.'//cdrun
      if (ntp > 0) fpname = 'potk2.'//cdrun
! energy time history
      if (ndw > 0) then
         allocate(wt((nloop-1)/ndw-(itime0/ndw)+1,4))
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
! dirichlet boundary conditions
      if (psolve==DIRICHLET_2D) then
         call fst_init(mixup,sct2,indx,indy)
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisd_init(ffd,ax,ay,affp,nx,ny)
! semi-periodic dirichlet boundary conditions
      else if (psolve==DIRICHLET_PERIODIC_2D) then
         call fst_init(mixup,sct2,indx,indy)
         call fft_init(mixup1,sct1,indx1,indy)
         call poism_init(ffb,ax,ay,affp,nx,ny)
! vacuum boundary conditions
      else if (psolve==VACUUM_3D) then
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisc3_init(ffg,q2,mixup2,sct2,ax,affp,indx,indy)
! neumann boundary conditions
      else if (psolve==NEUMANN_2D) then
         call fft_init(mixup2,sct2,indx1,indy1)
         call poisn_init(ffd,ax,ay,affp,nx,ny)
! semi-periodic neumann boundary conditions
      else if (psolve==NEUMANN_PERIODIC_2D) then
!        call fst_init(mixup,sct2,indx,indy)
         call fft_init(mixup1,sct1,indx1,indy)
         call poismn_init(ffb,ax,ay,affp,nx,ny)
! semi-periodic dirichlet-neumann boundary conditions
      else if (psolve==DIRICHLET_NEUMANN_PERIODIC_2D) then
         call fst_init(mixup,sct2,indx,indy)
         call fdt_init(sct2x,indx)
         call fft_init(mixup3,sct3,indx2,indy)
         call poismd_init(ffh,ax,ay,affp,nx,ny)
! semi-periodic vacuum boundary conditions
      else if (psolve==9) then
         call fft_init(mixupy,scty,indy)
         call poisb_init(ffv,bcd,affp,indx,ny)
      endif
!
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
            call vdistr(part,1,npxy,vtx,vty,vx0,vy0)
         endif
! beam electrons
         if (npxyb > 0) then
            call fdistr(part,npxy1,npxyb,ampdx,scaledx,shiftdx,ampdy,   &
     &scaledy,shiftdy,npxb,npyb,nx,ny,ipbc,ndprof,nsrand)
            call vdistr(part,npxy1,npxyb,vtdx,vtdy,vdx,vdy)
         endif
! initialize ions
         if (movion==1) then
            npxy1 = npxyi + 1
! background ions
            if (npxyi > 0) then
               call fdistr(parti,1,npxyi,ampdxi,scaledxi,shiftdxi,ampdyi&
     &,scaledyi,shiftdyi,npxi,npyi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,1,npxyi,vtxi,vtyi,vxi0,vyi0)
            endif
! beam ions
            if (npxybi > 0) then
               call fdistr(parti,npxy1,npxybi,ampdxi,scaledxi,shiftdxi, &
     &ampdyi,scaledyi,shiftdyi,npxbi,npybi,nx,ny,ipbc,ndprofi,nsrandi)
               call vdistr(parti,npxy1,npxybi,vtdxi,vtdyi,vdxi,vdyi)
            endif
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
! initialize ion charge density to zero
            call sguardp(qi,zero,nx,ny,ipbc,inorder)
! deposit ion charge
            call dpost(parti,qi,npi,qmi,tdposti,inorder,dopt)
! add guard cells for ion charge density
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
     &qi,irc,iuer)
            if (irc /= 0) cycle
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
     &nprec,fpname,irc,iuer)
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
! initialize charge density to background
      call sguardp(qe,zero,nx,ny,ipbc,inorder)
! deposit charge
      call dpost(part,qe,np,qme,tdpost,inorder,dopt)
! add guard cells for charge density
      call aguardp(qe,nx,ny,ipbc,inorder)
! add ion density
      if (movion==1) then
         call sguardp(qi,zero,nx,ny,ipbc,inorder)
         call dpost(parti,qi,npi,qmi,tdposti,inorder,dopt)
         call aguardp(qi,nx,ny,ipbc,inorder)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti)
            movion = 0
         endif
      endif
! add electron and ion densities
      call addqei(qe,qi,nx,ny,inorder)
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
!                 call bfopen(dent,modesxd,modesy2d,iud,ndrec,          &
!    &trim(fdname))
!              endif
!              call writebf(dent,modesxd,modesy2d,iud,ndrec,            &
!    &order=LINEAR)
!           endif
! transform ion density to real space
            if (jt==0) then
               isign = 1
!              call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,     &
!    &LINEAR)
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
!                 call bfopen(pott,modesxp,modesy2p,iup,nprec,          &
!    &trim(fpname))
!              endif
!              call writebf(pott,modesxp,modesy2p,iup,nprec,            &
!    &order=LINEAR)
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
! calculate force/charge in fourier space
!     call poisdx(q2,fxy2,ffd,we,nx,ny)
      call poisd(qe,fxye,ffd,we,nx,ny,inorder)
      isign = 1
      call fcst(fxye,isign,mixup,sct2,tfft,indx,indy,inorder)
!     call fft(fxy2,mixup2,sct2,tfft,indx1,indy1,LINEAR)
!     call hafdbl(fxye,fxy2,nx,ny,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
!
! semi-periodic dirichlet boundary conditions
!
      else if (psolve==DIRICHLET_PERIODIC_2D) then
! ion density diagnostic
      if ((ntd > 0) .or. (ndd > 0)) then
         it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
         jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
         if ((it==0) .or. (jt==0)) then
            qt = qi
! transform ion density to fourier space
!           call sglsin(qt,q1,nx,ny,inorder)
            isign = -1
!           call fft(q1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
            call fsft(qt,isign,mixup,sct2,tfft,indx,indy,inorder)
! calculate smoothing in fourier space
!           call poismx(q1,sfield1,ffb,nx,ny)
            call poism(qt,sfield,ffb,nx,ny,inorder)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(sfield,dent,nx,ny,modesxd,modesyd,inorder)
! write diagnostic output
!              if (ndrec==0) then
!                 ndrec = -1; iud = get_funit(12)
!                 call bfopen(dent,modesxd,modesy2d,iud,ndrec,          &
!    &trim(fdname))
!              endif
!              call writebf(dent,modesxd,modesy2d,iud,ndrec,            &
!    &order=LINEAR)
!           endif
! transform ion density to real space
            if (jt==0) then
               isign = 1
!              call fft(sfield1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
!              call hafsgl(sfield,sfield1,nx,ny,inorder)
               call fsft(sfield,isign,mixup,sct2,tfft,indx,indy,inorder)
               call cguardp(sfield,nx,ny,ipbc,inorder)
! display ion density
               call displays(sfield,' I DENSITY',ntime,999,2,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! transform charge to fourier space
!     call sglsin(qe,q1,nx,ny,inorder)
      isign = -1
!     call fft(q1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
      call fsft(qe,isign,mixup,sct2,tfft,indx,indy,inorder)
! potential diagnostic
      if ((ntp > 0) .or. (ndp > 0)) then
         it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
         jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
         if ((it==0) .or. (jt==0)) then
! solve for potential
!           call poismx(q1,sfield1,ffb,we,nx,ny)
            call poism(qe,sfield,ffb,we,nx,ny,inorder)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(sfield,pott,nx,ny,modesxp,modesyp,inorder)
! write diagnostic output
!              if (nprec==0) then
!                 nprec = -1; iup = get_funit(11)
!                 call bfopen(pott,modesxp,modesy2p,iup,nprec,          &
!    &trim(fpname))
!              endif
!              call writebf(pott,modesxp,modesy2p,iup,nprec,            &
!    &order=LINEAR)
!           endif
! transform potential to real space
            if (jt==0) then
               isign = 1
!              call fft(sfield1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
!              call hafsgl(sfield,sfield1,nx,ny,inorder)
               call fsft(sfield,isign,mixup,sct2,tfft,indx,indy,inorder)
               call cguardp(sfield,nx,ny,ipbc,inorder)
! display potential
               call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! calculate force/charge in fourier space
!     call poismx(q1,fxy1,ffb,we,nx,ny)
      call poism(qe,fxye,ffb,we,nx,ny,inorder)
      isign = 1
      call fcsft(fxye,isign,mixup,sct2,tfft,indx,indy,inorder)
!     call fft(fxy1,mixup1,sct1,tfft,indx1,indy,LINEAR)
!     call hafsgl(fxye,fxy1,nx,ny,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
!
! vacuum boundary conditions
!
      else if (psolve==VACUUM_3D) then
! ion density diagnostic
      if ((ntd > 0) .or. (ndd > 0)) then
         it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
         jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
         if ((it==0) .or. (jt==0)) then
            qt = qi
! transform ion density to fourier space
            call zdbl(qt,q2,nx,ny,inorder)
            isign = -1
            call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! calculate smoothing in fourier space
            call poisc(q2,sfield2,ffg,nx,ny)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(sfield,dent,nx,ny,modesxd,modesyd,inorder)
! write diagnostic output
!              if (ndrec==0) then
!                 ndrec = -1; iud = get_funit(12)
!                 call bfopen(dent,modesxd,modesy2d,iud,ndrec,          &
!    &trim(fdname))
!              endif
!              call writebf(dent,modesxd,modesy2d,iud,ndrec,            &
!    &order=LINEAR)
!           endif
! transform ion density to real space
            if (jt==0) then
               isign = 1
               call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,     &
     &LINEAR)
               call hafdbl(sfield,sfield2,nx,ny,inorder)
               call cguardp(sfield,nx,ny,ipbc,inorder)
! display ion density
               call displays(sfield,' I DENSITY',ntime,999,2,ndstyle,nx,&
     &ny,irc,inorder)
            if (irc==1) go to 2000
            endif
         endif
      endif
! transform charge to fourier space
      call zdbl(qe,q2,nx,ny,inorder)
      isign = -1
      call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! potential diagnostic
      if ((ntp > 0) .or. (ndp > 0)) then
         it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
         jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
         if ((it==0) .or. (jt==0)) then
! solve for potential
            call poisc(q2,sfield2,ffg,we,nx,ny)
! store selected fourier modes
!           if (it==0) then
!              call gtmodes(sfield,pott,nx,ny,modesxp,modesyp,inorder)
! write diagnostic output
!              if (nprec==0) then
!                 nprec = -1; iup = get_funit(11)
!                 call bfopen(pott,modesxp,modesy2p,iup,nprec,          &
!    &trim(fpname))
!              endif
!              call writebf(pott,modesxp,modesy2p,iup,nprec,            &
!    &order=LINEAR)
!           endif
            if (jt==0) then
               isign = 1
               call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,     &
     &LINEAR)
               call hafdbl(sfield,sfield2,nx,ny,inorder)
               call cguardp(sfield,nx,ny,ipbc,inorder)
! display potential
               call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,&
     &ny,irc,inorder)
               if (irc==1) go to 2000
            endif
         endif
      endif
! calculate force/charge in fourier space
      call poisc(q2,fxy2,ffg,we,nx,ny)
      call fft(fxy2,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call hafdbl(fxye,fxy2,nx,ny,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
!
! semi-periodic vacuum boundary conditions
!
      else if (psolve==9) then
! ion density diagnostic
      if (ntd > 0) then
         it = ntime/ntd
         if (ntime==ntd*it) then
            qt = qi
! transform electron density to fourier space
            isign = -1
            call fft(qt,isign,mixup,sct,tfft,indx,indy,inorder)
! calculate smoothing in fourier space
            sfield = 0.0
            call spois(qt,sfield,ffc,nx,ny,inorder)
! transform electron density to real space
            isign = 1
            call fft(sfield,isign,mixup,sct,tfft,indx,indy,inorder)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display electron density
            call displays(sfield,' E DENSITY',ntime,999,2,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,isign,mixup,sct,tfft,indx,indy,inorder)
      bv(:,1) = 0.
      bv(:,2) = 0.
      call bndryv(qe,ffc,bv,nx,ny,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! calculate potential in fourier space
            call pois(qe,sfield,ffc,we,nx,ny,inorder)
            call poisb(sfield,ffv,bv,bcd,wb,affp,indx,ny,inorder)
! transform potential to real space
            isign = 1
            call fft(sfield,isign,mixup,sct,tfft,indx,indy,inorder)
            call fft(bv(:,6),isign,mixupy,scty,tfft,indy,LINEAR)
            call ncguard(sfield,bv,nx,ny,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call pois(qe,fxye,ffc,tfield,we,nx,ny,inorder)
! solve for bounded part of electric fields
      call poisb(fxye,ffv,bv,bcd,wb,affp,indx,ny,inorder)
      we = we + wb
! transform force/charge to real space
      isign = 1
      call fft(fxye,mixup,sct,tfft,indx,indy,inorder)
      call fft(bv(:,5),isign,mixupy,scty,tfft,indy,LINEAR)
      call fft(bv(:,6),isign,mixupy,scty,tfft,indy,LINEAR)
      call ncguard(fxye,bv,nx,ny,inorder)
!
! neumann boundary conditions
!
      else if (psolve==NEUMANN_2D) then
! ion density diagnostic
      if (ntd > 0) then
         it = ntime/ntd
         if (ntime==ntd*it) then
            qt = qi
! transform electron density to fourier space
            call dblcos(qt,q2,nx,ny,inorder)
            isign = -1
            call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! calculate smoothing in fourier space
            call poisn(q2,sfield2,ffd,nx,ny)
! transform electron density to real space
            isign = 1
            call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,inorder)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display electron density
            call displays(sfield,' E DENSITY',ntime,999,2,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      call dblcos(qe,q2,nx,ny,inorder)
      isign = -1
      call fft(q2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! solve for potential
            call poisn(q2,sfield2,ffd,we,nx,ny)
            isign = 1
            call fft(sfield2,isign,mixup2,sct2,tfft,indx1,indy1,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,inorder)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poisn(q2,fxy2,ffd,we,nx,ny)
      call fft(fxy2,mixup2,sct2,tfft,indx1,indy1,LINEAR)
      call hafdbl(fxye,fxy2,nx,ny,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
!
! semi-periodic neumann boundary conditions
!
      else if (psolve==NEUMANN_PERIODIC_2D) then
! ion density diagnostic
      if (ntd > 0) then
         it = ntime/ntd
         if (ntime==ntd*it) then
            qt = qi
! transform electron density to fourier space
            call sglcos(qt,q1,nx,ny,inorder)
            isign = -1
            call fft(q1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
!           call fcft(qt,isign,mixup,sct2,tfft,indx,indy,inorder)
! calculate smoothing in fourier space
            call poismn(q1,sfield1,ffb,nx,ny)
!           call poismn(qt,sfield,ffb,nx,ny,inorder)
! transform electron density to real space
            isign = 1
            call fft(sfield1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
            call hafsgl(sfield,sfield1,nx,ny,inorder)
!           call fcft(sfield,isign,mixup,sct2,tfft,indx,indy,inorder)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display electron density
            call displays(sfield,' E DENSITY',ntime,999,2,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
      call sglcos(qe,q1,nx,ny,inorder)
      isign = -1
      call fft(q1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
!     call fcft(qe,isign,mixup,sct2,tfft,indx,indy,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! solve for potential
            call poismn(q1,sfield1,ffb,we,nx,ny)
!           call poismn(qe,sfield,ffb,we,nx,ny,inorder)
            isign = 1
            call fft(sfield1,isign,mixup1,sct1,tfft,indx1,indy,LINEAR)
            call hafsgl(sfield,sfield1,nx,ny,inorder)
!           call fcft(sfield,isign,mixup,sct2,tfft,indx,indy,inorder)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poismn(q1,fxy1,ffb,we,nx,ny)
!     call poismn(qe,fxye,ffb,we,nx,ny,inorder)
      isign = 1
!     call fscft(fxye,isign,mixup,sct2,tfft,indx,indy,inorder)
      call fft(fxy1,mixup1,sct1,tfft,indx1,indy,LINEAR)
      call hafsgl(fxye,fxy1,nx,ny,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
!
! semi-periodic dirichlet-neumann boundary conditions
!
      else if (psolve==DIRICHLET_NEUMANN_PERIODIC_2D) then
! ion density diagnostic
      if (ntd > 0) then
         it = ntime/ntd
         if (ntime==ntd*it) then
            qt = qi
! transform electron density to fourier space
!           call sgldsin(qt,q3,nx,ny,inorder)
            isign = -1
!           call fft(q3,isign,mixup3,sct3,tfft,indx2,indy,LINEAR)
            call fdsft(qt,isign,mixup,sct2,sct2x,tfft,indx,indy,inorder)
! calculate smoothing in fourier space
!           call poismdx(q3,sfield3,ffh,nx,ny)
            call poismd(qt,sfield,ffh,nx,ny,inorder)
! transform electron density to real space
            isign = 1
!           call fft(sfield3,isign,mixup3,sct3,tfft,indx2,indy,LINEAR)
!           call hafsgl(sfield1,sfield3,nx2,ny,LINEAR)
!           call hafsgl(sfield,sfield1,nx,ny,inorder)
            call fdsft(sfield,isign,mixup,sct2,sct2x,tfft,indx,indy,inor&
     &der)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display electron density
            call displays(sfield,' E DENSITY',ntime,999,2,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
         endif
      endif
! transform charge to fourier space
!     call sgldsin(qe,q3,nx,ny,inorder)
      isign = -1
!     call fft(q3,isign,mixup3,sct3,tfft,indx2,indy,LINEAR)
      call fdsft(qe,isign,mixup,sct2,sct2x,tfft,indx,indy,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = ntime/ntp
         if (ntime==ntp*it) then
! solve for potential
!           call poismdx(q3,sfield3,ffh,we,nx,ny)
            call poismd(qe,sfield,ffh,we,nx,ny,inorder)
            isign = 1
!           call fft(sfield3,isign,mixup3,sct3,tfft,indx2,indy,LINEAR)
!           call hafsgl(sfield1,sfield3,nx2,ny,LINEAR)
!           call hafsgl(sfield,sfield1,nx,ny,inorder)
            call fdsft(sfield,isign,mixup,sct2,sct2x,tfft,indx,indy,inor&
     &der)
            call cguardp(sfield,nx,ny,ipbc,inorder)
! display potential
            call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,nx,ny,&
     &irc,inorder)
            if (irc==1) go to 2000
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
!     call poismdx(q3,fxy3,ffh,we,nx,ny)
      call poismd(qe,fxye,ffh,we,nx,ny,inorder)
      isign = 1
      call fdcsft(fxye,isign,mixup,sct2,sct2x,tfft,indx,indy,inorder)
!     call fft(fxy3,mixup3,sct3,tfft,indx2,indy,LINEAR)
!     call hafsgl(fxy1,fxy3,nx2,ny,LINEAR)
!     call hafsgl(fxye,fxy1,nx,ny,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
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
! calculate force/charge in fourier space
      call pois(qe,fxye,ffc,we,tfield,nx,ny,inorder)
! transform force/charge to real space
      call fft(fxye,mixup,sct,tfft,indx,indy,inorder)
      call cguardp(fxye,nx,ny,ipbc,inorder)
!
      endif
!
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxye(1,:,:) = fxye(1,:,:) + etx
!     endif
! particle push and charge density update
      wke = 0.
! push electrons
      call pushg(part,fxye,np,qbme,dt,ci,wke,tpush,nx,ny,ipbc,relativity&
     &,inorder,popt)
! debug: zero force
!     call pushzfg(part,np,dt,ci,wke,tpush,nx,ny,ipbc,relativity)
! push ions
      if (movion==1) then
         wki = 0.
         call pushg(parti,fxye,npi,qbmi,dt,ci,wki,tpushi,nx,ny,ipbc,    &
     &relativity,inorder,popt)
! debug: zero force
!        call pushzfg(parti,npi,dt,ci,wki,tpushi,nx,ny,ipbc,relativity)
         wki = wki*rmass
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
      call esenergy(wt,we,wke,wki,ntw,ndw,itw,iuot,ntime)
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
! write primary file
            call restart_bwrite(iur1,itime,itime0,np,part,movion,npi,   &
     &parti,qi)
            call restart_dwrite(iur1,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname)
! write secondary file
            call restart_bwrite(iur2,itime,itime0,np,part,movion,npi,   &
     &parti,qi)
            call restart_dwrite(iur2,itime,itw,wt,ndrec,fdname,nprec,   &
     &fpname)
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
      write (iuot,*) 'bounded electrostatic code d0_beps2'
      write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
      write (iuot,*) 'main wall clock time = ', time, 'sec'
      totpush = tpush + tdpost
      write (iuot,*) 'electron push time = ', tpush, 'sec'
      write (iuot,*) 'electron charge deposit time = ', tdpost, 'sec'
      write (iuot,*) 'total electron push time = ', totpush, 'sec'
      write (iuot,*) 'electron sort time = ', tsort
      totpush = totpush + tsort
      write (iuot,*) 'total electron time = ', totpush, 'sec'
      if (movion==1) then
         totpushi = tpushi + tdposti
         write (iuot,*) 'ion push time = ', tpushi, 'sec'
         write (iuot,*) 'ion charge deposit time = ', tdposti, 'sec'
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
      end program d0_beps2
