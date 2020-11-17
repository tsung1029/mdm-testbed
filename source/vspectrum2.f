!-----------------------------------------------------------------------
! this program finds spectrum and autocorrelation for 2d periodic data
! interactively modifies input parameters
! written for macintosh - viktor k. decyk, ucla
! copyright 1990, regents of the university of california
! character*8 code(nc), where nc is the number of input parameters
! character*8 cp(ncc), where ncc is the number of character parameters
! dimension ip(nci), ap(ncr)
! where nci, ncr are the number of integer and real parameters
! update: march 2, 2012
      program vspectrum2
      use globals, only: LINEAR, QUADRATIC
      use fft2d
      use field2d
      use diag2d
      implicit none
      integer, parameter :: ncc=0, nci=11, ncr=3, nc=ncc+nci+ncr, ns=4
      integer :: idpal = 1
      integer :: nda, ndj, nde, nrec, modesx, modesy, modesxy
      integer :: lts, its, nts
      integer :: kxmax, kxmin, kymax, kymin, kxmx, kxmn, kymx, kymn
      integer :: ntq, ntc, nplot, nvf, istyle, itwo, iwmax, kymn0
      integer :: ltime, modesy2, modesxy2
      integer :: nx, nxh, ny, nyh, nxe, nye, nxeh, nxyh, nxhy, nx1, ny1
      integer :: it1, mta, nt2d, it, nt, nt2, inft, nft, nfth, isign
      integer :: i, j, k, n, ii, i1, j1, k1, irc, nts1, nts2, ntc1, ntc2
      integer :: j0, k0, j2, k2, jk, kk, iw, iw1, iw2, ncd, ntr
      integer :: inx, iny, nxd, nyd, nxdh, nvpot
      integer :: nf, nodesx, nodesy, nodesy2, mt
      real :: dti = 0.0
      real :: wmin, wmax, dw, pmin, anorm, bnorm
      real :: at1, at2, dtp, dnx, dny, dts, tt0, times, dkx, dky, ak, wk
      real :: ak2, dti2c, tfft
      double precision :: sum1, sum2
      complex :: zt1, zt2, zt3, zt4, zsc, zst
      logical :: ex
      integer, dimension(1) :: itm
      integer, dimension(nci) :: ip
      real, dimension(ncr) :: ap
      real, dimension(:,:), pointer :: vpotb, vpots, vpote, vpotr, vpk
      real, dimension(:,:,:), pointer :: vpotd
      complex, dimension(:,:,:), pointer :: vpotx, vpotp, vpoty
      real, dimension(:,:,:), pointer :: vpwk, vpck, vpak
      real, dimension(:), pointer :: potc, time
      complex, dimension(:,:,:,:), pointer ::  vpott
      complex, dimension(:,:), pointer ::  vpotts
      real, dimension(:), pointer :: g
      integer, dimension(:), pointer :: mixup, mixupt
      complex, dimension(:), pointer :: sct, sctt
      real, dimension(:), pointer :: wm, p, pc, ps, pcs
!
      character(len=60) :: prompt
      character(len=32) :: fname, ftname
      character(len=10) :: cdrun0
      character(len=20) :: cdrun, runid
      character(len=64) :: label, chr
      character(len=16), dimension(6) :: cnvf
      character(len=18), dimension(6) :: crnvf
      character(len=8), dimension(6) :: snvf, srnvf
      character(len=10), dimension(4) :: chrs, chws
      character(len=8), dimension(nc) :: code
      character(len=8), dimension(1) :: cp
!
      integer, external :: NDIAN, NDPREC
!
      equivalence (lts, ip(1)), (its, ip(2)), (nts, ip(3))
      equivalence (kxmin, ip(4)), (kxmax, ip(5))
      equivalence (kymin, ip(6)), (kymax, ip(7))
      equivalence (ntq, ip(8)), (ntc, ip(9)), (nplot, ip(10))
      equivalence (nvf, ip(11))
      equivalence (wmin, ap(1)), (wmax, ap(2)), (dw, ap(3))
   91 format (' SPECTRUM ANALYSIS FOR 2D DATA: ',A32)
   92 format (' RUNID= ',a20,' INDX=',i3,' INDY=',i3)
   93 format (' NTA=',i6,' MODESX=',i5,' MODESY=',i5,' PSOLVE=',i2)
   94 format (' NREC=',i10)
   95 format (' T0=',f8.1,' TEND=',f8.1,' DT=',f8.5)
  981 format (' ACTUAL LENGTH OF DATA = ',i6,' DATA EXPECTED = ',i6)
  982 format (' ERROR IN PARAMETERS, LTS,ITS,NTS,NT = ',4i7)
  983 format (' ERROR IN MODE PARAMETERS, KXMIN,KXMAX,MODESX = ',3i7)
  984 format (' ERROR IN MODE PARAMETERS, KYMIN,KYMAX,MODESY = ',3i7)
  985 format (' ERROR IN DISPLAY PARAMETERS, NTC,NTD,NTS = ',3i7)
  986 format (' ERROR IN PARAMETER, NVF = ',i7)
  987 format (' ERROR IN FREQUENCY PARAMETERS, IW,IWMAX = ',2i7)
  891 format (' KX, KY, AVERAGE MAGNETIC ENERGY')
  892 format (1x,2i4,1x,e14.7)
  893 format ('  PEAK FREQUENCY OF RAW DATA')
  894 format ('  PEAK FREQUENCY OF CORRRELATED DATA')
  895 format (' N ',i1,', KX ',i4,', KY ',i4,', AKX=',f7.4,', AKY=',f7.4&
     &,', AK=',f7.4)
  896 format (' MAX P(',f6.3,') = ',e14.7,' MAX P(',f7.3,') = ',e14.7)
!
  991 format (' <ENERGY> = ',e14.7,',')
  992 format (', SUM = ',e14.7)
  993 format (' INTEGRATED SPECTRUM VERSUS OMEGA')
  994 format (' INTEGRATED CORRELATION SPECTRUM VS OMEGA')
  995 format (1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)
  996 format (' N ',i1,', KX ',i4,', KY ',i4,', AKX=',f7.4,', AKY=',f7.4&
     &)
! 994 format ('  POWER SPECTRUM:',6x,' RAW DATA')
! 995 format (' P(',f6.3,') = ',e14.7,' P(',f6.3,') = ',e14.7)
! 996 format ('  POWER SPECTRUM:',6x,' RAW DATA',50x,' CORRELATED DATA')
! 997 format (' P(',f6.3,') = ',e14.7,' P(',f6.3,') = ',e14.7,9x,' P(',f&
!    &6.3,') = ',e14.7,' P(',f6.3,') = ',e14.7)
! istyle = (0,1) = use (brief,full) menu style
      data istyle /1/
      data cnvf /' VPOTENTIAL     ',' ELECTRIC FIELD ',' MAGNETIC FIELD &
     &',' JPERP        ',' EDOT           ',' JDOT         '/
      data crnvf /' VRPOTENTIAL      ',' RADIATIVE EFIELD ',' RADIATIVE &
     &BFIELD ',' JRPERP           ',' ERDOT            ',' JRDOT        &
     &    '/
      data snvf /' VPOT   ',' EFIELD ',' BFIELD ',' JPERP  ',' EDOT   ',&
     &' JDOT   '/
      data srnvf /' VRPOT  ',' ERFIELD',' BRFIELD',' JRPERP ',' ERDOT  '&
     &,' JRDOT  '/
      data chrs /'   REAL   ','IMAGINARY ',' POSITIVE ',' NEGATIVE '/
      data chws /' KY>0, W>0',' KY<0, W>0',' KY>0, W<0',' KY<0, W<0'/
      data code /'LTS   ','ITS   ','NTS   ','KXMIN ','KXMAX ','KYMIN ','&
     &KYMAX ','NTD   ','NTC   ','NPLOT ','NVF  ','WMIN  ','WMAX  ','DW  &
     &  '/
      data itwo /2/
      data pmin /1.0e-14/
! lts,its = initial data point and increment between
      data lts,its /1,1/
      data iwmax /5001/
! wmin,wmax,dw = initial, final frequency and increment
      data wmin,wmax,dw /0.,2.0,.01/
! nplot = number of plots per page
      data nplot /4/
! nvf = vector field (a=1,e=2,b=3,j=4,edot=5,jdot=6)
      data nvf /1/
! open graphics device
      call GROPEN0(1)
! get file name
      label = ' '
!   5 prompt = 'enter name of diagnostic metafile, q to quit:'
    5 prompt = 'enter runid, or q to quit:'
      call GTINPUT(label,prompt,chr,irc)
      if ((chr=='q').or.(chr=='Q')) go to 210
      nf = 1
      i = index(chr,'/')
      if (i > 1) nf = 2
!     fname = chr(i+1:)
      fname = 'diag2.'//chr(i+1:)
! open diagnostic metafile
      do j = 1, nf
         if (j==2) then
            write (cdrun0,'(i10)') idrun
            cdrun0 = adjustl(cdrun0)
            if (nda==0) then
               nodesx = modesxa; nodesy = modesya
               ftname = faname
            else if (ndj==0) then
               nodesx = modesxj; nodesy = modesyj
               ftname = fjname; cdrun0 = trim(cdrun0)//'j'
            else if (nde==0) then
               nodesx = modesxe; nodesy = modesye
               ftname = fename; cdrun0 = trim(cdrun0)//'r'
            endif
            nxh = 2**(indx-1); nyh = max(1,2**(indy-1))
            if (nodesx > nxh) nodesx = nxh
            if (nodesy > nyh) nodesy = nyh
            nodesy2 = 2*nodesy - 1
            allocate(vpoty(ndim,nodesx,nodesy2))
!           fname = chr(1:i-1)
            fname = 'diag2.'//chr(1:i-1)
            close(unit=19)
            chr = ftname
         endif
         open(unit=19,file=trim(fname),form='formatted',status='old',ios&
     &tat=irc)
         if (irc /= 0) then
            label = 'open error for '//trim(fname)
            go to 5
         endif
         read (19,vpot2d,iostat=nda)
         rewind 19
         read (19,vcur2d,iostat=ndj)
         rewind 19
         read (19,em2d,iostat=nde)
         rewind 19
         if (((nda==0).and.(ndj==0)).or.((nda==0).and.(nde==0)).or.((ndj&
     &==0).and.(nde==0))) then
            label = ' '
            if (nf==2) then
               label = 'for second file'
               if (j==2) label = 'for first file'
            endif
            if ((nda==0).and.(ndj==0).and.(nde==0)) then
               prompt = 'enter a, j or e for total, ion or radiative vpo&
     &t diagnostic:'
            else if ((nda==0).and.(ndj==0)) then
               prompt = 'enter a or j for total or ion vpot diagnostic:'
            else if ((nda==0).and.(nde==0)) then
               prompt = 'enter a or e for total or radiative vpot diagno&
     &stic:'
            else if ((ndj==0).and.(nde==0)) then
               prompt = 'enter j or e for ion or radiative vpot diagnost&
     &ic:'
            endif
    6       call GTINPUT(label,prompt,ftname,irc)
            if ((ftname=='q').or.(ftname=='Q')) then
               close(unit=19)
               label = ' '
               go to 5
            endif
            if (((ftname=='a').or.(ftname=='A')).and.(nda==0)) then
               ndj = -1; nde = -1
            else if (((ftname=='j').or.(ftname=='J')).and.(ndj==0)) then
               nda = -1; nde = -1
            else if (((ftname=='e').or.(ftname=='E')).and.(nde==0)) then
               nda = -1; ndj = -1
            else
               label = 'Invalid string'
               go to 6
            endif
! check file and number formats
            if (indian /= ndian()) then
               label = 'Endian mismatch for '//trim(fname)
               go to 5
            else if (rlprec /= NDPREC()) then
               label = 'Default real mismatch for '//trim(fname)
               go to 5
            endif
         else if ((nda /= 0).and.(ndj /= 0).and.(nde /= 0)) then
            label = 'No valid Namelists found'
            close(unit=19)
            go to 5
         endif
      enddo
! debug
!     ceng = 0.01
! end debug
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
      if (nda==0) then
         read (19,vpot2d,iostat=nda)
         modesx = modesxa; modesy = modesya
         nrec = narec
         fname = faname; ftname = 'vpot2out.'//cdrun
         bnorm = float((2**indx)*(2**indy))/ceng
      else if (ndj==0) then
         read (19,vcur2d,iostat=ndj)
         modesx = modesxj; modesy = modesyj
         nrec = njrec; nta = ntj
         cdrun = trim(cdrun)//'j'
         fname = fjname; ftname = 'vcur2out.'//cdrun
         bnorm = 1.0
      else if (nde==0) then
         read (19,em2d,iostat=nde)
         modesx = modesxe; modesy = modesye
         nrec = nerec; nta = nte
         cdrun = trim(cdrun)//'r'
         fname = fename; ftname = 'vpotr2out.'//cdrun
         bnorm = float((2**indx)*(2**indy))/ceng
      endif
      modesxy = modesx*modesy
      if (nf==2) cdrun = trim(cdrun)//'-'//trim(cdrun0)
      runid = 'vspectrum.'//cdrun
! text output file
      open(unit=18,file=trim(ftname),form='formatted',status='replace')
! echo input
      write (18,91) trim(adjustl(fname))
      write (18,92) runid, indx, indy
      write (18,93) nta, modesx, modesy, psolve
      write (18,94) nrec
      write (18,95) t0, tend, dt
! calculate initial parameters
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nyh = max(1,ny/2)
      nxe = nx + 2
      nye = ny + 1
      nxeh = nxe/2
      nx1 = nx + 1
      ny1 = ny + 1
      if (modesx > nxh) modesx = nxh
      if (modesy > nyh) modesy = nyh
      kxmin = 0
      kxmax = modesx - 1
      kymin = 0
      kymax = modesy - 1
      modesy2 = 2*modesy - 1
      modesxy2 = modesx*modesy2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dtp = dt*real(nta)
      if (dtp > 0) dti = 1.0/dtp
      dti2c = (dti*ci)**2
! time parameters
      if (nrec > 0) then
         mta = nrec
      else
         it1 = (tend - t0)/dt + .0001
         mta = (it1 - 1)/nta + 1
      endif
      nt2d = mta + mta
! allocate space data
      allocate(vpotx(ndim,modesx,modesy2),vpotp(ndim,modesx,modesy2))
      allocate(vpotb(modesx,modesy2),vpots(modesx,modesy2))
      allocate(vpote(modesx,modesy2),vpotr(modesx,modesy2))
      allocate(vpk(modesx,modesy2))
      allocate(vpwk(ndim,modesxy,ns),vpck(ndim,modesxy,ns))
      allocate(vpak(ndim,modesxy,ns))
! allocate time data
      allocate(potc(nt2d),time(mta))
      allocate(vpott(mta,ndim,modesx,modesy2),stat=nvpot)
      if (nvpot/=0) then
         write (2,*) 'attempting to allocate in sections'
         allocate(vpotts(mta,ndim),vpott(mta,ndim,modesx,1),stat=nvpot)
         if (nvpot==0) then
            nvpot = 1
         else
            nvpot = 2
         endif
      endif
! allocate data for non-periodic boundary conditions
      if (psolve==1) then
         inx = indx; iny = indy
         nxd = nx; nyd = ny + 1
      else if (psolve==2) then
         inx = indx+1; iny = indy+1
         nxd = nx + nx; nyd = ny + ny
         dnx = .5*dnx; dny = .5*dny
         kxmin = 1; kymin = 1
      else if (psolve==3) then
         inx = indx+1; iny = indy
         nxd = nx + nx; nyd = ny + 1
         dnx = .5*dnx
         kxmin = 1
      endif
      nxdh = nxd/2
! allocate spatial data
      allocate(vpotd(ndim,nxd,nyd))
! allocate data for ffts
      nxyh = max(nxd,nyd)/2
      nxhy = max(nxdh,nyd)
      allocate(mixup(nxhy),sct(nxyh))
! prepare fft tables
      call fft_init(mixup,sct,inx,iny)
! enter raw data display parameter
      label = ' '
      prompt = 'enter N to display raw data every N steps (0 for none)'
      call GTINPUT(label,prompt,ftname,irc)
      if (irc==1) go to 5
      read (ftname,*,iostat=irc) ntr
      if (irc /= 0) ntr = 1
! clear screen
      call CLRSCRN
      if (ntr > 0) then
! set number of plots per page
         call SETNPLT(4,irc)
         call STPALIT(idpal)
      endif
! get wave numbers
      call SAK2(vpk,nx,ny,modesx,modesy,modesx,modesy2)
      zsc = cmplx(999.0,999.0)
      zst = cmplx(2.0,2.0)
! debug
      zsc = cmplx(1.0,999.0)
      zst = cmplx(1.0,1.0)
! end debug
      if (nda==0) then
         label = ' MAGNETIC ENERGY VERSUS |K|'
         prompt = ' TRANSVERSE ELECTRIC ENERGY VERSUS |K|'
      else if (ndj==0) then
         label = ' ION CURRENT**2 VERSUS |K|'
      else if (nde==0) then
         label = ' RADIATIVE MAGNETIC ENERGY VERSUS |K|'
         prompt = ' TRANSVERSE RADIATIVE ELECTRIC ENERGY VERSUS |K|'
      endif
! read diagnostic input
      nt = 0
      call bfopen(vpotx,modesx,modesy2,15,nt,trim(fname))
      if (nt /= 1) go to 20
      if (nf==2) then
         mt = 0
         call bfopen(vpoty,nodesx,nodesy2,16,mt,trim(chr))
         if (mt /= 1) go to 20
      endif
   10 call readbf(vpotx,modesx,modesy2,15,nt,irc,order=LINEAR)
      if (irc /= 0) go to 20
      if (nf==2) then
         call readbf(vpoty,nodesx,nodesy2,16,mt,irc,order=LINEAR)
         if (irc /= 0) go to 20
         vpotx = vpotx - vpoty
      endif
      it = nt - 1
      if (nvpot==0) vpott(it,:,:,:) = vpotx
! create time array
      time(it) = t0 + dtp*real(it - 1)
! calculate magnetic energy
      if ((nda==0).or.(nde==0)) then
         if (ndim==3) then
            call WBR2(vpotx,vpotb,ceng,ci,potc(it),nx,ny,modesx,modesy,m&
     &odesx,modesy2)
         else if (ndim==2) then
            call WBR22(vpotx,vpotb,ceng,ci,potc(it),nx,ny,modesx,modesy,&
     &modesx,modesy2)
         endif
! calculate average of magnetic energy
         if (it==1) then
            sum1 = potc(it)
            vpots = vpotb
         else
            sum1 = sum1 + potc(it)
            vpots = vpots + vpotb
         endif
! calculate transverse electric energy
         if (it==1) vpotp = vpotx
         vpotp = (vpotp - vpotx)*dti
         if (ndim==3) then
            call WER2(vpotp,vpote,ceng,potc(mta+it),nx,ny,modesx,modesy,&
     &modesx,modesy2)
         else if (ndim==2) then
            call WER22(vpotp,vpote,ceng,potc(mta+it),nx,ny,modesx,modesy&
     &,modesx,modesy2)
         endif
! calculate average of transverse electric energy
         if (it==1) then
            sum2 = potc(mta+it)
            vpotr = vpote
        else
            sum2 = sum2 + potc(mta+it)
            vpotr = vpotr + vpote
         endif
      else if (ndj==0) then
         if (ndim==3) then
            call SQ2VPOT2(vpotx,vpotb,nx,ny,modesx,modesy,modesx,modesy2&
     &)
         else if (ndim==2) then
            call SQ2VPOT22(vpotx,vpotb,nx,ny,modesx,modesy,modesx,modesy&
     &2)
         endif
      endif
! display raw data
      if (ntr > 0) then
      it1 = (it-1)/ntr
      if ((it-1)==ntr*it1) then
         write (chr,*) 'NT = ', it
! display absolute value of vector potential
         call ptmodes(vpotd,vpotx,nx,ny,modesx,modesy,order=LINEAR)
         isign = 1
         call fft(vpotd,isign,mixup,sct,tfft,inx,iny,order=LINEAR)
         if (nda==0) then
            call displayv(vpotd,' VPOTENTIAL IN REAL SPACE',it,999,1,1,n&
     &x,ny,irc,LINEAR)
         else if (ndj==0) then
            call displayv(vpotd,' ION CURRENT IN REAL SPACE',it,999,1,1,&
     &nx,ny,irc,LINEAR)
         else if (nde==0) then
            call displayv(vpotd,' VRPOTENTIAL IN REAL SPACE',it,999,1,1,&
     &nx,ny,irc,LINEAR)
         endif
! display magnetic energy
         if (nda==0) then
            call displays(vpotb,' MAGNETIC ENERGY IN FOURIER SPACE',it,9&
     &99,1,1,modesx,modesy2,irc,LINEAR)
         else if (ndj==0) then
            call displays(vpotb,' ION CURRENT**2 IN FOURIER SPACE',it,99&
     &9,1,1,modesx,modesy2,irc,LINEAR)
         else if (nde==0) then
            call displays(vpotb,' MAGNETIC RADIATIVE ENERGY IN FOURIER S&
     &PACE',it,999,1,1,modesx,modesy2,irc,LINEAR)
         endif
!        call DISPC(vpotb,vpk,label,zsc,zst,2,modesxy2,modesxy2,1,chr,ch&
!    &rs,irc)
! display absolute value of time derivative of vector potential
         if ((nda==0).or.(nde==0)) then
            call ptmodes(vpotd,vpotp,nx,ny,modesx,modesy,order=LINEAR)
            isign = 1
            call fft(vpotd,isign,mixup,sct,tfft,inx,iny,order=LINEAR)
            if (nda==0) then
               call displayv(vpotd,' EPERP IN REAL SPACE',it,999,1,1,nx,&
     &ny,irc,LINEAR)
            else if (nde==0) then
               call displayv(vpotd,' ERPERP IN REAL SPACE',it,999,1,1,nx&
     &,ny,irc,LINEAR)
            endif
! display transverse electric energy
            if (nda==0) then
               call displays(vpote,' TRANSVERSE ELECTRIC ENERGY IN K SPA&
     &CE',it,999,1,1,modesx,modesy2,irc,LINEAR)
            else if (nde==0) then
               call displays(vpote,' TRANSVERSE RADIATIVE ELECTRIC ENERG&
     &Y IN K SPACE',it,999,1,1,modesx,modesy2,irc,LINEAR)
            endif
!           call DISPC(vpote,vpk,prompt,zsc,zst,2,modesxy2,modesxy2,1,ch&
!    &r,chrs,irc)
         endif
         if (irc > 127) then
            ntr = irc - 128
            if (ntr==0) call CLRSCRN
            irc = 0
         endif
         if (irc==1) go to 200
      endif
      endif
      vpotp = vpotx
      if (it < mta) go to 10
   20 call RSTSCRN
      nt = it
      if (nt /= mta) then
         write (2,981) nt, mta
      endif
      nt = min(nt,mta)
      call SETNPLT(1,irc)
! display magnetic energy
      if ((nda==0).or.(nde==0)) then
         if (nda==0) then
            label = ' MAGNETIC FIELD ENERGY VERSUS TIME'
         else if (nde==0) then
            label = ' RADIATIVE MAGNETIC FIELD ENERGY VERSUS TIME'
         endif
         chr = ' RUNID='
         chr = trim(chr)//cdrun
         call DISPS(potc,label,time(1),time(nt),999,2,nt,chr,irc)
         if (irc==1) go to 200
! display transverse electric energy
         if (nda==0) then
            label = ' TRANSVERSE ELECTRIC FIELD ENERGY VERSUS TIME'
         else if (nde==0) then
            label = ' TRANSVERSE RADIATIVE ELECTRIC FIELD ENERGY VERSUS &
     &TIME'
         endif
         chr = ' RUNID='//cdrun
         call DISPS(potc(mta+2),label,time(2),time(nt),999,2,nt-1,chr,ir&
     &c)
         if (irc==1) go to 200
! write out energies versus time
         ftname = 'energyt.'//cdrun
         open(unit=20,file=trim(ftname),form='formatted',status='replace&
     &')
         write (20,*) 'time wb wet'
         do j = 2, nt
            write (20,*) time(j), potc(j), potc(j+mta)
         enddo
         close(unit=20)
! display average magnetic field energy versus |k|
         vpots = vpots/real(nt)
         sum1 = sum1/real(nt)
         write (ftname,991) sum1
         chr = ' RUNID='//cdrun
         chr = trim(ftname)//chr
         if (nda==0) then
            label = ' TIME-AVERAGED MAGNETIC ENERGY VERSUS |K|'
         else if (nde==0) then
            label = ' TIME-AVERAGED RADIATIVE MAGNETIC ENERGY VERSUS |K|&
     &'
         endif
         call DISPC(vpots,vpk,label,zsc,zst,2,modesxy2,modesxy2,1,chr,ch&
     &rs,irc)
         if (irc==1) go to 200
! display average transverse field energy versus |k|
         vpotr = vpotr/real(nt)
         sum2 = sum2/real(nt)
         write (ftname,991) sum2
         chr = ' RUNID='//cdrun
         chr = trim(ftname)//chr
         if (nda==0) then
            label = ' TIME-AVERAGED TRANSVERSE ELECTRIC ENERGY VERSUS |K&
     &|'
         else if (nde==0) then
            label = ' TIME-AVERAGED TRANSVERSE RADIATIVE ELECTRIC ENERGY&
     &VERSUS |K|'
         endif
         call DISPC(vpotr,vpk,label,zsc,zst,2,modesxy2,modesxy2,1,chr,ch&
     &rs,irc)
         if (irc==1) go to 200
! write out energies versus abs(k)
         ftname = 'energyk.'//cdrun
         open(unit=21,file=trim(ftname),form='formatted',status='replace&
     &')
         write (21,*) 'k wb wet'
         do k = 1, modesy2
         do j = 1, modesx
            write (21,*) vpk(j,k), vpots(j,k), vpotr(j,k)
         enddo
         enddo
         close(unit=21)
      endif
! write out data
      write (18,891)
      k0 = 0
      k1 = modesy
      do j = 1, modesx
      j0 = j - 1
      write (18,892) j0, k0, vpots(j,k1)
      enddo
      do k = 2, modesy
      k0 = k - 1
      k1 = modesy + k0
      k2 = modesy - k0
      do j = 1, modesx
      j0 = j - 1
      write (18,892) j0, k0, vpots(j,k1)
      write (18,892) j0, -k0, vpots(j,k2)
      enddo
      enddo
      if (irc==1) go to 200
! quit if no memory available
      if (nvpot==2) then
         write (2,*) 'fatal vpott allocation error=', nvpot
         go to 210
      endif
! allocate correlation time data
      nt2 = nt + nt
      inft = 0
      nft = 1
   30 inft = inft + 1
      nft = 2*nft
      if (nt2 > nft) go to 30
      nfth = nft/2
      allocate(g(2*nft))  
      allocate(mixupt(nft),sctt(nfth))
      nts = nt
      ntq = nts
      ntc = nt/3
      if (ntq > nt) ntq = nt
! create fft table for correlations
      isign = 0
      call FFT1C(g,isign,mixupt,sctt,inft,nft,nfth)
! allocate frequency data 
      allocate(wm(iwmax),p(2*iwmax),pc(2*iwmax))
      allocate(ps(2*iwmax),pcs(2*iwmax))
! transpose data if it will not fit in memory
      if (nvpot==1) then
         fname = 't'//trim(fname)
         if (nf==2) fname = trim(fname)//'-'//trim(cdrun0)
         inquire(file=fname,exist=ex)
         if (ex) then
            it = 0
            call bfopen(vpotts,nt,ndim,17,it,trim(fname))
            if (it /= 1) go to 210
         else
            it = -1
            call bfopen(vpotts,nt,ndim,17,it,trim(fname))
            if (it /= 1) go to 210
            label = 'transposing data'
            do k = 1, modesy2
               write (prompt,'("k=",i6)') k
               call PTOTPUT(label,prompt)
               do i = 1, nt
                  call readbf(vpotx,modesx,modesy2,15,i,irc,order=LINEAR&
     &)
                  if (irc /= 0) go to 210
                  if (nf==2) then
                     call readbf(vpoty,nodesx,nodesy2,16,i,irc,order=LIN&
     &EAR)
                     if (irc /= 0) go to 210
                     vpotx = vpotx - vpoty
                  endif
                  vpott(i,:,:,1) = vpotx(:,:,k)
               enddo
               do j = 1, modesx
               vpotts = vpott(:,:,j,1)
               call writebf(vpotts,nt,ndim,17,it,order=LINEAR)
               enddo
            enddo
         endif
      endif
! open file for spectrum versus omega
      ftname = 'spectrumw.'//cdrun
      open(unit=22,file=trim(ftname),form='formatted',status='replace')
      write (22,*) 'spectrum versus omega'
! open file for  spectrum versus k
      ftname = 'spectrumk.'//cdrun
      open(unit=23,file=trim(ftname),form='formatted',status='replace')
      write (23,*) 'spectrum versus k'
! open file for  correlated spectrum versus omega
      ftname = 'cspectrumw.'//cdrun
      open(unit=24,file=trim(ftname),form='formatted',status='replace')
      write (24,*) 'correlated spectrum versus omega'
! open file for  correlated spectrum versus k
      ftname = 'cspectrumk.'//cdrun
      open(unit=25,file=trim(ftname),form='formatted',status='replace')
      write (25,*) 'correlated spectrum versus k'
! write out light wave dispersion relation
      ftname = 'lightw.'//cdrun
      open(unit=26,file=trim(ftname),form='formatted',status='replace')
      at1 = (vpk(modesx,modesy2)/ci)/100
      write (26,*) 'k w'
      do j = 1, 100
      ak = at1*real(j - 1)
      wk = sqrt(1.0 + ak*ak)
      write (26,*) ak*ci, wk
      enddo
      close(unit=26)
! enter analysis parameters
   50 call MENUCRV2(code,cp,ip,ap,nc,ncc,nci,ncr,istyle,irc)
      if (irc==1) go to 210
      if (irc==3) go to 50
      ntd = ntq
      prompt = ' Hit carrriage return or enter key to continue '
! make sure time display parameters make sense
      if ((lts < 1).or.(((nts - 1)*its + lts) > nt)) then
         write (label,982) lts, its, nts, nt
         call GTINPUT(label,prompt,ftname,irc)
         go to 50
      endif
      nts1 = nts + 1
      nts2 = nts + nts
      dts = dtp*real(its)
      tt0 = t0 + dtp*real(lts - 1)
! make sure mode number parameters make sense
      if ((kxmin < 0) .or. (kxmax >= modesx)) then
         write (label,983) kxmin, kxmax, modesx
         call GTINPUT(label,prompt,ftname,irc)
         go to 50
!        stop
      endif
      if ((kymin < 0) .or. (kymax >= modesy)) then
         write (label,984) kymin, kymax, modesy
         call GTINPUT(label,prompt,ftname,irc)
         go to 50
      endif
! make sure lag time and time display parameters make sense
      if ((ntd < 0).or.(ntd > nts).or.(ntc < 0).or.(ntc > nts)) then
         write (label,985) ntc, ntd, nts
         call GTINPUT(label,prompt,ftname,irc)
         go to 50
      endif
      ntc1 = ntc + 1
      ntc2 = 2*ntc
! make sure vector field parameter makes sense
      if ((nvf < 0).or.(nvf > 6)) then
         write (label,986) nvf
         call GTINPUT(label,prompt,ftname,irc)
         go to 50
      endif
! make sure frequency parameters make sense
      if (dw==0.0) then
         iw = 0
      else
         iw = (wmax - wmin)/dw + 1.001
      endif
      if (iw < 0) then
         write (label,987) iw, iwmax
         call GTINPUT(label,prompt,ftname,irc)
         go to 50
      else if (iw > iwmax) then
         deallocate(wm,p,pc,ps,pcs)
         iwmax = iw
         write (2,*) 'reallocating frequency data, iwmax=', iwmax
         allocate(wm(iwmax),p(2*iwmax),pc(2*iwmax))
         allocate(ps(2*iwmax),pcs(2*iwmax))
      endif
      iw2 = 2*iw
      iw1 = iw + 1
      do 60 k = 1, iw
      wm(k) = wmin + dw*real(k - 1)
   60 continue
! set number of plots per page
      call SETNPLT(nplot,irc)
! initialize timer
      call wtimer(times,ltime,-1)
! display analysis parameters
      call WPCORRV2(runid,indx,indy,nta,modesx,modesy,psolve,t0,tend,dt,&
     &lts,its,nts,kxmin,kxmax,kymin,kymax,ntd,ntc,nvf,wmin,wmax,dw,irc)
      if ((irc==1).or.(irc==3)) go to 50
      kxmn = kxmin + 1
      kxmx = kxmax + 1
      kymn = kymin + 1
      kymx = kymax + 1
! clear screen
      call CLRSCRN
      if (iw==0) go to 130
! begin main loop to find dispersion relation
      ps = 0.0; pcs = 0.0
      kymn0 = max(1,2*kymn-2)
      do 120 k = kymn0, 2*kymx-1
      k1 = k/2
      if ((k > 1) .and. (2*k1 /= k)) k1 = -k1
      k2 = abs(k1) - kymn + 2
      dky = dny*real(k1)
      kk = k
! unpack transposed data
      if (nvpot==1) then
         it = modesx*(k - 1) + 1
         do j = 1, modesx
            call readbf(vpotts,nt,ndim,17,it,irc,order=LINEAR)
            if (irc /= 0) go to 210
            vpott(:,:,j,1) = vpotts
         enddo
         kk = 1
      endif
      do 110 j = kxmn, kxmx
      j1 = j - 1
      j2 = j1 - kxmn + 2
      jk = j2 + (kxmx - kxmn + 1)*(k2 - 1)
      dkx = dnx*real(j1)
      ak2 = dkx*dkx + dky*dky
      ak = sqrt(ak2)
      do 100 n = 1, ndim
! load data in proper format
      do 70 i = 1, nts
      ii = lts + (i - 1)*its
! vector potential
      zt1 = vpott(ii,n,j,kk)
! electric field
      select case(nvf)
      case (2)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt1 = (zt2 - zt1)*dti
! magnetic field
      case (3)
         if (n==1) then
            zt1 = dky*vpott(ii,3,j,kk)
         else if (n==2) then
            zt1 = -dkx*vpott(ii,3,j,kk)
         else if (n==3) then
            zt1 = dkx*vpott(ii,2,j,kk) - dky*vpott(ii,1,j,kk)
         endif
! jperp
      case (4)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt3 = zt1
         if (ii < nt) zt3 = vpott(ii+1,n,j,kk)
         zt1 = ak2*zt1 + (zt3 - 2.0*zt1 + zt2)*dti2c
! edot
      case (5)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt3 = zt1
         if (ii < nt) zt3 = vpott(ii+1,n,j,kk)
         zt1 = -(zt3 - 2.0*zt1 + zt2)*dti2c
! jdot
      case (6)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt3 = zt1
         if (ii < nt) zt3 = vpott(ii+1,n,j,kk)
         zt4 = ak2*zt1 + (zt3 - 2.0*zt1 + zt2)*dti2c
         zt2 = zt1
         zt1 = zt3
         if ((ii+1) < nt) zt3 = vpott(ii+2,n,j,kk)
         zt1 = ak2*zt1 + (zt3 - 2.0*zt1 + zt2)*dti2c
         zt1 = (zt1 - zt4)*dti
      end select
      potc(i) = real(zt1)
      potc(i+nts) = aimag(zt1)
   70 continue
! perform frequency analysis of raw data
      call SPECT(potc,wm,p,tt0,dts,nts,iw,nts2,iw2)
! add contributions from different modes
      ps = ps + p
! extract maximum frequency value
      itm = maxloc(p(1:iw))
      at1 = wm(itm(1))
      itm = maxloc(p(iw1:iw2))
      at2 = wm(itm(1))
      if (k1 > 0) then
         vpwk(n,jk,1) = at1
         vpwk(n,jk,3) = at2
      else if (k1 < 0) then
         vpwk(n,jk,2) = at1
         vpwk(n,jk,4) = at2
      else
         vpwk(n,jk,1) = at1
         vpwk(n,jk,2) = at1
         vpwk(n,jk,3) = at2
         vpwk(n,jk,4) = at2
      endif
      vpak(n,jk,:) = ak
!     write (18,893)
!     write (18,895) n, j1, k1, dkx, dky, ak
!     if (k1 >= 0) then
!        write (18,896) at1, vpwk(n,jk,1), -at1, vpwk(n,jk,3)
!     else
!        write (18,896) at1, vpwk(n,jk,2), -at1, vpwk(n,jk,4)
!     endif
      if (ntc==0) then
! write out spectrum for raw data
!        write (18,994)
!        do 80 i = 1, iw
!        i1 = i + iw
!        at1 = wm(i)
!        at2 = -at1
!        write (18,995) at1, p(i), at2, p(i1)
!  80    continue
         go to 100
      endif
! auto-correlation function
      call CORRE(potc,g,g,potc,mixupt,sctt,inft,nts,ntc,nts2,nft,nfth,nt&
     &c2)
! perform frequency analysis of correlated data
      call SPECT(potc,wm,pc,tt0,dts,ntc,iw,ntc2,iw2)
! add contributions from different modes
      pcs = pcs + pc
! extract maximum frequency value
      itm = maxloc(pc(1:iw))
      at1 = wm(itm(1))
      itm = maxloc(pc(iw1:iw2))
      at2 = wm(itm(1))
      if (k1 > 0) then
         vpck(n,jk,1) = at1
         vpck(n,jk,3) = at2
      else if (k1 < 0) then
         vpck(n,jk,2) = at1
         vpck(n,jk,4) = at2
      else
         vpck(n,jk,1) = at1
         vpck(n,jk,2) = at1
         vpck(n,jk,3) = at2
         vpck(n,jk,4) = at2
      endif
!     write (18,894)
!     write (18,895) n, j1, k1, dkx, dky, ak
!     if (k1 >= 0) then
!        write (18,896) at1, vpck(n,jk,1), -at1, vpck(n,jk,3)
!     else
!        write (18,896) at1, vpck(n,jk,2), -at1, vpck(n,jk,4)
!     endif
! write out spectrum for raw and correlated data
!     write (18,996)
!     do 90 i = 1, iw
!     i1 = i + iw
!     at1 = wm(i)
!     at2 = -at1
!     write (18,997) at1, p(i), at2, p(i1), at1, pc(i), at2, pc(i1)
!  90 continue
  100 continue
  110 continue
  120 continue
! display frequency information
      call SETNPLT(1,irc)
! normalize and display integrated spectrum versus frequency
      if ((nvf==2).or.(ndj==0)) then
         anorm = bnorm
      else
         anorm = bnorm/(ci*ci)
      endif
      sum1 = 0.0d0
      sum2 = 0.0d0
      write (18,993)
      do i = 1, iw
         at1 = anorm*ps(i)
         at2 = anorm*ps(i+iw)
         sum1 = sum1 + at1
         sum2 = sum2 + at2
         ps(i) = at1
         ps(i+iw) = at2
         write (18,995) wm(i), at1, at2, at1+at2
      enddo
      sum1 = sum1 + sum2 - ps(iw+1)
      write (ftname,992) sum1
      label = ' INTEGRATED SPECTRUM VERSUS OMEGA'//ftname
      if (nda==0) then
         chr = trim(cnvf(nvf))//', RUNID='//cdrun
      else if (ndj==0) then
         chr = ' ION CURRENT, RUNID='//cdrun
      else if (nde==0) then
         chr = trim(crnvf(nvf))//', RUNID='//cdrun
      endif
      ncd = 58
      call DISPR(ps,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),chr&
     &s(3),irc)
      if ((irc==1).or.(irc==3)) go to 50
! write out spectrum versus omega
      write (22,*) 'nvf=', nvf
      write (22,*) 'omega ps(+) ps(-) ps log(ps)'
      do i = 1, iw
         at1 = ps(i)+ps(i+iw)
         write (22,*) wm(i), ps(i), ps(i+iw), at1, alog(at1)
      enddo
      end file 22
      backspace 22
! display maximum frequency
      zsc = cmplx(999.0,999.0)
      zst = cmplx(2.0,2.0)
      label = ' FREQUENCY W VERSUS |K|'
      it1 = ndim*(kxmx - kxmn + 1)*(kymx - kymn + 1)
      call DISPC(vpwk,vpak,label,zsc,zst,2,it1,ndim*modesxy,ns,chr,chws,&
     &irc)
      if ((irc==1).or.(irc==3)) go to 50
! write out spectrum versus k
      write (23,*) 'nvf=', nvf
      write (23,*) 'k ps'
      do k = 1, ns
      do j = 1, modesxy
      do i = 1, ndim
         write (23,*) vpak(i,j,k), vpwk(i,j,k)
      enddo
      enddo
      enddo
      end file 23
      backspace 23
      if (ntc > 0) then
! normalize and display integrated correlation spectrum versus frequency
         sum1 = 0.0d0
         sum2 = 0.0d0
         write (18,994)
         do i = 1, iw
            at1 = anorm*pcs(i)
            at2 = anorm*pcs(i+iw)
            sum1 = sum1 + at1
            sum2 = sum2 + at2
            pcs(i) = at1
            pcs(i+iw) = at2
            write (18,995) wm(i), at1, at2, at1+at2
         enddo
         sum1 = sum1 + sum2 - pcs(iw+1)
         write (ftname,992) sum1
         label = ' INTEGRATED CORRELATION SPECTRUM VS OMEGA'//ftname
         if (nda==0) then
            chr = trim(cnvf(nvf))//' CORRELATION, RUNID='//cdrun
         else if (ndj==0) then
            chr = 'ION CURRENT CORRELATION, RUNID='//cdrun
         else if (nde==0) then
            chr = trim(crnvf(nvf))//' CORRELATION, RUNID='//cdrun
         endif
         ncd = 58
         call DISPR(pcs,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd)&
     &,chrs(3),irc)
         if ((irc==1).or.(irc==3)) go to 50
! write out correlated spectrum versus omega
         write (24,*) 'nvf=', nvf
         write (24,*) 'omega ps(+) ps(-) ps log(ps)'
         do i = 1, iw
            at1 = ps(i)+ps(i+iw)
            write (24,*) wm(i), ps(i), ps(i+iw), at1, alog(at1)
         enddo
         end file 24
         backspace 24
! display maximum correlation frequency
         label = ' FREQUENCY W VERSUS |K|'
         call DISPC(vpck,vpak,label,zsc,zst,2,it1,ndim*modesxy,ns,chr,ch&
     &ws,irc)
         if ((irc==1).or.(irc==3)) go to 50
! write out correlated spectrum versus k
         write (25,*) 'nvf=', nvf
         write (25,*) 'k ps'
         do k = 1, ns
         do j = 1, modesxy
         do i = 1, ndim
            write (25,*) vpak(i,j,k), vpwk(i,j,k)
         enddo
         enddo
         enddo
         end file 25
         backspace 25
      endif
      call SETNPLT(nplot,irc)
! begin main loop to display individual modes
  130 kymn0 = max(1,2*kymn-2)
      do 190 k = kymn0, 2*kymx-1
      k1 = k/2
      if ((k > 1) .and. (2*k1 /= k)) k1 = -k1
      dky = dny*real(k1)
      kk = k
! unpack transposed data
      if (nvpot==1) then
         it = modesx*(k - 1) + 1
         do j = 1, modesx
            call readbf(vpotts,nt,ndim,17,it,irc,order=LINEAR)
            if (irc /= 0) go to 210
            vpott(:,:,j,1) = vpotts
         enddo
         kk = 1
      endif
      do 180 j = kxmn, kxmx
      j1 = j - 1
      dkx = dnx*real(j1)
      ak2 = dkx*dkx + dky*dky
      ak = sqrt(ak2)
      do 170 n = 1, ndim
! load data in proper format
      do 140 i = 1, nts
      ii = lts + (i - 1)*its
! vector potential
      zt1 = vpott(ii,n,j,kk)
! electric field
      select case(nvf)
      case (2)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt1 = (zt2 - zt1)*dti
! magnetic field
      case (3)
         if (n==1) then
            zt1 = dky*vpott(ii,3,j,kk)
         else if (n==2) then
            zt1 = -dkx*vpott(ii,3,j,kk)
         else if (n==3) then
            zt1 = dkx*vpott(ii,2,j,kk) - dky*vpott(ii,1,j,kk)
         endif
! jperp
      case (4)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt3 = zt1
         if (ii < nt) zt3 = vpott(ii+1,n,j,kk)
         zt1 = ak2*zt1 + (zt3 - 2.0*zt1 + zt2)*dti2c
! edot
      case (5)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt3 = zt1
         if (ii < nt) zt3 = vpott(ii+1,n,j,kk)
         zt1 = -(zt3 - 2.0*zt1 + zt2)*dti2c
! jdot
      case (6)
         zt2 = zt1
         if (ii > 1) zt2 = vpott(ii-1,n,j,kk)
         zt3 = zt1
         if (ii < nt) zt3 = vpott(ii+1,n,j,kk)
         zt4 = ak2*zt1 + (zt3 - 2.0*zt1 + zt2)*dti2c
         zt2 = zt1
         zt1 = zt3
         if ((ii+1) < nt) zt3 = vpott(ii+2,n,j,kk)
         zt1 = ak2*zt1 + (zt3 - 2.0*zt1 + zt2)*dti2c
         zt1 = (zt1 - zt4)*dti
      end select
      potc(i) = real(zt1)
      potc(i+nts) = aimag(zt1)
      time(i) = tt0 + dts*real(i - 1)
  140 continue
! time display of raw data
      if (ntd > 0) then
         if (nda==0) then
            label = trim(cnvf(nvf))//' VERSUS TIME, RUNID='//cdrun
         else if (ndj==0) then
            label = 'ION CURRENT VERSUS TIME, RUNID='//cdrun
         else if (nde==0) then
            label = trim(crnvf(nvf))//' VERSUS TIME, RUNID='//cdrun
         endif
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         call DISPR(potc,label,time(1),time(ntd),999,0,0,ntd,nts,itwo,ch&
     &r(1:ncd),chrs,irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
! plot complex phase
         if (nda==0) then
            label = ' REAL VERSUS IMAGINARY PART OF'//trim(cnvf(nvf))
         else if (ndj==0) then
            label = ' REAL VERSUS IMAGINARY PART OF ION CURRENT'
         else if (nde==0) then
            label = ' REAL VERSUS IMAGINARY PART OF'//trim(crnvf(nvf))
         endif
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         zsc = cmplx(999.,999.)
         zst = cmplx(0.,0.)
         call DISPC(potc,potc(nts1),label,zsc,zst,1,ntd,ntd,1,chr(1:ncd)&
     &,' ',irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
! log display of raw data
         do 150 i = 1, nts
         i1 = nts + i
         at1 = sqrt(potc(i)**2 + potc(i1)**2)
         if (at1 < pmin) at1 = pmin
         g(i) = alog(at1)
  150    continue
         if (nda==0) then
            label = ' LN |'//trim(adjustl(cnvf(nvf)))//'| VERSUS TIME, R&
     &UNID='
         else if (ndj==0) then
            label = ' LN |'//'ION CURRENT'//'| VERSUS TIME, RUNID='
         else if (nde==0) then
            label = ' LN |'//trim(adjustl(crnvf(nvf)))//'| VERSUS TIME, &
     &RUNID='
         endif
         label = trim(label)//cdrun
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         call DISPS(g,label,time(1),time(ntd),999,2,ntd,chr(1:ncd),irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
      endif
      if (iw > 0) then
! perform frequency analysis of raw data
         call SPECT(potc,wm,p,tt0,dts,nts,iw,nts2,iw2)
! power spectrum of raw data
         if (nda==0) then
            label = trim(cnvf(nvf))//' VERSUS OMEGA, RUNID='//cdrun
         else if (ndj==0) then
            label = 'ION CURRENT VERSUS OMEGA, RUNID='//cdrun
         else if (nde==0) then
            label = trim(crnvf(nvf))//' VERSUS OMEGA, RUNID='//cdrun
         endif
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         call DISPR(p,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),c&
     &hrs(3),irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
      endif
      if (ntc==0) go to 170
! auto-correlation function
      call CORRE(potc,g,g,potc,mixupt,sctt,inft,nts,ntc,nts2,nft,nfth,nt&
     &c2)
      if (ntd > 0) then
! time display of correlated data
         it1 = ntd
         if (it1 > ntc) it1 = ntc
         if (nda==0) then
            label = trim(snvf(nvf))//' CORRELATION  VERSUS TIME, RUNID='
         else if (ndj==0) then
            label = 'ION CURRENT CORRELATION  VERSUS TIME, RUNID='
         else if (nde==0) then
            label = trim(srnvf(nvf))//' CORRELATION  VERSUS TIME, RUNID=&
     &'
         endif
         label = trim(label)//cdrun
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         call DISPR(potc,label,time(1),time(it1),999,0,0,it1,ntc,itwo,ch&
     &r(1:ncd),chrs,irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
! plot complex phase
         if (nda==0) then
            label = ' REAL VERSUS IMAGINARY PART OF'//trim(snvf(nvf))
         else if (ndj==0) then
            label = ' REAL VERSUS IMAGINARY PART OF ION CURRENT'
         else if (nde==0) then
            label = ' REAL VERSUS IMAGINARY PART OF'//trim(srnvf(nvf))
         endif
         label = trim(label)//' CORRELATION'
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         zsc = cmplx(999.,999.)
         zst = cmplx(0.,0.)
         call DISPC(potc,potc(ntc1),label,zsc,zst,1,it1,it1,1,chr(1:ncd)&
     &,' ',irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
! log display of correlated data
         do 160 i = 1, ntc
         i1 = ntc + i
         at1 = sqrt(potc(i)**2 + potc(i1)**2)
         if (at1 < pmin) at1 = pmin
         g(i) = alog(at1)
  160    continue
         if (nda==0) then
            label = ' LN |'//trim(adjustl(snvf(nvf)))//' CORRELATION| VE&
     &RSUS TIME, RUNID='
         else if (ndj==0) then
            label = ' LN |'//'/ION CURRENT'//' CORRELATION| VERSUS TIME,&
     &RUNID='
         else if (nde==0) then
            label = ' LN |'//trim(adjustl(srnvf(nvf)))//' CORRELATION| V&
     &ERSUS TIME, RUNID='
         endif
         label = trim(label)//cdrun
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         call DISPS(g,label,time(1),time(it1),999,2,it1,chr(1:ncd),irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
      endif
      if (iw > 0) then
! perform frequency analysis of correlated data
         call SPECT(potc,wm,pc,tt0,dts,ntc,iw,ntc2,iw2)
! power spectrum of correlated data
         if (nda==0) then
            label = trim(snvf(nvf))//' CORRELATION VERSUS OMEGA, RUNID='
         else if (ndj==0) then
            label = 'ION CURRENT CORRELATION VERSUS OMEGA, RUNID='
         else if (nde==0) then
            label = trim(srnvf(nvf))//' CORRELATION VERSUS OMEGA, RUNID=&
     &'
         endif
         label = trim(label)//cdrun
         write (chr,996) n, j1, k1, dkx, dky
         ncd = 58
         call DISPR(pc,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),&
     &chrs(3),irc)
         if (irc==1) go to 200
         if (irc==3) go to 50
      endif
  170 continue
  180 continue
  190 continue
  200 call wtimer(times,ltime)
      go to 50
! close graphics device
  210 call GRCLOSE
      end program
