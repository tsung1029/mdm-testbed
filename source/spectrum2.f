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
      program spectrum2
      use globals, only: LINEAR, QUADRATIC
      use fft2d
      use field2d
      use diag2d
      implicit none
      integer, parameter :: ncc=0, nci=10, ncr=3, nc=ncc+nci+ncr, ns=4
      integer :: idpal = 1
      integer :: ndd, ndp, nrec, modesx, modesy, modesxy, lts, its, nts
      integer :: kxmax, kxmin, kymax, kymin, kxmx, kxmn, kymx, kymn
      integer :: ntq, ntc, nplot, istyle, itwo, iwmax, kymn0
      integer :: ltime, modesy2, modesxy2
      integer :: nx, nxh, ny, nyh, nxe, nye, nxeh, nxyh, nxhy, nx1, ny1
      integer :: it1, mtp, nt2d, it, nt, nt2, inft, nft, nfth, isign
      integer :: i, j, k, ii, i1, j1, k1, irc, nts1, nts2, ntc1, ntc2
      integer :: j0, k0, j2, k2, jk, kk, iw, iw1, iw2, ncd, ntr
      integer :: inx, iny, nxd, nyd, nxdh, npot
      integer :: nf, nodesx, nodesy, nodesy2, mt
      real :: wmin, wmax, dw, pmin, anorm
      real :: at1, at2, dtp, dnx, dny, dts, tt0, times, dkx, dky, ak, wk
      real :: ak2, vth, tfft
      double precision :: sum1, sum2
      complex :: zt1, zsc, zst
      logical :: ex
      integer, dimension(1) :: itm
      integer, dimension(nci) :: ip
      real, dimension(ncr) :: ap
      real, dimension(:,:), pointer :: potb, pots, potr, pk
      real, dimension(:,:), pointer :: potd
      complex, dimension(:,:), pointer :: potx, poty
      real, dimension(:,:), pointer :: pwk, pck, pak
      real, dimension(:), pointer :: potc, time
      complex, dimension(:,:,:), pointer ::  pott
      complex, dimension(:,:), pointer ::  potts
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
      equivalence (wmin, ap(1)), (wmax, ap(2)), (dw, ap(3))
   91 format (' SPECTRUM ANALYSIS FOR 2D DATA',A32)
   92 format (' RUNID= ',a20,' INDX=',i3,' INDY=',i3)
   93 format (' NTP=',i6,' MODESX=',i5,' MODESY=',i5,' PSOLVE=',i2)
   94 format (' NREC=',i10)
   95 format (' T0=',f8.1,' TEND=',f8.1,' DT=',f8.5)
  981 format (' ACTUAL LENGTH OF DATA = ',i6,' DATA EXPECTED = ',i6)
  982 format (' ERROR IN PARAMETERS, LTS,ITS,NTS,NT = ',4i7)
  983 format (' ERROR IN MODE PARAMETERS, KXMIN,KXMAX,MODESX = ',3i7)
  984 format (' ERROR IN MODE PARAMETERS, KYMIN,KYMAX,MODESY = ',3i7)
  985 format (' ERROR IN DISPLAY PARAMETERS, NTC,NTD,NTS = ',3i7)
  987 format (' ERROR IN FREQUENCY PARAMETERS, IW,IWMAX = ',2i7)
  891 format (' KX, KY, AVERAGE POTENTIAL ENERGY')
  892 format (1x,2i4,1x,e14.7)
  893 format ('  PEAK FREQUENCY OF RAW DATA')
  894 format ('  PEAK FREQUENCY OF CORRRELATED DATA')
  895 format (' KX ',i4,', KY ',i4,', AKX=',f7.4,', AKY=',f7.4,', AK=',f&
     &7.4)
  896 format (' MAX P(',f6.3,') = ',e14.7,' MAX P(',f7.3,') = ',e14.7)
!
  991 format (' <ENERGY> = ',e14.7,',')
  992 format (', SUM = ',e14.7)
  993 format (' INTEGRATED SPECTRUM VERSUS OMEGA')
  994 format (' INTEGRATED CORRELATION SPECTRUM VS OMEGA')
  995 format (1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)
  996 format (' KX ',i4,', KY ',i4,', AKX=',f7.4,', AKY=',f7.4,', AK=',f
     17.4)
! 994 format ('  POWER SPECTRUM:',6x,' RAW DATA')
! 995 format (' P(',f6.3,') = ',e14.7,' P(',f6.3,') = ',e14.7)
! 996 format ('  POWER SPECTRUM:',6x,' RAW DATA',50x,' CORRELATED DATA')
! 997 format (' P(',f6.3,') = ',e14.7,' P(',f6.3,') = ',e14.7,9x,' P(',f&
!    &6.3,') = ',e14.7,' P(',f6.3,') = ',e14.7)
! istyle = (0,1) = use (brief,full) menu style
      data istyle /1/
      data chrs /'   REAL   ','IMAGINARY ',' POSITIVE ',' NEGATIVE '/
      data chws /' KY>0, W>0',' KY<0, W>0',' KY>0, W<0',' KY<0, W<0'/
      data code /'LTS   ','ITS   ','NTS   ','KXMIN ','KXMAX ','KYMIN ','&
     &KYMAX ','NTD   ','NTC   ','NPLOT ','WMIN  ','WMAX  ','DW    '/
      data itwo /2/
      data pmin /1.0e-14/
! lts,its = initial data point and increment between
      data lts,its /1,1/
      data iwmax /5001/
! wmin,wmax,dw = initial, final frequency and increment
      data wmin,wmax,dw /0.,2.0,.01/
! nplot = number of plots per page
      data nplot /4/
! open graphics device
      call GROPEN0(1)
! get file name
      label = ' '
!   5 prompt = 'enter name of diagnostic metafile, q to quit:'
    5 prompt = 'enter runid, or q to quit:'
      call GTINPUT(label,prompt,chr,irc)
      if ((chr=='q').or.(chr=='Q')) go to 190
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
            if (ndp==0) then
               nodesx = modesxp; nodesy = modesyp
               ftname = fpname
            else if (ndd==0) then
               nodesx = modesxd; nodesy = modesyd
               ftname = fdname; cdrun0 = trim(cdrun0)//'d'
            endif
            nxh = 2**(indx-1); nyh = max(1,2**(indy-1))
            if (nodesx > nxh) nodesx = nxh
            if (nodesy > nyh) nodesy = nyh
            nodesy2 = 2*nodesy - 1
            allocate(poty(nodesx,nodesy2))
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
         read (19,den2d,iostat=ndd)
         rewind 19
         read (19,pot2d,iostat=ndp)
         rewind 19
         if ((ndp==0).and.(ndd==0)) then
            label = ' '
            if (nf==2) then
               label = 'for second file'
               if (j==2) label = 'for first file'
            endif
            prompt = 'enter p or d for potential or ion density diagnost&
     &ic:'
    6       call GTINPUT(label,prompt,ftname,irc)
            if ((ftname=='q').or.(ftname=='Q')) then
               close(unit=19)
               label = ' '
               go to 5
            endif
            if ((ftname=='p').or.(ftname=='P')) then
               ndd = -1
            else if ((ftname=='d').or.(ftname=='D')) then
               ndp = -1
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
         else if ((ndp /= 0).and.(ndd /=0)) then
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
      if (ndp==0) then
         read (19,pot2d,iostat=ndp)
         modesx = modesxp; modesy = modesyp
         nrec = nprec
         fname = fpname; ftname = 'pot2out.'//cdrun
         anorm = float((2**indx)*(2**indy))/ceng
      else if (ndd==0) then
         read (19,den2d,iostat=ndd)
         modesx = modesxd; modesy = modesyd
         nrec = ndrec; ntp = ntd
         cdrun = trim(cdrun)//'d'
         fname = fdname; ftname = 'den2out.'//cdrun
         anorm = 1.0
      endif
      modesxy = modesx*modesy
      if (nf==2) cdrun = trim(cdrun)//'-'//trim(cdrun0)
      runid = 'spectrum.'//cdrun
! text output file
      open(unit=18,file=trim(ftname),form='formatted',status='replace')
! echo input
      write (18,91) trim(adjustl(fname))
      write (18,92) runid, indx, indy
      write (18,93) ntp, modesx, modesy, psolve
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
      dtp = dt*real(ntp)
! time parameters
      if (nrec > 0) then
         mtp = nrec
      else
         it1 = (tend - t0)/dt + .0001
         mtp = (it1 - 1)/ntp + 1
      endif
      nt2d = mtp + mtp
! allocate space data
      allocate(potx(modesx,modesy2))
      allocate(potb(modesx,modesy2),pots(modesx,modesy2))
      allocate(potr(modesx,modesy2))
      allocate(pk(modesx,modesy2))
      allocate(pwk(modesxy,ns),pck(modesxy,ns),pak(modesxy,ns))
! allocate time data
      allocate(potc(nt2d),time(mtp))
      allocate(pott(mtp,modesx,modesy2),stat=npot)
      if (npot/=0) then
         write (2,*) 'attempting to allocate in sections'
         allocate(potts(mtp,1),pott(mtp,modesx,1),stat=npot)
         if (npot==0) then
            npot = 1
         else
            npot = 2
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
      allocate(potd(nxd,nyd))
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
      call SAK2(pk,nx,ny,modesx,modesy,modesx,modesy2)
      zsc = cmplx(999.0,999.0)
      zst = cmplx(2.0,2.0)
! debug
      zsc = cmplx(1.0,999.0)
      zst = cmplx(1.0,1.0)
! end debug
      if (ndp==0) then
         label = ' POTENTIAL ENERGY VERSUS |K|'
      else if (ndd==0) then
         label = ' ION DENSITY**2 VERSUS |K|'
      endif
! read diagnostic input
      nt = 0
      call bfopen(potx,modesx,modesy2,11,nt,trim(fname))
      if (nt /= 1) go to 20
      if (nf==2) then
         mt = 0
         call bfopen(poty,nodesx,nodesy2,12,mt,trim(chr))
         if (mt /= 1) go to 20
      endif
   10 call readbf(potx,modesx,modesy2,11,nt,irc,order=LINEAR)
      if (irc /= 0) go to 20
      if (nf==2) then
         call readbf(poty,nodesx,nodesy2,12,mt,irc,order=LINEAR)
         if (irc /= 0) go to 20
         potx = potx - poty
      endif
      it = nt - 1
      if (npot==0) pott(it,:,:) = potx
! create time array
      time(it) = t0 + dtp*real(it - 1)
! calculate potential energy
      if (ndp==0) then
         call WEL2(potx,potb,ceng,potc(it),nx,ny,modesx,modesy,modesx,mo&
     &desy2)
! calculate average of potential energy
         if (it==1) then
            sum1 = potc(it)
            pots = potb
         else
            sum1 = sum1 + potc(it)
            pots = pots + potb
         endif
      else if (ndd==0) then
         call SQ2POT2(potx,potr,nx,ny,modesx,modesy,modesx,modesy2)
      endif
! display raw data
      if (ntr > 0) then
      it1 = (it-1)/ntr
      if ((it-1)==ntr*it1) then
         write (chr,*) 'NT = ', it
         call ptmodes(potd,potx,nx,ny,modesx,modesy,order=LINEAR)
         isign = 1
         call fft(potd,isign,mixup,sct,tfft,inx,iny,order=LINEAR)
! display potential
         if (ndp==0) then
            call displays(potd,' POTENTIAL IN REAL SPACE',it,999,0,1,nx,&
     &ny,irc,LINEAR)
! display ion density
         else if (ndd==0) then
            call displays(potd,' ION DENSITY IN REAL SPACE',it,999,2,1,n&
     &x,ny,irc,LINEAR)
         endif
! display potential energy
         if (ndp==0) then
            call displays(potb,' POTENTIAL ENERGY IN FOURIER SPACE',it,9&
     &99,1,1,modesx,modesy2,irc,LINEAR)
! display ion density
         else if (ndd==0) then
            call displays(potr,' ION DENSITY**2 IN FOURIER SPACE',it,999&
     &,1,1,modesx,modesy2,irc,LINEAR)
         endif
!        call DISPC(potb,pk,label,zsc,zst,2,modesxy2,modesxy2,1,chr,chrs&
!    &,irc)
         if (irc > 127) then
            ntr = irc - 128
            if (ntr==0) call CLRSCRN
            irc = 0
         endif
         if (irc==1) go to 190
      endif
      endif
      if (it < mtp) go to 10
   20 call RSTSCRN
      nt = it
      if (nt /= mtp) then
         write (2,981) nt, mtp
      endif
      nt = min(nt,mtp)
      call SETNPLT(1,irc)
! display potential energy
      if (ndp==0) then
         label = ' POTENTIAL ENERGY VERSUS TIME'
         chr = ' RUNID='
         chr = trim(chr)//cdrun
         call DISPS(potc,label,time(1),time(nt),999,2,nt,chr,irc)
         if (irc==1) go to 180
! write out energies versus time
         ftname = 'penergyt.'//cdrun
         open(unit=20,file=trim(ftname),form='formatted',status='replace&
     &')
         write (20,*) 'time wel'
         do j = 2, nt
            write (20,*) time(j), potc(j)
         enddo
         close(unit=20)
! display average potential energy versus |k|
         pots = pots/real(nt)
         sum1 = sum1/real(nt)
         write (ftname,991) sum1
         chr = ' RUNID='//cdrun
         chr = trim(ftname)//chr
         label = ' TIME-AVERAGED POTENTIAL ENERGY VERSUS |K|'
         call DISPC(pots,pk,label,zsc,zst,2,modesxy2,modesxy2,1,chr,chrs&
     &,irc)
         if (irc==1) go to 180
! write out energies versus abs(k)
         ftname = 'penergyk.'//cdrun
         open(unit=21,file=trim(ftname),form='formatted',status='replace&
     &')
         write (21,*) 'k wel'
         do k = 1, modesy2
         do j = 1, modesx
            write (21,*) pk(j,k), pots(j,k)
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
      write (18,892) j0, k0, pots(j,k1)
      enddo
      do k = 2, modesy
      k0 = k - 1
      k1 = modesy + k0
      k2 = modesy - k0
      do j = 1, modesx
      j0 = j - 1
      write (18,892) j0, k0, pots(j,k1)
      write (18,892) j0, -k0, pots(j,k2)
      enddo
      enddo
      if (irc==1) go to 180
! quit if no memory available
      if (npot==2) then
         write (2,*) 'fatal pott allocation error=', npot
         go to 190
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
      if (npot==1) then
         fname = 't'//trim(fname)
         if (nf==2) fname = trim(fname)//'-'//trim(cdrun0)
         inquire(file=fname,exist=ex)
         if (ex) then
            it = 0
            call bfopen(potts,nt,1,17,it,trim(fname))
            if (it /= 1) go to 190
         else
            it = -1
            call bfopen(potts,nt,1,17,it,trim(fname))
            if (it /= 1) go to 190
            label = 'transposing data'
            do k = 1, modesy2
               write (prompt,'("k=",i6)') k
               call PTOTPUT(label,prompt)
               do i = 1, nt
                  call readbf(potx,modesx,modesy2,15,i,irc,order=LINEAR)
                  if (irc /= 0) go to 190
                  if (nf==2) then
                     call readbf(poty,nodesx,nodesy2,16,i,irc,order=LINE&
     &AR)
                     if (irc /= 0) go to 190
                     potx = potx - poty
                  endif
                  pott(i,:,1) = potx(:,k)
               enddo
               do j = 1, modesx
               potts(:,1) = pott(:,j,1)
               call writebf(potts,nt,1,17,it,order=LINEAR)
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
! write out plasma wave dispersion relation
      vth = 1.0
      ftname = 'plasmaw.'//cdrun
      open(unit=26,file=trim(ftname),form='formatted',status='replace')
      at1 = pk(modesx,modesy2)/100
      write (26,*) 'k w wfs'
      do j = 1, 100
      ak = at1*real(j - 1)
      wk = sqrt(1.0 + 3.0*(ak*vth)**2)
      at2 = sqrt(exp(-ak**2) + 3.0*(ak*vth)**2)
      write (26,*) ak, wk, at2
      enddo
      close(unit=26)
! enter analysis parameters
   50 call MENUCR2(code,cp,ip,ap,nc,ncc,nci,ncr,istyle,irc)
      if (irc==1) go to 190
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
      call WPCORR2(runid,indx,indy,ntp,modesx,modesy,psolve,t0,tend,dt,l&
     &ts,its,nts,kxmin,kxmax,kymin,kymax,ntd,ntc,wmin,wmax,dw,irc)
      if ((irc==1).or.(irc==3)) go to 50
      kxmn = kxmin + 1
      kxmx = kxmax + 1
      kymn = kymin + 1
      kymx = kymax + 1
! clear screen
      call CLRSCRN
      if (iw==0) go to 120
! begin main loop to find dispersion relation
      ps = 0.0; pcs = 0.0
      kymn0 = max(1,2*kymn-2)
      do 110 k = kymn0, 2*kymx-1
      k1 = k/2
      if ((k > 1) .and. (2*k1 /= k)) k1 = -k1
      k2 = abs(k1) - kymn + 2
      dky = dny*real(k1)
      kk = k
! unpack transposed data
      if (npot==1) then
         it = modesx*(k - 1) + 1
         do j = 1, modesx
            call readbf(potts,nt,1,17,it,irc,order=LINEAR)
            if (irc /= 0) go to 190
            pott(:,j,1) = potts(:,1)
         enddo
         kk = 1
      endif
      do 100 j = kxmn, kxmx
      j1 = j - 1
      j2 = j1 - kxmn + 2
      jk = j2 + (kxmx - kxmn + 1)*(k2 - 1)
      dkx = dnx*real(j1)
      ak2 = dkx*dkx + dky*dky
      ak = sqrt(ak2)
! load data in proper format
      do 70 i = 1, nts
      ii = lts + (i - 1)*its
! potential
      zt1 = pott(ii,j,kk)
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
         pwk(jk,1) = at1
         pwk(jk,3) = at2
      else if (k1 < 0) then
         pwk(jk,2) = at1
         pwk(jk,4) = at2
      else
         pwk(jk,1) = at1
         pwk(jk,2) = at1
         pwk(jk,3) = at2
         pwk(jk,4) = at2
      endif
      pak(jk,:) = ak
!     write (18,893)
!     write (18,895) j1, k1, dkx, dky, ak
!     if (k1 >= 0) then
!        write (18,896) at1, pwk(jk,1), -at1, pwk(jk,3)
!     else
!        write (18,896) at1, pwk(jk,2), -at1, pwk(jk,4)
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
         pck(jk,1) = at1
         pck(jk,3) = at2
      else if (k1 < 0) then
         pck(jk,2) = at1
         pck(jk,4) = at2
      else
         pck(jk,1) = at1
         pck(jk,2) = at1
         pck(jk,3) = at2
         pck(jk,4) = at2
      endif
!     write (18,894)
!     write (18,895) j1, k1, dkx, dky, ak
!     if (k1 >= 0) then
!        write (18,896) at1, pck(jk,1), -at1, pck(jk,3)
!     else
!        write (18,896) at1, pck(jk,2), -at1, pck(jk,4)
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
! display frequency information
      call SETNPLT(1,irc)
! normalize and display integrated spectrum versus frequency
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
      if (ndp==0) then
         chr = ' POTENTIAL, RUNID='//cdrun
      else if (ndd==0) then
         chr = ' ION DENSITY, RUNID='//cdrun
      endif
      ncd = 58
      call DISPR(ps,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),chr&
     &s(3),irc)
      if ((irc==1).or.(irc==3)) go to 50
! write out spectrum versus omega
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
      it1 = (kxmx - kxmn + 1)*(kymx - kymn + 1)
      call DISPC(pwk,pak,label,zsc,zst,2,it1,modesxy,ns,chr,chws,irc)
      if ((irc==1).or.(irc==3)) go to 50
! write out spectrum versus k
      write (23,*) 'k ps'
      do k = 1, ns
      do j = 1, modesxy
         write (23,*) pak(j,k), pwk(j,k)
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
         if (ndp==0) then
            chr = ' POTENTIAL CORRELATION, RUNID='//cdrun
         else if (ndd==0) then
            chr = ' ION DENSITY CORRELATION, RUNID='//cdrun
         endif
         ncd = 58
         call DISPR(pcs,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd)&
     &,chrs(3),irc)
         if ((irc==1).or.(irc==3)) go to 50
! write out correlated spectrum versus omega
         write (24,*) 'omega ps(+) ps(-) ps log(ps)'
         do i = 1, iw
            at1 = ps(i)+ps(i+iw)
            write (24,*) wm(i), ps(i), ps(i+iw), at1, alog(at1)
         enddo
         end file 24
         backspace 24
! display maximum correlation frequency
         label = ' FREQUENCY W VERSUS |K|'
         call DISPC(pck,pak,label,zsc,zst,2,it1,modesxy,ns,chr,chws,irc)
         if ((irc==1).or.(irc==3)) go to 50
! write out correlated spectrum versus k
         write (25,*) 'k ps'
         do k = 1, ns
         do j = 1, modesxy
            write (25,*) pak(j,k), pwk(j,k)
         enddo
         enddo
         end file 25
         backspace 25
      endif
      call SETNPLT(nplot,irc)
! begin main loop to display individual modes
  120 kymn0 = max(1,2*kymn-2)
      do 170 k = kymn0, 2*kymx-1
      k1 = k/2
      if ((k > 1) .and. (2*k1 /= k)) k1 = -k1
      dky = dny*real(k1)
      kk = k
! unpack transposed data
      if (npot==1) then
         it = modesx*(k - 1) + 1
         do j = 1, modesx
            call readbf(potts,nt,1,17,it,irc,order=LINEAR)
            if (irc /= 0) go to 190
            pott(:,j,1) = potts(:,1)
         enddo
         kk = 1
      endif
      do 160 j = kxmn, kxmx
      j1 = j - 1
      dkx = dnx*real(j1)
      ak2 = dkx*dkx + dky*dky
      ak = sqrt(ak2)
! load data in proper format
      do 130 i = 1, nts
      ii = lts + (i - 1)*its
! potential
      zt1 = pott(ii,j,kk)
      potc(i) = real(zt1)
      potc(i+nts) = aimag(zt1)
      time(i) = tt0 + dts*real(i - 1)
  130 continue
! time display of raw data
      if (ntd > 0) then
         if (ndp==0) then
            label = ' POTENTIAL VERSUS TIME, RUNID='//cdrun
         else if (ndd==0) then
            label = ' ION DENSITY VERSUS TIME, RUNID='//cdrun
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         call DISPR(potc,label,time(1),time(ntd),999,0,0,ntd,nts,itwo,ch&
     &r(1:ncd),chrs,irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
! plot complex phase
         if (ndp==0) then
            label = ' REAL VERSUS IMAGINARY PART OF POTENTIAL'
         else if (ndd==0) then
            label = ' REAL VERSUS IMAGINARY PART OF ION DENSITY'
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         zsc = cmplx(999.,999.)
         zst = cmplx(0.,0.)
         call DISPC(potc,potc(nts1),label,zsc,zst,1,ntd,ntd,1,chr(1:ncd)&
     &,' ',irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
! log display of raw data
         do 140 i = 1, nts
         i1 = nts + i
         at1 = sqrt(potc(i)**2 + potc(i1)**2)
         if (at1 < pmin) at1 = pmin
         g(i) = alog(at1)
  140    continue
         if (ndp==0) then
            label = ' LN |POTENTIAL| VERSUS TIME, RUNID='//cdrun
         else if (ndd==0) then
            label = ' LN |ION DENSITY| VERSUS TIME, RUNID='//cdrun
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         call DISPS(g,label,time(1),time(ntd),999,2,ntd,chr(1:ncd),irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
      endif
      if (iw > 0) then
! perform frequency analysis of raw data
         call SPECT(potc,wm,p,tt0,dts,nts,iw,nts2,iw2)
! power spectrum of raw data
         if (ndp==0) then
            label = ' POTENTIAL VERSUS OMEGA, RUNID='//cdrun
         else if (ndd==0) then
            label = ' ION DENSITY VERSUS OMEGA, RUNID='//cdrun
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         call DISPR(p,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),c&
     &hrs(3),irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
      endif
      if (ntc==0) go to 160
! auto-correlation function
      call CORRE(potc,g,g,potc,mixupt,sctt,inft,nts,ntc,nts2,nft,nfth,nt&
     &c2)
      if (ntd > 0) then
! time display of correlated data
         it1 = ntd
         if (it1 > ntc) it1 = ntc
         if (ndp==0) then
            label = ' POT CORRELATION VERSUS TIME, RUNID='//cdrun
         else if (ndd==0) then
            label = ' ION CORRELATION VERSUS TIME, RUNID='//cdrun
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         call DISPR(potc,label,time(1),time(it1),999,0,0,it1,ntc,itwo,ch&
     &r(1:ncd),chrs,irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
! plot complex phase
         if (ndp==0) then
            label = ' REAL VERSUS IMAGINARY PART OF POT CORRELATION'
         else if (ndd==0) then
            label = ' REAL VERSUS IMAGINARY PART OF ION CORRELATION'
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         zsc = cmplx(999.,999.)
         zst = cmplx(0.,0.)
         call DISPC(potc,potc(ntc1),label,zsc,zst,1,it1,it1,1,chr(1:ncd)&
     &,' ',irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
! log display of correlated data
         do 150 i = 1, ntc
         i1 = ntc + i
         at1 = sqrt(potc(i)**2 + potc(i1)**2)
         if (at1 < pmin) at1 = pmin
         g(i) = alog(at1)
  150    continue
         if (ndp==0) then
            label = ' LN |POT CORRELATION| VERSUS TIME, RUNID='//cdrun
         else if (ndd==0) then
            label = ' LN |ION CORRELATION| VERSUS TIME, RUNID='//cdrun
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         call DISPS(g,label,time(1),time(it1),999,2,it1,chr(1:ncd),irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
      endif
      if (iw > 0) then
! perform frequency analysis of correlated data
         call SPECT(potc,wm,pc,tt0,dts,ntc,iw,ntc2,iw2)
! power spectrum of correlated data
         if (ndp==0) then
            label = ' POT CORRELATION VERSUS OMEGA, RUNID='//cdrun
         else if (ndd==0) then
            label = ' ION CORRELATION VERSUS OMEGA, RUNID='//cdrun
         endif
         write (chr,996) j1, k1, dkx, dky, ak
         ncd = 58
         call DISPR(pc,label,wm(1),wm(iw),999,1,0,iw,iw,itwo,chr(1:ncd),&
     &chrs(3),irc)
         if (irc==1) go to 180
         if (irc==3) go to 50
      endif
  160 continue
  170 continue
  180 call wtimer(times,ltime)
      go to 50
! close graphics device
  190 call GRCLOSE
      end program
