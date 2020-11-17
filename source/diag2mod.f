!-----------------------------------------------------------------------
!
      module diag2d
!
! Fortran90 interface to 2d PIC Fortran77 library diag2lib.f, dlibgks2.f
! diag2mod.f contains diagnostic procedures:
!            defines module diag2d
! wtimer performs wall clock timing.
! get_funit returns an unconnected fortran unit number.
! bfopen => bfopen2 opens binary file for real 2d scalar data.
! bfopen => bfcopen2 opens binary file for complex 2d scalar data.
! bfopen => bfvcopen2 opens binary file for complex 2d vector data.
! writebf => ifwrite2 writes real 2d scalar data to a binary file.
!            calls FWRITE2
! writebf => ifcwrite2 writes complex 2d scalar data to a binary file.
!            calls FCWRITE2
! writebf => ifvcwrite2 writes complex 2d vector data to a binary file.
!            call FCWRITE2
! readbf => ifread2 reads real 2d scalar data from a binary file.
!           calls FREAD2
! readbf => ifcread2 reads complex 2d scalar data from a binary file.
!           calls FCREAD2
! readbf => ifvcread2 reads complex 2d vector data from a binary file.
!           calls FCREAD2
! vdist => ivdist2 calculates 2 or 3 component velocity distribution and
!          velocity moments.
!          calls VDIST2 or VDIST23
! adist => iadist2 calculates 2 or 3 component field part of canonical
!          momentum and field momentum moments.
!          calls ADIST2 or ADIST23
! grasp => igrasp23 displays particle in phase space.
!          calls GRASP23
! displayfv => idisplayv2 displays velocity distribution functions.
!              calls DISPR
! displayfa => idisplaya2 displays field momentum distribution functions.
!              calls DISPR
! displays => idscaler2 displays 2d scalar field in real space.
!             calls CARPET, or CONTUR
! displays => idscaler1 displays 1d scalar field in real space.
!             calls DISPS
! displayv => idvector2 displays amplitude of a 2d vector field in real
!             space.
!             calls CARPET, or CONTUR
! displayw => idisplayw displays time history of electric field, kinetic
!             and total energies.
!             calls DISPR
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: august 14, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use init2d, only: idrun, indx, indy, ntp, ntd, nta, ntj, nte,     &
     &psolve, tend, dt, ndim, omx, omy, omz, ci, t0, ceng, indian,      &
     &rlprec, inprec, den2d, modesxd, modesyd, ndrec, fdname, pot2d,    &
     &modesxp, modesyp, nprec, fpname, vpot2d, modesxa, modesya, narec, &
     &faname, vcur2d, modesxj, modesyj, njrec, fjname, em2d, modesxe,   &
     &modesye, nerec, nlsrerec, fename
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC
      public :: GROPEN, SETNPLT, STPALIT, GRCLOSE
      public :: idrun, indx, indy, ntp, ntd, nta, ntj, nte, psolve
      public :: tend, dt, ndim, omx, omy, omz, ci, t0, ceng
      public :: indian, rlprec, inprec
      public :: den2d, modesxd, modesyd, ndrec, fdname
      public :: pot2d, modesxp, modesyp, nprec, fpname
      public :: vpot2d, modesxa, modesya, narec, faname
      public :: vcur2d, modesxj, modesyj, njrec, fjname
      public :: em2d, modesxe, modesye, nerec, fename
      public :: nlsrerec
      public :: wtimer, get_funit, grasp, vdist, adist
      public :: displayfv, displayfa, displays, displayv, displayw
      public :: bfopen, writebf, readbf
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine GROPEN
         implicit none
         end subroutine
      end interface
      interface
         subroutine SETNPLT(nplt,irc)
         implicit none
         integer :: nplt, irc
         end subroutine
      end interface
      interface
         subroutine STPALIT(idpal)
         implicit none
         integer :: idpal
         end subroutine
      end interface
      interface
         subroutine GRCLOSE
         implicit none
         end subroutine
      end interface
      interface
         subroutine GRASP23(part,label,itime,isc,nx,ny,iyp,ixp,idimp,npx&
     &y,np,irc)
         implicit none
         integer :: itime, isc, nx, ny, iyp, ixp, idimp, npxy, np
         integer :: irc
         character(len=*) :: label
         real, dimension(idimp,np) :: part
         end subroutine
      end interface
      interface
         subroutine VDIST2(part,fv,fvm,idimp,np,nmv,nmvf)
         integer :: idimp, np, nmv, nmvf
         real, dimension(idimp,np) :: part
         real, dimension(nmvf,2) :: fv
         real, dimension(3,2) :: fvm
         end subroutine
      end interface
      interface
         subroutine VDIST23(part,fv,fvm,idimp,np,nmv,nmvf)
         integer :: idimp, np, nmv, nmvf
         real, dimension(idimp,np) :: part
         real, dimension(nmvf,3) :: fv
         real, dimension(3,3) :: fvm
         end subroutine
      end interface
      interface
         subroutine ADIST2(part,av,avm,idimp,np,nma,nmaf)
         integer :: idimp, np, nma, nmaf
         real, dimension(idimp,np) :: part
         real, dimension(nmaf,2) :: av
         real, dimension(2,2) :: avm
         end subroutine
      end interface
      interface
         subroutine ADIST23(part,av,avm,idimp,np,nma,nmaf)
         integer :: idimp, np, nma, nmaf
         real, dimension(idimp,np) :: part
         real, dimension(nmaf,3) :: av
         real, dimension(2,3) :: avm
         end subroutine
      end interface
      interface
         subroutine DISPR(f,label,xmin,xmax,isc,ist,mks,nx,nxv,ngs,chr,c&
     &hrs,irc)
         implicit none
         integer :: isc, ist, mks, nx, nxv, ngs, irc
         real :: xmin, xmax
         character(len=*) :: label, chr
         character(len=*), dimension(ngs) :: chrs
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine CARPET(f,label,isc,ist,nx,ny,nxv,chr,ntc,irc)
         implicit none
         integer :: isc, ist, nx, ny, nxv, ntc, irc
         character(len=*) :: label, chr
!        real, dimension(*) :: f
         real :: f
         end subroutine
      end interface
      interface
         subroutine CONTUR(f,lf,label,isc,ist,nx,ny,nxv,chr,nc,irc)
         implicit none
         integer :: isc, ist, nx, ny, nxv, nc, irc
         character(len=*) label, chr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxv,ny) :: lf
         end subroutine
      end interface
      interface
         subroutine FWRITE2(f,nx,ny,nxv,iunit,nrec,lrec,name)
         implicit none
         integer :: nx, ny, nxv, iunit, nrec, lrec
!        real, dimension(*) :: f
         real :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine FREAD2(f,nx,ny,nxv,iunit,nrec,lrec,name,ierr)
         implicit none
         integer :: nx, ny, nxv, iunit, nrec, lrec, ierr
!        real, dimension(*) :: f
         real :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine FCWRITE2(f,nx,ny,nxv,iunit,nrec,lrec,name)
         implicit none
         integer :: nx, ny, nxv, iunit, nrec, lrec
!        complex, dimension(*) :: f
         complex :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine FCREAD2(f,nx,ny,nxv,iunit,nrec,lrec,name,ierr)
         implicit none
         integer :: nx, ny, nxv, iunit, nrec, lrec, ierr
!        complex, dimension(*) :: f
         complex :: f
         character(len=*) :: name
         end subroutine
      end interface
      interface
         subroutine FWRITE1(f,nxp,iunit,nrec,name)
         implicit none
         integer :: nxp, iunit, nrec
         real, dimension(nxp) :: f
         character(len=*) :: name
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface grasp
         module procedure igrasp23
      end interface
!
      interface vdist
         module procedure ivdist2
      end interface
!
      interface adist
         module procedure iadist2
      end interface
!
      interface displayfv
         module procedure idisplayv2
      end interface
!
      interface displayfa
         module procedure idisplaya2
      end interface
!
      interface displays
         module procedure idscaler2
         module procedure idscaler1
      end interface
!
      interface displayv
         module procedure idvector2
      end interface
!
      interface displayw
         module procedure idisplayw
      end interface
!
      interface bfopen
         module procedure bfopen2
         module procedure bfcopen2
         module procedure bfvcopen2
      end interface
!
      interface writebf
         module procedure ifwrite2
         module procedure ifcwrite2
         module procedure ifvcwrite2
         module procedure ifwrite1
      end interface
!
      interface readbf
         module procedure ifread2
         module procedure ifcread2
         module procedure ifvcread2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine wtimer(time,itime,icntrl)
! this subroutine performs wall clock timing
! input: itime, icntrl, output: itime, time
! itime = initial time, gets updated on each call
! icntrl = (-1,0,1) = (initialize,ignore,read) clock
!          clock should be initialized before it is read!
! time = elapsed time, in seconds,  since itime was set
         implicit none
         real, intent(out) :: time
         integer, intent(inout) :: itime
         integer, intent(in), optional :: icntrl
         integer :: ltime, COUNT_RATE, COUNT_MAX
         ltime = 1
         if (present(icntrl)) ltime = icntrl
         if (ltime==0) return
! read clock and write time difference from last clock initialization
         if (ltime==1) then
            ltime = itime
            call system_clock(itime,COUNT_RATE,COUNT_MAX)
            ltime = itime - ltime
            if (ltime.lt.0) ltime = ltime + COUNT_MAX
            time = dble(ltime)/dble(COUNT_RATE)
! initialize clock
         else
           call system_clock(ltime,COUNT_RATE,COUNT_MAX)
           itime = ltime
         endif
         end subroutine wtimer
!
         function get_funit(start) result(funit)
! this function returns an unconnected fortran unit number,
! starting with unit = start.  returns -1 if none found
         integer, intent(in) :: start
         integer :: funit
! local data
         integer :: i
         logical :: connected
         funit = -1
! check connection status
         do i = start, 99
            inquire(unit=i,opened=connected)
            if (.not.connected) then
               funit = i
               exit
            endif
         enddo
         end function get_funit
!
         subroutine igrasp23(part,np,label,itime,isc,nx,ny,iyp,ixp,npxy,&
     &irc)
! displays phase space
         implicit none
         integer :: np, itime, isc, nx, ny, iyp, ixp, npxy, irc
         character(len=*) :: label
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call GRASP23(part,label,itime,isc,nx,ny,iyp,ixp,idimp,npxy,np,i&
     &rc)
         end subroutine igrasp23
!
         subroutine ivdist2(part,fv,fvm,np,nmv)
! calculates 2d velocity distribution, velocity moments, and entropy
         integer :: np, nmv
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: fv, fvm
! local data
         integer :: idimp, nmvf, idimv
         idimp = size(part,1)
         nmvf = size(fv,1); idimv = size(fv,2)
         if (idimv==2) then
            call VDIST2(part,fv,fvm,idimp,np,nmv,nmvf)
         else if (idimv==3) then
            call VDIST23(part,fv,fvm,idimp,np,nmv,nmvf)
         endif
         end subroutine ivdist2
!
         subroutine iadist2(part,av,avm,np,nma)
! calculates 2d  distribution of field part of canonical momentum
! and field momentum moments
         integer :: np, nma
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: av, avm
! local data
         integer :: idimp, nmaf, idima
         idimp = size(part,1)
         nmaf = size(av,1); idima = size(av,2)
         if (idima==2) then
            call ADIST2(part,av,avm,idimp,np,nma,nmaf)
         else if (idima==3) then
            call ADIST23(part,av,avm,idimp,np,nma,nmaf)
         endif
         end subroutine iadist2
!
         subroutine idisplayv2(fv,fvm,label,itime,nmv,idt,irc)
! displays velocity distribution functions
! fv = velocity distribution
! fvm = velocity moments
! itime = current time step
! nmv = number of velocity intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, nmv, idt, irc
         real, dimension(:,:), pointer :: fv, fvm
         character(len=*) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 1, mks = 0
         integer :: i, nmvf, nmv2, idimv
         real :: vmax, vmin
         character(len=12) :: c
         character(len=2) :: cs
         character(len=54) :: lbl
         character(len=45) :: chr
         character(len=10), dimension(3) :: chrs
   91    format(', T =',i7)
   92    format(' VD =',f9.6,' VTH =',f9.5)
   93    format(' VTX =',f9.5,' VTY =',f9.5)
   94    format(' VTX =',f9.5,' VTY =',f9.5,' VTZ =',f9.5)
! chrs = short array of characters to label individual line samples
         data chrs /'    VX    ','    VY    ','    VZ    '/
         idimv = size(fv,2)
         if ((idimv /= 2) .and. (idimv /= 3)) return
         nmvf = size(fv,1)
         nmv2 = 2*nmv + 1
         write (c,91) itime
! each velocity distributions on its own plot
         if (idt /= 2) then
            do i = 1, idimv
            cs = trim(adjustl(chrs(i)))
            lbl = trim(label)//' VELOCITY DISTRIBUTION VERSUS '//cs//c
            write (chr,92) fvm(1,i), fvm(2,i)
            vmax = fv(1,i)
            vmin = -vmax
            call DISPR(fv(2,i),lbl,vmin,vmax,isc,ist,mks,nmv2,nmvf,1,chr&
     &,chrs(i),irc)
            if (irc==1) return
            enddo
         endif
! all velocity distributions on common plot
         if (idt /= 1) then
            lbl = trim(label)//' VELOCITY DISTRIBUTIONS VERSUS '//'V'//c
            if (idimv==2) then
               write (chr,93) fvm(2,1), fvm(2,2)
               vmax = max(fv(1,1),fv(1,2))
               vmin = -vmax
            else if (idimv==3) then
               write (chr,94) fvm(2,1), fvm(2,2), fvm(2,3)
               vmax = max(fv(1,1),fv(1,2),fv(1,3))
               vmin = -vmax
            endif
            call DISPR(fv(2,1),lbl,vmin,vmax,isc,ist,mks,nmv2,nmvf,idimv&
     &,chr,chrs,irc)
         endif
         end subroutine idisplayv2
!
         subroutine idisplaya2(fa,fam,label,itime,nma,idt,irc)
! displays field momentum distribution functions
! fa = field momentum distribution
! fam = field momentum moments
! itime = current time step
! nma = number of field momentum intervals
! idt = (1,2,3) = display (individual,composite,both) functions
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, nma, idt, irc
         real, dimension(:,:), pointer :: fa, fam
         character(len=*) :: label
! isc = 999 = use default display scale
! ist = 1 = display positive values only
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 1, mks = 0
         integer :: i, nmaf, nma2, idima
         real :: amax, amin
         character(len=12) :: c
         character(len=2) :: cs
         character(len=54) :: lbl
         character(len=45) :: chr
         character(len=10), dimension(3) :: chrs
   91    format(', T =',i7)
   92    format(' AD =',f9.6,' ATH =',f9.5)
   93    format(' ATX =',f9.5,' ATY =',f9.5)
   94    format(' ATX =',f9.5,' ATY =',f9.5,' ATZ =',f9.5)
! chrs = short array of characters to label individual line samples
         data chrs /'    AX    ','    AY    ','    AZ    '/
         idima = size(fa,2)
         if ((idima /= 2) .and. (idima /= 3)) return
         nmaf = size(fa,1)
         nma2 = 2*nma + 1
         write (c,91) itime
! each field momentum distributions on its own plot
         if (idt /= 2) then
            do i = 1, idima
            cs = trim(adjustl(chrs(i)))
            lbl = trim(label)//' FIELD DISTRIBUTION VERSUS '//cs//c
            write (chr,92) fam(1,i), fam(2,i)
            amax = fa(1,i)
            amin = -amax
            call DISPR(fa(2,i),lbl,amin,amax,isc,ist,mks,nma2,nmaf,1,chr&
     &,chrs(i),irc)
            if (irc==1) return
            enddo
         endif
! all velocity distributions on common plot
         if (idt /= 1) then
            lbl = trim(label)//' FIELD DISTRIBUTIONS VERSUS '//'V'//c
            if (idima==2) then
               write (chr,93) fam(2,1), fam(2,2)
               amax = max(fa(1,1),fa(1,2))
               amin = -amax
            else if (idima==3) then
               write (chr,94) fam(2,1), fam(2,2), fam(2,3)
               amax = max(fa(1,1),fa(1,2),fa(1,3))
               amin = -amax
            endif
            call DISPR(fa(2,1),lbl,amin,amax,isc,ist,mks,nma2,nmaf,idima&
     &,chr,chrs,irc)
         endif
         end subroutine idisplaya2
!
         subroutine idscaler1(f,label,itime,isc,ist,nx,irc,inorder)
! displays 1d scalar field in real space
! f = 1d scalar field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of pot
! ist = flag for choosing positive and/or negative values
! nx = system length in x direction
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, isc, ist, nx, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:), pointer :: f
! local data
         integer :: nxv, mx, order
         real :: xmin, xmax
         character(len=12) :: lbl
   91    format(' T = ',i7)
         nxv = size(f)
         mx = nx
         if ((mx+1) <= nxv) mx = mx + 1
         xmin = 0.0; xmax = real(mx - 1)
         write (lbl,91) itime
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call DISPS(f(1),label,xmin,xmax,isc,ist,mx,lbl,irc)
         else if (order==CUBIC) then
            nxv = nxv - 2
            call DISPS(f(3),label,xmin,xmax,isc,ist,mx,lbl,irc)
         else
            nxv = nxv - 1
            call DISPS(f(2),label,xmin,xmax,isc,ist,mx,lbl,irc)
         endif
         end subroutine idscaler1
!
         subroutine idscaler2(pot,label,itime,isc,ist,idt,nx,ny,irc,inor&
     &der)
! displays 2d scalar field in real space
! pot = 2d scalar field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of pot
! ist = flag for choosing positive and/or negative values
! idt = (1,2,3) = display (color map,contour plot,both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, isc, ist, idt, nx, ny, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:,:) :: pot
! local data
         integer, dimension(size(pot,1),size(pot,2)) :: lf
! ntc = number of valid colors, should be power of 2, <= 256
! nc = number of contour lines
         integer :: ntc = 16, nc = 16
         integer :: nxv, mx, my, order
         character(len=12) :: lbl
   91    format(' T = ',i7)
         nxv = size(pot,1)
         mx = nx; my = ny
! plot guard cells if present
         if ((mx+1) <= nxv) mx = mx + 1
         if ((my+1) <= size(pot,2)) my = my + 1
         write (lbl,91) itime
         order = QUADRATIC
         if (present(inorder)) order = inorder
! color map plot for all values
         if (idt /= 2) then
            if (order==LINEAR) then
               call CARPET(pot(1,1),label,isc,ist,mx,my,nxv,lbl,ntc,irc)
            else if (order==CUBIC) then
               call CARPET(pot(3,3),label,isc,ist,mx,my,nxv,lbl,ntc,irc)
            else
               call CARPET(pot(2,2),label,isc,ist,mx,my,nxv,lbl,ntc,irc)
            endif
         endif
! contour map for all values
         if (idt /= 1) then
            if (order==LINEAR) then
               call CONTUR(pot(1,1),lf,label,isc,ist,mx,my,nxv,lbl,nc,ir&
     &c)
            else if (order==CUBIC) then
               call CONTUR(pot(3,3),lf,label,isc,ist,mx,my,nxv,lbl,nc,ir&
     &c)
            else
               call CONTUR(pot(2,2),lf,label,isc,ist,mx,my,nxv,lbl,nc,ir&
     &c)
            endif
         endif
         end subroutine idscaler2
!
         subroutine idvector2(vpot,label,itime,isc,ist,idt,nx,ny,irc,ino&
     &rder)
! displays amplitude of a 2d vector field in real space
! vpot = 2d vector field in real space
! label = field label
! itime = current time step
! isc = power of 2 scale of range of values of pot
! ist = flag for choosing positive and/or negative values
! idt = (1,2,3) = display (color map,contour plot,both)
! nx/ny = system length in x/y direction
! irc = return code (0 = normal return)
         implicit none
         integer :: itime, isc, ist, idt, nx, ny, irc
         integer, optional :: inorder
         character(len=*) label
         real, dimension(:,:,:), pointer :: vpot
! local data
         integer :: i, j, k, mx, my, order
         real :: sum1
         real, dimension(size(vpot,2),size(vpot,3)) :: pot
! calculate array size with guard cells
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            mx = nx + 1; my = ny + 1
         else if (order==CUBIC) then
            mx = nx + 5; my = ny + 5
         else
            mx = nx + 3; my = ny + 3
         endif
         mx = min(mx,size(vpot,2))
         my = min(my,size(vpot,3))
! find absolute value of vector field
         do k = 1, my
         do j = 1, mx
            sum1 = 0.0
            do i = 1, size(vpot,1)
            sum1 = sum1 + vpot(i,j,k)**2
            enddo
            pot(j,k) = sqrt(sum1)
         enddo
         enddo
! display amplitude
         call idscaler2(pot,label,itime,isc,ist,idt,nx,ny,irc,inorder)
         end subroutine idvector2
!
         subroutine idisplayw(wt,t0,dtw,nt,irc)
! displays time history of electric field, kinetic, and
! total energies
! wt = time history array for energies
! t0 = initial energy time
! dtw = time between energy values
! nt = number of energy values to be displayed
! irc = return code (0 = normal return)
         integer :: nt, irc
         real :: t0, dtw
         real, dimension(:,:), pointer :: wt
! isc = 999 = use default display scale
! ist = 2 = display minimum range
! mks = 0 = cycle through line styles
         integer :: isc = 999, ist = 2, mks = 0
         integer :: i, ntwd, ns
         real :: tmin, tmax
         character(len=36) :: lbl
         character(len=20), dimension(7) :: cs 
         character(len=10), dimension(7) :: chrs
! chrs = short array of characters to label individual line samples
         data cs /' TOTAL FIELD        ',' ELECTRON KINETIC   ',' ION KI&
     &NETIC     ',' TOTAL              ',' ES FIELD           ',' ET FIE&
     &LD        ',' MAGNETIC FIELD     '/
         data chrs /'TOT FIELD ','ELECT ENRG',' ION ENRG ','TOTAL ENRG',&
     &' EL FIELD ',' ET FIELD ',' B FIELD  '/
! quit if array is empty or incorrect
         if (nt <= 0) return
         ntwd = size(wt,1)
         ns = min(size(wt,2),7)
! tmin/tmax = range of time values in plot
         tmin = t0
         tmax = t0 + dtw*(nt - 1)
! display individual energies
         do i = 1, ns
         lbl = trim(cs(i))//' ENERGY VERSUS TIME'
         call DISPR(wt(1,i),lbl,tmin,tmax,isc,ist,mks,nt,ntwd,1,' ',chrs&
     &(i),irc)
         if (irc==1) return
         enddo
! all energies on common plot
         lbl = ' ENERGIES VERSUS TIME'
         call DISPR(wt(1,1),lbl,tmin,tmax,isc,ist,mks,nt,ntwd,ns,' ',chr&
     &s,irc)
         end subroutine idisplayw
!
         subroutine bfopen2(f,nx,ny,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for real 2d scalar data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, ny, iunit, nrec
         real, dimension(:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1,1)
         lrec = nx*ny*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfopen2
!
         subroutine ifwrite2(f,nx,ny,iunit,nrec,name,order)
! writes a subset of real 2d scalar data to a direct access binary file
         implicit none
         integer :: nx, ny, iunit, nrec
         real, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nx*ny*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         nxv = size(f,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FWRITE2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FWRITE2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,name)
            else
               call FWRITE2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FWRITE2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FWRITE2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,noname)
            else
               call FWRITE2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifwrite2
!
         subroutine ifread2(f,nx,ny,iunit,nrec,ierr,name,order)
! reads a subset of real 2d scalar data from a direct access binary file
         implicit none
         integer :: nx, ny, iunit, nrec, ierr
         real, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nx*ny*lrec
         endif
         nxv = size(f,1)
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FREAD2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,name,ierr)
            else if (inorder==CUBIC) then
               call FREAD2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,name,ierr)
            else
               call FREAD2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call FREAD2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,noname,ierr)
            else if (inorder==CUBIC) then
               call FREAD2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,noname,ierr)
            else
               call FREAD2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,noname,ierr)
            endif
         endif
         end subroutine ifread2
!
         subroutine bfcopen2(f,nx,ny,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 2d scalar data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, ny, iunit, nrec
         complex, dimension(:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, ierr
         if (nrec > 0) return
         inquire(iolength=lrec) f(1,1)
         lrec = nx*ny*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfcopen2
!
         subroutine ifcwrite2(f,nx,ny,iunit,nrec,name,order)
! writes a subset of complex 2d scalar data to a direct access binary
! file
         implicit none
         integer :: nx, ny, iunit, nrec
         complex, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         nxv = size(f,1)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nx*ny*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCWRITE2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FCWRITE2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,name)
            else
               call FCWRITE2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FCWRITE2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FCWRITE2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,noname)
            else
               call FCWRITE2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifcwrite2
!
         subroutine ifcread2(f,nx,ny,iunit,nrec,ierr,name,order)
! reads a subset of complex 2d scalar data from a direct access binary
! file
         implicit none
         integer :: nx, ny, iunit, nrec, ierr
         complex, dimension(:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: lrec, nxv, inorder
         character(len=1) :: noname = ' '
         nxv = size(f,1)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1)
            lrec = nx*ny*lrec
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCREAD2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,name,ierr)
            else if (inorder==CUBIC) then
               call FCREAD2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,name,ierr)
            else
               call FCREAD2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,name,ierr)
            endif
         else
            if (inorder==LINEAR) then
               call FCREAD2(f(1,1),nx,ny,nxv,iunit,nrec,lrec,noname,ierr&
     &)
            else if (inorder==CUBIC) then
               call FCREAD2(f(3,3),nx,ny,nxv,iunit,nrec,lrec,noname,ierr&
     &)
            else
               call FCREAD2(f(2,2),nx,ny,nxv,iunit,nrec,lrec,noname,ierr&
     &)
            endif
         endif
         end subroutine ifcread2
!
         subroutine bfvcopen2(f,nx,ny,iunit,nrec,fname)
! this subroutine opens direct access binary file
! for complex 2d vector data
! nrec = (0,-1) open (old, new file), reset to 1 if successful
         implicit none
         integer :: nx, ny, iunit, nrec
         complex, dimension(:,:,:), pointer :: f
         character(len=*) :: fname
! local data
         integer :: lrec, nnx, ierr
         if (nrec > 0) return
         nnx = size(f,1)*nx
         inquire(iolength=lrec) f(1,1,1)
         lrec = nnx*ny*lrec
         if (nrec==0) then
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='old',iostat=ierr)
         else
            open(unit=iunit,file=fname,form='unformatted',access='direct&
     &',recl=lrec,status='replace',iostat=ierr)
         endif
         if (ierr==0) nrec = 1
         end subroutine bfvcopen2
!
         subroutine ifvcwrite2(f,nx,ny,iunit,nrec,name,order)
! writes a subset of complex 2d vector data to a direct access binary
! file
         implicit none
         integer :: nx, ny, iunit, nrec
         complex, dimension(:,:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: nnx, nxv, lrec, inorder
         character(len=1) :: noname = ' '
         nnx = size(f,1)*nx
         nxv = size(f,1)*size(f,2)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = nnx*ny*lrec
            if (nrec < 0) iunit = get_funit(iunit)
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCWRITE2(f(1,1,1),nnx,ny,nxv,iunit,nrec,lrec,name)
            else if (inorder==CUBIC) then
               call FCWRITE2(f(1,3,3),nnx,ny,nxv,iunit,nrec,lrec,name)
            else
               call FCWRITE2(f(1,2,2),nnx,ny,nxv,iunit,nrec,lrec,name)
            endif
         else
            if (inorder==LINEAR) then
               call FCWRITE2(f(1,1,1),nnx,ny,nxv,iunit,nrec,lrec,noname)
            else if (inorder==CUBIC) then
               call FCWRITE2(f(1,3,3),nnx,ny,nxv,iunit,nrec,lrec,noname)
            else
               call FCWRITE2(f(1,2,2),nnx,ny,nxv,iunit,nrec,lrec,noname)
            endif
         endif
         end subroutine ifvcwrite2
!
         subroutine ifvcread2(f,nx,ny,iunit,nrec,ierr,name,order)
! reads a subset of complex 2d vector data from a direct access binary
! file
         implicit none
         integer :: nx, ny, iunit, nrec, ierr
         complex, dimension(:,:,:), pointer :: f
         character(len=*), optional :: name
         integer, optional :: order
! local data
         integer :: nnx, nxv, lrec, inorder
         character(len=1) :: noname = ' '
         nnx = size(f,1)*nx
         nxv = size(f,1)*size(f,2)
         if (nrec <= 0) then
            inquire(iolength=lrec) f(1,1,1)
            lrec = nnx*ny*lrec
         endif
         inorder = QUADRATIC
         if (present(order)) inorder = order
         if (present(name)) then
            if (inorder==LINEAR) then
               call FCREAD2(f(1,1,1),nnx,ny,nxv,iunit,nrec,lrec,name,ier&
     &r)
            else if (inorder==CUBIC) then
               call FCREAD2(f(1,3,3),nnx,ny,nxv,iunit,nrec,lrec,name,ier&
     &r)
            else
               call FCREAD2(f(1,2,2),nnx,ny,nxv,iunit,nrec,lrec,name,ier&
     &r)
            endif
         else
            if (inorder==LINEAR) then
               call FCREAD2(f(1,1,1),nnx,ny,nxv,iunit,nrec,lrec,noname,i&
     &err)
            else if (inorder==CUBIC) then
               call FCREAD2(f(1,3,3),nnx,ny,nxv,iunit,nrec,lrec,noname,i&
     &err)
            else
               call FCREAD2(f(1,2,2),nnx,ny,nxv,iunit,nrec,lrec,noname,i&
     &err)
            endif
         endif
         end subroutine ifvcread2
!
         subroutine ifwrite1(f,iunit,nrec,name)
! writes real 2d scalar data to a direct access binary file
         implicit none
         integer :: iunit, nrec
         real, dimension(:), pointer :: f
         character(len=*) :: name
! local data
         integer :: nxp
         nxp = size(f)
         call FWRITE1(f,nxp,iunit,nrec,name)
         end subroutine ifwrite1
!
      end module diag2d
