!-----------------------------------------------------------------------
!
      module simul2d
! Higher level subroutines for electrostatics
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: november 14, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use diag2d, only: get_funit, bfopen, vdist, writebf, displays,    &
     &displayfv, grasp
      use espush2d, only: dpost, push, rpush, pushzf, rpushzf, pushglx, &
     &rpushglx, pushgl, rpushgl, premoment2, primoment2, fft
      use field2d, only: sguard, aguard, cguard, spois, pois, gtmodes,  &
     &imoment
      implicit none
      private
      public :: restart_open, restart_bwrite, restart_dwrite
      public :: restart_bread, restart_dread
      public :: dpostg, pushg, bpushg, pushzfg, pushglg, pushglxg
      public :: initmodediag, initveldiag
      public :: dendiag, veldiag, phasediag, potdiag
      public :: emomtdiag, imomtdiag, esenergy
!
      contains
!
         subroutine restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
! open restart files
         implicit none
         integer, intent(in) :: nustrt, ntr, idrun
         integer :: iur1, iur2, iuer
! local data
         integer :: ierr
         character(len=10) :: cdrun
         character(len=32) :: fname
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
         iur1 = 0; iur2 = 0
! open old restart files
         if ((nustrt==2).or.(nustrt==0)) then
            iur1 = get_funit(16)
            fname = 'rstrt1.'//cdrun
            open(unit=iur1,file=trim(fname),form='unformatted',status='o&
     &ld',iostat=ierr)
            if (ierr /= 0)  then
               iur1 = -1
               write (iuer,*) 'Cannot open restart file1=',fname
            endif
            iur2 = get_funit(17)
            fname = 'rstrt2.'//cdrun
            open(unit=iur2,file=trim(fname),form='unformatted',status='o&
     &ld',iostat=ierr)
            if (ierr /= 0)  then
               iur2 = -1
               write (iuer,*) 'Cannot open restart file2=',fname
            endif
! open new restart files
         else if (ntr > 0) then
            iur1 = get_funit(16)
            fname = 'rstrt1.'//cdrun
            open(unit=iur1,file=trim(fname),form='unformatted',status='u&
     &nknown')
            iur2 = get_funit(17)
            fname = 'rstrt2.'//cdrun
            open(unit=iur2,file=trim(fname),form='unformatted',status='u&
     &nknown')
         endif
         end subroutine restart_open
!
         subroutine restart_bwrite(kunit,itime,itime0,np,part,movion,npi&
     &,parti,qi,ef,bf)
! write file for basic restart or continuation
         implicit none
         integer, intent(in) :: kunit, itime, itime0, np, movion, npi
         real, dimension(:,:), pointer :: part, parti
         real, dimension(:,:), pointer :: qi
         complex, dimension(:,:,:), pointer, optional :: ef, bf
! local data
         integer :: j, it
         character(len=1), dimension(8), save :: arch = ' '
         integer, external :: NDIAN, NDPREC, IDPREC
! determine architecture information
         if (arch(1)==' ') then
! determine if number format is big or little endian
            it = NDIAN()
            if (it==0) then
               arch(1) = 'L'
            else if (it==1) then
               arch(1) = 'B'
            endif
! determine if default reals are double precision
            it = NDPREC()
            if (it==0) then
               arch(2) = 'S'
            else if (it==1) then
               arch(2) = 'D'
            endif
! determine if default integers are double precision
            it = IDPREC()
            if (it==0) then
               arch(3) = 'S'
            else if (it==1) then
               arch(3) = 'D'
            endif
         endif
! write out architecture information
         write (kunit) (arch(j),j=1,8)
! write out current and initial time
         write (kunit) itime, itime0
! write out number of electrons, size of particle array
         write (kunit) np, size(part,1)
! write out electrons, if non-zero
         if (np > 0) write (kunit) (part(:,j),j=1,np)
! write out if ions are moving
         write (kunit) movion
         if (movion > 0) then
! write out number of ions, size of particle array
            write (kunit) npi, size(parti,1)
! write out ions, if non-zero
            if (npi > 0) write (kunit) (parti(:,j),j=1,npi)
         else if (movion==0) then
! write out ion density, if ions are not moving
            write (kunit) size(qi)
            write (kunit) qi
         endif
! write out electromagnetic fields, if present
         it = 0
         if (present(ef)) then
            write (kunit) size(ef)
            write (kunit) ef
         else
            write (kunit) it
         endif
         if (present(bf)) then
            write (kunit) size(bf)
            write (kunit) bf
         else
            write (kunit) it
         endif
! write out current time for later confirmation
         write (kunit) itime
         end subroutine restart_bwrite
!
         subroutine restart_dwrite(kunit,itime,itw,wt,ndrec,fdname,nprec&
     &,fpname,narec,faname,njrec,fjname,nerec,fename)
! write diagnostic file for restart or continuation
         implicit none
         integer, intent(in) :: kunit, itime
         integer :: itw, ndrec, nprec
         real, dimension(:,:), pointer :: wt
         character(len=*), intent(in) :: fdname, fpname
         integer, optional :: narec, njrec, nerec
         character(len=*), intent(in), optional :: faname, fjname
         character(len=*), intent(in), optional :: fename
! local data
         integer :: it, na, nj, ne, irc
         character(len=32) :: fname
         na = 0; nj = 0; ne = 0
         if ((present(narec)).and.(present(faname))) na = 1
         if ((present(njrec)).and.(present(fjname))) nj = 1
         if ((present(nerec)).and.(present(fename))) ne = 1
! write out number of diagnostics in file
         it = 6
         write (kunit) it
! write out current energy time step
         write (kunit) itw
         if (itw > 0) then
! write out size of energy array
            write (kunit) size(wt,2)
! write out energy values
            write (kunit) wt(1:itw,:)
         endif
! write out ion density diagnostic write location
         write (kunit) ndrec
! write out record length
         if (ndrec > 0) then
            inquire(file=fdname,recl=it,iostat=irc)
            if (irc /= 0) it = 0
            write (kunit) it
            if (it > 0) then
               fname = fdname
               write (kunit) fname
            endif
         endif
! write out potential diagnostic write location
         write (kunit) nprec
! write out record length
         if (nprec > 0) then
            inquire(file=fpname,recl=it,iostat=irc)
            if (irc /= 0) it = 0
            write (kunit) it
            if (it > 0) then
               fname = fpname
               write (kunit) fname
            endif
         endif
! write out vector potential diagnostic write location
         if (na==1) then
            write (kunit) narec
! write out record length
            if (narec > 0) then
               irc = 0
               inquire(file=faname,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = faname
                  write (kunit) fname
               endif
            endif
         else
            write (kunit) na
         endif
! write out ion current diagnostic write location
         if (nj==1) then
            write (kunit) njrec
! write out record length
            if (njrec > 0) then
               irc = 0
               inquire(file=fjname,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = fjname
                  write (kunit) fname
               endif
            endif
         else
            write (kunit) nj
         endif
! write out electromagnetic diagnostic write location
         if (ne==1) then
            write (kunit) nerec
! write out record length
            if (nerec > 0) then
               irc = 0
               inquire(file=fename,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = fename
                  write (kunit) fname
               endif
            endif
         else
            write (kunit) ne
         endif
! write current time for later confirmation
         write (kunit) itime
         end file kunit
         rewind kunit
         end subroutine restart_dwrite
!
         subroutine restart_bread(kunit,itime,itime0,np,part,movion,npi,&
     &parti,qi,irc,iuer,ef,bf)
! read file for basic restart or continuation
         implicit none
         integer, intent(in) :: kunit, movion, iuer
         integer :: itime, itime0, np, npi, irc
         real, dimension(:,:), pointer :: part, parti
         real, dimension(:,:), pointer :: qi
         complex, dimension(:,:,:), pointer, optional :: ef, bf
! local data
         integer :: j, it
         character(len=1), dimension(8), save :: arch = ' '
         integer, external :: NDIAN, NDPREC, IDPREC
! check if unit number is valid
         irc = 99; if (kunit < 0) return
! read in architecture information
         irc = 1; read (kunit,err=10) (arch(j),j=1,8)
! determine if number format is big or little endian
         it = NDIAN()
         if ((it==0) .and. (arch(1)/='L')) go to 10
         if ((it==1) .and. (arch(1)/='B')) go to 10
! determine if default reals are double precision
         it = NDPREC()
         if ((it==0) .and. (arch(2)/='S')) go to 10
         if ((it==1) .and. (arch(2)/='D')) go to 10
! determine if default integers are double precision
         it = IDPREC()
         if ((it==0) .and. (arch(3)/='S')) go to 10
         if ((it==1) .and. (arch(3)/='D')) go to 10
! read in current time
         irc = 2; read (kunit,err=10) itime, itime0
! read in number of electrons, size of particle array
         irc = 3; read (kunit,err=10) np, it
! read in electrons, if non-zero
         if (np > 0) then
            if ((it > size(part,1)) .or. (np > size(part,2))) go to 10
            irc = 4; read (kunit,err=10) (part(1:it,j),j=1,np)
         endif
! read in if ions are moving
         irc = 5; read (kunit,err=10) it
         if (it /= movion) go to 10
         if (movion > 0) then
! read in number of ions, size of particle array
            irc = 6; read (kunit,err=10) npi, it
! read in ions, if non-zero
            if (npi > 0) then
               if ((it>size(parti,1)) .or. (npi>size(parti,2))) go to 10
               irc = 7; read (kunit,err=10) (parti(1:it,j),j=1,npi)
            endif
! read in ion density, if ions are not moving
         else if (movion==0) then
            irc = 8; read (kunit,err=10) it
            if (it > size(qi)) go to 10
            irc = 9; read (kunit,err=10) qi
         endif
! read in first electromagnetic field, if present
         irc = 10; read (kunit,err=10) it
         if (it > 0) then
            if (.not.present(ef)) go to 10
            if (it > size(ef)) go to 10
            irc = 11; read (kunit,err=10) ef
         endif
! read in second electromagnetic field, if present
         irc = 12; read (kunit,err=10) it
         if (it > 0) then
            if (.not.present(bf)) go to 10
            if (it > size(bf)) go to 10
            irc = 13; read (kunit,err=10) bf
         endif
! read in current time for confirmation
         irc = 14; read (kunit,err=10) it
         if (it /= itime) go to 10
         irc = 0
         return
! write out errors
   10    write (iuer,*) 'Basic Restart Error, irc = ', irc
         if (irc==1) then
            write (iuer,*) 'Architecture=', (arch(j),j=1,8)
         else if (irc==3) then
            write (iuer,*) 'np,size(part,1)=', np, it
         else if (irc==5) then
            write (iuer,*) 'movion=', it
         else if (irc==6) then
            write (iuer,*) 'npi,size(parti,1)=', npi, it
         else if (irc==8) then
            write (iuer,*) 'size(qi)=', it
         else if (irc==10) then
            write (iuer,*) 'size(ef)=', it
         else if (irc==12) then
            write (iuer,*) 'size(bf)=', it
         else if (irc==14) then
            write (iuer,*) 'confirmation itime,it=', itime, it
         endif
         end subroutine restart_bread
!
         subroutine restart_dread(kunit,itime,itw,wt,iud,ndrec,fdname,iu&
     &p,nprec,fpname,irc,iuer,iua,narec,faname,iuj,njrec,fjname,iue,nere&
     &c,fename)
! read diagnostic file for restart or continuation
         implicit none
         integer, intent(in) :: kunit, itime
         integer :: itw, iud, ndrec, iup, nprec, irc, iuer
         real, dimension(:,:), pointer :: wt
         character(len=*) :: fdname, fpname
         integer, optional :: iua, narec, iuj, njrec, iue, nerec
         character(len=*), optional :: faname, fjname, fename
! local data
         integer :: it, nt, na, nj, ne, ierr
         character(len=32) :: fname
         na = 0; nj = 0; ne = 0
         if ((present(iua)).and.(present(narec)).and.(present(faname))) &
     &na = 1
         if ((present(iuj)).and.(present(njrec)).and.(present(fjname))) &
     &nj = 1
         if ((present(iue)).and.(present(nerec)).and.(present(fename))) &
     &ne = 1
! check if unit number is valid
         irc = -99; if (kunit < 0) return
! read in number of diagnostics in file
         irc = -98; read (kunit,err=20) nt
! read in current energy time step
         irc = -1; read (kunit,err=20) itw
         if (itw > 0) then
! read in size of energy array
            irc = -2; read (kunit,err=20) it
            if ((itw > size(wt,1)) .or. (it > size(wt,2))) go to 20
! read in energy values
            irc = -3; read (kunit,err=20) wt(1:itw,1:it)
         endif
         if (nt==1) go to 10
! read in density diagnostic write location
         irc = -4; read (kunit,err=20) ndrec
! read in record length and open file
         if (ndrec > 0) then
            irc = -5; read (kunit,err=20) it
            if (it < 1) go to 20
            irc = -6; read (kunit,err=20) fname
            fdname = fname
            irc = -7;  iud = get_funit(iud)
            open(unit=iud,file=fdname,form='unformatted',access='direct'&
     &,recl=it,status='old',iostat=ierr)
            if (ierr /= 0) go to 20
         endif
         if (nt==2) go to 10
! read in potential diagnostic write location
         irc = -8; read (kunit,err=20) nprec
! read in record length and open file
         if (nprec > 0) then
            irc = -9; read (kunit,err=20) it
            if (it < 1) go to 20
            irc = -10; read (kunit,err=20) fname
            fpname = fname
            irc = -11;  iup = get_funit(iup)
            open(unit=iup,file=fpname,form='unformatted',access='direct'&
     &,recl=it,status='old',iostat=ierr)
            if (ierr /= 0) go to 20
         endif
         if (nt==3) go to 10
! read in vector potential diagnostic write location
         irc = -12; read (kunit,err=20) it
         if (na==1) then
            narec = it
! read in record length and open file
            if (narec > 0) then
               irc = -13; read (kunit,err=20) it
               if (it < 1) go to 20
               irc = -14; read (kunit,err=20) fname
               faname = fname
               irc = -15;  iua = get_funit(iua)
               open(unit=iua,file=faname,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ierr)
               if (ierr /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
         if (nt==4) go to 10
! read in ion current diagnostic write location
         irc = -16; read (kunit,err=20) it
         if (nj==1) then
            njrec = it
! read in record length and open file
            if (njrec > 0) then
               irc = -17; read (kunit,err=20) it
               if (it < 1) go to 20
               irc = -18; read (kunit,err=20) fname
               fjname = fname
               irc = -19;  iuj = get_funit(iuj)
               open(unit=iuj,file=fjname,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ierr)
               if (ierr /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
         if (nt==5) go to 10
! read in electromagnetic diagnostic write location
         irc = -20; read (kunit,err=20) it
         if (ne==1) then
            nerec = it
! read in record length and open file
            if (nerec > 0) then
               irc = -21; read (kunit,err=20) it
               if (it < 1) go to 20
               irc = -22; read (kunit,err=20) fname
               fename = fname
               irc = -23;  iue = get_funit(iue)
               open(unit=iue,file=fename,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ierr)
               if (ierr /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
! read in current time for confirmation
   10    irc = -24; read (kunit,err=20) it
         if (it /= itime) go to 20
         irc = 0
         rewind kunit
         return
! write out errors
   20    write (iuer,*) 'Diagnostic Restart Error, irc = ', irc
         if (irc==(-1)) then
            write (iuer,*) 'itw=', itw
         else if (irc==(-2)) then
            write (iuer,*) 'itw,size(wt,2)=', itw, it
         else if ((irc==(-5)).or.(irc==(-9)).or.(irc==(-13)).or.(irc==(-&
     &17)).or.(irc==(-21))) then
            write (iuer,*) 'recl=', it
         else if ((irc==(-7)).or.(irc==(-11)).or.(irc==(-15)).or.(irc==(&
     &-19)).or.(irc==(-23))) then
            write (iuer,*) 'fname=', fname
         else if (irc==(-24)) then
            write (iuer,*) 'confirmation itime,it=', itime, it
         endif
         end subroutine restart_dread
!
!              call dpostg(part,qe,np,nx,ny,qme,tdpost,inorder,dopt)
         subroutine dpostg(part,q,nop,nx,ny,qm,tdpost,inorder,dopt)
! deposit charge
         implicit none
         integer :: nop, nx, ny
         integer, optional :: inorder, dopt
         real :: qm, tdpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: q
! local data
         real, save :: zero = 0.0
! initialize charge density to zero
         call sguard(q,zero,nx,ny,inorder)
! deposit charge density
         call dpost(part,q,nop,qm,tdpost,inorder,dopt)
! add guard cells for charge density
         call aguard(q,nx,ny,inorder)
         end subroutine dpostg
!
         subroutine pushg(part,fxy,nop,qbm,dt,ci,ek,tpush,nx,ny,ipbc,rel&
     &ativity,inorder,popt)
! push particles with 2d electrostatic fields
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! push particles
         if (relativity==1) then
            call rpush(part,fxy,nop,qbm,dt,ci,ek,tpush,nx,ny,ipbc,inorde&
     &r,popt)
         else
            call push(part,fxy,nop,qbm,dt,ek,tpush,nx,ny,ipbc,inorder,po&
     &pt)
         endif
         end subroutine pushg
!
         subroutine bpushg(part,fxy,bxy,nop,qbm,dt,ci,ek,tpush,nx,ny,ipb&
     &c,relativity,inorder,popt)
! push particles with 2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! push particles
         if (relativity==1) then
            call rpush(part,fxy,bxy,nop,qbm,dt,dt,ci,ek,tpush,nx,ny,ipbc&
     &,inorder,popt)
         else
            call push(part,fxy,bxy,nop,qbm,dt,dt,ek,tpush,nx,ny,ipbc,ino&
     &rder,popt)
         endif
         end subroutine bpushg
!
         subroutine pushzfg(part,nop,dt,ci,ek,tpush,nx,ny,ipbc,relativit&
     &y)
! push particles with no forces
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         real :: dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
! push particles
         if (relativity==1) then
            call rpushzf(part,nop,dt,ci,ek,tpush,nx,ny,ipbc)
         else
            call pushzf(part,nop,dt,ek,tpush,nx,ny,ipbc)
         endif
         end subroutine pushzfg
!
         subroutine pushglg(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpush,r&
     &elativity)
! push particles with 2d electrostatic fields using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! push particles
         if (relativity==1) then
            call rpushgl(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpush)
         else
            call pushgl(part,fxy,nop,qbm,dt,ek,nx,ny,ipbc,tpush)
         endif
         end subroutine pushglg
!
         subroutine pushglxg(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpush,&
     &relativity,popt)
! push particles with 2d electrostatic fields using gridless method
! optimized method
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy
! push particles
         if (relativity==1) then
            call rpushglx(part,fxy,nop,qbm,dt,ci,ek,nx,ny,ipbc,tpush,pop&
     &t)
         else
            call pushglx(part,fxy,nop,qbm,dt,ek,nx,ny,ipbc,tpush,popt)
         endif
         end subroutine pushglxg
!
         subroutine initmodediag(dent,ntd,nxh,nyh,modesxd,modesyd,iud,nd&
     &rec,fdname)
! initialize mode diagnostic
         implicit none
         integer :: ntd, nxh, nyh, modesxd, modesyd, iud, ndrec
         character(len=*) :: fdname
         complex, dimension(:,:), pointer :: dent
! local data
         integer :: modesy2d
         if (ntd <= 0) return
         if (modesxd > nxh) modesxd = nxh
         if (modesyd > nyh) modesyd = nyh
         modesy2d = 2*modesyd - 1
         allocate(dent(modesxd,modesy2d))
! open output file
         if (ndrec==0) then
            ndrec = -1; iud = get_funit(iud)
            call bfopen(dent,modesxd,modesy2d,iud,ndrec,trim(fdname))
         endif
         end subroutine initmodediag
!
         subroutine initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,nmv,ndim,iuv,&
     &fvname)
! initialize velocity diagnostic
         implicit none
         integer :: ntv, ndv, nmv, ndim, iuv
         real :: vtx, vty, vtz
         real, dimension(:,:), pointer :: fv, fvm
         character(len=*) :: fvname
         if ((ntv <= 0) .and. (ndv <= 0)) return
         allocate(fv(2*nmv+2,ndim),fvm(3,ndim))
! fix velocity range
         if (ndim==2) then
            fv(1,:) = 8.*max(vtx,vty)
         elseif (ndim==3) then
            fv(1,:) = 8.*max(vtx,vty,vtz)
         endif
         if (ntv > 0) then
            iuv = get_funit(iuv)
            open(unit=iuv,file=trim(fvname),form='formatted',status='unk&
     &nown')
! write captions
            if (ndim==2) then
               write (iuv,*) 'it vdx vdy vtx vty sk'
            else if (ndim==3) then
               write (iuv,*) 'it vdx vdy vdz vtx vty vtz sk'
            endif
         endif
         end subroutine initveldiag
!
         subroutine dendiag(qt,qi,sfield,dent,ffc,mixup,sct,tfft,ntd,ndd&
     &,nx,ny,modesxd,modesyd,iud,ndrec,indx,indy,ntime,ndstyle,iflg,irc,&
     &inorder)
! ion density diagnostic
         implicit none
         integer :: ntd, ndd, nx, ny, modesxd, modesyd, iud, ndrec
         integer :: indx, indy, ntime, ndstyle, iflg, irc
         real :: tfft
         integer, optional :: inorder
         real, dimension(:,:), pointer :: qt, qi, sfield
         complex, dimension(:,:), pointer :: ffc, dent
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2d, isign
         irc = 0
         if ((ntd > 0) .or. (ndd > 0)) then
            it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
            jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
            if ((it==0) .or. (jt==0)) then
               qt = qi
! transform ion density to fourier space
               if (iflg==1) then
                  isign = -1
                  call fft(qt,isign,mixup,sct,tfft,indx,indy,inorder)
               endif
! calculate smoothing in fourier space
               call spois(qt,sfield,ffc,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2d = 2*modesyd - 1
                  call gtmodes(sfield,dent,nx,ny,modesxd,modesyd,inorder&
     &)
! write diagnostic output
                  call writebf(dent,modesxd,modesy2d,iud,ndrec,order=LIN&
     &EAR)
               endif
! transform ion density to real space
               if (jt==0) then
                  isign = 1
                  call fft(sfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
                  call cguard(sfield,nx,ny,inorder)
! display ion density
                  call displays(sfield,' ION DENSITY',ntime,999,2,ndstyl&
     &e,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine dendiag
!
         subroutine veldiag(part,fv,fvm,ntv,ndv,np,nmv,iuv,ntime,label,i&
     &rc)
! velocity diagnostic
         implicit none
         integer :: ntv, ndv, np, nmv, iuv, ntime, irc
         character(len=*) :: label
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: fv, fvm
! local data
         integer :: it, jt, nt
         irc = 0
         if ((ntv > 0) .or. (ndv > 0)) then
            it = -1; if (ntv > 0) it = ntime - ntv*(ntime/ntv)
            jt = -1; if (ndv > 0) jt = ntime - ndv*(ntime/ndv)
            if ((it==0) .or. (jt==0)) then
! calculate paticle distribution function and moments
               call vdist(part,fv,fvm,np,nmv)
! print out velocity moments
               if (it==0) then
                  nt = ntime/ntv
                  write (iuv,*) nt, fvm(1,:), fvm(2,:), sum(fvm(3,:))
               endif
! display velocity distributions
               if (jt==0) then
                  call displayfv(fv,fvm,label,ntime,nmv,2,irc)
               endif
            endif
         endif
         end subroutine veldiag
!
         subroutine phasediag(part,nts,nds,nx,ny,np,npxy,ntime,label,irc&
     &)
! phase space diagnostic
         implicit none
         integer :: nts, nds, nx, ny, np, npxy, ntime, irc
         character(len=*) :: label
         real, dimension(:,:), pointer :: part
! local data
         integer :: it, jt, isc
         irc = 0
         if ((nts > 0) .or. (nds > 0)) then
            it = -1; if (nts > 0) it = ntime - nts*(ntime/nts)
            jt = -1; if (nds > 0) jt = ntime - nds*(ntime/nds)
            if ((it==0) .or. (jt==0)) then
               isc = 999
! plot particles vx versus x
               call grasp(part,np,label,ntime,isc,nx,ny,3,1,npxy,irc)
               if (irc /= 0) return
! plot particles vy versus y
               call grasp(part,np,label,ntime,isc,nx,ny,4,2,npxy,irc)
            endif
         endif
         end subroutine phasediag
!
         subroutine potdiag(qe,sfield,pott,ffc,mixup,sct,tfft,ntp,ndp,nx&
     &,ny,modesxp,modesyp,iup,nprec,indx,indy,ntime,ndstyle,irc,inorder)
! potential diagnostic
         implicit none
         integer :: ntp, ndp, nx, ny, modesxp, modesyp, iup, nprec
         integer :: indx, indy, ntime, ndstyle, irc
         integer :: iii,jjj
         real :: tfft
         integer, optional :: inorder
         real, dimension(:,:), pointer :: qe, sfield
         complex, dimension(:,:), pointer :: ffc, pott
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2p, isign
         real :: we
         irc = 0
         if ((ntp > 0) .or. (ndp > 0)) then
            it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
            jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
            if ((it==0) .or. (jt==0)) then
! calculate potential in fourier space
               call pois(qe,sfield,ffc,we,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2p = 2*modesyp - 1
                  call gtmodes(sfield,pott,nx,ny,modesxp,modesyp,inorder&
     &)
! write diagnostic output
                  call writebf(pott,modesxp,modesy2p,iup,nprec,order=LIN&
     &EAR)
               endif
               
               !print *, modesxp,' ',modesyp,' ',modesy2p
               
                !do jjj = 1, modesy2p
                !    do iii = 1,modesxp
                !        write(12345,'(ES24.14E3)'),pott(iii,jjj)
                !    enddo
                !enddo
                
! transform potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(sfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
                  call cguard(sfield,nx,ny,inorder)
! display potential
                  call displays(sfield,' POTENTIAL',ntime,999,0,ndstyle,&
     &nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine potdiag
!
         subroutine emomtdiag(part,pxe,pye,pze,sx,sy,sz,wx,wy,wz,ntm,np,&
     &ium,ntime,ndim)
! calculate electron momentum
         implicit none
         integer :: ntm, np, ium, ntime, ndim
         real :: pxe, pye, pze, sx, sy, sz, wx, wy, wz
         real, dimension(:,:), pointer :: part
! local data
         integer :: it
         if (ntm > 0) then
            it = ntime/ntm
! calculate the momentum in the electrons
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it >= 0) then
               call premoment2(part,ntime,np,ium,pxe,pye,pze,sx,sy,sz,wx&
     &,wy,wz,ndim,nprint=it)
            endif
         endif
         end subroutine emomtdiag
!
         subroutine imomtdiag(parti,qi,fxyze,dt,rmass,pxi,pyi,pzi,wx,wy,&
     &wz,ntm,movion,npi,ium,nx,ny,ntime,ndim,inorder)
! calculate ion and field momentum
         implicit none
         integer :: ntm, movion, npi, ium, nx, ny, ntime, ndim
         real :: dt, rmass, pxi, pyi, pzi, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:), pointer :: parti
         real, dimension(:,:), pointer :: qi
         real, dimension(:,:,:), pointer :: fxyze
! local data
         integer :: it
  996    format (' total momentum = ',3e14.7)
         if (ntm > 0) then
            it = ntime/ntm
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it >= 0) then
! calculate ion momentum
               if (movion==0) then
                  if (it==1) then
                     call imoment(qi,fxyze,ium,pxi,pyi,pzi,dt,wx,wy,wz,n&
     &x,ny,inorder)
                  endif
               else if (movion==1) then
                  call primoment2(parti,npi,ium,rmass,pxi,pyi,pzi,wx,wy,&
     &wz,ndim,nprint=it)
               endif
! print total momentum
               if (it==1) then
                  write (ium,996) wx, wy, wz
               endif
            endif
         endif
         end subroutine imomtdiag
!
         subroutine esenergy(wt,we,wke,wki,ntw,ndw,itw,iuot,ntime)
! electrostatic energy diagnostic
         implicit none
         integer :: ntw, ndw, itw, iuot, ntime
         real :: we, wke, wki
         real, dimension(:,:), pointer :: wt
! local data
         integer :: it, jt
         real :: wtot
  992    format (' field, kinetic, total energies = ',3e14.7)
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
               wtot = we + wke + wki
               if (it==0) then
                  write (iuot,992) we, wke, wtot
               endif
               if (jt==0) then
                  itw = itw + 1
                  wt(itw,:) = (/we,wke,wki,wtot/)
               endif
            endif
         endif
         end subroutine esenergy
!
      end module simul2d