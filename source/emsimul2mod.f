!-----------------------------------------------------------------------
!
      module emsimul2d
! Higher level subroutines for electromagnetics
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: november 16, 2009
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use diag2d, only: get_funit, bfopen, writebf, displayv
      use simul2d
      use empush2d, only: retard, djpost, pushzf, push3, rdjpost, rpush3&
     &, rpush3zf, dmjpost, rdmjpost, dcjpost, rdcjpost, premoment2,     &
     &djpostgl, rdjpostgl, push3gl, rpush3gl, fft, dmjpostgl, rdmjpostgl&
     &, dcjpostgl, rdcjpostgl, ordjpost
      use field2d, only: cguard, acguard, cuperp, sbpois, apois, avpot, &
     &avrpot, gtmodes, poynt
      implicit none
      private
      public :: restart_open, restart_bwrite, restart_dwrite
      public :: restart_bread, restart_dread
      public :: retardg, dpostg, pushg, djpostg, push3g, push3zfg
      public :: dmjpostg, dcjpostg, djpostglg, push3glg
      public :: dmjpostglg, dcjpostglg
      public :: initmodediag, initveldiag
      public :: dendiag, veldiag, phasediag, potdiag
      public :: initvmodediag, vpotdiagprep, vpotrdiagprep
      public :: vpotdiag, vcurdiag, avpotdiag, avcurdiag, vpotrdiag
      public :: EMdiag, ECharDendiag, ECurrdiag,ICharDendiag
      public :: dampfactorinit, dampEM, PartDiag
      public :: launchlaser, launchboundary
      public :: emomtdiag, fmomtdiag, imomtdiag, dmenergy, emenergy
      public :: EMdiagF
      public :: ECurdiag, TCurdiag, odjpostg, initpw2, partant, dipole_ant
      public :: pulsechop, cuant, pulse_fix
!

      contains
!
         subroutine retardg(part,nop,dtc,ci,nx,ny,ipbc,relativity,ndim)
! retards particle positions half time-step to deposit current
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: ndim
         real :: dtc, ci
         real, dimension(:,:), pointer :: part
         if (dtc==0.0) return
         if (relativity==1) then
            call retard(part,nop,dtc,ci,nx,ny,ipbc,ndim)
         else
            call retard(part,nop,dtc,nx,ny,ipbc,ndim)
         endif
         end subroutine retardg
!
!              call djpostg(part,cu,np,qme,dth,ci,tdjpost,nx,ny,ipbc,relativity, &
!    &inorder,djopt)
         subroutine djpostg(part,cu,nop,qm,dt,ci,tdjpost,nx,ny,ipbc,rela&
     &tivity,inorder,djopt)
! deposit current
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, djopt
         real :: qm, dt, ci, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: cu
! deposit current
         if (relativity==1) then
            call rdjpost(part,cu,nop,qm,dt,ci,tdjpost,nx,ny,ipbc,inorder&
     &,djopt)
         else
            call djpost(part,cu,nop,qm,dt,tdjpost,nx,ny,ipbc,inorder,djo&
     &pt)
         endif
         end subroutine djpostg
     
         subroutine odjpostg(part,vpart,cu,vcu,nop,qm,dth,ci,tdjpost,&
        &nx,ny,ipbc,relativity,inorder,djopt,dt)
! deposit charge conserving current following OSIRIS scheme
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, djopt
         real :: qm, dth, dt, ci, tdjpost
         real, dimension(:,:), pointer :: part, vpart
         real, dimension(:,:,:), pointer :: cu, vcu
! deposit current
         if (relativity==1) then
            call ordjpost(part,vpart,cu,vcu,nop,qm,dth,ci,tdjpost,&
           &nx,ny,ipbc,inorder,djopt,dt)
         else
     !       call djpost(part,cu,nop,qm,dt,tdjpost,nx,ny,ipbc,inorder,djo&
     !&pt)
         endif
         end subroutine odjpostg
!
         subroutine push3g(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,tpush,nx,ny&
     &,ipbc,relativity,inorder,popt)
! push particles with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
         if (relativity==1) then
            call rpush3(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,tpush,nx,ny,ip&
     &bc,inorder,popt)
         else
            call push3(part,fxy,bxy,nop,qbm,dt,dtc,ek,tpush,nx,ny,ipbc,i&
     &norder,popt)
         endif
         end subroutine push3g
!
         subroutine push3zfg(part,nop,dt,ci,ek,tpush,nx,ny,ipbc,ndim,rel&
     &ativity)
! push particles with no forces
         implicit none
         integer :: nop, nx, ny, ipbc, ndim, relativity
         real :: dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
! push particles
         if (relativity==1) then
            call rpush3zf(part,nop,dt,ci,ek,tpush,nx,ny,ipbc,ndim)
         else
            call pushzf(part,nop,dt,ek,tpush,nx,ny,ipbc)
         endif
         end subroutine push3zfg
!
         subroutine dmjpostg(part,amu,nop,qm,ci,tdcjpost,relativity,inor&
     &der,djopt)
! deposit momentum flux with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, relativity
         integer, optional :: inorder, djopt
         real :: qm, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: amu
! deposit momentum flux
         if (relativity==1) then
            call rdmjpost(part,amu,nop,qm,ci,tdcjpost,inorder,djopt)
         else
            call dmjpost(part,amu,nop,qm,tdcjpost,inorder,djopt)
         endif
         end subroutine dmjpostg
!
         subroutine dcjpostg(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,ci,td&
     &cjpost,relativity,inorder,djopt)
! deposit momentum flux, acceleration density, and current density
! with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, relativity
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy, cu, dcu, amu
         if (relativity==1) then
            call rdcjpost(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,ci,tdcjp&
     &ost,inorder,djopt)
         else
            call dcjpost(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,tdcjpost,&
     &inorder,djopt)
         endif
         end subroutine dcjpostg
!
         subroutine djpostglg(part,cu,nop,qm,dt,ci,nx,ny,ipbc,tdjpost,re&
     &lativity)
! deposit current using gridless method
         integer :: nop, nx, ny, ipbc, relativity
         real :: qm, dt, ci, tdjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: cu
         if (relativity==1) then
            call rdjpostgl(part,cu,nop,qm,dt,ci,nx,ny,ipbc,tdjpost)
         else
            call djpostgl(part,cu,nop,qm,dt,nx,ny,ipbc,tdjpost)
         endif
         end subroutine djpostglg
!
         subroutine push3glg(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,nx,ny,ipb&
     &c,tpush,relativity)
! push particles with 2-1/2d electromagnetic fields
! using gridless method
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         real :: qbm, dt, dtc, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy
! push particles
         if (relativity==1) then
            call rpush3gl(part,fxy,bxy,nop,qbm,dt,dtc,ci,ek,nx,ny,ipbc,t&
     &push)
         else
            call push3gl(part,fxy,bxy,nop,qbm,dt,dtc,ek,nx,ny,ipbc,tpush&
     &)
         endif
         end subroutine push3glg
!
         subroutine dmjpostglg(part,amu,nop,qm,ci,nx,ny,tdcjpost,relativ&
     &ity)
! deposit momentum flux using gridless method
! with 2-1/2d electromagnetic fields
         implicit none
         integer :: nop, nx, ny, relativity
         real :: qm, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: amu
! deposit momentum flux
         if (relativity==1) then
            call rdmjpostgl(part,amu,nop,qm,ci,nx,ny,tdcjpost)
         else
            call dmjpostgl(part,amu,nop,qm,nx,ny,tdcjpost)
! debug
!           call igmjpost2gl(part,amu,nop,qm,nx,ny,tdcjpost)
         endif
         end subroutine dmjpostglg
!
         subroutine dcjpostglg(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,ci,&
     &nx,ny,tdcjpost,relativity)
! deposit momentum flux, acceleration density, and curent density
! with 2-1/2d electromagnetic fields, using gridless method
         implicit none
         integer :: nop, nx, ny, relativity
         real :: qm, qbm, dt, ci, tdcjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, bxy, cu, dcu, amu
         if (relativity==1) then
            call rdcjpostgl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,ci,nx,&
     &ny,tdcjpost)
         else
            call dcjpostgl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,nx,ny,t&
     &dcjpost)
! debug
!           call igdcjpost2gl(part,fxy,bxy,cu,dcu,amu,nop,qm,qbm,dt,nx,n&
!    &y,tdcjpost)
         endif
         end subroutine dcjpostglg
!
         subroutine initvmodediag(vcurt,ntj,nxh,nyh,ndim,modesxj,modesyj&
     &,iuj,njrec,fjname)
! initialize vector mode diagnostic
         implicit none
         integer :: ntj, nxh, nyh, ndim, modesxj, modesyj, iuj, njrec
         character(len=*) :: fjname
         complex, dimension(:,:,:), pointer :: vcurt
! local data
         integer :: modesy2j
         if (ntj <= 0) return
         if (modesxj > nxh) modesxj = nxh
         if (modesyj > nyh) modesyj = nyh
         modesy2j = 2*modesyj - 1
         allocate(vcurt(ndim,modesxj,modesy2j))
! open output file
         if (njrec==0) then
            njrec = -1; iuj = get_funit(iuj)
            call bfopen(vcurt,modesxj,modesy2j,iuj,njrec,trim(fjname))
         endif
         end subroutine initvmodediag
!
!              call vpotdiagprep(cu,vfield,nte,nde,ntime-1) in newbbeps2.f
         subroutine vpotdiagprep(cu,vfield,nta,nda,ntime)
! vector potential diagnostic
         implicit none
         integer :: nta, nda, ntime
         real, dimension(:,:,:), pointer :: cu, vfield
! local data
         integer :: it, jt
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
! save old current
            if ((it==0) .or. (jt==0)) then
               vfield = cu
            endif
         endif
         end subroutine vpotdiagprep
!
         subroutine vpotrdiagprep(cu,vfield,nta,nda,ntime,diff)
! save data for electromagnetic diagnostic
         implicit none
         integer :: nta, nda, ntime
         real, dimension(:,:,:), pointer :: cu, vfield
         logical, optional :: diff
! local data
         integer :: is, js, it, jt
         logical :: ldiff
         if ((nta > 0) .or. (nda > 0)) then
            ldiff = .false.
            if (present(diff)) ldiff = diff
            is = -1; if (nta > 0) is = ntime - nta*(ntime/nta)
            js = -1; if (nda > 0) js = ntime - nda*(ntime/nda)
            it = -1; if (nta > 1) it = ntime - nta*((ntime-1)/nta) - 1
            jt = -1; if (nda > 1) jt = ntime - nda*((ntime-1)/nda) - 1
! save current
            if ((is==0) .or. (js==0) .or. (it==0) .or. (jt==0)) then
               if (ldiff) then
                  vfield = cu - vfield
               else
                  vfield = cu
               endif
            endif
         endif
         end subroutine vpotrdiagprep
!
         subroutine vpotdiag(cu,vfield,vpott,ffc,mixup,sct,ci,tfft,nta,n&
     &da,nx,ny,modesxa,modesya,iua,narec,indx,indy,ntime,ndstyle,irc,ino&
     &rder)
! static (darwin) vector potential diagnostic
         implicit none
         integer :: nta, nda, nx, ny, modesxa, modesya, iua, narec
         integer :: indx, indy, ntime, ndstyle, irc
         real :: ci, tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu, vfield
         complex, dimension(:,:,:), pointer :: vpott
         complex, dimension(:,:), pointer :: ffc
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2a, isign
         real :: wm
         irc = 0
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
            if ((it==0) .or. (jt==0)) then
! calculate vector potential in fourier space
               call apois(cu,vfield,ffc,ci,wm,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2a = 2*modesya - 1
                  call gtmodes(vfield,vpott,nx,ny,modesxa,modesya,inorde&
     &r)
! write diagnostic output
                  call writebf(vpott,modesxa,modesy2a,iua,narec,order=LI&
     &NEAR)
               endif
! transform vector potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
                  call cguard(vfield,nx,ny,inorder)
! display absolute value of vector potential
                  call displayv(vfield,' VECTOR POTENTIAL',ntime,999,1,n&
     &dstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine vpotdiag

! Particle diagnostics
         subroutine PartDiag(part,ntime,ntpart,idimp, np)
         
         implicit none
         integer :: ntime,ntpart,idimp,np
         real, dimension(:,:), pointer :: part
         
         !local data
         integer :: iii,jjj,it
         character(len=256) :: filestrpart = ''
         character(len=256) :: filestr = ''
         integer, parameter :: idpart = 12401
        
         filestrpart = '-'//idx_string(ntime,8)//'-0000.dat'; filestrpart = trim(filestrpart)
         filestr = 'MS//PART//electron//elec'//filestrpart;
         filestr = trim(filestr)
         !print *,'idimp,np', idimp, np
         if ((ntpart > 0)) then
            it = -1; if (ntpart > 0) it = ntime - ntpart*(ntime/ntpart)
            if ((it==0)) then
                open (unit = idpart, file = filestr)
                write (idpart,'(ES24.14E2)') ((part(iii,jjj),iii=1,idimp),jjj=1,np)
                close(unit = idpart)
            endif
         endif
         
         end subroutine PartDiag

! EM diagnostics in real space 
         subroutine EMdiag(fxyz,bxyz,ntime,nteme1,nteme2,nteme3, ntemb1,ntemb2, &
     &ntemb3, ntle1,ntle2,ntle3, ntlb1,ntlb2, ntlb3,nxe,nye,inorder)

         implicit none
         
         integer :: ntime,nxe,nye,inorder
         integer :: nteme1,nteme2,nteme3, ntemb1,ntemb2, ntemb3
         integer :: ntle1,ntle2,ntle3, ntlb1,ntlb2, ntlb3
         real, dimension(:,:,:), pointer :: fxyz,bxyz
         character(len=256) :: filestrpart = ''
         character(len=256) :: filestr = ''
         integer, parameter :: ide1 = 12341, ide2 = 12342, ide3 = 12343
         integer, parameter :: idb1 = 12344, idb2 = 12345, idb3 = 12346
         integer, parameter :: ide1l = 12381, ide2l = 12382, ide3l = 12383
         integer, parameter :: idb1l = 12384, idb2l = 12385, idb3l = 12386
         
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec, ycenter
! EM diagnostics in real space
! slice output
         iend = nxe - 3; istart = 2;  
         jend = nye - 2; jstart = 2; ycenter = (jend-1) / 2
         if (inorder==LINEAR) then
            iend = nxe - 2; istart = 1; 
            jend = nye - 1; jstart = 1; ycenter = jend / 2
         else if (inorder==CUBIC) then 
            iend = nxe - 4; istart = 3
            jend = nye - 3; jstart = 3; ycenter = (jend-2) / 2
         endif
         ycenter = max(1,ycenter)
         ! slice out
         filestrpart = '-'//idx_string(ntime,8)//'-0000.dat'; filestrpart = trim(filestrpart)
         if ((nteme1 > 0)) then
            it = -1; if (nteme1 > 0) it = ntime - nteme1*(ntime/nteme1)
            if ((it==0)) then
                filestr = 'MS//FLD//e1//e1'//filestrpart; filestr = trim(filestr)
                open (unit = ide1, file = filestr)
                write (ide1,'(ES24.14E2)') ((fxyz(1,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = ide1)
            endif
         endif
         if ((nteme2 > 0)) then
            it = -1; if (nteme2 > 0) it = ntime - nteme2*(ntime/nteme2)
            if ((it==0)) then
                filestr = 'MS//FLD//e2//e2'//filestrpart; filestr = trim(filestr)
                open (unit = ide2, file = filestr)
                write (ide2,'(ES24.14E3)') ((fxyz(2,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = ide2)
            endif
         endif
         if ((nteme3 > 0)) then
            it = -1; if (nteme3 > 0) it = ntime - nteme3*(ntime/nteme3)
            if ((it==0)) then
                filestr = 'MS//FLD//e3//e3'//filestrpart; filestr = trim(filestr)
                open (unit = ide3, file = filestr)
                write (ide3,'(ES24.14E2)') ((fxyz(3,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = ide3)
            endif
         endif
         ! b field
         if ((ntemb1 > 0)) then
            it = -1; if (ntemb1 > 0) it = ntime - ntemb1*(ntime/ntemb1)
            if ((it==0)) then
                filestr = 'MS//FLD//b1//b1'//filestrpart; filestr = trim(filestr)
                open (unit = idb1, file = filestr)
                write (idb1,'(ES24.14E2)') ((bxyz(1,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idb1)
            endif
         endif
         if ((ntemb2 > 0)) then
            it = -1; if (ntemb2 > 0) it = ntime - ntemb2*(ntime/ntemb2)
            if ((it==0)) then
                filestr = 'MS//FLD//b2//b2'//filestrpart; filestr = trim(filestr)
                open (unit = idb2, file = filestr)
                write (idb2,'(ES24.14E2)') ((bxyz(2,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idb2)
            endif
         endif
         if ((ntemb3 > 0)) then
            it = -1; if (ntemb3 > 0) it = ntime - ntemb3*(ntime/ntemb3)
            if ((it==0)) then
                filestr = 'MS//FLD//b3//b3'//filestrpart; filestr = trim(filestr)
                open (unit = idb3, file = filestr)
                write (idb3,'(ES24.14E2)') ((bxyz(3,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idb3)
            endif
         endif
! line output
         if ((ntle1 > 0)) then
            it = -1; if (ntle1 > 0) it = ntime - ntle1*(ntime/ntle1)
            if ((it==0)) then
                filestr = 'MS//FLD//e1line//e1line'//filestrpart; filestr = trim(filestr)
                open (unit = ide1l, file = filestr)
                write (ide1l,'(ES24.14E2)') (fxyz(1,iii,ycenter),iii=istart,iend)
                close(unit = ide1l)
            endif
         endif
         if ((ntle2 > 0)) then
            it = -1; if (ntle2 > 0) it = ntime - ntle2*(ntime/ntle2)
            if ((it==0)) then
                filestr = 'MS//FLD//e2line//e2line'//filestrpart; filestr = trim(filestr)
                open (unit = ide2l, file = filestr)
                write (ide2l,'(ES24.14E2)') (fxyz(2,iii,ycenter),iii=istart,iend)
                close(unit = ide2l)
            endif
         endif
         if ((ntle3 > 0)) then
            it = -1; if (ntle3 > 0) it = ntime - ntle3*(ntime/ntle3)
            if ((it==0)) then
                filestr = 'MS//FLD//e3line//e3line'//filestrpart; filestr = trim(filestr)
                open (unit = ide3l, file = filestr)
                write (ide3l,'(ES24.14E2)') (fxyz(3,iii,ycenter),iii=istart,iend)
                close(unit = ide3l)
            endif
         endif
         if ((ntlb1 > 0)) then
            it = -1; if (ntlb1 > 0) it = ntime - ntlb1*(ntime/ntlb1)
            if ((it==0)) then
                filestr = 'MS//FLD//b1line//b1line'//filestrpart; filestr = trim(filestr)
                open (unit = idb1l, file = filestr)
                write (idb1l,'(ES24.14E2)') (bxyz(1,iii,ycenter),iii=istart,iend)
                close(unit = idb1l)
            endif
         endif
         if ((ntlb2 > 0)) then
            it = -1; if (ntlb2 > 0) it = ntime - ntlb2*(ntime/ntlb2)
            if ((it==0)) then
                filestr = 'MS//FLD//b2line//b2line'//filestrpart; filestr = trim(filestr)
                open (unit = idb2l, file = filestr)
                write (idb2l,'(ES24.14E2)') (bxyz(2,iii,ycenter),iii=istart,iend)
                close(unit = idb2l)
            endif
         endif
         if ((ntlb3 > 0)) then
            it = -1; if (ntlb3 > 0) it = ntime - ntlb3*(ntime/ntlb3)
            if ((it==0)) then
                filestr = 'MS//FLD//b3line//b3line'//filestrpart; filestr = trim(filestr)
                open (unit = idb3l, file = filestr)
                write (idb3l,'(ES24.14E2)') (bxyz(3,iii,ycenter),iii=istart,iend)
                close(unit = idb3l)
            endif
         endif
         end subroutine EMdiag
         
! EM diagnostics in real space 
         subroutine EMdiagF(fxyz,bxyz,ntime,nteme1,nteme2,nteme3, ntemb1,ntemb2, &
     &ntemb3, ntle1,ntle2,ntle3, ntlb1,ntlb2, ntlb3,nxe,nye,inorder)

         implicit none
         
         integer :: ntime,nxe,nye,inorder
         integer :: nteme1,nteme2,nteme3, ntemb1,ntemb2, ntemb3
         integer :: ntle1,ntle2,ntle3, ntlb1,ntlb2, ntlb3
         complex, dimension(:,:,:), pointer :: fxyz,bxyz
         character(len=256) :: filestrpart = ''
         character(len=256) :: filestr = ''
         integer, parameter :: ide1 = 12441, ide2 = 12442, ide3 = 12443
         integer, parameter :: idb1 = 12444, idb2 = 12445, idb3 = 12446
         integer, parameter :: ide1l = 12481, ide2l = 12482, ide3l = 12483
         integer, parameter :: idb1l = 12484, idb2l = 12485, idb3l = 12486
         
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec, ycenter
! EM diagnostics in real space
! slice output
         
         istart = 1; iend = nxe / 2; jstart = 1; jend = nye
         ycenter = max(1,ycenter)
         ! slice out
         filestrpart = '-'//idx_string(ntime,8)//'-0000.dat'; filestrpart = trim(filestrpart)
         if ((nteme1 > 0)) then
            it = -1; if (nteme1 > 0) it = ntime - nteme1*(ntime/nteme1)
            if ((it==0)) then
                filestr = 'MS//FLD//e1r//e1r'//filestrpart; filestr = trim(filestr)
                open (unit = ide1, file = filestr)
                write (ide1,'(ES24.14E2)') ((real(fxyz(1,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                write (ide1,'(ES24.14E2)') ((aimag(fxyz(1,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                close(unit = ide1)
            endif
         endif
         if ((nteme2 > 0)) then
            it = -1; if (nteme2 > 0) it = ntime - nteme2*(ntime/nteme2)
            if ((it==0)) then
                filestr = 'MS//FLD//e2r//e2r'//filestrpart; filestr = trim(filestr)
                open (unit = ide2, file = filestr)
                write (ide2,'(ES24.14E2)') ((real(fxyz(2,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                write (ide2,'(ES24.14E2)') ((aimag(fxyz(2,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                close(unit = ide2)
            endif
         endif
         if ((nteme3 > 0)) then
            it = -1; if (nteme3 > 0) it = ntime - nteme3*(ntime/nteme3)
            if ((it==0)) then
                filestr = 'MS//FLD//e3r//e3r'//filestrpart; filestr = trim(filestr)
                open (unit = ide3, file = filestr)
                write (ide3,'(ES24.14E2)') ((real(fxyz(3,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                write (ide3,'(ES24.14E2)') ((aimag(fxyz(3,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                close(unit = ide3)
            endif
         endif
         ! b field
         if ((ntemb1 > 0)) then
            it = -1; if (ntemb1 > 0) it = ntime - ntemb1*(ntime/ntemb1)
            if ((it==0)) then
                filestr = 'MS//FLD//b1r//b1r'//filestrpart; filestr = trim(filestr)
                open (unit = idb1, file = filestr)
                write (idb1,'(ES24.14E2)') ((real(bxyz(1,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                close(unit = idb1)
            endif
         endif
         if ((ntemb2 > 0)) then
            it = -1; if (ntemb2 > 0) it = ntime - ntemb2*(ntime/ntemb2)
            if ((it==0)) then
                filestr = 'MS//FLD//b2r//b2r'//filestrpart; filestr = trim(filestr)
                open (unit = idb2, file = filestr)
                write (idb2,'(ES24.14E2)') ((real(bxyz(2,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                close(unit = idb2)
            endif
         endif
         if ((ntemb3 > 0)) then
            it = -1; if (ntemb3 > 0) it = ntime - ntemb3*(ntime/ntemb3)
            if ((it==0)) then
                filestr = 'MS//FLD//b3r//b3r'//filestrpart; filestr = trim(filestr)
                open (unit = idb3, file = filestr)
                write (idb3,'(ES24.14E2)') ((real(bxyz(3,iii,jjj)),iii=istart,iend),jjj=jstart,jend)
                close(unit = idb3)
            endif
         endif
         
         end subroutine EMdiagF

! Electron current diagnostics in real space 
         subroutine ECurdiag(cu,ntime,ntej1,ntej2,ntej3, ntlej1,ntlej2,ntlej3,nxe,nye,inorder)

         implicit none
         
         integer :: ntime,nxe,nye,inorder
         integer :: ntej1,ntej2,ntej3
         integer :: ntlej1,ntlej2,ntlej3
         real, dimension(:,:,:), pointer :: cu
         character(len=256) :: filestrpart = ''
         character(len=256) :: filestr = ''
         integer, parameter :: idej1 = 12311, idej2 = 12312, idej3 = 12313
         integer, parameter :: idej1l = 12321, idej2l = 12322, idej3l = 12323
         
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec, ycenter
! Electron current diagnostics in real space
! slice output
         iend = nxe - 3; istart = 2;  
         jend = nye - 2; jstart = 2; ycenter = (jend-1) / 2
         if (inorder==LINEAR) then
            iend = nxe - 2; istart = 1; 
            jend = nye - 1; jstart = 1; ycenter = jend / 2
         else if (inorder==CUBIC) then 
            iend = nxe - 4; istart = 3
            jend = nye - 3; jstart = 3; ycenter = (jend-2) / 2
         endif
         ycenter = max(1,ycenter)
         ! slice out
         filestrpart = '-'//idx_string(ntime,8)//'-0000.dat'; filestrpart = trim(filestrpart)
         if ((ntej1 > 0)) then
            it = -1; if (ntej1 > 0) it = ntime - ntej1*(ntime/ntej1)
            if ((it==0)) then
                filestr = 'MS//CURRENT//electron//j1//j1'//filestrpart; filestr = trim(filestr)
                open (unit = idej1, file = filestr)
                write (idej1,'(ES24.14E2)') ((cu(1,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idej1)
            endif
         endif
         if ((ntej2 > 0)) then
            it = -1; if (ntej2 > 0) it = ntime - ntej2*(ntime/ntej2)
            if ((it==0)) then
                filestr = 'MS//CURRENT//electron//j2//j2'//filestrpart; filestr = trim(filestr)
                open (unit = idej2, file = filestr)
                write (idej2,'(ES24.14E3)') ((cu(2,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idej2)
            endif
         endif
         if ((ntej3 > 0)) then
            it = -1; if (ntej3 > 0) it = ntime - ntej3*(ntime/ntej3)
            if ((it==0)) then
                filestr = 'MS//CURRENT//electron//j3//j3'//filestrpart; filestr = trim(filestr)
                open (unit = idej3, file = filestr)
                write (idej3,'(ES24.14E2)') ((cu(3,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idej3)
            endif
         endif
! line output
         if ((ntlej1 > 0)) then
            it = -1; if (ntlej1 > 0) it = ntime - ntlej1*(ntime/ntlej1)
            if ((it==0)) then
                filestr = 'MS//CURRENT//electron//j1line//j1line'//filestrpart; filestr = trim(filestr)
                open (unit = idej1l, file = filestr)
                write (idej1l,'(ES24.14E2)') (cu(1,iii,ycenter),iii=istart,iend)
                close(unit = idej1l)
            endif
         endif
         if ((ntlej2 > 0)) then
            it = -1; if (ntlej2 > 0) it = ntime - ntlej2*(ntime/ntlej2)
            if ((it==0)) then
                filestr = 'MS//CURRENT//electron//j2line//j2line'//filestrpart; filestr = trim(filestr)
                open (unit = idej2l, file = filestr)
                write (idej2l,'(ES24.14E2)') (cu(2,iii,ycenter),iii=istart,iend)
                close(unit = idej2l)
            endif
         endif
         if ((ntlej3 > 0)) then
            it = -1; if (ntlej3 > 0) it = ntime - ntlej3*(ntime/ntlej3)
            if ((it==0)) then
                filestr = 'MS//CURRENT//electron//j3line//j3line'//filestrpart; filestr = trim(filestr)
                open (unit = idej3l, file = filestr)
                write (idej3l,'(ES24.14E2)') (cu(3,iii,ycenter),iii=istart,iend)
                close(unit = idej3l)
            endif
         endif
         end subroutine ECurdiag !!!

! Total current diagnostics in real space
!Changed output format from '(ES24.14E2)' --> '(ES24.14E3)' to accomodate 3 digit exponent in Jy
         subroutine TCurdiag(cu,ntime,nttj1,nttj2,nttj3, ntltj1,ntltj2,ntltj3,nxe,nye,inorder)

         implicit none
         
         integer :: ntime,nxe,nye,inorder
         integer :: nttj1,nttj2,nttj3
         integer :: ntltj1,ntltj2,ntltj3,i,j
         real, dimension(:,:,:), pointer :: cu
         character(len=256) :: filestrpart = ''
         character(len=256) :: filestr = ''
         integer, parameter :: idtj1 = 12314, idtj2 = 12315, idtj3 = 12316
         integer, parameter :: idtj1l = 12324, idtj2l = 12325, idtj3l = 12326
         
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec, ycenter
!         print*, 'Jy =', cu(2,2:5,2)
!         print*, 'Jz =', cu(3,2:5,2)
! Total current diagnostics in real space
! slice output
         iend = nxe - 3; istart = 2;  
         jend = nye - 2; jstart = 2; ycenter = (jend-1) / 2
         if (inorder==LINEAR) then
            iend = nxe - 2; istart = 1; 
            jend = nye - 1; jstart = 1; ycenter = jend / 2
         else if (inorder==CUBIC) then 
            iend = nxe - 4; istart = 3
            jend = nye - 3; jstart = 3; ycenter = (jend-2) / 2
         endif
!         do 10 i = istart, iend
!         do 20 j = jstart, jend
!            if(abs(cu(1,j,2)).le.10**(-100)) cu(1,j,2) = 0.0
!      20 continue
!      10 continue
         ycenter = max(1,ycenter)
         ! slice out
         filestrpart = '-'//idx_string(ntime,8)//'-0000.dat'; filestrpart = trim(filestrpart)
         if ((nttj1 > 0)) then
            it = -1; if (nttj1 > 0) it = ntime - nttj1*(ntime/nttj1)
            if ((it==0)) then
                filestr = 'MS//CURRENT//total//j1//j1'//filestrpart; filestr = trim(filestr)
                open (unit = idtj1, file = filestr)
                write (idtj1,'(ES24.14E2)') ((cu(1,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idtj1)
            endif
         endif
         if ((nttj2 > 0)) then
            it = -1; if (nttj2 > 0) it = ntime - nttj2*(ntime/nttj2)
            if ((it==0)) then
                filestr = 'MS//CURRENT//total//j2//j2'//filestrpart; filestr = trim(filestr)
                open (unit = idtj2, file = filestr)
                write (idtj2,'(ES24.14E3)') ((cu(2,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idtj2)
            endif
         endif
         if ((nttj3 > 0)) then
            it = -1; if (nttj3 > 0) it = ntime - nttj3*(ntime/nttj3)
            if ((it==0)) then
                filestr = 'MS//CURRENT//total//j3//j3'//filestrpart; filestr = trim(filestr)
                open (unit = idtj3, file = filestr)
                write (idtj3,'(ES24.14E2)') ((cu(3,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idtj3)
            endif
         endif
! line output
         if ((ntltj1 > 0)) then
            it = -1; if (ntltj1 > 0) it = ntime - ntltj1*(ntime/ntltj1)
            if ((it==0)) then
                filestr = 'MS//CURRENT//total//j1line//j1line'//filestrpart; filestr = trim(filestr)
                open (unit = idtj1l, file = filestr)
                write (idtj1l,'(ES24.14E2)') (cu(1,iii,ycenter),iii=istart,iend)
                close(unit = idtj1l)
            endif
         endif
         if ((ntltj2 > 0)) then
            it = -1; if (ntltj2 > 0) it = ntime - ntltj2*(ntime/ntltj2)
            if ((it==0)) then
                filestr = 'MS//CURRENT//total//j2line//j2line'//filestrpart; filestr = trim(filestr)
                open (unit = idtj2l, file = filestr)
                write (idtj2l,'(ES24.14E2)') (cu(2,iii,ycenter),iii=istart,iend)
                close(unit = idtj2l)
            endif
         endif
         if ((ntltj3 > 0)) then
            it = -1; if (ntltj3 > 0) it = ntime - ntltj3*(ntime/ntltj3)
            if ((it==0)) then
                filestr = 'MS//CURRENT//total//j3line//j3line'//filestrpart; filestr = trim(filestr)
                open (unit = idtj3l, file = filestr)
                write (idtj3l,'(ES24.14E2)') (cu(3,iii,ycenter),iii=istart,iend)
                close(unit = idtj3l)
            endif
         endif
         end subroutine TCurdiag !!!
         
! function 
         function idx_string( i,length )
!-------------------------------------------------------------------------------
! this function generates a string representing the integer i using the number
! of digits specified (i.e. adding leading zeroes) 
!-------------------------------------------------------------------------------
         implicit none
  ! arguments and return value
         integer :: i,length  
         character(len = length) :: idx_string 
  ! local variables
         character(len = 18) :: d
  
         d = trim(tostring(length))
         write( idx_string, '(i'//d//'.'//d//')' ) i
         end function idx_string

         function tostring( i )
!-------------------------------------------------------------------------------
! Converts the given integer to a left justified string
!-------------------------------------------------------------------------------
         implicit none
         integer :: i
         character(len=15) :: tostring

         write( tostring, * ) i  
  ! discard leading blanks
         tostring = adjustl(tostring)
         end function tostring
         
         
! EM diagnostics in real space 
         subroutine ECurrdiag(cu,ntime,nte,nxe,nye,inorder)
! vector potential diagnostic
         implicit none
         
         integer :: ntime, nte,nxe,nye,inorder
         real, dimension(:,:,:), pointer :: cu
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec
! vector potential diagnostic
         if ((nte > 0)) then
            it = -1; if (nte > 0) it = ntime - nte*(ntime/nte)
            if ((it==0)) then
! EM diagnostics in real space
                iend = nxe - 3; istart = 2
                jend = nye - 2; jstart = 2
                if (inorder==LINEAR) then
                    iend = nxe - 2; istart = 1
                    jend = nye - 1; jstart = 1
                else if (inorder==CUBIC) then
                    iend = nxe - 4; istart = 3
                    jend = nye - 3; jstart = 3
                endif
                !write (12341,'(ES24.14E2)') ((fxyz(2,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                write (12361,'(ES24.14E2)') ((cu(1,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                write (12362,'(ES24.14E2)') ((cu(2,iii,jjj),iii=istart,iend),jjj=jstart,jend)
                write (12363,'(ES24.14E2)') ((cu(3,iii,jjj),iii=istart,iend),jjj=jstart,jend)
            endif   
         endif
         end subroutine ECurrdiag
         
! Electron charge density diagnostics in real space 
         subroutine ECharDendiag(qe,ntime,ntem,nxe,nye,inorder)
! vector potential diagnostic
         implicit none
         
         integer :: ntime, ntem,nxe,nye,inorder
         real, dimension(:,:), pointer :: qe
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec
         character(len=256) :: filestrpart = ''
         integer, parameter :: idqe = 12351
! vector potential diagnostic
         if ((ntem > 0)) then
            it = -1; if (ntem > 0) it = ntime - ntem*(ntime/ntem)
            if ((it==0)) then
! EM diagnostics in real space
                iend = nxe - 3; istart = 2
                jend = nye - 2; jstart = 2
                if (inorder==LINEAR) then
                    iend = nxe - 2; istart = 1
                    jend = nye - 1; jstart = 1
                else if (inorder==CUBIC) then
                    iend = nxe - 4; istart = 3
                    jend = nye - 3; jstart = 3
                endif
                filestrpart = 'MS//DENSITY//electron//echarge-'//idx_string(ntime,8)//'-0000.dat'; 
                filestrpart = trim(filestrpart)
                open (unit = idqe, file = filestrpart)
                write (idqe,'(ES24.14E2)') ((qe(iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idqe)
            endif   
         endif
         end subroutine ECharDendiag
         
         subroutine ICharDendiag(qi,ntime,nte,nxe,nye,inorder)
! vector potential diagnostic
         implicit none
         
         integer :: ntime, nte,nxe,nye,inorder
         real, dimension(:,:), pointer :: qi
         
! local data
         integer :: it,iii,jjj,istart,iend,jstart,jend
         integer :: lrec
         character(len=256) :: filestrpart = ''
         integer, parameter :: idqi = 12352
! vector potential diagnostic
         if ((nte > 0)) then
            it = -1; if (nte > 0) it = ntime - nte*(ntime/nte)
            if ((it==0)) then
! EM diagnostics in real space
                iend = nxe - 3; istart = 2
                jend = nye - 2; jstart = 2
                if (inorder==LINEAR) then
                    iend = nxe - 2; istart = 1
                    jend = nye - 1; jstart = 1
                else if (inorder==CUBIC) then
                    iend = nxe - 4; istart = 3
                    jend = nye - 3; jstart = 3
                endif
                filestrpart = 'MS//DENSITY//ion//icharge-'//idx_string(ntime,8)//'-0000.dat'; 
                filestrpart = trim(filestrpart)
                open (unit = idqi, file = filestrpart)
                write (idqi,'(ES24.14E2)') ((qi(iii,jjj),iii=istart,iend),jjj=jstart,jend)
                close(unit = idqi)
            endif   
         endif
         end subroutine ICharDendiag    
         
! Launch laser in real space 
         subroutine launchlaser(fxyze,bxyze,nxe,nye,nx,ny,ci,launchlas, &
     &lase0, laspol, lasstart,lasend,launchtype,lask,ntime,inorder,dt,&
     &ntlas,launchshape,lasw0,lasfocus,rayl, gamma,antlength,launchdir)
! launch laser to fxyze and bxyze     
! launchtype: 1 = box launch; 2 = antenna launch
! launchlas: 1 = launch laser; 2 = finishing launching laser; 0 = no laser
! launchshape: 1 = plane wave; 2 = Gaussian wave
         implicit none
         integer :: nxe,nye,nx,ny
         real, dimension(:,:,:), pointer :: fxyze,bxyze
         real :: ci, lask, lasw0, rayl, dt, lase0, gamma
         integer :: launchlas, laspol, lasstart,lasend, launchtype
         integer :: launchshape, lasfocus,ntime,inorder, ntlas,antlength
         integer :: launchdir
! local data
         real :: xcenter,laserk,omegaenv,a0,tt,xx,xx2,yy,RL,wz,y0,curv,beta
         integer :: ii, jj, ifant
         real :: per_env, lon_env, pi = 3.1415926536
! set EM field to fxyze, bxyze
         !print *, 'launching the laser'
         if (launchtype == 1) then ! box launch
            fxyze = 0.0; bxyze = 0.0
            if (launchshape == 1) then ! plane wave
                print *, 'plane box'
                call LaunchPlaneBox(fxyze,bxyze,lasstart,lasend,lask,lase0,&
                    &inorder,ci,laspol,nxe,nye,launchlas,launchdir)
            else if (launchshape == 2) then !
                print *, 'Gaussian box'
                call LaunchGaussianBox(fxyze,bxyze,lasstart,lasend,lask,lase0,&
                 &inorder,ci,laspol,nxe,nye,nx,ny,launchlas,rayl,gamma,lasw0,lasfocus)
            else if (launchshape == 3) then
                print*, 'radial box'
                call LaunchRadialBox(fxyze,bxyze,lasstart,lasend,lask,lase0,&
                    &inorder,ci,laspol,nxe,nye,launchlas)
            endif ! launchshape if
         else if (launchtype == 2) then ! antenna launch
            ifant = ntime - ntlas*(ntime/ntlas) ! use antenna every ntlas steps
            if (ifant == 0) then
                if (launchshape == 1) then
                    print *,'plane antenna'
                    call LaunchPlaneAntenna(fxyze,bxyze,lasstart,lasend,lask,lase0,&
                        &inorder,ci,laspol,nxe,nye,ntime,dt,launchlas,antlength) 
                else if (launchshape == 2) then
                    print *, 'Gaussian antenna'
                    call LaunchGaussianAntenna(fxyze,bxyze,lasstart,lasend,lask,   &
                        &lase0,inorder,ci,laspol,nxe,nye,nx,ny,ntime,dt,launchlas,   &
                        &rayl,gamma,lasw0,lasfocus,antlength)
                endif ! launchshape if
            endif ! ifant if
         endif ! launchtype if
         
         end subroutine launchlaser
!*****
! initialize a plane wave polarized in z traveling in x. Bare bones debugger.
      subroutine initpw2(fxy,bxy,x0,x1,nx,ny,nxe,nye,waves,ci,isign)
      real fxy, bxy
      integer x0, x1, nx, ny, nxe, nye, env, isign
      dimension fxy(3,nxe,nye), bxy(3,nxe,nye)
      real ci, waves
      !local vars
      real kx, ky, wavenum, lambdax, lambday, xmid
      integer dx, dy, j
      dx = x1 - x0
      dy = x1 - x0
      lambdax = dx/waves
      lambday = dy/waves
      kx = 6.283185307179586/lambdax
      ky = 6.283185307179586/lambday
      fxy = 0.0
      bxy = 0.0
      xmid = (x0+x1)/2.0
      do j = x0, x1
        fxy(1:2,j,:) = 0.0
        fxy(3,j,:) = sin(kx*(j-1))!*exp(-0.1*(j-xmid)**2)
        bxy(1,j,:) = 0.0
        bxy(2,j,:) = -ci*fxy(3,j,:)
        bxy(3,j,:) = 0.0
      enddo
        fxy(3,1,:) = 0.0
        fxy(3,nx,:) = 0.0
        bxy(2,1,:) = 0.0
        bxy(2,nx,:) = 0.0
      print*, fxy(3,1,:)
      print*, fxy(3,nx,:)
      print*, bxy(2,1,:)
      print*, bxy(2,nx,:)
!      if(isign < 0) then
!        xmid = (x0+x1)/2.0
!        do j = 1, nx
!          fxy(3,j,:) = fxy(3,j,:)*exp(-0.1*(j-xmid)**2)
!          bxy(2,j,:) = bxy(2,j,:)*exp(-0.1*(j-xmid)**2)
!        enddo
!      endif
!      print*, x0, x1, xmid
!        do j = 1, ny
!          fxy(1:2,:,j) = 0.0
!          fxy(3,:,j) = sin(ky*(j-1))
!          bxy(1,:,j) = fxy(3,:,j)
!          bxy(2,:,j) = 0.0
!          bxy(3,:,j) = 0.0
!        enddo
      return
      end
      subroutine dipole_ant(part,parti,omega,a0,dt,ntime,ntime_end,np,npi, &
     & nx,ny,idimp)
      real part, omega, a0, dt
      real parti
      integer ntime, ntime_end, nx, ny, idimp, np
      integer npi
      dimension part(idimp,np)
      dimension parti(idimp,npi)
      ! local vars
      real x0, y0, v, time

      x0 = float(nx/2)
      y0 = float(ny/2)
      if(ntime.le.ntime_end) then
        if(ntime==0) then
          part(1,1) = x0
          part(2,1) = y0
          parti(1,1) = x0
          parti(2,1) = y0
        endif
        time = ntime*dt
        v = a0*sin(omega*time)
      else
        v = 0.0
        part(1,1) = x0
        part(2,1) = y0
        parti(1,1) = x0
        parti(2,1) = y0
      endif
      !dipole in x
      part(3,1) = v
      parti(3,1) = -v
      print*, part(1,1)-x0, part(2,1)-y0
      !print*, "---"
      print*, parti(1,1)-x0, parti(2,1)-y0
      !dipole in y
      !part(4,2) = v
      !parti(4,2) = -v   
      return
      end
      subroutine partant(part,oscs,frac_num,dt,ci,ntime,np,nx,ny,idimp)
      implicit none
      real :: part, dt, ci
      integer :: oscs, frac_num, ntime, nx, ny, idimp, np
      dimension part(idimp,np)
      ! local vars
      real :: pi, x0, y0, v, time, t_end, env, omega_env, omega
      integer :: ntime_end

      x0 = float(nx/2)
      y0 = float(ny/2)

      pi = 3.141592653589793
      omega = 1!2*pi/float(frac_num)/ci
      t_end = 2*pi*float(oscs)/omega
      omega_env = 1.0!pi/t_end
      ntime_end = nint(t_end/dt)
      !print*, ntime_end

      if(ntime.le.ntime_end) then
        if(ntime==0) then
          part(1,1) = x0
          part(2,1) = y0
        endif
        time = ntime*dt
        env = sin(omega_env*time)**2
        v = env*cos(omega*time)
      else
        v = 0.0
        part(1,1) = x0
        part(2,1) = y0
      endif
      part(5,1) = v
      print*, 'radiator velocity =', part(5,1) 
      return
      end

      subroutine cuant(cu,oscs,frac_num,dt,ci,ntime,nxe,nye,inorder)
      implicit none
      real :: cu, dt, ci
      integer :: oscs, frac_num, ntime, nxe, nye, inorder
      dimension cu(3,nxe,nye)
      ! local vars
      real :: pi, v, time, t_end, env, omega_env, omega
      integer :: ntime_end
      integer :: istart, iend, jstart, jend, xmid, ymid

      iend = nxe - 3; istart = 2; jend = nye - 2; jstart = 2
      if (inorder==LINEAR) then
         iend = nxe - 2; istart = 1; jend = nye - 1; jstart = 1
      else if (inorder==CUBIC) then 
         iend = nxe - 4; istart = 3; jend = nye - 3; jstart = 3;
      endif

      xmid = (iend+istart)/2
      ymid = (jend+jstart)/2

      pi = 3.141592653589793
      omega = 2*pi/float(frac_num)/ci
      t_end = 2*pi*float(oscs)/omega
      omega_env = pi/t_end
      ntime_end = nint(t_end/dt)
      !print*, ntime_end

      if(ntime.le.ntime_end) then
        time = ntime*dt
        env = sin(omega_env*time)**2
        v = env*cos(omega*time)
      else
        v = 0.0
      endif
      cu(3,xmid,ymid) = v
      print*, 'radiator velocity =', v 
      return
      end

      subroutine pulsechop(fxyze,bxyze,x0,y0,nxe,nye,inorder)
      implicit none
      integer :: x0, y0, nxe, nye, inorder
      real, dimension(:,:,:), pointer :: fxyze, bxyze
      !local variables
      integer :: istart, iend, jstart, jend
      integer :: xlow, xup, ylow, yup
 
      iend = nxe - 3; istart = 2; jend = nye - 2; jstart = 2
      if (inorder==LINEAR) then
         iend = nxe - 2; istart = 1; jend = nye - 1; jstart = 1
      else if (inorder==CUBIC) then 
         iend = nxe - 4; istart = 3; jend = nye - 3; jstart = 3;
      endif

      xlow = istart + x0; xup = iend - x0
      ylow = jstart + y0; yup = jend - y0

      fxyze(:,istart:xlow,:) = 0.0
      fxyze(:,xup:iend,:) = 0.0
      fxyze(:,:,jstart:ylow) = 0.0
      fxyze(:,:,yup:jend) = 0.0

      bxyze(:,istart:xlow,:) = 0.0
      bxyze(:,xup:iend,:) = 0.0
      bxyze(:,:,jstart:ylow) = 0.0
      bxyze(:,:,yup:jend) = 0.0

      return
      end

     
      
!*****        
         subroutine dampEM(fxyze,factor,nxe,nye,nx,ny,inorder,dampe1,dampe2,dampe3)
! launch laser to fxyze and bxyze     
! launchtype: 1 = box launch; 2 = antenna launch
! launchlas: 1 = launch laser; 2 = finishing launching laser; 0 = no laser
! launchshape: 1 = plane wave; 2 = Gaussian wave
         implicit none
         integer :: nxe,nye,nx,ny
         real, dimension(:,:,:), pointer :: fxyze
         real, dimension(:,:), pointer :: factor
         integer :: inorder
         integer :: dampe1, dampe2, dampe3
! local data
         integer :: ii, jj
         integer :: istart, iend, jstart, jend
         
         iend = nxe - 3; istart = 2; jend = nye - 2; jstart = 2
         if (inorder==LINEAR) then
            iend = nxe - 2; istart = 1; jend = nye - 1; jstart = 1
         else if (inorder==CUBIC) then 
            iend = nxe - 4; istart = 3; jend = nye - 3; jstart = 3;
         endif
         
         do ii = istart, iend
            do jj = jstart, jend
                if (dampe1 == 1) fxyze(1,ii,jj) = fxyze(1,ii,jj) * factor(ii,jj)
                if (dampe2 == 1) fxyze(2,ii,jj) = fxyze(2,ii,jj) * factor(ii,jj)
                if (dampe3 == 1) fxyze(3,ii,jj) = fxyze(3,ii,jj) * factor(ii,jj)
            enddo
         enddo
         end subroutine dampEM
         
         subroutine dampfactorinit(factor,nxe,nye,nx,ny,inorder,&
      &dampystart,dampyend,dampflag,dampdrop)
! initialize damping form factor
         implicit none
         integer :: nxe,nye,nx,ny
         real, dimension(:,:), pointer :: factor
         integer :: dampystart, dampyend,inorder,dampflag
         real :: dampdrop
! local data
         integer :: istart, iend, jstart, jend
         integer :: ii,jj
         real :: pi = 3.1415926536
         
         if (dampflag == 0) then
            print *, 'Initializing damping form factor'
            iend = nxe - 3; istart = 2; jend = nye - 2; jstart = 2
            if (inorder==LINEAR) then
                iend = nxe - 2; istart = 1; jend = nye - 1; jstart = 1
            else if (inorder==CUBIC) then 
                iend = nxe - 4; istart = 3; jend = nye - 3; jstart = 3;
            endif
         
            factor = 1.0
            if ((dampystart .ne. 0) .and. (dampyend .ne. 0)) then
                do ii = istart,iend
                    do jj = jstart, jstart + dampystart 
                        factor(ii,jj) = 1 - &
                        & dampdrop * cos(0.5*pi * real(jj-jstart)/real(dampystart))**2
                    enddo
                    do jj = jend - dampyend, jend 
                        factor(ii,jj) = 1 - &
                        & dampdrop * cos(0.5*pi * real(jend-jj)/real(dampyend))**2
                    enddo
                enddo
            endif
            !print *, istart,iend,jstart,jend
            !write (12372,'(ES24.14E2)') ((factor(ii,jj),ii=istart,iend),jj=jstart,jend)
         
            dampflag = 1
         endif
         
         end subroutine dampfactorinit
         
         subroutine LaunchPlaneBox(fxyze,bxyze,lasstart,lasend,lask,lase0,&
      &inorder,ci,laspol,nxe,nye,launchlas,launchdir)
         implicit none
         real, dimension(:,:,:), pointer :: fxyze, bxyze
         real :: ci, lask, lase0
         integer :: laspol, lasstart,lasend,inorder, nxe,nye, launchlas
         integer :: launchdir
! local data
         integer :: xstart, xend, ystart, yend, xoffset, yoffset
         integer :: xstart_b, xend_b, ystart_b, yend_b
         real :: xlaunch, ylaunch, dt = 0.08
         real :: xcenter, ycenter, laserk, omegaenv, a0, per_env, lon_env
         real :: vphase, vgroup
         real ::  xcenter_bl, xcenter_bp, ycenter_bl, ycenter_bp
         integer :: ii, jj
         real :: pi = 3.1415926536

         ! laser parameters, and phase and group velocities for B offset.
         laserk = lask; omegaenv = pi / real(lasend-lasstart); a0 = lase0
         vphase = 2/dt/lask*asin(lask*dt/ci/2.0)
         vgroup = 1/ci/cos(asin(lask*dt/ci/2.0))

         if(launchdir==1.or.launchdir==-1) then

         xlaunch = float(launchdir)      
         yend = nye - 2; ystart = 2; xoffset = 1
         if (inorder==LINEAR) then
            yend = nye - 1; ystart = 1; xoffset = 0
         else if (inorder==CUBIC) then
            yend = nye - 3; ystart = 3; xoffset = 2
         endif
         xstart = lasstart + xoffset; xend = lasend + xoffset
         xstart_b = xstart - launchdir*nint(vgroup*dt/2.0);
         xend_b = xend - launchdir*nint(vgroup*dt/2.0) 
         xcenter = 0.5 * real(lasstart+lasend)
         xcenter_bl = xcenter - xlaunch*vgroup*dt/2.0
         xcenter_bp = xcenter - xlaunch*vphase*dt/2.0


         ! initialize E. can combine loop over y.      
         do ii = xstart,xend
            do jj = ystart, yend
                per_env = cos(laserk*(real(ii)-xcenter))
                lon_env = cos(omegaenv*(real(ii)-xcenter)) ** 2
                if (laspol == 2) then ! e2 polarized
                    fxyze(2,ii,jj) = a0 * lon_env * per_env
                 else if (laspol == 3) then ! e3 polarized
                    fxyze(3,ii,jj) = a0 * lon_env * per_env
                 endif
            enddo
         enddo

         !Initialize B
         do ii = xstart_b,xend_b
            do jj = ystart, yend
                per_env = cos(laserk*(real(ii)-xcenter_bp))
                lon_env = cos(omegaenv*(real(ii)-xcenter_bl)) ** 2
                if (laspol == 2) then ! e2 polarized
                    bxyze(3,ii,jj) = xlaunch * ci * a0 * lon_env * per_env
                 else if (laspol == 3) then ! e3 polarized
                    bxyze(2,ii,jj) = - xlaunch * ci * a0 * lon_env * per_env
                 endif
            enddo
         enddo

         elseif(launchdir==2.or.launchdir==-2) then

         launchdir = launchdir/2
         ylaunch = float(launchdir)
         xend = nxe - 2; xstart = 2; yoffset = 2
         if (inorder==LINEAR) then
            xend = nxe - 1; xstart = 1; yoffset = 0
         else if (inorder==CUBIC) then
            xend = nxe - 3; xstart = 3; yoffset = 2
         endif
         ystart = lasstart + yoffset; yend = lasend + yoffset
         ystart_b = ystart - launchdir*nint(vgroup*dt/2.0);
         yend_b = yend - launchdir*nint(vgroup*dt/2.0) 
         ycenter = 0.5 * real(lasstart+lasend)
         ycenter_bl = ycenter - ylaunch*vgroup*dt/2.0
         ycenter_bp = ycenter - ylaunch*vphase*dt/2.0
         !laserk = lask; omegaenv = pi / real(yend-ystart); a0 = lase0
         if(laspol==2) then
           print*, 're-polarizing to be along x'
           laspol = 1
         endif

         ! initialize E. can combine loop over x.       
         do ii = ystart,yend
            do jj = xstart, xend
                per_env = cos(laserk*(real(ii)-ycenter))
                lon_env = cos(omegaenv*(real(ii)-ycenter)) ** 2
                if (laspol == 1) then ! e1 polarized
                    fxyze(1,jj,ii) = a0 * lon_env * per_env
                 else if (laspol == 3) then ! e3 polarized
                    fxyze(3,jj,ii) = a0 * lon_env * per_env
                 endif
            enddo
         enddo
         ! initialize B.
         do ii = ystart_b,yend_b
            do jj = xstart, xend
                per_env = cos(laserk*(real(ii)-ycenter_bp))
                lon_env = cos(omegaenv*(real(ii)-ycenter_bl)) ** 2
                if (laspol == 1) then ! e1 polarized
                    bxyze(3,jj,ii) = - ylaunch * ci * a0 * lon_env * per_env
                 else if (laspol == 3) then ! e3 polarized
                    bxyze(1,jj,ii) = ylaunch * ci * a0 * lon_env * per_env
                 endif
            enddo
         enddo

         endif

         launchlas = 2 ! stop the laser launching
         
         end subroutine LaunchPlaneBox
!*****
         subroutine LaunchRadialBox(fxyze,bxyze,lasstart,lasend,lask,lase0,&
      &inorder,ci,laspol,nxe,nye,launchlas)
         implicit none
         real, dimension(:,:,:), pointer :: fxyze, bxyze
         real :: ci, dt = 0.08, lask, lase0
         integer :: laspol, lasstart,lasend,inorder, nxe,nye, launchlas
         integer :: launchdir
! local data
         integer :: xcenter, ycenter, xstart, xend, ystart, yend, xoffset, yoffset
         real :: r0, laserk, omegaenv, a0, b0, bmag, per_env, lon_env, x, y, theta
         real :: r0_bp, r0_bl, per_env_b, lon_env_b
         real :: vgroup, vphase
         integer :: lasstart_b, lasend_b
         integer :: ii, jj, ix, jy, r2
         real :: sintheta, costheta, pi = 3.1415926536
     
         ycenter = (nye - 3)/2 + 1; xcenter = (nxe - 4)/2 + 1
         if (inorder==LINEAR) then
            ycenter = (nye - 1)/2; xcenter = (nxe - 2)/2
         else if (inorder==CUBIC) then
            ycenter = (nye - 5)/2 + 2; xcenter = (nxe - 6)/2 + 2
         endif

         xstart = xcenter - lasend; xend = xcenter + lasend
         ystart = ycenter - lasend; yend = ycenter + lasend
         laserk = lask; omegaenv = pi / real(lasend-lasstart); a0 = lase0

         vphase = 2/dt/lask*asin(laserk*dt/ci/2.0)
         vgroup = 1/ci/cos(asin(laserk*dt/ci/2.0))
         lasend_b = lasend - nint(vgroup*dt/2.0)
         lasstart_b = lasstart - nint(vgroup*dt/2.0)

         r0 = real(lasend+lasstart)/2.0
         r0_bp = r0 - vphase*dt/2.0
         r0_bl = r0 - vgroup*dt/2.0
        
         do ii = xstart,xend
            do jj = ystart, yend

              ix = (ii-xcenter); jy = (jj-ycenter)
              r2 = ix**2 + jy**2

              if((r2.le.lasend**2).and.(lasstart**2).le.r2) then
                per_env = cos(laserk*(sqrt(real(r2))-r0))
                lon_env = cos(omegaenv*(sqrt(real(r2))-r0)) ** 2
                fxyze(3,ii,jj) = a0 * lon_env * per_env
              else
                fxyze(:,ii,jj) = 0.0
              endif          

              if((r2.le.lasend_b**2).and.(lasstart_b**2).le.r2) then
                per_env_b = cos(laserk*(sqrt(real(r2))-r0_bp))
                lon_env_b = cos(omegaenv*(sqrt(real(r2))-r0_bl)) ** 2
                b0 = a0! * sqrt(1+vgroup*dt/(2.0*sqrt(real(r2))))
                bmag = ci * b0 * lon_env_b * per_env_b
                sintheta = real(jy)/sqrt(real(r2))
                costheta = real(ix)/sqrt(real(r2))
                bxyze(1,ii,jj) = bmag*sintheta
                bxyze(2,ii,jj) = -bmag*costheta
              else
                bxyze(:,ii,jj) = 0.0
              endif

            enddo
         enddo

         launchlas = 2 ! stop the laser launching
         
         end subroutine LaunchRadialBox
!*****
         
         subroutine LaunchPlaneAntenna(fxyze,bxyze,lasxstart,lasxend,lask,lase0,&
      &inorder,ci,laspol,nxe,nye,ntime,dt,launchlas,antlength)
         implicit none
         real, dimension(:,:,:), pointer :: fxyze,bxyze
         real :: ci, lask, lase0, dt
         integer :: laspol,lasxstart,lasxend,inorder,nxe,nye,ntime,launchlas
         integer :: antlength
! local data
         integer :: xstart,xend,ystart,yend, ant_length,xoffset
         real :: xcenter,laserk, omegaenv, a0
         integer :: ii, jj
         real :: per_env, lon_env, tt, pi = 3.1415926536
         
         yend = nye - 2; ystart = 2; xoffset = 1
         if (inorder==LINEAR) then
            yend = nye - 1; ystart = 1; xoffset = 0
         else if (inorder==CUBIC) then
            yend = nye - 3; ystart = 3; xoffset = 2
         endif
         
         tt = real(ntime) * dt; ant_length = antlength
         if (ntime == 0) then
            xend = lasxend + xoffset; xstart = xend - ant_length
         else
            xend = lasxend + xoffset - int(tt/ci) - ant_length/2
            xstart = xend - ant_length/2
         endif
         xcenter = real(lasxstart+lasxend) * 0.5 + (tt/ci)
         laserk = lask; omegaenv = pi / real(lasxend-lasxstart); a0 = lase0
            
         do ii = xstart,xend
            do jj = ystart, yend
                per_env = cos(laserk*(real(ii)-xcenter))
                lon_env = cos(omegaenv*(real(ii)-xcenter)) ** 2
                if (laspol == 2) then ! e2 polarized
                    fxyze(2,ii,jj) = a0 * lon_env * per_env 
                    bxyze(3,ii,jj) = ci * a0 * lon_env * per_env 
                else if (laspol == 3) then ! e3 polarized
                    fxyze(3,ii,jj) = a0 * lon_env * per_env 
                    bxyze(2,ii,jj) = - ci * a0 * lon_env * per_env
                endif
            enddo
         enddo
         if ((xcenter-real(xstart))*2 .ge. real(lasxend-lasxstart)) then
            launchlas = 2 ! stop using antenna
         endif
         end subroutine LaunchPlaneAntenna
         
         subroutine LaunchGaussianBox(fxyze,bxyze,lasxstart,lasxend,lask,lase0,&
           &inorder,ci,laspol,nxe,nye,nx,ny,launchlas,rayl,gamma,lasw0,lasfocus)
         implicit none
         real, dimension(:,:,:), pointer :: fxyze,bxyze
         real :: ci, lask, lase0, rayl, gamma, lasw0
         integer :: laspol, lasxstart,lasxend, lasfocus
         integer :: inorder, nxe,nye, launchlas,nx,ny
! local data
         integer :: xstart,xend,ystart,yend, xoffset
         real :: xcenter,laserk,omegaenv,a0,wz,curv,xx,xx2,yy,y0,RL,beta
         integer :: ii, jj
         real :: per_env, lon_env, pi = 3.1415926536
         
         yend = nye - 2; ystart = 2; y0 = ny/2 + 1; xoffset = 1
         if (inorder==LINEAR) then
            yend = nye - 1; ystart = 1; y0 = ny/2; xoffset = 0
         else if (inorder==CUBIC) then
            yend = nye - 3; ystart = 3; y0 = ny/2 + 2; xoffset = 2
         endif
         
         beta = sqrt(1.0 - 1.0/gamma/gamma)
         xstart = lasxstart + xoffset; xend = lasxend + xoffset; 
         xcenter = real(lasxstart+lasxend) * 0.5
         laserk = lask; omegaenv = pi / real(xend-xstart); a0 = lase0; RL = rayl
         
         do ii = xstart,xend
            do jj = ystart, yend
                xx = real(ii-lasfocus); xx2 = real(ii)-xcenter
                yy = real(y0-jj) + 0.5
                wz = lasw0*sqrt(1+xx*xx/RL/RL); 
                if (ii==lasfocus) then
                    curv = 0.0
                else
                    curv = 0.5 * yy*yy/ (xx+RL*RL/xx) * (1+beta) 
                endif
                lon_env = cos(omegaenv*xx2) ** 2
                per_env = sqrt(lasw0/wz)*exp(-yy*yy/wz/wz) * &
                            & cos(laserk*xx2+laserk*curv-atan(xx/RL))
                if (laspol == 2) then ! e2 polarized
                    fxyze(2,ii,jj) = a0 * per_env * lon_env
                    bxyze(3,ii,jj) = ci * a0 * per_env * lon_env
                else if (laspol == 3) then ! e3 polarized
                    !print*, 'here'
                    fxyze(3,ii,jj) = a0 * per_env*lon_env
                    bxyze(2,ii,jj) = -ci * a0 * per_env * lon_env
                    print*, fxyze(3,ii,jj)
                endif
            enddo
         enddo
         print*, maxval(abs(fxyze(3,:,:)))
         launchlas = 2
         end subroutine LaunchGaussianBox        
         
         subroutine LaunchGaussianAntenna(fxyze,bxyze,lasxstart,lasxend,lask,lase0, &
      &inorder,ci,laspol,nxe,nye,nx,ny,ntime,dt,launchlas,rayl,gamma,lasw0,lasfocus,&
      &antlength)
         implicit none
         real, dimension(:,:,:), pointer :: fxyze,bxyze
         real :: ci, lask, lase0, dt,rayl,gamma,lasw0
         integer :: laspol, lasxstart,lasxend, lasfocus,antlength
         integer :: inorder, nxe,nye,ntime, launchlas, nx, ny
! local data
         integer :: xstart,xend,ystart,yend, ant_length,xoffset
         real :: laserk,omegaenv,a0,beta,xx,xx2,yy,y0,wz,curv,RL,ttovci
         integer :: ii, jj
         real :: per_env, lon_env, tt, xfocus,xcenter
         real :: pi = 3.1415926536
         
         yend = nye - 2; ystart = 2; y0 = ny/2 + 1; xoffset = 1
         if (inorder==LINEAR) then
            yend = nye - 1; ystart = 1; y0 = ny/2; xoffset = 0
         else if (inorder==CUBIC) then
            yend = nye - 3; ystart = 3; y0 = ny/2 + 2; xoffset = 2
         endif
         
         tt = real(ntime) * dt; ant_length = antlength
         beta = sqrt(1.0 - 1.0/gamma/gamma)
         ttovci = tt/ci; 
         ! for ntime == 0, length of antenna is ant_length
         ! for ntime > 0, length of antenna is ant_length/2
         if (ntime == 0) then
            if (ttovci - int(ttovci) > 0.999) then
                xend = lasxend+xoffset - ( int(ttovci)+1 )
            else 
                xend = lasxend+xoffset - int(ttovci)
            endif
            xstart = xend - ant_length
         else
            if (ttovci - int(ttovci) > 0.999) then
                xend = lasxend+xoffset - ( int(ttovci)+1 ) - ant_length/2
            else 
                xend = lasxend+xoffset - int(ttovci) - ant_length/2
            endif
            xstart = xend - ant_length/2
         endif
         xcenter = real(lasxstart+lasxend) * 0.5 + ttovci
         laserk = lask; omegaenv = pi / real(lasxend-lasxstart); a0 = lase0; RL = rayl
         xfocus = real(lasfocus+xoffset) - ttovci * beta
         
         do ii = xstart,xend
            do jj = ystart, yend
                xx = real(ii)-xfocus; xx2 = real(ii)-xcenter
                yy = real(y0-jj) + 0.5
                wz = lasw0 * sqrt(1+xx*xx/RL/RL); 
                if ( real(ii)-xfocus .le. 0.00001) then
                    curv = 0.0
                else
                    curv = 0.5 * yy*yy/ (xx+RL*RL/xx) * (1+beta)
                endif
                lon_env = cos(omegaenv*xx2) ** 2
                per_env = sqrt(lasw0/wz)*exp(-yy*yy/wz/wz) * &
                            & cos(laserk*xx2+laserk*curv-atan(xx/RL))
                if (laspol == 2) then ! e2 polarized
                    fxyze(2,ii,jj) = a0 * per_env * lon_env
                    bxyze(3,ii,jj) = ci * a0 * per_env * lon_env
                else if (laspol == 3) then ! e3 polarized
                    fxyze(3,ii,jj) = a0 * per_env * lon_env
                    bxyze(2,ii,jj) = -ci * a0 * per_env * lon_env
                endif
            enddo
         enddo
         if ((xcenter-real(xstart))*2 .ge. real(lasxend-lasxstart)) then
            launchlas = 2 ! stop using antenna
         endif
         end subroutine LaunchGaussianAntenna
         
         subroutine launchboundary(fxyze,bxyze,nxe,nye,nx,ny,ci,launchlas, &
     &lase0, laspol, lasxstart,lasxend,launchtype,lask,ntime)
! vector potential diagnostic
         implicit none
         
         integer :: nxe,nye,nx,ny
         real, dimension(:,:,:), pointer :: fxyze,bxyze
         real :: ci, lask
         integer :: launchlas, laspol, lasxstart,lasxend, launchtype
         integer :: ntime
         real :: lase0
! local data
         integer :: xstart,xend,xcenter,ystart,yend
         integer :: xstart2, xend2
         integer :: ii, jj
         real :: laserk, omegaenv, a0
! set EM field to fxyze, bxyze
         !print *, 'zeroing the boundary'
           xstart = 1; xend = 20; 
           xstart2 = nxe - 19; xend2 = nxe;
           ystart = 1; yend = nye
           do ii = xstart,xend
            do jj = ystart, yend
              fxyze(1,ii,jj) = 0.0; bxyze(1,ii,jj) = 0.0
              fxyze(2,ii,jj) = 0.0; bxyze(2,ii,jj) = 0.0
              fxyze(3,ii,jj) = 0.0; bxyze(3,ii,jj) = 0.0
            enddo
           enddo
           do ii = xstart2,xend2
            do jj = ystart, yend
              fxyze(1,ii,jj) = 0.0; bxyze(1,ii,jj) = 0.0
              fxyze(2,ii,jj) = 0.0; bxyze(2,ii,jj) = 0.0
              fxyze(3,ii,jj) = 0.0; bxyze(3,ii,jj) = 0.0
            enddo
           enddo
         
         end subroutine launchboundary
 !*****
         subroutine pulse_fix(f,lasstart,lasend,nxe,nye,ntime,inorder,ci,dt)
         implicit none
         real, dimension(:,:,:), pointer :: f
         integer :: lasstart, lasend, nxe, nye, ntime, ntime_end, inorder
         real :: ci, dt
! local data
         integer :: xcenter, ycenter, xstart, xend, ystart, yend, xoffset, yoffset
         integer :: ii, jj, ix, jy, r2
         real :: time
         integer :: ntime_zero
     
         ycenter = (nye - 3)/2 + 1; xcenter = (nxe - 4)/2 + 1
         if (inorder==LINEAR) then
            ycenter = (nye - 1)/2; xcenter = (nxe - 2)/2
         else if (inorder==CUBIC) then
            ycenter = (nye - 5)/2 + 2; xcenter = (nxe - 6)/2 + 2
         endif
         
         ntime_zero = nint(ci*(lasend-lasstart)/dt) + 5
         xstart = xcenter - lasstart; xend = xcenter + lasstart
         ystart = ycenter - lasstart; yend = ycenter + lasstart

         if(ntime==ntime_zero) then        
           do ii = xstart,xend
             do jj = ystart, yend           
               ix = (ii-xcenter); jy = (jj-ycenter)
               r2 = ix**2 + jy**2
               if(r2.le.lasstart**2) f(:,ii,jj) = 10**(-5.5)
             enddo
           enddo
         endif
         
         end subroutine pulse_fix        
!
         subroutine vcurdiag(cut,cu,vfield,vcurt,ffc,mixup,sct,tfft,ntj,&
     &ndj,nx,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,ndstyle,iflg,i&
     &rc,inorder)
! static (darwin) ion current diagnostic
         implicit none
         integer :: ntj, ndj, nx, ny, modesxj, modesyj, iuj, njrec
         integer :: indx, indy, ntime, ndstyle, iflg, irc
         real :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cut, cu, vfield
         complex, dimension(:,:,:), pointer :: vcurt
         complex, dimension(:,:), pointer :: ffc
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2j, isign
         irc = 0
         if ((ntj > 0) .or. (ndj > 0)) then
            it = -1; if (ntj > 0) it = ntime - ntj*(ntime/ntj)
            jt = -1; if (ndj > 0) jt = ntime - ndj*(ntime/ndj)
            if ((it==0) .or. (jt==0)) then
               cut = cu - vfield
! add guard cells for current
               call acguard(cut,nx,ny,inorder)
! transform ion current to fourier space
               if (iflg==1) then
                  isign = -1
                  call fft(cut,isign,mixup,sct,tfft,indx,indy,inorder)
               endif
! take transverse part of current
               call cuperp(cut,nx,ny,inorder)
! calculate smoothing in fourier space
               call sbpois(cut,vfield,ffc,nx,ny,inorder)
! calculate vector potential in fourier space
!              call apois(cut,vfield,ffc,ci,wm,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2j = 2*modesyj - 1
                  call gtmodes(vfield,vcurt,nx,ny,modesxj,modesyj,inorde&
     &r)
! write diagnostic output
                  call writebf(vcurt,modesxj,modesy2j,iuj,njrec,order=LI&
     &NEAR)
               endif
! transform ion current to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
                  call cguard(vfield,nx,ny,inorder)
! display absolute value of ion current
                  call displayv(vfield,' ION CURRENT',ntime,999,1,ndstyl&
     &e,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine vcurdiag
!
         subroutine avpotdiag(bxyz,vfield,vpott,mixup,sct,tfft,nta,nda,n&
     &x,ny,modesxa,modesya,iua,narec,indx,indy,ntime,ndstyle,irc,inorder&
     &)
! vector potential diagnostic
         implicit none
         integer :: nta, nda, nx, ny, modesxa, modesya, iua, narec
         integer :: indx, indy, ntime, ndstyle, irc
         real :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: vfield
         complex, dimension(:,:,:), pointer :: bxyz, vpott
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2a, isign
         irc = 0
! vector potential diagnostic
         if ((nta > 0) .or. (nda > 0)) then
            it = -1; if (nta > 0) it = ntime - nta*(ntime/nta)
            jt = -1; if (nda > 0) jt = ntime - nda*(ntime/nda)
            if ((it==0) .or. (jt==0)) then
! calculate vector potential in fourier space
               call avpot(bxyz,vfield,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2a = 2*modesya - 1
                  call gtmodes(vfield,vpott,nx,ny,modesxa,modesya,inorde&
     &r)
! write diagnostic output
                  call writebf(vpott,modesxa,modesy2a,iua,narec,order=LI&
     &NEAR)
               endif
! transform vector potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
                  call cguard(vfield,nx,ny,inorder)
! display absolute value of vector potential
                  call displayv(vfield,' VECTOR POTENTIAL',ntime,999,1,n&
     &dstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine avpotdiag
!
         subroutine avcurdiag(cut,cu,vfield,vcurt,ffc,mixup,sct,tfft,ntj&
     &,ndj,nx,ny,modesxj,modesyj,iuj,njrec,indx,indy,ntime,ndstyle,iflg,&
     &irc,inorder)
! ion current diagnostic for electromagnetic code
         implicit none
         integer :: ntj, ndj, nx, ny, modesxj, modesyj, iuj, njrec
         integer :: indx, indy, ntime, ndstyle, iflg, irc
         real :: tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cut, cu, vfield
         complex, dimension(:,:,:), pointer :: vcurt
         complex, dimension(:,:), pointer :: ffc
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: is, js, it, jt, modesy2j, isign
         irc = 0
         if ((ntj > 0) .or. (ndj > 0)) then
            is = -1; if (ntj > 0) is = ntime - ntj*((ntime)/ntj)
            js = -1; if (ndj > 0) js = ntime - ndj*((ntime)/ndj)
            it = -1; if (ntj > 0) it = ntime - ntj*((ntime-1)/ntj) - 1
            jt = -1; if (ndj > 0) jt = ntime - ndj*((ntime-1)/ndj) - 1
! save current if needed next time and not saved below
            if (((is==0) .or. (js==0)).and.((it/=0) .and. (jt/=0))) then
               cu = cut
            endif
            if ((it==0) .or. (jt==0)) then
! calculate averaged ion current
               vfield = 0.5*(cut + cu)
! save current if needed next time and not saved above
               if ((is==0) .or. (js==0)) then
                  cu = cut
               endif
! add guard cells for current
               call acguard(vfield,nx,ny,inorder)
! transform ion current to fourier space
               if (iflg==1) then
                  isign = -1
                  call fft(vfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
               endif
! take transverse part of current
               call cuperp(vfield,nx,ny,inorder)
! calculate smoothing in fourier space
               call sbpois(vfield,cut,ffc,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2j = 2*modesyj - 1
                  call gtmodes(cut,vcurt,nx,ny,modesxj,modesyj,inorder)
! write diagnostic output
                  call writebf(vcurt,modesxj,modesy2j,iuj,njrec,order=LI&
     &NEAR)
               endif
! transform ion current to real space
               if (jt==0) then
                  isign = 1
                  call fft(cut,isign,mixup,sct,tfft,indx,indy,inorder)
                  call cguard(cut,nx,ny,inorder)
! display absolute value of ion current
                  call displayv(cut,' ION CURRENT',ntime,999,1,ndstyle,n&
     &x,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine avcurdiag
!
         subroutine vpotrdiag(bxyz,cu,vfield,vpotr,ffc,mixup,sct,tfft,ci&
     &,nte,nde,nx,ny,modesxe,modesye,iue,nerec,indx,indy,ntime,ndstyle,i&
     &rc,inorder)
! electromagnetic diagnostic
         implicit none
         integer :: nte, nde, nx, ny, modesxe, modesye, iue, nerec, nlsrerec
         integer :: indx, indy, ntime, ndstyle, irc
         real :: ci, tfft
         integer, optional :: inorder
         real, dimension(:,:,:), pointer :: cu, vfield
         complex, dimension(:,:,:), pointer :: bxyz, vpotr
         complex, dimension(:,:), pointer :: ffc
         integer, dimension(:), pointer :: mixup
         integer :: iii,jjj,kkk
         complex, dimension(:), pointer :: sct
! local data
         integer :: it, jt, modesy2e, isign
         irc = 0
         if ((nte > 0) .or. (nde > 0)) then
            it = -1; if (nte > 0) it = ntime - nte*((ntime-1)/nte) - 1
            jt = -1; if (nde > 0) jt = ntime - nde*((ntime-1)/nde) - 1
            if ((it==0) .or. (jt==0)) then
! calculate averaged radiative vector potential
               vfield = 0.5*(vfield + cu)
               call avrpot(vfield,bxyz,ffc,ci,nx,ny,inorder)
! store selected fourier modes
               if (it==0) then
                  modesy2e = 2*modesye - 1
                  call gtmodes(vfield,vpotr,nx,ny,modesxe,modesye,inorde&
     &r)
! write diagnostic output
                  call writebf(vpotr,modesxe,modesy2e,iue,nerec,order=LI&
     &NEAR)
               endif
               
! transform radiative vector potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(vfield,isign,mixup,sct,tfft,indx,indy,inorder&
     &)
                  call cguard(vfield,nx,ny,inorder)
! display absolute value of radiative vector potential
                  call displayv(vfield,' RADIATIVE VPOTENTIAL',ntime,999&
     &,1,ndstyle,nx,ny,irc,inorder)
               endif
            endif
         endif
         end subroutine vpotrdiag
!
         subroutine fmomtdiag(part,qe,ffc,exyz,bxyz,pxe,pye,pze,sx,sy,sz&
     &,wx,wy,wz,ntm,np,ium,nx,ny,ntime,inorder)
! calculate electron and field momentum
         implicit none
         integer :: ntm, np, ium, nx, ny, ntime
         real :: pxe, pye, pze, sx, sy, sz, wx, wy, wz
         integer, optional :: inorder
         real, dimension(:,:), pointer :: part
         real, dimension(:,:), pointer :: qe
         complex, dimension(:,:,:), pointer :: exyz, bxyz
         complex, dimension(:,:), pointer :: ffc
! local data
         integer :: it
         if (ntm > 0) then
            it = ntime/ntm
! calculate the momentum in the electromagnetic field
            if (ntime==ntm*it) then
               call poynt(qe,exyz,bxyz,ffc,sx,sy,sz,nx,ny,inorder)
            endif
! calculate the momentum in the electrons
            it = ntime - ntm*it + 1
            if (it > 1) it = it - ntm
            if (it >= 0) then
               call premoment2(part,ntime,np,ium,pxe,pye,pze,sx,sy,sz,wx&
     &,wy,wz,nprint=it)
            endif
         endif
         end subroutine fmomtdiag
!
         subroutine dmenergy(wt,we,wf,wm,wke,wki,ntw,ndw,itw,iuot,ntime)
! darwin electromagnetic energy diagnostic
         implicit none
         integer :: ntw, ndw, itw, iuot, ntime
         real :: we, wf, wm, wke, wki
         real, dimension(:,:), pointer :: wt
! local data
         integer :: it, jt
         real :: wef, wtot
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
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
         end subroutine dmenergy
!
         subroutine emenergy(wt,we,wf,wm,wke,wki,ntw,ndw,itw,iuot,ntime)
! electromagnetic energy diagnostic
         implicit none
         integer :: ntw, ndw, itw, iuot, ntime
         real :: we, wf, wm, wke, wki
         real, dimension(:,:), pointer :: wt
! local data
         integer :: it, jt
         real :: wef, wtot
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
               wef = we + wf + wm
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
         end subroutine emenergy
!
      end module emsimul2d
