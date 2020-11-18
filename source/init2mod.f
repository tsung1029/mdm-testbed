!-----------------------------------------------------------------------
!
      module init2d
!
! Fortran90 interface to 2d PIC Fortran77 library init2lib.f
! init2mod.f contains interface procedures to initialize particle
!            co-ordinates:
!            defines module init2d
! distr => idistr2 initializes x, y and vx, vy co-ordinates for 2d code,
!          with uniform density and maxwellian velocity with drift
!          calls DISTR2
! distr => idistrh2 initializes x, y and vx, vy, vz co-ordinates for
!          magnetized 2-1/2d codes, with uniform density and maxwellian
!          velocity with drift.
!          calls DISTR2H
! distr => ibdistr2 calculates guiding centers for magnetized
!          2-1/2d codes.
!          calls GBDISTR2L, or GBZDISTR2L
! distr => irbdistr2 calculates guiding centers for relativistic,
!          magnetized 2-1/2d codes.
!          calls GRBDISTR2L, or GRBZDISTR2L
! fdistr => ifdistr2 initializes x, y co-ordinates for 2d code, with
!           various density profiles.
!           calls FDISTR2, and RDISTR2
! vdistr => ivdistr2 initializes vx, vy co-ordinates for 2d code, with
!           maxwellian velocity distribution with drift.
!           calls VDISTR2
! vdistr => ivdistrg2 initializes vx, vy co-ordinates for 2-1/2d code,
!           with maxwellian velocity distribution with drift.
!           calls VDISTR2, or VDISTR2H
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: march 1, 2012
!
      use input2d
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR, PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D,         &
     &VACUUM_2D, VACUUM_3D, NEUMANN_2D, NEUMANN_PERIODIC_2D,            &
     &DIRICHLET_NEUMANN_PERIODIC_2D
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D
      public :: VACUUM_2D, VACUUM_3D, NEUMANN_2D, NEUMANN_PERIODIC_2D
      public :: DIRICHLET_NEUMANN_PERIODIC_2D
      public :: input2
      public :: idrun, idrun0, idcode, indx, indy, npx, npy, npxb, npyb
      public :: inorder, popt, dopt, djopt, nustrt, ntr
      public :: ntw, ntp, nta, ntv, nts, ntm, nte
      public :: ndw, ndp, nda, ndv, nds, ndm, nde
      public :: nsv
      public :: tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0
      public :: vdx, vdy, vdz, vtdx, vtdy, vtdz
      public :: psolve, relativity
      public :: omx, omy, omz, ci, ax, ay
      public :: ndim, ndc, movion
      public :: sortime, nplot, idpal, ndstyle, sntasks
      public :: nsrand, ndprof
      public :: ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy
      public :: monitor, itpon, v0, w0
      public :: modesxp, modesyp, modesxa, modesya, modesxe, modesye
      public :: ions2
      public :: ntd, ntj, npxi, npyi, npxbi, npybi
      public :: ndd, ndj
      public :: qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0
      public :: vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, sortimi
      public :: nsrandi, ndprofi
      public :: ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi
      public :: ionoff
      public :: modesxd, modesyd, modesxj, modesyj
      public :: pot2d, vpot2d, em2d
      public :: t0, ceng, indian, rlprec, inprec
      public :: nprec, fpname, narec, faname, nerec, fename
      public :: nlsrerec
      public :: den2d, vcur2d
      public :: ndrec, fdname, njrec, fjname
      public :: distr, ldistr, fdistr, vdistr, revdistr
      public :: launchlas, lase0, laspol, lasstart,lasend,launchtype
      public :: launchdir
      public :: lask, ntlas, launchshape,lasw0,lasfocus,antlength
      public :: remass,affpin, laseronly, dampfield
      public :: rayl, gamma, smoothtype, smoothfactor
      public :: ntle1,ntle2,ntle3,ntlb1,ntlb2,ntlb3,ntqe,ntqi
      public :: ntcu,ntemb1,ntemb2,ntemb3,nteme1,nteme2,nteme3
      public :: edensity,smoothmask
      public :: plasmabstart,plasmabend,plasmabystart,plasmabyend
      public :: dampystart, dampyend, dampdrop,dampe1,dampe2,dampe3
      public :: npxbg, npybg, npxg, npyg
      public :: ntpart
      public :: beamstart,beamend,beamystart,beamyend
      public :: solvertype
      public :: xymaxwel
      public :: ntej1, ntej2, ntej3, ntlej1, ntlej2, ntlej3
      public :: nttj1, nttj2, nttj3, ntltj1, ntltj2, ntltj3
      public :: mzf, parmax, yee
      public :: vpml, xtrue, ytrue
!
! define interface to original Fortran77 procedures
      interface
         subroutine DISTR2(part,vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx,ny,&
     &ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
!        real, dimension(*) :: part
         real:: part
         end subroutine
      end interface
      interface
         subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,n&
     &op,nx,ny,ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(*) :: part
         real:: part
         end subroutine
      end interface
      interface
         subroutine LDISTR2(part,anlx,anly,npx,npy,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: anlx, anly
!        real, dimension(*) :: part
         real:: part
         end subroutine
      end interface
      interface
         subroutine FDISTR2(part,fnx,argx1,argx2,argx3,fny,argy1,argy2,a&
     &rgy3,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
!        real, dimension(*) :: part
         real:: part
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine FDISTR2B(part,fnx,argx1,argx2,argx3,fny,argy1,argy2,a&
     &rgy3,npx,npy,idimp,nop,nx,ny,ipbc,ierr,plasmabstart,plasmabend,plas&
     &mabystart,plasmabyend)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc, ierr
         integer :: plasmabstart,plasmabend,plasmabystart,plasmabyend
         real :: argx1, argx2, argx3, argy1, argy2, argy3
!        real, dimension(*) :: part
         real:: part
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine RDISTR2(part,fnx,argx1,argx2,argx3,fny,argy1,argy2,a&
     &rgy3,npx,npy,idimp,nop,nx,ny,ipbc)
         implicit none
         integer :: npx, npy, idimp, nop, nx, ny, ipbc
         real :: argx1, argx2, argx3, argy1, argy2, argy3
!        real, dimension(*) :: part
         real :: part
         real, external :: fnx, fny
         end subroutine
      end interface
      interface
         subroutine VDISTR2(part,vtx,vty,vdx,vdy,idimp,nop)
         implicit none
         integer :: idimp, nop
         real :: vtx, vty, vdx, vdy
!        real, dimension(*) :: part
         real :: part
         end subroutine
      end interface
      interface
         subroutine VDISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,idimp,nop)
         implicit none
         integer :: idimp, nop
         real :: vtx, vty, vtz, vdx, vdy, vdz
!        real, dimension(*) :: part
         real :: part
         end subroutine
      end interface
      interface
         subroutine BDISTR2L(part,bx,by,bz,qbm,idimp,nop,nx,ny,nxv)
         implicit none
         integer :: idimp, nop, nx, ny, nxv
         real :: qbm
         real, dimension(idimp,nop) :: part
         real, dimension(nxv,ny) :: bx, by, bz
         end subroutine
      end interface
      interface
         subroutine GBDISTR2L(part,bxy,qbm,idimp,nop,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm
         real, dimension(idimp,nop) :: part
         real :: bxy
         end subroutine
      end interface
      interface
         subroutine GBZDISTR2L(part,bz,qbm,idimp,nop,nx,ny,nxv,nyv,ipbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm
         real, dimension(idimp,nop) :: part
         real :: bz
         end subroutine
      end interface
      interface
         subroutine GRBDISTR2L(part,bxy,qbm,ci,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, ci
         real, dimension(idimp,nop) :: part
         real :: bxy
         end subroutine
      end interface
      interface
         subroutine GRBZDISTR2L(part,bz,qbm,ci,idimp,nop,nx,ny,nxv,nyv,i&
     &pbc)
         implicit none
         integer :: idimp, nop, nx, ny, nxv, nyv, ipbc
         real :: qbm, ci
         real, dimension(idimp,nop) :: part
         real :: bz
         end subroutine
      end interface
      interface
         subroutine REVDISTR2(part,vtx,vty,vdx,vdy,f,nstart,npxy,idimp,n&
     &op)
         implicit none
         integer :: nstart, npxy, idimp, nop
         real :: vtx, vty, vdx, vdy, f
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
      interface
         subroutine REVDISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,f,nstart,npx&
     &y,idimp,nop)
         implicit none
         integer :: nstart, npxy, idimp, nop
         real :: vtx, vty, vtz, vdx, vdy, vdz, f
         real, dimension(idimp,nop) :: part
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface distr
         module procedure idistr2
         module procedure idistrh2
         module procedure ibdistr2
         module procedure irbdistr2
      end interface
!
      interface ldistr
         module procedure ildistr2
      end interface
!
      interface fdistr
         module procedure ifdistr2
      end interface
!
      interface vdistr
         module procedure ivdistr2
!        module procedure ivdistrh2
         module procedure ivdistrg2
      end interface
!
      interface revdistr
         module procedure irevdistrg2
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine idistr2(part,nstart,nop,vtx,vty,vdx,vdy,npx,npy,nx,n&
     &y,ipbc)
! calculates initial particle co-ordinates and velocities in 2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: nstart, nop, npx, npy, nx, ny, ipbc
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call DISTR2(part(1,nstart),vtx,vty,vdx,vdy,npx,npy,idimp,nop,nx&
     &,ny,ipbc)
         end subroutine idistr2
!
         subroutine idistrh2(part,nstart,nop,vtx,vty,vtz,vdx,vdy,vdz,npx&
     &,npy,nx,ny,ipbc)
! calculates initial particle co-ordinates and velocities in 2-1/2d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: nstart, nop, npx, npy, nx, ny, ipbc
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call DISTR2H(part(1,nstart),vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idi&
     &mp,nop,nx,ny,ipbc)
         end subroutine idistrh2
!
         subroutine ildistr2(part,nstart,nop,anlx,anly,npx,npy,nx,ny,ipb&
     &c)
! calculates initial particle co-ordinates in 2d
! with bi-linear density profile
         implicit none
         integer :: nstart, nop, npx, npy, nx, ny, ipbc
         real :: anlx, anly
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call LDISTR2(part(1,nstart),anlx,anly,npx,npy,idimp,nop,nx,ny,i&
     &pbc)
         end subroutine ildistr2
!
         subroutine ifdistr2(part,nstart,nop,ampx,scalex,shiftx,ampy,sca&
     &ley,shifty,npx,npy,nx,ny,ipbc,ndpro,nsran,typeindex)
! calculates initial particle co-ordinates in 2d
! with various density profiles
         implicit none
         integer :: nstart, nop, npx, npy, nx, ny, ipbc, ndpro, nsran
         integer, optional :: typeindex
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, ierr
         real sxi, syi, zero
         real, external :: FLDISTR1, FLDISTR1B, FSDISTR1, FGDISTR1, FHDISTR1
         idimp = size(part,1)
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call FDISTR2(part(1,nstart),FLDISTR1,zero,zero,zero,FLDISTR1&
     &,zero,zero,zero,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
            if (nsran /= 0) then
               call RDISTR2(part(1,nstart),FLDISTR1,zero,zero,zero,FLDIS&
     &TR1,zero,zero,zero,npx,npy,idimp,nop,nx,ny,ipbc)
            endif
! linear density
         else if (ndpro==1) then
            call FDISTR2(part(1,nstart),FLDISTR1,ampx,sxi,shiftx,FLDISTR&
     &1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
            if (nsran /= 0) then
               call RDISTR2(part(1,nstart),FLDISTR1,ampx,sxi,shiftx,FLDI&
     &STR1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call FDISTR2(part(1,nstart),FSDISTR1,ampx,sxi,shiftx,FSDISTR&
     &1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
            if (nsran /= 0) then
               call RDISTR2(part(1,nstart),FSDISTR1,ampx,sxi,shiftx,FSDI&
     &STR1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc)
            endif
! gaussian density
         else if (ndpro==3) then
            call FDISTR2(part(1,nstart),FGDISTR1,ampx,sxi,shiftx,FGDISTR&
     &1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
            if (nsran /= 0) then
               call RDISTR2(part(1,nstart),FGDISTR1,ampx,sxi,shiftx,FGDI&
     &STR1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call FDISTR2(part(1,nstart),FHDISTR1,ampx,sxi,shiftx,FHDISTR&
     &1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc,ierr)
            if (nsran /= 0) then
               call RDISTR2(part(1,nstart),FHDISTR1,ampx,sxi,shiftx,FHDI&
     &STR1,ampy,syi,shifty,npx,npy,idimp,nop,nx,ny,ipbc)
            endif
         else if (ndpro==5) then
            if (typeindex == 1) then
            call FDISTR2B(part(1,nstart),FLDISTR1,zero,zero,zero,FLDISTR1&
     &,zero,zero,zero,npx,npy,idimp,nop,nx,ny,ipbc,ierr,plasmabstart,plasm&
     &abend,plasmabystart,plasmabyend)
            elseif (typeindex == 2) then
            call FDISTR2B(part(1,nstart),FLDISTR1,zero,zero,zero,FLDISTR1&
     &,zero,zero,zero,npx,npy,idimp,nop,nx,ny,ipbc,ierr,beamstart,beamend&
     ,beamystart,beamyend)
            else
            print *,'typeindex incorrect at fdistr'
            endif
            if (nsran /= 0) then
               call RDISTR2(part(1,nstart),FLDISTR1,zero,zero,zero,FLDIS&
     &TR1,zero,zero,zero,npx,npy,idimp,nop,nx,ny,ipbc)
            !print*, 'here'
            endif
         endif
         end subroutine ifdistr2
!
         subroutine ivdistr2(part,nstart,nop,vtx,vty,vdx,vdy)
! calculates initial particle velocities in 2d
! with maxwellian velocity with drift
         implicit none
         integer :: nstart, nop
         real :: vtx, vty, vdx, vdy
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call VDISTR2(part(1,nstart),vtx,vty,vdx,vdy,idimp,nop)
         end subroutine ivdistr2
!
         subroutine ivdistrh2(part,nstart,nop,vtx,vty,vtz,vdx,vdy,vdz)
! calculates initial particle velocities in 2-1/2d
! with maxwellian velocity with drift
         implicit none
         integer :: nstart, nop
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp
         idimp = size(part,1)
         call VDISTR2H(part(1,nstart),vtx,vty,vtz,vdx,vdy,vdz,idimp,nop)
         end subroutine ivdistrh2
!
         subroutine ivdistrg2(part,nstart,nop,vtx,vty,vtz,vdx,vdy,vdz,nd&
     &im)
! calculates initial particle velocities in 2d or 2-1/2d
! with maxwellian velocity with drift
         implicit none
         integer :: nstart, nop
         integer, optional :: ndim
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, nd
         idimp = size(part,1); nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call VDISTR2(part(1,nstart),vtx,vty,vdx,vdy,idimp,nop)
         case (3)
            call VDISTR2H(part(1,nstart),vtx,vty,vtz,vdx,vdy,vdz,idimp,n&
     &op)
         end select
         end subroutine ivdistrg2
!
         subroutine ibdistr2(part,bxy,nop,qbm,nx,ny,ipbc,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for 2d
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: bxy
! local data
         integer :: idimp, nxv, nyv, order
         idimp = size(part,1)
         nxv = size(bxy,2); nyv = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call GBZDISTR2L(part,bxy(1,1,1),qbm,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
            else if (order==CUBIC) then
               call GBZDISTR2L(part,bxy(1,3,3),qbm,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
            else
               call GBZDISTR2L(part,bxy(1,2,2),qbm,idimp,nop,nx,ny,nxv,n&
     &yv,ipbc)
            endif
         case (3)
            if (order==LINEAR) then
               call GBDISTR2L(part,bxy(1,1,1),qbm,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
            else if (order==CUBIC) then
               call GBDISTR2L(part,bxy(1,3,3),qbm,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
            else
               call GBDISTR2L(part,bxy(1,2,2),qbm,idimp,nop,nx,ny,nxv,ny&
     &v,ipbc)
            endif
         end select
         end subroutine ibdistr2
!
         subroutine irbdistr2(part,bxy,nop,qbm,ci,nx,ny,ipbc,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for relativistic 2d
         implicit none
         integer :: nop, nx, ny, ipbc
         integer, optional :: inorder
         real :: qbm, ci
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: bxy
! local data
         integer :: idimp, nxv, nyv, order
         idimp = size(part,1)
         nxv = size(bxy,2); nyv = size(bxy,3)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         select case(size(bxy,1))
         case (1)
            if (order==LINEAR) then
               call GRBZDISTR2L(part,bxy(1,1,1),qbm,ci,idimp,nop,nx,ny,n&
     &xv,nyv,ipbc)
            else if (order==CUBIC) then
               call GRBZDISTR2L(part,bxy(1,3,3),qbm,ci,idimp,nop,nx,ny,n&
     &xv,nyv,ipbc)
            else
               call GRBZDISTR2L(part,bxy(1,2,2),qbm,ci,idimp,nop,nx,ny,n&
     &xv,nyv,ipbc)
            endif
         case (3)
            if (order==LINEAR) then
               call GRBDISTR2L(part,bxy(1,1,1),qbm,ci,idimp,nop,nx,ny,nx&
     &v,nyv,ipbc)
            else if (order==CUBIC) then
               call GRBDISTR2L(part,bxy(1,3,3),qbm,ci,idimp,nop,nx,ny,nx&
     &v,nyv,ipbc)
            else
               call GRBDISTR2L(part,bxy(1,2,2),qbm,ci,idimp,nop,nx,ny,nx&
     &v,nyv,ipbc)
            endif
         end select
         end subroutine irbdistr2
!
         subroutine irevdistrg2(part,nstart,npxy,nop,vtx,vty,vtz,vdx,vdy&
     &,vdz,f,ndim)
! recalculates initial particle velocities in 2d or 2-1/2d
! with maxwellian velocity with drift
         implicit none
         integer :: nstart, npxy, nop
         integer, optional :: ndim
         real :: vtx, vty, vtz, vdx, vdy, vdz, f
         real, dimension(:,:), pointer :: part
! local data
         integer :: idimp, nd
         idimp = size(part,1); nd = 3
         if (present(ndim)) nd = ndim
         select case(nd)
         case (2)
            call REVDISTR2(part,vtx,vty,vdx,vdy,f,nstart,npxy,idimp,nop)
         case (3)
            call REVDISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,f,nstart,npxy,i&
     &dimp,nop)
         end select
         end subroutine irevdistrg2
!
      end module init2d
