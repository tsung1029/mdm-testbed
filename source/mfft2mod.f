!-----------------------------------------------------------------------
!
      module mfft2d
!
! Fortran90 interface to 2d PIC Fortran77 library mfft2lib.f
! mfft2mod.f contains multi-tasking interface procedures to perform
!            ffts:
!            defines module mfft2d
! fft => imfft2rx performs 2d real to complex fft and its inverse.
!        calls MFFT2RX
! fft => imfft2r2 performs multiple 2d complex to real fft for 2
!        component vector arrays.
!        calls MFFT2R2
! fft => imfft2r3 performs multiple 2d real to complex fft and its
!        inverse for 1, 2, or 3 component vector arrays.
!        calls MFFT2RX, MFFT2R2, or MFFT2R3
! fftn => imfft2rn performs multiple 2d real to complex fft and its
!         inverse for scalar or 2 and 3 component vector arrays.
!         calls MFFT2RN
! fsst => imfsst2rx performs 2d scalar real sine-sine transform.
!         calls MFSST2RX
! fcct => imfcct2rx performs 2d scalar real cosine-cosine transform.
!         calls MFCCT2RX
! fcst => imfcst2r3 performs multiple 2d vector real cosine-sine
!         transforms for 2 or 3 component vector arrays.
!         calls MFCST2R2, or MFCST2R3
! fsct => imfsct2r3 performs multiple 2d vector real sine-cosine
!         transforms for 2 or 3 component vector arrays.
!         calls MFSCT2R2, or MFSCT2R3
! fsft => imfsft2rx  performs 2d scalar real sine/periodic transform.
!         calls MFSFT2RX
! fcft => imfcft2rx  performs 2d scalar real cosine/periodic transform.
!         calls MFCFT2RX
! fcsft => imfcsft2r3 performs multiple 2d vector real cosine-sine/
!          periodic transforms for 2 or 3 component vector arrays.
!          calls  MFCSFT2R2, or MFCSFT2R3
! fscft => imfscft2r3 performs multiple 2d vector real sine-cosine/
!          periodic transforms for 2 or 3 component vector arrays.
!          calls  MFSCFT2R2, or MFSCFT2R3
! fdsft => imfdsft2rx performs 2d scalar real half sine/periodic
!          transform.
!          calls MFDSFT2RX
! fdcft => imfdcft2rx performs 2d scalar real half cosine/periodic
!          transform.
!          calls MFDCFT2RX
! fdcsft => imfdcsft2r3 performs multiple 2d vector real half
!           cosine-sine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls MFDCSFT2R2, or MFDCSFT2R3
! fdscft => imfdscft2r3 performs multiple 2d vector real half
!           sine-cosine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls MFDSCFT2R2, or MFDSCFT2R3
! written by viktor k. decyk, ucla
! copyright 2002, regents of the university of california
! update: december 4, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use fft2d, only: wtimer, fft_init, fftc_init, fst_init, fdt_init, &
     &iwfft2rn
      use mp0d, only: ntasks
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC
      public :: fft_init, fft, fftn, fftc_init
      public :: fst_init, fsst, fcct, fcst, fsct
      public :: fsft, fcft, fcsft, fscft
      public :: fdt_init, fdsft, fdcft, fdcsft, fdscft
      public :: iwfft2rn
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n&
     &xyhd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n&
     &xyhd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n&
     &xyhd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFFT2RN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,ndim&
     &,nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, ndim, nxhyd, nxyhd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         integer, dimension(nmt) :: nxyip, iftask
         complex, dimension(ndim,nxhd,nmt+1) :: ss
         end subroutine
      end interface
      interface
         subroutine MFSST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFSCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCST2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFSCT2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCST2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFSCT2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFSFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,n&
     &xhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,n&
     &xhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCSFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFSCFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFCSFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFSCFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFDSFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd&
     &,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFDCFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd&
     &,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFDCSFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFDSCFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFDCSFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
      interface
         subroutine MFDSCFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
         integer :: nmt, ierr
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         integer, dimension(nmt) :: nxyip, iftask
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft
         module procedure imfft2rx
         module procedure imfft2r2
         module procedure imfft2r3
      end interface
!
      interface fftn
         module procedure imfft2rn
      end interface
!
      interface fsst
         module procedure imfsst2rx
      end interface
!
      interface fcct
         module procedure imfcct2rx
      end interface
!
      interface fcst
         module procedure imfcst2r3
      end interface
!
      interface fsct
         module procedure imfsct2r3
      end interface
!
      interface fsft
         module procedure imfsft2rx
      end interface
!
      interface fcft
         module procedure imfcft2rx
      end interface
!
      interface fcsft
         module procedure imfcsft2r3
      end interface
!
      interface fscft
         module procedure imfscft2r3
      end interface
!
      interface fdsft
         module procedure imfdsft2rx
      end interface
!
      interface fdcft
         module procedure imfdcft2rx
      end interface
!
      interface fdcsft
         module procedure imfdcsft2r3
      end interface
!
      interface fdscft
         module procedure imfdscft2r3
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine imfft2rx(f,isign,mixup,sct,tfft,indx,indy,inorder)
! perform multi-tasking 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, nmt, ltime, order, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         if (order==LINEAR) then
            call MFFT2RX(f(1,1),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd,nxyip,iftask,nmt,ierr)
         else if (order==CUBIC) then
            call MFFT2RX(f(3,3),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd,nxyip,iftask,nmt,ierr)
         else
            call MFFT2RX(f(2,2),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfft2rx
!
         subroutine imfft2r2(f,mixup,sct,tfft,indx,indy,inorder)
! perform multi-tasking 2d vector real to complex fft
! for 2 component vectors
         implicit none
         integer :: indx, indy
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, nxhd, nyd, nxhyd, nxyhd, nmt, ltime
         integer :: order, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         if (order==LINEAR) then
            call MFFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,nxh&
     &yd,nxyhd,nxyip,iftask,nmt,ierr)
         else if (order==CUBIC) then
            call MFFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,nxh&
     &yd,nxyhd,nxyip,iftask,nmt,ierr)
         else
            call MFFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,nxh&
     &yd,nxyhd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfft2r2
!
         subroutine imfft2r3(f,isign,mixup,sct,tfft,indx,indy,inorder)
! perform multi-tasking 2d vector real to complex fft
         implicit none
         integer :: isign, indx, indy
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, nmt, ltime, order, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (order==LINEAR) then
               call MFFT2RX(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            else if (order==CUBIC) then
               call MFFT2RX(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            else
               call MFFT2RX(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            endif
         case (2)
            if (order==LINEAR) then
               call MFFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            else if (order==CUBIC) then
               call MFFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            else
               call MFFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            endif
         case (3)
            if (order==LINEAR) then
               call MFFT2R3(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            else if (order==CUBIC) then
               call MFFT2R3(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            else
               call MFFT2R3(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
            endif
         end select
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfft2r3
!
         subroutine imfft2rn(f,isign,mixup,sct,tfft,indx,indy,inorder)
! perform multi-tasking 2d vector real to complex fft
! for n component vectors
         implicit none
         integer :: isign, indx, indy
         integer, optional :: inorder
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: ndim, nxhd, nyd, nxhyd, nxyhd, nmt, ltime, order
         integer :: ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         complex, dimension(size(f,1),size(f,2)/2,ntasks+1) :: ss
         ndim = size(f,1); nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         nmt = ntasks
         order = QUADRATIC
         if (present(inorder)) order = inorder
! initialize timer
         call wtimer(tf,ltime,-1)
         if (order==LINEAR) then
            call MFFT2RN(f(1,1,1),ss,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &ndim,nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
         else if (order==CUBIC) then
            call MFFT2RN(f(1,3,3),ss,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &ndim,nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
         else
            call MFFT2RN(f(1,2,2),ss,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &ndim,nxhyd,nxyhd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfft2rn
!
         subroutine imfsst2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d scalar real sine-sine transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call MFSST2RX(f(1,1),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd,nxyip,iftask,nmt,ierr)
         else if (inorder==CUBIC) then
            call MFSST2RX(f(3,3),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd,nxyip,iftask,nmt,ierr)
         else
            call MFSST2RX(f(2,2),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfsst2rx
!
         subroutine imfcct2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d scalar real cosine-cosine transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call MFCCT2RX(f(1,1),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd,nxyip,iftask,nmt,ierr)
         else if (inorder==CUBIC) then
            call MFCCT2RX(f(3,3),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd,nxyip,iftask,nmt,ierr)
         else
            call MFCCT2RX(f(2,2),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfcct2rx
!
         subroutine imfcst2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d vector real cosine-sine transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call MFCST2R2(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFCST2R2(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFCST2R2(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call MFCST2R3(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFCST2R3(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFCST2R3(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfcst2r3
!
         subroutine imfsct2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d vector real sine-cosine transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call MFSCT2R2(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFSCT2R2(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFSCT2R2(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call MFSCT2R3(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFSCT2R3(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFSCT2R3(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfsct2r3
!
         subroutine imfsft2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d scalar real sine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call MFSFT2RX(f(1,1),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else if (inorder==CUBIC) then
            call MFSFT2RX(f(3,3),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else
            call MFSFT2RX(f(2,2),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfsft2rx
!
         subroutine imfcft2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d scalar real cosine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call MFCFT2RX(f(1,1),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else if (inorder==CUBIC) then
            call MFCFT2RX(f(3,3),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else
            call MFCFT2RX(f(2,2),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfcft2rx
!
         subroutine imfcsft2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d vector real cosine-sine/periodic transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call MFCSFT2R2(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFCSFT2R2(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFCSFT2R2(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call MFCSFT2R3(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFCSFT2R3(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFCSFT2R3(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfcsft2r3
!
         subroutine imfscft2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform multi-tasking 2d vector real sine-cosine/periodic transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nmt, ltime, inorder, ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call MFSCFT2R2(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFSCFT2R2(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFSCFT2R2(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call MFSCFT2R3(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFSCFT2R3(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFSCFT2R3(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfscft2r3
!
         subroutine imfdsft2rx(f,isign,mixup,sctd,sctdx,tfft,indx,indy,o&
     &rder)
! perform multi-tasking 2d scalar real half sine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, nmt, ltime, inorder
         integer :: ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call MFDSFT2RX(f(1,1),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else if (inorder==CUBIC) then
            call MFDSFT2RX(f(3,3),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else
            call MFDSFT2RX(f(2,2),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfdsft2rx
!
         subroutine imfdcft2rx(f,isign,mixup,sctd,sctdx,tfft,indx,indy,o&
     &rder)
! perform multi-tasking 2d scalar real half cosine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, nmt, ltime, inorder
         integer :: ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call MFDCFT2RX(f(1,1),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else if (inorder==CUBIC) then
            call MFDCFT2RX(f(3,3),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         else
            call MFDCFT2RX(f(2,2),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfdcft2rx
!
         subroutine imfdcsft2r3(f,isign,mixup,sctd,sctdx,tfft,indx,indy,&
     &order)
! perform multi-tasking 2d vector real reak cosine-sine/periodic
! transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, nmt, ltime, inorder
         integer :: ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call MFDCSFT2R2(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFDCSFT2R2(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFDCSFT2R2(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call MFDCSFT2R3(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFDCSFT2R3(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFDCSFT2R3(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfdcsft2r3
!
         subroutine imfdscft2r3(f,isign,mixup,sctd,sctdx,tfft,indx,indy,&
     &order)
! perform multi-tasking 2d vector real half sine-cosine/periodic
! transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, nmt, ltime, inorder
         integer :: ierr
         real :: tf
         integer, dimension(ntasks) :: nxyip, iftask
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         nmt = ntasks
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call MFDSCFT2R2(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFDSCFT2R2(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFDSCFT2R2(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call MFDSCFT2R3(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else if (inorder==CUBIC) then
               call MFDSCFT2R3(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            else
               call MFDSCFT2R3(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            stop
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine imfdscft2r3
!
      end module mfft2d
