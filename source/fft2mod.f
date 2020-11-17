!-----------------------------------------------------------------------
!
      module fft2d
!
! Fortran90 interface to 2d PIC Fortran77 library fft2lib.f
! fft2mod.f contains interface procedures to perform ffts:
!           defines module fft2d
! fft_init => iwfft2rinit initializes 2d real to complex fft tables.
!             calls WFFT2RINIT
! fft => iwfft2rx performs 2d real to complex fft and its inverse.
!        calls WFFT2RX
! fft => iwfft2r2 performs multiple 2d complex to real fft for 2
!        component vector arrays.
!        calls WFFT2R2
! fft => iwfft2r3 performs multiple 2d real to complex fft and its
!        inverse for 1, 2, or 3 component vector arrays.
!        calls WFFT2RX, WFFT2R2, or WFFT2R3
! fft => ifft2c performs 2d complex to complex fft and its inverse.
!        calls FFT2C
! fftn => iwfft2rn performs multiple 2d real to complex fft and its
!         inverse for scalar or 2 and 3 component vector arrays.
!         calls WFFT2RN
! fftc_init => ifft2cinit initializes 2d complex to complex fft tables.
!              calls FFT2C
! fst_init => iwfst2rinit initializes 2d real sine-cosine transform.
!             tables.
!             calls WFST2RINIT
! fsst => iwfsst2rx performs 2d scalar real sine-sine transform.
!         calls WFSST2RX
! fcct => iwfcct2rx performs 2d scalar real cosine-cosine transform.
!         calls WFCCT2RX
! fcst => iwfcst2r3 performs multiple 2d vector real cosine-sine
!         transforms for 2 or 3 component vector arrays.
!         calls  WFCST2R2, or WFCST2R3
! fsct => iwfsct2r3 performs multiple 2d vector real sine-cosine
!         transforms for 2 or 3 component vector arrays.
!         calls WFSCT2R2, or WFSCT2R3
! fsft => iwfsft2rx  performs 2d scalar real sine/periodic transform.
!         calls WFSFT2RX
! fcft => iwfcft2rx  performs 2d scalar real cosine/periodic transform.
!         calls WFCFT2RX
! fcsft => iwfcsft2r3 performs multiple 2d vector real cosine-sine/
!          periodic transforms for 2 or 3 component vector arrays.
!          calls  WFCSFT2R2, or WFCSFT2R3
! fscft => iwfscft2r3 performs multiple 2d vector real sine-cosine/
!          periodic transforms for 2 or 3 component vector arrays.
!          calls  WFSCFT2R2, or WFSCFT2R3
! fdt_init => iwfdt2rinit initializes 2d real half sine-cosine/periodic
!             transform tables.
!             calls WFDT2RINIT
! fdsft => iwfdsft2rx performs 2d scalar real half sine/periodic
!          transform.
!          calls WFDSFT2RX
! fdcft => iwfdcft2rx performs 2d scalar real half cosine/periodic
!          transform.
!          calls WFDCFT2RX
! fdcsft => iwfdcsft2r3 performs multiple 2d vector real half
!           cosine-sine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls WFDCSFT2R2, or WFDCSFT2R3
! fdscft => iwfdscft2r3 performs multiple 2d vector real half
!           sine-cosine/periodic transforms for 2 or 3 component vector
!           arrays.
!           calls WFDSCFT2R2, or WFDSCFT2R3
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: september 17, 2009
!
      use globals, only: LINEAR, QUADRATIC, CUBIC
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC
      public :: wtimer
      public :: fft_init, fft, fftn, fftc_init
      public :: fst_init, fsst, fcct, fcst, fsct
      public :: fsft, fcft, fcsft, fscft
      public :: fdt_init, fdsft, fdcft, fdcsft, fdscft
      public :: iwfft2rn
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine FFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nx&
     &yhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine FFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nx&
     &yhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine FFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nx&
     &yhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine FFT2C(f,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxyhd&
     &)
         integer :: isign, indx, indy, nxd, nyd, nxyd, nxyhd
!        complex, dimension(*) :: f
         complex :: f
         integer, dimension(nxyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         implicit none
         integer :: indx, indy
         integer :: nxhyd, nxyhd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      !interface
       !  subroutine WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n&
     !&xyhd)
      !   implicit none
      !   integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(*) :: f
      !   real :: f
      !   integer, dimension(nxhyd) :: mixup
      !   complex, dimension(nxyhd) :: sct
      !   end subroutine
      !end interface
      !interface
      !   subroutine WFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n&
     !&xyhd)
       !  implicit none
      !   integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(*) :: f
       !  real :: f
       !  integer, dimension(nxhyd) :: mixup
       !  complex, dimension(nxyhd) :: sct
       !  end subroutine
      !end interface
      !interface
       !  subroutine WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n&
     !&xyhd)
      !   implicit none
      !   integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
!        real, dimension(*) :: f
      !   real :: f
      !   integer, dimension(nxhyd) :: mixup
      !   complex, dimension(nxyhd) :: sct
      !   end subroutine
      !end interface
      interface
         subroutine WFFT2RN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,ndim&
     &,nxhyd,nxyhd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, ndim, nxhyd, nxyhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyhd) :: sct
         complex, dimension(ndim,nxhd) :: ss
         end subroutine
      end interface
      interface
         subroutine WFFT2CINIT(mixup,sct,indx,indy,nxyd,nxyhd)
         implicit none
         integer :: indx, indy
         integer :: nxyd, nxyhd
         integer, dimension(nxyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WFFT2C(f,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxyh&
     &d)
         integer :: isign, indx, indy, nxd, nyd, nxyd, nxyhd
!        complex, dimension(*) :: f
         complex :: f
         integer, dimension(nxyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WFFT2CT(f,g,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,n&
     &xyhd)
         integer :: isign, indx, indy, nxd, nyd, nxyd, nxyhd
!        complex, dimension(*) :: f
         complex :: f
         complex, dimension(nyd,nxd) :: g
         integer, dimension(nxyd) :: mixup
         complex, dimension(nxyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WFFT2RCX(f,fc,g,isign,mixup,sct,indx,indy,nxv,nyv,nx&
     &hd,nyd,nxh1d,nxhyd,nxhyhd)
         integer :: isign, indx, indy, nxv, nyv, nxhd, nyd, nxh1d
         integer :: nxhyd, nxhyhd
!        complex, dimension(*) :: f
         complex :: f
         complex, dimension(nxhd,nyd) :: fc
         complex, dimension(nyd,nxh1d) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxhyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WFFT2RCN(f,fc,g,isign,mixup,sct,indx,indy,ndim,nxv,n&
     &yv,nxhd,nyd,nxh1d,nxhyd,nxhyhd)
         integer :: isign, indx, indy, ndim, nxv, nyv, nxhd, nyd, nxh1d
         integer :: nxhyd, nxhyhd
!        complex, dimension(*) :: f
         complex :: f
         complex, dimension(nxhd,ndim,nyd) :: fc
         complex, dimension(nyd,ndim,nxh1d) :: g
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxhyhd) :: sct
         end subroutine
      end interface
      interface
         subroutine WFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
         implicit none
         integer :: indx, indy
         integer :: nxhyd, nxyd
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCST2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSCT2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCST2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSCT2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd&
     &,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,n&
     &xhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,n&
     &xhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCSFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSCFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFCSFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFSCFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,&
     &nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         end subroutine
      end interface
      interface
         subroutine WFDT2RINIT(sctdx,indx,nxd)
         implicit none
         integer :: indx, nxd
         complex, dimension(nxd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WFDSFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd&
     &,nyhd,nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WFDCFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd&
     &,nyhd,nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WFDCSFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WFDSCFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(2,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WFDCSFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         end subroutine
      end interface
      interface
         subroutine WFDSCFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
         implicit none
         integer :: isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
!        real, dimension(*) :: f
         real :: f
         complex, dimension(3,2*nxhd) :: ss
         integer, dimension(nxhyd) :: mixup
         complex, dimension(nxyd) :: sctd
         complex, dimension(2*nxhd) :: sctdx
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft_init
!        module procedure ifft2rxinit
         module procedure iwfft2rinit
      end interface
!
      interface fft
!        module procedure ifft2rx
!        module procedure ifft2r2
!        module procedure ifft2r3
         module procedure iwfft2rx
         module procedure iwfft2r2
         module procedure iwfft2r3
         module procedure iwfft2cr3
         module procedure ifft2c
      end interface
!
      interface fftn
         module procedure iwfft2rn
      end interface
!
      interface fftc_init
         module procedure ifft2cinit
      end interface
!
      interface fst_init
         module procedure iwfst2rinit
      end interface
!
      interface fsst
         module procedure iwfsst2rx
      end interface
!
      interface fcct
         module procedure iwfcct2rx
      end interface
!
      interface fcst
         module procedure iwfcst2r3
      end interface
!
      interface fsct
         module procedure iwfsct2r3
      end interface
!
      interface fsft
         module procedure iwfsft2rx
      end interface
!
      interface fcft
         module procedure iwfcft2rx
      end interface
!
      interface fcsft
         module procedure iwfcsft2r3
      end interface
!
      interface fscft
         module procedure iwfscft2r3
      end interface
!
      interface fdt_init
         module procedure iwfdt2rinit
      end interface
!
      interface fdsft
         module procedure iwfdsft2rx
      end interface
!
      interface fdcft
         module procedure iwfdcft2rx
      end interface
!
      interface fdcsft
         module procedure iwfdcsft2r3
      end interface
!
      interface fdscft
         module procedure iwfdscft2r3
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ifft2rxinit(mixup,sct,indx,indy)
! initialize 2d real to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, nxhd = 1, nyd = 1, nxhyd, nxyhd
         real :: f
         nxhyd = size(mixup); nxyhd = size(sct)
         call FFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd)
         end subroutine ifft2rxinit
!
         subroutine ifft2rx(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FFT2RX(f(1,1),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         else if (inorder==CUBIC) then
            call FFT2RX(f(3,3),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         else
            call FFT2RX(f(2,2),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,&
     &nxyhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft2rx
!
         subroutine ifft2r2(f,mixup,sct,tfft,indx,indy,order)
! perform 2d vector real to complex fft for 2 component vectors
         implicit none
         integer :: indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,nxhy&
     &d,nxyhd)
         else if (inorder==CUBIC) then
            call FFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,nxhy&
     &d,nxyhd)
         else
            call FFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,nxhy&
     &d,nxyhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft2r2
!
         subroutine ifft2r3(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d vector real to complex fft for 3 component vectors
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (inorder==LINEAR) then
               call FFT2RX(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            else if (inorder==CUBIC) then
               call FFT2RX(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            else
               call FFT2RX(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            endif
         case (2)
            if (inorder==LINEAR) then
               call FFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            else if (inorder==CUBIC) then
               call FFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            else
               call FFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            endif
         case (3)
            if (inorder==LINEAR) then
               call FFT2R3(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            else if (inorder==CUBIC) then
               call FFT2R3(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            else
               call FFT2R3(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,n&
     &xhyd,nxyhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft2r3
!
         subroutine iwfft2rinit(mixup,sct,indx,indy)
! initialize 2d real to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhyd, nxyhd
         nxhyd = size(mixup); nxyhd = size(sct)
         call WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
         end subroutine iwfft2rinit
!
         subroutine iwfft2rx(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d scalar real to complex fft
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFFT2RX(f(1,1),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         else if (inorder==CUBIC) then
            call WFFT2RX(f(3,3),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         else
            call WFFT2RX(f(2,2),isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd&
     &,nxyhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfft2rx
!
         subroutine iwfft2r2(f,mixup,sct,tfft,indx,indy,order)
! perform 2d vector complex to real fft for 2 component vectors
         implicit none
         integer :: indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 1, nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,nxh&
     &yd,nxyhd)
         else if (inorder==CUBIC) then
            call WFFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,nxh&
     &yd,nxyhd)
         else
            call WFFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,nxh&
     &yd,nxyhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfft2r2
!
         subroutine iwfft2r3(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d vector real to complex fft for 3 component vectors
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (inorder==LINEAR) then
               call WFFT2RX(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else if (inorder==CUBIC) then
               call WFFT2RX(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else
               call WFFT2RX(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            endif
         case (2)
            if (inorder==LINEAR) then
               call WFFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else if (inorder==CUBIC) then
               call WFFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else
               call WFFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            endif
         case (3)
            if (inorder==LINEAR) then
               call WFFT2R3(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else if (inorder==CUBIC) then
               call WFFT2R3(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else
               call WFFT2R3(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfft2r3
         
         subroutine iwfft2cr3(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d vector real to complex fft for 3 component vectors
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         complex, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         select case(size(f,1))
         case (1)
            if (inorder==LINEAR) then
               call WFFT2RX(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else if (inorder==CUBIC) then
               call WFFT2RX(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else
               call WFFT2RX(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            endif
         case (2)
            if (inorder==LINEAR) then
               call WFFT2R2(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else if (inorder==CUBIC) then
               call WFFT2R2(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else
               call WFFT2R2(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            endif
         case (3)
            if (inorder==LINEAR) then
               call WFFT2R3(f(1,1,1),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else if (inorder==CUBIC) then
               call WFFT2R3(f(1,3,3),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            else
               call WFFT2R3(f(1,2,2),isign,mixup,sct,indx,indy,nxhd,nyd,&
     &nxhyd,nxyhd)
            endif
         end select
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfft2cr3
!
         subroutine iwfft2rn(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d vector real to complex fft for n component vectors
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         complex, dimension(size(f,1),size(f,2)/2) :: ss
         integer :: ndim, nxhd, nyd, nxhyd, nxyhd, ltime, inorder
         real :: tf
         ndim = size(f,1); nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFFT2RN(f(1,1,1),ss,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &ndim,nxhyd,nxyhd)
         else if (inorder==CUBIC) then
            call WFFT2RN(f(1,3,3),ss,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &ndim,nxhyd,nxyhd)
         else
            call WFFT2RN(f(1,2,2),ss,isign,mixup,sct,indx,indy,nxhd,nyd,&
     &ndim,nxhyd,nxyhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfft2rn
!
         subroutine ifft2c(f,isign,mixup,sct,tfft,indx,indy,order)
! perform 2d scalar complex to complex fft
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         complex, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxd, nyd, nxyd, nxyhd, ltime, inorder
         real :: tf
         nxd = size(f,1); nyd = size(f,2)
         nxyd = size(mixup); nxyhd = size(sct)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FFT2C(f(1,1),isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxy&
     &hd)
         else if (inorder==CUBIC) then
            call FFT2C(f(3,3),isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxy&
     &hd)
         else
            call FFT2C(f(2,2),isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxy&
     &hd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft2c
!
         subroutine ifft2cinit(mixup,sct,indx,indy)
! initialize 2d complex to complex fft
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: isign = 0, nxd = 1, nyd = 1, nxyd, nxyhd
         complex :: f
         nxyd = size(mixup); nxyhd = size(sct)
         call FFT2C(f,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxyhd)
         end subroutine ifft2cinit
!
         subroutine iwfst2rinit(mixup,sctd,indx,indy)
! initialize 2d real sine-cosine transforms
         implicit none
         integer :: indx, indy
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhyd, nxyd
         nxhyd = size(mixup); nxyd = size(sctd)
         call WFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
         end subroutine iwfst2rinit
!
         subroutine iwfsst2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d scalar real sine-sine transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFSST2RX(f(1,1),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd)
         else if (inorder==CUBIC) then
            call WFSST2RX(f(3,3),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd)
         else
            call WFSST2RX(f(2,2),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfsst2rx
!
         subroutine iwfcct2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d scalar real cosine-cosine transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyd = size(f,2)
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFCCT2RX(f(1,1),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd)
         else if (inorder==CUBIC) then
            call WFCCT2RX(f(3,3),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd)
         else
            call WFCCT2RX(f(2,2),isign,mixup,sctd,indx,indy,nxhd,nyd,nxh&
     &yd,nxyd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfcct2rx
!
         subroutine iwfcst2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d vector real cosine-sine transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call WFCST2R2(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFCST2R2(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else
               call WFCST2R2(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call WFCST2R3(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFCST2R3(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else
               call WFCST2R3(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfcst2r3
!
         subroutine iwfsct2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d vector real sine-cosine transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         integer :: nxhd, nyd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyd = size(f,3)
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call WFSCT2R2(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFSCT2R2(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else
               call WFSCT2R2(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call WFSCT2R3(f(1,1,1),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFSCT2R3(f(1,3,3),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            else
               call WFSCT2R3(f(1,2,2),isign,mixup,sctd,indx,indy,nxhd,ny&
     &d,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfsct2r3
!
         subroutine iwfsft2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d scalar real sine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFSFT2RX(f(1,1),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd)
         else if (inorder==CUBIC) then
            call WFSFT2RX(f(3,3),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd)
         else
            call WFSFT2RX(f(2,2),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfsft2rx
!
         subroutine iwfcft2rx(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d scalar real cosine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFCFT2RX(f(1,1),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd)
         else if (inorder==CUBIC) then
            call WFCFT2RX(f(3,3),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd)
         else
            call WFCFT2RX(f(2,2),ss,isign,mixup,sctd,indx,indy,nxhd,nyhd&
     &,nxhyd,nxyd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfcft2rx
!
         subroutine iwfcsft2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d vector real cosine-sine/periodic transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call WFCSFT2R2(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFCSFT2R2(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else
               call WFCSFT2R2(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call WFCSFT2R3(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFCSFT2R3(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else
               call WFCSFT2R3(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfcsft2r3
!
         subroutine iwfscft2r3(f,isign,mixup,sctd,tfft,indx,indy,order)
! perform 2d vector real sine-cosine/periodic transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call WFSCFT2R2(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFSCFT2R2(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else
               call WFSCFT2R2(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call WFSCFT2R3(f(1,1,1),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFSCFT2R3(f(1,3,3),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            else
               call WFSCFT2R3(f(1,2,2),ss,isign,mixup,sctd,indx,indy,nxh&
     &d,nyhd,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfscft2r3
!
         subroutine iwfdt2rinit(sctdx,indx)
! initialize 2d real sine-cosine/periodic transforms
         implicit none
         integer :: indx
         complex, dimension(:), pointer :: sctdx
! local data
         integer :: nxd
         nxd = size(sctdx)
         call WFDT2RINIT(sctdx,indx,nxd)
         end subroutine iwfdt2rinit
!
         subroutine iwfdsft2rx(f,isign,mixup,sctd,sctdx,tfft,indx,indy,o&
     &rder)
! perform 2d scalar real half sine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFDSFT2RX(f(1,1),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd)
         else if (inorder==CUBIC) then
            call WFDSFT2RX(f(3,3),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd)
         else
            call WFDSFT2RX(f(2,2),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfdsft2rx
!
         subroutine iwfdcft2rx(f,isign,mixup,sctd,sctdx,tfft,indx,indy,o&
     &rder)
! perform 2d scalar real half cosine/periodic transform
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, ltime, inorder
         real :: tf
         nxhd = size(f,1)/2; nyhd = size(f,2)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call WFDCFT2RX(f(1,1),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd)
         else if (inorder==CUBIC) then
            call WFDCFT2RX(f(3,3),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd)
         else
            call WFDCFT2RX(f(2,2),ss,isign,mixup,sctd,sctdx,indx,indy,nx&
     &hd,nyhd,nxhyd,nxyd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfdcft2rx
!
         subroutine iwfdcsft2r3(f,isign,mixup,sctd,sctdx,tfft,indx,indy,&
     &order)
! perform 2d vector real half cosine-sine/periodic transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call WFDCSFT2R2(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFDCSFT2R2(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else
               call WFDCSFT2R2(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call WFDCSFT2R3(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFDCSFT2R3(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else
               call WFDCSFT2R3(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfdcsft2r3
!
         subroutine iwfdscft2r3(f,isign,mixup,sctd,sctdx,tfft,indx,indy,&
     &order)
! perform 2d vector real half sine-cosine/periodic transforms
         implicit none
         integer :: isign, indx, indy
         integer, optional :: order
         real :: tfft
         real, dimension(:,:,:), pointer :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sctd, sctdx
! local data
         complex, dimension(size(f,1),size(f,2)) :: ss
         integer :: nxhd, nyhd, nxhyd, nxyd, nxd, ltime, inorder
         real :: tf
         nxhd = size(f,2)/2; nyhd = size(f,3)/2
         nxhyd = size(mixup); nxyd = size(sctd); nxd = size(sctdx)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (size(f,1)==2) then
            if (inorder==LINEAR) then
               call WFDSCFT2R2(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFDSCFT2R2(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else
               call WFDSCFT2R2(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            endif
         else if (size(f,1)==3) then
            if (inorder==LINEAR) then
               call WFDSCFT2R3(f(1,1,1),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else if (inorder==CUBIC) then
               call WFDSCFT2R3(f(1,3,3),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            else
               call WFDSCFT2R3(f(1,2,2),ss,isign,mixup,sctd,sctdx,indx,i&
     &ndy,nxhd,nyhd,nxhyd,nxyd)
            endif
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine iwfdscft2r3
!
      end module fft2d
