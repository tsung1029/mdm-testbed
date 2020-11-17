!-----------------------------------------------------------------------
!
      module fft12d
!
! Fortran90 interface to 1d PIC Fortran77 library fft1lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: march 16, 2009
!
      use globals, only: LINEAR, QUADRATIC
      use diag2d, only: wtimer
      implicit none
      private
      public :: LINEAR, QUADRATIC
      public :: fft_init, fft
!
!
! define interface to original Fortran77 procedures
      interface
         subroutine FFT1RX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         implicit none
         integer :: isign, indx, nxd, nxhd
!        real, dimension(*) :: f
         real :: f
         integer, dimension(nxhd) :: mixup
         complex, dimension(nxhd) :: t, sct
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface fft_init
         module procedure ifft1rxinit
      end interface
!
      interface fft
         module procedure ifft1rx
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine ifft1rxinit(mixup,sct,indx)
! initialize 1d real to complex fft
         implicit none
         integer :: indx
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
         integer :: isign = 0, nxd = 1, nxhd
         real :: f
         complex, dimension(1) :: t
         nxhd = size(mixup)
         call FFT1RX(f,t,isign,mixup,sct,indx,nxd,nxhd)
         end subroutine ifft1rxinit
!
         subroutine ifft1rx(f,isign,mixup,sct,tfft,indx,order)
! perform 1d scalar real to complex fft
         implicit none
         integer :: isign, indx
         integer, optional :: order
         real :: tfft
         real, dimension(:) :: f
         integer, dimension(:), pointer :: mixup
         complex, dimension(:), pointer :: sct
! local data
         integer :: nxd, nxhd, ltime, inorder
         real :: tf
         complex, dimension(size(mixup)) :: t
         nxd = size(f); nxhd = size(mixup)
         inorder = QUADRATIC
         if (present(order)) inorder = order
! initialize timer
         call wtimer(tf,ltime,-1)
         if (inorder==LINEAR) then
            call FFT1RX(f(1),t,isign,mixup,sct,indx,nxd,nxhd)
         else
            call FFT1RX(f(2),t,isign,mixup,sct,indx,nxd,nxhd)
         endif
! record time
         call wtimer(tf,ltime)
         tfft = tfft + tf
         end subroutine ifft1rx
!
      end module fft12d
