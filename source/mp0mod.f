!-----------------------------------------------------------------------
!
      module mp0d
!
! Fortran90 interface for initializing Multi-tasking library MacMP.f
! mp0mod.f contains interface procedures to support multi-tasking:
!          defines module mp0d
! mpinit initializes for multiprocessing
! written by viktor k. decyk, ucla
! copyright 2000, regents of the university of california
! update: december 5, 2009
!
      implicit none
      private
      public :: mpinit, ncpus, ntasks
      public :: setnth, setnnth, setmth, npthreads, nfthreads
      public :: nnpthreads
!
! ncpus = number of cpus found
! ntasks = number of additional tasks for threaded programming
      integer, save :: ncpus = 1, ntasks = 0
! debug for GPUs
      integer, save :: npthreads = 1, nfthreads = 1, nnpthreads = 1
!
! define interface to original Fortran77 procedures
!
      interface
         subroutine MP_INIT(nproc)
         implicit none
         integer nproc
         end subroutine
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         function mpinit(sntasks)
! initialize for multiprocessing
         implicit none
         integer mpinit
         integer, optional :: sntasks
         call MP_INIT(ncpus)
         if (ncpus==0) then
            write (2,*) 'MPLibrary not installed'
            ntasks= 0
! return the number of processors on host computer
         else
            ntasks = ncpus - 1
! set specific number of tasks if requested
            if (present(sntasks)) then
               if (sntasks >= 0) ntasks = sntasks
            endif
         endif
         mpinit = ntasks
         end function
!
         subroutine setnth(sntasks)
! initialize for multithreading for particles
         implicit none
         integer :: sntasks
         npthreads = sntasks + 1
         end subroutine setnth
!
         subroutine setnnth(sntasks)
! initialize for multithreading for particles
         implicit none
         integer :: sntasks
         nnpthreads = sntasks + 1
         end subroutine setnnth
!
         subroutine setmth(sntasks)
! initialize for multithreading for fields
         implicit none
         integer :: sntasks
         nfthreads = sntasks + 1
         end subroutine setmth
!
      end module mp0d
