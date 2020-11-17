!-----------------------------------------------------------------------
!
      module hdsimul2d
! Higher level subroutines for darwin fields with hamiltonian
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: december 13, 2007
      use simul2d
      use emsimul2d
      use ehdpush2d
      implicit none
      private
      public :: restart_open, restart_bwrite, restart_dwrite
      public :: restart_bread, restart_dread
      public :: retardg, dpostg, djpostg
      public :: phasediag
      public :: cr8pcmg, dhjpostg, hpushg, hpushxhg, hpushxg
!
      contains
!
         subroutine cr8pcmg(part,axy,nop,qbm,ci,relativity,inorder)
! calculate canonical momentum for particles with 2-1/2d darwin fields
         implicit none
         integer :: nop, relativity
         integer, optional :: inorder
         real :: qbm, ci
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
         if (relativity==1) then
            call rcr8pcm(part,axy,nop,qbm,ci,inorder)
         else
            call cr8pcm(part,axy,nop,qbm,inorder)
         endif
         end subroutine cr8pcmg
!
         subroutine dhjpostg(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,ci,tdhjp&
     &ost,relativity,inorder,djopt)
! deposit current density with 2-1/2d darwin fields
         implicit none
         integer :: nop, relativity
         integer, optional :: inorder, djopt
         real :: qm, qbm, dt, ci, tdhjpost
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, axy, daxy, cu
         if (relativity==1) then
            call rdhjpost(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,ci,tdhjpost&
     &,inorder,djopt)
         else
            call dhjpost(part,fxy,axy,daxy,cu,nop,qm,qbm,dt,tdhjpost,ino&
     &rder,djopt)
         endif
         end subroutine dhjpostg
!
         subroutine hpushg(part,fxy,axy,daxy,nop,qbm,dt,ci,ek,tpush,nx,n&
     &y,ipbc,relativity,inorder,popt)
! update particle canonical momentum with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, ek, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fxy, axy, daxy
         if (relativity==1) then
            call rhpush(part,fxy,axy,daxy,nop,qbm,dt,ci,ek,tpush,nx,ny,i&
     &norder,popt)
         else
            call hpush(part,fxy,axy,daxy,nop,qbm,dt,ek,tpush,nx,ny,ipbc,&
     &inorder,popt)
         endif
         end subroutine hpushg
!
         subroutine hpushxhg(part,axy,nop,qbm,dt,ci,tpush,nx,ny,ipbc,rel&
     &ativity,inorder,popt)
! update particle position to t+dt/2 with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
         if (relativity==1) then
            call rhpushxh(part,axy,nop,qbm,dt,ci,tpush,nx,ny,inorder,pop&
     &t)
         else
            call hpushxh(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc,inorder,po&
     &pt)
         endif
         end subroutine hpushxhg
!
         subroutine hpushxg(part,axy,nop,qbm,dt,ci,tpush,nx,ny,ipbc,rela&
     &tivity,inorder,popt)
! update particle position to t+dt with 2-1/2d darwin fields
         implicit none
         integer :: nop, nx, ny, ipbc, relativity
         integer, optional :: inorder, popt
         real :: qbm, dt, ci, tpush
         real, dimension(:,:), pointer :: part
         real, dimension(:,:,:), pointer :: axy
         if (relativity==1) then
            call rhpushx(part,axy,nop,qbm,dt,ci,tpush,nx,ny,inorder,popt&
     &)
         else
            call hpushx(part,axy,nop,qbm,dt,tpush,nx,ny,ipbc,inorder,pop&
     &t)
         endif
         end subroutine hpushxg
!
      end module hdsimul2d