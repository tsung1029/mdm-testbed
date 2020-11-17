c 2d PIC multi-tasking library for relativistic pushing particles with
c darwin vector potential using hamiltonian formulation
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: september 22, 2007
c-----------------------------------------------------------------------
      subroutine MRCR8PCM23(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,nyv,s
     1p,idtask,nmt,ierr)
c multitasking particle push
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, ci, sx, sy, sz, sp
      dimension part(idimp,nop)
      dimension axy(3,nxv,nyv)
      dimension sp(3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RCR8PCM23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      sp(3,i) = 0.
      call MP_TASKSTART(idtask(i),RCR8PCM23,nargs,part(1,npo),axy,qbm,ci
     1,sp(1,i),sp(2,i),sp(3,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RCR8PCM23(part(1,npo),axy,qbm,ci,sx,sy,sz,idimp,npl,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      sx = sx + sp(1,i)
      sy = sy + sp(2,i)
      sz = sz + sp(3,i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRCR8PCM23L(part,axy,qbm,ci,sx,sy,sz,idimp,nop,nxv,nyv,
     1sp,idtask,nmt,ierr)
c multitasking particle push
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, ci,sx, sy, sz, sp
      dimension part(idimp,nop)
      dimension axy(3,nxv,nyv)
      dimension sp(3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RCR8PCM23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      sp(3,i) = 0.
      call MP_TASKSTART(idtask(i),RCR8PCM23L,nargs,part(1,npo),axy,qbm,c
     1i,sp(1,i),sp(2,i),sp(3,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RCR8PCM23L(part(1,npo),axy,qbm,ci,sx,sy,sz,idimp,npl,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      sx = sx + sp(1,i)
      sy = sy + sp(2,i)
      sz = sz + sp(3,i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRCR8PCM22(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv,sp,i
     1dtask,nmt,ierr)
c multitasking particle push
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, ci, sx, sy, sp
      dimension part(idimp,nop)
      dimension axy(2,nxv,nyv)
      dimension sp(2,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RCR8PCM22
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      call MP_TASKSTART(idtask(i),RCR8PCM22,nargs,part(1,npo),axy,qbm,ci
     1,sp(1,i),sp(2,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RCR8PCM22(part(1,npo),axy,qbm,ci,sx,sy,idimp,npl,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      sx = sx + sp(1,i)
      sy = sy + sp(2,i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRCR8PCM22L(part,axy,qbm,ci,sx,sy,idimp,nop,nxv,nyv,sp,
     1idtask,nmt,ierr)
c multitasking particle push
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, ci, sx, sy, sp
      dimension part(idimp,nop)
      dimension axy(2,nxv,nyv)
      dimension sp(2,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RCR8PCM22L
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      call MP_TASKSTART(idtask(i),RCR8PCM22L,nargs,part(1,npo),axy,qbm,c
     1i,sp(1,i),sp(2,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RCR8PCM22L(part(1,npo),axy,qbm,ci,sx,sy,idimp,npl,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      sx = sx + sp(1,i)
      sy = sy + sp(2,i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,nop,
     1nxv,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRHJPOST2
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 3
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRHJPOST2,nargs,part(1,npo),fxy,axy,da
     1xy,cup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRHJPOST2(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,npl,
     1nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,nop
     1,nxv,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRHJPOST2
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 3
      cup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRHJPOST2,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRHJPOST2(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,npl
     1,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 3
      cu(m,j) = cu(m,j) + cup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,nop
     1,nxv,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRHJPOST2L
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 3
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRHJPOST2L,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRHJPOST2L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,npl
     1,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,no
     1p,nxv,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRHJPOST2L
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 3
      cup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRHJPOST2L,nargs,part(1,npo),fxy,axy,
     1daxy,cup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRHJPOST2L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,np
     1l,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 3
      cu(m,j) = cu(m,j) + cup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,nop
     1,nxv,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRHJPOST22
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 2
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRHJPOST22,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRHJPOST22(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,npl
     1,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,no
     1p,nxv,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRHJPOST22
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 2
      cup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRHJPOST22,nargs,part(1,npo),fxy,axy,
     1daxy,cup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRHJPOST22(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,np
     1l,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      cu(m,j) = cu(m,j) + cup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,no
     1p,nxv,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRHJPOST22L
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 2
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRHJPOST22L,nargs,part(1,npo),fxy,axy,
     1daxy,cup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRHJPOST22L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,np
     1l,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,n
     1op,nxv,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, ci, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRHJPOST22L
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 2
      cup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRHJPOST22L,nargs,part(1,npo),fxy,axy
     1,daxy,cup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRHJPOST22L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,ci,idimp,n
     1pl,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      cu(m,j) = cu(m,j) + cup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,nx,
     1ny,nxv,nyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSH23
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRHPUSH23,nargs,part(1,npo),fxy,axy,da
     1xy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSH23(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,nx,
     1ny,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSH23(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,nx
     1,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSH23
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRHPUSH23,nargs,part(1,npo),fxy,axy,d
     1axy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSH23(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,nx
     1,ny,nxv,nxyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,nx
     1,ny,nxv,nyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSH23L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRHPUSH23L,nargs,part(1,npo),fxy,axy,d
     1axy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSH23L(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,nx
     1,ny,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSH23L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,n
     1x,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSH23L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRHPUSH23L,nargs,part(1,npo),fxy,axy,
     1daxy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSH23L(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,n
     1x,ny,nxv,nxyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,nx,
     1ny,nxv,nyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSH22
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRHPUSH22,nargs,part(1,npo),fxy,axy,da
     1xy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSH22(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,nx,
     1ny,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSH22(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,nx
     1,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSH22
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRHPUSH22,nargs,part(1,npo),fxy,axy,d
     1axy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSH22(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,nx
     1,ny,nxv,nxyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,nx
     1,ny,nxv,nyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSH22L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRHPUSH22L,nargs,part(1,npo),fxy,axy,d
     1axy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSH22L(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,nx
     1,ny,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSH22L(part,fxy,axy,daxy,qbm,dt,ci,ek,idimp,nop,n
     1x,ny,nxv,nxyv,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ci, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSH22L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRHPUSH22L,nargs,part(1,npo),fxy,axy,
     1daxy,qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSH22L(part(1,npo),fxy,axy,daxy,qbm,dt,ci,ek,idimp,npl,n
     1x,ny,nxv,nxyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nyv
     1,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHXH23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHXH23,nargs,part(1,npo),axy,qbm,
     1dt,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHXH23(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nyv
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHXH23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nx
     1yv,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHXH23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHXH23,nargs,part(1,npo),axy,qbm
     1,dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHXH23(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nx
     1yv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,ny
     1v,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHXH23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHXH23L,nargs,part(1,npo),axy,qbm
     1,dt,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHXH23L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,ny
     1v)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHXH23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n
     1xyv,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHXH23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHXH23L,nargs,part(1,npo),axy,qb
     1m,dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHXH23L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,n
     1xyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nyv
     1,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHXH22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHXH22,nargs,part(1,npo),axy,qbm,
     1dt,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHXH22(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nyv
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHXH22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nx
     1yv,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHXH22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHXH22,nargs,part(1,npo),axy,qbm
     1,dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHXH22(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nx
     1yv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,ny
     1v,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHXH22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHXH22L,nargs,part(1,npo),axy,qbm
     1,dt,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHXH22L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,ny
     1v)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHXH22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,n
     1xyv,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHXH22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHXH22L,nargs,part(1,npo),axy,qb
     1m,dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHXH22L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,n
     1xyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nyv,
     1idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHX23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHX23,nargs,part(1,npo),axy,qbm,d
     1t,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHX23(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHX23(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nxy
     1v,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHX23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHX23,nargs,part(1,npo),axy,qbm,
     1dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHX23(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nxy
     1v)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nyv
     1,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHX23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHX23L,nargs,part(1,npo),axy,qbm,
     1dt,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHX23L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nyv
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHX23L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nx
     1yv,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHX23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHX23L,nargs,part(1,npo),axy,qbm
     1,dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHX23L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nx
     1yv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nyv,
     1idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHX22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHX22,nargs,part(1,npo),axy,qbm,d
     1t,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHX22(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nyv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHX22(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nxy
     1v,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHX22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHX22,nargs,part(1,npo),axy,qbm,
     1dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHX22(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nxy
     1v)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nyv
     1,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRHPUSHX22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GRHPUSHX22L,nargs,part(1,npo),axy,qbm,
     1dt,ci,idimp,npp,nx,ny,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRHPUSHX22L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nyv
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRHPUSHX22L(part,axy,qbm,dt,ci,idimp,nop,nx,ny,nxv,nx
     1yv,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt, ci
      integer nop, idimp, nx, ny, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRHPUSHX22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSRHPUSHX22L,nargs,part(1,npo),axy,qbm
     1,dt,ci,idimp,npp,nx,ny,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRHPUSHX22L(part(1,npo),axy,qbm,dt,ci,idimp,npl,nx,ny,nxv,nx
     1yv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end