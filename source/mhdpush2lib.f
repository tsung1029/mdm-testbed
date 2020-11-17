c 2d PIC multi-tasking library for pushing particles with darwin
c vector potential using hamiltonian formulation
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: october 29, 2007
c-----------------------------------------------------------------------
      subroutine MCR8PCM23(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv,sp,id
     1task,nmt,ierr)
c multitasking particle canonical momentum
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, sx, sy, sz, sp
      dimension part(idimp,nop)
      dimension axy(3,nxv,nyv)
      dimension sp(3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external CR8PCM23
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      sp(3,i) = 0.
      call MP_TASKSTART(idtask(i),CR8PCM23,nargs,part(1,npo),axy,qbm,sp(
     11,i),sp(2,i),sp(3,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call CR8PCM23(part(1,npo),axy,qbm,sx,sy,sz,idimp,npl,nxv,nyv)
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
      subroutine MCR8PCM23L(part,axy,qbm,sx,sy,sz,idimp,nop,nxv,nyv,sp,i
     1dtask,nmt,ierr)
c multitasking particle canonical momentum
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, sx, sy, sz, sp
      dimension part(idimp,nop)
      dimension axy(3,nxv,nyv)
      dimension sp(3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external CR8PCM23L
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      sp(3,i) = 0.
      call MP_TASKSTART(idtask(i),CR8PCM23L,nargs,part(1,npo),axy,qbm,sp
     1(1,i),sp(2,i),sp(3,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call CR8PCM23L(part(1,npo),axy,qbm,sx,sy,sz,idimp,npl,nxv,nyv)
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
      subroutine MCR8PCM22(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv,sp,idtas
     1k,nmt,ierr)
c multitasking particle canonical momentum
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, sx, sy, sp
      dimension part(idimp,nop)
      dimension axy(2,nxv,nyv)
      dimension sp(2,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external CR8PCM22
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      call MP_TASKSTART(idtask(i),CR8PCM22,nargs,part(1,npo),axy,qbm,sp(
     11,i),sp(2,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call CR8PCM22(part(1,npo),axy,qbm,sx,sy,idimp,npl,nxv,nyv)
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
      subroutine MCR8PCM22L(part,axy,qbm,sx,sy,idimp,nop,nxv,nyv,sp,idta
     1sk,nmt,ierr)
c multitasking particle canonical momentum
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      real part, axy, qbm, sx, sy, sp
      dimension part(idimp,nop)
      dimension axy(2,nxv,nyv)
      dimension sp(2,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external CR8PCM22L
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      call MP_TASKSTART(idtask(i),CR8PCM22L,nargs,part(1,npo),axy,qbm,sp
     1(1,i),sp(2,i),idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call CR8PCM22L(part(1,npo),axy,qbm,sx,sy,idimp,npl,nxv,nyv)
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
      subroutine MGHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv,
     1nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GHJPOST2
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GHJPOST2,nargs,part(1,npo),fxy,axy,dax
     1y,cup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GHJPOST2(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nxv,
     1nyv)
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
      subroutine MGSHJPOST2(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv
     1,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSHJPOST2
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GSHJPOST2,nargs,part(1,npo),fxy,axy,da
     1xy,cup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSHJPOST2(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nxv
     1,nxyv)
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
      subroutine MGHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv
     1,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GHJPOST2L
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GHJPOST2L,nargs,part(1,npo),fxy,axy,da
     1xy,cup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GHJPOST2L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nxv
     1,nyv)
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
      subroutine MGSHJPOST2L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nx
     1v,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSHJPOST2L
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GSHJPOST2L,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSHJPOST2L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nx
     1v,nxyv)
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
      subroutine MGHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nxv
     1,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GHJPOST22
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GHJPOST22,nargs,part(1,npo),fxy,axy,da
     1xy,cup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GHJPOST22(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nxv
     1,nyv)
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
      subroutine MGSHJPOST22(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nx
     1v,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSHJPOST22
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GSHJPOST22,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSHJPOST22(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nx
     1v,nxyv)
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
      subroutine MGHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,nx
     1v,nyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GHJPOST22L
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GHJPOST22L,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GHJPOST22L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,nx
     1v,nyv)
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
      subroutine MGSHJPOST22L(part,fxy,axy,daxy,cu,qm,qbm,dt,idimp,nop,n
     1xv,nxyv,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, cu, qm, qbm, dt, cup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSHJPOST22L
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),GSHJPOST22L,nargs,part(1,npo),fxy,axy,
     1daxy,cup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSHJPOST22L(part(1,npo),fxy,axy,daxy,cu,qm,qbm,dt,idimp,npl,n
     1xv,nxyv)
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
      subroutine MGHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSH23
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GHPUSH23,nargs,part(1,npo),fxy,axy,dax
     1y,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSH23(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny,n
     1xv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSH23(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSH23
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSHPUSH23,nargs,part(1,npo),fxy,axy,da
     1xy,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSH23(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny,
     1nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), axy(3,nxv,nyv), daxy(5,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSH23L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GHPUSH23L,nargs,part(1,npo),fxy,axy,da
     1xy,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSH23L(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny,
     1nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSH23L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny
     1,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), axy(3,nxyv), daxy(5,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSH23L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSHPUSH23L,nargs,part(1,npo),fxy,axy,d
     1axy,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSH23L(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny
     1,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSH22
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GHPUSH22,nargs,part(1,npo),fxy,axy,dax
     1y,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSH22(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny,n
     1xv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSH22(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSH22
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSHPUSH22,nargs,part(1,npo),fxy,axy,da
     1xy,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSH22(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny,
     1nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), axy(2,nxv,nyv), daxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSH22L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GHPUSH22L,nargs,part(1,npo),fxy,axy,da
     1xy,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSH22L(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny,
     1nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSH22L(part,fxy,axy,daxy,qbm,dt,ek,idimp,nop,nx,ny
     1,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking momentum push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, axy, daxy, qbm, dt, ek, ekp
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), axy(2,nxyv), daxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSH22L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSHPUSH22L,nargs,part(1,npo),fxy,axy,d
     1axy,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSH22L(part(1,npo),fxy,axy,daxy,qbm,dt,ek,idimp,npl,nx,ny
     1,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipb
     1c,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHXH23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHXH23,nargs,part(1,npo),axy,qbm,d
     1t,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHXH23(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ipb
     1c)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHXH23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHXH23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHXH23,nargs,part(1,npo),axy,qbm,
     1dt,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHXH23(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,i
     1pbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ip
     1bc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHXH23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHXH23L,nargs,part(1,npo),axy,qbm,
     1dt,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHXH23L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ip
     1bc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHXH23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,
     1ipbc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHXH23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHXH23L,nargs,part(1,npo),axy,qbm
     1,dt,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHXH23L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,
     1ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipb
     1c,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHXH22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHXH22,nargs,part(1,npo),axy,qbm,d
     1t,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHXH22(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ipb
     1c)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHXH22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHXH22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHXH22,nargs,part(1,npo),axy,qbm,
     1dt,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHXH22(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,i
     1pbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ip
     1bc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHXH22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHXH22L,nargs,part(1,npo),axy,qbm,
     1dt,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHXH22L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ip
     1bc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHXH22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,
     1ipbc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHXH22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHXH22L,nargs,part(1,npo),axy,qbm
     1,dt,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHXH22L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,
     1ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc
     1,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHX23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHX23,nargs,part(1,npo),axy,qbm,dt
     1,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHX23(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ipbc
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHX23(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHX23
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHX23,nargs,part(1,npo),axy,qbm,d
     1t,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHX23(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,ip
     1bc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipb
     1c,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHX23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHX23L,nargs,part(1,npo),axy,qbm,d
     1t,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHX23L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ipb
     1c)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHX23L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHX23L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHX23L,nargs,part(1,npo),axy,qbm,
     1dt,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHX23L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,i
     1pbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipbc
     1,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHX22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHX22,nargs,part(1,npo),axy,qbm,dt
     1,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHX22(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ipbc
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHX22(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHX22
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHX22,nargs,part(1,npo),axy,qbm,d
     1t,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHX22(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,ip
     1bc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nyv,ipb
     1c,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxv,nyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHX22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHX22L,nargs,part(1,npo),axy,qbm,d
     1t,idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHX22L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nyv,ipb
     1c)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSHPUSHX22L(part,axy,qbm,dt,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc,idtask,nmt,ierr)
c multitasking co-ordinate push
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, axy, qbm, dt
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(2,nxyv)
      dimension idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSHPUSHX22L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GSHPUSHX22L,nargs,part(1,npo),axy,qbm,
     1dt,idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSHPUSHX22L(part(1,npo),axy,qbm,dt,idimp,npl,nx,ny,nxv,nxyv,i
     1pbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MCR8PCM23GL(part,axy,sctx,qbm,sx,sy,sz,idimp,nop,nx,ny,
     1nxvh,nyv,sctxp,sp,idtask,nmt,ierr)
c multitasking gridless particle canonical momentum
c sctxp = scratch arrays for sines and cosines
c sp = field momentum arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer idimp, nop, nx, ny, nxvh, nyv, idtask, nmt, ierr
      real part, qbm, sx, sy, sz, sp
      complex axy, sctx, sctxp
      dimension part(idimp,nop), axy(3,nxvh,nyv)
      dimension sctx(nxvh)
      dimension sctxp(nxvh,nmt), sp(3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external CR8PCM23GL
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      sp(1,i) = 0.
      sp(2,i) = 0.
      sp(3,i) = 0.
      call MP_TASKSTART(idtask(i),CR8PCM23GL,nargs,part(1,npo),axy,sctxp
     1(1,i),qbm,sp(1,i),sp(2,i),sp(3,i),idimp,npp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call CR8PCM23GL(part(1,npo),axy,sctx,qbm,sx,sy,sz,idimp,npl,nx,ny,
     1nxvh,nyv)
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
      subroutine MGHJPOST2GL(part,fxy,axy,daxy,cu,sctx,qm,qbm,dt,idimp,n
     1op,nx,ny,nxvh,nyv,cup,sctxp,idtask,nmt,ierr)
c multitasking gridless current deposition
c cup = current density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm, qbm, dt
      complex fxy, axy, daxy, cu, sctx, cup, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), axy(3,nxvh,nyv), daxy(5,nxvh,nyv)
      dimension cu(3,nxvh,nyv), sctx(nxvh)
      dimension cup(3,nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GHJPOST2GL
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxvh
      do 10 m = 1, 3
      cup(m,j,k,i) = cmplx(0.,0.)
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GHJPOST2GL,nargs,part(1,npo),fxy,axy,d
     1axy,cup(1,1,1,i),sctxp(1,i),qm,qbm,dt,idimp,npp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GHJPOST2GL(part(1,npo),fxy,axy,daxy,cu,sctx,qm,qbm,dt,idimp,n
     1pl,nx,ny,nxvh,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxvh
      do 50 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSH23GL(part,fxy,axy,daxy,sctx,qbm,dt,ek,idimp,nop,
     1nx,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
c multitasking gridless momentum push
c sctxp = scratch arrays for sines and cosines
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, ek, ekp
      complex fxy, axy, daxy, sctx, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), axy(3,nxvh,nyv), daxy(5,nxvh,nyv)
      dimension sctx(nxvh)
      dimension sctxp(nxvh,nmt), ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSH23GL
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GHPUSH23GL,nargs,part(1,npo),fxy,axy,d
     1axy,sctxp(1,i),qbm,dt,ekp(i),idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSH23GL(part(1,npo),fxy,axy,daxy,sctx,qbm,dt,ek,idimp,npl,
     1nx,ny,nxvh,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHXH23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxvh
     1,nyv,ipbc,sctxp,idtask,nmt,ierr)
c multitasking gridless co-ordinate push
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt
      complex axy, sctx, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxvh,nyv)
      dimension sctx(nxvh)
      dimension sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHXH23GL
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHXH23GL,nargs,part(1,npo),axy,sct
     1xp(1,i),qbm,dt,idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHXH23GL(part(1,npo),axy,sctx,qbm,dt,idimp,npl,nx,ny,nxvh
     1,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGHPUSHX23GL(part,axy,sctx,qbm,dt,idimp,nop,nx,ny,nxvh,
     1nyv,ipbc,sctxp,idtask,nmt,ierr)
c multitasking gridless co-ordinate push
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt
      complex axy, sctx, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), axy(3,nxvh,nyv)
      dimension sctx(nxvh)
      dimension sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GHPUSHX23GL
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),GHPUSHX23GL,nargs,part(1,npo),axy,sctx
     1p(1,i),qbm,dt,idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GHPUSHX23GL(part(1,npo),axy,sctx,qbm,dt,idimp,npl,nx,ny,nxvh,
     1nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end

