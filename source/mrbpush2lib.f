c-----------------------------------------------------------------------
c 2d PIC multi-tasking library for pushing relatistic particles with
c magnetic field and depositing current
c mrbpush2lib.f contains multi-tasking procedures to process
c relativistic particles with and without magnetic fields:
c MGRJPOST2 multi-tasking wrapper for GRJPOST2
c MGSRJPOST2 multi-tasking wrapper for GSRJPOST2
c MGSRJPOST2X multi-tasking wrapper for GSRJPOST2X
c MGRJPOST2L multi-tasking wrapper for GRJPOST2L
c MGSRJPOST2L multi-tasking wrapper for GSRJPOST2L
c MGSRJPOST2XL multi-tasking wrapper for GSRJPOST2XL
c MGRJPOST22 multi-tasking wrapper for GRJPOST22
c MGSRJPOST22 multi-tasking wrapper for GSRJPOST22
c MGSRJPOST22X multi-tasking wrapper for GSRJPOST22X
c MGRJPOST22L multi-tasking wrapper for GRJPOST22L
c MGSRJPOST22L multi-tasking wrapper for GSRJPOST22L
c MGSRJPOST22XL multi-tasking wrapper for GSRJPOST22XL
c MGRJPOST2C multi-tasking wrapper for GRJPOST2C
c MGRJPOST22C multi-tasking wrapper for GRJPOST22C
c MGRPUSH2 multi-tasking wrapper for GRPUSH2
c MGSRPUSH2 multi-tasking wrapper for GSRPUSH2
c MGRPUSH2L multi-tasking wrapper for GRPUSH2L
c MGSRPUSH2L multi-tasking wrapper for GSRPUSH2L
c MGRBPUSH2 multi-tasking wrapper for GRBPUSH2
c MGSRBPUSH2 multi-tasking wrapper for GSRBPUSH2
c MGRBPUSH2L multi-tasking wrapper for GRBPUSH2L
c MGSRBPUSH2L multi-tasking wrapper for GSRBPUSH2L
c MGRBPUSH23 multi-tasking wrapper for GRBPUSH23
c MGSRBPUSH23 multi-tasking wrapper for GSRBPUSH23
c MGRBPUSH23L multi-tasking wrapper for GRBPUSH23L
c MGSRBPUSH23L multi-tasking wrapper for GSRBPUSH23L
c MGRBPUSH22 multi-tasking wrapper for GRBPUSH22
c MGSRBPUSH22 multi-tasking wrapper for GSRBPUSH22
c MGRBPUSH22L multi-tasking wrapper for GRBPUSH22L
c MGSRBPUSH22L multi-tasking wrapper for GSRBPUSH22L
c MGRPUSH2C multi-tasking wrapper for GRPUSH2C
c MGRBPUSH2C multi-tasking wrapper for GRBPUSH2C
c MGRBPUSH23C multi-tasking wrapper for GRBPUSH23C
c MGRBPUSH22C multi-tasking wrapper for GRBPUSH22C
c MRPUSH2ZF multi-tasking wrapper for RPUSH2ZF
c MRPUSH23ZF multi-tasking wrapper for RPUSH23ZF
c MRDJPOST2GL multi-tasking wrapper for RDJPOST2GL
c MRPUSH2GL multi-tasking wrapper for RPUSH2GL
c MRPUSH2GLX multi-tasking wrapper for RPUSH2GLX
c MRBPUSH23GL multi-tasking wrapper for RBPUSH23GL
c MGRCJPOST2 multi-tasking wrapper for GRCJPOST2
c MGRCJPOST2L multi-tasking wrapper for GRCJPOST2L
c MGRCJPOST2C multi-tasking wrapper for GRCJPOST2C
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: november 25, 2009
c-----------------------------------------------------------------------
      subroutine MGRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc
     1,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRJPOST2
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
      call MP_TASKSTART(idtask(i),GRJPOST2,nargs,part(1,npo),cup(1,1,1,i
     1),qm,dt,ci,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRJPOST2(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nyv,ipbc
     1)
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
      subroutine MGSRJPOST2(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,ip
     1bc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRJPOST2
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
      call MP_TASKSTART(idtask(i),GSRJPOST2,nargs,part(1,npo),cup(1,1,i)
     1,qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST2(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,ip
     1bc)
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
      subroutine MGSRJPOST2X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,i
     1pbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3*nxyv)
      dimension cup(3*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSRJPOST2X
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 10 j = 1, 3*nxyv
      cup(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSRJPOST2X,nargs,part(1,npo),cup(1,i),
     1qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST2X(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,i
     1pbc)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 30 j = 1, 3*nxyv
      cu(j) = cu(j) + cup(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipb
     1c,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRJPOST2L
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
      call MP_TASKSTART(idtask(i),GRJPOST2L,nargs,part(1,npo),cup(1,1,1,
     1i),qm,dt,ci,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRJPOST2L(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nyv,ipb
     1c)
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
      subroutine MGSRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,i
     1pbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRJPOST2L
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
      call MP_TASKSTART(idtask(i),GSRJPOST2L,nargs,part(1,npo),cup(1,1,i
     1),qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST2L(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,i
     1pbc)
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
      subroutine MGSRJPOST2XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,
     1ipbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3*nxyv)
      dimension cup(3*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSRJPOST2XL
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 10 j = 1, 3*nxyv
      cup(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSRJPOST2XL,nargs,part(1,npo),cup(1,i)
     1,qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST2XL(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,
     1ipbc)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 30 j = 1, 3*nxyv
      cu(j) = cu(j) + cup(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipb
     1c,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRJPOST22
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
      call MP_TASKSTART(idtask(i),GRJPOST22,nargs,part(1,npo),cup(1,1,1,
     1i),qm,dt,ci,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRJPOST22(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nyv,ipb
     1c)
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
      subroutine MGSRJPOST22(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,i
     1pbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRJPOST22
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
      call MP_TASKSTART(idtask(i),GSRJPOST22,nargs,part(1,npo),cup(1,1,i
     1),qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST22(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,i
     1pbc)
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
      subroutine MGSRJPOST22X(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,
     1ipbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2*nxyv)
      dimension cup(2*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSRJPOST22X
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 10 j = 1, 2*nxyv
      cup(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSRJPOST22X,nargs,part(1,npo),cup(1,i)
     1,qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST22X(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,
     1ipbc)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 30 j = 1, 2*nxyv
      cu(j) = cu(j) + cup(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ip
     1bc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRJPOST22L
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
      call MP_TASKSTART(idtask(i),GRJPOST22L,nargs,part(1,npo),cup(1,1,1
     1,i),qm,dt,ci,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRJPOST22L(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nyv,ip
     1bc)
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
      subroutine MGSRJPOST22L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv,
     1ipbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRJPOST22L
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
      call MP_TASKSTART(idtask(i),GSRJPOST22L,nargs,part(1,npo),cup(1,1,
     1i),qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST22L(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv,
     1ipbc)
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
      subroutine MGSRJPOST22XL(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nxyv
     1,ipbc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2*nxyv)
      dimension cup(2*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSRJPOST22XL
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 10 j = 1, 2*nxyv
      cup(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSRJPOST22XL,nargs,part(1,npo),cup(1,i
     1),qm,dt,ci,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRJPOST22XL(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nxyv
     1,ipbc)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 30 j = 1, 2*nxyv
      cu(j) = cu(j) + cup(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRJPOST2C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipb
     1c,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRJPOST2C
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
      call MP_TASKSTART(idtask(i),GRJPOST2C,nargs,part(1,npo),cup(1,1,1,
     1i),qm,dt,ci,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRJPOST2C(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nyv,ipb
     1c)
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
      subroutine MGRJPOST22C(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ip
     1bc,cup,idtask,nmt,ierr)
c multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRJPOST22C
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
      call MP_TASKSTART(idtask(i),GRJPOST22C,nargs,part(1,npo),cup(1,1,1
     1,i),qm,dt,ci,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRJPOST22C(part(1,npo),cu,qm,dt,ci,npl,idimp,nx,ny,nxv,nyv,ip
     1bc)
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
      subroutine MGRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv,
     1ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRPUSH2
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRPUSH2,nargs,part(1,npo),fxy,qbm,dt,c
     1i,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRPUSH2(part(1,npo),fxy,qbm,dt,ci,ek,idimp,npl,nx,ny,nxv,nyv,
     1ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRPUSH2(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nxy
     1v,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRPUSH2
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRPUSH2,nargs,part(1,npo),fxy,qbm,dt,
     1ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRPUSH2(part(1,npo),fxy,qbm,dt,ci,ek,idimp,npl,nx,ny,nxv,nxy
     1v,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv
     1,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRPUSH2L
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRPUSH2L,nargs,part(1,npo),fxy,qbm,dt,
     1ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRPUSH2L(part(1,npo),fxy,qbm,dt,ci,ek,idimp,npl,nx,ny,nxv,nyv
     1,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRPUSH2L(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nx
     1yv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRPUSH2L
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRPUSH2L,nargs,part(1,npo),fxy,qbm,dt
     1,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRPUSH2L(part(1,npo),fxy,qbm,dt,ci,ek,idimp,npl,nx,ny,nxv,nx
     1yv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH2
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH2,nargs,part(1,npo),fxy,bxy,qbm
     1,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH2(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,ny
     1,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRBPUSH2(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRBPUSH2
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRBPUSH2,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRBPUSH2(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH2L
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH2L,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH2L(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,
     1ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRBPUSH2L
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRBPUSH2L,nargs,part(1,npo),fxy,bxy,q
     1bm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRBPUSH2L(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,
     1ny,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH23
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH23,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH23(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRBPUSH23(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,
     1ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRBPUSH23
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRBPUSH23,nargs,part(1,npo),fxy,bxy,q
     1bm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRBPUSH23(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,
     1ny,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,
     1ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH23L
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH23L,nargs,part(1,npo),fxy,bxy,q
     1bm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH23L(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,
     1ny,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx
     1,ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRBPUSH23L
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRBPUSH23L,nargs,part(1,npo),fxy,bxy,
     1qbm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRBPUSH23L(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx
     1,ny,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,ny
     1,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH22
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH22,nargs,part(1,npo),fxy,bz,qbm
     1,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH22(part(1,npo),fxy,bz,qbm,dt,dtc,ci,ek,idimp,npl,nx,ny
     1,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRBPUSH22(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRBPUSH22
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRBPUSH22,nargs,part(1,npo),fxy,bz,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRBPUSH22(part(1,npo),fxy,bz,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH22L
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH22L,nargs,part(1,npo),fxy,bz,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH22L(part(1,npo),fxy,bz,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRBPUSH22L(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,
     1ny,nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSRBPUSH22L
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSRBPUSH22L,nargs,part(1,npo),fxy,bz,q
     1bm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSRBPUSH22L(part(1,npo),fxy,bz,qbm,dt,dtc,ci,ek,idimp,npl,nx,
     1ny,nxv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRPUSH2C(part,fxy,qbm,dt,ci,ek,idimp,nop,nx,ny,nxv,nyv
     1,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRPUSH2C
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRPUSH2C,nargs,part(1,npo),fxy,qbm,dt,
     1ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRPUSH2C(part(1,npo),fxy,qbm,dt,ci,ek,idimp,npl,nx,ny,nxv,nyv
     1,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH2C
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH2C,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH2C(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,
     1ny,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH23C
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH23C,nargs,part(1,npo),fxy,bxy,q
     1bm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH23C(part(1,npo),fxy,bxy,qbm,dt,dtc,ci,ek,idimp,npl,nx,
     1ny,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRBPUSH22C(part,fxy,bz,qbm,dt,dtc,ci,ek,idimp,nop,nx,n
     1y,nxv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking relativistic magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ci, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GRBPUSH22C
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GRBPUSH22C,nargs,part(1,npo),fxy,bz,qb
     1m,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GRBPUSH22C(part(1,npo),fxy,bz,qbm,dt,dtc,ci,ek,idimp,npl,nx,n
     1y,nxv,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRPUSH2ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,idtask
     1,nmt,ierr)
c multitasking 2d zero force relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RPUSH2ZF
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),RPUSH2ZF,nargs,part(1,npo),dt,ci,ekp(i
     1),idimp,npp,nx,ny,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RPUSH2ZF(part(1,npo),dt,ci,ek,idimp,npl,nx,ny,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRPUSH23ZF(part,dt,ci,ek,idimp,nop,nx,ny,ipbc,ekp,idtas
     1k,nmt,ierr)
c multitasking 2-1/2d zero force relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, dt, ci, ek, ekp
      integer idimp, nop, nx, ny, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RPUSH23ZF
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),RPUSH23ZF,nargs,part(1,npo),dt,ci,ekp(
     1i),idimp,npp,nx,ny,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RPUSH23ZF(part(1,npo),dt,ci,ek,idimp,npl,nx,ny,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRDJPOST2GL(part,cu,sctx,qm,dt,ci,nop,idimp,nx,ny,nxvh,
     1nyv,ipbc,cup,sctxp,idtask,nmt,ierr)
c multitasking gridless current deposition for relativistic particles
c cup = current density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm, dt, ci
      complex cu, sctx, cup, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxvh,nyv), sctx(nxvh)
      dimension cup(3,nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external RDJPOST2GL
      data nargs /13/
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
      call MP_TASKSTART(idtask(i),RDJPOST2GL,nargs,part(1,npo),cup(1,1,1
     1,i),sctxp(1,i),qm,dt,ci,npp,idimp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call RDJPOST2GL(part(1,npo),cu,sctx,qm,dt,ci,npl,idimp,nx,ny,nxvh,
     1nyv,ipbc)
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
      subroutine MRPUSH2GL(part,fxy,sctx,qbm,dt,ci,ek,idimp,nop,nx,ny,nx
     1vh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
c multitasking gridless relativistic particle push
c sctxp = scratch arrays for sines and cosines
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, ci, ek, ekp
      complex fxy, sctx, sctxp
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxvh,nyv), sctx(nxvh)
      dimension sctxp(nxvh,nmt), ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RPUSH2GL
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),RPUSH2GL,nargs,part(1,npo),fxy,sctxp(1
     1,i),qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RPUSH2GL(part(1,npo),fxy,sctx,qbm,dt,ci,ek,idimp,npl,nx,ny,nx
     1vh,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,nop,n
     1x,ny,nxvh,nyv,nnp,ipbc,sctxpp,exypp,ekp,idtask,nmt,ierr)
c multitasking gridless relativistic particle push
c optimized version
c sctxpp = scratch arrays for sines and cosines
c exypp = scratch array for particle forces
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, ci, ek, ekp
      double precision exyp, exypp
      complex fxy, sctxp, sctxpp
      integer idimp, nop, nx, ny, nxvh, nyv, nnp, ipbc, idtask, nmt
      integer ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxvh,nyv), sctxp(nxvh,nnp), exyp(2,nnp)
      dimension sctxpp(nxvh,nnp,nmt), ekp(nmt), idtask(nmt)
      dimension exypp(2,nnp,nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RPUSH2GLX
      data nargs /16/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),RPUSH2GLX,nargs,part(1,npo),fxy,sctxpp
     1(1,1,i),exypp(1,1,i),qbm,dt,ci,ekp(i),idimp,npp,nx,ny,nxvh,nyv,nnp
     2,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RPUSH2GLX(part(1,npo),fxy,sctxp,exyp,qbm,dt,ci,ek,idimp,npl,n
     1x,ny,nxvh,nyv,nnp,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MRBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,no
     1p,nx,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
c multitasking gridless relativistic magnetized particle push
c sctxp = scratch arrays for sines and cosines
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, dtc, ci, ek, ekp
      complex fxy, bxy, sctx, sctxp
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), bxy(3,nxvh,nyv), sctx(nxvh)
      dimension sctxp(nxvh,nmt), ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external RBPUSH23GL
      data nargs /16/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),RBPUSH23GL,nargs,part(1,npo),fxy,bxy,s
     1ctxp(1,i),qbm,dt,dtc,ci,ekp(i),idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call RBPUSH23GL(part(1,npo),fxy,bxy,sctx,qbm,dt,dtc,ci,ek,idimp,np
     1l,nx,ny,nxvh,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRCJPOST2(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv,c
     1up,idtask,nmt,ierr)
c multitasking current density deposit for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, ci, cup
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRCJPOST2
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GRCJPOST2,nargs,part(1,npo),fxy,cup(1,
     11,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRCJPOST2(part(1,npo),fxy,cu,qm,qbm,dt,ci,idimp,npl,nxv,nyv)
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
      subroutine MGRCJPOST2L(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv,
     1cup,idtask,nmt,ierr)
c multitasking current density deposit for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, ci, cup
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRCJPOST2L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GRCJPOST2L,nargs,part(1,npo),fxy,cup(1
     1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRCJPOST2L(part(1,npo),fxy,cu,qm,qbm,dt,ci,idimp,npl,nxv,nyv)
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
      subroutine MGRCJPOST2C(part,fxy,cu,qm,qbm,dt,ci,idimp,nop,nxv,nyv,
     1cup,idtask,nmt,ierr)
c multitasking current density deposit for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, ci, cup
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRCJPOST2C
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GRCJPOST2C,nargs,part(1,npo),fxy,cup(1
     1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRCJPOST2C(part(1,npo),fxy,cu,qm,qbm,dt,ci,idimp,npl,nxv,nyv)
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
