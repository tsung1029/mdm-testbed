c-----------------------------------------------------------------------
c 2d PIC multi-tasking library for pushing particles and depositing
c charge
c mpush2lib.f contains multi-tasking procedures to process particles:
c MGPOST2 multi-tasking wrapper for GPOST2
c MGSPOST2 multi-tasking wrapper for GSPOST2
c MGSPOST2X multi-tasking wrapper for GSPOST2X
c MGPOST2L multi-tasking wrapper for GPOST2L
c MGSPOST2L multi-tasking wrapper for GSPOST2L
c MGSPOST2XL multi-tasking wrapper for GSPOST2XL
c MGPOST2C multi-tasking wrapper for GPOST2C
c MGPUSH2 multi-tasking wrapper for GPUSH2
c MGSPUSH2 multi-tasking wrapper for GSPUSH2
c MGPUSH2L multi-tasking wrapper for GPUSH2L
c MGSPUSH2L multi-tasking wrapper for GSPUSH2L
c MGPUSH2C multi-tasking wrapper for GPUSH2C
c MSORTP2Y multi-tasking wrapper for SORTP2Y
c MSORTP2YL multi-tasking wrapper for SORTP2YL
c MDSORTP2Y multi-tasking wrapper for DSORTP2Y
c MDSORTP2YL multi-tasking wrapper for DSORTP2YL
c MPUSH2ZF multi-tasking wrapper for PUSH2ZF
c MDPOST2GL multi-tasking wrapper for DPOST2GL
c MDPOST2GLX multi-tasking wrapper for DPOST2GLX
c MPUSH2GL multi-tasking wrapper for PUSH2GL
c MPUSH2GLX multi-tasking wrapper for PUSH2GLX
c MGCJPOST2 multi-tasking wrapper for GCJPOST2
c MGCJPOST2L multi-tasking wrapper for GCJPOST2L
c MGCJPOST2C multi-tasking wrapper for GCJPOST2C
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: november 24, 2009
c-----------------------------------------------------------------------
      subroutine MDPOST2(part,q,qm,nop,idimp,nx,ny,nxv,qp,idtask,nmt,ier
     1r)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv,ny)
      dimension qp(nxv,ny,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external DPOST2
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, ny
      do 10 j = 1, nx
      qp(j,k,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),DPOST2,nargs,part(1,npo),qp(1,1,i),qm,
     1npp,idimp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call DPOST2(part(1,npo),q,qm,npl,idimp,nx,ny,nxv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, ny
      do 40 j = 1, nx
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPOST2(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv,nyv)
      dimension qp(nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external GPOST2
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, nyv
      do 10 j = 1, nxv
      qp(j,k,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GPOST2,nargs,part(1,npo),qp(1,1,i),qm,
     1npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GPOST2(part(1,npo),q,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, nyv
      do 40 j = 1, nxv
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST2(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt,ier
     1r)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxyv)
      dimension qp(nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST2
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxyv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST2,nargs,part(1,npo),qp(1,i),qm,n
     1pp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST2(part(1,npo),q,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxyv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MSPOST2X(part,q,qm,nop,idimp,nx,ny,nxv,nxvy,qp,idtask,n
     1mt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nx, ny, nxv, nxvy, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxvy)
      dimension qp(nxvy,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external SPOST2X
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxvy
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),SPOST2X,nargs,part(1,npo),qp(1,i),qm,n
     1pp,idimp,nx,ny,nxv,nxvy)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call SPOST2X(part(1,npo),q,qm,npl,idimp,nx,ny,nxv,nxvy)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxvy
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST2X(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt,ie
     1rr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxyv)
      dimension qp(nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST2X
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxyv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST2X,nargs,part(1,npo),qp(1,i),qm,
     1npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST2X(part(1,npo),q,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxyv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDPOST2L(part,q,qm,nop,idimp,nx,ny,nxv,qp,idtask,nmt,ie
     1rr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv,ny)
      dimension qp(nxv,ny,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external DPOST2L
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, ny
      do 10 j = 1, nx
      qp(j,k,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),DPOST2L,nargs,part(1,npo),qp(1,1,i),qm
     1,npp,idimp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call DPOST2L(part(1,npo),q,qm,npl,idimp,nx,ny,nxv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, ny
      do 40 j = 1, nx
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPOST2L(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,ierr
     1)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv,nyv)
      dimension qp(nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external GPOST2L
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, nyv
      do 10 j = 1, nxv
      qp(j,k,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GPOST2L,nargs,part(1,npo),qp(1,1,i),qm
     1,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GPOST2L(part(1,npo),q,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, nyv
      do 40 j = 1, nxv
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST2L(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt,ie
     1rr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxyv)
      dimension qp(nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST2L
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxyv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST2L,nargs,part(1,npo),qp(1,i),qm,
     1npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST2L(part(1,npo),q,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxyv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MSPOST2XL(part,q,qm,nop,idimp,nx,ny,nxv,nxvy,qp,idtask,
     1nmt,ierr)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nx, ny, nxv, nxvy, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxvy)
      dimension qp(nxvy,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external SPOST2XL
      data nargs /9/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxvy
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),SPOST2XL,nargs,part(1,npo),qp(1,i),qm,
     1npp,idimp,nx,ny,nxv,nxvy)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call SPOST2XL(part(1,npo),q,qm,npl,idimp,nx,ny,nxv,nxvy)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxvy
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPOST2XL(part,q,qm,nop,idimp,nxv,nxyv,qp,idtask,nmt,i
     1err)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxyv)
      dimension qp(nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSPOST2XL
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 10 j = 1, nxyv
      qp(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),GSPOST2XL,nargs,part(1,npo),qp(1,i),qm
     1,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GSPOST2XL(part(1,npo),q,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 30 j = 1, nxyv
      q(j) = q(j) + qp(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPOST2C(part,q,qm,nop,idimp,nxv,nyv,qp,idtask,nmt,ierr
     1)
c multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxv,nyv)
      dimension qp(nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external GPOST2C
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, nyv
      do 10 j = 1, nxv
      qp(j,k,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GPOST2C,nargs,part(1,npo),qp(1,1,i),qm
     1,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call GPOST2C(part(1,npo),q,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, nyv
      do 40 j = 1, nxv
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH2(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ekp,idt
     1ask,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external PUSH2
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PUSH2,nargs,part(1,npo),fx,fy,qbm,dt,e
     1kp(i),idimp,npp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call PUSH2(part(1,npo),fx,fy,qbm,dt,ek,idimp,npl,nx,ny,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipbc
     1,ekp,idtask,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GPUSH2
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GPUSH2,nargs,part(1,npo),fxy,qbm,dt,ek
     1p(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GPUSH2(part(1,npo),fxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,nyv,ipbc
     1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPUSH2(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,ip
     1bc,ekp,idtask,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSPUSH2
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSPUSH2,nargs,part(1,npo),fxy,qbm,dt,e
     1kp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSPUSH2(part(1,npo),fxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,nxyv,ip
     1bc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH2L(part,fx,fy,qbm,dt,ek,idimp,nop,nx,ny,nxv,ekp,id
     1task,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external PUSH2L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PUSH2L,nargs,part(1,npo),fx,fy,qbm,dt,
     1ekp(i),idimp,npp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call PUSH2L(part(1,npo),fx,fy,qbm,dt,ek,idimp,npl,nx,ny,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipb
     1c,ekp,idtask,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GPUSH2L
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GPUSH2L,nargs,part(1,npo),fxy,qbm,dt,e
     1kp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GPUSH2L(part(1,npo),fxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,nyv,ipb
     1c)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSPUSH2L(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nxyv,i
     1pbc,ekp,idtask,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSPUSH2L
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSPUSH2L,nargs,part(1,npo),fxy,qbm,dt,
     1ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSPUSH2L(part(1,npo),fxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,nxyv,i
     1pbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGPUSH2C(part,fxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,nyv,ipb
     1c,ekp,idtask,nmt,ierr)
c multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GPUSH2C
      data nargs /12/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GPUSH2C,nargs,part(1,npo),fxy,qbm,dt,e
     1kp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GPUSH2C(part(1,npo),fxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,nyv,ipb
     1c)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MSORTP2Y(part,pt,ip,npic,idimp,nop,ny1,npicp,idtask,nmt
     1,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npicp
      integer idimp, nop, ny1, idtask, nmt, ierr
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(ny1)
      dimension npicp(ny1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external SORTP2Y
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),SORTP2Y,nargs,part(1,npo),pt(npo),ip(n
     1po),npicp(1,i),idimp,npp,ny1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call SORTP2Y(part(1,npo),pt(npo),ip(npo),npic,idimp,npl,ny1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MSORTP2YL(part,pt,ip,npic,idimp,nop,ny1,npicp,idtask,nm
     1t,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npicp
      integer idimp, nop, ny1, idtask, nmt, ierr
      dimension part(idimp,nop), pt(nop)
      dimension ip(nop), npic(ny1)
      dimension npicp(ny1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external SORTP2YL
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),SORTP2YL,nargs,part(1,npo),pt(npo),ip(
     1npo),npicp(1,i),idimp,npp,ny1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call SORTP2YL(part(1,npo),pt(npo),ip(npo),npic,idimp,npl,ny1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDSORTP2Y(parta,partb,npic,idimp,nop,ny1,npicp,idtask,n
     1mt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npicp
      integer idimp, nop, ny1, idtask, nmt, ierr
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
      dimension npicp(ny1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external DSORTP2Y
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),DSORTP2Y,nargs,parta(1,npo),partb(1,np
     1o),npicp(1,i),idimp,npp,ny1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call DSORTP2Y(parta(1,npo),partb(1,npo),npic,idimp,npl,ny1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDSORTP2YL(parta,partb,npic,idimp,nop,ny1,npicp,idtask,
     1nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npicp
      integer idimp, nop, ny1, idtask, nmt, ierr
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
      dimension npicp(ny1,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external DSORTP2YL
      data nargs /6/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle sorting tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      call MP_TASKSTART(idtask(i),DSORTP2YL,nargs,parta(1,npo),partb(1,n
     1po),npicp(1,i),idimp,npp,ny1)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c sort remaining particles
      npo = npl + 1
      npl = nop - npl
      call DSORTP2YL(parta(1,npo),partb(1,npo),npic,idimp,npl,ny1)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH2ZF(part,dt,ek,idimp,nop,nx,ny,ipbc,ekp,idtask,nmt
     1,ierr)
c multitasking zero force particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, dt, ek, ekp
      integer idimp, nop, nx, ny, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external PUSH2ZF
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PUSH2ZF,nargs,part(1,npo),dt,ekp(i),id
     1imp,npp,nx,ny,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call PUSH2ZF(part(1,npo),dt,ek,idimp,npl,nx,ny,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDPOST2GL(part,q,sctx,qm,nop,idimp,nx,ny,nxvh,nyv,qp,sc
     1txp,idtask,nmt,ierr)
c multitasking gridless charge deposition
c qp = charge density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm
      complex q, sctx, qp, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxvh,nyv), sctx(nxvh)
      dimension qp(nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external DPOST2GL
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, nyv
      do 10 j = 1, nxvh
      qp(j,k,i) = cmplx(0.,0.)
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),DPOST2GL,nargs,part(1,npo),qp(1,1,i),s
     1ctxp(1,i),qm,npp,idimp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call DPOST2GL(part(1,npo),q,sctx,qm,npl,idimp,nx,ny,nxvh,nyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, nyv
      do 40 j = 1, nxvh
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDPOST2GLX(part,q,sctxp,qm,nop,idimp,nx,ny,nxvh,nyv,nnp
     1,qp,sctxpp,idtask,nmt,ierr)
c multitasking gridless charge deposition
c optimized version
c qp = charge density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm
      complex q, sctxp, qp, sctxpp
      integer nop, idimp, nx, ny, nxvh, nyv, nnp, idtask, nmt, ierr
      dimension part(idimp,nop), q(nxvh,nyv), sctxp(nxvh,nnp)
      dimension qp(nxvh,nyv,nmt), sctxpp(nxvh,nnp,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k
      external DPOST2GLX
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start charge deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear charge arrays
      do 20 k = 1, nyv
      do 10 j = 1, nxvh
      qp(j,k,i) = cmplx(0.,0.)
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),DPOST2GLX,nargs,part(1,npo),qp(1,1,i),
     1sctxpp(1,1,i),qm,npp,idimp,nx,ny,nxvh,nyv,nnp)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining charge
      npo = npl + 1
      npl = nop - npl
      call DPOST2GLX(part(1,npo),q,sctxp,qm,npl,idimp,nx,ny,nxvh,nyv,nnp
     1)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 50 k = 1, nyv
      do 40 j = 1, nxvh
      q(j,k) = q(j,k) + qp(j,k,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPUSH2GL(part,fxy,sctx,qbm,dt,ek,idimp,nop,nx,ny,nxvh,n
     1yv,ipbc,sctxp,ekp,idtask,nmt,ierr)
c multitasking gridless particle push
c sctxp = scratch arrays for sines and cosines
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, ek, ekp
      complex fxy, sctx, sctxp
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxvh,nyv), sctx(nxvh)
      dimension sctxp(nxvh,nmt), ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external PUSH2GL
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PUSH2GL,nargs,part(1,npo),fxy,sctxp(1,
     1i),qbm,dt,ekp(i),idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call PUSH2GL(part(1,npo),fxy,sctx,qbm,dt,ek,idimp,npl,nx,ny,nxvh,n
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
      subroutine MPUSH2GLX(part,fxy,sctxp,exyp,qbm,dt,ek,idimp,nop,nx,ny
     1,nxvh,nyv,nnp,ipbc,sctxpp,exypp,ekp,idtask,nmt,ierr)
c multitasking gridless particle push
c optimized version
c sctxpp = scratch arrays for sines and cosines
c exypp = scratch array for particle forces
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, ek, ekp
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
      external PUSH2GLX
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PUSH2GLX,nargs,part(1,npo),fxy,sctxpp(
     11,1,i),exypp(1,1,i),qbm,dt,ekp(i),idimp,npp,nx,ny,nxvh,nyv,nnp,ipb
     1c)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call PUSH2GLX(part(1,npo),fxy,sctxp,exyp,qbm,dt,ek,idimp,npl,nx,ny
     1,nxvh,nyv,nnp,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGCJPOST2(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cup,i
     1dtask,nmt,ierr)
c multitasking current density deposit
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, cup
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GCJPOST2
      data nargs /10/
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
      call MP_TASKSTART(idtask(i),GCJPOST2,nargs,part(1,npo),fxy,cup(1,1
     1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GCJPOST2(part(1,npo),fxy,cu,qm,qbm,dt,idimp,npl,nxv,nyv)
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
      subroutine MGCJPOST2L(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cup,
     1idtask,nmt,ierr)
c multitasking current density deposit
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, cup
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GCJPOST2L
      data nargs /10/
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
      call MP_TASKSTART(idtask(i),GCJPOST2L,nargs,part(1,npo),fxy,cup(1,
     11,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GCJPOST2L(part(1,npo),fxy,cu,qm,qbm,dt,idimp,npl,nxv,nyv)
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
      subroutine MGCJPOST2C(part,fxy,cu,qm,qbm,dt,idimp,nop,nxv,nyv,cup,
     1idtask,nmt,ierr)
c multitasking current density deposit
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, cu, qm, qbm, dt, cup
      integer idimp, nop, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GCJPOST2C
      data nargs /10/
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
      call MP_TASKSTART(idtask(i),GCJPOST2C,nargs,part(1,npo),fxy,cup(1,
     11,1,i),qm,qbm,dt,idimp,npp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GCJPOST2C(part(1,npo),fxy,cu,qm,qbm,dt,idimp,npl,nxv,nyv)
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
