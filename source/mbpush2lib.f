c-----------------------------------------------------------------------
c 2d PIC multi-tasking library for pushing particles with magnetic field
c and depositing current
c mbpush2lib.f contains multi-tasking procedures to process particles
c              with magnetic fields:
c MGJPOST2 multi-tasking wrapper for GJPOST2
c MGSJPOST2 multi-tasking wrapper for GSJPOST2
c MGSJPOST2X multi-tasking wrapper for GSJPOST2X
c MGJPOST2L multi-tasking wrapper for GJPOST2L
c MGSJPOST2L multi-tasking wrapper for GSJPOST2
c MGSJPOST2XL multi-tasking wrapper for GSJPOST2XL
c MGJPOST22 multi-tasking wrapper for GJPOST22
c MGSJPOST22 multi-tasking wrapper for GSJPOST22
c MGSJPOST22X multi-tasking wrapper for GSJPOST22X
c MGJPOST22L multi-tasking wrapper for GJPOST22L
c MGSJPOST22L multi-tasking wrapper for GSJPOST22L
c MGSJPOST22XL multi-tasking wrapper for GSJPOST22XL
c MGJPOST2C multi-tasking wrapper for GJPOST2C
c MGJPOST22C multi-tasking wrapper for GJPOST22C
c MGBPUSH2 multi-tasking wrapper for GBPUSH2
c MGSBPUSH2 multi-tasking wrapper for GSBPUSH2
c MGBPUSH2L multi-tasking wrapper for GBPUSH2L
c MGSBPUSH2L multi-tasking wrapper for GSBPUSH2L
c MGBPUSH23 multi-tasking wrapper for GBPUSH23
c MGSBPUSH23 multi-tasking wrapper for GSBPUSH23
c MGBPUSH23L multi-tasking wrapper for GBPUSH23L
c MGSBPUSH23L multi-tasking wrapper for GSBPUSH23L
c MGBPUSH22 multi-tasking wrapper for GBPUSH22
c MGSBPUSH22 multi-tasking wrapper for GSBPUSH22
c MGBPUSH22L multi-tasking wrapper for GBPUSH22L
c MGSBPUSH22L multi-tasking wrapper for GSBPUSH22L
c MGBPUSH2C multi-tasking wrapper for GBPUSH2C
c MGBPUSH23C multi-tasking wrapper for GBPUSH23C
c MGBPUSH22C multi-tasking wrapper for GBPUSH22C
c MDJPOST2GL multi-tasking wrapper for DJPOST2GL
c MBPUSH23GL multi-tasking wrapper for BPUSH23GL
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: november 25, 2009
c-----------------------------------------------------------------------
      subroutine MDJPOST2(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv,cup
     1,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cux, cuy, cuz, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), cux(nxv,ny), cuy(nxv,ny), cuz(nxv,ny)
      dimension cup(nxv,ny,3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external DJPOST2
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 m = 1, 3
      do 20 k = 1, ny
      do 10 j = 1, nx
      cup(j,k,m,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),DJPOST2,nargs,part(1,npo),cup(1,1,1,i)
     1,cup(1,1,2,i),cup(1,1,3,i),qm,dt,npp,idimp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call DJPOST2(part(1,npo),cux,cuy,cuz,qm,dt,npl,idimp,nx,ny,nxv)
c wait for tasks to complete
      do 70 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 60 k = 1, ny
      do 50 j = 1, nx
      cux(j,k) = cux(j,k) + cup(j,k,1,i)
      cuy(j,k) = cuy(j,k) + cup(j,k,2,i)
      cuz(j,k) = cuz(j,k) + cup(j,k,3,i)
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,cup
     1,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GJPOST2
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
      do 10 m = 1, 3
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GJPOST2,nargs,part(1,npo),cup(1,1,1,i)
     1,qm,dt,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GJPOST2(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nyv,ipbc)
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
      subroutine MGSJPOST2(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc,c
     1up,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSJPOST2
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST2,nargs,part(1,npo),cup(1,1,i),
     1qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST2(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc)
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
      subroutine MSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy,cup,id
     1task,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxvy, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3*nxvy)
      dimension cup(3*nxvy,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external SJPOST2X
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 10 j = 1, 3*nxvy
      cup(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),SJPOST2X,nargs,part(1,npo),cup(1,i),qm
     1,dt,npp,idimp,nx,ny,nxv,nxvy)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call SJPOST2X(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxvy)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 30 j = 1, 3*nxvy
      cu(j) = cu(j) + cup(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSJPOST2X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc,
     1cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3*nxyv)
      dimension cup(3*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSJPOST2X
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST2X,nargs,part(1,npo),cup(1,i),q
     1m,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST2X(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc)
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
      subroutine MDJPOST2L(part,cux,cuy,cuz,qm,dt,nop,idimp,nx,ny,nxv,cu
     1p,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cux, cuy, cuz, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop), cux(nxv,ny), cuy(nxv,ny), cuz(nxv,ny)
      dimension cup(nxv,ny,3,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external DJPOST2L
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 m = 1, 3
      do 20 k = 1, ny
      do 10 j = 1, nx
      cup(j,k,m,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),DJPOST2L,nargs,part(1,npo),cup(1,1,1,i
     1),cup(1,1,2,i),cup(1,1,3,i),qm,dt,npp,idimp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call DJPOST2L(part(1,npo),cux,cuy,cuz,qm,dt,npl,idimp,nx,ny,nxv)
c wait for tasks to complete
      do 70 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 60 k = 1, ny
      do 50 j = 1, nx
      cux(j,k) = cux(j,k) + cup(j,k,1,i)
      cuy(j,k) = cuy(j,k) + cup(j,k,2,i)
      cuz(j,k) = cuz(j,k) + cup(j,k,3,i)
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,cu
     1p,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GJPOST2L
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
      do 10 m = 1, 3
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GJPOST2L,nargs,part(1,npo),cup(1,1,1,i
     1),qm,dt,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GJPOST2L(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nyv,ipbc)
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
      subroutine MGSJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc,
     1cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxyv)
      dimension cup(3,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSJPOST2L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST2L,nargs,part(1,npo),cup(1,1,i)
     1,qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST2L(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc)
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
      subroutine MSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxvy,cup,i
     1dtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxvy, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3*nxvy)
      dimension cup(3*nxvy,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external SJPOST2XL
      data nargs /10/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 20 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 10 j = 1, 3*nxvy
      cup(j,i) = 0.
   10 continue
      call MP_TASKSTART(idtask(i),SJPOST2XL,nargs,part(1,npo),cup(1,i),q
     1m,dt,npp,idimp,nx,ny,nxv,nxvy)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call SJPOST2XL(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxvy)
c wait for tasks to complete
      do 40 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 30 j = 1, 3*nxvy
      cu(j) = cu(j) + cup(j,i)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSJPOST2XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc
     1,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3*nxyv)
      dimension cup(3*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSJPOST2XL
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST2XL,nargs,part(1,npo),cup(1,i),
     1qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST2XL(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc
     1)
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
      subroutine MGJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,cu
     1p,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GJPOST22
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
      call MP_TASKSTART(idtask(i),GJPOST22,nargs,part(1,npo),cup(1,1,1,i
     1),qm,dt,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GJPOST22(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nyv,ipbc)
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
      subroutine MGSJPOST22(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc,
     1cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSJPOST22
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST22,nargs,part(1,npo),cup(1,1,i)
     1,qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST22(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc)
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
      subroutine MGSJPOST22X(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc
     1,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2*nxyv)
      dimension cup(2*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSJPOST22X
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST22X,nargs,part(1,npo),cup(1,i),
     1qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST22X(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc
     1)
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
      subroutine MGJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,c
     1up,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GJPOST22L
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
      call MP_TASKSTART(idtask(i),GJPOST22L,nargs,part(1,npo),cup(1,1,1,
     1i),qm,dt,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GJPOST22L(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nyv,ipbc)
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
      subroutine MGSJPOST22L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipbc
     1,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxyv)
      dimension cup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSJPOST22L
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST22L,nargs,part(1,npo),cup(1,1,i
     1),qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST22L(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipbc
     1)
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
      subroutine MGSJPOST22XL(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nxyv,ipb
     1c,cup,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2*nxyv)
      dimension cup(2*nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j
      external GSJPOST22XL
      data nargs /11/
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
      call MP_TASKSTART(idtask(i),GSJPOST22XL,nargs,part(1,npo),cup(1,i)
     1,qm,dt,npp,idimp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   20 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSJPOST22XL(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nxyv,ipb
     1c)
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
      subroutine MGJPOST2C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,cu
     1p,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GJPOST2C
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
      do 10 m = 1, 3
      cup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GJPOST2C,nargs,part(1,npo),cup(1,1,1,i
     1),qm,dt,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GJPOST2C(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nyv,ipbc)
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
      subroutine MGJPOST22C(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc,c
     1up,idtask,nmt,ierr)
c multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer nop, idimp, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GJPOST22C
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
      call MP_TASKSTART(idtask(i),GJPOST22C,nargs,part(1,npo),cup(1,1,1,
     1i),qm,dt,npp,idimp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GJPOST22C(part(1,npo),cu,qm,dt,npl,idimp,nx,ny,nxv,nyv,ipbc)
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
      subroutine MBPUSH2(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny,n
     1xv,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external BPUSH2
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),BPUSH2,nargs,part(1,npo),fx,fy,bx,by,b
     1z,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call BPUSH2(part(1,npo),fx,fy,bx,by,bz,qbm,dt,ek,idimp,npl,nx,ny,n
     1xv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH2
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH2,nargs,part(1,npo),fxy,bxy,qbm,
     1dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH2(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,nxv
     1,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSBPUSH2(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSBPUSH2
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSBPUSH2,nargs,part(1,npo),fxy,bxy,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSBPUSH2(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MBPUSH2L(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny,
     1nxv,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external BPUSH2L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),BPUSH2L,nargs,part(1,npo),fx,fy,bx,by,
     1bz,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call BPUSH2L(part(1,npo),fx,fy,bx,by,bz,qbm,dt,ek,idimp,npl,nx,ny,
     1nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH2L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH2L,nargs,part(1,npo),fxy,bxy,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH2L(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSBPUSH2L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,n
     1xv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSBPUSH2L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSBPUSH2L,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSBPUSH2L(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,n
     1xv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MBPUSH2CQ(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny
     1,nxv,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external BPUSH2CQ
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),BPUSH2CQ,nargs,part(1,npo),fx,fy,bx,by
     1,bz,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call BPUSH2CQ(part(1,npo),fx,fy,bx,by,bz,qbm,dt,ek,idimp,npl,nx,ny
     1,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH2CQ(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,n
     1yv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH2CQ
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH2CQ,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH2CQ(part(1,npo),fxy,bxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,n
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
      subroutine MBPUSH2CL(part,fx,fy,bx,by,bz,qbm,dt,ek,idimp,nop,nx,ny
     1,nxv,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, bx, by, bz, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fx(nxv,ny), fy(nxv,ny)
      dimension bx(nxv,ny), by(nxv,ny), bz(nxv,ny)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external BPUSH2CL
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),BPUSH2CL,nargs,part(1,npo),fx,fy,bx,by
     1,bz,qbm,dt,ekp(i),idimp,npp,nx,ny,nxv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call BPUSH2CL(part(1,npo),fx,fy,bx,by,bz,qbm,dt,ek,idimp,npl,nx,ny
     1,nxv)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH2CL(part,fxy,bxy,qbm,dt,ek,idimp,nop,nx,ny,nxv,n
     1yv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH2CL
      data nargs /13/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH2CL,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH2CL(part(1,npo),fxy,bxy,qbm,dt,ek,idimp,npl,nx,ny,nxv,n
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
      subroutine MGBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH23
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH23,nargs,part(1,npo),fxy,bxy,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH23(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSBPUSH23(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,n
     1xv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSBPUSH23
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSBPUSH23,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSBPUSH23(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,n
     1xv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,n
     1xv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH23L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH23L,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH23L(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,n
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
      subroutine MGSBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,
     1nxv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxyv), bxy(3,nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSBPUSH23L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSBPUSH23L,nargs,part(1,npo),fxy,bxy,q
     1bm,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSBPUSH23L(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,
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
      subroutine MGBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nxv
     1,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH22
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH22,nargs,part(1,npo),fxy,bz,qbm,
     1dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH22(part(1,npo),fxy,bz,qbm,dt,dtc,ek,idimp,npl,nx,ny,nxv
     1,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSBPUSH22(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSBPUSH22
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSBPUSH22,nargs,part(1,npo),fxy,bz,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSBPUSH22(part(1,npo),fxy,bz,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH22L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH22L,nargs,part(1,npo),fxy,bz,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH22L(part(1,npo),fxy,bz,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSBPUSH22L(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,n
     1xv,nxyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nxyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxyv), bz(nxyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GSBPUSH22L
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GSBPUSH22L,nargs,part(1,npo),fxy,bz,qb
     1m,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nxyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GSBPUSH22L(part(1,npo),fxy,bz,qbm,dt,dtc,ek,idimp,npl,nx,ny,n
     1xv,nxyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGBPUSH2C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH2C
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH2C,nargs,part(1,npo),fxy,bxy,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH2C(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
      subroutine MGBPUSH23C(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,n
     1xv,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH23C
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH23C,nargs,part(1,npo),fxy,bxy,qb
     1m,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH23C(part(1,npo),fxy,bxy,qbm,dt,dtc,ek,idimp,npl,nx,ny,n
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
      subroutine MGBPUSH22C(part,fxy,bz,qbm,dt,dtc,ek,idimp,nop,nx,ny,nx
     1v,nyv,ipbc,ekp,idtask,nmt,ierr)
c multitasking magnetized particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, qbm, dt, dtc, ek, ekp
      integer idimp, nop, nx, ny, nxv, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external GBPUSH22C
      data nargs /14/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),GBPUSH22C,nargs,part(1,npo),fxy,bz,qbm
     1,dt,dtc,ekp(i),idimp,npp,nx,ny,nxv,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call GBPUSH22C(part(1,npo),fxy,bz,qbm,dt,dtc,ek,idimp,npl,nx,ny,nx
     1v,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MDJPOST2GL(part,cu,sctx,qm,dt,nop,idimp,nx,ny,nxvh,nyv,
     1ipbc,cup,sctxp,idtask,nmt,ierr)
c multitasking gridless current deposition
c cup = current density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm, dt
      complex cu, sctx, cup, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop), cu(3,nxvh,nyv), sctx(nxvh)
      dimension cup(3,nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external DJPOST2GL
      data nargs /12/
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
      call MP_TASKSTART(idtask(i),DJPOST2GL,nargs,part(1,npo),cup(1,1,1,
     1i),sctxp(1,i),qm,dt,npp,idimp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call DJPOST2GL(part(1,npo),cu,sctx,qm,dt,npl,idimp,nx,ny,nxvh,nyv,
     1ipbc)
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
      subroutine MBPUSH23GL(part,fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,nop,nx
     1,ny,nxvh,nyv,ipbc,sctxp,ekp,idtask,nmt,ierr)
c multitasking gridless magnetized particle push
c sctxp = scratch arrays for sines and cosines
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qbm, dt, dtc, ek, ekp
      complex fxy, bxy, sctx, sctxp
      integer idimp, nop, nx, ny, nxvh, nyv, ipbc, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), bxy(3,nxvh,nyv), sctx(nxvh)
      dimension sctxp(nxvh,nmt), ekp(nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i
      external BPUSH23GL
      data nargs /15/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start particle push tasks
      do 10 i = 1, nmt
      npo = npp*(i - 1) + 1
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),BPUSH23GL,nargs,part(1,npo),fxy,bxy,sc
     1txp(1,i),qbm,dt,dtc,ekp(i),idimp,npp,nx,ny,nxvh,nyv,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   10 continue
c push remaining particles
      npo = npl + 1
      npl = nop - npl
      call BPUSH23GL(part(1,npo),fxy,bxy,sctx,qbm,dt,dtc,ek,idimp,npl,nx
     1,ny,nxvh,nyv,ipbc)
c wait for tasks to complete
      do 20 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   20 continue
      return
      end
