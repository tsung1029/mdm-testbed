c-----------------------------------------------------------------------
c 2d PIC multi-tasking library for pushing relativistic particles with
c darwin electric and magnetic fields and depositing current and
c derivative of current
c mrdpush2lib.f contains multi-tasking procedures to process
c               relativistic particles with darwin electric and magnetic
c               fields:
c MGRMJPOST2 multi-tasking wrapper for GRMJPOST2
c MGSRMJPOST2 multi-tasking wrapper for GSRMJPOST2
c MGRMJPOST2L multi-tasking wrapper for GRMJPOST2L
c MGSRMJPOST2L multi-tasking wrapper for GSRMJPOST2L
c MGRMJPOST22 multi-tasking wrapper for GRMJPOST22
c MGSRMJPOST22 multi-tasking wrapper for GSRMJPOST22
c MGRMJPOST22L multi-tasking wrapper for GRMJPOST22L
c MGSRMJPOST22L multi-tasking wrapper for GSRMJPOST22L
c MGRDCJPOST2 multi-tasking wrapper for GRDCJPOST2
c MGSRDCJPOST2 multi-tasking wrapper for GSRDCJPOST2
c MGRDCJPOST2L multi-tasking wrapper for GRDCJPOST2L
c MGSRDCJPOST2L multi-tasking wrapper for GSRDCJPOST2L
c MGRDCJPOST22 multi-tasking wrapper for GRDCJPOST22
c MGSRDCJPOST22 multi-tasking wrapper for GSRDCJPOST22
c MGRDCJPOST22L multi-tasking wrapper for GRDCJPOST22L
c MGSRDCJPOST22L multi-tasking wrapper for GSRDCJPOST22L
c MGRMJPOST2C multi-tasking wrapper for GRMJPOST2C
c MGRMJPOST22C multi-tasking wrapper for GRMJPOST22C
c MGRDCJPOST2C multi-tasking wrapper for GRDCJPOST2C
c MGRDCJPOST22C multi-tasking wrapper for GRDCJPOST22C
c MGRMJPOST2GL multi-tasking wrapper for GRMJPOST2GL
c MGRDCJPOST2GL multi-tasking wrapper for GRDCJPOST2GL
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: december 5, 2009
c-----------------------------------------------------------------------
      subroutine MGRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idtask
     1,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxv,nyv)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST2
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 4
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST2,nargs,part(1,npo),amup(1,1,1
     1,i),qm,ci,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST2(part(1,npo),amu,qm,ci,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 4
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRMJPOST2(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,idta
     1sk,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxyv)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRMJPOST2
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 4
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRMJPOST2,nargs,part(1,npo),amup(1,1,
     1i),qm,ci,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRMJPOST2(part(1,npo),amu,qm,ci,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 4
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idtas
     1k,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxv,nyv)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST2L
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 4
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST2L,nargs,part(1,npo),amup(1,1,
     11,i),qm,ci,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST2L(part(1,npo),amu,qm,ci,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 4
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRMJPOST2L(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,idt
     1ask,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxyv)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRMJPOST2L
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 4
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRMJPOST2L,nargs,part(1,npo),amup(1,1
     1,i),qm,ci,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRMJPOST2L(part(1,npo),amu,qm,ci,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 4
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idtas
     1k,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxv,nyv)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST22
      data nargs /8/
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
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST22,nargs,part(1,npo),amup(1,1,
     11,i),qm,ci,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST22(part(1,npo),amu,qm,ci,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRMJPOST22(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,idt
     1ask,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxyv)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRMJPOST22
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 2
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRMJPOST22,nargs,part(1,npo),amup(1,1
     1,i),qm,ci,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRMJPOST22(part(1,npo),amu,qm,ci,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idta
     1sk,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxv,nyv)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST22L
      data nargs /8/
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
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST22L,nargs,part(1,npo),amup(1,1
     1,1,i),qm,ci,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST22L(part(1,npo),amu,qm,ci,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRMJPOST22L(part,amu,qm,ci,nop,idimp,nxv,nxyv,amup,id
     1task,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxyv)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRMJPOST22L
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 2
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRMJPOST22L,nargs,part(1,npo),amup(1,
     11,i),qm,ci,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRMJPOST22L(part(1,npo),amu,qm,ci,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), dcup(3,nxv,nyv,nmt)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST2
      data nargs /14/
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
      amup(4,j,k,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST2,nargs,part(1,npo),fxy,bxy,c
     1up(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv
     2,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST2(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1npl,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
      amu(4,j,k) = amu(4,j,k) + amup(4,j,k,i)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      dimension cup(3,nxyv,nmt), dcup(3,nxyv,nmt)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRDCJPOST2
      data nargs /14/
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
      amup(4,j,i) = 0.
   20 continue
      call MP_TASKSTART(idtask(i),GSRDCJPOST2,nargs,part(1,npo),fxy,bxy,
     1cup(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxyv
     2)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRDCJPOST2(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,npl,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 3
      cu(m,j) = cu(m,j) + cup(m,j,i)
      dcu(m,j) = dcu(m,j) + dcup(m,j,i)
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
      amu(4,j) = amu(4,j) + amup(4,j,i)
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), dcup(3,nxv,nyv,nmt)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST2L
      data nargs /14/
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
      amup(4,j,k,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST2L,nargs,part(1,npo),fxy,bxy,
     1cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nx
     2v,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST2L(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,npl,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
      amu(4,j,k) = amu(4,j,k) + amup(4,j,k,i)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idim
     1p,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      dimension cup(3,nxyv,nmt), dcup(3,nxyv,nmt)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRDCJPOST2L
      data nargs /14/
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
      amup(4,j,i) = 0.
   20 continue
      call MP_TASKSTART(idtask(i),GSRDCJPOST2L,nargs,part(1,npo),fxy,bxy
     1,cup(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxy
     2v)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRDCJPOST2L(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idim
     1p,npl,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 3
      cu(m,j) = cu(m,j) + cup(m,j,i)
      dcu(m,j) = dcu(m,j) + dcup(m,j,i)
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
      amu(4,j) = amu(4,j) + amup(4,j,i)
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), dcup(2,nxv,nyv,nmt)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST22
      data nargs /14/
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST22,nargs,part(1,npo),fxy,bz,c
     1up(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nxv
     2,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST22(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp,
     1npl,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxyv), bz(nxyv)
      dimension cu(2,nxyv), dcu(2,nxyv), amu(2,nxyv)
      dimension cup(2,nxyv,nmt), dcup(2,nxyv,nmt)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRDCJPOST22
      data nargs /14/
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRDCJPOST22,nargs,part(1,npo),fxy,bz,
     1cup(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxyv
     2)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRDCJPOST22(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,npl,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      cu(m,j) = cu(m,j) + cup(m,j,i)
      dcu(m,j) = dcu(m,j) + dcup(m,j,i)
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), dcup(2,nxv,nyv,nmt)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST22L
      data nargs /14/
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST22L,nargs,part(1,npo),fxy,bz,
     1cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nx
     2v,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST22L(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,npl,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSRDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idim
     1p,nop,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxyv), bz(nxyv)
      dimension cu(2,nxyv), dcu(2,nxyv), amu(2,nxyv)
      dimension cup(2,nxyv,nmt), dcup(2,nxyv,nmt)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSRDCJPOST22L
      data nargs /14/
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSRDCJPOST22L,nargs,part(1,npo),fxy,bz
     1,cup(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,ci,idimp,npp,nxv,nxy
     2v)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSRDCJPOST22L(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idim
     1p,npl,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      cu(m,j) = cu(m,j) + cup(m,j,i)
      dcu(m,j) = dcu(m,j) + dcup(m,j,i)
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRMJPOST2C(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idtas
     1k,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxv,nyv)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST2C
      data nargs /8/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear current arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 4
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST2C,nargs,part(1,npo),amup(1,1,
     11,i),qm,ci,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST2C(part(1,npo),amu,qm,ci,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 4
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRMJPOST22C(part,amu,qm,ci,nop,idimp,nxv,nyv,amup,idta
     1sk,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, ci, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxv,nyv)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST22C
      data nargs /8/
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
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST22C,nargs,part(1,npo),amup(1,1
     1,1,i),qm,ci,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST22C(part(1,npo),amu,qm,ci,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), dcup(3,nxv,nyv,nmt)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST2C
      data nargs /14/
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
      amup(4,j,k,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST2C,nargs,part(1,npo),fxy,bxy,
     1cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nx
     2v,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST2C(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,npl,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
      amu(4,j,k) = amu(4,j,k) + amup(4,j,k,i)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,nop,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, ci, cup, dcup
      real amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), dcup(2,nxv,nyv,nmt)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST22C
      data nargs /14/
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST22C,nargs,part(1,npo),fxy,bz,
     1cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,ci,idimp,npp,nx
     2v,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST22C(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,ci,idimp
     1,npl,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxv
      do 50 m = 1, 2
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRMJPOST2GL(part,amu,sctx,qm,ci,nop,idimp,nx,ny,nxvh,n
     1yv,amup,sctxp,idtask,nmt,ierr)
c multitasking gridless momentum flux deposition
c for relativistic particles
c amup = momentum flux density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm, ci
      complex amu, sctx, amup, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxvh,nyv), sctx(nxvh)
      dimension amup(4,nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRMJPOST2GL
      data nargs /11/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxvh
      do 10 m = 1, 4
      amup(m,j,k,i) = cmplx(0.,0.)
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GRMJPOST2GL,nargs,part(1,npo),amup(1,1
     1,1,i),sctxp(1,i),qm,ci,npp,idimp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GRMJPOST2GL(part(1,npo),amu,sctx,qm,ci,npl,idimp,nx,ny,nxvh,n
     1yv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 70 k = 1, nyv
      do 60 j = 1, nxvh
      do 50 m = 1, 4
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGRDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,ci
     1,idimp,nop,nx,ny,nxvh,nyv,cup,dcup,amup,sctxp,idtask,nmt,ierr)
c multitasking gridless momentum flux, acceleration density, and current
c deposition for relativistic particles
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm, qbm, dt, ci
      complex fxy, bxy, cu, dcu, amu, sctx, cup, dcup, amup, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, idtask, nmt, ierr
      dimension part(idimp,nop)
      dimension fxy(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension cu(3,nxvh,nyv), dcu(3,nxvh,nyv), amu(4,nxvh,nyv)
      dimension sctx(nxvh)
      dimension cup(3,nxvh,nyv,nmt), dcup(3,nxvh,nyv,nmt)
      dimension amup(4,nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GRDCJPOST2GL
      data nargs /17/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 60 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 50 k = 1, nyv
      do 40 j = 1, nxvh
      do 10 m = 1, 3
      cup(m,j,k,i) = cmplx(0.,0.)
   10 continue
      do 20 m = 1, 3
      dcup(m,j,k,i) = cmplx(0.,0.)
   20 continue
      do 30 m = 1, 4
      amup(m,j,k,i) = cmplx(0.,0.)
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),GRDCJPOST2GL,nargs,part(1,npo),fxy,bxy
     1,cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),sctxp(1,i),qm,qbm,dt,ci,
     2idimp,npp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GRDCJPOST2GL(part(1,npo),fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,ci
     1,idimp,npl,nx,ny,nxvh,nyv)
c wait for tasks to complete
      do 120 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 110 k = 1, nyv
      do 100 j = 1, nxvh
      do 70 m = 1, 3
      cu(m,j,k) = cu(m,j,k) + cup(m,j,k,i)
   70 continue
      do 80 m = 1, 3
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
   80 continue
      do 90 m = 1, 4
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   90 continue
  100 continue
  110 continue
  120 continue
      return
      end
