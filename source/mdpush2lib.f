c-----------------------------------------------------------------------
c 2d PIC multi-tasking library for pushing particles with darwin
c electric and magnetic fields and depositing current and derivative
c of current
c mdpush2lib.f contains multi-tasking procedures to process particles
c             with darwin electric and magnetic fields:
c MGMJPOST2 multi-tasking wrapper for GMJPOST2
c MGSMJPOST2 multi-tasking wrapper for GSMJPOST2
c MGMJPOST2L multi-tasking wrapper for GMJPOST2L
c MGSMJPOST2L multi-tasking wrapper for GSMJPOST2L
c MGMJPOST22 multi-tasking wrapper for GMJPOST22
c MGSMJPOST22 multi-tasking wrapper for GSMJPOST22
c MGMJPOST22L multi-tasking wrapper for GMJPOST22L
c MGSMJPOST22L multi-tasking wrapper for GSMJPOST22L
c MGMJPOST2C multi-tasking wrapper for GMJPOST2C
c MGMJPOST22C multi-tasking wrapper for GMJPOST22C
c MGDCJPOST2 multi-tasking wrapper for GDCJPOST2
c MGSDCJPOST2 multi-tasking wrapper for GSDCJPOST2
c MGDCJPOST2L multi-tasking wrapper for GDCJPOST2L
c MGSDCJPOST2L multi-tasking wrapper for GSDCJPOST2L
c MGDCJPOST22 multi-tasking wrapper for GDCJPOST22
c MGSDCJPOST22 multi-tasking wrapper for GSDCJPOST22
c MGDCJPOST22L multi-tasking wrapper for GDCJPOST22L
c MGSDCJPOST22L multi-tasking wrapper for GSDCJPOST22L
c MGDCJPOST2C multi-tasking wrapper for GDCJPOST2C
c MGDCJPOST22C multi-tasking wrapper for GDCJPOST22C
c MGMJPOST2GL multi-tasking wrapper for GMJPOST2GL
c MGDCJPOST2GL multi-tasking wrapper for GDCJPOST2GL
c written by viktor k. decyk, ucla
c copyright 2006, regents of the university of california
c update: november 27, 2009
c-----------------------------------------------------------------------
      subroutine MGMJPOST2(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,nmt
     1,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxv,nyv)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST2
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 4
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GMJPOST2,nargs,part(1,npo),amup(1,1,1,
     1i),qm,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST2(part(1,npo),amu,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
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
      subroutine MGSMJPOST2(part,amu,qm,nop,idimp,nxv,nxyv,amup,idtask,n
     1mt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxyv)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSMJPOST2
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 4
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSMJPOST2,nargs,part(1,npo),amup(1,1,i
     1),qm,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GSMJPOST2(part(1,npo),amu,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 4
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGMJPOST2L(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,nm
     1t,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxv,nyv)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST2L
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 4
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GMJPOST2L,nargs,part(1,npo),amup(1,1,1
     1,i),qm,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST2L(part(1,npo),amu,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
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
      subroutine MGSMJPOST2L(part,amu,qm,nop,idimp,nxv,nxyv,amup,idtask,
     1nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxyv)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSMJPOST2L
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 4
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSMJPOST2L,nargs,part(1,npo),amup(1,1,
     1i),qm,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GSMJPOST2L(part(1,npo),amu,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 4
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGMJPOST22(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,nm
     1t,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxv,nyv)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST22
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 2
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GMJPOST22,nargs,part(1,npo),amup(1,1,1
     1,i),qm,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST22(part(1,npo),amu,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
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
      subroutine MGSMJPOST22(part,amu,qm,nop,idimp,nxv,nxyv,amup,idtask,
     1nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxyv)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSMJPOST22
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 2
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSMJPOST22,nargs,part(1,npo),amup(1,1,
     1i),qm,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GSMJPOST22(part(1,npo),amu,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGMJPOST22L(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,n
     1mt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxv,nyv)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST22L
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 2
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GMJPOST22L,nargs,part(1,npo),amup(1,1,
     11,i),qm,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST22L(part(1,npo),amu,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
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
      subroutine MGSMJPOST22L(part,amu,qm,nop,idimp,nxv,nxyv,amup,idtask
     1,nmt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxyv)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSMJPOST22L
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 30 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 20 j = 1, nxyv
      do 10 m = 1, 2
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSMJPOST22L,nargs,part(1,npo),amup(1,1
     1,i),qm,npp,idimp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GSMJPOST22L(part(1,npo),amu,qm,npl,idimp,nxv,nxyv)
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
      do 50 j = 1, nxyv
      do 40 m = 1, 2
      amu(m,j) = amu(m,j) + amup(m,j,i)
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGMJPOST2C(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,nm
     1t,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxv,nyv)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST2C
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 4
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GMJPOST2C,nargs,part(1,npo),amup(1,1,1
     1,i),qm,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST2C(part(1,npo),amu,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
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
      subroutine MGMJPOST22C(part,amu,qm,nop,idimp,nxv,nyv,amup,idtask,n
     1mt,ierr)
c multitasking momentum flux deposition
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, amu, qm, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(2,nxv,nyv)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST22C
      data nargs /7/
      npp = nop/(nmt + 1)
      npl = npp*nmt
      ierr = 0
c start momentum flux deposit tasks
      do 40 i = 1, nmt
      npo = npp*(i - 1) + 1
c clear momentum flux arrays
      do 30 k = 1, nyv
      do 20 j = 1, nxv
      do 10 m = 1, 2
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GMJPOST22C,nargs,part(1,npo),amup(1,1,
     11,i),qm,npp,idimp,nxv,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST22C(part(1,npo),amu,qm,npl,idimp,nxv,nyv)
c wait for tasks to complete
      do 80 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum momentum flux arrays
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
      subroutine MGDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,nop,
     1nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), dcup(3,nxv,nyv,nmt)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GDCJPOST2
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
      amup(4,j,k,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GDCJPOST2,nargs,part(1,npo),fxy,bxy,cu
     1p(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv
     2)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST2(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,npl,
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
      subroutine MGSDCJPOST2(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      dimension cup(3,nxyv,nmt), dcup(3,nxyv,nmt)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSDCJPOST2
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
      amup(4,j,i) = 0.
   20 continue
      call MP_TASKSTART(idtask(i),GSDCJPOST2,nargs,part(1,npo),fxy,bxy,c
     1up(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSDCJPOST2(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,npl
     1,nxv,nxyv)
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
      subroutine MGDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), dcup(3,nxv,nyv,nmt)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GDCJPOST2L
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
      amup(4,j,k,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GDCJPOST2L,nargs,part(1,npo),fxy,bxy,c
     1up(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,ny
     2v)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST2L(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,npl
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
      subroutine MGSDCJPOST2L(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,no
     1p,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxyv), bxy(3,nxyv)
      dimension cu(3,nxyv), dcu(3,nxyv), amu(4,nxyv)
      dimension cup(3,nxyv,nmt), dcup(3,nxyv,nmt)
      dimension amup(4,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSDCJPOST2L
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
      amup(4,j,i) = 0.
   20 continue
      call MP_TASKSTART(idtask(i),GSDCJPOST2L,nargs,part(1,npo),fxy,bxy,
     1cup(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSDCJPOST2L(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,np
     1l,nxv,nxyv)
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
      subroutine MGDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,nop,
     1nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), dcup(2,nxv,nyv,nmt)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GDCJPOST22
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GDCJPOST22,nargs,part(1,npo),fxy,bz,cu
     1p(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,nyv
     2)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST22(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,npl,
     1nxv,nyv)
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
      subroutine MGSDCJPOST22(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxyv), bz(nxyv)
      dimension cu(2,nxyv), dcu(2,nxyv), amu(2,nxyv)
      dimension cup(2,nxyv,nmt), dcup(2,nxyv,nmt)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSDCJPOST22
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSDCJPOST22,nargs,part(1,npo),fxy,bz,c
     1up(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSDCJPOST22(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,npl
     1,nxv,nxyv)
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
      subroutine MGDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), dcup(2,nxv,nyv,nmt)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GDCJPOST22L
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GDCJPOST22L,nargs,part(1,npo),fxy,bz,c
     1up(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,ny
     2v)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST22L(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,npl
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
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGSDCJPOST22L(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,no
     1p,nxv,nxyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nxyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxyv), bz(nxyv)
      dimension cu(2,nxyv), dcu(2,nxyv), amu(2,nxyv)
      dimension cup(2,nxyv,nmt), dcup(2,nxyv,nmt)
      dimension amup(2,nxyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, m
      external GSDCJPOST22L
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
      dcup(m,j,i) = 0.
      amup(m,j,i) = 0.
   10 continue
   20 continue
      call MP_TASKSTART(idtask(i),GSDCJPOST22L,nargs,part(1,npo),fxy,bz,
     1cup(1,1,i),dcup(1,1,i),amup(1,1,i),qm,qbm,dt,idimp,npp,nxv,nxyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GSDCJPOST22L(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,np
     1l,nxv,nxyv)
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
      subroutine MGDCJPOST2C(part,fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bxy, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(3,nxv,nyv), bxy(3,nxv,nyv)
      dimension cu(3,nxv,nyv), dcu(3,nxv,nyv), amu(4,nxv,nyv)
      dimension cup(3,nxv,nyv,nmt), dcup(3,nxv,nyv,nmt)
      dimension amup(4,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GDCJPOST2C
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
      amup(4,j,k,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GDCJPOST2C,nargs,part(1,npo),fxy,bxy,c
     1up(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,ny
     2v)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST2C(part(1,npo),fxy,bxy,cu,dcu,amu,qm,qbm,dt,idimp,npl
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
      subroutine MGDCJPOST22C(part,fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,nop
     1,nxv,nyv,cup,dcup,amup,idtask,nmt,ierr)
c multitasking momentum flux, acceleration, and current deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxy, bz, cu, dcu, amu, qm, qbm, dt, cup, dcup, amup
      integer nop, idimp, nxv, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), fxy(2,nxv,nyv), bz(nxv,nyv)
      dimension cu(2,nxv,nyv), dcu(2,nxv,nyv), amu(2,nxv,nyv)
      dimension cup(2,nxv,nyv,nmt), dcup(2,nxv,nyv,nmt)
      dimension amup(2,nxv,nyv,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GDCJPOST22C
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
      dcup(m,j,k,i) = 0.
      amup(m,j,k,i) = 0.
   10 continue
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),GDCJPOST22C,nargs,part(1,npo),fxy,bz,c
     1up(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),qm,qbm,dt,idimp,npp,nxv,ny
     2v)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST22C(part(1,npo),fxy,bz,cu,dcu,amu,qm,qbm,dt,idimp,npl
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
      dcu(m,j,k) = dcu(m,j,k) + dcup(m,j,k,i)
      amu(m,j,k) = amu(m,j,k) + amup(m,j,k,i)
   50 continue
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MGMJPOST2GL(part,amu,sctx,qm,nop,idimp,nx,ny,nxvh,nyv,a
     1mup,sctxp,idtask,nmt,ierr)
c multitasking gridless momentum flux deposition
c amup = momentum flux density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm
      complex amu, sctx, amup, sctxp
      integer nop, idimp, nx, ny, nxvh, nyv, idtask, nmt, ierr
      dimension part(idimp,nop), amu(4,nxvh,nyv), sctx(nxvh)
      dimension amup(4,nxvh,nyv,nmt), sctxp(nxvh,nmt), idtask(nmt)
c local data
      integer nargs, npp, npl, npo, i, j, k, m
      external GMJPOST2GL
      data nargs /10/
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
      call MP_TASKSTART(idtask(i),GMJPOST2GL,nargs,part(1,npo),amup(1,1,
     11,i),sctxp(1,i),qm,npp,idimp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GMJPOST2GL(part(1,npo),amu,sctx,qm,npl,idimp,nx,ny,nxvh,nyv)
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
      subroutine MGDCJPOST2GL(part,fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,idi
     1mp,nop,nx,ny,nxvh,nyv,cup,dcup,amup,sctxp,idtask,nmt,ierr)
c multitasking gridless momentum flux, acceleration density, and current
c deposition
c cup = current density arrays for tasks
c dcup = acceleration density arrays for tasks
c amup = momentum flux density arrays for tasks
c sctxp = scratch arrays for sines and cosines
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, qm, qbm, dt
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
      external GDCJPOST2GL
      data nargs /16/
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
      call MP_TASKSTART(idtask(i),GDCJPOST2GL,nargs,part(1,npo),fxy,bxy,
     1cup(1,1,1,i),dcup(1,1,1,i),amup(1,1,1,i),sctxp(1,i),qm,qbm,dt,idim
     2p,npp,nx,ny,nxvh,nyv)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining momentum flux
      npo = npl + 1
      npl = nop - npl
      call GDCJPOST2GL(part(1,npo),fxy,bxy,cu,dcu,amu,sctx,qm,qbm,dt,idi
     1mp,npl,nx,ny,nxvh,nyv)
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

