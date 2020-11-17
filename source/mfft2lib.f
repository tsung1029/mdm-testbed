c 2d PIC multi-tasking library for fast fourier transforms
c mfft2lib.f contains multi-tasking procedures to perform ffts:
c MFFT2RX multi-tasking wrapper for WFFT2RX
c MFFT2R2 multi-tasking wrapper for WFFT2R2
c MFFT2R3 multi-tasking wrapper for WFFT2R3
c MFFT2RN multi-tasking wrapper for WFFT2RN
c MFSST2RX multi-tasking wrapper for WFSST2RX
c MFSCT2RX multi-tasking wrapper for WFSCT2RX
c MFCST2RX multi-tasking wrapper for WFCST2RX
c MFCCT2RX multi-tasking wrapper for WFCCT2RX
c MFCST2R2 multi-tasking wrapper for WFCST2R2
c MFSCT2R2 multi-tasking wrapper for WFSCT2R2
c MFCST2R3 multi-tasking wrapper for WFCST2R3
c MFSCT2R3 multi-tasking wrapper for WFSCT2R3
c MFSFT2RX multi-tasking wrapper for WFSFT2RX
c MFCFT2RX multi-tasking wrapper for WFCFT2RX
c MFCSFT2R2 multi-tasking wrapper for WFCSFT2R2
c MFSCFT2R2 multi-tasking wrapper for WFSCFT2R2
c MFCSFT2R3 multi-tasking wrapper for WFCSFT2R3
c MFSCFT2R3 multi-tasking wrapper for WFSCFT2R3
c MFDSFT2RX multi-tasking wrapper for WFDSFT2RX
c MFDCFT2RX multi-tasking wrapper for WFDCFT2RX
c MFDCSFT2R2 multi-tasking wrapper for WFDCSFT2R2
c MFDSCFT2R2 multi-tasking wrapper for WFDSCFT2R2
c MFDCSFT2R3 multi-tasking wrapper for WFDCSFT2R3
c MFDSCFT2R3 multi-tasking wrapper for WFDSCFT2R3
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: december 3, 2009
c-----------------------------------------------------------------------
      subroutine MFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyh
     1d,nxyip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      integer nxyip, iftask, nmt, ierr
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nxh, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FFT2RXX, FFT2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nxh = 2**(indx - 1)
      ny = 2**indy
      nxp = nxh/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nxh - nxi
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RXX,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nypl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RXY,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RXY,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x fft tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RXX,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nypl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyh
     1d,nxyip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nxh, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FFT2R2X, FFT2R2Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nxh = 2**(indx - 1)
      ny = 2**indy
      nxp = nxh/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nxh - nxi
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R2X,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,nypl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R2Y,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R2Y,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x fft tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R2X,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,nypl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyh
     1d,nxyip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nxh, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FFT2R3X, FFT2R3Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nxh = 2**(indx - 1)
      ny = 2**indy
      nxp = nxh/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nxh - nxi
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R3X,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nypl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R3Y,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R3Y,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x fft tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2R3X,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nypl,nxhd,nyd,nxhy
     1d,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFFT2RN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,ndim,nx
     1hyd,nxyhd,nxyip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, ndim, nxhyd, nxyhd
      integer nxyip, iftask, nmt, ierr
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim,nxhd,nmt+1)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nmtx, nmty, nxh, ny, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl
      external FFT2RNX, FFT2RNY
      data nargs, margs  /13,14/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nxh = 2**(indx - 1)
      ny = 2**indy
      nmtt = nmt + 1
      nxp = nxh/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nxh - nxi
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RNX,margs,f,ss(1,1,i),isign,mix
     1up,sct,indx,indy,nxyip(i),nyp,nxhd,nyd,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call FFT2RNX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,nyi,nypl,
     1nxhd,nyd,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RNY,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,ndim
     1,nxhyd,nxyhd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RNY,nargs,f,isign,mixup,sct,ind
     1x,indy,nxyip(i),nxp,nxhd,nyd,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxpl,nxhd,nyd,ndim
     1,nxhyd,nxyhd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x fft tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FFT2RNX,margs,f,ss(1,1,i),isign,mix
     1up,sct,indx,indy,nxyip(i),nyp,nxhd,nyd,ndim,nxhyd,nxyhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x fft
         call FFT2RNX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,nyi,nypl,
     1nxhd,nyd,ndim,nxhyd,nxyhd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking real sine/sine transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FST2RXX, FST2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y sine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y sine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking real sine/cosine transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FST2RXX, FCT2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y cosine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y cosine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking real cosine/sine transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FCT2RXX, FST2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y sine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y sine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking real cosine/cosine transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FCT2RXX, FCT2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y cosine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y cosine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXY,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCST2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking 2 real sine-cosine transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FCST2R2X, FSCT2R2Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCST2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y sine-cosine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCT2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transforms
         call FSCT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y sine-cosine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCT2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call FSCT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCST2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSCT2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking 2 real sine-cosine transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FSCT2R2X, FCST2R2Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCT2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y sine-cosine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCST2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transforms
         call FCST2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y sine-cosine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCST2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call FCST2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCT2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCST2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking 3 real sine-cosine transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FCSST2R3X, FSCST2R3Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCSST2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y sine-cosine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCST2R3Y,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transforms
         call FSCST2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y sine-cosine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCST2R3Y,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call FSCST2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCSST2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSCT2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd,nxyip,iftask,nmt,ierr)
c multi-tasking 3 real sine-cosine transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nxi, nyi, nxp, nyp, nxpl, nypl
      integer i
      external FSCCT2R3X, FCSCT2R3Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nxp = nx/(nmt + 1)
      nyp = ny/(nmt + 1)
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi + 1
      nxi = nxi + 1
      nyi = nyi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCCT2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start y sine-cosine tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCSCT2R3Y,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y sine-cosine transforms
         call FCSCT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y sine-cosine tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCSCT2R3Y,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nxp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y sine-cosine transforms
         call FCSCT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCCT2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxhy
     1d,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking real sine/periodic transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension ss(2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FST2RXX, FDFT2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FST2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxhy
     1d,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking real cosine/periodic transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension ss(2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FCT2RXX, FDFT2RXY
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCT2RXX,nargs,f,isign,mixup,sctd,in
     1dx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nxh
     1yd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCSFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 2 real cosine-sine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension ss(2,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FCST2R2X, FDFT2R2Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine-sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCST2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x cosine-sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCST2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSCFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 2 real sine-cosine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension ss(2,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FSCT2R2X, FDFT2R2Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCT2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCT2R2X,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFCSFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 3 real cosine-sine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension ss(3,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FCSST2R3X, FDFT2R3Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine-sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCSST2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x cosine-sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FCSST2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFSCFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 3 real sine-cosine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension ss(3,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FSCCT2R3X, FDFT2R3Y
      data nargs /12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCCT2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FSCCT2R3X,nargs,f,isign,mixup,sctd,
     1indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nypl,nxhd,nyd,n
     1xhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFDSFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,ny
     1hd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking real sine/periodic transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd, sctdx
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      dimension ss(2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FDST2RXX, FDFT2RXY
      data nargs, margs /13, 12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDST2RXX,nargs,f,isign,mixup,sctd,s
     1ctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine transform
         call FDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd,
     1nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,nargs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDST2RXX,margs,f,isign,mixup,sctd,s
     1ctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call FDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd,
     1nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFDCFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,ny
     1hd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking real cosine/periodic transform
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd, sctdx
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      dimension ss(2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FDCT2RXX, FDFT2RXY
      data nargs, margs /13, 12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDCT2RXX,nargs,f,isign,mixup,sctd,s
     1ctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine transform
         call FDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd,
     1nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2RXY,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDCT2RXX,nargs,f,isign,mixup,sctd,s
     1ctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine transform
         call FDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd,
     1nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFDCSFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 2 real cosine-sine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd, sctdx
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(2,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FDCST2R2X, FDFT2R2Y
      data nargs, margs /13, 12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine-sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDCST2R2X,nargs,f,isign,mixup,sctd,
     1sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine-sine transforms
         call FDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd
     1,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x cosine-sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDCST2R2X,nargs,f,isign,mixup,sctd,
     1sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine-sine transforms
         call FDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd
     1,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFDSCFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 2 real sine-cosine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd, sctdx
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(2,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FDSCT2R2X, FDFT2R2Y
      data nargs, margs /13, 12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDSCT2R2X,nargs,f,isign,mixup,sctd,
     1sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd
     1,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R2Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDSCT2R2X,nargs,f,isign,mixup,sctd,
     1sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxhd
     1,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFDCSFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 3 real cosine-sine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd, sctdx
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(3,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FDCSST2R3X, FDFT2R3Y
      data nargs, margs /13, 12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x cosine-sine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDCSST2R3X,nargs,f,isign,mixup,sctd
     1,sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x cosine-sine transforms
         call FDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxh
     1d,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x cosine-sine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDCSST2R3X,nargs,f,isign,mixup,sctd
     1,sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x cosine-sine transforms
         call FDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxh
     1d,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MFDSCFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd,nxyip,iftask,nmt,ierr)
c multi-tasking for 3 real sine-cosine/periodic transforms
c nxyip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      complex f, ss, sctd, sctdx
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      integer nxyip, iftask, nmt, ierr
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(3,2*nxhd)
      dimension nxyip(nmt), iftask(nmt)
c local data
      integer nargs, margs, nmtx, nmty, nx, ny, nx1, nyh, nmtt, nxi, nyi
      integer nxp, nyp, nxpl, nypl, nxd, nyd, i
      external FDSCCT2R3X, FDFT2R3Y
      data nargs, margs /13, 12/
c calculate range of indices
      nmtx = nmt
      nmty = nmt
      nx = 2**indx
      ny = 2**indy
      nx1 = nx + 1
      nyh = ny/2
      nmtt = nmt + 1
      nxp = nx/nmtt
      nyp = ny/nmtt
      if (nxp.eq.0) nmtx = 0
      if (nyp.eq.0) nmty = 0
      nxi = nxp*nmt
      nyi = nyp*nmt
      nxpl = nx - nxi + 1
      nypl = ny - nyi
      nxi = nxi + 1
      nyi = nyi + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x sine-cosine tasks
         do 10 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDSCCT2R3X,nargs,f,isign,mixup,sctd
     1,sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x sine-cosine transforms
         call FDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxh
     1d,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 20 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start y fft tasks
         do 30 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 40 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c start y fft tasks
         do 50 i = 1, nmtx
         nxyip(i) = nxp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDFT2R3Y,margs,f,isign,mixup,sctd,i
     1ndx,indy,nxyip(i),nxp,nxd,nyhd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish y fft
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxpl,nxd,nyhd,nx
     1hyd,nxyd)
c wait for tasks to complete
         do 60 i = 1, nmtx
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   60    continue
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,1,nx1,nxd,nyhd)
c start x sine-cosine tasks
         do 70 i = 1, nmty
         nxyip(i) = nyp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),FDSCCT2R3X,nargs,f,isign,mixup,sctd
     1,sctdx,indx,indy,nxyip(i),nyp,nxhd,nyd,nxhyd,nxyd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish x sine-cosine transforms
         call FDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nypl,nxh
     1d,nyd,nxhyd,nxyd)
c wait for tasks to complete
         do 80 i = 1, nmty
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
      endif
      return
      end
