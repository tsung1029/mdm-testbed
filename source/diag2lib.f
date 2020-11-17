c-----------------------------------------------------------------------
c 2d PIC library for diagnostics
c diag2lib.f contains diagnostic procedures:
c VDIST2 calculates 2 component velocity distribution and velocity
c        moments.
c VDIST23 calculates 3 component velocity distribution and velocity
c         moments.
c ADIST2 calculates 2 component field part of canonical momentum
c        distribution and their moments.
c ADIST23 calculates 3 component field part of canonical momentum
c         distribution and their moments.
c FCWRITE2 writes 2d complex binary data to file.
c FCREAD2 writes 2d complex binary data to file.
c NDIAN determines whether number format is big or little endian.
c NDPREC determines whether default reals are double precision.
c IDPREC determines whether default integers are double precision.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: november 4, 2009
c-----------------------------------------------------------------------
      subroutine VDIST2(part,fv,fvm,idimp,np,nmv,nmvf)
c for 2d code, this subroutine calculates 2d velocity distribution
c and velocity moments
c input: all, output: fv
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift for i-th dimension is contained in fvm(1,i)
c vth for i-th dimension is contained in fvm(2,i)
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf,2), fvm(3,2)
      double precision sumvx, sumvy, sumvx2, sumvy2, anp
      real anmv, svx, svy
      integer j, nvx, nvy
      anmv = real(nmv)
      svx = anmv/fv(1,1)
      svy = anmv/fv(1,2)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j,1) = 0.
      fv(j,2) = 0.
   10 continue
c count particles in each velocity region
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 20 j = 1, np
      nvx = part(3,j)*svx + anmv
      sumvx = sumvx + part(3,j)
      sumvx2 = sumvx2 + part(3,j)**2
      nvy = part(4,j)*svy + anmv
      sumvy = sumvy + part(4,j)
      sumvy2 = sumvy2 + part(4,j)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx,1) = fv(nvx,1) + 1.
      if ((nvy.ge.2).and.(nvy.le.nmvf)) fv(nvy,2) = fv(nvy,2) + 1.
   20 continue
c calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(2,1) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(1,2) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j,1).gt.0.) then
         sumvx = sumvx + fv(j,1)
         sumvx2 = sumvx2 + fv(j,1)*dlog(dble(fv(j,1)*svx))
      endif
      if (fv(j,2).gt.0.) then
         sumvy = sumvy + fv(j,2)
         sumvy2 = sumvy2 + fv(j,2)*dlog(dble(fv(j,2)*svy))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      fvm(3,1) = sumvx
      fvm(3,2) = sumvy
      return
      end
      subroutine VDIST23(part,fv,fvm,idimp,np,nmv,nmvf)
c for 2-1/2d code, this subroutine calculates 3d velocity distribution
c and velocity moments
c input: all, output: fv
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift for i-th dimension is contained in fvm(1,i)
c vth for i-th dimension is contained in fvm(2,i)
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution
c is uniform in space and distributions in each dimension are
c independent.
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, np, nmv, nmvf
      real part, fv, fvm
      dimension part(idimp,np), fv(nmvf,3), fvm(3,3)
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      real anmv, svx, svy, svz
      integer j, nvx, nvy, nvz
      anmv = real(nmv)
      svx = anmv/fv(1,1)
      svy = anmv/fv(1,2)
      svz = anmv/fv(1,3)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j,1) = 0.
      fv(j,2) = 0.
      fv(j,3) = 0.
   10 continue
c count particles in each velocity region
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, np
      nvx = part(3,j)*svx + anmv
      sumvx = sumvx + part(3,j)
      sumvx2 = sumvx2 + part(3,j)**2
      nvy = part(4,j)*svy + anmv
      sumvy = sumvy + part(4,j)
      sumvy2 = sumvy2 + part(4,j)**2
      nvz = part(5,j)*svz + anmv
      sumvz = sumvz + part(5,j)
      sumvz2 = sumvz2 + part(5,j)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx,1) = fv(nvx,1) + 1.
      if ((nvy.ge.2).and.(nvy.le.nmvf)) fv(nvy,2) = fv(nvy,2) + 1.
      if ((nvz.ge.2).and.(nvz.le.nmvf)) fv(nvz,3) = fv(nvz,3) + 1.
   20 continue
c calculate velocity moments
      anp = 1.0d0/dble(np)
      sumvx = sumvx*anp
      fvm(1,1) = sumvx
      fvm(2,1) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(1,2) = sumvy
      fvm(2,2) = dsqrt(sumvy2*anp - sumvy**2)
      sumvz = sumvz*anp
      fvm(1,3) = sumvz
      fvm(2,3) = dsqrt(sumvz2*anp - sumvz**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j,1).gt.0.) then
         sumvx = sumvx + fv(j,1)
         sumvx2 = sumvx2 + fv(j,1)*dlog(dble(fv(j,1)*svx))
      endif
      if (fv(j,2).gt.0.) then
         sumvy = sumvy + fv(j,2)
         sumvy2 = sumvy2 + fv(j,2)*dlog(dble(fv(j,2)*svy))
      endif
      if (fv(j,3).gt.0.) then
         sumvz = sumvz + fv(j,3)
         sumvz2 = sumvz2 + fv(j,3)*dlog(dble(fv(j,3)*svz))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(3,1) = sumvx
      fvm(3,2) = sumvy
      fvm(3,3) = sumvz
      return
      end
      subroutine ADIST2(part,av,avm,idimp,np,nma,nmaf)
c for 2d code, this subroutine calculates 2d field part of canonical
c momentum distribution and their moments
c input: all, output: fv
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = canonical momentum px of particle n
c part(6,n) = canonical momentum py of particle n
c av = distribution function, number of particles in each momentum range
c maximum momentum (used for scaling) is contained in first element av.
c adrift for i-th dimension is contained in avm(1,i)
c ath for i-th dimension is contained in avm(2,i)
c the number of momentum bins used is 2*nma + 1, nmaf >= 2*nma+2
      implicit none
      integer idimp, np, nma, nmaf
      real part, av, avm
      dimension part(idimp,np), av(nmaf,2), avm(2,2)
      double precision sumax, sumay, sumax2, sumay2, anp
      real anma, sax, say, ax, ay
      integer j, nax, nay
      anma = real(nma)
      sax = anma/av(1,1)
      say = anma/av(1,2)
c zero out distribution
      do 10 j = 2, nmaf
      av(j,1) = 0.
      av(j,2) = 0.
   10 continue
c count particles in each momentum region
      anma = anma + 2.5
      sumax = 0.0d0
      sumay = 0.0d0
      sumax2 = 0.0d0
      sumay2 = 0.0d0
      do 20 j = 1, np
      ax = part(5,j) - part(3,j)
      nax = ax*sax + anma
      sumax = sumax + ax
      sumax2 = sumax2 + ax**2
      ay = part(6,j) - part(4,j)
      nay = ay*say + anma
      sumay = sumay + ay
      sumay2 = sumay2 + ay**2
      if ((nax.ge.2).and.(nax.le.nmaf)) av(nax,1) = av(nax,1) + 1.
      if ((nay.ge.2).and.(nay.le.nmaf)) av(nay,2) = av(nay,2) + 1.
   20 continue
c calculate field momentum moments
      anp = 1.0d0/dble(np)
      sumax = sumax*anp
      avm(1,1) = sumax
      avm(2,1) = dsqrt(sumax2*anp - sumax**2)
      sumay = sumay*anp
      avm(1,2) = sumay
      avm(2,2) = dsqrt(sumay2*anp - sumay**2)
      return
      end
      subroutine ADIST23(part,av,avm,idimp,np,nma,nmaf)
c for 2-1/2d code, this subroutine calculates 3d field part of canonical
c momentum distribution and their moments
c input: all, output: fv
c part(3,n) = velocity vx of particle n
c part(4,n) = velocity vy of particle n
c part(5,n) = velocity vz of particle n
c part(6,n) = canonical momentum px of particle n
c part(7,n) = canonical momentum py of particle n
c part(8,n) = canonical momentum pz of particle n
c av = distribution function, number of particles in each momentum range
c maximum momentum (used for scaling) is contained in first element av.
c adrift for i-th dimension is contained in avm(1,i)
c ath for i-th dimension is contained in avm(2,i)
c the number of momentum bins used is 2*nma + 1, nmaf >= 2*nma+2
      implicit none
      integer idimp, np, nma, nmaf
      real part, av, avm
      dimension part(idimp,np), av(nmaf,3), avm(2,3)
      double precision sumax, sumay, sumaz, sumax2, sumay2, sumaz2, anp
      real anma, sax, say, saz, ax, ay, az
      integer j, nax, nay, naz
      anma = real(nma)
      sax = anma/av(1,1)
      say = anma/av(1,2)
      saz = anma/av(1,3)
c zero out distribution
      do 10 j = 2, nmaf
      av(j,1) = 0.
      av(j,2) = 0.
      av(j,3) = 0.
   10 continue
c count particles in each momentum region
      anma = anma + 2.5
      sumax = 0.0d0
      sumay = 0.0d0
      sumaz = 0.0d0
      sumax2 = 0.0d0
      sumay2 = 0.0d0
      sumaz2 = 0.0d0
      do 20 j = 1, np
      ax = part(6,j) - part(3,j)
      nax = ax*sax + anma
      sumax = sumax + ax
      sumax2 = sumax2 + ax**2
      ay = part(7,j) - part(4,j)
      nay = ay*say + anma
      sumay = sumay + ay
      sumay2 = sumay2 + ay**2
      az = part(8,j) - part(5,j)
      naz = az*saz + anma
      sumaz = sumaz + az
      sumaz2 = sumaz2 + az**2
      if ((nax.ge.2).and.(nax.le.nmaf)) av(nax,1) = av(nax,1) + 1.
      if ((nay.ge.2).and.(nay.le.nmaf)) av(nay,2) = av(nay,2) + 1.
      if ((naz.ge.2).and.(naz.le.nmaf)) av(naz,3) = av(naz,3) + 1.
   20 continue
c calculate field momentum moments
      anp = 1.0d0/dble(np)
      sumax = sumax*anp
      avm(1,1) = sumax
      avm(2,1) = dsqrt(sumax2*anp - sumax**2)
      sumay = sumay*anp
      avm(1,2) = sumay
      avm(2,2) = dsqrt(sumay2*anp - sumay**2)
      sumaz = sumaz*anp
      avm(1,3) = sumaz
      avm(2,3) = dsqrt(sumaz2*anp - sumaz**2)
      return
      end
      subroutine FWRITE2(f,nx,ny,nxv,iunit,nrec,lrec,name)
c this subroutine writes real 2d data f to a direct access file
c f = input data to be written
c nx/ny = length of data f in x/y to write
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, ny, nxv, iunit, nrec, lrec, fname
c output: nrec
      implicit none
      integer nx, ny, nxv, iunit, nrec, lrec
      real f
      character*(*) name
      dimension f(nxv,ny)
c local data
      integer j, k
c open new file
      if (nrec.lt.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='replace')
         nrec = 1
c open old file
      else if (nrec.eq.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         return
      endif
      write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,ny)
      nrec = nrec + 1
      return
      end
      subroutine FREAD2(f,nx,ny,nxv,iunit,nrec,lrec,name,ierr)
c this subroutine reads real 2d data f from a file
c f = input data to be read
c nx/ny = length of data f in x/y to read
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for read, if nrec >  0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nren <= 0)
c input: nx, ny, nxv, iunit, nrec, lrec, fname
c output: f, nrec, ierr
      implicit none
      integer nx, ny, nxv, iunit, nrec, lrec, ierr
      real f
      character*(*) name
      dimension f(nxv,ny)
c local data
      integer j, k
      ierr = 0
      if (nrec.le.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         if (nrec.eq.0) return
         nrec = 1
      endif
      read (unit=iunit,rec=nrec,err=10) ((f(j,k),j=1,nx),k=1,ny)
      nrec = nrec + 1
      return
   10 ierr = 1
      return
      end
      subroutine FCWRITE2(f,nx,ny,nxv,iunit,nrec,lrec,name)
c this subroutine writes complex 2d data f to a direct access file
c f = input data to be written
c nx/ny = length of data f in x/y to write
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, ny, nxv, iunit, nrec, lrec, fname
c output: nrec
      implicit none
      integer nx, ny, nxv, iunit, nrec, lrec
      complex f
      character*(*) name
      dimension f(nxv,ny)
c local data
      integer j, k
c open new file
      if (nrec.lt.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='replace')
         nrec = 1
c open old file
      else if (nrec.eq.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         return
      endif
      write (unit=iunit,rec=nrec) ((f(j,k),j=1,nx),k=1,ny)
      nrec = nrec + 1
      return
      end
      subroutine FCREAD2(f,nx,ny,nxv,iunit,nrec,lrec,name,ierr)
c this subroutine reads complex 2d data f from a file
c f = input data to be read
c nx/ny = length of data f in x/y to read
c nxv = first dimension of data array f, must be >= nx
c iunit = fortran unit number
c nrec = current record number for read, if nrec >  0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nren <= 0)
c input: nx, ny, nxv, iunit, nrec, lrec, fname
c output: f, nrec, ierr
      implicit none
      integer nx, ny, nxv, iunit, nrec, lrec, ierr
      complex f
      character*(*) name
      dimension f(nxv,ny)
c local data
      integer j, k
      ierr = 0
      if (nrec.le.0) then
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='old')
         if (nrec.eq.0) return
         nrec = 1
      endif
      read (unit=iunit,rec=nrec,err=10) ((f(j,k),j=1,nx),k=1,ny)
      nrec = nrec + 1
      return
   10 ierr = 1
      return
      end
      subroutine FWRITE1(f,nxp,iunit,nrec,name)
c this subroutine write real data f to a direct access file
c f = input data to be written
c nxp = size of file f
c iunit = fortran unit number
c nrec = current record number for write (if negative, open file with
c recl=-nren)
c name = file name (used only if nren < 0)
c input: f, nxp, iunit, nrec, fname
c output: nrec
      implicit none
      integer nxp, iunit, nrec
      real f
      character*(*) name
      dimension f(nxp)
c local data
      integer lrec
      if (nrec.lt.0) then
         lrec = -nrec
         open(unit=iunit,file=name,form='unformatted',access='direct',re
     1cl=lrec,status='unknown')
         nrec = 1
      endif
      write (unit=iunit,rec=nrec) f
      nrec = nrec + 1
      return
      end
      function NDIAN()
c this function determines whether number format is big or little endian
c assumes ieee format for numbers
c ndian = (0,1) = architecture is (little-endian,big-endian)
      implicit none
      integer NDIAN
      integer*4 i
      real*8 d
      dimension i(2)
      equivalence (i,d)
      data d /1.0d0/
      save
c  big endian
      if (i(2).eq.0) then
         NDIAN = 1
c little endian
      else if (i(1).eq.0) then
         NDIAN = 0
c error
      else
         NDIAN = -1
      endif
      return
      end
      function NDPREC()
c this subroutine determines whether default reals are double precision
c ndprec = (0,1) = default reals are (normal,double-precision)
      implicit none
      integer NDPREC
      real small, prec, vresult
      data small /1.0e-12/
      save
      prec = 1.0 + small
      if (vresult(prec).gt.1.0) then
         NDPREC = 1
      else
         NDPREC = 0
      endif
      return
      end
      function IDPREC()
c this subroutine determines whether default integers are double
c precision
c idprec = (0,1) = default integers are (normal,double-precision)
      implicit none
      integer IDPREC
      integer ibig, iprec, iresult
      data ibig /2147483647/
      save
      iprec = ibig + 1
      if (iresult(iprec).gt.0.0) then
         IDPREC = 1
      else
         IDPREC = 0
      endif
      return
      end
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
      return
      end
      function iresult(iprec)
      implicit none
      integer iprec, iresult
      iresult = iprec
      return
      end
