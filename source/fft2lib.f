c-----------------------------------------------------------------------
c 2d PIC library for fast fourier transforms
c fft2lib.f contains procedures to perform ffts:
c FFT2C performs complex to complex fft and its inverse.
c WFFT2RINIT calculates tables needed by a two dimensional real to
c            complex fast fourier transform and its inverse.
c WFFT2RX performs real to complex fft and its inverse for scalar array,
c         with packed data.
c WFFT2R2 performs real to complex fft and its inverse for 2 component
c         vector array.
c WFFT2R3 performs real to complex fft and its inverse for 3 component
c         vector array.
c WFFT2RN performs real to complex fft and its inverse for n component
c         vector array.
c WFST2RINIT calculates tables needed by a two dimensional fast real
c            sine and cosine transforms and their inverses.
c WFSST2RX performs fast real sine/sine transform.
c WFSCT2RX performs fast real mixed sine/cosine transform.
c WFCST2RX performs fast real mixed cosine/sine transform.
c WFCCT2RX performs fast real cosine/cosine transform.
c WFCST2R2 performs fast real cosine/sine transform for 2 component
c          vector array for the electric field with dirichlet or
c          magnetic field with neumann boundary conditions.
c WFSCT2R2 performs fast real sine/cosine transform for 2 component
c          vector array for the magnetic field with dirichlet or
c          electric field with neumann boundary conditions.
c WFCST2R3 performs fast real cosine/sine transform for 3 component
c          vector array for the electric field with dirichlet or
c          magnetic field with neumann boundary conditions.
c WFSCT2R3 performs fast real sine/cosine transform for 3 component
c          vector array for the magnetic field with dirichlet or
c          electric field with neumann boundary conditions.
c WFSFT2RX performs fast real mixed sine/periodic transform.
c WFCFT2RX performs fast real mixed cosine/periodic transform.
c WFCSFT2R2 performs fast real cosine/sine/periodic transform for
c           2 component vector array for the electric field with
c           dirichlet or magnetic field with neumann boundary conditions
c WFSCFT2R2 performs fast real cosine/sine/periodic transform for
c           2 component vector array for the magnetic field with
c           dirichlet or electric field with neumann boundary conditions
c WFCSFT2R3 performs fast real cosine/sine/periodic transform for
c           3 component vector array for the electric field with
c           dirichlet or magnetic field with neumann boundary conditions
c WFSCFT2R3 performs fast real cosine/sine/periodic transform for
c           3 component vector array for the magnetic field with
c           dirichlet or electric field with neumann boundary conditions
c WFDT2RINIT calculates tables needed by a two dimensional fast real
c            sine DST-III/cosine DCT-III/periodic transforms and their
c            inverses.
c WFDSFT2RX performs fast real mixed sine DST-III/periodic transform.
c WFDCFT2RX performs fast real mixed cosine DCT-III/periodic transform
c WFDCSFT2R2 performs real cosine DCT-III/sine DST-III/periodic
c            transforms for the electric field with mixed
c            dirichlet-neumann or magnetic field with mixed
c            neumann-dirichlet boundary conditions.
c WFDCSFT2R2 performs real cosine DCT-III/sine DST-III/periodic
c            transforms for 2 component vector array for the electric
c            field with mixed dirichlet-neumann or magnetic field with
c            mixed neumann-dirichlet boundary conditions.
c WFDSCFT2R2 performs real sine DST-III/cosine DCT-III/periodic
c            transforms for 2 component vector array for the magnetic
c            field with mixed dirichlet-neumann or electric field with
c            mixed neumann-dirichlet boundary conditions.
c WFDCSFT2R3 performs real cosine DCT-III/sine DST-III/periodic
c            transforms for 3 component vector array for the electric
c            field with mixed dirichlet-neumann or magnetic field with
c            mixed neumann-dirichlet boundary conditions.
c WFDSCFT2R3 performs real sine DST-III/cosine DCT-III/periodic
c            transforms for 3 component vector array for the magnetic
c            field with mixed dirichlet-neumann or electric field with
c            mixed neumann-dirichlet boundary conditions.
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: february 19, 2010
c-----------------------------------------------------------------------
      subroutine FFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd
     1)
c this subroutine performs a two dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhd = first dimension of f
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      if (isign) 50, 10, 220
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
c bit-reverse array elements in x
   50 nrx = nxhy/nxh
      nry = nxhy/ny
      do 70 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 k = 1, ny
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 110 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 100 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 90 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 80 i = 1, ny
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 130 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 120 k = 1, ny
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(nxh2-j,k) = ani*conjg(t1 - t2)
  120 continue
  130 continue
      ani = 2.*ani
      do 140 k = 1, ny
      f(nxhh+1,k) = ani*conjg(f(nxhh+1,k))
      f(1,k) = ani*cmplx(real(f(1,k)) + aimag(f(1,k)),real(f(1,k)) - aim
     1ag(f(1,k)))
  140 continue
c bit-reverse array elements in y
      do 160 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 160
      do 150 j = 1, nxh
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
  150 continue
  160 continue
c then transform in y
      nry = nxy/ny
      do 200 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 170 i = 1, nxh
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  170 continue
  180 continue
  190 continue
  200 continue
c unscramble modes kx = 0, nx/2
      do 210 k = 2, nyh
      t1 = f(1,ny2-k)
      f(1,ny2-k) = .5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
      f(1,k) = .5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
  210 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  220 do 230 k = 2, nyh
      t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
      f(1,ny2-k) = conjg(f(1,k) - t1)
      f(1,k) = f(1,k) + t1
  230 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 250 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 j = 1, nxh
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
  240 continue
  250 continue
c first transform in y
      nry = nxy/ny
      do 290 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 280 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 270 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 260 i = 1, nxh
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  260 continue
  270 continue
  280 continue
  290 continue
c scramble coefficients
      kmr = nxy/nx
      do 310 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 300 k = 1, ny
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(nxh2-j,k) = conjg(t1 - t2)
  300 continue
  310 continue
      do 320 k = 1, ny
      f(nxhh+1,k) = 2.*conjg(f(nxhh+1,k))
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),real(f(1,k)) - aimag(f
     1(1,k)))
  320 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 340 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 340
      do 330 k = 1, ny
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
  330 continue
  340 continue
c then transform in x
      nrx = nxy/nxh
      do 380 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 370 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 360 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 350 i = 1, ny
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  350 continue
  360 continue
  370 continue
  380 continue
      return
      end
      subroutine FFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd
     1)
c this subroutine performs 2 two dimensional real to complex fast
c fourier transforms, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = 0, the fft tables are prepared
c if isign = -1, two inverse fourier transforms are performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      if (isign) 50, 10, 270
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
c swap complex components
   50 do 70 k = 1, ny
      do 60 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      nry = nxhy/ny
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 90
      do 80 k = 1, ny
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
   80 continue
   90 continue
c first transform in x
      nrx = nxy/nxh
      do 130 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 120 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 110 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 100 i = 1, ny
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 160 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 150 k = 1, ny
      do 140 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
  140 continue
  150 continue
  160 continue
      ani = 2.*ani
      do 180 k = 1, ny
      do 170 jj = 1, 2
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj
     1,1,k)) - aimag(f(jj,1,k)))
  170 continue
  180 continue
c bit-reverse array elements in y
      do 200 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 200
      do 190 j = 1, nxh
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  190 continue
  200 continue
c then transform in y
      nry = nxy/ny
      do 240 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 230 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 220 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 210 i = 1, nxh
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  210 continue
  220 continue
  230 continue
  240 continue
c unscramble modes kx = 0, nx/2
      do 260 k = 2, nyh
      do 250 jj = 1, 2
      t1 = f(jj,1,ny2-k)
      f(jj,1,ny2-k) = .5*cmplx(aimag(f(jj,1,k) + t1),real(f(jj,1,k) - t1
     1))
      f(jj,1,k) = .5*cmplx(real(f(jj,1,k) + t1),aimag(f(jj,1,k) - t1))
  250 continue
  260 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
cdir$ ivdep
  270 do 290 k = 2, nyh
      do 280 jj = 1, 2
      t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
      f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
      f(jj,1,k) = f(jj,1,k) + t1
  280 continue
  290 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 310 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 310
      do 300 j = 1, nxh
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  300 continue
  310 continue
c first transform in y
      nry = nxy/ny
      do 350 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 320 i = 1, nxh
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  320 continue
  330 continue
  340 continue
  350 continue
c scramble coefficients
      kmr = nxy/nx
      do 380 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 370 k = 1, ny
      do 360 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  360 continue
  370 continue
  380 continue
      do 400 k = 1, ny
      do 390 jj = 1, 2
      f(jj,nxhh+1,k) = 2.*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj,1,k
     1)) - aimag(f(jj,1,k)))
  390 continue
  400 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 420 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 420
      do 410 k = 1, ny
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  410 continue
  420 continue
c then transform in x
      nrx = nxy/nxh
      do 460 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 450 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 440 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 430 i = 1, ny
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
  430 continue
  440 continue
  450 continue
  460 continue
c swap complex components
      do 480 k = 1, ny
      do 470 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
  470 continue
  480 continue
      return
      end
      subroutine FFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyhd
     1)
c this subroutine performs 3 two dimensional real to complex fast
c fourier transforms, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = 0, the fft tables are prepared
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3, t4
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      if (isign) 50, 10, 270
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 20 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 40 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
c swap complex components
   50 do 70 i = 1, ny
      do 60 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),aimag(f(3,j,i)))
      at2 = aimag(f(2,j,i))
      f(2,j,i) = cmplx(aimag(f(1,j,i)),at1)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      nry = nxhy/ny
      do 90 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 90
      do 80 k = 1, ny
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   80 continue
   90 continue
c first transform in x
      nrx = nxy/nxh
      do 130 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 120 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 110 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 100 i = 1, ny
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
  100 continue
  110 continue
  120 continue
  130 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 160 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 150 k = 1, ny
      do 140 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
  140 continue
  150 continue
  160 continue
      ani = 2.*ani
      do 180 k = 1, ny
      do 170 jj = 1, 3
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj
     1,1,k)) - aimag(f(jj,1,k)))
  170 continue
  180 continue
c bit-reverse array elements in y
      do 200 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 200
      do 190 j = 1, nxh
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  190 continue
  200 continue
c then transform in y
      nry = nxy/ny
      do 240 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 230 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 220 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 210 i = 1, nxh
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  210 continue
  220 continue
  230 continue
  240 continue
c unscramble modes kx = 0, nx/2
      do 260 k = 2, nyh
      do 250 jj = 1, 3
      t1 = f(jj,1,ny2-k)
      f(jj,1,ny2-k) = .5*cmplx(aimag(f(jj,1,k) + t1),real(f(jj,1,k) - t1
     1))
      f(jj,1,k) = .5*cmplx(real(f(jj,1,k) + t1),aimag(f(jj,1,k) - t1))
  250 continue
  260 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
cdir$ ivdep
  270 do 290 k = 2, nyh
      do 280 jj = 1, 3
      t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
      f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
      f(jj,1,k) = f(jj,1,k) + t1
  280 continue
  290 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 310 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 310
      do 300 j = 1, nxh
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  300 continue
  310 continue
c first transform in y
      nry = nxy/ny
      do 350 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 320 i = 1, nxh
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  320 continue
  330 continue
  340 continue
  350 continue
c scramble coefficients
      kmr = nxy/nx
      do 380 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 370 k = 1, ny
      do 360 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  360 continue
  370 continue
  380 continue
      do 400 k = 1, ny
      do 390 jj = 1, 3
      f(jj,nxhh+1,k) = 2.*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj,1,k
     1)) - aimag(f(jj,1,k)))
  390 continue
  400 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 420 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 420
      do 410 k = 1, ny
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  410 continue
  420 continue
c then transform in x
      nrx = nxy/nxh
      do 460 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 450 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 440 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 430 i = 1, ny
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
  430 continue
  440 continue
  450 continue
  460 continue
c swap complex components
      do 480 i = 1, ny
      do 470 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(2,j,i)),aimag(f(3,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,aimag(f(1,j,i)))
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  470 continue
  480 continue
      return
      end
      subroutine FFT2C(f,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxyhd)
c this subroutine performs a two dimensional complex to complex fast
c fourier transform and its inverse, using complex arithmetic
c for isign = 0, input: all except f, output: mixup, sct
c for isign = (-1,1), input: all, output: f
c for isign = (-1,1), approximate flop count: 5*N*log2(N)
c where N = nx*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxd, nyd = first and second dimensions of f
c nxyd = maximum of (nx,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex f, sct, s, t
      dimension f(nxd,nyd), mixup(nxyd), sct(nxyhd)
      indxy = max0(indx,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nxy = 2**indxy
      if (isign.ne.0) go to 40
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxy
      lb = j - 1
      ll = 0
      do 10 k = 1, indxy
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
c bit-reverse array elements in x
   40 nrx = nxy/nx
      do 60 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 60
      do 50 i = 1, ny
      t = f(j1,i)
      f(j1,i) = f(j,i)
      f(j,i) = t
   50 continue
   60 continue
c bit-reverse array elements in y
      nry = nxy/ny
      do 80 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 80
      do 70 i = 1, nx
      t = f(i,k1)
      f(i,k1) = f(i,k)
      f(i,k) = t
   70 continue
   80 continue
      if (isign.gt.0) go to 190
c inverse fourier transform
c first transform in x
      do 120 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 110 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 100 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 90 i = 1, ny
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   90 continue
  100 continue
  110 continue
  120 continue
c then transform in y
      do 160 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 130 i = 1, nx
      t = s*f(i,j2)
      f(i,j2) = f(i,j1) - t
      f(i,j1) = f(i,j1) + t
  130 continue
  140 continue
  150 continue
  160 continue
c normalize result
      ani = 1./float(nx*ny)
      do 180 k = 1, ny
      do 170 j = 1, nx
      f(j,k) = f(j,k)*ani
  170 continue
  180 continue
      return
c forward fourier transform
c first transform in x
  190 do 230 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 200 i = 1, ny
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  200 continue
  210 continue
  220 continue
  230 continue
c then transform in y
      do 270 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 260 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 250 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 240 i = 1, nx
      t = s*f(i,j2)
      f(i,j2) = f(i,j1) - t
      f(i,j1) = f(i,j1) + t
  240 continue
  250 continue
  260 continue
  270 continue
      return
      end
      subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
c this subroutine calculates tables needed by a two dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, nxhyd, nxyhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyhd
      integer mixup
      complex sct
      dimension mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy, nxyh
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
      subroutine WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyh
     1d)
c wrapper function for real to complex fft, with packed data
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
      endif
      return
      end
      subroutine WFFT2RXX(f,isign,mixup,sct,indx,indy,nxh1d,nyd,nxhyd,nx
     1yhd)
c wrapper function for real to complex fft, with unpacked data
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxh1d, nyd, nxhyd, nxyhd
      dimension f(nxh1d,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh1, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh1 = 2**(indx - 1) + 1
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RXXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxh1d,nyd,nxhy
     1d,nxyhd)
c perform y fft
         call FFT2RXYX(f,isign,mixup,sct,indx,indy,nxi,nxh1,nxh1d,nyd,nx
     1hyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RXYX(f,isign,mixup,sct,indx,indy,nxi,nxh1,nxh1d,nyd,nx
     1hyd,nxyhd)
c perform x fft
         call FFT2RXXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxh1d,nyd,nxhy
     1d,nxyhd)
      endif
      return
      end
      subroutine WFFT2RXT(f,g,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,n
     1xyhd)
c wrapper function for real to complex fft
      implicit none
      complex f, g, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(nxhd,nyd), g(nyd,nxhd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2RXYT(f,g,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nx
     1hyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RXYT(f,g,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nx
     1hyd,nxyhd)
c perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
      endif
      return
      end
      subroutine WFFT2R2(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyh
     1d)
c wrapper function for 2 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
      endif
      return
      end
      subroutine WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,nxyh
     1d)
c wrapper function for 3 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
c perform y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
!      print*, 'fxF =', f(1,1:nxh,1)
!      print*, 'fxF =', f(1,1:nxh,2)
!      print*, 'fxF =', f(1,1:nxh,3)
!      print*, 'fxF =', f(1,1:nxh,4)
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd
     1,nxyhd)
c perform x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,
     1nxyhd)
      endif
      return
      end
      subroutine WFFT2RN(f,ss,isign,mixup,sct,indx,indy,nxhd,nyd,ndim,nx
     1hyd,nxyhd)
c wrapper function for 2d real to complex fft
      implicit none
      complex f, ss, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, ndim, nxhyd, nxyhd
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim,nxhd)
c local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2RNX(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,ndi
     1m,nxhyd,nxyhd)
c perform y fft
         call FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,
     1nxhyd,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,ndim,
     1nxhyd,nxyhd)
c perform x fft
         call FFT2RNX(f,ss,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,ndi
     1m,nxhyd,nxyhd)
      endif
      return
      end
      subroutine FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,nx
     1hyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic.
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
c first transform in x
      nrx = nxy/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = nyi, nyt
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(nxh2-j,k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.*ani
      do 90 k = nyi, nyt
      f(nxhh+1,k) = ani*conjg(f(nxhh+1,k))
      f(1,k) = ani*cmplx(real(f(1,k)) + aimag(f(1,k)),real(f(1,k)) - aim
     1ag(f(1,k)))
   90 continue
      return
c forward fourier transform
c scramble coefficients
  100 kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = nyi, nyt
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(nxh2-j,k) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 k = nyi, nyt
      f(nxhh+1,k) = 2.*conjg(f(nxhh+1,k))
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),real(f(1,k)) - aimag(f
     1(1,k)))
  130 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
  140 continue
  150 continue
c then transform in x
      nrx = nxy/nxh
      do 190 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 160 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
      subroutine FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,nx
     1hyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = first dimension of f >= nx/2
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 80
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      do 70 k = 2, nyh
      if (nxi.eq.1) then
         t1 = f(1,ny2-k)
         f(1,ny2-k) = .5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
         f(1,k) = .5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
      endif
   70 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   80 do 90 k = 2, nyh
      if (nxi.eq.1) then
         t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
         f(1,ny2-k) = conjg(f(1,k) - t1)
         f(1,k) = f(1,k) + t1
      endif
   90 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 110
      do 100 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
  100 continue
  110 continue
c first transform in y
      nry = nxy/ny
      do 150 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 120 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
      subroutine FFT2RXXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxh1d,nyd,
     1nxhyd,nxyhd)
c this subroutine performs the x part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic.  data is not packed
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2+1 and 1 <= k <= ny
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxh1d, nyd, nxhyd, nxyhd
      integer mixup
      complex f, sct
      dimension f(nxh1d,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer i, j, k, l, j1, j2, k1, k2, nrx, ns, ns2, km, kmr
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
c first transform in x
      nrx = nxy/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = nyi, nyt
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(nxh2-j,k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.*ani
      do 90 k = nyi, nyt
      f(nxhh+1,k) = ani*conjg(f(nxhh+1,k))
      f(nxh+1,k) = ani*cmplx(real(f(1,k)) - aimag(f(1,k)),0.0)
      f(1,k) = ani*cmplx(real(f(1,k)) + aimag(f(1,k)),0.0)
   90 continue
      return
c forward fourier transform
c scramble coefficients
  100 kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = nyi, nyt
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(nxh2-j,k) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 k = nyi, nyt
      f(nxhh+1,k) = 2.*conjg(f(nxhh+1,k))
      f(1,k) = cmplx(real(f(1,k)) + real(f(nxh+1,k)),real(f(1,k)) - real
     1(f(nxh+1,k)))
  130 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
  140 continue
  150 continue
c then transform in x
      nrx = nxy/nxh
      do 190 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 160 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
      subroutine FFT2RXYX(f,isign,mixup,sct,indx,indy,nxi,nxp,nxh1d,nyd,
     1nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic.  data is not packed
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxh1d = first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2+1 and 1 <= k <= ny
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxh1d, nyd, nxhyd, nxyhd
      integer mixup
      complex f, sct
      dimension f(nxh1d,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer i, j, k, l, j1, j2, k1, k2, nry, ns, ns2, km, kmr
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 70
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
      return
c forward fourier transform
c bit-reverse array elements in y
   70 nry = nxhy/ny
      do 90 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 90
      do 80 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
   80 continue
   90 continue
c first transform in y
      nry = nxy/ny
      do 130 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 120 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 110 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 100 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
      subroutine FFT2RXYT(f,g,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd
     1,nxhyd,nxyhd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic and transposing data for better cache usage
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c g = scratch array
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = first dimension of f
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, g, sct, t1, t2
      dimension f(nxhd,nyd), g(nyd,nxhd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
c transpose data
      do 20 k = 1, ny
      do 10 j = nxi, nxt
      g(k,j) = f(j,k)
   10 continue
   20 continue
      if (isign.gt.0) go to 90
c inverse fourier transform
      do 70 i = nxi, nxt
c bit-reverse array elements in y
      nry = nxhy/ny
      do 30 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 30
      t1 = g(k1,i)
      g(k1,i) = g(k,i)
      g(k,i) = t1
   30 continue
c then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      t2 = t1*g(j2,i)
      g(j2,i) = g(j1,i) - t2
      g(j1,i) = g(j1,i) + t2
   40 continue
   50 continue
   60 continue
   70 continue
c unscramble modes kx = 0, nx/2
      do 80 k = 2, nyh
      if (nxi.eq.1) then
         t1 = g(ny2-k,1)
         g(ny2-k,1) = .5*cmplx(aimag(g(k,1) + t1),real(g(k,1) - t1))
         g(k,1) = .5*cmplx(real(g(k,1) + t1),aimag(g(k,1) - t1))
      endif
   80 continue
      go to 160
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 do 100 k = 2, nyh
      if (nxi.eq.1) then
         t1 = cmplx(aimag(g(ny2-k,1)),real(g(ny2-k,1)))
         g(ny2-k,1) = conjg(g(k,1) - t1)
         g(k,1) = g(k,1) + t1
      endif
  100 continue
      do 150 i = nxi, nxt
c bit-reverse array elements in y
      nry = nxhy/ny
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 110
      t1 = g(k1,i)
      g(k1,i) = g(k,i)
      g(k,i) = t1
  110 continue
c first transform in y
      nry = nxy/ny
      do 140 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      t2 = t1*g(j2,i)
      g(j2,i) = g(j1,i) - t2
      g(j1,i) = g(j1,i) + t2
  120 continue
  130 continue
  140 continue
  150 continue
c transpose data back
  160 do 180 j = nxi, nxt
      do 170 k = 1, ny
      f(j,k) = g(k,j)
  170 continue
  180 continue
      return
      end
      subroutine FFT2R2X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,nx
     1hyd,nxyhd)
c this subroutine performs the x part of 2 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, two inverse fourier transforms are performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 140
c inverse fourier transform
c swap complex components
      do 20 k = nyi, nyt
      do 10 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
   10 continue
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = nyi, nyt
      do 90 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.*ani
      do 130 k = nyi, nyt
      do 120 jj = 1, 2
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj
     1,1,k)) - aimag(f(jj,1,k)))
  120 continue
  130 continue
      return
c forward fourier transform
c scramble coefficients
  140 kmr = nxy/nx
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = nyi, nyt
      do 150 jj = 1, 2
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, 2
      f(jj,nxhh+1,k) = 2.*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj,1,k
     1)) - aimag(f(jj,1,k)))
  180 continue
  190 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  200 continue
  210 continue
c then transform in x
      nrx = nxy/nxh
      do 250 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 k = nyi, nyt
      do 260 j = 1, nxh
      at1 = aimag(f(1,j,k))
      f(1,j,k) = cmplx(real(f(1,j,k)),real(f(2,j,k)))
      f(2,j,k) = cmplx(at1,aimag(f(2,j,k)))
  260 continue
  270 continue
      return
      end
      subroutine FFT2R2Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,nx
     1hyd,nxyhd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, two inverse fourier transforms are performed
c f(1:2,n,m) = (1/nx*ny)*sum(f(1:2,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:2,j,k) = sum(f(1:2,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3
      dimension f(2,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 90
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      do 80 k = 2, nyh
      if (nxi.eq.1) then
         do 70 jj = 1, 2
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = .5*cmplx(aimag(f(jj,1,k) + t1),real(f(jj,1,k) -
     1 t1))
         f(jj,1,k) = .5*cmplx(real(f(jj,1,k) + t1),aimag(f(jj,1,k) - t1)
     1)
   70    continue
      endif
   80 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 do 110 k = 2, nyh
      if (nxi.eq.1) then
         do 100 jj = 1, 2
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  100    continue
      endif
  110 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 130 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 130
      do 120 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  120 continue
  130 continue
c first transform in y
      nry = nxy/ny
      do 170 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 140 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
      subroutine FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,nx
     1hyd,nxyhd)
c this subroutine performs the x part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3, t4
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 140
c inverse fourier transform
c swap complex components
      do 20 i = nyi, nyt
      do 10 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),aimag(f(3,j,i)))
      at2 = aimag(f(2,j,i))
      f(2,j,i) = cmplx(aimag(f(1,j,i)),at1)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
   10 continue
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = nyi, nyt
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.*ani
      do 130 k = nyi, nyt
      do 120 jj = 1, 3
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj
     1,1,k)) - aimag(f(jj,1,k)))
  120 continue
  130 continue
      return
c forward fourier transform
c scramble coefficients
  140 kmr = nxy/nx
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = nyi, nyt
      do 150 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, 3
      f(jj,nxhh+1,k) = 2.*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj,1,k
     1)) - aimag(f(jj,1,k)))
  180 continue
  190 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  200 continue
  210 continue
c then transform in x
      nrx = nxy/nxh
      do 250 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
  220 continue
  230 continue
  240 continue
  250 continue
c swap complex components
      do 270 i = nyi, nyt
      do 260 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(2,j,i)),aimag(f(3,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,aimag(f(1,j,i)))
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  260 continue
  270 continue
      return
      end
      subroutine FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,nx
     1hyd,nxyhd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      complex f, sct, t1, t2, t3, t4
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 90
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble modes kx = 0, nx/2
      do 80 k = 2, nyh
      if (nxi.eq.1) then
         do 70 jj = 1, 3
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = .5*cmplx(aimag(f(jj,1,k) + t1),real(f(jj,1,k) -
     1 t1))
         f(jj,1,k) = .5*cmplx(real(f(jj,1,k) + t1),aimag(f(jj,1,k) - t1)
     1)
   70    continue
      endif
   80 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
   90 do 110 k = 2, nyh
      if (nxi.eq.1) then
         do 100 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  100    continue
      endif
  110 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 130 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 130
      do 120 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  120 continue
  130 continue
c first transform in y
      nry = nxy/ny
      do 170 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 140 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
      subroutine FFT2RNX(f,ss,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd
     1,ndim,nxhyd,nxyhd)
c this subroutine performs the x part of N two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c y, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
c for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
c where M = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:N,n,m) = (1/nx*ny)*sum(f(1:N,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:N,j,k) = sum(f(1:N,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c ss = scratch array
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c ndim = leading dimension of array f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, nyi, nyp, nxhd, nyd, ndim
      integer nxhyd, nxyhd
      complex f, ss, sct, t1, t2, t3
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
      dimension ss(ndim,nxhd)
c local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, ns, ns2, km, kmr, i, j, k, l, k1, k2, j1, j2, jj
      real ani
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 140
c inverse fourier transform
c swap complex components
      call SWAPC2N(f,ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 30
      do 20 k = nyi, nyt
      do 10 jj = 1, ndim
      t1 = f(jj,j1,k)
      f(jj,j1,k) = f(jj,j,k)
      f(jj,j,k) = t1
   10 continue
   20 continue
   30 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nyi, nyt
      do 40 jj = 1, ndim
      t2 = t1*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t2
      f(jj,j1,i) = f(jj,j1,i) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1./float(2*nx*ny)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = nyi, nyt
      do 90 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.*ani
      do 130 k = nyi, nyt
      do 120 jj = 1, ndim
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj
     1,1,k)) - aimag(f(jj,1,k)))
  120 continue
  130 continue
      return
c forward fourier transform
c scramble coefficients
  140 kmr = nxy/nx
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = nyi, nyt
      do 150 jj = 1, ndim
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, ndim
      f(jj,nxhh+1,k) = 2.*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),real(f(jj,1,k
     1)) - aimag(f(jj,1,k)))
  180 continue
  190 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 220 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 220
      do 210 k = nyi, nyt
      do 200 jj = 1, ndim
      t1 = f(jj,j1,k)
      f(jj,j1,k) = f(jj,j,k)
      f(jj,j,k) = t1
  200 continue
  210 continue
  220 continue
c then transform in x
      nrx = nxy/nxh
      do 270 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 260 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 250 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 240 i = nyi, nyt
      do 230 jj = 1, ndim
      t2 = t1*f(jj,j2,i)
      f(jj,j2,i) = f(jj,j1,i) - t2
      f(jj,j1,i) = f(jj,j1,i) + t2
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
c swap complex components
      call SWAPC2N(f,ss,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
      return
      end
      subroutine FFT2RNY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,nd
     1im,nxhyd,nxyhd)
c this subroutine performs the y part of N two dimensional real to
c complex fast fourier transforms, and their inverses, for a subset of
c x, using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: M*(5*log2(M) + 19/2)
c for isign = 1,  approximate flop count: M*(5*log2(M) + 15/2)
c where M = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, three inverse fourier transforms are performed
c f(1:N,n,m) = (1/nx*ny)*sum(f(1:N,j,k)*
c       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, two forward fourier transforms are performed
c f(1:N,j,k) = sum(f(1:N,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = second dimension of f
c nyd = third dimension of f
c ndim = leading dimension of array f
c nxhyd = maximum of (nx/2,ny)
c nxyhd = maximum of (nx,ny)/2
c fourier coefficients are stored as follows:
c f(2*j-1,k),f(2*j,k) = real, imaginary part of mode j-1,k-1, where
c 1 <= j <= nx/2 and 1 <= k <= ny, except for
c f(1,k),f(2,k) = real, imaginary part of mode nx/2,k-1, where
c ny/2+2 <= k <= ny, and
c f(2,1) = real part of mode nx/2,0 and
c f(2,ny/2+1) = real part of mode nx/2,ny/2
c written by viktor k. decyk, ucla
      implicit none
      integer isign, mixup, indx, indy, nxi, nxp, nxhd, nyd, ndim
      integer nxhyd, nxyhd
      complex f, sct, t1, t2
      dimension f(ndim,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
c local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, ns, ns2, km, kmr, i, j, k, l, k1, k2, j1, j2, jj
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 110
c inverse fourier transform
      nry = nxhy/ny
c bit-reverse array elements in y
      do 30 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 30
      do 20 j = nxi, nxt
      do 10 jj = 1, ndim
      t1 = f(jj,j,k1)
      f(jj,j,k1) = f(jj,j,k)
      f(jj,j,k) = t1
   10 continue
   20 continue
   30 continue
c then transform in y
      nry = nxy/ny
      do 80 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nxi, nxt
      do 40 jj = 1, ndim
      t2 = t1*f(jj,i,j2)
      f(jj,i,j2) = f(jj,i,j1) - t2
      f(jj,i,j1) = f(jj,i,j1) + t2
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 100 k = 2, nyh
      if (nxi.eq.1) then
         do 90 jj = 1, ndim
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = .5*cmplx(aimag(f(jj,1,k) + t1),real(f(jj,1,k) -
     1 t1))
         f(jj,1,k) = .5*cmplx(real(f(jj,1,k) + t1),aimag(f(jj,1,k) - t1)
     1)
   90    continue
      endif
  100 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  110 do 130 k = 2, nyh
      if (nxi.eq.1) then
         do 120 jj = 1, ndim
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  120    continue
      endif
  130 continue
c bit-reverse array elements in y
      nry = nxhy/ny
      do 160 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 160
      do 150 j = nxi, nxt
      do 140 jj = 1, ndim
      t1 = f(jj,j,k1)
      f(jj,j,k1) = f(jj,j,k)
      f(jj,j,k) = t1
  140 continue
  150 continue
  160 continue
c first transform in y
      nry = nxy/ny
      do 210 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 200 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 190 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 180 i = nxi, nxt
      do 170 jj = 1, ndim
      t2 = t1*f(jj,i,j2)
      f(jj,i,j2) = f(jj,i,j1) - t2
      f(jj,i,j1) = f(jj,i,j1) + t2
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
      return
      end
      subroutine SWAPC2N(f,s,isign,nxh,nyi,nyt,nxhd,nyd,ndim)
c this subroutine swaps components for multiple ffts
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nxh = complex dimension in x direction
c nyi/nyt = initial/final y index used
c nxhd = half of the second dimension of f
c nyd = third dimension of f
c ndim = leading dimension of array f
      implicit none
      integer isign, nxh, nyi, nyt, nxhd, nyd, ndim
      real f, s
      dimension f(ndim,2*nxhd,nyd), s(2*ndim*nxhd)
c local data
      integer i, j, k, ioff
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 60 k = nyi, nyt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1) = f(i,2*j-1,k)
         s(2*i+ioff) = f(i,2*j,k)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k) = s(i+ioff)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k) = s(i+ioff)
   40    continue
   50    continue
   60    continue
c complex to real
      else if (isign.gt.0) then
         do 120 k = nyi, nyt
         do 90 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 70 i = 1, ndim
         s(i+ioff) = f(i,2*j-1,k)
   70    continue
         ioff = ioff + ndim
         do 80 i = 1, ndim
         s(i+ioff) = f(i,2*j,k)
   80    continue
   90    continue
         do 110 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 100 i = 1, ndim
         f(i,2*j-1,k) = s(2*i+ioff-1)
         f(i,2*j,k) = s(2*i+ioff)
  100    continue
  110    continue
  120    continue
      endif
      return
      end
      subroutine WFFT2CINIT(mixup,sct,indx,indy,nxyd,nxyhd)
c this subroutine calculates tables needed by a two dimensional
c complex to complex fast fourier transform and its inverse.
c input: indx, indy, nxyd, nxyhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxyd = maximum of (nx,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxyd, nxyhd
      integer mixup
      complex sct
      dimension mixup(nxyd), sct(nxyhd)
c local data
      integer indxy, nxy, nxyh, j, k, lb, ll, jb, it
      real dnxy, arg
      indxy = max0(indx,indy)
      nxy = 2**indxy
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxy
      lb = j - 1
      ll = 0
      do 10 k = 1, indxy
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/float(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
      subroutine WFFT2C(f,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxyhd)
c wrapper function for complex to complex fft
      implicit none
      integer isign, indx, indy, nxd, nyd, nxyd, nxyhd
      integer mixup
      complex f, sct
      dimension f(nxd,nyd), mixup(nxyd), sct(nxyhd)
c local data
      integer nx, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2CX(f,isign,mixup,sct,indx,indy,nyi,ny,nxd,nyd,nxyd,nxy
     1hd)
c perform y fft
         call FFT2CY(f,isign,mixup,sct,indx,indy,nxi,nx,nxd,nyd,nxyd,nxy
     1hd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FFT2CY(f,isign,mixup,sct,indx,indy,nxi,nx,nxd,nyd,nxyd,nxy
     1hd)
c perform x fft
         call FFT2CX(f,isign,mixup,sct,indx,indy,nyi,ny,nxd,nyd,nxyd,nxy
     1hd)
      endif
      return
      end
      subroutine WFFT2CT(f,g,isign,mixup,sct,indx,indy,nxd,nyd,nxyd,nxyh
     1d)
c wrapper function for complex to complex fft
      implicit none
      integer isign, indx, indy, nxd, nyd, nxyd, nxyhd
      integer mixup
      complex f, g, sct
      dimension f(nxd,nyd), g(nyd,nxd), mixup(nxyd), sct(nxyhd)
c local data
      integer j, k, nx, ny, nxi, nyi
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call FFT2CX(f,isign,mixup,sct,indx,indy,nyi,ny,nxd,nyd,nxyd,nxy
     1hd)
c transpose data
         do 20 k = 1, ny
         do 10 j = 1, nx
         g(k,j) = f(j,k)
   10    continue
   20    continue
c perform y fft
         call FFT2CYT(g,isign,mixup,sct,indx,indy,nxi,nx,nxd,nyd,nxyd,nx
     1yhd)
c transpose data back
         do 40 j = 1, nx
         do 30 k = 1, ny
         f(j,k) = g(k,j)
   30    continue
   40    continue
c forward fourier transform
      else if (isign.gt.0) then
c transpose data
         do 60 k = 1, ny
         do 50 j = 1, nx
         g(k,j) = f(j,k)
   50    continue
   60    continue
c perform y fft
         call FFT2CYT(g,isign,mixup,sct,indx,indy,nxi,nx,nxd,nyd,nxyd,nx
     1yhd)
c transpose data back
         do 80 j = 1, nx
         do 70 k = 1, ny
         f(j,k) = g(k,j)
   70    continue
   80    continue
c perform x fft
         call FFT2CX(f,isign,mixup,sct,indx,indy,nyi,ny,nxd,nyd,nxyd,nxy
     1hd)
      endif
      return
      end
      subroutine WFFT2RC(f,f2c,g,isign,mixup,sct,indx,indy,nxv,nyv,nxd,n
     1yd,nxh1d,nxyd,nxyhd)
c wrapper function for real to complex fft
c for isign = -1, input is in real array f, output is in unpacked and
c transposed complex array g
c for isign = 1, input is in unpacked and transposed complex array g,
c output is in real array f
c f2c is a complex scratch array
c nxv >= nx, nyv >= ny, nxd = nx, nyd = ny, nxh1d = nx/2 + 1
      implicit none
      integer isign, indx, indy, nxv, nyv, nxd, nyd, nxh1d, nxyd, nxyhd
      integer mixup
      real f
      complex f2c, g, sct
      dimension f(nxv,nyv), f2c(nxd,nyd), g(nyd,nxh1d)
      dimension mixup(nxyd), sct(nxyhd)
c local data
      integer j, k, j1, nx, ny, nx2, nxh, nxh1, nxi, nyi
      real ani
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nxh1 = nxh + 1
      nx2 = nx + 2
c inverse fourier transform
      if (isign.lt.0) then
c copy real data to complex array
         do 20 k = 1, ny
         do 10 j = 1, nx
         f2c(j,k) = cmplx(f(j,k),0.0)
   10    continue
   20    continue
c perform x fft for some y
         call FFT2CXU(f2c,isign,mixup,sct,indx,indy,nyi,ny,nxd,nyd,nxyd,
     1nxyhd)
c transpose and normalize data
         ani = 1.0/real(nx*ny)
         do 40 k = 1, ny
         do 30 j = 1, nxh1
         g(k,j) = ani*f2c(j,k)
   30    continue
   40    continue
c perform y fft for some x
         call FFT2CYT(g,isign,mixup,sct,indx,indy,nxi,nxh1,nxh1d,nyd,nxy
     1d,nxyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft for some x
         call FFT2CYT(g,isign,mixup,sct,indx,indy,nxi,nxh1,nxh1d,nyd,nxy
     1d,nxyhd)
c transpose data
         do 60 j = 2, nxh
         j1 = nx2 - j
         do 50 k = 1, ny
         f2c(j,k) = g(k,j)
         f2c(j1,k) = conjg(g(k,j))
   50    continue
   60    continue
         do 70 k = 1, ny
         f2c(1,k) = g(k,1)
         f2c(nxh1,k) = g(k,nxh1)
   70    continue
c perform x fft for some y
         call FFT2CXU(f2c,isign,mixup,sct,indx,indy,nyi,ny,nxd,nyd,nxyd,
     1nxyhd)
c copy complex data to real array
         do 90 k = 1, ny
         do 80 j = 1, nx
         f(j,k) = real(f2c(j,k))
   80    continue
   90    continue
      endif
      return
      end
      subroutine WFFT2RCX(f,fc,g,isign,mixup,sct,indx,indy,nxv,nyv,nxhd,
     1nyd,nxh1d,nxhyd,nxhyhd)
c wrapper function for real to complex fft
c for isign = -1, input is in real array f, output is in unpacked and
c transposed complex array g
c for isign = 1, input is in unpacked and transposed complex array g,
c output is in real array f
c fc is a complex scratch array
c nxv >= nx, nyv >= ny, nxhd = nx/2, nyd = ny, nxh1d = nx/2 + 1
      implicit none
      integer isign, indx, indy, nxv, nyv, nxhd, nyd, nxh1d
      integer nxhyd, nxhyhd
      integer mixup
      real f
      complex fc, g, sct
      dimension f(nxv,nyv), fc(nxhd,nyd), g(nyd,nxh1d)
      dimension mixup(nxhyd), sct(nxhyhd)
c local data
      integer j, k, nx, ny, nx2, nxh, nxhh, nxh1, nxh2, nxi, nyi
      integer indx1
      real dnx, arg, ani
      complex t1, t2, t3
      data nxi, nyi /1,1/
c calculate range of indices
      indx1 = indx - 1
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nxhh = nx/4
      nxh1 = nxh + 1
      nxh2 = nxh + 2
      nx2 = nx + 2
      dnx = 6.28318530717959/real(nx)
c inverse fourier transform
      if (isign.lt.0) then
c copy real data to complex array
         do 20 k = 1, ny
         do 10 j = 1, nxh
         fc(j,k) = cmplx(f(2*j-1,k),f(2*j,k))
   10    continue
   20    continue
c perform x fft for some y
         call FFT2CXU(fc,isign,mixup,sct,indx1,indy,nyi,ny,nxhd,nyd,nxhy
     1d,nxhyhd)
c transpose data, unscramble coefficients, and normalize
         ani = 0.5/real(nx*ny)
         do 40 j = 2, nxhh
         arg = dnx*real(j-1)
         t3 = cmplx(sin(arg),cos(arg))
         do 30 k = 1, ny
         t2 = conjg(fc(nxh2-j,k))
         t1 = t2 + fc(j,k)
         t2 = (t2 - fc(j,k))*t3
         g(k,j) = ani*(t1 + t2)
         g(k,nxh2-j) = ani*conjg(t1 - t2)
   30    continue
   40    continue
         ani = 2.0*ani
         do 50 k = 1, ny
         g(k,nxhh+1) = ani*conjg(fc(nxhh+1,k))
         g(k,nxh1) = ani*cmplx(real(fc(1,k)) - aimag(fc(1,k)),0.0)
         g(k,1) = ani*cmplx(real(fc(1,k)) + aimag(fc(1,k)),0.0)
   50    continue
c perform y fft for some x
         call FFT2CYT(g,isign,mixup,sct,indx1,indy,nxi,nxh1,nxh1d,nyd,nx
     1hyd,nxhyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft for some x
         call FFT2CYT(g,isign,mixup,sct,indx1,indy,nxi,nxh1,nxh1d,nyd,nx
     1hyd,nxhyhd)
c transpose data and scramble coefficients
         do 70 j = 2, nxhh
         arg = dnx*real(j-1)
         t3 = cmplx(sin(arg),-cos(arg))
         do 60 k = 1, ny
         t2 = conjg(g(k,nxh2-j))
         t1 = t2 + g(k,j)
         t2 = (t2 - g(k,j))*t3
         fc(j,k) = t1 + t2
         fc(nxh2-j,k) = conjg(t1 - t2)
   60    continue
   70    continue
         do 80 k = 1, ny
         fc(nxhh+1,k) = 2.*conjg(g(k,nxhh+1))
         fc(1,k) = cmplx(real(g(k,1)) + real(g(k,nxh1)),real(g(k,1)) - r
     1eal(g(k,nxh1)))
   80    continue
c perform x fft for some y
         call FFT2CXU(fc,isign,mixup,sct,indx1,indy,nyi,ny,nxhd,nyd,nxhy
     1d,nxhyhd)
c copy complex data to real array
         do 100 k = 1, ny
         do 90 j = 1, nxh
         f(2*j-1,k) = real(fc(j,k))
         f(2*j,k) = aimag(fc(j,k))
   90    continue
  100    continue
      endif
      return
      end
      subroutine WFFT2RCN(f,fc,g,isign,mixup,sct,indx,indy,ndim,nxv,nyv,
     1nxhd,nyd,nxh1d,nxhyd,nxhyhd)
c wrapper function for multiple real to complex ffts
c for isign = -1, input is in real array f, output is in unpacked and
c transposed complex array g
c for isign = 1, input is in unpacked and transposed complex array g,
c output is in real array f
c fc is a complex scratch array
c ndim = leading dimension of array f
c nxv >= nx, nyv >= ny, nxhd = nx/2, nyd = ny, nxh1d = nx/2 + 1
      implicit none
      integer isign, indx, indy, ndim, nxv, nyv, nxhd, nyd, nxh1d
      integer nxhyd, nxhyhd
      integer mixup
      real f
      complex fc, g, sct
      dimension f(ndim,nxv,nyv), fc(nxhd,ndim,nyd), g(nyd,ndim,nxh1d)
      dimension mixup(nxhyd), sct(nxhyhd)
c local data
      integer i, j, k, nx, ny, nx2, nxh, nxhh, nxh1, nxh2, nxi, nyi
      integer indx1, nxp, nyp, nnxd, nnyd
      real dnx, arg, ani
      complex t1, t2, t3
      data nxi, nyi /1,1/
c calculate range of indices
      indx1 = indx - 1
      nx = 2**indx
      ny = 2**indy
      nxh = nx/2
      nxhh = nx/4
      nxh1 = nxh + 1
      nxh2 = nxh + 2
      nx2 = nx + 2
      nxp = ndim*nxh1
      nyp = ndim*ny
      nnxd = ndim*nxh1d
      nnyd = ndim*nyd
      dnx = 6.28318530717959/real(nx)
c inverse fourier transform
      if (isign.lt.0) then
c copy real data to complex array
         do 30 k = 1, ny
         do 20 j = 1, nxh
         do 10 i = 1, ndim
         fc(j,i,k) = cmplx(f(i,2*j-1,k),f(i,2*j,k))
   10    continue
   20    continue
   30    continue
c perform x fft for some y
         call FFT2CXU(fc,isign,mixup,sct,indx1,indy,nyi,nyp,nxhd,nnyd,nx
     1hyd,nxhyhd)
c transpose data, unscramble coefficients, and normalize
         ani = 0.5/real(nx*ny)
         do 60 j = 2, nxhh
         arg = dnx*real(j-1)
         t3 = cmplx(sin(arg),cos(arg))
         do 50 k = 1, ny
         do 40 i = 1, ndim
         t2 = conjg(fc(nxh2-j,i,k))
         t1 = t2 + fc(j,i,k)
         t2 = (t2 - fc(j,i,k))*t3
         g(k,i,j) = ani*(t1 + t2)
         g(k,i,nxh2-j) = ani*conjg(t1 - t2)
   40    continue
   50    continue
   60    continue
         ani = 2.0*ani
         do 80 k = 1, ny
         do 70 i = 1, ndim
         g(k,i,nxhh+1) = ani*conjg(fc(nxhh+1,i,k))
         g(k,i,nxh1) = ani*cmplx(real(fc(1,i,k)) - aimag(fc(1,i,k)),0.0)
         g(k,i,1) = ani*cmplx(real(fc(1,i,k)) + aimag(fc(1,i,k)),0.0)
   70    continue
   80    continue
c perform y fft for some x
         call FFT2CYT(g,isign,mixup,sct,indx1,indy,nxi,nxp,nnxd,nyd,nxhy
     1d,nxhyhd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft for some x
         call FFT2CYT(g,isign,mixup,sct,indx1,indy,nxi,nxp,nnxd,nyd,nxhy
     1d,nxhyhd)
c transpose data and scramble coefficients
         do 110 j = 2, nxhh
         arg = dnx*real(j-1)
         t3 = cmplx(sin(arg),-cos(arg))
         do 100 i = 1, ndim
         do 90 k = 1, ny
         t2 = conjg(g(k,i,nxh2-j))
         t1 = t2 + g(k,i,j)
         t2 = (t2 - g(k,i,j))*t3
         fc(j,i,k) = t1 + t2
         fc(nxh2-j,i,k) = conjg(t1 - t2)
   90    continue
  100    continue
  110    continue
         do 130 i = 1, ndim
         do 120 k = 1, ny
         fc(nxhh+1,i,k) = 2.*conjg(g(k,i,nxhh+1))
         fc(1,i,k) = cmplx(real(g(k,i,1)) + real(g(k,i,nxh1)),real(g(k,i
     1,1)) - real(g(k,i,nxh1)))
  120    continue
  130    continue
c perform x fft for some y
         call FFT2CXU(fc,isign,mixup,sct,indx1,indy,nyi,nyp,nxhd,nnyd,nx
     1hyd,nxhyhd)
c copy complex data to real array
         do 160 k = 1, ny
         do 150 i = 1, ndim
         do 140 j = 1, nxh
         f(i,2*j-1,k) = real(fc(j,i,k))
         f(i,2*j,k) = aimag(fc(j,i,k))
  140    continue
  150    continue
  160    continue
      endif
      return
      end
      subroutine FFT2CX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxd,nyd,nxyd
     1,nxyhd)
c this subroutine performs the x part of a two dimensional complex to
c complex fast fourier transform and its inverse,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = (-1,1), approximate flop count: 5*N*log2(N)
c where N = nx*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = (1/nx*ny)*sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxd, nyd = first and second dimensions of f
c nxyd = maximum of (nx,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxd, nyd, nxyd, nxyhd
      integer mixup
      complex f, sct
      dimension f(nxd,nyd), mixup(nxyd), sct(nxyhd)
      integer indxy, nx, nxh, ny, nxy, nyt, nrx, ns, ns2, km, kmr
      integer i, j, k, l, j1, j2, k1, k2
      real ani
      complex s, t
      if (isign.eq.0) return
      indxy = max0(indx,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nxy = 2**indxy
      nyt = nyi + nyp - 1
c bit-reverse array elements in x
      nrx = nxy/nx
      do 20 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 i = nyi, nyt
      t = f(j1,i)
      f(j1,i) = f(j,i)
      f(j,i) = t
   10 continue
   20 continue
      if (isign.gt.0) go to 90
c inverse fourier transform
c transform in x
      do 60 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = nyi, nyt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
c normalize result
      ani = 1./float(nx*ny)
      do 80 k = nyi, nyt
      do 70 j = 1, nx
      f(j,k) = f(j,k)*ani
   70 continue
   80 continue
      return
c forward fourier transform
c transform in x
   90 do 130 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 120 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 110 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 100 i = nyi, nyt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
      subroutine FFT2CY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxd,nyd,nxyd
     1,nxyhd)
c this subroutine performs the y part of a two dimensional complex to
c complex fast fourier transform and its inverse,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = (-1,1), approximate flop count: 5*N*log2(N)
c where N = nx*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd, nyd = first and second dimensions of f
c nxyd = maximum of (nx,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxd, nyd, nxyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxd,nyd), mixup(nxyd), sct(nxyhd)
      integer indxy, nx, ny, nyh, nxy, nxt, nry, ns, ns2, km, kmr
      integer i, j, k, l, j1, j2, k1, k2
      complex s, t
      if (isign.eq.0) return
      indxy = max0(indx,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxy = 2**indxy
      nxt = nxi + nxp - 1
c bit-reverse array elements in y
      nry = nxy/ny
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 i = nxi, nxt
      t = f(i,k1)
      f(i,k1) = f(i,k)
      f(i,k) = t
   10 continue
   20 continue
      if (isign.gt.0) go to 70
c inverse fourier transform
c transform in y
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t = s*f(i,j2)
      f(i,j2) = f(i,j1) - t
      f(i,j1) = f(i,j1) + t
   30 continue
   40 continue
   50 continue
   60 continue
      return
c forward fourier transform
c transform in y
   70 do 110 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 100 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 90 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 80 i = nxi, nxt
      t = s*f(i,j2)
      f(i,j2) = f(i,j1) - t
      f(i,j1) = f(i,j1) + t
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
      subroutine FFT2CYT(g,isign,mixup,sct,indx,indy,nxi,nxp,nxd,nyd,nxy
     1d,nxyhd)
c this subroutine performs the y part of a two dimensional complex to
c complex fast fourier transform and its inverse,
c using complex arithmetic
c for isign = (-1,1), input: all, output: g
c for isign = (-1,1), approximate flop count: 5*N*log2(N)
c where N = nx*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c g(m,n) = (sum(g(k,j)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c g(k,j) = sum(g(m,n)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd, nyd = first and second dimensions of f
c nxyd = maximum of (nx,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxd, nyd, nxyd, nxyhd
      complex g, sct
      integer mixup
      dimension g(nyd,nxd), mixup(nxyd), sct(nxyhd)
      integer indxy, nx, ny, nyh, nxy, nxt, nry, ns, ns2, km, kmr
      integer i, j, k, l, j1, j2, k1, k2
      complex s, t
      if (isign.eq.0) return
      indxy = max0(indx,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxy = 2**indxy
      nxt = nxi + nxp - 1
c bit-reverse array elements in y
      nry = nxy/ny
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 i = nxi, nxt
      t = g(k1,i)
      g(k1,i) = g(k,i)
      g(k,i) = t
   10 continue
   20 continue
      if (isign.gt.0) go to 70
c inverse fourier transform
c transform in y
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
      return
c forward fourier transform
c transform in y
   70 do 110 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 100 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 90 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 80 i = nxi, nxt
      t = s*g(j2,i)
      g(j2,i) = g(j1,i) - t
      g(j1,i) = g(j1,i) + t
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
      subroutine FFT2CXU(f,isign,mixup,sct,indx,indy,nyi,nyp,nxd,nyd,nxy
     1d,nxyhd)
c this subroutine performs the x part of a two dimensional complex to
c complex fast fourier transform and its inverse,
c using complex arithmetic.  The inverse fft does not normalize
c for isign = (-1,1), input: all, output: f
c for isign = (-1,1), approximate flop count: 5*N*log2(N)
c where N = nx*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(n,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxd, nyd = first and second dimensions of f
c nxyd = maximum of (nx,ny)
c nxyhd = one half of maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxd, nyd, nxyd, nxyhd
      integer mixup
      complex f, sct
      dimension f(nxd,nyd), mixup(nxyd), sct(nxyhd)
      integer indxy, nx, nxh, ny, nxy, nyt, nrx, ns, ns2, km, kmr
      integer i, j, k, l, j1, j2, k1, k2
      complex s, t
      if (isign.eq.0) return
      indxy = max0(indx,indy)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nxy = 2**indxy
      nyt = nyi + nyp - 1
c bit-reverse array elements in x
      nrx = nxy/nx
      do 20 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 i = nyi, nyt
      t = f(j1,i)
      f(j1,i) = f(j,i)
      f(j,i) = t
   10 continue
   20 continue
      if (isign.gt.0) go to 70
c inverse fourier transform
c transform in x
      do 60 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 30 i = nyi, nyt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   30 continue
   40 continue
   50 continue
   60 continue
      return
c forward fourier transform
c transform in x
   70 do 110 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 100 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 90 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 80 i = nyi, nyt
      t = s*f(j2,i)
      f(j2,i) = f(j1,i) - t
      f(j1,i) = f(j1,i) + t
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
      subroutine WFST2RINIT(mixup,sctd,indx,indy,nxhyd,nxyd)
c this subroutine calculates tables needed by a two dimensional
c fast real sine and cosine transforms and their inverses.
c input: indx, indy, nxhyd, nxyd
c output: mixup, sctd
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyd
      integer mixup
      complex sctd
      dimension mixup(nxhyd), sctd(nxyd)
c local data
      integer indx1, indx1y, nx, ny, nxy, nxhy
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles n*pi/nxy
      dnxy = 0.5*6.28318530717959/float(nxy)
      do 30 j = 1, nxy
      arg = dnxy*float(j - 1)
      sctd(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
      subroutine WFSST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for real sine/sine transform
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c perform y sine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y sine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
      subroutine WFSCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for real mixed sine/cosine transform
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c perform y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
      subroutine WFCST2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for real mixed cosine/sine transform
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c perform y cosine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y cosine transform
         call FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
      subroutine WFCCT2RX(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for real cosine/cosine transform
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c perform y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y cosine transform
         call FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxhy
     1d,nxyd)
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
      subroutine FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,n
     1xhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c sine transform and its inverse, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform is performed
c f(n,k) = (1/nx*ny)*sum(f(j,k)*sin(pi*n*j/nx))
c if isign = 1, a forward sine transform is performed
c f(j,k) = sum(f(n,k)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f => nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at1 = -aimag(sctd(j1))*at1
      at2 = .5*at2
      f(j,k) = at1 + at2
      f(nx+2-j,k) = at1 - at2
   10 continue
      f(1,k) = 0.0
      f(nxh+1,k) = 2.0*f(nxh+1,k)
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 100 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 90 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = ani*(t2 + t4)
         f(2*j,k) = ani*(t3 + t5)
         f(nx3-2*j,k) = ani*(t2 - t4)
         f(nx3-2*j+1,k) = ani*(t5 - t3)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 k = nyi, nyt
         f(nxh+1,k) = ani*f(nxh+1,k)
         f(nxh+2,k) = -ani*f(nxh+2,k)
         t2 = ani*(f(1,k) + f(2,k))
         f(2,k) = ani*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = ani*f(nx+1,k)
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = t2 + t4
         f(2*j,k) = t3 + t5
         f(nx3-2*j,k) = t2 - t4
         f(nx3-2*j+1,k) = t5 - t3
  120    continue
  130    continue
         do 140 k = nyi, nyt
         f(nxh+1,k) = 2.0*f(nxh+1,k)
         f(nxh+2,k) = -2.0*f(nxh+2,k)
         t2 = 2.0*(f(1,k) + f(2,k))
         f(2,k) = 2.0*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = 2.0*f(nx+1,k)
  140    continue
      endif
c perform recursion for sine transform
      do 160 k = nyi, nyt
      sum1 = .5*f(1,k)
      f(1,k) = 0.0
      f(2,k) = sum1
      do 150 j = 2, nxh
      sum1 = sum1 + f(2*j-1,k)
      f(2*j-1,k) = -f(2*j,k)
      f(2*j,k) = sum1
  150 continue
      f(nx+1,k) = 0.0
  160 continue
      return
      end
      subroutine FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,n
     1xhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c cosine transform and its inverse, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform is performed
c f(n,k) = (1/nx*ny)*(.5*f(1,k) + ((-1)**n)*f(nx+1,k) + sum(f(j,k)*
c       cos(pi*n*j/nx)))
c if isign = 1, a forward cosine transform is performed
c f(j,k) = 2*(.5*f(1,k) + ((-1)**j)*f(n+1,k) + sum(f(n,k)*
c       cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(1,k) - f(nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at2 = f(nx+2-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = -aimag(sctd(j1))*at2
      at1 = .5*at1
      f(j,k) = at1 - at2
      f(nx+2-j,k) = at1 + at2
   10 continue
      f(1,k) = .5*(f(1,k) + f(nx+1,k))
      f(nx+1,k) = sum1
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 100 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 90 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = ani*(t2 + t4)
         f(2*j,k) = ani*(t3 + t5)
         f(nx3-2*j,k) = ani*(t2 - t4)
         f(nx3-2*j+1,k) = ani*(t5 - t3)
   90    continue
  100    continue
         ani = 2.*ani
         do 110 k = nyi, nyt
         f(nxh+1,k) = ani*f(nxh+1,k)
         f(nxh+2,k) = -ani*f(nxh+2,k)
         t2 = ani*(f(1,k) + f(2,k))
         f(2,k) = ani*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = ani*f(nx+1,k)
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         t4 = f(nx3-2*j,k)
         t5 = -f(nx3-2*j+1,k)
         t2 = f(2*j-1,k) + t4
         t3 = f(2*j,k) + t5
         t6 = f(2*j-1,k) - t4
         t5 = f(2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(2*j-1,k) = t2 + t4
         f(2*j,k) = t3 + t5
         f(nx3-2*j,k) = t2 - t4
         f(nx3-2*j+1,k) = t5 - t3
  120    continue
  130    continue
         do 140 k = nyi, nyt
         f(nxh+1,k) = 2.0*f(nxh+1,k)
         f(nxh+2,k) = -2.0*f(nxh+2,k)
         t2 = 2.0*(f(1,k) + f(2,k))
         f(2,k) = 2.0*(f(1,k) - f(2,k))
         f(1,k) = t2
         f(nx+1,k) = 2.0*f(nx+1,k)
  140    continue
      endif
c perform recursion for cosine transform
      do 160 k = nyi, nyt
      sum1 = f(nx+1,k)
      f(nx+1,k) = f(2,k)
      f(2,k) = sum1
      do 150 j = 2, nxh
      sum1 = sum1 - f(2*j,k)
      f(2*j,k) = sum1
  150 continue
  160 continue
      return
      end
      subroutine FST2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxhd,nyd,n
     1xhyd,nxyd)
c this subroutine performs the y part of a two dimensional fast real
c sine transform and its inverse, for a subset of x,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 19)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform is performed
c f(n,m) = sum(f(n,k)*sin(pi*m*k/ny))
c if isign = 1, a forward sine transform is performed
c f(n,k) = 2*sum(f(n,m)*sin(pi*m*k/ny)
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny + 1
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      do 20 j = nxi, nxt
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = f(j,ny2-k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at1 = -aimag(sctd(k1))*at1
      at2 = .5*at2
      f(j,k) = at1 + at2
      f(j,ny2-k) = at1 - at2
   10 continue
      f(j,1) = 0.0
      f(j,nyh+1) = 2.0*f(j,nyh+1)
   20 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 40 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 40
      do 30 j = nxi, nxt
      t2 = f(j,2*k1-1)
      t3 = f(j,2*k1)
      f(j,2*k1-1) = f(j,2*k-1)
      f(j,2*k1) = f(j,2*k)
      f(j,2*k-1) = t2
      f(j,2*k) = t3
   30 continue
   40 continue
c then transform in y
      nry = nxy/nyh
      do 80 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 50 i = nxi, nxt
      t2 = real(t1)*f(i,2*j2-1) - aimag(t1)*f(i,2*j2)
      t3 = aimag(t1)*f(i,2*j2-1) + real(t1)*f(i,2*j2)
      f(i,2*j2-1) = f(i,2*j1-1) - t2
      f(i,2*j2) = f(i,2*j1) - t3
      f(i,2*j1-1) = f(i,2*j1-1) + t2
      f(i,2*j1) = f(i,2*j1) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         do 100 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 90 j = nxi, nxt
         t4 = f(j,ny+3-2*k)
         t5 = -f(j,ny+3-2*k+1)
         t2 = f(j,2*k-1) + t4
         t3 = f(j,2*k) + t5
         t6 = f(j,2*k-1) - t4
         t5 = f(j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(j,2*k-1) = ani*(t2 + t4)
         f(j,2*k) = ani*(t3 + t5)
         f(j,ny+3-2*k) = ani*(t2 - t4)
         f(j,ny+3-2*k+1) = ani*(t5 - t3)
   90    continue
  100    continue
         do 110 j = nxi, nxt
         f(j,nyh+2) = -f(j,nyh+2)
         t2 = f(j,1) + f(j,2)
         f(j,2) = f(j,1) - f(j,2)
         f(j,1) = t2
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = nxi, nxt
         t4 = f(j,ny+3-2*k)
         t5 = -f(j,ny+3-2*k+1)
         t2 = f(j,2*k-1) + t4
         t3 = f(j,2*k) + t5
         t6 = f(j,2*k-1) - t4
         t5 = f(j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(j,2*k-1) = t2 + t4
         f(j,2*k) = t3 + t5
         f(j,ny+3-2*k) = t2 - t4
         f(j,ny+3-2*k+1) = t5 - t3
  120    continue
  130    continue
         do 140 j = nxi, nxt
         f(j,nyh+1) = 2.0*f(j,nyh+1)
         f(j,nyh+2) = -2.0*f(j,nyh+2)
         t2 = 2.0*(f(j,1) + f(j,2))
         f(j,2) = 2.0*(f(j,1) - f(j,2))
         f(j,1) = t2
         f(j,ny+1) = 2.0*f(j,ny+1)
  140    continue
      endif
c perform recursion for sine transform
      do 160 j = nxi, nxt
      sum1 = .5*f(j,1)
      f(j,1) = 0.
      f(j,2) = sum1
      do 150 k = 2, nyh
      sum1 = sum1 + f(j,2*k-1)
      f(j,2*k-1) = -f(j,2*k)
      f(j,2*k) = sum1
  150 continue
      f(j,ny+1) = 0.0
  160 continue
      return
      end
      subroutine FCT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxhd,nyd,n
     1xhyd,nxyd)
c this subroutine performs the y part of a two dimensional fast real
c cosine transform and its inverse, for a subset of x,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 20)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform is performed
c f(n,m) = (.5*f(n,1) + ((-1)**m)*f(n,ny+1) + sum(f(n,k)*cos(pi*m*k/ny))
c if isign = 1, a forward cosine transform is performed
c f(n,k) = 2*(.5*f(n,1) + ((-1)**m)*f(n,ny+1) + sum(f(n,m)*
c       cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = half of the first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny + 1
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      do 20 j = nxi, nxt
      sum1 = .5*(f(j,1) - f(j,ny+1))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at2 = f(j,ny2-k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = -aimag(sctd(k1))*at2
      at1 = .5*at1
      f(j,k) = at1 - at2
      f(j,ny2-k) = at1 + at2
   10 continue
      f(j,1) = .5*(f(j,1) + f(j,ny+1))
      f(j,ny+1) = sum1
   20 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 40 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 40
      do 30 j = nxi, nxt
      t2 = f(j,2*k1-1)
      t3 = f(j,2*k1)
      f(j,2*k1-1) = f(j,2*k-1)
      f(j,2*k1) = f(j,2*k)
      f(j,2*k-1) = t2
      f(j,2*k) = t3
   30 continue
   40 continue
c then transform in y
      nry = nxy/nyh
      do 80 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 50 i = nxi, nxt
      t2 = real(t1)*f(i,2*j2-1) - aimag(t1)*f(i,2*j2)
      t3 = aimag(t1)*f(i,2*j2-1) + real(t1)*f(i,2*j2)
      f(i,2*j2-1) = f(i,2*j1-1) - t2
      f(i,2*j2) = f(i,2*j1) - t3
      f(i,2*j1-1) = f(i,2*j1-1) + t2
      f(i,2*j1) = f(i,2*j1) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         do 100 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 90 j = nxi, nxt
         t4 = f(j,ny+3-2*k)
         t5 = -f(j,ny+3-2*k+1)
         t2 = f(j,2*k-1) + t4
         t3 = f(j,2*k) + t5
         t6 = f(j,2*k-1) - t4
         t5 = f(j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(j,2*k-1) = ani*(t2 + t4)
         f(j,2*k) = ani*(t3 + t5)
         f(j,ny+3-2*k) = ani*(t2 - t4)
         f(j,ny+3-2*k+1) = ani*(t5 - t3)
   90    continue
  100    continue
         do 110 j = nxi, nxt
         f(j,nyh+2) = -f(j,nyh+2)
         t2 = f(j,1) + f(j,2)
         f(j,2) = f(j,1) - f(j,2)
         f(j,1) = t2
  110    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = nxi, nxt
         t4 = f(j,ny+3-2*k)
         t5 = -f(j,ny+3-2*k+1)
         t2 = f(j,2*k-1) + t4
         t3 = f(j,2*k) + t5
         t6 = f(j,2*k-1) - t4
         t5 = f(j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(j,2*k-1) = t2 + t4
         f(j,2*k) = t3 + t5
         f(j,ny+3-2*k) = t2 - t4
         f(j,ny+3-2*k+1) = t5 - t3
  120    continue
  130    continue
         do 140 j = nxi, nxt
         f(j,nyh+1) = 2.0*f(j,nyh+1)
         f(j,nyh+2) = -2.0*f(j,nyh+2)
         t2 = 2.0*(f(j,1) + f(j,2))
         f(j,2) = 2.0*(f(j,1) - f(j,2))
         f(j,1) = t2
         f(j,ny+1) = 2.0*f(j,ny+1)
  140    continue
      endif
c perform recursion for cosine transform
      do 160 j = nxi, nxt
      sum1 = f(j,ny+1)
      f(j,ny+1) = f(j,2)
      f(j,2) = sum1
      do 150 k = 2, nyh
      sum1 = sum1 - f(j,2*k)
      f(j,2*k) = sum1
  150 continue
  160 continue
      return
      end
      subroutine WFCST2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for 2 real cosine-sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
c perform y sine-cosine transforms
         call FSCT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y sine-cosine transforms
         call FSCT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxh
     1yd,nxyd)
c perform x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
      endif
      return
      end
      subroutine WFSCT2R2(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for 2 real sine-cosine transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
c perform y cosine-sine transforms
         call FCST2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y cosine-sine transforms
         call FCST2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nxh
     1yd,nxyd)
c perform x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
      endif
      return
      end
      subroutine FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,
     1nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
c             + sum(f(1,j,k)*cos(pi*n*j/nx)))
c f(2,n,k) = (1/nx*ny)*sum(f(2,j,k)*sin(pi*n*j/nx))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k) + sum(f(1,n,k)*
c       cos(pi*n*j/nx))
c f(2,j,k) = sum(f(2,n,k)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(1,1,k) - f(1,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,nx+2-j,k) = at1 + at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,nx+2-j,k) = at1 - at2
   10 continue
      f(1,1,k) = .5*(f(1,1,k) + f(1,nx+1,k))
      f(1,nx+1,k) = sum1
      f(2,1,k) = 0.0
      f(2,nxh+1,k) = 2.0*f(2,nxh+1,k)   
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 2
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 2
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine-sine transform
      do 220 k = nyi, nyt
      sum1 = f(1,nx+1,k)
      f(1,nx+1,k) = f(1,2,k)
      f(1,2,k) = sum1
      sum2 = .5*f(2,1,k)
      f(2,1,k) = 0.0
      f(2,2,k) = sum2
      do 210 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 + f(2,2*j-1,k)
      f(2,2*j-1,k) = -f(2,2*j,k)
      f(2,2*j,k) = sum2
  210 continue
      f(2,nx+1,k) = 0.0
  220 continue
      return
      end
      subroutine FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd,
     1nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
c f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + ((-1)**n)*f(2,nx+1,k)
c             + sum(f(2,j,k)*cos(pi*n*j/nx)))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
c f(2,j,k) = 2*(.5*f(2,1,k) + ((-1)**j)*f(2,n+1,k)
c             + sum(f(2,n,k)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(2,1,k) - f(2,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,nx+2-j,k) = at1 - at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,nx+2-j,k) = at1 + at2
   10 continue
      f(1,1,k) = 0.0
      f(1,nxh+1,k) = 2.0*f(1,nxh+1,k)
      f(2,1,k) = .5*(f(2,1,k) + f(2,nx+1,k))
      f(2,nx+1,k) = sum1
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 2
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 2
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 2
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for sine-cosine transform
      do 220 k = nyi, nyt
      sum1 = .5*f(1,1,k)
      f(1,1,k) = 0.0
      f(1,2,k) = sum1
      sum2 = f(2,nx+1,k)
      f(2,nx+1,k) = f(2,2,k)
      f(2,2,k) = sum2
      do 210 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k)
      f(1,2*j-1,k) = -f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 - f(2,2*j,k)
      f(2,2*j,k) = sum2
  210 continue
      f(1,nx+1,k) = 0.0
  220 continue
      return
      end
      subroutine FSCT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxhd,nyd,
     1nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 20)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,m) = sum(f(1,n,k)*sin(pi*m*k/ny))
c f(2,n,m) = (.5*f(2,n,1) + ((-1)**m)*f(2,n,ny+1) + sum(f(2,n,k)*
c            cos(pi*m*k/ny))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,n,k) = 2*sum(f(1,n,m)*sin(pi*m*k/ny)
c f(2,n,k) = 2*(.5*f(2,n,1) + ((-1)**m)*f(2,n,ny+1) + sum(f(2,n,m)*
c       cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = half of the first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny + 1
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      do 20 j = nxi, nxt
      sum1 = .5*(f(2,j,1) - f(2,j,ny+1))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = f(1,j,ny2-k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,j,ny2-k) = at1 - at2
      at2 = f(2,j,ny2-k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,j,ny2-k) = at1 + at2
   10 continue
      f(1,j,1) = 0.0
      f(1,j,nyh+1) = 2.0*f(1,j,nyh+1)
      f(2,j,1) = .5*(f(2,j,1) + f(2,j,ny+1))
      f(2,j,ny+1) = sum1
   20 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 50 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = nxi, nxt
      do 30 jj = 1, 2
      t2 = f(jj,j,2*k1-1)
      t3 = f(jj,j,2*k1)
      f(jj,j,2*k1-1) = f(jj,j,2*k-1)
      f(jj,j,2*k1) = f(jj,j,2*k)
      f(jj,j,2*k-1) = t2
      f(jj,j,2*k) = t3
   30 continue
   40 continue
   50 continue
c then transform in y
      nry = nxy/nyh
      do 100 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nxi, nxt
      do 60 jj = 1, 2
      t2 = real(t1)*f(jj,i,2*j2-1) - aimag(t1)*f(jj,i,2*j2)
      t3 = aimag(t1)*f(jj,i,2*j2-1) + real(t1)*f(jj,i,2*j2)
      f(jj,i,2*j2-1) = f(jj,i,2*j1-1) - t2
      f(jj,i,2*j2) = f(jj,i,2*j1) - t3
      f(jj,i,2*j1-1) = f(jj,i,2*j1-1) + t2
      f(jj,i,2*j1) = f(jj,i,2*j1) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = nxi, nxt
         do 110 jj = 1, 2
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = ani*(t2 + t4)
         f(jj,j,2*k) = ani*(t3 + t5)
         f(jj,j,ny+3-2*k) = ani*(t2 - t4)
         f(jj,j,ny+3-2*k+1) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         do 150 j = nxi, nxt
         do 140 jj = 1, 2
         f(jj,j,nyh+2) = -f(jj,j,nyh+2)
         t2 = f(jj,j,1) + f(jj,j,2)
         f(jj,j,2) = f(jj,j,1) - f(jj,j,2)
         f(jj,j,1) = t2
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 180 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 170 j = nxi, nxt
         do 160 jj = 1, 2
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = t2 + t4
         f(jj,j,2*k) = t3 + t5
         f(jj,j,ny+3-2*k) = t2 - t4
         f(jj,j,ny+3-2*k+1) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 j = nxi, nxt
         do 190 jj = 1, 2
         f(jj,j,nyh+1) = 2.0*f(jj,j,nyh+1)
         f(jj,j,nyh+2) = -2.0*f(jj,j,nyh+2)
         t2 = 2.0*(f(jj,j,1) + f(jj,j,2))
         f(jj,j,2) = 2.0*(f(jj,j,1) - f(jj,j,2))
         f(jj,j,1) = t2
         f(jj,j,ny+1) = 2.0*f(jj,j,ny+1)
  190    continue
  200    continue
      endif
c perform recursion for sine-cosine transform
      do 220 j = nxi, nxt
      sum1 = .5*f(1,j,1)
      f(1,j,1) = 0.
      f(1,j,2) = sum1
      sum2 = f(2,j,ny+1)
      f(2,j,ny+1) = f(2,j,2)
      f(2,j,2) = sum2
      do 210 k = 2, nyh
      sum1 = sum1 + f(1,j,2*k-1)
      f(1,j,2*k-1) = -f(1,j,2*k)
      f(1,j,2*k) = sum1
      sum2 = sum2 - f(2,j,2*k)
      f(2,j,2*k) = sum2
  210 continue
      f(1,j,ny+1) = 0.0
  220 continue
      return
      end
      subroutine FCST2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxhd,nyd,
     1nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 20)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,m) = (.5*f(1,n,1) + ((-1)**m)*f(1,n,ny+1) + sum(f(1,n,k)*
c            cos(pi*m*k/ny))
c f(2,n,m) = sum(f(2,n,k)*sin(pi*m*k/ny))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,n,k) = 2*(.5*f(1,n,1) + ((-1)**m)*f(1,n,ny+1) + sum(f(1,n,m)*
c       cos(pi*m*k/ny))
c f(3,n,k) = 2*sum(f(3,n,m)*sin(pi*m*k/ny)
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = half of the first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny + 1
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      do 20 j = nxi, nxt
      sum1 = .5*(f(1,j,1) - f(1,j,ny+1))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = f(1,j,ny2-k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,j,ny2-k) = at1 + at2
      at2 = f(2,j,ny2-k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,j,ny2-k) = at1 - at2
   10 continue
      f(1,j,1) = .5*(f(1,j,1) + f(1,j,ny+1))
      f(1,j,ny+1) = sum1
      f(2,j,1) = 0.0
      f(2,j,nyh+1) = 2.0*f(2,j,nyh+1)
   20 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 50 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = nxi, nxt
      do 30 jj = 1, 2
      t2 = f(jj,j,2*k1-1)
      t3 = f(jj,j,2*k1)
      f(jj,j,2*k1-1) = f(jj,j,2*k-1)
      f(jj,j,2*k1) = f(jj,j,2*k)
      f(jj,j,2*k-1) = t2
      f(jj,j,2*k) = t3
   30 continue
   40 continue
   50 continue
c then transform in y
      nry = nxy/nyh
      do 100 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nxi, nxt
      do 60 jj = 1, 2
      t2 = real(t1)*f(jj,i,2*j2-1) - aimag(t1)*f(jj,i,2*j2)
      t3 = aimag(t1)*f(jj,i,2*j2-1) + real(t1)*f(jj,i,2*j2)
      f(jj,i,2*j2-1) = f(jj,i,2*j1-1) - t2
      f(jj,i,2*j2) = f(jj,i,2*j1) - t3
      f(jj,i,2*j1-1) = f(jj,i,2*j1-1) + t2
      f(jj,i,2*j1) = f(jj,i,2*j1) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = nxi, nxt
         do 110 jj = 1, 2
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = ani*(t2 + t4)
         f(jj,j,2*k) = ani*(t3 + t5)
         f(jj,j,ny+3-2*k) = ani*(t2 - t4)
         f(jj,j,ny+3-2*k+1) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         do 150 j = nxi, nxt
         do 140 jj = 1, 2
         f(jj,j,nyh+2) = -f(jj,j,nyh+2)
         t2 = f(jj,j,1) + f(jj,j,2)
         f(jj,j,2) = f(jj,j,1) - f(jj,j,2)
         f(jj,j,1) = t2
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 180 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 170 j = nxi, nxt
         do 160 jj = 1, 2
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = t2 + t4
         f(jj,j,2*k) = t3 + t5
         f(jj,j,ny+3-2*k) = t2 - t4
         f(jj,j,ny+3-2*k+1) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 j = nxi, nxt
         do 190 jj = 1, 2
         f(jj,j,nyh+1) = 2.0*f(jj,j,nyh+1)
         f(jj,j,nyh+2) = -2.0*f(jj,j,nyh+2)
         t2 = 2.0*(f(jj,j,1) + f(jj,j,2))
         f(jj,j,2) = 2.0*(f(jj,j,1) - f(jj,j,2))
         f(jj,j,1) = t2
         f(jj,j,ny+1) = 2.0*f(jj,j,ny+1)
  190    continue
  200    continue
      endif
c perform recursion for cosine-sine transform
      do 220 j = nxi, nxt
      sum1 = f(1,j,ny+1)
      f(1,j,ny+1) = f(1,j,2)
      f(1,j,2) = sum1
      sum2 = .5*f(2,j,1)
      f(2,j,1) = 0.
      f(2,j,2) = sum2
      do 210 k = 2, nyh
      sum1 = sum1 - f(1,j,2*k)
      f(1,j,2*k) = sum1
      sum2 = sum2 + f(2,j,2*k-1)
      f(2,j,2*k-1) = -f(2,j,2*k)
      f(2,j,2*k) = sum2
  210 continue
      f(2,j,ny+1) = 0.0
  220 continue
      return
      end
      subroutine WFCST2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for 3 real cosine-sine transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
c perform y sine-cosine transforms
         call FSCST2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nx
     1hyd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y sine-cosine transforms
         call FSCST2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nx
     1hyd,nxyd)
c perform x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
      endif
      return
      end
      subroutine WFSCT2R3(f,isign,mixup,sctd,indx,indy,nxhd,nyd,nxhyd,nx
     1yd)
c wrapper function for 3 real sine-cosine transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      real f
      complex sctd
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
c local data
      integer nx, ny, nxi, nyi, nxt, nyt
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nxt = nx + 1
      nyt = ny + 1
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
c perform y cosine-sine transforms
         call FCSCT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nx
     1hyd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y cosineosine transforms
         call FCSCT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxhd,nyd,nx
     1hyd,nxyd)
c perform x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
      endif
      return
      end
      subroutine FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd
     1,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + ((-1)**n)*f(1,nx+1,k)
c             + sum(f(1,j,k)*cos(pi*n*j/nx)))
c f(2:3,n,k) = (1/nx*ny)*sum(f(2:3,j,k)*sin(pi*n*j/nx))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = 2*(.5*f(1,1,k) + ((-1)**j)*f(1,n+1,k) + sum(f(1,n,k)*
c       cos(pi*n*j/nx))
c f(2:3,j,k) = sum(f(2:3,n,k)*sin(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(1,1,k) - f(1,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      sum1 = sum1 + real(sctd(j1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,nx+2-j,k) = at1 + at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,nx+2-j,k) = at1 - at2
      at2 = f(3,nx+2-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(3,j,k) = at1 + at2
      f(3,nx+2-j,k) = at1 - at2
   10 continue
      f(1,1,k) = .5*(f(1,1,k) + f(1,nx+1,k))
      f(1,nx+1,k) = sum1
      f(2,1,k) = 0.0
      f(2,nxh+1,k) = 2.0*f(2,nxh+1,k)   
      f(3,1,k) = 0.0
      f(3,nxh+1,k) = 2.0*f(3,nxh+1,k)
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 3
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 3
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 k = nyi, nyt
      sum1 = f(1,nx+1,k)
      f(1,nx+1,k) = f(1,2,k)
      f(1,2,k) = sum1
      sum2 = .5*f(2,1,k)
      f(2,1,k) = 0.0
      f(2,2,k) = sum2
      sum3 = .5*f(3,1,k)
      f(3,1,k) = 0.0
      f(3,2,k) = sum3
      do 210 j = 2, nxh
      sum1 = sum1 - f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 + f(2,2*j-1,k)
      f(2,2*j-1,k) = -f(2,2*j,k)
      f(2,2*j,k) = sum2
      sum3 = sum3 + f(3,2*j-1,k)
      f(3,2*j-1,k) = -f(3,2*j,k)
      f(3,2*j,k) = sum3
  210 continue
      f(2,nx+1,k) = 0.0
      f(3,nx+1,k) = 0.0
  220 continue
      return
      end
      subroutine FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyp,nxhd,nyd
     1,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 18)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,k) = (1/nx*ny)*sum(f(1,j,k)*sin(pi*n*j/nx))
c f(2:3,n,k) = (1/nx*ny)*(.5*f(2:3,1,k) + ((-1)**n)*f(2:3,nx+1,k)
c             + sum(f(2:3,j,k)*cos(pi*n*j/nx)))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,j,k) = sum(f(1,n,k)*sin(pi*n*j/nx))
c f(2:3,j,k) = 2*(.5*f(2:3,1,k) + ((-1)**j)*f(2:3,n+1,k)
c             + sum(f(2:3,n,k)*cos(pi*n*j/nx))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = .5*(f(2,1,k) - f(2,nx+1,k))
      sum2 = .5*(f(3,1,k) - f(3,nx+1,k))
      do 10 j = 2, nxh
      j1 = 1 + kmr*(j - 1)
      at3 = -aimag(sctd(j1))
      at4 = real(sctd(j1))
      at2 = f(1,nx+2-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,nx+2-j,k) = at1 - at2
      at2 = f(2,nx+2-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      sum1 = sum1 + at4*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,nx+2-j,k) = at1 + at2
      at2 = f(3,nx+2-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      sum2 = sum2 + at4*at2
      at2 = at3*at2
      at1 = .5*at1
      f(3,j,k) = at1 - at2
      f(3,nx+2-j,k) = at1 + at2
   10 continue
      f(1,1,k) = 0.0
      f(1,nxh+1,k) = 2.0*f(1,nxh+1,k)
      f(2,1,k) = .5*(f(2,1,k) + f(2,nx+1,k))
      f(2,nx+1,k) = sum1
      f(3,1,k) = .5*(f(3,1,k) + f(3,nx+1,k))
      f(3,nx+1,k) = sum2
   20 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 50 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 50
      do 40 k = nyi, nyt
      do 30 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   30 continue
   40 continue
   50 continue
c first transform in x
      nrx = nxy/nxh
      do 100 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nyi, nyt
      do 60 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nxh
         ani = 1./float(2*nx*ny)
         do 130 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 120 k = nyi, nyt
         do 110 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = ani*(t2 + t4)
         f(jj,2*j,k) = ani*(t3 + t5)
         f(jj,nx3-2*j,k) = ani*(t2 - t4)
         f(jj,nx3-2*j+1,k) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         ani = 2.*ani
         do 150 k = nyi, nyt
         do 140 jj = 1, 3
         f(jj,nxh+1,k) = ani*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -ani*f(jj,nxh+2,k)
         t2 = ani*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = ani*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = ani*f(jj,nx+1,k)
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nxh
         do 180 j = 2, nxhh
         t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
         do 170 k = nyi, nyt
         do 160 jj = 1, 3
         t4 = f(jj,nx3-2*j,k)
         t5 = -f(jj,nx3-2*j+1,k)
         t2 = f(jj,2*j-1,k) + t4
         t3 = f(jj,2*j,k) + t5
         t6 = f(jj,2*j-1,k) - t4
         t5 = f(jj,2*j,k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,2*j-1,k) = t2 + t4
         f(jj,2*j,k) = t3 + t5
         f(jj,nx3-2*j,k) = t2 - t4
         f(jj,nx3-2*j+1,k) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 k = nyi, nyt
         do 190 jj = 1, 3
         f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
         f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
         t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
         f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
         f(jj,1,k) = t2
         f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 k = nyi, nyt
      sum1 = .5*f(1,1,k)
      f(1,1,k) = 0.0
      f(1,2,k) = sum1
      sum2 = f(2,nx+1,k)
      f(2,nx+1,k) = f(2,2,k)
      f(2,2,k) = sum2
      sum3 = f(3,nx+1,k)
      f(3,nx+1,k) = f(3,2,k)
      f(3,2,k) = sum3
      do 210 j = 2, nxh
      sum1 = sum1 + f(1,2*j-1,k)
      f(1,2*j-1,k) = -f(1,2*j,k)
      f(1,2*j,k) = sum1
      sum2 = sum2 - f(2,2*j,k)
      f(2,2*j,k) = sum2
      sum3 = sum3 - f(3,2*j,k)
      f(3,2*j,k) = sum3
  210 continue
      f(1,nx+1,k) = 0.0
  220 continue
      return
      end
      subroutine FSCST2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxhd,nyd
     1,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 20)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,m) = sum(f(1,n,k)*sin(pi*m*k/ny))
c f(2,n,m) = (.5*f(2,n,1) + ((-1)**m)*f(2,n,ny+1) + sum(f(2,n,k)*
c            cos(pi*m*k/ny))
c f(3,n,m) = sum(f(3,n,k)*sin(pi*m*k/ny))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,n,k) = 2*sum(f(1,n,m)*sin(pi*m*k/ny)
c f(2,n,k) = 2*(.5*f(2,n,1) + ((-1)**m)*f(2,n,ny+1) + sum(f(2,n,m)*
c       cos(pi*m*k/ny))
c f(3,n,k) = 2*sum(f(3,n,m)*sin(pi*m*k/ny)
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = half of the first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny + 1
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      do 20 j = nxi, nxt
      sum1 = .5*(f(2,j,1) - f(2,j,ny+1))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at2 = f(1,j,ny2-k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,j,ny2-k) = at1 - at2
      at2 = f(2,j,ny2-k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      sum1 = sum1 + real(sctd(k1))*at2
      at2 = at3*at2
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,j,ny2-k) = at1 + at2
      at2 = f(3,j,ny2-k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(3,j,k) = at1 + at2
      f(3,j,ny2-k) = at1 - at2
   10 continue
      f(1,j,1) = 0.0
      f(1,j,nyh+1) = 2.0*f(1,j,nyh+1)
      f(2,j,1) = .5*(f(2,j,1) + f(2,j,ny+1))
      f(2,j,ny+1) = sum1
      f(3,j,1) = 0.0
      f(3,j,nyh+1) = 2.0*f(3,j,nyh+1)
   20 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 50 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = nxi, nxt
      do 30 jj = 1, 3
      t2 = f(jj,j,2*k1-1)
      t3 = f(jj,j,2*k1)
      f(jj,j,2*k1-1) = f(jj,j,2*k-1)
      f(jj,j,2*k1) = f(jj,j,2*k)
      f(jj,j,2*k-1) = t2
      f(jj,j,2*k) = t3
   30 continue
   40 continue
   50 continue
c then transform in y
      nry = nxy/nyh
      do 100 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nxi, nxt
      do 60 jj = 1, 3
      t2 = real(t1)*f(jj,i,2*j2-1) - aimag(t1)*f(jj,i,2*j2)
      t3 = aimag(t1)*f(jj,i,2*j2-1) + real(t1)*f(jj,i,2*j2)
      f(jj,i,2*j2-1) = f(jj,i,2*j1-1) - t2
      f(jj,i,2*j2) = f(jj,i,2*j1) - t3
      f(jj,i,2*j1-1) = f(jj,i,2*j1-1) + t2
      f(jj,i,2*j1) = f(jj,i,2*j1) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = nxi, nxt
         do 110 jj = 1, 3
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = ani*(t2 + t4)
         f(jj,j,2*k) = ani*(t3 + t5)
         f(jj,j,ny+3-2*k) = ani*(t2 - t4)
         f(jj,j,ny+3-2*k+1) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         do 150 j = nxi, nxt
         do 140 jj = 1, 3
         f(jj,j,nyh+2) = -f(jj,j,nyh+2)
         t2 = f(jj,j,1) + f(jj,j,2)
         f(jj,j,2) = f(jj,j,1) - f(jj,j,2)
         f(jj,j,1) = t2
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 180 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 170 j = nxi, nxt
         do 160 jj = 1, 3
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = t2 + t4
         f(jj,j,2*k) = t3 + t5
         f(jj,j,ny+3-2*k) = t2 - t4
         f(jj,j,ny+3-2*k+1) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 j = nxi, nxt
         do 190 jj = 1, 3
         f(jj,j,nyh+1) = 2.0*f(jj,j,nyh+1)
         f(jj,j,nyh+2) = -2.0*f(jj,j,nyh+2)
         t2 = 2.0*(f(jj,j,1) + f(jj,j,2))
         f(jj,j,2) = 2.0*(f(jj,j,1) - f(jj,j,2))
         f(jj,j,1) = t2
         f(jj,j,ny+1) = 2.0*f(jj,j,ny+1)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 j = nxi, nxt
      sum1 = .5*f(1,j,1)
      f(1,j,1) = 0.
      f(1,j,2) = sum1
      sum2 = f(2,j,ny+1)
      f(2,j,ny+1) = f(2,j,2)
      f(2,j,2) = sum2
      sum3 = .5*f(3,j,1)
      f(3,j,1) = 0.
      f(3,j,2) = sum3
      do 210 k = 2, nyh
      sum1 = sum1 + f(1,j,2*k-1)
      f(1,j,2*k-1) = -f(1,j,2*k)
      f(1,j,2*k) = sum1
      sum2 = sum2 - f(2,j,2*k)
      f(2,j,2*k) = sum2
      sum3 = sum3 + f(3,j,2*k-1)
      f(3,j,2*k-1) = -f(3,j,2*k)
      f(3,j,2*k) = sum3
  210 continue
      f(1,j,ny+1) = 0.0
      f(3,j,ny+1) = 0.0
  220 continue
      return
      end
      subroutine FCSCT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxhd,nyd
     1,nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of x,
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions, using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 508.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 20)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms are performed
c f(1,n,m) = (.5*f(1,n,1) + ((-1)**m)*f(1,n,ny+1) + sum(f(1,n,k)*
c            cos(pi*m*k/ny))
c f(2,n,m) = sum(f(2,n,k)*sin(pi*m*k/ny))
c f(3,n,m) = (.5*f(3,n,1) + ((-1)**m)*f(3,n,ny+1) + sum(f(3,n,k)*
c            cos(pi*m*k/ny))
c if isign = 1, forward sine-cosine transforms are performed
c f(1,n,k) = 2*(.5*f(1,n,1) + ((-1)**m)*f(1,n,ny+1) + sum(f(1,n,m)*
c       cos(pi*m*k/ny))
c f(3,n,k) = 2*sum(f(3,n,m)*sin(pi*m*k/ny)
c f(3,n,k) = 2*(.5*f(3,n,1) + ((-1)**m)*f(3,n,ny+1) + sum(f(3,n,m)*
c       cos(pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxhd = half of the first dimension of f >= nx/2 + 1
c nyd = second dimension of f >= ny + 1
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.eq.0) return
c create auxiliary array in y
      kmr = nxy/ny
      do 20 j = nxi, nxt
      sum1 = .5*(f(1,j,1) - f(1,j,ny+1))
      sum2 = .5*(f(3,j,1) - f(3,j,ny+1))
      do 10 k = 2, nyh
      k1 = 1 + kmr*(k - 1)
      at3 = -aimag(sctd(k1))
      at4 = real(sctd(k1))
      at2 = f(1,j,ny2-k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      sum1 = sum1 + at4*at2
      at2 = at3*at2
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,j,ny2-k) = at1 + at2
      at2 = f(2,j,ny2-k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at3*at1
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,j,ny2-k) = at1 - at2
      at2 = f(3,j,ny2-k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      sum2 = sum2 + at4*at2
      at2 = at3*at2
      at1 = .5*at1
      f(3,j,k) = at1 - at2
      f(3,j,ny2-k) = at1 + at2
   10 continue
      f(1,j,1) = .5*(f(1,j,1) + f(1,j,ny+1))
      f(1,j,ny+1) = sum1
      f(2,j,1) = 0.0
      f(2,j,nyh+1) = 2.0*f(2,j,nyh+1)
      f(3,j,1) = .5*(f(3,j,1) + f(3,j,ny+1))
      f(3,j,ny+1) = sum2
   20 continue
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 50 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 50
      do 40 j = nxi, nxt
      do 30 jj = 1, 3
      t2 = f(jj,j,2*k1-1)
      t3 = f(jj,j,2*k1)
      f(jj,j,2*k1-1) = f(jj,j,2*k-1)
      f(jj,j,2*k1) = f(jj,j,2*k)
      f(jj,j,2*k-1) = t2
      f(jj,j,2*k) = t3
   30 continue
   40 continue
   50 continue
c then transform in y
      nry = nxy/nyh
      do 100 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 70 i = nxi, nxt
      do 60 jj = 1, 3
      t2 = real(t1)*f(jj,i,2*j2-1) - aimag(t1)*f(jj,i,2*j2)
      t3 = aimag(t1)*f(jj,i,2*j2-1) + real(t1)*f(jj,i,2*j2)
      f(jj,i,2*j2-1) = f(jj,i,2*j1-1) - t2
      f(jj,i,2*j2) = f(jj,i,2*j1) - t3
      f(jj,i,2*j1-1) = f(jj,i,2*j1-1) + t2
      f(jj,i,2*j1) = f(jj,i,2*j1) + t3
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble coefficients and normalize
c inverse fourier transform
      if (isign.lt.0) then
         kmr = nxy/nyh
         ani = 0.5
         do 130 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 120 j = nxi, nxt
         do 110 jj = 1, 3
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = ani*(t2 + t4)
         f(jj,j,2*k) = ani*(t3 + t5)
         f(jj,j,ny+3-2*k) = ani*(t2 - t4)
         f(jj,j,ny+3-2*k+1) = ani*(t5 - t3)
  110    continue
  120    continue
  130    continue
         do 150 j = nxi, nxt
         do 140 jj = 1, 3
         f(jj,j,nyh+2) = -f(jj,j,nyh+2)
         t2 = f(jj,j,1) + f(jj,j,2)
         f(jj,j,2) = f(jj,j,1) - f(jj,j,2)
         f(jj,j,1) = t2
  140    continue
  150    continue
c forward fourier transform
      else if (isign.gt.0) then
         kmr = nxy/nyh
         do 180 k = 2, nyhh
         t1 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
         do 170 j = nxi, nxt
         do 160 jj = 1, 3
         t4 = f(jj,j,ny+3-2*k)
         t5 = -f(jj,j,ny+3-2*k+1)
         t2 = f(jj,j,2*k-1) + t4
         t3 = f(jj,j,2*k) + t5
         t6 = f(jj,j,2*k-1) - t4
         t5 = f(jj,j,2*k) - t5
         t4 = t6*real(t1) - t5*aimag(t1)
         t5 = t6*aimag(t1) + t5*real(t1)
         f(jj,j,2*k-1) = t2 + t4
         f(jj,j,2*k) = t3 + t5
         f(jj,j,ny+3-2*k) = t2 - t4
         f(jj,j,ny+3-2*k+1) = t5 - t3
  160    continue
  170    continue
  180    continue
         do 200 j = nxi, nxt
         do 190 jj = 1, 3
         f(jj,j,nyh+1) = 2.0*f(jj,j,nyh+1)
         f(jj,j,nyh+2) = -2.0*f(jj,j,nyh+2)
         t2 = 2.0*(f(jj,j,1) + f(jj,j,2))
         f(jj,j,2) = 2.0*(f(jj,j,1) - f(jj,j,2))
         f(jj,j,1) = t2
         f(jj,j,ny+1) = 2.0*f(jj,j,ny+1)
  190    continue
  200    continue
      endif
c perform recursion for cosine transform
      do 220 j = nxi, nxt
      sum1 = f(1,j,ny+1)
      f(1,j,ny+1) = f(1,j,2)
      f(1,j,2) = sum1
      sum2 = .5*f(2,j,1)
      f(2,j,1) = 0.
      f(2,j,2) = sum2
      sum3 = f(3,j,ny+1)
      f(3,j,ny+1) = f(3,j,2)
      f(3,j,2) = sum3
      do 210 k = 2, nyh
      sum1 = sum1 - f(1,j,2*k)
      f(1,j,2*k) = sum1
      sum2 = sum2 + f(2,j,2*k-1)
      f(2,j,2*k-1) = -f(2,j,2*k)
      f(2,j,2*k) = sum2
      sum3 = sum3 - f(3,j,2*k)
      f(3,j,2*k) = sum3
  210 continue
      f(2,j,ny+1) = 0.0
  220 continue
      return
      end
      subroutine WFSFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxhy
     1d,nxyd)
c wrapper function for real sine/periodic transform
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine transform
         call FST2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
      subroutine WFCFT2RX(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxhy
     1d,nxyd)
c wrapper function for real cosine/periodic transform
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine transform
         call FCT2RXX(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxhy
     1d,nxyd)
      endif
      return
      end
      subroutine FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxd,nyhd,
     1nxhyd,nxyd)
c this subroutine performs the y part of a two dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c when the x transform is a real sine or cosine transform,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(j,m) = sum(f(j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(j,k) = sum(f(j,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd = first dimension of f, must be >= nx+1
c nyhd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c fourier coefficients are stored as follows:
c f(j,1) = real part of mode 0, f(j,2) = real part of mode ny/2
c f(2*j-1,k),f(2*j,k) = real,imaginary part of mode k-1, 0 < k < ny/2
c written by viktor k. decyk, ucla
      complex f, sctd, t1, t2, t3
      dimension f(nxd,nyhd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 100
c inverse fourier transform
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t2 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t2
   10 continue
   20 continue
c then transform in y
      nry = nxy/nyh
      do 60 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nyh
      ani = 0.5
      do 80 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 70 j = nxi, nxt
      t2 = conjg(f(j,nyh2-k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(j,nyh2-k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      do 90 j = nxi, nxt
      f(j,nyhh+1) = conjg(f(j,nyhh+1))
      f(j,1) = cmplx(real(f(j,1)) + aimag(f(j,1)),real(f(j,1)) - aimag(f
     1(j,1)))
   90 continue
      return
c forward fourier transform
c scramble coefficients
  100 kmr = nxy/nyh
      do 120 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 110 j = nxi, nxt
      t2 = conjg(f(j,nyh2-k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(j,nyh2-k) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 j = nxi, nxt
      f(j,nyhh+1) = 2.0*conjg(f(j,nyhh+1))
      f(j,1) = cmplx(real(f(j,1)) + aimag(f(j,1)),real(f(j,1)) - aimag(f
     1(j,1)))
  130 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 150 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 150
      do 140 j = nxi, nxt
      t2 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t2
  140 continue
  150 continue
c then transform in y
      nry = nxy/nyh
      do 190 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sctd(1+kmr*(j-1)))
      do 160 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
      subroutine RLTOCX2(f,s,isign,nyh,nxi,nxt,nxd,nyd)
c this subroutine transforms real to complex 2d arrays
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nyh = complex dimension in y direction
c nxi/nxt = initial/final x index used
c nxd = half of the first dimension of f
c nyd = second dimension of f
      implicit none
      integer isign, nyh, nxi, nxt, nxd, nyd
      real f, s
      dimension f(2*nxd*nyd), s(2*nxd)
c local data
      integer j, k, nx2d, joff
      nx2d = 2*nxd
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 30 k = 1, nyh
         joff = nx2d*(k - 1)
         do 10 j = nxi, nxt
         s(2*j-1) = f(j+joff)
         s(2*j) = f(j+joff+nxd)
   10    continue
         do 20 j = nxi, nxt
         f(2*j+joff-1) = s(2*j-1)
         f(2*j+joff) = s(2*j)
   20    continue
   30    continue
c complex to real
      else if (isign.gt.0) then
         do 60 k = 1, nyh
         joff = nx2d*(k - 1)
         do 40 j = nxi, nxt
         s(2*j-1) = f(2*j+joff-1)
         s(2*j) = f(2*j+joff)
   40    continue
         do 50 j = nxi, nxt
         f(j+joff) = s(2*j-1)
         f(j+joff+nxd) = s(2*j)
   50    continue
   60    continue
      endif
      return
      end
      subroutine WFCSFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 2 real cosine-sine/periodic transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine-sine transforms
         call FCST2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
      endif
      return
      end
      subroutine WFSCFT2R2(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 2 real sine-cosine/periodic transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(2,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine-cosine transforms
         call FSCT2R2X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nxh
     1yd,nxyd)
      endif
      return
      end
      subroutine FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxd,nyhd,
     1nxhyd,nxyd)
c this subroutine performs the y part of 2 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c when the x transform is a real sine or cosine transform,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:2,j,m) = sum(f(1:2,j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(1:2,j,k) = sum(f(1:2,j,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd = first dimension of f, must be >= nx+1
c nyhd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c fourier coefficients are stored as follows:
c f(j,1) = real part of mode 0, f(j,2) = real part of mode ny/2
c f(2*j-1,k),f(2*j,k) = real,imaginary part of mode k-1, 0 < k < ny/2
c written by viktor k. decyk, ucla
      complex f, sctd, t1, t2, t3
      dimension f(2,nxd,nyhd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
   10 continue
   20 continue
c then transform in y
      nry = nxy/nyh
      do 60 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nyh
      ani = 0.5
      do 90 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 80 j = nxi, nxt
      do 70 jj = 1, 2
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,j,nyh2-k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
   90 continue
      do 110 j = nxi, nxt
      do 100 jj = 1, 2
      f(jj,j,nyhh+1) = conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  100 continue
  110 continue
      return
c forward fourier transform
c scramble coefficients
  120 kmr = nxy/nyh
      do 150 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 140 j = nxi, nxt
      do 130 jj = 1, 2
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,j,nyh2-k) = conjg(t1 - t2)
  130 continue
  140 continue
  150 continue
      do 170 j = nxi, nxt
      do 160 jj = 1, 2
      f(jj,j,nyhh+1) = 2.0*conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  160 continue
  170 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 190 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 190
      do 180 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
  180 continue
  190 continue
c then transform in y
      nry = nxy/nyh
      do 230 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sctd(1+kmr*(j-1)))
      do 200 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
      subroutine RLTOCX22(f,s,isign,nyh,nxi,nxt,nxd,nyd)
c this subroutine transforms 2 real to complex 2d arrays
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nyh = complex dimension in y direction
c nxi/nxt = initial/final x index used
c ndim = first dimension of f
c nxd = half of the second dimension of f
c nyd = third dimension of f
      implicit none
      integer isign, nyh, nxi, nxt, nxd, nyd
      real f, s
      dimension f(2,2*nxd*nyd), s(2,2*nxd)
c local data
      integer i, j, k, nx2d, joff
      nx2d = 2*nxd
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 40 k = 1, nyh
         joff = nx2d*(k - 1)
         do 20 j = nxi, nxt
         do 10 i = 1, 2
         s(i,2*j-1) = f(i,j+joff)
         s(i,2*j) = f(i,j+joff+nxd)
   10    continue
   20    continue
         do 30 j = nxi, nxt
         f(1,2*j+joff-1) = s(1,2*j-1)
         f(2,2*j+joff-1) = s(1,2*j)
         f(1,2*j+joff) = s(2,2*j-1)
         f(2,2*j+joff) = s(2,2*j)
   30    continue
   40    continue
c complex to real
      else if (isign.gt.0) then
         do 80 k = 1, nyh
         joff = nx2d*(k - 1)
         do 50 j = nxi, nxt
         s(1,2*j-1) = f(1,2*j+joff-1)
         s(1,2*j) = f(2,2*j+joff-1)
         s(2,2*j-1) = f(1,2*j+joff)
         s(2,2*j)= f(2,2*j+joff)
   50    continue
         do 70 j = nxi, nxt
         do 60 i = 1, 2
         f(i,j+joff) = s(i,2*j-1)
         f(i,j+joff+nxd) = s(i,2*j)
   60    continue
   70    continue
   80    continue
      endif
      return
      end
      subroutine WFCSFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 3 real cosine-sine/periodic transforms
c for the electric field with dirichlet or magnetic field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(3,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine-sine transforms
         call FCSST2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
      endif
      return
      end
      subroutine WFSCFT2R3(f,ss,isign,mixup,sctd,indx,indy,nxhd,nyhd,nxh
     1yd,nxyd)
c wrapper function for 3 real sine-cosine/periodic transforms
c for the magnetic field with dirichlet or electric field with neumann
c boundary conditions
      implicit none
      complex f, sctd, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), ss(3,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nxt = nx + 1
      nyt = ny
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine-cosine transforms
         call FSCCT2R3X(f,isign,mixup,sctd,indx,indy,nyi,nyt,nxhd,nyd,nx
     1hyd,nxyd)
      endif
      return
      end
      subroutine FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxp,nxd,nyhd,
     1nxhyd,nxyd)
c this subroutine performs the y part of 3 two dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c when the x transform is a real sine or cosine transform,
c using complex arithmetic
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse fourier transform is performed
c f(1:3,j,m) = sum(f(1:3,j,k)*exp(-sqrt(-1)*2pi*m*k/ny))
c if isign = 1, a forward fourier transform is performed
c f(1:3,j,k) = sum(f(1:3,j,m)*exp(sqrt(-1)*2pi*m*k/ny))
c mixup = array of bit reversed addresses
c sctd = sine/cosine table
c nxi = initial x index used
c nxp = number of x indices used
c nxd = first dimension of f, must be >= nx+1
c nyhd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c fourier coefficients are stored as follows:
c f(j,1) = real part of mode 0, f(j,2) = real part of mode ny/2
c f(2*j-1,k),f(2*j,k) = real,imaginary part of mode k-1, 0 < k < ny/2
c written by viktor k. decyk, ucla
      complex f, sctd, t1, t2, t3, t4
      dimension f(3,nxd,nyhd), mixup(nxhyd), sctd(nxyd)
      indx1 = indx - 1
      indy1 = indy - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyhh = ny/4
      nyh2 = nyh + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 120
c inverse fourier transform
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 20 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxy/nyh
      do 60 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxy/nyh
      ani = 0.5
      do 90 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),-real(sctd(1+kmr*(k-1))))
      do 80 j = nxi, nxt
      do 70 jj = 1, 3
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,j,nyh2-k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
   90 continue
      do 110 j = nxi, nxt
      do 100 jj = 1, 3
      f(jj,j,nyhh+1) = conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  100 continue
  110 continue
      return
c forward fourier transform
c scramble coefficients
  120 kmr = nxy/nyh
      do 150 k = 2, nyhh
      t3 = cmplx(aimag(sctd(1+kmr*(k-1))),real(sctd(1+kmr*(k-1))))
      do 140 j = nxi, nxt
      do 130 jj = 1, 3
      t2 = conjg(f(jj,j,nyh2-k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,j,nyh2-k) = conjg(t1 - t2)
  130 continue
  140 continue
  150 continue
      do 170 j = nxi, nxt
      do 160 jj = 1, 3
      f(jj,j,nyhh+1) = 2.0*conjg(f(jj,j,nyhh+1))
      f(jj,j,1) = cmplx(real(f(jj,j,1)) + aimag(f(jj,j,1)),real(f(jj,j,1
     1)) - aimag(f(jj,j,1)))
  160 continue
  170 continue
      nry = nxhy/nyh
c bit-reverse array elements in y
      nry = nxhy/nyh
      do 190 k = 1, nyh
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 190
      do 180 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  180 continue
  190 continue
c then transform in y
      nry = nxy/nyh
      do 230 l = 1, indy1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyhh/ns
      kmr = 2*km*nry
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sctd(1+kmr*(j-1)))
      do 200 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
      subroutine RLTOCX23(f,s,isign,nyh,nxi,nxt,nxd,nyd)
c this subroutine transforms 2 real to complex 2d arrays
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nyh = complex dimension in y direction
c nxi/nxt = initial/final x index used
c ndim = first dimension of f
c nxd = half of the second dimension of f
c nyd = third dimension of f
      implicit none
      integer isign, nyh, nxi, nxt, nxd, nyd
      real f, s
      dimension f(3,2*nxd*nyd), s(3,2*nxd)
c local data
      integer i, j, k, nx2d, joff
      nx2d = 2*nxd
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 40 k = 1, nyh
         joff = nx2d*(k - 1)
         do 20 j = nxi, nxt
         do 10 i = 1, 3
         s(i,2*j-1) = f(i,j+joff)
         s(i,2*j) = f(i,j+joff+nxd)
   10    continue
   20    continue
         do 30 j = nxi, nxt
         f(1,2*j+joff-1) = s(1,2*j-1)
         f(2,2*j+joff-1) = s(1,2*j)
         f(3,2*j+joff-1) = s(2,2*j-1)
         f(1,2*j+joff) = s(2,2*j)
         f(2,2*j+joff) = s(3,2*j-1)
         f(3,2*j+joff) = s(3,2*j)
   30    continue
   40    continue
c complex to real
      else if (isign.gt.0) then
         do 80 k = 1, nyh
         joff = nx2d*(k - 1)
         do 50 j = nxi, nxt
         s(1,2*j-1) = f(1,2*j+joff-1)
         s(1,2*j) = f(2,2*j+joff-1)
         s(2,2*j-1) = f(3,2*j+joff-1)
         s(2,2*j) = f(1,2*j+joff)
         s(3,2*j-1) = f(2,2*j+joff)
         s(3,2*j) = f(3,2*j+joff)
   50    continue
         do 70 j = nxi, nxt
         do 60 i = 1, 3
         f(i,j+joff) = s(i,2*j-1)
         f(i,j+joff+nxd) = s(i,2*j)
   60    continue
   70    continue
   80    continue
      endif
      return
      end
      subroutine WFDT2RINIT(sctdx,indx,nxd)
c this subroutine calculates extra table needed by a one dimensional
c fast real sine and cosine transforms for mixed boundary conditions
c and their inverses.
c input: indx, nxd
c output: sctdx
c sctdx = sine/cosine table
c indx = exponent which determines length in x direction,
c where nx=2**indx
c nxd = must be >= nx
c written by viktor k. decyk, ucla
      implicit none
      integer indx, nxd
      complex sctdx
      dimension sctdx(nxd)
c local data
      integer j, nx
      real dnx, arg, at1, at2, at3, at4
      nx = 2**indx
c sine/cosine table for the angles n*pi/nx
      dnx = 6.28318530717959/float(nx + nx)
      arg = 0.5*dnx
      at3 = cos(arg)
      at4 = sin(arg)
      do 10 j = 1, nx
      arg = dnx*float(j - 1)
      at1 = cos(arg)
      at2 = sin(arg)
      at1 = at1*at4 + at2*at3
      sctdx(j) = cmplx(at1,1.0/at1)
   10 continue
      return
      end
      subroutine WFDSFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,ny
     1hd,nxhyd,nxyd)
c wrapper function for real sine DST-III/periodic transform
      implicit none
      complex f, sctd, sctdx, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      dimension ss(2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine transform
         call FDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,n
     1yd,nxhyd,nxyd)
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine transform
         call FDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,n
     1yd,nxhyd,nxyd)
      endif
      return
      end
      subroutine WFDCFT2RX(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,ny
     1hd,nxhyd,nxyd)
c wrapper function for real cosine/periodic transform
      implicit none
      complex f, sctd, sctdx, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      dimension ss(2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine transform
         call FDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,n
     1yd,nxhyd,nxyd)
c transform real to complex
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y fft
         call FDFT2RXY(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX2(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine transform
         call FDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,n
     1yd,nxhyd,nxyd)
      endif
      return
      end
      subroutine FDST2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyp,nxh
     1d,nyd,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c sine transform and its inverse, for a subset of y,
c using real arithmetic
c algorithm is similar to the one described in Numerical Recipies in
c Fortran, Second Ed.,by W. H. Press, B. P. Flannery, S. A. Teukolsky,
c and W. T. Vetterling, [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse sine transform DST-III is performed
c f(n,k) = (1/nx*ny)*(.5*f(nx+1,k)*(-1)**n + sum(f(j+1,k)*
c sin(pi*(n+1/2)*j/nx)))
c if isign = 1, a forward sine transform DST-II is performed
c f(j+1,k) = 2.0*sum(f(n,k)*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, sctdx, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
      if (isign.gt.0) go to 140
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      f(1,k) = 2.0*f(2,k)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(2*j1,k) - f(2*j1-2,k)
      at3 = f(2*j1-1,k)*at2 + at4*at1
      f(2*j1,k) = at4*at2 - f(2*j1-1,k)*at1
      f(2*j1-1,k) = at3
   10 continue
      f(2,k) = f(nx+1,k)
   20 continue
c scramble coefficients
      kmr = nxy/nxh
      do 40 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 30 k = nyi, nyt
      t4 = f(nx3-2*j,k)
      t5 = -f(nx3-2*j+1,k)
      t2 = f(2*j-1,k) + t4
      t3 = f(2*j,k) + t5
      t6 = f(2*j-1,k) - t4
      t5 = f(2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(2*j-1,k) = t2 + t4
      f(2*j,k) = t3 + t5
      f(nx3-2*j,k) = t2 - t4
      f(nx3-2*j+1,k) = t5 - t3
   30 continue
   40 continue
      do 50 k = nyi, nyt
      f(nxh+1,k) = 2.0*f(nxh+1,k)
      f(nxh+2,k) = -2.0*f(nxh+2,k)
      t2 = f(1,k) + f(2,k)
      f(2,k) = f(1,k) - f(2,k)
      f(1,k) = t2
   50 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 70 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 110 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 100 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 90 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 80 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) + aimag(t1)*f(2*j2,i)
      t3 = real(t1)*f(2*j2,i) - aimag(t1)*f(2*j2-1,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   80 continue
   90 continue
  100 continue
  110 continue
c perform recursion for sine transform and normalize
      ani = 1.0/float(4*nx*ny)
      do 130 k = nyi, nyt
      do 120 j = 1, nxh
      at2 = f(nx+1-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(j,k) = ani*(at1 + at2)
      f(nx+1-j,k) = ani*(at1 - at2)
  120 continue
      f(nx+1,k) = 0.0
  130 continue
      return
c forward fourier transform
c create auxiliary array in x
  140 kmr = nxy/nx
      do 160 k = nyi, nyt
      do 150 j = 1, nxh
      at2 = f(nx+1-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(j,k) = at1 + at2
      f(nx+1-j,k) = at1 - at2
  150 continue
  160 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 180 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 180
      do 170 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
  170 continue
  180 continue
c first transform in x
      nrx = nxy/nxh
      do 220 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 190 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
  190 continue
  200 continue
  210 continue
  220 continue
c unscramble coefficients
      kmr = nxy/nxh
      do 240 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 230 k = nyi, nyt
      t4 = f(nx3-2*j,k)
      t5 = -f(nx3-2*j+1,k)
      t2 = f(2*j-1,k) + t4
      t3 = f(2*j,k) + t5
      t6 = f(2*j-1,k) - t4
      t5 = f(2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(2*j-1,k) = t2 + t4
      f(2*j,k) = t3 + t5
      f(nx3-2*j,k) = t2 - t4
      f(nx3-2*j+1,k) = t5 - t3
  230 continue
  240 continue
      do 250 k = nyi, nyt
      f(nxh+1,k) = 2.0*f(nxh+1,k)
      f(nxh+2,k) = -2.0*f(nxh+2,k)
      t2 = 2.0*(f(1,k) + f(2,k))
      f(2,k) = 2.0*(f(1,k) - f(2,k))
      f(1,k) = t2
      f(nx+1,k) = 2.0*f(nx+1,k)
  250 continue
c perform recursion for sine transform
      kmr = nxy/nx
      do 270 k = nyi, nyt
      f(nx+1,k) = f(2,k)
      f(2,k) = 0.5*f(1,k)
      f(1,k) = 0.0
      sum1 = f(2,k)
      do 260 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(2*j-1,k)*at2 - f(2*j,k)*at1
      at2 = f(2*j-1,k)*at1 + f(2*j,k)*at2
      f(2*j-1,k) = at3
      sum1 = sum1 + at2
      f(2*j,k) = sum1
  260 continue
  270 continue
      return
      end
      subroutine FDCT2RXX(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyp,nxh
     1d,nyd,nxhyd,nxyd)
c this subroutine performs the x part of a two dimensional fast real
c cosine transform and its inverse, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, an inverse cosine transform DCT-III is performed
c f(n,k) = (1/nx*ny)*(.5*f(1,k) + sum(f(j,k)*cos(pi*(n+1/2)*j/nx)))
c if isign = 1, a forward cosine transform DCT-II is performed
c f(j,k) = 2.0*sum(f(n,k)*cos(pi*n*(j+1/2)/nx))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, sctdx, t1
      dimension f(2*nxhd,nyd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
      if (isign.gt.0) go to 140
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = f(nx,k)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(2*j1,k) - f(2*j1-2,k)
      at3 = f(2*j1-1,k)*at1 + at4*at2
      f(2*j1,k) = f(2*j1-1,k)*at2 - at4*at1
      f(2*j1-1,k) = at3
   10 continue
      f(2,k) = -2.0*sum1
   20 continue
c scramble coefficients
      kmr = nxy/nxh
      do 40 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 30 k = nyi, nyt
      t4 = f(nx3-2*j,k)
      t5 = -f(nx3-2*j+1,k)
      t2 = f(2*j-1,k) + t4
      t3 = f(2*j,k) + t5
      t6 = f(2*j-1,k) - t4
      t5 = f(2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(2*j-1,k) = t2 + t4
      f(2*j,k) = t3 + t5
      f(nx3-2*j,k) = t2 - t4
      f(nx3-2*j+1,k) = t5 - t3
   30 continue
   40 continue
      do 50 k = nyi, nyt
      f(nxh+1,k) = 2.0*f(nxh+1,k)
      f(nxh+2,k) = -2.0*f(nxh+2,k)
      t2 = f(1,k) + f(2,k)
      f(2,k) = f(1,k) - f(2,k)
      f(1,k) = t2
   50 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 70 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
   60 continue
   70 continue
c first transform in x
      nrx = nxy/nxh
      do 110 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 100 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 90 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 80 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) + aimag(t1)*f(2*j2,i)
      t3 = real(t1)*f(2*j2,i) - aimag(t1)*f(2*j2-1,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
   80 continue
   90 continue
  100 continue
  110 continue
c perform recursion for cosine transform and normalize
      ani = 1.0/float(4*nx*ny)
      do 130 k = nyi, nyt
      do 120 j = 1, nxh
      at2 = f(nx+1-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(j,k) = ani*(at1 - at2)
      f(nx+1-j,k) = ani*(at1 + at2)
  120 continue
      f(nx+1,k) = 0.0
  130 continue
      return
c forward fourier transform
c create auxiliary array in x
  140 kmr = nxy/nx
      do 160 k = nyi, nyt
      do 150 j = 1, nxh
      at2 = f(nx+1-j,k)
      at1 = f(j,k) + at2
      at2 = f(j,k) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(j,k) = at1 - at2
      f(nx+1-j,k) = at1 + at2
  150 continue
  160 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 180 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 180
      do 170 k = nyi, nyt
      t2 = f(2*j1-1,k)
      t3 = f(2*j1,k)
      f(2*j1-1,k) = f(2*j-1,k)
      f(2*j1,k) = f(2*j,k)
      f(2*j-1,k) = t2
      f(2*j,k) = t3
  170 continue
  180 continue
c first transform in x
      nrx = nxy/nxh
      do 220 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 190 i = nyi, nyt
      t2 = real(t1)*f(2*j2-1,i) - aimag(t1)*f(2*j2,i)
      t3 = aimag(t1)*f(2*j2-1,i) + real(t1)*f(2*j2,i)
      f(2*j2-1,i) = f(2*j1-1,i) - t2
      f(2*j2,i) = f(2*j1,i) - t3
      f(2*j1-1,i) = f(2*j1-1,i) + t2
      f(2*j1,i) = f(2*j1,i) + t3
  190 continue
  200 continue
  210 continue
  220 continue
c unscramble coefficients
      kmr = nxy/nxh
      do 240 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 230 k = nyi, nyt
      t4 = f(nx3-2*j,k)
      t5 = -f(nx3-2*j+1,k)
      t2 = f(2*j-1,k) + t4
      t3 = f(2*j,k) + t5
      t6 = f(2*j-1,k) - t4
      t5 = f(2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(2*j-1,k) = t2 + t4
      f(2*j,k) = t3 + t5
      f(nx3-2*j,k) = t2 - t4
      f(nx3-2*j+1,k) = t5 - t3
  230 continue
  240 continue
      do 250 k = nyi, nyt
      f(nxh+1,k) = 2.0*f(nxh+1,k)
      f(nxh+2,k) = -2.0*f(nxh+2,k)
      t2 = 2.0*(f(1,k) + f(2,k))
      f(2,k) = 2.0*(f(1,k) - f(2,k))
      f(1,k) = t2
      f(nx+1,k) = 2.0*f(nx+1,k)
  250 continue
c perform recursion for cosine transform
      kmr = nxy/nx
      do 270 k = nyi, nyt
      sum1 = -0.5*f(2,k)
      do 260 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(2*j1-1,k)*at1 + f(2*j1,k)*at2
      at2 = -f(2*j1-1,k)*at2 + f(2*j1,k)*at1
      f(2*j1-1,k) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2*j1,k) = at1
  260 continue
      f(2,k) = sum1
      f(nx+1,k) = 0.0
  270 continue
      return
      end
      subroutine WFDCSFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd)
c wrapper function for 2 real cosine-sine transforms
c for the electric field with mixed dirichlet-neumann or magnetic field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      complex f, sctd, sctdx, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(2,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call FDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,
     1nyd,nxhyd,nxyd)
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine-sine transform
         call FDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,
     1nyd,nxhyd,nxyd)
      endif
      return
      end
      subroutine WFDSCFT2R2(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd)
c wrapper function for 2 real sine-cosine transforms
c for the magnetic field with mixed dirichlet-neumann or electric field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      complex f, sctd, sctdx, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(2,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(2,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transform
         call FDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,
     1nyd,nxhyd,nxyd)
c transform real to complex
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R2Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX22(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine-cosine transform
         call FDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd,
     1nyd,nxhyd,nxyd)
      endif
      return
      end
      subroutine FDCST2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyp,nx
     1hd,nyd,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + sum(f(1,j,k)*cos(pi*(n+1/2)*j/nx)))
c f(2,n,k) = (1/nx*ny)*(.5*f(2,nx+1,k)*(-1)**n + sum(f(2,j+1,k)*
c sin(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j,k) = 2.0*sum(f(1,n,k)*cos(pi*n*(j+1/2)/nx))
c f(2,j+1,k) = 2.0*sum(f(2,n,k)*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, sctdx, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
      if (isign.gt.0) go to 180
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = f(1,nx,k)
      f(2,1,k) = 2.0*f(2,2,k)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k) - f(1,2*j1-2,k)
      at3 = f(1,2*j1-1,k)*at1 + at4*at2
      f(1,2*j1,k) = f(1,2*j1-1,k)*at2 - at4*at1
      f(1,2*j1-1,k) = at3
      at4 = f(2,2*j1,k) - f(2,2*j1-2,k)
      at3 = f(2,2*j1-1,k)*at2 + at4*at1
      f(2,2*j1,k) = at4*at2 - f(2,2*j1-1,k)*at1
      f(2,2*j1-1,k) = at3
   10 continue
      f(1,2,k) = -2.0*sum1
      f(2,2,k) = f(2,nx+1,k)
   20 continue
c scramble coefficients
      kmr = nxy/nxh
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 40 k = nyi, nyt
      do 30 jj = 1, 2
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
   30 continue
   40 continue
   50 continue
      do 70 k = nyi, nyt
      do 60 jj = 1, 2
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = f(jj,1,k) + f(jj,2,k)
      f(jj,2,k) = f(jj,1,k) - f(jj,2,k)
      f(jj,1,k) = t2
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = nyi, nyt
      do 80 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   80 continue
   90 continue
  100 continue
c first transform in x
      nrx = nxy/nxh
      do 150 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 120 i = nyi, nyt
      do 110 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) + aimag(t1)*f(jj,2*j2,i)
      t3 = real(t1)*f(jj,2*j2,i) - aimag(t1)*f(jj,2*j2-1,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      do 170 k = nyi, nyt
      do 160 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(1,j,k) = ani*(at1 - at2)
      f(1,nx+1-j,k) = ani*(at1 + at2)
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(2,j,k) = ani*(at1 + at2)
      f(2,nx+1-j,k) = ani*(at1 - at2)
  160 continue
      f(1,nx+1,k) = 0.0
      f(2,nx+1,k) = 0.0
  170 continue
      return
c forward fourier transform
c create auxiliary array in x
  180 kmr = nxy/nx
      do 200 k = nyi, nyt
      do 190 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,nx+1-j,k) = at1 + at2
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,nx+1-j,k) = at1 - at2
  190 continue
  200 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 230 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 230
      do 220 k = nyi, nyt
      do 210 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
  210 continue
  220 continue
  230 continue
c first transform in x
      nrx = nxy/nxh
      do 280 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 250 i = nyi, nyt
      do 240 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  240 continue
  250 continue
  260 continue
  270 continue
  280 continue
c unscramble coefficients
      kmr = nxy/nxh
      do 310 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 300 k = nyi, nyt
      do 290 jj = 1, 2
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
  290 continue
  300 continue
  310 continue
      do 330 k = nyi, nyt
      do 320 jj = 1, 2
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
      f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
      f(jj,1,k) = t2
      f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  320 continue
  330 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      do 360 k = nyi, nyt
      sum1 = -0.5*f(1,2,k)
      do 340 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(1,2*j1-1,k)*at1 + f(1,2*j1,k)*at2
      at2 = -f(1,2*j1-1,k)*at2 + f(1,2*j1,k)*at1
      f(1,2*j1-1,k) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(1,2*j1,k) = at1
  340 continue
      f(1,2,k) = sum1
      f(1,nx+1,k) = 0.0
      f(2,nx+1,k) = f(2,2,k)
      f(2,2,k) = 0.5*f(2,1,k)
      f(2,1,k) = 0.0
      sum1 = f(2,2,k)
      do 350 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(2,2*j-1,k)*at2 - f(2,2*j,k)*at1
      at2 = f(2,2*j-1,k)*at1 + f(2,2*j,k)*at2
      f(2,2*j-1,k) = at3
      sum1 = sum1 + at2
      f(2,2*j,k) = sum1
  350 continue
  360 continue
      return
      end
      subroutine FDSCT2R2X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyp,nx
     1hd,nyd,nxhyd,nxyd)
c this subroutine performs the x part of 2 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,nx+1,k)*(-1)**n + sum(f(1,j+1,k)*
c sin(pi*(n+1/2)*j/nx)))
c f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + sum(f(2,j,k)*cos(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j+1,k) = 2.0*sum(f(1,n,k)*sin(pi*n*(j+1/2)/nx)))
c f(2,j,k) = 2.0*sum(f(2,n,k)*cos(pi*n*(j+1/2)/nx))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, sctdx, t1
      dimension f(2,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
      if (isign.gt.0) go to 180
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      f(1,1,k) = 2.0*f(1,2,k)
      sum1 = f(2,nx,k)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k) - f(1,2*j1-2,k)
      at3 = f(1,2*j1-1,k)*at2 + at4*at1
      f(1,2*j1,k) = at4*at2 - f(1,2*j1-1,k)*at1
      f(1,2*j1-1,k) = at3
      at4 = f(2,2*j1,k) - f(2,2*j1-2,k)
      at3 = f(2,2*j1-1,k)*at1 + at4*at2
      f(2,2*j1,k) = f(2,2*j1-1,k)*at2 - at4*at1
      f(2,2*j1-1,k) = at3
   10 continue
      f(1,2,k) = f(1,nx+1,k)
      f(2,2,k) = -2.0*sum1
   20 continue
c scramble coefficients
      kmr = nxy/nxh
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 40 k = nyi, nyt
      do 30 jj = 1, 2
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
   30 continue
   40 continue
   50 continue
      do 70 k = nyi, nyt
      do 60 jj = 1, 2
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = f(jj,1,k) + f(jj,2,k)
      f(jj,2,k) = f(jj,1,k) - f(jj,2,k)
      f(jj,1,k) = t2
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = nyi, nyt
      do 80 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   80 continue
   90 continue
  100 continue
c first transform in x
      nrx = nxy/nxh
      do 150 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 120 i = nyi, nyt
      do 110 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) + aimag(t1)*f(jj,2*j2,i)
      t3 = real(t1)*f(jj,2*j2,i) - aimag(t1)*f(jj,2*j2-1,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      do 170 k = nyi, nyt
      do 160 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(1,j,k) = ani*(at1 + at2)
      f(1,nx+1-j,k) = ani*(at1 - at2)
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(2,j,k) = ani*(at1 - at2)
      f(2,nx+1-j,k) = ani*(at1 + at2)
  160 continue
      f(1,nx+1,k) = 0.0
      f(2,nx+1,k) = 0.0
  170 continue
      return
c forward fourier transform
c create auxiliary array in x
  180 kmr = nxy/nx
      do 200 k = nyi, nyt
      do 190 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,nx+1-j,k) = at1 - at2
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,nx+1-j,k) = at1 + at2
  190 continue
  200 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 230 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 230
      do 220 k = nyi, nyt
      do 210 jj = 1, 2
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
  210 continue
  220 continue
  230 continue
c first transform in x
      nrx = nxy/nxh
      do 280 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 250 i = nyi, nyt
      do 240 jj = 1, 2
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  240 continue
  250 continue
  260 continue
  270 continue
  280 continue
c unscramble coefficients
      kmr = nxy/nxh
      do 310 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 300 k = nyi, nyt
      do 290 jj = 1, 2
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
  290 continue
  300 continue
  310 continue
      do 330 k = nyi, nyt
      do 320 jj = 1, 2
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
      f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
      f(jj,1,k) = t2
      f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  320 continue
  330 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      do 360 k = nyi, nyt
      f(1,nx+1,k) = f(1,2,k)
      f(1,2,k) = 0.5*f(1,1,k)
      f(1,1,k) = 0.0
      sum1 = f(1,2,k)
      do 340 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(1,2*j-1,k)*at2 - f(1,2*j,k)*at1
      at2 = f(1,2*j-1,k)*at1 + f(1,2*j,k)*at2
      f(1,2*j-1,k) = at3
      sum1 = sum1 + at2
      f(1,2*j,k) = sum1
  340 continue
      sum1 = -0.5*f(2,2,k)
      do 350 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(2,2*j1-1,k)*at1 + f(2,2*j1,k)*at2
      at2 = -f(2,2*j1-1,k)*at2 + f(2,2*j1,k)*at1
      f(2,2*j1-1,k) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2,2*j1,k) = at1
  350 continue
      f(2,2,k) = sum1
      f(2,nx+1,k) = 0.0
  360 continue
      return
      end
      subroutine WFDCSFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd)
c wrapper function for 3 real cosine-sine transforms
c for the electric field with mixed dirichlet-neumann or magnetic field
c with mixed neumann-dirichelet boundary conditions
      implicit none
      complex f, sctd, sctdx, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(3,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x cosine-sine transform
         call FDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd
     1,nyd,nxhyd,nxyd)
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x cosine-sine transform
         call FDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd
     1,nyd,nxhyd,nxyd)
      endif
      return
      end
      subroutine WFDSCFT2R3(f,ss,isign,mixup,sctd,sctdx,indx,indy,nxhd,n
     1yhd,nxhyd,nxyd)
c wrapper function for 3 real sine-cosine transforms
c for the magnetic field with mixed dirichlet-neumann or electric field
c with mixed neumann-dirichlet boundary conditions
      implicit none
      complex f, sctd, sctdx, ss
      integer mixup
      integer isign, indx, indy, nxhd, nyhd, nxhyd, nxyd
      dimension f(3,2*nxhd,nyhd), mixup(nxhyd), sctd(nxyd)
      dimension sctdx(2*nxhd), ss(3,2*nxhd)
c local data
      integer nx, ny, nyh, nxi, nyi, nxt, nyt, nxd, nyd
      data nxi, nyi /1,1/
c calculate range of indices
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      nyt = ny
      nxt = nx + 1
      nxd = 2*nxhd
      nyd = 2*nxyd
c inverse fourier transform
      if (isign.lt.0) then
c perform x sine-cosine transform
         call FDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd
     1,nyd,nxhyd,nxyd)
c transform real to complex
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c forward fourier transform
      else if (isign.gt.0) then
c perform y ffts
         call FDFT2R3Y(f,isign,mixup,sctd,indx,indy,nxi,nxt,nxd,nyhd,nxh
     1yd,nxyd)
c transform complex to real
         call RLTOCX23(f,ss,isign,nyh,nxi,nxt,nxd,nyhd)
c perform x sine-cosine transform
         call FDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyt,nxhd
     1,nyd,nxhyd,nxyd)
      endif
      return
      end
      subroutine FDCSST2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyp,n
     1xhd,nyd,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,1,k) + sum(f(1,j,k)*cos(pi*(n+1/2)*j/nx)))
c f(2,n,k) = (1/nx*ny)*(.5*f(2,nx+1,k)*(-1)**n + sum(f(2,j+1,k)*
c sin(pi*(n+1/2)*j/nx)))
c f(3,n,k) = (1/nx*ny)*(.5*f(3,nx+1,k)*(-1)**n + sum(f(3,j+1,k)*
c sin(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j,k) = 2.0*sum(f(1,n,k)*cos(pi*n*(j+1/2)/nx))
c f(2,j+1,k) = 2.0*sum(f(2,n,k)*sin(pi*n*(j+1/2)/nx)))
c f(3,j+1,k) = 2.0*sum(f(3,n,k)*sin(pi*n*(j+1/2)/nx)))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, sctdx, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
      if (isign.gt.0) go to 180
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      sum1 = f(1,nx,k)
      f(2,1,k) = 2.0*f(2,2,k)
      f(3,1,k) = 2.0*f(3,2,k)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k) - f(1,2*j1-2,k)
      at3 = f(1,2*j1-1,k)*at1 + at4*at2
      f(1,2*j1,k) = f(1,2*j1-1,k)*at2 - at4*at1
      f(1,2*j1-1,k) = at3
      at4 = f(2,2*j1,k) - f(2,2*j1-2,k)
      at3 = f(2,2*j1-1,k)*at2 + at4*at1
      f(2,2*j1,k) = at4*at2 - f(2,2*j1-1,k)*at1
      f(2,2*j1-1,k) = at3
      at4 = f(3,2*j1,k) - f(3,2*j1-2,k)
      at3 = f(3,2*j1-1,k)*at2 + at4*at1
      f(3,2*j1,k) = at4*at2 - f(3,2*j1-1,k)*at1
      f(3,2*j1-1,k) = at3
   10 continue
      f(1,2,k) = -2.0*sum1
      f(2,2,k) = f(2,nx+1,k)
      f(3,2,k) = f(3,nx+1,k)
   20 continue
c scramble coefficients
      kmr = nxy/nxh
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 40 k = nyi, nyt
      do 30 jj = 1, 3
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
   30 continue
   40 continue
   50 continue
      do 70 k = nyi, nyt
      do 60 jj = 1, 3
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = f(jj,1,k) + f(jj,2,k)
      f(jj,2,k) = f(jj,1,k) - f(jj,2,k)
      f(jj,1,k) = t2
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = nyi, nyt
      do 80 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   80 continue
   90 continue
  100 continue
c first transform in x
      nrx = nxy/nxh
      do 150 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 120 i = nyi, nyt
      do 110 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) + aimag(t1)*f(jj,2*j2,i)
      t3 = real(t1)*f(jj,2*j2,i) - aimag(t1)*f(jj,2*j2-1,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      do 170 k = nyi, nyt
      do 160 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(1,j,k) = ani*(at1 - at2)
      f(1,nx+1-j,k) = ani*(at1 + at2)
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(2,j,k) = ani*(at1 + at2)
      f(2,nx+1-j,k) = ani*(at1 - at2)
      at2 = f(3,nx+1-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(3,j,k) = ani*(at1 + at2)
      f(3,nx+1-j,k) = ani*(at1 - at2)
  160 continue
      f(1,nx+1,k) = 0.0
      f(2,nx+1,k) = 0.0
      f(3,nx+1,k) = 0.0
  170 continue
      return
c forward fourier transform
c create auxiliary array in x
  180 kmr = nxy/nx
      do 200 k = nyi, nyt
      do 190 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(1,j,k) = at1 - at2
      f(1,nx+1-j,k) = at1 + at2
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(2,j,k) = at1 + at2
      f(2,nx+1-j,k) = at1 - at2
      at2 = f(3,nx+1-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(3,j,k) = at1 + at2
      f(3,nx+1-j,k) = at1 - at2
  190 continue
  200 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 230 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 230
      do 220 k = nyi, nyt
      do 210 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
  210 continue
  220 continue
  230 continue
c first transform in x
      nrx = nxy/nxh
      do 280 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 250 i = nyi, nyt
      do 240 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  240 continue
  250 continue
  260 continue
  270 continue
  280 continue
c unscramble coefficients
      kmr = nxy/nxh
      do 310 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 300 k = nyi, nyt
      do 290 jj = 1, 3
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
  290 continue
  300 continue
  310 continue
      do 330 k = nyi, nyt
      do 320 jj = 1, 3
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
      f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
      f(jj,1,k) = t2
      f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  320 continue
  330 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      do 370 k = nyi, nyt
      sum1 = -0.5*f(1,2,k)
      do 340 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(1,2*j1-1,k)*at1 + f(1,2*j1,k)*at2
      at2 = -f(1,2*j1-1,k)*at2 + f(1,2*j1,k)*at1
      f(1,2*j1-1,k) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(1,2*j1,k) = at1
  340 continue
      f(1,2,k) = sum1
      f(1,nx+1,k) = 0.0
      f(2,nx+1,k) = f(2,2,k)
      f(2,2,k) = 0.5*f(2,1,k)
      f(2,1,k) = 0.0
      sum1 = f(2,2,k)
      do 350 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(2,2*j-1,k)*at2 - f(2,2*j,k)*at1
      at2 = f(2,2*j-1,k)*at1 + f(2,2*j,k)*at2
      f(2,2*j-1,k) = at3
      sum1 = sum1 + at2
      f(2,2*j,k) = sum1
  350 continue
      f(3,nx+1,k) = f(3,2,k)
      f(3,2,k) = 0.5*f(3,1,k)
      f(3,1,k) = 0.0
      sum1 = f(3,2,k)
      do 360 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(3,2*j-1,k)*at2 - f(3,2*j,k)*at1
      at2 = f(3,2*j-1,k)*at1 + f(3,2*j,k)*at2
      f(3,2*j-1,k) = at3
      sum1 = sum1 + at2
      f(3,2*j,k) = sum1
  360 continue
  370 continue
      return
      end
      subroutine FDSCCT2R3X(f,isign,mixup,sctd,sctdx,indx,indy,nyi,nyp,n
     1xhd,nyd,nxhyd,nxyd)
c this subroutine performs the x part of 3 two dimensional fast real
c sine and cosine transforms and their inverses, for a subset of y,
c using real arithmetic
c algorithm is described in Numerical Recipies in Fortran, Second Ed.,
c by W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, 
c [Cambridge Univ. Press, 1992], p. 513.
c for isign = (-1,1), input: all, output: f
c approximate flop count: N*(5*log2(N) + 22)
c where N = (nx/2)*ny
c indx/indy = exponent which determines length in x/y direction,
c where nx=2**indx, ny=2**indy
c if isign = -1, inverse sine-cosine transforms DST-III/DCT-III are
c performed
c f(1,n,k) = (1/nx*ny)*(.5*f(1,nx+1,k)*(-1)**n + sum(f(1,j+1,k)*
c sin(pi*(n+1/2)*j/nx)))
c f(2,n,k) = (1/nx*ny)*(.5*f(2,1,k) + sum(f(2,j,k)*cos(pi*(n+1/2)*j/nx)))
c f(3,n,k) = (1/nx*ny)*(.5*f(3,1,k) + sum(f(3,j,k)*cos(pi*(n+1/2)*j/nx)))
c if isign = 1, forward sine-cosine transforms DST-II/DCT-II are
c performed
c f(1,j+1,k) = 2.0*sum(f(1,n,k)*sin(pi*n*(j+1/2)/nx)))
c f(2,j,k) = 2.0*sum(f(2,n,k)*cos(pi*n*(j+1/2)/nx))
c f(3,j,k) = 2.0*sum(f(3,n,k)*cos(pi*n*(j+1/2)/nx))
c mixup = array of bit reversed addresses
c sctd, sctdx = sine/cosine tables
c nyi = initial y index used
c nyp = number of y indices used
c nxhd = first dimension of f >= nx/2 + 1
c nyd = second dimension of f
c nxhyd = maximum of (nx/2,ny)
c nxyd = maximum of (nx,ny)
c written by viktor k. decyk, ucla
      complex sctd, sctdx, t1
      dimension f(3,2*nxhd,nyd), mixup(nxhyd), sctd(nxyd), sctdx(2*nxhd)
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nx3 = nx + 3
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.eq.0) return
      if (isign.gt.0) go to 180
c inverse fourier transform
c create auxiliary array in x
      kmr = nxy/nx
      do 20 k = nyi, nyt
      f(1,1,k) = 2.0*f(1,2,k)
      sum1 = f(2,nx,k)
      sum2 = f(3,nx,k)
      do 10 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at4 = f(1,2*j1,k) - f(1,2*j1-2,k)
      at3 = f(1,2*j1-1,k)*at2 + at4*at1
      f(1,2*j1,k) = at4*at2 - f(1,2*j1-1,k)*at1
      f(1,2*j1-1,k) = at3
      at4 = f(2,2*j1,k) - f(2,2*j1-2,k)
      at3 = f(2,2*j1-1,k)*at1 + at4*at2
      f(2,2*j1,k) = f(2,2*j1-1,k)*at2 - at4*at1
      f(2,2*j1-1,k) = at3
      at4 = f(3,2*j1,k) - f(3,2*j1-2,k)
      at3 = f(3,2*j1-1,k)*at1 + at4*at2
      f(3,2*j1,k) = f(3,2*j1-1,k)*at2 - at4*at1
      f(3,2*j1-1,k) = at3
   10 continue
      f(1,2,k) = f(1,nx+1,k)
      f(2,2,k) = -2.0*sum1
      f(3,2,k) = -2.0*sum2
   20 continue
c scramble coefficients
      kmr = nxy/nxh
      do 50 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),real(sctd(1+kmr*(j-1))))
      do 40 k = nyi, nyt
      do 30 jj = 1, 3
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
   30 continue
   40 continue
   50 continue
      do 70 k = nyi, nyt
      do 60 jj = 1, 3
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = f(jj,1,k) + f(jj,2,k)
      f(jj,2,k) = f(jj,1,k) - f(jj,2,k)
      f(jj,1,k) = t2
   60 continue
   70 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 100 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 100
      do 90 k = nyi, nyt
      do 80 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
   80 continue
   90 continue
  100 continue
c first transform in x
      nrx = nxy/nxh
      do 150 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 120 i = nyi, nyt
      do 110 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) + aimag(t1)*f(jj,2*j2,i)
      t3 = real(t1)*f(jj,2*j2,i) - aimag(t1)*f(jj,2*j2-1,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c perform recursion for cosine/sine transforms and normalize
      ani = 1.0/float(4*nx*ny)
      do 170 k = nyi, nyt
      do 160 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = .5*at1*aimag(sctdx(j))
      f(1,j,k) = ani*(at1 + at2)
      f(1,nx+1-j,k) = ani*(at1 - at2)
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(2,j,k) = ani*(at1 - at2)
      f(2,nx+1-j,k) = ani*(at1 + at2)
      at2 = f(3,nx+1-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at2 = .5*at2*aimag(sctdx(j))
      f(3,j,k) = ani*(at1 - at2)
      f(3,nx+1-j,k) = ani*(at1 + at2)
  160 continue
      f(1,nx+1,k) = 0.0
      f(2,nx+1,k) = 0.0
      f(3,nx+1,k) = 0.0
  170 continue
      return
c forward fourier transform
c create auxiliary array in x
  180 kmr = nxy/nx
      do 200 k = nyi, nyt
      do 190 j = 1, nxh
      at2 = f(1,nx+1-j,k)
      at1 = f(1,j,k) + at2
      at2 = f(1,j,k) - at2
      at1 = at1*real(sctdx(j))
      at2 = .5*at2
      f(1,j,k) = at1 + at2
      f(1,nx+1-j,k) = at1 - at2
      at2 = f(2,nx+1-j,k)
      at1 = f(2,j,k) + at2
      at2 = f(2,j,k) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(2,j,k) = at1 - at2
      f(2,nx+1-j,k) = at1 + at2
      at2 = f(3,nx+1-j,k)
      at1 = f(3,j,k) + at2
      at2 = f(3,j,k) - at2
      at2 = at2*real(sctdx(j))
      at1 = .5*at1
      f(3,j,k) = at1 - at2
      f(3,nx+1-j,k) = at1 + at2
  190 continue
  200 continue
c bit-reverse array elements in x
      nrx = nxhy/nxh
      do 230 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 230
      do 220 k = nyi, nyt
      do 210 jj = 1, 3
      t2 = f(jj,2*j1-1,k)
      t3 = f(jj,2*j1,k)
      f(jj,2*j1-1,k) = f(jj,2*j-1,k)
      f(jj,2*j1,k) = f(jj,2*j,k)
      f(jj,2*j-1,k) = t2
      f(jj,2*j,k) = t3
  210 continue
  220 continue
  230 continue
c first transform in x
      nrx = nxy/nxh
      do 280 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = 2*km*nrx
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sctd(1+kmr*(j-1))
      do 250 i = nyi, nyt
      do 240 jj = 1, 3
      t2 = real(t1)*f(jj,2*j2-1,i) - aimag(t1)*f(jj,2*j2,i)
      t3 = aimag(t1)*f(jj,2*j2-1,i) + real(t1)*f(jj,2*j2,i)
      f(jj,2*j2-1,i) = f(jj,2*j1-1,i) - t2
      f(jj,2*j2,i) = f(jj,2*j1,i) - t3
      f(jj,2*j1-1,i) = f(jj,2*j1-1,i) + t2
      f(jj,2*j1,i) = f(jj,2*j1,i) + t3
  240 continue
  250 continue
  260 continue
  270 continue
  280 continue
c unscramble coefficients
      kmr = nxy/nxh
      do 310 j = 2, nxhh
      t1 = cmplx(aimag(sctd(1+kmr*(j-1))),-real(sctd(1+kmr*(j-1))))
      do 300 k = nyi, nyt
      do 290 jj = 1, 3
      t4 = f(jj,nx3-2*j,k)
      t5 = -f(jj,nx3-2*j+1,k)
      t2 = f(jj,2*j-1,k) + t4
      t3 = f(jj,2*j,k) + t5
      t6 = f(jj,2*j-1,k) - t4
      t5 = f(jj,2*j,k) - t5
      t4 = t6*real(t1) - t5*aimag(t1)
      t5 = t6*aimag(t1) + t5*real(t1)
      f(jj,2*j-1,k) = t2 + t4
      f(jj,2*j,k) = t3 + t5
      f(jj,nx3-2*j,k) = t2 - t4
      f(jj,nx3-2*j+1,k) = t5 - t3
  290 continue
  300 continue
  310 continue
      do 330 k = nyi, nyt
      do 320 jj = 1, 3
      f(jj,nxh+1,k) = 2.0*f(jj,nxh+1,k)
      f(jj,nxh+2,k) = -2.0*f(jj,nxh+2,k)
      t2 = 2.0*(f(jj,1,k) + f(jj,2,k))
      f(jj,2,k) = 2.0*(f(jj,1,k) - f(jj,2,k))
      f(jj,1,k) = t2
      f(jj,nx+1,k) = 2.0*f(jj,nx+1,k)
  320 continue
  330 continue
c perform recursion for cosine/sine transform
      kmr = nxy/nx
      do 370 k = nyi, nyt
      f(1,nx+1,k) = f(1,2,k)
      f(1,2,k) = 0.5*f(1,1,k)
      f(1,1,k) = 0.0
      sum1 = f(1,2,k)
      do 340 j = 2, nxh
      at1 = real(sctd(1+kmr*(j-1)))
      at2 = -aimag(sctd(1+kmr*(j-1)))
      at3 = f(1,2*j-1,k)*at2 - f(1,2*j,k)*at1
      at2 = f(1,2*j-1,k)*at1 + f(1,2*j,k)*at2
      f(1,2*j-1,k) = at3
      sum1 = sum1 + at2
      f(1,2*j,k) = sum1
  340 continue
      sum1 = -0.5*f(2,2,k)
      do 350 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(2,2*j1-1,k)*at1 + f(2,2*j1,k)*at2
      at2 = -f(2,2*j1-1,k)*at2 + f(2,2*j1,k)*at1
      f(2,2*j1-1,k) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(2,2*j1,k) = at1
  350 continue
      f(2,2,k) = sum1
      f(2,nx+1,k) = 0.0
      sum1 = -0.5*f(3,2,k)
      do 360 j = 2, nxh
      j1 = nxh + 2 - j
      at1 = real(sctd(1+kmr*(j1-1)))
      at2 = -aimag(sctd(1+kmr*(j1-1)))
      at3 = f(3,2*j1-1,k)*at1 + f(3,2*j1,k)*at2
      at2 = -f(3,2*j1-1,k)*at2 + f(3,2*j1,k)*at1
      f(3,2*j1-1,k) = at3
      at1 = sum1
      sum1 = sum1 + at2
      f(3,2*j1,k) = at1
  360 continue
      f(3,2,k) = sum1
      f(3,nx+1,k) = 0.0
  370 continue
      return
      end
