!-----------------------------------------------------------------------
!
      module input2d
!
! input2mod.f defines namelists containing input and output variables:
!             defines module input2d
! written by viktor k. decyk, ucla
! copyright 2012, regents of the university of california
! update: march 1, 2012
!
      use globals, only: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, &
     &VECTOR, PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D,         &
     &VACUUM_2D, VACUUM_3D, NEUMANN_2D, NEUMANN_PERIODIC_2D,            &
     &DIRICHLET_NEUMANN_PERIODIC_2D
      implicit none
      private
      public :: LINEAR, QUADRATIC, CUBIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PERIODIC_2D, DIRICHLET_2D, DIRICHLET_PERIODIC_2D
      public :: VACUUM_2D, VACUUM_3D, NEUMANN_2D, NEUMANN_PERIODIC_2D
      public :: DIRICHLET_NEUMANN_PERIODIC_2D
      public :: input2
      public :: idrun, idrun0, idcode, indx, indy, npx, npy, npxb, npyb
      public :: inorder, popt, dopt, djopt, nustrt, ntr
      public :: ntw, ntp, nta, ntv, nts, ntm, nte
      public :: ndw, ndp, nda, ndv, nds, ndm, nde
      public :: nsv
      public :: tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0
      public :: vdx, vdy, vdz, vtdx, vtdy, vtdz
      public :: psolve, relativity
      public :: omx, omy, omz, ci, ax, ay
      public :: ndim, ndc, movion
      public :: sortime, nplot, idpal, ndstyle, sntasks
      public :: nsrand, ndprof
      public :: ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy
      public :: monitor, itpon, v0, w0
      public :: modesxp, modesyp, modesxa, modesya, modesxe, modesye
      public :: ions2
      public :: ntd, ntj, npxi, npyi, npxbi, npybi
      public :: ndd, ndj
      public :: qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0
      public :: vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, sortimi
      public :: nsrandi, ndprofi
      public :: ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi
      public :: ionoff
      public :: modesxd, modesyd, modesxj, modesyj
      public :: pot2d, vpot2d, em2d
      public :: t0, ceng, indian, rlprec, inprec
      public :: nprec, fpname, narec, faname, nerec, fename
      public :: nlsrerec
      public :: den2d, vcur2d
      public :: ndrec, fdname, njrec, fjname
      public :: launchlas, lase0, laspol, lasstart,lasend,antlength
      public :: launchdir
      public :: launchtype, lask, ntlas, launchshape,lasw0,lasfocus, rayl
      public :: remass,affpin, laseronly, dampfield, dampdrop
      public :: gamma,smoothtype,smoothfactor,smoothmask
      public :: ntle1,ntle2,ntle3,ntlb1,ntlb2,ntlb3,ntqe,ntqi
      public :: ntcu,ntemb1,ntemb2,ntemb3,nteme1,nteme2,nteme3,ntpart
      public :: edensity, dampystart,dampyend,dampe1,dampe2,dampe3
      public :: npxbg, npybg, npxg, npyg
      public :: plasmabstart,plasmabend,plasmabystart,plasmabyend
      public :: beamstart,beamend,beamystart,beamyend
      public :: solvertype
      public :: xymaxwel
      public :: ntej1, ntej2, ntej3, ntlej1, ntlej2, ntlej3
      public :: nttj1, nttj2, nttj3, ntltj1, ntltj2, ntltj3
      public :: mzf, parmax, yee
      public :: vpml, xtrue, ytrue
!
! Namelist Input
      save
! idrun/idrun0 = run identifier for current/old run
! idcode = code identifier
      integer :: idrun = 1, idrun0 = 0, idcode = 0
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! npx/npy = initial number of particles distributed in x/y direction
!     integer :: indx =   5, indy =   6, npx =      96, npy =     192
      integer :: indx =   6, indy =   7, npx =     384, npy =     768
!     integer :: indx =   7, indy =   8, npx =    1280, npy =    2560
! npxb/npyb = initial number of particles in beam in x/y direction
      integer :: npxb =   0, npyb =   0
!     integer :: npxb =  32, npyb =  64
!     integer :: npxb = 128, npyb = 256
!     integer :: npxb = 384, npyb = 768
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
      integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
      integer :: djopt = STANDARD
! nustrt = (0,1,2) = this is an (old start,new start,restart)
! ntr = number of time steps between restart routine
      integer :: nustrt = 1, ntr = 0
! ntw, ndw = number of time steps between energy diagnostic
! ntp, ndp = number of time steps between potential diagnostic
! nta, nda = number of time steps between vector potential diagnostic
! ntv, ndv = number of time steps between velocity-space diagnostic
! nts, nds = number of time steps between phase space diagnostic
      integer :: ntw = 0, ntp = 0, nta = 0, ntv = 0, nts = 0
      integer :: ndw = 0, ndp = 0, nda = 0, ndv = 0, nds = 0
! ntm, ndm = number of time steps between momentum diagnostic
! nte, nde = number of time steps between electromagnetic diagnostic
      integer :: ntm = 0, nte = 0
      integer :: ndm = 0, nde = 0
! nsv = velocity component(s) for phase-space display(s), if nts > 0
! 1 = vx, 2 = vy, 3 = vx and vy, 4 = vz, 5 = vx and vz, 6 = vy and vz,
! 7 = vx and vy and vz
      integer :: nsv = 1
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend =  65.000, dt = 0.2000000e+00
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! psolve = type of poisson solver = (1,2,3)
      integer :: psolve = PERIODIC_2D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! ci = reciprical of velocity of light
      real :: ci = 0.1
! ax/ay = half-width of particle in x/y direction
!     real :: ax = .816497, ay = .816497
!     real :: ax = .866025, ay = .866025
      real :: ax = .912871, ay = .912871
! ndim = number of velocity dimensions = 2 or 3
      integer :: ndim = 3
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
! movion = (0,1) = (no,yes) move the ions
      integer :: movion = 0
! sortime = number of time steps between electron sorting
      integer :: sortime = 50
! nplot = maximum number of plots per page
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
! ndstyle = (1,2,3) = display (color map,contour plot,both)
      integer :: nplot = 0, idpal = 1, ndstyle = 1
! sntasks = (-1,n) = set maximum number of tasks (-1 = number of cpus-1)
      integer :: sntasks = -1
! nsrand = (0,1) = (no,yes) randomize spatially positions locally
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
      integer :: nsrand = 0, ndprof = 0
! ampdx/ampdy = amplitude of density compared to uniform in x/y
! scaledx/scaledy = scale length for spatial coordinate in x/y
! shiftdx/shiftdy = shift of spatial coordinate in x/y
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0
! monitor = (0,1,2) = (disable,normal,extended) error processing
      integer :: monitor = 0
! itpon = time when external pump is turned on (-1=never)
      integer :: itpon = -1
! v0 = external pump strength, in units vos/vthermal
! w0 = external pump frequency, in units of wpe
      real :: v0 = 0.0, w0 = 0.0
! modesxp/modesyp = number of modes in x/y to keep for potential
! diagnostic
      integer :: modesxp = 11, modesyp = 11
! modesxa/modesya = number of modes in x/y to keep for vector potential
! diagnostic
      integer :: modesxa = 11, modesya = 11
! modesxe/modesye = number of modes in x/y to keep for electromagnetic
! diagnostic
      integer :: modesxe = 11, modesye = 11
! laser initialization
      integer :: launchlas = 0, lasstart = 1, lasend = 1400,laspol = 3
      integer :: lasfocus = 1400
      integer :: launchtype = 1, launchshape = 1
      integer :: launchdir = 1
      real :: lase0 = 0.5, lask = 1.0, lasw0 = 10.0, rayl = 0.0
! plasma block edge
      integer :: plasmabstart = 0, plasmabend = 0
      integer :: plasmabystart = 0, plasmabyend = 0
      integer :: beamstart = 0, beamend = 0
      integer :: beamystart = 0, beamyend = 0
      
      integer :: dampystart = 0, dampyend = 0
      real :: remass = 1.0, gamma = 1.0
      integer :: npxg = 1, npyg = 1, npxbg = 1, npybg = 1
      integer :: ntlas = 10, antlength = 50
      integer :: nteme1 = 0, nteme2 = 0, nteme3 = 0
      integer :: ntemb1 = 0, ntemb2 = 0, ntemb3 = 0
      integer :: ntqe = 0, ntqi = 0
      integer :: ntle1 = 0, ntle2 = 0, ntle3 = 0
      integer :: ntlb1 = 0, ntlb2 = 0, ntlb3 = 0
      integer :: ntcu = 0, ntpart = 0
      real :: edensity = 1.0, dampdrop = 1.0
      real :: affpin = 1.0,  smoothfactor = 1.0,smoothmask = 1.0
      integer :: laseronly = 0, smoothtype = 1, dampfield = 0
      integer :: dampe1=0,dampe2=0,dampe3=0
      integer :: solvertype = 1
      integer :: xymaxwel = 0
      integer ::  ntej1 = 0, ntej2 = 0, ntej3 = 0
      integer :: ntlej1 = 0, ntlej2 = 0, ntlej3 = 0
      integer ::  nttj1 = 0, nttj2 = 0, nttj3 = 0
      integer :: ntltj1 = 0, ntltj2 = 0, ntltj3 = 0
      integer :: mzf = 0, parmax = 0, yee = 0
      integer :: vpml = 0, xtrue = 1, ytrue = 1
! define namelist
      namelist /input2/ idrun, idrun0, idcode, indx, indy, npx, npy,    &
     &npxb, npyb, inorder, popt, dopt, djopt, nustrt, ntr, ntw, ntp,    &
     &nta, ntv, nts, ntm, nte, ndw, ndp, nda, ndv, nds, ndm, nde, nsv,  &
     &tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0, vdx, vdy, vdz, vtdx, &
     &vtdy, vtdz, psolve, relativity, omx, omy, omz, ci, ax, ay, ndim,  &
     &ndc, movion, sortime, nplot, idpal, ndstyle, sntasks, nsrand,     &
     &ndprof, ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy,         &
     &monitor, itpon, v0, w0, modesxp, modesyp, modesxa, modesya,       &
     &modesxe, modesye, launchlas, lase0,laspol,lasstart,lasend,        &
     &plasmabstart,plasmabend, launchtype, lask,launchdir,remass,       &
     &edensity, npxg, npyg, ntlas, launchshape,lasw0,lasfocus, rayl,    &
     &gamma,              antlength,affpin,laseronly,smoothtype,        &
     &smoothfactor,smoothmask,dampfield,plasmabystart,plasmabyend,      &
     &dampystart,dampyend,ntle1,ntle2,ntle3,ntlb1,ntlb2,ntlb3,dampdrop, &
     &ntcu,ntemb1,ntemb2,ntemb3,nteme1,nteme2,nteme3,ntqe,ntqi, dampe1, &
     &dampe2,dampe3,ntpart,beamstart,beamend,beamystart,beamyend,npxbg, &
     &npybg, solvertype, xymaxwel, ntej1, ntej2, ntej3, ntlej1, ntlej2, &
     &ntlej3, nttj1, nttj2, nttj3, ntltj1, ntltj2, ntltj3, mzf, parmax, &
     &yee, vpml, xtrue, ytrue
!
! Namelist for Ions
! ntd, ndd = number of time steps between ion density diagnostic
! ntj, ndj = number of time steps between ion current diagnostic
      integer :: ntd = 0, ntj = 0
      integer :: ndd = 0,  ndj = 0
! npxi/npyi = initial number of ions distributed in x/y/z direction
      integer :: npxi =  384, npyi =  768
! npxbi/npybi = initial number of ions in beam in x/y/z direction
      integer :: npxbi =   0, npybi =   0
! qmi = charge on ion, in units of e
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 100.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0
! sortimi = number of time steps between ion sorting
      integer :: sortimi = 250
! nsrandi = (0,1) = (no,yes) randomize spatially ion positions locally
! ndprofi = ion profile (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
      integer :: nsrandi = 0, ndprofi = 0
! ampdxi/ampdyi = amplitude of ion density compared to uniform in x/y
! scaledxi/scaledyi = scale length for spatial ion coordinate in x/y
! shiftdxi/shiftdyi = shift of spatial ion coordinate in x/y
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
      real :: ampdyi = 0.0, scaledyi = 0.0, shiftdyi = 0.0
! ionoff = time when ions are frozen and their charge saved (-1=never)
      integer :: ionoff = -1
! modesxd/modesyd = number of modes in x/y to keep for ion density
! diagnostic
      integer :: modesxd = 11, modesyd = 11
! modesxj/modesyj = number of modes in x/y to keep for ion current
! diagnostic
      integer :: modesxj = 11, modesyj = 11
! define namelist
      namelist /ions2/ ntd, ntj, ndd, ndj, npxi, npyi, npxbi, npybi,    &
     &qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0, vdxi,    &
     &vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, sortimi, nsrandi,       &
     &ndprofi, ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi,  &
     &ionoff, modesxd, modesyd, modesxj, modesyj
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
! indian = (0,1) = architecture is (little-endian,big-endian)
! rlprec = (0,1) = default reals are (normal,double-precision)
! inprec = (0,1) = default integers are (normal,double-precision)
      integer :: indian = 1, rlprec = 1, inprec = 0
!
! Namelist output for potential diagnostic
! nprec = current record number for potential writes
      integer :: nprec = 0
! fpname = file name for potential diagnostic
      character(len=32) :: fpname = 'potk2.0'
! define namelist
      namelist /pot2d/ idrun, indx, indy, ntp, modesxp, modesyp, psolve,&
     & omx, omy, omz, nprec, t0, tend, dt, ceng, indian, rlprec, inprec,&
     & fpname
!
! Namelist output for vector potential diagnostic
! narec = current record number for vector potential writes
      integer :: narec = 0
! faname = file name for vector potential diagnostic
      character(len=32) :: faname = 'vpotk2.0'
! define namelist
      namelist /vpot2d/ idrun, indx, indy, nta, modesxa, modesya,       &
     &psolve, ndim, omx, omy, omz, ci, narec, t0, tend, dt, ceng,       &
     &indian, rlprec, inprec, faname
!
! Namelist output for electromagnetic diagnostic
! nerec = current record number for electromagnetic writes
      integer :: nerec = 0
      integer :: nlsrerec = 0
! fename = file name for electromagnetic diagnostic
      character(len=32) :: fename = 'vpotrk2.0'
! define namelist
      namelist /em2d/ idrun, indx, indy, nte, modesxe, modesye, psolve, &
     &ndim, omx, omy, omz, ci, nerec, t0, tend, dt, ceng, indian,       &
     &rlprec, inprec, fename
!
! Namelist output for ion density diagnostic
! ndrec = current record number for ion density writes
      integer :: ndrec = 0
! fdname = file name for ion density diagnostic
      character(len=32) :: fdname = 'denk2.0'
! define namelist
      namelist /den2d/ idrun, indx, indy, ntd, modesxd, modesyd, psolve,&
     & ndrec, t0, tend, dt, ceng, indian, rlprec, inprec, fdname 
!
! Namelist output for ion current diagnostic
! njrec = current record number for ion current writes
      integer :: njrec = 0
! fjname = file name for ion current diagnostic
      character(len=32) :: fjname = 'vcurk2.0'
! define namelist
      namelist /vcur2d/ idrun, indx, indy, ntj, modesxj, modesyj,       &
     &psolve, ndim, omx, omy, omz, ci, njrec, t0, tend, dt, ceng,       &
     &indian, rlprec, inprec, fjname 
!
      end module input2d
