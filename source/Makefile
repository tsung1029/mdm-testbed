#Makefile for 2D PIC codes in new_beps2.source

GOBJS = nullgks2.o nullgks1.o
#CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile Absoft compiler with MacOS X

#FC90 = f90
#FC77 = f77
#CC = gcc

#OPTS90 = -O3 -N113 -YEXT_NAMES=LCS -YEXT_SFX=_
#OPTS77 = -O -N113 -f -N15
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s -N11
#LOPTS = -plainappl
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o

# Mac graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libmcX.o
#LIBS = $(CARBON)
# No graphics
#LIBS = $(CARBON)

# Makefile Nag compiler with MacOS X

#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS = -framework carbon 
#NOSAVE =
#LEGACY = -dusty

#MPOBJS = LnxMP.o

# No graphics
#LIBS =

# Makefile IBM compiler with MacOS X

#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qextname
#OPTS77 = -O3 -qautodbl=dbl4 -qextname -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#NOSAVE = -qnosave
#LEGACY =

#MPOBJS = LnxMP.o

# No graphics
#LIBS = $(CARBON)

# Makefile g95 compiler with MacOS X

#FC90 = g95
#FC77 = g95
#CC = gcc

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -fstatic
#MBOPTS = -fstatic
#LOPTS =
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS = $(CARBON) -lSystemStubs

# Makefile gfortran compiler with MacOS X

FC90 = gfortran
FC77 = gfortran
CC = gcc

OPTS90 = -O3 -fdefault-real-8 -ffree-form -fbounds-check
OPTS77 = -O3 -fdefault-real-8 -fbounds-check
CCOPTS = -O
MOPTS = -fno-automatic
MBOPTS = -fno-automatic
LOPTS =
NOSAVE = -frecursive
LEGACY =

MPOBJS = LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = dlibgks2.o libgks2.o libgks1.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# TIFF graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o librstr.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# ncar graphics (single precision only) from http://ngwww.ucar.edu/
#GOBJS = dlibgks2.o libgks2.o libgks1.o ncarstub.o
#LIBS = $(CARBON) -lSystemStubs -L/$(NCARG_ROOT)/lib -lncarg_gks -lncarg_c \
#-L/usr/X11R6/lib -lX11
# Tektronix graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs 
# No graphics
#LIBS = $(CARBON) -lSystemStubs

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#NOSAVE = -automatic
#LEGACY =

#MPOBJS = LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = $(CARBON) -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS = $(CARBON)

# Makefile Intel compiler with Linux

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#NOSAVE = -automatic
#LEGACY =

#MPOBJS = LnxMP.o LPProcessors.o

# No graphics
#LIBS =

# Makefile gfortran compiler with Linux

#FC90 = gfortran
#FC77 = gfortran
#CC = gcc

#OPTS90 = -O3 -fdefault-real-8 -ffree-form -fbounds-check
#OPTS77 = -O3 -fdefault-real-8 -fbounds-check
#CCOPTS = -O
#MOPTS = -fno-automatic
#MBOPTS = -fno-automatic
#LOPTS = -lpthread
#NOSAVE = -frecursive
#LEGACY =

#MPOBJS = LnxMP.o LPProcessors.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
#LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# Postcript printer
#GOBJS = dlibgks2.o libgks2.o libgks1.o libpsp.o libloc1.o
#LIBS =
# TIFF graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o librstr.o libloc1.o
#LIBS =
# Tektronix graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS = 
# No graphics
#LIBS =

# Makefile Pathscale compiler with AMD Opteron and Linux

#FC90 = pathf90
#FC77 = pathf90
#CC = gcc

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -static-data
#MBOPTS = -static-data
#LOPTS = -lpthread
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o LPProcessors.o

# No graphics
#LIBS =

# Makefile PGI compiler with AMD Opteron and Linux

#FC90 = pgf90
#FC77 = pgf90
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -Msave
#MBOPTS = -Msave
#LOPTS = -lpthread
#NOSAVE =
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

# No graphics
#LIBS =

#

GESOBJS = espush2mod.o simul2mod.o

MGESOBJS = mespush2mod.o msimul2mod.o

ESOBJS = globals.o input2mod.o init2mod.o diag2mod.o rbpush2mod.o \
bpush2mod.o push2mod.o fft2mod.o field2mod.o mp0mod.o init2lib.o \
rbpush2lib.o bpush2lib.o push2lib.o fft2lib.o field2lib.o diag2lib.o

MESOBJS = mrbpush2mod.o mbpush2mod.o mpush2mod.o mfft2mod.o \
mrbpush2lib.o mbpush2lib.o mpush2lib.o mfft2lib.o

GEMOBJS = empush2mod.o emsimul2mod.o simul2mod.o

MGEMOBJS = mempush2mod.o memsimul2mod.o msimul2mod.o

EMOBJS = rdpush2mod.o dpush2mod.o rdpush2lib.o dpush2lib.o

MEMOBJS = mrdpush2mod.o mdpush2mod.o mrdpush2lib.o mdpush2lib.o

GHDOBJS = ehdpush2mod.o hdsimul2mod.o

MGHDOBJS = mehdpush2mod.o mhdsimul2mod.o

HDOBJS = rhdpush2mod.o hdpush2mod.o rhdpush2lib.o hdpush2lib.o

MHDOBJS = mrhdpush2mod.o mrhdpush2lib.o mhdpush2mod.o mhdpush2lib.o

DESOBJS = dfield2mod.o bfield2mod.o cfield2mod.o nfield2mod.o \
mfield2mod.o hfield2mod.o npfield2mod.o vfield2mod.o fft12mod.o \
dfield2lib.o bfield2lib.o cfield2lib.o nfield2lib.o mfield2lib.o \
hfield2lib.o vfield2lib.o fft1lib.o

NPOBJS = nullMP.o

###
VPMLJS = vpml2mod.o
###

# Linkage rules

all : basic

basic : new_beps2 new_bbeps2_pei new_dbeps2 new_hdbeps2

bounded : new_d0_beps2 new_d0_bbeps2 new_d0_dbeps2

gridless: new_beps2gl new_bbeps2gl

threaded : threadbasic

threadbasic : new_mbeps2 new_mbbeps2_pei new_mdbeps2 new_mhdbeps2

threadbounded : new_d0_mbeps2 new_d0_mbbeps2 new_d0_mdbeps2

threadgridless : new_mbeps2gl new_mbbeps2gl

new_beps2 : new_beps2.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_beps2 \
        new_beps2.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

###
new_bbeps2_mdm : new_bbeps2_mdm.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(OPTS90) $(LOPTS) -o new_bbeps2_mdm \
        new_bbeps2_mdm.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS) $(LIBS)\
        $(VPMLJS)
###

new_bbeps2_pei : new_bbeps2_pei.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(OPTS90) $(LOPTS) -o new_bbeps2_pei \
        new_bbeps2_pei.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_dbeps2 : new_dbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_dbeps2 \
        new_dbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_hdbeps2 : new_hdbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(GHDOBJS) $(HDOBJS) \
                  $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(OPTS90) $(LOPTS) -o new_hdbeps2 \
        new_hdbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(GHDOBJS) $(HDOBJS) $(NPOBJS) \
        $(GOBJS) $(LIBS)

new_d0_beps2 : new_d0_beps2.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(NPOBJS) \
                   $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_beps2 \
        new_d0_beps2.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_d0_bbeps2 : new_d0_bbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) \
                    $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_bbeps2 \
        new_d0_bbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) $(NPOBJS) \
        $(GOBJS) $(LIBS)

new_d0_dbeps2 : new_d0_dbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) \
                    $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_dbeps2 \
        new_d0_dbeps2.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(DESOBJS) $(NPOBJS) \
        $(GOBJS) $(LIBS)

new_beps2gl : new_beps2gl.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_beps2gl \
        new_beps2gl.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_bbeps2gl : new_bbeps2gl.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(OPTS90) $(LOPTS) -o new_bbeps2gl \
        new_bbeps2gl.o $(ESOBJS) $(GEMOBJS) $(EMOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_mbeps2 : new_mbeps2.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mbeps2 \
        new_mbeps2.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_mbbeps2_pei : new_mbbeps2_pei.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                  $(MEMOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mbbeps2_pei \
        new_mbbeps2_pei.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(MPOBJS) $(GOBJS) $(LIBS)

new_mdbeps2 : new_mdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                  $(MEMOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mdbeps2 \
        new_mdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(MPOBJS) $(GOBJS) $(LIBS)

new_mhdbeps2 : new_mhdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                   $(MEMOBJS) $(MGHDOBJS) $(MHDOBJS) $(HDOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mhdbeps2 \
        new_mhdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(MGHDOBJS) $(MHDOBJS) $(HDOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_d0_mbeps2 : new_d0_mbeps2.o $(MGESOBJS) $(ESOBJS) $(DESOBJS) $(MESOBJS) \
                    $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_mbeps2 \
        new_d0_mbeps2.o $(MGESOBJS) $(ESOBJS) $(DESOBJS) $(MESOBJS) $(MPOBJS) \
	$(GOBJS) $(LIBS)

new_d0_mbbeps2 : new_d0_mbbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                     $(MEMOBJS) $(DESOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_mbbeps2 \
        new_d0_mbbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(DESOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_d0_mdbeps2 : new_d0_mdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                     $(MEMOBJS) $(DESOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_mdbeps2 \
        new_d0_mdbeps2.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(DESOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_mbeps2gl : new_mbeps2gl.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mbeps2gl \
        new_mbeps2gl.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPOBJS) $(GOBJS) $(LIBS)

new_mbbeps2gl : new_mbbeps2gl.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) \
                    $(MEMOBJS) $(MPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mbbeps2gl \
        new_mbbeps2gl.o $(ESOBJS) $(MESOBJS) $(MGEMOBJS) $(EMOBJS) $(MEMOBJS) \
        $(MPOBJS) $(GOBJS) $(LIBS)

# Compilation rules

LnxMP.o : LnxMP.c
	$(CC) $(CCOPTS) -c LnxMP.c

LPProcessors.o : LPProcessors.c
	$(CC) $(CCOPTS) -c LPProcessors.c

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f

libpsp.o : libpsp.f
	$(FC77) $(OPTS77) -c libpsp.f

librstr.o : librstr.f
	$(FC77) $(OPTS77) -c librstr.f

ncarstub.o : ncarstub.f
	$(FC77) $(OPTS77) -c ncarstub.f

libplt10.o : libplt10.f
	$(FC77) $(OPTS77) -c libplt10.f

plot10.o : plot10.f
	$(FC77) $(OPTS77) -c plot10.f

libt1.o : libt1.f
	$(FC77) $(OPTS77) -c libt1.f

libloc1.o : libloc1.f
	$(FC77) $(OPTS77) -c libloc1.f
	
dlibgks1.o : dlibgks1.f
	$(FC77) $(OPTS77) -c dlibgks1.f

dlibgks2.o : dlibgks2.f
	$(FC77) $(OPTS77) -c dlibgks2.f

libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

libgks2.o : libgks2.f
	$(FC77) $(OPTS77) -c libgks2.f

nullgks1.o : nullgks1.f
	$(FC77) $(OPTS77) -c nullgks1.f

nullgks2.o : nullgks2.f
	$(FC77) $(OPTS77) -c nullgks2.f

nullMP.o : nullMP.f
	$(FC77) $(OPTS77) -c nullMP.f

init2lib.o : init2lib.f
	$(FC77) $(OPTS77) -c init2lib.f

push2lib.o : push2lib.f
	$(FC77) $(OPTS77) -c push2lib.f

mpush2lib.o : mpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mpush2lib.f

bpush2lib.o : bpush2lib.f
	$(FC77) $(OPTS77) -c bpush2lib.f

mbpush2lib.o : mbpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mbpush2lib.f

dpush2lib.o : dpush2lib.f
	$(FC77) $(OPTS77) -c dpush2lib.f

mdpush2lib.o : mdpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mdpush2lib.f

hdpush2lib.o : hdpush2lib.f
	$(FC77) $(OPTS77) -c hdpush2lib.f

mhdpush2lib.o : mhdpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mhdpush2lib.f

rdpush2lib.o : rdpush2lib.f
	$(FC77) $(OPTS77) -c rdpush2lib.f

mrdpush2lib.o : mrdpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mrdpush2lib.f

rhdpush2lib.o : rhdpush2lib.f
	$(FC77) $(OPTS77) -c rhdpush2lib.f

mrhdpush2lib.o : mrhdpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mrhdpush2lib.f

rbpush2lib.o : rbpush2lib.f
	$(FC77) $(OPTS77) -c rbpush2lib.f

mrbpush2lib.o : mrbpush2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mrbpush2lib.f

fft2lib.o : fft2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c fft2lib.f

mfft2lib.o : mfft2lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mfft2lib.f

field2lib.o : field2lib.f
	$(FC77) $(OPTS77) -c field2lib.f

dfield2lib.o : dfield2lib.f
	$(FC77) $(OPTS77) -c dfield2lib.f

bfield2lib.o : bfield2lib.f
	$(FC77) $(OPTS77) -c bfield2lib.f

cfield2lib.o : cfield2lib.f
	$(FC77) $(OPTS77) -c cfield2lib.f

nfield2lib.o : nfield2lib.f
	$(FC77) $(OPTS77) -c nfield2lib.f

mfield2lib.o : mfield2lib.f
	$(FC77) $(OPTS77) -c mfield2lib.f

hfield2lib.o : hfield2lib.f
	$(FC77) $(OPTS77) -c hfield2lib.f

vfield2lib.o : vfield2lib.f
	$(FC77) $(OPTS77) $(LEGACY) -c vfield2lib.f

diag2lib.o : diag2lib.f
	$(FC77) $(OPTS77) -c diag2lib.f

fft1lib.o : fft1lib.f
	$(FC77) $(OPTS77) -c fft1lib.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

input2mod.o : input2mod.f globals.o
	$(FC90) $(OPTS90) -c input2mod.f

init2mod.o : init2mod.f input2mod.o
	$(FC90) $(OPTS90) -c init2mod.f

mpush2mod.o : mpush2mod.f push2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpush2mod.f

push2mod.o : push2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c push2mod.f

mbpush2mod.o : mbpush2mod.f bpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mbpush2mod.f

bpush2mod.o : bpush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c bpush2mod.f

mdpush2mod.o : mdpush2mod.f dpush2mod.o
	$(FC90) $(OPTS90) -c mdpush2mod.f

dpush2mod.o : dpush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c dpush2mod.f

mhdpush2mod.o : mhdpush2mod.f hdpush2mod.o
	$(FC90) $(OPTS90) -c mhdpush2mod.f

hdpush2mod.o : hdpush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c hdpush2mod.f

mrdpush2mod.o : mrdpush2mod.f rdpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mrdpush2mod.f

rdpush2mod.o : rdpush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c rdpush2mod.f

mrhdpush2mod.o : mrhdpush2mod.f rhdpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mrhdpush2mod.f

rhdpush2mod.o : rhdpush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c rhdpush2mod.f

mrbpush2mod.o : mrbpush2mod.f rbpush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mrbpush2mod.f

rbpush2mod.o : rbpush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c rbpush2mod.f

mspush2mod.o : mspush2mod.f spush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mspush2mod.f

msfield2mod.o : msfield2mod.f spush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c msfield2mod.f

mtpush2mod.o : mtpush2mod.f spush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mtpush2mod.f

mtfield2mod.o : mtfield2mod.f spush2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mtfield2mod.f

spush2mod.o : spush2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c spush2mod.f

mfft2mod.o : mfft2mod.f fft2mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mfft2mod.f

fft2mod.o : fft2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c fft2mod.f

espush2mod.o : espush2mod.f rbpush2mod.o bpush2mod.o push2mod.o fft2mod.o
	$(FC90) $(OPTS90) -c espush2mod.f

mespush2mod.o : mespush2mod.f mrbpush2mod.o mbpush2mod.o mpush2mod.o mfft2mod.o
	$(FC90) $(OPTS90) -c mespush2mod.f

empush2mod.o : empush2mod.f rbpush2mod.o rdpush2mod.o bpush2mod.o dpush2mod.o \
               push2mod.o fft2mod.o
	$(FC90) $(OPTS90) -c empush2mod.f

mempush2mod.o : mempush2mod.f mrbpush2mod.o mrdpush2mod.o mbpush2mod.o \
                mdpush2mod.o mpush2mod.o mfft2mod.o
	$(FC90) $(OPTS90) -c mempush2mod.f

ehdpush2mod.o : ehdpush2mod.f rhdpush2mod.o hdpush2mod.o
	$(FC90) $(OPTS90) -c ehdpush2mod.f

mehdpush2mod.o : mehdpush2mod.f mrhdpush2mod.o mhdpush2mod.o
	$(FC90) $(OPTS90) -c mehdpush2mod.f

field2mod.o : field2mod.f diag2mod.o
	$(FC90) $(OPTS90) $(LEGACY) -c field2mod.f

dfield2mod.o : dfield2mod.f globals.o
	$(FC90) $(OPTS90) -c dfield2mod.f

bfield2mod.o : bfield2mod.f globals.o
	$(FC90) $(OPTS90) -c bfield2mod.f

cfield2mod.o : cfield2mod.f globals.o
	$(FC90) $(OPTS90) -c cfield2mod.f

nfield2mod.o : nfield2mod.f globals.o
	$(FC90) $(OPTS90) -c nfield2mod.f

mfield2mod.o : mfield2mod.f globals.o
	$(FC90) $(OPTS90) -c mfield2mod.f

hfield2mod.o : hfield2mod.f globals.o
	$(FC90) $(OPTS90) -c hfield2mod.f

vfield2mod.o : vfield2mod.f globals.o
	$(FC90) $(OPTS90) -c vfield2mod.f

npfield2mod.o : npfield2mod.f field2mod.o dfield2mod.o bfield2mod.o \
                cfield2mod.o nfield2mod.o mfield2mod.o hfield2mod.o
	$(FC90) $(OPTS90) -c npfield2mod.f

diag2mod.o : diag2mod.f init2mod.o
	$(FC90) $(OPTS90) -c diag2mod.f

###
vpml2mod.o : vpml2mod.f
	$(FC90) $(OPTS90) -c vpml2mod.f
###

simul2mod.o : simul2mod.f diag2mod.o espush2mod.o field2mod.o
	$(FC90) $(OPTS90) -c simul2mod.f

msimul2mod.o : simul2mod.f diag2mod.o mespush2mod.o field2mod.o
	$(FC90) $(OPTS90) -o msimul2mod.o -c simul2mod.f

emsimul2mod.o : emsimul2mod.f simul2mod.o empush2mod.o
	$(FC90) $(OPTS90) -c emsimul2mod.f

memsimul2mod.o : emsimul2mod.f msimul2mod.o mempush2mod.o
	$(FC90) $(OPTS90) -o memsimul2mod.o -c emsimul2mod.f

hdsimul2mod.o : hdsimul2mod.f simul2mod.o emsimul2mod.o ehdpush2mod.o
	$(FC90) $(OPTS90) -c hdsimul2mod.f

mhdsimul2mod.o : hdsimul2mod.f msimul2mod.o memsimul2mod.o mehdpush2mod.o
	$(FC90) $(OPTS90) -o mhdsimul2mod.o -c hdsimul2mod.f

fft12mod.o : fft12mod.f globals.o diag2mod.o
	$(FC90) $(OPTS90) -c fft12mod.f

mp0mod.o : mp0mod.f
	$(FC90) $(OPTS90) -c mp0mod.f

new_beps2.o : new_beps2.f simul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_beps2.f

new_mbeps2.o : new_beps2.f msimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mbeps2.o -c new_beps2.f

new_dbeps2.o : new_dbeps2.f emsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_dbeps2.f

new_hdbeps2.o : new_hdbeps2.f hdsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_hdbeps2.f

new_mdbeps2.o : new_dbeps2.f memsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mdbeps2.o -c new_dbeps2.f

new_mhdbeps2.o : new_hdbeps2.f mhdsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mhdbeps2.o -c new_hdbeps2.f

####
new_bbeps2_mdm.o : new_bbeps2_mdm.f emsimul2mod.o mp0mod.o vpml2mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_bbeps2_mdm.f
###

new_bbeps2_pei.o : new_bbeps2_pei.f emsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_bbeps2_pei.f

new_mbbeps2_pei.o : new_bbeps2_pei.f memsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mbbeps2_pei.o -c new_bbeps2_pei.f

new_d0_beps2.o : new_d0_beps2.f simul2mod.o npfield2mod.o \
                 vfield2mod.o fft12mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_beps2.f

new_d0_mbeps2.o : new_d0_beps2.f msimul2mod.o npfield2mod.o \
                  vfield2mod.o fft12mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_d0_mbeps2.o -c new_d0_beps2.f

new_d0_bbeps2.o : new_d0_bbeps2.f emsimul2mod.o npfield2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -c new_d0_bbeps2.f

new_d0_mbbeps2.o : new_d0_bbeps2.f emsimul2mod.o npfield2mod.o \
                   mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -o new_d0_mbbeps2.o -c new_d0_bbeps2.f

new_d0_dbeps2.o : new_d0_dbeps2.f emsimul2mod.o npfield2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -c new_d0_dbeps2.f

new_d0_mdbeps2.o : new_d0_dbeps2.f emsimul2mod.o npfield2mod.o \
                   mp0mod.o
	$(FC90) $(OPTS90) $(MBOPTS) -o new_d0_mdbeps2.o -c new_d0_dbeps2.f

new_beps2gl.o : new_beps2gl.f simul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_beps2gl.f

new_mbeps2gl.o : new_beps2gl.f msimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mbeps2gl.o -c new_beps2gl.f

new_bbeps2gl.o : new_bbeps2gl.f emsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_bbeps2gl.f

new_mbbeps2gl.o : new_bbeps2gl.f memsimul2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mbbeps2gl.o -c new_bbeps2gl.f

new_sbeps2.o : new_sbeps2.f simul2mod.o mp0mod.o spush2mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_sbeps2.f

new_msbeps2.o : new_msbeps2.f msimul2mod.o mspush2mod.o msfield2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_msbeps2.f

new_mtbeps2.o : new_mtbeps2.f msimul2mod.o mtpush2mod.o mtfield2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_mtbeps2.f

new_mtbeps2_cu.o : new_mtbeps2_cu.f msimul2mod.o mtpush2mod.o mtfield2mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_mtbeps2_cu.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f all : new_beps2 new_bbeps2_pei new_dbeps2 new_hdbeps2 \
                new_d0_beps2 new_d0_bbeps2 new_d0_dbeps2 \
                new_beps2gl new_bbeps2gl new_bbeps2_mdm \
                new_mbeps2 new_mbbeps2_pei new_mdbeps2 new_mhdbeps2 \
                new_d0_mbeps2 new_d0_mbbeps2 new_d0_mdbeps2 \
                new_mbeps2gl new_mbbeps2gl
