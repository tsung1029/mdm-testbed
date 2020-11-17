# Makefile for spectrum2 and vspectrum2

GOBJS = nullgks2.o nullgks1.o gksnull.o
CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

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
#LEGACY =

# Mac graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libmcX.o
#LIBS = $(CARBON)
# No graphics
#LIBS =

# Makefile Nag compiler with MacOS X

#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY = -dusty

# No graphics
#LIBS =

# Makefile IBM compiler with MacOS X

#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

# No graphics
#LIBS =

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
#LEGACY =

# X11 graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS =  -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Tektronix
#GOBJS = dlibgks2.o libgks2.o libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS = -lSystemStubs
# No graphics
#LIBS = -lSystemStubs

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#FC77 = gfortran
#CC = gcc

#OPTS90 = -O3 -fdefault-real-8
#OPTS77 = -O3 -fdefault-real-8
#CCOPTS = -O
#MOPTS = -fno-automatic
#MBOPTS = -fno-automatic
#LOPTS =
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# X11 graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = dlibgks2.o libgks2.o libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# TIFF graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libgks2.o librstr.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# ncar graphics (single precision only) from http://ngwww.ucar.edu/
#GOBJS = dlibgks2.o libgks2.o libgks1.o libgks2.o ncarstub.o
#LIBS = $(CARBON) -lSystemStubs -L/$(NCARG_ROOT)/lib -lncarg_gks -lncarg_c \
#-L/usr/X11R6/lib -lX11
# Tektronix
#GOBJS = dlibgks2.o libgks2.o libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS = -lSystemStubs
# No graphics
#LIBS = -lSystemStubs

# Makefile Intel compiler with MacOS X

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS =
#LEGACY =

# X11 graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS =

# Makefile Intel compiler with Linux

#FC90 = ifc
#FC77 = ifc
#CC = gcc

#OPTS90 = -O3 -tpp7 -xW -r8
#OPTS77 = -O3 -tpp7 -xW -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS =
#LEGACY =

# No graphics
#LIBS =

# Makefile gfortran compiler with Linux

FC90 = gfortran
FC77 = gfortran
CC = gcc

OPTS90 = -O3 -fdefault-real-8
OPTS77 = -O3 -fdefault-real-8
CCOPTS = -O
MOPTS = -fno-automatic
MBOPTS = -fno-automatic
LOPTS = -lpthread
LEGACY =

MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
#LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# Postcript printer
#GOBJS = dlibgks2.o libgks2.o libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS =
# TIFF graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libgks2.o librstr.o libloc1.o
#LIBS =
# Tektronix
#GOBJS = dlibgks2.o libgks2.o libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS =
# No graphics
LIBS =

# Makefile Pathscale compiler with AMD Opteron and Linux

#FC90 = pathf90
#FC77 = pathf90
#CC = gcc

#OPTS90 = -O3 -r8 -static
#OPTS77 = -O3 -r8 -static
#CCOPTS = -O
#MOPTS = -static-data
#MBOPTS = -static-data
#LOPTS =
#LEGACY =

# No graphics
#LIBS =

#

OBJS = globals.o input2mod.o init2mod.o fft2mod.o field2mod.o \
diag2mod.o init2lib.o fft2lib.o fft1lib.o field2lib.o diag2lib.o \
libanls2.o nlibpars.o

# Linkage rule

all : spectrum2 vspectrum2

spectrum2 : spectrum2.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o spectrum2 spectrum2.o $(OBJS) \
        $(GOBJS) $(LIBS)

vspectrum2 : vspectrum2.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o vspectrum2 vspectrum2.o $(OBJS) \
        $(GOBJS) $(LIBS)

# Compilation rules

gksnull.o : gksnull.f
	$(FC77) $(OPTS77) -c gksnull.f

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f

libpsp.o : libpsp.f
	$(FC77) $(OPTS77) -c libpsp.f

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

libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

libgks2.o : libgks2.f
	$(FC77) $(OPTS77) -c libgks2.f

dlibgks2.o : dlibgks2.f
	$(FC77) $(OPTS77) -c dlibgks2.f

nullgks1.o : nullgks1.f
	$(FC77) $(OPTS77) -c nullgks1.f

nullgks2.o : nullgks2.f
	$(FC77) $(OPTS77) -c nullgks2.f

init2lib.o : init2lib.f
	$(FC90) $(OPTS77) -c init2lib.f

fft2lib.o : fft2lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c fft2lib.f

fft1lib.o : fft1lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c fft1lib.f

field2lib.o : field2lib.f
	$(FC90) $(OPTS77) -c field2lib.f

field1lib.o : field1lib.f
	$(FC90) $(OPTS77) -c field1lib.f

diag2lib.o : diag2lib.f
	$(FC90) $(OPTS77) $(LEGACY) -c diag2lib.f

libanls2.o : libanls2.f
	$(FC77) $(OPTS77) -c libanls2.f

nlibpars.o : nlibpars.f
	$(FC77) $(OPTS77) -c nlibpars.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

input2mod.o : input2mod.f globals.o
	$(FC90) $(OPTS90) -c input2mod.f

init2mod.o : init2mod.f input2mod.o
	$(FC90) $(OPTS90) -c init2mod.f

fft2mod.o : fft2mod.f diag2mod.o
	$(FC90) $(OPTS90) -c fft2mod.f

field2mod.o : field2mod.f globals.o
	$(FC90) $(OPTS90) $(LEGACY) -c field2mod.f

diag2mod.o : diag2mod.f init2mod.o
	$(FC90) $(OPTS90) -c diag2mod.f

spectrum2.o : spectrum2.f fft2mod.o field2mod.o diag2mod.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c spectrum2.f

vspectrum2.o : vspectrum2.f fft2mod.o field2mod.o diag2mod.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c vspectrum2.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f spectrum2 vspectrum2
