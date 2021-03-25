#SHELL = /bin/sh
 SYSTEM=gfortran
 

ifeq ($(SYSTEM),gfortran_imft)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   FFLAGS= -fdefault-real-8 -Wsurprising  -ffree-line-length-none -O3 -Wall 
   LAPACK = /usr/lib/liblapack.so.3 # -L/usr/lib -llapack #/usr/lib/liblapack-3.so
endif

ifeq ($(SYSTEM),gfortran)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   FFLAGS= -fdefault-real-8 -Wsurprising  -ffree-line-length-none -O3 -Wall 
   FFLAGS=  -Wsurprising  -ffree-line-length-none -O3 -Wall 
   #FFLAGS=  -Wsurprising  -ffree-line-length-none -O3 -Wall 
   #LAPACK = /usr/lib/liblapack-3.so
   LAPACK = /usr/lib/liblapack.so.3gf
endif


ifeq ($(SYSTEM),ifort)
#  this is for the intel  Fortran 90 compiler (
   FC= ifort
#   FFLAGS=   -i8 -r8  -WB  -assume buffered_io 
   FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel 
   LAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack -lguide
endif


BIN=$(HOME)/bin
EXENAME=$(BIN)/hypersonic_Kp              # Kp in hypersonic flow 
EXENAME1=$(BIN)/atmospheric_reentry       # atmospheric reentry in hypersonic flow
EXENAME2=$(BIN)/piston_method             # piston method in hypersonic flow

OBJS =  hypersonic_Kp.o
ALL: $(EXENAME)

OBJS1 = atmospheric_reentry.o
ALL: $(EXENAME1)

OBJS2 = piston_method.o
ALL: $(EXENAME2)

.SUFFIXES: .o .f90

.f90.o:
	$(FC)  -c $(FFLAGS)  $<

$(EXENAME): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS) 

$(EXENAME1): $(OBJS1)
	$(FC) $(FFLAGS) -o $@ $(OBJS1) 

#$(EXENAME2): $(OBJS2)
#	$(FC) $(FFLAGS) -o $@ $(OBJS2) 

clean:
	/bin/rm *.o *.mod *.f90~  Makefile~ *.*~ *.out *.bin beta  $(EXENAME)  $(EXENAME1) $(EXENAME2)

