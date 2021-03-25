#SHELL = /bin/sh
SYSTEM=gfortran
 

ifeq ($(SYSTEM),gfortran_imft)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   FFLAGS= -fdefault-real-8 -Wsurprising  -ffree-line-length-none -O3 -Wall 
   LAPACK = /usr/lib64/liblapack.so.3   # -L/usr/lib -llapack  
endif

ifeq ($(SYSTEM),gfortran)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   FFLAGS=  -Wsurprising  -ffree-line-length-none -O3 -Wall 
   #FFLAGS=  -Wsurprising  -ffree-line-length-none -O3 -Wall 
   #LAPACK = /usr/lib/liblapack-3.so
   LAPACK = /usr/lib/lapack/liblapack.so.3gf
endif

ifeq ($(SYSTEM),ifort_imft)
#  this is for the intel  Fortran 90 compiler (
   FC= ifort
#   FFLAGS=   -i8 -r8  -WB  -assume buffered_io 
#    FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel -check bound -g
#   ! option for hyperion
    FFLAGS=    -r8  -WB  -assume buffered_io -shared-intel
    LAPACK = -lmkl_intel -lmkl_sequential -lmkl_core -lmkl_lapack -lguide -lpthread
endif


BIN=$(HOME)/bin
EXENAME1=$(BIN)/cone

OBJS1 =  cone.o
ALL: $(EXENAME1)

.SUFFIXES: .o .f90


.f90.o:
	$(FC)  -c $(FFLAGS)  $<

$(EXENAME1):   $(OBJS1)
	$(FC) $(FFLAGS) -o $@ $(OBJS1) 

clean:
	/bin/rm *.o *.mod *.f90~  Makefile~ *.*~   $(EXENAME1)

 

