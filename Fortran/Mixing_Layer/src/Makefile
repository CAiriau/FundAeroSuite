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



ifeq ($(SYSTEM),pgf90)
#  this is for the Portland Group Fortran 90 compiler (tested with version 7.2-5)
   FC= pgf90
   #FFLAGS=  -i8 -r8
   FFLAGS=  -r8
   LAPACK = -L/PRODCOM/PGI/linux86/10.9/lib -llapack -lblas
endif

ifeq ($(SYSTEM),ifort)
#  this is for the intel  Fortran 90 compiler (
   FC= ifort
#   FFLAGS=   -i8 -r8  -WB  -assume buffered_io 
   FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel 
   # FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel -check bound -g
#   ! option for hyperion
#    FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel
   LAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack -lguide
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
EXENAME2=$(BIN)/mixing_layer  
EXENAME3=$(BIN)/mixing_layer_o1     # mixing layer with  Uref=(U1+U2)/2

OBJS2 = mixing_layer.o 
ALL: $(EXENAME2)

OBJS3 = mixing_layer_o1.o 
ALL: $(EXENAME3)

.SUFFIXES: .o .f90

.f90.o:
	$(FC)  -c $(FFLAGS)  $<

$(EXENAME2):   $(OBJS2)
	$(FC) $(FFLAGS) -o $@ $(OBJS2) 

$(EXENAME3):   $(OBJS3)
	$(FC) $(FFLAGS) -o $@ $(OBJS3) 

clean:
	/bin/rm *.o *.mod *.f90~  Makefile~ *.*~  $(EXENAME3)  $(EXENAME2) 
# /bin/rm  $(EXENAME3) 

