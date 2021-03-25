SYSTEM=gfortran

ifeq ($(SYSTEM),gfortran_laptop)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   FFLAGS= -fdefault-real-8 -Wsurprising  -ffree-line-length-none
# -O3 -Wall
   LAPACK = /usr/lib/liblapack.so

endif
ifeq ($(SYSTEM),gfortran_imft)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   FFLAGS= -fdefault-real-8 -Wsurprising  -ffree-line-length-none -O3 -Wall 
   LAPACK = /usr/lib/liblapack.so.3 # -L/usr/lib -llapack #/usr/lib/liblapack-3.so
endif
ifeq ($(SYSTEM),gfortran)
#  gfortran compiler (part of gcc 4.x.x)
   FC= gfortran
   #FFLAGS= -fdefault-real-8 -Wsurprising  -ffree-line-length-none -O3 -Wall 
   FFLAGS=   -ffree-line-length-none -O3 -Wall 
   LAPACK = /usr/lib/liblapack.so.3gf
endif
ifeq ($(SYSTEM),ifort)
#  this is for the intel  Fortran 90 compiler (
   FC= ifort
#   FFLAGS=   -i8 -r8  -WB  -assume buffered_io 
#    FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel -check bound -g
#   ! option for hyperion
   FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel
   LAPACK = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack -lguide
endif
ifeq ($(SYSTEM),ifort_imft)
#  this is for the intel  Fortran 90 compiler 
   FC= ifort
#   FFLAGS=   -i8 -r8  -WB  -assume buffered_io 
#    FFLAGS=    -r8  -WB  -assume buffered_io -mcmodel medium -shared-intel -check bound -g
#   
   FFLAGS=    -r8  -WB  -assume buffered_io -shared-intel
   LAPACK = -lmkl_intel -lmkl_sequential -lmkl_core -lmkl_lapack -lguide -lpthread
endif
BIN=$(HOME)/bin

EXENAME=$(BIN)/naca
EXENAME1=$(BIN)/naca45


OBJFS =naca.o
ALL:  $(EXENAME)

OBJFS1 =naca45_airfoil.o
ALL:  $(EXENAME1)

.SUFFIXES: .o .f90


.f90.o:  
	$(FC)  -c $(FFLAGS)  $<

$(EXENAME):  $(OBJFS)
	$(FC) $(FFLAGS) -o $@ $(OBJFS) 

$(EXENAME1):  $(OBJFS1)
	$(FC) $(FFLAGS) -o $@ $(OBJFS1) 

clean_mod:
	/bin/rm  $(OBJ_mod)

clean:
	rm -f   *~ *.*~  *.mod *.o $(EXENAME) $(EXENAME1)


