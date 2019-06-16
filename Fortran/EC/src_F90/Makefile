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

EXENAME=$(BIN)/EC


OBJFS = atmos.o mod_read.o  mod_EC.o  exercices_chap10.o exercices_chap11.o exercices_chap12.o mod_interaction_choc.o mod_Etude.o  main_EC.o
ALL:  $(EXENAME)

.SUFFIXES: .o .f90

file=atmos.o  mod_read.o mode_EC.o exercices_chap10.o exercices_chap11.o exercices_chap12.o  mod_interaction_choc.o mod_Etude.o 
OBJ_mod =  $(file).o $(file).mod 

.f90.o: atmos.o  mod_read.o mode_EC.o exercices_chap10.o exercices_chap11.o exercices_chap12.o  mod_interaction_choc.o mod_Etude.o 
	$(FC)  -c $(FFLAGS)  $<

$(EXENAME):  $(OBJFS)
	$(FC) $(FFLAGS) -o $@ $(OBJFS) 

exercices_chap10.o: 
	$(FC) $(FFLAGS) -c exercices_chap10.f90

exercices_chap11.o: 
	$(FC) $(FFLAGS) -c exercices_chap11.f90

exercices_chap12.o: 
	$(FC) $(FFLAGS) -c exercices_chap12.f90

mod_read.o: 
	$(FC) $(FFLAGS) -c mod_read.f90

atmos.o:
	$(FC) $(FFLAGS) -c atmos.f90

mod_EC.o: 
	$(FC) $(FFLAGS) -c mod_EC.f90

mod_interaction_choc.o: 
	$(FC) $(FFLAGS) -c mod_interaction_choc.f90

mod_Etude.o: 
	$(FC) $(FFLAGS) -c mod_Etude.f90

clean_mod:
	/bin/rm  $(OBJ_mod)

clean:
	rm -f   *~ *.*~  *.mod *.o $(EXENAME) 


