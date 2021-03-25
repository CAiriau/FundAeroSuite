# Naca airfoil

## Introduction

In many aerodynamic application an airfoil geometry is required. One can
find many web sites to get airfoil geometry. The best one is http://airfoiltools.com/

Ib this suite is provide two programs to calculate NACA 4 and 5 digits :

* the first one **naca.f90** is originated from Marc Drela f77 subroutines translated by me in f90
* the second one **naca45_airfoil.f90** is my personal application using fortran type and module.

One can choose its application !

Note:
    a matlab and a python versions can be found in the FundAeroSuite on GitHub.

## Compilation

* A makefile has been written. The compilation is carried out in a terminal :

\code{.sh}
make
\endcode

* `make clean` is used to clean the executable.

* With the make file, the both executable codes is named `naca` and `naca45` are moved  in `~/bin/` directory. It can be modified by the user.

* They can also be compiled using the shell script `x1.run.sh` with option 1 and 5

## Run

* The module can be used as this in a fortran program : 

\code{f90}
!=======================
PROGRAM Airfoil_Design
!=======================
use naca45_airfoil
implicit none
call run_airfoil
\endcode

or directly by using the module integrated function as :
\code{f90}
call NACA45(0012,101,0)
call SaveAirfoil
call deallocate_airfoil
\endcode

for a NACA 0012 airfoil for instance

* The user have to go in `NacaAirfoil/run` directory (or anywhere) and write in a terminal:

\code{.sh}
naca45           # or naca
\endcode

* A  shell script `x1.run.sh`, also located in this directory can be used, to run, to plot and to clean the run directory.

* script shell `gp` can be used to plot the airfoil. It is a wrapper of `gnuplot`. Run `gp` in a terminal to see the help