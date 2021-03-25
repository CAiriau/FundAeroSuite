# Busemann polar

## Introduction

The Busemann polar is the central curve to solve an oblique shock since it provides the shock angle depending on the 
deviation angle and the Mach number.

The curves are find in any aerodynamics or fluid mechanic books. The program describes here is able to generate the 
database (files).

The files can be used later with any plot tool (*xmgrace, python, matlab, ...,* ).

## Compilation

* A makefile has been written. The compilation is carried out in a terminal :

\code{.sh}
make
\endcode

* `make clean` is used to clean the executable.

* With the make file, the executable code is named `busemann` and is moves in `~/bin/` directory. It can be modified by the user.

* It can be compiled manually with:
\code{.sh}
gfortran busemann.f90 -o busemann.bin
\endcode

## Run

The user have to go in `Polar_oblique_shock/run` directory and write in a terminal:

\code{.sh}
busemann
\endcode

An help appears.

## Options

The user has to enter two parameters:

1. `Minit`, the Mach for the first generation of the polar.
It must be choosen in the following values :

    1.01      1.03      1.10      1.15      1.20      1.25      1.30      1.35      1.40      1.45
    1.50      1.60      1.70      1.80      2.00      2.20      2.40      2.60      2.80      3.00
    3.20      3.40      3.60      3.80      4.00      4.50      5.00      6.00     10.00    100.00

    This values can be changed in the source code

2. `n`,  the number of polar to plot.

* For instance with `Mint= 2` and `n=4`  the polar for Mach numbers 2, 2.2, 2.4 and 2.6 will be generated

## Outputs and database

As outputs we get:

* files named as `data_Mijk.dat` where **ijk= 100 x Mach**
* the file `unitary_downstream_Mach.dat` which provides for each Mach number, the deviation angle and the shock angle such that the downstream Mach number is equal to 1
* the file `theta_max.dat` which contains the limit of the plot in the (theta,sigma) plane.

## Caution