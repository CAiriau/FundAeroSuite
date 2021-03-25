# Hypersonic flow

## Introduction

The main exercises on the hypersonic regime are written in Python.
However, for the texbook "Fondamentale Aérodynamique" three short codes have been written in Fortran to
generate plots about :

* atmospheric reentry parameters
* Pressure coefficient in hypersonic flow
* piston method

## Compilation

* A makefile has been written. The compilation is carried out in a terminal :

\code{.sh}
make
\endcode

* `make clean` is used to clean the executable.

* With the make file, the three executable codes named `atmospheric_reentry` and `hypersonic_Kp` and
`piston_method` are moved  in `~/bin/` directory. It can be modified by the user.

## Run

* The user have to go in `Hypersonic/run` directory (or anywhere)  write in a terminal the code name
