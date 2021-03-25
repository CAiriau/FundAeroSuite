# Conical flow

## Introduction

Calculate the oblique shock on a sharp cone can be rather difficult. The main equations
are obtained from an ODE system which required a great accuracy, that means a very small step size
for integration.

Therefore it can be quite long to get the shock curves for conical flows.
The Fortran code has been used to generate data. Then it is possible
to extract data and interpolate to solve any conical supersonic flow. 

The data extraction and interpolation are performed in Python (including in the **FundAeroSuite**).

## Compilation

A makefile has been written. For the first time,  the following command as to be written in a terminal:

\code{.sh}
make clean; make
\endcode

 `make clean` is only necessary if the source of the module is modified. It due to the manner the Makefile is written.

For a normal compilation of the code, just write:
\code{.sh}
make
\endcode

With the make file, the executable code is named `cone` and is moves in `~/bin/` directory. It can be modified by the user.

It can be compiled manually with:
\code{.sh}
gfortran cone.f90 -o choc_conique.bin
\endcode

A `readme.md` file is also put in the main conical flow directory which provides some comments for compilation.

## Run

The user have to go in `Conical_flow/run` directory and write in a terminal:

\code{.sh}
cone
\endcode

## Options

The code options are found in the file `cone.in`. An example of the file can be found 
*doc/example/cone.in* :

\code{.sh}
# input file to solve conical shocks
2               # option : 1: the Mach numbers are in a table found in the source file,   2 : table below
1               # number of Mach number (option 2), if = 0 , the table is build with a Mach variation  delta_Mach
0.01            # delta_Mach, Mach number step
2.00d0          # Initial Mach number (for a single Mach calculus it is this value) 
1.02d0          # Final Mach number
1001            # Number of points on the plot for a given Mach number  
1.d-2           # Accuracy to determine the wall deviation angle, for a given shock wave angle (in degrees)  
1               # 1 : Newton algorithm + RK4, 0 : Euler order 1 integration, 2 : old version of Newton
10001           # number of discrete points  - angle - between the shock line and the wall
1.d-5           # abs(u_theta) on the cone wall
0               # velocity option =>  1: tangential velocity is saved in a file  
\endcode

## Database

The database is available in the `Python/conical_shock` directory

## Caution

Because the conical flow code has been translated in english, the name of the output files have benn changed accordingly.

However, if a new database is generated, the new file name will have to be changed in the Python code where conical flow
database is used.