# Falkner-skan flow

## Introduction

### Present code

* The boundary layer along a flat plate in incidence is similar to a boundary developed on a wedge.

* An autosimilar solution is possible and it is solution of the ODE (following Cousteix's notation ) :

\f$\displaystyle f'''+ f~f''+ \beta~(1-f'^2) = 0  \f$

with \f$ \displaystyle f'(\eta) = u(x,y)/u_e(x) \f$ 
and  \f$ \displaystyle \eta = y /\delta \f$
where
\f$ \displaystyle  \delta = \sqrt{\frac{\nu x}{Ue} }  \sqrt{\frac{2}{m+1}} \f$

\f$ \beta \f$ is associated to the streawise pressure gradient. 

* The boundary conditions are  \f$ f(0)=f'(0)=0 \f$ and \f$ \displaystyle \lim_{\eta \to \infty} f' =1 \f$

* This system is solved by a Newton algorithm with parameter \f$ s=f'(0) \f$ and the value of 
\f$ \eta_{\max} \f$  where \f$ u(\eta_{\max})= 1 \f$ is the target.

* The solution when \f$ \beta \f$ varies is given by an arc-length continuation method.

* \f$ \eta, f, f', f'' \f$ and boundary layer characteristics are given in output files.

* Many programs can be found in internet where FSK boundary layer is solved, but none (from my knowledge) takes into account of the variation of the 
\f$ \eta_{\max}\f$ parameter with respect to \f$ \beta \f$

### Textbook

- in our book "Aérodynamique fondamentale", Giovannini, Airiau, 2016 is solved:

\f$\displaystyle f'''+\frac{m+1}{2} f f''+ m (1-f'^2)=0  \f$ with \f$\displaystyle \delta = \sqrt{\frac{\nu x}{Ue} } \f$
where \f$ \displaystyle m = \frac{\beta}{2-\beta} \f$

* There is a scale factor between the solution of the both ODE
* In the program, the variable *xk* is this scale factor, and it is written in the file **save. characteristics.dat**
* A table with the main boundary layer parameters can be found in our text and exercises books.

## Compilation

* A makefile has been written. The compilation is carried out in a terminal :

\code{.sh}
make
\endcode

* `make clean` is used to clean the executable.

* With the make file, the  executable codes is named `fsk` and is moved  in `~/bin/` directory. It can be modified by the user.

* it can also be compiled using the shell script

\code{.sh}
x1.run.sh c
\endcode

## Run

\code{.sh}
fsk
\endcode

* A  shell script  also located in this directory can be used, to run and to clean the run directory:


\code{.sh}
x1.run.sh r
\endcode

## Options and parameters

Many parameters have to be given. 

* The first way is to run the code with the `read_option=1` as shown below 
and default parameters written in the main module are used.

\code{f90}
PROGRAM test_fsk
! read_option = 1 : parameters read in input file
!             = 0 : use default parameters defined in the module
use fsk_module
integer read_option

read_option = 1
call Main_Continuation(read_option)

END PROGRAM test_fsk
\endcode

* By setting `read_option=0` the input file given below and named `fsk_continuation` is read if it exists.
If it is not present, the default parameters are used.

An exemple of input file is  :

\code{.sh}
# Input file for the module fsk (Falkner-Skan boundary layer)
1               # flag  = 1 : the parameters are read in this file,  0 : default parameters
1               # flag_continuation = 1 :  FSK (single value), else continuation method to follow varying beta
0.0             # first beta value
0.00500d0       # step in beta =>  the number of point with beta is determined from L_max (see below)
1.00d-4         # deta : step in eta
0.10d-2         # dL  : arc length in the continuation method
0.5d0           # s_init =f''(0), initial guess for the Newton method 
101             # Max: number of beta for the continuation method
11.0d0          # eta_max approximate to 10 , close to separation set between 20 to 30
1               # integration scheme, =  O : Euler order 1, 1 : Runge-Kutta order 4
# end
\endcode

## Database

* To generate a database the continuation method has to be used.

* It can be tested just by setting `flag_continuation=0` in the previous file.

## Caution

Close to separation it can be difficult to get convergence. But the paramaters given here are efficient.
