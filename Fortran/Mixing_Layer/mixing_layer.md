# Mixing-layer flow

## Introduction

The  laminar or turbulent incompressible 2D mixing layer is solved with a autosimilar solution.


* The governing equation is an ODE :

  \f$ \displaystyle \frac{1}{2}~f~f''+ f'''~=~0 \f$

* The boundary conditions are given in the textbook. They depend on the reference length and velocity
which are different in the laminar and turbulent regime.

* The table below explains how to get turbulent profile from a laminar one. But it is not used in the code:

   \f$ \begin{array}{lll}
   {\rm variable-function | regime}  \qquad& {\rm   laminar} \qquad   & {\rm   turbulent} \vspace{0.3cm}  \\ \hline \hline
   \eta          &   \eta_{\ell}      &   \eta_t \vspace{0.3cm}\\
   f(\eta)       &  f_{\ell}          &   f_t    \vspace{0.3cm}\\
   U_{ref}       &   U_1              & \displaystyle \frac{U_1+U_2}{2}
  \end{array} \f$

 We get the same results by a change of variable and functions:

  \f$\eta_{\ell}~=~\alpha~\eta_t\f$ and  \f$ f_{\ell} = f_t~/~\alpha\f$ with
   \f$\displaystyle \alpha^2~=~\frac{2}{1+\lambda}\f$ and 
  \f$\displaystyle \lambda~=~\frac{U_2}{U_1}\f$

* There are two possibilities in this usual problem
    1. the x-axis correspond to the separation line (i.e v=0 => \f$f(0)= 0) \f$.

    2. on the x-axis, \f$ f(0) \neq 0 \f$, and the separation line must be calculated. It is the case here.

* Methodology:
    - Newton method : for case one, just with two parameters  \f$f'(0)=u_0\f$ and \f$f''(0)=g_0\f$  to define the shear layer boundaries,
    for the case two, there an additionnal Newton method with the parameter \f$f(0)=f_0\f$
    - first order Euler scheme (order 1) or Runge-Kutta order 4 scheme for the integration

* Heat transfer problem is also solved, but not validated for turbulent flow

* Main results can be found in a table in the textbook.

## Compilation

*  Two codes are provided, the first one **mixing_layer_o1** is very fast, it is the old version, with only an Euler scheme.
In the second version, **mixing_layer**, the Runge-Kutta scheme is implemented, and more parameters are included in the input file
**mixing_layer.in**.

* A makefile has been written. The compilation is carried out in a terminal :

\code{.sh}
make
\endcode

* `make clean` is used to clean the executable.


* With the make file, the  both executable codes `mixing_layer` and  `mixing_layer_o1` are  moved  in `~/bin/` directory. It can be modified by the user.


## Run
In the `run` directory, for instance :

\code{.sh}
mixing_layer 
\endcode

or 

\code{.sh}
./x1.run.sh 6
\endcode

* if 1 instead of 6 :  `mixing_layer_o1`
* A  shell script  also located in this directory can be used, to run and to clean the run directory:

\code{.sh}
x1.run.sh 0
\endcode

Other options are possible with this shell script.

## Options and parameters

Many parameters have to be given.

The input file, named **mixing_layer.in** is :


\code{.sh}
# input file for mixing_layer 
1                            # turbulence, O : Uref=U1 (laminar),  1 : Uref= (U1+U2)/2 (turbulent)
1.0                          # u1 
0.25                         # u2 or lambda =u2/u1 if u1 = 1
0.0004                       # deta, step in the eta direction
20.                          # eta_max,  il >= eta_max/deta
50000                        # il, half mesh size
1                            # 0: Euler O1, 1: RK4
\endcode


* for the code `mixing_layer_o1`, just the first four parameters are used, the other are in the source code.
* for the code `mixing_layer`, all the main parameters are read in the input file.

* There a relationship between **eta_max**, **deta** and **il** as written in the input file line 6.


* The first way is to run the code with the `read_option=1` as shown below 
and default parameters written in the main module are used.

\code{f90}
!> @details
!! program to test the module
program main_mixing_layer

use mod_mixing_layer
call solve_mixing_layer

end program main_mixing_layer
\endcode

## Database, output files

* There are several ouputts file :
    - **res.out** : contains the parameters and main ouputs (table)
    - **profile.m** : for Matlab used
    - **outputs.out** : profiles
    - **positive.dat** and **negative.dat** : upper and lower part of the profiles
    - **convergence.out** and **F0.out** : indications about convergence of the  Newton algorithm

## Caution

To get very high accuracy a grid step much be lower  than 0.0004