# Tutorial CAT - exercise list

# Compile and run 

## Compile and run from Python GUI

It is possible to compile and run the Fortran code with the Python GUI found in the Python directory.

Please see the Python documentation if it is your choice.

If you know Fortran and how to use Fortran it is better to follow the steps in the following.

## Compile in the **src** directory

* Go to that directory and write in a terminal :

\code{.sh}
make clean; make
\endcode

* On the main **CAT** directory there is a script shell to do the job as well :

\code{.sh}
./x3.compile_CAT.sh
\endcode

It makes the same thing.

* `make clean` is only necessary if the source of the module is modified. It due to the manner the Makefile is written.

## Run the CAT tool

There are several ways to run the **CAT tool** and are now described.

### 1 - Run and choose the options

\code{.sh} CAT \endcode

A menu let choose the options

**27** is the option number to run one or several exercises.

### 2 - Run when the option is written into a input file *CAT.in*

\code{.sh} CAT < CAT.in \endcode

The *CAT.in* file only contains the option number.

### 3 - Run with an input file and outputs are written into a file

\code{.sh} CAT < CAT.in > CAT.out \endcode

This latter solution can be simplified by using the a shell script:

\code{.sh} ./x1.run_CAT.sh  \endcode

### 4 - Run with options in command line

The following command lines are possible (after # it is a comment, not to write) :

\code{.sh}
CAT h                          # to get the help
CAT opt= n                     # n the option number
CAT opt= 14 M= 0.5 A/Ac= 3.5   # for subsonic solution
CAT opt= 14 M= 1.5 A/Ac= 3.5   # for supersonic solution
CAT opt= 17 M= 0.5 F= 0.4      # for subsonic solution
CAT opt= 17 M= 1.5 F= 0.4      # for supersonic solution
CAT opt= 10 M= 1               # Shock interaction, M is the testcase here
CAT opt= 1 M= 2.0 angle= 10.0  # angle is in degrees
\endcode

etc ...

Please respect the white space after CAT and the = signs.

# Choose one or several exercises / list

The choice is configured into the file  */run/exercise_choice* 
and an example can be found in *doc/example/exercise_choice.in*

\code{.sh}

# this file contains the choice of the exercises to solve
#  ligne 3 indicates the number of exercise to solve
#  In the next line indicate the chapter and the index of the exercise. Please follow the example below.
1                           # number of exercise
12 3                        # example of chapter and exercise
11 61
11 10

# end

# list of the exercises available :

          #**************
          # CHAPTER 10
          #**************
          #
          #               
10 1      #  Exercice_10_1    ! Normal shock wave / Calcul d'un choc droit
10 2      #  Exercice_10_2    ! Normal shock wave calculus / Calcul pratique d'un choc droit
10 3      #  Exercice_10_3    ! Subsonic Fanno's flow / Ecoulement de Fanno subsonique
10 4      #  Exercice_10_4    ! Supersonic Fanno's flow / Ecoulement de Fanno supersonique
10 5      #  Exercice_10_5    ! Rayleigh's flow / Ecoulement de Rayleigh
10 6      #  Exercice_10_6    ! Rayleigh's flow with combustion / Ecoulement de Rayleigh + combustion
10 7      #  Exercice_10_7    ! Pitot tube in subsonic regime / Pitot subsonique
10 8      #  Exercice_10_8    ! Pitot tube in supersonic regime / Pitot supersonique
10 9      #  Exercice_10_9    ! Regimes in nozzles  / Régimes dans une tuyere
10 10     #  Exercice_10_10   ! Laval Nozzle  / Tuyère de Laval
          #
          #**************
          # CHAPTER 11
          #**************
          #
11 1      #  Exercice_11_1    ! Flat plate in incidence / plaque en incidence
11 2      #  Exercice_11_2    ! Diamond airfoil / profil losangique
11 3      #  Exercice_11_3    ! Shock-expansion method / methode choc detente
11 4      #  Exercice_11_4    ! Shocks in a channel / Choc dans un canal
11 5      #  Exercice_11_5    ! Two obliques shock wave interactions / Interaction de 2 chocs obliques
11 6      #  Exercice_11_6    ! Backward facing step  (simple program) / Marche descendante (programmation basique)
11 61     #  Exercice_11_61   ! Backward facing step (2nd program) / Marche descendante (2nd programmation)
11 7      #  Exercice_11_7    ! Shock - isobar line interactions of a fluid at rest  / interaction choc - ligne isobare avec fluide au repos
11 8      #  Exercice_11_8    ! Expansion - isobar line interactions of a fluid at rest / interaction détente - ligne isobare avec fluide au repos
11 9      #  Exercice_11_9    ! Expansion - wall interaction / interaction détente - paroi
11 10     #  Exercice_11_10   ! Two upstream flow  separated by a isobar line  - oblique shock interaction / deux écoulements amont + Ligne de glissement + choc oblique
11 11     #  Exercice_11_11   ! Interaction of shock waves in a nozzle outlet / tuyère : interaction d'ondes de chocs en sortie
11 12     #  Exercice_11_12   ! Interaction of a centered expansion in a nozzle outlet / tuyère : interaction d'un faisceau de détente en sortie
          #
          #**************
          # CHAPTER 12
          #**************
          #
12 3      #  Exercice_12_3    ! Shock tube application / application du tube à choc.

the third  column is the Python reference and the last column is the name of the exercise.
\endcode

The outputs are in a file named with the exercise reference for instance
*Solution_10-1.out*

# Complex configurations

It it possible to run the code on a complex configurations, for instance
a flat plate with incidence, diamond airfoil in incidence.

To do this, you have to write the configuration the file `/run/case.in`.
An example is provided in *doc/example/case.in*:

\code{.sh}
#***********************************
# case description
#***********************************
3                       # number of zones (without the upstream zone indexed 0)
0.d0                    # incidence
# Zone 0
0                       # zone index
'Uniforme'              # flow type
2.2d0                   # Mach number
230.d0                  # Ti0
10000.d0                # P0
0.d0                    # flow deviation angle in °
0                       # upstream zone
# Zone 1
1                       # zone index
'Choc'                  # flow type
7.d0                    # flow deviation angle in °
0                       # upstream zone
# Zone 2
2                       # zone index
'Choc'                  # flow type
12.d0                   # flow deviation angle in °
1                       # upstream zone
# Zone 3
3                       # index de zone
'Choc'                  # type d'écoulement
0.d0                    # flow deviation angle in °
2                       # upstream zone

# end of the data

flow type can be :
'rien'          : nothing happen
'Uniforme'      : uniform flow
'Choc'          : shock wave
'Detente'       : expansion

Please keep the right format
\endcode

To run this case you have to select number  **99**

 \warning In case of airfoil in incidence, it is better to directly modify the flow
 deviation angle in the *case.in* file and let the "incidence" parameter to 0 (line 5).

# List of CAT options

 n°   |   Task                                              || n°   |   Task                                              |
 -----|-----------------------------------------------------|| -----|-----------------------------------------------------|
  0   | **Normal shock wave**                               || 12   |   Tables for the book                               |
  1   | **Oblique shock wave**                              ||  9   |   Shock table                                       |
  2   | **Isentropic compression and expansion**            ||  7   |   omega(M) table                                    |
  3   | **Prandt-Meyer function omega(Mach)**               || 23   |   personal shock tables                             |
  4   | **Mach(omega) : inverse of Prandt-Meyer function**  ||  6   |   Prandt-Meyer curve                                |
  5   | **Isentropic evolution (physical quantities)**      ||  8   |   Epicycloid curve                                  |
  .   | .                                                   ||  .   | .                                                   |
 10   | Shock interaction, (p,theta) plane                  ||      |                                                     |
 13   | **A / A_critical for a given Mach**                 || 15   |   A / A critique = f(Mach)  curve                   |
 14   | **Mach for a given A / A_critical**                 ||      |                                                     |
 16   | Fanno problem   F=f(Mach)                           || 18   |   Fanno curve                                       |
 17   | Inverse Fanno problem M=Fanno(F)                    ||      |                                                     |
 19   | Rayleigh problem                                    || 21   |   Rayleigh function table f(Mach)                   |
 20   | Rayleigh curve                                      || 22   |   Inverse Rayleigh function table                   |
  .   | .                                                   ||  .   | .                                                   |
 27   | **Exercise of the book**                            || 99   |  **case in the input file**                         |
 24   | **Standard atmosphere**                             ||      |                                                     |
  .   | .                                                   ||  .   | .                                                   |
 90   | May 2014 exam (example of application)              || 91   |   June 2014 exam (example of application)           |

The number indicates the value of the option to enter.

*Some options are hidden.*
