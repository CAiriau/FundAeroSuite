# FundAeroSuite

## Solutions of Fundamental Aerodynamics Exercices and Problems

The text of the problems and exercices are found into the textbook

---



**Exercices et Problèmes d'Aérodynamique Fondamentale**

Christophe Airiau, André Giovannini, Pierre Brancher

Cépadues Editions, France

2019,  I.S.B.N. : 9782364937253, 383 pages.

[Link](https://www.cepadues.com/livres/exercices-problemes-aerodynamique-fondamentale-accompagnes-des-codes-solutions-python-fortran-9782364937253.html)

---
The suite  contains some python, fortran, matlab and Maple sources.

### Warnings

* The master branch of the github site is build with  code written in French language (2019).

* The master_EN is build with code in English language (2021)

### Content

The numerical suite follows the chapters of the book  written in French:

1. *General equations of conservation and evolution* - Equations générales de conservation et d'évolution
2. *Mathematical models for aerodynamics* - Modèles mathématiques pour l'aérodynamique  
3. *Complex potential theory : application to airfoils* - Théorie des potentiels complexes : application aux profils d'aile
4. *Linearized theory of thin airfoils* - Théorie linéarisée des profils minces
5. *Lifting line Lanchester-Prandtl theory* - Théorie de la ligne portante de Lanchester-Prandtl
6. *Lifting surface theory and slender bodies* - Théorie de la surface portante et des corps élancés
7. *Panel method* - Aspects numériques : méthode des singularités
8. *Compressible subsonic and transonic flow* - Ecoulements compressible subsonique et transsonique
9. *Linearized supersonic flows* - Ecoulements supersoniques linéarisés
10. *One-dimensional compressible flows* - Ecoulements compressibles monodimensionnels
11. *Two-dimensional supersonic flows* - Ecoulements supersoniques bidimensionnels
12. *Method of characteristics for steady and unsteady flows* - Méthodes des caractéristiques en régimes stationnaire et instationnaire
13. *Slender bodies and wings in supersonic flows* - Ecoulements supersoniques : corps élancés et ailes
14. *Hypersonic flows* - Ecoulements hypersoniques
15. *Viscous flows and boundary layer* - Effets visqueux et couche limite

#### Codes

* The book contains the solutions of 120 problems and exercises and the numerical suite provides solutions of approximately 100 of them.
* Approximately 80 problems or exercises solutions are given in **Python**
* Approximately 22 problems or exercises solutions are given in **Fortran**
* Some solutions are also given in Matlab
* Some solutions are given in Maple (formal language)
* Any part of codes written in Matlab or Maple can also be found (in some way) in Fortran or Python.
* The codes are written in some educational way in order to be easily understood even if the reader is not a programmer.
* Some additional tools/codes/programs are also provided, as for instance to calculate :

  * NACA 4 airfoils geometry
  * Standard atmosphere quantities
  * ...

* Many comments are added inside the code, unfortunaly many are written in French. Some effort will be done to progressively translate them in English.

### Installation

* The Fortran codes can be compiled on any system (Windows, MacOs, Linus) as soon as a Fortran Compiler is installed
* The Python codes can be easily run with a minimun of knowledges :
  * a GUI can be used (validated only on Linux and Windows systems)
  * a main script can be used where  the number of the exercise and of the chapter have to be written
  * Python 3 can be installed for instance with anaconda, some interfaces as spyder3 can be used as well
* To use Maple files the Software software Maple has to be installed (it is not free).
 Some ascii files and pdf files are also given (*i n progress*)
* To use Matlab files, the Software Matlab has to be installed (it is not free). The Matlab script are in ascii files. Possibly the freeware Octave can be used.
* In each directory a README file can be found

### Documentation

* The book of exercices and problems provides the text and notations and some results.
* Documentation of Fortran codes can be generated with *Oxygen*
* Documentation of Python codes can be generated with *DOxygen*

## history of the sources

The initial fortran codes have been written in 2011 for a course in compressible aerodynamics at Toulouse III university. It has been improved for the content of the first book "Aerodynamique Fondamentale" published in 2016 (Cépadues ed. Toulouse)
and revised for the new book ("Exercices et problèmes d'aérodynamique fondamentale") which has been  pubished in July 2019.

The python sources have been written for the new book in 2017 and 2018.

The codes have been improved and translated in English in april 2020, during the covid19 French blackout.

 -------------------------------

Christophe Airiau, March 15, 2021

