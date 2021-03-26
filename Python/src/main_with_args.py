"""
Contents, references
====================

Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Main script to run in Python terminal
  ..

help
=====

1. Description of the task inside the script

* for **python** solution : solution( ( (chap_number, (exercice_number, ) ), ) )
    .. code-block:: python

        solution(( (2,(1,)),))   # exercise 1,  chapter 2

* for **fortran** solution :  run_fortran("compil+run", [[chap_number, [exercice_number]]])
    .. code-block:: python

        run_fortran("compil+run",  ( (10, (1,) ), )  # exercise 1,  chapter 10

2. Command line in a terminal:

* for **python** solution, to solve exercise 1,  chapter 2 :
    .. code-block:: bash

        python3 main_with_args -c 2 -e 1 -l p

    arguments are optional (it is the default configuration)

* for **fortran** solution, to solve exercise 1,  chapter 10 :
    .. code-block:: bash

        python3 main_with_args -c 10 -e 1 -l f

    all the arguments must be written!

3. Caution :

* The user has to enter the right value of the exercise depending of the language

* No verification is made in case of bad input parameters.

* for Fortran, the user has to verify once the default fortran path in "default_path.py"
"""

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from textwrap import dedent
from Book.Corrections import solution
from Book.Fortran.fortran import run_fortran

def main():
    """
    main program to solve a single exercice with Python
    """
    parser = ArgumentParser( \
        formatter_class = RawDescriptionHelpFormatter, \
        description = dedent('''\
            Script to solve a single exercise
            If no argument is provided, a demo is displayed.
            '''), \
        epilog = dedent('''\
            Examples:
                Get help
                    python {0} -h
                solve an exercise  (chapter: ch, exercise number : nu)
                    python {0} -c ch -e nu
                 solve an exercise Exercice2_1 (chapter 2, number 1)
                    python {0} -c 2 -e 1
            '''.format(os.path.basename(__file__))))
    parser.add_argument('-c','--chapter', type = int, default = 2, \
                        help = 'chapter index. Default is 2. Example: 2, 3, ..., 15')
    parser.add_argument('-e','--exercise', type = int, default = 1, \
                        help = 'exercise index. Default is 1. Example: 1, 2 , 3')
    parser.add_argument('-l', '--lang', type=str, default='p', \
                        help='language. Default is p, Python : p, Fortran : f, Example: p or f')
    args = parser.parse_args()

    if args.lang == "p":
        solution( ( (args.chapter, (args.exercise,) ),) )
    elif args.lang == "f":
        run_fortran("compil+run", ( (args.chapter, (args.exercise, ) ), ) )


if __name__ == "__main__":
    main()