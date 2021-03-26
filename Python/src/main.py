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

1. Python solutions:
  
* set language
    .. code-block:: python

        src = "python"

* define the list of exercises to solve :
    - a **single exercise** :  solution( ( (chap_number, (exercice_number, ) ), ) )
        .. code-block:: python

            solution(( (10,(5,)),))   # exercise 10,  chapter 5

    - **several exercises, several chapters**, write a list as
        .. code-block:: python

            solution(( (8, (6, 7)) , (9, (1, 2)) ))

2. Fortran solutions:

* set language
    .. code-block:: python
    
        src = "fortran"
  
* define the task with the task index
    0. "liste"  : exercise list with fortran coding
    1. "compilation" : to compile the source code
    2. "run" : to run the source code
    3. "compil+run" : make the both
    
    .. code-block:: python
    
        task = (3,)
    
* select the exercise reference into a list, similarly to python:
    .. code-block:: python
    
        choice = ( (10, (1, 2)), (11,(2,)), )
          
* Fortran tasks are not currently implemented with Windows. But it can be performed by-hand
          in the Fortran source directory.
* the Fortran source and run directory path must be implemented into the file :py:mod:`default_path`
* the list of exercises solved in Fortran can be found in :py:func:`Book.Fortran.fortran.exercise_list_fortran`

"""
import os
import numpy as np
from Book.Corrections import solution
from Book.Fortran.fortran import run_fortran

__licence__ = "AGPL-3.0"
__docformat__ = 'reStructuredText'


# ************************************
# main program
# ************************************

def main(src="python"):
    """
    main of the program
    """
    # choose language  solutions here
    src = "python"  # "fortran" or "python"
    # src = "fortran"
    # Fortran solutions
    task_list = ("listing", "compilation", "run", "compil+run")
    task = task_list[3]

    if src == "python":
        solution(((4, (7, )), ))
    elif src == "fortran":
        # example:
        # choice = ( (10, (1, 2)), (11,(2,)), )
        choice = ((10, (5, )), )
        run_fortran(task, choice)


if __name__ == "__main__":
    main(src="python")
