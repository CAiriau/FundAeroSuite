"""
Fundamental Aerodynamic Suite, FundAeroSuite
  ..

.. author: C. Airiau
.. date: 2020/04/10
.. licence: AGPL-3.0

Main script to run in Python terminal
"""


import os
import numpy             as np
import matplotlib.pyplot as plt
from Tools.misc                 import *
from CompressibleFlow.fonctions import *
from Book.Chapter2             import *
from Book.Chapter3             import *
from Book.Chapter4             import *
from Book.Chapter5             import *
from Book.Chapter6             import *
from Book.Chapter7             import *
from Book.Chapter7c            import *
from Book.Chapter8             import *
from Book.Chapter9             import *
from Book.Chapter10            import *
from Book.Chapter11            import *
from Book.Chapter12            import *
from Book.Chapter13            import *
from Book.Chapter14            import *
from Book.Corrections          import *
from Book.misc                 import *
from conical_shock.conical_shock  import *
from Book.panneaux_intro       import *
from Book.Fortran.fortran      import *

"""
# FONCTIONNEMENT :
 1) Codes en python :
     Mettre src="python"
     Définir la liste des exercices à résoudre par la liste solution[]
       suivre les exemples indiquées

 2) Codes en fortran :
     Mettre src="fortran"
     Définir la tâche par task=" "
     Définir la liste des exercices par choix=[]
     Attention faire d'abord task="liste" pour voir comment définir "choix"
     
     l'utilisation du fortran sous windows est à implémenter.
     
     Penser à mettre les bons chemins dans le fichier "default_path.py"
"""

#************************************
# main program
#************************************

# PYTHON
# première liste : liste des chapitres
# seconde liste : liste des exercices.
#

# utilisation :
# exemple 1 chapitre, 1 ou plusieurs exercices :
#           correction([8],[5,6])  ou
#           solution([[8,[5,6]]])
# exemple plusieurs chapitres, plusieurs exercices :
# solution([ [8,[6,7]] , [9,[1,2]] ])

src="python"  # "fortran" ou "python"

# FORTRAN
task_list = ["liste", "compilation", "run", "compil+run"]
task=task_list[3]
# ------------------------------------------------

if src=="python":
  solution([[10,[5]]])
elif src=="fortran":
    # example:    choix=[[10,[1,2]],[11,[2]]] 
    choice=[[10,[1,2]],[11,[2]]]  # liste des numéros d'exercices
    # la liste des exercices fortran se trouvent dans Book/Fortran/fortran.py
    run_fortran(task,choice)

