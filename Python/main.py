#!/usr/bin/env python
# -*- coding: utf-8 -*-


#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
#   Licence : AGPL-3.0

import os
import numpy             as np
import matplotlib.pyplot as plt
from Tools.misc                 import *
from CompressibleFlow.fonctions import *
from Livre.Chapter2             import *
from Livre.Chapter3             import *
from Livre.Chapter4             import *
from Livre.Chapter5             import *
from Livre.Chapter6             import *
from Livre.Chapter7             import *
from Livre.Chapter7c            import *
from Livre.Chapter8             import *
from Livre.Chapter9             import *
from Livre.Chapter10            import *
from Livre.Chapter11            import *
from Livre.Chapter12            import *
from Livre.Chapter13            import *
from Livre.Chapter14            import *
from Livre.Corrections          import *
from Livre.misc                 import *
from Choc_Conique.choc_conique  import *
from Livre.panneaux_intro       import *
from Livre.Fortran.fortran      import *

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
     
     Penser à mettre les bons chemins dans le fichier "CheminsParDefaut.py"
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
  solution([[2,[2]]])
elif src=="fortran":
    # example:    choix=[[10,[1,2]],[11,[2]]] 
    choice=[[10,[1,2]],[11,[2]]]  # liste des numéros d'exercices
    # la liste des exercices fortran se trouvent dans Livre/Fortran/fortran.py
    run_fortran(task,choice)

