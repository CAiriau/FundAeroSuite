#!/bin/python
"""
   résolution de problème d'aérodynamique
  @author{C. Airiau}
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
#   Licence : AGPL-3.0


import os
import platform
import numpy             as np
import matplotlib.pyplot as plt


#import tkinter
from tkinter.messagebox import * 
from tkinter.ttk import * 
from tkinter import filedialog  
from tkinter import * 

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
import CheminsParDefaut
from GUI.guiAero import *

par=set_params_gui()
par["ShowTree"]=True
par["ShowTitre"]=True
par["redirection"]=False
par["ShowCheck"]=False
par["ShowText"]=False
par["bugs"]="""
  Problèmes non résolus :
    - après la résolution d'un exercice, il faut appuyer 2 fois sur "Sortir"
        pour sortir
    - sous windows "Quitter" dans le menu "Général" ne fonctionne pas
    - le bouton "Sortir" n'apparaît pas sur les petits écrans
    - Sous windows il ne faut pas lancer le GUI sous spyder
    -              il faut utiliser windows shell et faire : python3 run_gui.py
        """
par["default"]  = {"tags":"chap","values":["Exercice5_12","py"]}  # exemple de cas par défaut
s= GuiWin(par)
