#!/bin/python
"""
   résolution de problème d'aérodynamique
  @author{C. Airiau}
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 

from pathlib import Path
import sys
print("plateforme                 : ",sys.platform)
print("home directory             : ",Path.home())
if  sys.platform=="linux":
    home = str(Path.home())
    FortranPath = home+"/ownCloud/Logiciels/EC/src_F90/"
    FortranRun  = home+"/ownCloud/Logiciels/EC/run/"
    PyRun       = home+"/ownCloud/Livre_Aerodynamique/Livre_Exercices_Aero/Codes/Python/Run/"

    FortranPath = home+"/Documents/Depot_FundAeroSuite/Fortran/EC/src_F90/"
    FortranRun  = home+"/Documents/Depot_FundAeroSuite/Fortran/EC/run/"
    PyRun       = home+"/Documents/Depot_FundAeroSuite/Python/Run/"

elif sys.platform=="win32":
    home="D:"
    FortranPath = home+"/ownCloud/Logiciels/EC/src_F90/"
    FortranRun  = home+"/ownCloud/Logiciels/EC/run/"
    #PyRun       = home+"/ownCloud/Livre_Aerodynamique/Livre_Exercices_Aero/Codes/Python/Run/"
    PyRun        = "F:/Christophe/Codes/Python/Run"
else:
    raise ValueError("pas de chemins implémentés dans le fichier CheminsParDefaut.py")

# 
