# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Default paths used for the gui and the fortran compiler for the gui for the Fundamental Aerodynamics Suite
    ..

"""
__licence__ = "AGPL-3.0"
__docformat__ = 'reStructuredText'

from pathlib import Path
import sys

print("plateform                  : ", sys.platform)
print("home directory             : ", Path.home())
if sys.platform == "linux":
    home = str(Path.home())
    FortranPath = home + "/ownCloud/Logiciels/CAT/src"
    FortranRun = home + "/ownCloud/Logiciels/CAT/run"
    PyRun = home + "/ownCloud/Livre_Aerodynamique/Livre_Exercices_Aero/Codes/Python/Run"

    FortranPath = home + "/Projects/FundAeroSuite/Fortran/CAT/src"
    FortranRun = home + "/Projects/FundAeroSuite/Fortran/CAT/run"
    PyRun = home + "/Projects/FundAeroSuite/Python/Run"

elif sys.platform == "win32":
    home = "D:"
    FortranPath = home + "/ownCloud/Logiciels/EC/src"
    FortranRun = home + "/ownCloud/Logiciels/EC/run"
    # PyRun       = home+"/ownCloud/Livre_Aerodynamique/Livre_Exercices_Aero/Codes/Python/Run"
    PyRun = "F:/Christophe/Codes/Python/Run"
else:
    raise ValueError("No path are implemented in the file 'default_path.py' ")
