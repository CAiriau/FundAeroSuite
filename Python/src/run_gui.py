"""
Contents, references
====================

Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Main script to run in Python GUI tkinter
    ..

Bugs
====

There are many bugs ...

- when an exercise is solved, we sometimes have to push twice the "Exit" button
- with Windows "Quit"  in the general menu does not work 
- "Exit" button is not displayed on too small screen
- with windows you must not run the GUI with spyder but you have to use the windows shell and write :  python3 run_gui.py

"""

# import os
# import platform
# import numpy             as np
# import matplotlib.pyplot as plt

# from tkinter import *
# import tkinter
# from tkinter.messagebox import *
# from tkinter.ttk import *
# from tkinter import filedialog


import Tools.misc

# from CompressibleFlow.fonctions import *
# from Book.Chapter2             import *
# from Book.Chapter3             import *
# from Book.Chapter4             import *
# from Book.Chapter5             import *
# from Book.Chapter6             import *
# from Book.Chapter7             import *
# from Book.Chapter7c            import *
# from Book.Chapter8             import *
# from Book.Chapter9             import *
# from Book.Chapter10            import *
# from Book.Chapter11            import *
# from Book.Chapter12            import *
# from Book.Chapter13            import *
# from Book.Chapter14            import *
# from Book.Corrections          import *
# from Book.misc                 import *
# from conical_shock.conical_shock  import *
# from Book.panneaux_intro       import *
# from Book.Fortran.fortran      import *
# import default_path 
from GUI.guiAero import GuiWin, set_params_gui


def main():
    par = set_params_gui()
    par["ShowTree"] = True
    par["ShowTitre"] = True
    par["redirection"] = False
    par["ShowCheck"] = False
    par["ShowText"] = False
    par["language"] = "en"  # "en" or "fr"

    par["default"] = {"tags": "chap", "values": ["Exercice5_12", "py"]}  # example of default case
    s = GuiWin(par)


if __name__ == "__main__":
    main()
