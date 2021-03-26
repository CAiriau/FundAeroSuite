# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Fortran package:
    ..
    
    * List of exercices with a fortran correction
    * function to run or compile Fortran source
    
"""

import glob
import os
import sys
import numpy             as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from default_path import home, PyRun, FortranRun, FortranPath


def exercise_list_fortran(display=False, lang="en"):
    exercise_list = []
    if lang == "fr":
        exercise_list.append({'ch': 10, 'num': 1, 'exo': 'Exercice_10_1', 'com': "Calcul d'un choc droit"})
        exercise_list.append({'ch': 10, 'num': 2, 'exo': 'Exercice_10_2', 'com': "Calcul pratique d'un choc droit"})
        exercise_list.append({'ch': 10, 'num': 30, 'exo': 'Exercice_10_10', 'com': "Tuyère de Laval"})
        exercise_list.append({'ch': 10, 'num': 3, 'exo': 'Exercice_10_3', 'com': "Ecoulement de Fanno subsonique"})
        exercise_list.append({'ch': 10, 'num': 4, 'exo': 'Exercice_10_4', 'com': "Ecoulement de Fanno supersonique"})
        exercise_list.append({'ch': 10, 'num': 5, 'exo': 'Exercice_10_5', 'com': "Ecoulement de Rayleigh"})
        exercise_list.append({'ch': 10, 'num': 6, 'exo': 'Exercice_10_6', 'com': "Ecoulement de Rayleigh + combustion"})
        exercise_list.append({'ch': 10, 'num': 7, 'exo': 'Exercice_10_7', 'com': "Pitot subsonique"})
        exercise_list.append({'ch': 10, 'num': 8, 'exo': 'Exercice_10_8', 'com': "Pitot supersonique"})
        exercise_list.append({'ch': 10, 'num': 9, 'exo': 'Exercice_10_9', 'com': "Régimes dans une tuyere"})

        exercise_list.append({'ch': 11, 'num': 1, 'exo': 'Exercice_11_1', 'com': "plaque en incidence"})
        exercise_list.append({'ch': 11, 'num': 2, 'exo': 'Exercice_11_2', 'com': "profil losangique"})
        exercise_list.append({'ch': 11, 'num': 3, 'exo': 'Exercice_11_3', 'com': "methode choc detente"})
        exercise_list.append({'ch': 11, 'num': 4, 'exo': 'Exercice_11_4', 'com': "Choc dans un canal"})
        exercise_list.append({'ch': 11, 'num': 5, 'exo': 'Exercice_11_5', 'com': "Interaction de 2 chocs obliques"})
        exercise_list.append(
            {'ch': 11, 'num': 6, 'exo': 'Exercice_11_6', 'com': "Marche descendante (programmation basique)"})
        exercise_list.append(
            {'ch': 11, 'num': 61, 'exo': 'Exercice_11_6bis', 'com': "Marche descendante (2nd programmation)"})
        exercise_list.append({'ch': 11, 'num': 7, 'exo': 'Exercice_11_7',
                          'com': "interaction choc - ligne isobare avec fluide au repos"})
        exercise_list.append({'ch': 11, 'num': 8, 'exo': 'Exercice_11_8',
                          'com': "interaction détente - ligne isobare avec fluide au repos"})
        exercise_list.append({'ch': 11, 'num': 9, 'exo': 'Exercice_11_9', 'com': "interaction détente - paroi"})
        exercise_list.append({'ch': 11, 'num': 10, 'exo': 'Exercice_11_10',
                          'com': "deux écoulements amont + Ligne de glissement + choc oblique"})
        exercise_list.append(
            {'ch': 11, 'num': 11, 'exo': 'Exercice_11_11', 'com': "tuyère : interaction d'ondes de chocs en sortie"})
        exercise_list.append({'ch': 11, 'num': 12, 'exo': 'Exercice_11_12',
                          'com': "tuyère : interaction d'un faisceau de détente en sortie"})

        exercise_list.append({'ch': 12, 'num': 3, 'exo': 'Exercice_12_3', 'com': " application du tube à choc."})

    else:

        exercise_list.append({'ch': 10, 'num': 1, 'exo': 'Exercice_10_1', 'com': "Normal shock wave"})
        exercise_list.append({'ch': 10, 'num': 2, 'exo': 'Exercice_10_2', 'com': "Normal shock wave calculus"})
        exercise_list.append({'ch': 10, 'num': 30, 'exo': 'Exercice_10_10', 'com': "Laval's nozzle"})
        exercise_list.append({'ch': 10, 'num': 3, 'exo': 'Exercice_10_3', 'com': "Subsonic Fanno's flow"})
        exercise_list.append({'ch': 10, 'num': 4, 'exo': 'Exercice_10_4', 'com': "Supersonic Fanno's flow"})
        exercise_list.append({'ch': 10, 'num': 5, 'exo': 'Exercice_10_5', 'com': "Rayleigh's flow"})
        exercise_list.append({'ch': 10, 'num': 6, 'exo': 'Exercice_10_6', 'com': "Rayleigh's flow with combustion"})
        exercise_list.append({'ch': 10, 'num': 7, 'exo': 'Exercice_10_7', 'com': "Pitot tube in subsonic regime"})
        exercise_list.append({'ch': 10, 'num': 8, 'exo': 'Exercice_10_8', 'com': "Pitot tube in supersonic regime"})
        exercise_list.append({'ch': 10, 'num': 9, 'exo': 'Exercice_10_9', 'com': "Regimes in nozzles"})

        exercise_list.append({'ch': 11, 'num': 1, 'exo': 'Exercice_11_1', 'com': "Flat plate in incidence"})
        exercise_list.append({'ch': 11, 'num': 2, 'exo': 'Exercice_11_2', 'com': "Diamond airfoil"})
        exercise_list.append({'ch': 11, 'num': 3, 'exo': 'Exercice_11_3', 'com': "Shock-expansion method"})
        exercise_list.append({'ch': 11, 'num': 4, 'exo': 'Exercice_11_4', 'com': "Shocks in a channel"})
        exercise_list.append({'ch': 11, 'num': 5, 'exo': 'Exercice_11_5', 'com': "Two obliques shock wave interactions"})
        exercise_list.append({'ch': 11, 'num': 6, 'exo': 'Exercice_11_6', 'com': "Backward facing step  (simple program)"})
        exercise_list.append({'ch': 11, 'num': 61, 'exo': 'Exercice_11_6bis', 'com': "Backward facing step (2nd program)"})
        exercise_list.append(
            {'ch': 11, 'num': 7, 'exo': 'Exercice_11_7', 'com': "Shock - isobar line interactions of a fluid at rest"})
        exercise_list.append({'ch': 11, 'num': 8, 'exo': 'Exercice_11_8',
                          'com': "Expansion - isobar line interactions of a fluid at rest"})
        exercise_list.append({'ch': 11, 'num': 9, 'exo': 'Exercice_11_9', 'com': "Expansion - wall interaction"})
        exercise_list.append({'ch': 11, 'num': 10, 'exo': 'Exercice_11_10',
                          'com': "Two upstream flow  separated by a isobar line  - oblique shock interaction"})
        exercise_list.append(
            {'ch': 11, 'num': 11, 'exo': 'Exercice_11_11', 'com': "Interaction of shock waves in a nozzle outlet"})
        exercise_list.append({'ch': 11, 'num': 12, 'exo': 'Exercice_11_12',
                          'com': "Interaction of a centered expansion in a nozzle outlet"})

        exercise_list.append({'ch': 12, 'num': 3, 'exo': 'Exercice_12_3', 'com': "Shock tube application"})

    if display:
        print("=" * 90)
        print("k  chap.    num.     exercise   \t comments")
        print("=" * 90)

        k = 0
        for dico in exercise_list:
            # for key,values in dico.items():
            print("%2i   %2i     %2i %20s  %s" % (k, dico['ch'], dico['num'], dico['exo'], dico['com']))
            k = k + 1

        print("to select an exercise enter the index given by the first column")
        print("in the main.py file")
    return exercise_list


def choice_exos(path, choice, liste):
    """
    Generate the Fortran input file
    """
    maintenant = datetime.now()
    fichier = open(path + "/exercise_choice.in", "w")
    # print('file to write in : ', fichier)
    fichier.write('# file generated with Python\n')
    fichier.write('# ' + str(maintenant) + '\n')
    fichier.write('# number of exercises, exercise list \n')
    N = 0
    for L in choice:
        N = N + len(L[1])
    fichier.write("%2i \t\t # number of exercises \n" % (N))
    for L in choice:
        ch = L[0]
        for num in L[1]:
            for n in range(len(liste)):
                # print("is it ch ", liste[n]["ch"], "number ", liste[n]["num"] , " ?" )
                if liste[n]["ch"] == ch and liste[n]["num"] == num:
                    # print("exercise found : n° ", liste[n])
                    exo="%2i %2i \t\t # %s \n" % (ch, num, liste[n]["com"])
                    print('exercise found  : ', exo)
                    fichier.write(exo)
                    break

    fichier.write("# end \n")
    fichier.close()
    print("exercise_choice.in has been generated")


def run_EC(SourcePath, RunPath, option="compilation"):
    """
    compilation and run Fortran source
    
    Params:
        * SourcePath (string) : the path to the Fortran directory source
        * RunPath (string) : the path to the directory where the executable code is run
        * option (string) : "compilation" or "run"
    """
    if option == "compilation":
        os.chdir(SourcePath)
        # print("files = ",os.listdir("./"))
        print("Compile of EC tool ...")
        os.system("make clean; make")
    elif option == "run":
        os.chdir(RunPath)
        print("run the EC tool ...")
        os.system("more exercise_choice.in")
        os.system("./x1.run_EC.sh")
        print("Solution file found : ")
        os.system("ls Solution*")
        print("files : ", os.listdir("/"))


def test_local(task):
    """
    local test of  compilation and run
    """
    print("directory lists   : ", os.listdir("/"))
    print("Current directory : ", os.getcwd())
    if task == 0:
        # compilation
        os.chdir("../Source")
        os.system("make clean; make")
    elif task == 1:
        os.chdir("../run")
        os.system("./tp3.bin")  # can be changed to another binary file


def run_fortran(task="listing", choice=( (10, 1, ), ), SrcDir=FortranPath, RunDir=FortranRun):
    """
    Definition of the different task to perform
    """
    Src = False
    Run = False

    if os.path.isdir(SrcDir):
        print("path of Fortran sources : ", SrcDir)
        Src = True
    else:
        raise NameError("Source file ' %s ' not found" % (SrcDir))

    if os.path.isdir(RunDir):
        print("path of Fortran run directory : ", RunDir)
        Run = True
    else:
        raise NameError(" Fortran execution directory ' %s 'not found" % (RunDir))

    if Src:
        print("exercise list     : ", choice)

    if Run:
        choice_exos(RunDir, choice, exercise_list_fortran(False))  # create  exercise.in

    if task == "listing" and Src:
        liste = exercise_list_fortran(display=True)
    elif task == "compilation" and Src:
        run_EC(SrcDir, RunDir, "compilation")
    elif task == "run" and Run:
        run_EC(SrcDir, RunDir, "run")
    elif task == "compil+run" and Src and Run:
        run_EC(SrcDir, RunDir, "compilation")
        run_EC(SrcDir, RunDir, "run")
