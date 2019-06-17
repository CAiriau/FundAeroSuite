# -*- coding: utf-8 -*-
"""
List of exercices with a fortran correction
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
 
import glob
import os
import numpy             as np
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from CheminsParDefaut import *



def liste_exo_fortran(display=False):

    liste_exo=[]
    liste_exo.append({'ch':10,'num':1,'exo':'Exercice_10_1','com':"Calcul d'un choc droit"})
    liste_exo.append({'ch':10,'num':2,'exo':'Exercice_10_2','com':"Calcul pratique d'un choc droit"})
    liste_exo.append({'ch':10,'num':30,'exo':'exo_Laval',    'com':"Tuyère de Laval"})
    liste_exo.append({'ch':10,'num':3,'exo':'Exercice_10_3','com':"Ecoulement de Fanno subsonique"})
    liste_exo.append({'ch':10,'num':4,'exo':'Exercice_10_4','com':"Ecoulement de Fanno supersonique"})
    liste_exo.append({'ch':10,'num':5,'exo':'Exercice_10_5','com':"Ecoulement de Rayleigh"})
    liste_exo.append({'ch':10,'num':6,'exo':'Exercice_10_6','com':"Ecoulement de Rayleigh + combustion"})
    liste_exo.append({'ch':10,'num':7,'exo':'Exercice_10_7','com':"Pitot subsonique"})
    liste_exo.append({'ch':10,'num':8,'exo':'Exercice_10_8','com':"Pitot supersonique"})
    liste_exo.append({'ch':10,'num':9,'exo':'Exercice_10_9','com':"Régimes dans une tuyere"})

    liste_exo.append({'ch':11,'num':1,'exo':'Exercice_11_1','com':"plaque en incidence"})
    liste_exo.append({'ch':11,'num':2,'exo':'Exercice_11_2','com':"profil losangique"})
    liste_exo.append({'ch':11,'num':3,'exo':'Exercice_11_3','com':"methode choc detente"})
    liste_exo.append({'ch':11,'num':4,'exo':'Exercice_11_4','com':"Choc dans un canal"})
    liste_exo.append({'ch':11,'num':5,'exo':'Exercice_11_5','com':"Interaction de 2 chocs obliques"})
    liste_exo.append({'ch':11,'num':6,'exo':'Exercice_11_6','com':"Marche descendante (programmation basique)"})
    liste_exo.append({'ch':11,'num':61,'exo':'Exercice_11_6bis','com':"Marche descendante (2nd programmation)"})
    liste_exo.append({'ch':11,'num':7,'exo':'Exercice_11_7','com':"interaction choc - ligne isobare avec fluide au repos"})
    liste_exo.append({'ch':11,'num':8,'exo':'Exercice_11_8','com':"interaction détente - ligne isobare avec fluide au repos"})
    liste_exo.append({'ch':11,'num':9,'exo':'Exercice_11_9','com':"interaction détente - paroi"})
    liste_exo.append({'ch':11,'num':10,'exo':'Exercice_11_10','com':"deux écoulements amont + Ligne de glissement + choc oblique"})
    liste_exo.append({'ch':11,'num':11,'exo':'Exercice_11_11','com':"tuyère : interaction d'ondes de chocs en sortie"})
    liste_exo.append({'ch':11,'num':12,'exo':'Exercice_11_12','com':"tuyère : interaction d'un faisceau de détente en sortie"})

    liste_exo.append({'ch':12,'num':3,'exo':'Exercice_12_3','com':" application du tube à choc."})
    if display:
        print("="*90)
        print("k  chap.    num.     exercice   \t commentaires")
        print("="*90)

        k=0
        for dico in liste_exo:
            #for key,values in dico.items():
            print("%2i   %2i     %2i %20s  %s"%(k,dico['ch'],dico['num'],dico['exo'],dico['com']))
            k=k+1

        print("Pour choisir le ou les exercices, il faut uniquement entrer l'index de la première colonne")
        print("dans le fichier main.py")
    return liste_exo

def choix_exos(path,choix,liste):
    """
    Génération du fichier input
    """
    maintenant=datetime.now()
    fichier=open(path+"choix_exercices.in","w")
    fichier.write('# fichier généré automatiquement en python\n')
    fichier.write('# '+ str(maintenant) +'\n')
    fichier.write('# number of exercices, exercice list \n')
    N=0
    for L in choix:
        N+=len(L[1])
    fichier.write("%2i \t\t # number of exercices \n"%(N))    

    for L in choix:
        ch=L[0]
        for num in L[1]:
            for n in range(len(liste)):
                if liste[n]["ch"]==ch and liste[n]["num"]==num :
                    print("exercise found : n° ",liste[n])
                    fichier.write("%2i %2i \t\t # %s \n"%(ch,num,liste[n]["com"]))
    fichier.write("# end \n")        
    fichier.close()
    print("liste_exercices.in has been generated")

def run_EC(SourcePath,RunPath,option="compilation"):
    """
    compilation
    run
    """
    if option=="compilation":
        os.chdir(SourcePath)
        #print("fichiers = ",os.listdir("./"))
        print("Compilation de l'outil EC ...")
        os.system("make clean; make")
    elif option=="run":
        os.chdir(RunPath)
        print("Lancement de l'outil EC ...")
        os.system("more choix_exercices.in")
        os.system("./x1.run_EC.sh")
        print("Solutions présentes : ")
        os.system("ls Solution*")
        print("fichiers = ",os.listdir("./"))



def test_local(task):
    """
    test en local de compilation et de lancement
    """
    print("dossiers = ",os.listdir("./"))
    print("Courant ",os.getcwd())
    if task==0:
        # compilation
        os.chdir("../Source")
        os.system("make clean; make")
    elif task==1:
        os.chdir("../RUN")
        os.system("./tp3.bin")


def run_fortran(task="liste",choix=[[10,1]],SrcDir=FortranPath,RunDir=FortranRun):
    """
    Définition des différentes tâches
    """
   
    Src=False
    Run=False

    if os.path.isdir(SrcDir) :
         print("path of Fortran sources : ", SrcDir)
         Src=True
    else:
        print("le dossier source ' %s ' n'existe pas"%(SrcDir))

    if os.path.isdir(RunDir) :
         print("path of Fortran run : ", RunDir)
         Run=True
    else:
        print("le dossier d'exécution du Fortran ' %s ' n'existe pas"%(RunDir))

    if Src:
        print("liste des exercices     : ",choix)
    
    if Run:    
        choix_exos(RunDir,choix,liste_exo_fortran(False)) # génération du fichier liste_exercices.in

    if task=="liste" and Src:
        liste=liste_exo_fortran(display=True)
    elif task=="compilation" and Src :
        run_EC(SrcDir,RunDir,"compilation")
    elif task=="run" and Run:
        run_EC(SrcDir,RunDir,"run")
    elif task=="compil+run" and Src and Run:
        run_EC(SrcDir,RunDir,"compilation")
        run_EC(SrcDir,RunDir,"run")




