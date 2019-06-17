#!/bin/python
"""
  Correction des exercices du chapitre 7 : test méthode des panneaux
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 



import os
import numpy as np
import matplotlib.pyplot as plt

from CompressibleFlow.fonctions import *
from IncompressibleFlow.PanelMethod import *
from scipy import integrate
from scipy.interpolate import interp1d
from IncompressibleFlow.Airfoil import *
from Tools.misc                     import *

#********************************************************************
def Exercice7_30():
#********************************************************************
    """ 
    Kp sur un profil naca , 
    """
    set_title("Solution de référence sur un profil Naca")

    # chargement des données pour le profil NACA0012
    reference_filepath = os.path.join('Livre/Data', 'NACA0012.dat')
    print('fichier de référence chargé : ',reference_filepath)
    with open(reference_filepath, 'r') as infile:
        x,v_2,v,Delta_v=np.loadtxt(infile, dtype=float, unpack=True,skiprows=1)
    print('x =' ,x)
    chord=1.0
    x=x/100.0
    
    alpha_ref=9.00              # incidence correspondant aux données du Kp
    alpha=alpha_ref                   # incidence de calcul
    c2=alpha/alpha_ref
    Kpe=1-(v +c2*Delta_v)**2
    Kpi=1-(v -c2*Delta_v)**2
    KpRef=1-v_2
    DeltaKp=Kpi-Kpe             #  chargement aérodynamique

    # Interpolation pour un calcul plus précis du coefficient de portance
    s=np.linspace(0,1,1001); f1 = interp1d(x, DeltaKp, kind='cubic')
    DeltaKp1=f1(s)
    Cl1=integrate.trapz(DeltaKp1,s)
    alpha1=np.rad2deg(Cl1/(2*np.pi))
    print("Interpolation  : Cl = %f , = alpha = %f"%(Cl1,alpha1))

    Cl=integrate.trapz(DeltaKp,x)
    alpha=np.rad2deg(Cl/(2*np.pi))
    print("Données brutes : Cl = %f , = alpha= %f"%(Cl,alpha))

    e=-1  # suivant qu'on trace le Kp ou -Kp
    width = 10;ech=0.1
    plt.figure(figsize=(width, width))
    plt.xlabel('x', fontsize=16);plt.ylabel('Kp', fontsize=16)
    plt.title('Coefficient de pression')
    plt.plot(x,e*Kpe,linestyle='-', linewidth=2, marker='o', markersize=9, color='red',label='Kp extrados')
    plt.plot(x,e*Kpi,linestyle='-', linewidth =2, marker='o', markersize=4, color='blue',label='Kp intrados')
    plt.plot(x,e*KpRef,linestyle='-', linewidth=2, marker='o', markersize=4, color='green',label=r'Kp ref $(\alpha=0)$')
    plt.plot(x,v,linestyle='--', linewidth=2, color='black',label=r'$v/V$')
    plt.plot(x,v_2,linestyle='-', linewidth=2, color='black',label=r'$v^2/V^2$ (épaisseur)')
    plt.plot(x,Delta_v,linestyle='--', linewidth=2,marker='o', markersize=3,  color='red',label=r'$\Delta v/V$ (incidence)')
    plt.plot(x,DeltaKp,linestyle='--', linewidth=2, marker='s', markersize=3, color='red',label=r'$\Delta Kp$')
    plt.plot(s,DeltaKp1,linestyle='--', linewidth=2,color='Green',label=r'$\Delta Kp$ (interpo.)')
    plt.legend(loc='upper right');plt.grid()
    plt.xlim(0-ech, chord+ech)
    #plt.ylim(-2.1, 1.1);
    plt.show()
