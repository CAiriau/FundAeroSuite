#!/bin/python
"""
    Correction des exercices du chapitre 12 suite
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from copy import deepcopy
from CompressibleFlow.fonctions import *
from CompressibleFlow.tables    import *
from Tools.misc                 import * 
 
 
# Variables globales

Largeur,Hauteur=10,8
Hequal=4

def Exercice12_3():
    """
    Tuyère avec une paroi sous forme polynomiale
    c'est une géométrie pour l'exercice 1.
    """
    opt_paroi=2
    n=11
    if opt_paroi==1:
    	y0=1.
    	L=5*y0
    	h=2*y0
    	hbar=h/L
    	y0bar=y0/L
    	x=np.linspace(0,L,n)
    	p0=np.deg2rad(10)
    else:
    	y0bar=1
    	p0=0.05
    	hbar=0
    	L=1.
    	x=np.linspace(1,8.5,n)**4/100000


    def paroi(opt,eta,y0bar,p0,hbar):
    	"""
    	Loi de paroi
    	eta = x/L, hbar=h/L
    	"""
    	if opt==1:
    		y= y0bar+p0*eta-(2*p0-3*hbar)*eta**2+(p0-2*hbar)*eta**3
    	else:
    		y=p0/2*eta**2+y0bar
    	return y

    fig = plt.figure(1,figsize=(Largeur,Hequal))
    fig.suptitle('Tuyère polynomiale',fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$y/L$',fontsize=14); ax.set_xlabel(r'$x/L$',fontsize=14)
    ax.plot(x/L, paroi(opt_paroi,x/L,y0bar,p0,hbar),'k-')
    plt.grid()
    plt.show()

    print(np.rad2deg(p0*x/L))
