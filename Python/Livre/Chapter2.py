#!/bin/python
"""
  Correction des exercices du chapitre 2
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 



import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import os
from CompressibleFlow.fonctions import *
from CompressibleFlow.tables    import *
from Tools.misc                 import *

twopi=2*np.pi

def streamfunction_source(r,theta,C,cas):
    """
    fonction de courant pour le point source
    retourne psi, phi
    Attention certaines fonctions ne marchent pas bien...
    il vaut mieux uniquement utiliser les cas 4 et 5 
    """

    x,y=r*np.cos(theta),r*np.sin(theta)
    if (cas==1):  # tourbillon ponctuel, validé mais ne fonctionne pas très bien.
        return C/2*np.log(r),C*theta,x,y
    elif (cas==2): # doublet source-puits, validé mais ne fonctionne pas très bien.
        return -C/twopi*np.sin(theta)/r,C/twopi*np.cos(theta)/r,x,y
    elif (cas==3):  # étirement pur, validé et fontionne correctement.
        return C*x*y,-C/2*np.log(x**2-y**2),x,y
    elif (cas==0): # source/puit, validé mais ne fonctionne pas très bien.
        return C*theta, C/2*np.log(r),x,y
    elif (cas==4):  # étirement pur,  validé et fonctionne très bien
        return C*x*y, C/2.*(x**2-y**2),x,y
    elif (cas==5):  # doublet source puits d'axe horizontal, validé et fonctionne très bien 
        return  C*np.sin(theta)/r, C*np.cos(theta)/r,x,y

def Exercice2_1():
    """ 
    fonctions de courant et iso potentielles
    """
    set_title("Tracé de fonctions de courant et d'isopotentielles")
    opt=2    # 1 : projection polaire, autre : projection cartésienne.

    cas=4
    Theta = np.radians(np.linspace(0, 360, 63))
    R = np.linspace(0.0, 4, 53)
    r, theta = np.meshgrid(R, Theta)
    psi,phi,x,y=streamfunction_source(r,theta,1.,cas)
    levels_psi=np.arange(-3.11,3.11,0.2)
    levels_phi=np.arange(-3.11,3.11,0.21)
    print('Min psi = %f, \t Max psi = %f '%(np.min(psi),np.max(psi)))
    print('Min phi = %f, \t Max phi = %f '%(np.min(phi),np.max(phi)))

    if opt==1:
        # dans un graphe polaire
        fig1, ax1 = plt.subplots(subplot_kw=dict(projection='polar'))
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        ax1.contour(theta, r, psi,levels_psi,colors='k',linewidths=2)
        #fig2, ax2 = plt.subplots(subplot_kw=dict(projection='polar'))
        matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
        ax1.contour(theta, r, phi,levels_phi,colors='r')
       
    else:
        # dans un graphe cartésien
        fig1, ax1 = plt.subplots()
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
        ax1.contour(x, y, psi,levels_psi,colors='k',linewidths=2)
        #matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
        ax1.contour(x, y, phi,levels_phi,colors='r')
        ax1.axis('equal')
    plt.show()
    return


def Exercice2_2():
    """ 
    Atmosphère Standard
    """
    set_title("Atmosphère Standard")
      
    Z    = np.array([10,10.5,18])
    Mach = np.array([0.8,0.82,2.05])
    p0   = 101325
    T0   = 288.15
    rho0 = 1.225 
  
    for (h,M) in zip(Z,Mach):
        sigma,delta,Theta=Atmosphere(h)
        T=Theta*T0
        a=sound_velocity(T)
        U0=a*M
        print("Z = %2.1f km, T = %3.2f K, p = %5.1f Pa, rho = %1.4f kg/m^3, a = %3.2f m/s"%(h,T,p0*delta,rho0*sigma,a))
        print('U0 = %3.2f m/s = %3.2f km/h'%(U0,U0*3.6))


    print("\n Données: \n")
    R=8.31432
    m_air=28.9644e-3
    r=R/m_air
    g0=9.80665
    g_sur_r=34.1631947

    print("r = %f, 1000 g0/r = %f = %f"%(r,1000*g0/r,g_sur_r))