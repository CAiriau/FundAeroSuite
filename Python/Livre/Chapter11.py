#!/bin/python
"""
    Correction des exercices du chapitre 11
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


import numpy as np
import matplotlib.pyplot as plt
from CompressibleFlow.fonctions import *
from CompressibleFlow.tables    import *
from Tools.misc                 import * 
 
#************************************
# main program
#********************* 

 
Largeur,Hauteur=10,8

    

def Exercice11_1():
    """ 
    Courbe p1/p0 en fonction d'une déviation pour plusieurs Mach
    """
    ifig=0;
    
    set_title("Courbes p1/p0 en fonction de la déviation")
     
    n=101
    Liste_Mach=np.array([1.3,1.6,2.0,2.3,2.5,2.7,3.0])
    #Liste_Mach=np.array([2.0,3.0])

    set_question('1 : Compression par choc')
     

    print('Mach  =', Liste_Mach)
     
    theta,r_p=[],[]
    for Mach in Liste_Mach:
        sigma=np.linspace(np.arcsin(1/Mach),np.pi/2,n)
        theta.append(theta_from_Mn0(sigma,Mach))
        r_p.append(P2_P1(Mach*np.sin(sigma)))
    
    fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
    fig.suptitle('Compression par choc', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    #fig.subplots_adjust(top=0.80)
    ax.set_ylabel(r'$P_2/P_1$',fontsize=20)
    ax.set_xlabel(r'$\theta(^\circ)$',fontsize=20)
    for k,Mach in enumerate(Liste_Mach):
        ax.plot(np.rad2deg(theta[k]),r_p[k],linewidth=2,label=r'$M_0=$%6.5f'%(Mach))
    ax.axis([0,35,1,11])
    xmajor_ticks,xminor_ticks = np.arange(0, 35, 10),np.arange(0, 35, 2)                                              
    ymajor_ticks,yminor_ticks = np.arange(1,12, 2),np.arange(1,12, 0.5)                                              
    ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True) 
    ax.tick_params(labelsize=16)
    ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
    ax.grid(which='both')                 
    ax.legend(loc='lower right')
   

    set_question('2 : Détente isentropique')
     
    theta=np.linspace(0,np.pi/4,n)
    
    r_p=[]
    for Mach in Liste_Mach:
        r_p.append(detente_isentropique(theta,Mach,gamma=gamma))
    ifig=ifig+1
    fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
    fig.suptitle('Détente isentropique', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    #fig.subplots_adjust(top=0.80)
    ax.set_ylabel(r'$P_2/P_1$',fontsize=20)
    ax.set_xlabel(r'$\theta(^\circ)$',fontsize=20)
    for k,Mach in enumerate(Liste_Mach):
        ax.plot(np.rad2deg(theta),r_p[k],linewidth=2,label=r'$M_0=$%6.5f'%(Mach))
    ax.axis([0,45,0,1])
    xmajor_ticks,xminor_ticks = np.arange(0, 45, 10),np.arange(0, 45, 2)                                              
    ymajor_ticks,yminor_ticks = np.arange(0,1.1, 0.25),np.arange(0,1, 0.05)                                              
    ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True) 
    ax.tick_params(labelsize=16)
    ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
    ax.grid(which='both')           
    ax.legend(loc='upper right')
    plt.show()

def Exercice11_2():
    """
    Réflexion d'un choc dans un canal
    """
    M0=1.5         # Mach amont des deux chocs
    theta0=5
    theta1=10
    ifig=0
    n=51
    plot_figure=True
    plot_grid=False


    thetaMin,thetaMax=-0.5,13
    dthetaMaj,dthetaMin=5,1
    rMax=2.5

    sigma0,Maval0,Paval0,omegaAmont,omegaAval=Choc_Oblique(M0,theta0)   
    sigma1,Maval1,Paval1,omegaAmont,omegaAval=Choc_Oblique(M0,theta1)  
   
    Liste_Mach=np.array([M0,Maval0,Maval1])
    #Liste_Mach=np.array([M0])
    #Liste_Mach=np.array([M0])

    # calcul des courbes en forme de coeur
    xc,yc=[],[]
    theta,r_p=[],[]
    for Mach in Liste_Mach:
        sigma=np.linspace(np.arcsin(1/Mach),np.pi/2,n)
        t=theta_from_Mn0(sigma,Mach)
        r=P2_P1(Mach*np.sin(sigma))
        theta.append(np.concatenate((-np.flipud(t),t)))
        r_p.append(np.concatenate((np.flipud(r),r)))

    # centre de ces courbes 
    xc.append(0); yc.append(1)                  # premier choc : pas de décalage
    xc.append(theta0);yc.append(Paval0)         # décalage du second choc, cas à 5°
    xc.append(theta1);yc.append(Paval1)         # décalage du second choc, cas à 10°
   
    if plot_figure :
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Interaction de chocs', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        #fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$P_2/P_1$',fontsize=20)
        ax.set_xlabel(r'$\theta(^\circ)$',fontsize=20)
        for k,Mach in enumerate(Liste_Mach):
            ax.plot(np.rad2deg(theta[k])+xc[k],r_p[k]*yc[k],linewidth=2,label=r'$M_0=$%6.5f'%(Mach))
        ax.axis([0,thetaMax,1,rMax])
        if plot_grid:
            xmajor_ticks,xminor_ticks = np.arange(thetaMin, thetaMax, dthetaMaj),np.arange(thetaMin, thetaMax, dthetaMin)                                              
            ymajor_ticks,yminor_ticks = np.arange(1,rMax, 0.5),np.arange(1,rMax, 0.1)                                              
            ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True) 
            ax.tick_params(labelsize=16)
            ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
            ax.grid(which='both')
        else:
            ax.grid()                 
        ax.legend(loc='lower right')
        plt.show()


def Exercice11_3():
    """
    Interaction de 2 chocs obliques
    """
    M0=2.0         # Mach amont des deux chocs
    theta0=5
    theta1=15
    thetag=0
    ifig=0
    n=51
    plot_figure=True
    plot_grid=True


    thetaMin,thetaMax=-0.,16
    dthetaMaj,dthetaMin=5,1
    rMax=2.4

    sigma0,Maval0,Paval0,omegaAmont,omegaAval=Choc_Oblique(M0,theta0,show=False,msg='Zone 1')   
    sigma1,Maval1,Paval1,omegaAmont,omegaAval=Choc_Oblique(Maval0,theta1-theta0,show=False,msg='Zone 2') 
    sigma2,Maval2,Paval2,omegaAmont,omegaAval=Choc_Oblique(M0,theta1-thetag,show=False,msg='Zone 4') 
   
    Liste_Mach=np.array([M0,Maval0])
    #Liste_Mach=np.array([M0])
    #Liste_Mach=np.array([M0])

    # calcul des courbes en forme de coeur
    xc,yc=[],[]
    theta,r_p=[],[]
    for Mach in Liste_Mach:
        sigma=np.linspace(np.arcsin(1/Mach),np.pi/2,n)
        t=theta_from_Mn0(sigma,Mach)
        r=P2_P1(Mach*np.sin(sigma))
        theta.append(np.concatenate((-np.flipud(t),t)))
        r_p.append(np.concatenate((np.flipud(r),r)))

    # centre de ces courbes 
    xc.append(0); yc.append(1)                  # premier choc : pas de décalage
    xc.append(theta0);yc.append(Paval0)         # décalage du second choc, cas à 5°
    #xc.append(theta1);yc.append(Paval1)         # décalage du second choc, cas à 10°
   
    if plot_figure :
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Interaction de chocs', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        #fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$P_2/P_1$',fontsize=20)
        ax.set_xlabel(r'$\theta(^\circ)$',fontsize=20)
        for k,Mach in enumerate(Liste_Mach):
            ax.plot(np.rad2deg(theta[k])+xc[k],r_p[k]*yc[k],linewidth=2,label=r'$M_0=$%6.5f'%(Mach))
        ax.axis([0,thetaMax,1,rMax])
        if plot_grid:
            xmajor_ticks,xminor_ticks = np.arange(thetaMin, thetaMax, dthetaMaj),np.arange(thetaMin, thetaMax, dthetaMin)                                              
            ymajor_ticks,yminor_ticks = np.arange(1,rMax, 0.5),np.arange(1,rMax, 0.1)                                              
            ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True) 
            ax.tick_params(labelsize=16)
            ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
            ax.grid(which='both')
        else:
            ax.grid()                 
        ax.legend(loc='lower right')
        plt.show()
