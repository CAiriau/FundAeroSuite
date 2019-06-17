#!/bin/python
"""
  Correction des exercices du chapitre 13
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
 

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Tools.misc                 import *
from CompressibleFlow.fonctions import * 
from scipy.interpolate import interp1d
from Choc_Conique.choc_conique  import * 

#************************************
# main program
#************************************

def Kp(M,theta_c):
    """
    calcul du Kp
    """
    beta=np.sqrt(M**2-1)
    Kp_exact=-theta_c**2*(
        2.0*np.log(beta*theta_c/(1.0+np.sqrt(1.0-(beta*theta_c)**2)))+1)
    Kp_approx=-theta_c**2*(2.0*np.log(beta*theta_c/2.0)+1)
    return Kp_exact,Kp_approx    


def Exercice13_1():
    """ 
    Obstacle conique en supersonique 
    """

    plot=True
         
    nMach=100
    set_title("Obstacle conique en supersonique avec la théorie des corps élancés")
    k=2        # choix de l'angle du cone
    theta_c=np.deg2rad(np.array([5.0, 10.0, 15.0]))
    print('Theta = ',theta_c)
      

    Kp_exact=np.zeros((len(theta_c),nMach))
    Kp_approx=np.zeros((len(theta_c),nMach))
    M0=np.zeros((len(theta_c),nMach))
    for k in np.arange(len(theta_c)):
        Kpe,Kpa=Kp(2.0,theta_c[k])
        Mach_maxi=np.sqrt(1.0/theta_c[k]**2-1)
        print('M0=2, Kp exact = %f, Kp_approx = %f for theta =  %f °, Mach Maxi = %f'%(Kpe,Kpa,np.rad2deg(theta_c[k]),Mach_maxi)) 
        M0[k,:]=np.linspace(1.1,Mach_maxi,nMach)             # Mach number 
        Kp_exact[k,:],Kp_approx[k,:] =Kp(M0[k,:],theta_c[k])    

    
        
    if plot:
        fig = plt.figure()
        fig.suptitle('Kp sur un cone', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('critère de Michel')
        ax.set_xlabel(r'$M_0$',fontsize=20)
        ax.set_ylabel(r'$Kp$',fontsize=20)
        ax.grid()
        ax.axis([1.0, 12, -0.05, 0.35])
        for k in np.arange(len(theta_c)):
            ax.plot(M0[k,:],Kp_exact[k,:],'-',label='Exact',linewidth=2,color='black') 
            ax.plot(M0[k,:],Kp_approx[k,:],'-',label='Approché',linewidth=2,color='blue') 
        ax.legend(loc='upper right')
        plt.show() 

    # calcul du Kp pour un choc conique, interpolation à partir d'un fichier

    reference_filepath = os.path.join('Livre/Data', 'Mach_2_choc_conique.dat')
    #reference_filepath = os.path.join('Livre/Data', 'test.dat')
    with open(reference_filepath, 'r') as infile:
        A=np.loadtxt(infile, dtype=float, unpack=True,skiprows=2)
        theta=A[0,1:200];kp_cc=A[2,1:200]
    print(A.shape)
    f1 = interp1d(theta,kp_cc, kind='cubic')
    print('pour les 3 angles', np.rad2deg(theta_c), '\t  Valeur des Kp  suivant les chocs coniques : = ',f1(np.rad2deg(theta_c)))


def Exercice13_2():
    """ 
    Ogive parabolique en supersonique 
    """
    plot=True
     
    hat_d=0.1   # d/L
    Mach=np.sqrt(2.0)      # Mach
    n=100
    eps=0.001
    x=np.linspace(eps,1.0-eps,n)
 
    def kp_ogive(X):
        """    
        Kp sur l'ogive
        """

        beta=np.sqrt(Mach**2-1)
        hat_R=2.*hat_d*X*(1.-X)
        var=beta*hat_R
        print(X/var)
        #f=np.log(X/var+np.sqrt((X/var)**2-1))
        f=np.arccosh(X/var)
        u_x=-2.*hat_d**2*((12*X*(X-1)+2+6*var**2)*f+ 6*(2-3*X)*np.sqrt(X**2-var**2))
        u_r=2*hat_d*(1-2*X)
        Kp=-2*u_x-u_r**2
        Kp1=-2*u_x 
        #delta=2*hat_d*(1-2*X)
        #Kp_ackeret=2*delta/beta
        return hat_R,u_x,u_r,Kp,Kp1  

    R,u_x,u_r,Kp,Kp1=kp_ogive(x)

    if plot:
        fig = plt.figure(1)
        fig.suptitle("Kp sur l'ogive", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('critère de Michel')
        ax.set_xlabel(r'$X$',fontsize=20)
        ax.set_ylabel(r'$Kp$',fontsize=20)
        ax.grid()

        #ax.axis([0, 1, -0.1, 0.4])
        ax.plot(x,Kp,'-',label='Exact',linewidth=2,color='black')    
        ax.plot(x,Kp1,'-',label=r'sans $u_r$',linewidth=2,color='blue')      
            
        ax.legend(loc='center')

        fig = plt.figure(2)
        fig.suptitle(r"vitesses sur l'ogive pour $\hat d$ = %3.2f"%(hat_d), fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('critère de Michel')
        ax.set_xlabel(r'$X$',fontsize=20)
        ax.set_ylabel(r'$u_r,u_x$',fontsize=20)
        ax.grid()
        #ax.axis([1.0, 12, -0.05, 0.35])
        ax.plot(x,u_x,'-',label=r'$u_x$',linewidth=2,color='black')    
        ax.plot(x,u_r,'-',label=r'$u_r$',linewidth=2,color='blue')      
            
        ax.legend(loc='upper right')

        fig = plt.figure(3)
        fig.suptitle("forme de l'ogive", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('critère de Michel')
        ax.set_xlabel(r'$X$',fontsize=20)
        ax.set_ylabel(r'$Y$',fontsize=20)
        ax.grid()
        ax.axis('equal')
        #ax.axis([1.0, 12, -0.05, 0.35])
        ax.plot(x,R,'-',label=r'$R$',linewidth=2,color='black')  
        ax.plot(x,-R,'-',label=r'$-R$',linewidth=2,color='black')     

        fig = plt.figure(4)
        fig.suptitle("Kp sur l'ogive + vitesses", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('critère de Michel')
        ax.set_xlabel(r'$X$',fontsize=20)
        ax.set_ylabel(r'$Kp$',fontsize=20)
        ax.grid()

        ax.axis([0, 1, -0.2, 0.4])
        ax.plot(x,Kp,'-',label='Exact',linewidth=2,color='black')    
        ax.plot(x,Kp1,'-',label=r'sans $u_r$',linewidth=2,color='blue') 
        #ax.plot(x,Kp_ackeret,'-',label=r'sans $u_r$',linewidth=2,color='green') 
        ax.plot(x,u_x,'-',label=r'$u_x$',linewidth=2,color='red')    
        ax.plot(x,u_r,'-',label=r'$u_r$',linewidth=2,color='magenta')   
        ax.plot(x,R,'-',label=r'$R$',linewidth=1,color='black')  
        ax.plot(x,-R,'-',label=r'$-R$',linewidth=1,color='black')       
            
        ax.legend(loc='upper center')        
        plt.show() 

def sears1(theta):
    """
    corps optimal avec volume fixé, S/S0
    """
    return np.sin(theta)**3
def sears2(theta):
    """
    corps optimal avec la section fixée
    """
    return np.sin(theta)-np.cos(theta)**2*np.log((1+np.sin(theta))/np.abs(np.cos(theta)))


def Exercice13_3():
    """ 
    Traînée minimale d'un corps élancé  supersonique 
    """
    set_title("Traînée minimale d'un corps élancé  supersonique")
    plot=True
    n=100
    ep=0.1
    t=np.linspace(0,np.pi,n)
    theta=np.flipud(t)
    x=(1+np.cos(theta))/2
   
    if n < 201:
        print('x= ',x)
        print(sears1(theta))
        print(sears2(theta))

    if plot:
        print('plot')
        fig = plt.figure(figsize=(10,3))
        fig.suptitle("Ogives optimales de Sears", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x/\ell$',fontsize=20)
        ax.set_ylabel(r'$S/S_0$',fontsize=20)
        ax.grid()
        ax.axis('equal')
        #ax.axis([1.0, 12, -0.05, 0.35])
        ax.plot(x,ep*np.sqrt(sears1(theta)),'-',label=r'$V$ fixé',linewidth=2,color='black')    
        ax.plot(x,ep*np.sqrt(sears2(theta)),'-',label=r'$S$ fixé',linewidth=2,color='blue')  
        ax.plot(x,-ep*np.sqrt(sears1(theta)),'-',label=r'$V$ fixé',linewidth=2,color='black')    
        ax.plot(x,-ep*np.sqrt(sears2(theta)),'-',label=r'$S$ fixé',linewidth=2,color='blue')   
        ax.legend(loc='upper right')        
        plt.show()

def Plot_iso_potential(Mach):
    """
    Tracé des isopotentielles en fonction du nombre de Mach
    """
    # paramètres

    Rmin    = 0.5              # rayon minimal pour les coordonnées polaires
    Rmax    = 10.0             # rayon maximal pour les coordonnées polaires
    m       = 15               # nombre de lignes de courant ou d'quipotentielles (environ)
    n       = 73               # nombre de points pour définir les lignes de courants ou équipotentielles
                               # un diviseur de 360+1 est préférable mais pas obligatoire :  37, 73, 181, 361 
       
    fig = plt.figure(figsize=(10,10))
    fig.suptitle(r"iso - $\Psi$,iso - $\phi$ pour $M=$ %3.2f "%Mach, fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    #fig.subplots_adjust(top=0.80)
    #ax.set_xlabel(r'$x$',fontsize=20)
    #ax.set_ylabel(r'$y$',fontsize=20)
    #ax.axis('equal')
    ech=1.1
   
    

    #****************************
    #  INCOMPRESSIBLE
    #****************************
    if Mach==0:
        beta=1.0
        theta=np.deg2rad(np.linspace(0,360,n))
        thetap=np.deg2rad(np.linspace(-90,90,m))
        R=np.linspace(Rmin,Rmax,m)
        for i in range(m):
            # lignes de courant
            ax.plot(np.cos(theta)*R[i],np.sin(theta)*R[i],'-',linewidth=1,color='black')  
        # équipotentielles
        a,b=np.cos(thetap)*Rmax,np.sin(thetap)*Rmax
        ax.plot([-a,a],[-b,b],'r--',linewidth=1,color='black')
        ax.axis("equal") 

    #****************************
    #  SUBSONIQUE COMPRESSIBLE
    #****************************
    elif Mach < 1:
        beta2=1.0-Mach**2
        beta=np.sqrt(beta2)
        print('beta=',beta)
        theta=np.deg2rad(np.linspace(0,360,n))                          # pour les lignes de courant
        # Vecteur t adapté au cas M=0.75
        t=np.deg2rad(np.array([25,40,50,57,65,71,75,78,80,82.5,85]))    # pour les iso potentielles
        
        thetap=np.concatenate((-np.flipud(t),t))
        p=len(thetap)
        R=np.linspace(Rmin,Rmax,m)                                  # nombre de ligne de courant
        
        # lignes de courant 
        for i in range(m):
            x=np.cos(theta)*R[i]
            y=np.sin(theta)*R[i]/beta
            ax.plot(x,y,'-',linewidth=1,color='black')  

        # équipotentielles
        eta=Rmax*np.cos(np.deg2rad(np.linspace(0,180,n)))                       # définitions des iso potentielles       
        for i in range(p): 
            phi=np.sign(eta)*np.tan(thetap[i])*pow(np.abs(eta),beta2)
            ax.plot(eta,phi,'--',linewidth=1,color='black')
        ax.plot([0,0],[-r,r],'--',linewidth=1,color='black')
        ax.plot([-Rmax,Rmax],[0,0],'--',linewidth=1,color='black')
        
        #ax.axis([-r,r,-r,r])
        ax.axis("equal")
        ax.set_ylim(-ech*Rmax/beta, ech*Rmax/beta)
        
    #****************************
    #  SUPERSONIQUE
    #****************************
    else:
        beta2=Mach**2-1.0
        beta=np.sqrt(beta2)
        pente=1/beta
        print("pente des lignes de courant sur l'asymptote ° = %4.2f"%(np.rad2deg(np.arctan(pente))))
        
        # Asymptote
        x_0=np.array([0 , Rmax])
        y_sup=np.array([0, pente*Rmax])
        y_inf=-y_sup
        hori=np.array([0,0])

        # hyperbole, ligne de courant
        R=np.linspace(Rmin,Rmax/beta,m)
        eps=0.00001   		# pour éviter les mauvais arrondis
        ax.plot(x_0,y_sup,'-',linewidth=1,color='black')
        ax.plot(x_0,y_inf,'-',linewidth=1,color='black')
        ax.plot(x_0,hori,'-',linewidth=1,color='black')
        for i in range(m):
        	if beta*R[i] < Rmax+eps:
        		x=np.linspace(beta*R[i]+eps,Rmax+eps,n)
        		y=np.sqrt(x**2/beta2-R[i]**2)
        		ax.plot(x,y,'-',linewidth=1,color='black')
        		ax.plot(x,-y,'-',linewidth=1,color='black')
        	else:
        		print("Problème")

        # Equipotentielles
        
        print("pente des équipotentielles sur l'asymptote en ° = %4.2f"%(-np.rad2deg(np.arctan(beta))))
        x0=np.linspace(Rmin,Rmax-0.2,m)
        # cône de Mach pour la définition de la zone de démarrage des équipotentielles sur l'asymptote)
        #ax.plot(x0, pente*x0,'-',linewidth=2,color='blue')
        #ax.plot(x0,-pente*x0,'-',linewidth=2,color='blue')
         
        for i in range(m):
            eta=np.linspace(x0[i],Rmax,n)
            phi=x0[i]/beta*pow(x0[i]/eta,beta2)
            ax.plot(eta, phi,'--',linewidth=1,color='black')
            ax.plot(eta,-phi,'--',linewidth=1,color='black')

        ax.axis("equal")     
     
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(True)
    #ax.axis('off')
    plt.show()
      

def Exercice13_4():
    """
    Tracé des iso-potentielles et lignes de courant
    figurer 13.2 du livre "Aérodynamique Fondamentale"
    """
    Mach = [0,0.75,2.0]           # on entre le Mach qu'on veut
    #Mach = [2]
    for Mach in Mach:
        Plot_iso_potential(Mach)
     

def Exercice13_5():
    """
    Corps conique élancé (ou pas)
    """    
    set_title("Choc conique sur un corps conique élancé (ou pas)")
    
    def Run_Choc_Conique(option):
    
        """
        Choc conique sur un corps conique élancé (ou pas)
        """ 
        if option==0:
            """
            Calcul d'un Mach unique, tant qu'il existe le fichier avec ce Mach.
            """
            main_single_Mach(M0=2.0,theta_c=np.float64(15.0),valeurs=True)

        elif option==1:
            """
            Mach aval, Kp, sigma en fonction de theta_c
            """
            main_multiple_Mach(np.array([1.05,1.1,1.3,1.5,1.6,2.0,5]),["sigma","Mach_aval","Kp"],Lees=False,
                t=[55,10,5],s=[70,10,5],Kp=[140,20,10],Mav=[6,1,0.5])
        elif option==2:
            """
             theta_x_Max, Kp correspondant en fonction du Mach
            """
            main_angle_deviation_maxi(MachMax=7.1)
        else:
            """
            Kp en fonction du Mach, pour un angle de cône donné
            """
            main_plot_Mach_Kp(np.array([1.055,7.]),10.)


    
    set_question("Theta Max")
    Run_Choc_Conique(2)
    set_question("Kp")
    Run_Choc_Conique(3)
    #Run_Choc_Conique(1)


