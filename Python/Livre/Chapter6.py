#!/bin/python
"""
  Correction des exercices du chapitre 5
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import os

from CompressibleFlow.fonctions import *
from Tools.misc                 import *

Hauteur=8; Largeur=10
 
def Helmbold(lamb,k=np.pi):
    """
    Loi de  Helmbold : d C_L/ d \alpha
    """
    x= 2*k/(np.pi*lamb)
    return 2*k/(x+np.sqrt(1+x**2))
 

def Exercice6_1():
    """ 
    Loi du coefficient de portance avec l'allongement
    """
    set_title("Loi du coefficient de portance avec l'allongement")

    print('Loi de  Helmbold')
    
    N          = 101                              # nombre de point sur l'axe lambda  
    lambda_Max = 25                              # valeur maximale pour le dessin
    lamb       = np.linspace(0.1,lambda_Max,N)     # allongement
    k          = np.array([2.5, 2.8, np.pi]) 
    
    # plot
    fig = plt.figure(1,figsize=(Largeur,Hauteur))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$\lambda$', fontsize=16)
    ax.set_ylabel(r'$C_{L_\alpha}$',fontsize=20)
    ax.set_title('Loi de Helmbold')
    lamb1=np.linspace(0,6)
    lamb2=np.linspace(5,lambda_Max)
    for k in k:
        ax.plot(lamb, Helmbold(lamb,k), linewidth=2,label=r'$k=$%2.1f'%(k))
        ax.plot([0,lambda_Max],[2*k, 2*k], linewidth=2,label=r'$k=$%2.1f'%(k))

    ax.plot(lamb1,np.pi*lamb1/2,       color='r', linestyle='-', linewidth=2,label=r'$\pi \lambda/2$')
    ax.plot(lamb2,np.pi*2*(1-2/lamb2),linewidth=2,label=r'approx')
     
    xmax,ymax=lambda_Max,8
    ax.axis([0, xmax, 0, ymax])                                   
    xmajor_ticks,xminor_ticks = np.arange(0, xmax, 2),np.arange(0, xmax, 0.5)                                              
    ymajor_ticks,yminor_ticks = np.arange(0,ymax, 1),np.arange(0,ymax, 1)                                              
    ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True)  
    ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
    ax.grid();
    ax.legend(loc='lower right')                                          
    plt.show()


def Kp_Schlichting(x, epsilon):
    """
    Coefficient de pression sur l'ellipsoide, loi de Schlichting
    """
    # cf. http://docs.desktop.aero/appliedaero/potential3d/Bodies.html
    c = np.sqrt(1-(2*epsilon)**2)           # epsilon = a/L (rayon/longueur)
    eta = 2*(1-c**2)/c**3*(np.arctanh(c)-c)
    U_max = 2./(2-eta)
    return 1 - U_max**2 * (1-(2*x-1)**2)/(1-c**2*(2*x-1)**2)

def Kp_VD(x,epsilon):
    """
    solution analytique
    """   
    return 4*epsilon**2*np.log(epsilon**2) +epsilon**2*(2-(1-2*x)**2)/(x*(1-x))

def compare_van_Driest():
    """
     Comparaison avec les valeurs expérimentales fournies par Van Driest
    """
    #x=np.array([0.05,0.17,0.36,0.58])
    #Kp=np.array([-0.010,-0.053,-0.061,-0.026])
    x =np.array([ 0.04944156, 0.09877804,  0.17424668,  0.26218079,  0.37486144])
    Kp = np.array([-0.01011674, -0.03771612, -0.05568793, -0.05826931, -0.06136363])
    leg=["A","B","C","D","E"]
    return x,Kp,leg

def Kp_min_analytique(epsilon):
    """
    - Kp minimal pour la solution asymptotique
    """
    return np.abs(8*epsilon**2*(1+np.log(epsilon)))

def Exercice6_3():
    """
    Ellipsoide
    """
    n=121;
    x=np.linspace(0,1,n)
    eps=np.array([0.05,0.1,0.2])

    set_question("2: Kp : solution exacte et asymptotique")
     
    
    for ifig in range(len(eps)):
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'$x/\ell$', fontsize=16)
        ax.set_ylabel(r'$Kp$',fontsize=20)
        ax.set_title(r'$Kp$ : analytique et Schlichting, $\varepsilon$ = %2.2f'%(eps[ifig]))
        ax.axis([0, 1, -0.4, 1])   
        ax.plot(x,Kp_Schlichting(x,eps[ifig]), color='k', linestyle='-', linewidth=2,label=r'$Schlichting$')  
        ax.plot(x,Kp_VD(x,eps[ifig]), color='k', linestyle='--', linewidth=2,label=r'$Asymptotique$')  
        ax.legend(loc='upper center')  
     

    set_question("3: Kp min: solution exacte et asymptotique")
    
    ifig+=1

    epsil=np.logspace(-3,-1,10)
    print(epsil)
    Kpmin=[]
    for eps in epsil:
        Kpmin.append(-Kp_Schlichting(0.5,eps))

    print(Kpmin)

    fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$ \varepsilon=h/\ell$', fontsize=16)
    ax.set_ylabel(r'$-Kp_{\min}$',fontsize=20)
    ax.set_title(r'$Kp$ : analytique et Schlichting')   
    ax.loglog(epsil,Kpmin, color='k', linestyle='-', linewidth=2,label=r'$Schlichting$')  
    ax.loglog(epsil,Kp_min_analytique(epsil), color='k', linestyle='--', linewidth=2,label=r'$Asymptotique$')  
    ax.grid(which='both')
    ax.legend(loc='upper center')  

    set_question("4: Kp, expérience de Van Driest")
    
    ifig+=1
    eps0=0.075
    x0,kp0,leg0=compare_van_Driest()
    print(x0,kp0)
    
    fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$x/\ell$', fontsize=16)
    ax.set_ylabel(r'$Kp$',fontsize=20)
    ax.set_title(r'$Kp$ : analytique  et expérience , $\varepsilon$ = %2.3f'%(eps0))
    for k in range(len(x0)):
       ax.plot(x0[k],kp0[k], linestyle='',marker='o', markersize=6,label=r'Exp. point :%s'%(leg0[k])) 
    ax.plot(x,Kp_Schlichting(x,eps0), color='k', linestyle='-', linewidth=2,label=r'$Schlichting$')
    ax.axis([0, 1, -0.1, 0.1])        
    ax.plot(x,Kp_VD(x,eps0), color='k', linestyle='--', linewidth=2,label=r'$Asymptotique$')  
    ax.legend(loc='upper center')  
    ax.grid()

    plt.show()
    
   
 
def Exercice6_4():
    """
    Analyse des données de la thèse de Van Driest  sur l'ellipsoide
    Ce n'est pas un exercice du livre.
    """
    set_info("Ce n'est pas un exercice du livre ! \n" +
        "Analyse de la thèse de Van Driest only !")
    ifig=0
    l=100
    h=15
    fichiers=["vanDriest_exp1.dat","vanDriest_exp2.dat","vanDriest_exp2_ogive.dat"]
    
    def lire_fichier(nom):
        """
        lecture des fichiers
        """
        filepath=os.path.join('Livre/Data', nom)
        print('Fichier à lire : ',filepath)
        with open(filepath, 'r') as infile:
            A=np.loadtxt(infile, dtype=float,unpack=False,skiprows=1)
        d=A.shape
        print(A.shape)
        print(A)
        return A[:,0],A[:,1]

    def affichage_courbe(nfig,x,y,titrex,titrey):
        """
        dessin des courbes
        """
        
        plt.figure(num=nfig,figsize=(Largeur,Hauteur))
        plt.xlabel(titrex, fontsize=16)
        plt.ylabel(titrey,fontsize=20)
        plt.plot(x,y, 'ko', linestyle='--',markersize=5, linewidth=2,label='van Driest exp.')  
        plt.grid(which='both')
        plt.legend(loc='best')  
        
    
    x0,Kp=lire_fichier(fichiers[0])
    affichage_courbe(ifig,x0,Kp,r"$X/\ell$",r"$Kp$")

    ifig+=1
    x1,y=lire_fichier(fichiers[1])
    x2,yo1=lire_fichier(fichiers[2])
    affichage_courbe(ifig,x1,y,r"$X$",r"$y_{ogive}$")
    affichage_courbe(ifig,x2,yo1,r"$X$",r"$y_{ellipsoide}$")
    
    

    def ogive_geometrie(x,l,h):
        """
        définition de la géométrie de l'ogive
        """
        y=h*np.sqrt(1-(x/l)**2)
        return y

    X=np.linspace(-l,l,101)
    Y=ogive_geometrie(X,l,h)
    
    plt.figure(ifig,figsize=(Largeur,Hauteur)) 
    plt.plot(X,Y, 'r-', linestyle='--',markersize=5, linewidth=2,label='géométrie de référence')  
    plt.plot(x0*l,ogive_geometrie(x0*l,l,h), 'bs', linestyle='--',markersize=7, linewidth=2,\
             label='Points expé.') 
    plt.legend(loc='lower center')
 
    xf=(1+x0)/2
    yf=Kp
    ifig+=1
    affichage_courbe(ifig,xf,yf,r"$x/c$",r"$Kp$")
    for x,kp in zip(xf,Kp):
        print("x = %2.3f,\t Kp = %2.4f"%(x,kp))
 
    plt.show()
    #data=np.zeros([d[0],d[1]+1])
