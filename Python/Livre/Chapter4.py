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
from Tools.misc                 import *
from scipy.optimize             import newton
from Livre.ConformalMapping.mappings  import *
from IncompressibleFlow.Airfoil       import *
from IncompressibleFlow.LinearTheoryProfile  import * 


twopi=2*np.pi


fz=16

def Exercice4_1():
    """
    Plaque plane en incidence
    """
    set_title("Plaque plane en incidence")
    par=set_parameters()
    s= AirfoilConfMap(par)
    s.map           = "Joukowski1"
    s.alpha         = 20 
    s.opt_levels    = "manuel"
    s.levels_iso    = [-1,1,0.2]
    s.adimChord     = True
    s.n_circle      = 101
    s.plot_velocity = False
    s.plot_circle   = False
    s.plot_airfoil  = False
    s.plot_Kp       = False
    s.run_airfoil()

def plot_profil(xe,xi,ye,yi,cas):
    """
    Dessin d'une liste de Kp
    """
    plt.figure(figsize=(16,4))
    plt.title("Comparaison des profils")
    for k in range(len(xe)):
        plt.plot(xe[k],ye[k],label="ye: "+cas[k])
        plt.plot(xi[k],yi[k],label="yi: "+cas[k])
    plt.grid()
    #plt.xlim(-2.05,2.05)

    plt.axis("equal")
    plt.legend(loc="best")
    plt.show()

 

def Exercice4_2():
    """ 
    profil de Joukowski
    """
    set_title("profil de Joukowski")

    xe,xi=[],[]
    ye,yi=[],[]
    cas=["Profil épais","Squelette"]
    par=set_parameters()
    set_question("a - Profil épais")
    s= AirfoilConfMap(par)
    # profil de Joukowski, le profil est toujours dessiné sous incidence nulle
    # dans ce cas la paramètre est Zc
    s.map           = "Joukowski1"
    s.a             = 1
    s.alpha         = 0
    s.Zc            = -0.1+1j*0
    s.plot_circle   = False
    s.plot_airfoil  = False
    s.opt_levels    = "manuel"
    s.levels_iso    = [-1,1,0.2]
    s.adimChord     = False
    s.plot_velocity = False
    s.plot_Kp       = True
    s.plot_psi      = False
    s.a             = 1
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    ye.append(s.ye)
    yi.append(s.yi)

    set_question("3 - Profil squelettique")

    s.Zc            = 0+1j*0.1     
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    ye.append(s.ye)
    yi.append(s.yi)

    plot_profil(xe,xi,ye,yi,cas)


def Exercice4_3():
    """ 
    profil famille_pq 
    """
    set_title("Profil famille_pq")
        
    set_question("1 à 4:  Profil sans incidence" )
    cas= 'analytique'
    par=set_parameters_airfoil()
    s= LinearTheoryAirfoil(par)
    s.airfoil       = ['famille_pq',"analytique","ys"]
    s.nFourier      = 20
    s.run_theory2D()

    set_question("5 à 7:  profil à l'incidence d'adaptation" )
    s.alpha         = 2.127744856489093
    s.run_theory2D()
    print("vérification p x epsilon = ",s.AirfoilParam[0]*s.AirfoilParam[1])
    print("Kp max, Kp min           = ",np.amax(s.Kp_s[2:40]),np.amin(s.Kp_s))
 

def Exercice4_4():
    """ 
    profil double cambrure +  becs et volets
    """
    set_title("profil double cambrure +  becs et volets")
        
    set_question("1- Profil" )
    cas= 'analytique'
    par=set_parameters_airfoil()
    s= LinearTheoryAirfoil(par)
    s.airfoil       = ['double courbure',"analytique","ys"]
    s.run_theory2D()

    set_question("3- Becs et Volets: volet" )

    # je rentre l'angle comme 1 rad, le problème étant linéaire je trouve les 
    # coefficients devant beta_v
    DeltaCL,DeltaCmF,DeltaAlpha_0,DeltaAlpha_a=EffetVoletBec(Volet=True,theta=0,longueur=0.25,angle=180/np.pi)

    set_question("4- Becs et Volets : volet, variation d'incidence" )
    delta_alpha=DeltaCL/(2*np.pi)
    print("variation de l'incidence $\delta_alpha / beta_v = ",delta_alpha)

    set_question("3- Becs et Volets: bec" )

    # je rentre l'angle comme 1 rad, le problème étant linéaire je trouve les 
    # coefficients devant beta_b
    DeltaCL,DeltaCmF,DeltaAlpha_0,DeltaAlpha_a=EffetVoletBec(Volet=False,theta=0,longueur=0.25,angle=180/np.pi)

def Exercice4_5():
    """
    Lois de squelette et d'épaisseur
    """
    set_title("Lois de squelette et d'épaisseur")
    print("Solution analytique et Matlab")

def Exercice4_6():
    """ 
    Exercice sur les becs et volets, Aile rectangulaire
    """
    set_title("Volets et becs : Aile rectangulaire")
    validation = False  # si on veut valider la fonction EffetVoletBec
    AileRectangle  = True   # exercice sur l'Aile Rectangle
    rho        = 1.2    # masse volumique en kg/m^3
    U0         = 50     # vitesse en m/s
    S          = 30     # surface de référence en m^2
    if validation : 
        EffetVoletBec(Volet=True,theta=53,longueur=0)
        EffetVoletBec(Volet=True,theta=0,longueur=0.2,angle=10)
        EffetVoletBec(Volet=False,theta=143)
        EffetVoletBec(Volet=False,longueur=0.2,angle=10)
  
    if AileRectangle:
        dCL,dCmF,dAlpha_0,dAlpha_a=EffetVoletBec(Volet=True,theta=0,longueur=0.2,angle=20)
        DeltaL = 1/2*rho*U0**2*S*dCL
         
        print("Delta L                           : %e N"%(DeltaL))
        dCL,dCmF,dAlpha_0,dAlpha_a=EffetVoletBec(Volet=False,longueur=0.05,angle=15)
        DeltaL = 1/2*rho*U0**2*S*dCL
        print("Delta L                           : %e N"%(DeltaL))

    # courbes pour la figure
    
    # volet :
    lv=np.linspace(0,0.3,301)
    lb=np.linspace(0,0.15,301)
    angle=1
    dCLv,dCmFv,dAlpha_0v,dAlpha_av=EffetVoletBec(Volet=True,theta=0,longueur=lv,angle=angle)
    dCLb,dCmFb,dAlpha_0b,dAlpha_ab=EffetVoletBec(Volet=False,theta=0,longueur=lb,angle=angle)
    c=100

   
    plt.figure(figsize=(10,8.5))
    plt.title(r"$C_L$ et $Cm_F$")
    plt.plot(c*lv,dCLv,label=r"volet:$\Delta C_L$")
    plt.plot(c*lb,dCLb,label=r"bec  :$\Delta C_L$")
    plt.plot(c*lv,dCmFv,label=r"volet:$\Delta Cm_F$")
    plt.plot(c*lb,dCmFb,label=r"bec  :$\Delta Cm_F$")
    plt.xlabel(r"$l/\ell$")
    plt.grid()
    #plt.xlim(-2.05,2.05)
    plt.legend(loc="best")
    
    plt.figure(figsize=(10,8.5))
    plt.title(r"$\alpha_0$ et $\alpha_a$")
    plt.plot(c*lv,np.rad2deg(dAlpha_0v),label=r"volet:$\Delta\alpha_0$")
    plt.plot(c*lb,np.rad2deg(dAlpha_0b),label=r"bec  :$\Delta\alpha_0$")
    plt.plot(c*lv,np.rad2deg(dAlpha_av),label=r"volet:$\Delta\alpha_a$")
    plt.plot(c*lb,np.rad2deg(dAlpha_ab),label=r"bec  :$\Delta\alpha_a$")
    plt.xlabel(r"$l/\ell$")
    plt.grid()
    #plt.xlim(-2.05,2.05)
    plt.legend(loc="best")

    plt.show()

    


def Exercice4_7():
    """
    Profil en théorie linéarisée : profil NACA
    """
    set_title("Profil NACA")
    cas             = "analytique" #'numérique'
    par             = set_parameters_airfoil()
    s               = LinearTheoryAirfoil(par)
    s.AirfoilParam = ['4412']
   

    nc          = 601
    if cas=="numérique":
    
        # test sur un profil NACA, 
        s.alpha     = 0
        s.airfoil   = ["naca4","numérique","y+-"]
        s.pnt       = nc
        
        set_question("1- Profil NACA "+s.AirfoilParam[0])

        # calcul du profil NACA
        # la corde est légèrement supérieure à 1 lorsque le profil est cambré
        # il faut donc le redimensionner et le définir entre  -1/2 <= x/c <= 1/2
        # cela crée une erreur numérique pas nécessairement négligeable
        # il faut mieux faire la solution analytique ci dessous.
       
        xNaca,yNaca,thetaNaca=naca(s.AirfoilParam[0],nc,half_cosine_spacing = True)
        chord = abs(np.amin(xNaca)-np.amax(xNaca))
        print("corde originale du profil NACA %s  : %f"%(s.AirfoilParam,chord))
        # correction sur la corde si elle est prise en compte (cela ne change pas grand chose): 
        xNaca,yNaca=(xNaca-np.amin(xNaca))/chord,yNaca/chord
        xNaca-=0.5
        indBA,indBF=np.argmin(xNaca),np.argmax(xNaca)
        print("x_BA = x[%i] = %f , \t x_BF = x[%i] = %f"%(indBA,xNaca[indBA],indBF,xNaca[indBF]))
        s.xe,s.ye=xNaca[indBF:indBA+1],yNaca[indBF:indBA+1]
        s.xi,s.yi=xNaca[indBA:],yNaca[indBA:]
        print('\t extrados : x = ',s.xe[:2],'...',s.xe[-3:])
        print('\t            y = ',s.ye[:2],'..', s.ye[-3:])
        print('\t intrados : x = ',s.xi[:2],'...',s.xi[-3:])
        print('\t            y = ',s.yi[:2],'...',s.yi[-3:])
    elif cas=="analytique":
        s.airfoil      = ["naca4","analytique","ys"]
        s.AirfoilParam = ['4412',]
        s.alpha        = 0
        s.nFourier     = 3 
         
    else:
        return

    
    set_question("2- Théorie linéaire")
    s.plot_Kp     = True
    s.plot_airfoil= True
    s.run_theory2D()
    #s.PlotDelta_s()

   
    par             = set_parameters_airfoil()
    s1               = LinearTheoryAirfoil(par)
     
    s1.airfoil      = ["naca4","analytique","ys"]
    s1.AirfoilParam = ['4412',]
    s1.alpha        = 0
    s1.nFourier     = 3

    s1.run_theory2D()

    par             = set_parameters_airfoil()
    s2               = LinearTheoryAirfoil(par)
    s2.airfoil      = ["naca4","analytique","ys"]
    s2.AirfoilParam = ['4412',]
    s2.alpha        = 0.514     # incidence d'adaptation.
    s2.nFourier     = 3
    s2.pnt          = nc
    s2.run_theory2D()


 

    
    
    
    fig = plt.figure()
    fig.suptitle('Comparaison des cambrures ', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$\eta$',fontsize=20)
    ax.grid()
    ax.plot(s.x,np.rad2deg(s.delta_s),'-',linewidth=1,color='black',label="numérique")   
    ax.plot(s1.x,np.rad2deg(s1.delta_s),'-',linewidth=1,color='red',label="analytique") 
    plt.legend(loc="best") 
    plt.show()   

    fig = plt.figure(figsize=(16,4))
    fig.suptitle('Comparaison des profils ', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$\eta$',fontsize=20)
    ax.grid()
    ax.plot(s.xe,s.ye,'-',linewidth=1,color='blue')
    ax.plot(s.xi,s.yi,'-',linewidth=1,color='blue')
    ax.plot(s.x,s.y_s,'-',linewidth=1,color='black')
    ax.plot(s1.x,s1.y_s,'-',linewidth=1,color='red')
    ax.axis("equal")


    fig = plt.figure()
    fig.suptitle('Comparaison des Kp', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$\eta$',fontsize=20)
    ax.grid()
    ax.plot(s.x[0:-1],s.Kp_s[0:-1],'-',linewidth=1,label=s.AirfoilParam[0]+"numérique")   
    ax.plot(s1.x[0:-1],s1.Kp_s[0:-1],'-',linewidth=1,label=s1.AirfoilParam[0]+"analytique")  
    ax.plot(s2.x[0:-1],s2.Kp_s[0:-1],'-',linewidth=1,label=s2.AirfoilParam[0]+"analytique")  
    plt.legend(loc="best") 
    plt.show()   



def Exercice4_8():
    """
    Profil en théorie linéarisée : profil NACA symmétrique
    """
    set_title("Profil NACA")
    
    par             = set_parameters_airfoil()
    s               = LinearTheoryAirfoil(par)
    s.AirfoilParam = ['2415',]
    nc          = 601
    s.alpha     = 5
    s.airfoil   = ["naca4","analytique","ys"]
    s.nFourier    = 3
    s.pnt         = nc
    s.plot_Kp     = True
    s.plot_airfoil= True
    s.run_theory2D()
