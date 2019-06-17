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
from Tools.misc                 import *
from scipy.optimize             import newton
from Livre.ConformalMapping.mappings  import *
from Livre.panneaux_intro  import *


twopi  = 2*np.pi
FigDim = [10,10]
fz    = 30



def Exercice3_1():
    """
    Cylindre
    """
    set_title("Coefficient de pression sur un cylindre")
    theta=np.linspace(0,360,361)
    set_question('1- cylindre, 2- sphère')
    plt.figure(figsize=(FigDim))
    plt.title(r'Kp pour le cylindre et la sphère', fontsize=fz, fontweight='bold')      
    plt.plot(theta,1-4*(np.sin(np.deg2rad(theta)))**2,linewidth=2, label='cylindre')
    plt.plot(theta,1-9/4*(np.sin(np.deg2rad(theta)))**2,linewidth=2,label='sphère')
    plt.xlabel(r"$\theta (^\circ)$",fontsize=fz);plt.ylabel(r'$Kp$',fontsize=fz)
    plt.legend(loc="best")
    plt.grid()
    plt.show()

def plot_Kp(xe,xi,Kpe,Kpi,cas):
    """
    Dessin d'une liste de Kp
    """
    plt.figure(figsize=(15,10))
    plt.title(r"$K_p$")
    for k in range(len(xe)):
        plt.plot(xe[k],Kpe[k],label="Kpe: "+cas[k])
        plt.plot(xi[k],Kpi[k],label="Kpi: "+cas[k])
    plt.grid()
    plt.xlim(-0.05,1.05)
    plt.legend(loc="best")
    plt.show()


def Exercice3_2():
    """
    Transformation de Joukowski
    """
    set_title("Transformation de Joukowski")
    xe,xi=[],[]
    Kpe,Kpi=[],[]
    cas=["Ellipse","Profil épais","Squelette","cambré"]
    set_question("1 - Ellipse")
    par=set_parameters()
    s= AirfoilConfMap(par)
    # profil de Joukowski, le profil est toujours dessiné sous incidence nulle
    # dans ce cas la paramètre est Zc
    s.map           = "Joukowski1"
    s.airfoil       = "Ellipse"
    s.alpha         = 0
    s.a             = 0.9
    s.plot_circle   = False
    s.plot_airfoil  = False
    s.opt_levels    = "manuel"
    s.levels_iso    = [-2,2,0.2]
    s.adimChord     = True
    s.plot_velocity = False
    s.plot_Kp       = True
   

    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    

    set_question("2 - Profil épais")

    par=set_parameters()
    s= AirfoilConfMap(par)
    # profil de Joukowski, le profil est toujours dessiné sous incidence nulle
    # dans ce cas la paramètre est Zc
    s.map           = "Joukowski1"
    s.a             = 1  
    s.R             = 0
    s.alpha         = 0
    s.Zc            = -0.1+1j*0
    s.plot_circle   = True
    s.plot_airfoil  = False
    s.opt_levels    = "manuel"
    s.levels_iso    = [-1,1,0.2]
    s.adimChord     = True
    s.plot_velocity = False
    s.plot_Kp       = True
    s.run_airfoil()

    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    
    set_question("3 - Profil squelettique")

    s.a             = 1  
    s.R             = 0
    s.Zc            = 0+1j*0.1     
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    
    set_question("4 - Profil cambré")

    s.a             = 1  
    s.R             = 0
    s.Zc            = -0.1+1j*0.1
    s.camberline    = True     
    s.plot_airfoil  = True

    s.run_airfoil()

    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    plot_Kp(xe,xi,Kpe,Kpi,cas)
     
                
def Exercice3_3():           
    """
    Transformation de Karman-Trefftz
    """
    set_title("Transformation de Karman-Trefftz")
    xe,xi=[],[]
    Kpe,Kpi=[],[]
    cas=["sym5","cambré5","sym10","cambré10"]
    set_question("3 - Profil symétrique")
    par=set_parameters()
    s= AirfoilConfMap(par)
    eps=[0.05,0.1]
    for epsi in eps:
        s.map           = "Karman-Trefftz"
        s.a             = 1  
        s.R             = 0
        s.alpha         = 0 
        s.beta          = 0
        s.Zc            =-epsi+1j*0
        s.camberline    = False    
        s.plot_airfoil  = False
        s.opt_levels    = "manuel"
        s.levels_iso    = [-2.,2.,0.1]
        s.k             = 1.9  
        s.plot_velocity = False 
        s.adimChord     = True
        s.plot_Kp       = True
        s.run_airfoil()

        xe.append(s.xe)
        xi.append(s.xi)
        Kpe.append(s.Kpe)
        Kpi.append(s.Kpi)

        set_question("3 - Profil cambré")
        s.Zc            =-epsi+1j*epsi
        s.camberline    = False   
        s.plot_airfoil  = True
        s.a             = 1  
        s.R             = 0
        s.run_airfoil()
        xe.append(s.xe)
        xi.append(s.xi)
        Kpe.append(s.Kpe)
        Kpi.append(s.Kpi)

    plot_Kp(xe,xi,Kpe,Kpi,cas)

def Exercice3_4():            
    """
    Transformation de Von Mises
    """
    set_title("Transformation de Von Mises")
    xe,xi=[],[]
    Kpe,Kpi=[],[]
    cas=["Zc=0;5","cambré5","cambré10"]
    set_question("9 - Profil eps=0.05")
    par=set_parameters()
    s= AirfoilConfMap(par)

   
    # Profil de von Mises symétrique
    # dans ce cas la paramètre est eps (epsilon) et k  
    # mu est calculé à partir de Zc
    # le profil est celui du livre d'exercice
    # le code ne fonctionne que si a=2 et mu_0+mu_1+a=0
    s.a             = 1  
    s.R             = 0
    s.map           = "von Mises"
    s.alpha         = 0.
    s.k             = 2
    s.eps           = 0.05
    s.Zc            = 0
    #s.Zc            = -s.eps*(1-1j)
    #s.camberline    = False    
    #s.plot_airfoil  = False
    s.opt_levels    = "manuel"
    s.levels_iso    = [-2.,2.,0.1] 
    s.plot_velocity = False 
    s.adimChord     = True
    s.plot_Kp       = True
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - Profil eps=0.05")
    s.a             = 1  
    s.R             = 0
    s.eps           = 0.05
    s.Zc            = -s.eps*(1-1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - Profil eps=0.10")
    s.a             = 1  
    s.R             = 0
    s.eps           = 0.10
    s.Zc            = -s.eps*(1-1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    plot_Kp(xe,xi,Kpe,Kpi,cas)


def Exercice3_5():           
    """
    Transformation de Van de Vooren et Jong
    """
    set_title("Transformation de Van de Vooren et Jong")
    xe,xi=[],[]
    Kpe,Kpi=[],[]
    cas=["sym,0.05","sym,0.10","cambré5","cambré10"]

    set_question("9 - Profil eps=0.05")
    par=set_parameters()
    s= AirfoilConfMap(par)

    s.map           = "van de Vooren"
    s.a             = 1  
    s.R             = 0
    s.alpha         = 0.
    s.k             = 2
    s.eps           = 0.05
    s.Zc            = 0+1j*0    # -s.eps*(1-1j)
    #s.camberline    = False    
    #s.plot_airfoil  = False
    s.opt_levels    = "manuel"
    s.levels_iso    = [-2.,2.,0.1] 
    s.plot_velocity = False 
    s.adimChord     = True
    s.plot_Kp       = True
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - Profil eps=0.10")
    s.a             = 1  
    s.R             = 0
    s.alpha         = 0
    s.eps           = 0.10
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - Profil eps=0.05")
    s.a             = 1  
    s.R             = 0
    s.alpha         = 0
    s.eps           = 0.05
    s.k             = 1.9
    s.Zc            = -s.eps*(1-1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - Profil eps=0.10")
    s.a             = 1  
    s.R             = 0
    s.alpha         = 0
    s.eps           = 0.10
    s.k             = 1.9
    s.Zc            = -s.eps*(1-1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    plot_Kp(xe,xi,Kpe,Kpi,cas)


def Exercice3_6():         
    """
    Profil à double courbure
    paramètres : a et b
     180° > argument de b > 90° 

    PROBLEME: suivant les valeurs des paramètres les lignes
    de courant peuvent traverser le profil !!! : problème de l'arctangente.
    """
    set_title("Profil à double courbure")

    set_alert("Problèmes possibles sur les lignes de courant")
    set_info("attention")
    cas = 'd'

    set_question("9 - Profil cas %s"%(cas))
    par=set_parameters()
    s= AirfoilConfMap(par)        
    s.map           = "double pointe"
    s.alpha         = 0.
    s.a             = 0  # il est calculé en fonction du rayon
    
    if cas =='c':
        """
        Cas à problème
        """
        beta            = 46
        s.b             = s.R/2*np.exp(1j*np.deg2rad(beta))
        s.Zc            = 0+1j*0

           #  case 'P2'
           #      b = a/2;
           #      Z_c = 0;

           #  case 'P3'
           #      beta = 20*pi/180;
           #      epsilon = 0.2;
           #      b = sqrt(epsilon) * a * exp(i*beta)
           #      Z_c = 0;
        

    elif cas=="d" : 
        """
        Cas à problème
        """
        s.eps           = 0.2 
        beta            = 20
        s.b             = np.sqrt(s.eps)*s.R*np.exp(1j*np.deg2rad(beta))
        s.Zc            = 0+1j*0
         
    elif cas=="b":
       
        beta            = 135
        s.b             = s.R/10*np.exp(1j*np.deg2rad(beta))
        s.Zc            = s.b  # -s.eps*(1-1j)

    elif cas=="a":
       
        beta            = 180
        s.b             = s.R/10*np.exp(1j*np.deg2rad(beta))
        s.Zc            = s.b  # -s.eps*(1-1j)
        

    s.run_airfoil()



def Exercice3_7():
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
    s.plot_velocity = True
    s.plot_Kp       = True
    s.run_airfoil()


 

def Exercice3_8():
    """
    Introduction à la méthode des panneaux
    """
    #******************************************
    #  MAIN
    #******************************************
    # Paramètres :
    U0              = 1                 #  vitesses infinie amont
    plot_panneaux   = False
    plot_gamma      = False
    plot_vitesse    = False             # pour dessiner la vitesse
    plot_psi        = False             # pour dessiner les iso psi
    plot_2D         = True              # pour dessiner les iso vitesses
    alpha           = np.deg2rad(5)     # incidence en degrés , doit être non nul
    nP              = 101                 # nombre de panneaux sur la plaque initiale
    corde           = 1.                # longueur de la plaque
    nv              = 5001               # nombre de points pour la vitesse locale (plot 1D)
    y_wall          = 0.01

    s=plaque_plane(nP,y_wall,nv=nv,alpha=alpha)
    dessiner_contours(s[0],s[2],alpha)
    nP1,nP2=21,3
    s1=plaque_plane(nP1,y_wall,nv=nv,alpha=alpha)
    s2=plaque_plane(nP2,y_wall,nv=nv,alpha=alpha)

    if plot_2D:
        option="uv"  # sinon mettre "uv" ou "Kp"
        vel=["v","u"]
        plt.figure(figsize=(12,12))
        if option=="Kp" :
             plt.plot(s[3].real,1-abs(s[4]),'-',label="U, nP = %i "%(nP))
             plt.plot(s1[3].real,1-abs(s1[4]),'-',label="U, nP = %i "%(nP1))
             plt.plot(s2[3].real,1-abs(s2[4]),'-',label="U, nP = %i "%(nP2))
        else:
            if "u" in vel:
                plt.plot(s[3].real,s[4].real,'-',label="u, nP = %i "%(nP))
                plt.plot(s1[3].real,s1[4].real,'-',label="u, nP = %i "%(nP1))
                plt.plot(s2[3].real,s2[4].real,'-',label="u, nP = %i "%(nP2))
                eps=0.01
                x=np.linspace(eps,1-eps,101)
                plt.plot(x,reference_plaque_plane(x,alpha),'r--',label="Analytique")
            if "v" in vel:
                plt.plot(s[3].real,-s[4].imag,'-',label="v, nP = %i "%(nP))
                plt.plot(s1[3].real,-s1[4].imag,'-',label="v, nP = %i "%(nP1))
                plt.plot(s2[3].real,-s2[4].imag,'-',label="v, nP = %i "%(nP2))
        plt.ylim(-1,2)
        plt.grid()
        plt.legend(loc="best")


        plt.figure(figsize=(10,8))
        plt.plot(s[0].real,s[2].real*len(s[2])/s[5],'--s',label=r"$\Gamma$, nP = %i "%(nP))
        plt.plot(s1[0].real,s1[2].real*len(s1[2])/s1[5],'-o',Markersize=10,label=r"$\Gamma$, nP = %i "%(nP1))
        plt.plot(s2[0].real,s2[2].real*len(s2[2])/s2[5],'-o',Markersize=10,label=r"$\Gamma$, nP = %i "%(nP2))
        plt.ylim(0,6)
        plt.grid()
        plt.legend(loc="best")

    

    s3=plaque_plane(nP,-y_wall,nv=nv,alpha=alpha)
    plt.figure(figsize=(12,12))
    plt.title(r"$\Delta u = u^+-u^-$ at $|y|$ = %f"%(y_wall))
    plt.plot(s[3].real,s[4].real-s3[4].real,'b-',label=r"$\Delta U$, nP = %i "%(nP))

    plt.grid()
    plt.legend(loc="best")
    
    if plot_2D or plot_vitesse or plot_psi or plot_panneaux :
        plt.show()

def rapport_circulation(eta,alpha,exacte=False):
    """
    Rapport des circulations
    exacte = False    : méthode approchée ou exacte
    """
    if exacte:
        return (1+eta**2/16*np.cos(alpha)**2)/(1+eta/4*np.sin(alpha))
    else:
        return 1+eta**2/16

def rapport_portance(eta,alpha,exacte=False):
    """
    Rapport des portances
    exacte = False    : méthode approchée ou exacte
    """
    K=rapport_circulation(eta,alpha,exacte)
    if exacte:
        return K*(1-K*eta*np.sin(alpha)/(4+eta*np.sin(alpha)))
    else:
        return K*(1-eta/4*np.sin(alpha)*K)

def Exercice3_9():
    """ 
    Effet de sol
    """
    set_title("Effet de sol")

    # courbe des maximuns
    al=np.linspace(2,15,27)
    eta_ref=np.linspace(0,20,1001)
    r_appro_max,eta_appro_max=[],[]
    r_exact_max,eta_exact_max,=[],[]
    for k in range(len(al)):
        r0=rapport_portance(eta_ref,np.deg2rad(al[k]),False)
        k0=np.argmax(r0)
        r_appro_max.append(np.amax(r0))
        eta_appro_max.append(eta_ref[k0])
        r0=rapport_portance(eta_ref,np.deg2rad(al[k]),True)
        k0=np.argmax(r0)
        r_exact_max.append(np.amax(r0))
        eta_exact_max.append(eta_ref[k0])
         
        

    option="eta"  # ou "1/eta"
    eta_lim=[17,17,9,11,9,12,12,14]
    alpha  =[3,3,5,5,10,10,15,15]
    eta=[]
    r_appro=[]
    r_exact=[]
    for k in range(len(eta_lim)):
        eta0=np.linspace(0.1,eta_lim[k],101)
        alpha0=np.deg2rad(alpha[k])
        if k%2==0:
            exact=False
        else:
            exact=True
        r_appro.append(rapport_portance(eta0,alpha0,exact))
        r_exact.append(rapport_portance(eta0,alpha0,exact))
        eta.append(eta0)

    fig = plt.figure(figsize=(10,8))
    fig.suptitle(r'Effet de sol, $C_L/C_{L_\infty}$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    
    ax.grid()
    for k in range(len(alpha)):
        if k%2==0:
            leg=r"appro., $\alpha=$ %2.2f"
        else:
            leg=r"exact,  $\alpha=$ %2.2f"
        if option=="eta":
            ax.plot(eta[k],r_appro[k],'-',linewidth=2,label=leg%(alpha[k]))
            ax.plot(eta_appro_max,r_appro_max,'k--')
            ax.plot(eta_exact_max,r_exact_max,'r--')
        else:
            ax.plot(1/eta[k],r_appro[k],'-',linewidth=2,label=leg%(alpha[k]))
            ax.plot([1/x for x in eta_appro_max],r_appro_max,'k--')
            ax.plot([1/x for x in eta_exact_max],r_exact_max,'r--')
    if option=="eta":
        plt.xlim(0,10)
        ax.set_xlabel(r'$\eta=\ell/h$',fontsize=20)
    else:
        plt.xlim(0,1)
        ax.set_xlabel(r'$\eta=h/\ell$',fontsize=20)
    plt.ylim(0,3)
    plt.legend(loc="best")  
    plt.show()  


def Exercice3_10():
    """
    Ailes en tandem
    """
    pass


def Exercice3_11():
    """
    Tansformation de Schwarz et Christoffel: expansion dans un canal, 
    """
    pass



def corps_Rankine(y,m):
    """
     dessin du corps de Rankine
     prendre des précautions à cause de y = 0
    """
    x2=[];
    for Y in y:
        if abs(Y) > 1e-6:
            x2.append(2*Y/np.tan(2*Y/m)+1-Y**2)
        else:
            x2.append(m+1-Y**2)
    # test pour éviter des soucis
    for i,X2 in enumerate(x2):
        if X2 < 0 : 
            # print("x^2 = %e, i = %i"%(X2,i))
            x2[i]=0
    return np.sqrt(x2)
 

def plot_func(eta,m):
    # dessin de la fonction pour la recherche du zéro    
    
    fig = plt.figure()
    fig.suptitle('recherche de h/a pour  m = %1.2f'%m, fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$eta$',fontsize=20)
    ax.grid()
    ax.plot(eta,np.tan(eta/m),'-',linewidth=2,color='black') 
    ax.plot(eta,1/eta,'-',linewidth=2,color='red')    
    plt.show()  

def Kp_Rankine(x,y,m):
    """
    Cacul du Kp sur le corps de Rankine
    """
    Kp=[]
    U=[]
    T=[]
    for X,Y in zip(x,y):
        z=complex(X,Y)
        w=1-m/(z**2-1)
        T.append(1+m-abs(w))
        U.append(abs(w))
        Kp.append(1-abs(w)**2)
    return Kp,U,T

def solve_Rankine(m,eta_init=1,display=False):
    """
    Résoud le problème pour une valeur du paramètre
    """
    # eta = h/a
    # ell = L/a
    ell = np.sqrt(1+m)
    def func(eta):
        return np.tan(eta/m)-1/eta
    #pour une fois j'utilise la méthode de Newton intrinsèque à scipy   
    eta=newton(func,eta_init)
    Uaxe=1+m/(1+eta**2)
    # recherche du max U, min Kp    
    y=np.linspace(0,eta,1000); x=corps_Rankine(y,m)
    Kp1,U1,T1=Kp_Rankine(x,y,m)
    Umax=np.amax(U1);Imax=np.argmax(U1)
    if display: print('Umax/U0 = %f, Kpmin = %f , x( %i ) = %f'%(Umax,Kp1[Imax],Imax,x[Imax]))    
    return ell,eta,Uaxe,Umax,Kp1[Imax],x[Imax]

def Exercice3_12():
    """ 
    Ovale de Rankine
    """
    set_title("Ovale de Rankine")

    # eta = h/a
    # ell = L/a
    #   
    #  ENTREE UTILISATEUR       
    #
    cas           = 12                # pour traiter un cas particulier, numéro du casn, cf m_table
    n             = 201              # nombre de points pour les courbes
    eps           = 0.4               # troncature du domaine si option=2
    plot_ovale    = True              # pour dessine la forme du corps
    plot_Kp       = True              # pour dessiner le Kp
    plot_U        = True              # pour dessiner U/U0
    plot_test     = False             # pour dessiner la fonction 1+m-(u^2+v^2)          
    overview      = True              # si option=1, calculer pour un grand nombre de valeurs de m
                                      # si faux, on peut un zoom au voisinage d'une valeur.
    option        = 4                # option principale :
    k_compare     = [1,15]            # indice des cas qu'on dessine sur une même figure pour option=4
    #   = 1 : tableau des valeurs
    #   = 2 : on trouve graphiquement la valeur de eta (pour m < 0.1)
    #   = 3 : tracer du corps et le Kp pour une valeur de m  particulière donnée par la variable "cas"
    #   =4  : tracés pour le livre.
    # pour les petites valeurs de m il faut donner quasiment la solution pour avoir convergence
    # vers la bonne valeur
        
    
    if overview :
         # tableau des valeurs de m
         m_table       = np.array([0.02,0.1,0.15,0.1589,0.2,0.3,0.4,0.5,0.6,0.639,0.65,0.70,0.8,0.88,0.99,1,2,10,20,50])
         # valeur initiale de eta avant convergence, pour chaque valeur de m         
         eta_table_init= np.array([0.0308,0.143,0.17,0.2,0.263,0.3,0.4,0.5,0.5,0.5,0.5,0.5,0.65,1,1,1,1,1,1,1])
    else:
        m_table=np.linspace(0.87,0.88,21)
        eta_table_init=np.ones(len(m_table))
    
    Nplot    = len(m_table)          
    m_cas    = m_table[cas] 
    eta_init = eta_table_init[cas]
    
    ar=[]          # aspect ratio, épaisserr relative
    eta_val=[]     # valeur de eta en fonction de m
    
    if option==1:
        # tableau des valeurs
        print(" m  \t\t L/a \t h/a \t L/h \t\t h/L \t U(x=0)  Kp_axe \t Umax \t Kpmin \t x_min ")
        for m,eta_init in zip(m_table,eta_table_init):
            ell,eta,Uaxe,Umax,Kpmin,x_max=solve_Rankine(m,eta_init)
            ar.append(eta/ell)
            eta_val.append(eta) 
            print('{:3.4f}   \t {:1.3f} \t {:1.3f} \t  {:2.3f} \t {:2.3f} \t {:1.3f} \t {:1.3f} \t {:1.3f} \t {:1.3f}  {:1.5f}'.format(m,ell,eta,ell/eta,eta/ell,Uaxe, 1-Uaxe**2,Umax,Kpmin,x_max))
 
        # plot the aspect ratio versus m
         
        fig = plt.figure(10)
        fig.suptitle('épaisseur relative', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$m$',fontsize=20)
        ax.set_ylabel(r'$h/\ell$',fontsize=20)
        ax.grid()
        #ax.axis([1.0, 12, -0.05, 0.35])
        ax.plot(m_table[:Nplot],ar[:Nplot],'-',linewidth=2,color='black',label=r"$h/\ell$") 
        ax.plot(m_table[:Nplot],eta_val[:Nplot],'-',linewidth=2,color='red',label=r"$\eta$") 
        ax.legend(loc='lower right') 
        plt.show()
       
    elif option==2:
        print('cas traité : m = %f '%(m_cas))
        eta=np.linspace(eps,1.5*m_cas,n)       
        #eta=m_cas*np.cos(np.linspace(eps,np.pi/2-eps,n))
        print('eta = ',eta)
        plot_func(eta,m_cas)

    elif option==3:
        # cas particulier
        ell,eta,Uaxe,Umax,Kpmin,x_max=solve_Rankine(m_cas,eta_init)
        print("m = %f, L/a = %f, h/a = %f, L/h = %f, h/L = %f"%(m_cas,ell,eta,ell/eta,eta/ell))
        print("U(x=0)/U0 = %f, Kp(x=0) = %f"%(Uaxe,1-Uaxe**2))
        y=np.linspace(-eta,eta,n)
        x=corps_Rankine(y,m_cas)
        Kp1,U1,T1=Kp_Rankine(x,y,m_cas)
        Kp2,U2,T2=Kp_Rankine(-x,y,m_cas)

        Um=np.amax(U1)
        Im=np.argmax(U1)
        
        print('Umax/U0 = %f, Kpmin = %f , x( %i ) = %f'%(Um,Kp1[Im],Im,x[Im]))
        if plot_ovale:
            fig = plt.figure(1)
            fig.suptitle('Ovale de rankine pour m = %1.2f'%m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$',fontsize=20)
            ax.set_ylabel(r'$y/a$',fontsize=20)
            ax.grid()
            #ax.axis([1.0, 12, -0.05, 0.35])
            ax.plot(x,y,'-',linewidth=2,color='black') 
            ax.plot(-x,y,'-',linewidth=2,color='black')
            ax.axis('equal') 
        if plot_Kp:
            fig = plt.figure(2)
            fig.suptitle('Kp: Ovale de rankine pour m = %1.2f'%m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$',fontsize=20)
            ax.set_ylabel(r'$Kp$',fontsize=20)
            ax.grid()
            #ax.axis([1.0, 12, -0.05, 0.35])
            ax.plot(x,Kp1,'-',linewidth=2,color='black') 
            ax.plot(-x,Kp2,'-',linewidth=2,color='red') 
        if plot_U:
            fig = plt.figure(3)
            fig.suptitle('U: Ovale de rankine pour m = %1.2f'%m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$',fontsize=20)
            ax.set_ylabel(r'$U/U_0$',fontsize=20)
            ax.grid()
            #ax.axis([1.0, 12, -0.05, 0.35])
            # Pour comparer à Katz et Plotkin p. 62
            #ax.plot(x,[U1**2 for U1 in U1],'-',linewidth=2,color='black') 
            #ax.plot(-x,[U2**2 for U2 in U2],'-',linewidth=2,color='red') 
            ax.plot(x,U1,'-',linewidth=2,color='black') 
            ax.plot(-x,U2,'-',linewidth=2,color='red') 
        if plot_test:
            fig = plt.figure(4)
            fig.suptitle(r'$1+m-z^2$: Ovale de rankine pour m = %1.2f'%m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$',fontsize=20)
            ax.set_ylabel(r'$T$',fontsize=20)
            ax.grid()
            #ax.axis([1.0, 12, -0.05, 0.35])
            ax.plot(x,T1,'-',linewidth=2,color='black') 
            ax.plot(-x,T2,'-',linewidth=2,color='red') 
        
    
        if plot_ovale or plot_Kp or plot_U or plot_test : plt.show()  
   
    elif option==4:
        fig = plt.figure(6)
        fig.suptitle("Forme de l'ovale de rankine", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x/a$',fontsize=20)
        ax.set_ylabel(r'$y/a$',fontsize=20)
        ax.grid()

        
        fig1 = plt.figure(7)
        fig1.suptitle('U: Ovale de rankine', fontsize=14, fontweight='bold')
        ax1 = fig1.add_subplot(111)
        fig1.subplots_adjust(top=0.80)
        ax1.set_xlabel(r'$x/a$',fontsize=20)
        ax1.set_ylabel(r'$U/U_0$',fontsize=20)
        ax1.grid()
        
        j=1
        p=[]
        for k in k_compare:
            m=m_table[k]
            ell,eta,Uaxe,Umax,Kpmin,x_max=solve_Rankine(m,eta_table_init[k])
            print("m = %f, L/a = %f, h/a = %f, L/h = %f, h/L = %f"%(m,ell,eta,ell/eta,eta/ell))
            print("U(x=0)/U0 = %f, Kp(x=0) = %f"%(Uaxe,1-Uaxe**2))
            y=np.linspace(-eta,eta,n)
            x=corps_Rankine(y,m)
            Kp1,U1,T1=Kp_Rankine(x,y,m)
            Kp2,U2,T2=Kp_Rankine(-x,y,m)
            p.append([x,y,U1])
            
            ax.plot(x,y,'-',linewidth=2) 
            ax.plot(-x,y,'-',linewidth=2)            
            #ax.axis([1.0, 12, -0.05, 0.35])
            ax.axis('equal') 
            
            ax1.plot(x,U1,'-',linewidth=2) 
            ax1.plot(-x,U2,'-',linewidth=2) 
        plt.show()
    return
      
 
 
def Exercice3_20():
    """ 
    Exercices sur les transformations conformes
    """
    set_title("Transformation conforme")

    par=set_parameters()
    s= AirfoilConfMap(par)

    # choix du cas : 
    cas=6
    if cas == 0:
        # profil de Joukowski, le profil est toujours dessiné sous incidence nulle
        # dans ce cas la paramètre est Zc
        s.map           = "Joukowski1"
        s.alpha         =  0
        s.Zc            =  -0.0895386257+1j*  0.0523359562  # -0.15+1j*0.1 #   
        s.opt_levels    = "manuel"
        s.levels_iso    = [-1,1,0.2]
        s.adimChord     = True
        s.n_circle      = 101
        s.plot_velocity = True
        s.plot_Kp       = True
    elif cas == 1:
        # profil de Joukowski, le profil bouge en fonction de l'incidence
        # la transformation conforme contient en terme en alpha
        # dans ce cas les paramètres sont lamb et beta
        s.map           = "Joukowski"
        s.alpha         = 0
        s.beta          = 3
        s.lamb          = 1.1
        s.n_circle      = 201
        s.plot_velocity = True
        s.plot_psi      = False  
        s.plot_Kp       = True
    elif cas == 2:
        # profil de Joukowski, le profil est toujours dessiné sous incidence nulle
        # dans ce cas la paramètre est Zc
        s.map           = "Karman-Trefftz"
        s.alpha         = 3 #-5.7391704773
        s.Zc            =-0.05+1j*0.05
        s.opt_levels    = "manuel"
        s.levels_iso    = [-2.,2.,0.1]
        s.k             = 1.9   # exposant de la transformation  
        s.plot_velocity = True  
        s.adimChord     = False

    elif cas == 3:
        # Profil de von Mises non symétrique
        # dans ce cas la paramètre est eps (epsilon) et k  
        # mu est calculé à partir de Zc
        # le profil est celui du livre d'exercice
        # le code ne fonctionne que si a=2 et mu_0+mu_1+a=0
        s.map           = "von Mises"
        s.alpha         = 0.
        s.k             = 1.9
        s.eps           = 0.05
        s.Zc            = -s.eps*(1-1j)
        s.levels_iso    = [-1,1,0.1]
         
    elif cas == 4:
        # Profil de von Mises symétrique
        # dans ce cas la paramètre est eps (epsilon) et k  
        # mu est calculé à partir de Zc
        # le profil est celui du livre d'exercice
        # le code ne fonctionne que si a=2 et mu_0+mu_1+a=0
        s.map           = "von Mises"
        s.alpha         = 0.
        s.k             = 2.0
        s.eps           = 0.1
        s.Zc            = 0
        s.levels_iso    = [-1,1,0.1]

    elif cas == 5:
        # Profil de van de Vooren
        # dans ce cas la paramètre est eps (epsilon) et k  
         
        s.map           = "van de Vooren"
        s.alpha         = 0.
        s.k             = 2
        s.eps           = 0.05
        s.Zc            = 0+1j*0    # -s.eps*(1-1j)
        s.levels_iso    = [-1,1,0.1]
        

    elif cas==6:
        # profil à double pointes,
        # paramètres : a et b
        s.n_circle      = 101
        s.map           = "double pointe"
        s.alpha         = 0.
        s.b             = s.R/2*np.exp(1j*np.pi/3)
        s.Zc            = 0+1j*0
        #s.eps           = 0.2 
        #beta            = 110
        #s.b             =  np.sqrt(s.eps)*s.R*np.exp(1j*np.deg2rad(beta))
        s.Psi_method    = "Psi"
        s.camberline    = True

       #  case 'P2'
       #      b = a/2;
       #      Z_c = 0;

       #  case 'P3'
       #      beta = 20*pi/180;
       #      epsilon = 0.2;
       #      b = sqrt(epsilon) * a * exp(i*beta)
       #      Z_c = 0;
            

    elif cas==7:
        # profil à double pointes,
        # paramètres : a et b
        # il y a un problème avec le bord d'attaque.

        s.map           = "double pointe"
        s.alpha         = 0.
        s.eps           = 0.2 
        beta            = np.rad2deg(np.pi/2+0.1)
        s.b             = np.sqrt(s.eps)*s.R*np.exp(1j*np.deg2rad(beta))
        s.Zc            = 0+1j*0
        s.levels_iso    = [-1,1,0.1]

    elif cas==8:
        # profil à double pointes,
        # paramètres : a et b
        # 180° > argument de b > 90° 
        s.camberline    = False 
        s.map           = "double pointe"
        s.alpha         = 0.
        s.eps           = 0.05 
        beta            = 135
        s.b             = np.sqrt(s.eps)*s.R*np.exp(1j*np.deg2rad(beta))
        s.Zc            = s.b  # -s.eps*(1-1j)
        s.levels_iso    = [-1,1,0.1]
        
       # case 'P1a'      
       #      beta = 180*pi/180;
       #      b = 0.1*exp(i*beta);
       #      Z_c = b;
            
       #  case 'P1b'  
       #      beta = 135*pi/180;
       #      b = 0.1*exp(i*beta);
       #      Z_c = b;
            

           

    s.run_airfoil()
    

