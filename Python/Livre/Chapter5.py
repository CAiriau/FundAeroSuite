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

from IncompressibleFlow.LiftingLine  import * 
from IncompressibleFlow.Prandtl  import * 



def Exercice5_1():
        """
        Aile optimale elliptique
        """
        set_title("Aile optimale elliptique")
        par=set_wing_parameters()
        s= WingAnalysis(par)
        s.plot_wing =True
        s.plot_circulation=False
        s.method    = 1
        s.S         = 22.48
        s.span      = 11.25
        s.lamb      = 0.
        # validated      
        s.run_analysis()

def Exercice5_2():
        """
        Vrillage d'une aile elliptique
        """
        set_title("Aile elliptique vrillée")
        set_question("1- cas 2" )
        par=set_wing_parameters()
        s= WingAnalysis(par)
        s.plot_wing = False
        s.plot_circulation=False
        s.method    = 1
        s.S         = 22.48
        s.span      = 11.25
        s.lamb      = 0.
        s.alphaV_law = 2 
        s.alphaV_par = np.deg2rad(1.)
        s.run_analysis()  

def Exercice5_3():
        """
        Braquage d'un volet d'une aile elliptique
        """
        set_title("Braquage d'un volet d'une aile elliptique")
        print("pas d'implémentation en python")

def Exercice5_4():
        """
        Calcul du foyer aérodynamique d'une aile
        """
        set_title("Calcul du foyer aérodynamique d'une aile")
        """
        Origine des moments au bord d'attaque sur l'axe de symétrie
        """

        CmF_profil=-0.05        # CmF du profil (on peut définir une loi en fonction de theta)
        Ny        = 101         # Nombre de points pour discrétiser l'aile
        Lambda    = np.deg2rad(5)# Angle pour définir la position du foyer du profil
        pente     = np.tan(Lambda)
        b         = 7            # demi envergure de l'aile
        x0        = -0.25         # position du foyer du profil sur l'axe de symétrie
        Ns        = 3            # nombre de mode de Fourier pour le calcul de la circulation de l'aile
        lamb      = 7            # allongement de l'aile
        alpha     =np.deg2rad(1) # incidence de l'aile alpah-alpha_0 constant
        #U0        = 50           # vitesse en m/s 
        #rho       = 1.3          # masse volumique
        xA        = -0.5         # moment au bord d'attaque de l'aile, sur l'axe de symétrie

        A=np.zeros(Ns)
        T=np.zeros(Ns)

        def x_F_profil(theta):
            """
            pour définir la position du foyer par profil
            """
            return x0+pente*b*abs(np.cos(theta))

        theta=np.linspace(0,np.pi,Ny) 

        # calcul avec l'aile elliptique uniquement
        S    = 4*b**2/lamb
        Lmoy = 2*b/lamb
        L0   = 4/np.pi*Lmoy
        x0   *= L0
        xA   *= L0

        A[0] = 2*alpha/(2+lamb)   # valeur pour l'aile elliptique
        print("b                    : ",b)
        print("lambda               : ",lamb)
        print("S                    : ",S)
        print("Lmoy                 : ",Lmoy)
        print("L0                   : ",L0)
        print("alpha (°)            : ",np.rad2deg(alpha))
        print("A1                   : ",A[0])
        print("Lambda (°)           : ",np.rad2deg(Lambda))
        print("xA/L0                : ",xA/L0)
        print("CmF_profil           : ",CmF_profil)
        print("x0/L                 : ",x0/L0)

         # ici soit je calcule A_1,A_3, ... soit je prends les valeurs de l'aile elliptique       
                                   # A_3, etc ... sont nuls
        # relations analytiques
        print("\nMETHODE ANALYTIQUE :\n")
        CmA = 32/(3*np.pi**2)*CmF_profil +(np.pi*(xA-x0)-4*b/3*pente)*lamb/Lmoy*A[0]
        CmA_T1=32/(3*np.pi**2)*CmF_profil
        xF  = x0+4*b/(3*np.pi)*pente

        print("xF/L0                : ",xF/L0)
        print("CmA/CmF_profil       : ",CmA/CmF_profil)
        print("CmA T1               : ",CmA_T1)
        
        print("\nMETHODE NUMERIQUE :\n")
        #  calcul par les intégrales

        def loi_corde(theta):
            return L0*np.sin(theta)

        def loi_xF(theta):
            return x0+b*abs(np.cos(theta))*pente

        def terme_profil_CmF(theta):
            return loi_corde(theta)**2*np.sin(theta)

        def terme_profil_xF(theta,n):
            return loi_xF(theta)*np.sin(theta)*np.sin(n*theta)

         
        T1,res1 = integrate.quad(terme_profil_CmF,0,np.pi)
        T2      = xA/Lmoy*np.pi*lamb*A[0]
        for n in range(Ns):
            T[n],res = integrate.quad(terme_profil_xF,0.,np.pi,args=(n+1,))
            #print("T[%i]              : %f"%(n,T[n]))
            T[n]*=-2*lamb/Lmoy*A[n]
            print("terme du point F     : ","n= ",n+1," : ",T[n],"résidu : ",res)

        T1*=CmF_profil/(2*Lmoy**2)
        print("terme du profil      : ",T1,"résidu : ",res1)
        print("terme du point A     : ",T2)
        print("CmA/CmF_profil       : ",(T1+T2+np.sum(T))/CmF_profil)



def Exercice5_5():
        """
        Sillage tourbillonnaire Airbus A380-800
        """
        set_title("Sillage tourbillonnaire Airbus A380-800")
        rho0    = 1.2       # masse volumique au niveau du sol [kg/m^3]
        rho     = 0.35      # masse volumique à l'altitude de vol [kg/m^3]
        Ud      = 80.       # vitesse de décollage [m/s]
        U       = 260.      # vitesse de croisière [m/s]
        S       = 840       # surface de l'aile [m^2]
        b       = 40        # demi envergure [m]
        W       = 4.e6      # Poids supporté par l'aile principale [N]

        set_question('1 - allongement, corde c_0')
        lamb    = 4 *b**2/S
        c_0     = 8*b/(np.pi*lamb)
        print("Allongement                              : %2.2f "%(lamb))
        print("Corde à l'emplanture                     : %2.2f m "%(c_0))

        set_question('2 - Cz au décollage ')
        print("w = ",W)
        Czdec=2*W/(rho0*Ud**2*S)
        print("Cz décollage                             : %2.3f "%(Czdec))
        alphadec = (1+2/lamb)*Czdec/(2*np.pi)
        print("Incidence equivalente                    : %2.3f "%(np.rad2deg(alphadec)))

        print("En vol de croisière : ")
        Cz=2*W/(rho*U**2*S)
        alpha = (1+2/lamb)*Cz/(2*np.pi)
        print("Cz                                       : %2.3f "%(Cz))
        print("Incidence equivalente                    : %2.3f "%(np.rad2deg(alpha)))

        set_question("3 - Circulation")
        
        def circulation(U,Cz):
            return (4*b*U*Cz)/(np.pi*lamb)
        Gamma0=circulation(Ud,Czdec)
        print("Circulation au décollage                 : %2.3f "%(Gamma0))
        Gamma=circulation(U,Cz)
        print("Circulation en croisière                 : %2.3f "%(Gamma))


        set_question("6 - Caractérisation d'un tourbillon marginal")
        y_T = np.pi*b/4
        R_T_sur_b =np.sqrt(2./3.-np.pi**2/16)
        print("y_T                                      : %2.3f m  "%(y_T))
        print("R_T / b                                  : %2.3f m  "%(R_T_sur_b))


        set_question("7 - vitesse de descente d'un tourbillon marginal")
        v_d=Gamma/(np.pi**2*b)
        print("Vitesse de descente d'un tourbillon      : %2.3f m/s "%(v_d))

        set_question("Cxi au décollage")
        Cxi     = Czdec**2/(np.pi*lamb)
        print("Cxi au décollage (sans sol) ,finesse     : %2.4f  %2.2f  "%(Cxi,Czdec/Cxi))
        def f_DeltaCxi(h):
            return Gamma0**2/(np.pi*2*Ud**2*S)*np.log(1+y_T**2/h**2)
        print("Delta Cxi, h/y_T                         : %2.4f %2.2f "%(f_DeltaCxi(y_T/2),1/2))
        print("Delta Cxi, h/y_T                         : %2.4f %2.2f "%(f_DeltaCxi(y_T/4),1/4))


def Exercice5_6():
        """
        Correction en soufflerie
        """
        set_title("Correction en soufflerie")
        rho0    = 1.2       # masse volumique au niveau du sol [kg/m^3]
        U0      = 30.      # vitesse dans la soufflerie [m/s]
        lamb    = 6         # allongement de l'aile
        b       = 0.5       # demi envergure [m]
        R       = 1.0       # rayon de la soufflerie

        Sv      = np.pi*R**2
        

        set_title("Aile optimale elliptique")
        par=set_wing_parameters()
        s= WingAnalysis(par)
        s.plot_wing =False
        s.plot_circulation=False
        s.method    = 1
        s.alpha     = np.deg2rad(4.)
        s.S         = 0.
        s.U_0       = U0
        s.span      = 2*b
        s.lamb      = lamb
        s.run_analysis()
        
        Sa=s.S
        CL=s.CL_elliptic
        CD=s.CD_elliptic
        x0=np.pi*b/4
        
        print("x0                           : %f"%(x0))
        print("R^2/x0                       : %f"%(R**2/x0))
        set_question("delta_Alpha_i")
        delta_Alpha_i = Sa/(8*Sv)*CL
        print("delta Alpha_i                : %f °"%(np.rad2deg(delta_Alpha_i)))
        print("delta Cx_i                   : %f"%(CL*delta_Alpha_i))
        print("delta Cx_i/Cx_i              : %f"%(CL*delta_Alpha_i/CD))
        print("Sa/Sv                        : %f"%(Sa/Sv))
        
def  Exercice5_7():
    """
    Tourbillon en fer à cheval
    """
    set_title("Tourbillon en fer à cheval")

    b=1
    Gamma= -np.pi
    eta=np.linspace(0.15,4,201)
    w=[]
    for x in eta:
        w.append(Vitesse_Tourbillon([x,0],0,-b,b,Gamma))

    plt.figure(figsize=(12, 8))
    plt.title(r'$\tilde w$', fontsize=fz, fontweight='bold')   
    plt.xlabel(r"$\eta=x/b$")
    plt.xlabel(r"$w/w_{ref}$")
    plt.plot(eta,w)
    plt.grid()
    plt.show()

    print("x/b = 4  : w/wref = ",Vitesse_Tourbillon([4*b,0],0,-b,b,Gamma))
    print("x/b = 0  : w/wref = ",Vitesse_Tourbillon([0,0],0,-b,b,Gamma))
          

def  Exercice5_8():
    """
    Vol en formation
    """
    set_title("Vol en formation")

    test=False
    b=5
    Gamma= 4*np.pi
    r=0.5
    d= r* b
    x=10

    if test:
        w= Vitesse_Tourbillon([x,0],0,-b,b,Gamma)
        print(" w = %f  "%(w))

    set_question("1- un oiseau isolé")
    w0=Vitesse_Tourbillon([0,0],0,-b,b,Gamma)
    print("vitesse au centre de l'aile         : ",w0)

    set_question("2- trois oiseaux volant de front")
    w_centre=Vitesse_Tourbillon([0,0],0,-3*b,3*b,Gamma)
    w_droite=Vitesse_Tourbillon([0,2*b],0,-3*b,3*b,Gamma)
    w_gauche=Vitesse_Tourbillon([0,-2*b],0,-3*b,3*b,Gamma)

    print("W_C/w0 = %f, W_D/w0 = %f W_G/w0 = %f "%(w_centre/w0,w_droite/w0,w_gauche/w0))

    set_question("3- trois oiseaux volant en chevron")
   
    r=np.linspace(0.4,0.43,21)

    r=[0.40,0.418768,0.42]
    Wc,Wd=[],[]
    for R in r:
        d=R*b
        # oiseau au centre
        w_oiseau={}
        w_centre=Vitesse_Tourbillon([0,0],0,-b,b,Gamma)
        w_droite=Vitesse_Tourbillon([0,0],d,b,3*b,Gamma)
        w_gauche=Vitesse_Tourbillon([0,0],d,-3*b,-b,Gamma)
        w_oiseau["centre"]=(w_centre+w_droite+w_gauche)/w0

        # oiseau  de droite
        w_centre=Vitesse_Tourbillon([d,2*b],0,-b,b,Gamma)
        w_droite=Vitesse_Tourbillon([d,2*b],d,b,3*b,Gamma)
        w_gauche=Vitesse_Tourbillon([d,2*b],d,-3*b,-b,Gamma)
        w_oiseau["droite"]=(w_centre+w_droite+w_gauche)/w0

        # oiseau  de gauche
        w_centre=Vitesse_Tourbillon([d,-2*b],0,-b,b,Gamma)
        w_droite=Vitesse_Tourbillon([d,-2*b],d,b,3*b,Gamma)
        w_gauche=Vitesse_Tourbillon([d,-2*b],d,-3*b,-b,Gamma)
        w_oiseau["gauche"]=(w_centre+w_droite+w_gauche)/w0

        wmean=(w_oiseau["droite"]+ w_oiseau["centre"]+ w_oiseau["gauche"] )/3
        Wc.append(w_oiseau["centre"])
        Wd.append(w_oiseau["droite"])
        
        print("%f %f %f %f %f "%(R,w_oiseau["gauche"],w_oiseau["centre"],w_oiseau["droite"], 45*wmean))

    print("Angle optimal  = %f °"%(np.rad2deg(np.arctan(1/0.418768))))
    # intersection = 0.418768
    plt.figure()
    plt.plot(r,Wc)
    plt.plot(r,Wd)
    plt.show()


def  Exercice5_9():
    """
    Tourbillons de sillage au voisinage du sol
    """
    set_title("Tourbillons de sillage au voisinage du sol")


    def fun(x,k):
        """
        trajectoires ou lignes de courant
        """
        return k*x/np.sqrt((k**2+1)*x**2-k**2)

    k=[0.5,1,2];
    npt=501
    a=3

    plt.figure(figsize=(12, 6))
    plt.title(r'Trajectoires', fontsize=fz, fontweight='bold')   
    for K in k:
        kmin=K/np.sqrt(K**2+1)
        print(" asymptote     : ",kmin)
        x=np.linspace(kmin*1.01,a,npt)
        y=fun(x,K)
        plt.plot(x[y<=a],y[y<=a], label=" k = %2.1f"%K)
        plt.plot(-x[y<=a],y[y<=a], label=" k = %2.1f"%K)
        plt.plot([-a,a],[kmin,kmin],'--')
        plt.plot([kmin,kmin],[0,a],'--')
        plt.plot([-kmin,-kmin],[0,a],'--')

    plt.xlim([-a,a])
    plt.axis('equal')
    plt.grid()
    plt.legend()
    plt.show()

def  Exercice5_10():
    """
    Résolution numérique du problème de Prandtl, cas général
    implémentation en utilisant la classe
    (sinon voir l'exercice 12)
    """

    set_title("Résolution numérique du problème de Prandtl, cas général")

    stop = 10       # numéro de question où il y a un arrêt prématuré
                   # pour aller jusqu'au bout mettre 10

    Ny          = 101               # nombre de points pour la méthode I
    plot_circulation = False
    gradient_Cz =True

    # pour avoir le gradient du Cz / alpha au lieu du Cz (tout est linéaire en alpha pour le Cz)
    if gradient_Cz:
        alpha_ref   = np.deg2rad(1)   # incidence de référence
    else:
        alpha_ref = np.deg2rad(4.1)

    lamb_ref    = 7                 # allongement de référence
    b_ref       = 7                 # demi-envergure de référence
    U0_ref      = 101               # vitesse de référence
    Ns          = 5                 # nombre de modes de Fourier
    k_ref       = np.pi             # demi gradient de portance
    plot_incidence = False           # effet d'incidence
    show_results=True               # affichage des valeurs CL, CD
    plot_allongement = True        # pour tester l'effet d'allongement

    set_question("1 - Cas de référence pour l'aile elliptique")

    Cz,Cxi,b,Gamma0=0.35,0.00557,7,45
    lamb = Cz**2/(np.pi*Cxi)
    A1   = Cz/(np.pi*lamb)
    alpha = A1*(2+lamb)/2
    U0 = Gamma0*Cz/(4*b*Cxi)
    U0bis = Gamma0/(4*b*A1)
    S = 4*b**2/lamb
    print("envergure                : %3.3f m"%(2*b))
    print("Allongement              : %3.3f "%(lamb))
    print("Surface de l'aile        : %3.3f m^2"%(S))
    print("vitesse                  : %3.3f m/s \t %3.3f m/s"%(U0, U0bis))
    print("indicence                : %3.3f °"%(np.rad2deg(alpha)))

             
    if stop==1 : return

   

    if plot_circulation:
        fig = plt.figure(10)
        plt.title(r'Loi de circulation', fontsize=fz, fontweight='bold') 
   

    set_question("2 - Méthode I : intégration en sin m theta")

    

    par=set_wing_parameters()
    s= WingAnalysis(par)
    s.plot_wing = False
    s.plot_circulation= False
    s.method    = 1      # sine approach
    s.S         = 0
    s.span      = 2*b_ref
    s.lamb      = lamb_ref
    s.alpha     = alpha_ref
    s.wing      ="elliptic"
    s.U_0       = U0_ref
    s.nFourier  = Ns
    s.ny        = Ny
    s.run_analysis()

    if stop==2 : return

    set_question("3 - Méthode II : Collocation")

    s.An        = 0
    s.S,s.l_0   = 0,0
    s.method    = 2 
    s.ny        = Ns+2     
    s.run_analysis()
  
    if stop==3 : return

    # =============================================================
    def wing_aerodynamics(lamb,b,alpha,U0,Angle=[np.deg2rad(0),np.deg2rad(0)],forme="rectangulaire",verif=False):
        """
        Calcul d'une aile rectangulaire, trapézoïdale ou delta
        """


        # **** Equations du livre : ***
        # c'est aussi dans la méthode n° 3 de la classe
        # On pourrait donc supprimer les 40 lignes ci dessous (jusqu'à la marque "possiblement effaçable") 
        def B_coef(n,m):
            return -np.float(8*(m*n)/( ((m+n)**2-1)*((m-n)**2-1)*np.pi))
        def C_coef(n,m):
            return np.float(4*(m**2+n**2-1)*(-1)**(int((m+n)/2))/( ((m+n)**2-1)*((m-n)**2-1)*np.pi))

        def calcul_An(mu0,mu1,alpha):
            """
            Matrice à inverser pour l'aile rectangulaire ou trapézoïdale
            Seulement pour les coefficients impairs
            """
            # construction de la matrice
            mat=np.zeros((3,3))
            mat[0,:] = [B_coef(1,1)+mu0+mu1*C_coef(1,1), B_coef(1,3)+3*mu1*C_coef(1,3), B_coef(1,5)+5*mu1*C_coef(1,5)]
            mat[1,:] = [B_coef(3,1)+mu1*C_coef(3,1), B_coef(3,3)+3*mu0+3*mu1*C_coef(3,3), B_coef(3,5)+5*mu1*C_coef(3,5)]
            mat[2,:] = [B_coef(5,1)+mu1*C_coef(5,1), B_coef(5,3)+3*mu1*C_coef(5,3), B_coef(5,5)+5*mu0+5*mu1*C_coef(5,5)]

            rhs=np.zeros((3))
            rhs[0]=alpha*(mu0+mu1*C_coef(1,1))
            rhs[1]=alpha*mu1*C_coef(1,3)
            rhs[2]=alpha*mu1*C_coef(1,5)
            print("mat  = ",mat)
            print("rhs  = ",rhs)
            An=np.linalg.solve(mat,rhs)
            return An


        if verif:
            print("Coefficients des matrices Bnm et Cnm :")
            for p in range(3):
                for q in range(3):
                    n,m=2*p+1,2*q+1
                    print("B[%i,%i] = %f, C[%i,%i] = %f %i"%(n,m,B_coef(n,m),n,m,C_coef(n,m),int((m+n)/2)))
        s.nFourier=Ns
        print("forme de l'aile   : ",forme)
       

        if forme=="elliptique":
            An=np.zeros(3)
            mu0= 2*k_ref/(np.pi*lamb)
            mu1=0
            An[0] =2/(lamb+2)*alpha
            s.wing="elliptic"
        else:
            mu0=k_ref/8*(4/lamb+np.sum(np.tan(Angle)))
            mu1= -k_ref/4*np.sum(np.tan(Angle))
            An=calcul_An(mu0,mu1,alpha)
            s.wing="tapered"
            

        Gamma_max=4*b*U0*(An[0]-An[1]+An[2])
        print("An     : ",An)
        print("lamb   : ",lamb)
        print("mu0    : ",mu0)
        print("mu1    : ",mu1)
        print("Gamma_max : ",Gamma_max)

        print("alpha  : ",np.rad2deg(alpha))
        CLa,CDa  = np.pi*lamb*An[0],np.pi*lamb*np.sum((2*np.arange(0,3)+1)*An**2)
        print("CL     : ",CLa)
        print("CD     : ",CDa)

        # *** possible effaçable ci-dessous car répétition de ce qu'il y a dans la classe.


        # METHOD I

        s.plot_wing = False
        s.WingAngle = Angle
        s.method    = 1      # sine approach
        s.S,s.l_0   = 0,0
        s.lamb,s.span  = lamb,2*b
        s.alpha,s.U_0   = alpha,U0
        s.nFourier  = Ns
        s.ny        = Ny
        s.run_analysis()
        
        CL,CD=s.CL,s.CD

        if plot_circulation: plt.plot(s.y,s.Gamma,label=forme)
           
        # METHOD II

        s.method    = 2      # collocation
        s.S,s.l_0   = 0,0
        s.lamb,s.span  = lamb,2*b
        s.alpha,s.U_0   = alpha,U0
        s.ny        = Ns+2
        s.run_analysis()

        # METHOD III

        if forme!="elliptique":
            s.method    = 3      
            s.S,s.l_0,s.b   = 0,0,0
            s.lamb,s.span  = lamb,2*b
            s.alpha,s.U_0   = alpha,U0
            s.WingAngle = Angle   
            s.nFourier  = Ns 
            s.ny        = Ny
        s.run_analysis()

        
        print("CL= %2.8f \t CL_a= %2.8f"%(CL,CLa))
        print("CD= %2.8f \t CD_a= %2.8f"%(CD,CDa))

        return CL,CD,CLa,CDa
    # =============================================================


    

    CL,CD,CLa,CDa=wing_aerodynamics(lamb_ref,b_ref,alpha_ref,U0_ref,Angle=[0,0],forme="elliptique")
    
    if stop==3 : return

    set_question("4 - Aile rectangulaire")

    lamb,U0,alpha=7,101,np.deg2rad(4.1)
    WingAngle= [0,0]
    CL,CD,CLa,CDa=wing_aerodynamics(lamb,b,alpha_ref,U0_ref,Angle=WingAngle,forme="rectangulaire")

    s.plot_wing = False
    s.WingAngle = WingAngle
    s.method    = 3      
    s.S,s.l_0   = 0,0
    s.lamb,s.span  = lamb,2*b
    s.alpha,s.U_0   = alpha,U0
    s.nFourier  = 5
    s.ny        = Ny
    s.run_analysis()
     
    if stop==4 : return   

    set_question("5 - Aile trapézoïdale")
    
    lamb,U0,alpha=7,101,np.deg2rad(4.1)
    WingAngle= [np.deg2rad(5),np.deg2rad(5)]
    CL,CD,CLa,CDa=wing_aerodynamics(lamb_ref,b_ref,alpha_ref,U0_ref,Angle=WingAngle,forme="trapézoïdale")
  
    if stop==5 : return

    set_question("6 - Aile delta")

    lamb,U0,alpha=7,101,np.deg2rad(4.1)
    WingAngle= [np.arctan(4/lamb),0]
    CL,CD,CLa,CDa=wing_aerodynamics(lamb_ref,b_ref,alpha_ref,U0_ref,Angle=WingAngle,forme="delta")
    
    if plot_circulation:
        plt.legend();plt.grid();plt.show()

    
        

    if stop==6 : return


    set_question("7 - Circulation")

    if stop==7 : return
    
    set_question("8 - effet d'incidence")  

   
    
    plot_circulation=False
    b=b_ref
    U0=U0_ref
    lamb=lamb_ref
    leg=["rectangulaire","trapézoïdale","delta","elliptique"]
    Na=28+1
    alpha_table=np.deg2rad(np.linspace(-4,10,Na))
    CL,CD=np.zeros((4,Na)),np.zeros((4,Na))
    CLa,CDa=np.zeros((4,Na)),np.zeros((4,Na))
    for i in range(Na):
        CL[0,i],CD[0,i],CLa[0,i],CDa[0,i]=wing_aerodynamics(lamb,b,alpha_table[i],U0,Angle=[0,0],forme="rectangulaire")
        CL[1,i],CD[1,i],CLa[1,i],CDa[1,i]=wing_aerodynamics(lamb,b,alpha_table[i],U0,Angle=[np.deg2rad(5),np.deg2rad(5)],forme="trapézoïdale")
        CL[2,i],CD[2,i],CLa[2,i],CDa[2,i]=wing_aerodynamics(lamb,b,alpha_table[i],U0,Angle=[np.arctan(4/lamb),0],forme="delta")
        CL[3,i],CD[3,i],CLa[3,i],CDa[3,i]=wing_aerodynamics(lamb,b,alpha_table[i],U0,Angle=[0,0],forme="elliptique")
        #CL[3,i],CD[3,i]=np.pi*lamb*2/(lamb+2)*alpha_table[i],np.pi*lamb*(2/(lamb+2)*alpha_table[i])**2
        #CLa[3,i],CDa[3,i]=CL[3,i],CD[3,i]

    if show_results:
        print("alpha = ",np.rad2deg(alpha_table))
        for i in range(4):
            print(leg[i],"\nCL= ",CL[i,:],"\nCD= ",CD[i,:])

    if plot_incidence:  
        
        fig = plt.figure(11)
        plt.title(r'$C_z$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(np.rad2deg(alpha_table),CL[i,:],label=leg[i])
            plt.plot(np.rad2deg(alpha_table),CLa[i,:],label=leg[i]+'_a')
        plt.ylabel(r"$C_z$")
        plt.xlabel(r"$\alpha (^\circ)$")
        plt.legend(loc='best'); plt.grid()    
        
        fig = plt.figure(12)
        plt.title(r'$C_x$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(np.rad2deg(alpha_table),CD[i,:],label=leg[i])
            plt.plot(np.rad2deg(alpha_table),CDa[i,:],label=leg[i]+'_a')
        plt.ylabel(r"$C_x$")
        plt.xlabel(r"$\alpha (^\circ)$")
        plt.legend(loc='best'); plt.grid()
        
        fig = plt.figure(13)
        plt.title(r'$C_x=f(C_z)$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(CL[i,:],CD[i,:],label=leg[i])
            plt.plot(CLa[i,:],CDa[i,:],label=leg[i]+'_a')
        plt.ylabel(r"$C_x$")
        plt.xlabel(r"$C_z$")
        plt.legend(loc='best'); plt.grid()

        fig = plt.figure(14)
        plt.title(r'$C_z/C_x$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(np.rad2deg(alpha_table),CL[i,:]/CD[i,:],label=leg[i])
            plt.plot(np.rad2deg(alpha_table),CLa[i,:]/CDa[i,:],label=leg[i]+'_a')
        plt.ylabel(r"$f$")
        plt.xlabel(r"$\alpha (^\circ)$")
        plt.legend(loc='best'); plt.grid()

        plt.show()
     
    if stop==8 : return
     
    set_question("9 - effet d'allongement")

    
    alpha=alpha_ref
    b=b_ref
    U0=U0_ref
    
    if plot_allongement:
        Na=21
        lamb_table=np.linspace(4,16,Na)
        CL,CD=np.zeros((4,Na)),np.zeros((4,Na))
        CLa,CDa=np.zeros((4,Na)),np.zeros((4,Na))
        for i in range(Na):
            CL[0,i],CD[0,i],CLa[0,i],CDa[0,i]=wing_aerodynamics(lamb_table[i],b,alpha,U0,Angle=[0,0],forme="rectangulaire")
            CL[1,i],CD[1,i],CLa[1,i],CDa[1,i]=wing_aerodynamics(lamb_table[i],b,alpha,U0,Angle=[np.deg2rad(1),np.deg2rad(1)],forme="trapézoïdale")
            CL[2,i],CD[2,i],CLa[2,i],CDa[2,i]=wing_aerodynamics(lamb_table[i],b,alpha,U0,Angle=[np.arctan(4/lamb_table[i]),0],forme="delta")
            CL[3,i],CD[3,i],CLa[3,i],CDa[3,i]=wing_aerodynamics(lamb_table[i],b,alpha,U0,Angle=[0,0],forme="elliptique")
            #CL[3,i],CD[3,i]=np.pi*lamb_table[i]*2/(lamb_table[i]+2)*alpha,np.pi*lamb_table[i]*(2/(lamb_table[i]+2)*alpha)**2
            #CLa[3,i],CDa[3,i]=CL[3,i],CD[3,i]
        
        if gradient_Cz:
            grad_Cz_elliptique=(2*lamb_table*np.pi)/(2+lamb_table)*np.pi/180
        else:
            grad_Cz_elliptique=1

        leg=["rectangulaire","trapézoïdale","delta","elliptique"]
        
        fig = plt.figure(15)
        plt.title(r'$C_z$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(lamb_table,CL[i,:]/grad_Cz_elliptique,label=leg[i])
            plt.plot(lamb_table,CLa[i,:]/grad_Cz_elliptique,label=leg[i]+'_a')
        plt.ylabel(r"$C_z$")
        plt.xlabel(r"$\lambda$")
        plt.legend(loc='best'); plt.grid()
       
        fig = plt.figure(16)
        plt.title(r'$C_x$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(lamb_table,CD[i,:],label=leg[i])
            plt.plot(lamb_table,CDa[i,:],label=leg[i]+'_a')
        plt.ylabel(r"$C_x$")
        plt.xlabel(r"$\lambda$")
        plt.legend(loc='best'); plt.grid()
         

        plt.show()

    print('Fin exercice ...')


def  Exercice5_11():
    """
    Résolution numérique du problème de Prandtl avec l'effet de sol.
    """

    set_title("Résolution numérique du problème de Prandtl, aile elliptique avec effet de sol")     
    
    gradient_Cz      = True

    # pour avoir le gradient du Cz / alpha au lieu du Cz (tout est linéaire en alpha pour le Cz)
    if gradient_Cz:
        alpha   = np.deg2rad(1)     # incidence de référence
    else:
        alpha   = np.deg2rad(8)

    lamb        = 7                 # allongement de référence
    b           = 7                 # demi-envergure de référence
    U0          = 50                # vitesse de référence
    Ns          = 3                 # nombre de modes de Fourier
    k           = np.pi             # demi gradient de portance
    variation_h = True              # effet d'altitude
    variation_lambda = True         # pour tester l'effet d'allongement
    M           = 51                # nombre de point sur l'aile en envergure
    k           = np.pi             # demi gradient de portance
    
    # 
    theta = np.linspace(0,np.pi,M)  # distribution angulaire
    y     = -b*np.cos(theta)        # coordonnées suivant l'envergure
    l0    = 8*b/(np.pi*lamb)        # corde moyenne
    S     = np.pi*l0*b/2            # Surface de l'aile

    print('Ns = %i, M=%i, alpha= %2.0f °, lambda= %f'%(Ns,M,np.rad2deg(alpha),lamb))
    print('lambda = %f, 2b = %f, l0= %f, S = %f m^2'%(lamb,2*b,l0,S))


    set_question("I - SANS EFFET DE SOL")
    
    A=Prandtl_E_OGE(Ns,lamb,b,l0,k,y,theta,alpha)
    CL=np.pi*lamb*A[0]             
    Cxi=np.pi*lamb*np.sum(np.arange(1,Ns+1)*A**2)
    A1_ref=alpha/(1+np.pi*lamb/(2*k))
    CL_ref=np.pi*lamb*A1_ref
    Cxi_ref=np.pi*lamb*A1_ref**2
    print('Calcul   : A_1 = %f, CL = %f, Cxi = %f '%(A[0],CL,Cxi))
    print('Reference: A_1 = %f, CL = %f, Cxi = %f '%(A1_ref,CL_ref,Cxi_ref))
    print("A1_bis = ",2*np.pi*lamb*alpha/(2+lamb))

    set_question("II - AVEC EFFET DE SOL, Lambda fixé")

    if variation_h:
        Nh        = 101
        h_sur_b   = np.linspace(0.0001,4.0001,Nh) 
        delta_CL  = np.zeros((Nh), dtype=float)
        delta_CDi = np.zeros((Nh), dtype=float)

        for i,h in enumerate(h_sur_b):
            Anew=Prandtl_E_IGE(Ns,lamb,b,l0,k,h,y,theta,alpha)
            #print('A = ',Anew)
            CL=np.pi*lamb*Anew[0]              # premier coefficient de la série
            Cxi=np.pi*lamb*np.sum(np.arange(1,Ns+1)*Anew**2)
            delta_CL[i]=(CL-CL_ref)/CL_ref
            delta_CDi[i]=(Cxi-Cxi_ref)/Cxi_ref
            if i <= 25 :
                print('h=%f,  A_1 = %2.6f, CL/CL_ref-1 = %f, Cxi/Cxi_ref-1 = %f'%
                    (h,Anew[0],CL/CL_ref-1,Cxi/Cxi_ref-1))
            
         
        # plot
        fig = plt.figure(1,figsize=(10,8))
        plt.title(r"Effet de sol, $\lambda$ = %f"%lamb, fontsize=fz, fontweight='bold')
        plt.xlabel(r'$h/b$', fontsize=16)
        plt.plot(h_sur_b,delta_CL, color='b', linestyle='-', linewidth=2,label=r'$\Delta C_L/C_L $')
        plt.plot(h_sur_b,delta_CDi, color='r', linestyle='-', linewidth=2,label=r'$\Delta C_{x_i}/C_{x_i} $')
        xmax,ymax=4,0.7
        plt.axis([0, xmax, 0, ymax])                                   
        plt.grid(which='both');plt.legend(loc="best")                                          
        

    set_question("II - AVEC EFFET DE SOL, Lambda variable, h fixé")
   
    if variation_lambda:
        Nl=101
        lamb_table = np.linspace(4,24,Nl)  
        delta_CL   = np.zeros((Nl), dtype=float)
        delta_CDi  = np.zeros((Nl), dtype=float)
        h          = 0.25  # h/b

        for i,lamb in enumerate(lamb_table):
            A1_ref  = alpha/(1+np.pi*lamb/(2*k))
            CL_ref  = np.pi*lamb*A1_ref
            Cxi_ref = np.pi*lamb*A1_ref**2
            l0   = 8*b/(np.pi*lamb) 
            Anew = Prandtl_E_IGE(Ns,lamb,b,l0,k,h,y,theta,alpha)
            CL   = np.pi*lamb*Anew[0]              # premier coefficient de la série
            Cxi  = np.pi*lamb*np.sum(np.arange(1,Ns+1)*Anew**2)
            delta_CL[i] = (CL-CL_ref)/CL_ref
            delta_CDi[i]= (Cxi-Cxi_ref)/Cxi_ref
            if i<= 25:
                print('lambda=%f,  A_1 = %2.6f, CL/CL_ref-1 = %f, Cxi/Cxi_ref-1 = %f'%
                    (lamb,Anew[0],CL/CL_ref-1,Cxi/Cxi_ref-1))
        
        # plot
        fig = plt.figure(2,figsize=(10,8))
        plt.title(r"Effet de sol, $h/b$ = %f"%(h), fontsize=fz, fontweight='bold')
        plt.xlabel(r'$\lambda$', fontsize=16)
        plt.plot(lamb_table,delta_CL, color='b', linestyle='-', linewidth=2,label=r'$\Delta C_L/C_L (\%)$')
        plt.plot(lamb_table,delta_CDi, color='r', linestyle='-', linewidth=2,label=r'$\Delta C_{x_i}/C_{x_i} (\%)$')
        plt.axis([4, 24, 0,0.4])    
        plt.xticks(np.arange(4,25,4))                               
        plt.grid();plt.legend(loc="best")                                          
        

    if variation_h or variation_lambda: plt.show()

 

def  Exercice5_12():
    """
    Aile elliptique du Spitfire
    """
    set_title("Aile d'un Spitfire")
    npt         = 101       # nombre de points suivant y
    b           = 18*0.3028 # demi envergure (pieds vers mètres)
    Lambda      = 5.5       # allongement de l'aile

    L0          = 8*b/(np.pi*Lambda) # corde à l'emplanture 
    S           = np.pi*L0*b/2       # surface de l'aile
    x0          = 0                  # ligne de position des foyers de chaque profil
    fleche      = 30                  # en degrés
    k           = np.tan(np.deg2rad(fleche)) # si on veut mettre un angle de flèche sur l'aile
     

    print("envergure                : %3.3f m"%(2*b))
    print("Allongement              : %3.3f "%(Lambda))
    print("Surface de l'aile        : %3.3f m^2"%(S))
    print("Corde à l'emplanture     : %3.3f m"%(L0))
    theta=np.linspace(0,np.pi/2,npt)
    y=b*np.cos(theta)
    x0=-k*y/b
    corde=L0*np.sqrt(1-(y/b)**2)
    xLE,xTE=x0+corde/4,x0-3*corde/4
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$y (m)$', fontsize=16)
    ax.set_ylabel(r'$x (m)$',fontsize=20)
    ax.set_title("Demi aile d'un Spitfire")
    ax.plot(y,xLE, color='b', linestyle='-', linewidth=2)
    ax.plot(y,xTE, color='r', linestyle='-', linewidth=2)
    ax.axis("equal")                              
    plt.show()

    

def Exercice5_13():
    """ 
    CL, CD : aile elliptique sans effet de sol
    """
    set_title("aile elliptique SANS effet de sol")

    print('Aile elliptique pour le moment')
     

    b=np.float64(1.0)           # demi-envergure
    lamb=np.float64(10.0)       # allongement
    alpha=np.deg2rad(np.float64(1)) # incidence en degrés
    N=2                         # nombre de termes de la série de Fourier
    M=21                        # nombre de point sur l'aile en envergure
    h=np.float(64)              # hauteur/demi-envergure ?
    Uinf=np.float(1.)           # vitesse infinie amont
    k=np.pi                     # demi gradient de portance
    theta=np.linspace(0,np.pi,M)# distribution angulaire
    y=-b*np.cos(theta)          # coordonnées suivant l'envergure
    l0=8*b/(np.pi*lamb)         # corde moyenne
    S=np.pi*l0*b/2               # Surface de l'aile
    

    print('N = %i, M=%i, alpha= %2.0f °, lambda= %f'%(N,M,np.rad2deg(alpha),lamb))
    print('lambda = %f, 2b = %f, l0= %f, S = %f m^2'%(lamb,2*b,l0,S))

       
    A=Prandtl_E_OGE(N,lamb,b,l0,k,y,theta,alpha)

    CL=np.pi*lamb*A[0]              # premier coefficient de la série
    Cxi=np.pi*lamb*np.sum(np.arange(1,N+1)*A**2)

    A1_ref=alpha/(1+np.pi*lamb/(2*k))
    CL_ref=np.pi*lamb*A1_ref
    Cxi_ref=np.pi*lamb*A1_ref**2

    print('Calcul   : A_1 = %f, CL = %f, Cxi = %f '%(A[0],CL,Cxi))
    print('Reference: A_1 = %f, CL = %f, Cxi = %f '%(A1_ref,CL_ref,Cxi_ref))
    
    set_title("aile elliptique AVEC effet de sol")

    Nh=31
    #h_sur_b=np.linspace(0.01,5.01,Nh) 
    h_sur_b=np.linspace(0.01,5.01,Nh) 
    err_CL = np.zeros((Nh), dtype=float)
    err_CDi = np.zeros((Nh), dtype=float)
   
    i=0

     

    for h in list(h_sur_b):
        Anew=Prandtl_E_IGE(N,lamb,b,l0,k,h,y,theta,alpha)
        CL=np.pi*lamb*Anew[0]              # premier coefficient de la série
        Cxi=np.pi*lamb*np.sum(np.arange(1,N+1)*Anew**2)
        err_CL[i]=CL/CL_ref-1
        err_CDi[i]=Cxi/Cxi_ref-1
        print('h=%f,  A_1 = %2.6f, CL/CL_ref-1 = %f, Cxi/Cxi_ref-1 = %f'%
            (h,Anew[0],CL/CL_ref-1,Cxi/Cxi_ref-1))
        i=i+1
     
     

    # plot
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$2h/b=h/L$', fontsize=16)
    #ax.set_ylabel(r'$Kp$',fontsize=20)
    ax.set_title('Rapport / à la référence en milieu infini')
    ax.plot(2*h_sur_b,err_CL, color='b', linestyle='-', linewidth=2,label=r'$e(C_L)$')
    ax.plot(2*h_sur_b,err_CDi, color='r', linestyle='-', linewidth=2,label=r'$e(C_{x_i})$')
    xmax,ymax=5,0.5
    ax.axis([0, xmax, 0, ymax])                                   
    xmajor_ticks,xminor_ticks = np.arange(0, xmax, 1),np.arange(0, xmax, 0.25)                                              
    ymajor_ticks,yminor_ticks = np.arange(0,ymax, 0.2),np.arange(0,ymax, 0.05)                                              
    ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True)  
    ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
    ax.grid(which='both');ax.legend()                                          
    plt.show()

   

def func(x,a,b,c):
    """
    fonction test
    """
    return a+b*x+c*x**2

def integralO2(a,b,c):
    """
    méthode des trapèzes
    """
    m=100
    x=np.linspace(0,1,m)
    dx=x[1]-x[0]
    for i in range(len(x)):
        val=func(x,a,b,c)
        s=dx*(np.sum(val[1:m-2])+0.5*(val[0]+val[m-1]))
    return s

def Exercice5_15():
    """ 
    test de l'intégrale quad avec des paramètres
    comparaison avec les trapèzes
    """
    A,B=1,2
    C=np.linspace(0,2,2)
    for c in list(C):
        G,res=integrate.quad(func,0,1.,args=(A,B,c,))
        G1=integralO2(A,B,c)
        print(G,G1,A+B/2+c/3)


def Exercice5_14():
    """
    Chargement d'une aile de P51 Mustang
    """
    set_title("Chargement d'une aile de P51 Mustang")
    f2m    = 0.3048
    l_0,l_t= 8.48*f2m,3.87*f2m
    b      = 5.64
    Lambda = np.rad2deg(np.arctan((l_0-l_t)/b))
    lamb   = 5.876
    U0     = 123.

    Ny     = 101
    alpha  = np.deg2rad(3.)
    Ns     = 5

    
    par=set_wing_parameters()
    s= WingAnalysis(par)
    s.plot_wing = True
    s.plot_circulation= True
    s.method    = 1      # sine approach
    s.S         = 0
    s.span      = 2*b
    s.lamb      = lamb
    s.alpha     = alpha
    s.wing      ="tapered"
    s.U_0       = U0
    s.nFourier  = Ns
    s.ny        = Ny
    s.WingAngle =[np.deg2rad(Lambda/2),np.deg2rad(Lambda/2)] 
    s.run_analysis()
    
