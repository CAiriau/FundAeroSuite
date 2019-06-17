#!/bin/python
"""
    Correction des exercices du chapitre 12
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


import numpy as np
from scipy import interpolate,integrate
import matplotlib.pyplot as plt
from copy import deepcopy
from CompressibleFlow.fonctions import *
from CompressibleFlow.tables    import *
from Tools.misc                 import * 
 
 
# Variables globales

Largeur,Hauteur=10,8
Hequal=4


def ParoiBasse(x,type_paroi):
    """
    Equation de la paroi inférieure
    type de paroi :
    1 : pente constante 
    2 : y=0
    3 : avec une équation
    """
    def f(opt,x,L):
        """
        Equation de la paroi basse
        """
        if opt==0:
            return (np.cos(np.pi*x/L)-1)/2
        elif opt==1:
            return -np.pi/(2*L) *np.sin(np.pi*x/L)
   
    y=0
    if type_paroi == 1:     # divergent, pente constante
        pente=-np.tan(5/180*np.pi)
        y=pente*x
    elif type_paroi == 2:   # convergent, pente constante
        y=0.*x
    else:
        L=10;
        lim=10;
        val=f(0,lim,L)
        pente=f(1,lim,L)
        if x <= lim:
            y=f(0,x,L)
        else:
            y= val+pente*(x-lim);
    return y

def ParoiHaute(x,type_paroi):
    """
    Equation de la paroi supérieure
    """
    def f1(opt,x,L):
        """
        Equation de la paroi haute
        """
        if opt==0:
            return (3-np.cos(np.pi*x/L))/2
        elif opt==1:
            return np.pi/(2*L) *np.sin(np.pi*x/L)

    y=0
    if type_paroi == 1: 
        pente=np.tan(5/180*np.pi)
        y=1+pente*x
    elif type_paroi == 2: 
        pente=np.tan(-1/180*np.pi)
        y=1+pente*x
    elif type_paroi == 3: 
        lim=6
        if x <= lim:
            y=1-0.05*np.sin(np.pi*x/6)**2
        else:
            y=1
    else:
        L=10
        lim=10
        val=f1(0,lim,L)
        pente=f1(1,lim,L)
        if x <= lim :
            y=f1(0,x,L);
        else:
            y= val+pente*(x-lim)   
    return y




def Exercice12_1():
    """ 
    Canal divergent
    Calcul de l'ecoulement dans une conduite bidimensionnelle par methode 
    des caracteristiques pour un gaz parfait compressible isentropique
    Parametres : 
        Mach0   Nombre de Mach en entree de la conduite
        N       Nombre de caracteristiques ascendantes/descendantes
        Ltot    Longueur de la conduite 
        paroi_basse
        paroi_haute
    """

    #****************************************************************
    #  PARAMETRES et DONNEES
    #****************************************************************
    cas=0

    test_omega_mu = False
    show_geometry = False
    show_Mach     = True
    show_Median   = False
    show_label    = False

    Mach_ref      = [2,1.8,1.3,2]
    pbasse        = [1,2,2,4]
    phaute        = [1,2,3,4]
    Mlim          = [[2.0,3.2],[1.55,1.8],[1.1,1.4],[2,3.5]]
    nlev          =[21,21,21,16]

    
    Nparoi        = 101     # nombre d'intervalles pour discrétiser les parois
    N             = 7      # nombre de caractéristiques à tirer, doit être impair (7 ou 41)
    Ltot          = 10.     # longueur totale pour le calcul
    Ls            = 10.      # longueur montrée sur le dessin
    nlevels       = nlev[cas]      # nombre de niveau pour les iso Mach


    type_paroi_basse,type_paroi_haute=pbasse[cas],phaute[cas]
    Mach0=Mach_ref[cas]     # Mach d'entrée
    Mmin,Mmax  =Mlim[cas][0],Mlim[cas][1]     # échelle des Machs pour les figures

    
    gam=gamma
    epsilon = 1e-6             # pas pour le calcul des dérivées
    Nmedian=int((N-1)/2)       # démarre sur l'axe pour la valeur médiane
    IterMax=10                 # nombre d'itérations maximales dans Newton
    Precision=1e-14            # précision pour Newton
    # print('Nmedian = ',Nmedian)
  
    #****************************************************************
    #  INITIALISATIONS
    #****************************************************************
    

    set_title("Méthodes des caractéristiques dans un canal divergent") 
    
    if test_omega_mu :
        set_question('Test des fonctions omega et mu')
        mu=np.arcsin(1./Mach0)
        omega=omega_mu(mu)
        mu1=mu_omega(omega)
        print('M0= ',Mach0,'mu = ', mu,'Omega = ',omega,'mu = ',mu1)
        return
    
    # Déclaration des tableaux initialisés à 0.
    XP                       = np.zeros((N),dtype=float)
    lambdaPlusP,lambdaMoinsP = np.zeros((N),dtype=float),np.zeros((N),dtype=float)
    thetaP                   = np.zeros((N),dtype=float)

    XI,YI                    = np.zeros((N-1),dtype=float),np.zeros((N-1),dtype=float)
    lambdaPlusI,lambdaMoinsI = np.zeros((N-1),dtype=float),np.zeros((N-1),dtype=float)
    omegaI,thetaI            = np.zeros((N-1),dtype=float),np.zeros((N-1),dtype=float)
    muI                      = np.zeros((N-1),dtype=float)


    # Déclaration des listes 
    MachP=[]                                # Mach sur l'axe
    XcP,YcP=[],[]                                  # position de ce Mach sur l'axe
    Xmesh,Ymesh,mumesh=[],[],[]
    Xbas,Xhaut,mubas,muhaut=[],[],[],[]

    # Initialisation : points principaux initiaux
    Ymin0,Ymax0 = ParoiBasse(0.,type_paroi_basse),ParoiHaute(0.,type_paroi_haute)
    YP=Ymin0+np.arange(N)*(Ymax0-Ymin0)/float(N-1)
    yMedian=(Ymin0+Ymax0)/2
    print('YP(0) = ',YP,'Ymedian = ',yMedian,'Nmedian = ',Nmedian,'Y(Nmedian) = ',YP[Nmedian])
     
    # print('XP = ',XP)
    mu0 = np.arcsin(1/Mach0)
    omega0 = omega_mu(mu0)

    # section d'entrée :
    muP = mu0*np.ones((N),dtype=float)
    omegaP = omega0*np.ones((N),dtype=float) 
    lambdaPlusP  = omegaP-thetaP
    lambdaMoinsP = omegaP+thetaP

    # sur l'axe :
    MachP.append(Mach0)
    XcP.append(XP[0])
    YcP.append(YP[Nmedian])

    
    # Preparation des tableaux pour presentation des resultats
    #Xmesh,Ymesh = XP,YP
    #mumesh = muP
    Xmesh.append(deepcopy(XP))
    Ymesh.append(deepcopy(YP))
    mumesh.append(deepcopy(muP))

    Xbas.append(XP[0])
    Xhaut.append(XP[-1]) 
    mubas.append(muP[0])
    muhaut.append(muP[-1])
    
    #****************************************************************
    #  FIGURE GENERALE, PAROIS
    #****************************************************************
    # pour le dessin de la paroi
    XL=np.arange(Nparoi)*Ltot/float(Nparoi-1)
    print("pas en x : ",XL[1]-XL[0])
    YBAS = [ParoiBasse(eta,type_paroi_basse) for eta in XL]
    YHAUT= [ParoiHaute(eta,type_paroi_haute) for eta in XL]

    fig = plt.figure(0,figsize=(10,Hequal))
    fig.suptitle('Tuyère',fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$y$',fontsize=14)
    ax.set_xlabel(r'$x$',fontsize=14)
    ax.plot(XL,YBAS,'k-',linewidth=1,label='bas')
    ax.plot(XL,YHAUT,'k-',linewidth=1,label='haut')
    #ax.axis([0,3,-0.5,2])
    ax.axis('equal')  
    ax.set_xlim(0, Ls)
    #ax.grid('')
    if show_geometry:         
        plt.show()
        return

       
    # Fonction interne
    def test_position(X,Y):
        """
        test des positions des points intermédiaires
        """
        m=len(Y)
        I1,I2=[],[]
        for k in range(m):
            if Y[k] < ParoiBasse(X[k],type_paroi_basse): 
                I1.append(k)
            elif Y[k] > ParoiHaute(X[k],type_paroi_haute): 
                I2.append(k)
        if I1 :
            print('Problème sur la paroi basse  pour les indices',I1)
        if I2 :
            print('Problème sur la paroi haute  pour les indices',I2)

    #****************************************************************
    #  BOUCLE POUR L'AVANCEMENT EN x
    #****************************************************************
    # 
    Xmin=0              # démarrage en x=0

    while Xmin < Ltot:
        
        # --------------------------------------------
        # Contruction des points intermédiaires
        # --------------------------------------------
        for j in range(N-1):
            XI[j] = (YP[j]-YP[j+1]-np.tan(thetaP[j]+muP[j])*XP[j]
                  + np.tan(thetaP[j+1]-muP[j+1])*XP[j+1])/(np.tan(thetaP[j+1]-muP[j+1])-np.tan(thetaP[j]+muP[j]))
            YI[j] = YP[j]+np.tan(thetaP[j]+muP[j])*(XI[j]-XP[j])
            lambdaPlusI[j]  = lambdaPlusP[j]
            lambdaMoinsI[j] = lambdaMoinsP[j+1]
            omegaI[j] = (lambdaPlusI[j]+lambdaMoinsI[j])/2
            thetaI[j] = (lambdaMoinsI[j]-lambdaPlusI[j])/2
            muI[j] = mu_omega(omegaI[j])
    
            # tracé des segments de caractéristiques       
            ax.plot([XP[j],XI[j]],[YP[j],YI[j]],'r-')
            ax.plot([XP[j+1],XI[j]],[YP[j+1],YI[j]],'b-')
    
        test_position(XI,YI)
        
        # ------------------------------------------------------ 
        # Construction des points principaux a l'étape suivante
        # points centraux
        # ------------------------------------------------------ 
        for j in np.arange(1,N-1):
            XP[j] = (YI[j-1]-YI[j]-np.tan(thetaI[j-1]+muI[j-1])*XI[j-1]
                  + np.tan(thetaI[j]-muI[j])*XI[j])/(np.tan(thetaI[j]-muI[j])-np.tan(thetaI[j-1]+muI[j-1]))
            YP[j] = YI[j-1]+np.tan(thetaI[j-1]+muI[j-1])*(XP[j]-XI[j-1])
            lambdaPlusP[j]  = lambdaPlusI[j-1]
            lambdaMoinsP[j] = lambdaMoinsI[j]
            omegaP[j] = (lambdaPlusP[j]+lambdaMoinsP[j])/2
            thetaP[j] = (lambdaMoinsP[j]-lambdaPlusP[j])/2
            muP[j] = mu_omega(omegaP[j])
            lambdaPlusP[j]  = omegaP[j]-thetaP[j]
            lambdaMoinsP[j] = omegaP[j]+thetaP[j]
            # tracé des segments de caractéristiques
            ax.plot([XP[j],XI[j]],[YP[j],YI[j]],'b-')       # C- intermédiaires
            ax.plot([XP[j],XI[j-1]],[YP[j],YI[j-1]],'r-')    # C+ intermidiaires
        
        # ------------------------------------------------------ 
        # premier point sur la paroi basse (méthode de Newton)
        # ------------------------------------------------------ 
        x = XI[0]
        dx = 1
        compteur=0
        while (np.abs(dx)>Precision) and (compteur<IterMax) :
            compteur+=1
            FF = YI[0]+np.tan(thetaI[0]-muI[0])*(x-XI[0])-ParoiBasse(x,type_paroi_basse)
            dFFdx = np.tan(thetaI[0]-muI[0]) - (ParoiBasse(x+epsilon,type_paroi_basse)-ParoiBasse(x-epsilon,type_paroi_basse))/(2*epsilon)
            dx=FF/dFFdx
            x -= dx
        if compteur == IterMax : print('ATTENTION : non-convergence de Newton pour point P1')
        #print("compteur (i=0) = %i, \t dx = %e"%(compteur, dx))
        XP[0] = x
        YP[0] = ParoiBasse(XP[0],type_paroi_basse)
        thetaP[0] = np.arctan((ParoiBasse(XP[0]+epsilon,type_paroi_basse)-ParoiBasse(XP[0]-epsilon,type_paroi_basse))/(2*epsilon))
        omegaP[0] = lambdaMoinsI[0]-thetaP[0]
        muP[0] = mu_omega(omegaP[0])
        lambdaPlusP[0] = omegaP[0]-thetaP[0]
        ax.plot([XP[0],XI[0]],[YP[0],YI[0]],'b')   # fin des C- 

        # ------------------------------------------------------ 
        # dernier point sur la paroi haute (methode de Newton)
        # ------------------------------------------------------ 
        NI=N-2
        x = XI[-1]
        dx = 1
        compteur=0
        while (np.abs(dx)>Precision) and (compteur<IterMax) :
            compteur+=1
            FF = YI[NI]+np.tan(thetaI[NI]+muI[NI])*(x-XI[NI])-ParoiHaute(x,type_paroi_haute)
            dFFdx = np.tan(thetaI[NI]+muI[NI]) - (ParoiHaute(x+epsilon,type_paroi_haute)-ParoiHaute(x-epsilon,type_paroi_haute))/(2*epsilon)
            dx=FF/dFFdx
            x -= dx
        if compteur == IterMax : print('ATTENTION : non-convergence de Newton pour point P_{N-1}')
        #print("compteur (i=N-1) = %i, \t dx = %e"%(compteur, dx))
        XP[-1] = x
        YP[-1] = ParoiHaute(XP[-1],type_paroi_haute)
        thetaP[-1] = np.arctan((ParoiHaute(XP[-1]+epsilon,type_paroi_haute)-ParoiHaute(XP[-1]-epsilon,type_paroi_haute))/(2*epsilon))
         
        omegaP[-1] = lambdaPlusI[-1]+thetaP[-1]
        muP[-1] = mu_omega(omegaP[-1])
        lambdaMoinsP[-1] = omegaP[-1]+thetaP[-1]
        ax.plot([XP[-1],XI[-1]],[YP[-1],YI[-1]],'r')   # fin des C+ 
         
        # -------------------------------------------------------- 
        # Remplissage des tableaux pour presentation des resultats
        # --------------------------------------------------------
        Xmesh.append(deepcopy(XP)); Ymesh.append(deepcopy(YP)); mumesh.append(deepcopy(muP))
        Xbas.append(XP[0]);   mubas.append(muP[0])
        Xhaut.append(XP[-1]); muhaut.append(muP[-1])

        #fin de la boucle
        Xmin = min(XP)
        XcP.append(XP[Nmedian])
        YcP.append(YP[Nmedian])
        MachP.append(mu2mach(muP[Nmedian]))
    
    ax.plot([0,Ltot],[yMedian,yMedian],'k-.')   # axe central
    if show_Median: 
        ax.plot(XcP,YcP,'ko')
    # je regarde si o� se trouve les points sur l'axe de sym�trie
    # et s'ils tombent sur le maillage

    #****************************************************************
    #  MACH au milieu du canal, comparaison avec la théorie 1D
    #****************************************************************
    if show_Mach:

        # solution 1D
        H=ParoiHaute(0,type_paroi_haute)-ParoiBasse(0,type_paroi_basse)
        print("len Xcp = ",len(XcP))
        Ycp,Mach_th,Yhaut,Ybas=[],[],[],[]
        for x in XcP:
            YcP.append((ParoiHaute(x,type_paroi_haute)+ParoiBasse(x,type_paroi_basse))/2)
            S=(ParoiHaute(x,type_paroi_haute)-ParoiBasse(x,type_paroi_basse))/H
            Mach_th.append(inv_S_sur_Scrit(S*S_sur_Scrit(Mach0,gamma=gam),Mach=3.,show=False,gamma=gam))


        fig = plt.figure(1,figsize=(Largeur,Hauteur))
        fig.suptitle("Mach sur l'axe",fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        ax.set_ylabel(r'$x$',fontsize=14)
        ax.set_xlabel(r'$M$',fontsize=14)
        ax.plot(XcP,MachP,'k-*',linewidth=1,label='Mach axe')
        ax.plot(XcP,Mach_th,'r-',linewidth=1,label='Mach 1D')
        ax.axis()  
        ax.set_xlim(0, Ls)

     
    #***********************************************
    #   Lignes iso Mach,
    #***********************************************

    print("len Xmesh : ",len(Xmesh))
    #Yhaut=[ParoiHaute(x,type_paroi_haute) for x in list(Xhaut)]
    #Ybas=[ParoiBasse(x,type_paroi_basse) for x in list(Xbas) ]

    m=len(Xmesh) 
    X,Y = np.zeros(([m,N]),dtype=float),np.zeros(([m,N]),dtype=float)
    mac  = np.zeros(([m,N]),dtype=float)
    for i,x,y,muj in zip(range(m),Xmesh,Ymesh,mumesh):
        for j,xl,yl,mul in  zip(range(N),x,y,muj):
            X[i,j],Y[i,j],mac[i,j]=xl,yl,mu2mach(mul)

 
    plt.figure(3,figsize=(10, Hequal))
    ax = fig.add_subplot(111)
    levels=np.linspace(Mmin,Mmax,nlevels) 
    print('Levels = ',levels) 
    cp = plt.contour(X,Y,mac,levels,colors='black', linestyles='-')
    #plt.scatter(X, Y,color='black', s=6, marker='o', linewidth=1);  
    if show_label : plt.clabel(cp, inline=True,fontsize=8)
    plt.title("Lignes Iso-Mach")
    plt.plot(XL,YBAS,'k-',linewidth=1,label='bas')
    plt.plot(XL,YHAUT,'k-',linewidth=1,label='haut')
    plt.xlabel("x",fontsize=14)
    plt.ylabel("y",fontsize=14)
    plt.axis('equal')
    plt.xlim(0, Ls)
     
     #***********************************************
    #   Lignes iso Mach en couleur, surfaces remplies
    #***********************************************

    plt.figure(4,figsize=(12, Hequal))
    plt.title(r'iso-Mach')
    plt.contourf(X,Y,mac,levels)
    plt.axis('equal')
    #plt.scatter(X, Y, c=mac)
    plt.axis('equal')
    plt.xlim(0, Ls)
    plt.colorbar()

    plt.show() 
    return
    

def Exercice12_2():
    """
    Tube à choc, seulement les tracés,
    les calculs sont fait en fortran avec le code EC.
    """
    set_title("Solution du tube à choc, uniquement le tracé ")
    
    p=np.array([1, 6.3922135814136576,25.862283608680318 ]) 
    t=np.array([0,9.88327061176361155E-003 ,1.78214998216048130E-002])*1000

    print("t = ",t)
    print("p = ",p)
    # rapport P2/P1 =P3/P1     =    6.3922135814136576
    # rapport P5/P1 article    =    25.862283608680318
    # rapport P5/P2 article    =    4.0459041737714898   


    # ##################################################
    #    Temps sur le capteur 
    # ##################################################

    #  axe en t en ms, echelle   =    1.0000000000000000     
    #  Temps A   ( s)           =   1.23540882647045144E-002
    #  Temps 2   ( s)           =   9.88327061176361155E-003
    #  Temps 3   ( s)           =   1.78214998216048130E-002
    #  Temps 2 / Temps A        =   0.80000000000000004     
    #  Temps 3 / Temps A        =    1.4425588873701534     


    fig = plt.figure(1,figsize=(Largeur,Hauteur))
    fig.suptitle('Tube à choc', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    #fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$P/P_0$',fontsize=20)
    ax.set_ylabel(r'$t (ms)$',fontsize=20)
    ax.plot(p,t,'s')
    ax.axis([0,30,0,30])
    ax.grid()           
    plt.show()

    return

#************************************************************
#
#  Tuyère de longueur minimale
#
#************************************************************

def zone_complexe1(I,theta_I,omega_I,mu_I,CL,gam=gamma):
    """
    Zone complexe, d'une caractéristique positive vers un axe
    en entrée : valeurs de point I sur une C+
    CL        : condition à la limite sur l'axe pour theta
    en sortie : positions des points, omega, theta et mu
    """
    n=len(theta_I)
    x,y=np.zeros((n,n),dtype=float),np.zeros((n,n),dtype=float)
    ome,the=np.zeros((n,n),dtype=float),np.zeros((n,n),dtype=float)
    mut=np.zeros((n,n),dtype=float)
    for l in range(n):
        [x[0,l],y[0,l]]=I[l]
        ome[0,l],the[0,l],mut[0,l]=omega_I[l],theta_I[l],mu_I[l]
        #print(' l = %i, x,y = %f %f , theta = %f '%(l,x[0,l],y[0,l],the[0,l]))
    for k in np.arange(1,n):
        for l in np.arange(k,n):
            if k==l:
                the[k,l]=CL                                     # Condition à la limite (sur l'axe)
                ome[k,l]=ome[k-1,l]+the[k-1,l]-the[k,l]         # C-
                mut[k,l]=mu_omega(ome[k,l],gam=gam)
                p1=0                                            # axe
                p2=(the[k-1,l]-mut[k-1,l]+the[k,l]-mut[k,l])/2  # C-
                [x[k,l],y[k,l]]=Intersection_droites([x[0,0],y[0,0]],[x[k-1,l],y[k-1,l]],p1,np.tan(p2))
            else:
                ome[k,l],the[k,l]=Intersection_complexe1(the[k-1,l],ome[k-1,l],the[k,l-1],ome[k,l-1])
                mut[k,l]=mu_omega(ome[k,l],gam=gam)
                p1=(the[k,l-1]+mut[k,l-1]+the[k,l]+mut[k,l])/2   # C+
                p2=(the[k-1,l]-mut[k-1,l]+the[k,l]-mut[k,l])/2   # C-
                [x[k,l],y[k,l]]=Intersection_droites([x[k,l-1],y[k,l-1]],[x[k-1,l],y[k-1,l]],np.tan(p1),np.tan(p2))
    return x,y,ome,the,mut

def zone_paroi(A0,J):
    """
    Calcul de la zone II, détermination de la paroi
    à partir des C+ issues de J
    démarrage de la paroi en A0
    A= [x,y,omega,theta,mu,M]
    J liste de type A.
    """
    n=len(J);print('len J = ',n)
    x,y=np.zeros((n),dtype=float),np.zeros((n),dtype=float)
    Pw=[]
    Pw.append(A0)
    for l in range(n):
        Pw.append(deepcopy(J[l]))
        p1=(Pw[l][3]+J[l][3])/2              # theta moyen à la paroi
        p2=J[l][3]+J[l][4]
        [Pw[l+1][0],Pw[l+1][1]]=Intersection_droites([Pw[l][0],Pw[l][1]],[J[l][0],J[l][1]],np.tan(p1),np.tan(p2))
    print('len Pw = ',len(Pw))
    return Pw

def affichage_A(A,init=False):
    """
    Méthode des caractéristiques
    Affichage des données contenues dans A
    A[[x,y,omega,theta,mu,M],...]
    """
    #print("x = %5.3f, y = %5.3f, omega = %5.3e °, theta = %5.3e °, mu = %5.2f °,M = %3.4f "
    #      %(A[0],A[1],np.rad2deg(A[2]),np.rad2deg(A[3]),np.rad2deg(A[4]),A[5]))
    if init:
        print("  x\t   y  \t   omega(°) theta(°)  mu(°)   M ")
    else:
        print("%7.3f  %7.3f  %7.3f  %7.3f  %5.2f  %3.4f "%(A[0],A[1],np.rad2deg(A[2]),
            np.rad2deg(A[3]),np.rad2deg(A[4]),A[5]))
    return

def distribution(dist_log,eps,N,theta_a):
    """
    distribution angulaire des caractéristiques C- issues du point (a)
    """
    if dist_log:
        # distribution logarithmique au niveau de l'axe
        theta_log=np.linspace(np.log(eps*theta_a),np.log(theta_a),N)
        theta_l=np.exp(theta_log)
    else:
        # distribution des angles des caractéristiques partant de (a)
        theta_l=np.zeros((N),dtype=float)
        Nw=N-3
        theta_l[0],theta_l[1],theta_l[2],theta_l[3]=theta_a/Nw*eps, theta_a/(Nw*100.),theta_a/(Nw*10.),theta_a/(Nw*2.)
        dt=np.linspace(theta_l[2],theta_a,Nw)
        theta_l[4:N]=dt[1:]
    print(len(theta_l),'theta (°) = ',np.rad2deg(theta_l))
    return theta_l

def TuyereMinimalResults(gam=gamma):
    """
    Résultats après simulations pour gamma=1.4, 53 caractéristiques
    """    
    # #   x            y     M     L approxAiriau 
    # 2.351813  1.176137   1.5    2.432995  
    # 4.830457  1.687388   2.0    4.654693  
    # 6.685007  2.096271   2.25   6.240735  
    # 9.167570  2.636511   2.5    8.332294  
    # 12.487851 3.337441   2.75   11.111386 
    # 16.905614 4.234424   3      14.805188 
    # 30.383045 6.790221   3.5    26.129196 
    # 53.068746 10.721976  4.00   45.399018 
 
    Mach=[1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4 ]
    L   =[2.351813,4.830457,6.685007,9.167570,12.487851,16.905614,30.383045,53.068746]
    # H: coordonnées en y du dernier point de la paroi S = 2 x H
    H   =[1.176137,1.687388,2.096271,2.636511,3.337441,4.234424,6.790221,10.721976]
    Lapprox=[2.432995,4.654693,6.240735 ,8.332294,11.111386,14.805188 ,26.129196,45.399018 ]
    fig = plt.figure(4,figsize=(10,8))
    fig.suptitle('Analyse - Tuyère de longueur minimale',fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$M$',fontsize=14); ax.set_ylabel(r'$5 H,L$',fontsize=14)
    ax.plot(Mach,L,label="L")  
    ax.plot(Mach,Lapprox,label="L appro.") 
    ax.plot(Mach,[5*h for h in H],label="5 x H") 
    ax.plot(Mach,[5*S_sur_Scrit(m,gam) for m in Mach],marker="o", markersize=5,label="H 1D")
    plt.legend()
    plt.grid()


def Propulsion_Reduite(Pw,Ms,hs,hc=1,pic=1e5,gam=gamma,show=True):
     

    """
    calcul du coefficient de propulsion CF
    en fonction des points sur la paroi
    Pw[[x,y,omega,theta,mu,M],...]
    """
    m=len(Pw)
    x=[]
    M=[]
    t=[]
    for P in Pw:
        x.append(P[0])
        t.append(np.tan(P[3])*p_pi(P[5],gam))   # tan theta x p/pi
        if show: print(' x = %f, \t M = %f, theta (°) = %f'%(P[0],P[5],np.rad2deg(P[3])))
    tmp=integrate.trapz(t,x=x)
    CF=2*tmp/hc-p_pi(Ms,gam)*(hs/hc-1)
    CFbis=(gam*Ms**2*hs/hc+1)*p_pi(Ms,gam)-(gam+1)*p_pi(1,gam)

    print("Ms                     = %f"%(Pw[-1][5]))
    print("CF  (int p-pa)         =%2.4f"%(CF))
    print("CF  (momentum)         =%2.4f"%(CFbis))
    print("error on CF (x100)     =%2.3f"%(100*abs(CFbis/CF-1)))
    print("F                      = %f kN "%(CFbis*pic/1000*hc)) 
    print("hs/hc                  =%f"%(hs/hc))
    print("L                      =%f"%(x[-1]))
   

  
    return 


def Exercice12_4():
    """
    Tuyère de longueur minimale
    connaissant le nombre de Mach en sortie Ms
    """
    set_title("Tuyère de longueur minimale") 


    print_A          = False     # affichage du contenu de A
    print_I          = False    # affichage du contenu de I
    print_J          = False    # affichage du contenu de J
    affichage_figure = True    # affichage de la figure de construction
    dist_log         = False    # choix de la distribution du faisceau centré au point (a)

    Ms= 4.5                       # nombre de Mach en sortie
    N = 53                      # nombre de caractérististiques tirées du point (a)
    gam=1.4                     # valeur de gamma pour ce gaz parfait
    eps=1e-3                    # petit paramètre pour démarrer le calcul au voisinage du col
    A=[]                        # A[[x,y,omega,theta,mu,M],...] : liste contenant tous les points du problèmes
    J=[]                        # points entre la zone I et la zone II
    I=[]                        # première C+ issue de l'axe
    K=[]                        # points sur l'axe
    a=[0.,1.]                   # position du point (a)
    xCmoins=[]                  # abscisse des points sur les C- 
    yCmoins=[]                  # ordonnées des points sur les C-
    xCplus=[]                   # abscisse des points sur les C+ 
    yCplus=[]                   # ordonnées des points sur les C+

    # les trois valeurs suivantes sont pour évaluer la poussée mais n'interviennent pas dans le calcul
    # de la forme de la tuyère.

    h_s  = 0                    # hauteur de la section de sortie: calcul plus loin
    h_c  = 2*a[1]               # un exemple de section au col  (! la théorie est en 2D ici, pas en cylindrique)
    pi_c = 1e5                 # un exemple de pi au col
    S_c  = h_c

    set_question('1 : Initialisation')

    if print_A : affichage_A(A,init=True)
    # premier point (a) :
    theta_a=omega_super(Ms,gamma=gam)/2.      # valeur initiale au point (a) :
    omega_a=theta_a;mu_a=mu_omega(omega_a,gam=gam); Mach_a=mu2mach(mu_a)
    A.append([a[0],a[1],omega_a,theta_a,mu_a,Mach_a])
    if print_A : affichage_A(A[0])
    print('theta max (°)        = %f'%(np.rad2deg(theta_a)))

    # C- issues de (a)
    theta_l=distribution(dist_log,eps,N,theta_a)
    # quantités omega et mu, lambdaMoins=Cmoins associées :
    omega_l=theta_l
    mu_l=[]
    for omega in omega_l:
        mu_l.append(mu_omega(omega,gam=gam))
    Cmoins_l=omega_l+theta_l;
    print('Cmoins_l = ',np.rad2deg(Cmoins_l))


    #*******************************************************************
    # résolution de la première C+ ( I_1 I_N), calcul des points I et J
    #*******************************************************************
    omega_I,theta_I,mu_I,Mach_I=[],[],[],[]
    
    set_question('2 : calcul du point I_1')

    # I ne contient que des coordonnées
    # A[1:N] contiendra tous les points I
    I.append(Intersection_droites(a,[0.,0.],np.tan(theta_l[0]-mu_l[0]),0.))
    omega_I.append(Cmoins_l[0]); theta_I.append(0.) ; mu_I.append(mu_omega(omega_I[0],gam=gam))
    Mach_I.append(mu2mach(mu_I[0]))
    A.append([I[0][0],I[0][1],omega_I[0],theta_I[0],mu_I[0],Mach_I[0]])
    K.append(A[1])
    print(K)
     
    if print_I: print("I : omega = %5.3e °  \t theta = %5.3e ° \t mu = %5.2f ° \t  M = %3.4f "%(np.rad2deg(omega_I[0]),np.rad2deg(theta_I[0]),np.rad2deg(mu_I[0]),Mach_I[0]))
    if print_A : affichage_A(A[1])

    # toutes les C- démarrent en (a)
    xCmoins.append([a[0],I[0][0]]); yCmoins.append([a[1],I[0][1]])
    
    # initialisation de la figure
    plt.figure(0,figsize=(10,Hequal))
    plt.title('Tuyère de longueur minimale',fontsize=14, fontweight='bold')
    plt.ylabel(r'$y$',fontsize=14)
    plt.xlabel(r'$x$',fontsize=14)
    plt.plot([a[0],I[0][0]],[a[1],I[0][1]],'k-')
    plt.axis("equal")

    set_question('3 : calcul du point I_2 à I_N')
    
    for l in np.arange(1,N):
        ome,the=Intersection_complexe1(theta_l[l],omega_l[l],theta_I[l-1],omega_I[l-1])
        mue=mu_omega(ome,gam=gam)
        p1=(the+theta_l[l]-mue-mu_l[l])/2                   # pente de la C- issue de (a)
        p2=(the+theta_I[l-1]+mu_I[l-1]+mue)/2               # pente de la C+ issue de I(l-1)
        I.append(Intersection_droites(a,I[l-1],np.tan(p1),np.tan(p2)))
        omega_I.append(ome);theta_I.append(the);mu_I.append(mue); Mach_I.append(mu2mach(mu_I[l]))
        A.append([I[l][0],I[l][1],omega_I[l],theta_I[l],mu_I[l],Mach_I[l]])
        if print_I: print("I : omega = %5.3e °  \t theta = %5.3e ° \t mu = %5.2f ° \t  M = %3.4f "%(np.rad2deg(omega_I[l]),np.rad2deg(theta_I[l]),np.rad2deg(mu_I[l]),Mach_I[l]))
        if print_A: affichage_A(A[l+1])
        plt.plot([I[l-1][0],I[l][0]],[I[l-1][1],I[l][1]],'r--')    # C+
        plt.plot([a[0],I[l][0]],[a[1],I[l][1]],'k-')               # C-
    #print('len A = %i, \t N =%i '%(len(A),N)) 

     

    #*************************************************
    # résolution de la ligne sous la première C+ :
    #*************************************************

    set_question('4 : calcul de la zone sous la C+ (I_1,I_N), points J')

    x,y,ome,the,mut=zone_complexe1(I,theta_I,omega_I,mu_I,CL=0,gam=gam)
    J.append(A[-1])                                  # points délimitant la première zone  J[0]=I[N-1]
    
    # pour l'extraction des C- et C+
    for j in np.arange(1,N):
        xtmp,ytmp=[],[]
        xtmp.append(a[0]);ytmp.append(a[1])
        for xl,yl in zip(x[:j+1,j],y[:j+1,j]):
            xtmp.append(xl); ytmp.append(yl)
        xCmoins.append( xtmp); yCmoins.append( ytmp)
    
    for k in np.arange(1,N):
        xp,yp=[],[]
        plt.plot(x[k,k],y[k,k],'ko')
        plt.plot(x[k,N-1],y[k,N-1],'bs')
        J.append([x[k,N-1],y[k,N-1],ome[k,N-1],the[k,N-1],mut[k,N-1],mu2mach(mut[k,N-1])])
        for l in np.arange(k,N):
            plt.plot([x[k-1,l],x[k,l]],[y[k-1,l],y[k,l]],'k-')
            plt.plot([x[k,l-1],x[k,l]],[y[k,l-1],y[k,l]],'k--')
            A.append([x[k,l],y[k,l],ome[k,l],the[k,l],mut[k,l],mu2mach(mut[k,l])])
            if the[k,l] == 0: K.append(A[-1])
            xp.append(x[k,l]); yp.append(y[k,l])
        xCplus.append(xp); yCplus.append(yp)

    for l in range(N):
        plt.plot(x[0,l],y[0,l],'ro')

    if print_J: 
        affichage_A(A[0],init=True)
        print('valeur de J : ')
        for l in range(N):
              affichage_A(J[l])


    #*************************************************
    # résolution de la zone 2 et de la paroi  :
    #*************************************************
    # points sur la paroi Pw, ce sont les derniers points de A.
    Pw=zone_paroi(A[0],J)           # calculs des positions des points Pw, les autres valeurs sont celles des points J
    for P in Pw: A.append(P)
    print('valeurs de Pw (paroi) : ')
    affichage_A(A[0],init=True)
    for l in range(len(Pw)): affichage_A(Pw[l])
    for l in range(N): plt.plot( [Pw[l][0],Pw[l+1][0]] ,[Pw[l][1],Pw[l+1][1]],'k-')
    for l in range(N): plt.plot([J[l][0],Pw[l+1][0]],[J[l][1],Pw[l+1][1]],'r--')
    

    #*****************************************************************
    # initialisation de la figure 2, figure finale pour les C+ et C-
    #*****************************************************************

    plt.figure(1,figsize=(10,Hequal))
    plt.title('Tuyère de longueur minimale',fontsize=14, fontweight='bold')
    plt.ylabel(r'$y$',fontsize=14); plt.xlabel(r'$x$',fontsize=14)
    m=len(A)
 
    for l in range(N): plt.plot( [Pw[l][0],Pw[l+1][0]] ,[Pw[l][1],Pw[l+1][1]],'k-')   # paroi
    xi,yi=[],[]
    # tracé des C+
    # première C+ passant par les points I  
    for l in np.arange(1,N+1):
        xi.append(A[l][0]); yi.append(A[l][1])    
    xi.append(Pw[1][0]); yi.append(Pw[1][1])
    plt.plot(xi,yi,'b-')
    # les autres C+
    for l in range(N-1):
        xi,yi=[],[]
        for x,y in zip(xCplus[l],yCplus[l]):
            xi.append(x); yi.append(y)
        xi.append(Pw[l+2][0]); yi.append(Pw[l+2][1])
        plt.plot(xi,yi,'b-')
    # la dernière C+ venant uniquement de l'axe et délimitant la zone III. 
    plt.plot( [J[-1][0],Pw[-1][0]] ,[J[-1][1],Pw[-1][1]],'b-')
    # C-
    for l in range(N):
        xi,yi=[],[]
        for x,y in zip(xCmoins[l],yCmoins[l]):
            xi.append(x); yi.append(y)
        plt.plot(xi,yi,'r-')
    plt.axis("equal")
    
    L,R_s=A[-1][0],A[-1][1]
    h_s=2*R_s
    Lapprox=(h_c+h_s)/2*np.sqrt(Ms**2-1)
    print('Longueur de la tuyère                 : %f'%(L))
    print('Longueur, formule approchée           : %f\n'%(Lapprox))
    print('1/2 hauteur de la section de sortie   : %f'%(R_s))
    print('Rapport des sections  hs/hc           : %f'%(h_s/h_c))
    print('hs/hc, théorie 1D                     : %f'%(S_sur_Scrit(Ms)))
    print("différence en pourcent sur hs/hc      : %f"%(100*(h_s/h_c/S_sur_Scrit(Ms)-1)))
    print('pa/pi_c                               : %f'%(p_pi(Ms,gam)))

    #*****************************************************************
    # figure 3, nombre de Mach
    #*****************************************************************

    plt.figure(3,figsize=(10,8))
    plt.title('Tuyère de longueur minimale',fontsize=14, fontweight='bold')
    plt.ylabel(r'$M$',fontsize=14); plt.xlabel(r'$x$',fontsize=14)
    plt.plot([Pw[l][0] for l in range(N+1)] ,[Pw[l][5] for l in range(N+1)],'ko-',label="paroi")   # paroi
    plt.plot([I[l][0] for l in range(N)] ,[Mach_I[l] for l in range(N)],label="I")           # points I
    plt.plot([J[l][0] for l in range(N)] ,[J[l][5] for l in range(N)],label="J")             # points J
    plt.plot([K[l][0] for l in range(N)] ,[K[l][5] for l in range(N)],label="K")             # points J

    plt.legend()
    plt.grid()


    Propulsion_Reduite(Pw,Ms,h_s,hc=h_c,pic=pi_c,gam=gam,show=False)

    TuyereMinimalResults(gam)
    
    if affichage_figure: plt.show()
    return
