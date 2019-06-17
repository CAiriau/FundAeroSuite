#!/bin/python
"""
  Correction des exercices du chapitre 14 : hypersonique
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

# fonctions générales pour les chocs obliques :

def Constante_Kc(K0):
    """
    angle de choc, Kc= M0 sigma approximation en hypersonique
    """
    tmp=(gamma+1.0)/4.0
    return K0*(tmp+np.sqrt(tmp**2+1.0/K0**2))

def Pc_P0(Kc):
    """
    pc/p0
    """
    return (2.0*gamma*Kc**2-(gamma-1.0))/(gamma+1.0)

def Kp_hyper(Mach,p_p0):
    """
    Kp en hypersonique
    """ 
    return 2/(gamma*Mach**2)*(p_p0-1)

def approx_omega(Mach,M0):
    """
    On teste la validité de l'approximation
    """
    f1=omega_super(M0)-omega_super(Mach)
    f2=omega_hyper(M0)-omega_hyper(Mach)
    return np.rad2deg(f1),np.rad2deg(f2)

def Solution_Lees():
    """
    Solution de LEES en supersonique pour les cones
    dans l'ordre Kp/theta et sigma/theta
    """
    return 2*(gamma+1)*(gamma+7)/(gamma+3)**2,2*(gamma+1)/(gamma+3)

def distance_detachement_choc():
    """
    Relation pour calculer le détachement du choc, 
    Delta/R
    """
    epsilon=(gamma-1.)/(gamma+1.)
    Gam=np.sqrt((gamma-1.)*(gamma+3.))/(gamma+1)
    return epsilon/(1+Gam-2*epsilon)


def Mach_aval_hyper(Mach,Kc=0,sigma=0,regime='hypersonique'):
    """
    Nombre de Mach en aval du choc hypersonique, en fonction de Kc (sigma) et Mach
    En "supersonique" Kc=sigma", l'angle de choc en radians et K=M0*sin(sigma)
    """
    if regime=='hypersonique':
        K=Kc
    else:
        K=Mach*np.sin(sigma)
    

    Num=((gamma+1)*K*Mach)**2-4*(K**2-1)*(gamma*K**2+1)
    Den=(2*gamma*K**2-gamma+1)*(2+(gamma-1)*K**2)
    return np.sqrt(Num/Den)

def Mach_aval_hyper_from_K0(Mach,K0,Kc,methode=1):
    """
    Nombre de Mach en aval du choc hypersonique
    seconde méthode, avec K0 (theta) et Kc, et Mach
    """
    K0New=2*(Kc**2-1)/(Kc *(gamma+1))
    print("Vérification du K0 par rapport au Kc: nouveau K0 = %f, \t ancien = %f"%(K0New,K0))
    if methode==1:
        Den=(2.+(gamma-1.)*Kc**2)*(2*gamma*Kc**2+1-gamma)
        Num=(gamma+1)*Mach*Kc
        return Num/np.sqrt(Den)
    else:
        Num=np.sqrt((2+(gamma-1.)*Kc**2)/(2*gamma*Kc**2+1-gamma))
        return Mach/(Kc-K0)*Num

def Mach_aval_asymptotique(theta):
    """
    Valeurs asymptotiques, theta en radians
    """
    c=np.sqrt(2/(gamma*(gamma-1)))
    # deux approximations suivant qu'on néglige pou pas gamma
    # la première est la meilleure, pour des petits angles il y a peu de différence
    return c*np.sqrt(1/theta**2-gamma),c/theta

def Pj_Pc(Mc,deltaTheta):
    """
    P_j/P_c
    """
    tmp=1+(gamma-1)/2*Mc*deltaTheta
    expo=2*gamma/(gamma-1)
    return np.power(tmp,expo)

def rapport_pression_de_Kp(Kp,Mach):
    """
    p/p0 connaissant le Kp
    """
    return 1+gamma/2*Mach**2*Kp

# autres fonctions pour ce chapitre

def losange_hyper_coef(alpha,tau):
    """
    Profil losangique en hypersonique, solution numérique
    """
    # la zone 0 ne sert à rien, c'est pour avoir les indices qui correspondent au zone.
    delta=np.array([alpha , tau, -tau, -tau, tau])-alpha
    Kp=2*delta**2
    if alpha < tau :
        Kp[3:]=0
    else:
        Kp[1:4:2]=0
    #print('Kp = ',Kp[1:])
    CL=0.5*(Kp[2]+Kp[4]-Kp[1]-Kp[3])
    #print('CL = ',CL)
    CD=0.5*(-Kp[2]*delta[2]-Kp[4]*delta[4]+Kp[1]*delta[1]+Kp[3]*delta[3])
    #print('CD = ',CD)
    F=CL/CD
    return Kp,CL,CD,F


def losange_hyper_analytique(alpha):
    """
    Profil losangique en hypersonique, solution analytique
    alpha est adimensionné : alpha/tau
    """
    n=len(alpha)
    CL,CD=np.zeros([n]),np.zeros([n])
    polaire=np.zeros([n])
    for k in range(n):
        if alpha[k] < 1:
            CL[k],CD[k]=4*alpha[k],2+6*alpha[k]**2
            polaire[k]=2+3/8*CL[k]
        else:
            CL[k],CD[k]=2+2*alpha[k]**2,6*alpha[k]+2*alpha[k]**3
            polaire[k]=(4+CL[k])*np.sqrt(0.5*CL[k]-1)
    fmax=1/np.sqrt(3)
    CL_opt,CD_opt=4/np.sqrt(3),4
    alpha_opt=np.rad2deg(1/np.sqrt(3))
    CD_min=2

    return CL,CD,CL/CD,polaire,fmax,CL_opt,CD_opt,alpha_opt,CD_min


def triangle_hyper_coef(alpha,tau):
    """
    Profil losangique en hypersonique,  solution numérique
    """
    # la zone 0 ne sert à rien, c'est pour avoir les indices qui correspondent au zone.
    delta=np.array([alpha , tau, 0, -tau ])-alpha
    Kp=2*delta**2
    if alpha < tau :
        Kp[3]=0
    else:
        Kp[1],Kp[3]=0,0
    CL=Kp[2]-(Kp[1]+Kp[3])/2
    CD=(Kp[1]*delta[1]+Kp[3]*delta[3])/2-Kp[2]*delta[2]
    F=CL/CD
    return Kp,CL,CD,F

def triangle_hyper_analytique(alpha):
    """
    Profil triangle en hypersonique, solution analytique
    alpha est adimensionné : alpha/tau
    """
    n=len(alpha)
    CL,CD=np.zeros([n]),np.zeros([n])
    polaire=np.zeros([n])
    for k in range(n):
        if alpha[k] < 1:
            CL[k],CD[k]=alpha[k]**2+2*alpha[k]-1,1-3*alpha[k]+3*alpha[k]**2+alpha[k]**3
            polaire[k]=6+np.sqrt(2+CL[k])*(CL[k]-4)
        else:
            CL[k],CD[k]=2*alpha[k]**2,2*alpha[k]**3
            polaire[k]=pow(CL[k],3./2.)/np.sqrt(2)
    fmax=1.248534555
    CL_opt,CD_opt=0.911027108,0.729677128
    alpha_opt=40.46074596
    CD_min=6-4*np.sqrt(2)

    return CL,CD,CL/CD,polaire,fmax,CL_opt,CD_opt,alpha_opt,CD_min


#**********************
def Exercice14_0():
#**********************
    """ 
    test du calcul du Mach aval 
    """
    set_title("Calcul du Mach aval pour un choc oblique sur un dièdre par différentes approximations")
    M0=100
    theta=5.0

    print("théorie exacte choc oblique ")
    sigma,Maval,Paval,omegaAmont,omegaAval=Choc_Oblique(M0,theta,show=False)
    print("theta en degrés            :",theta)
    print("theta en radians           :",np.deg2rad(theta))
    print("sigma en radians           :",sigma)
    print("sigma/sin sigma            :",sigma/np.sin(sigma))
    print("(sigma-theta)/sin (sigma)  :",(sigma-np.deg2rad(theta))/np.sin(sigma))
    K0,Kc=M0*np.deg2rad(theta),M0*sigma
    print("K0                         :",K0)
    print("Kc                         :",Kc)

    M1_hyper= Mach_aval_hyper(M0,Kc=Kc,regime='hypersonique')
    M1_super= Mach_aval_hyper(M0,sigma=sigma,regime='supersonique')
    print("M1 approx. hypersonique    :",M1_hyper)
    print("M1 approx. supersonique    :",M1_super)


    print('autre calcul')
    Kc1=Constante_Kc(K0)
    print('Kc approximation M0 infini : ',Kc1)
    print('sigma (M0 infini )(°)      : ',np.rad2deg(Kc1/M0))

    # M0 grand, theta et sigma petit
    M1_1=Mach_aval_hyper_from_K0(M0,K0,Kc1,methode=1)
    M1_2=Mach_aval_hyper_from_K0(M0,K0,Kc1,methode=2)
    print('M1 méthode 1               : ',M1_1)
    print('M1 méthode 2               : ',M1_2)

    print('M1 asymptotiques : ')
    print(Mach_aval_asymptotique(np.deg2rad(theta)))


#**********************
def Exercice14_1():
#**********************
    """ 
    profil losangique et triangulaire, méthode de Newton 
    """
    set_title("profil losangique et triangulaire, méthode de Newton ")
    plot=True               # pour sortir les graphiques
    plot_num=True           # visualiser le calcul numérique
    plot_ana=True           # visualiser le calcul analytique
    pas=100                  # pour la solution analytique 1 point tous les "pas" est tracé
    plot_max=True           # visualiser le point de finesse maximale
    plot_gmtry=True        # visualiser les deux profils étudiés
    tau_losange=0.1
    tau_triangle=0.1
    Msize=8                 # taille des symboles pour les points de finesse maxima
    lw=2                    # épaisseur des traits
    n=1001
    alpha=np.linspace(0,2*tau_losange,n) 
    CL=np.zeros([n,2]); CD=np.zeros([n,2]);F=np.zeros([n,2]);IndMax=np.zeros([2],dtype=int)

    for k in range(n):
        KP1,CL[k,0],CD[k,0],F[k,0]=losange_hyper_coef(alpha[k],tau_losange)
        KP2,CL[k,1],CD[k,1],F[k,1]=triangle_hyper_coef(alpha[k],tau_triangle)
        
        IndMax[0],IndMax[1]=np.argmax(F[:,0]),np.argmax(F[:,1])
 
    print("Finesse maximale : Los. : %f, \t Tri. : %f"%(F[IndMax[0],0],F[IndMax[1],1])) 
    print("Incidence en °   : Los. : %f, \t Tri. : %f"%(np.rad2deg(alpha[IndMax[0]]), np.rad2deg(alpha[IndMax[1]])))  
    print("CL               : Los. : %f, \t Tri. : %f"%(CL[IndMax[0],0], CL[IndMax[1],1]))  
    
     
    print("Solution analytique")        
    lCL,lCD,lF,lpolaire,lfmax,lCL_opt,lCD_opt,lalpha_opt,lCD_min=losange_hyper_analytique(alpha/tau_losange)
    tCL,tCD,tF,tpolaire,tfmax,tCL_opt,tCD_opt,talpha_opt,tCD_min=triangle_hyper_analytique(alpha/tau_triangle)
    
    print("Finesse maximale : Los. : %f, \t Tri. : %f"%(lfmax/tau_losange,tfmax/tau_triangle)) 
    print("  Incidence en ° : Los. : %f, \t Tri. : %f"%(lalpha_opt*tau_losange,talpha_opt*tau_triangle))  
    print("  CL             : Los. : %f, \t Tri. : %f"%(lCL_opt*tau_losange**2, tCL_opt*tau_triangle**2))  
    print("  CD             : Los. : %f, \t Tri. : %f"%(lCD_opt*tau_losange**2, tCD_opt*tau_triangle**2))  
    print("CD minimal       : Los. : %f, \t Tri. : %f"%(lCD_min*tau_losange**2, tCD_min*tau_triangle**2))  
    
   
    xl=alpha/tau_losange
    xt=alpha/tau_triangle
    xL=xl[0:n:pas]
    

    z=np.rad2deg(alpha)


    if plot: 
        # CD et CL vs alpha
        fig = plt.figure(1)
        fig.suptitle(r'$\tilde C_L$  et $\tilde C_D$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde \alpha$',fontsize=20)
        ax.set_ylabel(r'$\tilde C_D, \tilde C_L$',fontsize=20)
        ax.grid()
        
        if plot_num:
            ax.plot(xl,CL[:,0]/tau_losange**2,'-',label=r'Los.: $\tilde C_L$',linewidth=lw)
            ax.plot(xl,CD[:,0]/tau_losange**3,'--',label=r'Los.: $\tilde C_D$',linewidth=lw)
            ax.plot(xt,CL[:,1]/tau_triangle**2,'-',label=r'Tri.: $\tilde C_L$',linewidth=lw)
            ax.plot(xt,CD[:,1]/tau_triangle**3,'--',label=r'Tri.: $\tilde C_D$',linewidth=lw)
        if plot_ana:
            ax.plot(xL,lCL[0:n:pas],'s',label=r'Los. an.: $\tilde C_L$')
            ax.plot(xL,lCD[0:n:pas],'o',label=r'Los. an.: $\tilde C_D$')
            ax.plot(xL,tCL[0:n:pas],'s',label=r'Tri. an.: $\tilde C_L$')
            ax.plot(xL,tCD[0:n:pas],'o',label=r'Tri. an.: $\tilde C_D$')
        if plot_max:
            ax.plot(np.deg2rad(lalpha_opt),lCL_opt,'s',markersize=Msize,label=r'Los.: $\tilde C_L(f_{\max})$')
            ax.plot(np.deg2rad(lalpha_opt),lCD_opt,'s',markersize=Msize,label=r'Los.: $\tilde C_D(f_{\max})$')
            ax.plot(np.deg2rad(talpha_opt),tCL_opt,'s',markersize=Msize,label=r'Tri.: $\tilde C_L(f_{\max})$')
            ax.plot(np.deg2rad(talpha_opt),tCD_opt,'s',markersize=Msize,label=r'Tri.: $\tilde C_D(f_{\max})$')
        ax.legend(loc='upper left')

        # finesse vs CL
        fig = plt.figure(2)
        fig.suptitle(r'Finesse $\times \tau$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde C_L$',fontsize=20)
        ax.set_ylabel(r'$\tilde f$',fontsize=20)
        ax.grid()
        if plot_num:
            ax.plot(CL[:,0]/tau_losange**2,F[:,0]*tau_losange,'-',label="Los.",linewidth=lw)
            ax.plot(CL[:,1]/tau_triangle**2,F[:,1]*tau_triangle,'-',label="Tri.",linewidth=lw)
        if plot_ana:
            ax.plot(lCL[0:n:pas],lF[0:n:pas],'s',label="Los. an.")
            ax.plot(tCL[0:n:pas],tF[0:n:pas],'o',label="Tri. an.")
        if plot_max:
            ax.plot(lCL_opt,lfmax,'s',markersize=Msize,label=r'Los.: $f_{\max}$')
            ax.plot(tCL_opt,tfmax,'s',markersize=Msize,label=r'Tri.: $f_{\max}$')
        ax.legend(loc='upper right')

        # finesse vs alpha
        fig = plt.figure(3)
        fig.suptitle(r'Finesse $\times \tau$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde \alpha$',fontsize=20)
        ax.set_ylabel(r'$\tilde f$',fontsize=20)
        ax.grid()
        if plot_num:
            ax.plot(xl,F[:,0]*tau_losange,'-',label="Los.",linewidth=lw)
            ax.plot(xl,F[:,1]*tau_triangle,'-',label="Tri.",linewidth=lw)
        if plot_ana:
            ax.plot(xL,lF[0:n:pas],'s',label="Los. an.")
            ax.plot(xL,tF[0:n:pas],'o',label="Tri. an.")
        if plot_max:
            ax.plot(np.deg2rad(lalpha_opt),lfmax,'s',markersize=Msize,label=r'Los.: $f_{\max}$')
            ax.plot(np.deg2rad(talpha_opt),tfmax,'s',markersize=Msize,label=r'Tri.: $f_{\max}$')
        ax.legend(loc='lower right')

        # polaire CD vs CL
        fig = plt.figure(4)
        fig.suptitle(r'Polaire $\tilde C_D$=f(\tilde C_L)', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde C_L$',fontsize=20)
        ax.set_ylabel(r'$\tilde C_D$',fontsize=20)
        ax.grid()
        if plot_num:
            ax.plot(CD[:,0]/tau_losange**3,CL[:,0]/tau_losange**2,'--',label=r'Los.',linewidth=lw)
            ax.plot(CD[:,1]/tau_triangle**3,CL[:,1]/tau_triangle**2,'--',label=r'Tri.',linewidth=lw)
        if plot_ana:
            ax.plot(lCD[0:n:pas],lCL[0:n:pas],'o',label=r'los. an.')
            ax.plot(tCD[0:n:pas],tCL[0:n:pas],'o',label=r'tri. an.')
        if plot_max:
            ax.plot(lCD_opt,lCL_opt,'s',markersize=Msize,label=r'Los.: $f_{\max}$')
            ax.plot(tCD_opt,tCL_opt,'s',markersize=Msize,label=r'Tri.: $ f_{\max}$')
            
        ax.legend(loc='upper left')

        # forme des deux profils
        if plot_gmtry:
            fig = plt.figure(5)
            fig.suptitle('profils', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/\ell$',fontsize=20)
            ax.set_ylabel(r'$y/\ell$',fontsize=20)
            ax.axis('equal')
            ax.grid()
            ax.plot([0,1/2,1,1/2,0],[0,tau_losange,0,-tau_losange,0],'-',linewidth=lw)
            dec=0.3
            ax.plot([0,1/2,1,0],[0+dec,tau_triangle+dec,0+dec,0+dec],'-',linewidth=lw)


        plt.show() 





#**********************
def Exercice14_3():
#**********************
    """ 
    Kp pour des écoulements de dièdre 2D ou coniques, pour un nombre de Mach infini et d'autres
    """
    set_title(" Kp pour des écoulements de dièdre 2D ou coniques, pour un nombre de Mach infini ")

    
    M0=100
    theta=np.array([5 , 10, 15, 20 ])
    sigma=np.zeros([4]); Kp=np.zeros([4]);
    
    set_question('1 : dièdre')

    for k in range(4):
        sigma[k]=valeur_sigma(M0,theta[k])
        rap_p=P2_P1(M0*np.sin(np.deg2rad(sigma[k])))
        Kp[k]=Kp_hyper(M0,rap_p)

    print('Pour un nombre de Mach fini M=%f, choc oblique exact : '%(M0))
    print("theta    : ", theta)
    print("sigma    : ", sigma)
    print("100 Kp   : ", Kp*100)

    print('Pour un nombre de Mach fini M=%f, avec approximation :'%(M0))
    Kc=Constante_Kc(np.deg2rad(theta)*M0)
    Kp=Kp_hyper(M0,Pc_P0(Kc))
    print("theta    : ", theta)
    print("sigma    : ",np.rad2deg(Kc/M0))
    print("100 Kp   : ",Kp*100)

    print('Pour un nombre de Mach infini  ')
    print("theta    : ", theta)
    print("sigma    : ",(gamma+1)/2*theta)
    print("100 Kp   : ",(gamma+1)*(np.deg2rad(theta))**2*100)


    set_question('2 : choc conique')
    """ 
    Obstacle conique en supersonique 
    """
    # calcul du Kp pour un choc conique, interpolation à partir d'un fichier
    set_title("Interpolation dans les chocs coniques en super/hypersonique")
    reference_filepath = os.path.join('Livre/Data', 'Mach_50_choc_conique.dat')
    lw=2
    m=1
    with open(reference_filepath, 'r') as infile:
        A=np.loadtxt(infile, dtype=float,usecols=(0,1,2),unpack=True,skiprows=m)
        theta1=A[0,:];kp_cc=A[2,:] ; sigma1=A[1,:]
    print(A.shape)
    print(theta1.shape)
    f1 = interp1d(theta1,kp_cc, kind='cubic')
    f2 = interp1d(theta1,sigma1, kind='cubic')
    theta_c=np.array([5.0, 10.0, 15.0,20.0])

    print('Pour un nombre de Mach infini  ')
    print("theta    : ", theta_c)
    print("sigma    : ",f2(theta_c))
    print("100 Kp   : ",100*f1(theta_c))

    set_question('3 : Méthode de Newton + Tracés')
    
    print('Kp = 2 sin^2(theta)  ')
    print("theta            : ", theta)
    print("100 Kp (sin^2)   : ",100*2*np.sin(np.deg2rad(theta))**2)
    print("100 Kp (theta^2) : ",100*2*(np.deg2rad(theta))**2)

    n=21
    theta=np.linspace(2,20,n)
    #KpDiedre,KpCone,KpNewtonSin,KpNewtonApprox=np.zeros([n]),np.zeros([n]),np.zeros([n]),np.zeros([n])
    #sigmaDiedre,sigmaCone=np.zeros([n]),np.zeros([n])

    sigmaDiedre=(gamma+1)/2*theta
    sigmaCone=f2(theta)
    KpDiedre=(gamma+1)*(np.deg2rad(theta))**2*100
    KpCone=100*f1(theta)
    KpNewtonSin=100*2*np.sin(np.deg2rad(theta))**2
    KpNewtonApprox=100*2*(np.deg2rad(theta))**2
    
    Kp_Lees,sigma_Lees=Solution_Lees()
    fig = plt.figure(5)
    fig.suptitle('Coefficients de pression', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$\theta (^\circ)$',fontsize=20)
    ax.set_ylabel(r'$100\times Kp$ ou $\sigma(^\circ)$',fontsize=20)
    ax.grid()
    ax.plot(theta,sigmaDiedre,'--',linewidth=lw,label=r"$\sigma$ : Dièdre")
    ax.plot(theta,sigmaCone,'--',linewidth=lw,label=r"$\sigma$ : Cône")
    ax.plot(theta,KpDiedre,'-',linewidth=lw,label=r"$Kp$ : dièdre")
    ax.plot(theta,KpCone,'-',linewidth=lw,label=r"$Kp$ : cône")
    ax.plot(theta,KpNewtonSin,'-',linewidth=lw,label=r"$Kp : 2\sin^2\theta$")
    ax.plot(theta,KpNewtonApprox,'-',linewidth=lw,label=r"$Kp : 2 \theta^2$") 
    ax.plot(theta,100*Kp_Lees*(np.deg2rad(theta))**2,'-',linewidth=lw,label=r"$Kp : Lees$")  
    ax.plot(theta,sigma_Lees*theta,'--',linewidth=lw,label=r"$\sigma$ : Lees") 

    ax.legend(loc='upper left') 
    plt.show() 

    set_question('3 :Solution de LEES')
    theta=np.array([5 , 10, 15, 20 ])
    theta=np.linspace(5,20,16)
    print('Solution indépendante du nombre de Mach ')
    print("theta    : ", theta)
   
    print("sigma (°): ",sigma_Lees*theta)
    print("100 Kp   : ",100*Kp_Lees*(np.deg2rad(theta))**2)
    print("erreur Kp :",np.abs(Kp_Lees*(np.deg2rad(theta))**2/f1(theta)-1))
    print("erreur sigma :",np.abs(sigma_Lees*theta/f2(theta)-1))


def ogive_parabolique(e_sur_L,x,pente="angle"):
    """
    Profil, dérivée, angle, rayon de courbure
    """
    y=-e_sur_L*2.*x*(x-1)        # équation du profil
    yp=(-4*x+2)*e_sur_L         # dérivée première
    ypp=-4*e_sur_L              # dérivée seconde
    R=abs(pow(1+yp**2,3/2)/ypp)      # rayon de courbure
    # pente en radian       
    if pente=="angle":
        th=yp                        # on suppose les angles petits...
    else :
        th=np.arctan(yp)    
    #print('Theta en ° =',np.rad2deg(th))
    return y,th,R

def Kp_Newton(e_sur_L,C,x):
    """
    Coefficient de pression en hypersonique : formule de Newton
    """
    y,th,R=ogive_parabolique(e_sur_L,x)
    print("Newton")
    Kp=C*np.sin(th[th>0])**2
    for k in np.arange(len(Kp),len(x)):
        Kp=np.append(Kp,[0])
    return Kp

def Kp_Newton_Busemann(e_sur_L,C,x):
    """
    Coefficient de pression en hypersonique : formule de Newton +
    correction de Busemann
    """
    y,th,R=ogive_parabolique(e_sur_L,x)
    print("Newton Busemann")
    ind=th>0
    Kp=C*np.sin(th[ind])**2+2*y[ind]/R[ind]
    for k in np.arange(len(Kp),len(x)):
        Kp=np.append(Kp,[0])
    return Kp

#**********************
def Exercice14_4():
#**********************
    """ 
    Comparaison des méthodes de calcul autour d'un obstacle pointus
    """
    # 
    set_title("Comparaison des méthodes de calcul autour d'un obstacle pointus")
 
    plot_profil =False
    plot_courbes=True
    plot_Newton =True
     
    n=201   # nombre de points sur l'extrados du profil
    Mach_inf=np.array([3.5,10, 100.0]);m=len(Mach_inf)
    e_sur_L=0.10

    theta,x,Kp=np.zeros(n),np.zeros(n),np.zeros(n)
    ye,yi=np.zeros(n),np.zeros(n)
    x=np.linspace(0,1,n)

    #print('n = ',len(x))
  
    #---------------------------------------
    set_question('1 : profil lenticulaire')
    #---------------------------------------

    ye,theta,Re=ogive_parabolique(e_sur_L,x,pente="angle")
    yi=-ye
    # for i in range(n-1):
    #     theta[i]=np.arctan((ye[i+1]-ye[i])/(x[i+1]-x[i]))
    #     #print(i,theta[i])
    print('x(0)= %f, theta(0) = %f '%(x[0],theta[0]))
    

    if plot_profil: 
        fig = plt.figure()
        fig.suptitle('Profil', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        ax.axis('equal')
        ax.grid()
        ax.plot(x,ye,'-')
        ax.plot(x,yi,'-')
        ax.plot(x[:-1],theta[:-1],'o--')
        plt.show() 

    #-----------------------------------------------------    
    set_question("2: méthode choc-détente, Mach = infini") 
     #-----------------------------------------------------        
    
    theta_min=min(theta-theta[0])  
    print('Theta minimal du profil         : %f, %f °'%(theta_min,np.rad2deg(theta_min)))
  
    Kp=np.zeros([m,n])

    # traitement en théorie hypersonique


    print("\n ","-"*49,"\n RESULTAT POUR UN NOMBRE DE MACH INFINI \n","_"*50)

    sigma= theta[0]*(gamma+1)/2
    KpcInf= (gamma+1)*theta[0]**2
    Mc,tmp  = Mach_aval_asymptotique(theta[0])
    Kp_Infini=(gamma+1)*theta[0]**2*Pj_Pc(Mc,theta-theta[0])
    print("sigma infini                    : %7.4f °"%(np.rad2deg(sigma)))  
    print('100 x Kpc théorique  M_infini   : %7.4f '%(100*KpcInf))
    print("Mach en aval (M0 inf, théo.)    : %7.4f "%(Mc))
    print(" theta (°) \t 100 Kp")
    for thetap,kp in zip(theta,Kp_Infini):
        print(" %7.4f \t %7.4f"%(np.rad2deg(thetap),100*kp))
    print(100*Kp_Infini)
  
    #------------------------------------
    set_question("3: méthodes de Newton")  
    #------------------------------------  

     
    C=KpcInf/np.sin(theta[0])**2
    print('C = %f, C théorique = %f'%(C,gamma+1))
    #C=gamma+1
    Kp_N=Kp_Newton(e_sur_L,C,x)
    Kp_NB=Kp_Newton_Busemann(e_sur_L,C,x)
     
    print(" theta (°) \t 100 Kp_N \t 100 Kp_NB \t 100 Kp_N approché")
    for thetap,kp1,kp2 in zip(theta,Kp_N,Kp_NB):
        if thetap >=0:
            kp3=(gamma+1)*thetap**2
        else :
            kp3=0
        print(" %7.4f \t %7.4f \t %7.4f \t %7.4f"%(np.rad2deg(thetap),100*kp1,100*kp2,100*kp3))
 
    #------------------------------------------------------    
    set_question("4: méthodes pour un nombre de Mach FINI")  
    #------------------------------------------------------


    print("\n ","-"*49,"\n RESULTAT POUR UN NOMBRE DE MACH Amont FINI \n","_"*50)

    for k in range(m):      
        M0=Mach_inf[k]     
        Theta_limite=-2/(gamma-1)/M0
        print("\n","*"*49)
        print("Mach infini M0                  : %7.4f "%(M0))
        print("*"*50,"\n")
        print('Theta minimal autorisé          : %f, %f °'%(Theta_limite,np.rad2deg(Theta_limite)))
   
        # if theta_min < Theta_limite:
        #     x_bad=x[theta-theta[0]<Theta_limite]
        #     print("PROBLEME : LA THEORIE DE FONCTIONNE PLUS: DIMINUER M0 ou l'épaisseur du profil")
        #     print("A partir de x = \n",x_bad[0])
        print("Angle au bord d'attaque         : %7.4f ° , %7.4f rad "%(np.rad2deg(theta[0]),theta[0]))
        print('sin theta                       : %7.4f'%(np.sin(theta[0])))
        K0=M0*theta[0]
        Kc=Constante_Kc(K0)

        print("")
        print("THEORIE DU CHOC OBLIQUE :")
        sigma,Maval,Paval,omegaAmont,omegaAval=Choc_Oblique(M0,np.rad2deg(theta[0]),show=False)
        print("Kpc                             : %7.4f "%(Kp_hyper(M0,Paval)))
        print("sin(sigma-theta)                : %7.4f"%(np.sin(sigma-theta[0])))
        Mc=Mach_aval_hyper(M0,Kc)
        print()
        #print("sigma (°)                       : %7.4f °"%(np.rad2deg(sigma)) )
        #print("Mach aval                       : %7.4f  "%(Maval) )
        #print("Paval/Pamont                    : %7.4f  "%(Paval) )
       

        print("")
        print("THEORIE DE L'HYPERSONIQUE :")
        print("K0 au bord d'attaque            : %7.4f "%(K0))
        print('Kc après le choc                : %7.4f '%(Kc))
        print('(Kc-K0)/M0                      : %7.4f '%((Kc-K0)/M0))
        sigma=Kc/M0
        print("Angle de choc au bord d'attaque : %7.4f °"%(np.rad2deg(sigma)))
        Kpc=Kp_hyper(M0,Pc_P0(Kc))
        print("Kpc calcul                      : %7.4f "%(Kpc))
        Mc=Mach_aval_hyper(M0,Kc)
        print("Mach en aval du bord d'attaque  : %7.4f "%(Mc))
        
        
        Kp[k,0]=Kp_hyper(M0,Pc_P0(Kc))
        for i in np.arange(1,n):
            Kp[k,i]=Kp_hyper(M0,Pj_Pc(Mc,theta[i]-theta[0])*Pc_P0(Kc))
            #print('Kp[ %d ] = %f'%(i,Kp[i]))

        print("")
        print("THEORIE DE NEWTON :")
        C=Kpc/np.sin(theta[0])**2
        print('C = %f'%(C))
        Kp_N=Kp_Newton(e_sur_L,C,x)
        Kp_NB=Kp_Newton_Busemann(e_sur_L,C,x)
     
        print(" theta (°) \t 100 Kp_N  \t 100 Kp hypersonique")
        for thetap,kp1,kp2 in zip(theta,Kp_N,Kp[k,:]):
            print(" %7.4f \t %7.4f \t %7.4f "%(np.rad2deg(thetap),100*kp1,100*kp2))
 

    
    # data de la méthode des caractéristiques (A vérifier) 
    x_10=np.linspace(0,1,11)
    Kp_10=np.array([0.10812, 0.07437, 0.04750, 0.02656, 0.01125, 0.001562,
         -0.005, -0.009375, -0.012187, -0.013437, -0.014375])    
    Kp_inf=np.array([0.092813, 0.058125, 0.031875, 0.01625, 0.007812, 0.002187,
        0.0,-0.000625,-0.0009375,-0.001094,-0.00125])
    Kp_35=np.array([0.17477,0.13169,0.0923,0.0578,0.02769,0.00246,
            -0.01908,-0.03692,-0.05292,-0.06708,-0.07754 ])
   
    #------------------------------------------------------    
    set_question("5: méthode choc détente en supersonique")  
    #------------------------------------------------------

    # traitement en théorie supersonique
    print("")
    print('CAS SUPERSONIQUE (Rappel + Kp)\n')
    Kp_super=np.zeros([m,n])
    for k in range(m):      
        M0=Mach_inf[k]
        sigma=np.deg2rad(valeur_sigma(M0, np.rad2deg(theta[0]),sigma=0,show=False))
        print('-----')
        print("Mach infini M0                  : %7.4f "%(M0))
        print("Angle du choc                   : %7.4f °"%(np.rad2deg(sigma)))
        Mn=M0*np.sin(sigma)
        MnAval= Mach_Aval(Mn)
        Mc=MnAval/np.sin(sigma-theta[0])
        print("Mach aval Mc                    : %7.4f "%(Mc))
        Pc=P2_P1(Mn)
        print("Pression Pc                     : %7.4f "%(Pc))
        omega_c=omega_super(Mc)
        print("Angle omega_c                   : %7.4f °"%(np.rad2deg(omega_c)))

        Kp_super[k,0]=Kp_hyper(M0,Pc)
        for i in np.arange(1,n):
            Mach=invomega(omega_c+theta[0]-theta[i])
            Kp_super[k,i]=Kp_hyper(M0,p_pi(Mach)/p_pi(Mc)*Pc)

    #------------------------------------------------------    
    set_question("6: graphiques")  
    #------------------------------------------------------

    Msize=10
    if plot_courbes :
        # exemple de graphique
         
        fig = plt.figure()
        fig.suptitle('Kp en %', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.set_ylabel(r'$100 Kp$',fontsize=20)
        ax.grid()
        ax.plot(x_10,100*Kp_35,'ko',markersize=Msize,label=r'carac. : $M=3.5$')
        ax.plot(x_10,100*Kp_10,'ko',markersize=Msize,label=r'carac. : $M=10$')
        ax.plot(x_10,100*Kp_inf,'ko',markersize=Msize,label=r'carac.. : $M=\infty$')
        for k in range(m-1):
            ax.plot(x,100*Kp[k,:],'-',label=r'Thé. : $M = %2.1f$'%(Mach_inf[k]))
        ax.plot(x,100*Kp_super[0,:],marker='<',markersize=Msize,markevery=(0,0.1),label=r'Super : $M = %2.1f$'%(Mach_inf[0]))
        #ax.plot(x,100*Kp_super[1,:],marker='<',markersize=Msize,markevery=(0,0.1),label=r'Super : $M = %2.1f$'%(Mach_inf[1]))
        if plot_Newton:
            ax.plot(x,100*Kp_N,'-',label=r'The. : Newton')
            ax.plot(x,100*Kp_NB,'-',label=r'The. : Newton-Busemann')
        ax.plot(x,100*Kp_Infini,'-',label=r'$M=\infty$',linewidth=2)
        ax.legend(loc='upper right')
        plt.show() 

    print(' x \t theta (°) \t Kp (%) \t Kp (%) \t Kp (%) \t Kp(%) \t \t  Kp(%) \t Kp(%) \t \t theta (rad)')
    print("M0:\t\t\t %3.2f \t\t %3.2f \t\t %3.2f"%(Mach_inf[0],Mach_inf[1],Mach_inf[2]),"\t Infini \t Newon \t \t New. Bus.")
    for k in range(n):
        m=Kp[:,k]
        print( "%4.2f \t %+6.2f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f "%(x[k],np.rad2deg(theta[k]),
            100*Kp[0,k],100*Kp[1,k],100*Kp[2,k],100*Kp_Infini[k],100*Kp_N[k],100*Kp_NB[k],theta[k]))
        if k%5==0:
            print("")   

    # Kp de Newton et Newton+Busemann sont indépendant du Mach amont .
    # pour M0=3.5 il faut reprendre les formules en supersonique (pas implémenter ici)
    
    #print('test : ')
    #s=' '.join(['%f' % (Mach_inf[k]) for k in np.arange(len(Mach_inf))])
    #print(s)


def fx(p):
    """
    fonction principale qui donne la loi de x fonction de y'
    """
    return np.log(p)+1/p**2+3/(4*p**4)

def fr(p):
    """
    fonction principale qui donne la loi de r fonction de y'
    """
    return (1+p**2)**2/p**3

def r_approche(x,C1):
    """
    solution approchée du profil, à l'infini
    """
    return pow(C1,1/4)*pow(4./3.,3/4)*pow(x,3/4)

#**********************
def Exercice14_5():
#**********************
    """ 
    forme d'un corps de traînée minimale en hypersonique
    """
    # 
    set_title("forme d'un corps de traînée minimale en hypersonique")

    
    plot_profil =True
    plot_courbes=True
    displayTable=True
    r0,x0=.01,0.01
    p1,p0=0.1,np.sqrt(3)
    
    n=1001   # nombre de points sur l'extrados du profil

    
    r,x,p=np.zeros(n),np.zeros(n),np.zeros(n)
    erreur=np.zeros(n)
    p=np.linspace(p1,p0,n)

    X,R=fx(p),fr(p)
    C1=r0/fr(p0)
    C2=x0-C1*fx(p0)
    x=C1*X+C2
    r=C1*R
    r_a=r_approche(x,C1)
    erreur=np.abs(r_a/r-1)*100

    if displayTable:
        print("    p \t \t X \t R \t x \t     r \t\t r_a \t erreur")
        for k in range(n):
            print("%9.4f "*7 %(p[k],X[k],R[k],x[k],r[k],r_a[k],erreur[k]))
           
    print(" à x0 : fx = %f, fr=%f"%(fx(p0),fr(p0)))
    print(" à x1 : fx = %f, fr=%f"%(fx(p1),fr(p1)))
    print("  %f <= x <= %f et %f <= r <= %f"%(x[-1],x[0],r[-1],r[0]))
    print("C2 = %e, C1 = %e"%(C2,C1))
    print("1/C2 = %e, 1/C1 = %e"%(1/C2,1/C1))
    print("épaisseur relative  = %f "%(0.5*(r[1]-r[0])/(x[1]-x[0]) ))

    if plot_profil:
        fig = plt.figure(num=1)
        fig.suptitle('forme optimale', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\ln x$',fontsize=20)
        ax.set_ylabel(r'$\ln r$',fontsize=20)
        ax.grid()
        ax.plot(np.log(x),np.log(r),'b',linewidth=2,label=r'$r(x)$')
        c=np.log(r[0])-3/4*np.log(x[0])
        ax.plot(np.log(x),c+3/4*np.log(x),'r--',linewidth=2,label=r'$\frac{3}{4}\ln x + c$')
        ax.legend(loc='upper left')

        fig = plt.figure(num=2,figsize=(8, 3))
        fig.suptitle('forme optimale', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$r$',fontsize=20)
        ax.grid()
        ax.plot(x,r,'b',label=r'$r(x)$',linewidth=2)
        ax.plot(x,-r,'b',label=r'$r(x)$',linewidth=2)
        #ax.plot(x,r_a,'r',label=r'$r_a(x)$')
        #ax.legend(loc='upper left')
        ax.axis('equal')
        plt.show()
        
    if plot_courbes:
        fig = plt.figure(3)
        fig.suptitle('fonction X(p)  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r"$r'$",fontsize=20)
        ax.grid()
        ax.plot(X,p,'b',label=r'$X(x)$')
        #ax.legend(loc='upper left')

        fig = plt.figure(4)
        fig.suptitle('fonction R(p)  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$r$',fontsize=20)
        ax.set_ylabel(r"$r'$",fontsize=20)
        ax.grid()
        ax.plot(p,R,'b',label=r'$r(x)$')
        #ax.legend(loc='upper left')

        fig = plt.figure(5)
        fig.suptitle(r'fonction $x(p)$  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$p$',fontsize=20)
        ax.set_ylabel(r"$\ln x$",fontsize=20)
        ax.grid()
        ax.plot(p,np.log(x),'b',label=r'$r(x)$')
        ax.legend(loc='upper left')

        fig = plt.figure(6)
        fig.suptitle('fonction r(p)  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$p$',fontsize=20)
        ax.set_ylabel(r"$r$",fontsize=20)
        ax.grid()
        ax.plot(p,r,'b',label=r'$r(x)$')
        #ax.legend(loc='upper left')
        

        plt.show()

#**********************
def Exercice14_6():
#**********************
    """ 
    traînée d'un corps émoussé en hypersonique
    """ 
    set_title("traînée d'un corps émoussé en hypersonique")

    # Paramètres et options du problème :
    plot_courbes=True 
    displayTable=True
    M0=np.float64(20.)
    Lr_sur_Rq=np.float64(5.)
    na=11
    theta_d=np.deg2rad(np.float64(15.))
    R,Rq=np.float64(1.),np.float64(2.)


    set_question('1a : Calcul du Kp pour le corps axisymétrique')

    print("partie arrondie : ")
   
    Sref=np.pi*Rq**2


    theta=np.linspace(np.pi/2,theta_d,na)
    Kp_Na=2.0*np.sin(theta)**2
    print(" theta (°) \t Kp")
    for k in range(na):
        print(" %5.4f \t %5.4f "%(np.rad2deg(theta[k]),Kp_Na[k]))
    print("Partie dièdre 3D, cône tangent hypersonique : ")

    
    """ 
    Obstacle conique en supersonique 
    """
    # calcul du Kp pour un choc conique, interpolation à partir d'un fichier
    set_title("Méthode du cone tangent")
    reference_filepath = os.path.join('Livre/Data', 'Mach_20_choc_conique.dat')
    lw=2
    m=1
    with open(reference_filepath, 'r') as infile:
        A=np.loadtxt(infile, dtype=float,usecols=(0,1,2),unpack=True,skiprows=m)
        theta1=A[0,:];kp_cc=A[2,:] ; sigma1=A[1,:]
    print(A.shape)
    print(theta1.shape)
    f1 = interp1d(theta1,kp_cc, kind='cubic')
    f2 = interp1d(theta1,sigma1, kind='cubic')
    theta_c=np.linspace(50,np.rad2deg(theta_d),8)


    print('Pour un nombre de Mach M0       : %5.4f'%(M0))
    print(" theta (°) \t  Kp_n \t \t  sigma (°) \t Kp ")
    for the in theta_c:
        print(" %5.3f \t %5.3f \t \t  %5.3f \t %5.3f "%(the,2.0*np.sin(np.deg2rad(the))**2,f2(the),f1(the)))
    


    Kp_Lees,sigma_Lees=Solution_Lees()
    Kp_LEES=Kp_Lees*theta_d**2
    print("Solution de Lees, Kp            : %5.4f"%(Kp_LEES))
     

    print("distance de détachement Delta/R : %5.4f "%(distance_detachement_choc()))

    set_question('1b : Calcul du CD pour le corps axisymétrique')
    
     
    # partie émoussée
    CD_emoussee= np.pi/Sref*R**2*(1-np.sin(theta_d)**4)
    # partie cone
    Rc=R*np.cos(theta_d)
    CD_cone    = 2*np.pi/Sref*(Rq**2-Rc**2)*np.sin(theta_d)**2
    # partie culot
    Kp_q=0.
    CD_culot   = -Kp_q

    CD=CD_emoussee+CD_cone+CD_culot

    CD_ConeCulot=Kp_LEES *np.pi/Sref*(Rq**2-Rc**2)
    CD1=CD_emoussee+CD_ConeCulot

    print("calcul1 : CD = %5.4f\t CD_e = %5.4f\t CD_cone = %5.4f \t CD_culot = %5.4f "%
        (CD,CD_emoussee,CD_cone,CD_culot))

    print("calcul2 : CD = %5.4f\t CD_e = %5.4f\t CD_(cone+culot) = %5.4f"%
        (CD1,CD_emoussee,CD_ConeCulot))


    set_question('2a : Calcul du Kp pour la géométrie plane')
    theta_max_2D=np.deg2rad(45.293184)
    theta_2d=np.linspace(45,np.rad2deg(theta_d),7)
    print('Angles = ',theta_2d)
    print('Kp Méthode de Newton ')
    print(" theta (°) \t Kp_n \t Kp_dt (approx) \t K0 \t Kp_dt")
    for the in theta_2d:
        tm=np.deg2rad(the)
        print(" %2.2f \t %5.3f  \t %5.3f  \t %+05.3f \t %4.3f"%(the,2*np.sin(tm)**2,
         (gamma+1)*tm**2,M0*tm,2*tm*Constante_Kc(M0*tm)/M0))
       
  
    set_question('3 : Calcul du Kp au culot avec sillage')

    theta_s=-np.arctan(1./Lr_sur_Rq)
    Delta_theta=theta_s-theta_d
    print('Angle de déviation du sillage      : %6.3f'%(np.rad2deg(theta_s)))
    print('Angle de déviation dans la détente : %6.3f'%(np.rad2deg(Delta_theta)))
    #Kp1=Kp_LEES
    Kp1=Kp_Na[-1]
    P1_P0=1+gamma*M0**2/2*Kp1
    print('Kp1                                : %6.4f'%(Kp1))
    print('P1/P0                              : %6.4f'%(P1_P0))
    
    #print(p_pi(M0),inverse_p_pi(p_pi(M0)))
    P1_Pi1=P1_P0*p_pi(M0)/pi2_pi1(M0)
    print('rapport des pi à travers le choc   : %e'%(pi2_pi1(M0)))
    print('p1/pi1                             : %f'%(P1_Pi1))
    M1=inverse_p_pi(P1_Pi1)
    print('M1                                 : %6.4f'%(M1))
    print('omega(M1)                          : %6.4f'%(np.rad2deg(omega_super(M1))))
    omega2=omega_super(M1)+np.abs(Delta_theta)
    M2=invomega(omega2)
    print('omega(M2)                          : %6.4f'%(np.rad2deg(omega2)))
    print('M2                                 : %e'%(M2))
    P2_sur_P1=p_pi(M2)/p_pi(M1)
    print('P2_P1                              : %6.4f'%(P2_sur_P1))
    print('Kp2                                : %6.4f'%(Kp_hyper(M0,P2_sur_P1*P1_P0)))
    print('P2/P0                              : %e'%(P2_sur_P1*P1_P0)) 
    print('P0/P1                              : %e'%(1/P1_P0)) 

    opt=True
    if plot_courbes:
        fig = plt.figure(num=1)
        fig.suptitle('Kp sur le corps émoussé', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$Kp$',fontsize=20)
        ax.grid()
        na=75
        Lc=(Rq-Rc)/np.tan(theta_d)
        print('Lc = ',Lc)
        theta=np.linspace(np.pi/2,theta_d,na)
        Kp1=2.0*np.sin(theta)**2
        x1=(1-np.sin(theta))*R
        
        tc=np.deg2rad(theta_c)
        x2=(1-np.sin(tc))*R
        Kp2=f1(theta_c)
       
        tm=np.deg2rad(theta_2d)
        x3=(1-np.sin(tm))*R
        Kp3=2*tm*Constante_Kc(M0*tm)/M0
       
        
        if opt:
            Kp1=np.append(Kp1,Kp1[-1])
            x1=np.append(x1,x1[-1]+Lc)
            Kp2=np.append(Kp2,Kp2[-1])
            x2=np.append(x2,x2[-1]+Lc)
            Kp3=np.append(Kp3,Kp3[-1])
            x3=np.append(x3,x3[-1]+Lc)

        ax.plot(x1,Kp1,'b',linewidth=2,label=r'$Kp_n$')
        ax.plot(x2,Kp2,'r',linewidth=2,label=r'$Kp_{ct}$')
        ax.plot(x3,Kp3,'k',linewidth=2,label=r'$Kp_{dt}$')
        ax.legend(loc='upper right')
        ax.axis([0,1.,0,2])
        print("x_init : cone = %f, dièdre = %f"%(x2[0],x3[0]))
        print("corde= ",x1[-1])
        
        #ax.plot(np.log(x),c+3/4*np.log(x),'r--',linewidth=2,label=r'$\frac{3}{4}\ln x + c$')
        plt.savefig("fig14-exo6_solution.svg")

        plt.show()
#**********************
def Exercice14_100():
#**********************
    """
    Courbes des chocs coniques du rappel
    """    
    set_title("Courbes des chocs coniques du rappel")
    
    """
    Mach aval, Kp, sigma en fonction de theta_c
    """
    main_multiple_Mach(np.array([5,10,20,50]),["sigma","Mach_aval","Kp"],Lees=True,
        t=[20,2,1],s=[26,5,1],Kp=[30,5,1],Mav=[50,5,1]) 


#**********************
def Exercice14_102():
#**********************
    """
    Test de l'inversion de la fonction omega(M) en supersonique
    """
    set_title("Test inversion de la fonction  omega(Mach) en supersonique")
    M0=3
    theta=omega_super(M0)
    print("pour M0=%f on obtient theta = %f °"%(M0,np.rad2deg(theta)))
    print(invomega(theta,show=True))

    print(valeur_sigma(M0, 10,sigma=0,show=True))

#**********************
def test_approx_omega(cas):
#**********************
    """
    test de l'approximation sur la fonction de Prandtl-Meyer en hypersonique
    cas = 0 : plusieurs M0, cas =1 , 1 M0, cas = 2: fonction omega
    """
    
    n=101;
    M_init,M_end=10,15
    
     
    fig = plt.figure()
    
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$M$',fontsize=20)
   
    ax.grid()
    if cas==0:
        M0=np.array([3.5, 5, 10]);m=len(M0)
        erreur=np.zeros([m,n,2])
        M=np.zeros([m,n])
        for k in range(m):
            M[k,:]=np.linspace(M0[k],M_end,n)
            erreur[k,:,0],erreur[k,:,1]=approx_omega(M[k,:],M0[k])
        fig.suptitle(r'$\theta(M)$', fontsize=14, fontweight='bold')
        ax.set_ylabel(r'$\theta(^\circ)$',fontsize=20)
        for k in range(m):
            ax.plot(M[k,:],erreur[k,:,0],'-',label='exa., M0=%4.2f'%(M0[k]))
            ax.plot(M[k,:],erreur[k,:,1],'-',label='app., M0=%4.2f'%(M0[k]))
            ax.legend(loc='lower left')
    elif cas==1:
        
        Mach=np.linspace(M_init,M_end,n)
        fig.suptitle(r'$\theta(M), M0 $= %3.1f'%(M_init), fontsize=14, fontweight='bold')
        y1,y2=approx_omega(Mach,M_init)
        #ax.axis([M_init, M_end, -20, 0])
        ax.set_ylabel(r'$\theta(^\circ)$',fontsize=20)
        ax.plot(Mach,y1,label='Exact')
        ax.plot(Mach,y2,label='Approx')
        ax.legend(loc='upper right')
    else :
        Mach=np.linspace(M_init,M_end,n)
        fig.suptitle('erreurs M0 = %3.1f'%(M_init), fontsize=14, fontweight='bold')
        #ax.axis([M_init, M_end, 100,140])
        #ax.plot(Mach,np.rad2deg(omega_super(Mach)),label=r'$\omega$ exact')
        #ax.plot(Mach,np.rad2deg(omega_hyper(Mach)),label=r'$\omega$ approx.')
        ax.plot(Mach,(omega_hyper(Mach)-np.rad2deg(omega_super(Mach))),label=r'$\Delta\omega(^\circ)$')
        ax.legend(loc='lower right')
   
    plt.show() 