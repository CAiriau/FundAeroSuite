#!/bin/python
"""
    Correction des exercices du chapitre 9
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


 


import numpy as np
import matplotlib.pyplot as plt

from CompressibleFlow.fonctions import *
from Tools.misc                 import * 

#************************************
# main program
#************************************


def Exercice9_1():
    """ 
    Théorie linéaire : plaque plane en incidence
    """
    set_title("Plaque plane en incidence par la théorie linéarisée")

    # data, parameters
    alpha_deg=np.array([5, 10, 15 ])   # table of angle of attack
    M0=np.float64(2.);
    set_question('1 : Coefficients aérodynamiques')
    beta=np.sqrt(M0**2-1);
    alpha=alpha_deg*np.pi/180.;   
    CL=4/beta*alpha;
    CD=beta/4.*CL**2;
    f=CL/CD;

    print('# Mach  = %f '%(M0))
    print('# Alpha (°)\t CL \t\t CD \t\t f')
    for i in range(len(alpha_deg)):
        print(' %3.1f \t\t %3.5f \t %3.5f \t  %3.3f '%(alpha_deg[i],CL[i],CD[i],f[i]))


def CD_linear(CL,CD0,beta):
    """
    Polar in linearized theory
    """   
    return CD0+beta*CL**2/4     

def Exercice9_2():
    """ 
    Théorie linéaire :  profil losangique
    """
    set_title("Profil losangique par la théorie linéarisée")

    # data, parameters


    e_sur_L=np.float64(0.05);       # airfoil aspect ration
    alpha_deg=np.float64(3);          # Angle of attack  in degrees
    M0=np.float(1.4);               # Mach number
    
    mu0=np.arcsin(1/M0);  
    beta=np.sqrt(M0**2-1);
    print('Mach               = %f '%(M0));
    print('beta               = %f '%(beta));
    print('mu0                = %f '%(mu0/np.pi*180));
    theta=np.arctan(e_sur_L)
    print('theta in degrees   = %f '%(theta/np.pi*180));
    print('alpha in degrees   = %f '%(alpha_deg));
    alpha=alpha_deg/180*np.pi;

    set_question('1 : Kp')

    # take care of index (zone -1 ) of the exercice
    Kp=np.zeros((4))
    Kp[0]=2/beta*(theta-alpha);
    Kp[1]=-2/beta*(theta+alpha);
    Kp[2]=-Kp[1];
    Kp[3]=-Kp[0];
    np.set_printoptions(formatter={'float': '{: 2.3f}'.format})
    print(Kp, '\t Kp of the 4 zones');

    set_question('2 : Coefficients aérodynamiques');

    CL=4*alpha/beta;
    CD=4*(theta**2+alpha**2)/beta;
    CmA=2*alpha/beta;
    print('CL = %f, \t CD = %f CmA = %f'%(CL, CD,CmA));

    set_question('3: ');

    hstar=np.tan(mu0)/2;
    print('h*/L0                   = %f'%(hstar));

    set_question('4: ');

    alpha_deg=np.array(0);alpha=alpha_deg/180*np.pi;
    DeltaCD=-2*theta**2/beta;
    DeltaCL=2*theta/beta;
    print('alpha =  %f, \tDelta CL = %f, \t Delta CD = %f'%(alpha_deg, DeltaCL, DeltaCD));
    alpha_deg=np.array(3);alpha=alpha_deg/180*np.pi;
    DeltaCD=2*(alpha**2-theta**2)/beta;
    DeltaCL=2*(alpha+theta)/beta;
    print('alpha =  %f, \tDelta CL = %f, \t Delta CD = %f'%(alpha_deg, DeltaCL, DeltaCD));

def Exercice9_3():
    """ 
    Théorie linéaire : profil de traînée minimale
    """
    set_title("Profil de traînée minimale par la théorie linéarisée")

    # data, parameters
    plot=True 

    M0=2;
    M0=np.float(1.4);               # Mach number
    beta=np.sqrt(M0**2-1);
    CD0=np.float64(0.01);

    # Optimal aerodynamic efficiency
    CLopt=2*np.sqrt(CD0/beta);
    CDopt=2*CD0;
    fmax=1/np.sqrt(beta*CD0);
    alpha_opt=beta/4*CLopt/np.pi*180;   # in degrees

    print('CL = %f, CD= %f, fmax = %f  alpha = %f'%(CLopt,CDopt,fmax,alpha_opt)); 
    
    if plot:
        # exemple de graphique
        CL=np.linspace(-0.3, 0.3, 31)
        CD=CD_linear(CL,CD0,beta)
        f=CL/CD
        CLref=np.array([ 0, CLopt, 2*CLopt])
        CDref=np.array([ 0, CDopt, 2*CDopt])
        fig = plt.figure(1)
        fig.suptitle('Polaire', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$C_L$',fontsize=20)
        ax.set_ylabel(r'$C_D$',fontsize=20)
        ax.plot(CL,CD,label=r'$Polar$',color='red',linewidth=2)
        ax.plot(CLref,CDref,'--',label=r'$Tangente$',color='blue',linewidth=2)
        ax.plot(CLopt,CDopt,'s',label='Optimun',color='magenta',linewidth=2,markersize=6)
        ax.grid()
        ax.legend(loc='lower left')
        plt.show()


        CL=np.linspace(0, 0.3, 31)
        CD=CD_linear(CL,CD0,beta)
        f=CL/CD
        fig = plt.figure(2)
        fig.suptitle('Finesse', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\alpha (^\circ)$',fontsize=20)
        ax.set_ylabel(r'$f$',fontsize=20)
        ax.plot(CL*beta/4*180/np.pi,f,label=r'$C_L/C_D$',color='red',linewidth=2)
        ax.plot(alpha_opt,fmax,'s',label='Optimal ',color='black',linewidth=2,markersize=6)
        ax.grid()
        ax.legend(loc='lower right')
        plt.show()


 
def Exercice9_4():
    """ 
    Théorie linéaire : profil  lenticulaire
    """
    set_title("Profil lenticulaire par la théorie linéarisée")

    # data, parameters
    print('il faut regarder le fichier maple ou pdf');
 

def Exercice9_5():
    """ 
    Théorie linéaire : jet débouchant
    """
    set_title("Interaction jet supersonique - profil")

    # data, parameters
    plot=True 

    set_question('1 et 2');
    
    eps=np.float64(0.05);
    eta=np.linspace(-1,1,1001);   
    y=-eps*(1+eta)*(eta**2-1);
    dy=eps*(1-2*eta-3*eta**2);
 
    set_question('4 et 5');

    M0=np.sqrt(2);
    beta=np.sqrt(M0**2-1);
    mu0=np.arcsin(1/M0)*180/np.pi;
    N=len(eta);
    print('longueur de eta    = %d '%(N));
    print('Mach               = %f '%(M0));
    print('beta               = %f '%(beta));
    print('mu0 (°)            = %f '%(mu0));
    print('demi epaisseur     = %f '%(32/27*eps));

    rpression=np.float64(0.9);
    
    Kp1=Pressure_coefficient(rpression,M0)
    print('Kp1                = %f'%(Kp1));
    u2=-Kp1/2;
    print('u2/U0              = %f'%(u2));
    delta1=beta*u2/np.pi*180;
    print('delta1(°)          = %f'%(delta1));

    phi=np.arctan(beta*u2)/np.pi*180
    print('phi(°)             = %f'%(phi));

    set_question('6');

    d1=delta1;
    k=0.5
    d=dy[int(np.ceil(k*N-1))]*180/np.pi;
    print('delta (°) à %f de la corde  = %f'%(k,d));
    v=np.array([d1 ,d1,  d ,d, -d1,  d-d1, d ])
    u=np.array([d1, d1, 2*d1-d, d, d1, d+d1, d+2*d1])
    


    if plot:
        
        fig = plt.figure(1)
        fig.suptitle('Forme du profil', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$\eta$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        ax.plot(eta,y,color='red',linewidth=2)
        ax.grid()
         
        #plt.show()

        fig = plt.figure(2)
        fig.suptitle('Pente du profil en °', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$\eta$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        ax.plot(eta,dy/np.pi*180,color='red',linewidth=2)
        ax.plot(np.array([-1,1]),np.array([0, 0]),'--',color='black', linewidth=1)
        ax.grid()
        #ax.axis([-1, 1, 0, 12])
        #plt.show()

        fig = plt.figure(3)
        fig.suptitle("Plan de l'hodographe", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$u$',fontsize=20)
        ax.set_ylabel(r'$v$',fontsize=20)
        ax.plot(u,v,'o--',color='red',linewidth=2)
        ax.grid()
        #ax.axis([0.0, 1, 0, 12])
        # plt.savefig('fichier.svg')
        plt.show()

 
