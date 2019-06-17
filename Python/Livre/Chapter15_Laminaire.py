#!/bin/python
"""
  Correction des exercices du chapitre 15 : laminaire
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


 

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize             import fsolve
from Tools.misc                 import *
from CompressibleFlow.fonctions import * 
from scipy.interpolate          import interp1d
from scipy.special              import erfc
from scipy.special              import erf



 

def Ecoulement_Blasius():
    """
    Profil de Blasius, calcul fait en Fortran avec une très grande précision
    """
    # Solution de Blasius
    reference_filepath = os.path.join('Livre/Data', 'FSK_beta_0.dat')
    with open(reference_filepath, 'r') as infile:
        A=np.loadtxt(infile, dtype=float, usecols=(0,2),unpack=True,skiprows=2)
    eta_Blasius,U_Blasius=A[0,:],A[1,:]
    return eta_Blasius/eta_Blasius[-1],U_Blasius

def BL_caracteristics(eta,U):
    """
    caractéristiques de la couche limite
    """
    delta1,delta2=np.trapz(1-U,x=eta),np.trapz(U*(1-U),x=eta)
    return [delta1/delta2,delta1,delta2]

def Exo15_Lam_1():
    """ 
    Plaque plane : Loi en sinus et polynome
    cet correction va plus loin que le sujet car on va jusqu'à des polynomes de degrés 7.
    """

    set_title("Plaque plane : Loi polynomiale et sinus")
    plot=True
    zoom=False
    ref=False               # on dessine U-U_blasius
    n=1001                   # on calcule les grandeurs par la méthode des trapèzes, précision requise
    Ncas=14
    option_poly7=True       # pour avoir beaucoup plus de cas
    # les méthodes possibles :
    meth=['Blasius','Sinus','Poly 3','Poly 4','Poly 5 ','Poly 6 ','Poly 7','Poly 5H- 1','Poly 5H- 2',
    'Poly 6 +d1+d2','Poly 7 +d1+d2','Poly 8 +d1+d2','Poly 9 +d1+d2','Poly 9 + lsq']
    eta=np.linspace(0,1,n)
    CL,U=np.zeros([Ncas,3]),np.zeros([Ncas,n])
    
    set_question('1 et 2 ')
    # Solution de Blasius
    eta_Blasius,U_Blasius=Ecoulement_Blasius()
    CL[0,:]=BL_caracteristics(eta_Blasius,U_Blasius)

    f1 = interp1d(eta_Blasius,U_Blasius, kind='cubic')
    U[0,:]=f1(eta)
    # Solution en sinus
    k=1; U[k,:]= np.sin(np.pi/2*eta)
    CL[k,:]=BL_caracteristics(eta,U[k,:])
    # Solution Polynomiale de degré 3
    k=2; U[k,:]=3/2*eta-eta**3/2
    CL[k,:]=BL_caracteristics(eta,U[k,:])
    # Solution Polynomiale de degré 4
    k=3; U[k,:]=2*eta-2*eta**3+eta**4
    CL[k,:]=BL_caracteristics(eta,U[k,:])
  
    if option_poly7:
        # Solution Polynomiale de degré 5
        k=4;U[k,:]=5/2*eta-5*eta**3+5*eta**4-3/2*eta**5
        CL[k,:]=BL_caracteristics(eta,U[k,:])
        # Solution Polynomiale de degré 6
        k=5;U[k,:]=3*eta-10*eta**3+15*eta**4-9*eta**5+2*eta**6
        CL[k,:]=BL_caracteristics(eta,U[k,:])
        # Solution Polynomiale de degré 7   
        k=6;U[k,:]=7/2*eta-35/2*eta**3+35*eta**4-63/2*eta**5+14*eta**6-5/2 *eta**7
        CL[k,:]=BL_caracteristics(eta,U[k,:])
        
        # k= 7 et 8 sont réservés à la stratégie sur H, plus bas 

        # pour la stratégie suivante les calculs sont faits sous MAPLE 
        # stratégie sur delta_1 et delta_2, polynôme de degré 6   d1 := .21357; d2 := 0.8242e-1
        k=9;U[k,:] = -.972245073*eta**6+.5321734494*eta**5+4.442096167*eta**4-5.796878467*eta**3+2.794853923*eta
        CL[k,:]=BL_caracteristics(eta,U[k,:]) 
        # U2 := -.9722450731*eta^6+.5321734494*eta^5+4.442096167*eta^4-5.796878467*eta^3+2.794853923*eta
        #delta_1 := .2135700002  delta_2 := 0.8242000023e-1 H := 2.591239985
        # stratégie sur delta_1 et delta_2, polynôme de degré 7
        k=10;U[k,:]= 5.528385993*eta**7-21.17914117*eta**6+28.16630546*eta**5-12.44152997*eta**4-1.807621010*eta**3+2.733600700*eta
        CL[k,:]=BL_caracteristics(eta,U[k,:] )
        #U1 := 5.528385993*eta^7-21.17914117*eta^6+28.16630546*eta^5-12.44152997*eta^4-1.807621010*eta^3+2.733600700*eta
        # delta_1 := .2135699998 delta_2 := 0.8241999973e-1 H := 2.591239996
        # stratégie sur delta_1 et delta_2, polynôme de degré 8
        k=11;U[k,:] = -9.553848646*eta**8+49.47167446*eta**7-101.7101564*eta**6+101.8603675*eta**5-46.38869308*eta**4+4.661142338*eta**3+2.659513754*eta
        CL[k,:]=BL_caracteristics(eta,U[k,:])
        #U2 := -9.553848646*eta^8+49.47167446*eta^7-101.7101564*eta^6+101.8603675*eta^5-46.38869308*eta^4+4.661142338*eta^3+2.659513754*eta
        #delta_1 := .2135700116  delta_2 := 0.8242001159e-1 H := 2.591239767
        # stratégie sur delta_1 et delta_2, polynôme de degré 9
        k=12;U[k,:]= 9.378936756*eta**9-61.63551826*eta**8+169.4617881*eta**7-248.6453033*eta**6+203.0172616*eta**5-83.82909404*eta**4+10.64557881*eta**3+2.606350378*eta
        CL[k,:]=BL_caracteristics(eta,U[k,:])
        # U1 := 9.378936756*eta^9-61.63551826*eta^8+169.4617881*eta^7-248.6453033*eta^6+203.0172616*eta^5-83.82909404*eta^4+10.64557881*eta^3+2.606350378*eta
        # delta_1 := .2135699938  delta_2 := 0.8241999380e-1 H := 2.591240110
        C = np.polyfit(eta_Blasius,U_Blasius, 6)

        #le fitting ne marche pas, pas 1 en eta = 1
        print("fitting O9 : ",C)
        p9=np.poly1d(C)
        k=13;U[k,:]=p9(eta) 
        CL[k,:]=BL_caracteristics(eta,U[k,:])


    H=2.5911
    H=CL[0,0]

    print("Stragégie sur le facteur H  = %f"%(H))
    s=np.poly1d([104*H,2079-528*H,-7326*H+18711])
    print(s.r)

    print("Stragégie sur delta_1 et delta_2 : ")
    
    delta_1,delta_2=CL[0,1],CL[0,2]
    print('delta_1 = %f, \t delta_2 = %f'%(delta_1,delta_2))

    # Polynome de degré 7.    ce n'est pas absolument nécessaire, mais c'est fait ....
    print('Solution pour un polynome de degré 7')
    s7=np.poly1d([3,720*delta_1+360,862400*delta_1**2-582000*delta_1+514800*delta_2+39600])
    s7_sol=max(s7.r)
    print(s7.r,' la bonne  solution est ',s7_sol) 
    a0,a2 = 0,0 
    a1 = 6-(1/20)*s7_sol-14*delta_1
    a3 = -40+(3/2)*s7_sol+140*delta_1
    a4 = 75-5*s7_sol-280*delta_1
    a5 = -54+(27/4)*s7_sol+210*delta_1
    a6 = 14-(21/5)*s7_sol-56*delta_1
    a7 = s7_sol
    coef=[a7,a6,a5,a4,a3,a2,a1,a0]
    print("coefficients du polynôme : ",coef)
    p7=np.poly1d(coef)
    k=10
    U[k,:]=p7(eta)
    CL[k,:]=BL_caracteristics(eta,U[k,:])

    # Polynome de degré 9,  juste pour la forme
    s1=np.poly1d([62,55860*delta_1+16625, 157525200*delta_1**2-85186500*delta_1+77156625*delta_2+4375700])
    s1_sol=max(s1.r)
    print(s1.r,' la bonne  solution est ',s1_sol)   
    a1 = -(1/35)*s1_sol+8-24*delta_1
    a0,a2 = 0,0
    a3 = (8/5)*s1_sol-112+504*delta_1
    a4 = -8*s1_sol+350-1680*delta_1
    a5 = 18*s1_sol-504+2520*delta_1
    a6 = -(112/5)*s1_sol+392-2016*delta_1
    a7 = 16*s1_sol-160+840*delta_1
    a8 = -(216/35)*s1_sol+27-144*delta_1
    a9 = s1_sol
    coef=[a9,a8,a7,a6,a5,a4,a3,a2,a1,a0]
    print("coefficients du polynôme : ",coef)
    p9=np.poly1d(coef)
    k=12
    U[k,:]=p9(eta)
    CL[k,:]=BL_caracteristics(eta,U[k,:])

    # suivant qu'on dessine les profils ou l'erreur
    if ref:
        Uref=f1(eta)
    else:
        Uref=0

    if plot:
        if zoom:
            figsize=(6,6)
        else:
            figsize=(10,8)
        fig = plt.figure(figsize=figsize)
        fig.suptitle('Profil de couche limite', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('corrélations')
        ax.set_ylabel(r'$\eta$',fontsize=20)
        ax.set_xlabel(r'$\frac{u}{U_e}$',fontsize=20)
        ax.grid()
        
        ax.plot(U[2,:]-Uref,eta,'-',label='Polynôme de degré 3',linewidth=2,color='red')
        ax.plot(U[3,:]-Uref,eta,'-',label='Polynôme de degré 4',linewidth=2,color='green')
        ax.plot(U[1,:]-Uref,eta,'-',label='Sinus',linewidth=2,color='blue')
        if option_poly7:
            ax.plot(U[4,:]-Uref,eta,'--',label='Polynôme de degré 5',linewidth=2,color='red')
            ax.plot(U[5,:]-Uref,eta,'-.',label='Polynôme de degré 6',linewidth=2,color='black')
            ax.plot(U[6,:]-Uref,eta,':',label='Polynôme de degré 7',linewidth=2,color='blue')
            ax.plot(U[9,:]-Uref,eta,'--',label='Polynôme de degré 6 + d1+d2',linewidth=1,color='red')
            ax.plot(U[10,:]-Uref,eta,'-',label='Polynôme de degré 7 + d1+d2',linewidth=3,color='violet')
            ax.plot(U[11,:]-Uref,eta,':',label='Polynôme de degré 8 + d1+d2',linewidth=3,color='red')
            ax.plot(U[12,:]-Uref,eta,'--',label='Polynôme de degré 9 + d1+d2',linewidth=1,color='black')
            ax.plot(U[13,:]-Uref,eta,'--',label='Polynôme de degré 9 + lsq',linewidth=3,color='green')
        if ref:
            ax.plot(U[0,:]-Uref,eta,'-o',markersize=7,markevery=(0,50),label=r'Blasius',linewidth=1,color='black')
        else:
            ax.plot(U_Blasius,eta_Blasius/eta_Blasius[-1],'-o',markersize=7,markevery=(0,50),label=r'Blasius',linewidth=1,color='black')
        k=7
        for s in s.r:
            a0,a1,a2,a3,a4,a5 =0, -(1/3)*s+2,0,2* s-2,-(8/3)* s+1, s
            print('Poly 5-H',a0,a1,a2,a3,a4,a5)
            p=np.poly1d([a5,a4,a3,a2,a1,a0])
            U[k,:]=p(eta)
            CL[k,:]=BL_caracteristics(eta,U[k,:])
            print('s  = %f, H =%f  '%(s,CL[k,0]))
            ax.plot(U[k,:]-Uref,eta,'--',label='p5%f'%s,linewidth=1)
            k+=1
        ax.legend(loc='upper left')
        if ref==False:
            if zoom:
                ax.axis([0.98,1.02,0.6,1])
            else:
                ax.axis([0.,1.01,0.,1])


        plt.show() 

    if option_poly7==False:
        Ncas=4    
    
    print(" cas \t nom  \t\t\t H \t\t delta_1 \t delta_2     e(H) \t  e(delta_1) \t  e(delta_2) \t   norm(e)  \t norme (U-U_Blasius)")
    for k in range(Ncas): 
        erreur=np.abs(CL[k,:]/CL[0,:]-1)
        #print(erreur)
        print('%2i \t %15s   \t %7.5f  \t %7.5f  \t %7.5f  | %7.5f \t  %7.5f \t %7.5f \t %7.5f \t %7.5f '%(k,meth[k],CL[k,0],CL[k,1],CL[k,2],
            100*erreur[0],100*erreur[1],100*erreur[2],100*np.linalg.norm(erreur),100*np.linalg.norm(U[k,1:]-U[0,1:]) ))
         
   

def Exo15_Lam_2():
    """ 
    Plaque plane : Loi polynomiale
    """

    set_title("Plaque plane : Loi polynomiale")
    plot=True
    n=201   
    
    set_question('1 : Calcul de delta1 , delta2 et H, polynôme de degré 4')
    delta_1,delta_2=3/10,37/315
    print("Delta_1              : %5.3f"%(delta_1))
    print("Delta_2              : %5.3f"%(delta_2))
    print("H                    : %5.3f"%(delta_1/delta_2))
  
   
    set_question('5: plots')
    # Solution de Blasius
    eta_Blasius,U_Blasius=Ecoulement_Blasius()

    #Solution Polynomiale de degré 4 et 3
    eta=np.linspace(0,1,n)
    U_4=2*eta-2*eta**3+eta**4
    U_3=3/2*eta-eta**3/2
    
    if plot:

        fig = plt.figure()
        fig.suptitle('Profil de couche limite', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('corrélations')
        ax.set_xlabel(r'$\eta$',fontsize=20)
        ax.set_xlabel(r'$\frac{u}{U_e}$',fontsize=20)
        ax.grid()
        ax.plot(eta_Blasius/eta_Blasius[-1],U_Blasius,'-',label=r'Blasius',linewidth=2,color='black')
        ax.plot(eta,U_3,'-',label='Polyn\^ome de degr\'e 3',linewidth=2,color='red')
        ax.plot(eta,U_4,'-',label='Polyn\^ome de degr\'e 4',linewidth=2,color='blue')
         
        ax.legend(loc='lower right')
        plt.show() 

def Exo15_Lam_3():
    """
    couche limite avec aspiration constante
    """
    print('pas de correction numérique')


def Exo15_Lam_4():
    """ 
    Mise en mouvement d'une plaque, problème de Rayleigh
    """

    set_title("Plaque plane : Loi polynomiale")
    plot=True
    n=601
    eta=np.linspace(0,2,n)
    u=erfc(eta)

    def f(x):
        return erfc(x)-0.01
    x=fsolve(f, 2)
    print(x)
    x=2
    delta_1=(x*np.sqrt(np.pi)*erf(x)+np.exp(-x**2)-1)/np.sqrt(np.pi)
    delta_2= (-erf(x)**2*x*np.sqrt(np.pi)+x*np.sqrt(np.pi)*erf(x)+np.sqrt(2)*erf(np.sqrt(2)*x)-2*erf(x)*np.exp(-x**2)+np.exp(-x**2)-1)/np.sqrt(np.pi)
    print('solution analytique : H = %5.3f \t, delta_1 =  %5.3f, \t delta_2 =  %5.3f'%(delta_1/delta_2,delta_1,delta_2))
    H,delta_1,delta_2=BL_caracteristics(eta,f(eta))
    print('solution numerique  : H = %5.3f \t, delta_1 =  %5.3f, \t delta_2 =  %5.3f'%(delta_1/delta_2,delta_1,delta_2))

    if plot:

        fig = plt.figure()
        fig.suptitle('Profil de couche limite', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        #ax.set_title('corrélations')
        ax.set_ylabel(r'$\eta$',fontsize=20)
        ax.set_xlabel(r'$\frac{u}{U_0}$',fontsize=20)
        ax.grid()
        ax.plot(u,eta,'-',label=r'erfc',linewidth=2,color='black')      
        ax.legend(loc='upper right')
        plt.show() 

def Couette_instationnaire(Y,T=2,k=1,SingleMode=True):
    """
    Solution du profil sous forme adimentionnelle, mode k
    """  
    if SingleMode:
        return 2/np.pi*np.cos(k*np.pi)*np.exp(-k**2*np.pi**2*T)*np.sin(k*np.pi*Y)/k
    else:
        U=0
        for m in np.arange(1,k+1):
            amp=2/np.pi*np.cos(m*np.pi)*np.exp(-m**2*np.pi**2*T)*np.sin(m*np.pi*Y)/m
            U+=amp
            print('log amp du mode %2i =  %4.2f'%(m,np.log10(max(np.abs(amp/Y)))))
        return U
   


def Exo15_Lam_6():
    """ 
    Ecoulement de Couette instationnaire
    """

    set_title("Ecoulement de Couette instationnaire")
    plot=True
    ny=101
    T=np.array([0.01,0.02,0.05, 0.1, 0.2,0.3,0.5])
    nt=len(T)
    Y=np.linspace(0.0001,1,ny)

    # critère pour conserver de la précision pour les temps faibles
    limite=-np.log(0.001*np.abs(np.sin(np.pi*Y[0]))*np.pi/2)/np.pi**2
    #limite=np.log(2)/np.pi**2
    print('n^2 T = ',limite, ' n pour le premier temps : ',limite/T[0])
  
 
    if plot:
        fig = plt.figure(figsize=(10,8))
        fig.suptitle('Profil de couche limite', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        #ax.set_title('corrélations')
        ax.set_ylabel(r'$Y$',fontsize=20)
        ax.set_xlabel(r'$U$',fontsize=20)
        ax.grid()
        #for k in range(1):
        Nmode=50
        #for k in range(1):
        time=0.1
        for k in range(nt):
            Nmode=np.ceil(np.sqrt(limite/T[k]))
            print(' T = %5.3f, n minimal = %2i'%(T[k],Nmode))
            U=Couette_instationnaire(Y,T[k],Nmode,SingleMode=False)
            ax.plot(Y+U,Y,'-',label=r't=%f'%(T[k]),linewidth=2)      
        ax.axis([-0.05,1.05,0,1.05])
        ax.legend(loc='lower right')
        plt.show() 
    
def U(T,Y):
	"""
	Profil de vitesse
	"""
	eta=Y/np.sqrt(2)
	return -np.sin(T-eta) *np.exp(-eta)+np.sin(T)  

def Exo15_Lam_7():
    """ 
    Ecoulement instationnaire : second problème de Stokes
    """

    set_title("Second problème de Stokes")
    plot = True
    ny   = 101 		# nombre de points en Y
    ymax = 8        # limite supérieure du domaine
    nT   = 5        # nombre de points dans le temps
    epsilon=0.001
    print('Ymax =', -np.sqrt(2)*np.log(epsilon))
    print("erreur à l'infini   : ",np.exp(-ymax/np.sqrt(2)))

    # T=np.pi/2* np.linspace(-1,1,nT)
    T=np.pi/2*np.array([-1,-3/4,-1/2,-1/4, -1/8,0.,1/4,1/2,3/4,1])
    nt=len(T)
    Y=np.linspace(0,ymax,ny)

 	 
 
    if plot:
        fig = plt.figure(figsize=(11,8))
        fig.suptitle('Profil de couche limite, 2nd problème de Stokes', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_ylabel(r'$Y$',fontsize=20)
        ax.set_xlabel(r'$U$',fontsize=20)
        ax.grid()
       
        for t in T :
            ax.plot(U(t,Y),Y,'-',label=r't=%f'%(t*2/np.pi),linewidth=2)      
        #ax.axis([-0.05,1.05,0,1.05])
        ax.legend(loc='lower right')
        plt.show() 
    
def Exo15_Lam_8():
    """ 
    Ecoulement stationnaire : jet laminaire
    la solution est analytique
    """
    n=501
    xmax,ymax=5,1
    eta=np.linspace(-xmax,xmax,n)
    u=1/np.cosh(eta)**2
    v=2*eta*u-np.tanh(eta)
    fig,ax1=plt.subplots(figsize=(11,8))
    ax1.set_ylabel(r'$u$',fontsize=20)
    ax1.set_xlabel(r'$\eta$',fontsize=20)
    line1=ax1.plot(eta,u,'k-',linewidth=2,label='u')
    #plt.grid()
    ax2=ax1.twinx()
    ax2.set_ylabel(r'$v$',fontsize=20)
    line2=ax2.plot(eta,v,'r--',linewidth=2,label='v')
    #plt.legend(loc='lower right')

    xmajor_ticks,xminor_ticks = np.arange(-xmax, xmax+1, 1),np.arange(0, xmax, 0.5)                                              
    ymajor_ticks,yminor_ticks = np.arange(0,ymax+0.2, 0.2)   ,np.arange(0,ymax, 0.1)                                              
    ax1.set_xticks(xmajor_ticks);ax1.set_xticks(xminor_ticks, minor=True)  
    ymajor_ticks2,yminor_ticks2 = np.arange(-ymax,ymax+0.4, 0.4)   ,np.arange(-ymax,ymax, 0.1) 
    ax1.set_yticks(ymajor_ticks);ax1.set_yticks(yminor_ticks, minor=True)  
    ax2.set_yticks(ymajor_ticks2);ax2.set_yticks(yminor_ticks2, minor=True)  
    ax1.grid(b=True,which='major',color='#666666',linestyle='-');
    ax1.minorticks_on()
    ax1.grid(b=True,which='minor',color='#666666',linestyle='-',alpha=0.2);

    
    #ax1.legend();ax2.legend()
    first_legend=plt.legend( handles=line1,loc=(0.9,0.95))
    plt.gca().add_artist(first_legend)
    second_legend=plt.legend(handles=line2,loc=(0.9,0.90))
    ax=plt.gca().add_artist(second_legend)
    fig.tight_layout()
    plt.show() 