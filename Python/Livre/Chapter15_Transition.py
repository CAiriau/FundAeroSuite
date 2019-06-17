#!/bin/python
"""
  Correction des exercices du chapitre 15 : transition
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

#************************************
# main program
#************************************

def Michel(self):
    """
    Critère de Transition de Michel
    """
    # self : Re_x
    return 1.535*self**(0.444)

def Cebeci_Smith(self):
    """
    Critère de transition de Cebeci-Smith
    """
    # self : Re_x
    return 1.174*(1+22400./self)*self**0.46

def Re_crit(self):
    """
    Nombre de Reynolds Critique
    """  
    # self= H
    if (self>=2.591):
        Rec=54.2124/self/(self-2.48)+31.6/self
    else:
        #Rec=520./self+2.5e6/self*(1./self-1./2.591)**1.95
        Rec=520./self+2.5e6/self*(1./self-1./2.591)**1.95
    return Rec  

 
def Granville(m,H,k2):
    """
    Critère de Transition Granville
    """
    Lambda=m*k2**2
    return Re_crit(H)+375.0+np.exp(6.1+55*Lambda)

def Arnal(m,H,k2,Tu):
    """
    Critère de Transition d'Arnal
    """
    Lambda=m*k2**2
    return Re_crit(H)-206*np.exp(25.7*Lambda)*(np.log(16.8*Tu)-2.77*Lambda)

def FSK_correlation(self):
    """
    Corrélation valable pour   -0.14 < beta < 0.11
    pour une couche limite de Falkner-Skan
    self=beta
    Calculs fait par C. Airiau par une méthode de continuation.
    """
    Lambda= 4.5157e-5 + 0.2196 *self - 0.3872 *self**2 + 0.7967*self**3
    H     = 2.5917 - 1.404 *self + 3.885 *self**2 - 20.516 *self**3 + 104.06*self**4
    k_2   = 0.664 - 0.71153 *self + 0.9431 *self**2 - 1.9697 *self**3
    return k_2,H,Lambda


def print_correlations():
    """
    affichage des valeurs de référence pour la couche limite de Falkner-Skan
    """
    set_title("Corrélation pour la couche limite de Falner-Skan")
    n=51
    beta=np.linspace(-0.14,0.11,n)
    m=beta/(2-beta)
    print('beta \t\t  m \t\t k_2 \t\t H \t\t Lambda')
    for b in beta:
        k_2,H,Lambda=FSK_correlation(b)
        m_tmp=b/(2.-b)
        print('%2.5f \t'*5%(b,m_tmp,k_2,H,Lambda))
    return

def plot_correlations():
    """
    Calcul et tracé des corrélations pour FSK
    """
    # print_correlations()
    # lecture des fichiers de FSK
    reference_filepath = os.path.join('Livre/Data', 'BL_FSK_characteristics.dat')
    with open(reference_filepath, 'r') as infile:
        A=np.loadtxt(infile, dtype=float, unpack=True,skiprows=1)
    #     beta, m, H, s, d1/d, d2/d, Cf, xk, Eta Max, Lambda, End_BL
    print(A.shape)  
    beta,m,H,k_2,Lambda=A[[0,1,2,5,9],:]
    k2_coef=np.polyfit(beta,k_2,3); p_k2=np.poly1d(k2_coef)
    H_coef=np.polyfit(beta,H,4);p_H=np.poly1d(H_coef) 
    Lambda_coef=np.polyfit(beta,Lambda,3);p_Lambda=np.poly1d(Lambda_coef)
    eps=1e-10
    print('polynome k2 :', k2_coef)
    print('polynome H :', H_coef)
    print('polynome Lambda :',Lambda_coef)
    k2_1,H_1,Lambda_1=FSK_correlation(beta)   
    err_k2=(k2_1-p_k2(beta))/(p_k2(beta)+eps) 
    err_H=(H_1-p_H(beta))/(p_H(beta)+eps) 
    err_Lambda=(Lambda_1-p_Lambda(beta))/(p_Lambda(beta)+eps)   
   
    fig = plt.figure()
    fig.suptitle('erreurs sur les corrélations', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_title('corrélations')
    ax.set_xlabel(r'$\beta$',fontsize=20)
    ax.grid()
    ax.semilogy(beta,abs(err_k2),'-',label=r'k_2',linewidth=2,color='black')
    ax.semilogy(beta,abs(err_H),'-',label='H',linewidth=2,color='blue')
    ax.semilogy(beta,abs(err_Lambda),'-',label=r'\Lambda',linewidth=2,color='red')
    ax.legend(loc='upper left')
    plt.show() 
        
    m_test=0.
    beta_test=2.*m_test/(m_test+1.)
    print('beta_test = %1.3f \t m_test = %1.3f'%(beta_test,m_test))
    print('k2 = %1.4f, \t, H = %1.4f, \t, Lambda = %1.4f, %1.4f'%(p_k2(beta_test),p_H(beta_test),p_Lambda(beta_test),m_test*p_k2(beta_test)**2))
    return


         
def Exo15_Tra_1():
    """ 
    Comparaison des critères de transition 
    """

    plot=True
    flag_correlation=True
    flag_transition=True
    Tu=0.001        # Taux de turbulence
    n=51            # nombre de points pour Re_x
    k=4             # indice du tableau pour m

    if (flag_correlation) :plot_correlations()
       
    set_title("Comparaison des critères de transition")
    # données de Falkner-Skan 
    m=np.array([1.0,1./3.,1./7.,0.0,-1./21.,-0.0904])
    beta=2*m/(m+1)
    k2=np.array([0.2923, 0.4290, 0.5245, 0.6641, 0.7464, 0.8681])  # Rdelta2/Rdelta
    H=np.array([2.2162, 2.2970, 2.3843, 2.5911, 2.8011, 4.0292])
    
    if (k==4) :
        Rex=(np.linspace(600,800,n))**2    # m=-1/21
    else:   
        Rex=(np.linspace(1400,1900,n))**2  # m=0
    
    Id=np.ones(n)
   
    meth=['Michel', 'Cebeci-Smith', 'Granville','Arnal']
    print(' Valeur de m = %1.3f, beta = %1.3f \t k2= %1.4f , \t Lambda= m k2^2 = %1.4f'%(m[k],beta[k],k2[k], m[k]*k2[k]**2))
    print('Re delta  crit = %f, \t Re delta2  crit = %f'%(Re_crit(H[k])/k2[k],Re_crit(H[k])))
    print('Lambda = ',m[k]*k2[k]**2,' H = ',H[k])


    if plot:
        fig = plt.figure()
        fig.suptitle('Critère de transition', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_title('critère de Michel')
        ax.set_xlabel(r'$Re_x$',fontsize=20)
        #ax.set_ylabel(r'$\frac{P_i}{P_t}-1$',fontsize=20)
        ax.grid()
        #ax.axis([0., 1, 0, 20])
        ax.plot(np.sqrt(Rex),k2[k]*np.sqrt(Rex),'-',label='Falkner-Skan',linewidth=2,color='black')  # r'$Re_{\delta_2} $'
        ax.plot(np.sqrt(Rex),Michel(Rex),'-',label='Michel')
        ax.plot(np.sqrt(Rex),Cebeci_Smith(Rex),'-',label='Cebeci-Smith')
        ax.plot(np.sqrt(Rex),Id*Granville(m[k],H[k],k2[k]),'-',label='Granville')
        ax.plot(np.sqrt(Rex),Id*Arnal(m[k],H[k],k2[k],Tu),'-',label='Arnal')
        ax.legend(loc='upper left')
        plt.show() 

    Re=np.zeros((4,n))
    Rdeltk2=k2[k]*np.sqrt(Rex)
    Re[0,:]=Michel(Rex)-Rdeltk2
    Re[1,:]=Cebeci_Smith(Rex)- Rdeltk2
    Re[2,:]=Id*Granville(m[k],H[k],k2[k])- Rdeltk2
    Re[3,:]=Id*Arnal(m[k],H[k],k2[k],Tu)-Rdeltk2 

    # on a les solutions analytiques pour le critère de Granville et Arnal pour le cas de FSK.
    # je cherche les zéros de ces fonctions par interpolation linéaire 
    for i in np.arange(4):
        f1 = interp1d( Re[i,:],np.sqrt(Rex), kind='cubic')
        print('Re delta= %5.1f ,\t Méthode = %s, '%(f1(0),meth[i]))

    print("Critère de Granville : solution analytique : ")
    print(Granville(m[k],H[k],k2[k])/k2)
    print("Critère de Arnal     : solution analytique : ")
    print(Arnal(m[k],H[k],k2[k],Tu)/k2)
