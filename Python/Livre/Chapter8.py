#!/bin/python
"""
  Correction des exercices du chapitre 8
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 




import numpy as np
import matplotlib.pyplot as plt

from Tools.misc                 import *
from CompressibleFlow.fonctions import *
from Livre.Wavy.wavy_wall       import * 

#************************************
# main program
#************************************


def Exercice8_1():
    """ 
    Total pressure / Stagnation pressure = Isentropic Pressure
    """
    set_title("Comparaison Pression d'arrêt / Pression totale")

    # data, parameters

    M=np.array([0.3, 0.5, 0.7, 0.9 ])   # table of Mach numbers
    plot=True                        # to plot, set True

    set_question('2 : table des valeurs')
    M=np.linspace(0.1,2,20);
    pt=pt_p(M)
    pi=1/p_pi(M)
    rate=pi/pt
    val=np.log10(rate-1)

    print('# Mach \t Pi/P \t  Pt/P \t Pi/Pt-1 in % \t log10 error')
    for i in range(len(M)):
        print(' %3.1f \t %3.2f \t %3.2f \t  %3.3f \t %3.3f'%(M[i],pi[i],pt[i],100*(rate[i]-1),val[i]))
     
    set_question("3: graphe de l'erreur relative")
    if plot :
        # exemple de graphique
         
        plt.figure(figsize=(12,8))
        plt.title(r'$P_i / P_t$', fontsize=14, fontweight='bold')
        plt.title('relative error :'+ r'$\log_{10}(P_i / P_t-1)$')
        plt.xlabel(r'$M$',fontsize=20)
        plt.ylabel(r'$\frac{P_i}{P_t}-1$',fontsize=20)
        plt.yscale('log')
        plt.grid()
        #plt.legend(loc='upper left')
        #ax.axis([0., 1, 0, 20])
        plt.plot(M,rate-1,'o')
        plt.show() 


def Exercice8_2():
    """ 
    Pitot tube in subsonic flow
    """
    set_title(' Tube de Pitot double en écoulement subsonique')

    # data, parameters

    a=np.float64(340.0)                 # speed of sound of reference
    M=np.array([0.3, 0.5, 0.7, 0.9 ])   # table of Mach numbers
    Delta_p=np.float64(10000)           # variation of pressure given by the measure
    p=np.float64(41100)                 # static pressure given by the measure
    T=np.float64(243)                   # static temperature given by the measure
    plot=False                          # to plot, set True

    V=a*M
    Pi=1/p_pi(M);
    Rcomp=Pressure_coefficient(Pi,M)
    rapU=np.sqrt(Rcomp)-1
    rapU1=M**2/8

    set_question('5: tableau avec les erreurs relatives en %')


    print('# V \t Mach \t Pi/P \t  Rcomp \t Um/U-1 \t approx(Um/U_1)')
    for i in range(len(M)):
        print(' %3.1f \t %1.1f \t %1.3f \t %1.3f \t \t %2.2f \t \t %2.2f'%(V[i],M[i],Pi[i],Rcomp[i],
            100*rapU[i],100*rapU1[i]))
     
     
    if plot :
        # exemple de graphique
         
        fig = plt.figure()
        fig.suptitle('Approximate error vs accurate error ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.set_ylabel(r'$\varepsilon$ en $\%$',fontsize=20)
        ax.plot(M,100*(rapU-1),label=r'$\varepsilon$',color='black',linewidth=2)
        ax.plot(M,100*rapU1,'o',label=r'$\varepsilon1$',color='red')
        ax.grid()
        ax.legend(loc='upper left')
        ax.axis([0.2, 1, 0, 12])

        plt.show()

    set_question('6: application tube de Pitot')

    a=sound_velocity(T)
    print('c               = %3.2f m/s      '%(a))
    Mach=np.sqrt(2/(gamma-1)*( (Delta_p/p+1)**((gamma-1)/gamma)-1)) 
    print('M               = %2.3f          '%(Mach))
    print('V               = %3.2f km/h     '%(Mach*a*3.6))


def Exercice8_3():
    """ 
    Compressibility corrections (airfoil)
    """
    set_title(' Corrections de compressibilité en écoulement subsonique sur un profil')

    # data, parameters
    plot=True   

    set_question('1: tracé des Kp pour les 3 corrections')   
                      # to plot, set True
    M_exp=np.array(  [0., 0.1, 0.2,     0.3, 0.52,    0.6, 0.64,   0.68 ])   # table of Mach numbers
    Kp_exp=-np.array([0.6, 0.602, 0.610, 0.628, 0.786, 0.809, 0.858, 0.9])
    M=np.array(  [0., 0.1, 0.2, 0.3, 0.4,0.5, 0.6, 0.7 ])
    #M=np.linspace(0,0.7,101)
    Kp_res=np.zeros((len(M),3))
    
    
    Kp_res=compressibility_corrections(Kp_exp[0],M)


    print('# M \t Mach ')
    np.set_printoptions(formatter={'float': '{: 2.3f}'.format})
    print(M,'\t M')
    print(60*'-')
    #print(Kp_exp)
    print(Kp_res[:,0],'\t PG')
    print(Kp_res[:,1],'\t KT')
    print(Kp_res[:,2],'\t L')
    
    #    print(' %3.1f \t %1.1f \t %1.3f \t %1.3f \t \t %2.2f \t \t %2.2f'%(V[i],M[i],Pi[i],Rcomp[i],
    #        100*rapU[i],100*rapU1[i]))
     
    

    set_question('2:  Kp critique')   
    M0=np.linspace(0.6,0.75,7)
    Kp_crit=Kp_critical(M0)
    omega=p_pi(M0)
    print(M0,'\t M')
    print(Kp_crit,'\t Kp_crit')

    if plot :
        # exemple de graphique
        Mc=np.array([0.674, 0.688])
        Kpc=np.array([-0.894, -0.829])
        fig = plt.figure(1)
        fig.suptitle('Critical Kp ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.set_ylabel(r'$Kp$',fontsize=20)
        ax.plot(M,Kp_res [:,0],label=r'$Kp_{PG}$',color='red',linewidth=2)
        ax.plot(M,Kp_res [:,1],label=r'$Kp_{KT}$',color='blue',linewidth=2)
        ax.plot(M,Kp_res [:,2],label=r'$Kp_L$',color='black',linewidth=2)
        ax.plot(M_exp,Kp_exp [:],'o-',label=r'$Kp_{exp}$',color='green',linewidth=2)
        ax.plot(M0,Kp_crit[:],label=r'$Kp_{crit}$ (analytique)',color='magenta',linewidth=2)
        ax.plot(Mc,Kpc,'s',label=r'$Kp_{crit}$ (graphique)',color='black',linewidth=3,
              markersize=6)
        ax.grid()
        ax.legend(loc='lower left')
        #ax.axis([0.0, 1, 0, 12])

        plt.show()

    set_question('3:  Coefficient de portance')   

    M_exp=np.array( [0.,0.22,0.32,0.4,0.42,0.5, 0.54, 0.592, 0.608,0.665,0.73])   # table of Mach numbers
    CL_exp=np.array([6 ,6.2, 6.31,6.52,6.6, 7,7.28, 7.66,7.9,7.82,6.22])/100.0
    M=np.linspace(0.,0.7,15)
    CL_comp=CL_compressible(M,CL_exp[0])
    print('\t M = \t', M)
    print('\t CL= \t',CL_comp)

    if plot :
        # exemple de graphique
        fig = plt.figure(1)
        fig.suptitle('Lift coefficient', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.set_ylabel(r'$100 \times C_L$',fontsize=20)
        ax.plot(M_exp,100*CL_exp,'o-',label=r'$CL_{exp}$',color='red',linewidth=2)
        ax.plot(M,100*CL_comp,label=r'$CL_{PG}$',color='blue',linewidth=2)
        ax.grid()
        ax.legend(loc='upper left')
        #ax.axis([0.0, 1, 0, 12])

        plt.show()

def f_critical(self,Kp_inc,Lambda):
    """
        function to solve to find the critical Mach number
    """
    return PG_correction_wing(Kp_inc,self,Lambda)-Kp_critical(self)

def Newton_Mcrit(Kp_inc,Lambda):
    """
        Newton method to obtain the critical Mach number
    """
    iterMax=15
    errorMax=np.float64(1e-6)
    Mach=0.7
    dM=np.float64(0.001)
    f=f_critical(Mach,Kp_inc,Lambda)
    error=1.0
    i=0
    while (error > errorMax ) and (i<iterMax):
        fp=f_critical(Mach+dM,Kp_inc,Lambda)
        deltaM=-f/(fp-f)*dM
        error=np.abs(deltaM)
        Mach=Mach+deltaM
        f=f_critical(Mach,Kp_inc,Lambda)
        i=i+1
        #print( '%d \t %f \t %e \t %e '%(i,Mach,error,f))

    if (error>errorMax):
        print('Pas de convergence')
    else:
        Kp_crit=Kp_critical(Mach)
        print('Solution pour Lambda = %2.0f , Kp_inc= %f, \t M_crit  = %1.3f pour Kp_crit = %2.3f '%(Lambda,Kp_inc,Mach,Kp_crit))

    return Mach,Kp_critical(Mach)


def search_and_plot_Mcrit(Kp_inc,plot):
    """ core of exercice 8-4 : critical Mach number on a swept wing
    """
    M=np.linspace(0.,0.85,81)
    Lambda=np.array([0,15,30,45])
    M0=np.linspace(0.6,0.85,51)
    Kp_crit=Kp_critical(M0)
    Kp=np.zeros((len(M),len(Lambda)))
    Mc=np.zeros(len(Lambda));Kpc=np.zeros(len(Lambda))
    for i in range(len(Lambda)):
       Kp[:,i]=PG_correction_wing(Kp_inc,M,Lambda[i])
       Mc[i],Kpc[i]=Newton_Mcrit(Kp_inc,Lambda[i])
       

    print(Mc)
    print(Kpc)   

    if plot :
        # exemple de graphique
        #Mc=np.array([0.688,0.698,0.728,0.776])
        #Kpc=np.array([-0.829,-0.785,-0.669,-0.507])
        fig = plt.figure(1)
        fig.suptitle('Critical Kp in swept wing ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        #ax.set_title('relative error :'+ r'$(U_m-U)/U$')
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.set_ylabel(r'$Kp$',fontsize=20)

        for i in range(len(Lambda)):
            ax.plot(M,Kp[:,i],label=r'$Kp_{%s}$'%(Lambda[i]),linewidth=2)
         

        ax.plot(M0,Kp_crit[:],label=r'$Kp_{crit}$ (analytique)',color='magenta',linewidth=2)
        ax.plot(Mc,Kpc,'s',label=r'$Kp_{crit}$ (Newton)',color='black',linewidth=3,
               markersize=6)
        ax.grid()
        ax.legend(loc='lower left')
        #ax.axis([0.0, 1, 0, 12])

        plt.show()

    return

def Exercice8_4():
    """ 
    Compressibility corrections (wing) and critical Mach number
    """
    set_title(' Corrections de compressibilité en écoulement subsonique sur une aile en flèche')

    # data, parameters
    print('Solution analytique :')
    M0=np.linspace(0.1,0.85,16);print('M= ',M0)
    Kp_crit=Kp_critical(M0);print('Kp_critique= ',Kp_crit)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.20,top=0.90)
    ax.set_xlabel(r'$M$',fontsize=25)
    ax.set_ylabel(r'$Kp$',fontsize=25)
    ax.plot(M0,Kp_crit[:],'--',color='black',linewidth=2,label=r'$K_p$ analytique')
    xmajor_ticks,xminor_ticks = np.arange(0, 0.85, 0.1),np.arange(0, 0.85, 0.02)                                              
    ymajor_ticks,yminor_ticks = np.arange(-0.7,-0.35, 0.1),np.arange(-0.7,-0.35, 0.02)                                              
    ax.set_xticks(xmajor_ticks);ax.set_xticks(xminor_ticks, minor=True) 
    ax.tick_params(labelsize=25)
    ax.set_yticks(ymajor_ticks);ax.set_yticks(yminor_ticks, minor=True)  
    ax.grid(which='both');ax.legend(loc=3,fontsize=25)                                  
    ax.axis([0.0, 0.85, -0.7,-0.35 ])

    plt.show()

    print(M0, -0.4/np.sqrt(1-M0**2))


    Kp_inc=-0.6
    set_question('1 et 2:  Kp_inc=-0.6 ')
    search_and_plot_Mcrit(Kp_inc,True)

    set_question('3:  Kp_inc=-0.4 ')
    Kp_inc=-0.4
    search_and_plot_Mcrit(Kp_inc,True)

    M0=np.linspace(0.1,0.85,16);



def Exercice8_5():
    """ 
    Paroi ondulée
    """
    set_title('Paroi ondulée en écoulement subsonique')
    n=101

    lamb=1.
    y0=0.01*lamb
    M0=0.5
    beta=np.sqrt(1.-M0**2)
    U_0=1.
    A=2*np.pi*U_0*y0/lamb
    alpha=2*np.pi/lamb
    print("A = ",A/beta)
    print("beta = ",beta)
    print("coeff de l'exponentielle =",alpha*beta)
    print("2 pi beta = ",2*np.pi*beta )

    #x=np.linspace(0,2*lamb,n)
    x=0
    y=np.linspace(0,lamb,11)

    def u_ondulee(x,y):
        """
        Solution potentielle
        """
        return A/beta*np.exp(-alpha*beta*y)*np.sin(alpha*x)
     
    def v_ondulee(x,y):
        """
        Solution potentielle
        """
        return A*np.exp(-alpha*beta*y)*np.cos(alpha*x)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.20,top=0.90)
    ax.set_xlabel(r'$x$',fontsize=25)
    ax.set_ylabel(r'$y$',fontsize=25)
    ax.plot(y/lamb, A/beta*np.exp(-alpha*beta*y))
     
    plt.show()

    print('y = ',y)
    print('')
    print('u/U0 = ',100*A/beta*np.exp(-alpha*beta*y))

def Exercice8_6():
    """ 
    paroi ondulée en écoulement transsonique
    """
    set_title('Paroi ondulée en écoulement transsonique')
    print("Attention code complexe et méthodologie pointue, paramètres très sensibles")

    import Livre.Wavy.analytique       as an
    import Livre.Wavy.numerique        as nu
  
    # Utilisation :

    #
    # PARAMETRES 
    #
    """
    7 : plot des fonctions f
    8 : calcul des Kp
    9 : calcul du Mach critique
    Les questions sont indépendantes dans ce code, et traitables dans l'ordre qu'on choisit
    7 et 8 étant faites en même temps, choisir 78 dans la liste questions[]

    L'ensemble des paramètres se trouve dans Wavy.wavy_wall.py
    """

    questions=[78,9]                 # liste des questions, choix possibles: 78 ou 9 
    prm=set_parameters_wavy_wall()
    s=WavyWallTranssonique(prm)

    """
    Les questions 7 et 9 sont indépendantes
    """
    
    for q in questions:
        if q==78:
            set_question('7: Calcul des fonctions f_n, 8 : calcul des Kp')
            s.run_wavy_wall(sol_an=True,sol_nu=True,calcul_Mc=False,courbe_Mc=False,display=True)

        if q==9:
            # question 9 :
            set_question('9: Calcul du Mach critique et de k0_crit + courbe')
            # suivant solution analytique ou numérique mettre True or False à sol_an et sol_nu  
            s.run_wavy_wall(sol_an=True,sol_nu=False,calcul_Mc=True,courbe_Mc=True,display=False)
    return
     
     
def Exercice8_7():
    """ 
    Calcul de Kpinc en fonction du Mach critique
    """
    set_title('Calcul de Kp_inc en fonction du Mach critique')
    Mcr=0.62
    beta=np.sqrt(1-Mcr**2)

    Kp_inc=beta*Kp_critical(Mcr)
    print(" Mcr=   %f, \t Kp_inc = %f"%(Mcr,Kp_inc))


    Mcr1=Mcr+0.01;
    Kp_inc1=np.sqrt(1-Mcr1**2)*Kp_critical(Mcr1)
    print(" Mcr1=   %f, \t Kp_inc1 = %f"%(Mcr1,Kp_inc1))

    gradient=(Kp_inc1-Kp_inc)/(Mcr1-Mcr)
    print("gradient d Kp/d Mcr = %f"%(gradient))
    print("gradient d Mcr /d Kp  = %f"%(1/gradient))

    erreur_relative=gradient*Mcr/Kp_inc

    print('coefficient erreur relative Mcr = %f'%(erreur_relative))
    # 
    print('k = 1/coefficient erreur relative = %f'%(1/erreur_relative))

    print("approximation : ")

    grad=2/(gamma+1) * beta/Mcr*(1+2/Mcr**2)
    print('d Kp_i/ d Mc = ',grad)

    print("Mc^2 = %f, beta^2 = %f"%(Mcr**2,beta**2))
