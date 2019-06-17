# -*- coding: utf-8 -*-
"""
Created on sat May 25  2019

@author: cairiau
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Livre.Wavy.analytique       as an
import Livre.Wavy.numerique        as nu

 

#*******************************************
def set_parameters_wavy_wall():
#*******************************************
    """
    Set default parameters

    """
    prm={}  # Dictionnary of the parameters  

    prm["npt"]          = 20001           # nombre de points en y
    prm["N"]            = 6               # nombre de fonctions de la série
    prm["ymax"]         = 1.0             # dernier point en y
    prm["M0"]           = 0.8             # nombre de Mach amont
    prm["a"]            = 0.01            # amplitude de la perturbation de paroi (y_0)
    prm["lamb"]         = 1.0             # Lambda: longueur d'onde de la paroi
    prm["gam"]          = 1.4             # Cp/Cv
    prm["eps"]          = 0.0001          # coefficient pour le calcul du gradient numérique
    prm["display"]      = True            # pour afficher les graphiques
    prm["epsNL"]        = 1               # 1 : le terme non linéaire est pris à 100%
    prm["erreur"]       = 1e-4            # paramètre pour la convergence dans Newton
    prm["type_CI"]      = 0               # type de condition initiale sur s
    prm["cmesh"]        = 1               # dans Newton, calcul après convergence en cmesh x ylim
    prm["crit_conv"]    = 1               # 0: sur la dérivée, 1: évolution en exponentielle, 2 avec des poids
    # -------------------------------- 
    # il n'est pas conseillé de jouer avec les options suivantes sans conserver la version originale...
    # s : vecteur des paramètres: fonction f_k en y=0 (Y_2p)
    prm["test_ode"]     = False  # pour tester le système direct et son gradient, pour un Y donné
    prm["test_dir"]     = False  # intégration du système direct
    prm["test_grad"]    = False  # intégration du système direct et de son gradient
    prm["test_newton"]  = True   # utilisation de la méthode de Newton pour résoudre le problème (laisser à vrai)
    prm["test_grad_num"]= False  # calcul du gradient par différences finies, pour validation
    #
    # C'est ICI qu'on choisit de sortir les résultats uniquement
    # ou si on fait le calcul dans le détail.
    prm["calcul_Mc"]    = False  # calcul du nombre de Mach critique par une méthode de Newton-Raphson
    prm["courbe_Mc"]    = False   # Pour montrer la courbe Mc fonction de eta
    prm["msg"]          = False   # Pour afficher les détails des calculs (messages)
    #
    # -------------------------------- 
   
    prm["N_an"]         = 6               # nombre de modes pour la solution analytique
    prm["ylim"]         = 5               # distance maximale en y pour la solution analytique
    prm["npt_an"]       = 101             # nombre de points pour dessiner cette soltuion
    # -------------------------------- 
    prm["nx"]           = 201             # nombre de points en x pour calculer les Kp
    # -------------------------------- 
    prm["sol_an"]       = True            # vrai si on calcule la solution analytique
    prm["sol_nu"]       = True            # vrai si on calcule la solution numérique
    # ==============================================================
   
    return prm

 
class WavyWallTranssonique(object):
    """
    Classe principale pour effectuer les calculs en transsonique
    de l'écoulement autour d'une paroi ondulée
    """
    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#',60*'*')
        print('# %s'%("Transonic flow over a wavy wall"))
        print('#',60*'*','\n')
        
        self.npt        = prm["npt"]      # nombre de points en y
        self.N          = prm["N"]        # nombre de fonctions de la série
        self.ymax       = prm["ymax"]     # dernier point en y
        self.M0         = prm["M0"]       # nombre de Mach amont
        self.a          = prm["a"]        # amplitude de la perturbation de paroi (y_0)
        self.lamb       = prm["lamb"]     # Lambda: longueur d'onde de la paroi
        self.gam        = prm["gam"]      # Cp/Cv
        self.eps        = prm["eps"]      # coefficient pour le calcul du gradient numérique
        self.display    = prm["display"]  # pour afficher les graphiques
        self.epsNL      = prm["epsNL"]    # 1 : le terme non linéaire est pris à 100%
        self.erreur     = prm["erreur"]   # paramètre pour la convergence dans Newton
        self.type_CI    = prm["type_CI"]  # type de condition initiale sur s
        self.cmesh      = prm["cmesh"]    # dans Newton, calcul après convergence en cmesh x ylim
        self.crit_conv  = prm["crit_conv"]# 0: sur la dérivée, 1: évolution en exponentielle, 2 avec des poids
        # -------------------------------- 
        # il n'est pas conseillé de jouer avec les options suivantes sans conserver la version originale...
        # s : vecteur des paramètres: fonction f_k en y=0 (Y_2p)
        self.test_ode   = prm["test_ode"]  # pour tester le système direct et son gradient, pour un Y donné
        self.test_dir   = prm["test_dir"]  # intégration du système direct
        self.test_grad  = prm["test_grad"] # intégration du système direct et de son gradient
        self.test_newton= prm["test_newton"]    # utilisation de la méthode de Newton pour résoudre le problème (laisser à vrai)
        self.test_grad_num= prm["test_grad_num"]# calcul du gradient par différences finies, pour validation
        #
        # C'est ICI qu'on choisit de sortir les résultats uniquement
        # ou si on fait le calcul dans le détail.
        self.calcul_Mc  = prm["calcul_Mc"] # calcul du nombre de Mach critique par une méthode de Newton-Raphson
        self.courbe_Mc  = prm["courbe_Mc"] # Pour montrer la courbe Mc fonction de eta
        self.msg        = prm["msg"]       # Pour afficher les détails des calculs.
        #
        # -------------------------------- 
   
        self.N_an       = prm["N_an"]     # nombre de modes pour la solution analytique
        self.ylim       = prm["ylim"]     # distance maximale en y pour la solution analytique
        self.npt_an     = prm["npt_an"]   # nombre de points pour dessiner cette soltuion
        # -------------------------------- 
        self.nx         = prm["nx"]       # nombre de points en x pour calculer les Kp
        # -------------------------------- 
        self.sol_an     = prm["sol_an"]   # vrai si on calcule la solution analytique
        self.sol_nu     = prm["sol_nu"]          
        self.gam        = 1.4

    def run_wavy_wall(self,sol_an=True,sol_nu=False,calcul_Mc=False,courbe_Mc=False,display=False,msg=False):
        """
        main program
        """
        # some parameters :
        self.t       = 2 *self.a /self.lamb      # épaisseur relative
        self.tau     = 2*np.pi*self.a/self.lamb  # a*alpha
        self.k0      = (self.gam+1)*self.tau/(-self.M0**2+1)**(3/2)
        self.y_an    = np.linspace(0,self.ylim,self.npt_an)
        self.x       = np.linspace(0,4*np.pi,self.nx)
        self.beta2   = 1-self.M0**2

        # j'écrase certains paramètres initiaux qui sont des choix utilisateurs
        self.display,self.msg=display,msg
        self.sol_an,self.sol_nu=sol_an,sol_nu
        self.calcul_Mc,self.courbe_Mc=calcul_Mc,courbe_Mc

        # pour avoir moins d'arguments à passer dans les fonctions :
        config=[self.test_ode,self.test_dir,self.test_grad,self.test_newton,self.test_grad_num]
        real_params=[self.k0,self.ymax,self.epsNL,self.erreur,self.cmesh,self.eps]

        if sol_an :
            sol_type="ana"
        else:
            sol_type="num"
         
        s_init=np.zeros(self.N)
        if self.type_CI==2 and self.sol_an: 
            print("s is initiated from the analytical solution")
            s_init=f_an[0,:self.N]
        print("s_init = ", s_init)

        """
        attention : solution finale, les calculs ont déjà été faits et les valeurs récupérées
        C'est juste le graphique de la dernière question de l'exercice
        """
        if self.sol_an:
            """
            Génération de la solution analytique
            """
            self.print_info("=","Analytical solution")
            q_an,f_an=an.solution_analytique(self.N_an,self.npt_an,self.k0,self.x,self.y_an,display=False)
            champ_an=self.valeurs(q_an,self.beta2)
            print('Initial conditions (analytical solution): ')
            for k in range(self.N_an):
                print("f_%2i(0) ana     = %e "%(k,f_an[0,k]))           
        
        if self.sol_nu:
            """
            Génération de la solution analytique
            """
            self.print_info("=","Numerical solution")
            q_nu,f_nu,y_nu=nu.solution_numerique(self.N,self.npt,real_params,self.type_CI,config,self.crit_conv,self.x,s_init,display=False)
            champ_nu=self.valeurs(q_nu,self.beta2)
            s_init=f_nu[0,:self.N]          # mise à) jour de s_init pour accélérer le calcul du Mach critique, si nécessaire           
            print('Initial conditions (numerical solution): ')
            for k in np.arange(0,2*self.N,2):
                print("f_%2i(0) num     = %e "%(k,f_nu[0,k])) 

        if self.calcul_Mc:
            """
            CALCUL DU MACH CRITIQUE PAR UNE METHODE DE NEWTON-RAPSHON
            """
            self.print_info("=","Calcul du Mach critique en fonction de eta=y0/lambda")
            print("type of solution (ana. or num.) : ",sol_type)
            self.calcul_Mach_critique(s_init,real_params,self.M0,er0=0.0001,sol_type=sol_type,display=False)
            
        if self.courbe_Mc:
            self.courbe_finale_Mach_critique()


        intervalle=[self.x[0],self.x[-1],-0.4,0.2]

        if sol_an:
            if display:
                self.compare_champ(self.x,champ_nu,champ_an,limit=intervalle,option="Kp_an")

        if sol_nu and sol_an:
            if display:
                self.compare_champ(self.x,champ_nu,champ_an,limit=[self.x[0],self.x[-1],-2  ,0  ],option="S-1")
                # mis en commentaires car repris dans les figures plus bas.
                #self.compare_champ(self.x,champ_nu,champ_an,limit=intervalle                     ,option="u/U0")
                #self.compare_champ(self.x,champ_nu,champ_an,limit=intervalle                     ,option="Kp_lin")
                #self.compare_champ(self.x,champ_nu,champ_an,limit=intervalle                     ,option="Kp_nonlin")
                self.compare_f(self.y_an,f_an,y_nu,f_nu[:,np.arange(0,2*self.N,2)],[0,self.ymax,1e-6,1,],scale="semilogy")

        if sol_nu and (self.test_newton or self.test_dir) :
            if display:
                self.compare_champ(self.x,champ_nu,champ_an,limit=intervalle,option="Kp_nu")
            

        if display: plt.show() 
    
    """ 
        ======================================================
                figures pour les f(y), Kp, S(0)-1 
        ======================================================
    """      

    def compare_f(self,y_an,f_an,y_nu,f_nu,limit,scale="linear"):
        """
        Comparaison des fonctions f
        """
        print('nouvelle figure')
        plt.figure(figsize=(12,10))
        plt.title(r"Comparaison des fonctions $f_k$")
        plt.xlabel(r'$y$',fontsize=16)
        plt.grid()    
        plt.axis(limit)
            
        N_an,N_nu=f_an.shape,f_nu.shape
        line_an,line_nu=[],[]
      
        if scale=="linear":
            for k in range(N_an[1]):
                line_an[k],=plt.plot(y_an,f_an[:,k],'-', linewidth=2,marker="o",markersize=6,markevery=5,label=r"$k = $%i"%(k))
            for k in range(N_nu[1]):
                line_nu[k],=plt.plot(y_nu,f_nu[:,k],'-', linewidth=2,label=r"$k= $%i"%(k))

        elif scale=="semilogy":
            for k in range(N_an[1]):
                ligne=plt.semilogy(y_an,np.abs(f_an[:,k]),'-', linewidth=2,marker="o",markersize=6,markevery=5,label=r"$k = $%i"%(k))
                line_an.append(ligne[0])
            for k in range(N_nu[1]):
                ligne=plt.semilogy(y_nu,np.abs(f_nu[:,k]),'-', linewidth=2,label=r"$k= $%i"%(k))
                line_nu.append(ligne[0])
        first_legend=plt.legend( handles=line_an,loc=(0.2,0.08),prop={'size': 10},title='f analytique')
        plt.gca().add_artist(first_legend)
        second_legend=plt.legend(handles=line_nu,loc=(0.05,0.08),prop={'size': 10},title='f numerique')
        ax=plt.gca().add_artist(second_legend)
        #plt.show()


    def compare_champ(self,y,Q_nu,Q_an,limit=[0,1,-1,1],option="S-1"):
        """
        Comparaison des Kp
        Q[0] : U/U0
        Q[1] : Kp linéaire
        Q[2] : Kp non linéaire
        """
        option_list=["S-1","u/U0","Kp_lin","Kp_nonlin","Kp_nu","Kp_an"]
        k_opt=option_list.index(option)
            
        Titres=[r" $S(0)-1$",r" $u/U0$",r"$Kp_{lin}$",r"$Kp_{nonlin}$",r"$Kp_{num}$",r"$Kp_{ana}$"]
        

        plt.figure(figsize=(12,10))
        plt.title("Comparaison des %s"%Titres[k_opt])
        plt.ylabel(Titres[k_opt],fontsize=16) 
        plt.xlabel(r'$x$',fontsize=16)  

        plt.grid()
        plt.axis(limit)
        #plt.xlim(lim[0],lim[1])
        
        if k_opt < 4 :
            plt.plot(y,Q_an[k_opt],'k-',marker="o",markersize=6,markevery=5, linewidth=2,label="analytique")
            plt.plot(y,Q_nu[k_opt],'r--', linewidth=2,label="numérique")
        elif k_opt==4 :
            plt.plot(y,-2*Q_nu[1],'k-',marker="o",markersize=6,markevery=5, linewidth=2,label=r"$-2 u/U_0$")
            plt.plot(y,Q_nu[2],'g--',marker="o",markersize=6,markevery=4, linewidth=2,label=r"$Kp_{linear}$")
            plt.plot(y,Q_nu[3],'r-',linewidth=2,label=r"$Kp_{nonlinear}$")
            
        elif k_opt==5 :
            plt.plot(y,-2*Q_an[1],'k-',marker="o",markersize=6,markevery=5, linewidth=2,label=r"$-2 u/U_0$")
            plt.plot(y,Q_an[2],'g--',marker="o",markersize=6,markevery=4, linewidth=2,label=r"$Kp_{linear}$")
            plt.plot(y,Q_an[3],'r-',linewidth=2,label=r"$Kp_{nonlinear}$")
            
    #    plt.plot(y,Kp1,'b-', linewidth=2,marker="o",markersize=6,markevery=5,label="Kp analytique, i=1")
        plt.legend(loc="lower left")
        #plt.show()

    """ 
        ======================================================
                CALCUL des Kp et u et |u|_wall 
        ======================================================
    """  

    def rc2(self,beta2):
        """
        (ac/a0)^2
        """
        return 1-(self.gam-1)/(self.gam+1)*beta2

    def U_U0(self,qs,beta2):
        """
        formule de U/U0 
        """
        return np.sqrt(self.rc2(beta2)/(1-beta2))*(1+beta2/(self.gam+1)*qs)

    def Kp_lineaire(self,qs,beta2):
        """
        Calcul du Kp par la théorie linéarisée
        qs=-1+S
        """
        return -2*beta2/(self.gam+1)*(qs+1)

    def Kp_nonlineaire(self,qs,beta2):
        """
        Kp sans linéarisation
        """
        tmp=beta2+self.rc2(beta2)*(1+beta2/(self.gam+1)*qs)**2
        r_a2=(self.gam+1)/2-(self.gam-1)/2*tmp  # (a/a0)^2
        return 2/(self.gam*(1-beta2))*(pow(r_a2,self.gam/(self.gam-1))-1)

    def valeurs(self,qs,beta2):
        """
        Calcul des 3 fonctions U/U0, Kp linéaire et Kp nonlinéaire
        """
        return [qs,self.U_U0(qs,beta2)-1,self.Kp_lineaire(qs,beta2),self.Kp_nonlineaire(qs,beta2)]


    """ =============================================
                QUESTION 9 : CALCUL DU MACH CRITIQUE 
        =============================================
    """  

    def calcul_Mach_critique(self,s_init,params,M0_init=0.8,er0=0.00001,iterMax=25,sol_type="ana",display=False):
        """
        Calcul du Mach critique
        """ 
         
        i,erreur=1,1.0
        config=[False,False,False,True,False]
        eps=params[5]
        print("Newton début:  k0 = %f pour Mach = %f"%(params[0],M0_init)) 
        
        def fonc(Mach):
            """
            fonction définie pour Newton
            """
            params[0]=(self.gam+1)*self.tau/(1-Mach**2)**(3/2)               #params=[k0,ymax,epsNL,erreur,cmesh,eps]
            if display: print(" k0 = %f pour Mach = %f"%(params[0],Mach)) 
            if sol_type=="num":
                q,f,y=nu.solution_numerique(self.N,self.npt,params,self.type_CI,config,self.crit_conv,self.x,s_init,display=False)
                return max(q)
            else:
                return an.q(0,0,params[0])

        if display: print('iter = 1, target = %4.5e'%(fonc(M0_init)))
         
        Mcurrent=M0_init
        target0 = fonc(Mcurrent)
           
        while (erreur > er0) and (i <=iterMax) :
            dM=Mcurrent*eps  
            target1 = fonc(Mcurrent+dM)
            grad=(target1-target0)/dM
            dM0=-target0/grad
            erreur=abs(dM0/Mcurrent)
            Mcurrent+=dM0
            target0 = fonc(Mcurrent)
            if display:
                print('i =  %2i, \t erreur = %12.6f'%(i,erreur))
                print('M0 = ',Mcurrent," target = ",target0)
            i+=1
        if  i > iterMax :
            raise ValueError("pas de convergence")

        print("Newton Fin  :  k0 = %f pour Mach = %f"%(params[0],Mcurrent)) 
        return Mcurrent

    def fun_calcul_Mc(self,eta,k0):
        """
        Mach critique en fonction du facteur de forme eta
        """
        tmp=2*np.pi*(self.gam+1)/k0*eta
        return np.sqrt(1-pow(tmp,2/3))

    def plot_Mc(self,eta,Mc,Mc_ana,Titre="Mach critique",par=np.array([0,0.05])):
        """
        dessin de la solution  Mc= f(eta)
        """
        plt.figure(figsize=(12,10))
        plt.title(Titre)
        plt.xlabel(r'$\eta=a/\lambda$',fontsize=16)
        plt.grid()
        plt.xlim(par[0],par[1])
        plt.plot(eta,Mc,'-',linewidth=1,label="calcul")
        plt.plot(eta,Mc_ana,'o',markersize=7,markevery=4,linewidth=0,label="analytique")
        plt.legend(loc="best")
    

    def courbe_finale_Mach_critique(self):
        """
        Résultat final du problème
        Mc=f(eta)
        Las valeurs correspondent au cas Ymax=1, N=6
        """

        # résultats avec N=4 ?
        #eta=np.array([0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04])
        #M0=np.array([0.89461 ,0.82646,0.76463, 0.70487, 0.64509,0.58368,0.51897, 0.44866])
           
        print("résultats d'après les calculs effectués précisément:")
        eta=np.array([0.005, 0.01,0.015,0.04])
        M0=np.array([0.893867267649 , 0.825184961907,0.7628239872, 0.442712704326]) 
        K0=np.array([0.836687,0.836686, 0.836687, 0.836686])

        k0= (self.gam+1)*2*np.pi*eta/(1-M0**2)**(3/2)
        k0_ana=0.837684
        #plot_Mc(eta,M0,"Nombre de Mach critique",[0, 0.05])
        #plot_Mc(M0,k0,r"paramètre $k_0$",[0.4, 0.9])
        print("k0 = ",k0,)
        print("M0 = ",M0)
        k0_m=np.mean(k0)
        k0_m= 0.834717 
        e=np.linspace(eta[0],eta[-1],101)
        print(k0_m,k0_ana)
        self.plot_Mc(e,self.fun_calcul_Mc(e,k0_m),self.fun_calcul_Mc(e,k0_ana),"Nombre de Mach critique",par=[0,0.045])
        print("moyenne k0 = ",k0_m)
        plt.show()

    def print_info(self,symb,s):
        """
        print information on screen
        symb: symbol used
        s   : message
        """
        print('#',60*symb)
        print('# %s'%(s))
        print('#',60*symb,'\n')
