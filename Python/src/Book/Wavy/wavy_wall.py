# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Wavy package:
    ..
    
    * default parameters
    * main class
"""

#import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import Book.Wavy.analytique       as an
import Book.Wavy.numerique        as nu

 

#*******************************************
def set_parameters_wavy_wall():
#*******************************************
    """
    Set default parameters

    """
    prm={}  # Dictionnary of the parameters  

    prm["npt"]          = 20001           # point number in y-direction, must be large
    prm["N"]            = 6               # number of term in the serie (4 to 8)
    prm["ymax"]         = 1.0             # height of the domain in y, very sensitive parameter, low value
    prm["M0"]           = 0.8             # Upstream Mach number
    prm["a"]            = 0.01            # wall amplitude  (y_0)
    prm["lamb"]         = 1.0             # Lambda: wall wave length
    prm["gam"]          = 1.4             # Cp/Cv
    prm["eps"]          = 0.0001          # small parameter (variation) to calculate numerical gradient
    prm["display"]      = True            # to get plots
    prm["epsNL"]        = 1               #  0 <= eps <= 1,  1: non linear term is kept at 100%
    prm["erreur"]       = 1e-4            # convergence parameter in Newton method
    prm["type_CI"]      = 0               # type of initial conditions on  s
    prm["cmesh"]        = 1               # After Newton convergence, solve up to  cmesh x ylim
    prm["crit_conv"]    = 1               # convergence criteria : 
                                          # 0 : on the derivative, 1: exponential evolution, 2 : with weight
    #
    # Convergence is very difficult, please save all the parameters above and below  before changing them
    # s : vector of parameters for convergence method = functionn f_k @ y=0 (Y_2p)
    prm["test_ode"]     = False  # to test to solve the direct (primal) system and its gradient, for a given Y
    prm["test_dir"]     = False  # to test to solve the direct (primal) only
    prm["test_grad"]    = False  # integration of the primal systeme and its gradient
    prm["test_newton"]  = True   # use a Newton method to solve the problem, must always be true
    prm["test_grad_num"]= False  # to calculate the gradient by finite difference, for validation purpose
    #
    # The option to get results are below  
    # 
    prm["calcul_Mc"]    = False   # calculus of the critical Mach numbver by the Newton method
    prm["curve_Mc"]     = False   # to display the critical Mach number  curve Mc w.r.t eta
    prm["msg"]          = False   # to display on the terminal the details and comments (messages)
    #
    prm["N_an"]         = 6               # mode number for the analytical solution (must not be modified)
    prm["ylim"]         = 5               # y-direction upper domain limit for the analytical solution  
    prm["npt_an"]       = 101             # Number of points to plot the analytical solution
    # -------------------------------- 
    prm["nx"]           = 201             # number of point in  x to calculate the pressure coefficient Kp
    # -------------------------------- 
    prm["sol_an"]       = True            # True when calculatating analytical solution  
    prm["sol_nu"]       = True            # True when calculatating numerical solution
    return prm


class WavyWallTranssonique(object):
    """
    Main class Classe to calculate transonic flow on a wavy wall
    """
    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#',60*'*')
        print('# %s'%("Transonic flow over a wavy wall"))
        print('#',60*'*','\n')
        
        self.npt        = prm["npt"]      # point number in y-direction, must be large
        self.N          = prm["N"]        # number of term in the serie (4 to 8)
        self.ymax       = prm["ymax"]     # height of the domain in y, very sensitive parameter, low value
        self.M0         = prm["M0"]       # Upstream Mach number
        self.a          = prm["a"]        # wall amplitude  (y_0)
        self.lamb       = prm["lamb"]     # Lambda: wall wave length
        self.gam        = prm["gam"]      # Cp/Cv
        self.eps        = prm["eps"]      # small parameter (variation) to calculate numerical gradient
        self.display    = prm["display"]  # to get plots
        self.epsNL      = prm["epsNL"]    # 0 <= eps <= 1,  1: non linear term is kept at 100%
        self.erreur     = prm["erreur"]   # convergence parameter in Newton method
        self.type_CI    = prm["type_CI"]  # type of initial conditions on  s
        self.cmesh      = prm["cmesh"]    # After Newton convergence, solve up to  cmesh x ylim
        self.crit_conv  = prm["crit_conv"]# 0 : on the derivative, 1: exponential evolution, 2 : with weight

        #
        # Convergence is very difficult, please save all the parameters above and below  before changing them
        # s : vector of parameters for convergence method = functionn f_k @ y=0 (Y_2p)
        self.test_ode   = prm["test_ode"]  # to test to solve the direct (primal) system and its gradient, for a given Y
        self.test_dir   = prm["test_dir"]  # to test to solve the direct (primal) only
        self.test_grad  = prm["test_grad"] # integration of the primal systeme and its gradient
        self.test_newton= prm["test_newton"]    # use a Newton method to solve the problem, must always be true
        self.test_grad_num= prm["test_grad_num"]# to calculate the gradient by finite difference, for validation purpose
        #
        #  The option to get results are below  
        # 
        self.calcul_Mc  = prm["calcul_Mc"] #  calculus of the critical Mach numbver by the Newton method
        self.curve_Mc   = prm["curve_Mc"]  #  to display the critical Mach number  curve Mc w.r.t eta
        self.msg        = prm["msg"]       #  to display on the terminal the details and comments (messages)
        #
        self.N_an       = prm["N_an"]     # mode number for the analytical solution (must not be modified)
        self.ylim       = prm["ylim"]     # y-direction upper domain limit for the analytical solution  
        self.npt_an     = prm["npt_an"]   # Number of points to plot the analytical solution
        # -------------------------------- 
        self.nx         = prm["nx"]       # number of point in  x to calculate the pressure coefficient Kp
        # -------------------------------- 
        self.sol_an     = prm["sol_an"]   #  True when calculatating analytical solution  
        self.sol_nu     = prm["sol_nu"]   # True when calculatating numerical solution
        self.gam        = 1.4

    def run_wavy_wall(self, sol_an=True, sol_nu=False, calcul_Mc=False, curve_Mc=False, display=False, msg=False):
        """
        main program
        """
        # some parameters :
        self.t       = 2 *self.a /self.lamb      # relative thickness
        self.tau     = 2*np.pi*self.a/self.lamb  # a*alpha
        self.k0      = (self.gam+1)*self.tau/(-self.M0**2+1)**(3/2)
        self.y_an    = np.linspace(0,self.ylim,self.npt_an)
        self.x       = np.linspace(0,4*np.pi,self.nx)
        self.beta2   = 1-self.M0**2

        # some default parameters are modified and user dependent
        self.display, self.msg = display, msg
        self.sol_an , self.sol_nu = sol_an, sol_nu
        self.calcul_Mc, self.curve_Mc = calcul_Mc, curve_Mc

        # concatenation of the arguments to pass in a function, into a single list :
        config = [self.test_ode, self.test_dir, self.test_grad, self.test_newton, self.test_grad_num]
        real_params = [self.k0, self.ymax, self.epsNL, self.erreur, self.cmesh, self.eps]

        if sol_an :
            sol_type = "ana"
        else:
            sol_type = "num"
         
        s_init = np.zeros(self.N)
        if self.type_CI == 2 and self.sol_an: 
            print("s is initiated from the analytical solution")
            #s_init=f_an[0,:self.N]
            print('implementation problem of f_an: error')
        print("s_init = ", s_init)

        """
        Be careful :  the problem has been solved and outputs have been saved
        here it is just to plot solution from the database to get the critical Mach number
        """
        if self.sol_an:
            """
            to generate analytical solution
            """
            self._print_info_("=","Analytical solution")
            q_an, f_an = an.analytical_solution(self.N_an, self.npt_an, self.k0, self.x, self.y_an, display=False)
            field_an = self.valeurs(q_an, self.beta2)
            print('Initial conditions (analytical solution): ')
            for k in range(self.N_an):
                print("f_%2i(0) ana     = %e "%(k, f_an[0,k]))           
        
        if self.sol_nu:
            """
            to generate numerical solution
            """
            self._print_info_("=","Numerical solution")
            q_nu, f_nu, y_nu = nu.numerical_solution(self.N, self.npt, real_params, self.type_CI, config, self.crit_conv, self.x, s_init, display=False)
            field_nu = self.valeurs(q_nu,self.beta2)
            s_init = f_nu[0,:self.N]         #  s_init are updated to accelerate the critical Mach calculus, if necessary           
            print('Initial conditions (numerical solution): ')
            for k in np.arange(0, 2*self.N, 2):
                print("f_%2i(0) num     = %e "%(k, f_nu[0, k])) 

        if self.calcul_Mc:
            """
            get critical Mach number with a Newton method and a ODE system
            """
            self._print_info_("=","Critical Mach number w.r.t. eta=y0/lambda")
            print("type of solution (ana. or num.) : ", sol_type)
            self.calcul_critical_Mach(s_init, real_params, self.M0, er0=0.0001, sol_type=sol_type, display=False)
            
        if self.curve_Mc:
            self.curve_finale_critical_Mach()

        intervalle = [self.x[0], self.x[-1],-0.4, 0.2]

        if sol_an:
            if display:
                self.compare_field(self.x, field_nu, field_an, limit=intervalle, option="Kp_an")

        if sol_nu and sol_an:
            if display:
                self.compare_field(self.x, field_nu, field_an, limit=[self.x[0], self.x[-1], -2, 0], option="S-1")
                # next lines commented because written later in the figure
                #self.compare_field(self.x,field_nu,field_an,limit=intervalle                     ,option="u/U0")
                #self.compare_field(self.x,field_nu,field_an,limit=intervalle                     ,option="Kp_lin")
                #self.compare_field(self.x,field_nu,field_an,limit=intervalle                     ,option="Kp_nonlin")
                self.compare_f(self.y_an, f_an, y_nu, f_nu[:, np.arange(0, 2*self.N, 2)], [0, self.ymax, 1e-6, 1,], scale="semilogy")

        if sol_nu and (self.test_newton or self.test_dir) :
            if display:
                self.compare_field(self.x, field_nu, field_an, limit=intervalle, option="Kp_nu")    

        if display: plt.show() 
    
    """ 
        ======================================================
                figures for  f(y), Kp, S(0)-1 
        ======================================================
    """      

    def compare_f(self, y_an, f_an, y_nu, f_nu, limit, scale="linear"):
        """
        f analytical and numerical function comparisons
        """
        print('new figure')
        plt.figure(figsize=(12, 10))
        plt.title(r"Comparison of functions $f_k$")
        plt.xlabel(r'$y$', fontsize=16)
        plt.grid()    
        plt.axis(limit)
            
        N_an, N_nu = f_an.shape, f_nu.shape
        line_an, line_nu = [], []
      
        if scale=="linear":
            for k in range(N_an[1]):
                line_an[k], = plt.plot(y_an, f_an[:, k], '-', linewidth=2, marker="o", markersize=6, markevery=5, label=r"$k = $%i"%(k))
            for k in range(N_nu[1]):
                line_nu[k],= plt.plot(y_nu, f_nu[:, k], '-', linewidth=2, label=r"$k= $%i"%(k))

        elif scale=="semilogy":
            for k in range(N_an[1]):
                ligne = plt.semilogy(y_an, np.abs(f_an[:, k]), '-', linewidth=2, marker="o", markersize=6, markevery=5, label=r"$k = $%i"%(k))
                line_an.append(ligne[0])
            for k in range(N_nu[1]):
                ligne = plt.semilogy(y_nu, np.abs(f_nu[:, k]), '-', linewidth=2, label=r"$k= $%i"%(k))
                line_nu.append(ligne[0])
        first_legend = plt.legend( handles=line_an, loc=(0.2, 0.08), prop={'size': 10}, title='analytical f')
        plt.gca().add_artist(first_legend)
        second_legend = plt.legend(handles=line_nu, loc=(0.05, 0.08), prop={'size': 10}, title='numerical f')
        plt.gca().add_artist(second_legend)
        #plt.show()


    def compare_field(self, y, Q_nu, Q_an, limit=[0, 1, -1, 1], option="S-1"):
        """
        Kp comparison
        
        * Q[0] : U/U0
        * Q[1] : linear Kp
        * Q[2] : non linear Kp 
        """
        option_list=["S-1", "u/U0", "Kp_lin", "Kp_nonlin", "Kp_nu", "Kp_an"]
        k_opt = option_list.index(option)
        Titres = [r" $S(0)-1$", r" $u/U0$", r"$Kp_{lin}$", r"$Kp_{nonlin}$", r"$Kp_{num}$", r"$Kp_{ana}$"]
        plt.figure(figsize=(12, 10))
        plt.title("Comparison of %s"%Titres[k_opt])
        plt.ylabel(Titres[k_opt], fontsize=16) 
        plt.xlabel(r'$x$', fontsize=16)  
        plt.grid()
        plt.axis(limit)
        #plt.xlim(lim[0],lim[1])
        
        if k_opt < 4 :
            plt.plot(y, Q_an[k_opt], 'k-', marker="o", markersize=6, markevery=5, linewidth=2, label="analytical")
            plt.plot(y, Q_nu[k_opt], 'r--', linewidth=2, label="numerical")
        elif k_opt == 4 :
            plt.plot(y, -2*Q_nu[1], 'k-', marker="o", markersize=6, markevery=5, linewidth=2, label=r"$-2 u/U_0$")
            plt.plot(y, Q_nu[2], 'g--', marker="o", markersize=6, markevery=4, linewidth=2, label=r"$Kp_{linear}$")
            plt.plot(y, Q_nu[3], 'r-', linewidth=2, label=r"$Kp_{nonlinear}$")
        elif k_opt==5 :
            plt.plot(y, -2*Q_an[1], 'k-', marker="o", markersize=6, markevery=5, linewidth=2, label=r"$-2 u/U_0$")
            plt.plot(y, Q_an[2], 'g--', marker="o", markersize=6, markevery=4, linewidth=2, label=r"$Kp_{linear}$")
            plt.plot(y, Q_an[3], 'r-', linewidth=2, label=r"$Kp_{nonlinear}$")
            
    #    plt.plot(y,Kp1,'b-', linewidth=2,marker="o",markersize=6,markevery=5,label="analytical Kp, i=1")
        plt.legend(loc="lower left")
        #plt.show()

    """ 
        ======================================================
                get  Kp, u and|u|_wall 
        ======================================================
    """  

    def rc2(self, beta2):
        """
        (ac/a0)^2
        """
        return 1-(self.gam-1)/(self.gam+1)*beta2

    def U_U0(self, qs, beta2):
        """
        U/U0 formula 
        """
        return np.sqrt(self.rc2(beta2)/(1-beta2))*(1+beta2/(self.gam+1)*qs)

    def Kp_lineaire(self, qs, beta2):
        """
        Kp from linearized theory
        qs = -1+S
        """
        return -2*beta2/(self.gam+1)*(qs+1)

    def Kp_nonlineaire(self, qs, beta2):
        """
        Kp without  linearization
        """
        tmp = beta2+self.rc2(beta2)*(1+beta2/(self.gam+1)*qs)**2
        r_a2 = (self.gam+1)/2-(self.gam-1)/2*tmp  # (a/a0)^2
        return 2/(self.gam*(1-beta2))*(pow(r_a2, self.gam/(self.gam-1))-1)

    def valeurs(self, qs, beta2):
        """
        get the 3 functions U/U0, linear Kp and non linear Kp
        """
        return [qs, self.U_U0(qs,beta2)-1, self.Kp_lineaire(qs,beta2), self.Kp_nonlineaire(qs,beta2)]

    """ 
        =============================================
                Critical Mach number 
        =============================================
    """  

    def calcul_critical_Mach(self, s_init, params, M0_init=0.8, er0=0.00001, iterMax=25, sol_type="ana", display=False):
        """
        get critical Mach number 
        """ 
        i, erreur = 1, 1.0
        config = [False, False, False, True, False]
        eps = params[5]
        print("Newton start:  k0 = %f for Mach = %f"%(params[0], M0_init)) 
        
        def fonc(Mach):
            """
            defined function in Newton method
            """
            params[0] = (self.gam+1)*self.tau/(1-Mach**2)**(3/2)               #params=[k0,ymax,epsNL,erreur,cmesh,eps]
            if display: 
                print(" k0 = %f for Mach = %f"%(params[0],Mach)) 
            if sol_type=="num":
                q, f, y = nu.numerical_solution(self.N, self.npt ,params, self.type_CI, config, self.crit_conv, self.x, s_init, display=False)
                return max(q)
            else:
                return an.q(0, 0, params[0])

        if display: 
            print('iter = 1, target = %4.5e'%(fonc(M0_init)))
        
        Mcurrent = M0_init
        target0 = fonc(Mcurrent)
        
        while (erreur > er0) and (i <= iterMax) :
            dM = Mcurrent*eps  
            target1 = fonc(Mcurrent+dM)
            grad = (target1-target0)/dM
            dM0 = -target0/grad
            erreur = abs(dM0/Mcurrent)
            Mcurrent += dM0
            target0 = fonc(Mcurrent)
            if display:
                print('i =  %2i, \t errorr = %12.6f'%(i,erreur))
                print('M0 = ',Mcurrent," target = ",target0)
            i += 1
        if  i > iterMax :
            raise ValueError("no convergence in Newton's method")

        print("Newton end  :  k0 = %f for Mach = %f"%(params[0], Mcurrent)) 
        return Mcurrent

    def fun_calcul_Mc(self, eta, k0):
        """
        critical Mach w.r.t. the shape factor eta
        """
        tmp = 2*np.pi*(self.gam+1)/k0*eta
        return np.sqrt(1-pow(tmp, 2/3))

    def plot_Mc(self, eta, Mc, Mc_ana, Titre="Critical Mach number", par=np.array([0, 0.05])):
        """
        plot of   Mc= f(eta)
        """
        plt.figure(figsize=(12, 10))
        plt.title(Titre)
        plt.xlabel(r'$\eta=a/\lambda$', fontsize=16)
        plt.grid()
        plt.xlim(par[0], par[1])
        plt.plot(eta, Mc, '-', linewidth=1, label="numerical")
        plt.plot(eta, Mc_ana, 'o', markersize=7, markevery=4, linewidth=0, label="analytical")
        plt.legend(loc="best")
    

    def curve_finale_critical_Mach(self):
        """
        final results : Mc=f(eta)
        
        datas are for  Ymax=1, N=6
        """

        # résultats avec N=4 ?
        #eta=np.array([0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04])
        #M0=np.array([0.89461 ,0.82646,0.76463, 0.70487, 0.64509,0.58368,0.51897, 0.44866])
        
        print("result from database (accurate solutions):")
        eta = np.array([0.005, 0.01, 0.015, 0.04])
        M0 = np.array([0.893867267649 , 0.825184961907, 0.7628239872, 0.442712704326]) 
        K0 = np.array([0.836687, 0.836686, 0.836687, 0.836686])

        k0 = (self.gam+1)*2*np.pi*eta/(1-M0**2)**(3/2)
        k0_ana = 0.837684
        #plot_Mc(eta,M0,"Nombre de Mach critique",[0, 0.05])
        #plot_Mc(M0,k0,r"paramètre $k_0$",[0.4, 0.9])
        print("k0 = ", k0,)
        print("M0 = ", M0)
        k0_m = np.mean(k0)  # a test, value can be printed
        k0_m =  0.834717 
        e = np.linspace(eta[0], eta[-1], 101)
        print(k0_m, k0_ana)
        self.plot_Mc(e, self.fun_calcul_Mc(e, k0_m), self.fun_calcul_Mc(e, k0_ana), "Critical Mach number", par=[0, 0.045])
        print("k0 mean = ",k0_m)
        plt.show()

    def _print_info_(self, symb, s):
        """
        print information on screen
        symb: symbol used
        s   : message
        """
        print('#', 60*symb)
        print('# %s'%(s))
        print('#', 60*symb, '\n')
