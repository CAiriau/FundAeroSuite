#!/bin/py
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Turbulent boundary layer package
    ..
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import stats
from scipy.integrate import odeint


#*******************************************
def set_parameters():
#*******************************************
    """
    Class default parameter
    """
    prm={}  
    prm["name"]     = "initial velocity profile"
    prm["ypMin"]    = 0.001
    prm["ypMax"]    = 250
    prm["n"]        = 101
    prm["kappa"]    = 0.41
    prm["Utau"]     = 0
    prm["Rtau"]     = 1000
    prm["Cf"]       = 0.003
    prm["Model"]   = "Linear"
    #        "Laminar" : laminar boundary layer with suction
    #        "SCVonly"   : viscous sublayer approximation
    #        "CIonly"    : internal layer approximation  
    #        "SCVlinear" : u+=y+
    #        "Linear"    : L+ = kappa y+   pour tout y+
    #        "Mixed"     : "Linear" if  y+ <= ypScv else "Michel"
    #        "Michel"    : Michel's mixing lenght model
    #        "PowerLaw" : power law boundary layer velocity profile
    prm["Damping"]  = "No"
    prm["v0p"]      = 0
    prm["Rv0"]      = 0
    prm["ypScv"]    = 40
    prm["Ue"]       = 1
    prm["Uep"]      = 0
    prm["nu"]       = 1e-6
    prm["rho"]      = 1000
    prm["Aplus"]    = 26.0
    prm["AplusOpt"] = 0
    #                      A^+_0       |   A^+(v0+,p+)
    #                     ------------------------------
    #               = 0     fixé=A+    |  not calculated
    #               = 1     recalculé  |  not calculated
    #               = 2     fixé=A+    |  calculated
    #               = 3     calculated |  calculated
    prm["BplusOpt"] = 0  # same table as for  AplusOpt
    
    prm["AplusModel"]= ""
    prm["Bplus"]    = 36.0
    prm["Pplus"]    = 0.0
    prm["A"]        = 1/prm["kappa"]
    prm["B"]        = 0
    prm["Cu"]       = 0
    prm["Rex"]      = 0
    prm["N"]        = np.float(7.0)   # profile at power  1/N
    prm["ShowParameters"] = False 
    prm["com"]      = " "
    return prm


# to automatically calculate  A+ or B+, set initial values to 0.

class ProfilTurbulent(object):
    """
    properties and characteristics of boundary layer turbulent profile over a flat plate
    """

    def __init__(self, option, **kwargs):
        self.option = option
        self.name   = option['name']
        self.ypMin  = option['ypMin']               # y+, lower limit of the interval of y+
        self.ypMax  = option['ypMax']               # y+ upper limit of the interval of y+
        self.n      = option['n']                   # number of points (size of y+ vector)
        self.kappa  = option['kappa']               # Karman's constant    
        self.Utau   = option['Utau']                # turbulent skin friction velocity U_tau
        self.Rtau   = option['Rtau']                # Reynolds based on skin friction velocity
        self.Cf     = option['Cf']                  # skin friction coefficient
        self.Modele = option['Model']               # type of model
        self.Damping = option["Damping"]            # Application of tje damping correction in ["No",  "VanDriest" ,  "Kays"]
        self.Aplus  = option["Aplus"]               # constant in Van Driest law
        self.AplusOpt = option["AplusOpt"]          # option to calculate A+ wrt v0+ and p+  (True or False)
        self.Bplus  = option["Bplus"]               # constant in Kays damping law
        self.BplusOpt = option["BplusOpt"]          # option to calculate B+ wrt v0+ and p+  (True or False)
        self.AplusModel=option["AplusModel"]        # Model to calculate  A+ in ["Kays1", "Kays2", "Kays3", "Cebeci1"]
        self.pplus  = option["Pplus"]               # pressure gradient p+
        self.v0p    = option['v0p']                 # suction velocity / skin friction velocity (non dimensional value)
        self.v0     = ''                            # suction velocity (>0)
        self.Rv0    = option['Rv0']                 # Reynolds = v_0 delta/nu
        self.ypScv  = option['ypScv']               # y+ at the theoretical viscous sublayer boundary 
        self.Ue     = option['Ue']                  # external boundary layer velocity
        self.Uep    = option['Uep']                 # non dimentional external boundary layer velocit with U_tau : Ue+
        self.nu     = option['nu']                  # laminar kinematic viscosity
        self.rho    = option['rho']                 # density
        self.mu     = self.nu*self.rho              # dynamic viscosity        
        self.A      = option["A"]                   # log law slope
        self.B      = option["B"]                   # constant in log 
        self.N      = option["N"]                   # for power law, power coefficient of velocity in  1/N
        self.Cu     = option["Cu"]                  # multiplicative constant coefficient in the power law 
        self.yv     = ''                            # y_v  giving the boundary of the viscous sublayer with suction
        self.yp     = ''                            # y+
        self.up     = ''                            # u+
        
        self.y_vp   = ''                            # y_v^+, tests
        self.u_vp   = ''                            # u_v^+ tests
        self.y_sv   = 0
        self.up_sv  = 0
        self.plot_sv= False                         # for plotting  y_v^+ u_v^+
        self.dup_dyp= ''                            # du+/dy+ (output)
        self.Rdelta = ''                            # Reynolds based on the boundary layer thickness
        self.deltaCL= ''                            # boundary layer thickness
        self.deltav = ''                            # thickness nu/v0 (suction)
        self.Lplus  = ''                            # mixing length
        self.k1     = ''                            # delta1/delta
        self.k2     = ''                            # delta2/delta
        self.H      = ''                            # shape factor
        self.yp     = np.float(0.)
        self.Rex    = option["Rex"]                 # Re_x
        self.erreur = False
        self.ShowParameters = option['ShowParameters']
        self.ifig   = 0                             # figure index
        self.com    = option["com"]                 # comments

        self.set_title()

        if self.Modele == "Laminar" and self.Rv0 == 0 :
            self.erreur = True
            print("Laminar B.L. without suction not implemented  !")
        
        if self.v0p != 0:
            print("profile with suction")

        if self.pplus != 0:
            print("profile with pressure gradient")

        if self.Damping != "No":
            print("Need to calculate A+ or B+")
            if self.Damping == "Kays":
                print('B0+ (input)         = ',self.Bplus)
                if self.BplusOpt == 0 or self.BplusOpt == 2 :
                    print('B+_0 = B+ input')
                else :
                    self.Bplus = self.Bplus_function_kappa(self.kappa)
                    print('B0+ (calculated for kappa = %4.3f) = %4.3f'%(self.kappa, self.Bplus))
                if self.BplusOpt > 1 :
                    self.Bplus = self.B_Kays3(self.v0p,self.pplus)
                    print('B+ (v0+ = %4.3f, p+ = %4.3f)  = %4.3f'%(self.v0p, self.pplus, self.Bplus)) 

            if self.Damping == "VanDriest":
                print('A0+ (input)         = ',self.Aplus)
                if self.AplusOpt == 0 or self.AplusOpt == 2 :
                    print('A+_0 = A+ input')
                else :
                    self.Aplus = self.Aplus_function_kappa(self.kappa)
                    print('A0+ (calculated for kappa = %4.3f) = %4.3f'%(self.kappa, self.Aplus)) 
                if self.AplusOpt > 1 :
                    if self.AplusModel == 'Kays1':
                        s = self.A_vanDriest_Kays1(self.v0p, self.pplus)
                    elif self.AplusModel == 'Kays2':
                        s = self.A_vanDriest_Kays2(self.v0p, self.pplus)
                    elif self.AplusModel == 'Kays3':
                        s = self.A_vanDriest_Kays3(self.v0p, self.pplus)
                    elif self.AplusModel == 'Cebeci1':
                        s = self.A_vanDriest_Cebeci1(self.v0p, self.pplus)
                    self.Aplus = s
                print('A+ (v0+ = %4.3f, p+ = %4.3f)  = %4.3f'%(self.v0p, self.pplus, self.Aplus)) 
        if (self.Modele == "SCVonly" or self.Modele == "CIonly") and (self.v0p != 0 or self.pplus != 0): 
            self.com = "Faux"

    def get_bl_profile(self):
        """
        grid generation and get B.L. profile
        """
        if not self.erreur: 
            self.grid_generation()
            if self.Modele=="PowerLaw":
                self.power_profile()
            else:
                print('other profile than power law')
                self.evaluate_parameters()
                self.genere_profile()
            #print('U = ',self.up)
            print("Rdelta                         = %f "%(self.Rdelta))
            print("Max u+                         = %f "%(max(self.up)))

    # **********************************************
    #  POWER LAW for VELOCITY PROFILES
    # **********************************************

    def calcul_Rex(self):
        """
        Reynolds number based on x coordinate
        """
        self.Rex = pow(self.Rtau,1+3/self.N)*self.Cu**3*self.N/((self.N+2)*(self.N+3))

    def calcul_Cf(self):
        """
        skin friction coefficient with power law
        """
        self.Cf = 2/(self.Cu**2*pow(self.Rtau, 2/self.N))

    def set_Cu(self):
        """
        proportional coefficient on the power law  
        """
        if self.Cf != 0 :
            self.Cu = pow(self.Rtau,1/self.N)*np.sqrt(2/self.Cf)
        elif self.Rex != 0 :
            tmp = self.Rex*(self.N+2.)*(self.N+3.) / (self.N*pow(self.Rtau,1.+3./self.N))
            #tmp=(self.N+2.)*(self.N+3.) / (self.N*pow(self.Rtau,1.+3./self.N))
            #print("Rex = ",self.Rex)
            self.Cu = pow(tmp,1./3.)
        else:
            raise NameError("Problem with Cu calculus")

    def set_Cu_from_Cf(self):
        """
        from Rex, Cf, N get  Cu, Rtau
        """
        tmp = (self.N+2)*(self.N+3)*self.Rex/self.N*pow(self.Cf/2,(self.N+3)/2)
        self.Cu= pow(tmp,-1/self.N)

    def calcul_Rtau(self):
        """
        Rtau from Cf with the parameters N, Rex, Cu
        """
        self.Rtau = pow(self.Cf/2,-1*self.N/2)/pow(self.Cu,self.N)
        self.Rtau = pow(np.sqrt(2/self.Cf)/self.Cu,self.N)


    def power_profile(self):
        """
        get power law profile
        """

        def f_N(x):
            """
            Equation to solve with N as unknown 
            """
            return self.Rex-pow(self.Rtau,1+3/x)*self.Cu**3*x/((x+2)*(x+3))

        print("Power law characteristics")
        
        if self.N == 0 :
            print("solve N")
            self.N = fsolve(f_N,7)
            print('new value, N : ',self.N)
            self.calcul_Cf()
        elif self.Cf*self.Rtau != 0:   # Rtau and Cf known => calculate  Cu and Rex
            self.set_Cu()
            self.calcul_Rex()
        elif self.Rtau*self.Cu != 0:   # Rtau and Cu knwon > calculate  Cf and Rex
            self.calcul_Cf()
            self.calcul_Rex() 
        elif self.Rtau*self.Rex != 0: # Rtau and Rex known > calculate Cf and Cu
            print("Rtau and Rex given")
            self.set_Cu()
            self.calcul_Cf()
        else:
            print('Problem with this configuration!  verify your inputs : Cu, Rtau, Cf and Rex')
            print('Cu and Rtau arbitrarily given')
            self.Cu, self.Rtau = 7, 1000
            self.calcul_Cf()
            self.calcul_Rex() 

        self.k1, self.k2 = 1/(self.N+1), self.N/((self.N+1)*(self.N+2))
        self.Uep = self.Cu*pow(self.Rtau, 1/self.N)
        self.Rdelta = self.Rtau*np.sqrt(2/self.Cf)
        self.Rex = pow(self.Rtau, 1+3/self.N)*self.Cu**3*self.N/((self.N+2)*(self.N+3))
        self.H = self.k1/self.k2
        print("Rtau                           = %f "%(self.Rtau))
        print("Uep                            = %f "%(self.Uep))
        print("Cf                             = %f "%(self.Cf))
        print("Cu                             = %f "%(self.Cu))
        print("Rdelta                         = %f "%(self.Rdelta))
        print("Re_x                           = %e "%(self.Rex))
        print("delta_1/delta                  = %f "%(self.k1))
        print("delta_2/delta                  = %f "%(self.k2))
        print(" H                             = %f "%(self.H))

        self.up = self.Cu*pow(self.yp, 1/self.N)
        self.dup_dyp = self.up/self.N/self.yp
        # boundary layer characteristics
        C = pow(self.Cu,-2*self.N/(self.N+3))
        A = (self.N+2)*(self.N+3)/self.N
        N1 = (self.N+1)/(self.N+3)
        N2 = -2/(self.N+3)
        print("delta/x     = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A, N1), -2/(self.N+3), -2, self.N+3))
        print("Rdelta      = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A, N1), N1, self.N+1, self.N+3))
        print("Rdelta_1    = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A, N1)*self.k1, N1, self.N+1, self.N+3))
        print("Rdelta_2    = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A, N1)*self.k2, N1, self.N+1,  self.N+3))
        print("Cf/2        = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A, N2), N2, -2, self.N+3))
        print("Rtau        = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(pow(self.Cu, -3*self.N/(self.N+3))*pow(A, self.N/(self.N+3)), 
                                                                    self.N/(self.N+3), self.N, self.N+3))
        self.power_law()
        
        print("Michel's law                Cf = %f"%(self.Michel(self.Rex)))
        print("Schultz and Grunow law      Cf = %f"%(self.SchultzGrunow(self.Rex)))
        Rdt2 = self.k2*self.Rdelta
        print("Ludwieg et Tillman law      Cf = %f"%(self.LudwiegTillman(Rdt2,self.H)))
        print("Felsch et al                Cf = %f"%(self.Felsch(Rdt2,self.H)))

    def power_law(self,display=True):
        """
        power law for turbulent boundary layer, characteristics
        """
        C = pow(self.Cu, -2*self.N/(self.N+3))
        A = (self.N+2)*(self.N+3)/self.N
        N1 = (self.N+1)/(self.N+3)
        N2 = -2/(self.N+3)

        print("delta/x     = %f "%(C*pow(A,N1)* pow(self.Rex, -2/(self.N+3)) ))
        print("Rdelta      = %f "%(C*pow(A,N1)* pow(self.Rex, N1)            ))
        print("Rdelta_1    = %f "%(C*pow(A,N1)*self.k1*pow(self.Rex, N1)      ))
        print("Rdelta_2    = %f "%(C*pow(A,N1)*self.k2*pow(self.Rex, N1)      ))
        print("Cf          = %f "%(2*C*pow(A,N2)*pow(self.Rex, N2)           ))
        print("Rtau        = %f "%(pow(self.Cu, -3*self.N/(self.N+3))*pow(A, self.N/(self.N+3) )*pow(self.Rex, self.N/(self.N+3))))

    # *********************************************************
    #  GENERATION OF TURBULENT PROFILES FOR DIFFERENT MODELS
    # *********************************************************    

    def grid_generation(self):
        """
        Grid generation in y+
        """    
        if self.Modele == "Laminar" :
            self.yp = np.linspace(0, self.ypMax, self.n)
            #print("eta = ",self.yp)
        else:
            self.yp = np.exp(np.linspace(np.log(self.ypMin), np.log(self.ypMax), self.n))
            #print("y+ = ",self.yp)

    def  genere_profile(self):
        """
        Boundary layer profile generation
        """        
        if self.Modele == "Laminar":
            self.up = 1-np.exp(-1*self.yp)
            self.dup_dyp = np.exp(-1*self.yp)

        elif self.Modele == "SCVonly":
            print("BE CAREFUL : wrong model:  A+ assumed constant")  
            self.yv = self.calcul_yv()        # wrong since A+(v0+, p+)
            self.ypMax = self.yv
            self.grid_generation()           # need to generate a new grid, the first point has changed
            self.up = (1-np.exp(-self.v0p*self.yp))/self.v0p
            self.dup_dyp = np.exp(-self.v0p*self.yp)
        elif self.Modele == "CIonly":
            print("BE CAREFUL : wrong model:  A+ assumed constant")  
            self.yv=self.calcul_yv()        
            print(" yv = ",self.yv)
            self.ypMin=self.yv
            self.grid_generation() # need to generate a new grid, the first point has changed
            m = -self.v0p/(2*self.kappa)
            uv = (1-np.exp(-self.v0p*self.yv))/self.v0p
            tmp = m*np.log(self.yp/self.yv)+np.sqrt(1-self.v0p*uv)
            self.up = (1-tmp**2)/self.v0p
            self.dup_dyp = self.up  #It is wrong. just to fill the vector
            self.calcul_A_B() 
        elif self.Modele=="SCVlinear":
            self.up = self.yp    
            self.dup_dyp = np.ones(len(self.yp))
        elif self.Modele == "loi_log":
            self.loi_log()
        else :
            self.turbulent_profile_with_suction()

        print("Max u+                         = %f "%(max(self.up)))
        
        self.y_vp = 0
        if self.Rtau == self.ypMax:
            print("Profile defined in the whole B.L. thickness")
            print("Characteristics are calculated")
            self.Uep = max(self.up)
            self.Rdelta = self.Uep*self.Rtau
            self.Cf = 2/self.Uep**2
            print('NEW OUTPUTS')
            print('')    
            print('Ue+                 = %f'%(self.Uep))
            print('Rdelta              = %f'%(self.Rdelta))
            print('Cf                  = %f'%(self.Cf))

            if self.v0p != 0:
                self.plot_sv = True
                self.y_vp = self.yp
                self.u_vp = self.solve_yv_uv(self.yp)      #-self.up[:,0]
                #print(np.where(self.yp>10))
                v = self.u_vp-self.up[:, 0]
                xloc = np.where(self.yp > 10)
                vloc, yloc = v[xloc], self.yp[xloc]
                #print("yloc = ",yloc)
                #print("vloc = ",vloc)
                i = 0
                while vloc[i] > 0 and i < len(vloc):
                    i = i+1
                #print("i = ",i,vloc[i],yloc[i-1],yloc[i])
                self.yp_sv = yloc[i-1]-(yloc[i]-yloc[i-1])/(vloc[i]-vloc[i-1])*vloc[i-1]
                self.up_sv = self.solve_yv_uv(self.yp_sv)
                print("Viscous sublayer boundary   : y+ = %f, u+= %f"%(self.yp_sv, self.up_sv))

        else:
            print("The velocity profile is not defined in the whole B.L. height")


    # **********************************************
    #  VISCOUS SUBLAYER BOUNDARY
    # **********************************************

    def solve_yv_uv(self,y):
        """
        Viscous sublayer boundary
        plot the function to see ...
        """
        tmp = np.sqrt(1-self.v0p*self.Uep)+self.v0p/(2*self.kappa)*np.log(self.Rtau/y)
        return (1-tmp**2)/self.v0p

    # **********************************************
    #  DAMPING LAWS
    # **********************************************

    def regression_lineaire(self, filename, skiprows=0):
        """
        linear regression of a curve found in a file
        """
        reference_filepath = os.path.join('Book/Data',filename)
        with open(reference_filepath, 'r') as infile:
            Data = np.loadtxt(infile, dtype=float, unpack=True, skiprows=skiprows)
        slope, intercept, r_value, p_value, std_err = stats.linregress(Data[0], Data[1])
        print('y = %f + %f x , correlation = %f, sigma = %f '%(intercept, slope, r_value, std_err))

        # second way to do it, with less outputs
        coef = np.polyfit(Data[0], Data[1], 1)
        print("Polynomial : ",coef)
        # p=np.poly1d(coef)   
        return [intercept,slope]


    def A_B_regression(self):
        """
        linear regression of  A+ and B+, Kays's plots
        Int J. Heat Mass transfer. Vol. 15. pp. 1023-1044. W.M. Kays, 1972
        """
        A = self.regression_lineaire('A_function_K.dat', skiprows=2)
        B = self.regression_lineaire('B_function_K.dat', skiprows=2)
        return A,B

    def plot_loi_Van_Driest(self,xaxis="Log"):
        """
        plot  Van Driest wall
        """
        self.grid_generation()
        self.up = 1-np.exp(-1*self.yp/self.Aplus)
        self.plot_profile(mode=xaxis, legx=r"$y^+$", legy=r"$\chi(y^+)$")

    def Aplus_function_kappa(self,kappa):
        """
        Constant A+=f(kappa) (kappa = Karman's constant)
        after linear regression from Kays's data
        """
        return -13.047+90.171*kappa

    def Bplus_function_kappa(self,kappa):
        """
        Constant B+=g(kappa) (kappa = Karman's constant)
        after linear regression from Kays's data
        """
        return -8.0503+102.24*kappa

    def B_Kays3(self,v0plus,Pplus):
        """
        constant for damping law w.r.t  presssure gradient and transpiration velocity (Kays)
        Int J. Heat Mass transfer. Vol. 15. pp. 1023-1044. W.M. Kays, 1972
        p. 1035
        """
        if isinstance(Pplus, float):
            return self.Bplus/( 9*v0plus+3.35*Pplus/(1+4.0*v0plus)+1.0 )
        else:
            s = [] 
            for p in Pplus:
                sol = self.Bplus/( 9*v0plus+3.35*p/(1+4.0*v0plus)+1.0 )
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
            return s        
    

    def A_vanDriest_Kays3(self, v0plus, Pplus):
        """
        van Driest's constant  w.r.t  presssure gradient and transpiration velocity (Kays)
        Int J. Heat Mass transfer. Vol. 15. pp. 1023-1044. W.M. Kays, 1972
        p. 1035
        """
        if isinstance(Pplus, float):
            return self.Aplus/( 5.15*(v0plus+5.86*Pplus/(1+5.0*v0plus))+1.0 )
        else:
            s = [] 
            for p in Pplus:
                sol = self.Aplus/( 5.15*(v0plus+5.86*p/(1+5.0*v0plus))+1.0 )
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
            return s        


    def A_vanDriest_Kays2(self, v0plus, Pplus):
        """
        van Driest's constant  w.r.t  presssure gradient and transpiration velocity (Kays), 
        approach with order 2 surface,
        limited by  v0+ and p+ values,
        Andersen, Kays et Moffat, JFM 69,part 2, 1975, p.364
        """
        x = np.log(v0plus+0.08)
        a = [[-6.71751, -1.50414], [-4.81589,-1.24311], [-1.27827, -0.388216]]

        if isinstance(Pplus, float):
            y = np.log(Pplus+0.04)
            u = [0, 0, 0]
            for i in range(3):
                u[i] = a[i][0] + a[i][1]*y
                #print('a[ %i , 0] = %f \t a[ %i , 1] = %f '%(i,a[i][0],i,a[i][1]))
            tmp = u[0]+x*u[1]+x**2*u[2]
            return 24*np.exp(tmp)

        else:
            s = []
            for p in Pplus:
                y = np.log(p+0.04)
                u = [0, 0, 0]
                for i in range(3):
                    u[i] = a[i][0]+a[i][1]*y
                    #print('a[ %i , 0] = %f \t a[ %i , 1] = %f '%(i,a[i][0],i,a[i][1]))
                tmp = u[0]+x*u[1]+x**2*u[2]
                sol = 24*np.exp(tmp)
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
            return s

    def A_vanDriest_Kays1(self,v0plus,Pplus):
        """
        van Driest's constant  w.r.t  presssure gradient and transpiration velocity (Kays), 
        Andersen, Kays et Moffat, JFM 69,part 2, 1975, P.367
        """
        if v0plus < 0 : 
            a = 9.0
        else :
            a = 7.1
        
        if isinstance(Pplus, float):
            b, c = 4.25, 10.0
            if Pplus > 0  :
                b, c = 2.9, 0.0     # b= 2.9 in the paper and  = 2 in the report
            return self.Aplus/(a*(v0plus+b*Pplus/(1+c*v0plus))+1.0)

        else:
            s = []
            for p in Pplus:
                b, c = 4.25, 10.0
                if p > 0  : 
                    b, c = 2.9, 0.0    
                # can the  denominator be null ? 
                # alpha,beta,gam=a*c,a+c,1+a*b*p
                # delta=beta**2-4*alpha*gam
                # if delta >= 0:
                #     print("p+ = %f D= 0 pour v0 = %f et   %f MAIS v0 =  %f"%(p,(-beta+np.sqrt(delta))/(2*alpha),(-beta-np.sqrt(delta))/(2*alpha),v0plus))
                sol = self.Aplus/(a*(v0plus+b*p/(1+c*v0plus))+1.0)
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
            return s

    def A_vanDriest_Cebeci1(self, v0plus, Pplus):
        """
        van Driest's constant  w.r.t  presssure gradient and transpiration velocity.
        Cebeci, 1969
        """
        yplus = 11.8
        eps = 1.0 
        # in the reference it is  -1, but trends are opposite to Kays, we have to set the value to +1,
        # Kays's definition is the opposite of  Cebeci's one ( dp/dx for Kays, -dp/dx for Cebeci)
        if isinstance(Pplus, float):
            if v0plus == 0 : 
                return self.Aplus/np.sqrt(1+eps*Pplus*yplus)
            else :
                tmp = np.exp(v0plus*yplus)
                return self.Aplus/np.sqrt(tmp+eps*Pplus/v0plus*(tmp-1))
        else:
            s = []
            for p in Pplus:
                if v0plus == 0 : 
                    sol = self.Aplus/np.sqrt(1+eps*p*yplus)
                else :
                    tmp = np.exp(v0plus*yplus)
                    sol = self.Aplus/np.sqrt(tmp+eps*p/v0plus*(tmp-1))
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
            return s

    # **********************************************
    #  viscous sublayer and layer   
    # **********************************************

    def calcul_A_B(self):
        """
        set log law coefficients
        """
        n_inf, n_sup = np.argmin(np.abs(self.yp-100)), np.argmin(np.abs(self.yp-1000))
        print("Log law constants: A and B, n_inf = %i, n_sup = %i"%(n_inf, n_sup))
        self.A = (self.up[n_sup]-self.up[n_inf])/np.log(self.yp[n_sup]/self.yp[n_inf])
        self.B = self.up[n_inf]-self.A*np.log(self.yp[n_inf])
        print('A  + %f, B = %f, 1/A = %f'%(self.A,self.B,1/self.A))

    def loi_log(self):
        """
        asymptotic log law with suction
        """
        ul = self.B+self.A*np.log(self.yp)
        print('A  + %f, B = %f, 1/A = %f'%(self.A,self.B,1/self.A))
        self.up = ul[ul > 0]
        self.yp = self.yp[ul > 0]                
        self.dup_dyp = self.A/self.yp

    def SCV(self,yp,v0):
        """
        viscous sublayer with approximation
        """
        return (1-np.exp(-v0*yp))/v0

    def calcul_yv(self):
        """
        viscous sublayer thickness w.r.t. various criteria 
        v0 is v0+ , return y_v^+

        .. warning::
            it is assumed that A+ does not depend on v0+ and p+, but it is wrong
        """
        def f1(y):
            """
            Continuous derivative , exact value obtained from the ODE
            """
            return 1-(self.kappa*y)**2*np.exp(-self.v0p*y)

        def f3(y):
            """
            Continuous derivative calculate with the solution
            """
            return -2*np.exp(-self.v0p*y)*self.kappa*y+2*np.sqrt(np.exp(-self.v0p*y))
        def s2(v0):
            """
            Continuous derivative, first order expansion of f1  
            """
            # instead of looking for the zeros of this function I'm looking for the polynomial roots of
            # kappa**2*v0*y**3-kappa**2*y**2+1
            s = np.poly1d([self.kappa**2*self.v0p,-self.kappa**2,0,1])
            I = np.argmin((s.r-2)**2)
            return s.r[I]

        def s4(v0):
            """
            Continuous derivative, first order expansion of f2 and roots of a 2nd order polynomial
            """
            eta = v0/self.kappa
            if eta <= 6-4*np.sqrt(2):
                return  (2+eta-np.sqrt(eta**2-12*eta+4))/(4*self.v0p)
            else:
                print('limit value : %f'%( (6-4*np.sqrt(2))*self.kappa))
                print('error on the approximated value of y_v')
                return 1e10
            
        print("Problem in y_v calculus")
        y0 = 0.1
        s1 = fsolve(f1, y0)
        s3 = fsolve(f3, y0)
        if self.ShowParameters:
            print('viscous sublayer boundary')
            print('prb 1 : exact ', s1[0])
            print('prb 3 : exact ', s3[0])
            print('prb 2 : appro.', s2(self.v0p))
            print('prb 4 : appro.', s4(self.v0p))
        return s1[0]        

    
            # modele =
    #        "SCVonly" : approximated viscous sublayer
    #        "CIonly"  : approximated internal layer
    #        "Linear"  : L+ = kappa y+   for all y+
    #        "Mixed"   : "Linear" if  y+ <= ypScv else  "Michel"
    #        "Michel"  : Michel's mixing length
    #        "PowerLaw" : power law


    def calculate_mixing_length(self):
        """
        calculate mixing lenght, for comparisons
        """
        Lplus = []
        for yp in self.yp:
            Lplus.append(self.mixing_length(yp))
        self.Lplus = Lplus

    def mixing_length(self,yp):
        """
        choose the mixing length model
        
        Params:
            yp (real) :y+
        
        Returns:
            real : L
        """
        c = 0.085
        k = self.Rtau*c
        L = 0.0
        
        if self.Modele == "Linear":
            # print('Linear',yp)
            L = self.kappa*yp
        elif self.Modele == "Michel":
            L = k*np.tanh(self.kappa*yp/k)
        elif  yp <= self.ypScv:
            L=self.kappa*yp
        else:
            # print('Mixed autre ',self.Modele)
            L = k*np.tanh(self.kappa*yp/k)
        if self.Damping == "VanDriest":
                L = L*(1-np.exp(-yp/self.Aplus))
        elif self.Damping == "Kays":   
                if yp <= self.Bplus:
                    L = L*yp/self.Bplus
        return L


    def dudy(self, u, y):
        """
        du+/dy+
        """
        L2 = self.mixing_length(y)**2
        dudy = (-1+np.sqrt(1+4*L2*(1-u*self.v0p)))/(2*L2)
        dudy = 2*(1-u*self.v0p)/(1+np.sqrt(1+4*L2*(1-u*self.v0p)))
        return dudy


    def turbulent_profile_with_suction(self):
        """
        profile with suction
        """
        self.up = odeint(self.dudy,self.ypMin*(1-self.v0*self.ypMin),self.yp)
        #n_inf,n_sup=np.argmin(np.abs(yp-100)), np.argmin(np.abs(yp-1000))
        #A,B=calcul_AB(yp,up,[n_inf,n_sup])
        #print('A  + %f, B = %f, 1/A = %f'%(A,B,1/A))
        #u_log,y_log=loi_log(yp,A,B)
        self.dup_dyp = self.dudy(self.up[:,0],self.yp) # take care of the dimentions for odeint outputs!


    def evaluate_parameters(self):
        """
        Test to calculate  different parameters of the T B L
        """
        eps = 1e-7
        if self.Cf != 0 :
            self.Uep = np.sqrt(2/self.Cf)
        else:
            self.Cf = 2/self.Uep
        if self.Rv0 == 0:
            self.Rv0 = self.v0p*self.Rtau
        else:
            self.v0p = self.Rv0/self.Rtau

        self.Rdelta  = self.Uep*self.Rtau
        self.deltaCL = self.nu*self.Rdelta/self.Ue
        self.Utau    = self.Ue/self.Uep
        self.v0      = self.Rv0*self.nu/self.deltaCL
        if np.abs(self.v0) > eps:
            self.deltav  = self.nu/self.v0
        else:
            self.deltav = 1e10

        if self.ShowParameters:
            print('INPUTS')
            print('Cf                  = ', self.Cf)
            print('R_tau               = ', self.Rtau)
            print('U_e                 = ', self.Ue)
            print('nu                  = ', self.nu)
            print('Re_v_0              = ', self.Rv0)
            print('A+                  = ', self.Aplus)
            print('B+                  = %4.3f'%(self.Bplus_function_kappa(self.kappa)))
            print('A+ (calculus)       = %4.3f'%(self.Aplus_function_kappa(self.kappa)))
            print('')
            print('OUTPUTS')
            print('')    
            print('Ue+                 = ', self.Uep)
            print('Rdelta              = ', self.Rdelta)
            print('delta_CL            = ', self.deltaCL)
            print('U_tau               = ', self.Utau)
            print('Rtau (verification) = ', self.Utau*self.deltaCL/self.nu)
            print('v0+                 = ', self.v0p)
            print('v0/Ue               = ', self.Rv0/self.Rdelta)
            print('v0                  = ', self.v0)
            if self.deltav != 1e10:
                print('delta=nu/v0         = ', self.deltav)
                print('delta_CL/delta      = ', self.deltaCL/self.deltav)
                print('v_0/u_tau           = ', self.v0/self.Utau)
            else:
                print("no blowing or suction")


    # **********************************************
    #  SKIN FRICTION LAWS ON A FLAT PLATE
    # **********************************************

    def LudwiegTillman(self,Rdelta2,H):
        """
        flat plat turbulent skin friction
        """
        return 0.246*pow(10, -0.678*H)*pow(Rdelta2, -0.268)

    def Felsch(self,Rdelta2,H):
        """
        flat plat turbulent skin friction
        """
        return 0.058*pow(0.93-1.95*np.log10(H),1.705)*pow(Rdelta2,-0.268)

    def SchultzGrunow(self,Rex):
        """
        flat plat turbulent skin friction
        """
        return 0.370*pow(np.log10(Rex),-2.584)

    def Michel(self,Rex):
        """
        flat plat turbulent skin friction
        """
        return 0.0368*pow(Rex,-1./6.)

    def Cf_Power_Law(self,Rex):
        """
        Cf  determined from the power law parameters
        """
        C = pow(self.Cu, -2*self.N/(self.N+3))
        A = (self.N+2)*(self.N+3)/self.N
        N2 = -2/(self.N+3)
        return  2*C*pow(A, N2)*pow(Rex, N2)

    # **********************************************
    #   VARIOUS PLOTS 
    # **********************************************

    def plot_profile(self, mode="Log", legx=r"$y^+$", legy=r"$u^+$"):
        """
        simple plots  for turbulent profile
        
        Params:
            * mode (string): x axis scale : Log or Linear or Laminar
            * legx (string): legend on x axis
            * legy (string): legend on y axis
        """
        fig = plt.figure()
        fig.suptitle('T B L profile', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(legx,fontsize=20)
        ax.set_ylabel(legy,fontsize=20)
        if self.Modele == "Laminar" or mode == "Linear":    
            ax.plot(self.yp, self.up, '-', linewidth=2, color='black')  
        else:
            ax.semilogx(self.yp, self.up, '-', linewidth=2, color='black')
        
        ax.grid()
        plt.show()

    def plot_A_vanDriest(self, v0p, pp, mode="Log", modele='Kays1'):
        """
        A+ law plot w.r.t. v0+ and p+
        """
        fig = plt.figure(self.ifig, figsize=(10,8))
        fig.suptitle(r"Model  %s : $A^+$ w.r.t $v_0^+$ a,d $p^+$"%modele, fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r"$p_0^+$", fontsize=20)
        ax.set_ylabel(r"$A^+/A_0^+$", fontsize=20)
        for v0 in list(v0p):
            if modele == 'Kays1':
                s = self.A_vanDriest_Kays1(v0, pp)
            elif modele == 'Kays2':
                s = self.A_vanDriest_Kays2(v0,pp)
            elif modele == 'Kays3':
                s=self.A_vanDriest_Kays3(v0,pp)
            elif modele == 'Cebeci1':
                s=self.A_vanDriest_Cebeci1(v0,pp)
            else :
                print("model not implemented")
                s = self.Aplus
            if mode == "Linear":    
                ax.plot(pp, [s/self.Aplus for s in s], '-', linewidth=2, label=r"v_0^+ = %f"%(v0))  
            else:
                ax.semilogx(pp, [s/self.Aplus for s in s], '-', linewidth=2, label=r"v_0^+ = %f"%(v0))  
        ax.legend()  
        ax.grid()
        self.ifig += 1

    def plot_B_Kays(self, v0p, pp, mode="Log"):
        """
        d A+ law plot w.r.t. v0+ and p+ in Kays, 1971
        B_Kays3(self, v0plus, Pplus)
        """
        fig = plt.figure(self.ifig, figsize=(10,8))
        fig.suptitle(r"Kays model : $B^+$ w.r.t $v_0^+$ and $p^+$", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r"$p_0^+$", fontsize=20)
        ax.set_ylabel(r"$B^+$", fontsize=20)
        for v0 in list(v0p):
            s = self.B_Kays3(v0,pp)
            if mode == "Linear":    
                ax.plot(pp, [s for s in s], '-', linewidth=2, label=r"v_0^+ = %f"%(v0))  
            else:
                ax.semilogx(pp, [s for s in s], '-', linewidth=2, label=r"v_0^+ = %f"%(v0))  
        ax.legend()  
        ax.grid()
        self.ifig += 1

    def plot_Constante_A_B(self):
        """
        constants A+ and B+ plot w.r.t kappa
        """
        fig = plt.figure(self.ifig, figsize=(10, 8))
        fig.suptitle(r"constants $A_0^+$ and $B_0^+$ w.r.t Karman's constant $\kappa$", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r"$\kappa$", fontsize=20)
        ax.set_ylabel(r"$A^+,B^+$", fontsize=20)
        kappa=np.linspace(0.39, 0.45, 51)
        ax.plot(kappa, self.Aplus_function_kappa(kappa), 'r-', linewidth=2, label=r"A^+")
        ax.plot(kappa, self.Bplus_function_kappa(kappa), 'k--', linewidth=2, label=r"B^+")
        ax.legend()  
        ax.grid()   
        self.ifig += 1

    # **********************************************
    #   MISC
    # **********************************************

    def set_title(self):
        """Set a title of a case"""
        print('#',60*'*')
        print('# %s'%(self.name))
        print('#',60*'*','\n')
