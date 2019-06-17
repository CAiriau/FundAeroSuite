#!/bin/python
# -*- coding: utf-8 -*-
"""
  Main program for  conformal mapping for airfoil, Joukowski transform
  mappings file, with the class
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 

import os
import numpy                as np
import numpy.ma             as ma
import matplotlib.pyplot    as plt
from scipy                  import integrate
from scipy.interpolate      import interp1d
from scipy.interpolate      import InterpolatedUnivariateSpline
 

fz=16   # taille des fontes pour les légendes

#*******************************************
def set_parameters_airfoil():
#*******************************************
    """
    Set default parameters

    """
    prm={}  # Dictionnary of the parameters  

    # physical parameters
    prm["probleme"]     = "portant"     # "portant" "portant+épais"
    prm["Uinf"]         = 1.            # vitesse à l'infini amont en m/s
    prm["rho"]          = 1.3           # masse volumique en kg/m^3
    prm["alpha"]        = 0.0           # incidence en degrés
    prm["Umax"]         = 10            # filtrage de la vitesse sur les points singuliers
    prm["KpLim"]        = [-3,1]            # Pour la visualisation des Kp

    prm["airfoil"]      = ["parabole","analytique","h(x)"]  # cas de l'étude, nom du profil : "parabole", "famille_pq" "double courbure"
                                        # airfoil[0]: "parabole" ou "double pointe" ou "naca4"
                                        #
    prm["npt"]          = 201           # nombre de points pour définir la cambrure ou l'extrados
    prm["eps"]          = 0.1           # petit paramètre pour les définitions des profils
    prm["chord"]        = 1.            # corde du profil    
    prm["nFourier"]     = 5             # nombre de coefficients de la série de Fourier
    # plot parameters
    prm["plot_squelette"] = True        # pour tracer le squelette
    prm["plot_Kp"]      = True          # pour tracer le coefficient de pression Kp sur le corps
    prm["plot_airfoil"] = True          # pour dessiner le profil
    prm["plot_velocity"]= False         # pour dessiner U/Uinf sur le profil
    prm["adimChord"]    = True          # vrai : profil sans dimension
    prm["ordre"]        = 2             # ordre du schéma d'intégration (2 ou 4)
    prm["xsize"]        = 10            # largeur des figures  
    prm["ysize"]        = 8.5           # hauteur des figures  
    eps,q,CL=0.1,7/8.,0.35
    #prm["AirfoilParam"] = [eps, -CL/(np.pi*eps*(3./4-q)) ,q] #[epsilon, p, q]
    p = CL*8/(np.pi*3*eps)
    prm["AirfoilParam"] = [eps, p ,q]
    # options diverses

    return prm


class LinearTheoryAirfoil(object):
    """
    airfoil designed by a conformal mapping
    """
    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#',60*'*')
        print('# %s'%("Linear Theory for 2D section"))
        print('#',60*'*','\n')

        self.probleme   = prm["probleme"] # "portant" "portant+épais"
        # physical parameters
         
        self.Uinf       = prm["Uinf"]   # vitesse à l'infini amont en m/s
        self.rho        = prm["rho"]    # masse volumique en kg/m^3
        self.alpha      = prm["alpha"]  # incidence en degrés
        self.Gamma      = 0.            # circulation à calculée 
        self.Kp         = 0.            # coefficient de pression
        self.Kp_s       = 0.            # coefficient de pression du profil squelettique
        self.chord      = prm["chord"]  # corde du profil

        self.camberMax  = 0.            # cambrure relative maximale
        self.thicknessMax  = 0.         # épaisseur relative maximale
        self.alpha_0    = 0.            # incidence de portance nulle
        self.alpha_ad   = 0.            # incidence de portance nulle
        self.CL         = 0.            # coefficient de portante
        self.xF         = 0.            # position du foyer aérodynamique / corde
        self.CmF        = 0.            # coefficient du moment de tangage au foyer
        self.CmBA       = 0.            # coefficient du moment de tangage au bord d'attaque
        self.delta_s    = 0.            # loi de cambrure du profil squelettique
        self.y_s        = 0.            # équation du squelette
        self.y_e        = 0.            # équation pour l'épaisseur
        self.xe         = 0.            # extrados si données du profil numériques
        self.ye         = 0.            # extrados si données du profil numériques
        self.xi         = 0.            # intrados si données du profil numériques
        self.yi         = 0.            # intrados si données du profil numériques
        self.Lref       = 1.            # longueur de référence pour les graphiques

        self.AirfoilParam=  prm["AirfoilParam"] #[epsilon, p, q]

        self.Umax       = prm["Umax"]    # filtrage de la vitesse sur les points singuliers
        self.KpLim      = prm["KpLim"]  # Pour la visualisation des Kp
        

        self.airfoil    = prm["airfoil"]  # cas de l'étude, nom du profil, "analytique" ou "numérique", 3eme paramètre : "h(x)","y+-",""
        self.npt        = prm["npt"]       # nombre de points pour définir la cambrure ou l'extrados
        self.eps        = prm["eps"]       # petit paramètre pour les définitions des profils
        self.nFourier   = prm["nFourier"]  # nombre de coefficients de la série de Fourier
        self.ordre      = prm["ordre"]     # ordre du schéma d'intégration (2 ou 4)         
                 

        # plot parameters
        self.plot_squelette=prm["plot_squelette"] # pour tracer le squelette

        self.plot_Kp      = prm["plot_Kp"]        # pour tracer le coefficient de pression Kp sur le corps
        self.plot_airfoil = prm["plot_airfoil"] # pour dessiner le profil
        self.plot_velocity= prm["plot_velocity"]# pour dessiner U/Uinf sur le profil
        self.xsize        = prm["xsize"]        # largeur des figures  
        self.ysize        = prm["ysize"]        # hauteur des figures  
        self.adimChord    = prm["adimChord"]    # vrai : profil sans dimension
        self.xlabel       = "x"
        print("initialization done")

    


    def run_theory2D(self):
        """
        On applique la théorie sur un profil donné
        """
        if self.adimChord : 
            self.Lref=self.chord
            self.xlabel=r"$x/c$"


        self.theta      = np.linspace(0,np.pi,self.npt)
        self.x          = self.chord/2*np.cos(self.theta)
        
        if self.airfoil[0]=="naca4":
            self.naca4_analytique()
        self.solution_analytique()
      
        if self.airfoil[2]=="y+-": 
            self.airfoil_geometry()
            self.ordre=2
        
       
            
        
        if self.probleme=="portant":
            if  self.airfoil[1]!="numérique": 
                self.delta_s    = self.camber(self.x)
                self.camberline()
            self.SerieFourier()
            for i in range(self.nFourier):
                print("a[%i]                                     : %f"%(i,self.a[i]))
            self.Aerodynamics()
            self.Kp_portant()
        self.graphiques()

        
    def naca4_analytique(self):
        """
        Initialisation des paramètres pour le NACA 4
        """
        number=self.AirfoilParam[0]
        m = float(number[0])/100.0
        p = float(number[1])/10.0
        t = float(number[2:])/100.0
        self.AirfoilParam[1:]=[m,p,t]
        print("Paramètres du profil NACA 4               : m = %f, p = %f, t = %f"%(m,p,t))
        

    def solution_analytique(self):
        """
        solution analytique si elle est implémentée
        """
        print(self.airfoil[0])
        if self.airfoil[0]=='famille_pq':
            [eps,p,q] = self.AirfoilParam
            print("pour validation, résultats analytiques :")
            print ("eps = %f, p = %f, q= %f"%(eps,p,q))
            print("a0 = %f, a1 = %f, a2 = %f Cl = %f alpha_0 = %f° "%(p*eps/2,p*eps*(1-2*q),3*p*eps/4,-np.pi*p*eps*(3/4-q),np.rad2deg(p*eps/2*(3/4-q))))
        elif self.airfoil[0]=='double courbure':
            print("pour validation, résultats analytiques :")
            print("a0 = %f, a1 = %f, a2 = %f Cl = %f alpha_0 = %f° "%(-2*self.eps,0,-3*self.eps,np.pi*self.eps,np.rad2deg(-self.eps/2)))
        elif self.airfoil[0]=='naca4':
            [m,p]=[self.AirfoilParam[1],self.AirfoilParam[2]]
            if m!=0:
                theta_c=np.arccos(2*p-1)
                a_0 =4* m/(np.pi* p**2* (1-p)**2)* (1-2*p)* (2 *np.sqrt(p*(1-p))+ (1-2*p)*theta_c)  -4 *m/p**2* (1-2*p)-4*np.deg2rad(self.alpha)  
                a_1 = 2*m/(np.pi* p**2* (1-p)**2)* (1-2*p)* (2 *(1-2*p)*np.sqrt(p*(1-p))+theta_c)   - 2* m /p**2 
                a_2 = 32* m* (1-2*p)*np.sqrt(p*(1-p)) / (3*p*np.pi*(1-p))
                print("a0 = %f, a1 = %f, a2 = %f Cl = %f alpha_0 = %f° "%(a_0,a_1,a_2,-np.pi/2*(a_0+a_1), np.rad2deg(a_0+a_1)/4))
    
    def thicknessLaw(self):
        """
        Loi d'épaisseur, EN COURS ...
        """
        if self.airfoil[0]=="double courbure":
            eta=self.x/self.chord
            mu=self.AirfoilParam[0]
            self.y_e=mu*self.eps*self.chord/2*(1-eta)*np.sqrt(1-eta**2) 

        else:
            self.y_e= np.zeros(len(self.x))  

    def thicknessSlope(self,x):
        """
        pente de la loi d'épaisseur,  EN COURS ...
        """
        if self.airfoil[0]=="double courbure":
            eta=self.x/self.chord
            mu=self.AirfoilParam[0]
            return -mu*self.eps*(1+2*eta)*np.sqrt((1-eta)/(1+eta)) 
 
    def camberline(self):
        """
        Equation du squelette analytique
        """    
        if self.airfoil[0]=='parabole':
            a = self.chord/4
            self.y_s= self.eps/a**2/2*(2*a-self.x)*(2*a+self.x)
        elif self.airfoil[0]=="famille_pq":
            [eps,p,q] = self.AirfoilParam
            eta=self.x/self.chord
            self.y_s = eps*p*self.chord*(eta+1/2)*(1/2-eta)*(q-eta-1/2)
        elif self.airfoil[0]=="double courbure":
            eta=self.x/self.chord
            self.y_s=self.eps*self.chord/2*eta*(1-eta**2) 
        elif self.airfoil[0]=="naca4":
            
            self.y_s=np.zeros(self.npt)
            [m,p,t]=self.AirfoilParam[1:]
            if p*m != 0 :
                for k in range(self.npt):
                    x=self.x[k]+0.5*self.chord
                    if x/self.chord <= p:
                        self.y_s[k]=m/p**2*x*(2*p-x/self.chord)
                    else:
                        self.y_s[k]=m/(1-p)**2*(self.chord-x)*(1-2*p+x/self.chord)

        else:
            self.y_s= np.zeros(self.npt)  

    def camber(self,x):
        """
        Equation de la cambrure
        """
        if self.airfoil[0]=='parabole':
            return -self.eps*x/(self.chord/4)**2-np.deg2rad(self.alpha)
        elif self.airfoil[0]=="famille_pq":
            [eps,p,q] = self.AirfoilParam
            eta=self.x/self.chord
            return eps*p*(eta*(1-2*q) -1/4+3*eta**2)-np.deg2rad(self.alpha)
        elif self.airfoil[0]=="double courbure":
            eta=self.x/self.chord
            return self.eps*(1-3*eta**2/(self.chord/2)**2) -np.deg2rad(self.alpha)
        elif self.airfoil[0]=="naca4":
            delta=np.zeros(self.npt)
            [m,p,t]=self.AirfoilParam[1:]
            if p*m != 0 :
                for k in range(self.npt):
                    x=self.x[k]+0.5*self.chord
                    if x/self.chord <= p:
                        delta[k]=2*m/p**2* (p-x/self.chord)-np.deg2rad(self.alpha)
                    else:
                        delta[k]=2*m/(1-p)**2*(p-x/self.chord)-np.deg2rad(self.alpha)
            else:
                print("profil symétrique")
                delta=-np.deg2rad(self.alpha)*np.ones(self.npt)
            return delta
        else:
            return np.zeros(len(x))    
    
    def integrand(self,theta,n=0):
        """
        intégrand 
        """
        delta=camber(self,self.chord/2*np.cos(theta))
        return delta*np.cos(n*theta)

    def SerieFourier(self):
        """
        n premiers coefficients a_k = 4/pi int_0^pi delta(theta)cos(n theta) dtheta
        """
        self.a=np.zeros(self.nFourier)
        if self.ordre==2:
            print('ordre 2')
            for i in range(self.nFourier):
                self.a[i]=4/np.pi*np.trapz(self.delta_s*np.cos(i*self.theta),self.theta)
        elif ordre==4:
            print('ordre 4')
            for i in range(self.nFourier):
                self.a[i],res=4/np.pi*integrate.quad(self.integrand,0,np.pi,args=(i,))
        else:
            raise ValueError('Mauvais ordre dans SerieFourier')
        
    def Aerodynamics(self):
        """
        Calcul des coefficient aérodynamiques
        """
        self.CL         = -np.pi/2 *(self.a[0]+self.a[1])
        self.alpha_0    = (self.a[0]+self.a[1]) / 4
        self.CmF        = -np.pi/8 *(self.a[2]+self.a[1])
        self.alpha_ad   = self.a[0]/ 4
        self.CmBA       = -np.pi/8 *(self.a[0]+2*self.a[1]+self.a[2])
        
        if self.alpha != 0 :
            print("ATTENTION :  \t le calcul de alpha_0 et alpha_d n'est pas correct")
            print("\t \t car l'incidence n'est pas nulle")

        print("alpha                                     : %f °"%(self.alpha))
        print("CL                                        : %f"%(self.CL))
        if self.alpha==0:
            print("CL = 2 pi (alpha - alpha_0)               : %f"%(2*np.pi*(np.deg2rad(self.alpha)-self.alpha_0)))
            print("alpha_0                                   : %f °"%(np.rad2deg(self.alpha_0)))
            print("alpha_adaptation                          : %2.15f °"%(np.rad2deg(self.alpha_ad)))
        print("Cm au foyer                               : %f"%(self.CmF))
        print("Cm au bord d'attaque                      : %f"%(self.CmBA))


    def Kp_portant(self):
        """
        Calcul du Kp pour le problème portant
        """
        self.Kp_s=np.zeros(self.npt)
        for i in np.arange(1,self.nFourier):
            self.Kp_s+=self.a[i]*np.sin(i*self.theta)
            print("a %i %f"%(i,self.a[i]))

        if abs(self.a[0]) >= 1e-7 :
            print("prise en compte de a_0                   : %e"%(self.a[0]))
            self.Kp_s += self.a[0]/2*np.tan(self.theta/2)
            self.Kp_s[-1]==float("inf")
     
        #print("Kp_s = ",self.Kp_s)

    def graphiques(self):
        """
        différents graphiques
        """

        if self.plot_squelette:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Profil du problème '+self.probleme, fontsize=fz, fontweight='bold')   
            plt.xlabel(self.xlabel,fontsize=fz)

            if self.probleme=="portant":
                print("len x, y_s :",len(self.x),len(self.y_s),self.Lref)
                plt.plot(self.x/self.Lref,self.y_s/self.Lref,color='red',label=r'$y_s$')
                plt.ylabel(r'$y_s$',fontsize=fz)
            plt.axis("equal")    
            plt.grid()
            plt.legend()
            plt.xlim(-0.55*self.chord/self.Lref,0.55*self.chord/self.Lref)
            
        
        if self.plot_Kp:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Coefficient de pression du problème '+self.probleme, fontsize=fz, fontweight='bold')   
            plt.xlabel(self.xlabel,fontsize=fz)

            if self.probleme=="portant":
                plt.plot(self.x/self.Lref,self.Kp_s,color='red',label=r'$Kp_s$')
                plt.ylabel(r'$Kp_s$',fontsize=fz)

            plt.grid()
            plt.legend()
            plt.ylim(self.KpLim[0],self.KpLim[1])
        plt.show()

    def airfoil_geometry(self):
        """
        Définir la cambrure à partir des données du profil
        en entrée xplus, yplus,xmoins, ymoins (+ et - pour extrados et intrados)
        """
        # il faut interpoler le profil dans le nouveau maillage x
        ye = interp1d(self.xe, self.ye, kind='linear')
        yi = interp1d(self.xi, self.yi, kind='linear')
         
        print("extrados: min et max de x (profil)        : %f et %f"%(np.min(self.xe),np.max(self.xe)))
        print("intrados: min et max de x (profil)        : %f et %f"%(np.min(self.xi),np.max(self.xi)))
        print("min et max de x (nouveau)                 : %f et %f"%(np.min(self.x),np.max(self.x)))
        # je vérifie l'incidence géométrique du profil
        incidence_profil       = np.arctan(self.ye[0]-self.ye[-1])
        print("incidence initiale                        : %f °, yBF = %f yBA = %f x_BF = %f x_BA = %f"%(np.rad2deg(incidence_profil),self.ye[0],self.ye[-1],self.xe[0],self.xe[-1] ))
        print("corde                                     : %f"%(self.chord))
        # profil avec le nouveau vecteur x :
        y_plus,y_moins = ye(self.x),yi(self.x)
        self.y_s,self.y_e=(y_plus+y_moins)/2,(y_plus-y_moins)/2
        # petit problème pour les profils non symétriques, il faut mieux extrapoler le y_s au bord de fuite
        #ys_tmp= InterpolatedUnivariateSpline(np.flipud(self.x[-4:-1]),np.flipud(self.y_s[-4:-1]),k=2)
        ys_tmp= InterpolatedUnivariateSpline(self.x[-4:-1],self.y_s[-4:-1],k=2)
        self.y_s[-1]=ys_tmp(self.x[-1])    
    

        self.epaisseur,ind_e = np.amax(self.y_e),np.argmax(self.y_e)
        self.cambrure,ind_c  = np.amax(self.y_s),np.argmax(self.y_s)
        print(" cambrure / c                             : %5.4f à x(%4d)/c = %2.4f soit %f pourcents de corde"%(self.cambrure/self.Lref,ind_c,\
                                                                   self.x[ind_c]/self.Lref,100*(self.x[ind_c]/self.Lref+0.5)))
        print(" epaisseur / c                            : %5.4f à x(%4d)/c = %2.4f soit %f  pourcents de corde"%(2*self.epaisseur/self.Lref,\
                                                                   ind_e,self.x[ind_e]/self.Lref,100*(self.x[ind_e]/self.Lref+0.5)))

        # calcul des pentes :
        dy=self.y_s[1:]-self.y_s[0:-1]
        dx=self.x[1:]-self.x[0:-1]
        print("len y = ",len(dy), "len ys = ",len(self.y_s),len(self.y_s[:-1]))
        print("dy[0] =", dy[0],self.y_s[1]-self.y_s[0])
       
         
        self.delta_s=np.zeros(self.npt)
        self.delta_s[0:-1]=dy/dx
        #self.delta_s[-2:]= self.delta_s[-3]
        #self.delta_s[-1]= self.delta_s[-2] 
        # petit problème pour les profils non symétriques, il faut mieux extrapoler le delta_s au bord de fuite
        #ys_tmp= InterpolatedUnivariateSpline(np.flipud(self.x[-4:-1]),np.flipud(self.y_s[-4:-1]),k=2)
        deltas_tmp= InterpolatedUnivariateSpline(self.x[-4:-1],self.delta_s[-4:-1],k=2)
        self.delta_s[-1]=deltas_tmp(self.x[-1])    

        self.delta_s-=np.deg2rad(self.alpha)
        
        self.PlotDelta_s()

        if self.plot_airfoil:
            self.PlotAirfoil()
        #print(self.delta_s)
        #print(dx[-2:],dy[-2:])

    def PlotAirfoil(self):
        """
        Dessin du profil
        """
        symbole=False

        plt.figure(figsize=(18, 4))
        plt.title("Profil")
        plt.xlabel("x",fontsize=fz)
        plt.ylabel("y",fontsize=fz)
        if symbole:
            plt.plot(self.x,self.y_s,'bo',label=r'$y_s$')
            plt.plot(self.x,self.y_e,'b--',label=r'$y_e$')
            plt.plot(self.xe,self.ye,'ko',label=r'ext.')
            plt.plot(self.xi,self.yi,'rs',label=r'int.')
        else:
            plt.plot(self.x,self.y_s,'k-',label=r'$y_s$')
            plt.plot(self.x,self.y_e,'r-',label=r'$y_e$')
            plt.plot(self.xe,self.ye,'k--',label=r'ext.')
            plt.plot(self.xi,self.yi,'r--',label=r'int.')
        plt.grid()
        plt.axis("equal")
        print("x = ",self.x[:2],self.x[-3:])
        print("ys = ",self.y_s[:2],self.y_s[-3:])
        print("deltas = ",self.delta_s[:2],self.delta_s[-3:])
        plt.legend()
        plt.show()        

    def PlotDelta_s(self):
        """
        Dessin de delta_s
        """
        plt.figure(figsize=(14, 6))
        plt.title("loi de pente du profil squelettique")
        plt.xlabel("x",fontsize=fz)
        plt.ylabel(r"$\delta_s$",fontsize=fz)
        plt.plot(self.x,np.rad2deg(self.delta_s),'bo')
        #plt.plot(self.x,self.y_s,'rs')
        plt.grid()
        plt.xlim(-0.55,0.55)
        plt.ylim(-20,20)
        plt.show()        
