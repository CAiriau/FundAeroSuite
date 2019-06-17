#!/bin/python
# -*- coding: utf-8 -*-
"""
  Profil de la couche limite turbulente
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize             import fsolve
from scipy                      import stats
from Tools.misc                 import *
from CompressibleFlow.fonctions import * 
from scipy.interpolate          import interp1d
from scipy.special              import erfc
from scipy.special              import erf
from scipy.integrate            import odeint

 

#*******************************************
def set_parameters():
#*******************************************
    """
     Fixer les paramètres par défaut
    """
    prm={}  # Dictionnaire des paramètres  
    prm["nom"]      = "Profil initial"
    prm["ypMin"]    = 0.001
    prm["ypMax"]    = 250
    prm["n"]        = 101
    prm["kappa"]    = 0.41
    prm["Utau"]     = 0
    prm["Rtau"]     = 1000
    prm["Cf"]       = 0.003
    prm["Modele"]   = "Linear"
    #        "Laminaire" : Ecoulement laminaire avec aspiration
    #        "SCVonly"   : approximation pour la sous couche visqueuse
    #        "CIonly"    : approximation pour la couche interne
    #        "SCVlinear" : u+=y+
    #        "Linear"    : L+ = kappa y+   pour tout y+
    #        "Mixed"     : "Linear" si  y+ <= ypScv sinon "Michel"
    #        "Michel"    : Longueur de Mélange de Michel
    #        "Puissance" : Loi puissance
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
    #               = 0     fixé=A+    |  pas calculé
    #               = 1     recalculé  |  pas calculé
    #               = 2     fixé=A+    |  recalculé
    #               = 3     recalculé  |  recalculé
    prm["BplusOpt"] = 0  # même tableau que pour AplusOpt
           
    prm["AplusModel"]= ""
    prm["Bplus"]    = 36.0
    prm["Pplus"]    = 0.0
    prm["A"]        = 1/prm["kappa"]
    prm["B"]        = 0
    prm["Cu"]       = 0
    prm["Rex"]      = 0
    prm["N"]        = 7.0   # profil à la puissance 1/N
    prm["ShowParameters"] = False 
    prm["com"]      =" "
    return prm


# Pour calculer automatiquement A+ ou B+, mettre les valeurs initiales à 0.

class ProfilTurbulent(object):
    """
    propriétés et caractéristiques du profil
    """


    def __init__(self, option, **kwargs):
        self.option = option
        self.name   = option['nom']
        self.ypMin  = option['ypMin']               # y+ début
        self.ypMax  = option['ypMax']               # y+ fin
        self.n      = option['n']                   # nombre de points
        self.kappa  = option['kappa']               # constante de Karman    
        self.Utau   = option['Utau']                # vitesse de frottement U_tau
        self.Rtau   = option['Rtau']                # Reynolds basé sur la vitesse de frottement
        self.Cf     = option['Cf']                  # Coefficient de frottement
        self.Modele = option['Modele']              # type de modèle
        self.Damping = option["Damping"]            # Application de la Correction pour l'amortissmeent : No,  VanDriest ,  Kays
        self.Aplus  = option["Aplus"]               # constante dans la loi de Van Driest
        self.AplusOpt = option["AplusOpt"]          # option pour calculer A+ en fonction de v0+ et p+ ou pas.
        self.Bplus  = option["Bplus"]               # constante dans la loi d'amortissement de Kays
        self.BplusOpt = option["BplusOpt"]          # option pour calculer B+ en fonction de v0+ et p+ ou pas.
        self.AplusModel=option["AplusModel"]        # Modèle pour le calcul de A+ : Kays1, Kays2, Kays3, Cebeci1
        self.pplus  = option["Pplus"]               # gradient de pression
        self.v0p    = option['v0p']                 # vitesse aspiration / vitesse de frottement
        self.v0     = ''                            # vitesse d'aspiration (>0)
        self.Rv0    = option['Rv0']                 # Reynolds  v_0 delta/nu
        self.ypScv  = option['ypScv']               # y+ limite de la sous couche visqueuse théorique
        self.Ue     = option['Ue']                  # vitesse extérieure
        self.Uep    = option['Uep']                 # vitesse extérieure réduite
        self.nu     = option['nu']                  # viscosité cinématique
        self.rho    = option['rho']                 # masse volumique
        self.mu     = self.nu*self.rho              # viscosité dynamique        
        self.A      = option["A"]                   # pente de la loi log
        self.B      = option["B"]                   # ordonnée à l 'origine de la loi log
        self.N      = option["N"]                   # pour la loi puissance, profil à la puissance 1/N
        self.Cu     = option["Cu"]                  # coefficient de la loi puissance
        self.yv     = ''                            # y_v délimitant la sous couche visqueuse avec aspiration
        self.yp     = ''                            # y+
        self.up     = ''                            # u+
        # pour les tests
        self.y_vp   = ''                            # y_v^+
        self.u_vp   = ''                            # u_v^+
        self.y_sv   = 0
        self.up_sv  = 0
        self.plot_sv= False                         # pour dessiner y_v^+ u_v^+
        self.dup_dyp= ''                            # du+/dy+ (sortie)
        self.Rdelta = ''                            # Reynolds basé sur l'épaisseur de la couche limite
        self.deltaCL= ''                            # épaisseur de la couche limite
        self.deltav = ''                            # épaisseur nu/v0 (aspiration)
        self.Lplus  = ''                            # longueur de mélange
        self.k1     = ''                            # delta1/delta
        self.k2     = ''                            # delta2/delta
        self.H      = ''                            # facteur de forme
        self.Rex    = option["Rex"]                 # Re_x
        self.erreur = False
        self.ShowParameters = option['ShowParameters']
        self.ifig   = 0                             # Numérotation des figures
        self.com    = option["com"]                 # commentaire

        self.set_title()

        if self.Modele=="Laminaire" and self.Rv0 == 0:
            self.erreur=True
            print("Pas d'implémentation de profil laminaire sans aspiration !")
             
        
        if self.v0p!= 0:
            print("profil avec aspiration")

        if self.pplus!= 0:
            print("profil avec gradient de pression")

        if self.Damping != "No":
            print("Besoin de recalculer A+ ou B+")
            if self.Damping=="Kays":
                print('B0+ (entré)         = ',self.Bplus)
                if self.BplusOpt==0 or self.BplusOpt==2 :
                    print('B+_0 = B+ entré')
                else :
                    self.Bplus=self.Bplus_function_kappa(self.kappa)
                    print('B0+ (calculé pour kappa = %4.3f) = %4.3f'%(self.kappa,self.Bplus))
                if self.BplusOpt > 1 :
                    self.Bplus=self.B_Kays3(self.v0p,self.pplus)
                    print('B+ (v0+ = %4.3f, p+ = %4.3f)  = %4.3f'%(self.v0p,self.pplus,self.Bplus)) 

            if self.Damping=="VanDriest":
                print('A0+ (entré)         = ',self.Aplus)
                if self.AplusOpt==0 or self.AplusOpt==2 :
                    print('A+_0 = A+ entré')
                else :
                    self.Aplus=self.Aplus_function_kappa(self.kappa)
                    print('A0+ (calculé pour kappa = %4.3f) = %4.3f'%(self.kappa,self.Aplus)) 
                if self.AplusOpt > 1 :
                    if self.AplusModel=='Kays1':
                        s=self.A_vanDriest_Kays1(self.v0p,self.pplus)
                    elif self.AplusModel=='Kays2':
                        s=self.A_vanDriest_Kays2(self.v0p,self.pplus)
                    elif self.AplusModel=='Kays3':
                        s=self.A_vanDriest_Kays3(self.v0p,self.pplus)
                    elif self.AplusModel=='Cebeci1':
                        s=self.A_vanDriest_Cebeci1(self.v0p,self.pplus)
                    self.Aplus=s
                print('A+ (v0+ = %4.3f, p+ = %4.3f)  = %4.3f'%(self.v0p,self.pplus,self.Aplus)) 
        if (self.Modele=="SCVonly" or self.Modele=="CIonly") and (self.v0p!=0 or self.pplus!=0): self.com="Faux"

    def calcul_profil(self):
        """
        Génération de la grille et calcul du profil
        """
        if not self.erreur: 
            self.genere_grille()
            if self.Modele=='Puissance':
                self.profil_puissance()
            else:
                print('Autre profil que la loi Puissance')
                self.calcul_parametres()
                self.genere_profil()
            #print('U = ',self.up)
            print("Rdelta                         = %f "%(self.Rdelta))
            print("Max u+                         = %f "%(max(self.up)))

    # **********************************************
    #  LOI PUISSANCE pour les PROFILS DE VITESSE
    # **********************************************

    def calcul_Rex(self):
        """
        Nombre de Reynolds basé sur l'abcisse x
        """
        self.Rex=pow(self.Rtau,1+3/self.N)*self.Cu**3*self.N/((self.N+2)*(self.N+3))

    def calcul_Cf(self):
        """
        Coefficient de frottement
        """
        self.Cf=2/(self.Cu**2*pow(self.Rtau,2/self.N))

    def calcul_Cu(self):
        """
        Coefficient de proportionnalité de la loi puissance
        """
        if self.Cf != 0:
            self.Cu=pow(self.Rtau,1/self.N)*np.sqrt(2/self.Cf)
        elif self.Rex != 0:
            tmp=self.Rex*(self.N+2.)*(self.N+3.) / (self.N*pow(self.Rtau,1.+3./self.N))
            #tmp=(self.N+2.)*(self.N+3.) / (self.N*pow(self.Rtau,1.+3./self.N))
            #print("Rex = ",self.Rex)
            self.Cu=pow(tmp,1./3.)
        else:
            print("Problème dans le calcul de Cu")

    def calcul_Cu_from_Cf(self):
        """
        on connaît : Rex, Cf, N => Cu, Rtau
        """
        tmp=(self.N+2)*(self.N+3)*self.Rex/self.N*pow(self.Cf/2,(self.N+3)/2)
        self.Cu= pow(tmp,-1/self.N)

    def calcul_Rtau(self):
        """
        Calcul du Rtau à partir du Cf, connaissant N, Rex, Cu
        """
        self.Rtau=pow(self.Cf/2,-self.N/2)/pow(self.Cu,self.N)
        self.Rtau=pow(np.sqrt(2/self.Cf)/self.Cu,self.N)


    def profil_puissance(self):
        """
        Calcul avec la loi en puissance
        """

        def f_N(x):
            """
            Equations à résoudre si N est inconnu
            """
            return self.Rex-pow(self.Rtau,1+3/x)*self.Cu**3*x/((x+2)*(x+3))

        print("Calcul des caractéristiques de la loi puissance")
        
        if self.N==0:
            print("Calcul de N")
            self.N=fsolve(f_N,7)
            print('nouvelle valeur de N : ',self.N)
            self.calcul_Cf()
        elif self.Cf*self.Rtau != 0:                # Rtau et Cf connus > calcul de Cu et Rex
            self.calcul_Cu()
            self.calcul_Rex()
        elif self.Rtau*self.Cu != 0:                # Rtau et Cu connus > calcul de Cf et Rex
            self.calcul_Cf()
            self.calcul_Rex() 
        elif self.Rtau*self.Rex != 0:               # Rtau et Rex connus > calcul de Cf et Cu
            print("Rtau et Rex donnés")
            self.calcul_Cu()
            self.calcul_Cf()
        else:
            print('Configuration à problème, vérifier les entrées pour Cu, Rtau, Cf et Rex')
            print('Cu et Rtau fixés arbitrairement')
            self.Cu,self.Rtau=7,1000
            self.calcul_Cf()
            self.calcul_Rex() 

        self.k1,self.k2=1/(self.N+1),self.N/((self.N+1)*(self.N+2))
        self.Uep=self.Cu*pow(self.Rtau,1/self.N)
        self.Rdelta=self.Rtau*np.sqrt(2/self.Cf)
        self.Rex=pow(self.Rtau,1+3/self.N)*self.Cu**3*self.N/((self.N+2)*(self.N+3))
        self.H=self.k1/self.k2
        print("Rtau                           = %f "%(self.Rtau))
        print("Uep                            = %f "%(self.Uep))
        print("Cf                             = %f "%(self.Cf))
        print("Cu                             = %f "%(self.Cu))
        print("Rdelta                         = %f "%(self.Rdelta))
        print("Re_x                           = %e "%(self.Rex))
        print("delta_1/delta                  = %f "%(self.k1))
        print("delta_2/delta                  = %f "%(self.k2))
        print(" H                             = %f "%(self.H))
       

        self.up=self.Cu*pow(self.yp,1/self.N)
        self.dup_dyp=self.up/self.N/self.yp
        # caractéristiques de la couche limite
        C=pow(self.Cu,-2*self.N/(self.N+3))
        A=(self.N+2)*(self.N+3)/self.N
        N1=(self.N+1)/(self.N+3)

        N2=-2/(self.N+3)

        print("delta/x     = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A,N1), -2/(self.N+3),-2, self.N+3 ))
        print("Rdelta      = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A,N1), N1      ,self.N+1,  self.N+3     ))
        print("Rdelta_1    = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A,N1)*self.k1, N1 , self.N+1,  self.N+3    ))
        print("Rdelta_2    = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A,N1)*self.k2, N1 , self.N+1,  self.N+3   ))
        print("Cf/2        = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(C*pow(A,N2)        , N2 ,-2, self.N+3   ))
        print("Rtau        = %f  Re_x^ %+5.3f \t ou %+2.1f / %+2.1f"%(pow(self.Cu,-3*self.N/(self.N+3))*pow(A, self.N/(self.N+3) ), self.N/(self.N+3) , self.N,  self.N+3 ))

        
        self.lois_puissance()
        
        print("Loi de Michel               Cf = %f"%(self.Michel(self.Rex)))
        print("Loi de Schultz et Grunow    Cf = %f"%(self.SchultzGrunow(self.Rex)))
        Rdt2=self.k2*self.Rdelta
        print("Loi de Ludwieg et Tillman   Cf = %f"%(self.LudwiegTillman(Rdt2,self.H)))
        print("Loi de Felsch et al         Cf = %f"%(self.Felsch(Rdt2,self.H)))

    def lois_puissance(self,display=True):
        """
        Lois qui résultent de la loi puissance
        """
        C=pow(self.Cu,-2*self.N/(self.N+3))
        A=(self.N+2)*(self.N+3)/self.N
        N1=(self.N+1)/(self.N+3)
        N2=-2/(self.N+3)

        print("delta/x     = %f "%(C*pow(A,N1)* pow(self.Rex, -2/(self.N+3)) ))
        print("Rdelta      = %f "%(C*pow(A,N1)* pow(self.Rex, N1)            ))
        print("Rdelta_1    = %f "%(C*pow(A,N1)*self.k1*pow(self.Rex,N1)      ))
        print("Rdelta_2    = %f "%(C*pow(A,N1)*self.k2*pow(self.Rex,N1)      ))
        print("Cf          = %f "%(2*C*pow(A,N2)*pow(self.Rex, N2)           ))
        print("Rtau        = %f "%(pow(self.Cu,-3*self.N/(self.N+3))*pow(A, self.N/(self.N+3) )*pow(self.Rex,self.N/(self.N+3))))

    # *********************************************************
    #  GENERATION D'UN PROFIL TURBULENT POUR DIFFERENTS MODELES
    # *********************************************************    

    def genere_grille(self):
        """
        Génération de la grille en y
        """    
        if self.Modele=="Laminaire":
            self.yp=np.linspace(0,self.ypMax,self.n)
            #print("eta = ",self.yp)
        else:
            self.yp=np.exp(np.linspace(np.log(self.ypMin),np.log(self.ypMax),self.n))
            #print("y+ = ",self.yp)

    def genere_profil(self):
        """
        Génération du profil
        """        
        if self.Modele == "Laminaire":
            self.up=1-np.exp(-self.yp)
            self.dup_dyp=np.exp(-self.yp)

        elif self.Modele == "SCVonly":
            print("ATTENTION : ce modèle est faux : on suppose A+ constante")    # alors pourquoi le conserver ?
            self.yv=self.calcul_yv()        # calcul faux car A+(v0+, p+)
            self.ypMax=self.yv
            self.genere_grille() # besoin de le refaire car j'ai changé le minimun
            self.up=(1-np.exp(-self.v0p*self.yp))/self.v0p
            self.dup_dyp=np.exp(-self.v0p*self.yp)
        elif self.Modele == "CIonly":
            print("ATTENTION : ce modèle est faux : on suppose A+ constante")
            self.yv=self.calcul_yv()        # calcul faux car A+(v0+, p+)
            print(" yv = ",self.yv)
            self.ypMin=self.yv
            self.genere_grille() # besoin de le refaire car j'ai changé le minimun
            m=-self.v0p/(2*self.kappa)
            uv=(1-np.exp(-self.v0p*self.yv))/self.v0p
            tmp=m*np.log(self.yp/self.yv)+np.sqrt(1-self.v0p*uv)
            self.up=(1-tmp**2)/self.v0p
            self.dup_dyp=self.up  # c'est faux, c'est pour remplir le vecteur uniquement
            self.calcul_A_B() # constante pour la loi log approchée
        elif self.Modele=="SCVlinear":
            self.up=self.yp    
            self.dup_dyp=np.ones(len(self.yp))
        elif self.Modele=="loi_log":
            self.loi_log()
        else :
            self.profil_turbulent_aspiration()

        print("Max u+                         = %f "%(max(self.up)))
        
        self.y_vp=0
        if self.Rtau==self.ypMax:
            print("Le profil est défini sur toute la couche limite")
            print("On recalcule quelques grandeurs")
            self.Uep=max(self.up)
            self.Rdelta=self.Uep*self.Rtau
            self.Cf=2/self.Uep**2
            print('NOUVELLES SORTIES')
            print('')    
            print('Ue+                 = %f'%(self.Uep))
            print('Rdelta              = %f'%(self.Rdelta))
            print('Cf                  = %f'%(self.Cf))

            if self.v0p!= 0:
                self.plot_sv=True
                self.y_vp=self.yp
                self.u_vp=self.solve_yv_uv(self.yp)      #-self.up[:,0]
                #print(np.where(self.yp>10))
                v=self.u_vp-self.up[:,0]
                xloc=np.where(self.yp>10)
                vloc,yloc=v[xloc],self.yp[xloc]
                #print("yloc = ",yloc)
                #print("vloc = ",vloc)
                i=0
                while vloc[i]>0 and i < len(vloc):
                    i=i+1
                #print("i = ",i,vloc[i],yloc[i-1],yloc[i])
                self.yp_sv=yloc[i-1]-(yloc[i]-yloc[i-1])/(vloc[i]-vloc[i-1])*vloc[i-1]
                self.up_sv=self.solve_yv_uv(self.yp_sv)
                print("limite de la sous couche visqueuse   : y+ = %f, u+= %f"%(self.yp_sv,self.up_sv))


        else:
            print("Le profil n'est pas défini sur toute la couche limite")


    # **********************************************
    #  LIMITE DE LA SOUS COUCHE VISQUEUSE
    # **********************************************
 
    def solve_yv_uv(self,y):
        """
        Délimitation de la sous couche visqueuse
        on trace la fonction pour voir...
        """
        #self.y_vp=np.exp(np.linspace(np.log(self.ypMin),np.log(self.ypMax),self.n))
       
        tmp=np.sqrt(1-self.v0p*self.Uep)+self.v0p/(2*self.kappa)*np.log(self.Rtau/y)
        return (1-tmp**2)/self.v0p
        
         
    # **********************************************
    #  LOIS D'AMORTISSEMENTS
    # **********************************************


    def regression_lineaire(self,filename,skiprows=0):
        """
        Régression linéaire d'une courbe située dans un fichier
        """
        reference_filepath = os.path.join('Livre/Data',filename)
        with open(reference_filepath, 'r') as infile:
            Data=np.loadtxt(infile, dtype=float, unpack=True,skiprows=skiprows)
        slope, intercept, r_value, p_value, std_err = stats.linregress(Data[0],Data[1])
        print('y = %f + %f x , correlation = %f, sigma = %f '%(intercept,slope,r_value,std_err))

        # autre façon de faire, moins d'information
        coef=np.polyfit(Data[0],Data[1],1)
        print("Polynôme : ",coef)
        # p=np.poly1d(coef)   # polynome fonction
        return [intercept,slope]


    def A_B_regression(self):
        """
        Régression linéaire des constantes A+ et B+, courbes de Kays
        Int J. Heat Mass transfer. Vol. 15. pp. 1023-1044. W.M. Kays, 1972

        """
        A=self.regression_lineaire('A_function_K.dat',skiprows=2)
        B=self.regression_lineaire('B_function_K.dat',skiprows=2)

        #return [A,B]

    def plot_loi_Van_Driest(self,xaxis="Log"):
        """
        dessin de la correction de Van Driest
        """
        self.genere_grille()
        self.up=1-np.exp(-self.yp/self.Aplus)
        self.plot_profile(mode=xaxis,legx=r"$y^+$",legy=r"$\chi(y^+)$")

    def Aplus_function_kappa(self,kappa):
        """
        Constante A+=f(kappa) (kappa = constante de Karman)
        après régression de la courbe donnée par Kays
        """
        return -13.047+90.171*kappa

    def Bplus_function_kappa(self,kappa):
        """
        Constante B+=g(kappa) (kappa = constante de Karman)
        après régression de la courbe donnée par Kays
        """
        return -8.0503+102.24*kappa

    def B_Kays3(self,v0plus,Pplus):
        """
        constante pour l'amortissement en fonction du gradient de pression et 
        de la vitesse de transpiration (Kays)
        Int J. Heat Mass transfer. Vol. 15. pp. 1023-1044. W.M. Kays, 1972
        p. 1035
        """
        if isinstance(Pplus,float):
            return self.Bplus/( 9*v0plus+3.35*Pplus/(1+4.0*v0plus)+1.0 )
        else:
            s=[] 
            for p in Pplus:
                sol=self.Bplus/( 9*v0plus+3.35*p/(1+4.0*v0plus)+1.0 )
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
               
            return s        
    

    def A_vanDriest_Kays3(self,v0plus,Pplus):
        """
        constante de van Driest en fonction du gradient de pression et 
        de la vitesse de transpiration (Kays)
        Int J. Heat Mass transfer. Vol. 15. pp. 1023-1044. W.M. Kays, 1972
        p. 1035
        """
        if isinstance(Pplus,float):
            return self.Aplus/( 5.15*(v0plus+5.86*Pplus/(1+5.0*v0plus))+1.0 )
        else:
            s=[] 
            for p in Pplus:
                sol=self.Aplus/( 5.15*(v0plus+5.86*p/(1+5.0*v0plus))+1.0 )
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
               
            return s        


    def A_vanDriest_Kays2(self,v0plus,Pplus):
        """
        constante de van Driest en fonction du gradient de pression et 
        de la vitesse de transpiration (Kays), approche par surface d'ordre 2
        Il marche beaucoup moins bien, limité dans les valeurs de v0+ et p+
        Andersen, Kays et Moffat, JFM 69,part 2, 1975, p.364
        """
        x=np.log(v0plus+0.08)
        a=[[-6.71751, -1.50414],[-4.81589,-1.24311],[-1.27827, -0.388216]]

        if isinstance(Pplus,float):
            y=np.log(Pplus+0.04)
            u=[0,0,0]
            for i in range(3):
                u[i]=a[i][0]+a[i][1]*y
                #print('a[ %i , 0] = %f \t a[ %i , 1] = %f '%(i,a[i][0],i,a[i][1]))
            tmp= u[0]+x*u[1]+x**2*u[2]
            return 24*np.exp(tmp)

        else:

            s=[]
            for p in Pplus:
                y=np.log(p+0.04)
                u=[0,0,0]
                for i in range(3):
                    u[i]=a[i][0]+a[i][1]*y
                    #print('a[ %i , 0] = %f \t a[ %i , 1] = %f '%(i,a[i][0],i,a[i][1]))
                tmp= u[0]+x*u[1]+x**2*u[2]
                sol=24*np.exp(tmp)
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
            return s

    def A_vanDriest_Kays1(self,v0plus,Pplus):
        """
        constante de van Driest en fonction du gradient de pression et 
        de la vitesse de transpiration (Kays)
        Andersen, Kays et Moffat, JFM 69,part 2, 1975, P.367
        """
        if v0plus < 0 : 
            a=9.0
        else :
            a=7.1
        
        if isinstance(Pplus,float):
            b,c=4.25,10.0
            if Pplus > 0  : b,c=2.9,0.0     # b= 2.9 dans l'article, 2 dans le rapport
            return self.Aplus/(a*(v0plus+b*Pplus/(1+c*v0plus))+1.0)

        else:

            s=[]
            for p in Pplus:
                b,c=4.25,10.0
                if p > 0  : b,c=2.9,0.0     # b= 2.9 dans l'article, 2 dans le rapport
                # je cherche si le dénominateur peut s'annuler :
                # alpha,beta,gam=a*c,a+c,1+a*b*p
                # delta=beta**2-4*alpha*gam
                # if delta >= 0:
                #     print("p+ = %f D= 0 pour v0 = %f et   %f MAIS v0 =  %f"%(p,(-beta+np.sqrt(delta))/(2*alpha),(-beta-np.sqrt(delta))/(2*alpha),v0plus))
                sol=self.Aplus/(a*(v0plus+b*p/(1+c*v0plus))+1.0)
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
               
            return s

    def A_vanDriest_Cebeci1(self,v0plus,Pplus):
        """
        constante de van Driest en fonction du gradient de pression et 
        de la vitesse de transpiration. Relation de Cebeci, 1969
        """
        yplus=11.8
        eps=1.0                # originalement c'est -1, comme on trouve une tendance opposée à Kays, on doit mettre 1,
                               # la définition de Kays est l'opposée de celle de Cebeci (en dp/dx pour Kays, -dp/dx pour Cebeci)
        if isinstance(Pplus,float):
            if v0plus == 0 : 
                return self.Aplus/np.sqrt(1+eps*Pplus*yplus)
            else :
                tmp=np.exp(v0plus*yplus)
                return self.Aplus/np.sqrt(tmp+eps*Pplus/v0plus*(tmp-1))
        else:
            s=[]
            
            for p in Pplus:
                if v0plus == 0 : 
                    sol=self.Aplus/np.sqrt(1+eps*p*yplus)
                else :
                    tmp=np.exp(v0plus*yplus)
                    sol=self.Aplus/np.sqrt(tmp+eps*p/v0plus*(tmp-1))
                if  (sol > 100) or (sol < 0) :
                    s.append(np.nan)
                else:
                    s.append(sol)
               
            return s

    # **********************************************
    #  SOUS-COUCHE et COUCHE VISQUEUSE
    # **********************************************

    def calcul_A_B(self):
        """
        Calcul des constantes de la loi log
        """
        n_inf,n_sup=np.argmin(np.abs(self.yp-100)), np.argmin(np.abs(self.yp-1000))
        print("Calcul des constantes A et B, n_inf = %i, n_sup = %i"%(n_inf,n_sup))
        self.A=(self.up[n_sup]-self.up[n_inf])/np.log(self.yp[n_sup]/self.yp[n_inf])
        self.B=self.up[n_inf]-self.A*np.log(self.yp[n_inf])
        print('A  + %f, B = %f, 1/A = %f'%(self.A,self.B,1/self.A))

    def loi_log(self):
        """
        Loi log asymptotique en présence d'aspiration
        """
        ul=self.B+self.A*np.log(self.yp)
        print('A  + %f, B = %f, 1/A = %f'%(self.A,self.B,1/self.A))
        self.up=ul[ul>0]
        self.yp=self.yp[ul>0]                
        self.dup_dyp=self.A/self.yp

    def SCV(yp,v0):
        """
        Sous couche visqueuse avec approximation
        """
        return (1-np.exp(-v0*yp))/v0

    def calcul_yv(self):
        """
        Plusieurs calculs de l'épaisseur en fonction de différents critères
        v0 est v0+ ici, on retourne y_v^+

        ATTENTION : ce calcul repose sur le fait que A+ est indépendant de v0+ et de p+
                    Mais c'est faux.
        """
        def f1(y):
            """
            Dérivée continue, version exacte à partir de l'équatd$ion différentielle
            """
            return 1-(self.kappa*y)**2*np.exp(-self.v0p*y)

        def f3(y):
            """
            Dérivée continue, version exacte à partir de la dérivée des solutions
            """
            return -2*np.exp(-self.v0p*y)*self.kappa*y+2*np.sqrt(np.exp(-self.v0p*y))
        def s2(v0):
            """
            Dérivée continue, version DL ordre 1 de la fonction f1  
            """
            # au lieu de chercher les zéros de cette fonction, je cherche les racines du polynome
            # return kappa**2*v0*y**3-kappa**2*y**2+1
            s=np.poly1d([self.kappa**2*self.v0p,-self.kappa**2,0,1])
            I=np.argmin((s.r-2)**2)
            return s.r[I]

        def s4(v0):
            """
            Dérivée continue, version DL à partir de la dérivée des solutions (fonction f2) 
            directement la solution (racine d'un polynome de degré 2)
            """
            eta=v0/self.kappa
            if eta <= 6-4*np.sqrt(2):
                return  (2+eta-np.sqrt(eta**2-12*eta+4))/(4*self.v0p)
            else:
                print('valeur limite = %f'%( (6-4*np.sqrt(2))*self.kappa))
                print('erreur sur le calcul approché de y_v')
                return 1e10
            
        print("Problème sur le calcul de y_v")
        y0=0.1
        s1=fsolve(f1,y0)
        s3=fsolve(f3,y0)
        if self.ShowParameters:
            print('Délimitation de la sous couche visqueuse')
            print('prb 1 : exact ',s1[0])
            print('prb 3 : exact ',s3[0])
            print('prb 2 : appro.',s2(self.v0p))
            print('prb 4 : appro.',s4(self.v0p))
        return s1[0]        

    

            # modele =
    #        "SCVonly" : approximation pour la sous couche visqueuse
    #        "CIonly"  : approximation pour la couche interne
    #        "Linear"  : L+ = kappa y+   pour tout y+
    #        "Mixed"   : "Linear" si  y+ <= ypScv sinon "Michel"
    #        "Michel"  : Longueur de Mélange de Michel
    #        "Puissance" : Loi puissance


    def calcul_longueur_melange(self):
        """
        Calcul de la longueur de mélange pour comparer
        """
        Lplus=[]
        for yp in self.yp:
            Lplus.append(self.longueur_melange(yp))
        self.Lplus=Lplus

    def longueur_melange(self,yp):
        """
        choix du modèle de longeur de mélange
        """
        c=0.085
        k=self.Rtau*c
        L=0.0
         
        if self.Modele == "Linear":
           # print('Linear',yp)
            L=self.kappa*yp
        elif self.Modele == "Michel":
            L=k*np.tanh(self.kappa*yp/k)
        elif  yp <= self.ypScv:
             L=self.kappa*yp
        else:
            # print('Mixed autre ',self.Modele)
            L=k*np.tanh(self.kappa*yp/k)
        if self.Damping =="VanDriest":
                L=L*(1-np.exp(-yp/self.Aplus))
        elif self.Damping=="Kays":   
                if yp<=self.Bplus:
                    L=L*yp/self.Bplus

        return L


    def dudy(self,u,y):
        """
        du+/dy+=
        """
        L2=self.longueur_melange(y)**2
        dudy=(-1+np.sqrt(1+4*L2*(1-u*self.v0p)))/(2*L2)
        dudy=2*(1-u*self.v0p)/(1+np.sqrt(1+4*L2*(1-u*self.v0p)))
        return dudy


    def profil_turbulent_aspiration(self):
        """
        Solution avec aspiration
        """
        self.up=odeint(self.dudy,self.ypMin*(1-self.v0*self.ypMin),self.yp)
        #n_inf,n_sup=np.argmin(np.abs(yp-100)), np.argmin(np.abs(yp-1000))
        #A,B=calcul_AB(yp,up,[n_inf,n_sup])
        #print('A  + %f, B = %f, 1/A = %f'%(A,B,1/A))
        #u_log,y_log=loi_log(yp,A,B)
        self.dup_dyp=self.dudy(self.up[:,0],self.yp) # attention aux dimensions en sortie de odeint


    def calcul_parametres(self):
        """
        Test de calcul de différents paramètres de la CL turbulente
        """
        eps=1e-7
        if self.Cf != 0 :
            self.Uep = np.sqrt(2/self.Cf)
        else:
            self.Cf = 2/self.Uep
        if self.Rv0==0:
            self.Rv0=self.v0p*self.Rtau
        else:
            self.v0p = self.Rv0/self.Rtau

        self.Rdelta  = self.Uep*self.Rtau
        self.deltaCL = self.nu*self.Rdelta/self.Ue
        self.Utau    = self.Ue/self.Uep
        self.v0      = self.Rv0*self.nu/self.deltaCL
        if np.abs(self.v0) > eps:
            self.deltav  = self.nu/self.v0
        else:
            self.deltav=1e10

        if self.ShowParameters:
            print('ENTREES')
            print('Cf                  = ',self.Cf)
            print('R_tau               = ',self.Rtau)
            print('U_e                 = ',self.Ue)
            print('nu                  = ',self.nu)
            print('Re_v_0              = ',self.Rv0)
            print('A+                  = ',self.Aplus)
            print('B+                  = %4.3f'%(self.Bplus_function_kappa(self.kappa)))
            print('A+ (calcul)         = %4.3f'%(self.Aplus_function_kappa(self.kappa)))
            print('')
            print('SORTIES')
            print('')    
            print('Ue+                 = ',self.Uep)
            print('Rdelta              = ',self.Rdelta)
            print('delta_CL            = ',self.deltaCL)
            print('U_tau               = ',self.Utau)
            print('Rtau (vérification) = ',self.Utau*self.deltaCL/self.nu)
            print('v0+                 = ',self.v0p)
            print('v0/Ue               = ',self.Rv0/self.Rdelta)
            print('v0                  = ',self.v0)
            if self.deltav != 1e10:
                print('delta=nu/v0         = ',self.deltav)
                print('delta_CL/delta      = ',self.deltaCL/self.deltav)
                print('v_0/u_tau           = ',self.v0/self.Utau)
            else:
                print("pas d'aspiration ou de soufflage")


    # **********************************************
    #  LOIS DE FROTTEMENT SUR LA PLAQUE PLANE
    # **********************************************

    def LudwiegTillman(self,Rdelta2,H):
        """
        Loi qui donne le coefficient de frottement sur une plaque plane
        """
        return 0.246*pow(10,-0.678*H)*pow(Rdelta2,-0.268)

    def Felsch(self,Rdelta2,H):
        """
        Loi qui donne le coefficient de frottement sur une plaque plane
        """
        return 0.058*pow(0.93-1.95*np.log10(H),1.705)*pow(Rdelta2,-0.268)

    def SchultzGrunow(self,Rex):
        """
        Loi qui donne le coefficient de frottement sur une plaque plane
        """
        return 0.370*pow(np.log10(Rex),-2.584)

    def Michel(self,Rex):
        """
        Loi qui donne le coefficient de frottement sur une plaque plane
        """
        return 0.0368*pow(Rex,-1./6.)

    def Cf_Loi_Puissance(self,Rex):
        """
        # calcul du Cf à partir des paramètres de la loi puissance 
        """
        C=pow(self.Cu,-2*self.N/(self.N+3))
        A=(self.N+2)*(self.N+3)/self.N
        N2=-2/(self.N+3)
        return  2*C*pow(A,N2)*pow(Rex,N2)

    # **********************************************
    #   DIVERS PLOTS 
    # **********************************************

    def plot_profile(self, mode="Log",legx=r"$y^+$",legy=r"$u^+$"):
        """
        dessin du profil
        """
        fig = plt.figure()
        fig.suptitle('Profil de la couche limite', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(legx,fontsize=20)
        ax.set_ylabel(legy,fontsize=20)
        if self.Modele=="Laminaire" or mode=="Linear":    
            ax.plot(self.yp,self.up,'-',linewidth=2,color='black')  
        else:
            ax.semilogx(self.yp,self.up,'-',linewidth=2,color='black')
           
        ax.grid()
        plt.show()

    def plot_A_vanDriest(self,v0p,pp,mode="Log",modele='Kays1'):
        """
        dessin de la loi de A+ en fonction de v0+ et de p+
        """
        fig = plt.figure(self.ifig, figsize=(10,8))
        fig.suptitle(r"Modèle  %s : $A^+$ en fonction de $v_0^+$ et de $p^+$"%modele, fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r"$p_0^+$",fontsize=20)
        ax.set_ylabel(r"$A^+/A_0^+$",fontsize=20)
        for v0 in list(v0p):
            if modele=='Kays1':
                s=self.A_vanDriest_Kays1(v0,pp)
            elif modele=='Kays2':
                s=self.A_vanDriest_Kays2(v0,pp)
            elif modele=='Kays3':
                s=self.A_vanDriest_Kays3(v0,pp)
            elif modele=='Cebeci1':
                s=self.A_vanDriest_Cebeci1(v0,pp)
            else :
                print("le modèle choisi n'est pas implémenté")
                s=self.Aplus
            if mode=="Linear":    
                ax.plot(pp,[s/self.Aplus for s in s],'-',linewidth=2,label=r"v_0^+ = %f"%(v0))  
            else:
                ax.semilogx(pp,[s/self.Aplus for s in s],'-',linewidth=2,label=r"v_0^+ = %f"%(v0))  
        ax.legend()  
        ax.grid()
        self.ifig+=1

    def plot_B_Kays(self,v0p,pp,mode="Log"):
        """
        dessin de la loi de B+ en fonction de v0+ et de p+ donnée dans Kays, 1971
         B_Kays3(self,v0plus,Pplus)
        """
        fig = plt.figure(self.ifig, figsize=(10,8))
        fig.suptitle(r"Modèle de Kays : $B^+$ en fonction de $v_0^+$ et de $p^+$", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r"$p_0^+$",fontsize=20)
        ax.set_ylabel(r"$B^+$",fontsize=20)
        for v0 in list(v0p):
            s=self.B_Kays3(v0,pp)
            if mode=="Linear":    
                ax.plot(pp,[s for s in s],'-',linewidth=2,label=r"v_0^+ = %f"%(v0))  
            else:
                ax.semilogx(pp,[s for s in s],'-',linewidth=2,label=r"v_0^+ = %f"%(v0))  
        ax.legend()  
        ax.grid()
        self.ifig+=1

    def plot_Constante_A_B(self):
        """
        graphique des constantes A+ et B+
        """
        fig = plt.figure(self.ifig, figsize=(10,8))
        fig.suptitle(r"constantes $A_0^+$ et $B_0^+$ en fonction de la constante de Karman $\kappa$", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r"$\kappa$",fontsize=20)
        ax.set_ylabel(r"$A^+,B^+$",fontsize=20)
        kappa=np.linspace(0.39,0.45,51)
        ax.plot(kappa,self.Aplus_function_kappa(kappa),'r-',linewidth=2,label=r"A^+")
        ax.plot(kappa,self.Bplus_function_kappa(kappa),'k--',linewidth=2,label=r"B^+")
        ax.legend()  
        ax.grid()   
        self.ifig+=1

    # **********************************************
    #   DIVERS  
    # **********************************************

    def set_title(self):
        """Set a title of a case"""
        print('#',60*'*')
        print('# %s'%(self.name))
        print('#',60*'*','\n')
