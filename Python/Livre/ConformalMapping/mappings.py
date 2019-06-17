#!/bin/python
# -*- coding: utf-8 -*-
"""
  Main program for  conformal mapping for airfoil, Joukowski transform
  mappings file, with the class
  C. Airiau, november, 2018
"""
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.interpolate              import interp1d
from scipy.interpolate      import InterpolatedUnivariateSpline

fz=16   # taille des fontes pour les légendes

#*******************************************
def set_parameters():
#*******************************************
    """
    Set default parameters

    """
    prm={}  # Dictionnary of the parameters  

    # physical parameters

    prm["Uinf"]         = 1.            # vitesse à l'infini amont en m/s
    prm["rho"]          = 1.3           # masse volumique en kg/m^3
    prm["alpha"]        = 0.0           # incidence en degrés
    prm["Umax"]         = 10            # filtrage de la vitesse sur les points singuliers
    prm["Kpmin"]        = -5            # Pour la visualisation des Kp

    # mapping parameter
    prm["airfoil"]      = "Jouko1"      # cas de l'étude, nom du profil
    prm["map"]          = "Joukowski"   # nom de la transformation
    prm["Zc"]           = 0.0+1j*0.0    # position du centre du cercle
    prm["R"]            = 1.            # rayon du cercle
    prm["beta"]         = 0.0           # beta en degrés, angle pour le calcul de la circulation
    prm["Xgrille"]      = np.arange(-3,3, 0.05)  # pour définir la grille dans le plan complexe Z suivant X
    prm["Ygrille"]      = np.arange(-3,3, 0.05)  # pour définir la grille dans le plan complexe Z suivant Y
    prm["lamb"]         = 1.2           # pour le calcul de a = lamb x R
    prm["k"]            = 1.9           # exposant dans la transformation de Karman-Trefftz
    prm["n_circle"]     = 361           # nombre de points sur le cercle dans le plan conforme
    prm["Zs"]           = 1+1j*0        # point singulier du bord de fuite dans le plan conforme
    prm["zs"]           = 1+1j*0        # point singulier du bord de fuite dans le plan physique
    prm["eps"]          = 0.05          # petit paramètre pour les transformations, quand c'est utile

    # pour les iso psi :
    prm["n_iso"]        = 51            # nombres de courbes de niveau visé
    prm["levels_iso"]   = [-2.4,3.4,0.2]    # niveaux [psi_min,psi_max,delta_psi]
    prm["opt_levels"]   = "manuel"      # option pour le calcul des courbes de niveau "manuel" ou "nombre"
    prm["eps_levels"]   = 0.1           # troncature pour dessiner les iso psi [psi_min+eps, psi_max-eps]
    prm["Psi_limit"]    = [-1,1]        # Psi_min,Psi_max
    prm["Psi_method"]   = "Psi"         # "Psi" or "U", méthode pour tracer les isolines de courant
    
    # plot parameters
    prm["plot_Kp"]      = True          # pour tracer le coefficient de pression Kp sur le corps
    prm["calcul_Kp"]    = True         # pour calculer le coefficient de pression Kp sur le corps
    prm["plot_airfoil"] = True          # pour dessiner le profil
    prm["plot_velocity"]= False         # pour dessiner U/Uinf sur le profil
    prm["plot_circle"]  = True          # pour dessiner le cercle dans le plan conforme Z
    prm["plot_psi"]     = True          # pour dessiner les lignes de courant dans le plan physique
    prm["adimChord"]    = True          # vrai : profil sans dimension
    prm["camberline"]   = True          # pour calculer la loi de cambrure et d'épaisseur
    
    prm["xsize"]        = 10            # largeur des figures  
    prm["ysize"]        = 8.5           # hauteur des figures  
    # options diverses


    return prm

# Plan I    : plan du cercle centré  , plan conforme      Z'= R exp(i theta)
# Plan II   : plan du cercle décentré, plan conforme      Z = Zc + Z'
# Plan III  ; plan du profil,          plan physique      z = H(Z') = H(Z-Zc)


class AirfoilConfMap(object):
    """
    airfoil designed by a conformal mapping
    """
    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#',60*'*')
        print('# %s'%("Conformal mappings"))
        print('#',60*'*','\n')
        zeroc=0+1j*0
        # physical parameters
         
        self.Uinf       = prm["Uinf"]   # vitesse à l'infini amont en m/s
        self.rho        = prm["rho"]    # masse volumique en kg/m^3
        self.alpha      = prm["alpha"]  # incidence en degrés
        self.Gamma      = 0.            # circulation à calculée 
        self.w          = 0.            # vitesse complexe sur le profil
        self.Kp         = 0.            # coefficient de pression
        self.Psi        = 0.            # iso Psi, lignes de courant
        self.chord      = 0.            # corde du profil
        self.camber     = 0.            # cambrure relative
        self.thickness  = 0.            # épaisseur relative
        self.A0         = 0.            # premier coefficient de la série de Laurent
        self.A1         = 0.            # second coefficient de la série de Laurent
        self.alpha0     = 0.            # incidence de portance nulle
        self.CL         = 0.            # coefficient de portante
        self.zF         = zeroc         # position du foyer aérodynamique / corde
        self.CmF        = 0.            # coefficient du moment de tangage au foyer

        self.Umax       = prm["Umax"]   # filtrage de la vitesse sur les points singuliers
        self.Kpmin      = prm["Kpmin"]  # Pour la visualisation des Kp
        self.camberline = prm["camberline"] # pour calculer la loi de cambrure et d'épaisseur
        # mapping parameter
        
        self.airfoil    = prm["airfoil"]   # cas de l'étude, nom du profil
        self.map        = prm["map"]       # nom de la transformation
        self.Zc         = prm["Zc"]        # position du centre du cercle
        self.R          = prm["R"]         # rayon du cercle
        self.beta       = prm["beta"]      # beta en degrés, angle pour le calcul de la circulation
        self.eps        = prm["eps"]       # petit paramètre pour les transformations, quand c'est utile
        self.Xgrille    = prm["Xgrille"]   # pour définir la grille dans le plan complexe Z suivant X
        self.Ygrille    = prm["Ygrille"]   # pour définir la grille dans le plan complexe Z suivant Y
        self.lamb       = prm["lamb"]      # pour le calcul de a = R / lamb 
        self.k          = prm["k"]         # exposant dans la transformation de Karman-Trefftz
        self.n_circle   = prm["n_circle"]  # nombre de points sur le cercle dans le plan conforme
        self.Zs         = prm["Zs"]        # point singulier du bord de fuite dans le plan conforme
        self.z_profil   = zeroc            # profil dans le plan physique à calculer
        self.Z_cercle   = zeroc            # cercle décentré dans le plan complexe (plan II)
        self.Z_BA       = zeroc            # position du bord d'attaque sur le cercle
        self.z_BA       = zeroc            # position du bord d'attaque sur le profil
        self.Z_arret    = zeroc            # position du point d'arret sur le cercle
        self.z_arret    = zeroc            # position du point d'arret sur le profil
        self.mu         = [0,0]            # paramètres pour la transformation de von Mises
        self.b          = zeroc            # paramètres pour la transformation double pointe
        self.a          = 0                # paramètre de la transformation de Joukowski
        self.xe         = 0
        self.xi         = 0
        self.ye         = 0
        self.yi         = 0
        self.Kpe        = 0
        self.Kpi        = 0

        # pour les iso psi :
        
        self.n_iso      = prm["n_iso"]  # nombres de courbes de niveau visé
        self.levels_iso = prm["levels_iso"]  # niveaux [psi_min,psi_max,delta_psi]
        self.opt_levels = prm["opt_levels"]  # option pour le calcul des courbes de niveau "manuel" ou "nombre"
        self.eps_levels = prm["eps_levels"]  # troncature pour dessiner les iso psi [psi_min+eps, psi_max-eps]
        self.Psi_limit  = prm["Psi_limit"]   # Psi_min,Psi_max
        self.Psi_method = prm["Psi_method"]  # "Psi" or "U", méthode pour tracer les isolines de courant
        
        # plot parameters
        self.calcul_Kp    = prm["calcul_Kp"]    # pour calculer le coefficient de pression Kp sur le corps
        self.plot_Kp    = prm["plot_Kp"]        # pour tracer le coefficient de pression Kp sur le corps
        self.plot_airfoil = prm["plot_airfoil"] # pour dessiner le profil
        self.plot_velocity= prm["plot_velocity"]# pour dessiner U/Uinf sur le profil
        self.plot_circle  = prm["plot_circle"]  # pour dessiner le cercle dans le plan conforme Z
        self.xsize        = prm["xsize"]        # largeur des figures  
        self.ysize        = prm["ysize"]        # hauteur des figures  
        self.plot_psi     = prm["plot_psi"]     # pour dessiner les lignes de courant dans le plan physique
        self.adimChord    = prm["adimChord"]    # vrai : profil sans dimension
        print("initialization done")

        if self.lamb<=1:                #R/a doit être  supérieur à 1
            raise ValueError('R/a  doit être  supérieur à 1')

        

    def foyer(self,A1,corde):
        """
        Position du foyer en fonction de A1 et alpha_0
        notons : beta = -alpha_0
        """
        z_F=-A1/self.R*np.exp(1j*self.beta)
        print("z_F/corde                                 : %f + i %f"%(z_F.real/corde,z_F.imag/corde))
        print("z_F                                       : %f + i %f"%(z_F.real,z_F.imag))


    def run_airfoil(self):
        """
        main function to run the code
        """

        print("attention : si     a=0  a est calculé en fonction de R (R=1 par défaut)")
        print("            sinon       R est calculé en fonction de a ")

        self.set_title("Type de transformation :" + self.map)
        self.beta        = np.deg2rad(self.beta)
        self.alpha       = np.deg2rad(self.alpha)

        if self.map=="Joukowski":
        
            self.a=self.R/self.lamb
            self.Zs=self.a*np.exp(-1j*self.alpha)
            self.Zc=np.exp(-1j*self.alpha)*(self.a-self.R*np.exp(-1j*self.beta))
            self.Z_BA=self.Zc+self.R*np.exp(1j*(np.pi+self.beta-self.alpha))
            self.Z_arret=self.Zc+self.R*np.exp(1j*(self.beta+np.pi+self.alpha ))


        
        #elif (self.map=="Joukowski1") or (self.map=="Karman-Trefftz")  \
        #      or (self.map=="von Mises") or  self.map=="van de Vooren" or self.map=="double pointe":

        elif self.map in ["Joukowski1","Karman-Trefftz","von Mises","van de Vooren","double pointe"]:
        
            if self.a==0 :
                print(" a is determined")
                self.a=np.sqrt(self.R**2-self.Zc.imag**2)+self.Zc.real
            else : 
                self.R=np.abs(self.a-self.Zc)
                print("==> Le paramètre   'a' est  fixé")

            self.Zs=self.a
            #self.beta=np.arccos((self.a-self.Zc.real)/self.R)
            self.beta=np.arcsin((self.Zc.imag)/self.R)
            self.Z_BA=2*self.Zc.real-self.a               # ce doit être une approximation... faux avec une cambrure élevée
            self.Z_arret=self.Zc+self.R*np.exp(1j*(self.beta+np.pi+2*self.alpha ))

                

        if self.map=="Joukowski1":   
            if self.airfoil=="Ellipse":
                print("cas de l'ellipse")
                # pour l'ellipse il faut fixer uniquement a
                if self.a >= self.R:
                    self.a=0.9*self.R
                    print('on diminue a ',self.a)
                self.Zs=self.R
                self.beta=0
                self.Z_BA=-self.R
                self.beta=0
                self.Z_arret= self.Zc+self.R*np.exp(1j*(self.beta+np.pi+2*self.alpha ))
                print("Corde analytique                          : %f"%(2*(self.R+self.a**2/self.R)))     
            corde=3*self.a-2*self.Zc.real+self.a**2/(self.a-2*self.Zc.real)
            print("Corde analytique (Joukowski)              : %f"%(corde))

        elif self.map=="von Mises":
            self.mu=[-self.a+self.eps*(1+1j),-self.eps*(1+1j)]
            if self.k!=2 : 
                print('changement de la valeur de k à 2')
                self.k=2
            Xc=self.Zc.real
            Corde = (4*(self.a**2+(-2*Xc+1/2*self.mu[0])*self.a+1/2*self.mu[0]**2))*(self.a-Xc)**2/(self.a*(-2*Xc+self.a)**2)
            corde = Corde.real
            print("Corde analytique                          : %f et | %f |"%(corde,abs(Corde)))
            A1 = self.a**2+self.a*self.mu[0]+self.mu[0]**2
            print("A1                                        : %f + i %f"%(A1.real,A1.imag))
            self.foyer(A1,corde)

        elif self.map=="van de Vooren":
            self.L=self.a*2*pow(2/(1+self.eps),self.k-1)
            print("corde L                                   : %f"%(self.L))
            
        elif self.map=="double pointe":
            # pas très clair de changer la position du BA en fonction de alpha non nul, sinon y'a un souci
            # de précision géométrique
            if self.alpha != 0: self.Z_BA=self.Zc+self.R*np.exp(1j*(2*self.beta+np.pi-self.alpha))
            print(self.Z_BA)
        
        elif self.map=="Karman-Trefftz":
            if abs(self.Zc.real) <= 1e-10 :
                raise ValueError('X_c doit être non nul dans Karman-Trefftz')
            A=(1-self.a/self.Zc.real)
            print("A = ",A)
            corde = 2*self.k*self.a*pow(A,self.k)/(pow(A,self.k)-1)
            A1    = (self.k**2-1)*self.a**2/3
            print("Corde analytique                          : %f"%(corde))
            print("A1                                        : %f + i %f"%(A1.real,A1.imag))
            self.foyer(A1,corde)
        

        self.z_arret=self.H(self.Z_arret)
        self.z_BA=self.H(self.Z_BA) 
        self.zs=self.H(self.Zs)
        self.Gamma=-4*np.pi*self.Uinf*self.R*np.sin(self.alpha+self.beta)          #circulation autour du profil
        
        print("a                                         : %f"%(self.a))
        print("R                                         : %f"%(self.R))
        print("beta                                      : %3.10f °"%(np.rad2deg(self.beta)))
        print("Circulation                               : %f"%(self.Gamma))
        print("Plan du cercle                            :")
        print("Centre du cercle Zc                       : %2.10f+ i %2.10f"%(self.Zc.real,self.Zc.imag))
        print("Point singulier Zs                        : %f+ i %f"%(self.Zs.real,self.Zs.imag))
        print("Bord d'attaque Z_BA                       : %f+ i %f"%(self.Z_BA.real,self.Z_BA.imag))
        print("Bord d'attaque Z_arrêt                    : %f+ i %f"%(self.Z_arret.real,self.Z_arret.imag))
        print("Plan du profil                            :")
        print("Bord de fuite  z_s                        : %f+ i %f"%(self.zs.real,self.zs.imag))
        print("Bord d'attaque z_BA                       : %f+ i %f"%(self.z_BA.real,self.z_BA.imag))
        print("Point d'arrêt z_arrêt                     : %f+ i %f"%(self.z_arret.real,self.z_arret.imag))

        self.Z_cercle=self.cercle()            # cercle dans le plan II
        self.z_profil=self.H(self.Z_cercle)    # profil dans le plan III
        
        self.profil_geometrie() 

        CL=8*np.pi/self.chord*np.sin(self.alpha+self.beta)
        print("CL                                        : %f "%(CL))
         
        self.SerieLaurent()
        # Kp ou U  
        if self.calcul_Kp:
            self.dessiner_Kp(self)

        
        #print("Methode pour Psi                          : %s "%(self.Psi_method))
        if self.Psi_method=="Psi":
            # lignes de courant :
            self.LignesCourant()
            if self.plot_psi: self.dessiner_fonction_courant()
        elif self.Psi_method=="U":
            self.iso_vitesse()
        else:
            raise ValueError('Mauvaise option pour tracer les lignes de courant')
       


    def iso_vitesse(self):
        """
        Calcul des vitesses en dehors du profil
        """
        # Définition de la grille dans le plan complexe du cercle
        X,Y=np.meshgrid(self.Xgrille,self.Ygrille)
        Zgrid=X+1j*Y
        Zgrid=ma.masked_where(np.absolute(Zgrid-self.Zc)<=self.R, Zgrid)
        zgrid=self.H(Zgrid)
         
        W=self.vitesse_W(Zgrid-self.Zc)
        w_g  = W /self.H(Zgrid,opt=1)
        u_g,v_g=w_g.real,-w_g.imag
        
        if self.plot_psi:
            planZ=False
            print('stream plot')
                # Nx,Ny=len(self.Xgrille),len(self.Ygrille)
            #y_start,y_end=self.zgrid[0,0].imag,self.zgrid[0,-1].imag
            #x_start,x_end=self.zgrid[0,0].real,self.zgrid[-1,0].real
            size = 15
                #plt.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
            plt.figure(figsize=(size,size*0.8))
            plt.xlabel('x', fontsize=16);plt.ylabel('y', fontsize=16)
            #plt.xlim(x_start,x_end);plt.ylim(y_start,y_end)
            #print(len(self.zgrid.real),len(u_g))
            #print(np.shape(self.zgrid.real))
            #print(np.shape(u_g))
            

            #plt.streamplot(self.zgrid.real, self.zgrid.imag, self.w.real, -self.w.imag, density=1.0, linewidth=1, arrowsize=1, arrowstyle='->')
            if planZ:
                plt.streamplot(Zgrid.real, Zgrid.imag, W.real, -W.imag, density=1.0, linewidth=1, arrowsize=1, arrowstyle='->')
                plt.scatter(Zgrid.real, Zgrid.imag,color='green', s=4, marker='o', linewidth=0); 
                plt.plot(self.Z_cercle.real,self.Z_cercle.imag,label='profil')
            else:
                #for k in range(len(zgrid.real)):
                #    print(len(zgrid.imag[k,:]),len(w_g.real[k,:]),max(zgrid.imag[k,:]),min(zgrid.imag[k,:]))
                #plt.streamplot(zgrid.real, zgrid.imag, w_g.real, -w_g.imag, density=1.0, linewidth=1, arrowsize=1, arrowstyle='->')
                plt.scatter(zgrid.real, zgrid.imag,color='green', s=4, marker='o', linewidth=0); 
                plt.plot(self.z_profil.real,self.z_profil.imag,label='profil')
                plt.plot(self.zs.real,self.zs.imag,'ks',label='B.F.')
                plt.plot(self.z_BA.real,self.z_BA.imag,'ko',label='B.A.')
                plt.plot(self.z_arret.real,self.z_arret.imag,'bs',label='Arrêt')
            #plt.scatter(xs, np.zeros((Ns), dtype=float),color='red', s=30, marker='o', linewidth=0); 
            # pour afficher les points servant au calcul des lignes de courant
               
            #plt.plot(xp,yp,'r-')
            #
            # 
            #plt.legend()
            plt.show()



    def SerieLaurent(self):
        """
        Décomposition en série de Laurent des transformations
        """
        serie=True
        if self.map=="Joukowski1":
            self.A0,self.A1=-self.Zc,self.a**2
        elif self.map=="Karman-Trefftz":
            self.A0,self.A1=-self.Zc,self.a**2*(self.k**2-1)/3
            print("Dièdre au bord de fuite                  : %f °"%((2-self.k)*180))
        elif self.map=="von Mises":
            self.A0,self.A1=-self.Zc,self.a**2-self.mu[0]*self.mu[1]
        elif self.map=="van de Vooren":
            self.A0,self.A1=self.L-self.k*self.a,self.k*(self.k-1)*(1-2*self.eps)*self.a**2/2
            print("Dièdre au bord de fuite                  : %f °"%((2-self.k)*180))
        elif self.map=="double pointe":
            self.A0,self.A1=self.Zc,self.a**2+self.b**2
        else:
            serie=False
            print('pas de série de Laurent disponible pour cette transformation')
     
        if serie:
            self.alpha0=-self.beta  # en radian
            d2=abs(self.A1)
            gamma=np.angle(self.A1,deg=False)/2
            self.CmF=4*np.pi*d2/self.chord**2*np.sin(2*(gamma-self.alpha0))
            self.zF=-self.A1/self.a*np.exp(-1j*self.alpha0)/self.chord
             
            print("Pente du premier axe Delta_1              : %f °"%(np.rad2deg(self.alpha0)))
            print("d^2                                       : %f"%(d2))
            print("gamma pour le second axe Delta_2          : %f °"%(np.rad2deg(gamma)))
            print("CmF                                       : %f "%(self.CmF))
            print("Position du foyer zF /chord               : %f + i %f "%(self.zF.real,self.zF.imag))
        
        return serie

    def H(self,Z,opt=0):
        """
        Transformation conforme de Joukowski
        a : paramètre de la transformation
        Z  : plan complexe du cercle
        opt= 0 : z=H(Z)
        opt= 1 : H'(Z)=dz/dZ
        """
        if self.map=="Joukowski":

            if opt==0:
                return Z+(np.exp(-1j*2*self.alpha)*self.a**2)/Z
            else:
                return 1-(np.exp(-1j*2*self.alpha)*self.a**2)/Z**2
        
        elif self.map=="Joukowski1":
        
            if opt==0:
                return Z+self.a**2/Z
            else:
                return 1-self.a**2/Z**2
        
        elif self.map=="Karman-Trefftz":
        
            G  = (Z-self.a)/(Z+self.a)
            if opt==0:
                return self.k*self.a*(1+pow(G,self.k))/(1-pow(G,self.k))
            else:
                Gp =2*self.a/(Z+self.a)**2
                return 2*self.k**2*self.a*Gp*pow(G,self.k-1)/(1-pow(G,self.k))**2
        
        elif self.map=="von Mises":

            if opt==0:
                return Z+(self.a**2-self.mu[0]*self.mu[1])/Z+0.5*self.a*self.mu[0]*self.mu[1]/Z**2
            else:
                return (1-self.a/Z)*(1-self.mu[0]/Z)*(1-self.mu[1]/Z)

        elif self.map=="van de Vooren":
            
            Z=Z+1j*0
            if opt==0:
                return  self.L+pow(Z-self.a,self.k)/pow(Z-self.eps*self.a,self.k-1)
            else:
                return (Z+(self.a*(self.k-1-self.eps*self.k)))*pow(Z-self.a,self.k-1)/pow(Z-self.eps*self.a,self.k)

        elif self.map=="double pointe":
            
            Z=Z+1j*0

            if opt==0:
                return  Z+(self.a**2+self.b**2)/Z - (self.a*self.b)**2/3/Z**3
            else:
                return  1-(self.a**2+self.b**2)/Z**2 + (self.a*self.b)**2/Z**4      

        else:
            raise ValueError('Mauvaise définition du mapping (self.map)')


    def cercle(self):
        """
        Tracer d'un cercle de centre Z_c, de rayon R dans le plan complexe Z=(X,Y)
        """
        theta=np.linspace(0,2*np.pi,self.n_circle)
        Z=self.Zc+self.R*np.exp(1j*theta)

        alpha_Zs=np.rad2deg(np.arcsin((self.Zc-self.Zs).imag/self.R))
        print("Angle du point singulier sur le cercle    : %f ° "%(alpha_Zs))

        
        if self.plot_circle:
            k=1.2
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Cercle dans le plan complexe', fontsize=fz, fontweight='bold')   
            plt.plot(Z.real,Z.imag,'k-')
            plt.plot([-k*self.R,k*self.R],[0,0],color='blue')
            plt.plot([0,0],[-k*self.R,k*self.R],color='blue')
            plt.plot([-k*self.R,k*self.R],[self.Zc.imag,self.Zc.imag],'r--')
            plt.plot([self.Zc.real,self.Zc.real],[-k*self.R,k*self.R],'r--')
            plt.plot(self.Zc.real,self.Zc.imag,'ro',label='centre')
            plt.plot(self.Zs.real,self.Zs.imag,'ks',label='B.F.')
            plt.plot(self.Z_BA.real,self.Z_BA.imag,'ko',label='B.A.')
            plt.plot(self.Z_arret.real,self.Z_arret.imag,'bs',label='Arrêt')
            plt.axis("equal")
            plt.legend()
            plt.show()
        return Z

    def profil_geometrie(self):
        """
        Caractéristiques géométrique du profil
        """
        print("Caractéristique géométrique du profil : ")
        npt=self.n_circle
        Z_tmp=self.zs-self.z_BA
        incidence = -np.rad2deg(np.arctan(Z_tmp.imag/Z_tmp.real))  
        print("inclinaison du profil en °                : %f  %f "%(incidence,- np.angle(Z_tmp,deg=True)))
        self.chord=(Z_tmp*np.exp(1j*np.deg2rad(incidence))).real
        print("corde du profil                           : %f  %f "%(self.chord,Z_tmp.real/np.cos(np.deg2rad(incidence))))
        if self.adimChord: 
            self.Lref=self.chord
        else:
            self.Lref=1.0


        if self.map=="Joukowski":
            # extrados

            theta_e=np.linspace(-self.alpha-self.beta+0.01,np.pi+self.beta-self.alpha,npt)
            self.Ze=np.flipud(self.Zc+self.R*np.exp(1j*theta_e))
            ze=(self.H(self.Ze)-self.z_BA)*np.exp(1j*np.deg2rad(incidence))
            #print(np.rad2deg(theta_e))
            # intrados
            theta_i=np.linspace(np.pi+self.beta-self.alpha,2*np.pi-self.alpha-self.beta-0.01,npt)
            self.Zi=self.Zc+self.R*np.exp(1j*theta_i)
            zi=(self.H(self.Zi)-self.z_BA)*np.exp(1j*np.deg2rad(incidence))
            #print(np.rad2deg(theta_i))
        
        else:
            # on commence les coordonnées du profil au bord de d'attaque : translation de -z_BA
            theta_e=np.linspace(-self.beta,np.pi+self.beta,npt)
            theta_i=np.linspace(np.pi+self.beta,2*np.pi-self.beta,npt)
            self.Ze=self.Zc+self.R*np.exp(1j*theta_e)
            # il faut rajouter le bord d'attaque et de fuite pour avoir le profil complet 
            ze=np.concatenate(([0+1j*0], np.flipud((self.H(self.Ze,opt=0)-self.z_BA)*np.exp(1j*np.deg2rad(incidence)))))
            ze[-1]=self.chord
            self.Zi=self.Zc+self.R*np.exp(1j*theta_i)
            zi=np.concatenate(( [0+1j*0],(self.H(self.Zi,opt=0)-self.z_BA)*np.exp(1j*np.deg2rad(incidence))))
            zi[-1]=self.chord

            #print("ze",ze[0:5],ze[-5:])
             
        # Recherche du bord de fuite réel du profil :
        # en fait il y a des soucis au bord d'attaque
        Min=max([np.min(ze.real),np.min(zi.real)])
        Max=min([np.max(ze.real),np.max(zi.real)])
        #x=(1-np.cos(np.linspace(0,np.pi-0.1,npt)))/2   # répartition en cosinus
        x=np.linspace(Min,Max,npt)

        # on recherche le bord d'attaque réel :
        Imin=[np.argmin(ze.real),np.argmin(zi.real)]
        print('Imin                                      : ',Imin)
        #print("BA: Ze,Zi                                 : %3.4e + i %3.4e, %3.4e + i %3.4e"%(self.Ze[Imin[0]].real,self.Ze[Imin[0]].imag ,self.Zi[Imin[1]].real,self.Zi[Imin[1]].imag))
        print("BA: ze,zi                                 : %3.4e + i %3.4e, %3.4e + i %3.4e"%(ze[Imin[0]].real,ze[Imin[0]].imag ,zi[Imin[1]].real,zi[Imin[1]].imag)) 
        print("BA : thetae, thetai                       : %3.2f, %3.2f"%(np.rad2deg(theta_e[npt-Imin[0]-1]),np.rad2deg(theta_i[Imin[1]])))
        c=abs(ze[Imin[0]]-ze[-1])
        print("corde corrigée  + erreur relative         : %f %f "%(c,abs(c/self.chord-1)))
        print("incidence corrigée                        : %f °"%(np.angle(-ze[Imin[0]]+ze[-1],deg=True)))


        if self.adimChord: 
            ze,zi = ze/self.chord,zi/self.chord
            x_lim = 1.0
            x     = x/self.chord
            xlabel,ylabel= r'$x/c$',r'$y/c$'
        else:
            x_lim = self.chord
            
            xlabel,ylabel= r'$x$',r'$y$'
        #N=len(ze)
        #print(ze[0:5],ze[N-5:N])

        print("ze: min,max                               : %5.15e  %5.15e "%(np.min(ze.real),np.max(ze.real)))
        print("zi: min,max                               : %5.15e  %5.15e "%(np.min(zi.real),np.max(zi.real)))
        print("x : min,max                               : %5.15e  %5.15e "%(np.min(x),np.max(x)))
      
        if self.camberline:
            interpol=True
            # attention la seconde méthode ne marche pas toujours ...,
            # laisser interpol=True sauf en cas de problème sur les bords.
            if interpol:
             
                fe = interp1d(ze.real, ze.imag, kind='linear')
                fi = interp1d(zi.real, zi.imag, kind='linear')
            else:
                fe = InterpolatedUnivariateSpline(ze.real, ze.imag, k=2)
                fi = InterpolatedUnivariateSpline(zi.real, zi.imag, k=2)
           
            
            # pour éviter les problèmes d'interpolation dues aux arrondis
            ValMax=min([np.max(ze.real),np.max(zi.real)])
            ValMin=max([np.min(ze.real),np.min(zi.real)])
            x[0]=ValMin;x[-1]=ValMax
            ye,yi = fe(x),fi(x)
            ys,yt=(ye+yi)/2,(ye-yi)/2
            #print("yt",yt[0:5],yt[-5:])

            self.thickness,ind_t= np.amax(yt),np.argmax(yt)
            self.camber,ind_c   = np.amax(ys),np.argmax(ys)
            print(" cambrure / c                             : %5.4f à x(%4d)/c = %2.4f "%(self.camber/self.Lref,ind_c,x[ind_c]/self.Lref))
            print(" epaisseur / c                            : %5.4f à x(%4d)/c = %2.4f "%(2*self.thickness/self.Lref,ind_t,x[ind_t]/self.Lref))

        if self.plot_airfoil:
            plt.figure(figsize=(16, 4))
            plt.title(r'Profil sous incidence nulle', fontsize=fz, fontweight='bold')   
            
            plt.plot(ze.real,ze.imag,'k-',linewidth=1,label=r'$y_e$ (upper wall)')
            plt.plot(zi.real,zi.imag,'b-',linewidth=1,label=r'$y_i$ (lower wall)')
            #plt.plot((self.z_profil-self.z_BA).real/self.Lref,(self.z_profil-self.z_BA).imag/self.Lref,'r',linewidth=1,label='profil')
            #plt.plot(x,ye,'rs')
            #plt.plot(x,yi,'rs')
            if self.camberline:
                plt.plot(x,ys,'b--',linewidth=1,label=r'$y_s$ (camber)')
                plt.plot(x,yt,'r--',linewidth=1,label=r'$y_t$ (thickness)')
            plt.xlabel(xlabel,fontsize=fz)
            plt.ylabel(ylabel,fontsize=fz)
           
            plt.axis("equal")
            plt.legend()
            plt.grid()
            plt.xlim(-0.05,x_lim+0.05)
            plt.show()


    
    def w_Kp(self,Z,display=False):
        """
        Calcul de la vitesse et du Kp
        """
        W=self.vitesse_W(Z-self.Zc)
        if display:
            print("Kp \t \t W \t \t z \t \t \t \t theta \t \t H'  : ")
            s="{%3.4e} \t {%3.4e} \t {%3.4e} \t {%3.4e} \t {%4.2f} \t {%3.4e}"
            for Wtmp,Ztmp in zip(W,Z):
                ztmp=(self.H(Ztmp)-self.z_BA)/self.Lref 
                Kp=1-abs(Wtmp/self.H(Ztmp,opt=1))**2
                print(s%(Kp,abs(Wtmp),ztmp.real,ztmp.imag,np.angle(Ztmp-self.Zc,deg=True),abs(self.H(Ztmp,opt=1)))) 
        
        
        Hp=self.H(Z,opt=1)
        #I=np.where(abs(Hp)<1e-10)
        #print('I = ',I)
        #Hp[I]=1e-10
        w=W*1e-10 #pour l'initialiser
        for i in range(len(W)):
            #hp=self.H(Z[i],opt=1)
            hp=Hp[i]
            if abs(hp)>=1e-10:
                w[i]=W[i]/hp
            else:
                print("exception :  Z, H(Z)' = 0                 : ", Z[i],hp)
                w[i]=0


        for i in range(len(w)):
            if abs(w[i])> self.Umax or (np.isnan(w[i])):
                ztmp=self.H(Z[i])
                print("Filtrage de la vitesse complexe en z = %2.5f + i %2.5f pour w = %3.3e + i %3.3e"%(ztmp.real,ztmp.imag,w[i].real,w[i].imag))
                w[i]=float("-inf")
        Kp = 1-np.abs(w/self.Uinf)**2
        print('Kp max                                    : %f '%abs(np.max(Kp)))
        print('Kp min                                    : %f '%abs(np.min(Kp)))
         

        return w,Kp

    def dessiner_Kp(self,x):
        """
        Dessin du Kp fonction de x
        """
        option_Kp=False
        if option_Kp:
            self.w,self.Kp= self.w_Kp(self.Z_cercle)  
        else:
            self.w,self.Kp= self.w_Kp(self.Z_cercle)  
            print("extrados :")
            self.we,self.Kpe=self.w_Kp(self.Ze,display=False)
            print("intrados :")
            self.wi,self.Kpi=self.w_Kp(self.Zi,display=False)   

         
        if self.plot_Kp:
           
            plt.figure(1,figsize=(self.xsize, self.xsize))
            plt.title(r'Coefficient de pression', fontsize=fz, fontweight='bold')   
            if option_Kp:
                plt.plot(self.z_profil.real/self.Lref,self.Kp,color='black')
            else:
                plt.plot((self.z_profil-self.z_BA).real/self.Lref,self.Kp,color='black',label='Kp')
                
                ze_tmp=self.H(self.Ze)-self.z_BA
                zi_tmp=self.H(self.Zi)-self.z_BA
               
                self.xe=ze_tmp.real/self.Lref
                self.xi=zi_tmp.real/self.Lref

                self.ye=ze_tmp.imag/self.Lref
                self.yi=zi_tmp.imag/self.Lref

                 
                plt.plot(self.xe,self.Kpe,color='red',label=r'$Kp^+$')
                plt.plot(self.xi,self.Kpi,color='blue',label=r'$Kp^-$')
                plt.plot(ze_tmp.real/self.Lref,self.H(self.Ze).imag/self.Lref ,color='red',label=r'$y^+$')
                plt.plot(zi_tmp.real/self.Lref,self.H(self.Zi).imag/self.Lref ,color='blue',label=r'$y^-$')
           
            
            #plt.axis("equal")
            plt.grid()
            plt.legend()
            plt.ylim(self.Kpmin,-self.Kpmin)
            plt.show()

        if self.plot_velocity:
            
            plt.figure(2,figsize=(self.xsize, self.xsize))
            plt.title(r'vitesse', fontsize=fz, fontweight='bold')   
            if option_Kp:
                plt.plot(self.z_profil.real/self.Lref,self.w,color='black')
            else:
                ze_tmp=self.H(self.Ze)-self.z_BA
                zi_tmp=self.H(self.Zi)-self.z_BA
                plt.plot(ze_tmp.real/self.Lref,abs(self.we),color='red',label=r'$U^+$')
                plt.plot(zi_tmp.real/self.Lref,abs(self.wi),color='blue',label=r'$U^-$')
                plt.plot(ze_tmp.real/self.Lref,self.H(self.Ze).imag/self.Lref ,color='red',label=r'$y^+$')
                plt.plot(zi_tmp.real/self.Lref,self.H(self.Zi).imag/self.Lref ,color='blue',label=r'$y^-$')
            #plt.axis("equal")
            plt.grid()
            plt.legend()
            plt.show()

            print('w max                                     : %f '%np.max(np.abs(self.w)))
            print('w min                                     : %f '%np.min(np.abs(self.w)))
            print(self.w)
    def F(self,Z):
        """
        potentiel complexe dans le plan I
        """
        U=np.zeros(Z.shape, dtype=np.complex)

        if self.map=="Joukowski":
            with np.errstate(divide='ignore'):
                for m in range(Z.shape[0]):
                    for n in range(Z.shape[1]):             # attention au  numpy bug https://github.com/numpy/numpy/issues/8516
                                                        # on doit le faire élément par élement
                        U[m,n]=self.Gamma*np.log((Z[m,n])/self.R)/(2*np.pi)
            return self.Uinf*(Z + self.R**2/Z) - 1j*U             # potentiel complexe de l'écoulement (f(z)=F(Z) proriété utilisée ici) 
        else:
            with np.errstate(divide='ignore'):
                for m in range(Z.shape[0]):
                    for n in range(Z.shape[1]):             # attention au  numpy bug https://github.com/numpy/numpy/issues/8516
                                                        # on doit le faire élément par élement
                        U[m,n]=self.Gamma*np.log((Z[m,n])/self.R)/(2*np.pi)
            return self.Uinf*(Z *np.exp(-1j*self.alpha) + self.R**2*np.exp(1j*self.alpha)/Z) - 1j*U
              # potentiel complexe de l'écoulement (f(z)=F(Z) proriété utilisée ici) 

    def vitesse_W(self,Z):
        """
        vitesse W dans le plan complexe
        """    
        if self.map=="Joukowski":
            return self.Uinf*( 1 - self.R**2/Z**2) - 1j*self.Gamma/(2*np.pi*Z)
        else: 
            return self.Uinf*( np.exp(-1j*self.alpha) - self.R**2*np.exp(1j*self.alpha)/Z**2) - 1j*self.Gamma/(2*np.pi*Z)


    def LignesCourant(self):
        """
        calcul des lignes de courant, psi(x,y) dans le plan conforme puis physique
        """
        #alpha, beta sont en degrés
         
        # Définition de la grille dans le plan complexe du cercle
        X,Y=np.meshgrid(self.Xgrille,self.Ygrille)
        Zgrid=X+1j*Y
        Zgrid=ma.masked_where(np.absolute(Zgrid-self.Zc)<=self.R, Zgrid)
        self.zgrid=self.H(Zgrid)
        f=self.F(Zgrid-self.Zc)   # f(z) = f(Z)
        self.Psi=f.imag
        Psi_min=np.min(self.Psi)
        Psi_max=np.max(self.Psi)
        #print('Psi_min                                   : %f '%(Psi_min))
        #print('Psi_max                                   : %f '%(Psi_max))
        self.Psi_limit[0:1]=[Psi_min,Psi_max]
        print("Psi, min,max                              :",self.Psi_limit)      

    def get_contours(self,mplcont):
        """
        Récupération les segments de lignes calculés par plt.contour
        """
        conts=mplcont.allsegs                            
        xline=[]
        yline=[]

        for  cont in conts:
            if len(cont)!=0:
                for arr in cont: 
                    
                    xline+=arr[:,0].tolist()
                    yline+=arr[:,1].tolist()
                    xline.append(None) 
                    yline.append(None)

        return xline, yline


    def dessiner_fonction_courant(self):
        """
        Tracer des lignes de courant pour le profil de Joukowski
        alpha en degrés
        opt_levels: 'manuel' ou 'nombre'
        """ 
      
        if self.opt_levels=='manuel':
            iso_psi=np.arange(self.levels_iso[0], self.levels_iso[1], self.levels_iso[2]).tolist()
        elif self.opt_levels=='nombre':
            #DeltaPsi=(Psi_max-Psi_min)/(n_levels-1)
            #iso_psi=np.arange(Psi_min, Psi_max, DeltaPsi).tolist() 
            iso_psi=np.linspace(self.Psi_limit[0]+self.eps_levels,self.Psi_limit[1]-self.eps_levels,self.n_iso).tolist()
        else:
            raise ValueError("opt_levels : 'manuel' ou 'nombre' ")
 
     
        plt.figure(figsize=(self.xsize, self.ysize))
        plt.title(r'iso-$\Psi$, $\alpha = $ %4.2f °'%(np.rad2deg(self.alpha)), fontsize=fz, fontweight='bold')   
        psi= plt.contour(self.zgrid.real, self.zgrid.imag, self.Psi, levels=iso_psi, colors='blue')
        psi_x, psi_y=self.get_contours(psi)
         
        plt.plot(psi_x,psi_y) 
        plt.plot(self.z_profil.real,self.z_profil.imag,label='profil')
        plt.plot(self.zs.real,self.zs.imag,'ks',label='B.F.')
        plt.plot(self.z_BA.real,self.z_BA.imag,'ko',label='B.A.')
        plt.plot(self.z_arret.real,self.z_arret.imag,'bs',label='Arrêt')
        plt.plot([self.z_BA.real,self.zs.real],[self.z_BA.imag,self.zs.imag],'r-',linewidth=1)
        plt.legend()
        plt.show()

    def set_title(self,texte):
        """Set a title """
        print('#',60*'*')
        print('# %s'%(texte))
        print('#',60*'*','\n')
