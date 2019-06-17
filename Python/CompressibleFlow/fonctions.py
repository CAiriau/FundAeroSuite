#!/bin/python
'''
 General functions for compressible flows programs
''' 
#    FundAeroSuite
#   @date   : Novembre 2018 
#   @author : Christophe.airiau@imft.fr 
 

import numpy as np
import matplotlib.pyplot as plt
 
R_gaz_parfait=np.float64(8.314462)
M_molaire_Air=28.976e-3
r_Air=R_gaz_parfait/M_molaire_Air
print("r de l'air très précis = ",r_Air) 
# constantes
gamma=np.float64(1.4)       # Cp/Cv
r=np.float64(287)           # r perfect gas


def sound_velocity(self,gamma=gamma):
    """ Sound velocity of a gas"""
    return np.sqrt(gamma*r*self)

def Pressure_coefficient(self,M,gamma=gamma):
    """ Pressure coefficient function of p/p0 and M, p0 :static reference pressure """
    return 2/gamma/M**2*(self-1);

def p_pi(Mach,gamma=gamma):
    """
    ratio pressure/isentropic pressure function of Mach number
    """
    return (T_Ti(Mach,gamma=gamma))**(gamma/(gamma-1));

def rho_rhoi(Mach,gamma=gamma):
    """
    ratio rho/isentropic rho function of Mach number
    """
    return (T_Ti(Mach,gamma=gamma))**(1./(gamma-1));

def T_Ti(Mach,gamma=gamma):
    """
    ratio Temperature/isentropic Temperature function of Mach number
    """
    return 1./(1+(gamma-1)/2*Mach**2) ;

def pt_p(Mach,gamma=gamma):
    """
    ratio Total pressure / static pressure function of Mach number
    """
    return 1+gamma/2*Mach**2 ;

def inverse_p_pi(rapport,gamma=gamma):
    """
    Calcul du nombre de Mach connaissant la valeur de p/pi
    """
    return np.sqrt(2./(gamma-1)*(pow(rapport,(1.-gamma)/gamma)-1.))


def compressibility_corrections(self,M):
    """
    3 compressibility corrections: Prandtl-Glauert (PG), Karman-Tsien (KT),
    Laitone (L) as a function of Mach number
    self = Kp incompressible flow, M : Mach number
    """
    beta=np.sqrt(1-M**2)
    Kp_out=np.zeros((len(M),3))
    Kp_PG=self/beta
    Kp_KT=self/(beta+M**2/2*self/(1+beta))
    Kp_L =self/(beta+M**2*self/(2*beta)*(1+(gamma-1)/2*M**2))
    Kp_out[:,0]=Kp_PG
    Kp_out[:,1]=Kp_KT
    Kp_out[:,2]=Kp_L
    return Kp_out

def Kp_critical(self,gamma=gamma):
    """ Critical Kp function of the Mach number
    """
    return 2/(gamma*self**2)*(p_pi(1,gamma=gamma)/p_pi(self,gamma=gamma)-1)

def CL_compressible(self,CL):
    """
        Prandtl-Glauert Correction on the C_L
    """
    return CL/np.sqrt(1-self**2)

def PG_correction_wing(self,M,Lambda):
    """
    3 compressibility corrections of a swept wing with Prandtl-Glauert (PG),  as a function of Mach number
    self = Kp incompressible flow, M : Mach number, Lambda : swept angle
    """
    Lamb=Lambda/180*np.pi   # angle in radians
    cos=np.cos(Lamb); beta=np.sqrt(1-(M*cos)**2)
    return self/beta*cos

def omega_super(Mach,gamma=gamma):
    """
    fonction de Prandtl-Meyer en radians, pour les écoulements supersoniques
    """
    a=(gamma+1)/(gamma-1);x=Mach**2-1.0
    return np.sqrt(a)*np.arctan(np.sqrt(x/a))-np.arctan(np.sqrt(x))
    

def omega_hyper(Mach,gamma=gamma):
    """
    fonction de Prandtl-Meyer en radians, pour les écoulements hypersoniques
    """
    omega_inf=np.pi/2*(np.sqrt((gamma+1)/(gamma-1))-1.0);
    return  omega_inf-2.0/((gamma-1)*Mach)


def invomega(angle,show=False,gamma=gamma):
    """
    Inversion de la fonction omega en écoulement supersonique
    angle en radians
    Méthode de Newton
    """
    er0,dmach,mach,iterMax=1e-6,1e-3,2.0,20
    i,erreur=1,1.0

    while (erreur > er0) and (i <=iterMax) :
        omega0=omega_super(mach,gamma=gamma)-angle
        omega1=omega_super(mach+dmach,gamma=gamma)-angle
        dm=- dmach * omega0/(omega1-omega0)
        erreur=abs(dm/mach)
        if show:
            print('i =  %2i, \t erreur = %12.6f, \t dm = %12.6f, \t M = %12.6f'%(i,erreur,dm,mach))
        mach+=dm
        i=i+1

    if  i > iterMax :
        print("pas de convergence")
        mach=-1

    return mach

def detente_isentropique(theta,Mach,gamma=gamma):
    """
    rapport des pressions à travers une détente isentropique
    """
    omega_aval=omega_super(Mach,gamma=gamma)+np.abs(theta)
    rpi=[]
    for omega in omega_aval:
        rpi.append(p_pi(invomega(omega),gamma=gamma)/p_pi(Mach,gamma=gamma))
    return rpi       

#********   CHOC OBLIQUE    *********

def Angle_Choc(sigma,teta,M0):
    """     
    fonction donnant l'angle du choc en radian
    """
    return np.tan(sigma-teta)/np.tan(sigma)-2.0/(gamma+1.0)/M0**2/(np.sin(sigma))**2 - (gamma-1.0)/(gamma+1.0)  
       

def valeur_sigma(M0, theta_d,sigma=0,show=False):
    """
    méthode de Newton pour déterminer l'angle de choc oblique en supersonique
    angle de déviation en degrés
    """
    teta=np.deg2rad(theta_d)
    if sigma == 0: sigma=np.arcsin(1.0/M0)
    er0,ds,erreur=1e-6,np.deg2rad(1e-2),1.0
    iterMax,i=20,1
    if show:
        print("Iterations \t sigma(°) \t dsigma  \t erreur")
            
    while (abs(erreur) > er0) and (i<=iterMax):
        f0=Angle_Choc(sigma,teta,M0)
        f1=Angle_Choc(sigma+ds,teta,M0)
        dsigma=- ds * f0/(f1-f0)
        erreur=dsigma/sigma
        if show: 
            print(' %2i \t %12.6f \t %12.6f \t %12.6f'%(i,np.rad2deg(sigma),np.rad2deg(dsigma),erreur))
        sigma+=dsigma
        i=i+1
    if i > iterMax: 
        print('pas de convergence,angle de déviation  est supérieur à celui détachant le choc' )
        sigma=-1
    elif sigma < 0:
        print('Non convergence, angle négatif : Il faut augmenter le Mach amont ou diminuer angle de déviation')
        sigma=-1
    elif np.rad2deg(sigma)>90: 
        print('Non convergence, angle > 90 :Il faut augmenter le Mach amont ou diminuer angle de déviation')
        sigma=-1
    elif (show): print('convergence sur sigma')
                   
    return np.rad2deg(sigma)
       
def P2_P1(mach,gamma=gamma):      
    """ 
    saut de pression à travers un choc
    """      
    return 2.0*gamma/(gamma+1.0)*mach**2-(gamma-1.0)/(gamma+1.0)
          
def Inverse_P2_P1(rapport):
    """     
    Fonction inverse, connaissant le saut de pression à travers un choc
        on retrouve le Mach normal amont
    """      
    return np.sqrt(1.0/(2.0*gamma)* ((gamma+1.0)* rapport+ (gamma-1.0)))
          
def rho2_rho1(mach,gamma=gamma):
    """    
     saut de masse volumique à travers un choc   
    """  
    return 1.0/( 2.0/((gamma+1.0)* mach**2)+ (gamma-1.0)/(gamma+1.0) )

def pi2_pi1(mach,gamma=gamma):
    """
      rapport des pressions isentropiques à travers un choc
    """
    return pow(P2_P1(mach,gamma=gamma),-1/(gamma-1))*pow(rho2_rho1(mach,gamma=gamma),gamma/(gamma-1))
      
def Mach_Aval(mach,gamma=gamma):
    """       
    Mach normal aval à travers un choc
    """      
    return np.sqrt((1.0+ 0.5*(gamma-1.0)* mach**2)/(gamma*mach**2-0.5*(gamma-1.0)))

def Choc_Oblique(Mach,theta,show=True,msg='Choc oblique'):
    """
    Résolution du choc oblique
    Angle en degrés
    """
    print("\n %s "%(msg))
    sigma=np.deg2rad(valeur_sigma(Mach, theta,sigma=0,show=show))
    print("Mach amont                      : %7.4f "%(Mach))
    print("Angle du choc                   : %7.4f °"%(np.rad2deg(sigma)))
    Mn=Mach*np.sin(sigma)
    MnAval= Mach_Aval(Mn)
    print("Mach normal amont               : %7.4f"%(Mn))
    Maval=MnAval/np.sin(sigma-np.deg2rad(theta))
    print("Mach aval                       : %7.4f "%(Maval))
    Paval=P2_P1(Mn)
    print("Pression Pamont/Paval           : %7.4f "%(Paval))

    omegaAmont=omega_super(Mach)
    print("Angle omega Amont               : %7.4f °"%(np.rad2deg(omegaAmont)))
    omegaAval=omega_super(Maval)
    print("Angle omega Aval                : %7.4f °"%(np.rad2deg(omegaAval)))
    # il faut faire le rapport PiAval/PiAmont
    return sigma,Maval,Paval,omegaAmont,omegaAval

def theta_from_Mn0(sigma,Mach,gamma=gamma):
    """
    theta fonction du Mach et de l'angle de choc pour un choc oblique
    """
    Mn0_2=(Mach*np.sin(sigma))**2
    t=(2.+(gamma-1)*Mn0_2)/((gamma+1)*Mn0_2)
    return sigma-np.arctan(t*np.tan(sigma))

def inv_pi2_sur_pi1(r,Mach=2.,show=False,gamma=gamma):
    """
    Calcul du Mach connaissant pi2/pi1 à travers un choc droit
    """
    er0,ds,erreur=1e-6,1e-3,1.0
    iterMax,i=20,1
    if show:
        print("Iterations \t Mach \t dMach  \t erreur")
            
    while (abs(erreur) > er0) and (i<=iterMax):
        f0=pi2_pi1(Mach,gamma=gamma)
        dM=Mach*ds
        f1=pi2_pi1(Mach+dM,gamma=gamma)
        dMach=- dM * (f0-r)/(f1-f0)
        erreur=dMach/Mach
        if show: 
            print(' %2i \t %12.6f \t %12.6f \t %12.6f'%(i,Mach,dMach,erreur))
        Mach+=dMach
        i=i+1
    if i > iterMax: 
        print('pas de convergence' )
        Mach=-1
    elif Mach < 0:
        print('Non convergence, Mach négatif')
        Mach=-1
    elif (show): print('convergence sur le Mach')
    return Mach

#********   CHOC OBLIQUE  (FIN)  *********  
      
#********   ECOULEMENT SUPERSONIQUE 1D   *********

def S_sur_Scrit(M,gamma=gamma):
    """
    S/S_critique en fonction du Mach
    """
    omega=1+(gamma-1.)/2.*M**2
    n=(gamma+1)/(2*(gamma-1))
    return pow(2*omega/(gamma+1),n)/M

def inv_S_sur_Scrit(r,Mach=2.,show=False,gamma=gamma):
    """
    Calcul du Mach connaissant S/Scritique
    """
    er0,ds,erreur=1e-6,1e-3,1.0
    iterMax,i=20,1
    if show:
        print("Iterations \t Mach \t dMach  \t erreur")
            
    while (abs(erreur) > er0) and (i<=iterMax):
        f0=S_sur_Scrit(Mach,gamma=gamma)
        dM=Mach*ds
        f1=S_sur_Scrit(Mach+dM,gamma=gamma)
        dMach=- dM * (f0-r)/(f1-f0)
        erreur=dMach/Mach
        if show: 
            print(' %2i \t %12.6f \t %12.6f \t %12.6f'%(i,Mach,dMach,erreur))
        Mach+=dMach
        i=i+1
    if i > iterMax: 
        print('pas de convergence' )
        Mach=-1
    elif Mach < 0:
        print('Non convergence, Mach négatif')
        Mach=-1
    elif (show): print('convergence sur le Mach')
    return Mach
def qm_1D(S,Mach,pi,Ti,gamma=gamma,r=r_Air):
    """
    débit massique dans une tuyère 1D
    """
    omega=pow(1+(gamma-1)/2*Mach**2,-0.5*(gamma+1)/(gamma-1))
    return S*np.sqrt(gamma/r)*pi/np.sqrt(Ti)*omega*Mach
#********   METHODE DES CARACTERISTIQUES *********

def mu2mach(mu):
    """
    calcul du Mach à partir de mu
    """ 
    #return np.sqrt(1/np.tan(mu)**2+1)
    return 1/np.sin(mu)

def omega_mu(mu,gam=gamma):
    """
    omega fonction de mu, angles en radians
    """
    c=np.sqrt((gam+1)/(gam-1))
    return mu+c*np.arctan(1./(np.tan(mu)*c))-np.pi/2

def mu_omega(omega,gam=gamma):
    """
    Inversion de la relation omega_mu par méthode de Newton
    """
    # Initialisation 
    F = 1
    it = 0
    eps = 1e-6
    # (Point de départ : formule approchee) 
    mu = np.pi/2*(1- pow((omega/2.27),0.25))
    
    while (abs(F)>1e-10) and (it<20) :     #  Boucle du Newton
        F = (omega_mu(mu,gam=gam)-omega)           #  fonction a annuler
        dF = (omega_mu(mu+eps,gam=gam)-omega_mu(mu-eps,gam=gam))/(2*eps)  # derivée de la fonction
        mu = mu-F/dF                         # iteration de Newton
        it += 1
    return mu

def Intersection_droites(M1,M2,p1,p2):
    """
    Intersection de deux droites :
        D1 passant par M1 et de pente p1
        D2 passant par M2 et de pente p2
    """
    x=(M2[1]-M1[1]+p1*M1[0]-p2*M2[0])/(p1-p2)
    y=M1[1]+(x-M1[0])*p1
    return [x,y]

def Intersection_complexe1(thetaMoins,omegaMoins,thetaPlus,omegaPlus):
    """
    Intersection de deux caractéristiques
    C^-  : thetaMoins+omegaMoins = lambdaMoins
    C^+  : omegaPlus -thetaPlus  = lambdaPlus
    """
    lambdaPlus,lambdaMoins=omegaPlus-thetaPlus,omegaMoins+thetaMoins
    ome,the=(lambdaPlus+lambdaMoins)/2,(lambdaMoins-lambdaPlus)/2
    return ome,the

