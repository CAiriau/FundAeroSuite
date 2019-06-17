#!/bin/python
"""
  Correction des exercices du chapitre 15 partie turbulence, ancienne version ?
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
from scipy.integrate            import odeint

kappa=0.41
global ModeleTurbulence
ModeleTurbulence=0

def u_aspiration_laminaire(y):
    """
    Loi de u/u0 fonction de y/delta pour l'écoulement laminaire
    """ 
    return 1-np.exp(-y)


def u_aspiration_sous_couche_visqueuse(y,v0):
    """
    Loi de u+ fonction de y+ et v0+ dans la sous couche visqueuse
    """ 
    return (1-np.exp(-v0*y))/v0

def u_aspiration_couche_interne(y,yv,uv,v0):
    """
    Loi de u+ fonction de y+ et v0+,uv+,yv+ dans la couche interne
    """ 
    m=-v0/(2*kappa)
    tmp=m*np.log(y/yv)+np.sqrt(1-v0*uv)
    return (1-tmp**2)/v0
 
def Longueur_Melange(yp,Rtau,yp_scv):
    """
    choix du modèle de longeur de mélange
    """
    c,Aplus=0.085,26.0
    k=Rtau*c
    L=k*np.tanh(kappa*yp/k)
    if ModeleTurbulence==3:
        L=L*(1-np.exp(-yp/Aplus))
    elif ModeleTurbulence==0:
        return kappa*yp
    elif (ModeleTurbulence==1) and (yp <= yp_scv):
        L=kappa*yp
    return L
    
def dudy(u,y,v0,Rtau,yp_scv):
    """
    du+/dy+=
    """
    L2=Longueur_Melange(y,Rtau,yp_scv)**2
    dudy=(-1+np.sqrt(1+4*L2*(1-u*v0)))/(2*L2)
    return dudy

def loi_log(y,A,B):
    """
    Loi log asymptotique en présence d'aspiration
    """
    ul=B+A*np.log(y)
    return ul[ul>0],y[ul>0]

def calcul_AB(y,u,ind):
    """
    calcul des constantes en fonction de la solution
    en deux points définis par les index
    """
    A=(u[ind[1]]-u[ind[0]])/np.log(y[ind[1]]/y[ind[0]])
    B=u[ind[0]]-A*np.log(y[ind[0]])
    return A,B

def calcul_parametres(Cf,Rtau,Re_v0,U0,nu,display=True):
    """
    Test de calcul de différents paramètres de la CL turbulente
    """
    Ue_plus=np.sqrt(2/Cf)
    Rdelta=Ue_plus*Rtau
    delta_CL=nu*Rdelta/U0
    u_tau=U0/Ue_plus
    v0_plus=Re_v0/Rtau
    v_0=Re_v0*nu/delta_CL
    delta=nu/v_0
    if display:
        print('ENTREES')
        print('Cf                  = ',Cf)
        print('R_tau               = ',Rtau)
        print('U_e                 = ',U0)
        print('nu                  = ',nu)
        print('Re_v_0              = ',Re_v0)
        print('')
        print('SORTIES')
        print('')    
        print('Ue+                 = ',Ue_plus)
        print('Rdelta              = ',Rdelta)
        print('delta_CL            = ', delta_CL)
        print('U_tau               = ',u_tau)
        print('Rtau (vérification) = ',u_tau*delta_CL/nu)
        print('v0+                 = ',v0_plus)
        print('v0/U0               = ',Re_v0/Rdelta)
        print('v0                  = ',v_0)
        print('delta=nu/v0         = ',delta)
        print('delta+              = ',1/v0_plus)
        print('delta_CL/delta      = ', delta_CL/delta)
        print('v_0/u_tau           = ',v_0/u_tau)
    return v0_plus,u_tau,v_0,delta_CL,delta

def profil_turbulent_aspiration(Cf,Rtau,Re_v0,yp0,yp_scv,n):
    """
    Solution avec aspiration
    eps : premier point en y+
    """

    print('Cf                  = ',Cf)
    print('R_tau               = ',Rtau)
    print('Re_v_0              = ',Re_v0) 
    Ue_plus=np.sqrt(2/Cf)
    print('Ue+                 = ',Ue_plus)
    v0_plus=Re_v0/Rtau
    print('v0+                 = ',v0_plus)
    print('v0/Ue               = ',v0_plus/Ue_plus)
    print('Modele de longueur de mélange : ',ModeleTurbulence)
    yp=np.exp(np.linspace(np.log(yp0),np.log(Rtau),n))
    up=odeint(dudy,yp0*(1-v0_plus*yp0),yp,args=(v0_plus,Rtau,yp_scv))
    n_inf,n_sup=np.argmin(np.abs(yp-100)), np.argmin(np.abs(yp-1000))
    A,B=calcul_AB(yp,up,[n_inf,n_sup])
    print('A  + %f, B = %f, 1/A = %f'%(A,B,1/A))
    u_log,y_log=loi_log(yp,A,B)

    return up,yp,u_log,y_log,v0_plus 

def calcul_limite_sous_couche_visqueuse(v0):
    """
    Plusieurs calculs de l'épaisseur en fonction de différents critères
    v0 est v0+ ici, on retourne y_v^+
    """
    def f1(y,v0):
        """
        Dérivée continue, version exacte à partir de l'équatd$ion différentielle
        """
        return 1-(kappa*y)**2*np.exp(-v0*y)

    def f3(y,v0):
        """
        Dérivée continue, version exacte à partir de la dérivée des solutions
        """
        return -2*np.exp(-v0*y)*kappa*y+2*np.sqrt(np.exp(-v0*y))
    def s2(v0):
        """
        Dérivée continue, version DL ordre 1 de la fonction f1  
        """
        # au lieu de chercher les zéros de cette fonction, je cherche les racines du polynome
        # return kappa**2*v0*y**3-kappa**2*y**2+1
        s=np.poly1d([kappa**2*v0,-kappa**2,0,1])
        I=np.argmin((s.r-2)**2)
        return s.r[I]

    def s4(v0):
        """
        Dérivée continue, version DL à partir de la dérivée des solutions (fonction f2) 
        directement la solution (racine d'un polynome de degré 2)
        """
        eta=v0/kappa
        if eta <= 6-4*np.sqrt(2):
            return  (2+eta-np.sqrt(eta**2-12*eta+4))/(4*v0)
        else:
            print('valeur limite = %f'%( (6-4*np.sqrt(2))*kappa))
            print('erreur sur le calcul approché de y_v')
            return 1e10
        
   
    y0=0.1
    s1=fsolve(f1,y0,args=(v0,))
    s3=fsolve(f3,y0,args=(v0,))
    print('prb 1 : exact ',s1[0])
    print('prb 3 : exact ',s3[0])
    print('prb 2 : appro.',s2(v0))
    print('prb 4 : appro.',s4(v0))
    return s1[0]    

def Exo15_Tur_1():
    """ 
    profil turbulent avec aspiration (ancienne version sans classe)
    """
    global ModeleTurbulence

    set_title("Loi d'aspiration en laminaire et turbulent")
    plot=True
    nl=101              # nombre de points dans le cas laminaire
    nv=101              # nombre de points dans la sous-couche visqueuse     
    nvi=101             # nombre de points dans la couche interne  
    nt=201 
    plot_laminaire = False
    plot_ssv       = True     # sous-couche visqueuse
    plot_i         = True    # couche interne
    plot_sol_exacte= True     # la figure finale
    plot_ssv_i     = True    # pour  ajouter les approximations  sur la figure finale
    plot_Lm        = True   # pour voir les différents modèles de longueur de mélange appliqués
    Etude_yv       = False   # pour regarder le calcul de y_v
    Un_yv          = False  # choix de l'utilisateur sinon il y a une liste
    Re_tau   = 1000
    nu       = 1e-6           # viscosité cinématique
    rho      = 1e3            # masse volumique
    y_v_plus = 0             # limite de la sous couche visqueuse en +, k% de l'épaisseur de la couche limite
    Cf       = 3.e-3          # coefficient de frottement
    U_0      = 1              # vitesse infinie amont
    Re_lam   = 500             #  v_0 delta/nu  
    yp_0     = 0.001          # premier point dans la CL
    yp_end   = 250            # dernier point de la couche interne
    yp_limite = 40            # yp limite de la sous couche visqueuse pour  le modèle de turbulence uniquement

    
    #
    # ETUDE SUR LA LIMITE DE LA SOUS COUCHE VISQUEUSE
    # 
    if Etude_yv:
        set_question('0 - Etude de y_v')
        if Un_yv:
            ans=float(input('Entrer le nombre Re_lam: \n '))
            Re=[ans]
        else:
            Re=[1e-6,0.01,0.1,1,10,100,150,200,300]
        for Re_lam in Re:
            v_0_plus,u_tau,v_0,delta_CL,delta=calcul_parametres(Cf,Re_tau,Re_lam,U_0,nu,display=False)
            print('v_0_+                        : %5.2e '%(v_0_plus))
            y_v_plus = calcul_limite_sous_couche_visqueuse(v_0_plus)
            print('yv_plus                      : %5.2e '%(y_v_plus))
    else:
        set_question('0 - calcul de quantités')
        v_0_plus,u_tau,v_0,delta_CL,delta=calcul_parametres(Cf,Re_tau,Re_lam,U_0,nu,display=True)
        print('v_0_+                        : %5.2e '%(v_0_plus))
        y_v_plus = calcul_limite_sous_couche_visqueuse(v_0_plus)
        print('yv_plus                      : %5.2e '%(y_v_plus))
        if y_v_plus >= 6:
            print('probleme avec y_v')
            yv_plus=10

 
    set_question('1 - Ecoulement laminaire')

    eta=np.linspace(0,6,nl)    # l'épaisseur est de l'ordre de 6 unités y/delta
    if plot_laminaire:
        fig = plt.figure()
        fig.suptitle('Ecoulement laminaire', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$\eta=y/\delta$',fontsize=20)
        ax.set_ylabel(r'$\frac{u}{U_0}$',fontsize=20)
        ax.grid()
        ax.plot(eta,u_aspiration_laminaire(eta),'-',linewidth=2,color='black')      
        
    
    set_question('2 - Ecoulement turbulent, sous-couche visqueuse')
    print('Re tau                       : %5.2e '%(Re_tau))
    print('U_tau                        : %5.2e '%(u_tau))
    
    y_plus=np.exp(np.linspace(np.log(yp_0),np.log(y_v_plus),nv))
    u_plus=u_aspiration_sous_couche_visqueuse(y_plus,v_0_plus)
    u_v_plus=u_plus[-1]

    print('v_0_+                        : %5.2e '%(v_0_plus))
    print('u_v_+                        : %5.2e '%(u_v_plus))
    print('y_v_+                        : %5.2e '%(y_v_plus))

    if plot_ssv:
        fig = plt.figure()
        fig.suptitle(r'Ecoulement turbulent, sous-couche visqueuse, $v_0^+$ = %1.4f'%(v_0_plus), fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+ v_0^+$',fontsize=20)
        ax.grid()
        ax.plot(y_plus,u_plus*v_0_plus,'-',linewidth=2,color='blue')      
        

    set_question('3 - Ecoulement turbulent, couche interne')

    print('U_tau                        : %5.2e '%(u_tau))
    print('Exposant m                   : %5.2e '%(-kappa*v_0_plus/2))
    y_plus_i=np.exp(np.linspace(np.log(y_v_plus),np.log(yp_end),nvi))
    u_plus_i=u_aspiration_couche_interne(y_plus_i,y_v_plus,u_v_plus,v_0_plus)

    if plot_i:
        fig = plt.figure()
        fig.suptitle(r'Ecoulement turbulent, couche interne, $v_0^+$ = %1.4f'%(v_0_plus), fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+ v_0^+$',fontsize=20)
        ax.grid()
        ax.semilogx(y_plus,u_plus*v_0_plus,'-',linewidth=1,color='blue',label='sous-couche')      
        ax.semilogx(y_plus_i,u_plus_i*v_0_plus,'o',linewidth=2,color='black',label='interne')    
        ax.legend(loc='lower right')  
         

    set_question("4 - Ecoulement turbulent, intégration de l'équation complète")
    print('modele = ',ModeleTurbulence)
    # variable globale "modele"  qui permet de choisir le modèle de turbulence 
    # modele = 0 : L+ = kappa y+   pour tout y+
    #        = 1 :  le même si  y+ <= yp_limite sinon le modèle de Michel
    #        = 2 : modèle de Michel
    #        = 3 : modèle de Michel avec la fonction de Van Driest
    
    # loi log de référence : 
    Yref= np.exp(np.linspace(np.log(yp_0),np.log(Re_tau),nt))
    uref,yref=loi_log(Yref,1/kappa,5.28)
    ModeleTurbulence=3
    Re=np.array([0.001,100,Re_lam,Re_lam,Re_lam])
    leg=[3,3,3,2,0]
    Ncas=len(Re)
    u_t,y_plus_t,L=np.zeros([Ncas,nt]),np.zeros([Ncas,nt]),np.zeros([Ncas,nt])

    v0Plus=np.zeros([Ncas])
    #v0+=1e-6
    for k in range(Ncas):
        ModeleTurbulence=leg[k]
        u_tmp,y_plus_tmp,u_log,y_log,v0Plus[k]=profil_turbulent_aspiration(Cf,Re_tau,Re[k],yp_0,yp_limite,nt)
        for i in range(nt):
            u_t[k,i],y_plus_t[k,i],L[k,i]=u_tmp[i],y_plus_tmp[i],Longueur_Melange(Yref[i],Re_tau,yp_limite)
        #u_t[k,:],y_plus_t[k,:],u_log,y_log,v0Plus[k]=profil_turbulent_aspiration(Cf,Re_tau,Re[k],yp_0,yp_limite,nt)
    
   
    if plot_sol_exacte:
        fig = plt.figure()
        fig.suptitle('Ecoulement turbulent', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas-2):
            ax.semilogx(y_plus_t[k,:],u_t[k,:],'-',linewidth=2,label=r'$v_0^+$ = %1.4f, mod. n° %i'%(v0Plus[k],leg[k]))   
        ax.semilogx(y_plus_t[Ncas-2,:],u_t[Ncas-2,:],'s',linewidth=2,markevery=(0,5),label=r'$v_0^+$ = %1.4f, mod. n° %i'%(v0Plus[Ncas-2],leg[Ncas-2])) 
        ax.semilogx(y_plus_t[Ncas-1,:],u_t[Ncas-1,:],'-o',linewidth=2,markevery=(0,10),label=r'$v_0^+$ = %1.4f, mod. n° %i'%(v0Plus[Ncas-1],leg[Ncas-1]))      
        #ax.semilogx(y_log,u_log,'--',linewidth=2,color='red' ,label='loi log')
        ax.semilogx(yref,uref,'--',linewidth=2,color='green' ,label='loi log Ref')
        if plot_ssv_i:
            ax.semilogx(y_plus,u_plus,'-',linewidth=3,color='blue',label=r'sous-couche, $v_0^+$ %1.4f'%(v_0_plus))      
            ax.semilogx(y_plus_i,u_plus_i,'-',linewidth=3,color='orange',label=r'interne, $v_0^+$ %1.4f'%(v_0_plus))
        ax.axis([0.1,200,0,30])
        ax.legend(loc='upper left')         
        
    #
    # ETUDE SUR LE MODELE DE LONGUEUR DE MELANGE
    #     
    if plot_Lm:
        fig = plt.figure()
        fig.suptitle('Longueur de Mélange', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$L^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas-2):
            #ax.semilogx(Yref[:],L[k,:],'-',linewidth=2,label='Modele v0+ = %1.4f, mod. n° %i'%(v0Plus  
            ax.plot(Yref[:],L[k,:],'-',linewidth=2,label=r'$v_0^+$ = %1.4f, mod. n° %i'%(v0Plus[k],leg[k])) 
        ax.plot(Yref[:],L[Ncas-2,:],'-s',linewidth=2,markevery=(0,10),label=r'$v_0^+$ %1.4f, mod. n° %i'%(v0Plus[Ncas-2],leg[Ncas-2]))  
        ax.plot(Yref[:],L[Ncas-1,:],'-o',linewidth=2,markevery=(0,0.1),label=r'$v_0^+$ = %1.4f, mod. n° %i'%(v0Plus[Ncas-1],leg[Ncas-1]))
        #ax.semilogx(Yref[:],L[k,:],'s',linewidth=2,markevery=(0,5),label='Modele v0+ plot mod. n° %i'%(v0Plus[Ncas-2],leg[Ncas-2])) 
        #ax.semilogx(Yref[:],L[k,:],'-o',linewidth=2,markevery=(0,10),label='Modele v0+ = %1.4f, mod. n° %i'%(v0Plus[Ncas-1],leg[Ncas-1]))      
        ax.axis([0,200,0,100])
        ax.legend(loc='upper left')         
    plt.show() 


    return


