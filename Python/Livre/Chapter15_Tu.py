#!/bin/python
"""
  Correction des exercices du chapitre 15 - Turbulence
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
import Livre.Turbulence.Turbulence    as Tu
from Livre.Turbulence.misc            import *

Largeur=10
Hauteur=8 
 
def Exo15_Tur_Aspiration():
    """ 
    Etude des profils de turbulence avec ou sans aspiration
    """
    set_title("Etude des profils de la couche limite turbulente")
    ifig=0
    display       = True
    display_perso = True
    displayLm     = False
    option_semilog= False   # seulement pour tracer la longueur de mélange
    etude_aplus   = False   # étude de la loi A+ en fonction de v0⁺
    cas           = 1       # pour le livre d'exercice prendre cas=1


    Profil=[]
    if cas == 0:
        L=[0,1,2,3,4]            # premiers cas
    elif cas == 1:
        L=[6,7,8,9,10,11,12,13]  # figure du dernier exercice
    else:
        L=[5,6,7,8]              # on met ce qu'on veut
   
    Ncas=len(L)
    m=0 
    for k in L:
        set=Tu.set_parameters()
        if k==0:                    # Laminaire
            set["nom"]  = 'CAS LAMINAIRE'
            set["ypMax"]= 6
            set["v0p"]  = 0.01
             
        if k==1:                    # Sous couche visqueuse approchée
            set["nom"]  = 'CAS SOUS COUCHE VISQUEUSE APPROCHEE'
            set["Modele"] = "SCVonly"  #uniquement valable si v0p = 0
            set["ypMax"]  = 5
            set["v0p"]    = 0.0
            
        if k==2:                   # couche interne approchée 
            set["nom"]  = 'CAS COUCHE INTERNE APPROCHEE'
            set["Modele"] = "CIonly" #uniquement valable si v0p = 0
            set["v0p"]    = 0.0
            set['Damping'] = "VanDriest"
             
        if k==3:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.0
           
        if k==4:                    # Modèle de Michel
            set["nom"]  = 'CAS MODELE DE MICHEL SANS CORRECTION'
            set["Modele"] = "Michel" 
            set["v0p"]    = 0
             
        if k==5:                    # Modèle de Michel
            set["nom"]  = 'CAS MODELE DE MICHEL AVEC CORRECTION'
            set["Modele"] = "Michel" 
            set["v0p"]    = 0
            set['Damping'] = "VanDriest"
        
        if k==6:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.0
            #set["Aplus"]  = 0
            #set["Bplus"]  = 0
            set['Damping'] = "VanDriest"
            
           
        if k==7:                    # Loi Log
            set["nom"]  = 'CAS MODELE AVEC LA LOI LOG'
            set["Modele"] = "loi_log" 
            set["A"]=1/set["kappa"]
            set["B"]=5.28

        if k==8:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.03
            set["Aplus"]  = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"

        if k==9:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.03

            set["Aplus"]  = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Cebeci1"
            
        if k==10:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = 0.
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"

        if k==11:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = 0.0
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Cebeci1"

        if k==12:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = -0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"

        if k==13:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = 0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"
                
            

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].calcul_profil()
        Profil[m].calcul_longueur_melange()
        print(Profil[m].Modele,Profil[m].v0p)
        m+=1
         
    if displayLm:
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Longueur de mélange', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if option_semilog: 
                ax.semilogx(Profil[k].yp, Profil[k].Lplus,label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                    Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
            else:
                ax.plot(Profil[k].yp, Profil[k].Lplus,label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                    Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        if option_semilog:
            k=0;ax.semilogx(Profil[k].yp, Profil[k].Lplus,'bo',label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        else:
            k=0;ax.plot(Profil[k].yp, Profil[k].Lplus,'bo',label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        ax.axis([1,200,0,100])
        ax.legend(loc='upper left')    
        ifig+=1
     
    if display:
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Ecoulement turbulent 1', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        #fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if Profil[k].com!="Faux" :
                ax.semilogx(Profil[k].yp,Profil[k].up,label=r'$v_0^+$ = %1.4f, mod.  %s, D: %s, $A^+$ = %4.2f, %s'%(
                    Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus,Profil[k].com))
        ax.axis([1,200,0,25])
        ax.legend(loc='upper left')       
        ifig+=1

    if display_perso and cas==1:    
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Ecoulement turbulent 2', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        k=0;ax.semilogx(Profil[k].yp,Profil[k].up,'s',markevery=(0,5),label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
            Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        k=1;ax.semilogx(Profil[k].yp,Profil[k].up,'o',markevery=(0,7),label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
            Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))  
        k=2;ax.semilogx(Profil[k].yp,Profil[k].up,'-o',label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
            Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        k=3;ax.semilogx(Profil[k].yp,Profil[k].up,'*',label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
            Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))  
        k=4;ax.semilogx(Profil[k].yp,Profil[k].up,'-o',markevery=(0,11),label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
            Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        k=5;ax.semilogx(Profil[k].yp,Profil[k].up,'--',label=r'$v_0^+$ = %1.4f, mod.  %s , D : %s, $A^+$ = %4.2f'%(
            Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        ax.axis([1,200,0,30])
        ax.legend(loc='upper left')   

    if etude_aplus:
        Profil[-1].ifig=ifig+1
        print('influence sur A+ de la transpiration et du gradient de pression')
        Profil[-1].Aplus=26          # A+ donnée par Kays (pas A+ = 26 pour p_0^+=0 et v_0^+=0 )
        v0p=[-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.2,0.3]
        pp =np.linspace(-0.039,0.06,51)
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays1')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays2')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays3')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Cebeci1')
        Profil[-1].plot_B_Kays(v0p,pp,mode="Linear")
        Profil[-1].plot_Constante_A_B() 
        Profil[-1].A_B_regression()

    if display or display_perso or displayLm or etude_aplus:
        plt.show()      


def Exo15_Tur_Loi_Puissance():
    """ 
    Etude des profils de turbulence, loi en puissance
    """
    set_title("Etude des profils de la couche limite turbulente, loi en puissance")

    display       = True
    display_perso = True
    displayLm     = False
    option_semilog= False 
 

    Profil=[]
    L=[0]
    Ncas=len(L)
    m=0 
    for k in L:
        set=Tu.set_parameters()
        if k==0:                     
            set["nom"]  = 'Loi en puissance, Rtau et Cu, et N donnés'
            set["Modele"] = "Puissance" 
            set["ypMax"]= 1
            set["v0p"]  = 0
            set["Rtau"] = 2000
            set["Cu"]   = 8.90
            set["Cf"]   = 0
            set["ShowParameters"] = True

        if k==1:                    
            set["nom"]  = 'Loi en puissance, Rtau et Rex et N  donnés'
            set["Modele"] = "Puissance" 
            set["ypMax"]= 1
            set["v0p"]  = 0
            set["Rtau"] = 2000
            set["Cu"]   = 0
            set["Rex"]  = 2.85e6
            set["N"]    = 9
            set["Cf"]   = 0.00309
            set["ShowParameters"] = True  

        if k==2:                    
            set["nom"]  = 'Loi en puissance, Rtau,Rex et Cu donnés, calculer N'
            set["Modele"] = "Puissance" 
            set["ypMax"]= 1
            set["v0p"]  = 0
            set["Rtau"] = 2000
            set["Cu"]   = 8.75
            set["Rex"]  = 2.85e6
            set["Cf"]   = 0
            set["N"]    = 0.
            set["ShowParameters"] = True  
              
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].calcul_profil()
        Profil[m].plot_profile(mode="Linear")
        m+=1      

    set_question("Dernière question ")
    set=Tu.set_parameters()
    set["nom"]  = 'Loi en puissance, Rtau et Cu inconnus, loi de Michel'
    set["Modele"] = "Puissance" 
    set["ypMax"]= 1
    set["v0p"]  = 0
    set["Rtau"] = 0
    set["Cu"]   = 0
    set["Cf"]   = 0.00309
    set["N"]    = 9
    set["Rex"]  = 2.84959e6
    set["ShowParameters"] = True
    Profil.append(Tu.ProfilTurbulent(set))
    Profil[0].calcul_Cu_from_Cf()
    print('Cu =',Profil[0].Cu)
    Profil[0].calcul_Rtau()
    print('Rtau =',Profil[0].Rtau)   


    set_question("Tracé des lois")
    RexMin,RexMax,nR=1e6,1e7,51
    Rex=np.linspace(RexMin,RexMax,nR)
    print(Profil[0].SchultzGrunow(Rex))
    print(Profil[0].Michel(Rex))
    print(Profil[0].Cf_Loi_Puissance(Rex))
    for Prof in Profil:
        print("Nom = %s , \t N = %i" %(Prof.name,Prof.N))
    Rex_sol=Profil[0].Rex
    print('Rex solution = ',Rex_sol)

    if display_perso:   
        CfMin,CfMax= 0.0024,0.0036
        fig = plt.figure(figsize=(10,8))
        fig.suptitle('Loi de frottement', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$Rex$',fontsize=20)
        ax.set_ylabel(r'$Cf$',fontsize=20) 
        if option_semilog:
            ax.semilogx(Rex,Profil[0].SchultzGrunow(Rex),label='Schultz-Grunow')
            ax.semilogx(Rex,Profil[0].Michel(Rex), label='Michel')
            ax.semilogx(Rex,Profil[0].Cf_Loi_Puissance(Rex), label='présent')
            ax.semilogx([Rex_sol, Rex_sol], [CfMin,CfMax], label='ref')
        else:
            ax.plot(Rex,Profil[0].SchultzGrunow(Rex),label='Schultz-Grunow')
            ax.plot(Rex,Profil[0].Michel(Rex), label='Michel')
            ax.plot(Rex,Profil[0].Cf_Loi_Puissance(Rex), label='présent')
            ax.plot([Rex_sol, Rex_sol], [CfMin,CfMax], label='ref')

        ax.axis([RexMin,RexMax,CfMin,CfMax])
        ax.grid()
        ax.legend(loc='upper right')   
        plt.show()
    


def Exo15_Tur_VanDriest():
    """ 
    Etude des profils de turbulence, effet de la correction de Van Driest
    """
    set_title("Etude des profils de la couche limite turbulente, effet de la correction de Van Driest")

    plot_CorrectionVanDriest = False
    display_profile          = True
    analyse_lois_amortissement = False
    ifig=0
    Profil=[];
    m=0         # initialisation du nombre de cas 
    if plot_CorrectionVanDriest:
        set_question("0 - tracer de la correction de Van Driest")
        set=Tu.set_parameters()
        set["nom"]  = 'Loi avec correction de Van Driest'
        set["Modele"] = "Linear"; set["ypMax"]= 140
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[0].plot_loi_Van_Driest(xaxis="Linear")  # ou "Linear" (échelle linéaire ou semi log)
        m+=1  

    set_question("6- Comparaison des modèles")
    Profil=[]; 
    m=0;
    L=[1,2,3,4]   # choix des différents cas traités
    Ncas=len(L)
   
    for k in L:
        set=Tu.set_parameters()
        set["ypMax"]     = 250
        set["Rtau"]      = 10000 
        set["ShowParameters"] = True
        if k==0:                    # modèle linéaire
            set["nom"]       = 'CAS MODELE LINEAIRE, sans la correction de Van Driest'
            set["Modele"]    = "Linear" 
            set['Damping']   = "No"    
        if k==1:                    # modèle linéaire
            set["nom"]       = 'CAS MODELE LINEAIRE, avec la correction de Van Driest'
            set["Modele"]    = "Linear" 
            set['Damping']   = "VanDriest"
        if k==2:                    # Loi Log
            set["nom"]       = 'CAS MODELE AVEC LA LOI LOG'
            set["Modele"]    = "loi_log" 
            set['Damping']   = "No"
            set["ypMin"]     = 1.
            set["A"]         = 1/set["kappa"]
            set["B"]         = 5.28
        if k==3:                     # u+=y+
            set["nom"]       = 'CAS MODELE DE SOUS COUCHE VISQUEUSE LINEAIRE'
            set["Modele"]    = "SCVlinear"
            set["ypMax"]     = 11.5
            set['Damping']   = "No"
        if k==4:                    # Modèle de Michel
            set["nom"]       = 'CAS MODELE DE MICHEL AVEC CORRECTION'
            set["Modele"]    = "Michel" 
            set['Damping']   = "VanDriest"
         
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].calcul_profil()
        #Profil[m].plot_profile(mode="Log")
        m+=1   


    # dessin automatique des profils   

    if display_profile:
        fig = plt.figure(ifig,figsize=(10,8))
        fig.suptitle('Ecoulement turbulent', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas):
            ax.semilogx(Profil[k].yp,Profil[k].up,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        ax.axis([0.1,300,0,20])
        ax.legend(loc='upper left')       
        plt.show()
        ifig+=1
    
    if analyse_lois_amortissement:
        Profil[-1].ifig=ifig
        print('influence sur A+ de la transpiration et du gradient de pression')
        Profil[-1].Aplus=26          # A+ donnée par Kays (pas A+ = 26 pour p_0^+=0 et v_0^+=0 )
        v0p=[-0.1,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.2,0.3]
        pp =np.linspace(-0.039,0.06,21)
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays1')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays2')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays3')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Cebeci1')
        Profil[-1].plot_B_Kays(v0p,pp,mode="Linear")
        Profil[-1].plot_Constante_A_B()
        plt.show()      
        Profil[-1].A_B_regression()

def Exo15_Tur_tau():
    """ 
    Etude des profils de turbulence, effet de la correction de Van Driest
    tracer de tau laminaire et tau turbulent
    """
    set_title("Etude des profils de la couche limite turbulente, effet de la correction de Van Driest")

    plot_CorrectionVanDriest = False
    display_profile          = True
    analyse_lois_amortissement = False
    ifig=0
    Profil=[];
    m=0         # initialisation du nombre de cas 
    if plot_CorrectionVanDriest:
        set_question("0 - tracer de la correction de Van Driest")
        set=Tu.set_parameters()
        set["nom"]  = 'Loi avec correction de Van Driest'
        set["Modele"] = "Linear"; set["ypMax"]= 140
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[0].plot_loi_Van_Driest(xaxis="Linear")  # ou "Linear" (échelle linéaire ou semi log)
        m+=1  

    set_question("6- Comparaison des modèles")
    Profil=[]; 
    m=0;
    L=[1,2,3,4]   # choix des différents cas traités
    Ncas=len(L)
   
    for k in L:
        set=Tu.set_parameters()
        set["ypMax"]     =  250
        set["Rtau"]      = 10000 
        set["ShowParameters"] = True
        set["n"]         = 501
        if k==0:                    # modèle linéaire
            set["nom"]       = 'CAS MODELE LINEAIRE, sans la correction de Van Driest'
            set["Modele"]    = "Linear" 
            set['Damping']   = "No"    
        if k==1:                    # modèle linéaire
            set["nom"]       = 'CAS MODELE LINEAIRE, avec la correction de Van Driest'
            set["Modele"]    = "Linear" 
            set['Damping']   = "VanDriest"
        if k==2:                    # Loi Log
            set["nom"]       = 'CAS MODELE AVEC LA LOI LOG'
            set["Modele"]    = "loi_log" 
            set['Damping']   = "No"
            set["ypMin"]     = 1.
            set["A"]         = 1/set["kappa"]
            set["B"]         = 5.28
        if k==3:                     # u+=y+
            set["nom"]       = 'CAS MODELE DE SOUS COUCHE VISQUEUSE LINEAIRE'
            set["Modele"]    = "SCVlinear"
            set["ypMax"]     = 11.5
            set['Damping']   = "No"
        if k==4:                    # Modèle de Michel
            set["nom"]       = 'CAS MODELE DE MICHEL AVEC CORRECTION'
            set["Modele"]    = "Michel" 
            set['Damping']   = "VanDriest"
            set["ypMax"]     =  set["Rtau"] 
            
         
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].calcul_profil()
        #Profil[m].plot_profile(mode="Log")
        m+=1   


    # dessin automatique des profils   

    k=Ncas-1
    print(Profil[k].dup_dyp.shape)

    if display_profile:
        fig = plt.figure(ifig,figsize=(10,8))
        fig.suptitle('Ecoulement turbulent', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas):
            ax.semilogx(Profil[k].yp,Profil[k].up,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        #ax.axis([0.1,300,0,20])
        #ax.axis([0.1,1000,0,100])
        ax.legend(loc='upper left')       
        plt.show()
        ifig+=1
    
        fig = plt.figure(ifig,figsize=(10,8))
        fig.suptitle('Ecoulement turbulent', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y/\delta$',fontsize=20)
        ax.set_ylabel(r'$u/Ue$',fontsize=20)
        ax.grid()
        Rtau=set["Rtau"]
        Uep=Profil[-1].Uep

        for k in range(Ncas):
            ax.plot(Profil[k].yp/Rtau,Profil[k].up/Uep,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        #ax.axis([0.1,300,0,20])
        #ax.axis([0.1,1000,0,100])
        ax.legend(loc='upper left')       
        plt.show()
        ifig+=1


        fig = plt.figure(ifig,figsize=(10,8))
        fig.suptitle('Ecoulement turbulent', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_ylabel(r'$y^+$',fontsize=20)
        ax.set_xlabel(r'$\dfrac{du^+}{dy^+}$',fontsize=20)
        ax.grid()
        Rtau=set["Rtau"]
        Uep=Profil[-1].Uep

        k=Ncas-1
        ax.semilogy(Profil[k].dup_dyp,Profil[k].yp,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        #ax.plot(Profil[k].dup_dyp,Profil[k].yp,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        #k=0
        #ax.plot(Profil[k].yp,Profil[k].dup_dyp,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        #ax.axis([0.1,300,0,1])
        ax.axis([0,1,0.1,1000])
        ax.legend(loc='upper left')       
        plt.show()
        ifig+=1
    
def Exo15_Tur_Aspiration_Test():
    """ 
    Etude des profils de turbulence avec ou sans aspiration
    """
    set_title("Etude des profils de la couche limite turbulente")
    ifig=0
    display       = True
    display_perso = False
    displayLm     = False
    option_semilog= False  # seulement pour tracer la longueur de mélange
    etude_aplus   = False   # étude de la loi A+ en fonction de v0⁺
    cas           = 0


    Profil=[]
    if cas == 0:
        L=[5,6,7,8]             
    elif cas == 1:
        L=[6,7,8,9,10,11,12,13]  # figure du dernier exercice
     
   
    Ncas=len(L)
    m=0 
    for k in L:
        set=Tu.set_parameters() 
        set["ShowParameters"] = True    
        if k==5:                    # Modèle de Michel
            set["nom"]  = 'CAS MODELE DE MICHEL AVEC CORRECTION'
            set["Modele"] = "Michel" 
            set["v0p"]    = -0.03
            set['Damping'] = "VanDriest"
            set["ypMax"]     =  set["Rtau"] 
        
        if k==6:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.03
            #set["Aplus"]  = 0
            #set["Bplus"]  = 0
            set['Damping'] = "VanDriest"
            set["ypMax"]     =  set["Rtau"]
            
           
        if k==7:                    # Loi Log
            set["nom"]  = 'CAS MODELE AVEC LA LOI LOG'
            set["Modele"] = "loi_log" 
            set["A"]=1/set["kappa"]
            set["B"]=5.28

        if k==8:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.03
            set["Aplus"]  = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"
            set["ypMax"]     =  set["Rtau"] 

        if k==9:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = -0.03

            set["Aplus"]  = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Cebeci1"
            
        if k==10:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = 0.
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"

        if k==11:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = 0.0
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Cebeci1"

        if k==12:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = -0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"

        if k==13:                    # modèle linéaire
            set["nom"]  = 'CAS MODELE LINEAIRE'
            set["Modele"] = "Linear" 
            set["v0p"]    = 0.03
            set["Aplus"]  = 26
            set["Pplus"]    = 0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel']= "Kays1"
                
            

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].calcul_profil()
        Profil[m].calcul_longueur_melange()
        print(Profil[m].Modele,Profil[m].v0p)
        m+=1
         
    if displayLm:
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Longueur de mélange', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if option_semilog: 
                ax.semilogx(Profil[k].yp, Profil[k].Lplus,label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                    Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
            else:
                ax.plot(Profil[k].yp, Profil[k].Lplus,label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                    Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        if option_semilog:
            k=0;ax.semilogx(Profil[k].yp, Profil[k].Lplus,'bo',label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        else:
            k=0;ax.plot(Profil[k].yp, Profil[k].Lplus,'bo',label=r'mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
        ax.axis([1,200,0,100])
        ax.legend(loc='upper left')    
        ifig+=1
     
    if display:
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle('Ecoulement turbulent 1', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        #fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$',fontsize=20)
        ax.set_ylabel(r'$u^+$',fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if Profil[k].com!="Faux" :
                ax.semilogx(Profil[k].yp,Profil[k].up,label=r'$v_0^+$ = %1.4f, mod.  %s, D: %s, $A^+$ = %4.2f, %s'%(
                    Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus,Profil[k].com))
                if Profil[k].plot_sv:
                    ax.semilogx(Profil[k].y_vp,Profil[k].u_vp,label=r'$u_v^+$')
                    ax.semilogx(Profil[k].yp_sv,Profil[k].up_sv,"ko",markersize=7)
        ax.axis([1,set["Rtau"],0,37])
        ax.legend(loc='upper left')       
        ifig+=1

    if display_perso:    
        if len(Profil)>=6:
            fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
            fig.suptitle('Ecoulement turbulent 2', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80) 
            ax.set_xlabel(r'$y^+$',fontsize=20)
            ax.set_ylabel(r'$u^+$',fontsize=20)
            ax.grid()
            k=0;ax.semilogx(Profil[k].yp,Profil[k].up,'s',markevery=(0,5),label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
            k=1;ax.semilogx(Profil[k].yp,Profil[k].up,'o',markevery=(0,7),label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))  
            k=2;ax.semilogx(Profil[k].yp,Profil[k].up,'-o',label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
            k=3;ax.semilogx(Profil[k].yp,Profil[k].up,'*',label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))  
            k=4;ax.semilogx(Profil[k].yp,Profil[k].up,'-o',markevery=(0,11),label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f'%(
                Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
            k=5;ax.semilogx(Profil[k].yp,Profil[k].up,'--',label=r'$v_0^+$ = %1.4f, mod.  %s , D : %s, $A^+$ = %4.2f'%(
                Profil[k].v0p,Profil[k].Modele,Profil[k].Damping,Profil[k].Aplus))
            ax.axis([1,200,0,30])
            ax.legend(loc='upper left')  
        else:
            print('len (Profil)      :',len(Profil)) 

    if etude_aplus:
        Profil[-1].ifig=ifig+1
        print('influence sur A+ de la transpiration et du gradient de pression')
        Profil[-1].Aplus=26          # A+ donnée par Kays (pas A+ = 26 pour p_0^+=0 et v_0^+=0 )
        v0p=[-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.2,0.3]
        pp =np.linspace(-0.039,0.06,51)
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays1')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays2')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Kays3')
        Profil[-1].plot_A_vanDriest(v0p,pp,mode="Linear",modele='Cebeci1')
        Profil[-1].plot_B_Kays(v0p,pp,mode="Linear")
        Profil[-1].plot_Constante_A_B() 
        Profil[-1].A_B_regression()

    if display or display_perso or displayLm or etude_aplus:
        plt.show()      



