#!/bin/python
"""
    Correction des exercices du chapitre 10
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 


import numpy as np
import matplotlib.pyplot as plt
from CompressibleFlow.fonctions import *
from CompressibleFlow.tables    import *
from Tools.misc                 import * 
 
Largeur,Hauteur=10,8
#************************************
# main program
#************************************


    

def valX(M,gamma=gamma):
    """
    X en fonction du Mach
    """
    omega=1.+(gamma-1.)/2.*M**2
    return 1./omega

def u_sur_ac(M,gamma=gamma):
    """
    u/ac en fonction du Mach
    """
    omega=1+(gamma-1.)/2.*M**2
    return M*np.sqrt((gamma+1)/(2.*omega))


def S_sur_Sc(x,gamma=gamma):
    """
    S / Sc en fonction de x=p/pi0
    """
    z=pow(x,-1./(gamma-1.))
    tmp=np.sqrt((gamma-1.0)/(gamma+1))*pow(2./((gamma+1.)),1./(gamma-1.))
    return tmp*z/np.sqrt(1-z)

def CF(M,gamma=gamma):
    """
    Coefficient de propulsion
    """
    omega=1+(gamma-1.)/2.*M**2
    n=(gamma+1.)/(2.*(gamma-1.))
    return gamma*pow((gamma+1.)/2.,-n)*M/np.sqrt(omega)

def CF_m(M,gamma=gamma):
    """
    Coefficient de propulsion donné dans le livre de Shapiro après correction de l'erreur
    """
    tmp1=1-pow(p_pi(M),(gamma-1)/gamma)
    n=(gamma+1.)/(gamma-1.)
    tmp2= 2./(gamma-1)*pow(2/(gamma+1.),n)
    return gamma*np.sqrt(tmp1*tmp2)

def CF_add(M,rpa,gamma=gamma):
    """
    Partie supplémentaire en pression du CF
    """
    return (p_pi(M,gamma=gamma)-rpa)*S_sur_Scrit(M,gamma=gamma)

def CFmax(r,gamma=gamma):
    """
    Coefficient de Poussée Maximale  en fonction de r = pa/pi0
    """
    z=pow(r,(gamma-1.)/gamma)
    tmp=pow((gamma+1.)/2.,-(gamma+1)/(2*(gamma-1.)))
    return gamma*tmp*np.sqrt(2./(gamma-1.)) *np.sqrt(1-z)

def plotCF(rpa,borne,n,option,ifig,gamma=gamma):   
    """
    tracer du Coefficient de propulsion
    """
    
    fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
    #fig.suptitle(r'Poussée $F/F_0$ pour $p_a/p_{i_0}=$%e'%(rpa), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_ylabel(r'$F/F_0$',fontsize=20)
    if option=='M':
        ax.set_xlabel(r'$M$',fontsize=20)
        for k in range(len(rpa)):
            r=rpa[k];
            Mtmp=np.linspace(borne[k][0],borne[k][1],n)
            y=CF(Mtmp,gamma=gamma)+CF_add(Mtmp,r,gamma=gamma)
            ax.plot(Mtmp,y,linewidth=2,label=r'$\frac{p_a}{p_{i_0}}=$%6.5f'%(r))
    elif option=='P':
        ax.set_xlabel(r'$p_s/p_{i_0}$',fontsize=20)
        for k in range(len(rpa)):
            r=rpa[k];
            Mtmp=np.linspace(borne[k][0],borne[k][1],n)
            y=CF(Mtmp,gamma=gamma)+CF_add(Mtmp,r,gamma=gamma)
            ax.plot(p_pi(Mtmp),y,linewidth=2,label=r'$\frac{p_a}{p_{i_0}}=$%6.5f'%(r))
    elif option=='Pa':
        ax.set_xlabel(r'$p_s/p_a$',fontsize=20)
        for k in range(len(rpa)):
            r=rpa[k];
            Mtmp=np.linspace(borne[k][0],borne[k][1],n)
            y=CF(Mtmp,gamma=gamma)+CF_add(Mtmp,r,gamma=gamma) 
            ax.plot(p_pi(Mtmp)/r,y,linewidth=2,label=r'$\frac{p_a}{p_{i_0}}=$%6.5f'%(r))
    elif option=='S':
        ax.set_xlabel(r'$S_s/S_c$',fontsize=20)
        for k in range(len(rpa)):
            r=rpa[k];
            Mtmp=np.linspace(borne[k][0],borne[k][1],n)
            y=CF(Mtmp,gamma=gamma)+CF_add(Mtmp,r,gamma=gamma)
            ax.plot(S_sur_Scrit(Mtmp),y,linewidth=2,label=r'$\frac{p_a}{p_{i_0}}=$%6.5f'%(r))
    else :
        ax.set_xlabel(r'$X$',fontsize=20)
        for k in range(len(rpa)):
            r=rpa[k];
            Mtmp=np.linspace(borne[k][0],borne[k][1],n)
            y=CF(Mtmp,gamma=gamma)+CF_add(Mtmp,r,gamma=gamma)
            ax.plot(valX(Mtmp),y,linewidth=2,label=r'$\frac{p_a}{p_{i_0}}=$%6.5f'%(r))
    ax.grid()
    ax.legend(loc='upper right')
    Stmp=S_sur_Scrit(Mtmp)
    Xtmp=valX(Mtmp)
    Ptmp=p_pi(Mtmp)     
    print('range M              = %f \t %f '%(Mtmp[0],Mtmp[-1]))   
    print('range ps/pi0         = %f \t %f '%(Ptmp[0],Ptmp[-1]))  
    for r in rpa:
        print('range ps/pa      = %f \t %f pour rpa= %f  :  '%(Ptmp[0]/r,Ptmp[-1]/r,r))   
        print('CF_max                              = %f '%(CFmax(r)));
    
    print('range S              = %f \t %f '%(Stmp[0],Stmp[-1]))
    print('range X              = %f \t %f '%(Xtmp[0],Xtmp[-1]))       

def Exercice10_1():
    """ 
    Poussée d'un moteur fusée
    """
    ifig=0;
    
    set_title("Poussée d'un moteur fusée")
    plot_ps=False # quelques courbes classiques
    plot_CF=True
    legM='M'        # 'M', 'P' ou 'S' ou 'Pa' ou 'X'
    Sc=10e-4;
    pi0=10e5;
    rpa=[0.005,0.01]
    rapport_section=10.
    n=31
    Ss=Sc*rapport_section

    
    if plot_ps:
        x=np.linspace(0.1,4,n)
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle(r' $p/p_{i_0}=f(M)$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$p/p_{i_0}$',fontsize=20)
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.plot(x,p_pi(x),linewidth=2)
        ax.grid()
        ifig+=1
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle(r' $S/S_{crit}=f(M)$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$S/S_{crit}$',fontsize=20)
        ax.set_xlabel(r'$M$',fontsize=20)
        ax.plot(x,S_sur_Scrit(x),linewidth=2)
        ax.grid()
        x=np.linspace(1.,6,n)
        ifig+=1
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle(r' $p_s/p_a=f(M)$ pour choc droit', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$p_s/p_a$',fontsize=20)
        ax.set_xlabel(r'$M_s$',fontsize=20)
        ax.plot(x,P2_P1(x),linewidth=2)
        ax.grid()

        ifig+=1
        rapport=p_pi(Mach_Aval(x))*pi2_pi1(x)
        rapport1=p_pi(x)*P2_P1(x)
        fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
        fig.suptitle(r' $p_s/p_{i_0}=f(M)$ pour choc droit', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$p_s/p_{i_0}$',fontsize=20)
        ax.set_xlabel(r'$M_s$',fontsize=20)
        ax.plot(x,rapport,linewidth=2)
        ax.plot(x,rapport1,linewidth=2)
        ax.grid()

        # Mach_Aval(mach)   pi2_pi1(self)  P2_P1(mach) 
     
    set_question('4 : Application')
    Ms=inv_S_sur_Scrit(rapport_section,Mach=2.,show=False)
    print('Mach Ms            = %f '%(Ms));
    print('S / Sc  (verif.)   = %f '%(S_sur_Scrit(Ms)));  
    ps_pi0=p_pi(Ms)
    print('ps/pi0             = %f '%(ps_pi0));

   
    rpa.append(ps_pi0)
    print('liste des rapports p_a/pi_0',rpa) 
    

    us_sur_ac=u_sur_ac(Ms)
    print('us/ac              = %f '%(us_sur_ac));
    CF1=gamma*p_pi(1.)*us_sur_ac  # qm us/(Sc pi0)
    print('qm us/(Sc pi0)     = %f '%(CF1));
    F1=CF1*pi0*Sc
    print('F1                 = %f N'%(F1));
    CF1bis=gamma*rapport_section*us_sur_ac**2*rho_rhoi(Ms)*T_Ti(1.)
    print('qm us/(Sc pi0)     = %f '%(CF1bis));
    print('qm us/(Sc pi0)     = %f '%(CF(Ms)));
    print('CF_m (Shapiro)     = %f '%(CF_m(Ms,gamma=gamma)))
    print('CF_m (Shapiro) M=1 = %f '%(CF_m(1.,gamma=gamma)))
    print('CFT (M=1) pa=0     = %f '%(CF_m(1.,gamma=gamma)+p_pi(1)))
    print('den.(Shapirp),pa=0 = %f '%(2*pow(2/(gamma+1),1/(gamma-1))))
        
    
    for r in rpa:
        print()
        print('pa/pi0             = %f'%(r));
        Mopt=inverse_p_pi(r)
        print('Ms optimal         = %f'%(Mopt))
        print('S/Sc optimam       = %f'%(S_sur_Scrit(Mopt)))
        CF2=CF_add(Ms,r) 
        F2=CF2*pi0*Sc
        print('CF2                = %f N'%(CF2));
        print('F2                 = %f N'%(F2));
        print('F1+F2              = %f N'%(F1+F2));
        print('CF1+CF2            = %f '%(CF1+CF2));
        print('ps/pa              = %f '%(ps_pi0/r));
    
    print('CF_max             = %f '%(CFmax(ps_pi0)));
    print('F_max              = %f N'%(CFmax(ps_pi0)*pi0*Sc));
    

    
    #set_question('1 : F sur F0')
     
   
    if plot_CF:
        table_r=[0.004, 0.005,0.0073,0.01]
        print('pi0/pa = ',[1/r for r in list(table_r)])
        
        for leg in ['M','S','P','Pa','X']:
            ifig+=1
            plotCF(table_r,[[3.5,5.5],[3.5,5.5],[3.2,5],[3.,4.5]],n,leg,ifig,gamma=1.4)

    
     
    if plot_CF or plot_ps:
        plt.show()    
     


    set_question('5 : Moteur Vulcain')
    rapportS=58.5 
    D=2.15
    S=np.pi*D**2/4
    Sc=S/rapportS
    gam=1.2
    pi0=115e5
    Ms=inv_S_sur_Scrit(rapportS,Mach=5.,show=False,gamma=gam)
    print('gamma              = %f '%(gam));
    print('Mach Ms            = %f '%(Ms));
    ps_pi0=p_pi(Ms,gamma=gam)
    print('ps/pi0             = %f '%(ps_pi0));
    print('ps                 = %f '%(pi0*ps_pi0));
    print('CF_max             = %f '%(CFmax(ps_pi0,gamma=gam)));
    print('F_max              = %f N'%(CFmax(ps_pi0,gamma=gam)*pi0*Sc));
    
    p0=101325
    altitude=np.linspace(13.7,13.71,11)
    for alt in altitude:
        sigma, delta, theta= Atmosphere(alt)
        print(alt,delta*p0)

def Exercice10_2():
    """ 
    Poussée d'un moteur fusée, formule de Mattingly
    """
    ifig=0;
    verb=False       # verbose
    figure=True
    set_title("Poussée d'un moteur fusée, formules de Mattingly")

    # Tc = température de la chambre de combustion,  Ti0 pour nous
    # Pc = pression dans la chambre de combustion,   Pi0 pour nous

    gam=1.2             # on fixe la valeur de gamma
    m=(gam+1)/(gam-1)   # puissance utile dans les formules 
    R=283               # constante du gaz de combustion   
    Tc=500              # température de la chambre de combustion en K  
    n=21
    Pc_sur_Pe=np.linspace(5,500,n)
    Pe_sur_Pc=1/Pc_sur_Pe     
    Gamma=np.sqrt(gam*pow((gam+1)/2,-m))         # coefficient utile
    print('gam                = %f '%(gam));
    print('Gamma              = %f '%(Gamma));
  
    Cstar=np.sqrt(R*Tc)/Gamma
    print('C^*                = %f m/s'%(Cstar));  # vitesse caractéristique

    #*****************************
    # fig. 3.28 p. 199
    #*****************************
    def vitesse_exit(r,Tc):
        """
        vitesse d 'éjection
        r= Pe/Pc
        """
        tmp,ac=1-pow(r,(gam-1)/gam),gam*R*Tc
        Vmax=np.sqrt(2/(gam-1)*ac)
        Ve=np.sqrt(tmp)*Vmax
        return Ve,Vmax

    Ve,Vmax=vitesse_exit(Pe_sur_Pc,Tc)
    if verb:
        print('Vmax               = %f m/s'%(Vmax));
        print('Ve = ',Ve/Vmax)
    if figure: SimplePlot(0,'$Ve$ réduit',Pc_sur_Pe,Ve/Vmax,[r'$P_c/P_e$',r'$Ve/Vmax$'])
    # pas besoin de dessin
  
    #*****************************
    # fig. 3.30 p. 202
    #*****************************
    def ratio_Ae_sur_At(r):
        """
        Ae/At en fonction de r=Pe/Pc
        """
        tmp=2*gam/(gam-1)*(pow(r,2/gam)-pow(r,1.+1./gam))
        return Gamma/np.sqrt(tmp)

    Ae_sur_At=ratio_Ae_sur_At(Pe_sur_Pc)
    if verb: print('Ae/At = ',Ae_sur_At)  
    if figure: SimplePlot(1,'rapport des sections',Pc_sur_Pe,Ae_sur_At,[r'$P_c/P_e$',r'$\varepsilon=Ae/At$'])
    
    def CFi_opt(r):
        """
        Coefficient de propulsion idéal optimal,  r=Pe/Pc
        """
        tmp=2*gam/(gam-1)*(1-pow(r,gam-1/gam))
        return Gamma*np.sqrt(tmp)
    CFiMax= CFi_opt(Pe_sur_Pc)
    if verb: print('CFi_opt = ',CFiMax)
    if figure: SimplePlot(2,r'$CF_{opt}$',Pc_sur_Pe,CFiMax,[r'$P_c/P_e$',r'$C_{F_\max}$'])
    
    def CFi_vac(r):
        """
        CFi dans le vide, r=Pe/Pc
        """
        return CFi_opt(r)+r*ratio_Ae_sur_At(r)

    CFi_vide=CFi_vac(Pe_sur_Pc)
    if verb: print('CFi vac = ',CFi_vide)  
    if figure: SimplePlot(3,r'$Cf_{vac}$',Pc_sur_Pe,CFi_vide,[r'$P_c/P_e$',r'$C_{F_{vac}}$'])
  
    def CFi_opt2(r,eps):
        """
        Coefficient de propulsion idéal optimal
        r=Pe/Pc, eps=Ae/At
        """
        return 2*gam/(gam-1)*eps*(1-pow(r,1-1/gam))*pow(r,1/gam)
    tabx,taby=[],[]
    epsilon=[5,10,30]
    for eps in epsilon: 
        print('eps = ',eps)
        CFiMax1= CFi_opt2(Pe_sur_Pc,eps)
        if verb: print('CFi_opt = ', CFiMax1)
        taby.append(CFiMax1)
        taby.append(CFiMax)
        tabx.append(Pc_sur_Pe)
        tabx.append(Pc_sur_Pe)
     
    if figure: SimplePlot(4,r'$CF2_{opt}$ pour $\varepsilon$= %4.1f'%(eps),tabx,taby,[r'$P_c/P_e$',r'$C_{F_\max}$'],len(tabx))
    #SimplePlot(4,r'$CF2_{opt}$ pour $\varepsilon$= %4.1f'%(eps),[Pc_sur_Pe,Pc_sur_Pe],[CFiMax,CFiMax1],[r'$P_c/P_e$',r'$C_{F_\max}$'])

    def CFi(r,ra,eps):
        """
        CFi fonction de r=Pe/Pc, eps=Ae/At et de ra=Pe/Pc
        """
        return CFi_opt(r)+(r-ra)*ratio_Ae_sur_At(r)



def Exercice10_3():
    """ 
    Tube de Pitot dans une tuyère amorcée
    """
    
    set_title("Tube de Pitot dans une tuyère amorcée")
    pi0    = 10.         # pression dans la chambre en bar   
    Ti0    = 500.        # température dans la chambre en K
    pi2    = 3.283       # pression totale fournie par le Pitot en bar.
    gam    = 1.4         # Cp/Cv
    S_c    = 0.2         # section au col en m^2
   
    set_question('1 : Mach en sortie')

    # je pars de la solution
    Ms     = 3.0
    print("pi2/pi0              = %f"%(pi2_pi1(3,gamma=gam)))
    # je calcule réellement Ms
    Ms_calcul=inv_pi2_sur_pi1(pi2/pi0,Mach=2.,show=False,gamma=gam)
    print("Ms calculé (vérif)   = %f"%(Ms_calcul))

    set_question('2 : M2 et p2')
    M2= Mach_Aval(Ms,gamma=gam)
    p2=pi2* p_pi(M2,gamma=gam)
    print("M2                   = %f"%(M2))
    print("p2                   = %f bar"%(p2))

    set_question('3 : Section de sortie')

    #rapport S_s/S_c 
    r=S_sur_Scrit(Ms,gamma=gam)
    print("p2                   = %f bar"%(p2))
    print("S_s/S_c              = %f "%(r))
    S_s=S_c*r
    print("S_s                  = %f m^2 "%(S_s))

    set_question('4 : ps,Ts, us')
    ps=pi0*p_pi(Ms)
    print("ps                   = %f bar"%(ps))
    Ts=Ti0*T_Ti(Ms)
    print("Ts                   = %f K"%(Ts))
    us=Ms*np.sqrt(gam*r_Air*Ts)
    print("us                   = %f m/s"%(us))

def Exercice10_4():
    """
    Moteur oxygène-hydrogène
    """
    set_title("Moteur oxygène-hydrogène")
    pi0    = 15e5         # pression totale dans la chambre  de combustion en Pa
    Ti0    = 4000.       # température dans la chambre en K
    ps     = 1174        # pression en sortie
    gam    = 1.22        # Cp/Cv
    r      = 519.6       # r= Cp-Cv
    Ds     = 3.0         # diamètre de la section de sortie.

    Ms=inverse_p_pi(ps/pi0,gamma=gam)
    print("Ms                   = %f"%(Ms))
    Ts=Ti0*T_Ti(Ms,gamma=gam)
    print("Ts                   = %f K"%(Ts))
    rapport_section=S_sur_Scrit(Ms,gamma=gam)
    print("S_s/S_c              = %f"%(rapport_section))
    S_s=np.pi*Ds**2/4
    print("S_s                  = %f m^2"%(S_s))
    qm=qm_1D(S_s,Ms,pi0,Ti0,gamma=gam,r=r)
    print("qm                   = %f kg/s"%(qm))
    a_s=np.sqrt(gam*r*Ts)
    us=Ms*a_s
    Fm=qm*us
    S_c=S_s/rapport_section
    C_F=Fm/(pi0*S_c)
    print("as                   = %f m/s"%(a_s))
    print("us                   = %f m/s"%(us))
    print("S_c                  = %f m^2"%(S_c))
    print("F_m                  = %f N"%(Fm))
    print("C_F                  = %f "%(C_F))