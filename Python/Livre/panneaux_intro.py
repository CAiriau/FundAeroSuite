# -*- coding: utf-8 -*-
"""
Exercice : Introduction à la méthode des panneaux :
champ des vitesses généré par un panneau
"""
#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 

 

import numpy    as np
import numpy.ma as ma
import matplotlib.pyplot as plt

plot_panneaux   = False
plot_gamma      = False
plot_vitesse    = False             # pour dessiner la vitesse
plot_psi        = False             # pour dessiner les iso psi


def reference_plaque_plane(x,alpha):
    """
    Solution de référence
     0< x < 1
     alpha en radians
     en sortie : vitesse à la paroi
    """
    return np.cos(alpha)+np.sin(alpha)*np.sqrt((1-x)/x) 


def potentiel_complexe(z,xP,Gamma,alpha,U0=1):
    """
    Potentiel complexe de n tourbillons
    """

    if np.isscalar(z):
        return np.sum(-1J * Gamma/(2*np.pi)*np.log(z-xP))+U0*z*np.exp(-1j*alpha)
    else:
        print("size z in pot. compl. : ",z.shape)
        f=np.zeros(z.shape)+1j*np.zeros(z.shape)
        for i in range(len(Gamma)):
            f+=-1J * Gamma[i]/(2*np.pi)*np.log(z-xP[i])
        return f+U0*z*np.exp(-1j*alpha)

def vitesse_complexe(z,xP,Gamma,alpha,U0=1):
    """
    Calcul de la vitesse complexe
    """
    if np.isscalar(z):
        return np.sum(-1J * Gamma/(2*np.pi)/(z-xP))+U0*np.exp(-1j*alpha)
    else:

        w=np.zeros(z.shape)+1j*np.zeros(z.shape)
        for i in range(len(Gamma)):
            w+=-1J * Gamma[i]/(2*np.pi)/(z-xP[i])
        return w+U0*np.exp(-1j*alpha)

def get_contours(mplcont):
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

def dessiner_fonction_courant(zgrid,Psi,levels,titre,corde=1):
        """
        Tracer des lignes de courant
        opt_levels: 'manuel' ou 'nombre'
        """
        fz=20
        iso_psi=levels.tolist()
        plt.figure(figsize=(10, 10))
        plt.title(r'iso-%s, $\alpha = $ %4.2f °'%(titre,np.rad2deg(alpha)), fontsize=fz, fontweight='bold')
        psi= plt.contour(zgrid.real, zgrid.imag, Psi, levels=iso_psi, colors='blue')
        plt.plot([0,corde],[0,0],'r-',linewidth=5)
        psi_x, psi_y=get_contours(psi)
        plt.plot(psi_x,psi_y)
        plt.legend()

def montrer_gamma(xP,Gamma,corde=1):
    """
    dessiner la distribution du Gamma
    """
    plt.figure(figsize=(15,15))
    plt.title("panneaux")
    plt.plot([0,corde],[0,0],'k--',linewidth=1)
    plt.plot(xP,Gamma,'rs',markersize=5)
    #plt.ylim(-0.01,0.05)


def montrer_panneaux(A,xP,xC,corde=1):
    """
    figure pour montrer le découpage
    """
    plt.figure(figsize=(15,3))
    plt.title("panneaux")
    plt.plot([0,corde],[0,0],'k--',linewidth=1,label='Plaque')
    plt.plot(A,np.zeros(len(A)),'ko',markersize=10,label='limite panneau')
    plt.plot(xP,np.zeros(len(xP)),'bs',markersize=5,label='tourbillons')
    plt.plot(xC,np.zeros(len(xC)),'rs',markersize=5,label='contrôle')
    plt.ylim(-0.01,0.05)
    plt.legend(loc="best")


def plaque_plane(nP,y_wall=0.01,nv=1001,corde=1,U0=1,alpha=0.1):
    """
    Calcul de la plaque plane
    nP  : nombre de panneaux
    """
    A=np.linspace(0,corde,nP+1)         # position des premiers points de chaque panneau
    print("position des panneaux                    :\n",A)

    Lp=A[1]-A[0]                        # Longueur de chaque panneau
    xP=A[0:-1]+1/4*Lp                   # Position des tourbillons
    xC=A[0:-1]+3/4*Lp                   # Position des points de contrôle
    print("xP                                       :\n",xP)

    # Calcul de la répartition de circulation
    Mat=np.zeros((nP,nP))
    for j in range(nP):
        Mat[:,j]=1/(xC-xP[j])           # matrice d'influence
    Rhs=-np.ones(nP)*U0*np.sin(alpha)*2*np.pi
    
    Gamma=np.linalg.solve(Mat,Rhs)      #  résolution du système linéaire
    print("Gamma                                    :\n",Gamma)
    Gamma0 = -np.pi*U0*corde*np.sin(alpha);
    print("Gamma 0                                  :\n",Gamma0)
    if plot_panneaux:
        montrer_panneaux(A,xP,xC)

    if plot_gamma:
        montrer_gamma(xP,Gamma)
    # on vérifie la résolution, on doit trouver -1 au résultat ci-dessous:
    print("sum(Gamma)/(pi sin alpha) =-1            : ",np.sum(Gamma)/(np.pi*np.sin(alpha)))

    wC=vitesse_complexe(xC,xP,Gamma,alpha)
    print("vitesse aux points de collocation : ",wC)

    z_wall=np.linspace(-0.5,1.5,nv)+1j*y_wall*np.ones(nv)
    w_wall=vitesse_complexe(z_wall,xP,Gamma,alpha)
    return [xP,xC,Gamma,z_wall,w_wall,Gamma0]

def dessiner_contours(xP,Gamma,alpha=0.1):
    """
    Pour dessiner les iso-vitesses ou les lignes de courant
    """
    npt=101
    n_iso=101
    Xgrille=np.linspace(-1,2,npt)
    Ygrille=np.linspace(-0.5,0.5,npt)

    X,Y=np.meshgrid(Xgrille,Ygrille)
    Zgrid=X+1j*Y
    Zgrid=ma.masked_where(np.absolute(Zgrid.imag)==0.00, Zgrid)
    f=potentiel_complexe(Zgrid,xP,Gamma,alpha)
    w=vitesse_complexe(Zgrid,xP,Gamma,alpha)

    Psi=f.imag
    Psi_limit=[]
    Psi_limit.append(np.min(Psi))
    Psi_limit.append(np.max(Psi))
    print("Psi, min,max                              :",Psi_limit)
    if plot_psi:
        dessiner_fonction_courant(Zgrid,Psi,np.linspace(Psi_limit[0],Psi_limit[1],n_iso),'r$\Psi$')
    v_limit=[]
    v_limit.append(np.min(-w.imag))
    v_limit.append(np.max(-w.imag))

    if plot_vitesse:
        dessiner_fonction_courant(Zgrid,-w.imag,np.linspace(v_limit[0],Psi_limit[1],n_iso),'v')



# 