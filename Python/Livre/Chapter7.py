#!/bin/python
"""
  Correction des exercices du chapitre 7
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 



import os
import numpy as np
import matplotlib.pyplot as plt

from CompressibleFlow.fonctions     import *
from IncompressibleFlow.PanelMethod import *
from IncompressibleFlow.Airfoil     import *
from scipy.interpolate              import interp1d
from Tools.misc                     import *

def Exercice7_0():
    """
    calcul des paramètres du profil de Van-Hooren
    suivant les options on peut :
    - calculer t si on connait k ou l'angle du bord de fuite pour une épaisseur donnée
    - calculet k si on connait t pour une épaisseur donnée
    - tracer dans le plan (t,k) la valeur de l'épaisseur
    """
    height,width=8,10
    set_title("Profil de Van-Hooren : définition des paramètres")
    chord=np.float64(1.)                # chord length
    eps=0.1                             # to avoid 0 division at theta=0 or 180 °
    npt=360/2+1                             # number of points on the airfoil
    t=np.float64(0.055)                 # thickness parameter (epsilon in the exercise)
    k=1.900                             # trailing edge angle parameter
    delta_TE=-10.0                      # half Trailing edge  angle (negative value)
                                        # put 0 if k is set, otherwhise  is calculated
    thickness_ref= 0.15                 # la valeur de l'épaisseur visée
    
    if delta_TE != 0:
        k=delta_TE/90.0+2
    
    n=int((npt-1)/2)
    plot_profil=False

    option=0                            # = 0 : calcul de k en fonction t, Newton
                                        # = 1 : calcul de t en fonction k, Newton
                                        # = 2 : carte t,k
    

    def Thickness(self,option=0):
        """
        Calcul de l'épaisseur du profil
        """
        if option==0:
            x,y,theta=Van_Hooren_gmtry(chord,t,self,npt,plot_profil,origin=False,dtheta=0.0)
        else:
            x,y,theta=Van_Hooren_gmtry(chord,self,k,npt,plot_profil,origin=False,dtheta=0.0)
        Imax=np.argmax(y)
        return 2*y[Imax],x[Imax]

    def ThicknessTable(t,k):
        """
        Calcul de l'épaisseur du profil
        """
        x,y,theta=Van_Hooren_gmtry(chord,t,k,npt,plot_profil,origin=False,dtheta=0.0)
        Imax=np.argmax(y)
        # print(Imax, np.rad2deg(theta))
        return 2*y[Imax],x[Imax],np.rad2deg(theta[Imax]),
        

    if option<2:
        if option==0:
            par=k
        else:
            par=t
        errorMax=np.float64(0.0001)
        iterMax=35
        dpar=np.float64(0.0001)*par
        error=1.0
        i=0
        thickness,xpos=Thickness(par,option)
        print(' iter \t parameter \t error \t \t thickness \t  delta parameter')
        while (error > errorMax ) and (i<iterMax):
            thickness0,xpos=Thickness(par+dpar,option)
            deltaPar=(thickness_ref-thickness)/(thickness0-thickness)*dpar
            error=np.abs(deltaPar)
            par=par+deltaPar
            thickness,xpos=Thickness(par,option)
            i=i+1
            print( '%d \t %f \t %e \t %f \t %f  '%(i,par,error,thickness,deltaPar))

        if (error>errorMax):
            print('Pas de convergence')
        else:
            if option==0:
                k=par
            else:
                t=par
            x,y,theta=Van_Hooren_gmtry(chord,t,k,npt,plot_profil,origin=False,dtheta=0.0)
            print('Solution pour k = %f, t = %f, delta_TE = %f  , thickness= %f,  x_pos  = %1.3f '%(k,t,(k-2)*90,thickness,xpos))
            angle=np.rad2deg( np.arctan2(y[n]-y[n-1],x[n]-x[n-1]) ) 
            print('theta = %f, Bord de fuite = %f = %f'%(np.rad2deg(theta[n]),angle,(k-2)*90))

    elif option==2:
        N=41
        t_min,t_max=0,0.1
        k_min,k_max=1.7,2
        t,k=definition_grille(N,t_min,t_max,k_min,k_max)
        thickness,xpos = np.zeros((N,N), dtype=float),np.zeros((N,N), dtype=float)
        thetaMax = np.zeros((N,N), dtype=float)
        for i in range(N):
            for j in range(N):
                thickness[i,j],xpos[i,j],thetaMax[i,j]= ThicknessTable(t[i,j],k[i,j])

        plot_option=1        
        
        plt.figure(figsize=(width, height))
        if plot_option==0:
            plt.imshow(thickness, vmin=thickness.min(), vmax=thickness.max(), origin='lower',
                    extent=[t.min(), t.max(), k.min(), k.max()])
            #plt.scatter(t, k, c=thickness)
            plt.colorbar()
        elif plot_option==1: 
            levels=[0.05,0.10,0.15,0.20,0.25,0.30,0.35]   
            cp = plt.contour(t,k,thickness,levels,colors='black', linestyles='-')
            #plt.clabel(cp, inline=True,fontsize=10)
            plt.title(r'iso-thickness')
            plt.xlabel(r't, $\varepsilon$')
            plt.ylabel('k')
            plt.grid()  

            

          
        #plt.figure(figsize=(10, 10))
        if plot_option==0:
            plt.figure(figsize=(width, height))
            plt.imshow(xpos, vmin=xpos.min(), vmax=xpos.max(), origin='lower',
                    extent=[t.min(), t.max(), k.min(), k.max()])
            #plt.scatter(t, k, c=thickness)
            plt.colorbar()
        elif plot_option==1: 
            
            #levels=[-25,-20,-16,-14,-12,-10,-8,-6,-4,-2]   
            # iso contour de la position du maximun suivant x
            # cp = plt.contour(t,k,xpos,colors='black', linestyles='--')
            #plt.clabel(cp, inline=True,fontsize=10)
            plt.xlabel(r't,$\varepsilon$')
            plt.ylabel('k')  
            #plt.scatter(t, k,color='black', s=6, marker='o', linewidth=1);   
            plt.grid() 



            plt.figure(3,figsize=(width, height))
            # pour récupérer le contour correspondant à l'épaisseur de référence
            cs = plt.contour(t,k,thickness,[thickness_ref])
            p = cs.collections[0].get_paths()[0]
            v = p.vertices

            plt.plot([w[0] for w in v],[w[1] for w in v],
                    linestyle='-', linewidth=1, color='black')
            plt.xlabel(r't,$\varepsilon$')
            plt.ylabel('k') 
            titre='thickness = %f'%(thickness_ref)
            plt.grid()
            plt.axis([0,0.1,1.8,2])
            plt.title(titre)
            x,y = v[:,0],v[:,1]
            f1 = interp1d(x, y, kind='cubic')
            t1 =np.linspace(t_min,t_max,N)
            f2 = interp1d(y,x, kind='cubic')
            print('Interpolation pour une épaisseur de %f'%(thickness_ref))
            print(' t  \t \t k')
            for i in range(N):
                print('%f \t %f '%(t1[i],f1(t1[i])))

            k_min=np.amin(f1(t1))
            n=200
            k_min=int(n*k_min+1)/n
            k_max=int(n*np.amax(f1(t1))-1)/n
            print(k_min,k_max)
            print(' t  \t \t k')
            N1=int(np.ceil((k_max-k_min)*n+1./n))
            print('N1 =',N1 )
            k2 =np.linspace(k_min,k_max,N1)
            for i in range(N1):
                print('%f \t %f '%(f2(k2[i]),k2[i]))


        plt.show()

                 


def Exercice7_1():
    """ 
    Profil de Van-Hooren
    """
    set_title("Profil de Van-Hooren : solution analytique")

    # data, parameters
    U0=np.float64(1)                    # reference velocity
    chord=np.float64(1.)                # chord length
    ell=chord
    aoa=np.array([0,5.0,10.0]);             # angle of attack
    epsilon=np.float64(0.02)             # thickness parameter
    # Paramètre pour avoir un profil proche du NACA 0015
    t=np.float64(0.055)                 # thickness parameter (epsilon in the exercise)
    k=1.900
    k_naca0015,eps_naca0015=k,t  # np.float64(1.851),np.float64(0.1)            #1.851,0.1
    k=k_naca0015                               # trailing edge angle parameter
    epsilon=eps_naca0015                             # to avoid 0 division at theta=0 or 180 °
    eps=0.01
    npt=361                             # number of points on the airfoil
    #Kp_calculus_option=False;


   
    theta=np.linspace(-np.pi-eps,np.pi-eps,npt)
    theta=np.linspace(0+eps,2*np.pi-eps,npt)

    
    beta=1-k+epsilon*k                # parameter for perturbed velocity and Kp
    plot=True                           # to plot, set True
    flag_arctan=2                       # 3 types of arctan functions
    angle_flag=False                   # to plot the angles instead of the airfoil geometry

    set_question('1 : Profil')

    a=chord/2.0**k*(1.0+epsilon)**(k-1)   # parameter of the profil chord -> a 
    print('k = %f, epsilon = %f, a = %f '%(k,epsilon,a))

    r1=np.sqrt(np.sin(theta)**2+(np.cos(theta)-1)**2)
    r2=np.sqrt(np.sin(theta)**2+(np.cos(theta)-epsilon)**2)
    

    n=len(theta);
    print('n = %d '%n)

    t1=np.arctan2(np.sin(theta),(np.cos(theta)-1),)
    t2=np.arctan2(np.sin(theta),(np.cos(theta)-epsilon))
    

    theta1=np.zeros([n])
    theta2=np.zeros([n])

    for i in range(n):
        sint=np.sin(theta[i])
        cost=np.cos(theta[i])

        if theta[i]==0: 
            theta1[i]=np.pi/2
        else:
            theta1[i]=np.arctan(sint/(cost-1))+np.pi

        if flag_arctan==1:
            if cost-epsilon < 0 and  sint > 0 :
                theta2[i]=np.arctan(sint/(cost-epsilon))+np.pi     
            elif cost-epsilon<0 and  sint<0:
                theta2[i]=np.arctan(sint/(cost-epsilon))+np.pi 
            elif cost-epsilon >0 and sint<0: 
                theta2[i]=np.arctan(sint/(cost-epsilon))+2*np.pi 
            else:
                theta2[i]=np.arctan(sint/(cost-epsilon))
        elif flag_arctan==2:
            theta2[i]=arctan3(sint,cost-epsilon)

    if flag_arctan==3:
        theta2=t2
        theta1=t1

    
    phi=k*theta1-(k-1)*theta2
    r=a*r1**k/r2**(k-1)
    x=r*np.cos(phi)+ell;y=r*np.sin(phi)
    print('epaisseur maximale = %f'%(2*np.amax(y)))
    Imax=np.argmax(y)
    print('Imax = %i, x_max = %f'%(Imax,x[Imax]))
    if plot :
        # exemple de graphique
        nfig=1 
        plt.figure(nfig,figsize=(10,3))
        plt.title('profil', fontsize=14, fontweight='bold')
        #ax.set_title('relative error :'+ r'$\log_{10}(P_i / P_t-1)$')
        if angle_flag:
            plt.plot(theta,theta1,'-')
            plt.plot(theta,theta2,'--')
            plt.plot(theta,t1,'o')
            plt.plot(theta,t2,'s')
        else:
            if k==k_naca0015 and epsilon==eps_naca0015 :  
                xNaca,yNaca,thetaNaca=naca('0015',101,half_cosine_spacing = True)
                plt.plot(xNaca,yNaca,'-',linewidth=2,label='NACA 0015')
            #reference_filepath = os.path.join('Livre/Data', 'NACA0015.dat')
            #with open(reference_filepath, 'r') as infile:
            #   xNaca,yNaca=np.loadtxt(infile, dtype=float, unpack=True)
            plt.plot(x,y,'-',linewidth=2,label='Van de Vooren')
            plt.axis('equal')
            plt.xlim(-0.05,1.05)
            #ax.plot(theta,r)
        plt.xlabel(r'$x$',fontsize=20)
        plt.ylabel(r'$y$',fontsize=20)
        plt.grid()
        plt.legend(loc='upper left')
        plt.legend()
        
   
    set_question('2 : u,v and Kp')


    r3=np.sqrt(np.sin(theta)**2+(np.cos(theta)-beta)**2)
    t3=np.arctan2(np.sin(theta),np.cos(theta)-beta)
    THETA=theta+t3+(k-1)*theta1-k*theta2

    for alpha in list(aoa):
        Alpha=alpha/180.*np.pi
        
        #if Kp_calculus_option:     
        U=2*U0*(np.sin(Alpha)-np.sin(Alpha-theta))*r2**k/(r3*r1**(k-1))
        u=U*np.sin(THETA);v=-U*np.cos(THETA)
        # else:
        #     a1=np.cos((k-1)*theta1)*np.cos(k*theta2)+np.sin((k-1)*theta1)*np.sin(k*theta2)
        #     b1=np.sin((k-1)*theta1)*np.cos(k*theta2)-np.cos((k-1)*theta1)*np.sin(k*theta2)
        #     c1=(np.cos(k*theta2))**2+(np.sin(k*theta2))**2
        #     p=a*(1-k+k*epsilon)
        #     d1=a1*(a*np.cos(theta)-p)-b1*a*np.sin(theta)
        #     d2=a1*a*np.sin(theta)+b1*(a*np.cos(theta)-p)
        #     temp=2*c1*(np.sin(Alpha)-np.sin(Alpha-theta))/(d1**2+d2**2)
        #     com2=temp*(r2**k)/(r1**(k-1))
        #     u=d1*np.sin(theta)+d2*np.cos(theta)
        #     v=-(d1*np.cos(theta)-d2*np.sin(theta))
        #     U=U0*com2*np.sqrt(u**2+v**2)

             
        Kp=1-(U/U0)**2
        print('Kp minimale = %f'%(np.amin(Kp)))
        Cl=8*np.pi*a/ell*np.sin(Alpha)
        print('alpha = %2.0f, Cl = %1.3f'%(alpha,Cl))
        np.savetxt('Livre/Data/Kp_VanVooren_%i_deg.dat'%(alpha),np.c_[x,Kp],fmt='%f',delimiter=' ',newline='\n')
        if plot :
              
            # dessin des perturbations de vitesse
            nfig=nfig+1; plt.figure(nfig)
            plt.title('perturbation (u,v), alpha = %2.1f'%alpha, fontsize=14, fontweight='bold')
            plt.xlabel(r'$x/\ell$',fontsize=20)
            plt.ylabel(r'$u,v$',fontsize=20)
            plt.grid()
            plt.plot(x,u,'-',label=r'$\frac{u}{U_0}$',color='black',linewidth=2)
            plt.plot(x,v,'--',label=r'$\frac{v}{U_0}$',color='red',linewidth=2)
            plt.legend(loc='lower right')

            # dessin des perturbations de vitesse
            nfig=nfig+1; plt.figure(nfig)
            plt.title(r'Kp, alpha = %2.1f'%alpha, fontsize=14, fontweight='bold')
            plt.xlabel(r'$x$',fontsize=20)
            plt.ylabel(r'$Kp$',fontsize=20)
            plt.grid()
            #ax.legend(loc='upper left')
            plt.plot(x,Kp,'-')

            # dessin des perturbations de vitesse
            plt.figure(10)
            plt.title(r'Kp', fontsize=14, fontweight='bold')
            plt.xlabel(r'$x$',fontsize=20)
            plt.ylabel(r'$Kp$',fontsize=20)
            plt.grid()
            plt.axis([-0.05,1.05, -5, 1.2])
            plt.plot(x,Kp,'-',label=r'$\alpha =$ %3.0f'%(alpha))
            plt.legend(loc='lower right')
        

    plt.show() 

def Kp_profil_parabolique(theta,alpha,t):
    """
    Kp du profil parabolique
    t = cambrer/chord
    """    
    a0=-4*alpha
    a1=-8*t
    print('a0 = %2.3f, \t a1= %2.3f '%(a0,a1))
    Cl=-np.pi/2*(a0+a1)
    #Kp=-2*np.sin(theta)*(4*t+alpha/(1+np.cos(theta)))
    Kp=np.sin(theta)*(a1+a0/2/(1+np.cos(theta)))

    return Cl,Kp
    
def Kp_epaisseur(theta,t):
    """
    Kp du profil symmétrique épais
    """    
    b0=16/9*np.sqrt(3)*t
    b1=-b0
    print('b0 = %2.3f, \t b1= %2.3f '%(b0,b1))
    Kp=-(b0/2 + b1*np.cos(theta))
    return Kp
    

def Exercice7_2():
    """ 
    Profil parabolique
    """
    set_title("Profil parabolique : solution analytique")

    # data, parameters
    epaisseur=True
    U0=np.float64(1)                    # reference velocity
    chord=np.float64(1.)                # chord length
    ell=chord
    aoa=np.array([0.,5.0,10.0])/180*np.pi;             # angle of attack
    c=np.float64(0.1)             # camber parameter e/L
    eps=1/180*np.pi                             # to avoid 0 division at theta=0 or 180 ° 
    npt=181                             # number of points on the airfoil
    plot=True                           # to plot, set True
    t=0.1                              # thickness parameter          

    theta=np.linspace(0,np.pi-eps,npt)
    x=np.cos(theta)/2              # x/ chord
    yc=c*(1-4*x**2)
    if epaisseur:
        ye=2*t/np.sqrt(27)*np.sqrt(-16*x**4+16*x**3-4*x+1)
        KP_epais=Kp_epaisseur(theta,t)
    n=len(theta);              
    m=len(aoa)   
    KP=np.zeros([n,m])
    Cl=np.zeros([m])
    nfig=0

    if plot :
        # exemple de graphique
        fig = plt.figure(0)
        fig.suptitle('profil', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        ax.grid()    
        ax.plot(x,yc,'-',linewidth=2)
        if epaisseur :
            ax.plot(x,yc+ye,'-',linewidth=2)
            ax.plot(x,yc-ye,'-',linewidth=2)
            ax.plot(x,ye,'-',linewidth=2)
      
        ax.axis([-0.5, 0.5, -0.25, 0.25])
        ax.axis('equal')
             

    for i in range(m):
        Cl[i],KP[:,i]=Kp_profil_parabolique(theta,aoa[i],c)
        print(' alpha = %2.2f ° \t CL = %2.2f '%(aoa[i]/np.pi*180,Cl[i]))

        if plot :  
            fig = plt.figure(1)
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            #ax.set_title('Kp')
            ax.set_xlabel(r'$x/L$',fontsize=20)
            ax.set_ylabel(r'$Kp$',fontsize=20)
            ax.grid()
            ax.axis([-0.5, 0.5, -4, 0])
            ax.plot(x,KP[:,i])
            if epaisseur:
                fig1 = plt.figure(2)
                ax1 = fig1.add_subplot(111)
                fig.subplots_adjust(top=0.80)
                #ax.set_title('Kp')
                ax1.set_xlabel(r'$x/L$',fontsize=20)
                ax1.set_ylabel(r'$Kp$',fontsize=20)
                ax1.grid()
                ax1.axis([-0.5, 0.5, -4, 4])
                ax1.plot(x,KP[:,i]+KP_epais[:])
                ax1.plot(x,-KP[:,i]+KP_epais[:])
             

    if plot :        
        plt.show()

#********************************************************************
def Exercice7_3():
#********************************************************************
    """
    Discrete Vortex Method on the parabolic mean line airfoil
    Le code vient de la traduction du code original de Katz et Plotkins qui
    ont des circulations positives pour une portance positive
    Dans nos conventions c'est le contraire ...
    """

   


    def mean_line(x):
        """
        equation of the mean line
        t : thickness parameter
        """
        if TestCase == 1 :
            return 4.*t*x/chord*(1.-x/chord)
        else:
            return t*(1-2*x/chord)*(1.+2*x/chord)

    def slope(x):
        """
        slope of the mean line dy/dx
        """
        if TestCase == 1 :
            return 4.*t/chord*(1.-2.*x/chord);
        else:
            return -8*t/chord*x/chord
        
    def normal(x):
        """
        normal vector components from the mean line equation
        """
        detadx=slope(x);sq=np.sqrt(1+detadx**2);
        nx,ny= -detadx/sq, 1./sq
        return nx,ny

    def Aero_Reference(x,alpha):
        """
        Pressure coefficient and lift coefficient, reference case
        """
        if TestCase == 1 :
            temp=16*t/chord*np.sqrt(x/chord*(1-x/chord));
            Kp_ref=2.*np.sqrt((chord-x)/x)*alpha+temp;
            CL_ref=2.*np.pi*(alpha+2*t/chord)
        else:
            Kp_ref=2*np.sqrt(1-4*(x/chord)**2)*(4*t/chord+alpha/(1+2*x/chord))
            CL_ref=2*np.pi*(alpha+2*t/chord)
            Cm0=np.pi/2*(alpha+4*t/chord)
        return Kp_ref,CL_ref,Cm0
     
    def solve_circulation(nx,ny):
        """
        Circulation for a vortex point distribution
        """
        A = np.zeros((N, N), dtype=float)  
        RHS = np.zeros((N), dtype=float) 
        for i in range(N):
            for j in range(N):
                x0,y0=xc[i]-xv[j],yc[i]-yv[j]
                r=np.sqrt(x0**2+y0**2)
                if  r >= 0.001:
                    tmp=1/(2*np.pi*r)
                    u,v=-tmp*y0/r,tmp*x0/r
                else:
                    u,v=0,0
                    print('u=v=0')
                A[i,j]=u*nx[i]+v*ny[i];
                
            #   the rhs vector 
            RHS[i]=-uinf*nx[i]-vinf*ny[i];
        gamma = np.linalg.solve(A, RHS)
        return gamma

    def point_location_local(N,x0,dx):
        """
        the points are defined locally from the mean line equation and its local slope
        """
        # initialization
        xc = np.zeros(N, dtype=float); yc = np.zeros(N, dtype=float) 
        xv  = np.zeros(N, dtype=float); yv  = np.zeros(N, dtype=float) 
        nx = np.zeros(N, dtype=float); ny = np.zeros(N, dtype=float) 

        xc = x0+dx*(np.float64(range(N))+0.75)           # control point
        xv = x0+dx*(np.float64(range(N))+0.25)           # vortex point
        yc,yv = mean_line(xc),mean_line(xv)
        nx,ny=normal(xc)
        print('local: len(xc)',len(xc))
        if display:
            print('Local: xc =' ,xc);print('xv =' ,xv)

        return xc,yc,xv,yv,nx,ny


    def point_location_panel(N,x0,dx,distribution):
        """
        the points are defined locally from the panels 
        """
        # initialization
        xa = np.zeros(N, dtype=float); ya = np.zeros(N, dtype=float) 
        xb = np.zeros(N, dtype=float); yb = np.zeros(N, dtype=float) 
        xc = np.zeros(N, dtype=float); yc = np.zeros(N, dtype=float) 
        xv  = np.zeros(N, dtype=float); yv  = np.zeros(N, dtype=float) 
        nx = np.zeros(N, dtype=float); ny = np.zeros(N, dtype=float)
        if distribution == 'cos':
            theta= np.zeros(N+2, dtype=float)
            eps=5/180*np.pi
            theta=np.linspace(np.pi-eps,eps,N+2)
            if TestCase==1:
                xa=(np.cos(theta[0:N])+1)/2
                xb=(np.cos(theta[1:N+1])+1)/2
            else:
                xa=np.cos(theta[0:N])/2
                xb=np.cos(theta[1:N+1])/2
        else:
            xa,xb=x0+dx*np.float64(range(N)),x0+dx*(np.float64(range(N))+1)
        
        ya=mean_line(xa);yb=mean_line(xb)
        dx_out=xb-xa
        xc,xv=xa+0.75*dx_out,xa+0.25*dx_out  # control point,   vortex point         
        yc,yv=ya+0.75*(yb-ya),ya+0.25*(yb-ya)  # control point,   vortex point 
        #yc,yv = mean_line(xc),mean_line(xv)

        beta=np.arctan2(yb-ya,xb-xa)
        nx,ny=-np.sin(beta),np.cos(beta)
        print('Panel: len(xc)',len(xc))
        if display:
            print('Panel: xc =' ,xc);print('xv =' ,xv)

        return xc,yc,xv,yv,nx,ny,dx_out    

    def Aerodynamic_load(gamma,dx,x,alpha):
        """
        Aerodynamic coefficient K_p and C_L
        """ 
        dL=-rho*Vinf*gamma;             # elementary lift per panel
        print('x0',x0,alpha)
        Cm0=np.sum((x-x0)/chord*dL/(q*chord))*np.cos(alpha)
        Kp=dL/(2*dx*q);                 # lower surface Cp^- per span unit
        Lift=np.sum(dL);
        CL=Lift/(q*chord);
        return Kp,CL,Cm0   

  

    def plot_geometry():
        """
        plot the geometry
        """
       
        # one windows
        fig, (ay1, ay2) = plt.subplots(1,2, figsize=(20, 10))
        lines=ay1.plot(xv, yv, color='r', linestyle='-', marker='s', markersize=6,linewidth=2,label='local')
        lines=ay1.plot(xvp,yvp,linestyle='-', linewidth=1, marker='o', markersize=6, color='blue',label='panel')
        ay1.set_ylabel(r'$y_v$', fontsize=16)
        ay1.set_xlabel(r'$x_v$', fontsize=16)
        #lgd = ay1.legend((lines), ('local','panel'), loc=3,ncol=3, bbox_to_anchor=(0, 1), fontsize=16)
        ay1.legend()
        ay1.set_xlim([xmin, xmax])
        ay1.set_ylim([0, 1.1*t])
        ay1.set_title('airfoil mean line')
        ay1.grid()

        lines = ay2.plot(xv, slope(xv), color='r', linestyle='-', marker='s', markersize=6,linewidth=2)
        lines = ay2.plot(xvp,slope(xvp),linestyle='-', linewidth=1, marker='o', markersize=6, color='blue')
        #ay2.set_ylabel(r'$\frac{dy}{dx}$', fontsize=16)
        ay2.set_xlabel(r'$x_v$', fontsize=16) 
        ay2.set_title('airfoil mean line slope')
        ay2.grid()
         
        ay2.set_xlim([xmin, xmax])
        plt.show(lines)

    #****************************************************
    # Main program
    #****************************************************
    set_title("Profil parabolique : tourbillons ponctuels")

    TestCase=2                      # = 1 : parabolic profile x in [0 1]
                                    # = 2 : parabolic profile x in [-1/2 1/2]
    width,height = 10,8                      # plot size
    # upstream flow
    Vinf=np.float64(1)                    # reference velocity
    rho=1.4
    alpha_1=10                            # angle of attack in degrees
    print('Angle of attack  = %2.1f °'%(alpha_1))
    # Airfoil
    chord=np.float64(1.)                # chord length
    t=np.float64(0.1)*chord             # camber parameter  (sometimes called c)
    if TestCase == 1:
        x0=0
    else:
        x0=-chord/2
    # Disctretization with panels
    N=90                             # number of points on the airfoil 
    dx=chord/N                         # constant step size
    #                                   
    print('dx = %f '%(dx))

    if N < 15 :
        display=True
    else:
        display=False
    alpha=np.pi/180*alpha_1;
    # velocity projection
    uinf=np.cos(alpha)*Vinf;
    vinf=np.sin(alpha)*Vinf;
    q=0.5*rho*Vinf**2;                  # dynamic pressure

    # first option : slope are calculated directly from the analytical equation
    xc,yc,xv,yv,nx,ny=point_location_local(N,x0,dx)
    # second option the normal is defined from panels, 
    # all the points are differents in the panel
    # possible distribution "regular (reg) or cosine (cos)
    # with cosine it does not work properly. Why ?
    xcp,ycp,xvp,yvp,nxp,nyp,dxp=point_location_panel(N,x0,dx,'reg')

    if xv[N-1] >= 0.5:
        xmin,xmax=0,chord
    else:
        xmin,xmax=-chord/2,chord/2


    #
    plot_geometry()

    gamma=solve_circulation(nx,ny)
    gammap=solve_circulation(nxp,nyp)

    

    # Aerodynamic load

    Kp,CL,Cm=Aerodynamic_load(gamma,dx,xv,alpha)
    Kpp,CLp,Cmp=Aerodynamic_load(gammap,dxp,xvp,alpha)
    Kp_ref,CL_ref,Cm_ref=Aero_Reference(xv,alpha) # analytic solution

    print('CL Local                        = %12.8f '%(CL))
    print('CL Panel                        = %12.8f '%(CLp)) 
    print('CL   (exact)                    = %12.8f '%(CL_ref))
    print('Cm Local                        = %12.8f '%(Cm))
    print('Cm Panel                        = %12.8f '%(Cmp))
    print('Cm   (exact)                    = %12.8f '%(Cm_ref))

     

    # plot circulation
    plt.figure(figsize=(width, height))
    plt.grid()
    plt.title(r'Circulation $\gamma(x)$')
    plt.xlabel(r'$x_v$', fontsize=16)
    plt.ylabel(r'$\gamma$', fontsize=16)
    plt.plot(xv,gamma, color='b', linestyle='-', marker='o', markersize=6,linewidth=2,label='local')
    plt.plot(xvp,gammap, color='r', linestyle='--', linewidth=2,label='panel')
    plt.plot(xv,-Kp_ref*Vinf*dx, color='k', linestyle='-', linewidth=2,label='reference')
    plt.legend(loc='lower right')
    plt.xlim(xmin, xmax)


    # plot Kp
    plt.figure(figsize=(width, height))
    plt.grid()
    plt.title(r'K_p$')
    plt.xlabel(r'$x_v$', fontsize=16)
    plt.ylabel(r'$Kp_i$', fontsize=16)
    plt.plot(xv, Kp_ref, color='k', linestyle='-', linewidth=3,label='reference')
    plt.plot(xv,Kp,linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305',label='local')
    plt.plot(xvp,Kpp,linestyle='--', linewidth=1,  color='blue',label='panel')
    plt.legend()
    plt.xlim(xmin, xmax)
    plt.show()

#********************************************************************
def Exercice7_4():
#********************************************************************
    """ 
    Profil de Van-Hooren
    """
    set_title("Profil de Van-Hooren : solution analytique, 2nd version")
    set_title("Distribution ponctuelle de sources")

    # data, parameters
    U0=np.float64(1)                    # reference velocity
    chord=np.float64(1.)                # chord length
    rho0=1.4                            # air densitu (not used here)
    ell=chord           
    aoa=np.array([0.]);                 # angle of attack = 0 here, thickness effect
    #t,k=np.float64(0.02),1.85                   # thickness parameter 
    t,k=0.055,1.9                         # trailing edge angle parameter
    npt=181                              # number of points on the airfoil and of control points
    
    """
     Options de méthodologie, de calcul et de visualisations :
    """

    show_control_points     = False             # to have some intermediate plot  set True
    symmetrie               = False             # solve the control point only on the upper wall, otherwise on the full airfoil surface
    plot_streamlines        = True
    plot_vectors            = False             # to plot the slip velocity on the airfoil surface
    plot_Kp                 = True              # plot the pressure coefficient
    plot_profil             = False             # geometry    
    N_display               = 31                # a parameter to avoid to many data on screen
    zero_source_sum         = True              # to impose sum s_k = 0
    row_source              = 0                 # ou int(np.ceil(npt/3)), ce paramètre n'a pas trop d'importance en fait, on peut mettre 0 aussi 
    plot_sources            = True              # to plot the source distribution
    improved_control_point  = False             # change the control point at the leading and trailing edge if symmetry is used
    more_details            = False 
    translation_xc          = np.float(0.1)     # translate the first and last source points at the L.E. and T.E.


    # pour bien fonctionner mettre vrai le test_somme_source et ne pas utiliser la symétrie


    fst = Freestream(U_inf=U0, alpha=aoa,rho=rho0)  # écoulement amont

    # Geometry
    x,y,theta=Van_Hooren_gmtry(chord,t,k,npt,plot_profil,dtheta=0)
    

    n=int((npt-1)/2)
    xe,ye=x[0:n+1],y[0:n+1]
    xi,yi=xe,-ye                # extrados, intrados
    print("longueur de xe = %i, n= %i "%(len(xe),int(n)))
    u=(0,n);print(" [xe_BA, xe_BF], [ye_BA,ye_BF] ", xe[list(u)],ye[list(u)])


    if plot_profil:
        fig = plt.figure(figsize=(10,5)); fig.suptitle('profil', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111); fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$',fontsize=20); ax.set_ylabel(r'$y$',fontsize=20)
        ax.plot(xe,ye,'r-',linewidth=2); ax.plot(xi,yi,'b-',linewidth=2)
        ax.grid();ax.axis('equal')
       

    # définition des panneaux

    if symmetrie:
        # create panels
        print("On utilise la symmétrie, calcul sur l'extrados")
        N=len(xe)-1;print('nombre de panneaux = %i'%(N))
        panels = np.zeros(N, dtype=object)
        for i in range(N):
            panels[i] = Panel(xe[i], ye[i], xe[i+1], ye[i+1])
        if more_details:
            print("{%s} \t \t {%s} \t  {%s} "%('i', 'delta (°)', 'x_control'))
            for i in range(N):
                print('{%4i} \t \t {%06.2f} \t {%+06.4f} '%(i,np.rad2deg(panels[i].delta),panels[i].xc))
        if improved_control_point:
            print('Improved BC')
            i=0;
            print('Before change at L.E.       : xc = %f'%(panels[i].xc))
            panels[i] = Panel(xe[i], ye[i], xe[i+1], ye[i+1],coef=0.25)
            print('After change at L.E.        : xc = %f'%(panels[i].xc))
            
            i=N-1;
            print('before change at T.E.       : xc = %f'%(panels[i].xc)) 
            panels[i] = Panel(xe[i], ye[i], xe[i+1], ye[i+1],coef=0.95)
            print('After change at T.E.        : xc = %f'%(panels[i].xc)) 
            
    else:
        print("On n'utilise pas la symmétrie, calcul sur l'extrados et l'intrados")
        # create panels
        N=len(x)-1;print('nombre de panneaux = %i'%(N))
        panels = np.zeros(N, dtype=object)
        for i in range(N):
            panels[i] = Panel(x[i], y[i], x[i+1], y[i+1])
            if N<N_display:
                print('i = %3i, delta = %3.3f °, x = %f '%(i,panels[i].delta*180/np.pi,panels[i].xc))

    if show_control_points:
        # plot discretized geometry
        width = 10
        plt.figure(figsize=(width, 0.5*width))
        plt.grid()
        plt.title("Points de contrôle sur les panneaux")
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.plot(xe, ye, color='k', linestyle='-', linewidth=2)
        plt.plot([panel.xc for panel in panels],[panel.yc for panel in panels],
                    linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305')
        plt.axis('scaled', adjustable='box')
        ech=0.1
        plt.xlim(0-ech, chord+ech)
        ymax=np.max(y)+ech
        plt.ylim(-ymax, ymax);

    # position des sources sur l'axe

    Ns=len(xe)-1
    xs = np.zeros(Ns, dtype=float); ys = np.zeros(N, dtype=float) 
    
    for i in range(Ns):
        xs[i]=np.copy(panels[i].xc)
    if more_details: print('Position des sources xs=',xs)
    
    
    if symmetrie and improved_control_point:
        #dx_tmp=translation_xc*(xs[1]-xs[0]);
        dx_tmp=chord/(npt-1)*translation_xc
        print('translation de xs at L.E. =%f  at xs = %f'%(dx_tmp,xs[0]))
        print('xs = %f %f'%(xs[0],xs[Ns-1]))
        xs[0]+=dx_tmp
        #xs[Ns-1]-=dx_tmp*0.5   # possibilité
        print('xs = %f %f'%(xs[0],xs[Ns-1]))

    # la dernière ligne de A permet d'imposer la somme des sources = 1
    # idem pour la dernière ligne de RHS
    A,B = np.zeros((N, Ns), dtype=float),np.zeros((N, Ns), dtype=float)
    RHS,RHSt = np.zeros((N), dtype=float),np.zeros((N), dtype=float)
    Ut,Kp=np.zeros((N),dtype=float),np.zeros((N),dtype=float)
    xp = np.zeros((N), dtype=float)        # position du centre des panneaux

    # source contribution on a panel from others
    for i, panel_i in enumerate(panels):
        for j in range(Ns):
            x0,y0=panel_i.xc-xs[j],panel_i.yc-ys[j]
            r2=x0**2+y0**2
            r=np.sqrt(r2)
            if  r >= 0.00001:
                tmp=1/(2*np.pi*r2)
                u,v=tmp*x0,tmp*y0
            else:
                u,v=0,0;print('u=v=0')
            A[i,j]=u*panel_i.nx+v*panel_i.ny;
            B[i,j]=u*panel_i.tx+v*panel_i.ty;
            # print(i,j,x0,y0,panel_i.nx,panel_i.ny,A[i,j])
            
        xp[i]=panel_i.xc   
        RHS[i]=-fst.u_inf*panel_i.nx-fst.v_inf*panel_i.ny
        RHSt[i]=fst.u_inf*panel_i.tx+fst.v_inf*panel_i.ty
        #print('RHS = ',RHS[i])
        #print(fst.u_inf,panel_i.nx,fst.v_inf,panel_i.ny)
    if zero_source_sum:
        RHS[row_source]=0;A[row_source,:]=1
        print('On impose sum sigma_k = 0 sur la ligne %i, à x = %f'%(row_source,xs[row_source]))
    if Ns==N:
        print('inversion directe')
        S=np.linalg.solve(A, RHS)
    else:
        print('Résolution par les moindres carrés')
        S,res,rank,s = np.linalg.lstsq(A, RHS)

    if more_details: print('vecteur des sources : S =',S )
     
    print("somme ses S_k           = %3.10f"%(np.sum(S)))
    error=np.dot(A,S)-RHS
    print("erreur d'inversion du système N_2 = %12.10e"%(np.sqrt(np.dot(error,error))))
    if plot_sources:
        ymax=np.max(np.abs(S))
        ymax=np.max(S)
        ymin=np.min(S)
        width = 10
        plt.figure(figsize=(width, 0.8*width))
        plt.grid()
        plt.xlabel('x', fontsize=16)
        plt.ylabel('q', fontsize=16)
        plt.title('Sources')
        plt.plot(xs,S,linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305')
        ech=0.05
        plt.xlim(0-ech, chord+ech)
        #ymax=0.05
        plt.ylim(ymin, ymax);

    # calcul de la vitesse tangentielle au point de contrôle 
    Ut=np.dot(B,S)+RHSt    
    # coefficient de pression
    Kp=1-Ut**2                  # coefficient de pression
    for i,panel in enumerate(panels):
        panel.Kp,panel.Ut=Kp[i],Ut[i]
    
    if plot_Kp:
         # load geometry from data file
        reference_filepath = os.path.join('Livre/Data', 'Constant_S_Source.dat')
        with open(reference_filepath, 'r') as infile:
            xKatz,CpKatz,ind=np.loadtxt(infile, dtype=float, unpack=True,skiprows=2)
        filepath_ref = os.path.join('Livre/Data', 'Kp_VanVooren_0_deg.dat')
        with open(filepath_ref, 'r') as infile:
            xref,Cpref=np.loadtxt(infile, dtype=float, unpack=True)    
        width = 10
        plt.figure(figsize=(width, 0.8*width))
        plt.grid()
        plt.xlabel('x', fontsize=16)
        plt.ylabel('Kp', fontsize=16)
        plt.title('Coefficient de pression')
        #plt.plot(xp,Kp,linestyle='-', linewidth=1, marker='x', markersize=6, color='#CD2305')
        plt.plot(xp,Kp,linestyle='-', linewidth=1,  color='#CD2305',label='Distribution ponctuelle')
        plt.plot((xKatz+1)/2,CpKatz,linestyle='', linewidth=1, marker='o', markersize=4,markevery=2, color='black',label='Distribution uniforme (Katz)')
        plt.plot(xref,Cpref,linestyle='-', linewidth=1, color='blue',label='ref. analytique')
        plt.legend()
        ech=0.05
        plt.xlim(0-ech, chord+ech)
        plt.ylim(-1.1, 1.1);


    if plot_streamlines:
        X,Y=definition_grille(100,-0.2,1.2,-0.25,0.25)
        champ_vitesse(fst.u_inf,fst.v_inf,X,Y,xs,S,plot_streamlines,x,y)

    if plot_vectors:    
        # Points de contrôle :
        xc,yc = np.zeros((N), dtype=float),np.zeros((N), dtype=float)
        for i, panel_i in enumerate(panels):
            xc[i],yc[i]=panel_i.xc,panel_i.yc 
        champ_vitesse_vecteur(fst.u_inf,fst.v_inf,xc,yc,xs,S,plot_vectors,x,y)
    
       
    plt.show()



#********************************************************************
def Exercice7_10():
#********************************************************************
    """ 
    Profil NACA 
    """
    set_title("Profil NACA")
    number='4412'
    npt=41
    x,y,theta=naca(number,npt,half_cosine_spacing = True)
    plt.figure(figsize=(10, 10))
    plt.plot(x,y,linestyle='-', linewidth=2, marker='o', markersize=3, color='#CD2305')
    plt.xlabel('x')
    plt.ylabel('y') 
    plt.title('NACA 4  %s'%(number))
    plt.grid()
    plt.axis('equal')
    plt.xlim(-0.05,1.05)
    plt.show()
