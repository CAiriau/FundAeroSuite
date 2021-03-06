
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:40:21 2018

@author: cairiau
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

 

def titre(s):
    print(50*"=","\n","  %s\n"%(s),49*"=","\n")

def plot_u(y,Ysol,Titre="Solution pour u",lim=np.array([0,5]),scale="linear",display=False):
    """
    dessin de la solution
    """
    if display:
        plt.figure(figsize=(12,10))
        plt.title(Titre)
        plt.xlabel(r'$x$',fontsize=16)
        plt.grid()
        #ax.axis([-1,1,0,1])
        plt.xlim(lim[0],lim[1])
        
        if scale=="linear":
            plt.plot(y,Ysol,'-', linewidth=2,label="u")
        elif scale=="semilogy":
            plt.semilogy(y,Ysol,'-', linewidth=2,label="u")
        
        plt.legend(loc="lower left")
        #plt.show()
        

def plot_solution(y,Ysol,Titre="Solution pour Ys",lim=np.array([0,5]),scale="linear",display=False):
    """
    dessin de la solution
    """

    if display:
        t=Ysol.shape
        print("dimension des données                    : %i x %i "%(t[0],t[1]))
        N=int(t[1]/2) 
        plt.figure(figsize=(12,10))
        plt.title(Titre)
        plt.xlabel(r'$y$',fontsize=16)
        #ax.set_ylabel(r'$f_k$',fontsize=16)
        plt.grid()
        #ax.axis([-1,1,0,1])
        plt.xlim(lim[0],lim[1])
        
        if scale=="linear":
            for k in range(2*N): 
                plt.plot(y,Ysol[:,k],'-', linewidth=2,label="k = %i"%(k))
        elif scale=="semilogy":
            for k in range(2*N): 
                plt.semilogy(y,Ysol[:,k],'-', linewidth=2,label="k = %i"%(k))
        plt.legend(loc="best")
        #plt.show()
        
def calcul_system_coeff(N,display=False,eps=1):
    """
    calcul de la série pour le système ODE
    """
    testnum=False
    if display: titre("Calcul des coefficients de l'ODE")

    s=[]
    
    if display: print ("N = ",N)
    Lin=[]
    TNL=[]
    for n in np.arange(1,N+1):
        c=[]
        if display : print("ligne: %i, n= %i, linéaire : \t  %f f_%i "%(n-1,n,n*n,n-1))
        Lin.append(n*n)   # attention : terme linéaire pas quadratique
        for m in range(n+1):
            if testnum:
                c1=0
            else:
                c1=-np.float(n/4*m*(n-m))*eps
            if np.abs(m*c1)>0:  
                if display: print("ligne: %i, n= %i, m = %i \t  %f x f_%i f_%i"%(n-1,n,m,c1,m-1,n-m-1))
                c.append([m-1,n-m-1,c1])
            else:
                #if display: print("n= %i, m = %i : out of limit  "%(n,m))
                pass
            
        if display: print("")
        for m in range(N+1):
            if testnum:
                c2=0
            else:
                c2=-np.float(n/2* (m+1)*(n+m+1))*eps
            test=m+n > N-1
            if c2 < 0 and test==False:
                if display : print("ligne: %i,n= %i, m = %i : \t  %f x f_%i f_%i"%(n-1,n,m,-n/2* (m+1)*(n+m+1),m+1-1,n+m+1-1))
                c.append([m+1-1,n+m+1-1,c2])
            else:
                #if display: print("n= %i, m = %i : out of limit  "%(n,m))
                pass
        s.append(len(c))
        TNL.append(c)  
        if display: print(50*"-")
        
    if display :   
        print('nombre de coefficients par valeur de n   :',s)
        print('nombre total de  coefficients            : %i'%(sum(s)))
        print ("termes linéaires                        :",Lin)
    return TNL,Lin,s
    



def derivative_sys(TNL,s,display=False):
    """
    Gradient of the system with respect to f_k
    """
    if display : titre("calcul du gradient du système direct")
    N=len(s);
    delta_TNL=[]
    for n in range(N):
        if display: print(' \n Ligne n° %i :\n'%(n+1))
        coef=TNL[n]
        c=np.zeros((N,N))
        for L in coef:
            #print(L)
            i,j,a=L[0],L[1],L[2]
            if i==j==n:
                if display: print("terme          f_%i^2"%(i))
                c[i,j] =2*a
            else:
                if display: print("terme          %f  f_%i x f_%i "%(a,i,j))
                c[i,j]+=a
                c[j,i]+=a
        if display: print("Matrice c = \n",c)
        delta_TNL.append(c)
        if display: 
            print("delta_TNL = ")
            for L in delta_TNL: print(L)
    return delta_TNL

def sys_direct(Y,y,TNL,Lin,display=False):
    """
    Système direct
    """
    N=len(Lin)
    dY=np.zeros(2*N)
    dY[np.arange(0,2*N,2)]=Y[np.arange(1,2*N+1,2)]  # composante paire : fonction
    if display :
        print("Index dérivée : Y",np.arange(0,2*N,2),"' =  Y", np.arange(1,2*N+1,2))
    for n in range(N):
        dY[2*n+1]= Lin[n]*Y[2*n]          # partie linéaire de la dérivée de la composante impaire
    for n in range(N):
        coef=TNL[n]
        if display : print("TNL_%i      ="%(n),coef)
        for L in coef:
            i,j,a=L[0],L[1],L[2]
            dY[2*n+1]+=a*Y[2*i]*Y[2*j]   # f"_n = somme (a f_i f_j)
    return dY    

def grad_sys(Y,y,TNL,delta_TNL,Lin,display=False):
    """
    système linéaire, gradient du système précédent par
    rapport à une variable s
    et simultanément système direct
    """
    
    N=len(Lin); Nt=N*(N+1)
    dY=np.zeros(2*Nt)
    dY[np.arange(0,2*Nt,2)]=Y[np.arange(1,2*Nt+1,2)]  # composante paire : fonction
    if display :
        print("Index dérivée : Y",np.arange(0,2*Nt,2),"' =  Y", np.arange(1,2*Nt+1,2))
        print("dY, IC: ",dY)
        print("Lin :",Lin)
    for n in range(N):
        u=np.arange(2*n+1,2*Nt,2*N)
        dY[u]= Lin[n]*Y[u-1]          # partie linéaire de la dérivée de la composante impaire
    
        if display: print('Index termes linéaires : ',u)
    dY1=np.copy(dY);
    if display: print("dY linéaire     = ",dY)
    for n in range(N):
        coef=TNL[n]
        for L in coef:
            i,j,a=L[0],L[1],L[2]
            dY[2*n+1]+=a*Y[2*i]*Y[2*j]   # f"_n = somme (a f_i f_j)
    # fin du système direct
    # Y(0:2*N)     : champ direct Y
    # Y(2N+1,4N)  : dY/ds_1
    # Y(4N+1,6N)   : dY/ds_2
    if display: print("dY non linéaire = ",dY)
    # cette partie est très compliquée pour trouver les bons indices...
    for n in range(N):
        c=delta_TNL[n].dot(Y[np.arange(0,2*N,2)])
        w=np.arange(2*(N+n)+1,2*Nt,2*N)        # position de la ligne pour Y'
        for k in range(N):
            v=2*N*(k+1)+np.arange(0,2*N,2)
            dY[w[k]]+=c.dot(Y[v]) 
            if display: print("c(%i).Y("%(n),v,")= Y'[%i]"%w[k])        
    if display: print('delta dY         = ',dY-dY1,' end')
    return dY


def solve_syst_direct(s,y,TNL,delta_TNL,Lin,parameters=[0,0],display=False):        
    """
    résolution du problème : relation entre les sorties (target) 
    et les entrées (s)
    """
    if display: print(50*"=","\n"," Système direct seul : \n",50*"=","\n")
    k0,N=parameters[0],len(Lin)
    if len(s) == N :
        Y0=np.zeros(2*N)
        Y0[1]=-k0
        Y0[np.arange(0,2*N,2)]=s
        if display:
            print('IC :  Y0                                 : ',Y0)
            print('IC :  s0                                 : ',s)
    else:
        raise ValueError("Mauvaise dimension du paramètre s")

    Ys=odeint(sys_direct,Y0,y,args=(TNL,Lin,),Dfun=None,full_output=False)
    if parameters[1]== 0:
        target=Ys[-1,np.arange(1,2*N,2)]
    elif parameters[1]== 1:
        target=Ys[-1,np.arange(1,2*N,2)]+np.arange(1,N+1)*Ys[-1,np.arange(0,2*N,2)]
        if display: print("target  indices : ",np.arange(1,2*N,2),np.arange(1,N+1),np.arange(0,2*N,2) )
    else:
        w=np.exp(np.arange(1,N+1)*y[-1])
        #w=10**(np.arange(0,N))
        target=w*(Ys[-1,np.arange(1,2*N,2)]+np.arange(1,N+1)*Ys[-1,np.arange(0,2*N,2)])
        if display: print("target  indices : ",np.arange(1,2*N,2),np.arange(1,N+1),np.arange(0,2*N,2) )
    
    
    if display:
        print("size of the solution Ys                  :",Ys.shape)
        print("f_k(y_max)                               :",target)
        print("Y(y_max)                                 :",Ys[-1,:])
        print("Y(0)                                     :",Ys[ 0,:])
    return target,Ys

def solve_syst_direct_gradient(s,y,TNL,delta_TNL,Lin,parameters=[0,0],display=False):        
    """
    résolution du problème : relation entre les sorties (target) 
    et les entrées (s),
    on résoud en même temps les gradients
    """
    if display: print(50*"=","\n"," Système direct + gradient : \n",50*"=","\n")
    k0,N=parameters[0],len(Lin)
    Nt=N*(N+1)
    
    if len(s) == N :
        Y0=np.zeros(2*Nt)
        Y0[1]=-k0
        Y0[np.arange(0,2*N,2)]=s
        ind= 2*(np.arange(0,N)*(N+1)+N)
        Y0[ind]=1
        if display:
            print("index pour delta Y(0)=1                  : ",ind)
            print('IC :  Y0                                 : ',Y0)
            print('IC :  s0                                 : ',s)
    else:
        raise ValueError("Mauvaise dimension du paramètre s")
    
    
    Ys=odeint(grad_sys,Y0,y,args=(TNL,delta_TNL,Lin),Dfun=None,full_output=False)
    if parameters[1]==0:
        target=Ys[-1,np.arange(1,2*N,2)]
        grad_target=Ys[-1,np.arange(2*N+1,2*Nt+1,2)] 
    elif parameters[1]==1:
        target=Ys[-1,np.arange(1,2*N,2)]+np.arange(1,N+1)*Ys[-1,np.arange(0,2*N,2)]
        grad_target=Ys[-1,np.arange(2*N+1,2*Nt+1,2)]+np.tile(np.arange(1,N+1),N)*Ys[-1,np.arange(2*N,2*Nt,2)]
        if display:        
            print("target  indices : ",np.arange(1,2*N,2),+np.arange(1,N+1),np.arange(0,2*N,2) )
            print("grad target ind : ",np.arange(2*N+1,2*Nt+1,2),np.tile(np.arange(1,N+1),N),np.arange(2*N,2*Nt,2) )
    else:
        w=np.exp(np.arange(1,N+1)*y[-1])
        #w=10**(np.arange(0,N))
        target=w*(Ys[-1,np.arange(1,2*N,2)]+np.arange(1,N+1)*Ys[-1,np.arange(0,2*N,2)])
        grad_target=np.tile(w,N)*(Ys[-1,np.arange(2*N+1,2*Nt+1,2)]+np.tile(np.arange(1,N+1),N)*Ys[-1,np.arange(2*N,2*Nt,2)])
       
    
    if display:
        print("size of the solution Ys                  :",Ys.shape)
        print("f_k(y_max)                               :",target)
        print("grad (f_k(y_max))                        :",grad_target)
        for k in range(N+1):
            print("Y(0)                                     :",Ys[0,k*2*N+np.arange(0,2*N)])
        for k in range(N+1):    
            print("Y(N)                                     :",Ys[-1,k*2*N+np.arange(0,2*N)])
    matrix_grad_target=np.reshape(grad_target,(N,N))
    if display: print("grad matrix                              :\n",matrix_grad_target)
    return matrix_grad_target,target,Ys



def run_newton(s,y,TNL,delta_TNL,Lin,parameters=[1.0,0],er0=0.001,display=False):
    """
    Méthode de Newton pour résoudre le système ODE non linéaire
    avec des conditions aux limites sur les deux bords
    """
    if display: print(50*"=","\n"," Méthode de Newton : \n",50*"=","\n")
   
    s0=np.copy(s)

    iterMax = 25
    i,erreur=1,1.0
    while (erreur > er0) and (i <=iterMax) :
        A,b,Ys=solve_syst_direct_gradient(s0,y,TNL,delta_TNL,Lin,parameters=parameters,display=False)
        ds0=np.linalg.solve(A,-b)
        erreur=np.linalg.norm(ds0)/np.linalg.norm(s0)
        s0+=ds0
        if display:
            print('i =  %2i, \t erreur = %12.6f'%(i,erreur))
            print('s = ',s0," t = ",b)
         
        #res=int(input(" entrer 0 pour continuer"))
        res=0
        if res==0:
            pass
        else:
            break
        i+=1
    if  i > iterMax :
        #raise ValueError("pas de convergence")
        print("pas de convergence")

    target,Ys=solve_syst_direct(s0,y,TNL,delta_TNL,Lin,parameters=parameters)
    return s0,target,Ys

     
def calcul_des_pentes_log(y,Ys,display=False):
    """
    Calcul des pentes logarithmiques
    """
    npt=len(y)
    i=npt-3
    p=[]
    p.append(np.log(abs(Ys[i,:])))
    p.append( np.log(abs(Ys[npt-1,:])))
    pente=(p[1]-p[0])/(y[npt-1]-y[int(i)])
    if display:    
        print("Evaluation des pentes log en ymax        :", np.round(pente,2))

def q(f,x,indy,N):
    """
    Calcul de la vitesse sur la paroi
    """
    return -1+np.sum(np.arange(1,N+1)*np.cos(np.arange(1,N+1)*x)*f[indy,np.arange(0,2*N,2)])


"""
    **************************************************************
    PROGRAMME PRINCIPAL
    y : variable espace
    Y : la fonction : Y_0=f_0, Y_1=f'_0, Y_2=f_1, Y_3=f'_1, etc ..
    **************************************************************
"""
 
 
def solution_numerique(N,npt,real_params,type_CI,config,crit_conv,x,sinit,display=False,msg=False): 
    """
    solution numérique en résolvant un système d'ODE non linéaires fortement couplées 
    msg     : pour afficher les messages
    display : pour afficher les courbes
    """ 
    [test_ode,test_dir,test_grad,test_newton,test_grad_num]=config
    [k0,ymax,epsNL,erreur,cmesh,eps]=real_params
    
    Nt=N*(N+1)
    disp_opt=False                   # on pourrait mettre msg ... 
    if msg: print("critère de convergence = ", crit_conv)  
        # paramètre du problème
    # y= np.logspace(-7,np.log10(ymax),npt)     # pas utile
    y=np.linspace(0,ymax,npt)

     
    # initialisation, calcul des coefficients du système ODE et de son gradient
    TNL,Lin,s=calcul_system_coeff(N,display=disp_opt,eps=epsNL)
    delta_TNL=derivative_sys(TNL,s,display=disp_opt)

    if type_CI==1: 
        s0=np.sin(2*np.pi*np.random.rand(N))            # random initial conditions
    elif type_CI==0:
        s0=np.zeros(N)
        s0[0]=k0
    else:
        s0=sinit
     
    #s0  =np.array([-0.95307878 , 0.0718775,  -0.01646209])
    #s0  =np.array([-0.85467194,  0.04825353, -0.00916035])
    s_init=np.copy(s0)
    u=np.zeros(len(x))
     
    if msg: print("IC pour s = ",s0)

    if test_ode: 
        titre("test du système direct pour une valeur de y")
        Y=np.random.rand(2*N)
        print("Valeur initiale de Y = ",Y)
        Ys=sys_direct(Y,y,TNL,Lin,display=True)
        print(" dY / dy = ",Ys)
        
        titre("test du système direct + son gradient pour une valeur de y")
        # différents types de conditions initiales  : (à décommenter)
        #Y=np.random.rand(2*N*(N+1))
        #Y0=np.arange(1,2*Nt+1)
        Y0=np.zeros(2*Nt)
        Y0[np.arange(1,2*N,2)]=k0
        Y0[np.arange(0,2*N,2)]=s0
        ind= 2*(np.arange(0,N)*(N+1)+N)
        Y0[ind]=1
        Y=Y0
        print("Valeur initiale de Y = ",Y)
        Ys=grad_sys(Y,y,TNL,delta_TNL,Lin,display=True)
        print(" dY / dy = ",Ys)
        
    if test_dir:
        titre("test du système direct pour y in [0 ymax]")
        target,Ys=solve_syst_direct(s0,y,TNL,delta_TNL,Lin,[k0,crit_conv],display=True)  
        plot_solution(y,abs(Ys),"Solution directe",np.array([0.,1])*ymax,scale="linear",display=True)
        calcul_des_pentes_log(y,Ys) 
        

    if test_grad: 
        titre("test du système direct + gradient pour y in [0 ymax]")
        grad_target_g,target_g,Ys_g=solve_syst_direct_gradient(s0,y,TNL,delta_TNL,Lin,[k0,crit_conv],display=True) 
        # choix des sorties pour l'affichage (décommenter)
        #ind=np.arange(2*N,2*N*(N+1),2)   # uniquement les fonctions
        #ind=np.arange(0,2*N)             # la solution du système direct               
        ind=np.arange(2*N+1,2*Nt,2)       # uniquement la solution du système avec les gradients
        ind=[5,11]                        # choix manuel
        print("index pour les sorties graphiques        : ",ind)
        plot_solution(y,abs(Ys_g[:,ind]),"Solution directe + gradient",np.array([0.,1])*ymax,scale="semilogy",display=False)
        calcul_des_pentes_log(y,Ys_g[:,ind]) 
        Ys=Ys_g
        
        

    if test_newton:
        s,target,Ys=run_newton(s0,y,TNL,delta_TNL,Lin,[k0,crit_conv],er0=erreur,display=False)
        #plot_solution(y,Ys[:,np.arange(1,2*N+1,2)],"Solution directe",np.array([0.,1])*ymax,display=True)
        #plot_solution(y,np.abs(Ys[:,:]),"Solution directe (dans test_newton)",np.array([0.,1])*ymax,scale="semilogy",display=display)
        if msg:        
            print("delta s / |s| (%)= ",(s-s_init)/np.linalg.norm(s_init)*100)
            #print("delta s / s (%)= ",(s-s_init)/s_init*100)
            print(" s  = ",s)
            print(" s0 = ",s0)
            print(Ys.shape,np.arange(1,2*N+1,2))
        y=np.linspace(0,cmesh*ymax,npt)
        target,Ys=solve_syst_direct(s,y,TNL,delta_TNL,Lin,[k0,crit_conv],display=msg)  
        plot_solution(y,abs(Ys),"Solution directe convergée (test_newton)",np.array([0.,cmesh])*ymax,scale="semilogy",display=display)
        calcul_des_pentes_log(y,Ys) 

        u=np.zeros(len(x))
        for i,xloc in enumerate(x):
            u[i]=q(Ys,xloc,0,N)
        plot_u(x,u,Titre=r"$U_{num}$ (dans test_newton)",lim=np.array([x[0],x[-1]]),scale="linear",display=display)     
        
                
        #print(u.shape,x.shape)
        #print(Ys[-5:-1,np.arange(1,2*N+1,2)])
        #print(Ys[0:10,2*N-1],Ys[-10:-1,2*N-1])
        if msg : print("critère pour le Mach critique  0 =",q(Ys,0,0,N))
            
    if test_grad_num:
        Ys=[]
        target,Ys=solve_syst_direct(s0,y,TNL,delta_TNL,Lin,[k0,crit_conv])
       
        ds=eps*np.linalg.norm(s0)
        for i in range(N): 
            s1=np.copy(s0); s1[i]+=ds
            target_new,Ys=solve_syst_direct(s1,y,TNL,delta_TNL,Lin,[k0,crit_conv])
            dt=target_new-target   
            print("i= %i"%(i), " dt/ds = ", dt/ds," dt = ",dt)
            
    return u,Ys,y