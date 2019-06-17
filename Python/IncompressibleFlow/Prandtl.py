#!/bin/python
"""
  Correction des exercices du chapitre 5
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
 

def Vitesse_Tourbillon(X,x0,ya,yb,Gamma):
	"""
	Vitesse induite par un tourbillon en fer à cheval
	X =(x,y)
	A=(x0,y_A)
	B=(x0,y_B)
	Gamma : circulation
	"""
	x,y=X[0],X[1]
	g=Gamma/(4*np.pi)

	if x==x0 and y==(ya+yb)/2 :
		return -Gamma/np.pi/(yb-ya)
	elif x==x0:
		return -g*(1/(y-ya)-1/(y-yb))
	else:
		db=np.sqrt((x-x0)**2+(y-yb)**2)
		da=np.sqrt((x-x0)**2+(y-ya)**2)
		
		w1 = g / (x-x0)*( (y-yb)/db -(y-ya)/da)
		w2 = g / (y-yb)*(1 +(x-x0)/db)
		w3 = -g / (y-ya)*(1 +(x-x0)/da)

		w=w1+w2+w3
		w4=-g / (y-ya)*(1 +da/(x-x0)) + g / (y-yb)*(1 + db /(x-x0))
		# w4 est bien égale à w

	return w

def Prandtl_E_OGE(Nc,lamb,b,l0,k,y,theta,alpha):
        """
        Calcul des coefficients de la série de Prandtl pour l'aile elliptique
        sans effet de sol
        """
        l=l0*np.sin(theta)          # évolution de la corde en envergure
        My=len(theta)
        mat = np.zeros((My-2,Nc), dtype=float)  
        RHS = np.zeros((My-2), dtype=float) 
        u=np.arange(1,My-1)
        for n in range(Nc):
            mat[:,n]=4*b*np.sin((n+1)*theta[u])+k*l[u]*(n+1)*np.sin((n+1)*theta[u])/np.sin(theta[u]);
        RHS=k*alpha*l[u]
        A,res,rank,s = np.linalg.lstsq(mat, RHS)
       
        if Nc==My-2:
            print(np.linalg.solve(mat, RHS))
        print('A_n = ',A)
        return A

def integrand(x,theta,n,h_sur_b):
    """
    intégrand pour le calcul IGE
    """
    tmp=np.cos(x)-np.cos(theta)
    return (tmp*np.cos(n*x))/(tmp**2+(2*h_sur_b)**2)
    
def solve_integralO2(theta,n,h):
    """
    Intégrale d'ordre 2, méthode des trapèzes
    """
    s=np.zeros(len(theta), dtype=float)
    m=720*8
    phi=np.linspace(0,np.pi,m)
    dphi=phi[1]-phi[0]
    for i in range(len(theta)):
        val=integrand(phi,theta[i],n,h)
        s[i]=dphi*(np.sum(val[1:m-2])+(val[0]+val[m-1])/2)
    return s

def solve_integralO4(theta,n,h_sur_b):
    """
    Intégrale avec des quadratures
    """
    s=np.zeros(len(theta), dtype=float)
    for i in range(len(theta)):
        s[i],res=integrate.quad(integrand,0,np.pi,args=(theta[i],n,h_sur_b,))
    return s

def Prandtl_E_IGE(Nc,lamb,b,l0,k,h,y,theta,alpha):
        """
        Calcul des coefficients de la série de Prandtl pour l'aile elliptique
        avec l'effet de sol
        h : h / b
        """
        #l=l0*np.sin(theta)          # évolution de la corde en envergure
        My=len(theta)
        mat = np.zeros((My-2,Nc), dtype=float)  
        RHS = np.zeros((My-2), dtype=float) 
        A = np.zeros((Nc), dtype=float) 
        G = np.zeros((My-2), dtype=float) 
        u=np.arange(1,My-1)
        for n in range(Nc):
            coef=2*k*(n+1)/(np.pi*lamb)
            G=solve_integralO4(theta[u],n+1,h)
            #G1=solve_integralO2(theta[u],n+1,h)
            #t1,t2=np.linalg.norm(G),np.linalg.norm(G1)
            #print('erreur = ',t2/t1-1)
            mat[:,n]=(1+coef)*np.sin((n+1)*theta[u])-coef/np.pi*G*np.sin(theta[u]) 

        RHS=2*k*alpha/(np.pi*lamb)*np.sin(theta[u]) 
        
        A,res,rank,s = np.linalg.lstsq(mat, RHS)
        
        if Nc==My-2:
            print(np.linalg.solve(mat, RHS))
        #print('A_1 = %3.15f'%(A[0]))
        return A


 
