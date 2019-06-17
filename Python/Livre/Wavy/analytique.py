# -*- coding: utf-8 -*-
"""
Created on Tue May 22 11:40:21 2018

@author: cairiau
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def f_values(y,k):
    """
    calcul de la fonction vectorielle f
    """
   
    f=np.zeros(6)
    a=np.zeros(3)
    a[0] = (1+(11/256)*k**2+1861*k**4/(3*256**2))*k;
    a[1] = -(5/256+2765*k**2/(9*256**2)+(1/64+33*k**2/(64*256))*y)*k**3;
    a[2] = -(-65/(576*256)+7*y/(72*256)+9*y**2/(32*256)+y**3/(8*256))*k**5;
    f[0] = a[0]*np.exp(-y)+a[1]*np.exp(-3*y)+a[2]*np.exp(-5*y)

    b=np.zeros(3)
    b[0] = (1/16+419*k**2/(72*256)+234215*k**4/(432*256**2)+(1/8+(11/1024)*k**2+4085*k**4/(24*256**2))*y)*k**2;
    b[1] = -k**4*(125/(36*256)+45895*k**2/(108*256**2)+(5/256+815*k**2/(288*512))*y+y**2*(1/128+11*k**2*(1/8192)));
    b[2] = -k**6*(-12245/(144*256**2)-6245*y/(72*256**2)+79*y**2/(3*256**2)+47*y**3/(192*256)+y**4/(12*256));
    f[1] = b[0]*np.exp(-2*y)+b[1]*np.exp(-4*y)+b[2]*np.exp(-6*y)
 

    c=np.zeros(2)
    c[0] = k**3*(1/72+23603*k**2/(27*256**2)+y*(1/24+259*k**2/(72*256))+y**2*(1/32+33*k**2*(1/8192)));
    c[1] = -k**5*(1765/(3*256**2)+155*y/(32*256)+117*y**2/(32*256)+(1/256)*y**3);
    f[2] = c[0]*np.exp(-3*y)+c[1]*np.exp(-5*y)

    d=np.zeros(2)
    d[0] = k**4*(7/1536+390547*k**2/(720*256**2)+y*(7/384+32947*k**2/(36*256**2))+y**2*(3/128+617*k**2/(288*256))+y**3*(1/96+11*k**2/(24*256)));
    d[1] = -k**6*(35251/(90*256**2)+2465*y/(576*256)+455*y**2/(96*256)+119*y**3/(48*256)+25*y**4/(48*256));
    f[3] = d[0]*np.exp(-4*y)+d[1]*np.exp(-6*y)

    f[4] = k**5*np.exp(-5*y)*(7/3840+7*y*(1/768)+25*y**2*(1/1536)+5*y**3*(1/384)+25*y**4/(24*256));
    f[5] = k**6*np.exp(-6*y)*(91/(54*2048)+91*y/(72*256)+23*y**2*(1/2048)+13*y**3*(1/1024)+15*y**4*(1/2048)+9*y**5*(1/5120));
     
    return f

def q(x,y,k):
    """
    Calcul de la vitesse sur la paroi
    """
    return -1+np.sum(np.arange(1,7)*np.cos(np.arange(1,7)*x)*f_values(y,k))


def plot_solution(y,Ysol,Titre="Solution",lim=np.array([0,5]),scale="linear",display=False):
    """
    dessin de la solution
    """
    if display:
        t=Ysol.shape
        print("dimension des données                    : %i x %i "%(t[0],t[1]))
        N=t[1] 
        plt.figure(1,figsize=(12,10))
        plt.title(Titre)
        plt.xlabel(r'$y$',fontsize=16)
        plt.grid()
        #ax.axis([-1,1,0,1])
        plt.xlim(lim[0],lim[1])
        
        if scale=="linear":
            for k in range(N): 
                plt.plot(y,Ysol[:,k],'-', linewidth=2,label="k = %i"%(k))
        elif scale=="semilogy":
            for k in range(N): 
                plt.semilogy(y,Ysol[:,k],'-', linewidth=2,label="k = %i"%(k))
        plt.legend(loc="best")
         
def plot_u(y,Ysol,Titre="Solution",lim=np.array([0,5]),scale="linear",display=False):
    """
    dessin de la solution
    """
    if display:
        plt.figure(2,figsize=(12,10))
        plt.title(Titre)
        plt.xlabel(r'$x$',fontsize=16)
        plt.grid()
        #ax.axis([-1,1,0,1])
        plt.xlim(lim[0],lim[1])
        
        if scale=="linear":
            plt.plot(y,Ysol,'-', linewidth=2,label="u")
        elif scale=="semilogy":
            plt.semilogy(y,Ysol,'-', linewidth=2,label="u")
        
        #plt.legend(loc="lower left")
         

def solution_analytique(N,npt,k0,x,y,display="False"):
    """
    calcul de la solution analytique
    """
    f=np.zeros((npt,N))
    for  i in range(npt):
        f[i,:]=f_values(y[i],k0)
    plot_solution(y,f,Titre=r"Solution analytique $f$ ",lim=np.array([0,y[-1]]),scale="semilogy",display=display)
    #print("f(0) = " ,f[0,:])

    u=np.zeros(len(x))
    for i,xloc in enumerate(x):
        u[i]=q(xloc,y[0],k0)
    plot_u(x,u,Titre=r"$|{\bf u}|_{wall}$ (analytique)",lim=np.array([ x[0], x[-1] ]),scale="linear",display=display)       
        
    return u,f 
 

#plot({log10(abs(f[1])), log10(abs(f[2])), log10(abs(f[3])), log10(abs(f[4])), log10(abs(f[5])), log10(abs(f[6]))}, y = 0 .. 10);
# Conditionx aux limites
# 
# for i to 6 do f_b[i, 1] = evalf(subs(y = 0, f[i])); f_b[i, 2] = evalf(subs(y = 0, diff(f[i], y))); printf("i= %d, f(0) = %e  f0'=%e, \n", i, f_b[i, 1], f_b[i, 2]) end do;
# for i to 6 do erreur[i] = simplify(f[i]/(f_b[i, 1]*exp(-i*y))-1) end do;
# plot({log10(abs(erreur[1])), log10(abs(erreur[2])), log10(abs(erreur[3])), log10(abs(erreur[4])), log10(abs(erreur[5])), log10(abs(erreur[6]))}, y = 0 .. 1);
# eq[1] = simplify(diff(f[1], y, y)-f[1]+f[1]*f[2]+3*f[2]*f[3]+6*f[3]*f[4]+10*f[4]*f[5]+15*f[5]*f[6]);
# eq[2] = simplify(diff(f[2], y, y)-4*f[2]+(1/2)*f[1]*f[1]+3*f[1]*f[3]+8*f[2]*f[4]+15*f[3]*f[5]+24*f[4]*f[6]);
# eq[3] = simplify(diff(f[3], y, y)-9*f[3]+3*f[1]*f[2]+6*f[1]*f[4]+15*f[2]*f[5]+27*f[3]*f[6]);
# eq[4] = simplify(diff(f[4], y, y)-16*f[4]+6*f[1]*f[3]+4*f[2]**2+10*f[1]*f[5]+24*f[2]*f[6]);
# eq[5] = simplify(diff(f[5], y, y)-25*f[5]+10*f[1]*f[4]+15*f[2]*f[3]+15*f[1]*f[6]);
# eq[6] = simplify(diff(f[6], y, y)-36*f[4]+15*f[1]*f[5]+(27*(1/3))*f[3]**2+24*f[2]*f[4]);
# plot({log10(abs(eq[1])), log10(abs(eq[2])), log10(abs(eq[3])), log10(abs(eq[4])), log10(abs(eq[5])), log10(abs(eq[6]))}, y = 0 .. 10);
# for l to 6 do printf(" eq en y= 0, log10(|eq[ %d ]|) = %e \n", l, subs(y = 0, log10(abs(eq[l])))) end do;
# u = -1+sum(j*f[j]*cos(j*x), j = 1 .. 6);
 
 
