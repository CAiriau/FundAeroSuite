"""
 implémentation des profils NACA ou autres ici
"""

import numpy as np
import matplotlib.pyplot as plt

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 
 


def naca4(number, n, finite_TE = False, half_cosine_spacing = False):
    """
    Returns 2*n+1 points in [0 1] for the given 4 digit NACA number string
    """

    m = float(number[0])/100.0
    p = float(number[1])/10.0
    t = float(number[2:])/100.0

    a0,a1,a2,a3 = +0.2969,-0.1260,-0.3516, +0.2843
     

    if finite_TE:
        a4 = -0.1015 # For finite thick TE
    else:
        a4 = -0.1036 # For zero thick TE

    if half_cosine_spacing:
        beta = np.linspace(0.0,np.pi,n+1)
        #x = [(0.5*(1.0-np.cos(xx))) for xx in beta]  # Half cosine based spacing
        x = 0.5*(1.0-np.cos(beta))  # Half cosine based spacing
    else:
        x = np.linspace(0.0,1.0,n+1)

    yt = 5*t*(a0*np.sqrt(x)+a1*x+a2*np.power(x,2)+a3*np.power(x,3)+a4*np.power(x,4)) 

    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]

    if p == 0:
        xu,yu = x,yt
        xl,yl = x,-yt
        xc = xc1 + xc2
        print('len  =',len(xc1),len(xc2),len(xc))
        zc = [0]*len(xc)
        theta=[0]*len(xc)
    else:
        yc1 = [m/np.power(p,2)*xx*(2*p-xx) for xx in xc1]
        yc2 = [m/np.power(1-p,2)*(1-2*p+xx)*(1-xx) for xx in xc2]
        zc = yc1 + yc2

        dyc1_dx = [m/np.power(p,2)*(2*p-2*xx) for xx in xc1]
        dyc2_dx = [m/np.power(1-p,2)*(2*p-2*xx) for xx in xc2]
        dyc_dx = dyc1_dx + dyc2_dx

        theta = [np.arctan(xx) for xx in dyc_dx]

        xu = [xx - yy * np.sin(zz) for xx,yy,zz in zip(x,yt,theta)]
        yu = [xx + yy * np.cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

        xl = [xx + yy * np.sin(zz) for xx,yy,zz in zip(x,yt,theta)]
        yl = [xx - yy * np.cos(zz) for xx,yy,zz in zip(zc,yt,theta)]

    X = np.append(xu[::-1] , xl[1:])
    Z = np.append(yu[::-1] , yl[1:])

    return X,Z,theta

def naca5(number, n, finite_TE = False, half_cosine_spacing = False):
    """
    Returns 2*n+1 points in [0 1] for the given 5 digit NACA number string
    """

    

def naca(number, n, finite_TE = False, half_cosine_spacing = False):
    if len(number)==4:
        return naca4(number, n, finite_TE, half_cosine_spacing)
    elif len(number)==5:
        return naca5(number, n, finite_TE, half_cosine_spacing)
    else:
        raise Exception

def EffetVoletBec(Volet=True,theta=0,longueur=0.0,angle=0):
    """
    Effet des volets et bec par la théorie linéarisée
    On entre soit  theta (°), soit longueur/corde
    angle = angle de braquage en degrés   
    """
  
    if isinstance(longueur,float):
        if (theta==0) and (longueur==0):
                raise ValueError('un des deux paramètres doit être non nul dans EffetVolet')

        if Volet:
            print("Braquage d'un volet : ")
            if theta==0:
                theta=np.rad2deg(np.arccos(1-2*longueur))
            elif longueur==0:
                longueur= (1-np.cos(np.deg2rad(theta)))/2


            print("longueur volet/corde              : %f"%(longueur))
            print("angle theta                       : %f"%(theta))
            
            # theta est en degré.

            # variations des coefficients A[]
            #A[0:2]       = []
            t            = np.deg2rad(theta)
            beta         = np.deg2rad(angle)
            DeltaCL      = 2*beta*(t+np.sin(t)) 
            DeltaCmF     = beta*np.sin(t)*(1+np.cos(t))/2
            DeltaAlpha_0 = -beta*(t+np.sin(t))/np.pi
            DeltaAlpha_a = -beta*t/np.pi

        else :
            print("Braquage d'un bec : ")
            # un bec
            if theta==0:
                theta=np.rad2deg(np.arccos(2*longueur-1))
            elif longueur==0:
                longueur= (1+np.cos(np.deg2rad(theta)))/2

            print("longueur bec/corde                : %f"%(longueur))
            print("angle theta                       : %f"%(theta))

            # variations des coefficients A[]
            #A[0:2]       = []
            t            = np.deg2rad(theta)
            beta         = np.deg2rad(angle)
            DeltaCL      = -2*beta*(np.pi-t-np.sin(t)) 
            DeltaCmF     = beta*np.sin(t)*(1+np.cos(t))/2
            DeltaAlpha_0 = beta*(np.pi-t-np.sin(t))/np.pi
            DeltaAlpha_a = beta*(1-t/np.pi)


        if angle != 0 :

            print("angle beta                        : %f °"%(angle))
            print("angle beta                        : %f rad"%(beta))
            print("Delta CL                          : %f"%(DeltaCL))
            print("Delta CmF                         : %f"%(DeltaCmF))
            print("Delta alpha_0                     : %f rad, %f°"%(DeltaAlpha_0,np.rad2deg(DeltaAlpha_0)))
            print("Delta alpha_a                     : %f rad, %f°"%(DeltaAlpha_a,np.rad2deg(DeltaAlpha_a)))
    else:
        
        if Volet:
            print("Braquage d'un volet : ")
            theta=np.rad2deg(np.arccos(1-2*longueur))
            # theta est en degré.
            # variations des coefficients A[]
            #A[0:2]       = []
            t            = np.deg2rad(theta)
            beta         = np.deg2rad(angle)
            DeltaCL      = 2*beta*(t+np.sin(t)) 
            DeltaCmF     = beta*np.sin(t)*(1+np.cos(t))/2
            DeltaAlpha_0 = -beta*(t+np.sin(t))/np.pi
            DeltaAlpha_a = -beta*t/np.pi
        else :
            print("Braquage d'un bec : ")
            # un bec
            theta=np.rad2deg(np.arccos(2*longueur-1))
            # variations des coefficients A[]
            #A[0:2]       = []
            t            = np.deg2rad(theta)
            beta         = np.deg2rad(angle)
            DeltaCL      = -2*beta*(np.pi-t-np.sin(t)) 
            DeltaCmF     = beta*np.sin(t)*(1+np.cos(t))/2
            DeltaAlpha_0 = beta*(np.pi-t-np.sin(t))/np.pi
            DeltaAlpha_a = beta*(1-t/np.pi)

    return DeltaCL,DeltaCmF,DeltaAlpha_0,DeltaAlpha_a

 
 
