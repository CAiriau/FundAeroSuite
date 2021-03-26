#!/bin/py
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Module : incompressibleFlow.Airfoils ,   to get NACA airfoil geometry
    ..
    
    * naca4
    * naca5
    * naca
    * slats_flaps_effects
"""

import numpy as np
import matplotlib.pyplot as plt


def naca4(number, n, finite_TE=False, half_cosine_spacing=False):
    """
    get NACA 4 series geometry
    
    Args:
        number (str) : name of the NACA 4 airfoil
        n (int) : number of point to define the x-coordinate
        finite_TE (bool) : True if the trailing thickness is not null
        half_cosine_spacing (bool) : True if the x-vector is defined from arccos 
    
    Returns:
        list of float: X,Z,theta

        * X : 2 x n + 1  points in [0 1] for the given 4 digit NACA number string
        * Z : upper and lower wall coordinates
        * theta : slope of the camberline in radians
    """
    m = float(number[0]) / 100.0
    p = float(number[1]) / 10.0
    t = float(number[2:]) / 100.0
    a0, a1, a2, a3 = +0.2969, -0.1260, -0.3516, +0.2843

    if finite_TE:
        a4 = -0.1015  # For finite thick TE
    else:
        a4 = -0.1036  # For zero thick TE

    if half_cosine_spacing:
        beta = np.linspace(0.0, np.pi, n + 1)
        # x = [(0.5*(1.0-np.cos(xx))) for xx in beta]  # Half cosine based spacing
        x = 0.5 * (1.0 - np.cos(beta))  # Half cosine based spacing
    else:
        x = np.linspace(0.0, 1.0, n + 1)

    yt = 5 * t * (a0 * np.sqrt(x) + a1 * x + a2 * np.power(x, 2) + a3 * np.power(x, 3) + a4 * np.power(x, 4))
    xc1 = [xx for xx in x if xx <= p]
    xc2 = [xx for xx in x if xx > p]

    if p == 0:
        # symmetrical airfoil
        xu, yu = x, yt
        xl, yl = x, -yt
        xc = xc1 + xc2
        print('NACA4 : len  =', len(xc1), len(xc2), len(xc))
        zc = [0] * len(xc)
        theta = [0] * len(xc)
    else:
        # airfoil with a camberline
        yc1 = [m / np.power(p, 2) * xx * (2 * p - xx) for xx in xc1]
        yc2 = [m / np.power(1 - p, 2) * (1 - 2 * p + xx) * (1 - xx) for xx in xc2]
        zc = yc1 + yc2
        dyc1_dx = [m / np.power(p, 2) * (2 * p - 2 * xx) for xx in xc1]
        dyc2_dx = [m / np.power(1 - p, 2) * (2 * p - 2 * xx) for xx in xc2]
        dyc_dx = dyc1_dx + dyc2_dx
        theta = [np.arctan(xx) for xx in dyc_dx]
        xu = [xx - yy * np.sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yu = [xx + yy * np.cos(zz) for xx, yy, zz in zip(zc, yt, theta)]
        xl = [xx + yy * np.sin(zz) for xx, yy, zz in zip(x, yt, theta)]
        yl = [xx - yy * np.cos(zz) for xx, yy, zz in zip(zc, yt, theta)]

    X = np.append(xu[::-1], xl[1:])
    Z = np.append(yu[::-1], yl[1:])

    return X, Z, theta


def naca5(number, n, finite_TE=False, half_cosine_spacing=False):
    """
    get NACA 5 series geometry 
    
    .. warning:
        Not implemented yet in Python
    
    Args:
        number (int) : name of the NACA 4 airfoil
        n (int) : number of point to define the x-coordinate
        finite_TE (bool) : True if the trailing thickness is not null
        half_cosine_spacing (bool) : True if the x-vector is defined from arccos 
    
    Returns:
        list of float: X,Z,theta 

        * X : 2 x n + 1  points in [0 1] for the given 4 digit NACA number string
        * Z : upper and lower wall coordinates
        * theta : slope of the camberline in radians
    """

    raise NameError("NAVA 5 airfoil not implemented yet in Python, found in Fortran")


def naca(number, n, finite_TE=False, half_cosine_spacing=False):
    """
    get NACA 4 series geometry 
    
    .. warning:
        Not implemented yet in Python
    
    Args:
        number (str) : name of the NACA 4 airfoil
        n (int) : number of point to define the x-coordinate
        finite_TE (bool) : True if the trailing thickness is not null
        half_cosine_spacing (bool) : True if the x-vector is defined from arccos 
        
    Returns:
        list of float: X,Z,theta

        * X : 2 x n + 1  points in [0 1] for the given 4 digit NACA number string
        * Z : upper and lower wall coordinates
        * theta : slope of the camberline in radians
    """

    if len(number) == 4:
        return naca4(number, n, finite_TE, half_cosine_spacing)
    elif len(number) == 5:
        return naca5(number, n, finite_TE, half_cosine_spacing)
    else:
        raise NameError("NACA : bad choice of airfoil name")

def flap_effect(theta, angle):
    """
    Effect of flap with a given angle, length given by theta
    """
    t = np.deg2rad(theta)
    beta = np.deg2rad(angle)
    DeltaCL = 2 * beta * (t + np.sin(t))
    DeltaCmF = beta * np.sin(t) * (1 + np.cos(t)) / 2
    DeltaAlpha_0 = -beta * (t + np.sin(t)) / np.pi
    DeltaAlpha_a = -beta * t / np.pi
    return DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a

def slat_effect(theta, angle):
    """
    Effect of slat with a given angle, length given by theta
    """
    t = np.deg2rad(theta)
    beta = np.deg2rad(angle)
    DeltaCL = -2 * beta * (np.pi - t - np.sin(t))
    DeltaCmF = beta * np.sin(t) * (1 + np.cos(t)) / 2
    DeltaAlpha_0 = beta * (np.pi - t - np.sin(t)) / np.pi
    DeltaAlpha_a = beta * (1 - t / np.pi)
    return DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a



def slats_flaps_effects(Flaps=True, theta=0, length=0.0, angle=0):
    """
    Effect of flaps and slat with the linearized airfoil theory
    
    Args:
        Flaps (bool): if True use of flaps depending on theta value, else slats if lenght non equal to zero 
        theta (float) : angle in degree defining the size of the flaps or slats, along the chord
        length (float) : size of the flaps or slats with respect to chord (in percents)
        angle (float) : angle in degrees of inclination of the flaps or slats
    
    Returns:  
        float : DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a
        
        * DeltaCL : CL variation
        * DeltaCmF : CmF variation (pitching coefficient at the aerodynamic center)
        * DeltaAlpha_0 : variation of the zero-lift angle in radians
        * DeltaAlpha_a : variation of the adaptation angle in radians
    """
    if isinstance(length, float):
        if (theta == 0) and (length == 0):
            raise ValueError('theta or length parameter must be non null in slats_flaps_effects function')

        if Flaps:
            print("use of a flap : ")
            if theta == 0:
                theta = np.rad2deg(np.arccos(1 - 2 * length))
            elif length == 0:
                length = (1 - np.cos(np.deg2rad(theta))) / 2
            print("length flap/chord                 : %f" % (length))
            print("theta angle                       : %f °" % (theta))
            DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a = flap_effect(theta, angle)


        else:
            print("use of slat : ")
            if theta == 0:
                theta = np.rad2deg(np.arccos(2 * length - 1))
            elif length == 0:
                length = (1 + np.cos(np.deg2rad(theta))) / 2
            print("length slat/chord                : %f" % (length))
            print("angle theta                       : %f" % (theta))
            DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a = slat_effect(theta, angle)

        if angle != 0:
            print("angle beta                        : %f °" % (angle))
            print("angle beta                        : %f rad" % (np.deg2rad(angle)))
            print("Delta CL                          : %f" % (DeltaCL))
            print("Delta CmF                         : %f" % (DeltaCmF))
            print("Delta alpha_0                     : %f rad, %f°" % (DeltaAlpha_0, np.rad2deg(DeltaAlpha_0)))
            print("Delta alpha_a                     : %f rad, %f°" % (DeltaAlpha_a, np.rad2deg(DeltaAlpha_a)))

    else:

        if Flaps:
            print("use of flap : ")
            theta = np.rad2deg(np.arccos(1 - 2 * length))
            # theta est en degree.
            DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a = flap_effect(theta, angle)

        else:
            print("use of slat : ")
            theta = np.rad2deg(np.arccos(2 * length - 1))
            DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a = slat_effect(theta, angle)


    return DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a