#!/bin/py
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Module : CompressibleFlow.fonctions:  a collection of function to study compressible flow:
    ..
    
    * isentropic state
    * transonic flow
    * supersonic 1D 
    * supersonic 2D (oblique shock and Prandtl-Meyer expansion)
"""

import numpy as np
import matplotlib.pyplot as plt

__R_ideal_gas__ = np.float64(8.314462)
__molar_mass_Air__ = 28.976e-3
r_Air = __R_ideal_gas__ / __molar_mass_Air__
# print("accurate Air r= R / M  = ", r_Air)     # 286.9430563224738
# constants
gamma = np.float64(1.4)  # :math:`C_p/C_v`
r = np.float64(287)  # r ideal gas (can replace accurate value r_Air)


def sound_velocity(T, gamma=gamma):
    """ 
    Sound velocity of a gas
    
    Args:
        * T (real) : temperature [K]
        * gamma (real) : :math:`C_p/C_v`
        
    Returns :
        real : sound velocity [m/s]
    """
    return np.sqrt(gamma * r * T)


def Pressure_coefficient(ratio, M, gamma=gamma):
    """
    Pressure coefficient function of p/p0 and M, 
    
    Args:
        * ratio (real) : p/p_0, p0 :static reference pressure
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : Kp, pressure coefficient
    """
    return 2 / gamma / M ** 2 * (ratio - 1)


def p_pi(Mach, gamma=gamma):
    """
    ratio pressure/isentropic pressure function of Mach number

    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : p/p_i
    """
    return (T_Ti(Mach, gamma=gamma)) ** (gamma / (gamma - 1))


def rho_rhoi(Mach, gamma=gamma):
    """
    ratio rho/isentropic rho function of Mach number
    
    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : rho/rho_i
    """
    return (T_Ti(Mach, gamma=gamma)) ** (1. / (gamma - 1))


def T_Ti(Mach, gamma=gamma):
    """
    ratio Temperature/isentropic Temperature function of Mach number
    
    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : T/T_i
    """
    return 1. / (1 + (gamma - 1) / 2 * Mach ** 2)


def pt_p(Mach, gamma=gamma):
    """
    ratio Total pressure / static pressure function of Mach number with :math:`p_t = p + \\rho u^2/2`
    
    Args:
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns: 
        real : p_t/p
    """
    return 1 + gamma / 2 * Mach ** 2


def inverse_p_pi(ratio, gamma=gamma):
    """
    Calcul du nombre de Mach connaissant la valeur de p/pi
    
    Args:
        * ratio (real) : p/p_i, p_i : isentropic pressure
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : M, Mach number
    """
    return np.sqrt(2. / (gamma - 1) * (pow(ratio, (1. - gamma) / gamma) - 1.))


def compressibility_corrections(Kpinc, M):
    """
    3 compressibility corrections: Prandtl-Glauert (PG), Karman-Tsien (KT),
    Laitone (L) as a function of Mach number
    Kpinc = Kp incompressible flow, M : Mach number
    
    Args:
        * Kpinc (real): Kp for incompressible flow
        * M (real): Mach number
    
    Returns:
        list of real: Kp_out, Pressure coefficient for compressible flow, with 3 corrections 
    """
    beta = np.sqrt(1 - M ** 2)
    Kp_out = np.zeros((len(M), 3))
    Kp_PG = Kpinc / beta
    Kp_KT = Kpinc / (beta + M ** 2 / 2 * Kpinc / (1 + beta))
    Kp_L = Kpinc / (beta + M ** 2 * Kpinc / (2 * beta) * (1 + (gamma - 1) / 2 * M ** 2))
    Kp_out[:, 0] = Kp_PG
    Kp_out[:, 1] = Kp_KT
    Kp_out[:, 2] = Kp_L
    return Kp_out


def Kp_critical(Mach, gamma=gamma):
    """ 
    Critical Kp function of the Mach number
    
    Args:
        * Mach (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`

    Returns:
        real : Kp_c, critical pressure coefficient
    """
    return 2 / (gamma * Mach ** 2) * (p_pi(1, gamma=gamma) / p_pi(Mach, gamma=gamma) - 1)


def CL_compressible(Mach, CL):
    """
    Prandtl-Glauert Correction on the C_L
    
    Args:
        * Mach (real) : Mach number
        * CL (real) : lift coefficient in incompressible flow
        
    Returns:
        real : lift coefficient in compressible flow w.r.t. Mach
        
    """
    return CL / np.sqrt(1 - Mach ** 2)


def PG_correction_wing(Kpinc, M, Lambda):
    """
    3 compressibility corrections of a swept wing with Prandtl-Glauert (PG),  as a function of Mach number
    
    Args:
        * Kp_inc (real) : Kp incompressible flow
        * M (real) : Mach number
        * Lambda (real) : swept angle in degrees
        
    Returns:
        real: Kp compressible (with correction)
    """
    cos = np.cos(np.deg2rad(Lambda))
    beta = np.sqrt(1 - (M * cos) ** 2)
    return Kpinc *cos / beta


def omega_super(Mach, gamma=gamma):
    """
    Prandtl-Meyer function en radians for supersonic flows
    
    Args: 
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`
        
    Return :
        real : omega, Prandtl-Meyer angle in radians
    """
    a = (gamma + 1) / (gamma - 1);
    x = Mach ** 2 - 1.0
    return np.sqrt(a) * np.arctan(np.sqrt(x / a)) - np.arctan(np.sqrt(x))


def omega_hyper(Mach, gamma=gamma):
    """
    Prandtl-Meyer function en radians for hypersonic flows
    
    Args: 
        * M (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`
        
    Return :
        real : omega, Prandtl-Meyer angle in radians
    """
    omega_inf = np.pi / 2 * (np.sqrt((gamma + 1) / (gamma - 1)) - 1.0)
    return omega_inf - 2.0 / ((gamma - 1) * Mach)


def invomega(angle, show=False, gamma=gamma):
    """
    Inversse of Prandtl-Meyer function  for superponic flows, solve by a Newton method
    
    Args:
        * angle (real) : angle of the PM function in radians
        * show (boolean) : True , details of Newton method are displayed
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : Mach number
    """
    er0, dmach, mach, iterMax = 1e-6, 1e-3, 2.0, 20
    i, erreur = 1, 1.0

    while (erreur > er0) and (i <= iterMax):
        omega0 = omega_super(mach, gamma=gamma) - angle
        omega1 = omega_super(mach + dmach, gamma=gamma) - angle
        dm = -dmach * omega0 / (omega1 - omega0)
        erreur = abs(dm / mach)
        if show:
            print('i =  %2i, \t error = %12.6f, \t dm = %12.6f, \t M = %12.6f' % (i, erreur, dm, mach))
        mach += dm
        i += 1

    if i > iterMax:
        print("no convergence")
        mach = -1
    return mach


def isentropic_expansion(theta, Mach, gamma=gamma):
    """
    pressure ratio in a isentropic expansion
    
    Args:
        * theta (real) : deviation angle in radians
        * M (real) : Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream p / upstream p
    """
    downstream_omega = omega_super(Mach, gamma=gamma) + np.abs(theta)
    rpi = []
    for omega in downstream_omega:
        rpi.append(p_pi(invomega(omega), gamma=gamma) / p_pi(Mach, gamma=gamma))
    return rpi


# ********   OBLIQUE SHOCK WAVE    *********

def shock_angle(sigma, teta, M0, gamma=gamma):
    """     
    ratio of density equation  (f=0) across a shock, to get the shock angle in radian
    used in Newton method (f = 0)
    
    Args:
        * sigma (real) : shock angle in radians
        * theta (real) : deviation angle in radians
        * M0 (real) : upstream Mach number
        * gamma (real) : :math:`C_p/C_v`
    
    Returns :
        real : f (is = 0 at convergence)
    """
    return np.tan(sigma - teta) / np.tan(sigma) - 2.0 / (gamma + 1.0) / M0 ** 2 / (np.sin(sigma)) ** 2 - (
                gamma - 1.0) / (gamma + 1.0)


def valeur_sigma(M0, theta_d, sigma=0, show=False, gamma=gamma):
    """
    Newton metho to solve the oblique shock angle in supersonic flow 

    Args:
        * M0 (real) : upstream Mach number
        * theta_d (real) : deviation angle in degrees
        * sigma (real) : initial guess of the shock angle in radians
        * show (boolean) : to display details on screen
    
    Returns :
     * sigma (real) : shock angle in radians
    """
    teta = np.deg2rad(theta_d)
    if sigma == 0:
        sigma = np.arcsin(1.0 / M0)
    er0, ds, erreur = 1e-6, np.deg2rad(1e-2), 1.0
    iterMax, i = 20, 1
    if show:
        print("Iterations \t sigma(°) \t dsigma  \t error")

    while (abs(erreur) > er0) and (i <= iterMax):
        f0 = shock_angle(sigma, teta, M0, gamma=gamma)
        f1 = shock_angle(sigma + ds, teta, M0, gamma=gamma)
        dsigma = -ds * f0 / (f1 - f0)
        erreur = dsigma / sigma
        if show:
            print(' %2i \t %12.6f \t %12.6f \t %12.6f' % (i, np.rad2deg(sigma), np.rad2deg(dsigma), erreur))
        sigma += dsigma
        i += 1
    if i > iterMax:
        print('no convergence, deviation angle greater than the one leading to detached shock')
        sigma = -1
    elif sigma < 0:
        print('No convergence, negative angle  : increase the upstream Mach number or decrease the deviation angle')
        sigma = -1
    elif np.rad2deg(sigma) > 90:
        print('No convergence , angle > 90 : increase the upstream Mach number or decrease the deviation angle')
        sigma = -1
    elif show:
        print('convergence on sigma reached')
    return np.rad2deg(sigma)


def P2_P1(mach, gamma=gamma):
    """ 
    pressure ratio across the shock wave, oblique or normal shock wave
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream p / upstream p
    """
    return 2.0 * gamma / (gamma + 1.0) * mach ** 2 - (gamma - 1.0) / (gamma + 1.0)


def Inverse_P2_P1(ratio, gamma=gamma):
    """     
    inverse function, get the Mach number from a pressure ratio across a shock
    
    Args:
        * ratio (real) : downstream p / upstream p 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real :  normal Mach number 
    """
    return np.sqrt(1.0 / (2.0 * gamma) * ((gamma + 1.0) * ratio + (gamma - 1.0)))


def rho2_rho1(mach, gamma=gamma):
    """    
    density ratio across the shock wave, oblique or normal shock wave
    
    Args:
        * mach (real) : Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream rho / upstream rho   
    """
    return 1.0 / (2.0 / ((gamma + 1.0) * mach ** 2) + (gamma - 1.0) / (gamma + 1.0))


def pi2_pi1(mach, gamma=gamma):
    """
    isentropic pressure ratio across the shock wave, oblique or normal shock wave
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream p_i / upstream p_i
    """
    return pow(P2_P1(mach, gamma=gamma), -1 / (gamma - 1)) * pow(rho2_rho1(mach, gamma=gamma), gamma / (gamma - 1))


def downstream_Mach(mach, gamma=gamma):
    """       
    downstream normal Mach across a shock
    
    Args:
        * mach (real) : normal Mach number 
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : downstream normal Mach number
    """
    return np.sqrt((1.0 + 0.5 * (gamma - 1.0) * mach ** 2) / (gamma * mach ** 2 - 0.5 * (gamma - 1.0)))


def oblique_shock(Mach, theta, show=True, msg='oblique shock', gamma=gamma):
    """
    solve oblique shock wave
    
    Args:
        * mach (real) : Mach number
        * theta (real): deviation angle (positive) in degrees
        * show (boolean) : to show details
        * msg (string) : message on screen
        * gamma (real) : :math:`C_p/C_v`
    
    Returns: 
    
        real: sigma, down_Mach, p_ratio, omega_upstream, omega_downstream
        
        * sigma : shock angle in radians
        * down_Mach : downstream Mach number across the shock
        * p_ratio : ratio of static pressure across the shock
        * omega_upstream : upstream Prandtl-Meyger function in radians
        * omega_downstream : downstream Prandtl-Meyger function in radians
    """
    print("\n %s " % (msg))
    sigma = np.deg2rad(valeur_sigma(Mach, theta, sigma=0, show=show, gamma=gamma))
    print("Upstream Mach number            : %7.4f " % (Mach))
    print("Shock angle                     : %7.4f °" % (np.rad2deg(sigma)))
    upstream_normal_Mach = Mach * np.sin(sigma)
    downstream_normal_Mach = downstream_Mach(upstream_normal_Mach)
    print("Upstream normal Mach number     : %7.4f" % (upstream_normal_Mach))
    down_Mach = downstream_normal_Mach / np.sin(sigma - np.deg2rad(theta))
    print("Downstream Mach number          : %7.4f " % (down_Mach))
    p_ratio = P2_P1(upstream_normal_Mach)
    print("Static pressure ratio           : %7.4f " % (p_ratio))
    omega_upstream = omega_super(Mach)
    print("Upstream angle omega            : %7.4f °" % (np.rad2deg(omega_upstream)))
    omega_downstream = omega_super(down_Mach)
    print("Downstream angle omega          : %7.4f °" % (np.rad2deg(omega_downstream)))
    # the isentropic ratio across the shock should be added
    return sigma, down_Mach, p_ratio, omega_upstream, omega_downstream


def theta_from_Mn0(sigma, Mach, gamma=gamma):
    """
    deviation angle theta function of the  Mach number and of the shock angle for an oblique shock wave
    
    Args:
        * sigma (real) : shock angle in radians
        * Mach (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real: deviation angle theta in radians
    """
    Mn0_2 = (Mach * np.sin(sigma)) ** 2
    t = (2. + (gamma - 1) * Mn0_2) / ((gamma + 1) * Mn0_2)
    return sigma - np.arctan(t * np.tan(sigma))


def inv_pi2_over_pi1(ratio, Mach=2., show=False, gamma=gamma):
    """
    get  the  Mach number from the ratio  pi2/pi1 across a normal shock, using a Newton method
    
    Args:
        * ratio (real) : ratio of isentropic pressure across the shock
        * Mach (real) : initial guess of the Mach number
        * show (boolean) : to display details on screen
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : Mach number
    """
    er0, ds, erreur = 1e-6, 1e-3, 1.0
    iterMax, i = 20, 1
    if show:
        print("Iterations \t Mach \t dMach  \t errorr")

    while (abs(erreur) > er0) and (i <= iterMax):
        f0 = pi2_pi1(Mach, gamma=gamma)
        dM = Mach * ds
        f1 = pi2_pi1(Mach + dM, gamma=gamma)
        dMach = - dM * (f0 - ratio) / (f1 - f0)
        erreur = dMach / Mach
        if show:
            print(' %2i \t %12.6f \t %12.6f \t %12.6f' % (i, Mach, dMach, erreur))
        Mach += dMach
        i += 1
    if i > iterMax:
        print('no  convergence')
        Mach = -1
    elif Mach < 0:
        print('No convergence, Mach < 0')
        Mach = -1
    elif (show):
        print('convergence reached for the  Mach number')
    return Mach


# ********   OBLIQUE SHOCK WAVE (end)   *********


# ********   1D SUPERSONIC FLOW   *********

def S_over_Scrit(M, gamma=gamma):
    """
    S/S_critical
    
    Args:
        * Mach (real) : Mach number
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : ratio :math:`S/S_{crit}`
    """
    omega = 1 + (gamma - 1.) / 2. * M ** 2
    n = (gamma + 1) / (2 * (gamma - 1))
    return pow(2 * omega / (gamma + 1), n) / M


def inv_S_over_Scrit(ratio, Mach=2., show=False, gamma=gamma):
    """
    get Mach number from the section ration S/S_critical with a Newton method
    
    Args:
        * ratio (real) : section ratio :math:`S/S_{crit}`
        * Mach (real) : initial guess of the Mach number
        * show (boolean) : to show details on screen
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : Mach number
    """
    er0, ds, erreur = 1e-6, 1e-3, 1.0
    iterMax, i = 20, 1
    if show:
        print("Iterations \t Mach \t dMach  \t errorr")

    while (abs(erreur) > er0) and (i <= iterMax):
        f0 = S_over_Scrit(Mach, gamma=gamma)
        dM = Mach * ds
        f1 = S_over_Scrit(Mach + dM, gamma=gamma)
        dMach = -dM * (f0 - ratio) / (f1 - f0)
        erreur = dMach / Mach
        if show:
            print(' %2i \t %12.6f \t %12.6f \t %12.6f' % (i, Mach, dMach, erreur))
        Mach += dMach
        i += 1
    if i > iterMax:
        print('no convergence')
        Mach = -1
    elif Mach < 0:
        print('No convergence, Mach < 0')
        Mach = -1
    elif (show):
        print('convergence reached on the Mach number')
    return Mach


def qm_1D(S, Mach, pi, Ti, gamma=gamma, r=r_Air):
    """
    mass flow rate in a 2D nozzle :math:`q_m`
    
    Args:
        * S (real) : section [m^2]
        * Mach (real) : Mach number
        * pi (real) : isentropic pressure [Pa]
        * Ti (real) : isentropic temperature [K]
        * gamma (real) : :math:`C_p/C_v`
        * r (real) : :math:`Cp - Cv`
        
    Returns:
        real :  mass flow rate [kg/s]
    """
    omega = pow(1 + (gamma - 1) / 2 * Mach ** 2, -0.5 * (gamma + 1) / (gamma - 1))
    return S * np.sqrt(gamma / r) * pi / np.sqrt(Ti) * omega * Mach


# ********   CHARACTERISTIC METHOD *********

def mu2mach(mu):
    """
    get Mach number from mu angle = arcsin (1/Mach)
    
    Args:
        mu (real) : Mach angle in radians
        
    Returns:
        real : Mach number
    """
    # return np.sqrt(1/np.tan(mu)**2+1)
    return 1 / np.sin(mu)


def omega_mu(mu, gam=gamma):
    """
    Prandtl-Meyer (omega) function of mu, angles in radians
    
    Args:
        * mu (real) : Mach angle in radians
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : Prandtl-Meyer angle = omega in radians
    """
    c = np.sqrt((gam + 1) / (gam - 1))
    return mu + c * np.arctan(1. / (np.tan(mu) * c)) - np.pi / 2


def mu_omega(omega, gam=gamma):
    """
    get Mach angle (mu) from Prandtl-Meyer function (omega) with a Newton method
    
    Args:
        * omega (real) :  Prandtl-Meyer angle = omega in radians
        * gamma (real) : :math:`C_p/C_v`
    
    Returns:
        real : mu, Mach angle in radians 
    
    """
    F = 1
    it = 0
    eps = 1e-6
    #  stating point : approximate formula
    mu = np.pi / 2 * (1 - pow((omega / 2.27), 0.25))

    while (abs(F) > 1e-10) and (it < 20):
        F = (omega_mu(mu, gam=gam) - omega)
        dF = (omega_mu(mu + eps, gam=gam) - omega_mu(mu - eps, gam=gam)) / (2 * eps)
        mu += -F / dF
        it += 1
    return mu


def line_intersection(M1, M2, p1, p2):
    """
    Intersection of two lines :
        * D1 ( point M1 and slope  p1 )
        * D2 ( point M2 and slope  p2 )
    """
    x = (M2[1] - M1[1] + p1 * M1[0] - p2 * M2[0]) / (p1 - p2)
    y = M1[1] + (x - M1[0]) * p1
    return [x, y]


def complex_line_intersection(thetaMinus, omegaMinus, thetaPlus, omegaPlus):
    """
    Intersection of two  characteristic lines
    C^-  : thetaMinus + omegaMinus = lambdaMinus
    C^+  : omegaPlus - thetaPlus  = lambdaPlus
    
    Args:
        * thetaMinus (real) : theta on the negative characteristic
        * omegaMinus (real) : theta on the negative characteristic
        * thetaPlus (real) : theta on the positive characteristic
        * omegaPlus (real) : theta on the positive characteristic
    
    Returns:
        real : ome, the
        * omega_PM : Prandtl-Meyer angle, omega in radians
        * deviation : deviation angle, theta in radians
    """
    lambdaPlus, lambdaMinus = omegaPlus - thetaPlus, omegaMinus + thetaMinus
    omega_PM, deviation = (lambdaPlus + lambdaMinus) / 2, (lambdaMinus - lambdaPlus) / 2
    return omega_PM, deviation
