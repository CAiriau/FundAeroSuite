#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 14 exercises
    ..
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Tools.misc import set_title, set_question
from CompressibleFlow.fonctions import gamma, oblique_shock, valeur_sigma, P2_P1, downstream_Mach
from CompressibleFlow.fonctions import omega_hyper, invomega, p_pi, pi2_pi1, inverse_p_pi, omega_super
from scipy.interpolate import interp1d
from conical_shock.conical_shock import main_multiple_Mach


def Constant_Kc(K0):
    """
    Oblique shock wave, shock angle, Kc= M0 sigma,  hypersonic approximation
    """
    tmp = (gamma + 1.0) / 4.0
    return K0 * (tmp + np.sqrt(tmp ** 2 + 1.0 / K0 ** 2))


def Pc_P0(Kc):
    """
    Oblique shock wave, pc/p0, hypersonic approximation
    """
    return (2.0 * gamma * Kc ** 2 - (gamma - 1.0)) / (gamma + 1.0)


def Kp_hyper(Mach, p_p0):
    """
    Kp in hypersonic regime
    """
    return 2 / (gamma * Mach ** 2) * (p_p0 - 1)


def approx_omega(Mach, M0):
    """
    the approximation validity is tested
    """
    f1 = omega_super(M0) - omega_super(Mach)
    f2 = omega_hyper(M0) - omega_hyper(Mach)
    return np.rad2deg(f1), np.rad2deg(f2)


def Solution_Lees():
    """
    LEES  solution in supersonic flow for cones 
    
    Returns:
        real : Kp/theta, sigma/theta
    """
    return 2 * (gamma + 1) * (gamma + 7) / (gamma + 3) ** 2, 2 * (gamma + 1) / (gamma + 3)


def detached_shock_length():
    """
    detached shock value, Delta/R
    """
    epsilon = (gamma - 1.) / (gamma + 1.)
    Gam = np.sqrt((gamma - 1.) * (gamma + 3.)) / (gamma + 1)
    return epsilon / (1 + Gam - 2 * epsilon)


def downstream_Mach_hyper(Mach, Kc=0, sigma=0, regime='hypersonic'):
    """
    Downstream Mach for a hypersonic shock wave, as a function od  Kc (sigma) and Mach numbers
    In "supersonic" Kc=sigma", where the shock angle is in radians and K = M0*sin(sigma)
    """
    if regime == 'hypersonic':
        K = Kc
    else:
        K = Mach * np.sin(sigma)

    Num = ((gamma + 1) * K * Mach) ** 2 - 4 * (K ** 2 - 1) * (gamma * K ** 2 + 1)
    Den = (2 * gamma * K ** 2 - gamma + 1) * (2 + (gamma - 1) * K ** 2)
    return np.sqrt(Num / Den)


def downstream_Mach_hyper_from_K0(Mach, K0, Kc, method=1):
    """
    Downstream Mach for a hypersonic shock wave
    second method, with K0 (theta) and Kc, et Mach
    """
    K0New = 2 * (Kc ** 2 - 1) / (Kc * (gamma + 1))
    print("Verification of K0 wrt to Kc: new K0 = %f, \t old one = %f" % (K0New, K0))
    if method == 1:
        Den = (2. + (gamma - 1.) * Kc ** 2) * (2 * gamma * Kc ** 2 + 1 - gamma)
        Num = (gamma + 1) * Mach * Kc
        return Num / np.sqrt(Den)
    else:
        Num = np.sqrt((2 + (gamma - 1.) * Kc ** 2) / (2 * gamma * Kc ** 2 + 1 - gamma))
        return Mach / (Kc - K0) * Num


def downstream_Mach_asymptotic(theta):
    """
    asymptotic values of downstream Mach for an hypersonic shock wave, theta en radians
    """
    c = np.sqrt(2 / (gamma * (gamma - 1)))
    # there are two approximations, depending of neglecting gamma or not
    # the first one is the best, for small deviation angles the differences are very small
    return c * np.sqrt(1 / theta ** 2 - gamma), c / theta


def Pj_Pc(Mc, deltaTheta):
    """
    P_j/P_c
    """
    tmp = 1 + (gamma - 1) / 2 * Mc * deltaTheta
    expo = 2 * gamma / (gamma - 1)
    return np.power(tmp, expo)


def pressure_ratio_from_Kp(Kp, Mach):
    """
    p/p0 knowing  Kp
    """
    return 1 + gamma / 2 * Mach ** 2 * Kp


# additional used functions

def diamond_hyper_coef(alpha, tau):
    """
    Diamond airfoil in hypersonic, numerical solution.
    """
    # 0 index is not used, the other index numbered the zone
    delta = np.array([alpha, tau, -tau, -tau, tau]) - alpha
    Kp = 2 * delta ** 2
    if alpha < tau:
        Kp[3:] = 0
    else:
        Kp[1:4:2] = 0
    # print('Kp = ',Kp[1:])
    CL = 0.5 * (Kp[2] + Kp[4] - Kp[1] - Kp[3])
    # print('CL = ',CL)
    CD = 0.5 * (-Kp[2] * delta[2] - Kp[4] * delta[4] + Kp[1] * delta[1] + Kp[3] * delta[3])
    # print('CD = ',CD)
    F = CL / CD  # Aerodynamic efficient
    return Kp, CL, CD, F


def diamond_hyper_analytical(alpha):
    """
    Diamond airfoil in hypersonic, analytical solution.
    alpha is divided  : alpha/tau
    """
    n = len(alpha)
    CL, CD = np.zeros([n]), np.zeros([n])
    polar = np.zeros([n])
    for k in range(n):
        if alpha[k] < 1:
            CL[k], CD[k] = 4 * alpha[k], 2 + 6 * alpha[k] ** 2
            polar[k] = 2 + 3 / 8 * CL[k]
        else:
            CL[k], CD[k] = 2 + 2 * alpha[k] ** 2, 6 * alpha[k] + 2 * alpha[k] ** 3
            polar[k] = (4 + CL[k]) * np.sqrt(0.5 * CL[k] - 1)
    fmax = 1 / np.sqrt(3)
    CL_opt, CD_opt = 4 / np.sqrt(3), 4
    alpha_opt = np.rad2deg(1 / np.sqrt(3))
    CD_min = 2

    return CL, CD, CL / CD, polar, fmax, CL_opt, CD_opt, alpha_opt, CD_min


def triangle_hyper_coef(alpha, tau):
    """
    Diamond airfoil in hypersonic,  numerical solution.
    """
    # # 0 index is not used, the other index numbered the zone
    delta = np.array([alpha, tau, 0, -tau]) - alpha
    Kp = 2 * delta ** 2
    if alpha < tau:
        Kp[3] = 0
    else:
        Kp[1], Kp[3] = 0, 0
    CL = Kp[2] - (Kp[1] + Kp[3]) / 2
    CD = (Kp[1] * delta[1] + Kp[3] * delta[3]) / 2 - Kp[2] * delta[2]
    f = CL / CD
    return Kp, CL, CD, f


def triangle_hyper_analytical(alpha):
    """
    Triangle airfoil in hypersonic, analytical solution.
    alpha is divided  : alpha/tau
    """
    n = len(alpha)
    CL, CD = np.zeros([n]), np.zeros([n])
    polar = np.zeros([n])
    for k in range(n):
        if alpha[k] < 1:
            CL[k], CD[k] = alpha[k] ** 2 + 2 * alpha[k] - 1, 1 - 3 * alpha[k] + 3 * alpha[k] ** 2 + alpha[k] ** 3
            polar[k] = 6 + np.sqrt(2 + CL[k]) * (CL[k] - 4)
        else:
            CL[k], CD[k] = 2 * alpha[k] ** 2, 2 * alpha[k] ** 3
            polar[k] = pow(CL[k], 3. / 2.) / np.sqrt(2)
    fmax = 1.248534555
    CL_opt, CD_opt = 0.911027108, 0.729677128
    alpha_opt = 40.46074596
    CD_min = 6 - 4 * np.sqrt(2)
    return CL, CD, CL / CD, polar, fmax, CL_opt, CD_opt, alpha_opt, CD_min


def Exercice14_0():
    """ 
    Downstream Mach of an oblique shock wave over a plane wall with  different approximations
    """
    set_title("Downstream Mach of an oblique shock wave over a plane wall with  different approximations")
    M0 = 100
    theta = 5.0

    print("Exact theory of oblique shock wave")
    sigma, Mdown, Pdown, omegaupst, omegadown = oblique_shock(M0, theta, show=False)
    print("theta in degrees           :", theta)
    print("theta in radians           :", np.deg2rad(theta))
    print("sigma en radians           :", sigma)
    print("sigma/sin sigma            :", sigma / np.sin(sigma))
    print("(sigma-theta)/sin (sigma)  :", (sigma - np.deg2rad(theta)) / np.sin(sigma))
    K0, Kc = M0 * np.deg2rad(theta), M0 * sigma
    print("K0                         :", K0)
    print("Kc                         :", Kc)

    M1_hyper = downstream_Mach_hyper(M0, Kc=Kc, regime='hypersonic')
    M1_super = downstream_Mach_hyper(M0, sigma=sigma, regime='supersonic')
    print("M1 approx. hypersonic      :", M1_hyper)
    print("M1 approx. supersonic      :", M1_super)

    print('other calculus')
    Kc1 = Constant_Kc(K0)
    print('Kc approximation, M0 infinite : ', Kc1)
    print('sigma (M0 infinite )(°)      : ', np.rad2deg(Kc1 / M0))

    # M0 large, theta et sigma small
    M1_1 = downstream_Mach_hyper_from_K0(M0, K0, Kc1, method=1)
    M1_2 = downstream_Mach_hyper_from_K0(M0, K0, Kc1, method=2)
    print('M1 method  1                : ', M1_1)
    print('M1 method  2                : ', M1_2)
    print('M1 asymptotic : ')
    print(downstream_Mach_asymptotic(np.deg2rad(theta)))


def Exercice14_1():
    """ 
    Diamond and triangular airfoils Kp: Newton's method
    """
    set_title("Diamond and triangular airfoils Kp: Newton's method")
    plot = True
    plot_num = True     # to plot numerical solution
    plot_ana = True     # to plot analytical solution
    pas = 100           # to plot analytical solution with a given step in the points drawn
    plot_max = True     # to plot the point of maximal aerodynamic efficiency
    plot_gmtry = True   # to plot the both airfoils
    tau_diamond = 0.1
    tau_triangle = 0.1
    Msize = 8           # symbol size  of point of maximal aerodynamic efficiency
    lw = 2              # line with
    n = 1001
    alpha = np.linspace(0, 2 * tau_diamond, n)  # angle of attack 
    CL = np.zeros([n, 2])
    CD = np.zeros([n, 2])
    F = np.zeros([n, 2])
    IndMax = np.zeros([2], dtype=int)

    for k in range(n):
        KP1, CL[k, 0], CD[k, 0], F[k, 0] = diamond_hyper_coef(alpha[k], tau_diamond)
        KP2, CL[k, 1], CD[k, 1], F[k, 1] = triangle_hyper_coef(alpha[k], tau_triangle)

        IndMax[0], IndMax[1] = np.argmax(F[:, 0]), np.argmax(F[:, 1])

    print("max. aerodynamic efficiency : Los. : %f, \t Tri. : %f" % (F[IndMax[0], 0], F[IndMax[1], 1]))
    print("Incidence in °   : Los. : %f, \t Tri. : %f" % (np.rad2deg(alpha[IndMax[0]]), np.rad2deg(alpha[IndMax[1]])))
    print("CL               : Los. : %f, \t Tri. : %f" % (CL[IndMax[0], 0], CL[IndMax[1], 1]))

    print("Analytical solution")
    lCL, lCD, lF, lpolar, lfmax, lCL_opt, lCD_opt, lalpha_opt, lCD_min = diamond_hyper_analytical(alpha / tau_diamond)
    tCL, tCD, tF, tpolar, tfmax, tCL_opt, tCD_opt, talpha_opt, tCD_min = triangle_hyper_analytical(
        alpha / tau_triangle)

    print("max. aerodynamic efficiency : Los. : %f, \t Tri. : %f" % (lfmax / tau_diamond, tfmax / tau_triangle))
    print("  Incidence in ° : Los. : %f, \t Tri. : %f" % (lalpha_opt * tau_diamond, talpha_opt * tau_triangle))
    print("  CL             : Los. : %f, \t Tri. : %f" % (lCL_opt * tau_diamond ** 2, tCL_opt * tau_triangle ** 2))
    print("  CD             : Los. : %f, \t Tri. : %f" % (lCD_opt * tau_diamond ** 2, tCD_opt * tau_triangle ** 2))
    print("CD minimal       : Los. : %f, \t Tri. : %f" % (lCD_min * tau_diamond ** 2, tCD_min * tau_triangle ** 2))

    xl = alpha / tau_diamond
    xt = alpha / tau_triangle
    xL = xl[0:n:pas]

    if plot:
        # CD et CL vs alpha
        fig = plt.figure(1)
        fig.suptitle(r'$\tilde C_L$  and $\tilde C_D$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde \alpha$', fontsize=20)
        ax.set_ylabel(r'$\tilde C_D, \tilde C_L$', fontsize=20)
        ax.grid()

        if plot_num:
            ax.plot(xl, CL[:, 0] / tau_diamond ** 2, '-', label=r'Los.: $\tilde C_L$', linewidth=lw)
            ax.plot(xl, CD[:, 0] / tau_diamond ** 3, '--', label=r'Los.: $\tilde C_D$', linewidth=lw)
            ax.plot(xt, CL[:, 1] / tau_triangle ** 2, '-', label=r'Tri.: $\tilde C_L$', linewidth=lw)
            ax.plot(xt, CD[:, 1] / tau_triangle ** 3, '--', label=r'Tri.: $\tilde C_D$', linewidth=lw)
        if plot_ana:
            ax.plot(xL, lCL[0:n:pas], 's', label=r'Los. an.: $\tilde C_L$')
            ax.plot(xL, lCD[0:n:pas], 'o', label=r'Los. an.: $\tilde C_D$')
            ax.plot(xL, tCL[0:n:pas], 's', label=r'Tri. an.: $\tilde C_L$')
            ax.plot(xL, tCD[0:n:pas], 'o', label=r'Tri. an.: $\tilde C_D$')
        if plot_max:
            ax.plot(np.deg2rad(lalpha_opt), lCL_opt, 's', markersize=Msize, label=r'Los.: $\tilde C_L(f_{\max})$')
            ax.plot(np.deg2rad(lalpha_opt), lCD_opt, 's', markersize=Msize, label=r'Los.: $\tilde C_D(f_{\max})$')
            ax.plot(np.deg2rad(talpha_opt), tCL_opt, 's', markersize=Msize, label=r'Tri.: $\tilde C_L(f_{\max})$')
            ax.plot(np.deg2rad(talpha_opt), tCD_opt, 's', markersize=Msize, label=r'Tri.: $\tilde C_D(f_{\max})$')
        ax.legend(loc='upper left')

        # aerodynamic efficiency f  vs CL
        fig = plt.figure(2)
        fig.suptitle(r'Finesse $\times \tau$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde C_L$', fontsize=20)
        ax.set_ylabel(r'$\tilde f$', fontsize=20)
        ax.grid()
        if plot_num:
            ax.plot(CL[:, 0] / tau_diamond ** 2, F[:, 0] * tau_diamond, '-', label="Los.", linewidth=lw)
            ax.plot(CL[:, 1] / tau_triangle ** 2, F[:, 1] * tau_triangle, '-', label="Tri.", linewidth=lw)
        if plot_ana:
            ax.plot(lCL[0:n:pas], lF[0:n:pas], 's', label="Los. an.")
            ax.plot(tCL[0:n:pas], tF[0:n:pas], 'o', label="Tri. an.")
        if plot_max:
            ax.plot(lCL_opt, lfmax, 's', markersize=Msize, label=r'Los.: $f_{\max}$')
            ax.plot(tCL_opt, tfmax, 's', markersize=Msize, label=r'Tri.: $f_{\max}$')
        ax.legend(loc='upper right')

        # aerodynamic efficiency f  vs alpha
        fig = plt.figure(3)
        fig.suptitle(r'aerodynamic efficiency $\times \tau$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde \alpha$', fontsize=20)
        ax.set_ylabel(r'$\tilde f$', fontsize=20)
        ax.grid()
        if plot_num:
            ax.plot(xl, F[:, 0] * tau_diamond, '-', label="Los.", linewidth=lw)
            ax.plot(xl, F[:, 1] * tau_triangle, '-', label="Tri.", linewidth=lw)
        if plot_ana:
            ax.plot(xL, lF[0:n:pas], 's', label="Los. an.")
            ax.plot(xL, tF[0:n:pas], 'o', label="Tri. an.")
        if plot_max:
            ax.plot(np.deg2rad(lalpha_opt), lfmax, 's', markersize=Msize, label=r'Los.: $f_{\max}$')
            ax.plot(np.deg2rad(talpha_opt), tfmax, 's', markersize=Msize, label=r'Tri.: $f_{\max}$')
        ax.legend(loc='lower right')

        # polar CD vs CL
        fig = plt.figure(4)
        fig.suptitle(r'Polar $\tilde C_D$=f(\tilde C_L)', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\tilde C_L$', fontsize=20)
        ax.set_ylabel(r'$\tilde C_D$', fontsize=20)
        ax.grid()
        if plot_num:
            ax.plot(CD[:, 0] / tau_diamond ** 3, CL[:, 0] / tau_diamond ** 2, '--', label=r'Los.', linewidth=lw)
            ax.plot(CD[:, 1] / tau_triangle ** 3, CL[:, 1] / tau_triangle ** 2, '--', label=r'Tri.', linewidth=lw)
        if plot_ana:
            ax.plot(lCD[0:n:pas], lCL[0:n:pas], 'o', label=r'los. an.')
            ax.plot(tCD[0:n:pas], tCL[0:n:pas], 'o', label=r'tri. an.')
        if plot_max:
            ax.plot(lCD_opt, lCL_opt, 's', markersize=Msize, label=r'Los.: $f_{\max}$')
            ax.plot(tCD_opt, tCL_opt, 's', markersize=Msize, label=r'Tri.: $ f_{\max}$')
        ax.legend(loc='upper left')

        # airfoil geometry
        if plot_gmtry:
            fig = plt.figure(5)
            fig.suptitle('Airfoils', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/\ell$', fontsize=20)
            ax.set_ylabel(r'$y/\ell$', fontsize=20)
            ax.axis('equal')
            ax.grid()
            ax.plot([0, 1 / 2, 1, 1 / 2, 0], [0, tau_diamond, 0, -tau_diamond, 0], '-', linewidth=lw)
            dec = 0.3
            ax.plot([0, 1 / 2, 1, 0], [0 + dec, tau_triangle + dec, 0 + dec, 0 + dec], '-', linewidth=lw)

        plt.show()


def Exercice14_3():
    
    """ 
    Kp of a 2D wedge or conical wall, for various Mach number
    """
    set_title("Kp of a 2D wedge or conical wall, for various Mach number")

    M0 = 100
    theta = np.array([5, 10, 15, 20])
    sigma = np.zeros([4])
    Kp = np.zeros([4])

    set_question('1 : wedge')

    for k in range(4):
        sigma[k] = valeur_sigma(M0, theta[k])
        rap_p = P2_P1(M0 * np.sin(np.deg2rad(sigma[k])))
        Kp[k] = Kp_hyper(M0, rap_p)

    print('for a finite Mach number M=%f, oblique shock exact solution : ' % (M0))
    print("theta    : ", theta)
    print("sigma    : ", sigma)
    print("100 Kp   : ", Kp * 100)

    print('for a finite Mach number M=%f, approximated Kp :' % (M0))
    Kc = Constant_Kc(np.deg2rad(theta) * M0)
    Kp = Kp_hyper(M0, Pc_P0(Kc))
    print("theta    : ", theta)
    print("sigma    : ", np.rad2deg(Kc / M0))
    print("100 Kp   : ", Kp * 100)

    print('for a infinite Mach number ')
    print("theta    : ", theta)
    print("sigma    : ", (gamma + 1) / 2 * theta)
    print("100 Kp   : ", (gamma + 1) * (np.deg2rad(theta)) ** 2 * 100)

    set_question('2 : conical shock wave')
    """ 
    Conical body in supersonic flow  
    """
    # interpolation in the database
    set_title("Interpolation in the conical shock wave tables for super/hyper-sonic flows")
    reference_filepath = os.path.join('Book/Data', 'Mach_50_choc_conique.dat')
    lw = 2
    m = 1
    with open(reference_filepath, 'r') as infile:
        A = np.loadtxt(infile, dtype=float, usecols=(0, 1, 2), unpack=True, skiprows=m)
        theta1 = A[0, :]
        kp_cc = A[2, :]
        sigma1 = A[1, :]
    print(A.shape)
    print(theta1.shape)
    f1 = interp1d(theta1, kp_cc, kind='cubic')
    f2 = interp1d(theta1, sigma1, kind='cubic')
    theta_c = np.array([5.0, 10.0, 15.0, 20.0])

    print('for an infinite Mach number')
    print("theta    : ", theta_c)
    print("sigma    : ", f2(theta_c))
    print("100 Kp   : ", 100 * f1(theta_c))

    set_question("3 : Newton's method + plots")

    print('Kp = 2 sin^2(theta)')
    print("theta            : ", theta)
    print("100 Kp (sin^2)   : ", 100 * 2 * np.sin(np.deg2rad(theta)) ** 2)
    print("100 Kp (theta^2) : ", 100 * 2 * (np.deg2rad(theta)) ** 2)
    n = 21
    theta = np.linspace(2, 20, n)
    # KpWedge,KpCone,KpNewtonSin,KpNewtonApprox=np.zeros([n]),np.zeros([n]),np.zeros([n]),np.zeros([n])
    # sigmaWedge,sigmaCone=np.zeros([n]),np.zeros([n])
    sigmaWedge = (gamma + 1) / 2 * theta
    sigmaCone = f2(theta)
    KpWedge = (gamma + 1) * (np.deg2rad(theta)) ** 2 * 100
    KpCone = 100 * f1(theta)
    KpNewtonSin = 100 * 2 * np.sin(np.deg2rad(theta)) ** 2
    KpNewtonApprox = 100 * 2 * (np.deg2rad(theta)) ** 2

    Kp_Lees, sigma_Lees = Solution_Lees()
    fig = plt.figure(5)
    fig.suptitle('Pressure coefficients', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
    ax.set_ylabel(r'$100\times Kp$ or $\sigma(^\circ)$', fontsize=20)
    ax.grid()
    ax.plot(theta, sigmaWedge, '--', linewidth=lw, label=r"$\sigma$ : wedge")
    ax.plot(theta, sigmaCone, '--', linewidth=lw, label=r"$\sigma$ : cone")
    ax.plot(theta, KpWedge, '-', linewidth=lw, label=r"$Kp$ : wedge")
    ax.plot(theta, KpCone, '-', linewidth=lw, label=r"$Kp$ : cone")
    ax.plot(theta, KpNewtonSin, '-', linewidth=lw, label=r"$Kp : 2\sin^2\theta$")
    ax.plot(theta, KpNewtonApprox, '-', linewidth=lw, label=r"$Kp : 2 \theta^2$")
    ax.plot(theta, 100 * Kp_Lees * (np.deg2rad(theta)) ** 2, '-', linewidth=lw, label=r"$Kp : Lees$")
    ax.plot(theta, sigma_Lees * theta, '--', linewidth=lw, label=r"$\sigma$ : Lees")

    ax.legend(loc='upper left')
    plt.show()

    set_question('3 :LEES solution')
    theta = np.array([5, 10, 15, 20])
    theta = np.linspace(5, 20, 16)
    print('Solution independant of the Mach number')
    print("theta    : ", theta)

    print("sigma (°): ", sigma_Lees * theta)
    print("100 Kp   : ", 100 * Kp_Lees * (np.deg2rad(theta)) ** 2)
    print("error Kp :", np.abs(Kp_Lees * (np.deg2rad(theta)) ** 2 / f1(theta) - 1))
    print("error sigma :", np.abs(sigma_Lees * theta / f2(theta) - 1))


def ogive_parabolique(e_over_L, x, pente="angle"):
    """
    Airfoil, slope (camber derivative), wall deviation, radius of curvature
    """
    y = -e_over_L * 2. * x * (x - 1)  # airfoil camber line
    yp = (-4 * x + 2) * e_over_L  # first derivative
    ypp = -4 * e_over_L  # second derivative
    R = abs(pow(1 + yp ** 2, 3 / 2) / ypp)  # radius of curvature
    # slope in  radian       
    if pente == "angle":
        th = yp  # small deviation are assumed
    else:
        th = np.arctan(yp)
        # print('Theta en ° =',np.rad2deg(th))
    return y, th, R


def Kp_Newton(e_over_L, C, x):
    """
    Hypersonic pressure coefficient : Newton's formula
    """
    y, th, R = ogive_parabolique(e_over_L, x)
    print("Newton")
    Kp = C * np.sin(th[th > 0]) ** 2
    for k in np.arange(len(Kp), len(x)):
        Kp = np.append(Kp, [0])
    return Kp


def Kp_Newton_Busemann(e_over_L, C, x):
    """
    Hypersonic pressure coefficient : Newton's formula +
    Busemann's correction
    """
    y, th, R = ogive_parabolique(e_over_L, x)
    print("Newton Busemann")
    ind = th > 0
    Kp = C * np.sin(th[ind]) ** 2 + 2 * y[ind] / R[ind]
    for k in np.arange(len(Kp), len(x)):
        Kp = np.append(Kp, [0])
    return Kp


def Exercice14_4():
    """ 
    Sharp obstacle Kp with different approaches
    """
    set_title("Sharp obstacle Kp with different approaches")
    plot_profil = False
    plot_curves = True
    plot_Newton = True
    n = 201  # number of points on the airfoil upper surface 
    Mach_inf = np.array([3.5, 10, 100.0])
    m = len(Mach_inf)
    e_over_L = 0.10

    theta, x, Kp = np.zeros(n), np.zeros(n), np.zeros(n)
    ye, yi = np.zeros(n), np.zeros(n)
    x = np.linspace(0, 1, n)

    # ---------------------------------------
    set_question('1 : lenticular airfoil')
    # ---------------------------------------

    ye, theta, Re = ogive_parabolique(e_over_L, x, pente="angle")
    yi = -ye
    # for i in range(n-1):
    #     theta[i]=np.arctan((ye[i+1]-ye[i])/(x[i+1]-x[i]))
    #     #print(i,theta[i])
    print('x(0)= %f, theta(0) = %f ' % (x[0], theta[0]))

    if plot_profil:
        fig = plt.figure()
        fig.suptitle('Airfoil', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$', fontsize=20)
        ax.set_ylabel(r'$y$', fontsize=20)
        ax.axis('equal')
        ax.grid()
        ax.plot(x, ye, '-')
        ax.plot(x, yi, '-')
        ax.plot(x[:-1], theta[:-1], 'o--')
        plt.show()

    # -----------------------------------------------------    
    set_question("2: shock-expansion method, infinite Mach number")
    # -----------------------------------------------------        

    theta_min = min(theta - theta[0])
    print('Min Theta on the airfoil         : %f, %f °' % (theta_min, np.rad2deg(theta_min)))
    Kp = np.zeros([m, n])
    # Hypersonic flow theory:
    print("\n ", "-" * 49, "\nINFINITE MACH NUMBER result \n", "_" * 50)
    sigma = theta[0] * (gamma + 1) / 2
    KpcInf = (gamma + 1) * theta[0] ** 2
    Mc, tmp = downstream_Mach_asymptotic(theta[0])
    Kp_infinite = (gamma + 1) * theta[0] ** 2 * Pj_Pc(Mc, theta - theta[0])
    print("sigma infinite case             : %7.4f °" % (np.rad2deg(sigma)))
    print('100 x Kpc theo.  M_infinite     : %7.4f ' % (100 * KpcInf))
    print("Downstream Mach  (M0 inf, theo.): %7.4f " % Mc)
    print(" theta (°) \t 100 Kp")
    for thetap, kp in zip(theta, Kp_infinite):
        print(" %7.4f \t %7.4f" % (np.rad2deg(thetap), 100 * kp))
    print(100 * Kp_infinite)

    # ------------------------------------
    set_question("3: Newton's method")
    # ------------------------------------  
    C = KpcInf / np.sin(theta[0]) ** 2
    print('C = %f, C theo.  = %f' % (C, gamma + 1))
    # C=gamma+1
    Kp_N = Kp_Newton(e_over_L, C, x)
    Kp_NB = Kp_Newton_Busemann(e_over_L, C, x)

    print(" theta (°) \t 100 Kp_N \t 100 Kp_NB \t 100 Kp_N approx.")
    for thetap, kp1, kp2 in zip(theta, Kp_N, Kp_NB):
        if thetap >= 0:
            kp3 = (gamma + 1) * thetap ** 2
        else:
            kp3 = 0
        print(" %7.4f \t %7.4f \t %7.4f \t %7.4f" % (np.rad2deg(thetap), 100 * kp1, 100 * kp2, 100 * kp3))

    # ------------------------------------------------------    
    set_question("4: methods for FINITE Mach number")
    # ------------------------------------------------------

    print("\n ", "-" * 49, "\n RESULTS for FINITE UPSTREAM MACH NUMBER  \n", "_" * 50)
    for k in range(m):
        M0 = Mach_inf[k]
        Theta_limit = -2 / (gamma - 1) / M0
        print("\n", "*" * 49)
        print("Upstream Mach M0                  : %7.4f " % M0)
        print("*" * 50, "\n")
        print('Minimal Theta allowed             : %f, %f °' % (Theta_limit, np.rad2deg(Theta_limit)))

        # if theta_min < Theta_limit:
        #     x_bad=x[theta-theta[0]<Theta_limit]
        #     print("PROBLEM : THE THEORY NO LONGER WORKS : DECREASE  M0 or airfoil thickness")
        #     print("from  x = \n",x_bad[0])
        print("Angle at the leading edge         : %7.4f ° , %7.4f rad " % (np.rad2deg(theta[0]), theta[0]))
        print('sin theta                         : %7.4f' % (np.sin(theta[0])))
        K0 = M0 * theta[0]
        Kc = Constant_Kc(K0)

        print("")
        print("OBLIQUE SHOCK WAVE THEORY :")
        sigma, Mdown, Pdown, omegaupst, omegadown = oblique_shock(M0, np.rad2deg(theta[0]), show=False)
        print("Kpc                             : %7.4f " % (Kp_hyper(M0, Pdown)))
        print("sin(sigma-theta)                : %7.4f" % (np.sin(sigma - theta[0])))
        Mc = downstream_Mach_hyper(M0, Kc)
        print()
        # print("sigma (°)                       : %7.4f °"%(np.rad2deg(sigma)) )
        # print("Mach downstream                 : %7.4f  "%(Mdown) )
        # print("Pdown/Pupstream                 : %7.4f  "%(Pdown) )

        print("")
        print("HYPERSONIC THEORY :")
        print("K0 at the leading edge          : %7.4f " % K0)
        print('Kc downstream the shock         : %7.4f ' % Kc)
        print('(Kc-K0)/M0                      : %7.4f ' % ((Kc - K0) / M0))
        sigma = Kc / M0
        print("Leading edge shock angle        : %7.4f °" % (np.rad2deg(sigma)))
        Kpc = Kp_hyper(M0, Pc_P0(Kc))
        print("Kpc calculus                    : %7.4f " % Kpc)
        Mc = downstream_Mach_hyper(M0, Kc)
        print("Downtream Mach at the l.e.      : %7.4f " % Mc)

        Kp[k, 0] = Kp_hyper(M0, Pc_P0(Kc))
        for i in np.arange(1, n):
            Kp[k, i] = Kp_hyper(M0, Pj_Pc(Mc, theta[i] - theta[0]) * Pc_P0(Kc))
            # print('Kp[ %d ] = %f'%(i,Kp[i]))
        print("")
        print("NEWTON's THEORY :")
        C = Kpc / np.sin(theta[0]) ** 2
        print('C = %f' % (C))
        Kp_N = Kp_Newton(e_over_L, C, x)
        Kp_NB = Kp_Newton_Busemann(e_over_L, C, x)

        print(" theta (°) \t 100 Kp_N  \t 100 x hypersonic Kp")
        for thetap, kp1, kp2 in zip(theta, Kp_N, Kp[k, :]):
            print(" %7.4f \t %7.4f \t %7.4f " % (np.rad2deg(thetap), 100 * kp1, 100 * kp2))

    # data from MoC (to verify) 
    x_10 = np.linspace(0, 1, 11)
    Kp_10 = np.array([0.10812, 0.07437, 0.04750, 0.02656, 0.01125, 0.001562,
                      -0.005, -0.009375, -0.012187, -0.013437, -0.014375])
    Kp_inf = np.array([0.092813, 0.058125, 0.031875, 0.01625, 0.007812, 0.002187,
                       0.0, -0.000625, -0.0009375, -0.001094, -0.00125])
    Kp_35 = np.array([0.17477, 0.13169, 0.0923, 0.0578, 0.02769, 0.00246,
                      -0.01908, -0.03692, -0.05292, -0.06708, -0.07754])

    # ------------------------------------------------------    
    set_question("5: shock-expansion method in  supersonic flow")
    # ------------------------------------------------------
    print("")
    print('SUPERSONIC CASE (Reminder + Kp)\n')
    Kp_super = np.zeros([m, n])
    for k in range(m):
        M0 = Mach_inf[k]
        sigma = np.deg2rad(valeur_sigma(M0, np.rad2deg(theta[0]), sigma=0, show=False))
        print('-----')
        print("Upstream Mach M0                : %7.4f " % (M0))
        print("Shock angle                     : %7.4f °" % (np.rad2deg(sigma)))
        Mn = M0 * np.sin(sigma)
        Mndown = downstream_Mach(Mn)
        Mc = Mndown / np.sin(sigma - theta[0])
        print("Downstream Mach  Mc             : %7.4f " % (Mc))
        Pc = P2_P1(Mn)
        print("Pressure Pc                     : %7.4f " % (Pc))
        omega_c = omega_super(Mc)
        print("Angle omega_c                   : %7.4f °" % (np.rad2deg(omega_c)))

        Kp_super[k, 0] = Kp_hyper(M0, Pc)
        for i in np.arange(1, n):
            Mach = invomega(omega_c + theta[0] - theta[i])
            Kp_super[k, i] = Kp_hyper(M0, p_pi(Mach) / p_pi(Mc) * Pc)

    # ------------------------------------------------------    
    set_question("6: plots")
    # ------------------------------------------------------
    Msize = 10
    if plot_curves:

        fig = plt.figure()
        fig.suptitle('Kp en %', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$M$', fontsize=20)
        ax.set_ylabel(r'$100 Kp$', fontsize=20)
        ax.grid()
        ax.plot(x_10, 100 * Kp_35, 'ko', markersize=Msize, label=r'carac. : $M=3.5$')
        ax.plot(x_10, 100 * Kp_10, 'ko', markersize=Msize, label=r'carac. : $M=10$')
        ax.plot(x_10, 100 * Kp_inf, 'ko', markersize=Msize, label=r'carac.. : $M=\infty$')
        for k in range(m - 1):
            ax.plot(x, 100 * Kp[k, :], '-', label=r'Thé. : $M = %2.1f$' % (Mach_inf[k]))
        ax.plot(x, 100 * Kp_super[0, :], marker='<', markersize=Msize, markevery=(0, 0.1),
                label=r'Super : $M = %2.1f$' % (Mach_inf[0]))
        # ax.plot(x,100*Kp_super[1,:],marker='<',markersize=Msize,markevery=(0,0.1),label=r'Super : $M = %2.1f$'%(
        # Mach_inf[1])) 
        if plot_Newton:
            ax.plot(x, 100 * Kp_N, '-', label=r'The. : Newton')
            ax.plot(x, 100 * Kp_NB, '-', label=r'The. : Newton-Busemann')
        ax.plot(x, 100 * Kp_infinite, '-', label=r'$M=\infty$', linewidth=2)
        ax.legend(loc='upper right')
        plt.show()

    print(' x \t theta (°) \t Kp (%) \t Kp (%) \t Kp (%) \t Kp(%) \t \t  Kp(%) \t Kp(%) \t \t theta (rad)')
    print("M0:\t\t\t %3.2f \t\t %3.2f \t\t %3.2f" % (Mach_inf[0], Mach_inf[1], Mach_inf[2]),
          "\t Infini \t Newon \t \t New. Bus.")
    for k in range(n):
        m = Kp[:, k]
        print("%4.2f \t %+6.2f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f \t %+6.4f " % (
        x[k], np.rad2deg(theta[k]),
        100 * Kp[0, k], 100 * Kp[1, k], 100 * Kp[2, k], 100 * Kp_infinite[k], 100 * Kp_N[k], 100 * Kp_NB[k], theta[k]))
        if k % 5 == 0:
            print("")

            # Kp from Newton and Newton+Busemann are independent of upstream Mach number


def fx(p):
    """
    function providing  x as a function of p=dy/dx=y'
    """
    return np.log(p) + 1 / p ** 2 + 3 / (4 * p ** 4)


def fr(p):
    """
    radius distribution as a function of r p=y'
    """
    return (1 + p ** 2) ** 2 / p ** 3


def r_approx(x, C1):
    """
    approximated formula to get the body geometry, for x to infinite 
    """
    return pow(C1, 1 / 4) * pow(4. / 3., 3 / 4) * pow(x, 3 / 4)


# **********************
def Exercice14_5():
    # **********************
    """ 
    Minimal drag body in hypersonic flow
    """
    # 
    set_title("Minimal drag body in hypersonic flow")
    plot_body = True
    plot_curves = True
    r0, x0 = .01, 0.01
    p1, p0 = 0.1, np.sqrt(3)
    n = 1001            # number of point on the body upper wall 
    displayTable = True

    r, x, p = np.zeros(n), np.zeros(n), np.zeros(n)
    error = np.zeros(n)
    p = np.linspace(p1, p0, n)
    X, R = fx(p), fr(p)
    C1 = r0 / fr(p0)
    C2 = x0 - C1 * fx(p0)
    x = C1 * X + C2
    r = C1 * R
    r_a = r_approx(x, C1)
    error = np.abs(r_a / r - 1) * 100

    if displayTable:
        print("    p \t \t X \t R \t x \t     r \t\t r_a \t error")
        for k in range(n):
            print("%9.4f " * 7 % (p[k], X[k], R[k], x[k], r[k], r_a[k], error[k]))

    print(" @ x0 : fx = %f, fr=%f" % (fx(p0), fr(p0)))
    print(" @ x1 : fx = %f, fr=%f" % (fx(p1), fr(p1)))
    print("  %f <= x <= %f et %f <= r <= %f" % (x[-1], x[0], r[-1], r[0]))
    print("C2 = %e, C1 = %e" % (C2, C1))
    print("1/C2 = %e, 1/C1 = %e" % (1 / C2, 1 / C1))
    print("relative thickness  = %f " % (0.5 * (r[1] - r[0]) / (x[1] - x[0])))

    if plot_body:
        fig = plt.figure(num=1)
        fig.suptitle('optimal shape', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\ln x$', fontsize=20)
        ax.set_ylabel(r'$\ln r$', fontsize=20)
        ax.grid()
        ax.plot(np.log(x), np.log(r), 'b', linewidth=2, label=r'$r(x)$')
        c = np.log(r[0]) - 3 / 4 * np.log(x[0])
        ax.plot(np.log(x), c + 3 / 4 * np.log(x), 'r--', linewidth=2, label=r'$\frac{3}{4}\ln x + c$')
        ax.legend(loc='upper left')

        fig = plt.figure(num=2, figsize=(8, 3))
        fig.suptitle('optimal shape', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$', fontsize=20)
        ax.set_ylabel(r'$r$', fontsize=20)
        ax.grid()
        ax.plot(x, r, 'b', label=r'$r(x)$', linewidth=2)
        ax.plot(x, -r, 'b', label=r'$r(x)$', linewidth=2)
        # ax.plot(x,r_a,'r',label=r'$r_a(x)$')
        # ax.legend(loc='upper left')
        ax.axis('equal')
        plt.show()

    if plot_curves:
        fig = plt.figure(3)
        fig.suptitle('X(p)  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$', fontsize=20)
        ax.set_ylabel(r"$r'$", fontsize=20)
        ax.grid()
        ax.plot(X, p, 'b', label=r'$X(x)$')
        # ax.legend(loc='upper left')

        fig = plt.figure(4)
        fig.suptitle('R(p)  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$r$', fontsize=20)
        ax.set_ylabel(r"$r'$", fontsize=20)
        ax.grid()
        ax.plot(p, R, 'b', label=r'$r(x)$')
        # ax.legend(loc='upper left')

        fig = plt.figure(5)
        fig.suptitle(r'$x(p)$  ', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$p$', fontsize=20)
        ax.set_ylabel(r"$\ln x$", fontsize=20)
        ax.grid()
        ax.plot(p, np.log(x), 'b', label=r'$r(x)$')
        ax.legend(loc='upper left')

        fig = plt.figure(6)
        fig.suptitle('r(p)', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$p$', fontsize=20)
        ax.set_ylabel(r"$r$", fontsize=20)
        ax.grid()
        ax.plot(p, r, 'b', label=r'$r(x)$')
        # ax.legend(loc='upper left')
        plt.show()



def Exercice14_6():
    """ 
    Drag of blunt body in hypersonic flow
    """
    set_title("Drag of blunt body in hypersonic flow")

    plot_curves = True
    M0 = np.float64(20.)
    Lr_over_Rq = np.float64(5.)
    na = 11
    theta_d = np.deg2rad(np.float64(15.))
    R, Rq = np.float64(1.), np.float64(2.)

    set_question('1a : Kp for axisymmetrical body')

    print("Blunt part : ")
    Sref = np.pi * Rq ** 2
    theta = np.linspace(np.pi / 2, theta_d, na)
    Kp_Na = 2.0 * np.sin(theta) ** 2
    print(" theta (°) \t Kp")
    for k in range(na):
        print(" %5.4f \t %5.4f " % (np.rad2deg(theta[k]), Kp_Na[k]))
    print("Wedge  3D part, tangent hypersonic angle : ")

    """ 
    conical shock wave  in supersonic flow 
    """
    # interpolation from database
    set_title("tangential cone method")
    reference_filepath = os.path.join('Book/Data', 'Mach_20_choc_conique.dat')
    lw = 2
    m = 1
    with open(reference_filepath, 'r') as infile:
        A = np.loadtxt(infile, dtype=float, usecols=(0, 1, 2), unpack=True, skiprows=m)
        theta1 = A[0, :]
        kp_cc = A[2, :]
        sigma1 = A[1, :]
    print(A.shape)
    print(theta1.shape)
    f1 = interp1d(theta1, kp_cc, kind='cubic')
    f2 = interp1d(theta1, sigma1, kind='cubic')
    theta_c = np.linspace(50, np.rad2deg(theta_d), 8)

    print('for M0  : %5.4f' % (M0))
    print(" theta (°) \t  Kp_n \t \t  sigma (°) \t Kp ")
    for the in theta_c:
        print(" %5.3f \t %5.3f \t \t  %5.3f \t %5.3f " % (the, 2.0 * np.sin(np.deg2rad(the)) ** 2, f2(the), f1(the)))

    Kp_Lees, sigma_Lees = Solution_Lees()
    Kp_LEES = Kp_Lees * theta_d ** 2
    print("Lees solution, Kp             : %5.4f" % (Kp_LEES))

    print("detached shock length Delta/R : %5.4f " % (detached_shock_length()))

    set_question('1b : Drag coefficient CD for axi-symmetrical body')

    # blunt part
    CD_emoussee = np.pi / Sref * R ** 2 * (1 - np.sin(theta_d) ** 4)
    # conical part
    Rc = R * np.cos(theta_d)
    CD_cone = 2 * np.pi / Sref * (Rq ** 2 - Rc ** 2) * np.sin(theta_d) ** 2
    # rear  part (culot)
    Kp_q = 0.
    CD_culot = -Kp_q
    CD = CD_emoussee + CD_cone + CD_culot
    CD_ConeCulot = Kp_LEES * np.pi / Sref * (Rq ** 2 - Rc ** 2)
    CD1 = CD_emoussee + CD_ConeCulot
    print("calcul1 : CD = %5.4f\t CD_e = %5.4f\t CD_cone = %5.4f \t CD_culot = %5.4f " %
          (CD, CD_emoussee, CD_cone, CD_culot))
    print("calcul2 : CD = %5.4f\t CD_e = %5.4f\t CD_(cone+base) = %5.4f" %
          (CD1, CD_emoussee, CD_ConeCulot))

    set_question('2a : Kp for 2D plane geometry')
    theta_max_2D = np.deg2rad(45.293184)
    theta_2d = np.linspace(45, np.rad2deg(theta_d), 7)
    print('Angles = ', theta_2d)
    print("Kp (Newton's method")
    print(" theta (°) \t Kp_n \t Kp_dt (approx) \t K0 \t Kp_dt")
    for the in theta_2d:
        tm = np.deg2rad(the)
        print(" %2.2f \t %5.3f  \t %5.3f  \t %+05.3f \t %4.3f" % (the, 2 * np.sin(tm) ** 2,
                                                                  (gamma + 1) * tm ** 2, M0 * tm,
                                                                  2 * tm * Constant_Kc(M0 * tm) / M0))

    set_question('3 : Rear Kp with wake')

    theta_s = -np.arctan(1. / Lr_over_Rq)
    Delta_theta = theta_s - theta_d
    print('Wake deviation angle               : %6.3f' % (np.rad2deg(theta_s)))
    print('deviation angle in the expansion   : %6.3f' % (np.rad2deg(Delta_theta)))
    # Kp1=Kp_LEES
    Kp1 = Kp_Na[-1]
    P1_P0 = 1 + gamma * M0 ** 2 / 2 * Kp1
    print('Kp1                                : %6.4f' % Kp1)
    print('P1/P0                              : %6.4f' % P1_P0)
    # print(p_pi(M0),inverse_p_pi(p_pi(M0)))
    P1_Pi1 = P1_P0 * p_pi(M0) / pi2_pi1(M0)
    print('pi ration through the shock        : %e' % (pi2_pi1(M0)))
    print('p1/pi1                             : %f' % P1_Pi1)
    M1 = inverse_p_pi(P1_Pi1)
    print('M1                                 : %6.4f' % M1)
    print('omega(M1)                          : %6.4f' % (np.rad2deg(omega_super(M1))))
    omega2 = omega_super(M1) + np.abs(Delta_theta)
    M2 = invomega(omega2)
    print('omega(M2)                          : %6.4f' % (np.rad2deg(omega2)))
    print('M2                                 : %e' % M2)
    P2_over_P1 = p_pi(M2) / p_pi(M1)
    print('P2_P1                              : %6.4f' % P2_over_P1)
    print('Kp2                                : %6.4f' % (Kp_hyper(M0, P2_over_P1 * P1_P0)))
    print('P2/P0                              : %e' % (P2_over_P1 * P1_P0))
    print('P0/P1                              : %e' % (1 / P1_P0))

    opt = True
    if plot_curves:
        fig = plt.figure(num=1)
        fig.suptitle('Kp on the blunt body', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$', fontsize=20)
        ax.set_ylabel(r'$Kp$', fontsize=20)
        ax.grid()
        na = 75
        Lc = (Rq - Rc) / np.tan(theta_d)
        print('Lc = ', Lc)
        theta = np.linspace(np.pi / 2, theta_d, na)
        Kp1 = 2.0 * np.sin(theta) ** 2
        x1 = (1 - np.sin(theta)) * R

        tc = np.deg2rad(theta_c)
        x2 = (1 - np.sin(tc)) * R
        Kp2 = f1(theta_c)

        tm = np.deg2rad(theta_2d)
        x3 = (1 - np.sin(tm)) * R
        Kp3 = 2 * tm * Constant_Kc(M0 * tm) / M0

        if opt:
            Kp1 = np.append(Kp1, Kp1[-1])
            x1 = np.append(x1, x1[-1] + Lc)
            Kp2 = np.append(Kp2, Kp2[-1])
            x2 = np.append(x2, x2[-1] + Lc)
            Kp3 = np.append(Kp3, Kp3[-1])
            x3 = np.append(x3, x3[-1] + Lc)

        ax.plot(x1, Kp1, 'b', linewidth=2, label=r'$Kp_n$')
        ax.plot(x2, Kp2, 'r', linewidth=2, label=r'$Kp_{ct}$')
        ax.plot(x3, Kp3, 'k', linewidth=2, label=r'$Kp_{dt}$')
        ax.legend(loc='upper right')
        ax.axis([0, 1., 0, 2])
        print("x_init : cone = %f, wedge = %f" % (x2[0], x3[0]))
        print("chord= ", x1[-1])

        # ax.plot(np.log(x),c+3/4*np.log(x),'r--',linewidth=2,label=r'$\frac{3}{4}\ln x + c$')
        plt.savefig("fig14-exo6_solution.svg")

        plt.show()


def Exercice14_100():
    """
    Conical shock angle curves
    """
    set_title("Conical shock angle curves")
    """
    Downstream Mach , wall pressure coefficient Kp, shock angel sigma as a function of the cone angle theta_c
    """
    main_multiple_Mach((5, 10, 20, 50), ("sigma", "downstream_Mach", "Kp"), Lees=True,
                       t=(20, 2, 1), s=(26, 5, 1), Kp=(30, 5, 1), Mav=(50, 5, 1))


def Exercice14_102():
    """
    Test : inverted omega(Mach) function  in supersonic flow
    """
    set_title("Test : inverted omega(Mach) function  in supersonic flow")
    M0 = 3
    theta = omega_super(M0)
    print("for M0=%f we get theta = %f °" % (M0, np.rad2deg(theta)))
    print(invomega(theta, show=True))

    print(valeur_sigma(M0, 10, sigma=0, show=True))


def Exercice14_101(cas=0):
    """
    Test: Prandtl-Meyer approximated function in   hypersonic flow
    cas = 0 : several M0, cas = 1 , single M0, cas = 2: omega (P.M.) function
    """
    n = 101
    M_init, M_end = 10, 15
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$M$', fontsize=20)
    ax.grid()
    if cas == 0:
        M0 = np.array([3.5, 5, 10])
        m = len(M0)
        error = np.zeros([m, n, 2])
        M = np.zeros([m, n])
        for k in range(m):
            M[k, :] = np.linspace(M0[k], M_end, n)
            error[k, :, 0], error[k, :, 1] = approx_omega(M[k, :], M0[k])
        fig.suptitle(r'$\theta(M)$', fontsize=14, fontweight='bold')
        ax.set_ylabel(r'$\theta(^\circ)$', fontsize=20)
        for k in range(m):
            ax.plot(M[k, :], error[k, :, 0], '-', label='exa., M0=%4.2f' % (M0[k]))
            ax.plot(M[k, :], error[k, :, 1], '-', label='app., M0=%4.2f' % (M0[k]))
            ax.legend(loc='lower left')
    elif cas == 1:
        Mach = np.linspace(M_init, M_end, n)
        fig.suptitle(r'$\theta(M), M0 $= %3.1f' % M_init, fontsize=14, fontweight='bold')
        y1, y2 = approx_omega(Mach, M_init)
        # ax.axis([M_init, M_end, -20, 0])
        ax.set_ylabel(r'$\theta(^\circ)$', fontsize=20)
        ax.plot(Mach, y1, label='Exact')
        ax.plot(Mach, y2, label='Approx')
        ax.legend(loc='upper right')
    else:
        Mach = np.linspace(M_init, M_end, n)
        fig.suptitle('errors M0 = %3.1f' % M_init, fontsize=14, fontweight='bold')
        # ax.axis([M_init, M_end, 100,140])
        # ax.plot(Mach,np.rad2deg(omega_super(Mach)),label=r'$\omega$ exact')
        # ax.plot(Mach,np.rad2deg(omega_hyper(Mach)),label=r'$\omega$ approx.')
        ax.plot(Mach, (omega_hyper(Mach) - np.rad2deg(omega_super(Mach))), label=r'$\Delta\omega(^\circ)$')
        ax.legend(loc='lower right')
    plt.show()
