#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 11 exercises 
    ..
"""

import numpy as np
import matplotlib.pyplot as plt
from CompressibleFlow.fonctions import gamma, P2_P1, oblique_shock, theta_from_Mn0, isentropic_expansion
from CompressibleFlow.tables import Atmosphere
from Tools.misc import set_question, set_title

width, height = 10, 8


def Exercice11_1():
    """ 
    plots of  p1/p0 with respect to flow deviation angle for different Mach number
    """
    ifig = 0
    set_title("plots of  p1/p0 with respect to flow deviation angle for different Mach number")
    n = 101
    Mach_list = np.array([1.3, 1.6, 2.0, 2.3, 2.5, 2.7, 3.0])

    set_question('1 : Compression by shock ')

    print('Mach  = ', Mach_list)

    theta, r_p = [], []
    for Mach in Mach_list:
        sigma = np.linspace(np.arcsin(1 / Mach), np.pi / 2, n)
        theta.append(theta_from_Mn0(sigma, Mach))
        r_p.append(P2_P1(Mach * np.sin(sigma)))

    fig = plt.figure(ifig, figsize=(width, height))
    fig.suptitle('Compression by shock', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    # fig.subplots_adjust(top=0.80)
    ax.set_ylabel(r'$P_2/P_1$', fontsize=20)
    ax.set_xlabel(r'$\theta(^\circ)$', fontsize=20)
    for k, Mach in enumerate(Mach_list):
        ax.plot(np.rad2deg(theta[k]), r_p[k], linewidth=2, label=r'$M_0=$%6.5f' % Mach)
    ax.axis([0, 35, 1, 11])
    xmajor_ticks, xminor_ticks = np.arange(0, 35, 10), np.arange(0, 35, 2)
    ymajor_ticks, yminor_ticks = np.arange(1, 12, 2), np.arange(1, 12, 0.5)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.tick_params(labelsize=16)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.grid(which='both')
    ax.legend(loc='lower right')

    set_question('2 : Isentropic expansion')

    theta = np.linspace(0, np.pi / 4, n)
    r_p = []
    for Mach in Mach_list:
        r_p.append(isentropic_expansion(theta, Mach, gamma=gamma))
    ifig += 1
    fig = plt.figure(ifig, figsize=(width, height))
    fig.suptitle('Isentropic expansion', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    # fig.subplots_adjust(top=0.80)
    ax.set_ylabel(r'$P_2/P_1$', fontsize=20)
    ax.set_xlabel(r'$\theta(^\circ)$', fontsize=20)
    for k, Mach in enumerate(Mach_list):
        ax.plot(np.rad2deg(theta), r_p[k], linewidth=2, label=r'$M_0=$%6.5f' % Mach)
    ax.axis([0, 45, 0, 1])
    xmajor_ticks, xminor_ticks = np.arange(0, 45, 10), np.arange(0, 45, 2)
    ymajor_ticks, yminor_ticks = np.arange(0, 1.1, 0.25), np.arange(0, 1, 0.05)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.tick_params(labelsize=16)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.grid(which='both')
    ax.legend(loc='upper right')
    plt.show()


def Exercice11_2():
    """
    Oblique shock wave reflection of a channel wall
    """
    M0 = 1.5  # upstream of the two shocks
    theta0 = 5  #  angle in degrees of the first wall deviation
    theta1 = 10  # angle in degrees of the second wall deviation
    ifig = 0
    n = 51
    plot_figure = True
    plot_grid = False

    thetaMin, thetaMax = -0.5, 13  # in degrees
    dthetaMaj, dthetaMin = 5, 1
    rMax = 2.5

    sigma0, M_downstream0, P_downstream0, omg_up, omg_down = oblique_shock(M0, theta0)
    sigma1, M_downstream1, P_downstream1, omg_up, omg_down = oblique_shock(M0, theta1)

    Mach_list = np.array([M0, M_downstream0, M_downstream1])

    # calculus of the heart-shape curves
    xc, yc = [], []
    theta, r_p = [], []
    for Mach in Mach_list:
        sigma = np.linspace(np.arcsin(1 / Mach), np.pi / 2, n)
        t = theta_from_Mn0(sigma, Mach)
        r = P2_P1(Mach * np.sin(sigma))
        theta.append(np.concatenate((-np.flipud(t), t)))
        r_p.append(np.concatenate((np.flipud(r), r)))

    # center of these curves 
    xc.append(0)
    yc.append(1)        # first shock : so y-shift
    xc.append(theta0)
    yc.append(P_downstream0)   # second shock shift, case at 5°
    xc.append(theta1)
    yc.append(P_downstream1)   # second shock shift, case at 10°

    if plot_figure:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('shock interaction', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        # fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$P_2/P_1$', fontsize=20)
        ax.set_xlabel(r'$\theta(^\circ)$', fontsize=20)
        for k, Mach in enumerate(Mach_list):
            ax.plot(np.rad2deg(theta[k]) + xc[k], r_p[k] * yc[k], linewidth=2, label=r'$M_0=$%6.5f' % Mach)
        ax.axis([0, thetaMax, 1, rMax])
        if plot_grid:
            xmajor_ticks, xminor_ticks = np.arange(thetaMin, thetaMax, dthetaMaj), np.arange(thetaMin, thetaMax,
                                                                                             dthetaMin)
            ymajor_ticks, yminor_ticks = np.arange(1, rMax, 0.5), np.arange(1, rMax, 0.1)
            ax.set_xticks(xmajor_ticks)
            ax.set_xticks(xminor_ticks, minor=True)
            ax.tick_params(labelsize=16)
            ax.set_yticks(ymajor_ticks)
            ax.set_yticks(yminor_ticks, minor=True)
            ax.grid(which='both')
        else:
            ax.grid()
        ax.legend(loc='lower right')
        plt.show()


def Exercice11_3():
    """
    Interaction of 2  oblique shock waves
    """
    M0 = 2.0        # upstream Mach number
    theta0 = 5      # deviation angle in degrees
    theta1 = 15     # deviation angle in degrees
    thetag = 0      # angle of the pressure line in degrees
    ifig = 0
    n = 51
    plot_figure = True
    plot_grid = True

    thetaMin, thetaMax = -0., 16
    dthetaMaj, dthetaMin = 5, 1
    rMax = 2.4

    sigma0, M_downstream0, P_downstream0, omg_up, omg_down = oblique_shock(M0, theta0, show=False, msg='Zone 1')
    sigma1, M_downstream1, P_downstream1, omg_up, omg_down = oblique_shock(M_downstream0, theta1 - theta0,
                                                                           show=False, msg='Zone 2')
    sigma2, M_downstream2, P_downstream2, omg_up, omg_down = oblique_shock(M0, theta1 - thetag, show=False, msg='Zone 4')

    Mach_list = np.array([M0, M_downstream0])
    # Mach_list=np.array([M0])
    # Mach_list=np.array([M0])

    # calculus of heart-shape curve
    xc, yc = [], []
    theta, r_p = [], []
    for Mach in Mach_list:
        sigma = np.linspace(np.arcsin(1 / Mach), np.pi / 2, n)
        t = theta_from_Mn0(sigma, Mach)
        r = P2_P1(Mach * np.sin(sigma))
        theta.append(np.concatenate((-np.flipud(t), t)))
        r_p.append(np.concatenate((np.flipud(r), r)))

    # curve  center 
    xc.append(0)
    yc.append(1)                # first shock : no shift
    xc.append(theta0)
    yc.append(P_downstream0)    # Shift of the second shock, case at 5°
    # xc.append(theta1)
    # yc.append(P_downstream1)  # Shift of the second shock, case at 10°

    if plot_figure:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('Interaction of 2  oblique shock waves', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        # fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$P_2/P_1$', fontsize=20)
        ax.set_xlabel(r'$\theta(^\circ)$', fontsize=20)
        for k, Mach in enumerate(Mach_list):
            ax.plot(np.rad2deg(theta[k]) + xc[k], r_p[k] * yc[k], linewidth=2, label=r'$M_0=$%6.5f' % (Mach))
        ax.axis([0, thetaMax, 1, rMax])
        if plot_grid:
            xmajor_ticks, xminor_ticks = np.arange(thetaMin, thetaMax, dthetaMaj), np.arange(thetaMin, thetaMax,
                                                                                             dthetaMin)
            ymajor_ticks, yminor_ticks = np.arange(1, rMax, 0.5), np.arange(1, rMax, 0.1)
            ax.set_xticks(xmajor_ticks)
            ax.set_xticks(xminor_ticks, minor=True)
            ax.tick_params(labelsize=16)
            ax.set_yticks(ymajor_ticks)
            ax.set_yticks(yminor_ticks, minor=True)
            ax.grid(which='both')
        else:
            ax.grid()
        ax.legend(loc='lower right')
        plt.show()
