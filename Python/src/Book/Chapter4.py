#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2017 - 2020*

*license : AGPL-3.0*

Correction of the chapter 4 exercises
    ..
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import os
from Tools.misc import set_title, set_question
from scipy.optimize import newton
from Book.ConformalMapping.mappings import set_parameters, AirfoilConfMap
from IncompressibleFlow.Airfoil import slats_flaps_effects, naca
from IncompressibleFlow.LinearTheoryProfile import set_parameters_airfoil, LinearTheoryAirfoil
from Tools.misc import set_title, set_question, set_info, set_alert

def Exercice4_1():
    """
    Flate plate with incidence
    """
    set_title("Flate plate with incidence")
    par = set_parameters()
    s = AirfoilConfMap(par)
    s.map = "Joukowski1"
    s.alpha = 20
    s.opt_levels = "manuel"
    s.levels_iso = [-1, 1, 0.2]
    s.adimChord = True
    s.n_circle = 101
    s.plot_velocity = False
    s.plot_circle = False
    s.plot_airfoil = False
    s.plot_Kp = False
    s.run_airfoil()
    plt.show()

def plot_profil(xe, xi, ye, yi, cas):
    """
    plot of a Kp
    """
    plt.figure(figsize=(16, 4))
    plt.title("Airfoil comparisons")
    for k in range(len(xe)):
        plt.plot(xe[k], ye[k], label="ye: " + cas[k])
        plt.plot(xi[k], yi[k], label="yi: " + cas[k])
    plt.grid()
    # plt.xlim(-2.05,2.05)

    plt.axis("equal")
    plt.legend(loc="best")



def Exercice4_2():
    """ 
    Joukowski's airfoil
    """
    set_title(" Joukowski's airfoil")

    xe, xi = [], []
    ye, yi = [], []
    cas = ["Thick airfoil", "Camber line"]  # thick airfoil, cambered line
    par = set_parameters()
    set_question("a - thick airfoil")
    s = AirfoilConfMap(par)
    # zero incidence plot of the airfoil
    # pameter:  Zc
    s.map = "Joukowski1"
    s.a = 1
    s.alpha = 0
    s.Zc = -0.1 + 1j * 0
    s.plot_circle = False
    s.plot_airfoil = False
    s.opt_levels = "manuel"
    s.levels_iso = [-1, 1, 0.2]
    s.adimChord = False
    s.plot_velocity = False
    s.plot_Kp = True
    s.plot_psi = False
    s.a = 1
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    ye.append(s.ye)
    yi.append(s.yi)

    set_question("3 - Cambered line")

    s.Zc = 0 + 1j * 0.1
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    ye.append(s.ye)
    yi.append(s.yi)

    plot_profil(xe, xi, ye, yi, cas)
    plt.show()

def Exercice4_3():
    """ 
    pq family airfoils 
    """
    set_title("pq family airfoils")

    set_question("1 à 4:  airfoil with zero incidence")
    # cas= 'analytical'
    par = set_parameters_airfoil()
    s = LinearTheoryAirfoil(par)
    s.airfoil = ['family_pq', "analytical", "ys"]
    s.nFourier = 20
    s.run_theory2D()

    set_question("5 à 7:  airfoil at adaptation incidence")
    s.alpha = 2.127744856489093
    s.run_theory2D()
    print("verification p x epsilon = ", s.AirfoilParam[0] * s.AirfoilParam[1])
    print("Kp max, Kp min           = ", np.amax(s.Kp_s[2:40]), np.amin(s.Kp_s))
    plt.show()

def Exercice4_4():
    """ 
    Double cambered airfoil +  slats and flaps
    """
    set_title("Double cambered airfoil +  slats and flaps")

    set_question("1- Airfoil")
    # cas= 'analytical'
    par = set_parameters_airfoil()
    s = LinearTheoryAirfoil(par)
    s.airfoil = ['double-cambered', "analytical", "ys"]
    s.run_theory2D()

    set_question("3- slats and flaps: flaps")

    # The angle is  1 rad, since the problem is linear I can find the coefficient in front of beta_v 

    DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a = slats_flaps_effects(Flaps=True, theta=0, length=0.25,
                                                                        angle=180 / np.pi)

    set_question("4- slats and flaps: flap, variation of incidence")
    delta_alpha = DeltaCL / (2 * np.pi)
    print("variation of incidence  delta_alpha / beta_v = ", delta_alpha)

    set_question("3- slats and flaps: slats")

    # The angle is  1 rad, since the problem is linear I can find the coefficient in front of beta_v 
    DeltaCL, DeltaCmF, DeltaAlpha_0, DeltaAlpha_a = slats_flaps_effects(Flaps=False, theta=0, length=0.25,
                                                                        angle=180 / np.pi)
    plt.show()

def Exercice4_5():
    """
    Camber and thickness laws
    """
    set_title("Camber and thickness laws")
    print("Analytical and Matlab solutions only")


def Exercice4_6():
    """ 
    Slots and flaps on a rectangular wing
    """
    set_title("Slots and flaps on a rectangular wing")
    validation = False  # to validate the 'slats_flaps_effects' function
    AileRectangle = True  # exercise with the  rectangular wing
    rho = 1.2  # air density  kg/m^3
    U0 = 50  # velocity in m/s
    S = 30  # reference surface in  m^2
    if validation:
        slats_flaps_effects(Flaps=True, theta=53, length=0)
        slats_flaps_effects(Flaps=True, theta=0, length=0.2, angle=10)
        slats_flaps_effects(Flaps=False, theta=143)
        slats_flaps_effects(Flaps=False, length=0.2, angle=10)

    if AileRectangle:
        dCL, dCmF, dAlpha_0, dAlpha_a = slats_flaps_effects(Flaps=True, theta=0, length=0.2, angle=20)
        DeltaL = 1 / 2 * rho * U0 ** 2 * S * dCL
        print("Delta L                           : %e N" % (DeltaL))
        dCL, dCmF, dAlpha_0, dAlpha_a = slats_flaps_effects(Flaps=False, length=0.05, angle=15)
        DeltaL = 1 / 2 * rho * U0 ** 2 * S * dCL
        print("Delta L                           : %e N" % (DeltaL))

    # curves for the plot:

    # slot :
    npt = 301
    lv = np.linspace(0, 0.3, npt)
    lb = np.linspace(0, 0.15, npt)
    angle = 1
    dCLv, dCmFv, dAlpha_0v, dAlpha_av = slats_flaps_effects(Flaps=True, theta=0, length=lv, angle=angle)
    dCLb, dCmFb, dAlpha_0b, dAlpha_ab = slats_flaps_effects(Flaps=False, theta=0, length=lb, angle=angle)
    c = 100

    plt.figure(figsize=(10, 8.5))
    plt.title(r"$C_L$ et $Cm_F$")
    plt.plot(c * lv, dCLv, label=r"flap:$\Delta C_L$")
    plt.plot(c * lb, dCLb, label=r"slat  :$\Delta C_L$")
    plt.plot(c * lv, dCmFv, label=r"flap:$\Delta Cm_F$")
    plt.plot(c * lb, dCmFb, label=r"slat  :$\Delta Cm_F$")
    plt.xlabel(r"$l/\ell$")
    plt.grid()
    # plt.xlim(-2.05,2.05)
    plt.legend(loc="best")

    plt.figure(figsize=(10, 8.5))
    plt.title(r"$\alpha_0$ et $\alpha_a$")
    plt.plot(c * lv, np.rad2deg(dAlpha_0v), label=r"flap:$\Delta\alpha_0$")
    plt.plot(c * lb, np.rad2deg(dAlpha_0b), label=r"slat  :$\Delta\alpha_0$")
    plt.plot(c * lv, np.rad2deg(dAlpha_av), label=r"flap:$\Delta\alpha_a$")
    plt.plot(c * lb, np.rad2deg(dAlpha_ab), label=r"slat  :$\Delta\alpha_a$")
    plt.xlabel(r"$l/\ell$")
    plt.grid()
    # plt.xlim(-2.05,2.05)
    plt.legend(loc="best")

    plt.show()


def Exercice4_7():
    """
    Linearized airfoil theory : NACA airfoil
    """
    plot_option = False
    
    set_title("NACA airfoil")
    cas = 'analytical'         # "analytical"  # or 'numerical'  for the camberline
    par = set_parameters_airfoil()
    s = LinearTheoryAirfoil(par)
    s.AirfoilParam = ['4412']

    nc = 601
    if cas == "numerical":

        # test fo  NACA airfoil, 
        s.alpha = 0
        s.airfoil = ["naca4", "numerical", "y+-"]
        s.pnt = nc

        set_question("1- NACA  airfoil " + s.AirfoilParam[0])

        # NACA airfoil geometry
        # the chord is slightly greater than  1 for a cambered airfoil
        # so we have to dimensionalize the geometry again and to do make a change of coordinates   -1/2 <= x/c <= 1/2
        # It creates a numerical error which is not always negligible
        # It is better to use the analytical solution below

        xNaca, yNaca, thetaNaca = naca(s.AirfoilParam[0], nc, half_cosine_spacing=True)
        chord = abs(np.amin(xNaca) - np.amax(xNaca))
        print("original chord of NACA airfoil %s  : %f" % (s.AirfoilParam, chord))
        # Chord correction if it is taken into account (not absolutely necessary) 
        xNaca, yNaca = (xNaca - np.amin(xNaca)) / chord, yNaca / chord
        xNaca -= 0.5
        indLE, indTE = np.argmin(xNaca), np.argmax(xNaca)
        print("x_LE = x[%i] = %f , \t x_TE = x[%i] = %f" % (indLE, xNaca[indLE], indTE, xNaca[indTE]))
        s.xe, s.ye = xNaca[indTE:indLE + 1], yNaca[indTE:indLE + 1]
        s.xi, s.yi = xNaca[indLE:], yNaca[indLE:]
        print('\t upper wall : x = ', s.xe[:2], '...', s.xe[-3:])
        print('\t            y = ', s.ye[:2], '..', s.ye[-3:])
        print('\t lower wall : x = ', s.xi[:2], '...', s.xi[-3:])
        print('\t            y = ', s.yi[:2], '...', s.yi[-3:])
    elif cas == "analytical":
        s.airfoil = ["naca4", "analytical", "ys"]
        s.AirfoilParam = ['4412', ]
        s.alpha = 5
        s.nFourier = 5

    else:
        return

    set_question("2- Linear airfoil theory")
    
    set_alert("alpha = %2.1f" % s.alpha)
    
    s.plot_Kp = True
    s.plot_airfoil = True
    s.run_theory2D()
    # s.PlotDelta_s()

    set_alert("alpha = 0 " )
    par = set_parameters_airfoil()
    s1 = LinearTheoryAirfoil(par)

    s1.airfoil = ["naca4", "analytical", "ys"]
    s1.AirfoilParam = ['4412', ]
    s1.alpha = 0
    s1.nFourier = 5

    s1.run_theory2D()

    set_alert("alpha = alpha_adaptation " )
    par = set_parameters_airfoil()
    s2 = LinearTheoryAirfoil(par)
    s2.airfoil = ["naca4", "analytical", "ys"]
    s2.AirfoilParam = ['4412', ]
    s2.alpha = 0.514  # adaptation incidence in degrees.
    s2.nFourier = 5
    s2.pnt = nc
    s2.run_theory2D()

    if plot_option:
        fig = plt.figure()
        fig.suptitle('Cambered line comparison', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\eta$', fontsize=20)
        ax.grid()
        ax.plot(s.x, np.rad2deg(s.delta_s), '-', linewidth=1, color='black', label="numerical")
        ax.plot(s1.x, np.rad2deg(s1.delta_s), '-', linewidth=1, color='red', label="analytical")
        plt.legend(loc="best")
        plt.show()

        fig = plt.figure(figsize=(16, 4))
        fig.suptitle('Airfoil comparison', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\eta$', fontsize=20)
        ax.grid()
        ax.plot(s.xe, s.ye, '-', linewidth=1, color='blue')
        ax.plot(s.xi, s.yi, '-', linewidth=1, color='blue')
        ax.plot(s.x, s.y_s, '-', linewidth=1, color='black')
        ax.plot(s1.x, s1.y_s, '-', linewidth=1, color='red')
        ax.axis("equal")

        fig = plt.figure()
        fig.suptitle('Kp comparison', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$\eta$', fontsize=20)
        ax.grid()
        ax.plot(s.x[0:-1], s.Kp_s[0:-1], '-', linewidth=1, label=s.AirfoilParam[0] + "numerical")
        ax.plot(s1.x[0:-1], s1.Kp_s[0:-1], '-', linewidth=1, label=s1.AirfoilParam[0] + "analytical")
        ax.plot(s2.x[0:-1], s2.Kp_s[0:-1], '-', linewidth=1, label=s2.AirfoilParam[0] + "analytical")
        plt.legend(loc="best")
        plt.show()


def Exercice4_8():
    """
    Linearized airfoil theory : symmetrical NACA airfoil 
    """
    set_title("NACA airfoil")

    par = set_parameters_airfoil()
    s = LinearTheoryAirfoil(par)
    s.AirfoilParam = ['2415', ]
    nc = 601
    s.alpha = 5
    s.airfoil = ["naca4", "analytical", "ys"]
    s.nFourier = 3
    s.pnt = nc
    s.plot_Kp = True
    s.plot_airfoil = True
    s.run_theory2D()
