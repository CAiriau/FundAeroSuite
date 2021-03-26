#!/bin/python
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Module : incompressibleFlow.LinearTheoryProfile
    ..

Python module for conformal mapping for airfoil, Joukowski transform
mappings file, with the class
* set_parameters_airfoil
* LinearTheoryAirfoil

Main program for  conformal
"""

import os
import numpy                as np
import numpy.ma             as ma
import matplotlib.pyplot    as plt
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from Tools.misc import set_title, set_question, set_info, set_alert

fz = 16  # legend font size (option)


def set_parameters_airfoil():
    """
    Set default parameters
    """
    prm = {}

    # physical parameters
    prm["probleme"] = "portant"         # "portant" "portant+épais"  ["camber" or "camber+thickness"]
    prm["Uinf"] = 1.                    # upstream velocity in m/s
    prm["rho"] = 1.3                    # density in kg/m^3
    prm["alpha"] = 0.0                  # angle of attack in degrees
    prm["Umax"] = 10                    # velocity filter for singular points
    prm["KpLim"] = [-3, 1]              # y-axis range for Kp plot

    prm["airfoil"] = ["parabole", "analytical",
                      "h(x)"]           # case name  : "parabole", "family_pq" "double-cambered"
    #                                   airfoil[0]: "parabole" or "double sharp" or "naca4"
    #
    prm["npt"] = 201                    # number of points to define the camberline or the airfoil upper wall
    prm["eps"] = 0.1                    # small parameter for airfoil definition
    prm["chord"] = 1.                   # airfoil chord
    prm["nFourier"] = 5                 # number of coefficients to troncate the Fourier's serie
    # plot parameters
    prm["plot_squelette"] = True        # True to plot camberline
    prm["plot_Kp"] = True               # True to plot pressure coefficient Kp on the airfoil
    prm["plot_airfoil"] = True          # True to plot the airfoil
    prm["plot_velocity"] = False        # True to plot U/Uinf on the airfoil
    prm["adimChord"] = True             # True to get a nondimensional airfoil w.r.t. chord
    prm["order"] = 2                    # order of the integration scheme (2 or 4)
    prm["xsize"] = 10                   # figure width
    prm["ysize"] = 8.5                  # figure length
    eps, q, CL = 0.1, 7 / 8., 0.35
    # prm["AirfoilParam"] = [eps, -CL/(np.pi*eps*(3./4-q)) ,q] #[epsilon, p, q]
    p = CL * 8 / (np.pi * 3 * eps)
    prm["AirfoilParam"] = [eps, p, q]

    return prm

class LinearTheoryAirfoil(object):
    """
    airfoil designed by a conformal mapping
    """
    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#', 60 * '*')
        print('# %s' % ("Linear Theory for 2D section"))
        print('#', 60 * '*', '\n')

        self.probleme = prm["probleme"]  # "portant" "portant+épais",  ["camber" or "camber+thickness"]
        # physical parameters

        self.Uinf = prm["Uinf"]         # upstream velocity in m/s
        self.rho = prm["rho"]           # density in kg/m^3
        self.alpha = prm["alpha"]       # angle of attack in degrees (AoA)
        self.Gamma = 0.                 # circulation to calculate
        self.Kp = 0.                    # pressure coefficient
        self.Kp_s = 0.                  # pressure coefficient due to camberline
        self.chord = prm["chord"]       # airfoil chord

        self.camberMax = 0.             # maximal relative camberline
        self.thicknessMax = 0.          # maximal relative thickness
        self.alpha_0 = 0.               # AoA at zero lift
        self.alpha_ad = 0.              # AoA for adapted airfoil
        self.CL = 0.                    # lift coefficient
        self.xF = 0.                    # aerodynamic center x-coordinate / chord
        self.CmF = 0.                   # pitching momentum coefficient at the aerodynamic center
        self.Cm_LE = 0.                  # pitching momentum coefficient at the leading edge
        self.delta_s = 0.               # camberline slope law
        self.y_s = 0.                   # camberline equation
        self.y_e = 0                    # thickness law
        self.xe = [0.]                  # upper wall x-coordinate if airfoil coordinate given  from numerical data
        self.ye = [0.]                  # upper wall y-coordinate if airfoil coordinate given  from numerical data
        self.xi = 0.                    # lower wall x-coordinate if airfoil coordinate given  from numerical data
        self.yi = [0.]                  # lower wall y-coordinate if airfoil coordinate given  from numerical data
        self.Lref = 1.                  # reference length of the x-coordinates in some plots

        self.AirfoilParam = prm["AirfoilParam"]  # [epsilon, p, q]

        self.Umax = prm["Umax"]         # velocity filter for singular points
        self.KpLim = prm["KpLim"]       # y-axis range for Kp plot
        self.airfoil = prm["airfoil"]   # [testcase type, airfoil name, type of data]
                                        # testcase:"analytical" or "numerical"
                                        # type of date : "h(x)"  or "y+-" or ""
        self.npt = prm["npt"]           # number of points to define the camberline or the airfoil upper wall
        self.eps = prm["eps"]           # small parameter for airfoil definition
        self.nFourier = prm["nFourier"] # number of coefficients to troncate the Fourier's serie
        self.order = prm["order"]       # order of the integration scheme (2 or 4)

        # plot parameters
        self.plot_squelette = prm["plot_squelette"]  #  True to plot camberline

        self.plot_Kp = prm["plot_Kp"]  # True to plot pressure coefficient Kp on the airfoil
        self.plot_airfoil = prm["plot_airfoil"]     # True to plot the airfoil
        self.plot_velocity = prm["plot_velocity"]   # True to plot U/Uinf on the airfoil
        self.xsize = prm["xsize"]       # figure width
        self.ysize = prm["ysize"]       # figure length
        self.adimChord = prm["adimChord"]  # True to get a nondimensional airfoil w.r.t. chord
        self.xlabel = "x"
        print("initialization done")

    def run_theory2D(self):
        """
        application of the linearized aifoil theory to a given airfoil
        """
        if self.adimChord:
            self.Lref = self.chord
            self.xlabel = r"$x/c$"

        self.theta = np.linspace(0, np.pi, self.npt)
        self.x = self.chord / 2 * np.cos(self.theta)

        if self.airfoil[0] == "naca4":
            self.naca4_analytical()
        self.solution_analytical()
        
        set_info("Numerical solution : \n")

        if self.airfoil[2] == "y+-":
            self.airfoil_geometry()
            self.order = 2

        if self.probleme == "portant":
            if self.airfoil[1] != "numerical":
                self.delta_s = self.camber(self.x)
                self.camberline()
            self.SerieFourier()
            for i in range(self.nFourier):
                print("a[%i]                                     : %f" % (i, self.a[i]))
            self.Aerodynamics()
            self.Kp_portant()
        self.graphiques()

    def naca4_analytical(self):
        """
        Initialization of the  NACA 4 parameters
        """
        number = self.AirfoilParam[0]
        m = float(number[0]) / 100.0
        p = float(number[1]) / 10.0
        t = float(number[2:]) / 100.0
        self.AirfoilParam[1:] = [m, p, t]
        print("NACA 4 parameters                         : m = %f, p = %f, t = %f" % (m, p, t))

    def solution_analytical(self):
        """
        analytical solution if it is implemented
        """
        print(self.airfoil[0],  ": analytical solution")
        if self.airfoil[0] == 'family_pq':
            [eps, p, q] = self.AirfoilParam
            print("for validation, analytical results :")
            print("eps = %f, p = %f, q= %f" % (eps, p, q))
            print("a0 = %f, a1 = %f, a2 = %f Cl = %f alpha_0 = %f° " % (
            p * eps / 2, p * eps * (1 - 2 * q), 3 * p * eps / 4, -np.pi * p * eps * (3 / 4 - q),
            np.rad2deg(p * eps / 2 * (3 / 4 - q))))
        elif self.airfoil[0] == 'double-cambered':
            print("for validation, analytical results :")
            print("a0 = %f, a1 = %f, a2 = %f Cl = %f alpha_0 = %f° " % (
            -2 * self.eps, 0, -3 * self.eps, np.pi * self.eps, np.rad2deg(-self.eps / 2)))
        elif self.airfoil[0] == 'naca4':
            [m, p] = [self.AirfoilParam[1], self.AirfoilParam[2]]
            if m != 0:
                theta_c = np.arccos(2 * p - 1)
                a_0 = 4 * m / (np.pi * p ** 2 * (1 - p) ** 2) * (1 - 2 * p) * (
                            2 * np.sqrt(p * (1 - p)) + (1 - 2 * p) * theta_c) - 4 * m / p ** 2 * (
                                  1 - 2 * p) 
                a_0p = a_0 - 4 * np.deg2rad(self.alpha)
                a_1 = 2 * m / (np.pi * p ** 2 * (1 - p) ** 2) * (1 - 2 * p) * (
                            2 * (1 - 2 * p) * np.sqrt(p * (1 - p)) + theta_c) - 2 * m / p ** 2
                a_2 = 32 * m * (1 - 2 * p) * np.sqrt(p * (1 - p)) / (3 * p * np.pi * (1 - p))
                print("a0 (alpha = 0) = %f, \t a0 (alpha = %2.1f) = %f" % (a_0,  self.alpha, a_0p))
                print("a1 = %f, \t a2 = %f " % (a_1, a_2))
                print("Cl = %f,  \t alpha_0 = %2.3f° \n" % (-np.pi / 2 * (a_0p + a_1), np.rad2deg(a_0 + a_1) / 4))

    def thicknessLaw(self):
        """
        thickness law, just for "double-cambered", additional implementation required
        """
        if self.airfoil[0] == "double-cambered":
            eta = self.x / self.chord
            mu = self.AirfoilParam[0]
            self.y_e = mu * self.eps * self.chord / 2 * (1 - eta) * np.sqrt(1 - eta ** 2)

        else:
            self.y_e = np.zeros(len(self.x))

    def thicknessSlope(self, x):
        """
        slope of the thickness law, just for "double-cambered", additional implementation required
        """
        if self.airfoil[0] == "double-cambered":
            eta = self.x / self.chord
            mu = self.AirfoilParam[0]
            return -mu * self.eps * (1 + 2 * eta) * np.sqrt((1 - eta) / (1 + eta))
        else:
            np.zeros(len(self.x))

    def camberline(self):
        """
        cambered airfoil with analytical equations
        """
        if self.airfoil[0] == 'parabole':
            a = self.chord / 4
            self.y_s = self.eps / a ** 2 / 2 * (2 * a - self.x) * (2 * a + self.x)
        elif self.airfoil[0] == "family_pq":
            [eps, p, q] = self.AirfoilParam
            eta = self.x / self.chord
            self.y_s = eps * p * self.chord * (eta + 1 / 2) * (1 / 2 - eta) * (q - eta - 1 / 2)
        elif self.airfoil[0] == "double-cambered":
            eta = self.x / self.chord
            self.y_s = self.eps * self.chord / 2 * eta * (1 - eta ** 2)
        elif self.airfoil[0] == "naca4":
            self.y_s = np.zeros(self.npt)
            [m, p, t] = self.AirfoilParam[1:]
            if p * m != 0:
                for k in range(self.npt):
                    x = self.x[k] + 0.5 * self.chord
                    if x / self.chord <= p:
                        self.y_s[k] = m / p ** 2 * x * (2 * p - x / self.chord)
                    else:
                        self.y_s[k] = m / (1 - p) ** 2 * (self.chord - x) * (1 - 2 * p + x / self.chord)
        else:
            self.y_s = np.zeros(self.npt)

    def camber(self, x):
        """
        camberline equation
        """
        if self.airfoil[0] == 'parabole':
            return -self.eps * x / (self.chord / 4) ** 2 - np.deg2rad(self.alpha)
        elif self.airfoil[0] == "family_pq":
            [eps, p, q] = self.AirfoilParam
            eta = self.x / self.chord
            return eps * p * (eta * (1 - 2 * q) - 1 / 4 + 3 * eta ** 2) - np.deg2rad(self.alpha)
        elif self.airfoil[0] == "double-cambered":
            eta = self.x / self.chord
            return self.eps * (1 - 3 * eta ** 2 / (self.chord / 2) ** 2) - np.deg2rad(self.alpha)
        elif self.airfoil[0] == "naca4":
            delta = np.zeros(self.npt)
            [m, p, t] = self.AirfoilParam[1:]
            if p * m != 0:
                for k in range(self.npt):
                    x = self.x[k] + 0.5 * self.chord
                    if x / self.chord <= p:
                        delta[k] = 2 * m / p ** 2 * (p - x / self.chord) - np.deg2rad(self.alpha)
                    else:
                        delta[k] = 2 * m / (1 - p) ** 2 * (p - x / self.chord) - np.deg2rad(self.alpha)
            else:
                print("symmetrical airfoil")
                delta = -np.deg2rad(self.alpha) * np.ones(self.npt)
            return delta
        else:
            return np.zeros(len(x))

    def integrand(self, theta, n=0):
        """
        integrand
        """
        delta = self.camber(self.chord / 2 * np.cos(theta))
        return delta * np.cos(n * theta)

    def SerieFourier(self):
        """
        n first coefficient a_k = 4/pi int_0^pi delta(theta)cos(n theta) dtheta
        """
        self.a = np.zeros(self.nFourier)
        if self.order == 2:
            print('order 2')
            for i in range(self.nFourier):
                self.a[i] = 4 / np.pi * np.trapz(self.delta_s * np.cos(i * self.theta), self.theta)
        elif self.order == 4:
            print('order 4')
            for i in range(self.nFourier):
                self.a[i], res = 4 / np.pi * integrate.quad(self.integrand, 0, np.pi, args=(i,))
        else:
            raise ValueError('bad order for the SerieFourier')

    def Aerodynamics(self):
        """
        Calculus of the aerodynamic coefficients
        """
        self.CL = -np.pi / 2 * (self.a[0] + self.a[1])
        self.alpha_0 = (self.a[0] + self.a[1]) / 4
        self.CmF = -np.pi / 8 * (self.a[2] + self.a[1])
        self.alpha_ad = self.a[0] / 4
        self.Cm_LE = -np.pi / 8 * (self.a[0] + 2 * self.a[1] + self.a[2])

        if self.alpha != 0:
            print("BE CAREFUL :  \t the calculus of alpha_0 and alpha_d are wrong")
            print("\t \t because of Angle of Attack is not null")

        print("alpha                                     : %f °" % (self.alpha))
        print("CL                                        : %f" % (self.CL))
        if self.alpha == 0:
            print("CL = 2 pi (alpha - alpha_0)               : %f" % (
                        2 * np.pi * (np.deg2rad(self.alpha) - self.alpha_0)))
            print("alpha_0                                   : %f °" % (np.rad2deg(self.alpha_0)))
            print("alpha_adaptation                          : %2.15f °" % (np.rad2deg(self.alpha_ad)))
        print("Cm at the aerodynamic center              : %f" % (self.CmF))
        print("Cm at the leading edge                    : %f" % (self.Cm_LE))

    def Kp_portant(self):
        """
        Calculus of Kp for the lifting problem
        """
        self.Kp_s = np.zeros(self.npt)
        for i in np.arange(1, self.nFourier):
            self.Kp_s += self.a[i] * np.sin(i * self.theta)
            print("a %i %f" % (i, self.a[i]))

        if abs(self.a[0]) >= 1e-7:
            print("taking account of a_0                   : %e" % (self.a[0]))
            self.Kp_s += self.a[0] / 2 * np.tan(self.theta / 2)
            self.Kp_s[-1] == float("inf")

        # print("Kp_s = ",self.Kp_s)

    def graphiques(self):
        """
        various plots
        """

        if self.plot_squelette:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Airfoil used in the problem ' + self.probleme, fontsize=fz, fontweight='bold')
            plt.xlabel(self.xlabel, fontsize=fz)

            if self.probleme == "portant":
                print("len x, y_s :", len(self.x), len(self.y_s), self.Lref)
                plt.plot(self.x / self.Lref, self.y_s / self.Lref, color='red', label=r'$y_s$')
                plt.ylabel(r'$y_s$', fontsize=fz)
            plt.axis("equal")
            plt.grid()
            plt.legend()
            plt.xlim(-0.55 * self.chord / self.Lref, 0.55 * self.chord / self.Lref)

        if self.plot_Kp:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Pressure coefficient of the problem  ' + self.probleme, fontsize=fz, fontweight='bold')
            plt.xlabel(self.xlabel, fontsize=fz)

            if self.probleme == "portant":
                plt.plot(self.x / self.Lref, self.Kp_s, color='red', label=r'$Kp_s$')
                plt.ylabel(r'$Kp_s$', fontsize=fz)

            plt.grid()
            plt.legend()
            plt.ylim(self.KpLim[0], self.KpLim[1])
        plt.show()

    def airfoil_geometry(self):
        """
        Definition of the camberline from airfoil data
        Inputs:
        x+, y+, x-, y- (+ and - for upper and lower wall)
        """
        # The airfoil must be interpolated in a new mesh x
        ye = interp1d(self.xe, self.ye, kind='linear')
        yi = interp1d(self.xi, self.yi, kind='linear')

        print("upper wall: min and max of x (airfoil)    : %f et %f" % (np.min(self.xe), np.max(self.xe)))
        print("lower wall: min and max of x (airfoil)    : %f et %f" % (np.min(self.xi), np.max(self.xi)))
        print("min and max of x (new)                    : %f et %f" % (np.min(self.x), np.max(self.x)))
        # The geometrical incidence of the airfoil is verified
        incidence_profil = np.arctan(self.ye[0] - self.ye[-1])
        print("initial incidence                         : %f °, yTE = %f yLE = %f x_TE = %f x_LE = %f" % (
        np.rad2deg(incidence_profil), self.ye[0], self.ye[-1], self.xe[0], self.xe[-1]))
        print("chord                                     : %f" % (self.chord))
        # Airfoil with the new x-vector :
        y_plus, y_minus = ye(self.x), yi(self.x)
        self.y_s, self.y_e = (y_plus + y_minus) / 2, (y_plus - y_minus) / 2
        # small problems with non symmetrical airfoils, it is better to extrapolate y_s at the trailing edge
        # ys_tmp= InterpolatedUnivariateSpline(np.flipud(self.x[-4:-1]),np.flipud(self.y_s[-4:-1]),k=2)
        ys_tmp = InterpolatedUnivariateSpline(self.x[-4:-1], self.y_s[-4:-1], k=2)
        self.y_s[-1] = ys_tmp(self.x[-1])

        thickness, ind_e = np.amax(self.y_e), np.argmax(self.y_e)
        self.cambrure, ind_c = np.amax(self.y_s), np.argmax(self.y_s)
        print(" camberline / c                             : %5.4f à x(%4d)/c = %2.4f soit %f pourcents de chord" % (
        self.cambrure / self.Lref, ind_c, \
        self.x[ind_c] / self.Lref, 100 * (self.x[ind_c] / self.Lref + 0.5)))
        print(" thickness / c                            : %5.4f à x(%4d)/c = %2.4f soit %f  pourcents de chord" % (
        2 * thickness / self.Lref, \
        ind_e, self.x[ind_e] / self.Lref, 100 * (self.x[ind_e] / self.Lref + 0.5)))

        # slopes :
        dy = self.y_s[1:] - self.y_s[0:-1]
        dx = self.x[1:] - self.x[0:-1]
        print("len y = ", len(dy), "len ys = ", len(self.y_s), len(self.y_s[:-1]))
        print("dy[0] =", dy[0], self.y_s[1] - self.y_s[0])

        self.delta_s = np.zeros(self.npt)
        self.delta_s[0:-1] = dy / dx
        # self.delta_s[-2:]= self.delta_s[-3]
        # self.delta_s[-1]= self.delta_s[-2]
        #  small problems with non symmetrical airfoils, it is better to extrapolate y_s at the trailing edge
        # ys_tmp= InterpolatedUnivariateSpline(np.flipud(self.x[-4:-1]),np.flipud(self.y_s[-4:-1]),k=2)
        deltas_tmp = InterpolatedUnivariateSpline(self.x[-4:-1], self.delta_s[-4:-1], k=2)
        self.delta_s[-1] = deltas_tmp(self.x[-1])
        self.delta_s -= np.deg2rad(self.alpha)

        self.PlotDelta_s()

        if self.plot_airfoil:
            self.PlotAirfoil()


    def PlotAirfoil(self):
        """
        Airfoil plot
        """
        symbole = False

        plt.figure(figsize=(18, 4))
        plt.title("Profil")
        plt.xlabel("x", fontsize=fz)
        plt.ylabel("y", fontsize=fz)
        if symbole:
            plt.plot(self.x, self.y_s, 'bo', label=r'$y_s$')
            plt.plot(self.x, self.y_e, 'b--', label=r'$y_e$')
            plt.plot(self.xe, self.ye, 'ko', label=r'ext.')
            plt.plot(self.xi, self.yi, 'rs', label=r'int.')
        else:
            plt.plot(self.x, self.y_s, 'k-', label=r'$y_s$')
            plt.plot(self.x, self.y_e, 'r-', label=r'$y_e$')
            plt.plot(self.xe, self.ye, 'k--', label=r'ext.')
            plt.plot(self.xi, self.yi, 'r--', label=r'int.')
        plt.grid()
        plt.axis("equal")
        print("x = ", self.x[:2], self.x[-3:])
        print("ys = ", self.y_s[:2], self.y_s[-3:])
        print("deltas = ", self.delta_s[:2], self.delta_s[-3:])
        plt.legend()
        plt.show()

    def PlotDelta_s(self):
        """
        camberline plot (delta_s)
        """
        plt.figure(figsize=(14, 6))
        plt.title("slope of the camberline")
        plt.xlabel("x", fontsize=fz)
        plt.ylabel(r"$\delta_s$", fontsize=fz)
        plt.plot(self.x, np.rad2deg(self.delta_s), 'bo')
        # plt.plot(self.x,self.y_s,'rs')
        plt.grid()
        plt.xlim(-0.55, 0.55)
        plt.ylim(-20, 20)
        plt.show()
