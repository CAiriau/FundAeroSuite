#!/bin/python
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Conformal mapping package

for airfoil, Joukowski transform, mappings file, with the main class
"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline

fz = 16  # font size for legend


# *******************************************
def set_parameters():
    # *******************************************
    """
    Set default parameters
    """
    prm = {}
    # physical parameters
    prm["Uinf"] = 1.        # upstream velocity m/s
    prm["rho"] = 1.3        # gas density kg/m^3
    prm["alpha"] = 0.0      # angle of attack in degrees
    prm["Umax"] = 10        # velocity filter for singular points
    prm["Kpmin"] = -5       # y-axis range for Kp plot  

    # mapping parameter
    prm["airfoil"] = "Jouko1"   # testcase name, airfoil name
    prm["map"] = "Joukowski"    # transformation name
    prm["Zc"] = 0.0 + 1j * 0.0  # circle center position 
    prm["R"] = 1.               # circle radius
    prm["beta"] = 0.0           # beta in degrees, angle in the circulation formula
    prm["Xgrille"] = np.arange(-3, 3, 0.05)  # define the X-grid in the Z complex plane  
    prm["Ygrille"] = np.arange(-3, 3, 0.05)  # define the Y-grid in the Z complex plane
    prm["lamb"] = 1.2           # for a = lamb x R
    prm["k"] = 1.9              # Karman-Trefftz transformation exponent
    prm["n_circle"] = 361       # In the conformal plane, number of points on the circle
    prm["Zs"] = 1 + 1j * 0      # In the conformal plane, trailing edge singular point affix
    prm["zs"] = 1 + 1j * 0      # In the physical plane, trailing edge singular point affix
    prm["eps"] = 0.05           # small parameter in the transformation when it is used

    # for iso psi :
    prm["n_iso"] = 51               # number of level curve targeted 
    prm["levels_iso"] = [-2.4, 3.4, 0.2]  # levels [psi_min,psi_max,delta_psi]
    prm["opt_levels"] = "manuel"    # option to calculate the levels :  "manuel" or "nombre"
    prm["eps_levels"] = 0.1         # truncated value to draw iso psi [psi_min+eps, psi_max-eps]
    prm["Psi_limit"] = [-1, 1]      # Psi_min,Psi_max
    prm["Psi_method"] = "Psi"       # "Psi" or "U", method to plot streamlines

    # plot parameters
    prm["plot_Kp"] = True           # to plot pressure coefficient on the body
    prm["calcul_Kp"] = True         # to calculate pressure coefficient on the body
    prm["plot_airfoil"] = True      # to plot airfoil
    prm["plot_velocity"] = False    # to plot  U/Uinf on the body
    prm["plot_circle"] = True       # to plot the circle in the conformal plane Z  
    prm["plot_psi"] = True          # to plot streamlines in the physical plane
    prm["adimChord"] = True         # true : non dimensional airfoil
    prm["camberline"] = True        # to calculate the camber line and thickness laws

    prm["xsize"] = 10               # figure width 
    prm["ysize"] = 8.5              # figure height  

    return prm


# ----------|----------------------------------|-----------------|---------------------| 
# Plane I   | plane of a centered circle       | conformal plane | Z'= R exp(i theta)  |
# Plane II  | plane with a non centered circle | conformal plane | Z = Zc + Z'         |
# Plane III | airfoil plane                    | plan physique   | z = H(Z') = H(Z-Zc) |
# ----------|----------------------------------|-----------------|---------------------| 

class AirfoilConfMap(object):
    """
    airfoil designed by a conformal mapping
    """

    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#', 60 * '*')
        print('# %s' % "Conformal mappings")
        print('#', 60 * '*', '\n')
        zeroc = 0 + 1j * 0
        # physical parameters

        self.Uinf = prm["Uinf"]     # upstream velocity m/s
        self.rho = prm["rho"]       # gas density kg/m^3
        self.alpha = prm["alpha"]   # angle of attack in degrees
        self.Gamma = 0.             # circulation to calculate 
        self.w = 0.                 # complex velocity on the airfoil
        self.Kp = 0.                # pressure coefficient
        self.Psi = 0.               # iso Psi, streamlines
        self.chord = 0.             # airfoil chord
        self.camber = 0.            # relative camber
        self.thickness = 0.         # relative thickness
        self.A0 = 0.                # first coefficient of the Laurent's expansion
        self.A1 = 0.                # second coefficient of the Laurent's expansion
        self.alpha0 = 0.            # zero-lift angle of attack
        self.CL = 0.                # lift coefficient
        self.zF = zeroc             # aerodynamic center affix /  chord
        self.CmF = 0.               # pitching coefficient at the aerodynamic center

        self.Umax = prm["Umax"]     # velocity filter for singular points
        self.Kpmin = prm["Kpmin"]   # y-axis range for Kp plot  
        self.camberline = prm["camberline"]  # to calculate the camber line and thickness laws
        # mapping parameter

        self.airfoil = prm["airfoil"]   # testcase name, airfoil name
        self.map = prm["map"]           # transformation name
        self.Zc = prm["Zc"]             # circle center position 
        self.R = prm["R"]               # circle radius
        self.beta = prm["beta"]         # beta in degrees, angle in the circulation formula
        self.eps = prm["eps"]           # small parameter in the transformation when it is used
        self.Xgrille = prm["Xgrille"]   # define the X-grid in the Z complex plane
        self.Ygrille = prm["Ygrille"]   # define the Y-grid in the Z complex plane
        self.lamb = prm["lamb"]         # for a = lamb x R 
        self.k = prm["k"]               # Karman-Trefftz transformation exponent
        self.n_circle = prm["n_circle"]  # In the conformal plane, number of points on the circle
        self.Zs = prm["Zs"]             # In the conformal plane, trailing edge singular point affix
        self.z_airfoil = zeroc          # airfoil affix in the physical plane  
        self.Z_cercle = zeroc           # non centered circle in the conformal plane  (plane II)
        self.Z_le = zeroc               # leading edge airfoil affix on the circle (conformal map)
        self.z_le = zeroc               # leading edge airfoil affix  (physical map)
        self.Z_stagnation = zeroc       # affix of the stagnation point on the circle
        self.z_stagnation = zeroc       # affix of the stagnation point on the airfoil
        self.mu = [0, 0]                # von Mises's transformation parameter
        self.b = zeroc                  # parameter for the transformation with sharp L.E and T.E airfoil
        self.a = 0                      # Joukowski's transformation parameter
        self.xe = 0
        self.xi = 0
        self.ye = 0
        self.yi = 0
        self.Kpe = 0
        self.Kpi = 0

        # for iso psi :

        self.n_iso = prm["n_iso"]  # number of level curve targeted
        self.levels_iso = prm["levels_iso"]  # levels [psi_min,psi_max,delta_psi]
        self.opt_levels = prm["opt_levels"]  # option to calculate the levels :  "manuel" or "number"
        self.eps_levels = prm["eps_levels"]  # truncated value to draw iso psi [psi_min+eps, psi_max-eps]
        self.Psi_limit = prm["Psi_limit"]    # Psi_min,Psi_max
        self.Psi_method = prm["Psi_method"]  # "Psi" or "U", method to plot streamlines
        # plot parameters
        self.calcul_Kp = prm["calcul_Kp"]   # to calculate pressure coefficient on the body
        self.plot_Kp = prm["plot_Kp"]       # to plot pressure coefficient on the body
        self.plot_airfoil = prm["plot_airfoil"]     # to plot airfoil
        self.plot_velocity = prm["plot_velocity"]   # to plot  U/Uinf on the body
        self.plot_circle = prm["plot_circle"]       # to plot the circle in the conformal plane Z
        self.xsize = prm["xsize"]            # figure width
        self.ysize = prm["ysize"]            # figure height
        self.plot_psi = prm["plot_psi"]      # to plot streamlines in the physical plane
        self.adimChord = prm["adimChord"]    # true : non dimensional airfoil
        print("initialization done")

        if self.lamb <= 1:
            raise ValueError('R/a must be greater than 1')

    def aero_center(self, A1, chord):
        """
        Aerodynamic center position given from A1 and alpha_0,
        here  beta = -alpha_0
        """
        z_F = -A1 / self.R * np.exp(1j * self.beta)
        print("z_F/chord                                 : %f + i %f" % (z_F.real / chord, z_F.imag / chord))
        print("z_F                                       : %f + i %f" % (z_F.real, z_F.imag))

    def run_airfoil(self):
        """
        main function to run the code
        """

        print("be careful: if     a=0  a is calculated from R (R=1 by default)")
        print("            else        R is calculated from en a ")

        self.set_title("Type of transformation :" + self.map)
        self.beta = np.deg2rad(self.beta)
        self.alpha = np.deg2rad(self.alpha)

        if self.map == "Joukowski":
            self.a = self.R / self.lamb
            self.Zs = self.a * np.exp(-1j * self.alpha)
            self.Zc = np.exp(-1j * self.alpha) * (self.a - self.R * np.exp(-1j * self.beta))
            self.Z_le = self.Zc + self.R * np.exp(1j * (np.pi + self.beta - self.alpha))
            self.Z_stagnation = self.Zc + self.R * np.exp(1j * (self.beta + np.pi + self.alpha))



        # elif (self.map=="Joukowski1") or (self.map=="Karman-Trefftz")  \
        #      or (self.map=="von Mises") or  self.map=="van de Vooren" or self.map=="double sharp":

        elif self.map in ["Joukowski1", "Karman-Trefftz", "von Mises", "van de Vooren", "double sharp"]:

            if self.a == 0:
                print(" a is determined")
                self.a = np.sqrt(self.R ** 2 - self.Zc.imag ** 2) + self.Zc.real
            else:
                self.R = np.abs(self.a - self.Zc)
                print("==> parameter 'a' is given")

            self.Zs = self.a
            # self.beta=np.arccos((self.a-self.Zc.real)/self.R)
            self.beta = np.arcsin((self.Zc.imag) / self.R)
            self.Z_le = 2 * self.Zc.real - self.a  #   approximation, wrong for large camberline 
            self.Z_stagnation = self.Zc + self.R * np.exp(1j * (self.beta + np.pi + 2 * self.alpha))

        if self.map == "Joukowski1":
            if self.airfoil == "Ellipse":
                print("ellipse case")
                # for an ellipse  a is only given
                if self.a >= self.R:
                    self.a = 0.9 * self.R
                    print(' a is decreased ', self.a)
                self.Zs = self.R
                self.beta = 0
                self.Z_le = -self.R
                self.beta = 0
                self.Z_stagnation = self.Zc + self.R * np.exp(1j * (self.beta + np.pi + 2 * self.alpha))
                print("Chord (analytical)                         : %f" % (2 * (self.R + self.a ** 2 / self.R)))
            chord = 3 * self.a - 2 * self.Zc.real + self.a ** 2 / (self.a - 2 * self.Zc.real)
            print("Analytical chord (Joukowski)              : %f" % (chord))

        elif self.map == "von Mises":
            self.mu = [-self.a + self.eps * (1 + 1j), -self.eps * (1 + 1j)]
            if self.k != 2:
                print(' k value is changed to 2')
                self.k = 2
            Xc = self.Zc.real
            Chord = (4 * (self.a ** 2 + (-2 * Xc + 1 / 2 * self.mu[0]) * self.a + 1 / 2 * self.mu[0] ** 2)) * (
                        self.a - Xc) ** 2 / (self.a * (-2 * Xc + self.a) ** 2)
            chord = Chord.real
            print("Analytical chord                          : %f et | %f |" % (chord, abs(Chord)))
            A1 = self.a ** 2 + self.a * self.mu[0] + self.mu[0] ** 2
            print("A1                                        : %f + i %f" % (A1.real, A1.imag))
            self.aero_center(A1, chord)

        elif self.map == "van de Vooren":
            self.L = self.a * 2 * pow(2 / (1 + self.eps), self.k - 1)
            print("chord L                                   : %f" % (self.L))

        elif self.map == "double sharp":
            # not clear to change the leading edge position w.r.t. non zero alpha, 
            # otherwise there is a problem with the geometrical accuracy 
            if self.alpha != 0:
                self.Z_le = self.Zc + self.R * np.exp(1j * (2 * self.beta + np.pi - self.alpha))
            print(self.Z_le)

        elif self.map == "Karman-Trefftz":
            if abs(self.Zc.real) <= 1e-10:
                raise ValueError('X_c must be different from O in Karman-Trefftz')
            A = 1 - self.a / self.Zc.real
            print("A = ", A)
            chord = 2 * self.k * self.a * pow(A, self.k) / (pow(A, self.k) - 1)
            A1 = (self.k ** 2 - 1) * self.a ** 2 / 3
            print("analytical chord                          : %f" % (chord))
            print("A1                                        : %f + i %f" % (A1.real, A1.imag))
            self.aero_center(A1, chord)

        self.z_stagnation = self.H(self.Z_stagnation)
        self.z_le = self.H(self.Z_le)
        self.zs = self.H(self.Zs)
        self.Gamma = -4 * np.pi * self.Uinf * self.R * np.sin(self.alpha + self.beta)  # airfoil circulation

        print("a                                         : %f" % (self.a))
        print("R                                         : %f" % (self.R))
        print("beta                                      : %3.10f °" % (np.rad2deg(self.beta)))
        print("Circulation                               : %f" % (self.Gamma))
        print("Circle plane                              :")
        print("Circle center Zc                          : %2.10f+ i %2.10f" % (self.Zc.real, self.Zc.imag))
        print("Singular point Zs                         : %f+ i %f" % (self.Zs.real, self.Zs.imag))
        print("Leading edge Z_le                         : %f+ i %f" % (self.Z_le.real, self.Z_le.imag))
        print("Stagnation point close to l.e. Z_stag.    : %f+ i %f" % (self.Z_stagnation.real, self.Z_stagnation.imag))
        print("Airfoil plane                             :")
        print("Trailing edge  z_s                        : %f+ i %f" % (self.zs.real, self.zs.imag))
        print("leading edge z_le                         : %f+ i %f" % (self.z_le.real, self.z_le.imag))
        print("stagnation point  z_stagnation            : %f+ i %f" % (self.z_stagnation.real, self.z_stagnation.imag))

        self.Z_cercle = self.cercle()                   # circle in  plane II
        self.z_airfoil = self.H(self.Z_cercle)          # circle in plane III
        self.airfoil_geometry()

        CL = 8 * np.pi / self.chord * np.sin(self.alpha + self.beta)
        print("CL                                        : %f " % (CL))

        self.SerieLaurent()
        #  Kp ou U  
        if self.calcul_Kp:
            self.display_Kp(self)

        if self.Psi_method == "Psi":
            #  streamlines :
            self.Streamlines()
            if self.plot_psi: self.streamlines_plot()
        elif self.Psi_method == "U":
            self.iso_velocity()
        else:
            raise ValueError('bad option to plot streamlines')

    def iso_velocity(self):
        """
        complex velocity outside airfoil
        """
        # Grid definition in the circle plane
        X, Y = np.meshgrid(self.Xgrille, self.Ygrille)
        Zgrid = X + 1j * Y
        Zgrid = ma.masked_where(np.absolute(Zgrid - self.Zc) <= self.R, Zgrid)
        zgrid = self.H(Zgrid)

        W = self.velocity_W(Zgrid - self.Zc)
        # w_g  = W /self.H(Zgrid,opt=1)
        # u_g,v_g=w_g.real,-w_g.imag   # unused

        if self.plot_psi:
            planZ = False
            print('streamlines plot')
            # Nx,Ny=len(self.Xgrille),len(self.Ygrille)
            # y_start,y_end=self.zgrid[0,0].imag,self.zgrid[0,-1].imag
            # x_start,x_end=self.zgrid[0,0].real,self.zgrid[-1,0].real
            size = 15
            # plt.figure(figsize=(size, (y_end-y_start)/(x_end-x_start)*size))
            plt.figure(figsize=(size, size * 0.8))
            plt.xlabel('x', fontsize=16)
            plt.ylabel('y', fontsize=16)
            # plt.xlim(x_start,x_end);plt.ylim(y_start,y_end) print(len(self.zgrid.real),len(u_g)) print(np.shape(
            # self.zgrid.real)) print(np.shape(u_g)) plt.streamplot(self.zgrid.real, self.zgrid.imag, self.w.real,
            # -self.w.imag, density=1.0, linewidth=1, arrowsize=1, arrowstyle='->')
            if planZ:
                plt.streamplot(Zgrid.real, Zgrid.imag, W.real, -W.imag, density=1.0, linewidth=1, arrowsize=1,
                               arrowstyle='->')
                plt.scatter(Zgrid.real, Zgrid.imag, color='green', s=4, marker='o', linewidth=0)
                plt.plot(self.Z_cercle.real, self.Z_cercle.imag, label='airfoil')
            else:
                # for k in range(len(zgrid.real)): print(len(zgrid.imag[k,:]),len(w_g.real[k,:]),max(zgrid.imag[k,
                # :]),min(zgrid.imag[k,:])) plt.streamplot(zgrid.real, zgrid.imag, w_g.real, -w_g.imag, density=1.0,
                # linewidth=1, arrowsize=1, arrowstyle='->')
                plt.scatter(zgrid.real, zgrid.imag, color='green', s=4, marker='o', linewidth=0)
                plt.plot(self.z_airfoil.real, self.z_airfoil.imag, label='airfoil')
                plt.plot(self.zs.real, self.zs.imag, 'ks', label='T.E.')
                plt.plot(self.z_le.real, self.z_le.imag, 'ko', label='L.E.')
                plt.plot(self.z_stagnation.real, self.z_stagnation.imag, 'bs', label='Stagn.')
            # plt.scatter(xs, np.zeros((Ns), dtype=float),color='red', s=30, marker='o', linewidth=0); 
            # to plot points used for streamlining
            # plt.legend()
            #plt.show()

    def SerieLaurent(self):
        """
        transformation written as Laurent's series
        """
        serie = True
        if self.map == "Joukowski1":
            self.A0, self.A1 = -self.Zc, self.a ** 2
        elif self.map == "Karman-Trefftz":
            self.A0, self.A1 = -self.Zc, self.a ** 2 * (self.k ** 2 - 1) / 3
            print("trailing edge wedge angle                : %f °" % ((2 - self.k) * 180))
        elif self.map == "von Mises":
            self.A0, self.A1 = -self.Zc, self.a ** 2 - self.mu[0] * self.mu[1]
        elif self.map == "van de Vooren":
            self.A0, self.A1 = self.L - self.k * self.a, self.k * (self.k - 1) * (1 - 2 * self.eps) * self.a ** 2 / 2
            print("trailing edge wedge angle                : %f °" % ((2 - self.k) * 180))
        elif self.map == "double sharp":
            self.A0, self.A1 = self.Zc, self.a ** 2 + self.b ** 2
        else:
            serie = False
            print('No Laurent series implemented for this transformation')

        if serie:
            self.alpha0 = -self.beta  # radian
            d2 = abs(self.A1)
            gamma = np.angle(self.A1, deg=False) / 2
            self.CmF = 4 * np.pi * d2 / self.chord ** 2 * np.sin(2 * (gamma - self.alpha0))
            self.zF = -self.A1 / self.a * np.exp(-1j * self.alpha0) / self.chord

            print("First axis slope Delta_1                  : %f °" % (np.rad2deg(self.alpha0)))
            print("d^2                                       : %f" % (d2))
            print("gamma for the second axis Delta_2         : %f °" % (np.rad2deg(gamma)))
            print("CmF                                       : %f " % self.CmF)
            print("aero_center zF /chord affix              : %f + i %f " % (self.zF.real, self.zF.imag))

        return serie

    def H(self, Z, opt=0):
        """
        Joukowski conformal transformation
        
        Params:
            * a (complex) : transformation parameter
            * Z (complex) : affix in the conformal plane (of the circle)
            * opt (integer) ! 0 : z=H(Z),  1 : H'(Z)=dz/dZ
            
        Returns:
            complex : z = H(Z)
        """
        if self.map == "Joukowski":
            if opt == 0:
                return Z + (np.exp(-1j * 2 * self.alpha) * self.a ** 2) / Z
            else:
                return 1 - (np.exp(-1j * 2 * self.alpha) * self.a ** 2) / Z ** 2

        elif self.map == "Joukowski1":
            if opt == 0:
                return Z + self.a ** 2 / Z
            else:
                return 1 - self.a ** 2 / Z ** 2

        elif self.map == "Karman-Trefftz":
            G = (Z - self.a) / (Z + self.a)
            if opt == 0:
                return self.k * self.a * (1 + pow(G, self.k)) / (1 - pow(G, self.k))
            else:
                Gp = 2 * self.a / (Z + self.a) ** 2
                return 2 * self.k ** 2 * self.a * Gp * pow(G, self.k - 1) / (1 - pow(G, self.k)) ** 2

        elif self.map == "von Mises":
            if opt == 0:
                return Z + (self.a ** 2 - self.mu[0] * self.mu[1]) / Z + 0.5 * self.a * self.mu[0] * self.mu[1] / Z ** 2
            else:
                return (1 - self.a / Z) * (1 - self.mu[0] / Z) * (1 - self.mu[1] / Z)

        elif self.map == "van de Vooren":
            Z = Z + 1j * 0
            if opt == 0:
                return self.L + pow(Z - self.a, self.k) / pow(Z - self.eps * self.a, self.k - 1)
            else:
                return (Z + (self.a * (self.k - 1 - self.eps * self.k))) * pow(Z - self.a, self.k - 1) / pow(
                    Z - self.eps * self.a, self.k)

        elif self.map == "double sharp":
            Z = Z + 1j * 0
            if opt == 0:
                return Z + (self.a ** 2 + self.b ** 2) / Z - (self.a * self.b) ** 2 / 3 / Z ** 3
            else:
                return 1 - (self.a ** 2 + self.b ** 2) / Z ** 2 + (self.a * self.b) ** 2 / Z ** 4

        else:
            raise ValueError('bad conformal mapping option (self.map)')

    def cercle(self):
        """
        Plot of the circle of center Z_c, of radius R in the complex plane Z=(X,Y)
        """
        theta = np.linspace(0, 2 * np.pi, self.n_circle)
        Z = self.Zc + self.R * np.exp(1j * theta)
        alpha_Zs = np.rad2deg(np.arcsin((self.Zc - self.Zs).imag / self.R))
        print("Angle of singular point on the circle     : %f ° " % alpha_Zs)

        if self.plot_circle:
            k = 1.2
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Circle in the complex plane Z', fontsize=fz, fontweight='bold')
            plt.plot(Z.real, Z.imag, 'k-')
            plt.plot([-k * self.R, k * self.R], [0, 0], color='blue')
            plt.plot([0, 0], [-k * self.R, k * self.R], color='blue')
            plt.plot([-k * self.R, k * self.R], [self.Zc.imag, self.Zc.imag], 'r--')
            plt.plot([self.Zc.real, self.Zc.real], [-k * self.R, k * self.R], 'r--')
            plt.plot(self.Zc.real, self.Zc.imag, 'ro', label='center')
            plt.plot(self.Zs.real, self.Zs.imag, 'ks', label='T.E.')
            plt.plot(self.Z_le.real, self.Z_le.imag, 'ko', label='L.E.')
            plt.plot(self.Z_stagnation.real, self.Z_stagnation.imag, 'bs', label='Stagn.')
            plt.axis("equal")
            plt.legend()
            # plt.show()
        return Z

    def airfoil_geometry(self):
        """
        geometry characteristics of the airfoil
        """
        print("geometry characteristics of the airfoil : ")
        npt = self.n_circle
        Z_tmp = self.zs - self.z_le
        incidence = -np.rad2deg(np.arctan(Z_tmp.imag / Z_tmp.real))
        print("airfoil tilt in ° (inclination)           : %f  %f " % (incidence, -np.angle(Z_tmp, deg=True)))
        self.chord = (Z_tmp * np.exp(1j * np.deg2rad(incidence))).real
        print("airfoil chord                             : %f  %f " % (
        self.chord, Z_tmp.real / np.cos(np.deg2rad(incidence))))
        if self.adimChord:
            self.Lref = self.chord
        else:
            self.Lref = 1.0

        if self.map == "Joukowski":
            # upper wall
            theta_e = np.linspace(-self.alpha - self.beta + 0.01, np.pi + self.beta - self.alpha, npt)
            self.Ze = np.flipud(self.Zc + self.R * np.exp(1j * theta_e))
            ze = (self.H(self.Ze) - self.z_le) * np.exp(1j * np.deg2rad(incidence))
            # print(np.rad2deg(theta_e))
            # lower wall
            theta_i = np.linspace(np.pi + self.beta - self.alpha, 2 * np.pi - self.alpha - self.beta - 0.01, npt)
            self.Zi = self.Zc + self.R * np.exp(1j * theta_i)
            zi = (self.H(self.Zi) - self.z_le) * np.exp(1j * np.deg2rad(incidence))
            # print(np.rad2deg(theta_i))

        else:
            # change of coordinate : origin at the leading edge : translation of -z_le
            theta_e = np.linspace(-self.beta, np.pi + self.beta, npt)
            theta_i = np.linspace(np.pi + self.beta, 2 * np.pi - self.beta, npt)
            self.Ze = self.Zc + self.R * np.exp(1j * theta_e)
            # to get the complete airfoil :
            ze = np.concatenate(
                ([0 + 1j * 0], np.flipud((self.H(self.Ze, opt=0) - self.z_le) * np.exp(1j * np.deg2rad(incidence)))))
            ze[-1] = self.chord
            self.Zi = self.Zc + self.R * np.exp(1j * theta_i)
            zi = np.concatenate(
                ([0 + 1j * 0], (self.H(self.Zi, opt=0) - self.z_le) * np.exp(1j * np.deg2rad(incidence))))
            zi[-1] = self.chord
            # print("ze",ze[0:5],ze[-5:])

        # Look for the real trailing edge of the airfoil :
        # and we can have some small geometrical problem in the leading edge
        Min = max([np.min(ze.real), np.min(zi.real)])
        Max = min([np.max(ze.real), np.max(zi.real)])
        # x=(1-np.cos(np.linspace(0,np.pi-0.1,npt)))/2    
        x = np.linspace(Min, Max, npt)

        # Look for the real leading edge of the airfoil :
        Imin = [np.argmin(ze.real), np.argmin(zi.real)]
        print('Imin                                      : ', Imin)
        # print("le: Ze,Zi                                 : %3.4e + i %3.4e, %3.4e + i %3.4e"%(self.Ze[Imin[
        # 0]].real,self.Ze[Imin[0]].imag ,self.Zi[Imin[1]].real,self.Zi[Imin[1]].imag)) 
        print("le: ze,zi                                 : %3.4e + i %3.4e, %3.4e + i %3.4e" % (
        ze[Imin[0]].real, ze[Imin[0]].imag, zi[Imin[1]].real, zi[Imin[1]].imag))
        print("le : thetae, thetai                       : %3.2f, %3.2f" % (
        np.rad2deg(theta_e[npt - Imin[0] - 1]), np.rad2deg(theta_i[Imin[1]])))
        c = abs(ze[Imin[0]] - ze[-1])
        print("corrected chord + relative error          : %f %f " % (c, abs(c / self.chord - 1)))
        print("corrected angle of attack                 : %f °" % (np.angle(-ze[Imin[0]] + ze[-1], deg=True)))

        if self.adimChord:
            ze, zi = ze / self.chord, zi / self.chord
            x_lim = 1.0
            x = x / self.chord
            xlabel, ylabel = r'$x/c$', r'$y/c$'
        else:
            x_lim = self.chord
            xlabel, ylabel = r'$x$', r'$y$'

        print("ze: min,max                               : %5.15e  %5.15e " % (np.min(ze.real), np.max(ze.real)))
        print("zi: min,max                               : %5.15e  %5.15e " % (np.min(zi.real), np.max(zi.real)))
        print("x : min,max                               : %5.15e  %5.15e " % (np.min(x), np.max(x)))

        if self.camberline:
            interpol = True
            # be careful, sometimes the second method does not work
            # let interpol=True except in case of problem on the vector boundary (LE, TE).
            if interpol:
                fe = interp1d(ze.real, ze.imag, kind='linear')
                fi = interp1d(zi.real, zi.imag, kind='linear')
            else:
                fe = InterpolatedUnivariateSpline(ze.real, ze.imag, k=2)
                fi = InterpolatedUnivariateSpline(zi.real, zi.imag, k=2)

            # To avoid interpolation problem with round-off error :   
            ValMax = min([np.max(ze.real), np.max(zi.real)])
            ValMin = max([np.min(ze.real), np.min(zi.real)])
            x[0] = ValMin
            x[-1] = ValMax
            ye, yi = fe(x), fi(x)
            ys, yt = (ye + yi) / 2, (ye - yi) / 2
            # print("yt",yt[0:5],yt[-5:])

            self.thickness, ind_t = np.amax(yt), np.argmax(yt)
            self.camber, ind_c = np.amax(ys), np.argmax(ys)
            print(" cambrer / c                              : %5.4f à x(%4d)/c = %2.4f " % (
            self.camber / self.Lref, ind_c, x[ind_c] / self.Lref))
            print(" thickness / c                            : %5.4f à x(%4d)/c = %2.4f " % (
            2 * self.thickness / self.Lref, ind_t, x[ind_t] / self.Lref))

        if self.plot_airfoil:
            plt.figure(figsize=(16, 4))
            plt.title(r'airfoil at 0 AoA', fontsize=fz, fontweight='bold')
            plt.plot(ze.real, ze.imag, 'k-', linewidth=1, label=r'$y_e$ (upper wall)')
            plt.plot(zi.real, zi.imag, 'b-', linewidth=1, label=r'$y_i$ (lower wall)')
            # plt.plot((self.z_airfoil-self.z_le).real/self.Lref,(self.z_airfoil-self.z_le).imag/self.Lref,
            # 'r',linewidth=1,label='profil')
            # plt.plot(x,ye,'rs')
            # plt.plot(x,yi,'rs')
            if self.camberline:
                plt.plot(x, ys, 'b--', linewidth=1, label=r'$y_s$ (camber)')
                plt.plot(x, yt, 'r--', linewidth=1, label=r'$y_t$ (thickness)')
            plt.xlabel(xlabel, fontsize=fz)
            plt.ylabel(ylabel, fontsize=fz)
            plt.axis("equal")
            plt.legend()
            plt.grid()
            plt.xlim(-0.05, x_lim + 0.05)
            # plt.show()

    def w_Kp(self, Z, display=False):
        """
        return velocity and Kp on the airfoil walls
        """
        W = self.velocity_W(Z - self.Zc)
        if display:
            print("Kp \t \t W \t \t z \t \t \t \t theta \t \t H'  : ")
            s = "{%3.4e} \t {%3.4e} \t {%3.4e} \t {%3.4e} \t {%4.2f} \t {%3.4e}"
            for Wtmp, Ztmp in zip(W, Z):
                ztmp = (self.H(Ztmp) - self.z_le) / self.Lref
                Kp = 1 - abs(Wtmp / self.H(Ztmp, opt=1)) ** 2
                print(s % (
                Kp, abs(Wtmp), ztmp.real, ztmp.imag, np.angle(Ztmp - self.Zc, deg=True), abs(self.H(Ztmp, opt=1))))

        Hp = self.H(Z, opt=1)
        w = W * 1e-10  # for initialization
        for i in range(len(W)):
            hp = Hp[i]
            if abs(hp) >= 1e-10:
                w[i] = W[i] / hp
            else:
                print("exception :  Z, H(Z)' = 0                 : ", Z[i], hp)
                w[i] = 0

        for i in range(len(w)):
            if abs(w[i]) > self.Umax or (np.isnan(w[i])):
                ztmp = self.H(Z[i])
                print("Filtered complex velocity in z = %2.5f + i %2.5f for w = %3.3e + i %3.3e" % (ztmp.real,
                                                                                                     ztmp.imag,
                                                                                                     w[i].real,
                                                                                                     w[i].imag))
                w[i] = float("-inf")
        Kp = 1 - np.abs(w / self.Uinf) ** 2
        print('Kp max                                    : %f ' % abs(np.max(Kp)))
        print('Kp min                                    : %f ' % abs(np.min(Kp)))

        return w, Kp

    def display_Kp(self, x):
        """
        Kp(x) plot
        """
        option_Kp = False
        if option_Kp:
            self.w, self.Kp = self.w_Kp(self.Z_cercle)
        else:
            self.w, self.Kp = self.w_Kp(self.Z_cercle)
            print("upper wall :")
            self.we, self.Kpe = self.w_Kp(self.Ze, display=False)
            print("lower wall :")
            self.wi, self.Kpi = self.w_Kp(self.Zi, display=False)

        if self.plot_Kp:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Pressure coefficient', fontsize=fz, fontweight='bold')
            if option_Kp:
                plt.plot(self.z_airfoil.real / self.Lref, self.Kp, color='black')
            else:
                plt.plot((self.z_airfoil - self.z_le).real / self.Lref, self.Kp, color='black', label='Kp')
                ze_tmp = self.H(self.Ze) - self.z_le
                zi_tmp = self.H(self.Zi) - self.z_le
                self.xe = ze_tmp.real / self.Lref
                self.xi = zi_tmp.real / self.Lref
                self.ye = ze_tmp.imag / self.Lref
                self.yi = zi_tmp.imag / self.Lref

                plt.plot(self.xe, self.Kpe, color='red', label=r'$Kp^+$')
                plt.plot(self.xi, self.Kpi, color='blue', label=r'$Kp^-$')
                plt.plot(ze_tmp.real / self.Lref, self.H(self.Ze).imag / self.Lref, color='red', label=r'$y^+$')
                plt.plot(zi_tmp.real / self.Lref, self.H(self.Zi).imag / self.Lref, color='blue', label=r'$y^-$')

            # plt.axis("equal")
            plt.grid()
            plt.legend()
            plt.ylim(self.Kpmin, -self.Kpmin)
            # plt.show()

        if self.plot_velocity:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'velocity', fontsize=fz, fontweight='bold')
            if option_Kp:
                plt.plot(self.z_airfoil.real / self.Lref, self.w, color='black')
            else:
                ze_tmp = self.H(self.Ze) - self.z_le
                zi_tmp = self.H(self.Zi) - self.z_le
                plt.plot(ze_tmp.real / self.Lref, abs(self.we), color='red', label=r'$U^+$')
                plt.plot(zi_tmp.real / self.Lref, abs(self.wi), color='blue', label=r'$U^-$')
                plt.plot(ze_tmp.real / self.Lref, self.H(self.Ze).imag / self.Lref, color='red', label=r'$y^+$')
                plt.plot(zi_tmp.real / self.Lref, self.H(self.Zi).imag / self.Lref, color='blue', label=r'$y^-$')
            # plt.axis("equal")
            plt.grid()
            plt.legend()
            # plt.show()

            print('w max                                     : %f ' % np.max(np.abs(self.w)))
            print('w min                                     : %f ' % np.min(np.abs(self.w)))
            print(self.w)

    def F(self, Z):
        """
        complex velocity potential in the plane I, f(z)=F(Z) 
        """
        U = np.zeros(Z.shape, dtype=np.complex)

        if self.map == "Joukowski":
            with np.errstate(divide='ignore'):
                for m in range(Z.shape[0]):
                    for n in range(Z.shape[1]):  # be careful numpy bug https://github.com/numpy/numpy/issues/8516
                        # element by element
                        U[m, n] = self.Gamma * np.log((Z[m, n]) / self.R) / (2 * np.pi)
            return self.Uinf * (Z + self.R ** 2 / Z) - 1j * U
        else:
            with np.errstate(divide='ignore'):
                for m in range(Z.shape[0]):
                    for n in range(Z.shape[1]):  # be careful numpy bug https://github.com/numpy/numpy/issues/8516
                        # element by element
                        U[m, n] = self.Gamma * np.log((Z[m, n]) / self.R) / (2 * np.pi)
            return self.Uinf * (Z * np.exp(-1j * self.alpha) + self.R ** 2 * np.exp(1j * self.alpha) / Z) - 1j * U

    def velocity_W(self, Z):
        """
        complex velocity W in the complex plane I (circle)
        """
        if self.map == "Joukowski":
            return self.Uinf * (1 - self.R ** 2 / Z ** 2) - 1j * self.Gamma / (2 * np.pi * Z)
        else:
            return self.Uinf * (
                        np.exp(-1j * self.alpha) - self.R ** 2 * np.exp(1j * self.alpha) / Z ** 2) - 1j * self.Gamma / (
                               2 * np.pi * Z)

    def Streamlines(self):
        """
        calculus of streamlines, psi(x,y) in the conformal and physical planes
        """
        # alpha, beta are in degrees

        # Grid in the circle plane
        X, Y = np.meshgrid(self.Xgrille, self.Ygrille)
        Zgrid = X + 1j * Y
        Zgrid = ma.masked_where(np.absolute(Zgrid - self.Zc) <= self.R, Zgrid)
        self.zgrid = self.H(Zgrid)
        f = self.F(Zgrid - self.Zc)  # f(z) = f(Z)
        self.Psi = f.imag
        Psi_min = np.min(self.Psi)
        Psi_max = np.max(self.Psi)
        # print('Psi_min                                   : %f '%(Psi_min))
        # print('Psi_max                                   : %f '%(Psi_max))
        self.Psi_limit[0:1] = [Psi_min, Psi_max]
        print("Psi, min,max                              :", self.Psi_limit)

    def get_contours(self, mplcont):
        """
        get  the piece of lines to calculate plt.contour
        """
        conts = mplcont.allsegs
        xline = []
        yline = []
        for cont in conts:
            if len(cont) != 0:
                for arr in cont:
                    xline += arr[:, 0].tolist()
                    yline += arr[:, 1].tolist()
                    xline.append(None)
                    yline.append(None)
        return xline, yline

    def streamlines_plot(self):
        """
        Streamline plot for Joukowski airfoil,
        alpha in degrees,
        opt_levels: 'manuel' or 'nombre'
        """
        if self.opt_levels == 'manuel':
            iso_psi = np.arange(self.levels_iso[0], self.levels_iso[1], self.levels_iso[2]).tolist()
        elif self.opt_levels == 'nombre':
            # DeltaPsi=(Psi_max-Psi_min)/(n_levels-1)
            # iso_psi=np.arange(Psi_min, Psi_max, DeltaPsi).tolist() 
            iso_psi = np.linspace(self.Psi_limit[0] + self.eps_levels, self.Psi_limit[1] - self.eps_levels,
                                  self.n_iso).tolist()
        else:
            raise ValueError("opt_levels : 'manuel' or 'nombre' ")

        plt.figure(figsize=(self.xsize, self.ysize))
        plt.title(r'iso-$\Psi$, $\alpha = $ %4.2f °' % (np.rad2deg(self.alpha)), fontsize=fz, fontweight='bold')
        psi = plt.contour(self.zgrid.real, self.zgrid.imag, self.Psi, levels=iso_psi, colors='blue')
        psi_x, psi_y = self.get_contours(psi)
        plt.plot(psi_x, psi_y)
        plt.plot(self.z_airfoil.real, self.z_airfoil.imag, label='airfoil')
        plt.plot(self.zs.real, self.zs.imag, 'ks', label='T.E.')
        plt.plot(self.z_le.real, self.z_le.imag, 'ko', label='L.E.')
        plt.plot(self.z_stagnation.real, self.z_stagnation.imag, 'bs', label='Stagnation')
        plt.plot([self.z_le.real, self.zs.real], [self.z_le.imag, self.zs.imag], 'r-', linewidth=1)
        plt.legend()
        # plt.show()

    def set_title(self, texte):
        """Set a title """
        print('#', 60 * '*')
        print('# %s' % (texte))
        print('#', 60 * '*', '\n')
