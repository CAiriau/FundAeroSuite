#!/bin/python
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Module : incompressibleFlow.LiftingLine
    ..

Python class for wing analysis
    * set_wing_parameters
    * WingAnalysis
"""
import os
import numpy as np
import matplotlib.pyplot as plt

fz = 16  # font size for captions


def set_wing_parameters():
    """
    Set default parameters
    for main parameters:
    Parameters:
    * lamb
    * S
    * b
    * l_0

    if :
        l_0 =0, one between S,b and lamb must equal to zero
    else :
        set l_0 and b
    """
    prm = {}  # Dictionnary of the parameters

    prm["l_0"] = 0              # chord at y = 0
    prm["span"] = 7             # span of the wing, be careful : span = 2 b
    prm["ny"] = 201             # number of points in y
    prm["WingAngle"] = [np.deg2rad(10), np.deg2rad(5)]  #
    #                           gamma_le: leading edge angle
    #                           gamma_te: trailing edge angle
    prm["lamb"] = 6             # lambda : wing aspect ratio
    prm["wing"] = "elliptical"  # wing shape or lift distribution, "elliptical","tapered",
    prm["U_0"] = 100.           # upstream velocity [m/s]
    prm["k_alpha"] = np.pi      # lift gradient at alpha_0
    prm["alpha"] = np.deg2rad(1.)  # geometric angle of attac
    prm["x_A"] = 0.             # x position on the wing
    prm["nFourier"] = 5         # number of elements of the Fourier's serie
    prm["S"] = 0.               # wing surface
    prm["method"] = 1           # type of method to get the coefficients
    prm["alpha0_law"] = "constante"  # alpha_0 law for wing section
    #                               "constante", "linear", "other"
    prm["alpha0_par"] = [np.deg2rad(0.), np.deg2rad(-0.)]
    prm["alphaV_law"] = 1        # alpha_V (wing twisting) law for wing section
    #  1 : constant,
    #  2 : alpha epsilon (z/b)^2
    #  3 : alpha epsilon |z/b|
    #   4 : alpha epsilon z/b
    prm["alphaV_par"] = np.deg2rad(0.)

    prm["plot_wing"] = True         # plot the wing plantform
    prm["plot_circulation"] = True

    return prm


class WingAnalysis(object):
    """
    Prandtl theory,
    Numerical approach to calculate the circulation of the wing
    """

    def __init__(self, prm, **kwargs):
        """
        Initial values of the parameters
        """
        print('#', 60 * '*')
        print('# %s' % ("Prandtl theory"))
        print('#', 60 * '*', '\n')

        self.l_0 = prm["l_0"]           # chord at y = 0
        self.span = prm["span"]         # span of the wing, be careful : span = 2 b
        self.ny = prm["ny"]             # number of points in y
        self.WingAngle = prm["WingAngle"]
        #                               gamma_le: leading edge angle
        #                               gamma_te: trailing edge angle
        self.lamb = prm["lamb"]         # lambda
        self.wing = prm["wing"]         # wing shape or lift distribution
        self.U_0 = prm["U_0"]           # upstream velocity [m/s]
        self.k_alpha = prm["k_alpha"]   # lift gradient at alpha_0
        self.alpha = prm["alpha"]       # geometric angle of attack
        self.x_A = prm["x_A"]           # x position on the wing
        self.nFourier = prm["nFourier"] # number of elements of the Fourier's serie
        self.S = prm["S"]               # wing surface
        self.b = self.span / 2          # half span

        self.method = prm["method"]     # 1 : sine integral, 2 : collocation points, 3 : analytical for tapered wing

        self.alpha0_law = prm["alpha0_law"]  # alpha_0 law for wing section  "constante", "linear", "other"
        self.alpha0_par = prm["alpha0_par"]
        self.alphaV_law = prm["alphaV_law"]
        self.alphaV_par = prm["alphaV_par"]

        self.alpha_0_elliptic = np.deg2rad(0.)

        self.plot_wing = prm["plot_wing"]  # plot the wing plantform
        self.plot_circulation = prm["plot_circulation"]
        self.xsize = 10
        self.ysize = 8
        self.pmax = 0                       # to calculate odd coefficients n = 2 p + 1

    def run_analysis(self):
        """
        To calculate the wing aerodynamic properties
        1 : elliptic wing
        2 : tapered wing (and rectangular and delta shapes)
        3 : chord law, elliptic wing : not implemented
        4 : load data : not implemented
        """
        print("\n", "*" * 80, "\n")
        self.Lambda_LE = self.WingAngle[0]
        self.Lambda_TE = self.WingAngle[1]
        self.b = self.span / 2  # half span

        if self.nFourier % 2 == 0:
            self.pmax = int(self.nFourier / 2 - 1)
        else:
            self.pmax = int((self.nFourier - 1) / 2)

        self.wing_plantform()   # define the wing shape
        self.alpha0_airfoil()   # define the zero lift incidence
        self.alphaV_wing()      # define the wing twisting
        self.elliptical_wing()  # reference values for elliptical wing
        print("p max, n Fourier                 : %3i %3i" % (self.pmax, self.nFourier))
        self.wing_solve_An()    # An coefficients
        self.wing_aerodynamics()  # Aerodynamic coefficients and circulation

        if self.plot_wing or self.plot_circulation: plt.show()

    def wing_solve_An(self):
        """
        to solve the An coefficients of Prandtl equation
        """
        print("AERODYNAMICS:")

        mu_theta = self.k_alpha * self.chord / (4 * self.b)
        b_theta = mu_theta * (self.alpha - self.alpha_0 + self.alpha_V)
        sint = np.sin(self.theta)

        if self.method == 2:
            print('*** collocation method ***')
            Mat = np.zeros((self.ny - 2, self.nFourier))
            u = np.arange(1, self.ny - 1)
            for i in range(self.nFourier):
                Mat[:, i] = (1 + (mu_theta[u] * (i + 1)) / np.sin(self.theta[u])) * np.sin((i + 1) * self.theta[u])
                # Mat[:,i] = (sint[u]+mu_theta[u]*(i+1))*np.sin((i+1)*self.theta[u])
            # print('size of M : ',M.shape)
            # print('size of b : ',b_theta[u].shape)

            if self.nFourier == int(self.ny - 2):
                print("Ns = Nc-2, square matrix : direct inversion")
                # self.An=np.linalg.solve(Mat,b_theta[u]*sint[u])
                self.An = np.linalg.solve(Mat, b_theta[u])

                # E=max(np.abs(np.matmul(Mat,self.An)-b_theta[u]*sint[u]))
                E = max(np.abs(np.matmul(Mat, self.An) - b_theta[u]))
                print("Inversion error                  : %2.8e " % (E))

            else:
                print("Ns != Nc-2 , least square approach : Ns = %i , Nc = %i " % (self.nFourier, self.ny))
                # self.An,res,rank,s = np.linalg.lstsq(Mat,b_theta[u]*sint[u],rcond=-1)
                self.An, res, rank, s = np.linalg.lstsq(Mat, b_theta[u], rcond=-1)
                print("res  = ", res)


        elif self.method == 1:
            print("*** Sine integral method ***")
            RHS = np.zeros((self.nFourier))
            Mat = np.zeros((self.nFourier, self.nFourier))
            for i in range(self.nFourier):
                b1_theta = b_theta * sint * np.sin((i + 1) * self.theta)
                b2_theta = -np.trapz(b1_theta, x=self.theta)  # numerical integration
                RHS[i] = b2_theta  # build the vector V
                for j in np.arange(self.nFourier):
                    Fj = (sint + mu_theta * (j + 1)) * np.sin((j + 1) * self.theta)
                    Si = np.sin((i + 1) * self.theta)
                    Vec_ij = Fj * Si
                    Mat[i, j] = -np.trapz(Vec_ij, x=self.theta)  # numerical integration and build the Matrix M

            self.An = np.linalg.solve(Mat, RHS)
            E = max(np.abs(np.matmul(Mat, self.An) - RHS))
            print("Inversion error                  : %2.8e " % (E))

        elif self.method == 3 and self.wing == "tapered":
            print("*** Analytical method for tapered wing ***")
            print("True if alpha = constant, no twisting")
            self.An = np.zeros(self.nFourier)
            self.calcul_An()


        else:
            raise NameError("bad choice of method to solve An")

        for i in range(self.nFourier):
            print("A[%i]                             : %2.6e" % (i + 1, self.An[i]))

    def wing_aerodynamics(self):
        """
        Aerodynamics coefficients  and circulation
        """
        # resolution of the circulation

        S = np.zeros((self.ny, self.nFourier))
        for j in range(self.nFourier):
            S[:, j] = np.sin((j + 1) * self.theta)
        self.Gamma = 4 * self.b * self.U_0 * np.matmul(S, self.An)
        if self.plot_circulation:
            plt.figure(figsize=(self.xsize, self.ysize))
            plt.title(r'Circulation $\Gamma$', fontsize=fz, fontweight='bold')
            plt.plot(self.y, self.Gamma, 'b-')

        # determine wing coefficients
        self.CL = np.pi * self.lamb * self.An[0]
        u = np.arange(0, self.nFourier)
        v = np.arange(1, self.nFourier + 1)
        # print("u = ",u,self.nFourier)
        # print("v = ",v)
        self.CD = np.pi * self.lamb * np.sum(v * self.An[u] ** 2)
        CD = self.CL ** 2 / (self.lamb * np.pi) + np.pi * self.lamb * np.sum(
            np.arange(2, self.nFourier + 1) * self.An[np.arange(1, self.nFourier)] ** 2)
        self.Cmx = -np.pi / 4 * self.lamb * self.An[1]
        self.Cmx = -np.pi / 4 * self.lamb * self.An[1]
        u = np.arange(1, self.nFourier - 1)
        self.Cmy = np.pi / 4 * self.lamb * sum((2 * u + 1) * self.An[u] * self.An[u + 1])
        if abs(self.CD) <= 1e-6:
            print('CD = 0')
            self.f = float("inf")
        else:
            self.f = self.CL / self.CD

        print("alpha                            : %2.7f ° " % (np.rad2deg(self.alpha)))
        print("Max Gamma                        : %2.7f m^2/s" % (max(self.Gamma)))
        print("4 b U0 A1                        : %2.7f m^2/s" % (4 * self.b * self.U_0 * self.An[0]))
        print("CL                               : %2.7f " % (self.CL))
        print("CD                               : %1.7e " % (self.CD))
        print("CD                               : %1.7e " % (CD))

        print("f                                : %3.8f" % (self.f))
        print("Cmx                              : %2.6e" % (self.Cmx))
        print("Cmy                              : %2.6e" % (self.Cmy))

    def elliptical_wing(self):
        """
        Reference value for an elliptic wing
        """
        self.A1_elliptic = 2 * (self.alpha - self.alpha_0_elliptic) / (self.lamb + 2)
        self.CL_elliptic = np.pi * self.lamb * self.A1_elliptic
        self.CD_elliptic = self.A1_elliptic * self.CL_elliptic
        print("A1_elliptic = ", self.A1_elliptic, "for alpha =  ", np.rad2deg(self.alpha))
        if abs(self.A1_elliptic) <= 1e-9:
            print('A1 = 0')
            self.f_elliptic= float("inf")
        else:
            self.f_elliptic = 1 / self.A1_elliptic
        self.alpha_i_elliptic = np.rad2deg(-self.A1_elliptic)
        self.Gamma0_elliptic = 4 * self.b * self.U_0 * self.A1_elliptic
        print("-" * 80)
        print("EQUIVALENT ELLIPTICAL WING (ANALYTICAL) :")
        print("A1 elliptical                : %2.6e " % (self.A1_elliptic))
        print("alpha_i                      : %2.6e " % (np.rad2deg(-self.A1_elliptic)))
        print("CL elliptical                : %2.4f " % (self.CL_elliptic))
        print("CD elliptical                : %1.7f " % (self.CD_elliptic))
        print("f  elliptical                : %2.8f " % (self.f_elliptic))
        print("lambda/(2+lambda)            : %2.6f " % (self.lamb / (2 + self.lamb)))
        print("Gamma0 elliptical            : %2.6f " % (self.Gamma0_elliptic))
        print("-" * 80)

    def alphaV_wing(self):
        """
        Wing twisting law
        """
        if self.alphaV_law == 1:
            self.alpha_V = np.ones(self.ny) * self.alphaV_par
        elif self.alphaV_law == 2:
            self.alpha_V = self.y ** 2 / self.b ** 2 * self.alphaV_par
        elif self.alphaV_law == 3:
            self.alpha_V = abs(self.y) / self.b * self.alphaV_par
        elif self.alphaV_law == 4:
            self.alpha_V = self.y / self.b * self.alphaV_par
        else:
            self.alpha_V = np.zeros(self.ny)

    def alpha0_airfoil(self):
        """
        Angle of zero lift law with respect to y
        """
        if self.alpha0_law == "constante":
            self.alpha_0 = np.ones(self.ny) * self.alpha0_par[0]
        else:
            self.alpha_0 = abs(self.y) * self.alpha0_par[1]

    def wing_plantform(self):
        """
        Wing plantform definition
        """
        print("WING DATA : ")

        if self.wing == "elliptical":
            """ Elliptic wing """
            print("=" * 5 + "Elliptical wing " + "=" * 5)

            if self.span * self.l_0 != 0:
                print("l_0 and b known")
                self.S = np.pi / 2 * self.l_0 * self.b
                self.lamb = self.span ** 2 / self.S

            elif self.lamb * self.b != 0:
                print("lambda and b are known")
                self.S = self.span ** 2 / self.lamb

            elif self.b * self.S != 0:
                print("b and S are known")
                self.lamb = self.span ** 2 / self.S

            elif self.lamb * self.S != 0:
                print("lambda and S are known")
                self.span = np.sqrt(self.lamb * self.S)
                self.b = self.span / 2
            else:
                raise NameError("error : verify the geometrical parameters")

            self.l_0 = 2 * self.S / (np.pi * self.b)

            self.span_coordinates()

            print("S (analytic)                     : %4.4f m^2" % self.S)
            print("b (analytic)                     : %2.4f m" % self.b)
            print("c_mean (analytic)                : %2.4f m" % (self.S / self.span))
            print("lambda (analytic)                : %2.4f " % self.lamb)

            self.chord = self.l_0 * np.sqrt(1 - (self.y / self.b) ** 2)
            self.LE, self.TE = self.chord / 4, -3 * self.chord / 4
            self.S_num = np.trapz(self.chord, x=self.y)
            self.lamb_num = self.span ** 2 / self.S_num
            meanchord = self.S_num / self.span

            print("S num                            : %4.4f m^2" % self.S_num)
            print("c_mean num                       : %2.4f m" % meanchord)
            print("lambda num                       : %2.4f " % self.lamb_num)

        if self.wing == "tapered":
            """
            Tapered wing
            """
            print("=" * 5 + "Tapered wing " + "=" * 5)

            tmp = np.tan(self.Lambda_LE) + np.tan(self.Lambda_TE)

            if self.span * self.l_0 != 0:
                print("l_0 and b known")
                chordtip = self.l_0 - tmp * self.b
                self.S = self.b * (self.l_0 + chordtip)
                self.lamb = self.span ** 2 / self.S

            elif self.lamb * self.b != 0:
                print("lambda and b are known")
                self.S = self.span ** 2 / self.lamb

            elif self.b * self.S != 0:
                print("b and S are known")
                self.lamb = self.span ** 2 / self.S

            elif self.lamb * self.S != 0:
                print("lambda and S are known")
                self.span = np.sqrt(self.lamb * self.S)
                self.b = self.span / 2
            else:
                raise ("error : verify the geometrical parameters")

            self.span_coordinates()
            self.l_0 = (self.S / self.b + self.b * tmp) / 2
            chordtip = self.l_0 - tmp * self.b  # chord tip
            self.LE = self.l_0 / 2 - abs(np.tan(self.Lambda_LE)) * abs(self.y)
            self.TE = -self.l_0 / 2 + abs(np.tan(self.Lambda_TE)) * abs(self.y)
            self.chord = self.LE - self.TE
            # self.S=self.b*(self.l_0+chordtip)
            self.lamb_num = self.span ** 2 / self.S
            meanchord = self.S / self.span
            self.S_num = np.trapz(self.chord, x=self.y)

            print("tip chord                        : %4.7f m" % chordtip)
            print("root chord                       : %4.7f m" % self.l_0)
            print("tip/root chord                   : %4.4f " % (chordtip / self.l_0))
            print("lambda num                       : %2.4f " % self.lamb_num)
            print("S num, S an.                     : %4.4f m^2 , %4.4f m^2" % (self.S_num, self.S))
            print("mean chord num                   : %4.4f m" % meanchord)
            print("mean chord num                   : %4.4f m" % ((self.l_0 + chordtip) / 2))

            self.mu0 = self.k_alpha * self.l_0 / (4 * self.b)
            self.mu0 = self.k_alpha * (1 / (2 * self.lamb) + tmp / 8)
            self.mu1 = -self.k_alpha / 4 * tmp
            print("mu0                              : %2.8e " % self.mu0)
            print("mu1                              : %2.8e " % self.mu1)

        print("root chord                       : %2.4f m" % self.l_0)
        print("b                                : %4.7f m" % self.b)
        print("lambda                           : %4.7f " % self.lamb)
        print("U0                               : %4.7f m/s" % self.U_0)

        if self.plot_wing:
            plt.figure(figsize=(self.xsize, self.xsize))
            plt.title(r'Wing plantform', fontsize=fz, fontweight='bold')
            plt.plot(self.y, self.LE, 'b-')
            plt.plot(self.y, self.TE, 'b-')
            if abs(self.TE[-1] - self.LE[-1]) >= 1e-10:
                plt.plot([self.b, self.b], [self.TE[-1], self.LE[-1]], 'b-')
                plt.plot([-self.b, -self.b], [self.TE[-1], self.LE[-1]], 'b-')
            plt.axis("equal")

    def B_coef(self, n, m):
        """
        coefficients for a rectangular wing
        """
        return -(8 * (m * n) / (((m + n) ** 2 - 1) * ((m - n) ** 2 - 1) * np.pi))

    def C_coef(self, n, m):
        """
        coefficients for a tapered wing
        """
        return (4 * (m ** 2 + n ** 2 - 1) * (-1) ** ((m + n) / 2) / (((m + n) ** 2 - 1) * ((m - n) ** 2 - 1) * np.pi))

    def calcul_An(self):
        """
        tapered or rectangular or delta wing 
        Only odd coefficients are solved
        output A_{2p+1} (A[2p] in python)
        """
        # construction de la matrice
        m = 2 * np.arange(self.pmax + 1) + 1

        Mat = np.zeros((self.pmax + 1, self.pmax + 1))
        Mat = self.mu0 * np.diag(m)

        for i in range(self.pmax + 1):
            n = 2 * i + 1
            Mat[i, :] += self.B_coef(n, m) + self.mu1 * m * self.C_coef(n, m)
        # mat[0,:] = [B_coef(1,1)+mu0+mu1*C_coef(1,1), B_coef(1,3)+3*mu1*C_coef(1,3), B_coef(1,5)+5*mu1*C_coef(1,5)]
        # mat[1,:] = [B_coef(3,1)+mu1*C_coef(3,1), B_coef(3,3)+3*mu0+3*mu1*C_coef(3,3), B_coef(3,5)+5*mu1*C_coef(3,5)]
        # mat[2,:] = [B_coef(5,1)+mu1*C_coef(5,1), B_coef(5,3)+3*mu1*C_coef(5,3), B_coef(5,5)+5*mu0+5*mu1*C_coef(5,5)]

        rhs = np.zeros((self.pmax + 1))
        # print("m odd =",m )
        rhs[0] = self.mu0
        rhs += self.mu1 * self.C_coef(1, m)
        rhs *= self.alpha
        # print("mat  = ",Mat)
        # print("rhs  = ",rhs)
        An = np.linalg.solve(Mat, rhs)
        self.An[m - 1] = An
        E = max(np.abs(np.matmul(Mat, An) - rhs))
        print("Inversion error                  : %2.8e " % (E))

    def span_coordinates(self, distribution="cosine"):
        """
        points distribution in the spanwise direction
        """
        if distribution == "cosine":
            self.theta = np.linspace(0, np.pi, self.ny)
            self.y = -self.b * np.cos(self.theta)
        else:
            self.y = np.linspace(-self.b, self.b, self.ny)
            self.theta = np.arccos(-self.y / self.b)
        self.y_half = self.y[int((len(self.y) - 1) / 2):]

        # print('len : ',len(self.y),len(self.y_half))
