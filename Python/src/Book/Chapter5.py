#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 5 exercises
    ..
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import os

# from CompressibleFlow.fonctions import *
from Tools.misc import set_title, set_question

from IncompressibleFlow.LiftingLine import set_wing_parameters, WingAnalysis
from IncompressibleFlow.Prandtl import Prandtl_E_IGE, Prandtl_E_OGE, vortex_velocity

fz = 16

def Exercice5_1():
    """
    Elliptical optimal wing
    """
    set_title("Elliptical optimal wing")
    par = set_wing_parameters()
    s = WingAnalysis(par)
    s.plot_wing = True
    s.plot_circulation = False
    s.method = 1
    s.S = 22.48
    s.span = 11.25
    s.lamb = 0.
    # validated
    s.run_analysis()


def Exercice5_2():
    """
    Twisted elliptical wing
    """
    set_title("Twisted elliptical wing")
    set_question("1- case 2")
    par = set_wing_parameters()
    s = WingAnalysis(par)
    s.plot_wing = False
    s.plot_circulation = False
    s.method = 1
    s.S = 22.48
    s.span = 11.25
    s.lamb = 0.
    s.alphaV_law = 2
    s.alphaV_par = np.deg2rad(1.)
    s.run_analysis()


def Exercice5_3():
    """
    Flaps on a elliptical wing
    """
    set_title("Flaps on a elliptical wing")
    print("not implemented in python yet")


def Exercice5_4():
    """
    Wing aerodynamic center
    """
    set_title(" Wing aerodynamic center")
    """
    The momentum are calculated at the leading edge point, on the symmetry axis
    """

    CmF_profile = -0.05      # CmF of the airfoil (we can define a law with respect to  theta)
    Ny = 101                # Number of points to discretize the wing
    Lambda = np.deg2rad(5)  # Angle to define the airfoil aerodynamic center on the wing
    slope = np.tan(Lambda)
    b = 7                   # wing half-span
    x0 = -0.25              # position of the aerodynamic center of the airfoil on the symmetry axis of the wing
    #                       (wing root)
    Ns = 3                  # number of Fourier mode for the circulation serie
    lamb = 7                # wing aspect ration
    alpha = np.deg2rad(1)   # constant wing incident alpha-alpha_0
    # U0        = 50        # velocity in  m/s
    # rho       = 1.3       # air density
    xA = -0.5               # position /chord of the aerodynamic center at the root wing leading edge

    A = np.zeros(Ns)
    T = np.zeros(Ns)

    def x_F_profile(theta):
        """
        to define the aerodynamic center position for each airfoil
        not used today
        """
        return x0 + slope * b * abs(np.cos(theta))

    theta = np.linspace(0, np.pi, Ny)

    # only calculated on a elliptical wing
    S = 4 * b ** 2 / lamb
    Lmoy = 2 * b / lamb
    L0 = 4 / np.pi * Lmoy
    x0 *= L0
    xA *= L0

    A[0] = 2 * alpha / (2 + lamb)  # elliptical wing value
    print("b                    : ", b)
    print("lambda               : ", lamb)
    print("S                    : ", S)
    print("Lmoy                 : ", Lmoy)
    print("L0                   : ", L0)
    print("alpha (°)            : ", np.rad2deg(alpha))
    print("A1                   : ", A[0])
    print("Lambda (°)           : ", np.rad2deg(Lambda))
    print("xA/L0                : ", xA / L0)
    print("CmF_profile           : ", CmF_profile)
    print("x0/L                 : ", x0 / L0)

    # here either  A_1,A_3, ... are calculated  or I'am taking the elliptical wing values
    # A_3, etc ... are null
    # analytical relationships
    print("\nANALYCAL METHOD :\n")
    CmA = 32 / (3 * np.pi ** 2) * CmF_profile + (np.pi * (xA - x0) - 4 * b / 3 * slope) * lamb / Lmoy * A[0]
    CmA_T1 = 32 / (3 * np.pi ** 2) * CmF_profile
    xF = x0 + 4 * b / (3 * np.pi) * slope

    print("xF/L0                : ", xF / L0)
    print("CmA/CmF_airfoil      : ", CmA / CmF_profile)
    print("CmA T1               : ", CmA_T1)

    print("\nNUMERICAL METHOD :\n")

    #  Integral calculus

    def law_chord(theta):
        """ chord law """
        return L0 * np.sin(theta)

    def law_xF(theta):
        """aerodynamic center position """
        return x0 + b * abs(np.cos(theta)) * slope

    def term_profile_CmF(theta):
        return law_chord(theta) ** 2 * np.sin(theta)

    def term_profile_xF(theta, n):
        return law_xF(theta) * np.sin(theta) * np.sin(n * theta)

    T1, res1 = integrate.quad(term_profile_CmF, 0, np.pi)
    T2 = xA / Lmoy * np.pi * lamb * A[0]
    for n in range(Ns):
        T[n], res = integrate.quad(term_profile_xF, 0., np.pi, args=(n + 1,))
        # print("T[%i]              : %f"%(n,T[n]))
        T[n] *= -2 * lamb / Lmoy * A[n]
        print("term with  F point      : ", "n= ", n + 1, " : ", T[n], "residu : ", res)

    T1 *= CmF_profile / (2 * Lmoy ** 2)
    print("airfoil term         : ", T1, "residu : ", res1)
    print("Point A term         : ", T2)
    print("CmA/CmF_airfoil      : ", (T1 + T2 + np.sum(T)) / CmF_profile)


def Exercice5_5():
    """
    Airbus A380-800 vortex wake
    """
    set_title("Airbus A380-800 vortex wake")
    rho0 = 1.2      # air density at the ground level [kg/m^3]
    rho = 0.35      # air density at the flight level [kg/m^3]
    Ud = 80.        # take-off speed [m/s]
    U = 260.        # cruise velocity [m/s]
    S = 840         # wing surface [m^2]
    b = 40          # halp-span [m]
    W = 4.e6        # weight supported by the mean wing (wing load) [N]

    set_question('1 - wing aspect ratio, chord c_0')
    lamb = 4 * b ** 2 / S
    c_0 = 8 * b / (np.pi * lamb)
    print("wing aspect ratio                        : %2.2f " % (lamb))
    print("Wing root chord                          : %2.2f m " % (c_0))

    set_question('2 - Cz at take-off ')
    print("w = ", W)
    Czdec = 2 * W / (rho0 * Ud ** 2 * S)
    print("Cz take-off                              : %2.3f " % (Czdec))
    alphadec = (1 + 2 / lamb) * Czdec / (2 * np.pi)
    print("Equivalent incidence                     : %2.3f " % (np.rad2deg(alphadec)))

    print("Cruise flight : ")
    Cz = 2 * W / (rho * U ** 2 * S)
    alpha = (1 + 2 / lamb) * Cz / (2 * np.pi)
    print("Cz                                       : %2.3f " % (Cz))
    print("Equivalent incidence                     : %2.3f " % (np.rad2deg(alpha)))

    set_question("3 - Circulation")

    def circulation(U, Cz):
        return (4 * b * U * Cz) / (np.pi * lamb)

    Gamma0 = circulation(Ud, Czdec)
    print("Circulation at take-off                  : %2.3f " % (Gamma0))
    Gamma = circulation(U, Cz)
    print("Circulation at cruise                    : %2.3f " % (Gamma))

    set_question("6 - vortex wing characterization")
    y_T = np.pi * b / 4
    R_T_over_b = np.sqrt(2. / 3. - np.pi ** 2 / 16)
    print("y_T                                      : %2.3f m  " % (y_T))
    print("R_T / b                                  : %2.3f m  " % (R_T_over_b))

    set_question("7 - descent speed of the vortex wing")
    v_d = Gamma / (np.pi ** 2 * b)
    print("descent speed of a vortex wing          : %2.3f m/s " % (v_d))

    set_question("Cxi at take-off")
    Cxi = Czdec ** 2 / (np.pi * lamb)
    print("Cxi at take-off (without ground), lift-to-drag ratio   : %2.4f  %2.2f  " % (Cxi, Czdec / Cxi))

    def f_DeltaCxi(h):
        return Gamma0 ** 2 / (np.pi * 2 * Ud ** 2 * S) * np.log(1 + y_T ** 2 / h ** 2)

    print("Delta Cxi, h/y_T                         : %2.4f %2.2f " % (f_DeltaCxi(y_T / 2), 1 / 2))
    print("Delta Cxi, h/y_T                         : %2.4f %2.2f " % (f_DeltaCxi(y_T / 4), 1 / 4))


def Exercice5_6():
    """
    Wind tunnel corrections
    """
    set_title("Wind tunnel corrections")
    # rho0    = 1.2         # air density at the flight level [kg/m^3]
    U0 = 30.                # speed in the wind tunnel [m/s]
    lamb = 6                # wing aspect ratio
    b = 0.5                 # half-span [m]
    R = 1.0                 # wind tunnel radius [m]

    Sv = np.pi * R ** 2

    set_title("optimal elliptical wing")
    par = set_wing_parameters()
    s = WingAnalysis(par)
    s.plot_wing = False
    s.plot_circulation = False
    s.method = 1
    s.alpha = np.deg2rad(4.)
    s.S = 0.
    s.U_0 = U0
    s.span = 2 * b
    s.lamb = lamb
    s.run_analysis()

    Sa = s.S
    CL = s.CL_elliptic
    CD = s.CD_elliptic
    x0 = np.pi * b / 4

    print("x0                           : %f" % (x0))
    print("R^2/x0                       : %f" % (R ** 2 / x0))
    set_question("delta_alpha_i")
    delta_alpha_i = Sa / (8 * Sv) * CL
    print("delta Alpha_i                : %f °" % (np.rad2deg(delta_alpha_i)))
    print("delta Cx_i                   : %f" % (CL * delta_alpha_i))
    print("delta Cx_i/Cx_i              : %f" % (CL * delta_alpha_i / CD))
    print("Sa/Sv                        : %f" % (Sa / Sv))


def Exercice5_7():
    """
    Horse-shoe vortex
    """
    set_title("Horse-shoe vortex")

    b = 1
    Gamma = -np.pi
    eta = np.linspace(0.15, 4, 201)
    w = []
    for x in eta:
        w.append(vortex_velocity([x, 0], 0, -b, b, Gamma))

    plt.figure(figsize=(12, 8))
    plt.title(r'$\tilde w$', fontsize=fz, fontweight='bold')
    plt.xlabel(r"$\eta=x/b$")
    plt.xlabel(r"$w/w_{ref}$")
    plt.plot(eta, w)
    plt.grid()
    plt.show()

    print("x/b = 4  : w/wref = ", vortex_velocity([4 * b, 0], 0, -b, b, Gamma))
    print("x/b = 0  : w/wref = ", vortex_velocity([0, 0], 0, -b, b, Gamma))


def Exercice5_8():
    """
    Formation flying
    """
    set_title("Formation flying")

    test = False
    b = 5
    Gamma = 4 * np.pi
    r = 0.5
    d = r * b
    x = 10

    if test:
        w = vortex_velocity([x, 0], 0, -b, b, Gamma)
        print(" w = %f  " % (w))

    set_question("1- isolated bird")
    w0 = vortex_velocity([0, 0], 0, -b, b, Gamma)
    print("velocity at the wing center         : ", w0)

    set_question("2- three bird flying on a front line")
    w_centre = vortex_velocity([0, 0], 0, -3 * b, 3 * b, Gamma)
    w_right = vortex_velocity([0, 2 * b], 0, -3 * b, 3 * b, Gamma)
    w_left = vortex_velocity([0, -2 * b], 0, -3 * b, 3 * b, Gamma)

    print("W_C/w0 = %f, W_D/w0 = %f W_G/w0 = %f " % (w_centre / w0, w_right / w0, w_left / w0))

    set_question("3- three birds flying in V")

    r = np.linspace(0.4, 0.43, 21)

    r = [0.40, 0.418768, 0.42]
    Wc, Wd = [], []
    for R in r:
        d = R * b
        # bird at the center
        w_bird = {}
        w_centre = vortex_velocity([0, 0], 0, -b, b, Gamma)
        w_right = vortex_velocity([0, 0], d, b, 3 * b, Gamma)
        w_left = vortex_velocity([0, 0], d, -3 * b, -b, Gamma)
        w_bird["centre"] = (w_centre + w_right + w_left) / w0

        # bird on the right
        w_centre = vortex_velocity([d, 2 * b], 0, -b, b, Gamma)
        w_right = vortex_velocity([d, 2 * b], d, b, 3 * b, Gamma)
        w_left = vortex_velocity([d, 2 * b], d, -3 * b, -b, Gamma)
        w_bird["right"] = (w_centre + w_right + w_left) / w0

        # bird on the left
        w_centre = vortex_velocity([d, -2 * b], 0, -b, b, Gamma)
        w_right = vortex_velocity([d, -2 * b], d, b, 3 * b, Gamma)
        w_left = vortex_velocity([d, -2 * b], d, -3 * b, -b, Gamma)
        w_bird["left"] = (w_centre + w_right + w_left) / w0

        wmean = (w_bird["right"] + w_bird["centre"] + w_bird["left"]) / 3
        Wc.append(w_bird["centre"])
        Wd.append(w_bird["right"])

        print("%f %f %f %f %f " % (R, w_bird["left"], w_bird["centre"], w_bird["right"], 45 * wmean))

    print("Optimal angle  = %f °" % (np.rad2deg(np.arctan(1 / 0.418768))))
    # intersection = 0.418768
    plt.figure()
    plt.plot(r, Wc)
    plt.plot(r, Wd)
    plt.show()


def Exercice5_9():
    """
    Wake vortex with ground effect
    """
    set_title(" Wake vortex with ground effect")

    def fun(x, k):
        """
        trajectory or streamline
        """
        return k * x / np.sqrt((k ** 2 + 1) * x ** 2 - k ** 2)

    k = [0.5, 1, 2]
    npt = 501
    a = 3

    plt.figure(figsize=(12, 6))
    plt.title(r'Trajectory', fontsize=fz, fontweight='bold')
    for K in k:
        kmin = K / np.sqrt(K ** 2 + 1)
        print(" asymptote     : ", kmin)
        x = np.linspace(kmin * 1.01, a, npt)
        y = fun(x, K)
        plt.plot(x[y <= a], y[y <= a], label=" k = %2.1f" % K)
        plt.plot(-x[y <= a], y[y <= a], label=" k = %2.1f" % K)
        plt.plot([-a, a], [kmin, kmin], '--')
        plt.plot([kmin, kmin], [0, a], '--')
        plt.plot([-kmin, -kmin], [0, a], '--')

    plt.xlim([-a, a])
    plt.axis('equal')
    plt.grid()
    plt.legend()
    plt.show()


def Exercice5_10():
    """
    Numerical solution of the Prandt equation : general case
    implementation using a  class
    (otherwise see exercise 12)
    """

    set_title("Numerical solution of the Prandt equation : general case")

    stop = 10  # question number where the function is stopped
    # to run all set to 10

    Ny = 101  # number of point for the method I
    plot_circulation = False
    gradient_Cz = True

    # to get the  gradient of Cz / alpha  instead of  this of  Cz (Cz is linear with alpha)
    if gradient_Cz:
        alpha_ref = np.deg2rad(1)  # incidence of reference (1 degree)
    else:
        alpha_ref = np.deg2rad(4.1)

    lamb_ref = 7            # wing aspect ratio
    b_ref = 7               # half-span
    U0_ref = 101            # reference velocity
    Ns = 5                  # number of Fourier mode
    k_ref = np.pi           # half lift  gradient
    plot_incidence = False  # incidence effect
    show_results = True     # to display CL, CD
    plot_wing_aspect_ratio = True # to test wing aspect ratio effect

    set_question("1 - reference case : elliptical wing")

    Cz, Cxi, b, Gamma0 = 0.35, 0.00557, 7, 45
    lamb = Cz ** 2 / (np.pi * Cxi)
    A1 = Cz / (np.pi * lamb)
    alpha = A1 * (2 + lamb) / 2
    U0 = Gamma0 * Cz / (4 * b * Cxi)
    U0bis = Gamma0 / (4 * b * A1)
    S = 4 * b ** 2 / lamb
    print("span                     : %3.3f m" % (2 * b))
    print("wing aspect ratio        : %3.3f " % (lamb))
    print("wing surface             : %3.3f m^2" % (S))
    print("velocity                 : %3.3f m/s \t %3.3f m/s" % (U0, U0bis))
    print("incidence                : %3.3f °" % (np.rad2deg(alpha)))

    if stop == 1: return

    if plot_circulation:
        fig = plt.figure(10)
        plt.title(r'circulation law', fontsize=fz, fontweight='bold')

    set_question("2 - Method I : integration of  sin m theta")

    par = set_wing_parameters()
    s = WingAnalysis(par)
    s.plot_wing = False
    s.plot_circulation = False
    s.method = 1  # sine approach
    s.S = 0
    s.span = 2 * b_ref
    s.lamb = lamb_ref
    s.alpha = alpha_ref
    s.wing = "elliptical"
    s.U_0 = U0_ref
    s.nFourier = Ns
    s.ny = Ny
    s.run_analysis()

    if stop == 2: return

    set_question("3 - Method II : collocation points")

    s.An = 0
    s.S, s.l_0 = 0, 0
    s.method = 2
    s.ny = Ns + 2
    s.run_analysis()

    if stop == 3: return

    # =============================================================
    def wing_aerodynamics(lamb, b, alpha, U0, Angle=[np.deg2rad(0), np.deg2rad(0)], plantform="rectangular", verif=False):
        """
        Calculus of rectangular, tappered or delta wing
        """

        # **** see book equations ***
        # it is also in the method n° 3 of the class
        # the next 40 lines could be deleted  (up to  "possibly deleted") 

        def B_coef(n, m):
            return -np.float(8 * (m * n) / (((m + n) ** 2 - 1) * ((m - n) ** 2 - 1) * np.pi))

        def C_coef(n, m):
            return np.float(4 * (m ** 2 + n ** 2 - 1) * (-1) ** (int((m + n) / 2)) / (
                        ((m + n) ** 2 - 1) * ((m - n) ** 2 - 1) * np.pi))

        def calcul_An(mu0, mu1, alpha):
            """
            Matrix to invert for rectangular or tapered wing
            only odd coefficients
            """
            # build the matrix
            mat = np.zeros((3, 3))
            mat[0, :] = [B_coef(1, 1) + mu0 + mu1 * C_coef(1, 1), B_coef(1, 3) + 3 * mu1 * C_coef(1, 3),
                         B_coef(1, 5) + 5 * mu1 * C_coef(1, 5)]
            mat[1, :] = [B_coef(3, 1) + mu1 * C_coef(3, 1), B_coef(3, 3) + 3 * mu0 + 3 * mu1 * C_coef(3, 3),
                         B_coef(3, 5) + 5 * mu1 * C_coef(3, 5)]
            mat[2, :] = [B_coef(5, 1) + mu1 * C_coef(5, 1), B_coef(5, 3) + 3 * mu1 * C_coef(5, 3),
                         B_coef(5, 5) + 5 * mu0 + 5 * mu1 * C_coef(5, 5)]

            rhs = np.zeros((3))
            rhs[0] = alpha * (mu0 + mu1 * C_coef(1, 1))
            rhs[1] = alpha * mu1 * C_coef(1, 3)
            rhs[2] = alpha * mu1 * C_coef(1, 5)
            print("mat  = ", mat)
            print("rhs  = ", rhs)
            An = np.linalg.solve(mat, rhs)
            return An

        if verif:
            print("Coefficient of matrices Bnm and Cnm :")
            for p in range(3):
                for q in range(3):
                    n, m = 2 * p + 1, 2 * q + 1
                    print(
                        "B[%i,%i] = %f, C[%i,%i] = %f %i" % (n, m, B_coef(n, m), n, m, C_coef(n, m), int((m + n) / 2)))
        s.nFourier = Ns
        print("wing shape   : ", plantform)

        if plantform == "elliptical":
            An = np.zeros(3)
            mu0 = 2 * k_ref / (np.pi * lamb)
            mu1 = 0
            An[0] = 2 / (lamb + 2) * alpha
            s.wing = "elliptical"
        else:
            mu0 = k_ref / 8 * (4 / lamb + np.sum(np.tan(Angle)))
            mu1 = -k_ref / 4 * np.sum(np.tan(Angle))
            An = calcul_An(mu0, mu1, alpha)
            s.wing = "tapered"

        Gamma_max = 4 * b * U0 * (An[0] - An[1] + An[2])
        print("An     : ", An)
        print("lamb   : ", lamb)
        print("mu0    : ", mu0)
        print("mu1    : ", mu1)
        print("Gamma_max : ", Gamma_max)

        print("alpha  : ", np.rad2deg(alpha))
        CLa, CDa = np.pi * lamb * An[0], np.pi * lamb * np.sum((2 * np.arange(0, 3) + 1) * An ** 2)
        print("CL     : ", CLa)
        print("CD     : ", CDa)

        # *** possibly deleted above 

        #  METHOD I

        s.plot_wing = False
        s.WingAngle = Angle
        s.method = 1  # sine approach
        s.S, s.l_0 = 0, 0
        s.lamb, s.span = lamb, 2 * b
        s.alpha, s.U_0 = alpha, U0
        s.nFourier = Ns
        s.ny = Ny
        s.run_analysis()

        CL, CD = s.CL, s.CD

        if plot_circulation: plt.plot(s.y, s.Gamma, label=plantform)

        #  METHOD II

        s.method = 2  # collocation
        s.S, s.l_0 = 0, 0
        s.lamb, s.span = lamb, 2 * b
        s.alpha, s.U_0 = alpha, U0
        s.ny = Ns + 2
        s.run_analysis()

        #  METHOD III

        if plantform != "elliptical":
            s.method = 3
            s.S, s.l_0, s.b = 0, 0, 0
            s.lamb, s.span = lamb, 2 * b
            s.alpha, s.U_0 = alpha, U0
            s.WingAngle = Angle
            s.nFourier = Ns
            s.ny = Ny
        s.run_analysis()

        print("CL= %2.8f \t CL_a= %2.8f" % (CL, CLa))
        print("CD= %2.8f \t CD_a= %2.8f" % (CD, CDa))

        return CL, CD, CLa, CDa

    # =============================================================

    CL, CD, CLa, CDa = wing_aerodynamics(lamb_ref, b_ref, alpha_ref, U0_ref, Angle=[0, 0], plantform="elliptical")

    if stop == 3: return

    set_question("4 - rectangulat wing")

    lamb, U0, alpha = 7, 101, np.deg2rad(4.1)
    WingAngle = [0, 0]
    CL, CD, CLa, CDa = wing_aerodynamics(lamb, b, alpha_ref, U0_ref, Angle=WingAngle, plantform="rectangular")

    s.plot_wing = False
    s.WingAngle = WingAngle
    s.method = 3
    s.S, s.l_0 = 0, 0
    s.lamb, s.span = lamb, 2 * b
    s.alpha, s.U_0 = alpha, U0
    s.nFourier = 5
    s.ny = Ny
    s.run_analysis()

    if stop == 4: return

    set_question("5 - tapered wing")

    lamb, U0, alpha = 7, 101, np.deg2rad(4.1)
    WingAngle = [np.deg2rad(5), np.deg2rad(5)]
    CL, CD, CLa, CDa = wing_aerodynamics(lamb_ref, b_ref, alpha_ref, U0_ref, Angle=WingAngle, plantform="tapered")

    if stop == 5: return

    set_question("6 - delta wing")

    lamb, U0, alpha = 7, 101, np.deg2rad(4.1)
    WingAngle = [np.arctan(4 / lamb), 0]
    CL, CD, CLa, CDa = wing_aerodynamics(lamb_ref, b_ref, alpha_ref, U0_ref, Angle=WingAngle, plantform="delta")

    if plot_circulation:
        plt.legend();
        plt.grid();
        plt.show()

    if stop == 6: return

    set_question("7 - Circulation")

    if stop == 7: return

    set_question("8 - incidence effect")

    plot_circulation = False
    b = b_ref
    U0 = U0_ref
    lamb = lamb_ref
    leg = ["rectangular", "tapered", "delta", "elliptical"]
    Na = 28 + 1
    alpha_table = np.deg2rad(np.linspace(-4, 10, Na))
    CL, CD = np.zeros((4, Na)), np.zeros((4, Na))
    CLa, CDa = np.zeros((4, Na)), np.zeros((4, Na))
    for i in range(Na):
        CL[0, i], CD[0, i], CLa[0, i], CDa[0, i] = wing_aerodynamics(lamb, b, alpha_table[i], U0, Angle=[0, 0],
                                                                     plantform="rectangular")
        CL[1, i], CD[1, i], CLa[1, i], CDa[1, i] = wing_aerodynamics(lamb, b, alpha_table[i], U0,
                                                                     Angle=[np.deg2rad(5), np.deg2rad(5)],
                                                                     plantform="tapered")
        CL[2, i], CD[2, i], CLa[2, i], CDa[2, i] = wing_aerodynamics(lamb, b, alpha_table[i], U0,
                                                                     Angle=[np.arctan(4 / lamb), 0], plantform="delta")
        CL[3, i], CD[3, i], CLa[3, i], CDa[3, i] = wing_aerodynamics(lamb, b, alpha_table[i], U0, Angle=[0, 0],
                                                                     plantform="elliptical")
        # CL[3,i],CD[3,i]=np.pi*lamb*2/(lamb+2)*alpha_table[i],np.pi*lamb*(2/(lamb+2)*alpha_table[i])**2
        # CLa[3,i],CDa[3,i]=CL[3,i],CD[3,i]

    if show_results:
        print("alpha = ", np.rad2deg(alpha_table))
        for i in range(4):
            print(leg[i], "\nCL= ", CL[i, :], "\nCD= ", CD[i, :])

    if plot_incidence:
        fig = plt.figure(11)
        plt.title(r'$C_z$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(np.rad2deg(alpha_table), CL[i, :], label=leg[i])
            plt.plot(np.rad2deg(alpha_table), CLa[i, :], label=leg[i] + '_a')
        plt.ylabel(r"$C_z$")
        plt.xlabel(r"$\alpha (^\circ)$")
        plt.legend(loc='best')
        plt.grid()

        fig = plt.figure(12)
        plt.title(r'$C_x$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(np.rad2deg(alpha_table), CD[i, :], label=leg[i])
            plt.plot(np.rad2deg(alpha_table), CDa[i, :], label=leg[i] + '_a')
        plt.ylabel(r"$C_x$")
        plt.xlabel(r"$\alpha (^\circ)$")
        plt.legend(loc='best')
        plt.grid()

        fig = plt.figure(13)
        plt.title(r'$C_x=f(C_z)$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(CL[i, :], CD[i, :], label=leg[i])
            plt.plot(CLa[i, :], CDa[i, :], label=leg[i] + '_a')
        plt.ylabel(r"$C_x$")
        plt.xlabel(r"$C_z$")
        plt.legend(loc='best')
        plt.grid()

        fig = plt.figure(14)
        plt.title(r'$C_z/C_x$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(np.rad2deg(alpha_table), CL[i, :] / CD[i, :], label=leg[i])
            plt.plot(np.rad2deg(alpha_table), CLa[i, :] / CDa[i, :], label=leg[i] + '_a')
        plt.ylabel(r"$f$")
        plt.xlabel(r"$\alpha (^\circ)$")
        plt.legend(loc='best')
        plt.grid()

        plt.show()

    if stop == 8: return

    set_question("9 - wing aspect ratio effect")

    alpha = alpha_ref
    b = b_ref
    U0 = U0_ref

    if plot_wing_aspect_ratio:
        Na = 21
        lamb_table = np.linspace(4, 16, Na)
        CL, CD = np.zeros((4, Na)), np.zeros((4, Na))
        CLa, CDa = np.zeros((4, Na)), np.zeros((4, Na))
        for i in range(Na):
            CL[0, i], CD[0, i], CLa[0, i], CDa[0, i] = wing_aerodynamics(lamb_table[i], b, alpha, U0, Angle=[0, 0],
                                                                         plantform="rectangular")
            CL[1, i], CD[1, i], CLa[1, i], CDa[1, i] = wing_aerodynamics(lamb_table[i], b, alpha, U0,
                                                                         Angle=[np.deg2rad(1), np.deg2rad(1)],
                                                                         plantform="tapered")
            CL[2, i], CD[2, i], CLa[2, i], CDa[2, i] = wing_aerodynamics(lamb_table[i], b, alpha, U0,
                                                                         Angle=[np.arctan(4 / lamb_table[i]), 0],
                                                                         plantform="delta")
            CL[3, i], CD[3, i], CLa[3, i], CDa[3, i] = wing_aerodynamics(lamb_table[i], b, alpha, U0, Angle=[0, 0],
                                                                         plantform="elliptical")
            # CL[3,i],CD[3,i]=np.pi*lamb_table[i]*2/(lamb_table[i]+2)*alpha,np.pi*lamb_table[i]*(2/(lamb_table[i]+2)*alpha)**2
            # CLa[3,i],CDa[3,i]=CL[3,i],CD[3,i]

        if gradient_Cz:
            grad_Cz_elliptical = (2 * lamb_table * np.pi) / (2 + lamb_table) * np.pi / 180
        else:
            grad_Cz_elliptical = 1

        leg = ["rectangular", "tapered", "delta", "elliptic"]

        fig = plt.figure(15)
        plt.title(r'$C_z$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(lamb_table, CL[i, :] / grad_Cz_elliptical, label=leg[i])
            plt.plot(lamb_table, CLa[i, :] / grad_Cz_elliptical, label=leg[i] + '_a')
        plt.ylabel(r"$C_z$")
        plt.xlabel(r"$\lambda$")
        plt.legend(loc='best');
        plt.grid()

        fig = plt.figure(16)
        plt.title(r'$C_x$', fontsize=fz, fontweight='bold')
        for i in range(4):
            plt.plot(lamb_table, CD[i, :], label=leg[i])
            plt.plot(lamb_table, CDa[i, :], label=leg[i] + '_a')
        plt.ylabel(r"$C_x$")
        plt.xlabel(r"$\lambda$")
        plt.legend(loc='best');
        plt.grid()

        plt.show()

    print('end of  exercise ...')

    if gradient_Cz:
        alpha_ref = np.deg2rad(1)
    else:
        alpha_ref = np.deg2rad(4.1)

    lamb_ref = 7
    b_ref = 7
    U0_ref = 101
    Ns = 5
    k_ref = np.pi
    plot_incidence = False
    show_results = True             # to display CL, CD
    plot_wing_aspect_ratio = True


def Exercice5_11():
    """
    Numerical solution of the Prandt equation with ground effect
    """

    set_title("Numerical solution of the Prandt equation with ground effect, elliptical wing")

    gradient_Cz = True

    # to get the  gradient of Cz / alpha  instead of  this of  Cz (Cz is linear with alpha)
    if gradient_Cz:
        alpha = np.deg2rad(1)  # incidence of reference (1 degree)
    else:
        alpha = np.deg2rad(8)

    lamb = 7  # wing aspect ratio
    b = 7  # half-span
    # U0          = 50                # reference velocity
    Ns = 3  # number of Fourier mode
    k = np.pi  # half lift  gradient
    variation_h = True  # altitude effect
    variation_lambda = True  # to test wing aspect ratio effect
    M = 51  # number of points along the span

    # 
    theta = np.linspace(0, np.pi, M)  # angle distribution 
    y = -b * np.cos(theta)  # span coordinate vector
    l0 = 8 * b / (np.pi * lamb)  # mean chord
    S = np.pi * l0 * b / 2  # wing surface

    print('Ns = %i, M=%i, alpha= %2.0f °, lambda= %f' % (Ns, M, np.rad2deg(alpha), lamb))
    print('lambda = %f, 2b = %f, l0= %f, S = %f m^2' % (lamb, 2 * b, l0, S))

    set_question("I - WITHOUT GROUND EFFECT")

    A = Prandtl_E_OGE(Ns, lamb, b, l0, k, y, theta, alpha)
    CL = np.pi * lamb * A[0]
    Cxi = np.pi * lamb * np.sum(np.arange(1, Ns + 1) * A ** 2)
    A1_ref = alpha / (1 + np.pi * lamb / (2 * k))
    CL_ref = np.pi * lamb * A1_ref
    Cxi_ref = np.pi * lamb * A1_ref ** 2
    print('Calcul   : A_1 = %f, CL = %f, Cxi = %f ' % (A[0], CL, Cxi))
    print('Reference: A_1 = %f, CL = %f, Cxi = %f ' % (A1_ref, CL_ref, Cxi_ref))
    print("A1_bis = ", 2 * np.pi * lamb * alpha / (2 + lamb))

    set_question("II - WITH GROUND EFFECT, Lambda is given")

    if variation_h:
        Nh = 101
        h_over_b = np.linspace(0.0001, 4.0001, Nh)
        delta_CL = np.zeros((Nh), dtype=float)
        delta_CDi = np.zeros((Nh), dtype=float)

        for i, h in enumerate(h_over_b):
            Anew = Prandtl_E_IGE(Ns, lamb, b, l0, k, h, y, theta, alpha)
            # print('A = ',Anew)
            CL = np.pi * lamb * Anew[0]  # first coefficient of the Fourier serie
            Cxi = np.pi * lamb * np.sum(np.arange(1, Ns + 1) * Anew ** 2)
            delta_CL[i] = (CL - CL_ref) / CL_ref
            delta_CDi[i] = (Cxi - Cxi_ref) / Cxi_ref
            if i <= 25:
                print('h=%f,  A_1 = %2.6f, CL/CL_ref-1 = %f, Cxi/Cxi_ref-1 = %f' %
                      (h, Anew[0], CL / CL_ref - 1, Cxi / Cxi_ref - 1))

        # plot
        fig = plt.figure(1, figsize=(10, 8))
        plt.title(r"Ground effect, $\lambda$ = %f" % lamb, fontsize=fz, fontweight='bold')
        plt.xlabel(r'$h/b$', fontsize=16)
        plt.plot(h_over_b, delta_CL, color='b', linestyle='-', linewidth=2, label=r'$\Delta C_L/C_L $')
        plt.plot(h_over_b, delta_CDi, color='r', linestyle='-', linewidth=2, label=r'$\Delta C_{x_i}/C_{x_i} $')
        xmax, ymax = 4, 0.7
        plt.axis([0, xmax, 0, ymax])
        plt.grid(which='both');
        plt.legend(loc="best")

    set_question("II - WITH GROUND EFFECT, Lambda is variable, h is given")

    if variation_lambda:
        Nl = 101
        lamb_table = np.linspace(4, 24, Nl)
        delta_CL = np.zeros((Nl), dtype=float)
        delta_CDi = np.zeros((Nl), dtype=float)
        h = 0.25  #  h/b

        for i, lamb in enumerate(lamb_table):
            A1_ref = alpha / (1 + np.pi * lamb / (2 * k))
            CL_ref = np.pi * lamb * A1_ref
            Cxi_ref = np.pi * lamb * A1_ref ** 2
            l0 = 8 * b / (np.pi * lamb)
            Anew = Prandtl_E_IGE(Ns, lamb, b, l0, k, h, y, theta, alpha)
            CL = np.pi * lamb * Anew[0]  # first coefficient of the Fourier serie
            Cxi = np.pi * lamb * np.sum(np.arange(1, Ns + 1) * Anew ** 2)
            delta_CL[i] = (CL - CL_ref) / CL_ref
            delta_CDi[i] = (Cxi - Cxi_ref) / Cxi_ref
            if i <= 25:
                print('lambda=%f,  A_1 = %2.6f, CL/CL_ref-1 = %f, Cxi/Cxi_ref-1 = %f' %
                      (lamb, Anew[0], CL / CL_ref - 1, Cxi / Cxi_ref - 1))

        # plot
        fig = plt.figure(2, figsize=(10, 8))
        plt.title(r"Effet de sol, $h/b$ = %f" % (h), fontsize=fz, fontweight='bold')
        plt.xlabel(r'$\lambda$', fontsize=16)
        plt.plot(lamb_table, delta_CL, color='b', linestyle='-', linewidth=2, label=r'$\Delta C_L/C_L (\%)$')
        plt.plot(lamb_table, delta_CDi, color='r', linestyle='-', linewidth=2, label=r'$\Delta C_{x_i}/C_{x_i} (\%)$')
        plt.axis([4, 24, 0, 0.4])
        plt.xticks(np.arange(4, 25, 4))
        plt.grid();
        plt.legend(loc="best")

    if variation_h or variation_lambda: plt.show()


def Exercice5_12():
    """
    Spitfire wing
    """
    set_title("Spitfire wing")
    npt = 101  # number of point along span  y
    b = 18 * 0.3028  # half span (foot to meters )
    Lambda = 5.5  # wing aspect ration

    L0 = 8 * b / (np.pi * Lambda)  # root wing chord
    S = np.pi * L0 * b / 2  # wing surface
    x0 = 0  # line giving the aerodynamic center of each airfoil along the span
    fleche = 30  # en degrees
    k = np.tan(np.deg2rad(fleche))  # if we want to set a swept wing angle

    print("Span                     : %3.3f m" % (2 * b))
    print("Wing aspect ratio        : %3.3f " % (Lambda))
    print("Wing surface             : %3.3f m^2" % (S))
    print("root wing chord          : %3.3f m" % (L0))
    theta = np.linspace(0, np.pi / 2, npt)
    y = b * np.cos(theta)
    x0 = -k * y / b
    chord = L0 * np.sqrt(1 - (y / b) ** 2)
    xLE, xTE = x0 + chord / 4, x0 - 3 * chord / 4
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$y (m)$', fontsize=16)
    ax.set_ylabel(r'$x (m)$', fontsize=20)
    ax.set_title("Spitfire half wing")
    ax.plot(y, xLE, color='b', linestyle='-', linewidth=2)
    ax.plot(y, xTE, color='r', linestyle='-', linewidth=2)
    ax.axis("equal")
    plt.show()


def Exercice5_13():
    """ 
    CL, CD : elliptical wing without ground effect
    """
    set_title("elliptical wing without ground effect")

    print('elliptical wing only')

    b = np.float64(1.0)  # half-span
    lamb = np.float64(10.0)  # wing aspect ratio
    alpha = np.deg2rad(np.float64(1))  # incidence in degrees
    N = 2  # number of Fourier mode
    M = 21  # number of points along the span
    h = np.float(64)  # height / half-span ?
    #  Uinf=np.float(1.)           # upstream velocity
    k = np.pi  # half lift gradient of the airfoil
    theta = np.linspace(0, np.pi, M)  # angle distribution
    y = -b * np.cos(theta)  # span coordinate vector
    l0 = 8 * b / (np.pi * lamb)  # mean chord
    S = np.pi * l0 * b / 2  # wing surface

    print('N = %i, M=%i, alpha= %2.0f °, lambda= %f' % (N, M, np.rad2deg(alpha), lamb))
    print('lambda = %f, 2b = %f, l0= %f, S = %f m^2' % (lamb, 2 * b, l0, S))

    A = Prandtl_E_OGE(N, lamb, b, l0, k, y, theta, alpha)

    CL = np.pi * lamb * A[0]
    Cxi = np.pi * lamb * np.sum(np.arange(1, N + 1) * A ** 2)

    A1_ref = alpha / (1 + np.pi * lamb / (2 * k))
    CL_ref = np.pi * lamb * A1_ref
    Cxi_ref = np.pi * lamb * A1_ref ** 2

    print('Calculus : A_1 = %f, CL = %f, Cxi = %f ' % (A[0], CL, Cxi))
    print('Reference: A_1 = %f, CL = %f, Cxi = %f ' % (A1_ref, CL_ref, Cxi_ref))

    set_title("elliptical wing with ground effect")

    Nh = 31
    # h_over_b=np.linspace(0.01,5.01,Nh)
    h_over_b = np.linspace(0.01, 5.01, Nh)
    err_CL = np.zeros((Nh), dtype=float)
    err_CDi = np.zeros((Nh), dtype=float)
    i = 0

    for h in list(h_over_b):
        Anew = Prandtl_E_IGE(N, lamb, b, l0, k, h, y, theta, alpha)
        CL = np.pi * lamb * Anew[0]
        Cxi = np.pi * lamb * np.sum(np.arange(1, N + 1) * Anew ** 2)
        err_CL[i] = CL / CL_ref - 1
        err_CDi[i] = Cxi / Cxi_ref - 1
        print('h=%f,  A_1 = %2.6f, CL/CL_ref-1 = %f, Cxi/Cxi_ref-1 = %f' %
              (h, Anew[0], CL / CL_ref - 1, Cxi / Cxi_ref - 1))
        i = i + 1

    # plot
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$2h/b=h/L$', fontsize=16)
    # ax.set_ylabel(r'$Kp$',fontsize=20)
    ax.set_title('Ratio / to the infinite height')
    ax.plot(2 * h_over_b, err_CL, color='b', linestyle='-', linewidth=2, label=r'$e(C_L)$')
    ax.plot(2 * h_over_b, err_CDi, color='r', linestyle='-', linewidth=2, label=r'$e(C_{x_i})$')
    xmax, ymax = 5, 0.5
    ax.axis([0, xmax, 0, ymax])
    xmajor_ticks, xminor_ticks = np.arange(0, xmax, 1), np.arange(0, xmax, 0.25)
    ymajor_ticks, yminor_ticks = np.arange(0, ymax, 0.2), np.arange(0, ymax, 0.05)
    ax.set_xticks(xmajor_ticks);
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks);
    ax.set_yticks(yminor_ticks, minor=True)
    ax.grid(which='both');
    ax.legend()
    plt.show()


def func(x, a, b, c):
    """
    test function
    """
    return a + b * x + c * x ** 2


def integralO2(a, b, c):
    """
    trapeze integral rule
    
    why do not use numpy function ?
    """
    m = 100
    x = np.linspace(0, 1, m)
    dx = x[1] - x[0]
    val = func(x, a, b, c)
    s = dx * (np.sum(val[1:m - 2]) + 0.5 * (val[0] + val[m - 1]))
    return s


def Exercice5_15():
    """ 
    test  of the quad integral with parameters
    comparison with trapeze rule
    """
    A, B = 1, 2
    C = np.linspace(0, 2, 2)
    for c in list(C):
        G, res = integrate.quad(func, 0, 1., args=(A, B, c,))
        G1 = integralO2(A, B, c)
        print(G, G1, A + B / 2 + c / 3)


def Exercice5_14():
    """
    Tapered wing load of the  P51 Mustang
    """
    set_title("Tapered wing load of the  P51 Mustang")
    f2m = 0.3048
    l_0, l_t = 8.48 * f2m, 3.87 * f2m
    b = 5.64
    Lambda = np.rad2deg(np.arctan((l_0 - l_t) / b))
    lamb = 5.876
    U0 = 123.

    Ny = 101
    alpha = np.deg2rad(3.)
    Ns = 5

    par = set_wing_parameters()
    s = WingAnalysis(par)
    s.plot_wing = True
    s.plot_circulation = True
    s.method = 1  # sine approach
    s.S = 0
    s.span = 2 * b
    s.lamb = lamb
    s.alpha = alpha
    s.wing = "tapered"
    s.U_0 = U0
    s.nFourier = Ns
    s.ny = Ny
    s.WingAngle = [np.deg2rad(Lambda / 2), np.deg2rad(Lambda / 2)]
    s.run_analysis()
