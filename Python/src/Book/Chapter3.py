#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2017 - 2020*

*license : AGPL-3.0*

Correction of the chapter 3 exercises
    ..
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import os
from Tools.misc import set_title, set_question, set_info, set_alert
from scipy.optimize import newton
from Book.ConformalMapping.mappings import set_parameters, AirfoilConfMap
from Book.panneaux_intro import plaque_plane, dessiner_contours, reference_plaque_plane


FigDim = [10, 10]  # figure size
fz = 16  # fontsize


def Exercice3_1():
    """
    Cylinder flow
    """
    set_title("Pressure coefficient on a  cylinder or a sphere")
    theta = np.linspace(0, 360, 361)
    set_question('1- cylinder, 2- sphere')
    plt.figure(figsize=(FigDim))
    plt.title(r'Kp on a  cylinder or a sphere', fontsize=fz, fontweight='bold')
    plt.plot(theta, 1 - 4 * (np.sin(np.deg2rad(theta))) ** 2, linewidth=2, label='cylindre')
    plt.plot(theta, 1 - 9 / 4 * (np.sin(np.deg2rad(theta))) ** 2, linewidth=2, label='sphère')
    plt.xlabel(r"$\theta (^\circ)$", fontsize=fz);
    plt.ylabel(r'$Kp$', fontsize=fz)
    plt.legend(loc="best")
    plt.grid()
    plt.show()


def plot_Kp(xe, xi, Kpe, Kpi, cas):
    """
    plot of a Kp list
    
    Pars:
        xe (real) : x coordinate on the upper wall
        xi (real) : x coordinate on the lower wall
        Kpe (real) : Kp on the upper wall
        Kpi (real) : Kp on the lower wall
    """
    plt.figure(figsize=(15, 10))
    plt.title(r"$K_p$")
    for k in range(len(xe)):
        plt.plot(xe[k], Kpe[k], label="Kpe: " + cas[k])
        plt.plot(xi[k], Kpi[k], label="Kpi: " + cas[k])
    plt.grid()
    plt.xlim(-0.05, 1.05)
    plt.legend(loc="best")



def Exercice3_2():
    """
    Joukowski's transformation
    """
    set_title("Joukowski's transformation")
    xe, xi = [], []
    Kpe, Kpi = [], []
    cas = ["Ellipse", "Thick airfoil", "Cambered line", "Cambered airfoil"]  # title for plot
    set_question("1 - Ellipse")
    par = set_parameters()
    s = AirfoilConfMap(par)
    # Joukowski's airfoi, it is always plotted under zero angle of attack
    #  the parameter is Zc, the circle center
    s.map = "Joukowski1"
    s.airfoil = "Ellipse"
    s.alpha = 0
    s.a = 0.9
    s.plot_circle = False
    s.plot_airfoil = False
    s.opt_levels = "manuel"
    s.levels_iso = [-2, 2, 0.2]
    s.adimChord = True
    s.plot_velocity = False
    s.plot_Kp = True

    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("2 - Thick airfoil")

    par = set_parameters()
    s = AirfoilConfMap(par)
    # Joukowski's airfoil, it is always plotted under zero angle of attack
    #  the parameter is Zc, the circle center
    s.map = "Joukowski1"
    s.a = 1
    s.R = 0
    s.alpha = 0
    s.Zc = -0.1 + 1j * 0
    s.plot_circle = True
    s.plot_airfoil = False
    s.opt_levels = "manuel"
    s.levels_iso = [-1, 1, 0.2]
    s.adimChord = True
    s.plot_velocity = False
    s.plot_Kp = True
    s.run_airfoil()

    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("3 - Cambered line")

    s.a = 1
    s.R = 0
    s.Zc = 0 + 1j * 0.1
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("4 - Cambered airfoil")

    s.a = 1
    s.R = 0
    s.Zc = -0.1 + 1j * 0.1
    s.camberline = True
    s.plot_airfoil = True

    s.run_airfoil()

    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    plot_Kp(xe, xi, Kpe, Kpi, cas)
    plt.show()

def Exercice3_3():
    """
    Karman-Trefftz's transformation
    """
    set_title("Karman-Trefftz's transformation")
    xe, xi = [], []
    Kpe, Kpi = [], []
    cas = ["sym5", "cambered-5", "sym-10", "cambered-10"]
    set_question("3 - symmetrical airfoil")
    par = set_parameters()
    s = AirfoilConfMap(par)
    eps = [0.05, 0.1]
    for epsi in eps:
        s.map = "Karman-Trefftz"
        s.a = 1
        s.R = 0
        s.alpha = 0
        s.beta = 0
        s.Zc = -epsi + 1j * 0
        s.camberline = False
        s.plot_airfoil = False
        s.opt_levels = "manuel"
        s.levels_iso = [-2., 2., 0.1]
        s.k = 1.9
        s.plot_velocity = False
        s.adimChord = True
        s.plot_Kp = True
        s.run_airfoil()

        xe.append(s.xe)
        xi.append(s.xi)
        Kpe.append(s.Kpe)
        Kpi.append(s.Kpi)

        set_question("3 - cambered airfoil")
        s.Zc = -epsi + 1j * epsi
        s.camberline = False
        s.plot_airfoil = True
        s.a = 1
        s.R = 0
        s.run_airfoil()
        xe.append(s.xe)
        xi.append(s.xi)
        Kpe.append(s.Kpe)
        Kpi.append(s.Kpi)

    plot_Kp(xe, xi, Kpe, Kpi, cas)
    plt.show()

def Exercice3_4():
    """
    Von Mises's transformation
    """
    set_title("Von Mises's transformation")
    xe, xi = [], []
    Kpe, Kpi = [], []
    cas = ["Zc=0;5", "cambered-5", "cambered-10"]
    set_question("9 - Airfoil eps=0.05")
    par = set_parameters()
    s = AirfoilConfMap(par)

    # von Mises symmetrical airfoil
    # the parameters are eps (epsilon) a,d k  
    # mu is determined from Zc
    # The airfoil is the one given in the book
    # it works when a=2 and mu_0+mu_1+a=0
    s.a = 1
    s.R = 0
    s.map = "von Mises"
    s.alpha = 0.
    s.k = 2
    s.eps = 0.05
    s.Zc = 0
    # s.Zc            = -s.eps*(1-1j)
    # s.camberline    = False
    # s.plot_airfoil  = False
    s.opt_levels = "manuel"
    s.levels_iso = [-2., 2., 0.1]
    s.plot_velocity = False
    s.adimChord = True
    s.plot_Kp = True
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - airfoil eps=0.05")
    s.a = 1
    s.R = 0
    s.eps = 0.05
    s.Zc = -s.eps * (1 - 1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - airfoil eps=0.10")
    s.a = 1
    s.R = 0
    s.eps = 0.10
    s.Zc = -s.eps * (1 - 1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    plot_Kp(xe, xi, Kpe, Kpi, cas)
    plt.show()

def Exercice3_5():
    """
    Van de Vooren et Jong's transformation
    """
    set_title("Van de Vooren et Jong's transformation")
    xe, xi = [], []
    Kpe, Kpi = [], []
    cas = ["sym, 0.05", "sym, 0.10", "cambered-5", "cambered-10"]

    set_question("9 - airfoil eps=0.05")
    par = set_parameters()
    s = AirfoilConfMap(par)

    s.map = "van de Vooren"
    s.a = 1
    s.R = 0
    s.alpha = 0.
    s.k = 2
    s.eps = 0.05
    s.Zc = 0 + 1j * 0  # -s.eps*(1-1j)
    # s.camberline    = False
    # s.plot_airfoil  = False
    s.opt_levels = "manuel"
    s.levels_iso = [-2., 2., 0.1]
    s.plot_velocity = False
    s.adimChord = True
    s.plot_Kp = True
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - airfoil eps=0.10")
    s.a = 1
    s.R = 0
    s.alpha = 0
    s.eps = 0.10
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - airfoil eps=0.05")
    s.a = 1
    s.R = 0
    s.alpha = 0
    s.eps = 0.05
    s.k = 1.9
    s.Zc = -s.eps * (1 - 1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)

    set_question("9 - airfoil eps=0.10")
    s.a = 1
    s.R = 0
    s.alpha = 0
    s.eps = 0.10
    s.k = 1.9
    s.Zc = -s.eps * (1 - 1j)
    s.run_airfoil()
    xe.append(s.xe)
    xi.append(s.xi)
    Kpe.append(s.Kpe)
    Kpi.append(s.Kpi)
    plot_Kp(xe, xi, Kpe, Kpi, cas)
    plt.show()

def Exercice3_6():
    """
    Double-cambered airfoil
    
    Parameters :
        a and  b (complex),  180° > arg(b) > 90° 

    Problem:
        Depending on the parameter choice, the streamlines can enter into the airfoil, it is related to arctan
    
    """
    set_title("Double-cambered airfoil")

    set_alert(" Depending on the parameter choice, the streamlines can enter into the airfoil, it is related to arctan")
    set_info("Be careful")
    cas = 'd'

    set_question("9 - airfoil case %s" % (cas))
    par = set_parameters()
    s = AirfoilConfMap(par)
    s.map = "double sharp"  # double sharp leading and trailing edge
    s.alpha = 0.
    s.a = 0  # It is determined from the circle radius

    if cas == 'c':
        """
        Case with problems
        """
        beta = 46
        s.b = s.R / 2 * np.exp(1j * np.deg2rad(beta))
        s.Zc = 0 + 1j * 0

        #  case 'P2'
        #      b = a/2;
        #      Z_c = 0;

        #  case 'P3'
        #      beta = 20*pi/180;
        #      epsilon = 0.2;
        #      b = sqrt(epsilon) * a * exp(i*beta)
        #      Z_c = 0;


    elif cas == "d":
        """
        Case with problems
        """
        s.eps = 0.2
        beta = 20
        s.b = np.sqrt(s.eps) * s.R * np.exp(1j * np.deg2rad(beta))
        s.Zc = 0 + 1j * 0

    elif cas == "b":
        beta = 135
        s.b = s.R / 10 * np.exp(1j * np.deg2rad(beta))
        s.Zc = s.b  # -s.eps*(1-1j)

    elif cas == "a":
        beta = 180
        s.b = s.R / 10 * np.exp(1j * np.deg2rad(beta))
        s.Zc = s.b  # -s.eps*(1-1j)

    s.run_airfoil()
    plt.show()

def Exercice3_7():
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
    s.plot_velocity = True
    s.plot_Kp = True
    s.run_airfoil()
    plt.show()

def Exercice3_8():
    """
    Flate plate : introduction to panels
    """

    # Parameter :
    # U0              = 1       #  upstream velocity, unused
    plot_panels = False
    # plot_gamma      = False
    plot_vitesse = False        # to plot velocity
    plot_psi = False            # to plot the streamlines
    plot_2D = True              # to plot the iso-velocity
    alpha = np.deg2rad(5)       # incidence in degrees , must be different of 0 !
    nP = 101                    # panel number of the initial plate
    # chord           = 1.      # plate length (chord)
    nv = 5001                   # number of points to define the grid and to calculate the local velocity (plot 1D)
    y_wall = 0.01               # first point from the wall height

    s = plaque_plane(nP, y_wall, nv=nv, alpha=alpha)
    dessiner_contours(s[0], s[2], alpha)
    nP1, nP2 = 21, 3
    s1 = plaque_plane(nP1, y_wall, nv=nv, alpha=alpha)
    s2 = plaque_plane(nP2, y_wall, nv=nv, alpha=alpha)

    if plot_2D:
        option = "uv"           # possibly  "uv" ou "Kp"
        vel = ["v", "u"]
        plt.figure(figsize=(12, 12))
        if option == "Kp":
            plt.plot(s[3].real, 1 - abs(s[4]), '-', label="U, nP = %i " % (nP))
            plt.plot(s1[3].real, 1 - abs(s1[4]), '-', label="U, nP = %i " % (nP1))
            plt.plot(s2[3].real, 1 - abs(s2[4]), '-', label="U, nP = %i " % (nP2))
        else:
            if "u" in vel:
                plt.plot(s[3].real, s[4].real, '-', label="u, nP = %i " % (nP))
                plt.plot(s1[3].real, s1[4].real, '-', label="u, nP = %i " % (nP1))
                plt.plot(s2[3].real, s2[4].real, '-', label="u, nP = %i " % (nP2))
                eps = 0.01
                x = np.linspace(eps, 1 - eps, 101)
                plt.plot(x, reference_plaque_plane(x, alpha), 'r--', label="Analytical")
            if "v" in vel:
                plt.plot(s[3].real, -s[4].imag, '-', label="v, nP = %i " % (nP))
                plt.plot(s1[3].real, -s1[4].imag, '-', label="v, nP = %i " % (nP1))
                plt.plot(s2[3].real, -s2[4].imag, '-', label="v, nP = %i " % (nP2))
        plt.ylim(-1, 2)
        plt.grid()
        plt.legend(loc="best")

        plt.figure(figsize=(10, 8))
        plt.plot(s[0].real, s[2].real * len(s[2]) / s[5], '--s', label=r"$\Gamma$, nP = %i " % (nP))
        plt.plot(s1[0].real, s1[2].real * len(s1[2]) / s1[5], '-o', Markersize=10, label=r"$\Gamma$, nP = %i " % (nP1))
        plt.plot(s2[0].real, s2[2].real * len(s2[2]) / s2[5], '-o', Markersize=10, label=r"$\Gamma$, nP = %i " % (nP2))
        plt.ylim(0, 6)
        plt.grid()
        plt.legend(loc="best")

    s3 = plaque_plane(nP, -y_wall, nv=nv, alpha=alpha)
    plt.figure(figsize=(12, 12))
    plt.title(r"$\Delta u = u^+-u^-$ at $|y|$ = %f" % (y_wall))
    plt.plot(s[3].real, s[4].real - s3[4].real, 'b-', label=r"$\Delta U$, nP = %i " % (nP))

    plt.grid()
    plt.legend(loc="best")

    if plot_2D or plot_vitesse or plot_psi or plot_panels:
        plt.show()


def circulation_ratio(eta, alpha, exacte=False):
    """
    Circulation ratio (ground effect)
    exacte = False    : exact or approximated method
    """
    if exacte:
        return (1 + eta ** 2 / 16 * np.cos(alpha) ** 2) / (1 + eta / 4 * np.sin(alpha))
    else:
        return 1 + eta ** 2 / 16


def lift_ratio(eta, alpha, exacte=False):
    """
    Lift ratio  (ground effect)
    exacte = False    : exact or approximated method
    """
    K = circulation_ratio(eta, alpha, exacte)
    if exacte:
        return K * (1 - K * eta * np.sin(alpha) / (4 + eta * np.sin(alpha)))
    else:
        return K * (1 - eta / 4 * np.sin(alpha) * K)


def Exercice3_9():
    """ 
    Ground effect
    """
    set_title("Ground effect")

    # curve of the maximal values
    al = np.linspace(2, 15, 27)
    eta_ref = np.linspace(0, 20, 1001)
    r_appro_max, eta_appro_max = [], []
    r_exact_max, eta_exact_max, = [], []
    for k in range(len(al)):
        r0 = lift_ratio(eta_ref, np.deg2rad(al[k]), False)
        k0 = np.argmax(r0)
        r_appro_max.append(np.amax(r0))
        eta_appro_max.append(eta_ref[k0])
        r0 = lift_ratio(eta_ref, np.deg2rad(al[k]), True)
        k0 = np.argmax(r0)
        r_exact_max.append(np.amax(r0))
        eta_exact_max.append(eta_ref[k0])

    option = "eta"  # ou "1/eta"
    eta_lim = [17, 17, 9, 11, 9, 12, 12, 14]
    alpha = [3, 3, 5, 5, 10, 10, 15, 15]
    eta = []
    r_appro = []
    r_exact = []
    for k in range(len(eta_lim)):
        eta0 = np.linspace(0.1, eta_lim[k], 101)
        alpha0 = np.deg2rad(alpha[k])
        if k % 2 == 0:
            exact = False
        else:
            exact = True
        r_appro.append(lift_ratio(eta0, alpha0, exact))
        r_exact.append(lift_ratio(eta0, alpha0, exact))
        eta.append(eta0)

    fig = plt.figure(figsize=(10, 8))
    fig.suptitle(r'Ground effect, $C_L/C_{L_\infty}$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)

    ax.grid()
    for k in range(len(alpha)):
        if k % 2 == 0:
            leg = r"appro., $\alpha=$ %2.2f"
        else:
            leg = r"exact,  $\alpha=$ %2.2f"
        if option == "eta":
            ax.plot(eta[k], r_appro[k], '-', linewidth=2, label=leg % (alpha[k]))
            ax.plot(eta_appro_max, r_appro_max, 'k--')
            ax.plot(eta_exact_max, r_exact_max, 'r--')
        else:
            ax.plot(1 / eta[k], r_appro[k], '-', linewidth=2, label=leg % (alpha[k]))
            ax.plot([1 / x for x in eta_appro_max], r_appro_max, 'k--')
            ax.plot([1 / x for x in eta_exact_max], r_exact_max, 'r--')
    if option == "eta":
        plt.xlim(0, 10)
        ax.set_xlabel(r'$\eta=\ell/h$', fontsize=20)
    else:
        plt.xlim(0, 1)
        ax.set_xlabel(r'$\eta=h/\ell$', fontsize=20)
    plt.ylim(0, 3)
    plt.legend(loc="best")
    plt.show()


def Exercice3_10():
    """
    Tandem wings
    """
    pass


def Exercice3_11():
    """
    Schwarz-Christoffel 's tranformation: channel expansion 
    """
    pass


def corps_Rankine(y, m):
    """
    Rankine's oval plot
    
    Problems:
        take care at y=0
          
    """
    x2 = []
    for Y in y:
        if abs(Y) > 1e-6:
            x2.append(2 * Y / np.tan(2 * Y / m) + 1 - Y ** 2)
        else:
            x2.append(m + 1 - Y ** 2)
    # test to avoid problems
    for i, X2 in enumerate(x2):
        if X2 < 0:
            # print("x^2 = %e, i = %i"%(X2,i))
            x2[i] = 0
    return np.sqrt(x2)


def plot_func(eta, m):
    """
    plot of the function to research zeros (Rankine's oval)
    """
    fig = plt.figure()
    fig.suptitle('search of  h/a for  m = %1.2f' % m, fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$eta$', fontsize=20)
    ax.grid()
    ax.plot(eta, np.tan(eta / m), '-', linewidth=2, color='black')
    ax.plot(eta, 1 / eta, '-', linewidth=2, color='red')
    plt.show()


def Kp_Rankine(x, y, m):
    """
    Kp for the Rankine's body
    """
    Kp = []
    U = []
    T = []
    for X, Y in zip(x, y):
        z = complex(X, Y)
        w = 1 - m / (z ** 2 - 1)
        T.append(1 + m - abs(w))
        U.append(abs(w))
        Kp.append(1 - abs(w) ** 2)
    return Kp, U, T


def solve_Rankine(m, eta_init=1, display=False):
    """
    Solve the Rankine's geometry for a given parameter m
    """
    # eta = h/a
    # ell = L/a
    ell = np.sqrt(1 + m)

    def func(eta):
        return np.tan(eta / m) - 1 / eta

    # use of the Newton method intrinsic in scipy
    eta = newton(func, eta_init)
    Uaxe = 1 + m / (1 + eta ** 2)
    # research of  max U, min Kp    
    y = np.linspace(0, eta, 1000);
    x = corps_Rankine(y, m)
    Kp1, U1, T1 = Kp_Rankine(x, y, m)
    Umax = np.amax(U1);
    Imax = np.argmax(U1)
    if display: print('Umax/U0 = %f, Kpmin = %f , x( %i ) = %f' % (Umax, Kp1[Imax], Imax, x[Imax]))
    return ell, np.float(eta), Uaxe, Umax, Kp1[Imax], x[Imax]


def Exercice3_12():
    """ 
    Rankine's oval
    """

    # a use of a class could be better

    set_title("Rankine's oval")

    # eta = h/a
    # ell = L/a
    #   
    # USER INPUT      
    #
    cas = 12            # to solve a given test case, case number: n, see m_table
    n = 201             # number of points for the plot
    eps = 0.4           # truncation of the domain if  option=2
    plot_oval = True   # to plot the body
    plot_Kp = True      # to plot the Kp
    plot_U = True       # to plot U/U0
    plot_test = False   # to plot the function 1+m-(u^2+v^2)
    overview = True     # if option=1, calculate the bodies for a large number of  m values
    #                     if False, we can have a  zoom near a given value of m.
    option = 4          # main option :
    #                       = 1 : table of values
    #                       = 2 : the value of eta (for m < 0.1) can be found graphically
    #                       = 3 : plot of the body and of Kp for a given  m  by the  variable "cas"
    #                       = 4  : plots found in the book
    # For the small  m  value the initial guess (in Newton method) must be very close to the converged value 

    k_compare = [1, 15]  # table of the case "cas" which are plotted on the same figure (option=4)

    if overview:
        # table for  m
        m_table = np.array(
            [0.02, 0.1, 0.15, 0.1589, 0.2, 0.3, 0.4, 0.5, 0.6, 0.639, 0.65, 0.70, 0.8, 0.88, 0.99, 1, 2, 10, 20, 50])
        # table of initial value for  eta before convergence, for the table of  m above         
        eta_table_init = np.array(
            [0.0308, 0.143, 0.17, 0.2, 0.263, 0.3, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.65, 1, 1, 1, 1, 1, 1, 1])
    else:
        m_table = np.linspace(0.87, 0.88, 21)
        eta_table_init = np.ones(len(m_table))

    Nplot = len(m_table)
    m_cas = m_table[cas]
    eta_init = eta_table_init[cas]

    ar = []  # aspect ratio
    eta_val = []  # eta value as a function of  m

    if option == 1:
        # table
        print(" m  \t\t L/a \t h/a \t L/h \t\t h/L \t U(x=0)  Kp_axe \t Umax \t Kpmin \t x_min ")
        for m, eta_init in zip(m_table, eta_table_init):
            ell, eta, Uaxe, Umax, Kpmin, x_max = solve_Rankine(m, eta_init)
            ar.append(eta / ell)
            eta_val.append(eta)
            print(
                '{:3.4f}   \t {:1.3f} \t {:1.3f} \t  {:2.3f} \t {:2.3f} \t {:1.3f} \t {:1.3f} \t {:1.3f} \t {:1.3f}  '
                '{:1.5f}'.format(m, ell, eta, ell / eta, eta / ell, Uaxe, 1 - Uaxe ** 2, Umax, Kpmin, x_max))

        # plot the aspect ratio versus m

        fig = plt.figure(10)
        fig.suptitle('relative thickness', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$m$', fontsize=20)
        ax.set_ylabel(r'$h/\ell$', fontsize=20)
        ax.grid()
        # ax.axis([1.0, 12, -0.05, 0.35])
        ax.plot(m_table[:Nplot], ar[:Nplot], '-', linewidth=2, color='black', label=r"$h/\ell$")
        ax.plot(m_table[:Nplot], eta_val[:Nplot], '-', linewidth=2, color='red', label=r"$\eta$")
        ax.legend(loc='lower right')
        plt.show()

    elif option == 2:
        print('case : m = %f ' % (m_cas))
        eta = np.linspace(eps, 1.5 * m_cas, n)
        # eta=m_cas*np.cos(np.linspace(eps,np.pi/2-eps,n))
        print('eta = ', eta)
        plot_func(eta, m_cas)

    elif option == 3:
        # particular case
        ell, eta, Uaxe, Umax, Kpmin, x_max = solve_Rankine(m_cas, eta_init)
        print("m = %f, L/a = %f, h/a = %f, L/h = %f, h/L = %f" % (m_cas, ell, eta, ell / eta, eta / ell))
        print("U(x=0)/U0 = %f, Kp(x=0) = %f" % (Uaxe, 1 - Uaxe ** 2))
        y = np.linspace(-eta, eta, n)
        x = corps_Rankine(y, m_cas)
        Kp1, U1, T1 = Kp_Rankine(x, y, m_cas)
        Kp2, U2, T2 = Kp_Rankine(-x, y, m_cas)

        Um = np.amax(U1)
        Im = np.argmax(U1)

        print('Umax/U0 = %f, Kpmin = %f , x( %i ) = %f' % (Um, Kp1[Im], Im, x[Im]))
        if plot_oval:
            fig = plt.figure(1)
            fig.suptitle("Rankine's oval for m = %1.2f" % m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$', fontsize=20)
            ax.set_ylabel(r'$y/a$', fontsize=20)
            ax.grid()
            # ax.axis([1.0, 12, -0.05, 0.35])
            ax.plot(x, y, '-', linewidth=2, color='black')
            ax.plot(-x, y, '-', linewidth=2, color='black')
            ax.axis('equal')
        if plot_Kp:
            fig = plt.figure(2)
            fig.suptitle('Kp: Rankine oval for m = %1.2f' % m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$', fontsize=20)
            ax.set_ylabel(r'$Kp$', fontsize=20)
            ax.grid()
            # ax.axis([1.0, 12, -0.05, 0.35])
            ax.plot(x, Kp1, '-', linewidth=2, color='black')
            ax.plot(-x, Kp2, '-', linewidth=2, color='red')
        if plot_U:
            fig = plt.figure(3)
            fig.suptitle('U: Rankine oval for m = %1.2f' % m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$', fontsize=20)
            ax.set_ylabel(r'$U/U_0$', fontsize=20)
            ax.grid()
            # ax.axis([1.0, 12, -0.05, 0.35])
            # Pour comparer à Katz et Plotkin p. 62
            # ax.plot(x,[U1**2 for U1 in U1],'-',linewidth=2,color='black')
            # ax.plot(-x,[U2**2 for U2 in U2],'-',linewidth=2,color='red')
            ax.plot(x, U1, '-', linewidth=2, color='black')
            ax.plot(-x, U2, '-', linewidth=2, color='red')
        if plot_test:
            fig = plt.figure(4)
            fig.suptitle(r'$1+m-z^2$: Rankine oval for m = %1.2f' % m_cas, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$x/a$', fontsize=20)
            ax.set_ylabel(r'$T$', fontsize=20)
            ax.grid()
            # ax.axis([1.0, 12, -0.05, 0.35])
            ax.plot(x, T1, '-', linewidth=2, color='black')
            ax.plot(-x, T2, '-', linewidth=2, color='red')

        if plot_oval or plot_Kp or plot_U or plot_test: plt.show()

    elif option == 4:
        fig = plt.figure(6)
        fig.suptitle("Rankine's ovals", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x/a$', fontsize=20)
        ax.set_ylabel(r'$y/a$', fontsize=20)
        ax.grid()

        fig1 = plt.figure(7)
        fig1.suptitle("U: Rankine's ovals", fontsize=14, fontweight='bold')
        ax1 = fig1.add_subplot(111)
        fig1.subplots_adjust(top=0.80)
        ax1.set_xlabel(r'$x/a$', fontsize=20)
        ax1.set_ylabel(r'$U/U_0$', fontsize=20)
        ax1.grid()

        p = []
        for k in k_compare:
            m = m_table[k]
            ell, eta, Uaxe, Umax, Kpmin, x_max = solve_Rankine(m, eta_table_init[k])
            print("m = %f, L/a = %f, h/a = %f, L/h = %f, h/L = %f" % (m, ell, eta, ell / eta, eta / ell))
            print("U(x=0)/U0 = %f, Kp(x=0) = %f" % (Uaxe, 1 - Uaxe ** 2))
            y = np.linspace(-eta, eta, n)
            x = corps_Rankine(y, m)
            Kp1, U1, T1 = Kp_Rankine(x, y, m)
            Kp2, U2, T2 = Kp_Rankine(-x, y, m)
            p.append([x, y, U1])

            ax.plot(x, y, '-', linewidth=2)
            ax.plot(-x, y, '-', linewidth=2)
            # ax.axis([1.0, 12, -0.05, 0.35])
            ax.axis('equal')

            ax1.plot(x, U1, '-', linewidth=2)
            ax1.plot(-x, U2, '-', linewidth=2)
        plt.show()
    return


def Exercice3_20():
    """ 
    Conformal transformation applications (old version)
    """
    set_title("Conformal transformations")

    par = set_parameters()
    s = AirfoilConfMap(par)

    # choose a case : 
    cas = 6

    if cas == 0:
        # Joukowski's airfoil, airfoil plotted under zero incidence 
        # parameter : Zc
        s.map = "Joukowski1"
        s.alpha = 0
        s.Zc = -0.0895386257 + 1j * 0.0523359562  # -0.15+1j*0.1 #
        s.opt_levels = "manuel"
        s.levels_iso = [-1, 1, 0.2]
        s.adimChord = True
        s.n_circle = 101
        s.plot_velocity = True
        s.plot_Kp = True
    elif cas == 1:
        # Joukowski's airfoil, the airfoil plot move as a function of the incidence
        # The transformation contains  alpha
        #  the parameters are  lamb and beta
        s.map = "Joukowski"
        s.alpha = 0
        s.beta = 3
        s.lamb = 1.1
        s.n_circle = 201
        s.plot_velocity = True
        s.plot_psi = False
        s.plot_Kp = True
    elif cas == 2:
        # Joukowski's airfoil, airfoil plotted under zero incidence 
        # parameter : Zc
        s.map = "Karman-Trefftz"
        s.alpha = 3  # -5.7391704773
        s.Zc = -0.05 + 1j * 0.05
        s.opt_levels = "manuel"
        s.levels_iso = [-2., 2., 0.1]
        s.k = 1.9  # exponent of the transformation
        s.plot_velocity = True
        s.adimChord = False

    elif cas == 3:
        # Non symmetrical von Mises's airfoil
        # parameters: eps (epsilon) and k  
        # mu if found from  Zc
        # The airfoil is the one found in the book
        # it works only with a=2 and mu_0+mu_1+a=0
        s.map = "von Mises"
        s.alpha = 0.
        s.k = 1.9
        s.eps = 0.05
        s.Zc = -s.eps * (1 - 1j)
        s.levels_iso = [-1, 1, 0.1]

    elif cas == 4:
        # Symmetrical von Mises's airfoil
        # parameters: eps (epsilon) and k   
        # mu if found from  Zc
        # The airfoil is the one found in the boo
        # it works only with a=2 and mu_0+mu_1+a=0
        s.map = "von Mises"
        s.alpha = 0.
        s.k = 2.0
        s.eps = 0.1
        s.Zc = 0
        s.levels_iso = [-1, 1, 0.1]

    elif cas == 5:
        # van de Vooren's airfoil
        # parameters: eps (epsilon) and k    

        s.map = "van de Vooren"
        s.alpha = 0.
        s.k = 2
        s.eps = 0.05
        s.Zc = 0 + 1j * 0  # -s.eps*(1-1j)
        s.levels_iso = [-1, 1, 0.1]


    elif cas == 6:
        # double sharp airfoil,
        # parameters : a et b
        s.n_circle = 101
        s.map = "double sharp"
        s.alpha = 0.
        s.b = s.R / 2 * np.exp(1j * np.pi / 3)
        s.Zc = 0 + 1j * 0
        # s.eps           = 0.2
        # beta            = 110
        # s.b             =  np.sqrt(s.eps)*s.R*np.exp(1j*np.deg2rad(beta))
        s.Psi_method = "Psi"
        s.camberline = True

    #  case 'P2'
    #      b = a/2;
    #      Z_c = 0;

    #  case 'P3'
    #      beta = 20*pi/180;
    #      epsilon = 0.2;
    #      b = sqrt(epsilon) * a * exp(i*beta)
    #      Z_c = 0;

    elif cas == 7:
        # double sharp leading and trailing edge
        # parameters : a and b
        # Problem at the leading edge
        s.map = "double sharp"
        s.alpha = 0.
        s.eps = 0.2
        beta = np.rad2deg(np.pi / 2 + 0.1)
        s.b = np.sqrt(s.eps) * s.R * np.exp(1j * np.deg2rad(beta))
        s.Zc = 0 + 1j * 0
        s.levels_iso = [-1, 1, 0.1]

    elif cas == 8:
        # double sharp leading and trailing edge
        # parameters : a and b
        # 180° > arg(b)  > 90° 
        s.camberline = False
        s.map = "double sharp"
        s.alpha = 0.
        s.eps = 0.05
        beta = 135
        s.b = np.sqrt(s.eps) * s.R * np.exp(1j * np.deg2rad(beta))
        s.Zc = s.b  # -s.eps*(1-1j)
        s.levels_iso = [-1, 1, 0.1]

        # case 'P1a'      
        #      beta = 180*pi/180;
        #      b = 0.1*exp(i*beta);
        #      Z_c = b;

        #  case 'P1b'  
        #      beta = 135*pi/180;
        #      b = 0.1*exp(i*beta);
        #      Z_c = b;

    s.run_airfoil()
    plt.show()