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
from Tools.misc import set_title, set_question, set_info

Hauteur = 8
Largeur = 10


def Helmbold(lamb, k=np.pi):
    """
    Helmbold law : d C_L/ d \alpha
    """
    x = 2 * k / (np.pi * lamb)
    return 2 * k / (x + np.sqrt(1 + x ** 2))


def Exercice6_1():
    """ 
    Lift coefficient law with  the wing aspect ratio
    """
    set_title("Lift coefficient law with  the wing aspect ratio")

    print("Helmbold's law")

    N = 101                                 # number of points for lambda
    lambda_Max = 25                         # maximal value on the plot
    lamb = np.linspace(0.1, lambda_Max, N)  # wing aspect ratio
    k = np.array([2.5, 2.8, np.pi])

    # plot
    fig = plt.figure(1, figsize=(Largeur, Hauteur))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$\lambda$', fontsize=16)
    ax.set_ylabel(r'$C_{L_\alpha}$', fontsize=20)
    ax.set_title("Helmbold's law")
    lamb1 = np.linspace(0, 6)
    lamb2 = np.linspace(5, lambda_Max)
    for k in k:
        ax.plot(lamb, Helmbold(lamb, k), linewidth=2, label=r'$k=$%2.1f' % (k))
        ax.plot([0, lambda_Max], [2 * k, 2 * k], linewidth=2, label=r'$k=$%2.1f' % (k))

    ax.plot(lamb1, np.pi * lamb1 / 2, color='r', linestyle='-', linewidth=2, label=r'$\pi \lambda/2$')
    ax.plot(lamb2, np.pi * 2 * (1 - 2 / lamb2), linewidth=2, label=r'approx')

    xmax, ymax = lambda_Max, 8
    ax.axis([0, xmax, 0, ymax])
    xmajor_ticks, xminor_ticks = np.arange(0, xmax, 2), np.arange(0, xmax, 0.5)
    ymajor_ticks, yminor_ticks = np.arange(0, ymax, 1), np.arange(0, ymax, 1)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.grid()
    ax.legend(loc='lower right')
    plt.show()


def Kp_Schlichting(x, epsilon):
    """
    Kp on a  ellipsoide, Schlichting's law
    `see <http://docs.desktop.aero/appliedaero/potential3d/Bodies.html>`_
    """

    c = np.sqrt(1 - (2 * epsilon) ** 2)  # epsilon = a/L (radius/length)
    eta = 2 * (1 - c ** 2) / c ** 3 * (np.arctanh(c) - c)
    U_max = 2. / (2 - eta)
    return 1 - U_max ** 2 * (1 - (2 * x - 1) ** 2) / (1 - c ** 2 * (2 * x - 1) ** 2)


def Kp_VD(x, epsilon):
    """
    analycal solution
    """
    return 4 * epsilon ** 2 * np.log(epsilon ** 2) + epsilon ** 2 * (2 - (1 - 2 * x) ** 2) / (x * (1 - x))


def compare_van_Driest():
    """
    Comparison with experimental values provided by Van Driest
    """
    # x=np.array([0.05,0.17,0.36,0.58])
    # Kp=np.array([-0.010,-0.053,-0.061,-0.026])
    x = np.array([0.04944156, 0.09877804, 0.17424668, 0.26218079, 0.37486144])
    Kp = np.array([-0.01011674, -0.03771612, -0.05568793, -0.05826931, -0.06136363])
    leg = ["A", "B", "C", "D", "E"]
    return x, Kp, leg


def analytical_Kp_min(epsilon):
    """
    -Kp minimal for the asymptotic solution
    """
    return np.abs(8 * epsilon ** 2 * (1 + np.log(epsilon)))


def Exercice6_3():
    """
    Ellipsoide
    """
    n = 121
    epsi = 1e-3
    x = np.linspace(0.+epsi, 1-epsi, n)
    eps = np.array([0.05, 0.1, 0.2])

    set_question("2: Kp : exact and asymptotic solutions")

    for ifig in range(len(eps)):
        fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
        ax = fig.add_subplot(111)
        ax.set_xlabel(r'$x/\ell$', fontsize=16)
        ax.set_ylabel(r'$Kp$', fontsize=20)
        ax.set_title(r'$Kp$ : analytical and Schlichting, $\varepsilon$ = %2.2f' % (eps[ifig]))
        ax.axis([0, 1, -0.4, 1])
        ax.plot(x, Kp_Schlichting(x, eps[ifig]), color='k', linestyle='-', linewidth=2, label='Schlichting')
        ax.plot(x, Kp_VD(x, eps[ifig]), color='k', linestyle='--', linewidth=2, label='Asymptotic')
        ax.legend(loc='upper center')

    set_question("3: Kp min: exact and asymptotic solutions")

    ifig += 1

    epsil = np.logspace(-3, -1, 10)
    print(epsil)
    Kpmin = []
    for eps in epsil:
        Kpmin.append(-Kp_Schlichting(0.5, eps))

    print(Kpmin)

    fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$ \varepsilon=h/\ell$', fontsize=16)
    ax.set_ylabel(r'$-Kp_{\min}$', fontsize=20)
    ax.set_title(r'$Kp$ : analytical and Schlichting')
    ax.loglog(epsil, Kpmin, color='k', linestyle='-', linewidth=2, label='Schlichting')
    ax.loglog(epsil, analytical_Kp_min(epsil), color='k', linestyle='--', linewidth=2, label='Asymptotic')
    ax.grid(which='both')
    ax.legend(loc='upper center')

    set_question("4: Kp, Van Driest's experiment data")

    ifig += 1
    eps0 = 0.075
    x0, kp0, leg0 = compare_van_Driest()
    print(x0, kp0)

    fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$x/\ell$', fontsize=16)
    ax.set_ylabel(r'$Kp$', fontsize=20)
    ax.set_title(r'$Kp$ : analytical  and experiment , $\varepsilon$ = %2.3f' % (eps0))
    for k in range(len(x0)):
        ax.plot(x0[k], kp0[k], linestyle='', marker='o', markersize=6, label=r'Exp. point :%s' % (leg0[k]))
    ax.plot(x, Kp_Schlichting(x, eps0), color='k', linestyle='-', linewidth=2, label='Schlichting')
    ax.axis([0, 1, -0.1, 0.1])
    ax.plot(x, Kp_VD(x, eps0), color='k', linestyle='--', linewidth=2, label='Asymptotic')
    ax.legend(loc='upper center')
    ax.grid()

    plt.show()


def Exercice6_4():
    """
    Data analysis of  Van Driest phd about the ellipsoidal flo
    
    *it is not an exercise of the book !*
    """
    set_info("it is not an exercise of the book ! \n" +
             "Only an analysis of Van Driest's thesis !")
    ifig = 0
    l = 100
    h = 15
    fichiers = ["vanDriest_exp1.dat", "vanDriest_exp2.dat", "vanDriest_exp2_ogive.dat"]

    def lire_fichier(nom):
        """
        read the file
        """
        filepath = os.path.join('Book/Data', nom)
        print('File to read : ', filepath)
        with open(filepath, 'r') as infile:
            A = np.loadtxt(infile, dtype=float, unpack=False, skiprows=1)
        print(A.shape)
        print(A)
        return A[:, 0], A[:, 1]

    def display_curve(nfig, x, y, titrex, titrey):
        """
        plots
        """
        plt.figure(num=nfig, figsize=(Largeur, Hauteur))
        plt.xlabel(titrex, fontsize=16)
        plt.ylabel(titrey, fontsize=20)
        plt.plot(x, y, 'ko', linestyle='--', markersize=5, linewidth=2, label='van Driest exp.')
        plt.grid(which='both')
        plt.legend(loc='best')

    x0, Kp = lire_fichier(fichiers[0])
    display_curve(ifig, x0, Kp, r"$X/\ell$", r"$Kp$")

    ifig += 1
    x1, y = lire_fichier(fichiers[1])
    x2, yo1 = lire_fichier(fichiers[2])
    display_curve(ifig, x1, y, r"$X$", r"$y_{ogive}$")
    display_curve(ifig, x2, yo1, r"$X$", r"$y_{ellipsoide}$")

    def ogive_geometry(x, l, h):
        """
        ogive geometry
        """
        y = h * np.sqrt(1 - (x / l) ** 2)
        return y

    X = np.linspace(-l, l, 101)
    Y = ogive_geometry(X, l, h)

    plt.figure(ifig, figsize=(Largeur, Hauteur))
    plt.plot(X, Y, 'r-', linestyle='--', markersize=5, linewidth=2, label='reference geometry')
    plt.plot(x0 * l, ogive_geometry(x0 * l, l, h), 'bs', linestyle='--', markersize=7, linewidth=2, \
             label='exp.')
    plt.legend(loc='lower center')

    xf = (1 + x0) / 2
    yf = Kp
    ifig += 1
    display_curve(ifig, xf, yf, r"$x/c$", r"$Kp$")
    for x, kp in zip(xf, Kp):
        print("x = %2.3f,\t Kp = %2.4f" % (x, kp))

    plt.show()

