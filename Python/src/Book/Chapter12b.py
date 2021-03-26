#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 12 exercises  (cont')
    ..
"""
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from copy import deepcopy

width, height = 10, 8
Hequal = 4


def Exercice12_3():
    """
    Nozzle with a wall given by a polynomials
    exercise 1 geometry
    """
    opt_wall = 2
    n = 11
    if opt_wall == 1:
        y0 = 1.
        L = 5 * y0
        h = 2 * y0
        hbar = h / L
        y0bar = y0 / L
        x = np.linspace(0, L, n)
        p0 = np.deg2rad(10)
    else:
        y0bar = 1
        p0 = 0.05
        hbar = 0
        L = 1.
        x = np.linspace(1, 8.5, n) ** 4 / 100000

    def wall(opt, eta, y0bar, p0, hbar):
        """
        wall equation
        eta = x/L, hbar = h/L
        """
        if opt == 1:
            y = y0bar + p0 * eta - (2 * p0 - 3 * hbar) * eta ** 2 + (p0 - 2 * hbar) * eta ** 3
        else:
            y = p0 / 2 * eta ** 2 + y0bar
        return y

    fig = plt.figure(1, figsize=(width, Hequal))
    fig.suptitle('Polynomial nozzle', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$y/L$', fontsize=14)
    ax.set_xlabel(r'$x/L$', fontsize=14)
    ax.plot(x / L, wall(opt_wall, x / L, y0bar, p0, hbar), 'k-')
    plt.grid()
    plt.show()

    print(np.rad2deg(p0 * x / L))
