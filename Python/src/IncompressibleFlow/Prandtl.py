#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Module : incompressibleFlow.Prandtl , Prandtl'wing theory and applications
    ..

* functions
* intregration of Prandtl' equation
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


def vortex_velocity(X, x0, ya, yb, Gamma):
    """
    induced velocity from a horseshoe vortex with A=(x0,y_A) and B=(x0,y_B)

    Args:
	    X (float) : =(x,y)
	    x0 (float) : x-coordinate
	    ya (float) : y-coordinate
	    yb (float) : y-coordinate
	    Gamma (float) : circulation

    Returns:
        real : w , induced velocity

	"""

    x, y = X[0], X[1]
    g = Gamma / (4 * np.pi)

    if x == x0 and y == (ya + yb) / 2:
        return -Gamma / np.pi / (yb - ya)
    elif x == x0:
        return -g * (1 / (y - ya) - 1 / (y - yb))
    else:
        db = np.sqrt((x - x0) ** 2 + (y - yb) ** 2)
        da = np.sqrt((x - x0) ** 2 + (y - ya) ** 2)
        w1 = g / (x - x0) * ((y - yb) / db - (y - ya) / da)
        w2 = g / (y - yb) * (1 + (x - x0) / db)
        w3 = -g / (y - ya) * (1 + (x - x0) / da)
        w = w1 + w2 + w3
    # w4=-g / (y-ya)*(1 +da/(x-x0)) + g / (y-yb)*(1 + db /(x-x0))
    # w4 is equal to  w, verification
    return w


def Prandtl_E_OGE(Nc, lamb, b, l0, k, y, theta, alpha):
    """
    Laurent serie coefficients for en elliptic wing, without ground effect

    OGE: out ground effect

    Args:
         Nc (int) : number of collocation points
         lamb (float) : wing aspect ratio
         b (float) : span = 2 x b
         l0 (float) : root wing chord
         k (float) : half lift gradient of the airfoil
         theta (float) : angle to define the spanwise coordinate y
         alpha (float) : angle of attack in radians, assumed constant

    Returns:
        float : A_n coefficients
    """
    l = l0 * np.sin(theta)  # chord law in the spanwise direction
    My = len(theta)
    mat = np.zeros((My - 2, Nc), dtype=float)
    RHS = np.zeros((My - 2), dtype=float)
    u = np.arange(1, My - 1)
    for n in range(Nc):
        mat[:, n] = 4 * b * np.sin((n + 1) * theta[u]) + k * l[u] * (n + 1) * np.sin((n + 1) * theta[u]) / np.sin(
            theta[u])
    RHS = k * alpha * l[u]
    A, res, rank, s = np.linalg.lstsq(mat, RHS, rcond=-1)

    if Nc == My - 2:
        print(np.linalg.solve(mat, RHS))
    print('A_n = ', A)
    return A


def integrand(x, theta, n, h_over_b):
    """
    integrand of the integrel with ground effect  (IGE)
    """
    tmp = np.cos(x) - np.cos(theta)
    return (tmp * np.cos(n * x)) / (tmp ** 2 + (2 * h_over_b) ** 2)


def solve_integralO2(theta, n, h):
    """
    order 2 integral, trapeze rule
    """
    s = np.zeros(len(theta), dtype=float)
    m = 720 * 8
    phi = np.linspace(0, np.pi, m)
    dphi = phi[1] - phi[0]
    for i in range(len(theta)):
        val = integrand(phi, theta[i], n, h)
        s[i] = dphi * (np.sum(val[1:m - 2]) + (val[0] + val[m - 1]) / 2)
    return s


def solve_integralO4(theta, n, h_over_b):
    """
    Integral with quadrature rule
    """
    s = np.zeros(len(theta), dtype=float)
    for i in range(len(theta)):
        s[i], res = integrate.quad(integrand, 0, np.pi, args=(theta[i], n, h_over_b,))
    return s


def Prandtl_E_IGE(Nc, lamb, b, l0, k, h, y, theta, alpha):
    """
    Laurent serie coefficients for en elliptic wing, in ground effect

    IGE: in ground effect

    Args:
         Nc (int) : number of collocation points
         lamb (float) : wing aspect ratio
         b (float) : span = 2 x b
         l0 (float) : root wing chord
         h (float) : distance of the wing from the ground  / b
         k (float) : half lift gradient of the airfoil
         theta (float) : angle to define the spanwise coordinate y
         alpha (float) : angle of attack in radians, assumed constant

    Returns:
        float: A_n coefficients
        """
    # l=l0*np.sin(theta)
    My = len(theta)
    mat = np.zeros((My - 2, Nc), dtype=float)
    RHS = np.zeros((My - 2), dtype=float)
    A = np.zeros((Nc), dtype=float)
    G = np.zeros((My - 2), dtype=float)
    u = np.arange(1, My - 1)
    for n in range(Nc):
        coef = 2 * k * (n + 1) / (np.pi * lamb)
        G = solve_integralO4(theta[u], n + 1, h)
        # test with different integral method
        # G1=solve_integralO2(theta[u],n+1,h)
        # t1,t2=np.linalg.norm(G),np.linalg.norm(G1)
        # print('error = ',t2/t1-1)
        mat[:, n] = (1 + coef) * np.sin((n + 1) * theta[u]) - coef / np.pi * G * np.sin(theta[u])

    RHS = 2 * k * alpha / (np.pi * lamb) * np.sin(theta[u])

    A, res, rank, s = np.linalg.lstsq(mat, RHS, rcond=-1)

    if Nc == My - 2:
        print(np.linalg.solve(mat, RHS))
    # print('A_1 = %3.15f'%(A[0]))
    return A
