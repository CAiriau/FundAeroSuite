# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: `Pr. L.A. Barba <https://github.com/barbagroup/AeroPytho>`_

.. codeauthor:: modidied by C. Airiau

*date : 2017*

*license : AGPL-3.0*

This module contains various functions related to potential flow and panel methods
  ..
  
It is a small introduction to panel method

"""

import numpy    as np
import numpy.ma as ma
import matplotlib.pyplot as plt

plot_panels = False         # to show panels
plot_gamma = False          # to plot circulation
plot_vitesse = False        # to plot velocity
plot_psi = False            # to plot iso streamlines


def reference_plaque_plane(x, alpha):
    """
    
    reference solution on a flat plate

    Args:
        x (real): :math:`0< x < 1`
        alpha (real): angle of attack in radians
    
    Returns:
        real: wall velocity
    """
    return np.cos(alpha) + np.sin(alpha) * np.sqrt((1 - x) / x)


def potentiel_complexe(z, xP, Gamma, alpha, U0=1):
    """
    Complex potential of n vortex
    
    Args:
        z (complex): location
        xP (real): Point vortex position
        Gamma (real): vortex circulation
        alpha (real): angle of attack in radians
        U0 (real): upstream velocity
        
    Returns:
        complex: :math:`f` complex potential
    """

    if np.isscalar(z):
        return np.sum(-1J * Gamma / (2 * np.pi) * np.log(z - xP)) + U0 * z * np.exp(-1j * alpha)
    else:
        print("size z in pot. compl. : ", z.shape)
        f = np.zeros(z.shape) + 1j * np.zeros(z.shape)
        for i in range(len(Gamma)):
            f += -1J * Gamma[i] / (2 * np.pi) * np.log(z - xP[i])
        return f + U0 * z * np.exp(-1j * alpha)


def complex_velocity(z, xP, Gamma, alpha, U0=1):
    """
    Potential complex velocity
    
    Args:
        z (complex): location
        xP (real): Point vortex position
        Gamma (real): vortex circulation
        alpha (real): angle of attack in radians
        U0 (real): upstream velocity
        
    Returns:
        complex: :math:`w` complex potential velocity
        
    """

    if np.isscalar(z):
        return np.sum(-1J * Gamma / (2 * np.pi) / (z - xP)) + U0 * np.exp(-1j * alpha)
    else:

        w = np.zeros(z.shape) + 1j * np.zeros(z.shape)
        for i in range(len(Gamma)):
            w += -1J * Gamma[i] / (2 * np.pi) / (z - xP[i])
        return w + U0 * np.exp(-1j * alpha)


def get_contours(mplcont):
    """
    To get the lines calculated with plt.contour
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


def dessiner_fonction_courant(zgrid, Psi, levels, titre, chord=1):
    """
    to plot streamlines
    """
    fz = 20
    iso_psi = levels.tolist()
    plt.figure(figsize=(10, 10))
    plt.title(r'iso-%s, $\alpha = $ %4.2f °' % (titre, np.rad2deg(alpha)), fontsize=fz, fontweight='bold')
    psi = plt.contour(zgrid.real, zgrid.imag, Psi, levels=iso_psi, colors='blue')
    plt.plot([0, chord], [0, 0], 'r-', linewidth=5)
    psi_x, psi_y = get_contours(psi)
    plt.plot(psi_x, psi_y)
    plt.legend()


def montrer_gamma(xP, Gamma, chord=1):
    """
    to plot :math:`\Gamma` circulation distribution
    """
    plt.figure(figsize=(15, 15))
    plt.title("panel")
    plt.plot([0, chord], [0, 0], 'k--', linewidth=1)
    plt.plot(xP, Gamma, 'rs', markersize=5)
    # plt.ylim(-0.01,0.05)


def montrer_panneaux(A, xP, xC, chord=1):
    """
    to show panels
    """
    plt.figure(figsize=(15, 3))
    plt.title("panneaux")
    plt.plot([0, chord], [0, 0], 'k--', linewidth=1, label='Plaque')
    plt.plot(A, np.zeros(len(A)), 'ko', markersize=10, label='limite panneau')
    plt.plot(xP, np.zeros(len(xP)), 'bs', markersize=5, label='tourbillons')
    plt.plot(xC, np.zeros(len(xC)), 'rs', markersize=5, label='contrôle')
    plt.ylim(-0.01, 0.05)
    plt.legend(loc="best")


def plaque_plane(nP, y_wall=0.01, nv=1001, chord=1, U0=1, alpha=0.1):
    """
    Flat plate flow
    
    Args:
        nP (integer): panel number
        y_wall (real): distance from the wall
        nv (integer): number of discrete value in x for z
        chord (real): chord 
        U0 (real): upstream velocity
        alpha (real): angle of attack
    
    Outputs:
        * real array - xP : panel 1/4 abscisse 
        * real array - xC : panel 3/4 absicsse
        * real array - Gamma : panel circulation
        * complex array - z_wall : positions where are calculated velocity
        * complex - w_wall : complexe veloicty
        * real - Gamma0 : total circulation
        
    """
    A = np.linspace(0, chord, nP + 1)  # location of the first panel point
    print("Panel positions                          :\n", A)

    Lp = A[1] - A[0]  # panel length
    xP = A[0:-1] + 1 / 4 * Lp  # vortex location (1/4 panel length)
    xC = A[0:-1] + 3 / 4 * Lp  # control point location (3/4 panel length)
    print("xP                                       :\n", xP)

    # Circulation distribution :
    Mat = np.zeros((nP, nP))
    for j in range(nP):
        Mat[:, j] = 1 / (xC - xP[j])  # influence matrix
    Rhs = -np.ones(nP) * U0 * np.sin(alpha) * 2 * np.pi

    Gamma = np.linalg.solve(Mat, Rhs)  # linear system solution
    print("Gamma                                    :\n", Gamma)
    Gamma0 = -np.pi * U0 * chord * np.sin(alpha);
    print("Gamma 0                                  :\n", Gamma0)
    if plot_panels:
        montrer_panneaux(A, xP, xC)

    if plot_gamma:
        montrer_gamma(xP, Gamma)
    # Verification of the solution , we must find -1 below:
    print("sum(Gamma)/(pi sin alpha) =-1            : ", np.sum(Gamma) / (np.pi * np.sin(alpha)))

    wC = complex_velocity(xC, xP, Gamma, alpha)
    print("Collocation point velocity : ", wC)

    z_wall = np.linspace(-0.5, 1.5, nv) + 1j * y_wall * np.ones(nv)
    w_wall = complex_velocity(z_wall, xP, Gamma, alpha)
    return [xP, xC, Gamma, z_wall, w_wall, Gamma0]


def dessiner_contours(xP, Gamma, alpha=0.1):
    """
    To plot iso velocity lines and streamlines 
    """
    npt = 101
    n_iso = 101
    Xgrille = np.linspace(-1, 2, npt)
    Ygrille = np.linspace(-0.5, 0.5, npt)

    X, Y = np.meshgrid(Xgrille, Ygrille)
    Zgrid = X + 1j * Y
    Zgrid = ma.masked_where(np.absolute(Zgrid.imag) == 0.00, Zgrid)
    f = potentiel_complexe(Zgrid, xP, Gamma, alpha)
    w = complex_velocity(Zgrid, xP, Gamma, alpha)

    Psi = f.imag
    Psi_limit = []
    Psi_limit.append(np.min(Psi))
    Psi_limit.append(np.max(Psi))
    print("Psi, min,max                              :", Psi_limit)
    if plot_psi:
        dessiner_fonction_courant(Zgrid, Psi, np.linspace(Psi_limit[0], Psi_limit[1], n_iso), 'r$\Psi$')
    v_limit = []
    v_limit.append(np.min(-w.imag))
    v_limit.append(np.max(-w.imag))

    if plot_vitesse:
        dessiner_fonction_courant(Zgrid, -w.imag, np.linspace(v_limit[0], Psi_limit[1], n_iso), 'v')
