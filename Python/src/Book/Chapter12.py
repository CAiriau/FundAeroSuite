#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 12 exercises 
    ..
"""

import numpy as np
from scipy import interpolate, integrate
import matplotlib.pyplot as plt
from copy import deepcopy
from CompressibleFlow.fonctions import gamma, p_pi, omega_mu, mu_omega, mu2mach, inv_S_over_Scrit
from CompressibleFlow.fonctions import S_over_Scrit, line_intersection, complex_line_intersection, omega_super
from Tools.misc import set_title, set_question

Largeur, Hauteur = 10, 8
Hequal = 4


def lower_wall(x, type_wall):
    """
    Equation  of the lower wall
    wall type  :
    1 : constante slope 
    2 : y=0
    3 : with a given equation
    """

    def f(opt, x, L):
        """
        Equation  of the lower wall
        """
        if opt == 0:
            return (np.cos(np.pi * x / L) - 1) / 2
        elif opt == 1:
            return -np.pi / (2 * L) * np.sin(np.pi * x / L)

    y = 0
    if type_wall == 1:      #  diverging,  constant slope
        slope = -np.tan(5 / 180 * np.pi)
        y = slope * x
    elif type_wall == 2:    # converging,  constant slope
        y = 0. * x
    else:
        L = 10
        lim = 10
        val = f(0, lim, L)
        slope = f(1, lim, L)
        if x <= lim:
            y = f(0, x, L)
        else:
            y = val + slope * (x - lim)
    return y


def upper_wall(x, type_wall):
    """
    Equation  of the upper wall
    """

    def f1(opt, x, L):
        """
        Equation  of the upper wall
        """
        if opt == 0:
            return (3 - np.cos(np.pi * x / L)) / 2
        elif opt == 1:
            return np.pi / (2 * L) * np.sin(np.pi * x / L)

    y = 0
    if type_wall == 1:
        slope = np.tan(5 / 180 * np.pi)
        y = 1 + slope * x
    elif type_wall == 2:
        slope = np.tan(-1 / 180 * np.pi)
        y = 1 + slope * x
    elif type_wall == 3:
        lim = 6
        if x <= lim:
            y = 1 - 0.05 * np.sin(np.pi * x / 6) ** 2
        else:
            y = 1
    else:
        L = 10
        lim = 10
        val = f1(0, lim, L)
        slope = f1(1, lim, L)
        if x <= lim:
            y = f1(0, x, L)
        else:
            y = val + slope * (x - lim)
    return y


def Exercice12_1():
    """ 
    Diverging channel
    
    2D channel compressible isentropic flow with the method of caracteristics (MoC) 
    for an ideal gas   
    
    Args :
     
        Mach0 :  channel entrance Mach number 
        N     :  number of positive and negative slope characteristics 
        Ltot  :  channel length
        lower_wall: lower wall geometry option
        upper_wall: upper wall geometry option
    """

    # ****************************************************************
    #  PARAMETERS et DATA
    # ****************************************************************
    cas = 0

    test_omega_mu = False
    show_geometry = False
    show_Mach = True
    show_Median = False
    show_label = False

    Mach_ref = [2, 1.8, 1.3, 2]
    lower_wall_data = [1, 2, 2, 4]
    upper_wall_data = [1, 2, 3, 4]
    Mlim = [[2.0, 3.2], [1.55, 1.8], [1.1, 1.4], [2, 3.5]]
    nlev = [21, 21, 21, 16]

    Nwall = 101             # number of intervals to discretize the wall
    N = 7                   # number of characteristic  to draw, need to be odd (7 or 41 is adviced)
    Ltot = 10.              # total channel length for calculations
    Ls = 10.                # length displayed
    nlevels = nlev[cas]     # level number for iso Mach lines

    lower_wall_type, upper_wall_type = lower_wall_data[cas], upper_wall_data[cas]
    Mach0 = Mach_ref[cas]   # Entrance Mach number
    Mmin, Mmax = Mlim[cas][0], Mlim[cas][1]     # Mach scaling on the figures

    gam = gamma
    epsilon = 1e-6              # step for derivatives
    Nmedian = int((N - 1) / 2)  # calculus will start on the axis for median value
    IterMax = 10                # maximal number of iterations in the Newton method
    Precision = 1e-14           # accuracy in  Newton function
    # print('Nmedian = ',Nmedian)

    # ****************************************************************
    #  INITIALIZATIONS
    # ****************************************************************

    set_title("Characteristic method in a diverging channel")

    if test_omega_mu:
        set_question('Test of  omega and mu functions')
        mu = np.arcsin(1. / Mach0)
        omega = omega_mu(mu)
        mu1 = mu_omega(omega)
        print('M0= ', Mach0, 'mu = ', mu, 'Omega = ', omega, 'mu = ', mu1)
        return

    # Table declaration
    XP = np.zeros((N), dtype=float)
    lambdaPlusP, lambdaMinusP = np.zeros((N), dtype=float), np.zeros((N), dtype=float)
    thetaP = np.zeros((N), dtype=float)
    XI, YI = np.zeros((N - 1), dtype=float), np.zeros((N - 1), dtype=float)
    lambdaPlusI, lambdaMinusI = np.zeros((N - 1), dtype=float), np.zeros((N - 1), dtype=float)
    omegaI, thetaI = np.zeros((N - 1), dtype=float), np.zeros((N - 1), dtype=float)
    muI = np.zeros((N - 1), dtype=float)

    # Lists declaration 
    MachP = []          # Mach on the main axis
    XcP, YcP = [], []   # location of this Mach on the axis
    Xmesh, Ymesh, mumesh = [], [], []
    X_lower, X_upper, mulower, muupper = [], [], [], []

    #  Initialization : main initial  points
    Ymin0, Ymax0 = lower_wall(0., lower_wall_type), upper_wall(0., upper_wall_type)
    YP = Ymin0 + np.arange(N) * (Ymax0 - Ymin0) / float(N - 1)
    yMedian = (Ymin0 + Ymax0) / 2
    print('YP(0) = ', YP, 'Ymedian = ', yMedian, 'Nmedian = ', Nmedian, 'Y(Nmedian) = ', YP[Nmedian])

    # print('XP = ',XP)
    mu0 = np.arcsin(1 / Mach0)
    omega0 = omega_mu(mu0)

    # Entrance section :
    muP = mu0 * np.ones((N), dtype=float)
    omegaP = omega0 * np.ones((N), dtype=float)
    lambdaPlusP = omegaP - thetaP
    lambdaMinusP = omegaP + thetaP

    # on the main axis :
    MachP.append(Mach0)
    XcP.append(XP[0])
    YcP.append(YP[Nmedian])

    #  table preparation to present results
    # Xmesh,Ymesh = XP,YP
    # mumesh = muP
    Xmesh.append(deepcopy(XP))
    Ymesh.append(deepcopy(YP))
    mumesh.append(deepcopy(muP))

    X_lower.append(XP[0])
    X_upper.append(XP[-1])
    mulower.append(muP[0])
    muupper.append(muP[-1])

    # ****************************************************************
    #  FIGURES, WALLS
    # ****************************************************************
    # to draw walls
    XL = np.arange(Nwall) * Ltot / float(Nwall - 1)
    print("x - step  : ", XL[1] - XL[0])
    Y_lower = [lower_wall(eta, lower_wall_type) for eta in XL]
    Y_upper = [upper_wall(eta, upper_wall_type) for eta in XL]

    fig = plt.figure(0, figsize=(10, Hequal))
    fig.suptitle('Nozzle', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$y$', fontsize=14)
    ax.set_xlabel(r'$x$', fontsize=14)
    ax.plot(XL, Y_lower, 'k-', linewidth=1, label='lower')
    ax.plot(XL, Y_upper, 'k-', linewidth=1, label='upper')
    # ax.axis([0,3,-0.5,2])
    ax.axis('equal')
    ax.set_xlim(0, Ls)
    # ax.grid('')
    if show_geometry:
        plt.show()
        return

    # internal functions
    def test_position(X, Y):
        """
        to test the location of intermediate points
        """
        m = len(Y)
        I1, I2 = [], []
        for k in range(m):
            if Y[k] < lower_wall(X[k], lower_wall_type):
                I1.append(k)
            elif Y[k] > upper_wall(X[k], upper_wall_type):
                I2.append(k)
        if I1:
            print('Problem with lower wall and indexes', I1)
        if I2:
            print('Problem with upper wall and indexes', I2)

    # ****************************************************************
    #  LOOP for x MARCHING  
    # ****************************************************************
    # 
    Xmin = 0  # starting on x=0

    while Xmin < Ltot:

        # --------------------------------------------
        # Build intermediate points
        # --------------------------------------------
        for j in range(N - 1):
            XI[j] = (YP[j] - YP[j + 1] - np.tan(thetaP[j] + muP[j]) * XP[j]
                     + np.tan(thetaP[j + 1] - muP[j + 1]) * XP[j + 1]) / (
                                np.tan(thetaP[j + 1] - muP[j + 1]) - np.tan(thetaP[j] + muP[j]))
            YI[j] = YP[j] + np.tan(thetaP[j] + muP[j]) * (XI[j] - XP[j])
            lambdaPlusI[j] = lambdaPlusP[j]
            lambdaMinusI[j] = lambdaMinusP[j + 1]
            omegaI[j] = (lambdaPlusI[j] + lambdaMinusI[j]) / 2
            thetaI[j] = (lambdaMinusI[j] - lambdaPlusI[j]) / 2
            muI[j] = mu_omega(omegaI[j])

            #  plot of the characteristic lines
            ax.plot([XP[j], XI[j]], [YP[j], YI[j]], 'r-')
            ax.plot([XP[j + 1], XI[j]], [YP[j + 1], YI[j]], 'b-')

        test_position(XI, YI)

        # ------------------------------------------------------ 
        # Build the main point at the next step   
        # Central points (on the axis)  
        # ------------------------------------------------------ 
        for j in np.arange(1, N - 1):
            XP[j] = (YI[j - 1] - YI[j] - np.tan(thetaI[j - 1] + muI[j - 1]) * XI[j - 1]
                     + np.tan(thetaI[j] - muI[j]) * XI[j]) / (
                                np.tan(thetaI[j] - muI[j]) - np.tan(thetaI[j - 1] + muI[j - 1]))
            YP[j] = YI[j - 1] + np.tan(thetaI[j - 1] + muI[j - 1]) * (XP[j] - XI[j - 1])
            lambdaPlusP[j] = lambdaPlusI[j - 1]
            lambdaMinusP[j] = lambdaMinusI[j]
            omegaP[j] = (lambdaPlusP[j] + lambdaMinusP[j]) / 2
            thetaP[j] = (lambdaMinusP[j] - lambdaPlusP[j]) / 2
            muP[j] = mu_omega(omegaP[j])
            lambdaPlusP[j] = omegaP[j] - thetaP[j]
            lambdaMinusP[j] = omegaP[j] + thetaP[j]
            # plot of the characteristic lines
            ax.plot([XP[j], XI[j]], [YP[j], YI[j]], 'b-')  # C- intermediate
            ax.plot([XP[j], XI[j - 1]], [YP[j], YI[j - 1]], 'r-')  # C+ intermediate

        # ------------------------------------------------------ 
        # First point of the lower wall (Newton's method)
        # ------------------------------------------------------ 
        x = XI[0]
        dx = 1
        counter = 0
        while (np.abs(dx) > Precision) and (counter < IterMax):
            counter += 1
            FF = YI[0] + np.tan(thetaI[0] - muI[0]) * (x - XI[0]) - lower_wall(x, lower_wall_type)
            dFFdx = np.tan(thetaI[0] - muI[0]) - (
                        lower_wall(x + epsilon, lower_wall_type) - lower_wall(x - epsilon, lower_wall_type)) / (
                                2 * epsilon)
            dx = FF / dFFdx
            x -= dx
        if counter == IterMax: print('BE CAREFUL : non-convergence of  Newton for point P1')
        # print("counter (i=0) = %i, \t dx = %e"%(counter, dx))
        XP[0] = x
        YP[0] = lower_wall(XP[0], lower_wall_type)
        thetaP[0] = np.arctan(
            (lower_wall(XP[0] + epsilon, lower_wall_type) - lower_wall(XP[0] - epsilon, lower_wall_type)) / (
                        2 * epsilon))
        omegaP[0] = lambdaMinusI[0] - thetaP[0]
        muP[0] = mu_omega(omegaP[0])
        lambdaPlusP[0] = omegaP[0] - thetaP[0]
        ax.plot([XP[0], XI[0]], [YP[0], YI[0]], 'b')  # end of C-

        # ------------------------------------------------------ 
        # last point on the upper wall (Newton's method)
        # ------------------------------------------------------ 
        NI = N - 2
        x = XI[-1]
        dx = 1
        counter = 0
        while (np.abs(dx) > Precision) and (counter < IterMax):
            counter += 1
            FF = YI[NI] + np.tan(thetaI[NI] + muI[NI]) * (x - XI[NI]) - upper_wall(x, upper_wall_type)
            dFFdx = np.tan(thetaI[NI] + muI[NI]) - (
                        upper_wall(x + epsilon, upper_wall_type) - upper_wall(x - epsilon, upper_wall_type)) / (
                                2 * epsilon)
            dx = FF / dFFdx
            x -= dx
        if counter == IterMax: print('BE CAREFUL : non-convergence of  Newton for point P_{N-1}')
        # print("counter (i=N-1) = %i, \t dx = %e"%(counter, dx))
        XP[-1] = x
        YP[-1] = upper_wall(XP[-1], upper_wall_type)
        thetaP[-1] = np.arctan(
            (upper_wall(XP[-1] + epsilon, upper_wall_type) - upper_wall(XP[-1] - epsilon, upper_wall_type)) / (
                        2 * epsilon))

        omegaP[-1] = lambdaPlusI[-1] + thetaP[-1]
        muP[-1] = mu_omega(omegaP[-1])
        lambdaMinusP[-1] = omegaP[-1] + thetaP[-1]
        ax.plot([XP[-1], XI[-1]], [YP[-1], YI[-1]], 'r')  # end of C+

        # -------------------------------------------------------- 
        # filling tables for results 
        # --------------------------------------------------------
        Xmesh.append(deepcopy(XP))
        Ymesh.append(deepcopy(YP))
        mumesh.append(deepcopy(muP))
        X_lower.append(XP[0])
        mulower.append(muP[0])
        X_upper.append(XP[-1])
        muupper.append(muP[-1])

        # end of the loop
        Xmin = min(XP)
        XcP.append(XP[Nmedian])
        YcP.append(YP[Nmedian])
        MachP.append(mu2mach(muP[Nmedian]))

    ax.plot([0, Ltot], [yMedian, yMedian], 'k-.')  # central axis
    if show_Median:
        ax.plot(XcP, YcP, 'ko')
    # I am looking the positions of the points on the central (symmetry) axis   
    # and if they belong to the mesh

    # ****************************************************************
    #  MACH number on the symmetry axis, comparison with 1D theory
    # ****************************************************************
    if show_Mach:

        # solution 1D
        H = upper_wall(0, upper_wall_type) - lower_wall(0, lower_wall_type)
        print("len Xcp = ", len(XcP))
        YcP, Mach_th = [], []  #  possible problem with YcP
        for x in XcP:
            YcP.append((upper_wall(x, upper_wall_type) + lower_wall(x, lower_wall_type)) / 2)
            S = (upper_wall(x, upper_wall_type) - lower_wall(x, lower_wall_type)) / H
            Mach_th.append(inv_S_over_Scrit(S * S_over_Scrit(Mach0, gamma=gam), Mach=3., show=False, gamma=gam))

        fig = plt.figure(1, figsize=(Largeur, Hauteur))
        fig.suptitle("Mach number on the symmetry axis", fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        ax.set_ylabel(r'$x$', fontsize=14)
        ax.set_xlabel(r'$M$', fontsize=14)
        ax.plot(XcP, MachP, 'k-*', linewidth=1, label='Mach on the axis')
        ax.plot(XcP, Mach_th, 'r-', linewidth=1, label='Mach 1D')
        ax.axis()
        ax.set_xlim(0, Ls)

    # ***********************************************
    #   iso Mach lines
    # ***********************************************

    print("len Xmesh : ", len(Xmesh))
    # Y_upper=[upper_wall(x,upper_wall_type) for x in list(X_upper)]
    # Y_lower=[lower_wall(x,lower_wall_type) for x in list(X_lower) ]

    m = len(Xmesh)
    X, Y = np.zeros(([m, N]), dtype=float), np.zeros(([m, N]), dtype=float)
    mac = np.zeros(([m, N]), dtype=float)
    for i, x, y, muj in zip(range(m), Xmesh, Ymesh, mumesh):
        for j, xl, yl, mul in zip(range(N), x, y, muj):
            X[i, j], Y[i, j], mac[i, j] = xl, yl, mu2mach(mul)

    plt.figure(3, figsize=(10, Hequal))
    levels = np.linspace(Mmin, Mmax, nlevels)
    print('Levels = ', levels)
    cp = plt.contour(X, Y, mac, levels, colors='black', linestyles='-')
    # plt.scatter(X, Y,color='black', s=6, marker='o', linewidth=1);
    if show_label:
        plt.clabel(cp, inline=True, fontsize=8)
    plt.title("Iso-Mach lines")
    plt.plot(XL, Y_lower, 'k-', linewidth=1, label='lower wall')
    plt.plot(XL, Y_upper, 'k-', linewidth=1, label='upper wall')
    plt.xlabel("x", fontsize=14)
    plt.ylabel("y", fontsize=14)
    plt.axis('equal')
    plt.xlim(0, Ls)

    # ***********************************************
    #  Iso Mach colored lines ,filled  surfaces
    # ***********************************************

    plt.figure(4, figsize=(12, Hequal))
    plt.title(r'iso-Mach')
    plt.contourf(X, Y, mac, levels)
    plt.axis('equal')
    # plt.scatter(X, Y, c=mac)
    plt.axis('equal')
    plt.xlim(0, Ls)
    plt.colorbar()
    plt.show()
    return


def Exercice12_2():
    """
    Shock tube
    
    Plots only,
    
    Calculations are performed in Fortran with the  EC code.
    """
    set_title("Shock tube solution, plots only ")

    p = np.array([1, 6.3922135814136576, 25.862283608680318])
    t = np.array([0, 9.88327061176361155E-003, 1.78214998216048130E-002]) * 1000

    print("t = ", t)
    print("p = ", p)
    # ratio P2/P1 =P3/P1     =    6.3922135814136576
    # ratio P5/P1 article    =    25.862283608680318
    # ratio P5/P2 article    =    4.0459041737714898   

    # ##################################################
    #   sensor times 
    # ##################################################

    #  axe en t en ms, echelle   =    1.0000000000000000     
    #  Temps A   ( s)           =   1.23540882647045144E-002
    #  Temps 2   ( s)           =   9.88327061176361155E-003
    #  Temps 3   ( s)           =   1.78214998216048130E-002
    #  Temps 2 / Temps A        =   0.80000000000000004     
    #  Temps 3 / Temps A        =    1.4425588873701534     

    fig = plt.figure(1, figsize=(Largeur, Hauteur))
    fig.suptitle('Shock tube', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    # fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$P/P_0$', fontsize=20)
    ax.set_ylabel(r'$t (ms)$', fontsize=20)
    ax.plot(p, t, 's')
    ax.axis([0, 30, 0, 30])
    ax.grid()
    plt.show()

    return


# ************************************************************
#
#  Nozzle of minimal length
#
# ************************************************************

def zone_complexe1(I, theta_I, omega_I, mu_I, CL, gam=gamma):
    """
    Complex zone , from a positive characteristic position towards axis  
    
    Inputs    : values of points I on a  C+
    CL        : theta boundary condition on the axis 
    
    Outputs   : positions of points, omega, theta et mu
    """
    n = len(theta_I)
    x, y = np.zeros((n, n), dtype=float), np.zeros((n, n), dtype=float)
    ome, the = np.zeros((n, n), dtype=float), np.zeros((n, n), dtype=float)
    mut = np.zeros((n, n), dtype=float)
    for l in range(n):
        [x[0, l], y[0, l]] = I[l]
        ome[0, l], the[0, l], mut[0, l] = omega_I[l], theta_I[l], mu_I[l]
        # print(' l = %i, x,y = %f %f , theta = %f '%(l,x[0,l],y[0,l],the[0,l]))
    for k in np.arange(1, n):
        for l in np.arange(k, n):
            if k == l:
                the[k, l] = CL  # axis boundary condition
                ome[k, l] = ome[k - 1, l] + the[k - 1, l] - the[k, l]  # C-
                mut[k, l] = mu_omega(ome[k, l], gam=gam)
                p1 = 0  # axis
                p2 = (the[k - 1, l] - mut[k - 1, l] + the[k, l] - mut[k, l]) / 2  # C-
                [x[k, l], y[k, l]] = line_intersection([x[0, 0], y[0, 0]], [x[k - 1, l], y[k - 1, l]], p1, np.tan(p2))
            else:
                ome[k, l], the[k, l] = complex_line_intersection(the[k - 1, l], ome[k - 1, l], the[k, l - 1],
                                                                 ome[k, l - 1])
                mut[k, l] = mu_omega(ome[k, l], gam=gam)
                p1 = (the[k, l - 1] + mut[k, l - 1] + the[k, l] + mut[k, l]) / 2  # C+
                p2 = (the[k - 1, l] - mut[k - 1, l] + the[k, l] - mut[k, l]) / 2  # C-
                [x[k, l], y[k, l]] = line_intersection([x[k, l - 1], y[k, l - 1]], [x[k - 1, l], y[k - 1, l]],
                                                       np.tan(p1), np.tan(p2))
    return x, y, ome, the, mut


def zone_wall(A0, J):
    """
    Determination of  zone II, wall geometry calculus
    from  C+ coming from  J points
    
    start at the wall at A_0
    
    A = [x, y, omega, theta, mu, M]
    
    J list of type A
    """
    n = len(J)
    print('len J = ', n)
    Pw = [A0]
    for l in range(n):
        Pw.append(deepcopy(J[l]))
        p1 = (Pw[l][3] + J[l][3]) / 2  # mean theta at the wall
        p2 = J[l][3] + J[l][4]
        [Pw[l + 1][0], Pw[l + 1][1]] = line_intersection([Pw[l][0], Pw[l][1]], [J[l][0], J[l][1]], np.tan(p1),
                                                         np.tan(p2))
    print('len Pw = ', len(Pw))
    return Pw


def display_A(A, init=False):
    """
    Characteristic method: 
    display data A vector
    A[[x, y, omega, theta, mu, M],...]
    """
    # print("x = %5.3f, y = %5.3f, omega = %5.3e °, theta = %5.3e °, mu = %5.2f °,M = %3.4f "
    #      %(A[0],A[1],np.rad2deg(A[2]),np.rad2deg(A[3]),np.rad2deg(A[4]),A[5]))
    if init:
        print("  x\t   y  \t   omega(°) theta(°)  mu(°)   M ")
    else:
        print("%7.3f  %7.3f  %7.3f  %7.3f  %5.2f  %3.4f " % (A[0], A[1], np.rad2deg(A[2]),
                                                             np.rad2deg(A[3]), np.rad2deg(A[4]), A[5]))
    return


def distribution(dist_log, eps, N, theta_a):
    """
    Angular distribution (slopes) of the C- from (a) point
    """
    if dist_log:
        # logarithm distribution on the symmetry axis
        theta_log = np.linspace(np.log(eps * theta_a), np.log(theta_a), N)
        theta_l = np.exp(theta_log)
    else:
        # characteristic angles from (a) point
        theta_l = np.zeros((N), dtype=float)
        Nw = N - 3
        theta_l[0], theta_l[1], theta_l[2], theta_l[3] = theta_a / Nw * eps, theta_a / (Nw * 100.), theta_a / (
                    Nw * 10.), theta_a / (Nw * 2.)
        dt = np.linspace(theta_l[2], theta_a, Nw)
        theta_l[4:N] = dt[1:]
    print(len(theta_l), 'theta (°) = ', np.rad2deg(theta_l))
    return theta_l


def TuyereMinimalResults(gam=gamma):
    """
    Nozzle of minimal length
    Results after  simulations with gamma=1.4, 53 characteristics
    Post-treatment
    """
    # #   x            y     M     L approxAiriau 
    # 2.351813  1.176137   1.5    2.432995  
    # 4.830457  1.687388   2.0    4.654693  
    # 6.685007  2.096271   2.25   6.240735  
    # 9.167570  2.636511   2.5    8.332294  
    # 12.487851 3.337441   2.75   11.111386 
    # 16.905614 4.234424   3      14.805188 
    # 30.383045 6.790221   3.5    26.129196 
    # 53.068746 10.721976  4.00   45.399018 
    Mach = [1.5, 2, 2.25, 2.5, 2.75, 3, 3.5, 4]
    L = [2.351813, 4.830457, 6.685007, 9.167570, 12.487851, 16.905614, 30.383045, 53.068746]
    # H:  y coordinates of the last point of the wall. The section is  S = 2 x H
    H = [1.176137, 1.687388, 2.096271, 2.636511, 3.337441, 4.234424, 6.790221, 10.721976]
    Lapprox = [2.432995, 4.654693, 6.240735, 8.332294, 11.111386, 14.805188, 26.129196, 45.399018]
    fig = plt.figure(4, figsize=(10, 8))
    fig.suptitle('Analysys - Nozzle of minimal length', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$M$', fontsize=14)
    ax.set_ylabel(r'$5 H,L$', fontsize=14)
    ax.plot(Mach, L, label="L")
    ax.plot(Mach, Lapprox, label="L appro.")
    ax.plot(Mach, [5 * h for h in H], label="5 x H")
    ax.plot(Mach, [5 * S_over_Scrit(m, gam) for m in Mach], marker="o", markersize=5, label="H 1D")
    plt.legend()
    plt.grid()


def propulsion_coefficient(Pw, Ms, hs, hc=1., pic=1e5, gam=gamma, show=True):
    """
    Propulsion coefficient CF with respect ot the wall points
    Pw[[x, y, omega, theta, mu, M],...]
    """
    x = []
    t = []
    for P in Pw:
        x.append(P[0])
        t.append(np.tan(P[3]) * p_pi(P[5], gam))  # tan theta x p/pi
        if show:
            print(' x = %f, \t M = %f, theta (°) = %f' % (P[0], P[5], np.rad2deg(P[3])))
    tmp = integrate.trapz(t, x=x)
    CF = 2 * tmp / hc - p_pi(Ms, gam) * (hs / hc - 1)
    CFbis = (gam * Ms ** 2 * hs / hc + 1) * p_pi(Ms, gam) - (gam + 1) * p_pi(1, gam)

    print("Ms                     = %f" % (Pw[-1][5]))
    print("CF  (int p-pa)         =%2.4f" % CF)
    print("CF  (momentum)         =%2.4f" % CFbis)
    print("error on CF (x100)     =%2.3f" % (100 * abs(CFbis / CF - 1)))
    print("F                      = %f kN " % (CFbis * pic / 1000 * hc))
    print("hs/hc                  =%f" % (hs / hc))
    print("L                      =%f" % (x[-1]))

    return


def Exercice12_4():
    """
    Nozzle of minimal length when the exit Mach number Ms is known
    by the MoC
    """
    set_title("Nozzle of minimal length")

    print_A = False         # To plot point A-type  data
    print_I = False         # To plot point I-type  data
    print_J = False         # To plot point J-type  data
    display_figure = True   # to plot the detail of constructions
    dist_log = False        # to choose the logarithm point distribution on the symmetry axis

    Ms = 4.5        # Exit Mach number
    N = 53          # number of characteristics originating from  point (a)
    gam = 1.4       # gamma value for this ideal gas
    eps = 1e-3      # small parameter to start calculations close to the throat
    A = []          # A[[x, y, omega, theta, mu, M],...] : list of all the points of the problem
    J = []          # points between  zone I and zone II
    I = []          # first C+ from symmetry axis
    K = []          # points on the main axis
    a = [0., 1.]    # point (a) location
    xCminus = []    # abscissas of the points on C-
    yCminus = []    # ordinates of points on C-
    xCplus = []     # abscissas of points C+
    yCplus = []     # ordinates of points on C+

    # The next three parameters are used to evaluate the thrust but are not used here to determine the shape of the 
    # nozzle 

    h_s = 0         # Section height at exit section
    h_c = 2 * a[1]  # example of throat section (! 2D cartesian theory)
    pi_c = 1e5      # isentropic pressure at the throat (critical condition)
    # S_c  = h_c                # unused
    set_question('1 : Initialization')

    if print_A: display_A(A, init=True)
    # first point (a) :
    theta_a = omega_super(Ms, gamma=gam) / 2.  # initial value at point (a) :
    omega_a = theta_a
    mu_a = mu_omega(omega_a, gam=gam)
    Mach_a = mu2mach(mu_a)
    A.append([a[0], a[1], omega_a, theta_a, mu_a, Mach_a])
    if print_A: display_A(A[0])
    print('theta max (°)        = %f' % (np.rad2deg(theta_a)))

    # C- from point (a)
    theta_l = distribution(dist_log, eps, N, theta_a)
    # omega and mu, lambda Minus = associated C Minus  :
    omega_l = theta_l
    mu_l = []
    for omega in omega_l:
        mu_l.append(mu_omega(omega, gam=gam))
    Cminus_l = omega_l + theta_l
    print('Cminus_l = ', np.rad2deg(Cminus_l))

    # *******************************************************************
    # resolution of the first C+ ( I_1 I_N), solve points I and J
    # *******************************************************************
    omega_I, theta_I, mu_I, Mach_I = [], [], [], []

    set_question('2 : calculus of point I_1')

    # just coordinates are found in I 
    # A[1:N] will get all points I with all quantities
    I.append(line_intersection(a, [0., 0.], np.tan(theta_l[0] - mu_l[0]), 0.))
    omega_I.append(Cminus_l[0])
    theta_I.append(0.)
    mu_I.append(mu_omega(omega_I[0], gam=gam))
    Mach_I.append(mu2mach(mu_I[0]))
    A.append([I[0][0], I[0][1], omega_I[0], theta_I[0], mu_I[0], Mach_I[0]])
    K.append(A[1])
    print(K)

    if print_I:
        print("I : omega = %5.3e °  \t theta = %5.3e ° \t mu = %5.2f ° \t  M = %3.4f " % (
        np.rad2deg(omega_I[0]), np.rad2deg(theta_I[0]), np.rad2deg(mu_I[0]), Mach_I[0]))
    if print_A: display_A(A[1])

    # all the C- start in (a)
    xCminus.append([a[0], I[0][0]])
    yCminus.append([a[1], I[0][1]])

    # plot initialization
    plt.figure(0, figsize=(10, Hequal))
    plt.title('Nozzle of minimal length', fontsize=14, fontweight='bold')
    plt.ylabel(r'$y$', fontsize=14)
    plt.xlabel(r'$x$', fontsize=14)
    plt.plot([a[0], I[0][0]], [a[1], I[0][1]], 'k-')
    plt.axis("equal")

    set_question('3 : points from I_2 to I_N')

    for l in np.arange(1, N):
        ome, the = complex_line_intersection(theta_l[l], omega_l[l], theta_I[l - 1], omega_I[l - 1])
        mue = mu_omega(ome, gam=gam)
        p1 = (the + theta_l[l] - mue - mu_l[l]) / 2             # C- slope from (a)
        p2 = (the + theta_I[l - 1] + mu_I[l - 1] + mue) / 2     # C+ slope from I(l-1)
        I.append(line_intersection(a, I[l - 1], np.tan(p1), np.tan(p2)))
        omega_I.append(ome)
        theta_I.append(the)
        mu_I.append(mue)
        Mach_I.append(mu2mach(mu_I[l]))
        A.append([I[l][0], I[l][1], omega_I[l], theta_I[l], mu_I[l], Mach_I[l]])
        if print_I:
            print("I : omega = %5.3e °  \t theta = %5.3e ° \t mu = %5.2f ° \t  M = %3.4f " % (
            np.rad2deg(omega_I[l]), np.rad2deg(theta_I[l]), np.rad2deg(mu_I[l]), Mach_I[l]))
        if print_A: display_A(A[l + 1])
        plt.plot([I[l - 1][0], I[l][0]], [I[l - 1][1], I[l][1]], 'r--')     # C+
        plt.plot([a[0], I[l][0]], [a[1], I[l][1]], 'k-')                    # C-
    # print('len A = %i, \t N =%i '%(len(A),N))

    # *************************************************
    # resolution of characteristics under the first C+ :
    # *************************************************

    set_question('4 : Zone under  C+ (I_1,I_N), points J')

    x, y, ome, the, mut = zone_complexe1(I, theta_I, omega_I, mu_I, CL=0, gam=gam)
    J.append(A[-1])  # points bounded the first zone  J[0]=I[N-1]

    # to get  C- and C+
    for j in np.arange(1, N):
        xtmp, ytmp = [], []
        xtmp.append(a[0])
        ytmp.append(a[1])
        for xl, yl in zip(x[:j + 1, j], y[:j + 1, j]):
            xtmp.append(xl)
            ytmp.append(yl)
        xCminus.append(xtmp)
        yCminus.append(ytmp)

    for k in np.arange(1, N):
        xp, yp = [], []
        plt.plot(x[k, k], y[k, k], 'ko')
        plt.plot(x[k, N - 1], y[k, N - 1], 'bs')
        J.append([x[k, N - 1], y[k, N - 1], ome[k, N - 1], the[k, N - 1], mut[k, N - 1], mu2mach(mut[k, N - 1])])
        for l in np.arange(k, N):
            plt.plot([x[k - 1, l], x[k, l]], [y[k - 1, l], y[k, l]], 'k-')
            plt.plot([x[k, l - 1], x[k, l]], [y[k, l - 1], y[k, l]], 'k--')
            A.append([x[k, l], y[k, l], ome[k, l], the[k, l], mut[k, l], mu2mach(mut[k, l])])
            if the[k, l] == 0:
                K.append(A[-1])
            xp.append(x[k, l])
            yp.append(y[k, l])
        xCplus.append(xp)
        yCplus.append(yp)

    for l in range(N):
        plt.plot(x[0, l], y[0, l], 'ro')

    if print_J:
        display_A(A[0], init=True)
        print('J values : ')
        for l in range(N):
            display_A(J[l])

    # *************************************************
    # resolution of zone 2 and of the wall :
    # *************************************************
    # wall point Pw, the last points of  A list.
    Pw = zone_wall(A[0], J)  # to calculate the  positions of points Pw, the other values are  points J
    for P in Pw:
        A.append(P)
    print('Pw (wall) : ')
    display_A(A[0], init=True)
    for l in range(len(Pw)):
        display_A(Pw[l])
    for l in range(N):
        plt.plot([Pw[l][0], Pw[l + 1][0]], [Pw[l][1], Pw[l + 1][1]], 'k-')
    for l in range(N):
        plt.plot([J[l][0], Pw[l + 1][0]], [J[l][1], Pw[l + 1][1]], 'r--')

    # *****************************************************************
    # figure 2 initialization, final figure for C+ and C-
    # *****************************************************************

    plt.figure(1, figsize=(10, Hequal))
    plt.title('Nozzle of minimal length', fontsize=14, fontweight='bold')
    plt.ylabel(r'$y$', fontsize=14)
    plt.xlabel(r'$x$', fontsize=14)
    # m = len(A)

    for l in range(N):
        plt.plot([Pw[l][0], Pw[l + 1][0]], [Pw[l][1], Pw[l + 1][1]], 'k-')  # wall
    xi, yi = [], []
    # plot of  C+
    # first C+ from I points  
    for l in np.arange(1, N + 1):
        xi.append(A[l][0])
        yi.append(A[l][1])
    xi.append(Pw[1][0])
    yi.append(Pw[1][1])
    plt.plot(xi, yi, 'b-')
    # the other C+
    for l in range(N - 1):
        xi, yi = [], []
        for x, y in zip(xCplus[l], yCplus[l]):
            xi.append(x)
            yi.append(y)
        xi.append(Pw[l + 2][0])
        yi.append(Pw[l + 2][1])
        plt.plot(xi, yi, 'b-')
    # the last C+ coming from the symmetry axe and delimiting zone III. 
    plt.plot([J[-1][0], Pw[-1][0]], [J[-1][1], Pw[-1][1]], 'b-')
    # C-
    for l in range(N):
        xi, yi = [], []
        for x, y in zip(xCminus[l], yCminus[l]):
            xi.append(x)
            yi.append(y)
        plt.plot(xi, yi, 'r-')
    plt.axis("equal")

    L, R_s = A[-1][0], A[-1][1]
    h_s = 2 * R_s
    Lapprox = (h_c + h_s) / 2 * np.sqrt(Ms ** 2 - 1)
    print('Nozzle length                         : %f' % L)
    print('Lenght, approximation                 : %f\n' % Lapprox)
    print('1/2 section heigth at exit            : %f' % R_s)
    print('Section ratio  hs/hc                  : %f' % (h_s / h_c))
    print('hs/hc, 1D theory                      : %f' % (S_over_Scrit(Ms)))
    print("difference in percents on hs/hc       : %f" % (100 * (h_s / h_c / S_over_Scrit(Ms) - 1)))
    print('pa/pi_c                               : %f' % (p_pi(Ms, gam)))

    # *****************************************************************
    # figure 3, Mach number
    # *****************************************************************

    plt.figure(3, figsize=(10, 8))
    plt.title('Nozzle of minimal length', fontsize=14, fontweight='bold')
    plt.ylabel(r'$M$', fontsize=14)
    plt.xlabel(r'$x$', fontsize=14)
    plt.plot([Pw[l][0] for l in range(N + 1)], [Pw[l][5] for l in range(N + 1)], 'ko-', label="wall")  # wall
    plt.plot([I[l][0] for l in range(N)], [Mach_I[l] for l in range(N)], label="I")  # points I
    plt.plot([J[l][0] for l in range(N)], [J[l][5] for l in range(N)], label="J")  # points J
    plt.plot([K[l][0] for l in range(N)], [K[l][5] for l in range(N)], label="K")  # points J

    plt.legend()
    plt.grid()

    propulsion_coefficient(Pw, Ms, h_s, hc=h_c, pic=pi_c, gam=gam, show=False)

    TuyereMinimalResults(gam)

    if display_figure:
        plt.show()
    return
