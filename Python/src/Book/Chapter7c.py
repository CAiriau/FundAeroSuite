#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 7 exercises (cont')
    ..
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from IncompressibleFlow.PanelMethod import Van_Hooren_gmtry, grid_definition, Panel, Freestream, arctan3, \
    velocity_field, velocity_field_vector
from IncompressibleFlow.PanelMethod import plot_profile_directions, plot_profile_with_panels, \
    plot_doublet_derivative_panels, plot_source_panels, aerodynamics, plot_doublet_panels

from IncompressibleFlow.Airfoil import naca
from scipy.interpolate import interp1d
from Tools.misc import set_title, set_question


def distance_rotation(panel_k, panel_j):
    """
    to calculate influence coefficients
    convert collocation points to local panel coordinates.
    k : control point, j : singularity panel
    """
    xt, yt = panel_k.xc - panel_j.xa, panel_k.yc - panel_j.ya
    x2t, y2t = panel_j.xb - panel_j.xa, panel_j.yb - panel_j.ya
    Cos, Sin = np.cos(panel_j.theta), np.sin(panel_j.theta)
    # (x_p,y_p) -> (x,y)
    # x_k,y_k    =  xt*Cos-yt *Sin, xt*Sin+yt*Cos
    # x_2,y_2    = x2t*Cos-y2t*Sin, 0.
    # (x,y) -> (x_p,y_p)
    x_k, y_k = xt * Cos + yt * Sin, -xt * Sin + yt * Cos
    # x_2,y_2    = x2t*Cos+y2t*Sin, -x2t*Sin+y2t*Cos; print('y2 = ',y_2)  # vérfié : y_2 = 0
    x_2, y_2 = x2t * Cos + y2t * Sin, 0
    # x_2  = simply panel length
    # print(x_k,y_k,x_2)

    r_k, r_2 = np.sqrt(x_k ** 2 + y_k ** 2), np.sqrt((x_k - x_2) ** 2 + y_k ** 2)
    theta_k, theta_2 = np.arctan2(y_k, x_k), np.arctan2(y_k, x_k - x_2)

    return x_k, y_k, x_2, y_2, r_k, r_2, theta_k, theta_2


def Exercice7_5():
    """ 
    Uniform distribution of  sources and of  doublets for the  Van-Hooren airfoil
    """
    # very difficult problem, depending on the numerous parameters and options. One of the problem is the arctan 
    # functions - sometimes atan or sometimes atan2 -  which create phase discontinuities des discontinuités on the 
    # angle We also need to take care of the rotation direction (anti-clockwise (here) anti-clockwise (usual 
    # english/american publication)) when panels are described around the  airfoil We have to take care when changing
    # the parameters. When it works, save your parameters... 

    set_title("Uniform distribution of  sources and of  doublets for the  Van-Hooren airfoil")

    # data, parameters

    U0 = np.float64(1)          # reference velocity
    chord = np.float64(1.)      # chord length
    rho0 = 1.4                  # air density (not used here)
    aoa = np.array([5.000])     # angle of attack = 0 here, thickness effect
    t = np.float64(0.055)       # thickness parameter
    k = 1.9                     # trailing edge angle parameter
    npt = 91                    # number of points on the airfoil and of control points
    number = '4412'             # '4412' Naca 4 digit number, in case the option is selected
    Alpha = np.deg2rad(aoa[0])

    def_profile = 0             # 0: Van Vooren, 1 : NACA
    stop = 0                    # stop  at a given point in the program
    plot_option = False         # to have some intermediate plot  set True
    plot_streamlines = False    # to implement
    plot_vectors = False        # to implement  to plot the slip velocity on the airfoil surface
    plot_Kp = True              # plot the pressure coefficient
    plot_profile = False        # geometry (only for Van de Vooren airfoil)
    plot_panels = False         # to plot the panels and their center
    plot_sources = True         # to plot the source distribution
    plot_doublet = True         # to plot the doublet distribution

    N_interior = False          # True: the normal direction remains interior to the airfoil
    Ut_from_dmudx = False       # True: Ut from velocity, False: from potential line
    #                               otherwise  Ut is calculated  from the  potential phi
    algo_doublet = True         # Personal programming for the matrices

    # option to print data on screen :
    print_profile = False
    print_sigma = False
    print_RHS = False
    print_A = False
    print_mu = False
    print_Kp = False

    signe = -1  # related to the rotation direction used (trigonometry or clockwise) along the airfoil
    delta_is_theta = False  # delta (panel slop) no longer in [-90° +90°] but in [-180° +180°], and  equal to theta
    # 
    if N_interior is False:
        print("N interior is False")

    # Geometry
    if def_profile == 0:
        x, y, theta = Van_Hooren_gmtry(chord, t, k, npt, plot_profile, origin=True, dtheta=0.0)
    else:
        x, y, theta = naca(number, (npt - 1) / 2, half_cosine_spacing=True)

    n = int((npt - 1) / 2)
    # to verify the location of the leading and trailing edge
    # this 4 next lines can be deleted.
    xe, ye = x[0:n + 1], y[0:n + 1]
    print("length of xe = %i, n= %i " % (len(xe), int(n)))
    print(" [xe_LE, xe_TE], [ye_LE,ye_TE] ", xe[list((0, -1))], ye[list((0, -1))])

    # Stream flow conditions
    fst = Freestream(U_inf=U0, alpha=aoa, rho=rho0)
    print("Panel definition")
    N = len(x) - 1
    print('Number of panels on the airfoil = %i' % (N))
    # Ns, Nd = N, N+1                         # number of sources, number of doublets
    panels = np.zeros(N, dtype=object)
    xc = np.zeros(N, dtype=float)
    for i in range(N):
        panels[i] = Panel(x[i], y[i], x[i + 1], y[i + 1], n_interior=N_interior, delta_is_theta=delta_is_theta)
        xc[i] = panels[i].xc  # panel center x-coordinate

    plot_profile_directions(panels)
    # print('xc',xc)
    # print('theta = ',theta)

    # display : airfoil, collocation points, angles - slope 
    if print_profile:
        print("{%s} \t {%s} \t \t  {%s} \t \t {%s} \t \t {%s} \t \t {%s} \t \t {%s} " % (
            'theta (°)', 'x', 'y', 'i', 'xc ', 'delta(°)', 'theta(°)'))
        for i in range(N):
            print('{%06.2f} \t {%+06.4f} \t {%+06.4f} \t {%4i} \t {%+06.4f} \t {%+06.4f} \t {%+06.4f}' % (
                np.rad2deg(theta[i]),
                x[i], y[i], i, panels[i].xc, np.rad2deg(panels[i].delta), np.rad2deg(panels[i].theta)))

    if plot_panels:
        plot_profile_with_panels(panels)

        # calculation of the influence coefficients

    sigma = np.zeros(N, dtype=float)

    # source contribution on a panel from others
    # j : source or  doublet panels
    # k : panel with the collocation point where the influence is calculated
    # page 151 of the book "Aérodynamique Instationnaire", section 7.3.2.2
    # Aprime : source, A : doublet
    # coef = 0.5/np.pi
    for j in range(N):
        #  panels[j].sigma = (fst.u_inf*panels[j].nx+fst.v_inf*panels[j].ny)  # *panels[i].epsilon
        panels[j].sigma = np.sin(panels[j].theta - fst.alpha)
        sigma[j] = panels[j].sigma

        # print('Source vector : S =',sigma)
    print("sum_k  S_k           = %f" % (np.sum(sigma)))
    somme = 0.
    for i in range(N):
        somme += panels[i].sigma * panels[i].length
    print("integral of S_k      = %f" % (somme))

    dl = np.zeros((N), dtype=float)


    def calcul_distribution_doublet1():
        """
        Method found in the book
        """
        print('Book algorithm')
        A, Aprime = np.zeros((N + 1, N + 1), dtype=float), np.zeros((N, N), dtype=float)
        RHS = np.zeros(N, dtype=float)

        for k in range(N):
            for j in range(N):
                x_k, y_k, x_2, y_2, r_k, r_2, theta_k, theta_2 = distance_rotation(panels[k], panels[j])
                if k == 1:
                    dl[j] = x_2
                if k == j:
                    # The next two followings expressions are equivalent
                    # A[k,j],Aprime[k,j]=signe*0.5,0.5/np.pi*panels[k].length*np.log(panels[k].length/2.0)
                    A[k, j], Aprime[k, j] = signe * 0.5, 1. / np.pi * x_k * np.log(r_k)

                else:
                    A[k, j] = signe * 0.5 / np.pi * (theta_2 - theta_k)
                    Aprime[k, j] = 1. / (2 * np.pi) * (
                            x_k * np.log(r_k) - (x_k - x_2) * np.log(r_2) + y_k * (theta_2 - theta_k))

        A[:, N] = 0
        # Kutta's trailing edge condition
        A[N, N] = 1  # mu_w
        A[N, 0], A[N, N - 1] = 1, -1  # mu_1 -mu_N  
        # wake influence :

        for k in range(N):
            # theta_k=np.arctan2(panels[k].yc-panels[N-1].yb,panels[k].xc-panels[N-1].xb)
            theta_k = np.arctan((panels[k].yc - panels[N - 1].yb) / (panels[k].xc - panels[N - 1].xb))
            A[k, N] = -0.5 / np.pi * (0.0 - theta_k)

        RHS = -1 * np.dot(Aprime, sigma)
        RHS = np.append(RHS, 0)  # last row (or component)
        return A, RHS, dl


    def calcul_distribution_doublet_Katz():
        """
        This code is a translation in python of the original fortran77 code provided by Katz and Plotkins
        """
        print("Algorithm found in Katz' book")
        RHS = np.zeros(N, dtype=float);
        temp = np.zeros(N, dtype=float)
        A = np.zeros((N + 1, N + 1), dtype=float)
        for k in range(N):
            temp[:] = 0.0
            for j in range(N):
                x_k, y_k, x_2, y_2, r_k, r_2, theta_k, theta_2 = distance_rotation(panels[k], panels[j])
                if k == 1:
                    dl[j] = x_2
                if k == j:
                    A[k, j] = signe * 0.5
                    temp[j] = panels[j].sigma / np.pi * x_k * np.log(r_k)
                    # print(x_k,r_k,1/np.pi*x_k*np.log(r_k),k)
                else:
                    A[k, j] = -0.5 / np.pi * (theta_2 - theta_k)
                    temp[j] = panels[j].sigma / (2 * np.pi) * (
                            x_k * np.log(r_k) - (x_k - x_2) * np.log(r_2) + y_k * (theta_2 - theta_k))
                    #   add wake influence
            xw_tmp = panels[k].xc - panels[N - 1].xb
            yw_tmp = panels[k].yc - panels[N - 1].yb
            dthw = -np.arctan(yw_tmp / xw_tmp)
            A[k, N] = -1 / (2 * np.pi) * dthw
            RHS[k] = -sum(temp)

            #     add an explicit Kutta condition
        A[N, :] = 0
        A[N, 0] = 1
        A[N, N - 1] = -1
        A[N, N] = 1
        RHS = np.append(RHS, 0)
        return A, RHS, dl

    # The both approach must give the same results ! it is a validation 

    if stop == 1:
        return

    if algo_doublet:
        A, RHS, dl = calcul_distribution_doublet1()
    else:
        A, RHS, dl = calcul_distribution_doublet_Katz()

    if print_sigma:
        print("sigma = ", sigma)
    if print_RHS:
        print('RHS = ', RHS)
    if print_A:
        print('A = ', A)

    print('Direct inversion')
    mu = np.linalg.solve(A, RHS)
    if print_mu:
        print('mu = ', mu)
    print('len(mu) =  %i, \t N = %i' % (len(mu), N))
    print('sum(mu[:]) = %f, sum(mu[1:-2] = %f' % (np.sum(mu), np.sum(mu[0: -1])))
    if stop == 2:
        return

    for i in range(N):
        panels[i].mu = mu[i]
    print('mu_w = %f' % (mu[len(mu) - 1]))
    print('mu_1 = %f, mu_N = %f, mu_w = %f' % (mu[0], mu[N - 1], mu[N]))

    Cl = 2 * mu[N] / (chord * U0)
    print('Cl = %f, 2 pi alpha = %f' % (Cl, 2 * np.pi * Alpha))
    if Alpha == 0:
        print('estimation of alpha_0 = %f °' % (np.rad2deg(-Cl / (2 * np.pi))))

    if plot_doublet:
        plot_doublet_panels(panels)
    if plot_sources:
        plot_source_panels(panels)

    # tangential velocity  : 2  possible approaches
    x_kp = np.zeros(N, dtype=float)
    Ut = np.zeros(N, dtype=float)
    print(N, len(panels), len(mu))

    ds = []
    for i in range(N - 1):
        ds.append(-(panels[i].length + panels[i + 1].length) / 2 * panels[i].epsilon)
        x_kp[i] = panels[i].xb
        panels[i].dmudx = (mu[i + 1] - mu[i]) / ds[i]
        #
    # print('ds = ',ds[0:4],ds[-4:-1])
    # 
    # Method 2, by derivative of the normal doublet
    #
    if Ut_from_dmudx:
        print('tangential velocity  with d mu / dx')
        for i in range(N - 1):
            delta = (panels[i].delta + panels[i + 1].delta) / 2
            Ut[i] = U0 * np.cos(delta - fst.alpha) - panels[i].dmudx
    else:
        # 
        # Method 1, with the velocity potential
        #
        print('tangential velocity  d phi / dx')
        phi = np.zeros(N, dtype=float)
        for i in range(N):
            phi[i] = U0 * panels[i].xc * np.cos(fst.alpha) + U0 * panels[i].yc * np.sin(fst.alpha)
        phi -= mu[:len(mu) - 1]
        dr = []
        for i in range(N - 1):
            dr.append(-(dl[i + 1] + dl[i]) / 2 * panels[i].epsilon)
            Ut[i] = (phi[i + 1] - phi[i]) / dr[i]
        #  print('dr = ',dr[0:4],dr[-4:-1])     

    # in fact  dr=ds: we could simplify the coding

    if plot_doublet:
        plot_doublet_derivative_panels(panels)

    Kp = 1 - Ut ** 2 / U0 ** 2  # pressure coefficient 
    print(" N                     : %i" % (N))
    print("Upper wall point index : %i \t %i" % (0, n - 1))
    print("Lower wall point index : %i \t %i" % (n, N - 1))
    upper_wall = np.arange(0, n - 1)
    lower_wall = np.arange(n, N - 1)

    xe_kp = np.flipud(np.concatenate(([1], x_kp[upper_wall], [0])))
    Kpe = np.flipud(np.concatenate(([1], Kp[upper_wall], [0])))
    Ute = np.flipud(np.concatenate(([0], Ut[upper_wall], [0])))
    xi_kp = np.concatenate(([0], x_kp[lower_wall], [1]))
    Kpi = np.concatenate(([0], Kp[lower_wall], [1]))
    Uti = np.concatenate(([0], Ut[lower_wall], [0]))

    def plot_tangential_velocity(xe, xi, Ute, Uti):
        """
        plot the tangentiel velocity 
        """
        width = 10
        ech = 0.05
        plt.figure(figsize=(width, 0.8 * width))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('Vt', fontsize=16)
        plt.title('Tangential velocity')
        plt.plot(xe, Ute, linestyle='-', linewidth=1, marker='o', markersize=8, color='red', label='upper_wall')
        plt.plot(xi, Uti, linestyle='-', linewidth=1, marker='o', markersize=8, color='blue', label='lower_wall')
        # plt.plot(np.append(0,x_kp[lower_wall]),np.append(0,Kp[lower_wall]),linestyle='-', linewidth=1, marker='o',
        # markersize=8, color='blue',label='lower_wall')
        plt.legend(loc='lower right')
        plt.grid()
        plt.xlim(0 - ech, chord + ech)

    plot_tangential_velocity(xe_kp, xi_kp, Ute, Uti)
    aerodynamics(xe_kp, xi_kp, Kpe, Kpi)

    if print_Kp:
        print('Upper wall')
        for i in upper_wall:
            print(x_kp[i], Kp[i], panels[i].dmudx)
        print('Lower wall')
        for i in lower_wall:
            print(x_kp[i], Kp[i], panels[i].dmudx)
    print('x at N-1 = %f \t , Upper wall Kp = %f, \t Lower wall Kp  %f' % (x_kp[0], Kp[0], Kp[N - 2]))
    print('Upper wall point number = %i, Lower wall point number = %i' % (len(upper_wall), len(lower_wall)))
    if plot_Kp:
        # load geometry from data file
        reference_filepath = os.path.join('Book/Data', 'Kp_VanVooren_5_deg.dat')
        with open(reference_filepath, 'r') as infile:
            xKatz, CpKatz = np.loadtxt(infile, dtype=float, unpack=True, skiprows=0)
        reference_filepath = os.path.join('Book/Data', 'Constant_Source_Doublet_Potential.dat')
        with open(reference_filepath, 'r') as infile:
            xKatz_num, CpKatz_num, ind = np.loadtxt(infile, dtype=float, unpack=True, skiprows=2)
        reference_filepath = os.path.join('Book/Data', 'Kp_VanVooren_0_deg.dat')
        with open(reference_filepath, 'r') as infile:
            xKatz0, CpKatz0 = np.loadtxt(infile, dtype=float, unpack=True, skiprows=0)

        width = 10
        ech = 0.05
        plt.figure(figsize=(width, 0.8 * width))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('Kp', fontsize=16)
        plt.title('Pressure coefficient')
        plt.plot(np.concatenate(([1], x_kp[upper_wall], [0])), np.concatenate(([1], Kp[upper_wall], [0])),
                 linestyle='-',
                 linewidth=1, marker='o', markersize=8, color='red', label='Upper wall')
        plt.plot(np.append(0, x_kp[lower_wall]), np.append(0, Kp[lower_wall]), linestyle='-', linewidth=1,
                 marker='o', markersize=8, color='blue', label='Lower wall')
        plt.plot(xKatz, CpKatz, linestyle='', linewidth=1, marker='o', markersize=4, markevery=2,
                 color='red', label='ici Th.')
        plt.plot((xKatz_num + 1) / 2, CpKatz_num, linestyle='--', linewidth=3, marker='', markersize=4, markevery=2,
                 color='green', label='Katz Num.')
        plt.plot(xKatz0, CpKatz0, linestyle='', linewidth=1, marker='o', markersize=4, markevery=2, color='black',
                 label='ici Th. alpha=0°')
        plt.legend(loc='lower right')
        plt.grid()
        plt.xlim(0 - ech, chord + ech)
        # plt.ylim(-2.1, 1.1);

    # if plot_streamlines:
    #     X,Y=grid_definition(50,-0.2,1.2,-0.25,0.25)
    #     velocity_field(fst.u_inf,fst.v_inf,X,Y,xs,S,plot_streamlines,x,y)

    # if plot_vectors:    
    #     # Points de contrôle :
    #     xc,yc = np.zeros((N), dtype=float),np.zeros((N), dtype=float)
    #     for i, panel_i in enumerate(panels):
    #         xc[i],yc[i]=panel_i.xc,panel_i.yc 
    #     velocity_field_vector(fst.u_inf,fst.v_inf,xc,yc,xs,S,plot_vectors,x,y)

    if plot_option or plot_streamlines or plot_vectors or plot_Kp or plot_profile:
        plt.show()
    elif plot_panels or plot_sources or plot_doublet:
        plt.show()

    print('Normal end of execution')
