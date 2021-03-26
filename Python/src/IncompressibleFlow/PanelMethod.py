#!/bin/python

"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Module : incompressibleFlow.PanelMethod
    ..

Python module providing functions for panel method
    * Freestream class
    * Panel class
    * Van Hooren airfoil
    * various functions related to the panel method and plots

"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


class Freestream:
    """
    Freestream conditions.
    """

    def __init__(self, U_inf=1.0, alpha=0.0, rho=1.0):
        """
        Sets the freestream speed and angle (in degrees).
        
        Parameters
        ----------
        U_inf: float, optional
            Freestream speed;
            default: 1.0.
        alpha: float, optional
            Angle of attack in degrees;
            default 0.0.
        rho : density
        """
        self.U_inf = U_inf
        self.alpha = alpha * np.pi / 180.0  # degrees to radians
        self.rho = rho
        self.q = 1 / 2 * rho * U_inf ** 2
        # velocity components
        self.u_inf, self.v_inf = U_inf * np.cos(self.alpha), U_inf * np.sin(self.alpha)


class Panel:
    """
    Contains information related to a panel.
    """
    def __init__(self, xa, ya, xb, yb, coef=0.5, n_interior=True, delta_is_theta=False):
        """
        Initializes the panel.

        * The panel is defined by two points A and B.
        * Sets the end-points and calculates the center-point, length and angle (with the x-axis) of the panel.
        * Defines if the panel is located on the upper or lower surface of the geometry.
        * Initializes the source-strength, tangential velocity, and pressure coefficient of the panel to zero.
        
        Parameters
        ---------_
        xa: float,
            x-coordinate of the first end-point.
        ya: float,
            y-coordinate of the first end-point.
        xb: float,
            x-coordinate of the second end-point.
        yb: float,
            y-coordinate of the second end-point.
        delta: float
            slope of the panel ( pi/2< delta < pi/2)
        theta: float
            angle of the panel ( pi  < delta < pi  )
        """
        self.xa, self.ya = xa, ya  # panel starting-point
        self.xb, self.yb = xb, yb  # panel ending-point
        self.xc, self.yc = xa + (xb - xa) * coef, ya + (yb - ya) * coef  # panel center
        self.length = np.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)  # panel length

        # orientation of panel (angle between x-axis and panel's tangent)
        self.theta = np.arctan2((yb - ya), (xb - xa))

        if xb - xa <= 0.0:
            self.delta = np.arcsin(-(yb - ya) / self.length)
            # self.xa,self.ya=xb,yb
            # self.xb,self.yb=xa,ya
        elif xb - xa > 0.0:
            self.delta = np.arcsin((yb - ya) / self.length)

        if delta_is_theta:
            self.delta = self.theta

        # panel location
        self.epsilon = 0
        if xb - xa >= 0:
            self.loc = 'lower_wall'  # upper surface
            self.epsilon = -1
        else:
            self.loc = 'upper_wall'  # lower surface
            self.epsilon = 1

        self.sigma = 0.0    # source strength
        self.mu = 0.0       # doublet intensity
        self.dmudx = 0.0    # doublet intensity
        self.Ut = 0.0       # tangential velocity
        self.Kp = 0.0       # pressure coefficient

        # tangential and normal directions
        if n_interior:
            self.tx, self.ty = np.cos(self.theta), np.sin(self.theta)
            self.nx, self.ny = -np.sin(self.theta), np.cos(self.theta)
        else:
            # normal exterior, t always downstream direction
            self.tx, self.ty = np.cos(self.delta), np.sin(self.delta)
            self.nx, self.ny = -self.ty * self.epsilon, self.tx * self.epsilon


def arctan3(sint, cost):
    """ 
    Angle in radians from a continuous arc-tangent (to verify)
    """
    if sint == 0:
        angle = np.pi / 2
    else:
        angle = np.arctan(sint / cost)
    if cost < 0 and sint > 0:
        theta = angle + np.pi
    elif cost < 0 and sint < 0:
        theta = angle + np.pi
    elif cost > 0 and sint < 0:
        theta = angle + 2 * np.pi
    else:
        theta = angle
    return theta


def Van_Hooren_gmtry(chord, epsilon, k, npt, plot_option, origin=False, dtheta=0):
    """
    Van-Hooren Airfoil geometry

    Args:
        chord (float) : chord
        npt (int) : numbers of points 16, 31, 61, 91, 181, 361, etc ...
        k (float) : airfoil parameter
        epsilon (float) : airfoil thickness parameter
        origin (bool):  True: first table value at the trailing edge, else at the leading edge
        dtheta (float) : angle in degrees to avoid calculation just at the leading edge

    Returns:
        list of float : x, y, theta=arctan(y/x)
    """
    eps = dtheta / 180.0 * np.pi                            # to avoid 0 division at theta=0 or 180 °
    theta = np.linspace(-np.pi + eps, np.pi - eps, npt)     # theta=0 => x= 0 : Leading Edge
    if origin:
        theta = theta + np.pi                               # theta=0 => x= 1 : Trailing Edge

    # eta=1-k+epsilon*k                                     # parameter for perturbed velocity and Kp (unused)
    flag_arctan = 3                                         # 3 types of arctan functions
    angle_flag = False                                      # to plot the angles instead of the airfoil geometry
    ell = chord
    a = chord / 2.0 ** k * (1.0 + epsilon) ** (k - 1)       # parameter of the profil chord -> a
    print('k = %f, epsilon = %f, a = %f ' % (k, epsilon, a))
    r1 = np.sqrt(np.sin(theta) ** 2 + (np.cos(theta) - 1) ** 2)
    r2 = np.sqrt(np.sin(theta) ** 2 + (np.cos(theta) - epsilon) ** 2)
    n = len(theta)
    # print('n = %d '%n)
    t1 = np.arctan2(np.sin(theta), (np.cos(theta) - 1), )
    t2 = np.arctan2(np.sin(theta), (np.cos(theta) - epsilon))
    theta1, theta2 = np.zeros([n]), np.zeros([n])
    for i in range(n):
        sint = np.sin(theta[i])
        cost = np.cos(theta[i])
        if theta[i] == 0:
            theta1[i] = np.pi / 2
        else:
            if np.abs(cost - 1) <= 1e-13:
                # print('error with  theta1 or i= %i, theta= %f'%(i,theta[i]))
                theta1[i] = 0
            else:
                theta1[i] = np.arctan(sint / (cost - 1)) + np.pi

        if flag_arctan == 1:
            if cost - epsilon < 0 and sint > 0:
                theta2[i] = np.arctan(sint / (cost - epsilon)) + np.pi
            elif cost - epsilon < 0 and sint < 0:
                theta2[i] = np.arctan(sint / (cost - epsilon)) + np.pi
            elif cost - epsilon > 0 and sint < 0:
                theta2[i] = np.arctan(sint / (cost - epsilon)) + 2 * np.pi
            else:
                theta2[i] = np.arctan(sint / (cost - epsilon))
        elif flag_arctan == 2:
            theta2[i] = arctan3(sint, cost - epsilon)

    if flag_arctan == 3:
        theta2 = t2
        theta1 = t1

    phi = k * theta1 - (k - 1) * theta2
    r = a * r1 ** k / r2 ** (k - 1)
    x = r * np.cos(phi) + ell;
    y = r * np.sin(phi)

    if plot_option:
        print('Van de Vooren airfoil')
        nfig = 1
        fig = plt.figure(nfig)
        fig.suptitle('Airfoil', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$x$', fontsize=20)
        ax.set_ylabel(r'$y$', fontsize=20)
        ax.grid()
        # ax.legend(loc='upper left')
        if angle_flag:
            ax.plot(theta, theta1, '-')
            ax.plot(theta, theta2, '--')
            ax.plot(theta, t1, 'o')
            ax.plot(theta, t2, 's')
        else:
            ax.plot(x, y, '-', linewidth=2)
            ax.axis('equal')
            # ax.plot(theta,r)
        plt.show()

    if origin is False:
        return np.flipud(x), np.flipud(y), np.flipud(theta)
    else:
        return x, y, theta


def grid_definition(N, xmin, xmax, ymin, ymax):
    """
    Grid with constant step size
    """
    X, Y = np.meshgrid(np.linspace(xmin, xmax, N), np.linspace(ymin, ymax, N))
    return X, Y


def velocity_field(uinf, vinf, X, Y, xs, S, plot_option, x_airfoil, y_airfoil):
    """
    Velocity field and streamlines on a grid, point sources are distributed on the x-axis
    """
    Nx, Ny = X.shape
    print('Nx = %i, Ny = %i' % (X.shape))
    Ns = len(xs)
    eps = 1.e-5
    u, v = uinf * np.ones((X.shape), dtype=float), vinf * np.ones(X.shape, dtype=float)
    for i in range(Ns):
        r = (X - xs[i]) ** 2 + Y ** 2
        u = u + S[i] / (2 * np.pi) * (X - xs[i]) / (r + eps)
        v = v + S[i] / (2 * np.pi) * Y / (r + eps)

    if plot_option:
        y_start, y_end = Y[0, 0], Y[Ny - 1, 0]
        x_start, x_end = X[0, 0], X[0, Nx - 1]
        # print('X=',X[:,:]); print('Y=',Y[:,:])
        # print(x_start,x_end,y_start,y_end)

        plt.figure()
        plt.xlabel('x', fontsize=16);
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end);
        plt.ylim(y_start, y_end)
        plt.streamplot(X, Y, u, v, density=8.0, linewidth=1, arrowsize=1, arrowstyle='->')
        plt.scatter(xs, np.zeros((Ns), dtype=float), color='red', s=30, marker='o', linewidth=0)
        # To plot points use to calculate streamlines uncomment next line.
        # plt.scatter(X, Y,color='green', s=4, marker='o', linewidth=0);    
        plt.plot(x_airfoil, y_airfoil, 'r-')


def velocity_field_vector(uinf, vinf, X, Y, xs, S, plot_option, x_airfoil, y_airfoil, plt_show=False):
    """
    Velocity field  on a grid, point sources are distributed on the x-axis

    Plot of velocity vector at the panel control points
    """
    Nx = X.shape
    print('Nx = %i' % (Nx))
    Ns = len(xs)
    eps = 1.e-5
    u, v = uinf * np.ones(X.shape, dtype=float), vinf * np.ones(Y.shape, dtype=float)
    for i in range(Ns):
        r = (X - xs[i]) ** 2 + Y ** 2
        u = u + S[i] / (2 * np.pi) * (X - xs[i]) / (r + eps)
        v = v + S[i] / (2 * np.pi) * Y / (r + eps)
    # print('X:',len(X),X); print('Y:',len(Y),Y)
    umax = np.max(np.sqrt(u ** 2 + v ** 2))
    print('umax = %f' % (umax))
    u, v = u / umax, v / umax
    # print('u=',u); print('v=',v)
    if plot_option:
        x_start, x_end = -0.1, 1.5
        y_start, y_end = -0.2, 0.5
        size = 10
        plt.figure(figsize=(size, (y_end - y_start) / (x_end - x_start) * size))
        plt.xlabel('x', fontsize=16)
        plt.ylabel('y', fontsize=16)
        plt.xlim(x_start, x_end)
        plt.ylim(y_start, y_end)
        plt.quiver(X, Y, u, v, angles='xy', scale=5, width=0.003)  # velocity vector
        plt.scatter(xs, np.zeros((Ns), dtype=float), color='#CD2305', s=80, marker='o',
                    linewidth=0)  # position des sources
        plt.scatter(X, Y, color='green', s=2, marker='o', linewidth=1)  # control points
        plt.plot(x_airfoil, y_airfoil, 'r-')
        if plt_show:
            plt.show()


def plot_profile_with_panels(panels):
    """
    panel plot with their control points
    """
    print(len(panels))
    chord = panels[len(panels) - 2].xb
    print('panel number = ', len(panels))
    y = [panel.ya for panel in panels]
    width = 10
    plt.figure(figsize=(width, 0.8 * width))
    plt.grid()
    plt.title("Panels and control points")
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    # plt.plot(xe, ye, color='k', linestyle='-', linewidth=2)
    plt.plot([panel.xc for panel in panels], [panel.yc for panel in panels],
             linestyle='', linewidth=0, marker='s', markersize=5, color='blue', label='Control points')
    plt.plot(np.append([panel.xa for panel in panels], panels[0].xa),
             np.append([panel.ya for panel in panels], panels[0].ya),
             linestyle='-', linewidth=1, marker='o', markersize=6, color='#CD2305', label='Control points')
    plt.axis('scaled', adjustable='box')
    ech = 0.05
    plt.xlim(0 - ech, chord + ech)
    ymax = np.max(y) + ech
    plt.ylim(-ymax, ymax)
    plt.legend()
    return


def plot_source_panels(panels, plt_show=False):
    """
    plot of the source on the panels
    """
    chord = panels[len(panels) - 1].xb
    sigma = [panel.sigma for panel in panels]
    ymin, ymax = np.min(sigma), np.max(sigma)
    xc = [panel.xc for panel in panels]
    width = 10
    n = -1
    for i in range(len(panels)):
        if panels[i].epsilon == 1:
            n += 1
    print('n upper wall = ', n)
    N = len(panels)

    plt.figure(figsize=(width, 0.8 * width))
    plt.grid()
    plt.xlabel('x', fontsize=16)
    plt.ylabel('q', fontsize=16)
    plt.title('Sources')
    plt.plot([panel.xc for panel in panels], [panel.sigma for panel in panels], linestyle='-', linewidth=1, marker='o',
             markersize=6, color='#CD2305')
    plt.plot(xc[0:n], sigma[0:n], linestyle='-', linewidth=1, marker='o', markersize=6, color='red', label='upper w.')
    plt.plot(xc[n:N], sigma[n:N], linestyle='-', linewidth=1, marker='o', markersize=6, color='blue', label='lower w.')
    ech = 0.05
    plt.xlim(0 - ech, chord + ech)
    plt.ylim(ymin - 0.01, ymax + 0.01)
    if plt_show:
        plt.show()

def plot_doublet_panels(panels, plt_show=False):
    """
    plot of the doublets from the panels
    """
    chord = panels[len(panels) - 1].xb
    mu = [panel.mu for panel in panels]
    xc = [panel.xc for panel in panels]
    ymax = np.max(mu)
    ymin = np.min(mu)
    print('mu_min= %f \t mu_max = %f' % (ymin, ymax))
    width = 10
    plt.figure(figsize=(width, 0.8 * width))
    plt.grid()
    plt.xlabel('x', fontsize=16)
    plt.ylabel(r'$\mu$', fontsize=16)
    plt.title('Doublets')
    n = -1
    for i in range(len(panels)):
        if panels[i].epsilon == 1:
            n += 1
    print('n upper wall = ', n)
    N = len(panels)
    plt.plot([panel.xc for panel in panels], [panel.mu for panel in panels], linestyle='-', linewidth=1, marker='o',
             markersize=6, color='#CD2305')
    plt.plot(xc[0:n], mu[0:n], linestyle='-', linewidth=1, marker='o', markersize=6, color='red', label='upper w.')
    plt.plot(xc[n:N], mu[n:N], linestyle='-', linewidth=1, marker='o', markersize=6, color='blue', label='lower w.')
    ech = 0.05
    plt.xlim(0 - ech, chord + ech)
    plt.ylim(ymin - 0.01, ymax + 0.01)
    plt.legend()
    if plt_show: plt.show()
    S = 0
    for i in range(len(mu)):
        S += mu[i]
    print('sum(mu) for all panels  :', S)


def plot_doublet_derivative_panels(panels, plt_show=False):
    """
    plot of the doublets on the panels
    """
    chord = panels[len(panels) - 1].xb
    dmudx = [panel.dmudx for panel in panels]
    for i in range(len(dmudx)):
        if np.isinf(np.abs(dmudx[i])):
            dmudx[i] = 0
            print('infinite d mu/ dx  at %i' % (i))
    xb = [panel.xb for panel in panels]
    width = 10;
    ech = 0.05
    plt.figure(figsize=(width, 0.8 * width))
    plt.xlabel('x', fontsize=16)
    plt.ylabel(r'$d\mu/dx$', fontsize=16)
    plt.title('Doublet derivative')
    n = -1
    for i in range(len(panels)):
        if panels[i].epsilon == 1:
            n += 1
    N = len(panels)
    print('n upper wall= %i, N = %i' % (n, N))
    ymax = np.max(np.abs(dmudx[2:N - 2]))
    print('ymax = ', ymax)
    plt.plot(xb[0:n + 1], dmudx[0:n + 1], linestyle='-', linewidth=1, marker='o', markersize=6, color='red',
             label='upper w.')
    plt.plot(xb[n + 1:N - 1], dmudx[n + 1:N - 1], linestyle='-', linewidth=1, marker='o', markersize=6, color='blue',
             label='lower w.')
    plt.xlim(0 - ech, chord + ech)
    # plt.ylim(-ymax, ymax);
    plt.legend()
    plt.grid()
    if plt_show:
        plt.show()

def plot_profile_directions(panels, plt_show=False):
    """
    plot of the airfoil with the normal and tangential directions of the panels
    """
    x = [panel.xc for panel in panels]
    y = [panel.yc for panel in panels]
    x_start, x_end = -0.6, 1.5
    y_start, y_end = -0.5, 0.5
    size = 10
    s = 0.5         # vector length on the plot
    plt.figure(figsize=(size, (y_end - y_start) / (x_end - x_start) * size))
    plt.xlabel('x', fontsize=16)
    plt.ylabel('y', fontsize=16)
    plt.xlim(x_start, x_end)
    plt.ylim(y_start, y_end)
    plt.quiver(x, y, [panel.tx * s for panel in panels], [panel.ty * s for panel in panels], angles='xy', scale=5,
               width=0.003, color='red')
    plt.quiver(x, y, [panel.nx * s for panel in panels], [panel.ny * s for panel in panels], angles='xy', scale=5,
               width=0.003, color='blue')
    plt.scatter(x, y, color='green', s=2, marker='o', linewidth=1)  # control points
    plt.plot([panel.xa for panel in panels], [panel.ya for panel in panels], 'k-')
    plt.title('tangential and normal directions')
    # plt.axis('scaled', adjustable='box')
    if plt_show:
        plt.show()


def aerodynamic_coefficients(panels, c=1.0, alpha=0.0):
    """
    Evaluation of the aerodynamic coefficients from the pressure coefficient given on each panel
    """
    Alpha = np.deg2rad(alpha)
    # Ry=np.sum(panels[:].Kp*panels[:].length*panels[:].ny)/c
    # Rx=np.sum(panels[:].Kp*panels[:].length*panels[:].ny)/c
    Rx, Ry = 0, 0
    for i, panel in enumerate(panels):
        Ry -= panel.Kp * panel.length * panel.ny / c
        Rx -= panel.Kp * panel.length * panel.nx / c
        print('i = %i, Kp = %f, l = %f, nx = %f, ny = %f' % (i, panel.Kp, panel.length, panel.nx, panel.ny))
    print('Rx = %f, \t Ry = %f' % (Rx, Ry))
    CD = Rx * np.cos(Alpha) + Ry * np.sin(Alpha)
    CL = Ry * np.cos(Alpha) - Rx * np.sin(Alpha)
    print('CL = %f, \t CD = %f' % (CL, CD))
    return CL, CD


def aerodynamics(xe, xi, Kpe, Kpi, plt_show=False):
    """ 
    Aerodynamic coefficients
    """
    npt0 = 2 * len(Kpe) + 1
    # x=np.linspace(0.,1.,npt0)
    beta = np.linspace(0.0, np.pi, npt0)
    x = 0.5 * (1.0 - np.cos(beta))
    fe = interp1d(xe, Kpe, kind='linear')
    fi = interp1d(xi, Kpi, kind='linear')
    kpe, kpi = fe(x), fi(x)
    DeltaKp = kpi - kpe
    plt.figure()
    plt.xlabel('x/c', fontsize=16)
    plt.ylabel(r'$\Delta K_p$', fontsize=16)
    plt.plot(x, DeltaKp, linestyle='-', linewidth=2, color='red')
    plt.plot(xe, Kpi - Kpe, linestyle='', linewidth=1, marker='o', markersize=4, markevery=2, color='black')
    plt.plot(xi, Kpi - Kpe, linestyle='', linewidth=1, marker='s', markersize=4, markevery=2, color='blue')
    plt.title('Load')

    Cle, Cli = np.trapz(Kpe, x=xe), np.trapz(Kpi, x=xi)
    Cl = np.trapz(DeltaKp, x)
    Cm = np.trapz(x * DeltaKp, x)
    print('Aerodynamic:\n CL = %f, Cle = %f, Cli = %f, CL= Cli-CLe = %f ' % (Cl, Cle, Cli, (Cli - Cle)))
    print('Cm0 = ', Cm, ', x center of lift = ', Cm / Cl)

    width = 10
    ech = 0.05
    plt.figure(figsize=(width, width))
    plt.xlabel('x', fontsize=16)
    plt.ylabel('Kp', fontsize=16)
    plt.title('Pressure coefficient')
    plt.plot(xe, Kpe, linestyle='-', linewidth=1, marker='o', markersize=8, color='red', label='upper_wall')
    plt.plot(xi, Kpi, linestyle='-', linewidth=1, marker='o', markersize=8, color='blue', label='lower_wall')
    plt.plot(x, kpe, linestyle='-', color='black', label='upper_wall')
    plt.plot(x, kpi, linestyle='--', color='black', label='lower_wall')
    plt.legend(loc='lower right')
    plt.grid()
    plt.xlim(0 - ech, 1 + ech)
    if plt_show: plt.show()

    print('xe \t  Kpe \t\t xi \t  Kpi \t\t\t x \t\t  kpe \t  kpi')
    for i in np.arange(0, 5):
        print('%f %f \t %f %f \t | \t %f \t %f  %f ' % (xe[i], Kpe[i], xi[i], Kpi[i], x[i], kpe[i], kpi[i]))
    ne = len(xe)
    print('xe \t  Kpe \t\t xi \t  Kpi')
    for i in np.arange(ne - 6, ne):
        print('%f %f \t %f %f ' % (xe[i], Kpe[i], xi[i], Kpi[i]))
