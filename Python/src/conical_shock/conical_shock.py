# -*- coding: utf-8 -*-
"""

Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2020/04/10*

*licence : AGPL-3.0*

Module to extract solution of conical shocks from a database
  ..

The database has been build with a Fortran code also available in FundAeroSuite
"""
import glob
import os
import numpy             as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


__fileroot__ = "conical_shock/MACH/Mach_"
first_Mach = len(__fileroot__)
last_Mach = first_Mach + 5


def Mach_list(display=False):
    """
    List of Mach number available in the database
    """
    File_List = glob.glob('conical_shock/MACH/*.dat')
    # print(File_List)
    Mach = []
    for file in File_List:
        Mach.append((int(file[first_Mach:last_Mach]) / 1000.))
    Mach = sorted(Mach)
    print("Number of Mach found : ", len(Mach))
    if display:
        print(Mach)
    return Mach


def find_value(Mach, M0):
    """
    search if a Mach number is in the list of Mach available
    
    Params:
        * Mach (list of real) : list of Mach available
        * M0 (real) : Mach number define by the user
    
    Returns:
        int : m , index of the list to find the input Mach M0
    """
    try:
        m = Mach.index(M0)
        print('index of Mach list = %f  : %d' % (M0, m))
        return m
    except:
        print("The value of  M=%f is not in the database" % M0)
        mn, idx = min((abs(Mach[i] - M0), i) for i in range(len(Mach)))
        if Mach[idx] < M0:
            print('Nearest Mach value available : ', Mach[idx:idx + 2])
        else:
            print('Nearest Mach value available : ', Mach[idx - 1:idx + 1])
        return -1


def get_data_from_database(M0):
    """
    for a given Mach number, get the data from the database
    """
    filepath = os.path.join('conical_shock/MACH', 'Mach_%05d.dat' % (M0 * 1000))
    print('File to read : ', filepath)

    # content of a file, by column
    # 0 :  Theta at wall in deg 
    # 1 :  shock angle in deg      
    # 2 :  Kp           
    # 3 :  Inverse Mach number      
    # 4 :  upstream Mach number
    # 5 :  downstream Mach number

    with open(filepath, 'r') as infile:
        A = np.loadtxt(infile, dtype=float, unpack=True, skiprows=2)
        # theta1=A[0,:];kp_cc=A[2,:] ; sigma1=A[1,:]
    d = A.shape
    data = np.zeros([d[0], d[1] + 1])
    print('data table size : ', data.shape)
    # I have to build the first column when  theta=0
    data[:, 0] = [0.0, np.rad2deg(np.arcsin(1 / M0)), 0, 1 / M0 - 1, M0, M0]
    data[:, 1:] = A
    return data


def Solution_Lees(gamma=1.4):
    """
    LEES's solution on supersonic cone flow  
    
    Returns:
        real: Kp/theta and sigma/theta
    """
    return 2 * (gamma + 1) * (gamma + 7) / (gamma + 3) ** 2, 2 * (gamma + 1) / (gamma + 3)


def main_single_Mach(M0=2.0, theta_c=np.float64(15.0), display_values=True):
    """
    main function to get data for a single Mach as input
    
    Params:
        * M0 (real) : upstream Mach number
        * theta_c (real) : deviation angle in deg
        * display_values (boolean) : True  print values in the terminal
    """
    # Main options :
    plot_sigma = False
    plot_Kp = False
    plot_downstream_Mach = True
    display = False
    interpolate_option = 'linear'  # linear or cubic, linear in case of problem
    print("interpolation option: ", interpolate_option)
    lw = 2  # line width in the plots

    Mach = Mach_list()
    m = find_value(Mach, M0)
    if m == -1:
        return
    data = get_data_from_database(M0)

    if display:
        if plot_sigma:
            fig = plt.figure(1)
            fig.suptitle(r'conical shock solution for M0 = %f' % M0, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
            ax.set_ylabel(r'$\sigma (^\circ)$', fontsize=20)
            ax.grid()
            ax.plot(data[0, :], data[1, :], '-', label=r"$\sigma$", linewidth=lw)
        if plot_Kp:
            fig = plt.figure(2)
            fig.suptitle(r'conical shock solution for M0 = %f' % M0, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
            ax.set_ylabel(r'$100\times Kp$', fontsize=20)
            ax.grid()
            ax.plot(data[0, :], 100 * data[2, :], '-', label=r"$100 \times K_p$", linewidth=lw)
        if plot_downstream_Mach:
            fig = plt.figure(3)
            fig.suptitle(r'conical shock solution for M0 = %f' % M0, fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
            ax.set_ylabel(r'$M_{down.}$', fontsize=20)
            ax.grid()
            ax.plot(data[0, :], data[5, :], '-', label=r"$M_{downstream}$", linewidth=lw)

            # ax.legend(loc='upper right')
        plt.show()

    if display_values:
        Kp_func = interp1d(data[0, :], data[2, :], kind=interpolate_option)
        sigma_func = interp1d(data[0, :], data[1, :], kind=interpolate_option)
        downstream_Mach_func = interp1d(data[0, :], data[5, :], kind=interpolate_option)
        print('Kp = %f, sigma= %f °, downstream Mach number = %f' % (
        Kp_func(theta_c), sigma_func(theta_c), downstream_Mach_func(theta_c)))


def main_multiple_Mach(M0, liste_plot=("sigma"), Lees=False, t=(20, 2, 1), s=(26, 5, 1), Kp=(30, 5, 1), Mav=(50, 5, 1),
                       gam=1.4):
    """
    main function to get data with multiple Mach as input (last version)
    
    Params:
    
    * M0 (real table) : upstream Mach number list
    * liste_plot (list of string) : the list of what the user want to plot,"sigma", "downstream_Mach", "Kp"
    * Lees (boolean) : True: Lees's solution is given
    * t (real list) : theta angle parameters in deg: theta_max, grille
    * s (real list) : sigma angle parameters in deg: sigma_max, grille    
    * Kp (real list) : 100 x Kp : max, grille
    * Mav (real list) : downstream Mach : max, grille
    """
    if Lees:
        c_Kp, c_sigma = Solution_Lees(gamma=gam)
        print("Lees's coefficients : c_Kp = %f, c_sigma = %f" % (c_Kp, c_sigma))
        npt = 21

    plot_sigma, plot_Kp, plot_downstream_Mach = False, False, False
    for l in liste_plot:
        if l == "sigma":
            plot_sigma = True
        elif l == "downstream_Mach":
            plot_downstream_Mach = True
        elif l == "Kp":
            plot_Kp = True

    print("plot options: ", plot_sigma, plot_downstream_Mach, plot_Kp)
    display = True
    interpolate_option = 'linear'  # linear or cubic, choose linear is case of problem
    print("Interpolation option : ", interpolate_option)
    lw = 2
    print("input Mach table", M0)

    Mach = Mach_list(display=False)
    # test if the values are in the database
    DATA = []
    for Ma in M0:
        m = find_value(Mach, Ma)
        if m == -1:
            raise NameError('Modify the value of the Mach number: STOP')
        DATA.append([Ma, get_data_from_database(Ma)])
        # DATA[n][k][j][i]
    # n : Mach index in the table
    # k : 0 : M0, 1: data table
    # j : column of the data file
    # i : row of the data file
    N = len(M0)

    """
    The y-axis detail have to be modified depending the option choosen
    """

    xmax = 20
    xmax = t[0]
    if display:
        xmajor_ticks, xminor_ticks = np.arange(0, xmax, t[1]), np.arange(0, xmax, t[2])
        if plot_sigma:
            ymax = s[0]
            ymajor_ticks, yminor_ticks = np.arange(0, ymax, s[1]), np.arange(0, ymax, s[2])
            fig = plt.figure(1)
            fig.suptitle(r'Conical shock: shock angle in deg', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
            ax.set_ylabel(r'$\sigma (^\circ)$', fontsize=20)
            ax.grid(which='both')
            for k in range(N):
                ax.plot(DATA[k][1][0][:], DATA[k][1][1][:], '-', label=r"$M_0=$%4.2f" % (M0[k]), linewidth=lw)
            if Lees:
                theta = np.linspace(0, xmax, npt)
                ax.plot(theta, theta * c_sigma, 'ko-', markersize=5, label=r"Lees")
            ax.legend(loc='lower right', fontsize=12)
            ax.set_xticks(xmajor_ticks)
            ax.set_xticks(xminor_ticks, minor=True)
            ax.set_yticks(ymajor_ticks)
            ax.set_yticks(yminor_ticks, minor=True)
            ax.axis([0.0, xmax, 0, ymax])
            ax.tick_params(labelsize=12)

        if plot_Kp:
            ymax = Kp[0]
            ymajor_ticks, yminor_ticks = np.arange(0, ymax, Kp[1]), np.arange(0, ymax, Kp[2])
            fig = plt.figure(2)
            fig.suptitle(r'Conical shock : $Kp$', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
            ax.set_ylabel(r'$100 \times Kp$', fontsize=20)
            ax.grid(which='both')
            for k in range(N):
                ax.plot(DATA[k][1][0][:], 100 * DATA[k][1][2][:], '-', label=r"$M_0=$%4.2f" % (M0[k]), linewidth=lw)
            if Lees:
                theta = np.linspace(0, xmax, npt)
                ax.plot(theta, np.deg2rad(theta) ** 2 * 100 * c_Kp, 'ko-', markersize=5, label=r"Lees")
            ax.legend(loc='upper left', fontsize=12)
            ax.set_xticks(xmajor_ticks)
            ax.set_xticks(xminor_ticks, minor=True)
            ax.set_yticks(ymajor_ticks)
            ax.set_yticks(yminor_ticks, minor=True)
            ax.axis([0.0, xmax, 0, ymax])
            ax.tick_params(labelsize=12)

        if plot_downstream_Mach:
            # ymax=DATA[N-1][0]   # max is calculated
            ymax = Mav[0]
            ymajor_ticks, yminor_ticks = np.arange(1, ymax, Mav[1]), np.arange(0, ymax, Mav[2])
            fig = plt.figure(3)
            fig.suptitle(r'Conical shock : downstream Mach number ', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$\theta (^\circ)$', fontsize=20)
            ax.set_ylabel(r'$M_{aval}$', fontsize=20)
            ax.grid(which='both')
            for k in range(N):
                ax.plot(DATA[k][1][0][:], DATA[k][1][5][:], '-', label=r"$M_0=$%4.2f" % (M0[k]), linewidth=lw)
            ax.legend(loc='upper right')
            ax.set_xticks(xmajor_ticks)
            ax.set_xticks(xminor_ticks, minor=True)
            ax.set_yticks(ymajor_ticks)
            ax.set_yticks(yminor_ticks, minor=True)
            ax.axis([0.0, xmax, 0, ymax])
            ax.tick_params(labelsize=12)
        plt.show()


def main_deviation_angle_maxi(MachMax=7.1):
    """
    get the maximal of deviation angle  
    and the pressure coefficient Kp for this angle
    """
    interpolate_option = 'linear'
    thetaM = 10.0
    filepath = os.path.join('conical_shock/COURBES', 'theta_cone_max.dat')
    print('File to read : ', filepath)
    # content of a file, by column
    # 0 : upstream Mach number  
    # 1 : shock angle in deg      
    # 2 : theta_max in deg          
    # 3 : Inverse Mach number on the cone wall       
    # 4 : Kp
    with open(filepath, 'r') as infile:
        A = np.loadtxt(infile, dtype=float, unpack=True, skiprows=1)
    print(A.shape)
    Mach = A[0, :]
    thetac = A[1, :]
    thetaMax = A[2, :]
    Kpc = A[4, :]
    f = interp1d(thetaMax, Mach, kind=interpolate_option)
    print('minimal Mach  = ', f(thetaM))

    # Figure Mach - Kpc
    xmax = MachMax
    ymax = 1.65
    xmajor_ticks, xminor_ticks = np.arange(1, xmax, 1), np.arange(1, xmax, 0.2)
    ymajor_ticks, yminor_ticks = np.arange(0, ymax, 0.4), np.arange(0, ymax, 0.1)
    fig = plt.figure(0, figsize=(10, 8))
    fig.suptitle(r'Conical shock : $Kp(\theta_{\max})$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$Mach$', fontsize=20)
    ax.set_ylabel(r'$Kp$', fontsize=20)
    ax.plot(Mach, Kpc)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.axis([1.0, xmax, 0, ymax])
    ax.grid(which='both')
    ax.tick_params(labelsize=12)

    # Figure Mach - Kpc ZOOM
    xmax = 1.61
    ymax = 1.05
    xmajor_ticks, xminor_ticks = np.arange(1, xmax, 0.2), np.arange(1, xmax, 0.05)
    ymajor_ticks, yminor_ticks = np.arange(0, ymax, 0.2), np.arange(0, ymax, 0.05)
    fig = plt.figure(1, figsize=(10, 8))
    fig.suptitle(r'Conical shock : $Kp(\theta_{\max})$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$Mach$', fontsize=20)
    ax.set_ylabel(r'$Kp$', fontsize=20)
    ax.plot(Mach, Kpc)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.axis([1.0, xmax, 0, ymax])
    ax.grid(which='both')
    ax.tick_params(labelsize=12)

    # Figure Mach - theta_max
    xmax = 7.1
    ymax = 61
    xmajor_ticks, xminor_ticks = np.arange(1, xmax, 1), np.arange(1, xmax, 0.2)
    ymajor_ticks, yminor_ticks = np.arange(0, ymax, 10), np.arange(0, ymax, 5)
    fig = plt.figure(2, figsize=(10, 8))
    fig.suptitle(r'Conical shock : $\theta_{\max}$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$Mach$', fontsize=20)
    ax.set_ylabel(r'$\theta_{\max}$', fontsize=20)
    ax.plot(Mach, thetaMax)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.axis([1.0, xmax, 0, ymax])
    ax.grid(which='both')
    ax.tick_params(labelsize=12)

    # Figure Mach - theta_max ZOOM
    xmax = 1.61
    ymax = 36
    xmajor_ticks, xminor_ticks = np.arange(1, xmax, 0.2), np.arange(1, xmax, 0.05)
    ymajor_ticks, yminor_ticks = np.arange(0, ymax, 10), np.arange(0, ymax, 2)
    fig = plt.figure(3, figsize=(10, 8))
    fig.suptitle(r'Conical shock : $\theta_{\max}$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'$Mach$', fontsize=20)
    ax.set_ylabel(r'$\theta_{\max}$', fontsize=20)
    ax.plot(Mach, thetaMax)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.axis([1.0, xmax, 0, ymax])
    ax.grid(which='both')
    ax.tick_params(labelsize=12)
    plt.show()


def main_plot_Mach_Kp(Mach_intervalle, theta_c):
    """
    for a list of Mach number, get the pressure coefficient for a given cone angle
    from the database   
    """
    gam = 1.4
    Mmax = Mach_intervalle[1]
    interpolate_option = 'linear'  # linear or cubic, choose linear is case of problem
    File_List = glob.glob('conical_shock/MACH/*.dat')
    Mach = []
    MachList = []
    for file in File_List:
        M = int(file[first_Mach:last_Mach]) / 1000.
        if Mach_intervalle[0] <= M <= Mach_intervalle[1]:
            Mach.append(int(file[first_Mach:last_Mach]) / 1000.)
            MachList.append(file)
    Mach = sorted(Mach)
    print("length of the Mach list : ", len(Mach))
    print("Mach list : ", Mach)
    Mach, Kp = [], []
    MachList = sorted(MachList)
    for file in MachList:
        with open(file, 'r') as infile:
            data = np.loadtxt(infile, dtype=float, unpack=True, skiprows=2)
        Kp_func = interp1d(data[0, :], data[2, :], kind=interpolate_option)
        Mach.append(int(file[first_Mach:last_Mach]) / 1000.)
        Kp.append(Kp_func(theta_c))
    M = np.linspace(5, Mmax, 51)
    # ymax=170
    # ymajor_ticks,yminor_ticks = np.arange(0, ymax, 20),np.arange(0, ymax, 4)
    ymin, ymax = 0.05, 0.23
    xmajor_ticks, xminor_ticks = np.arange(1, Mmax, 1), np.arange(1, Mmax, 0.2)
    ymajor_ticks, yminor_ticks = np.arange(ymin, ymax, 0.05), np.arange(ymin, ymax, 0.01)
    fig = plt.figure(1, figsize=(10, 8))
    fig.suptitle(r'Conical shock : $Kp$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    # fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$M$', fontsize=20)
    ax.set_ylabel(r'$Kp$', fontsize=20)
    ax.grid(which='both')
    ax.plot(Mach, Kp, linewidth=2)
    ax.plot(M, 2 * np.sin(theta_c * np.pi / 180) ** 2 * np.ones(M.shape))
    ax.plot(M, 2 * (gam + 1) * (gam + 7) / (gam + 3) ** 2 * (theta_c / 180 * np.pi) ** 2 * np.ones(M.shape))
    print(2 * (gam + 1) * (gam + 7) / (gam + 3) ** 2 * (theta_c / 180 * np.pi) ** 2,
          2 * np.sin(theta_c * np.pi / 180) ** 2)
    # ax.legend(loc='upper left',fontsize=12)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.axis([1.0, Mmax, ymin, ymax])
    ax.tick_params(labelsize=12)

    ymin, ymax = 0.11, 0.23
    ymajor_ticks, yminor_ticks = np.arange(ymin, ymax, 0.02), np.arange(0, ymax, 0.01)
    xmajor_ticks, xminor_ticks = np.arange(1, 1.6, 0.2), np.arange(1, 1.6, 0.05)
    fig = plt.figure(2, figsize=(10, 8))
    fig.suptitle(r'Conical shock : $Kp$', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    # fig.subplots_adjust(top=0.80)
    ax.set_xlabel(r'$M$', fontsize=20)
    ax.set_ylabel(r'$Kp$', fontsize=20)
    ax.plot(Mach, Kp, linewidth=2)
    # ax.legend(loc='upper left',fontsize=12)
    ax.set_xticks(xmajor_ticks)
    ax.set_xticks(xminor_ticks, minor=True)
    ax.set_yticks(ymajor_ticks)
    ax.set_yticks(yminor_ticks, minor=True)
    ax.axis([1.0, 1.6, ymin, ymax])
    ax.grid(which='both')
    ax.tick_params(labelsize=12)
    plt.show()


# *********************************************
#  main programm
# *********************************************

def run_conical_shock(option):
    """
    main program to get conical shock results
    """
    if option == 0:
        """
        single Mach as input, must be in the database.
        """
        main_single_Mach(M0=2.0, theta_c=np.float64(15.0), display_values=True)

    elif option == 1:
        """
        Downstream Mach number, Kp, sigma as a function of theta_c, the cone angle
        """
        main_multiple_Mach((1.05, 1.1, 1.3, 1.5, 1.6, 2.0, 5), ("sigma", "downstream_Mach", "Kp"))
    elif option == 2:
        """
        theta_x_Max and Kp  as a function of a upstream Mach number
        """
        main_deviation_angle_maxi(MachMax=7.1)
    else:
        """
        Kp as a function of the upstream Mach, for a given cone angle
        """
        main_plot_Mach_Kp((1.055, 7.), 10.)


if __name__ == "__main__":
    run_conical_shock(0)
