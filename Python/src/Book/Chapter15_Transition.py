#!/bin/py
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 15 exercises : transitional flows
    ..
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Tools.misc import set_title, set_question
from scipy.interpolate import interp1d


def Michel(Re_x):
    """
    Michel's transition criteria
    """
    return 1.535 * Re_x ** 0.444


def Cebeci_Smith(Re_x):
    """
    Cebeci-Smith's transition criteria
    """
    return 1.174 * (1 + 22400. / Re_x) * Re_x ** 0.46


def Re_crit(H):
    """
    Critical Reynolds number
    """
    if H >= 2.591:
        Rec = 54.2124 / H / (H - 2.48) + 31.6 / H
    else:
        Rec = 520. / H + 2.5e6 / H * (1. / H - 1. / 2.591) ** 1.95
    return Rec


def Granville(m, H, k2):
    """
    Granville's transition criteria
    """
    Lambda = m * k2 ** 2
    return Re_crit(H) + 375.0 + np.exp(6.1 + 55 * Lambda)


def Arnal(m, H, k2, Tu):
    """
    Arnal's transition criteria
    """
    Lambda = m * k2 ** 2
    return Re_crit(H) - 206 * np.exp(25.7 * Lambda) * (np.log(16.8 * Tu) - 2.77 * Lambda)


def FSK_correlation(beta):
    """
    Correlation in the range   -0.14 < beta < 0.11
    fo a Falkner-Skan boundary layer
    Calculation perfomed by C.  Airiau with a continuation method
    """
    Lambda = 4.5157e-5 + 0.2196 * beta - 0.3872 * beta ** 2 + 0.7967 * beta ** 3
    H = 2.5917 - 1.404 * beta + 3.885 * beta ** 2 - 20.516 * beta ** 3 + 104.06 * beta ** 4
    k_2 = 0.664 - 0.71153 * beta + 0.9431 * beta ** 2 - 1.9697 * beta ** 3
    return k_2, H, Lambda


def print_correlations():
    """
    display the reference data for a  Falkner-Skan boundary layer
    """
    set_title("Correlation for a  Falkner-Skan boundary layer")
    n = 51
    beta = np.linspace(-0.14, 0.11, n)
    m = beta / (2 - beta)
    print('beta \t\t  m \t\t k_2 \t\t H \t\t Lambda')
    for b in beta:
        k_2, H, Lambda = FSK_correlation(b)
        m_tmp = b / (2. - b)
        print('%2.5f \t' * 5 % (b, m_tmp, k_2, H, Lambda))
    return


def plot_correlations():
    """
    Calculus and plot of the correlations for FSK
    """
    # print_correlations()
    # read FSK files
    reference_filepath = os.path.join('Book/Data', 'BL_FSK_characteristics.dat')
    with open(reference_filepath, 'r') as infile:
        A = np.loadtxt(infile, dtype=float, unpack=True, skiprows=1)
    #     beta, m, H, s, d1/d, d2/d, Cf, xk, Eta Max, Lambda, End_BL
    print(A.shape)
    beta, m, H, k_2, Lambda = A[[0, 1, 2, 5, 9], :]
    k2_coef = np.polyfit(beta, k_2, 3)
    p_k2 = np.poly1d(k2_coef)
    H_coef = np.polyfit(beta, H, 4)
    p_H = np.poly1d(H_coef)
    Lambda_coef = np.polyfit(beta, Lambda, 3)
    p_Lambda = np.poly1d(Lambda_coef)
    eps = 1e-10
    print('polynomial k2 :', k2_coef)
    print('polynomial H :', H_coef)
    print('polynomial Lambda :', Lambda_coef)
    k2_1, H_1, Lambda_1 = FSK_correlation(beta)
    err_k2 = (k2_1 - p_k2(beta)) / (p_k2(beta) + eps)
    err_H = (H_1 - p_H(beta)) / (p_H(beta) + eps)
    err_Lambda = (Lambda_1 - p_Lambda(beta)) / (p_Lambda(beta) + eps)

    fig = plt.figure()
    fig.suptitle('correlation errors', fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_title('correlations')
    ax.set_xlabel(r'$\beta$', fontsize=20)
    ax.grid()
    ax.semilogy(beta, abs(err_k2), '-', label=r"$k_2$", linewidth=2, color='black')
    ax.semilogy(beta, abs(err_H), '-', label='H', linewidth=2, color='blue')
    ax.semilogy(beta, abs(err_Lambda), '-', label=r"$\Lambda$", linewidth=2, color='red')
    ax.legend(loc='upper left')
    plt.show()

    m_test = 0.
    beta_test = 2. * m_test / (m_test + 1.)
    print('beta_test = %1.3f \t m_test = %1.3f' % (beta_test, m_test))
    print('k2 = %1.4f, \t, H = %1.4f, \t, Lambda = %1.4f, %1.4f' % (
    p_k2(beta_test), p_H(beta_test), p_Lambda(beta_test), m_test * p_k2(beta_test) ** 2))
    return


def Exercice15_9():
    """ 
    Comparison of criteria transition 
    """

    plot = True
    flag_correlation = True
    Tu = 0.001      # turbulence rate
    n = 51          # size of  Re_x vector
    k = 4           # index of table for m

    if flag_correlation:
        plot_correlations()

    set_title("Comparison of criteria transition")
    # Falkner-Skan data
    m = np.array([1.0, 1. / 3., 1. / 7., 0.0, -1. / 21., -0.0904])
    beta = 2 * m / (m + 1)
    k2 = np.array([0.2923, 0.4290, 0.5245, 0.6641, 0.7464, 0.8681])  # Rdelta2/Rdelta
    H = np.array([2.2162, 2.2970, 2.3843, 2.5911, 2.8011, 4.0292])

    if k == 4:
        Rex = (np.linspace(600, 800, n)) ** 2  # m=-1/21
    else:
        Rex = (np.linspace(1400, 1900, n)) ** 2  # m=0

    Id = np.ones(n)

    meth = ['Michel', 'Cebeci-Smith', 'Granville', 'Arnal']
    print(
        ' m = %1.3f, beta = %1.3f \t k2= %1.4f , \t Lambda= m k2^2 = %1.4f' % (m[k], beta[k], k2[k], m[k] * k2[k] ** 2))
    print('Re delta  crit = %f, \t Re delta2  crit = %f' % (Re_crit(H[k]) / k2[k], Re_crit(H[k])))
    print('Lambda = ', m[k] * k2[k] ** 2, ' H = ', H[k])

    if plot:
        fig = plt.figure()
        fig.suptitle('transition criteria', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_title("Michel's criteria")
        ax.set_xlabel(r'$Re_x$', fontsize=20)
        # ax.set_ylabel(r'$\frac{P_i}{P_t}-1$',fontsize=20)
        ax.grid()
        # ax.axis([0., 1, 0, 20])
        ax.plot(np.sqrt(Rex), k2[k] * np.sqrt(Rex), '-', label='Falkner-Skan', linewidth=2,
                color='black')  # r'$Re_{\delta_2} $'
        ax.plot(np.sqrt(Rex), Michel(Rex), '-', label='Michel')
        ax.plot(np.sqrt(Rex), Cebeci_Smith(Rex), '-', label='Cebeci-Smith')
        ax.plot(np.sqrt(Rex), Id * Granville(m[k], H[k], k2[k]), '-', label='Granville')
        ax.plot(np.sqrt(Rex), Id * Arnal(m[k], H[k], k2[k], Tu), '-', label='Arnal')
        ax.legend(loc='upper left')
        plt.show()

    Re = np.zeros((4, n))
    Rdeltk2 = k2[k] * np.sqrt(Rex)
    Re[0, :] = Michel(Rex) - Rdeltk2
    Re[1, :] = Cebeci_Smith(Rex) - Rdeltk2
    Re[2, :] = Id * Granville(m[k], H[k], k2[k]) - Rdeltk2
    Re[3, :] = Id * Arnal(m[k], H[k], k2[k], Tu) - Rdeltk2

    # We got analytical solution for Granville and Arnal's criteria for FSK's case.
    # I'am looking for  zeros the functions by linear interpolation 
    for i in np.arange(4):
        f1 = interp1d(Re[i, :], np.sqrt(Rex), kind='cubic')
        print('Re delta= %5.1f ,\t Method = %s, ' % (f1(0), meth[i]))

    print("Granville's criteria : analytical condition : ")
    print(Granville(m[k], H[k], k2[k]) / k2)
    print("Arnal's criteria : analytical solution")
    print(Arnal(m[k], H[k], k2[k], Tu) / k2)
