#!/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 10 exercises 
    ..
"""

import numpy as np
import matplotlib.pyplot as plt
from CompressibleFlow.fonctions import gamma, p_pi, S_over_Scrit, P2_P1, downstream_Mach, pi2_pi1, rho_rhoi, T_Ti
from CompressibleFlow.fonctions import inv_S_over_Scrit, inverse_p_pi, r_Air, qm_1D, sound_velocity, rho2_rho1, \
    inv_pi2_over_pi1
from CompressibleFlow.tables import Atmosphere
from Tools.misc import set_title, set_question, simple_plot

Largeur, Hauteur = 10, 8


def valX(M, gamma=gamma):
    """
    X as a function of  Mach number,
    X is a function of the ratio pressure
    """
    omega = 1. + (gamma - 1.) / 2. * M ** 2
    return 1. / omega


def u_over_ac(M, gamma=gamma):
    """
    u/a_c as a function of Mach number
    """
    omega = 1 + (gamma - 1.) / 2. * M ** 2
    return M * np.sqrt((gamma + 1) / (2. * omega))


def S_over_Sc(x, gamma=gamma):
    """
    S / Sc as a function of  x=p/pi_0
    """
    z = pow(x, -1. / (gamma - 1.))
    tmp = np.sqrt((gamma - 1.0) / (gamma + 1)) * pow(2. / ((gamma + 1.)), 1. / (gamma - 1.))
    return tmp * z / np.sqrt(1 - z)


def CF(M, gamma=gamma):
    """
    propulsion coefficient
    """
    omega = 1 + (gamma - 1.) / 2. * M ** 2
    n = (gamma + 1.) / (2. * (gamma - 1.))
    return gamma * pow((gamma + 1.) / 2., -n) * M / np.sqrt(omega)


def CF_m(M, gamma=gamma):
    """
    Propulsion coefficient given in Shapiro' book after error correction
    """
    tmp1 = 1 - pow(p_pi(M), (gamma - 1) / gamma)
    n = (gamma + 1.) / (gamma - 1.)
    tmp2 = 2. / (gamma - 1) * pow(2 / (gamma + 1.), n)
    return gamma * np.sqrt(tmp1 * tmp2)


def CF_add(M, rpa, gamma=gamma):
    """
    additional part of the propulsion coefficient due to pression
    """
    return (p_pi(M, gamma=gamma) - rpa) * S_over_Scrit(M, gamma=gamma)


def CFmax(r, gamma=gamma):
    """
    Maximal thrust coefficient as a function of  r = pa/pi0
    """
    z = pow(r, (gamma - 1.) / gamma)
    tmp = pow((gamma + 1.) / 2., -(gamma + 1) / (2 * (gamma - 1.)))
    return gamma * tmp * np.sqrt(2. / (gamma - 1.)) * np.sqrt(1 - z)


def plotCF(rpa, borne, n, option, ifig, gamma=gamma):
    """
    plot of the propulsion coefficient CF
    """

    fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
    # fig.suptitle(r'Thrust $F/F_0$ for $p_a/p_{i_0}=$%e'%(rpa), fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.80)
    ax.set_ylabel(r'$F/F_0$', fontsize=20)
    if option == 'M':
        ax.set_xlabel(r'$M$', fontsize=20)
        for k in range(len(rpa)):
            r = rpa[k]
            Mtmp = np.linspace(borne[k][0], borne[k][1], n)
            y = CF(Mtmp, gamma=gamma) + CF_add(Mtmp, r, gamma=gamma)
            ax.plot(Mtmp, y, linewidth=2, label=r'$\frac{p_a}{p_{i_0}}=$%6.5f' % r)
    elif option == 'P':
        ax.set_xlabel(r'$p_s/p_{i_0}$', fontsize=20)
        for k in range(len(rpa)):
            r = rpa[k]
            Mtmp = np.linspace(borne[k][0], borne[k][1], n)
            y = CF(Mtmp, gamma=gamma) + CF_add(Mtmp, r, gamma=gamma)
            ax.plot(p_pi(Mtmp), y, linewidth=2, label=r'$\frac{p_a}{p_{i_0}}=$%6.5f' % r)
    elif option == 'Pa':
        ax.set_xlabel(r'$p_s/p_a$', fontsize=20)
        for k in range(len(rpa)):
            r = rpa[k]
            Mtmp = np.linspace(borne[k][0], borne[k][1], n)
            y = CF(Mtmp, gamma=gamma) + CF_add(Mtmp, r, gamma=gamma)
            ax.plot(p_pi(Mtmp) / r, y, linewidth=2, label=r'$\frac{p_a}{p_{i_0}}=$%6.5f' % r)
    elif option == 'S':
        ax.set_xlabel(r'$S_s/S_c$', fontsize=20)
        for k in range(len(rpa)):
            r = rpa[k]
            Mtmp = np.linspace(borne[k][0], borne[k][1], n)
            y = CF(Mtmp, gamma=gamma) + CF_add(Mtmp, r, gamma=gamma)
            ax.plot(S_over_Scrit(Mtmp), y, linewidth=2, label=r'$\frac{p_a}{p_{i_0}}=$%6.5f' % r)
    else:
        ax.set_xlabel(r'$X$', fontsize=20)
        for k in range(len(rpa)):
            r = rpa[k]
            Mtmp = np.linspace(borne[k][0], borne[k][1], n)
            y = CF(Mtmp, gamma=gamma) + CF_add(Mtmp, r, gamma=gamma)
            ax.plot(valX(Mtmp), y, linewidth=2, label=r'$\frac{p_a}{p_{i_0}}=$%6.5f' % r)
    ax.grid()
    ax.legend(loc='upper right')
    Stmp = S_over_Scrit(Mtmp)
    Xtmp = valX(Mtmp)
    Ptmp = p_pi(Mtmp)
    print('range M              = %f \t %f ' % (Mtmp[0], Mtmp[-1]))
    print('range ps/pi0         = %f \t %f ' % (Ptmp[0], Ptmp[-1]))
    for r in rpa:
        print('range ps/pa      = %f \t %f pour rpa= %f  :  ' % (Ptmp[0] / r, Ptmp[-1] / r, r))
        print('CF_max                              = %f ' % (CFmax(r)))

    print('range S              = %f \t %f ' % (Stmp[0], Stmp[-1]))
    print('range X              = %f \t %f ' % (Xtmp[0], Xtmp[-1]))


def Exercice10_1():
    """ 
    Rocket nozzle thrust
    """
    ifig = 0
    set_title("Rocket nozzle thrust")
    plot_ps = False  # True : usual plots
    plot_CF = True
    # legM = 'M'        # 'M', 'P' ou 'S' ou 'Pa' ou 'X'
    Sc = 10e-4  #  m^2
    pi0 = 10e5  # Pa
    rpa = [0.005, 0.01]
    section_ratio = 10.  # section ratio
    n = 31
    Ss = Sc * section_ratio

    if plot_ps:
        x = np.linspace(0.1, 4, n)
        fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
        fig.suptitle(r' $p/p_{i_0}=f(M)$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$p/p_{i_0}$', fontsize=20)
        ax.set_xlabel(r'$M$', fontsize=20)
        ax.plot(x, p_pi(x), linewidth=2)
        ax.grid()
        ifig += 1
        fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
        fig.suptitle(r' $S/S_{crit}=f(M)$', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$S/S_{crit}$', fontsize=20)
        ax.set_xlabel(r'$M$', fontsize=20)
        ax.plot(x, S_over_Scrit(x), linewidth=2)
        ax.grid()
        x = np.linspace(1., 6, n)
        ifig += 1
        fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
        fig.suptitle(r' $p_s/p_a=f(M)$ for normal shock waves', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$p_s/p_a$', fontsize=20)
        ax.set_xlabel(r'$M_s$', fontsize=20)
        ax.plot(x, P2_P1(x), linewidth=2)
        ax.grid()

        ifig += 1
        ratio = p_pi(downstream_Mach(x)) * pi2_pi1(x)
        ratio1 = p_pi(x) * P2_P1(x)
        fig = plt.figure(ifig, figsize=(Largeur, Hauteur))
        fig.suptitle(r' $p_s/p_{i_0}=f(M)$ for normal shock waves', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$p_s/p_{i_0}$', fontsize=20)
        ax.set_xlabel(r'$M_s$', fontsize=20)
        ax.plot(x, ratio, linewidth=2)
        ax.plot(x, ratio1, linewidth=2)
        ax.grid()

    set_question('4 : Application')
    Ms = inv_S_over_Scrit(section_ratio, Mach=2., show=False)
    print('Mach Ms            = %f ' % Ms)
    print('S / Sc  (verif.)   = %f ' % (S_over_Scrit(Ms)))
    ps_pi0 = p_pi(Ms)
    print('ps/pi0             = %f ' % ps_pi0)
    rpa.append(ps_pi0)
    print('list of ratios  p_a/pi_0', rpa)
    us_over_ac = u_over_ac(Ms)
    print('us/ac              = %f ' % us_over_ac)
    CF1 = gamma * p_pi(1.) *  us_over_ac  # qm us/(Sc pi0)
    print('qm us/(Sc pi0)     = %f ' % CF1)
    F1 = CF1 * pi0 * Sc
    print('F1                 = %f N' % F1)
    CF1bis = gamma * section_ratio *  us_over_ac ** 2 * rho_rhoi(Ms) * T_Ti(1.)
    print('qm us/(Sc pi0)     = %f ' % CF1bis)
    print('qm us/(Sc pi0)     = %f ' % (CF(Ms)))
    print('CF_m (Shapiro)     = %f ' % (CF_m(Ms, gamma=gamma)))
    print('CF_m (Shapiro) M=1 = %f ' % (CF_m(1., gamma=gamma)))
    print('CFT (M=1) pa=0     = %f ' % (CF_m(1., gamma=gamma) + p_pi(1)))
    print('den.(Shapirp),pa=0 = %f ' % (2 * pow(2 / (gamma + 1), 1 / (gamma - 1))))

    for r in rpa:
        print()
        print('pa/pi0             = %f' % r)
        Mopt = inverse_p_pi(r)
        print('Ms optimal         = %f' % Mopt)
        print('S/Sc optimal       = %f' % (S_over_Scrit(Mopt)))
        CF2 = CF_add(Ms, r)
        F2 = CF2 * pi0 * Sc
        print('CF2                = %f  ' % CF2)
        print('F2                 = %f N' % F2)
        print('F1+F2              = %f N' % (F1 + F2))
        print('CF1+CF2            = %f ' % (CF1 + CF2))
        print('ps/pa              = %f ' % (ps_pi0 / r))

    print('CF_max             = %f ' % (CFmax(ps_pi0)))
    print('F_max              = %f N' % (CFmax(ps_pi0) * pi0 * Sc))

    if plot_CF:
        table_r = [0.004, 0.005, 0.0073, 0.01]
        print('pi0/pa = ', [1 / r for r in list(table_r)])

        for leg in ['M', 'S', 'P', 'Pa', 'X']:
            ifig += 1
            plotCF(table_r, [[3.5, 5.5], [3.5, 5.5], [3.2, 5], [3., 4.5]], n, leg, ifig, gamma=1.4)

    if plot_CF or plot_ps:
        plt.show()

    set_question('5 : Vulcain engine (Ariane rocket)')
    ratioS = 58.5
    D = 2.15
    S = np.pi * D ** 2 / 4
    Sc = S / ratioS
    gam = 1.2
    pi0 = 115e5
    Ms = inv_S_over_Scrit(ratioS, Mach=5., show=False, gamma=gam)
    print('gamma              = %f ' % gam)
    print('Mach Ms            = %f ' % Ms)
    ps_pi0 = p_pi(Ms, gamma=gam)
    print('ps/pi0             = %f ' % ps_pi0)
    print('ps                 = %f ' % (pi0 * ps_pi0))
    print('CF_max             = %f ' % (CFmax(ps_pi0, gamma=gam)))
    print('F_max              = %f N' % (CFmax(ps_pi0, gamma=gam) * pi0 * Sc))

    p0 = 101325
    altitude = np.linspace(13.7, 13.71, 11)
    for alt in altitude:
        sigma, delta, theta = Atmosphere(alt)
        print(alt, delta * p0)


def Exercice10_2():
    """ 
    Rocket nozzle thrust, Mattingly's formules
    """
    verb = False  # verbose
    plot_figure = True
    set_title(" Rocket nozzle thrust, Mattingly's formules")

    # Tc = combustion chamber temperature,   Ti0 for us
    # Pc = combustion chamber  pressure,   Pi0 for us

    gam = 1.2  # gamma is set
    m = (gam + 1) / (gam - 1)  # power coefficient used later  
    R = 283  # combustion gas constant  
    Tc = 500  # combustion chamber temperature in K  
    n = 21
    Pc_over_Pe = np.linspace(5, 500, n)
    Pe_over_Pc = 1 / Pc_over_Pe
    Gamma = np.sqrt(gam * pow((gam + 1) / 2, -m))  # useful coefficient
    print('gam                = %f ' % gam)
    print('Gamma              = %f ' % Gamma)

    Cstar = np.sqrt(R * Tc) / Gamma
    print('C^*                = %f m/s' % Cstar)  # characteristic velocity

    # *****************************
    # fig. 3.28 p. 199
    # *****************************
    def exit_velocity(r, Tc):
        """
        ejection velocity of gas
        as a function of r= Pe/Pc
        """
        tmp, ac = 1 - pow(r, (gam - 1) / gam), gam * R * Tc
        Vmax = np.sqrt(2 / (gam - 1) * ac)
        Ve = np.sqrt(tmp) * Vmax
        return Ve, Vmax

    Ve, Vmax = exit_velocity(Pe_over_Pc, Tc)
    if verb:
        print('Vmax               = %f m/s' % Vmax)
        print('Ve = ', Ve / Vmax)
    if plot_figure:
        simple_plot( '$Ve$ réduit', Pc_over_Pe, Ve / Vmax, [r'$P_c/P_e$', r'$Ve/Vmax$'])

    # *****************************
    # fig. 3.30 p. 202
    # *****************************
    def ratio_Ae_over_At(r):
        """
        Ae/At as a function of  r=Pe/Pc
        """
        tmp = 2 * gam / (gam - 1) * (pow(r, 2 / gam) - pow(r, 1. + 1. / gam))
        return Gamma / np.sqrt(tmp)

    Ae_over_At = ratio_Ae_over_At(Pe_over_Pc)
    if verb:
        print('Ae/At = ', Ae_over_At)
    if plot_figure:
        simple_plot( 'Section ratio', Pc_over_Pe, Ae_over_At, [r'$P_c/P_e$', r'$\varepsilon=Ae/At$'])

    def CFi_opt(r):
        """
        optimal ideal propulsion coefficient,  r=Pe/Pc
        """
        tmp = 2 * gam / (gam - 1) * (1 - pow(r, gam - 1 / gam))
        return Gamma * np.sqrt(tmp)

    CFiMax = CFi_opt(Pe_over_Pc)
    if verb:
        print('CFi_opt = ', CFiMax)
    if plot_figure:
        simple_plot( r'$CF_{opt}$', Pc_over_Pe, CFiMax, [r'$P_c/P_e$', r'$C_{F_\max}$'])

    def CFi_vac(r):
        """
        CFi in the vacuum, r=Pe/Pc
        """
        return CFi_opt(r) + r * ratio_Ae_over_At(r)

    CFi_vacuum = CFi_vac(Pe_over_Pc)
    if verb:
        print('CFi vac = ', CFi_vacuum)
    if plot_figure:
        simple_plot( r'$Cf_{vac}$', Pc_over_Pe, CFi_vacuum, [r'$P_c/P_e$', r'$C_{F_{vac}}$'])

    def CFi_opt2(r, eps):
        """
        optimal ideal propulsion coefficient (second expression)
        r=Pe/Pc, eps=Ae/At
        """
        return 2 * gam / (gam - 1) * eps * (1 - pow(r, 1 - 1 / gam)) * pow(r, 1 / gam)

    tabx, taby = [], []
    epsilon = [5, 10, 30]
    for eps in epsilon:
        print('eps = ', eps)
        CFiMax1 = CFi_opt2(Pe_over_Pc, eps)
        if verb:
            print('CFi_opt = ', CFiMax1)
        taby.append(CFiMax1)
        taby.append(CFiMax)
        tabx.append(Pc_over_Pe)
        tabx.append(Pc_over_Pe)

        if plot_figure:
            simple_plot( r'$CF2_{opt}$ for $\varepsilon$= %4.1f' % eps, tabx, taby,
                        [r'$P_c/P_e$', r'$C_{F_\max}$'], len(tabx))

    # simple_plot(4,r'$CF2_{opt}$ pour $\varepsilon$= %4.1f'%(eps),[Pc_over_Pe,Pc_over_Pe],[CFiMax,CFiMax1],[r'$P_c/P_e$',r'$C_{F_\max}$'])

    def CFi(r, ra):
        """
        CFi with respect to r=Pe/Pc, eps=Ae/At et de ra=Pe/Pc,
        not used
        """
        return CFi_opt(r) + (r - ra) * ratio_Ae_over_At(r)

    plt.show()

def Exercice10_3():
    """ 
    Pitot tube into a supersonic nozzle
    """
    # index s means 'sortie' ie exit in english

    set_title("Pitot tube into a supersonic nozzle")
    pi0 = 10.       # pressure chamber in bar
    Ti0 = 500.      # temperature chamber in K
    pi2 = 3.283     # total pressure given by the Pitot tube in  bar.
    gam = 1.4       # Cp/Cv
    S_c = 0.2       # throat section in m^2

    set_question('1 : exit Mach')

    # My departure is the solution
    Ms = 3.0
    print("pi2/pi0              = %f" % (pi2_pi1(3, gamma=gam)))
    # Then I recalculate Ms
    Ms_calcul = inv_pi2_over_pi1(pi2 / pi0, Mach=2., show=False, gamma=gam)
    print("Ms solved (verif)   = %f" % Ms_calcul)

    set_question('2 : M2 and p2')
    M2 = downstream_Mach(Ms, gamma=gam)
    p2 = pi2 * p_pi(M2, gamma=gam)
    print("M2                   = %f" % M2)
    print("p2                   = %f bar" % p2)

    set_question('3 : exit section')

    # ratio S_s/S_c 
    r = S_over_Scrit(Ms, gamma=gam)
    print("p2                   = %f bar" % p2)
    print("S_s/S_c              = %f " % r)
    S_s = S_c * r
    print("S_s                  = %f m^2 " % S_s)

    set_question('4 : ps, Ts, us')
    ps = pi0 * p_pi(Ms)
    print("ps                   = %f bar" % ps)
    Ts = Ti0 * T_Ti(Ms)
    print("Ts                   = %f K" % Ts)
    us = Ms * np.sqrt(gam * r_Air * Ts)
    print("us                   = %f m/s" % us)


def Exercice10_4():
    """
    Oxygen - hydrogen rocket engine
    """
    set_title(" Oxygen - hydrogen rocket engine")
    pi0 = 15e5      # total pressure given by the Pitot tube in  Pa
    Ti0 = 4000.     # temperature chamber in K
    ps = 1174       # exit pressure
    gam = 1.22      # Cp/Cv
    r = 519.6       # r= Cp-Cv
    Ds = 3.0        # exit section diameter

    Ms = inverse_p_pi(ps / pi0, gamma=gam)
    print("Ms                   = %f" % Ms)
    Ts = Ti0 * T_Ti(Ms, gamma=gam)
    print("Ts                   = %f K" % Ts)
    section_ratio = S_over_Scrit(Ms, gamma=gam)
    print("S_s/S_c              = %f" % section_ratio)
    S_s = np.pi * Ds ** 2 / 4
    print("S_s                  = %f m^2" % S_s)
    qm = qm_1D(S_s, Ms, pi0, Ti0, gamma=gam, r=r)
    print("qm                   = %f kg/s" % qm)
    a_s = np.sqrt(gam * r * Ts)
    us = Ms * a_s
    Fm = qm * us
    S_c = S_s / section_ratio
    C_F = Fm / (pi0 * S_c)
    print("as                   = %f m/s" % a_s)
    print("us                   = %f m/s" % us)
    print("S_c                  = %f m^2" % S_c)
    print("F_m                  = %f N" % Fm)
    print("C_F                  = %f " % C_F)


def Exercice10_5():
    """
    normal shock wave example
    """
    set_title("normal shock wave example")
    M0 = 3
    p0 = 1e5
    rho0 = 1.25
    set_question("1 - upstream quantities")
    T0 = p0 / (rho0 * r_Air)
    Ti0 = T0 / T_Ti(M0)
    pi0 = p0 / p_pi(M0)
    rhoi0 = pi0 / (r_Air * Ti0)
    a0 = sound_velocity(T0)
    u0 = M0 * a0
    print("Ti0 = %f K, \t pi0 = %e Pa,  rhoi0 = %f kg/m^3" % (Ti0, pi0, rhoi0))
    print("T0 = %f K, \t  a0 = %f m/s, \t u0 = %f m/s" % (T0, a0, u0))
    set_question("2 - normal shock ; downstream quantities")
    M1 = downstream_Mach(M0)
    pressure_ratio = P2_P1(M0)
    density_ratio = rho2_rho1(M0)
    isentropic_pressure_ratio = pi2_pi1(M0)
    rho1 = density_ratio * rho0
    p1 = pressure_ratio * p0
    pi1 = isentropic_pressure_ratio * pi0
    T1 = p1 / (rho1 * r_Air)
    a1 = sound_velocity(T1)
    u1 = a1 * M1
    print("p0  = %e" % p0)
    print("M1 = %f, \t\t  p1/p0 = %f, \t rho1/rho0 = %f" % (M1, pressure_ratio, density_ratio))
    print("T1 = %f K, \t p1 = %e Pa, \t u1 = %f m/s " % (T1, p1, u1))
    print("pi1/pi0 = %f , \t pi1 = %e Pa " % (isentropic_pressure_ratio, pi1))
    set_question("3 - entropy variation")
    print("s1-s0 = ", -r_Air * np.log(isentropic_pressure_ratio))
