#!/bin/py
# -*- coding: utf-8 -*-
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2018 - 2020*

*license : AGPL-3.0*

Correction of the chapter 15 exercises : Turbulent flow (with Turbulence class)
    ..
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.special import erfc, erf
from scipy.integrate import odeint
import Book.Turbulence.Turbulence as Tu
from Tools.misc import set_title, set_question

# from Book.Turbulence.misc            import

width = 10
height = 8


def Exercice15_11():
    """ 
    Turbulent analytical velocity profile with suction
    """
    set_title("Turbulent analytical velocity profile with suction")
    ifig = 0
    display = True
    display_perso = True
    displayLm = False
    option_semilog = False  # only to plot mixing length
    etude_aplus = False     # study of the A+ law wrt v0+
    cas = 1                 # for the exercise book set cas=1

    Profil = []
    if cas == 0:
        L = [0, 1, 2, 3, 4]                 # first testcase
    elif cas == 1:
        L = [6, 7, 8, 9, 10, 11, 12, 13]    # figures of the last exercise of the book
    else:
        L = [5, 6, 7, 8]                    # an example, user dependent

    Ncas = len(L)
    m = 0
    for k in L:
        set = Tu.set_parameters()
        if k == 0:  # Laminar flow
            set["name"] = 'LAMINAR CASE'
            set["ypMax"] = 6
            set["v0p"] = 0.01

        if k == 1:  # approximated viscous sublayer
            set["name"] = 'APPROXIMATED VISCOUS LAYER'
            set["Model"] = "SCVonly"  # only true if  v0p = 0
            set["ypMax"] = 5
            set["v0p"] = 0.0

        if k == 2:  # approximated internal layer 
            set["name"] = 'APPROXIMATED INTERNAL LAYER'
            set["Model"] = "CIonly"  # only true if  v0p = 0
            set["v0p"] = 0.0
            set['Damping'] = "VanDriest"

        if k == 3:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.0

        if k == 4:  # Michel's mixing length model
            set["name"] = "MICHEL's MODEL WITHOUT CORRECTION"
            set["Model"] = "Michel"
            set["v0p"] = 0

        if k == 5:  # Michel's mixing length model
            set["name"] = "MICHEL's MODEL WITHOUT CORRECTION"
            set["Model"] = "Michel"
            set["v0p"] = 0
            set['Damping'] = "VanDriest"

        if k == 6:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.0
            # set["Aplus"]  = 0
            # set["Bplus"]  = 0
            set['Damping'] = "VanDriest"

        if k == 7:  # Log law
            set["name"] = 'LOG LAW MODEL'
            set["Model"] = "loi_log"
            set["A"] = 1 / set["kappa"]
            set["B"] = 5.28

        if k == 8:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.03
            set["Aplus"] = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        if k == 9:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.03
            set["Aplus"] = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Cebeci1"

        if k == 10:  # linear mode
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = 0.
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        if k == 11:  # linear mode
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = 0.0
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Cebeci1"

        if k == 12:  # linear mode
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = -0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        if k == 13:  # linear mode
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = 0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].get_bl_profile()
        Profil[m].calculate_mixing_length()
        print(Profil[m].Modele, Profil[m].v0p)
        m += 1

    if displayLm:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('Mixing length', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if option_semilog:
                ax.semilogx(Profil[k].yp, Profil[k].Lplus, label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                    Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            else:
                ax.plot(Profil[k].yp, Profil[k].Lplus, label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                    Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        if option_semilog:
            k = 0
            ax.semilogx(Profil[k].yp, Profil[k].Lplus, 'bo', label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        else:
            k = 0
            ax.plot(Profil[k].yp, Profil[k].Lplus, 'bo', label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        ax.axis([1, 200, 0, 100])
        ax.legend(loc='upper left')
        ifig += 1

    if display:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('Turbulent flow 1', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        # fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if Profil[k].com != "wrong":
                ax.semilogx(Profil[k].yp, Profil[k].up, label=r'$v_0^+$ = %1.4f, mod.  %s, D: %s, $A^+$ = %4.2f, %s' % (
                    Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus, Profil[k].com))
        ax.axis([1, 200, 0, 25])
        ax.legend(loc='upper left')
        ifig += 1

    if display_perso and cas == 1:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('Turbulent flow 2', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        k = 0
        ax.semilogx(Profil[k].yp, Profil[k].up, 's', markevery=(0, 5),
                    label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                        Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        k = 1
        ax.semilogx(Profil[k].yp, Profil[k].up, 'o', markevery=(0, 7),
                    label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                        Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        k = 2
        ax.semilogx(Profil[k].yp, Profil[k].up, '-o', label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
            Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        k = 3
        ax.semilogx(Profil[k].yp, Profil[k].up, '*', label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
            Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        k = 4
        ax.semilogx(Profil[k].yp, Profil[k].up, '-o', markevery=(0, 11),
                    label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                        Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        k = 5
        ax.semilogx(Profil[k].yp, Profil[k].up, '--', label=r'$v_0^+$ = %1.4f, mod.  %s , D : %s, $A^+$ = %4.2f' % (
            Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        ax.axis([1, 200, 0, 30])
        ax.legend(loc='upper left')

    if etude_aplus:
        Profil[-1].ifig = ifig + 1
        print('A+ influence on the transpiration velocity and  on the pressure gradient')
        Profil[-1].Aplus = 26  # A+ given by  Kays's law ( A+ != to 26 for p_0^+=0 and v_0^+=0 )
        v0p = [-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3]
        pp = np.linspace(-0.039, 0.06, 51)
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays1')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays2')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays3')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Cebeci1')
        Profil[-1].plot_B_Kays(v0p, pp, mode="Linear")
        Profil[-1].plot_Constante_A_B()
        Profil[-1].A_B_regression()

    if display or display_perso or displayLm or etude_aplus:
        plt.show()


def Exercice15_12():
    """ 
    Turbulent analytical velocity profile with the power law
    """
    set_title("Turbulent analytical velocity profile with the power law")

    display_perso = True
    option_semilog = False

    Profil = []
    L = [0]
    m = 0
    for k in L:
        set = Tu.set_parameters()
        if k == 0:
            set["name"] = 'Power law, Rtau and Cu, and N given'
            set["Model"] = "PowerLaw"
            set["ypMax"] = 1
            set["v0p"] = 0
            set["Rtau"] = 2000
            set["Cu"] = 8.90
            set["Cf"] = 0
            set["ShowParameters"] = True

        if k == 1:
            set["name"] = 'Power law, Rtau and Rex and N given'
            set["Model"] = "PowerLaw"
            set["ypMax"] = 1
            set["v0p"] = 0
            set["Rtau"] = 2000
            set["Cu"] = 0
            set["Rex"] = 2.85e6
            set["N"] = 9
            set["Cf"] = 0.00309
            set["ShowParameters"] = True

        if k == 2:
            set["name"] = 'Power law, Rtau,Rex and Cu given, solve N'
            set["Model"] = "PowerLaw"
            set["ypMax"] = 1
            set["v0p"] = 0
            set["Rtau"] = 2000
            set["Cu"] = 8.75
            set["Rex"] = 2.85e6
            set["Cf"] = 0
            set["N"] = 0.
            set["ShowParameters"] = True

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].get_bl_profile()
        Profil[m].plot_profile(mode="Linear")
        m += 1

    set_question("Last question")
    set = Tu.set_parameters()
    set["name"] = "Power law, Rtau and Cu inconnus, loi de Michel's law"
    set["Model"] = "PowerLaw"
    set["ypMax"] = 1
    set["v0p"] = 0
    set["Rtau"] = 0
    set["Cu"] = 0
    set["Cf"] = 0.00309
    set["N"] = 9
    set["Rex"] = 2.84959e6
    set["ShowParameters"] = True
    Profil.append(Tu.ProfilTurbulent(set))
    Profil[0].set_Cu_from_Cf()
    print('Cu =', Profil[0].Cu)
    Profil[0].calcul_Rtau()
    print('Rtau =', Profil[0].Rtau)

    set_question("law plots")
    RexMin, RexMax, nR = 1e6, 1e7, 51
    Rex = np.linspace(RexMin, RexMax, nR)
    print(Profil[0].SchultzGrunow(Rex))
    print(Profil[0].Michel(Rex))
    print(Profil[0].Cf_Power_Law(Rex))
    for Prof in Profil:
        print("Name = %s , \t N = %i" % (Prof.name, Prof.N))
    Rex_sol = Profil[0].Rex
    print('Rex solution = ', Rex_sol)

    if display_perso:
        CfMin, CfMax = 0.0024, 0.0036
        fig = plt.figure(figsize=(10, 8))
        fig.suptitle('skin-friction law', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$Rex$', fontsize=20)
        ax.set_ylabel(r'$Cf$', fontsize=20)
        if option_semilog:
            ax.semilogx(Rex, Profil[0].SchultzGrunow(Rex), label='Schultz-Grunow')
            ax.semilogx(Rex, Profil[0].Michel(Rex), label='Michel')
            ax.semilogx(Rex, Profil[0].Cf_Power_Law(Rex), label='présent')
            ax.semilogx([Rex_sol, Rex_sol], [CfMin, CfMax], label='ref')
        else:
            ax.plot(Rex, Profil[0].SchultzGrunow(Rex), label='Schultz-Grunow')
            ax.plot(Rex, Profil[0].Michel(Rex), label='Michel')
            ax.plot(Rex, Profil[0].Cf_Power_Law(Rex), label='present')
            ax.plot([Rex_sol, Rex_sol], [CfMin, CfMax], label='ref')

        ax.axis([RexMin, RexMax, CfMin, CfMax])
        ax.grid()
        ax.legend(loc='upper right')
        plt.show()


def Exercice15_13():
    """ 
    Turbulent analytical velocity profile with the  Van Driest's correction
    """
    set_title("Turbulent analytical velocity profile with the  Van Driest's correction")

    plot_CorrectionVanDriest = False
    display_profile = True
    damping_law_analysis = False
    ifig = 0
    Profil = []
    m = 0  # initialization of the testcase number 
    if plot_CorrectionVanDriest:
        set_question("0 - plot of Van Driest's correction")
        set = Tu.set_parameters()
        set["name"] = "Van Driest's correction"
        set["Model"] = "Linear"
        set["ypMax"] = 140
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[0].plot_loi_Van_Driest(xaxis="Linear")  # ou "Linear" (linear or semilog scale)
        m += 1

    set_question("6- Model comparisons")
    Profil = []
    m = 0
    L = [1, 2, 3, 4]  # list of the testcase
    Ncas = len(L)

    for k in L:
        set = Tu.set_parameters()
        set["ypMax"] = 250
        set["Rtau"] = 10000
        set["ShowParameters"] = True
        if k == 0:
            set["name"] = 'LINEAR MODEL, without Van Driest correction'
            set["Model"] = "Linear"
            set['Damping'] = "No"
        if k == 1:
            set["name"] = 'LINEAR MODEL, with Van Driest correction'
            set["Model"] = "Linear"
            set['Damping'] = "VanDriest"
        if k == 2:  # Loi Log
            set["name"] = 'with LOG LAW MODEL'
            set["Model"] = "loi_log"
            set['Damping'] = "No"
            set["ypMin"] = 1.
            set["A"] = 1 / set["kappa"]
            set["B"] = 5.28
        if k == 3:  # u+=y+
            set["name"] = 'LINEAR VISCOUS SUBLAYER'
            set["Model"] = "SCVlinear"
            set["ypMax"] = 11.5
            set['Damping'] = "No"
        if k == 4:  # Michel's mixing length
            set["name"] = "MICHEL's MODEL with CORRECTION"
            set["Model"] = "Michel"
            set['Damping'] = "VanDriest"

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].get_bl_profile()
        # Profil[m].plot_profile(mode="Log")
        m += 1

    if display_profile:
        fig = plt.figure(ifig, figsize=(10, 8))
        fig.suptitle('turbulent flow', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        for k in range(Ncas):
            ax.semilogx(Profil[k].yp, Profil[k].up, label=r'mod.  %s, D : : %s' % (Profil[k].Modele, Profil[k].Damping))
        ax.axis([0.1, 300, 0, 20])
        ax.legend(loc='upper left')
        plt.show()
        ifig += 1

    if damping_law_analysis:  # damping laws analysis
        Profil[-1].ifig = ifig
        print('A+ influence on the transpiration velocity and pressure gradient')
        Profil[-1].Aplus = 26  # A+ given Kays ( A+ != 26 for p_0^+=0 and v_0^+=0 )
        v0p = [-0.1, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3]
        pp = np.linspace(-0.039, 0.06, 21)
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays1')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays2')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays3')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Cebeci1')
        Profil[-1].plot_B_Kays(v0p, pp, mode="Linear")
        Profil[-1].plot_Constante_A_B()
        plt.show()
        Profil[-1].A_B_regression()


def Exercice15_14():
    """ 
    Turbulent analytical velocity profile with the  Van Driest's correction
    plot of laminar and turbulent tau (stress)
    """
    set_title("Turbulent analytical velocity profile with the  Van Driest's correction")

    plot_CorrectionVanDriest = False
    display_profile = True
    ifig = 0
    Profil = []
    m = 0
    if plot_CorrectionVanDriest:
        set_question("0 - plot of Van Driest's correction")
        set = Tu.set_parameters()
        set["name"] = "Van Driest's correction"
        set["Model"] = "Linear"
        set["ypMax"] = 140
        Profil.append(Tu.ProfilTurbulent(set))
        Profil[0].plot_loi_Van_Driest(xaxis="Linear")  # ou "Linear" (linear or semilog scale)
        m += 1

    set_question("6- Model comparison")
    Profil = []
    m = 0
    L = [1, 2, 3, 4]
    Ncas = len(L)

    for k in L:
        set = Tu.set_parameters()
        set["ypMax"] = 250
        set["Rtau"] = 10000
        set["ShowParameters"] = True
        set["n"] = 501

        if k == 0:
            set["name"] = 'LINEAR MODEL, without Van Driest correction'
            set["Model"] = "Linear"
            set['Damping'] = "No"
        if k == 1:
            set["name"] = 'LINEAR MODEL, with Van Driest correction'
            set["Model"] = "Linear"
            set['Damping'] = "VanDriest"
        if k == 2:  # Loi Log
            set["name"] = 'with LOG LAW MODEL'
            set["Model"] = "loi_log"
            set['Damping'] = "No"
            set["ypMin"] = 1.
            set["A"] = 1 / set["kappa"]
            set["B"] = 5.28
        if k == 3:  # u+=y+
            set["name"] = 'LINEAR VISCOUS SUBLAYER'
            set["Model"] = "SCVlinear"
            set["ypMax"] = 11.5
            set['Damping'] = "No"
        if k == 4:  # Michel's mixing length
            set["name"] = "MICHEL's MODEL with CORRECTION"
            set["Model"] = "Michel"
            set['Damping'] = "VanDriest"
            set["ypMax"] = set["Rtau"]

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].get_bl_profile()
        # Profil[m].plot_profile(mode="Log")
        m += 1

    k = Ncas - 1
    print(Profil[k].dup_dyp.shape)

    if display_profile:
        fig = plt.figure(ifig, figsize=(10, 8))
        fig.suptitle('turbulent flow', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        for k in range(Ncas):
            ax.semilogx(Profil[k].yp, Profil[k].up, label=r'mod.  %s, D : : %s' % (Profil[k].Modele, Profil[k].Damping))
        # ax.axis([0.1,300,0,20])
        # ax.axis([0.1,1000,0,100])
        ax.legend(loc='upper left')
        plt.show()
        ifig += 1

        fig = plt.figure(ifig, figsize=(10, 8))
        fig.suptitle('Turbulent flow', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$y/\delta$', fontsize=20)
        ax.set_ylabel(r'$u/Ue$', fontsize=20)
        ax.grid()
        Rtau = set["Rtau"]
        Uep = Profil[-1].Uep

        for k in range(Ncas):
            ax.plot(Profil[k].yp / Rtau, Profil[k].up / Uep,
                    label=r'mod.  %s, D : : %s' % (Profil[k].Modele, Profil[k].Damping))
        # ax.axis([0.1,300,0,20])
        # ax.axis([0.1,1000,0,100])
        ax.legend(loc='upper left')
        plt.show()
        ifig += 1

        fig = plt.figure(ifig, figsize=(10, 8))
        fig.suptitle('Turbulent flow', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_ylabel(r'$y^+$', fontsize=20)
        ax.set_xlabel(r'$\dfrac{du^+}{dy^+}$', fontsize=20)
        ax.grid()
        Rtau = set["Rtau"]
        Uep = Profil[-1].Uep

        k = Ncas - 1
        ax.semilogy(Profil[k].dup_dyp, Profil[k].yp,
                    label=r'mod.  %s, D : : %s' % (Profil[k].Modele, Profil[k].Damping))
        # ax.plot(Profil[k].dup_dyp,Profil[k].yp,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        # k=0
        # ax.plot(Profil[k].yp,Profil[k].dup_dyp,label=r'mod.  %s, D : : %s'%(Profil[k].Modele,Profil[k].Damping))
        # ax.axis([0.1,300,0,1])
        ax.axis([0, 1, 0.1, 1000])
        ax.legend(loc='upper left')
        plt.show()
        ifig += 1


def Exercice15_15():
    """ 
    Study of turbulence profile with ou without suction
    """
    set_title("turbulence profile with ou without suction")
    ifig = 0
    display = True
    display_perso = False
    displayLm = False
    option_semilog = False
    etude_aplus = False     # study of the  A+  law w.r.t. v0⁺
    cas = 0
    L = [1]
    Profil = []
    if cas == 0:
        L = [5, 6, 7, 8]
    elif cas == 1:
        L = [6, 7, 8, 9, 10, 11, 12, 13]

    Ncas = len(L)
    m = 0
    for k in L:
        set = Tu.set_parameters()
        set["ShowParameters"] = True
        if k == 5:  # Michel's model
            set["name"] = "MICHEL's MODEL with CORRECTION"
            set["Model"] = "Michel"
            set["v0p"] = -0.03
            set['Damping'] = "VanDriest"
            set["ypMax"] = set["Rtau"]

        if k == 6:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.03
            # set["Aplus"] = 0
            # set["Bplus"] = 0
            set['Damping'] = "VanDriest"
            set["ypMax"] = set["Rtau"]

        if k == 7:  # Loi Log
            set["name"] = 'MODEL WITH LOG LAW'
            set["Model"] = "loi_log"
            set["A"] = 1 / set["kappa"]
            set["B"] = 5.28

        if k == 8:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.03
            set["Aplus"] = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"
            set["ypMax"] = set["Rtau"]

        if k == 9:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = -0.03
            set["Aplus"] = 26
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Cebeci1"

        if k == 10:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = 0.
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        if k == 11:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = 0.0
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Cebeci1"

        if k == 12:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = -0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        if k == 13:  # linear model
            set["name"] = 'LINEAR MODEL'
            set["Model"] = "Linear"
            set["v0p"] = 0.03
            set["Aplus"] = 26
            set["Pplus"] = 0.03
            set["AplusOpt"] = 2
            set['Damping'] = "VanDriest"
            set['AplusModel'] = "Kays1"

        Profil.append(Tu.ProfilTurbulent(set))
        Profil[m].get_bl_profile()
        Profil[m].calculate_mixing_length()
        print(Profil[m].Modele, Profil[m].v0p)
        m += 1

    if displayLm:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('Mixing length', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.80)
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if option_semilog:
                ax.semilogx(Profil[k].yp, Profil[k].Lplus, label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                    Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            else:
                ax.plot(Profil[k].yp, Profil[k].Lplus, label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                    Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        if option_semilog:
            k = 0
            ax.semilogx(Profil[k].yp, Profil[k].Lplus, 'bo', label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        else:
            k = 0
            ax.plot(Profil[k].yp, Profil[k].Lplus, 'bo', label=r'mod.  %s, D : %s, $A^+$ = %4.2f' % (
                Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
        ax.axis([1, 200, 0, 100])
        ax.legend(loc='upper left')
        ifig += 1

    if display:
        fig = plt.figure(ifig, figsize=(width, height))
        fig.suptitle('Turbulent flow 1', fontsize=14, fontweight='bold')
        ax = fig.add_subplot(111)
        # fig.subplots_adjust(top=0.80) 
        ax.set_xlabel(r'$y^+$', fontsize=20)
        ax.set_ylabel(r'$u^+$', fontsize=20)
        ax.grid()
        for k in range(Ncas):
            if Profil[k].com != "wrong":
                ax.semilogx(Profil[k].yp, Profil[k].up, label=r'$v_0^+$ = %1.4f, mod.  %s, D: %s, $A^+$ = %4.2f, %s' % (
                    Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus, Profil[k].com))
                if Profil[k].plot_sv:
                    ax.semilogx(Profil[k].y_vp, Profil[k].u_vp, label=r'$u_v^+$')
                    ax.semilogx(Profil[k].yp_sv, Profil[k].up_sv, "ko", markersize=7)
        ax.axis([1, set["Rtau"], 0, 37])
        ax.legend(loc='upper left')
        ifig += 1

    if display_perso:
        if len(Profil) >= 6:
            fig = plt.figure(ifig, figsize=(width, height))
            fig.suptitle('Ecoulement turbulent 2', fontsize=14, fontweight='bold')
            ax = fig.add_subplot(111)
            fig.subplots_adjust(top=0.80)
            ax.set_xlabel(r'$y^+$', fontsize=20)
            ax.set_ylabel(r'$u^+$', fontsize=20)
            ax.grid()
            k = 0
            ax.semilogx(Profil[k].yp, Profil[k].up, 's', markevery=(0, 5),
                        label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                            Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            k = 1
            ax.semilogx(Profil[k].yp, Profil[k].up, 'o', markevery=(0, 7),
                        label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                            Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            k = 2
            ax.semilogx(Profil[k].yp, Profil[k].up, '-o', label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            k = 3
            ax.semilogx(Profil[k].yp, Profil[k].up, '*', label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            k = 4
            ax.semilogx(Profil[k].yp, Profil[k].up, '-o', markevery=(0, 11),
                        label=r'$v_0^+$ = %1.4f, mod.  %s, D : %s, $A^+$ = %4.2f' % (
                            Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            k = 5
            ax.semilogx(Profil[k].yp, Profil[k].up, '--', label=r'$v_0^+$ = %1.4f, mod.  %s , D : %s, $A^+$ = %4.2f' % (
                Profil[k].v0p, Profil[k].Modele, Profil[k].Damping, Profil[k].Aplus))
            ax.axis([1, 200, 0, 30])
            ax.legend(loc='upper left')
        else:
            print('len (Profil)      :', len(Profil))

    if etude_aplus:
        Profil[-1].ifig = ifig + 1
        print('A+ influence on the transpiration velocity and  on the pressure gradient')
        Profil[-1].Aplus = 26  # A+ given by  Kays's law ( A+ != to 26 for p_0^+=0 and v_0^+=0 )
        v0p = [-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.3]
        pp = np.linspace(-0.039, 0.06, 51)
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays1')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays2')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Kays3')
        Profil[-1].plot_A_vanDriest(v0p, pp, mode="Linear", modele='Cebeci1')
        Profil[-1].plot_B_Kays(v0p, pp, mode="Linear")
        Profil[-1].plot_Constante_A_B()
        Profil[-1].A_B_regression()

    if display or display_perso or displayLm or etude_aplus:
        plt.show()
