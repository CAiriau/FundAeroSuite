# !/bin/python
"""
Fundamental Aerodynamic Suite, FundAeroSuite

.. codeauthor:: C. Airiau

*date : 2017 - 2020*

*license : AGPL-3.0*

Python module which makes the relationship between the exercises number, their names and the pyhton function names
    ..

The list of the exercises is also given in French and English

"""

import os

from Book.Chapter2 import Exercice2_1, Exercice2_2
from Book.Chapter3 import Exercice3_1, Exercice3_2, Exercice3_3, Exercice3_4, Exercice3_5, Exercice3_20
from Book.Chapter3 import Exercice3_6, Exercice3_7, Exercice3_7, Exercice3_8, Exercice3_9, Exercice3_12
from Book.Chapter4 import Exercice4_1, Exercice4_2, Exercice4_3, Exercice4_4, Exercice4_5
from Book.Chapter4 import Exercice4_6, Exercice4_7, Exercice4_8
from Book.Chapter5 import Exercice5_1, Exercice5_2, Exercice5_3, Exercice5_4, Exercice5_5
from Book.Chapter5 import Exercice5_6, Exercice5_7, Exercice5_8, Exercice5_9, Exercice5_10
from Book.Chapter5 import Exercice5_11, Exercice5_12, Exercice5_13, Exercice5_14, Exercice5_15
from Book.Chapter6 import Exercice6_1, Exercice6_3, Exercice6_4
from Book.Chapter7 import Exercice7_0, Exercice7_1, Exercice7_2, Exercice7_3, Exercice7_4, Exercice7_10
from Book.Chapter7c import Exercice7_5
from Book.Chapter8 import Exercice8_1, Exercice8_2, Exercice8_3, Exercice8_4, Exercice8_5
from Book.Chapter8 import Exercice8_6, Exercice8_7
from Book.Chapter9 import Exercice9_1, Exercice9_2, Exercice9_3, Exercice9_4, Exercice9_5
from Book.Chapter10 import Exercice10_1, Exercice10_2, Exercice10_3, Exercice10_4, Exercice10_5
from Book.Chapter11 import Exercice11_1, Exercice11_2, Exercice11_3
from Book.Chapter12 import Exercice12_1, Exercice12_2, Exercice12_4
from Book.Chapter12b import Exercice12_3
from Book.Chapter13 import Exercice13_1, Exercice13_2, Exercice13_3, Exercice13_4, Exercice13_5
from Book.Chapter14 import Exercice14_0, Exercice14_1, Exercice14_3, Exercice14_4, Exercice14_5
from Book.Chapter14 import Exercice14_6, Exercice14_100, Exercice14_102, Exercice14_101
from Book.Chapter15_Laminaire import Exercice15_1, Exercice15_2, Exercice15_4, Exercice15_6
from Book.Chapter15_Laminaire import Exercice15_7, Exercice15_8
from Book.Chapter15_Transition import Exercice15_9
from Book.Chapter15_Turbulent import Exercice15_10
from Book.Chapter15_Tu import Exercice15_11, Exercice15_12, Exercice15_13
from Book.Chapter15_Tu import Exercice15_14, Exercice15_15
from Book.misc import Exercice7_30


# from conical_shock.conical_shock  import


def liste_chapters(lang="fr"):
    """
    list of the chapters, in french and english

    Args:
        lang (string): langage "fr" or "en"
    """
    liste_chapter = []
    if lang == "fr":
        liste_chapter.append("Equations de conservation et d'évolution")
        liste_chapter.append("Modèles mathématiques pour l'aérodynamique")
        liste_chapter.append("Théorie des potentiels complexes: application aux profils d'aile")
        liste_chapter.append("Théorie linéarisée des profils minces")
        liste_chapter.append("Théorie de la ligne portante de Lanchester-Prandtl")
        liste_chapter.append("Théorie de la surface portante et des corps élancés")
        liste_chapter.append("Aspects numériques : méthodes des singularités")
        liste_chapter.append("Ecoulements compressibles subsonique et transsonique")
        liste_chapter.append("Ecoulements supersoniques linéarisés")
        liste_chapter.append("Ecoulements compressibles monodimensionnels")
        liste_chapter.append("Ecoulements supersoniques bidimensionnels")
        liste_chapter.append("Méthodes des caractéristiques en régimes stationnaire et instationnaire")
        liste_chapter.append("Ecoulement supersonique: corps élancés et aile")
        liste_chapter.append("Ecoulements hypersoniques")
        liste_chapter.append("Effets visqueux et couche limite")
    else:
        liste_chapter.append("General equations of conservation and evolution")
        liste_chapter.append("Mathematical models for aerodynamics")
        liste_chapter.append("Complex potential theory : application to airfoils")
        liste_chapter.append("Linearized theory of thin airfoils")
        liste_chapter.append("Lanchester-Prandtl lifting line theory")
        liste_chapter.append("Lifting surface theory and slender bodies")
        liste_chapter.append("Numerical aspects : méthod of singularities")
        liste_chapter.append("Subsonic and transonic compressible flows")
        liste_chapter.append("Linearized supersonic flows")
        liste_chapter.append("One-dimensional compressible flows")
        liste_chapter.append("Two-dimensional supersonic flows")
        liste_chapter.append("Method of characteristics for steady and unsteady flows")
        liste_chapter.append("Slender bodies and wings in supersonic flows")
        liste_chapter.append("Hypersonic flows")
        liste_chapter.append("Viscous flows and boundary layer")
    return liste_chapter


def liste_exo_python(display=False, lang="fr"):
    """
    list of the exercises by chapter, in french and english

    Args:
        lang (string): langage fr or en
    """

    liste_exo = []
    if lang == "fr":
        liste_exo.append({'ch': 2, 'exo': 'Exercice2_1', 'com': "Fonction de courant et potentiel des vitesses"})
        liste_exo.append({'ch': 2, 'exo': 'Exercice2_2', 'com': "Atmosphère standard"})

        liste_exo.append({'ch': 3, 'exo': 'Exercice3_1', 'com': "Cylindre "})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_2', 'com': "Transformation de Joukowski"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_3', 'com': "Transformation de Karman-Trefftz"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_4', 'com': "Transformation de Von Mises"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_5', 'com': "Transformation de Van de Vooren et Jong"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_6', 'com': "Profil à double courbure"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_7', 'com': "Plaque plane en incidence "})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_8', 'com': "Plaque plane : introduction panneaux "})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_9', 'com': "Effet de sol"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_12', 'com': "Ovale de Rankine"})

        liste_exo.append({'ch': 4, 'exo': 'Exercice4_1', 'com': "Plaque plane en incidence"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_2', 'com': "profil de Joukowski"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_3', 'com': "profil famille pq "})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_4', 'com': "profil double cambrure +  becs et volets "})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_5', 'com': "Lois de squelette et d'épaisseur"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_6', 'com': "Exercice sur les becs et volets, Aile rectangulaire "})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_7', 'com': "Profil en théorie linéarisée : profil NACA"})

        liste_exo.append({'ch': 5, 'exo': 'Exercice5_1', 'com': "Aile optimale elliptique "})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_2', 'com': "Vrillage d'une aile elliptique"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_3', 'com': "Braquage d'un volet d'une aile elliptique"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_4', 'com': "Calcul du foyer aérodynamique d'une aile"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_5', 'com': "Sillage tourbillonnaire d'un Airbus A380-800"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_6', 'com': "Corrections de soufflerie"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_7', 'com': "Tourbillons en fer à cheval"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_8', 'com': "Vol en formation"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_9', 'com': "Tourbillons de sillage au voisinage du sol"})
        liste_exo.append(
            {'ch': 5, 'exo': 'Exercice5_10', 'com': "Résolution numérique du problème de Prandtl, cas général"})
        liste_exo.append(
            {'ch': 5, 'exo': 'Exercice5_11', 'com': "Résolution numérique du problème de Prandtl avec effet "
                                                    "de sol"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_12', 'com': "demi-aile d'un spitfire"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_14', 'com': "Chargement d'une aile trapézoïdale de P51 Mustang"})
        # liste_exo.append({'ch':5,'exo':'Exercice5_13','com':"Résolution numérique du problème de Prandtl, cas général (
        # OLD)"})

        liste_exo.append({'ch': 6, 'exo': 'Exercice6_1', 'com': "Loi du coefficient de portance avec l'allongement "})
        liste_exo.append({'ch': 6, 'exo': 'Exercice6_3', 'com': "Ellipsoide"})
        liste_exo.append({'ch': 6, 'exo': 'Exercice6_4', 'com': "Analyse des données de la thèse de Van Driest sur "
                                                                "l'ellipsoide"})

        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_0', 'com': "calcul des paramètres du profil de Van-Hooren {7.1.1}"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_1', 'com': "Profil de Van-Hooren, carte des paramètres dans le plan "
                                                   "k-t (épaisseur) {7.1.1}"})
        liste_exo.append({'ch': 7, 'exo': 'Exercice7_2', 'com': "Profil parabolique {7.1.2}"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_3', 'com': "Discrete Vortex Method on the parabolic mean line airfoil "
                                                   "{7.2.1}"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_4', 'com': "Profil de Van-Hooren, distribution de sources, Approche "
                                                   "de Neumann {7.2.2}"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_5', 'com': "Distribution uniforme de sources et de doublets sur le "
                                                   "profil de Van-Hooren {7.2.3} "})
        liste_exo.append({'ch': 7, 'exo': 'Exercice7_10', 'com': " Profil NACA : dessin"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_30', 'com': "Solution de référence sur un profil Naca depuis fichiers "
                                                    "NACA"})

        liste_exo.append({'ch': 8, 'exo': 'Exercice8_1', 'com': "Total pressure / Stagnation pressure = Isentropic "
                                                                "Pressure {8.1.1}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_2', 'com': "Pitot tube in subsonic flow {8.1.2} "})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_3', 'com': "Compressibility corrections (airfoil) {8.2.1}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_4', 'com': "Compressibility corrections (wing) and critical Mach "
                                                                "number {8.2.2}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_5', 'com': "paroi ondulée {8.3}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_6', 'com': "paroi ondulée en transsonique {8.4} "})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_7', 'com': "Calcul de Kpinc en fonction du Mach critique {8.5}"})

        liste_exo.append(
            {'ch': 9, 'exo': 'Exercice9_1', 'com': "Théorie linéaire : plaque plane en incidence {9.2.1} "})
        liste_exo.append({'ch': 9, 'exo': 'Exercice9_2', 'com': "Théorie linéaire :  profil losangique {9.2.3}"})
        liste_exo.append(
            {'ch': 9, 'exo': 'Exercice9_3', 'com': "Théorie linéaire : profil de traînée minimale {9.2.4} "})
        liste_exo.append({'ch': 9, 'exo': 'Exercice9_4', 'com': "Théorie linéaire : profil  lenticulaire {9.2.5}"})
        liste_exo.append({'ch': 9, 'exo': 'Exercice9_5', 'com': "Interaction jet supersonique - profil {9.2.6}"})

        liste_exo.append({'ch': 10, 'exo': 'Exercice10_1', 'com': "Poussée d'une tuyère de fusée"})
        liste_exo.append(
            {'ch': 10, 'exo': 'Exercice10_2', 'com': "Poussée d'une tuyère de fusée, formules du Mattingly."})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_3', 'com': "Tube de Pitot dans une tuyère amorcée"})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_4', 'com': "Moteur oxygène - hydrogène"})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_5', 'com': "Exemple de calcul d'un choc droit"})

        liste_exo.append({'ch': 11, 'exo': 'Exercice11_1', 'com': "courbes p_1/p_0 en fonction de la déviation"})
        liste_exo.append({'ch': 11, 'exo': 'Exercice11_2', 'com': "Réflexion d'un choc dans un canal"})
        liste_exo.append({'ch': 11, 'exo': 'Exercice11_3', 'com': "Interaction de 2 chocs obliques "})

        liste_exo.append({'ch': 12, 'exo': 'Exercice12_1', 'com': "Canal divergent"})
        liste_exo.append({'ch': 12, 'exo': 'Exercice12_2', 'com': "Tube à choc"})
        liste_exo.append({'ch': 12, 'exo': 'Exercice12_3', 'com': "Tuyère d'une fusée"})
        liste_exo.append({'ch': 12, 'exo': 'Exercice12_4', 'com': "Tuyère de longueur minimale"})

        liste_exo.append({'ch': 13, 'exo': 'Exercice13_1', 'com': "Obstacle conique en supersonique "})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_2', 'com': "ogive parabolique en supersonique"})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_3', 'com': "Traînée minimale d'un corps élancé  supersonique "})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_4', 'com': "Tracé des iso-potentielles et lignes de courant "})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_5', 'com': "Corps conique élancé "})

        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_0', 'com': "Calcul du Mach aval pour un choc oblique sur un dièdre "
                                                     "par différentes approximations "})
        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_1', 'com': "profil losangique et triangulaire, méthode de Newton"})
        # liste_exo.append({'ch':14,'exo':'Exercice14_2','com':""})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_3', 'com': "Kp pour des écoulements de dièdre 2D ou coniques, "
                                                                  "pour un nombre de Mach infini et d'autres "})
        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_4', 'com': "Comparaison des méthodes de calcul autour d'un obstacle "
                                                     "pointus"})
        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_5', 'com': "forme d'un corps de traînée minimale en hypersonique "})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_6', 'com': "Traînée d'un corps émoussé en hypersonique "})
        liste_exo.append({'ch': 13, 'exo': 'Exercice14_100', 'com': "Courbes des chocs coniques  "})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_101', 'com': "test de l'approximation sur la fonction de "
                                                                    "Prandtl-Meyer en hypersonique"})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_102', 'com': "Test de l'inversion de la fonction omega(M) en "
                                                                    "supersonique"})

        liste_exo.append({'ch': 15, 'exo': 'Exercice15_1', 'com': "Plaque plane : Loi en sinus et polynome "})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_2', 'com': "Plaque plane : Loi polynomiale"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_4', 'com': "Mise en mouvement d'une plaque, problème de Rayleigh"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_6', 'com': "Ecoulement de Couette instationnaire"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_7', 'com': "Ecoulement instationnaire : second problème de Stokes"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_8', 'com': "Jet laminaire"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_9', 'com': "Calcul de la position de la transition"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_10', 'com': "Ancienne version, mais fonctionne, pas de Class"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_11', 'com': "profil turbulent analytique avec aspiration"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_12', 'com': "profil turbulent analytique avec une loi puissance"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_13', 'com': "profil turbulent analytique avec la correction de Van "
                                                      "Driest"})
    else:
        liste_exo.append({'ch': 2, 'exo': 'Exercice2_1', 'com': "Stream functions and potential of velocity"})
        liste_exo.append({'ch': 2, 'exo': 'Exercice2_2', 'com': "Standard atmosphere"})

        liste_exo.append({'ch': 3, 'exo': 'Exercice3_1', 'com': "Cylinder "})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_2', 'com': "Joukowski's transformation"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_3', 'com': "Karman-Trefftz's transformation"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_4', 'com': "Von Mises's transformation"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_5', 'com': "Van de Vooren and Jong transformation"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_6', 'com': "Double-cambered airfoil"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_7', 'com': "Flate plate with incidence"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_8', 'com': "Flate plate : introduction to panels "})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_9', 'com': "Ground effect"})
        # liste_exo.append({'ch':3,'exo':'Exercice3_10','com':"Tandem wings"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_12', 'com': "Rankine's oval"})
        liste_exo.append({'ch': 3, 'exo': 'Exercice3_20', 'com': "Conformal transformation applications"})

        liste_exo.append({'ch': 4, 'exo': 'Exercice4_1', 'com': "Flate plate with incidence"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_2', 'com': "Joukowski airfoil"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_3', 'com': "pq family airfoils "})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_4', 'com': "Double cambered airfoil +  slats and flaps "})
        # liste_exo.append({'ch':4,'exo':'Exercice4_5','com':"Camber and thickness laws"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_6', 'com': "Slots and flaps on a rectangular wing"})
        liste_exo.append({'ch': 4, 'exo': 'Exercice4_7', 'com': "Linearized airfoil theory : NACA airfoil"})

        liste_exo.append({'ch': 5, 'exo': 'Exercice5_1', 'com': "Elliptic optimal wing"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_2', 'com': "Twisted elliptic wing"})
        # liste_exo.append({'ch':5,'exo':'Exercice5_3','com':"Flaps on a elliptic wing"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_4', 'com': "Wing aerodynamic center"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_5', 'com': "Airbus A380-800 vortex wake"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_6', 'com': "Wind tunnel corrections"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_7', 'com': "Horse-shoe vortex"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_8', 'com': "Formation flying"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_9', 'com': "Wake vortex with ground effect"})
        liste_exo.append(
            {'ch': 5, 'exo': 'Exercice5_10', 'com': "Numerical solution of the Prandt equation : general case"})
        liste_exo.append(
            {'ch': 5, 'exo': 'Exercice5_11', 'com': "Numerical solution of the Prandt equation with ground "
                                                    "effect"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_12', 'com': "Spitfire wing"})
        liste_exo.append({'ch': 5, 'exo': 'Exercice5_14', 'com': "Tappered wing load of the  P51 Mustang"})
        # liste_exo.append({'ch':5,'exo':'Exercice5_13','com':"Résolution numérique du problème de Prandtl, cas général (
        # OLD)"})

        liste_exo.append({'ch': 6, 'exo': 'Exercice6_1', 'com': "Lift coefficient law with  the wing aspect ratio"})
        liste_exo.append({'ch': 6, 'exo': 'Exercice6_3', 'com': "Ellipsoide"})
        liste_exo.append(
            {'ch': 6, 'exo': 'Exercice6_4', 'com': "Data analysis of  Van Driest phd about the ellipsoidal flow"})

        liste_exo.append({'ch': 7, 'exo': 'Exercice7_0', 'com': "Van-Hooren airfoil parameters {7.2.1}"})
        liste_exo.append({'ch': 7, 'exo': 'Exercice7_1', 'com': "Van-Hooren's airfoil, parameter map in the k-t plane ("
                                                                "thickness) {7.1.1}"})
        liste_exo.append({'ch': 7, 'exo': 'Exercice7_2', 'com': "Parabolic airfoil {7.1.2}"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_3', 'com': "Discrete Vortex Method on the parabolic mean line airfoil "
                                                   "{7.2.1}"})
        liste_exo.append({'ch': 7, 'exo': 'Exercice7_4', 'com': "Van-Hooren's airfoil, sources distribution, Neumann's "
                                                                "approach {7.2.2}"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_5', 'com': "Uniform distribution of  sources and of  doublets for the "
                                                   " Van-Hooren airfoil {7.2.3} "})
        liste_exo.append({'ch': 7, 'exo': 'Exercice7_10', 'com': "Airfoil NACA : drawings"})
        liste_exo.append(
            {'ch': 7, 'exo': 'Exercice7_30', 'com': "Reference solution of a Naca airfoil from NACA files"})

        liste_exo.append({'ch': 8, 'exo': 'Exercice8_1', 'com': "Total pressure / Stagnation pressure = Isentropic "
                                                                "Pressure {8.1.1}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_2', 'com': "Pitot tube in subsonic flow {8.1.2} "})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_3', 'com': "Compressibility corrections (airfoil) {8.2.1}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_4', 'com': "Compressibility corrections (wing) and critical Mach "
                                                                "number {8.2.2}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_5', 'com': "Wavy wall flow {8.3}"})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_6', 'com': "Wavy wall flow on transonic regime {8.4} "})
        liste_exo.append({'ch': 8, 'exo': 'Exercice8_7', 'com': "Kpinc with respect to the critical Mach number {8.4}"})

        liste_exo.append({'ch': 9, 'exo': 'Exercice9_1', 'com': "Linearized theory : flat plate in incidence {9.2.1} "})
        liste_exo.append({'ch': 9, 'exo': 'Exercice9_2', 'com': "Linearized theory : diamond airfoil {9.2.3}"})
        liste_exo.append({'ch': 9, 'exo': 'Exercice9_3', 'com': "Linearized theory : airfoil of minimal drag {9.2.4} "})
        # liste_exo.append({'ch':9,'exo':'Exercice9_4','com':"Linearized theory : lenticular airfoil {9.2.5}"})
        liste_exo.append({'ch': 9, 'exo': 'Exercice9_5', 'com': "Supersonic jet - airfoil interaction {9.2.6}"})

        liste_exo.append({'ch': 10, 'exo': 'Exercice10_1', 'com': "Rocket nozzle thrust"})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_2', 'com': "Rocket nozzle thrust, Mattingly's formules"})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_3', 'com': "Pitot tube into a supersonic nozzle"})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_4', 'com': "Oxygen - hydrogen rocket engine"})
        liste_exo.append({'ch': 10, 'exo': 'Exercice10_5', 'com': "normal shock wave example"})

        liste_exo.append(
            {'ch': 11, 'exo': 'Exercice11_1', 'com': "p_1/p_0 plots with respect to velocity deviation angle"})
        liste_exo.append({'ch': 11, 'exo': 'Exercice11_2', 'com': "Oblique shock wave reflection of a channel wall"})
        liste_exo.append({'ch': 11, 'exo': 'Exercice11_3', 'com': "Interaction between two oblique shock waves"})

        liste_exo.append({'ch': 12, 'exo': 'Exercice12_1', 'com': "Diverging channel"})
        liste_exo.append({'ch': 12, 'exo': 'Exercice12_2', 'com': "Shock tube"})
        liste_exo.append({'ch': 12, 'exo': 'Exercice12_3', 'com': "Rocket nozzle"})
        liste_exo.append({'ch': 12, 'exo': 'Exercice12_4', 'com': "Nozzle of minimal length"})

        liste_exo.append({'ch': 13, 'exo': 'Exercice13_1', 'com': "Conical obstacle in supersonic flow"})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_2', 'com': "Parabolic ogive in supersonic flow"})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_3', 'com': "Minimal drag of a slender body in supersonic flow"})
        liste_exo.append(
            {'ch': 13, 'exo': 'Exercice13_4', 'com': "Plot of the iso-potential lines and  of the stream lines"})
        liste_exo.append({'ch': 13, 'exo': 'Exercice13_5', 'com': "Slender conical body"})

        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_0', 'com': "Downstream Mach of an oblique shock wave over a plane "
                                                     "wall with  different approximations"})
        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_1', 'com': "Diamond and triangular airfoils Kp: Newton's method"})
        # liste_exo.append({'ch':14,'exo':'Exercice14_2','com':""})
        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_3', 'com': "Kp of a 2D wedge or conical wall, for various Mach number"})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_4', 'com': "Sharp obstacle Kp with different approaches"})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_5', 'com': "Minimal drag body in hypersonic flow"})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_6', 'com': "Drag of blunt body in hypersonic flow"})
        liste_exo.append({'ch': 13, 'exo': 'Exercice14_100', 'com': "Conical shock angle curves"})
        liste_exo.append({'ch': 14, 'exo': 'Exercice14_101', 'com': "Test: Prandtl-Meyer approximated function in   "
                                                                    "hypersonic flow"})
        liste_exo.append(
            {'ch': 14, 'exo': 'Exercice14_102', 'com': "Test : inverted omega(Mach) function  in supersonic "
                                                       "flow"})

        liste_exo.append({'ch': 15, 'exo': 'Exercice15_1', 'com': "Flat plate : sine and polynomial laws "})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_2', 'com': "Flat plate : polynomial laws"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_4', 'com': "Instantaneous movement of a plate, Rayleigh's problem"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_6', 'com': "Unsteady Couette's flow"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_7', 'com': "Unsteady flow : second Stokes problem"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_8', 'com': "Steady laminar jet"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_9', 'com': "Transition to turbulence onset location"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_10', 'com': "old version, it works without  python class"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_11', 'com': "Turbulent analytical velocity profile with suction"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_12', 'com': "Turbulent analytical velocity profile with the power law"})
        liste_exo.append(
            {'ch': 15, 'exo': 'Exercice15_13', 'com': "Turbulent analytical velocity profile with the  Van "
                                                      "Driest's correction"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_14', 'com': "analytical turbulent profile with the Van Driest "
                                                                   "correction  (for tests)"})
        liste_exo.append({'ch': 15, 'exo': 'Exercice15_15', 'com': "analytical turbulent profile with suction ("
                                                                   "for tests) "})


    if display:
        print("=" * 90)
        if lang == "fr":
            print("k  chap.     exercice   \t commentaires")
        else:
            print("k  chap.     exercise   \t comments")
        print("=" * 90)

        k = 0
        for dico in liste_exo:
            # for key,values in dico.items():
            print("%2i   %2i %20s  %s" % (k, dico['ch'], dico['exo'], dico['com']))
            k = k + 1
        if lang == "fr":
            print("Pour choisir un  ou plusieurs exercices, il faut uniquement entrer l'index de la première colonne")
            print("dans le fichier main.py")
        else:
            print("To choose one or several exercises you have to enter the index of the first column")
            print("in the main.py file")

    return liste_exo


def print_exo(chap, nbre):
    """
    format of the exercise references
    """
    ch = " " * 10
    print("_" * 75, '\n')
    print('%s \tCh : %d \t EXERCISE n° %d' % (ch, chap, nbre))
    print("_" * 75)


def correction(chapter, exos):
    """ 
    Selection of  the exercises and of the chapter 

    Args:
        chapter (int): chapter of the exercise
        nbre (int): list of the exercise inside a chapter
    """
    ch = '*' * 25
    i = chapter
    print('\n%s \tCHAPTER n° %d \t  %s\n' % (ch, i, ch))
    if i == 2:
        '''
        Modèles mathématiques pour l'aérodynamique 
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """
                Fonction de courant et potentiel des vitesses
                *Stream functions and potential of velocity*
                """
                Exercice2_1()
            if j == 2:
                """
                Atmosphère Standard
                *Standard atmosphere*
                """
                Exercice2_2()
    if i == 3:
        '''
        Théorie des potentiels complexes : applications aux profils d'aile
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """
                Cylindre
                """
                Exercice3_1()
            elif j == 2:
                """
                Transformation de Joukowski
                """
                Exercice3_2()
            elif j == 3:
                """
                Transformation de Karman-Trefftz
                """
                Exercice3_3()
            elif j == 4:
                """
                Transformation de Von Mises
                """
                Exercice3_4()
            elif j == 5:
                """
                Transformation de Van de Vooren et Jong
                """
                Exercice3_5()
            elif j == 6:
                """
                Profil à double courbure
                *Double-cambered airfoil*
                """
                Exercice3_6()
            elif j == 7:
                """
                Plaque plane en incidence
                *Flate plate with incidence*
                """
                Exercice3_7()
            elif j == 8:
                """
                Plaque plane : introduction panneaux
                *Flate plate : introduction to panels*
                """
                Exercice3_8()

            elif j == 9:
                """
                Effet de sol
                *Ground effect*
                """
                Exercice3_9()
            elif j == 12:
                """
                Ovale de Rankine
                """
                Exercice3_12()
            elif j == 20:
                """
                application of complex mappings (old version)
                """
                Exercice3_20()

    if i == 4:
        '''
        Théorie linéarisée des profils minces
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """
                Plaque plane en incidence
                *Flate plate with incidence*
                """
                Exercice4_1()
            if j == 2:
                """
                profil de Joukowski
                """
                Exercice4_2()
            if j == 3:
                """
                profil family_pq 
                """
                Exercice4_3()
            if j == 4:
                """
                profil double cambrure +  becs et volets 
                *Double cambered airfoil +  slats and flaps*
                """
                Exercice4_4()
            if j == 5:
                """
                Lois de squelette et d'épaisseur
                *Camber and thickness laws*
                """
                Exercice4_5()
            if j == 6:
                """
                Exercice sur les becs et volets, Aile rectangulaire
                *Slots and flaps on a rectangular wing*
                """
                Exercice4_6()
            if j == 7:
                """
                Profil en théorie linéarisée : profil NACA
                *Linearized airfoil theory : NACA airfoil*
                """
                Exercice4_7()
            if j == 8:
                """
                Profil en théorie linéarisée : profil NACA
                *Linearized airfoil theory : NACA airfoil*
                """
                Exercice4_8()

    if i == 5:
        '''
        Théorie de la ligne portante de Lanchester-Prandtl
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """
                Aile optimale elliptique 
                *Elliptic optimal wing*
                """
                Exercice5_1()
            if j == 2:
                """
                Vrillage d'une aile elliptique
                *Twisted elliptic wing*
                """
                Exercice5_2()

            if j == 3:
                """
                Braquage d'un volet d'une aile elliptique
                *Flaps on a elliptic wing*
                """
                Exercice5_3()

            if j == 4:
                """
                Calcul du foyer aérodynamique d'une aile
                *Wing aerodynamic center*
                """
                Exercice5_4()

            if j == 5:
                """
                Sillage tourbillonnaire d'un Airbus A380-800
                *Airbus A380-800 vortex wake*
                """
                Exercice5_5()

            if j == 6:
                """
                Correction en soufflerie
                *Wind tunnel corrections*
                """
                Exercice5_6()

            if j == 7:
                """
                Tourbillon en fer à cheval
                *Horse-shoe vortex*
                """
                Exercice5_7()

            if j == 8:
                """
                Vol en formation
                *Flight in formation*
                """
                Exercice5_8()

            if j == 9:
                """
                Tourbillons de sillage au voisinage du sol
                *Wake vortex with ground effect*
                """
                Exercice5_9()

            if j == 10:
                """
                Résolution numérique du problème de Prandtl, cas général
                *Numerical solution of the Prandt equation : general case*
                """
                Exercice5_10()

            if j == 11:
                """
                Résolution numérique du problème de Prandtl avec effet de sol
                *Numerical solution of the Prandt equation with ground effect*
                """
                Exercice5_11()

            if j == 12:
                """
                demi-aile d'un spitfire
                *Spitfire wing*
                """
                Exercice5_12()

            if j == 13:
                """
                aile elliptique sans effet de sol
                elliptic wing without ground effect
                """
                Exercice5_13()

            if j == 14:
                """
                chargement de l'aile trapézoidale du P51 Mustang
                *Tappered wing load of the  P51 Mustang*
                """
                Exercice5_14()

            if j == 15:
                """
                test de l'intégrale avec quad 
                """
                Exercice5_15()

    if i == 6:
        '''
        Théorie de la surface portante et des corps élancés
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """
                Loi du coefficient de portance avec l'allongement
                *Lift coefficient law with  the wing aspect ratio*
                """
                Exercice6_1()
            if j == 3:
                """
                Ellipsoide
                """
                Exercice6_3()
            if j == 4:
                """
                Analyse des données de la thèse de Van Driest sur l'ellipsoide
                *Data analysis of  Van Driest phd about the ellipsoidal flow*
                """
                Exercice6_4()

    if i == 7:
        '''
        Aspects numériques : méthode des singularités
        '''
        for j in exos:
            print_exo(i, j)
            if j == 0:
                """
                calcul des paramètres du profil de Van-Hooren {7.2.1} 
                *Van-Hooren airfoil parameters*
                """
                Exercice7_0()
            if j == 1:
                """ 
                Profil de Van-Hooren, carte des paramètres dans le plan k-t (épaisseur) {7.1.1} 
                *Van-Hooren's airfoi, parameter map in the k-t plane (thickness)*
                """
                Exercice7_1()
            if j == 2:
                """ 
                Profil parabolique {7.2.2}
                """
                Exercice7_2()
            if j == 3:
                """
                Discrete Vortex Method on the parabolic mean line airfoil {7.3.1}
                """
                Exercice7_3()
            if j == 4:
                """ 
                Profil de Van-Hooren, distribution de sources, Approche de Neumann {7.3.2}
                """
                Exercice7_4()
            if j == 5:
                """ 
                Distribution uniforme de sources et de doublets sur le profil de Van-Hooren {7.3.3}
                """
                Exercice7_5()
            if j == 10:
                """ 
                Profil NACA : dessin
                """
                Exercice7_10()
            if j == 30:
                """ 
                Solution de référence sur un profil Naca depuis fichiers NACA
                dans misc.py
                """

                Exercice7_30()

    if i == 8:
        '''
        Ecoulements compressible subsonique et transsonique 
        '''
        for j in exos:
            print_exo(i, j)

            if j == 1:
                """ 
                *Total pressure / Stagnation pressure = Isentropic Pressure {8.2.1}*
                """
                Exercice8_1()
            if j == 2:
                """ 
                *Pitot tube in subsonic flow {8.2.2}*
                """
                Exercice8_2()
            if j == 3:
                """ 
                *Compressibility corrections (airfoil) {8.3.1}*
                """
                Exercice8_3()
            if j == 4:
                """ 
                *Compressibility corrections (wing) and critical Mach number {8.3.2}*
                """
                Exercice8_4()
            if j == 5:
                """ 
                paroi ondulée {8.3}
                *Wavy wall flow {8.3}*
                """
                Exercice8_5()
            if j == 6:
                """ 
                paroi ondulée en transsonique {8.4}
                *Wavy wall flow on transonic regime {8.4}* 
                
                """
                Exercice8_6()
            if j == 7:
                """ 
                Calcul de Kpinc en fonction du Mach critique {8.5}
                *Kpinc with respect to the critical Mach number {8.5}*
                """
                Exercice8_7()

    if i == 9:
        '''
        Ecoulements supersoniques linéarisés  
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """ 
                Théorie linéaire : plaque plane en incidence {9.2.1}
                *Linearized theory : flat plate in incidence {9.2.1}*
                """
                Exercice9_1()
            if j == 2:
                """ 
                Théorie linéaire :  profil losangique {9.2.3}
                *Linearized theory : diamond airfoil {9.2.3}*
                """
                Exercice9_2()
            if j == 3:
                """ 
                Théorie linéaire : profil de traînée minimale {9.2.4}
                *Linearized theory : airfoil of minimal drag {9.2.4}*
                """
                Exercice9_3()
            if j == 4:
                """ 
                Théorie linéaire : profil  lenticulaire {9.2.5}
                *Linearized theory : lenticular airfoil {9.2.5}*
                """
                Exercice9_4()
            if j == 5:
                """ 
                Interaction jet supersonique - profil {9.2.6}
                *Supersonic jet - airfoil interaction {9.2.6}*
                """
                Exercice9_5()

    if i == 10:
        '''
        Ecoulements supersoniques monodimensionnels 
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """ 
                Poussée d'une tuyère de fusée
                *Rocket nozzle thrust*
                """
                Exercice10_1()
            if j == 2:
                """ 
                Poussée d'une tuyère de fusée, formules du Mattingly.
                *Rocket nozzle thrust, Mattingly's formules*
                """
                Exercice10_2()
            if j == 3:
                """
                Tube de Pitot dans une tuyère amorcée
                *Pitot tube into a supersonic nozzle*
                """
                Exercice10_3()
            if j == 4:
                """
                Moteur oxygène - hydrogène
                *Oxygen - hydrogen rocket engine*
                """
                Exercice10_4()
            if j == 5:
                """
                Calcul élémentaire d'un choc droit
                *normal shock wave example*
                """
                Exercice10_5()
    if i == 11:
        '''
        Ecoulements supersoniques bidimensionnels 
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """ 
                courbes p_1/p_0 en fonction de la déviation
                *p_1/p_0 plots with respect to velocity deviation angle*
                """
                Exercice11_1()
            if j == 2:
                """ 
                Réflexion d'un choc dans un canal
                *Oblique shock wave reflection of a channel wall*
                """
                Exercice11_2()
            if j == 3:
                """ 
                Interaction de 2 chocs obliques
                *Interaction between two oblique shock waves*
                """
                Exercice11_3()

    if i == 12:
        '''
        Méthode des caractéristiques 
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """ 
                Canal divergent
                *Diverging channel*
                """
                Exercice12_1()
            if j == 2:
                """ 
                Tube à choc
                *Shock tube*
                """
                Exercice12_2()
            if j == 3:
                """ 
                Tuyère d'une fusée
                Rocket nozzle
                """
                Exercice12_3()
            if j == 4:
                """ 
                Tuyère de longueur minimale
                *Nozzle of minimal length*
                """
                Exercice12_4()


    elif i == 13:
        '''
        Ecoulements supersoniques : corps élancés et ailes
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """ 
                Obstacle conique en supersonique 
                *Conical obstacle in supersonic flow*
                """
                Exercice13_1()
            if j == 2:
                """ 
                ogive parabolique en supersonique 
                *Parabolic ogive in supersonic flow*
                """
                Exercice13_2()
            if j == 3:
                """ 
                Traînée minimale d'un corps élancé  supersonique 
                *Minimal drag of a slender body in supersonic flow*
                """
                Exercice13_3()
            if j == 4:
                """
                Tracé des iso-potentielles et lignes de courant
                *Plot of the iso-potential lines and  of the stream lines*
                """
                Exercice13_4()
            if j == 5:
                """
                Corps conique élancé
                *Slender conical body*
                """
                Exercice13_5()
    elif i == 14:
        '''
        Ecoulements hypersoniques
        '''
        for j in exos:
            print_exo(i, j)
            if j == 0:
                '''
                Calcul du Mach aval pour un choc oblique sur un dièdre par différentes approximations
                *Downstream Mach of an oblique shock wave over a plane wall with  different approximations*
                '''
                Exercice14_0()
            if j == 1:
                Exercice14_1()
                """
                profil losangique et triangulaire, méthode de Newton 
                *Diamond and triangular airfoils Kp: Newton's method*
                """
            # if j == 2:         il n'existe pas, pas besoin de python
            #   Exercice14_2()
            if j == 3:
                Exercice14_3()
                '''
                Kp pour des écoulements de dièdre 2D ou coniques, pour un nombre de Mach infini et d'autres
                *Kp of a 2D plane or conical wall, for various Mach number*
                '''
            if j == 4:
                '''
                Comparaison des méthodes de calcul autour d'un obstacle pointus
                *Sharp obstacle Kp with different approaches*
                '''
                Exercice14_4()
            if j == 5:
                '''
                forme d'un corps de traînée minimale en hypersonique
                *Minimal drag body in hypersonic flow*
                '''
                Exercice14_5()
            if j == 6:
                '''
                traînée d'un corps émoussé en hypersonique
                *Drag of blunt body in hypersonic flow*
                '''
                Exercice14_6()
            if j == 100:
                """ 
                Courbes des chocs coniques pour les rappels
                *Conical shock angle curves*
                """
                Exercice14_100()

            if j == 101:
                """
                test de l'approximation sur la fonction de Prandtl-Meyer en hypersonique
                cas = 0 : plusieurs M0, cas =1 , 1 M0, cas = 2: fonction omega

                *Test: Prandtl-Meyer approximated function in hypersonic flow*
                """
                Exercice14_101(cas=0)

            if j == 102:
                """
                Test de l'inversion de la fonction omega(M) en supersonique
                *Test : inverted omega(Mach) function  in supersonic flow*
                """
                Exercice14_102()
    elif i == 15:
        '''
        Effets visqueux et couche limite 
        '''
        for j in exos:
            print_exo(i, j)
            if j == 1:
                """ 
                Plaque plane : Loi en sinus et polynome
                cette correction va plus loin que le sujet car on va jusqu'à
                    des polynomes de degrés 7.
                *Flat plate : sine and polynomial laws (up to 7th degrees)*
                """
                Exercice15_1()
            if j == 2:
                """ 
                Plaque plane : Loi polynomiale
                *Flat plate : polynomial laws*
                """
                Exercice15_2()
            if j == 4:
                """ 
                Mise en mouvement d'une plaque, problème de Rayleigh
                *Instantaneous movement of a plate, Rayleigh's problem*
                """
                Exercice15_4()
            if j == 6:
                """ 
                Ecoulement de Couette instationnaire
                *Unsteady Couette's flow*
                """
                Exercice15_6()
            if j == 7:
                """ 
                Ecoulement instationnaire : second problème de Stokes
                *Unsteady flow : second Stokes problem*
                """
                Exercice15_7()
            if j == 8:
                """ 
                Ecoulement stationnaire : jet laminaire
                *Steady Laminar jet*
                """
                Exercice15_8()
            if j == 9:
                """
                Calcul de la position de la transition
                *Transition to turbulence  onset location*
                """
                Exercice15_9()
            if j == 10:
                '''
                profil turbulent avec aspiration (ancienne version sans classe)
                *old version, it works without  python class*
                '''
                Exercice15_10()
            if j == 11:
                """
                profil turbulent analytique avec aspiration
                Utilisation de la classe Turbulence
                *Turbulent analytical velocity profile with suction*
                """
                Exercice15_11()
            if j == 12:
                """
                profil turbulent analytique avec une loi puissance
                Utilisation de la classe Turbulence
                *Turbulent analytical velocity profile with the power law*
                """
                Exercice15_12()
            if j == 13:
                """
                profil turbulent analytique avec la correction de Van Driest
                Utilisation de la classe Turbulence
                *Turbulent analytical velocity profile with the  Van Driest's correction*
                """
                Exercice15_13()
            if j == 14:
                """
                profil turbulent analytique avec la correction de Van Driest
                tau fonction de y/delta
                Utilisation de la classe Turbulence
                *analytical turbulent profile with the Van Driest correction  (with the turbulence class)*
                """
                Exercice15_14()
            if j == 15:
                """
                profil turbulent analytique avec aspiration
                Utilisation de la classe Turbulence
                analytical turbulent profile with suction (turbulence class)
                """
                Exercice15_15()


def solution(liste):
    """
    choose the exercises in the input list

    Args:
        liste (list): list of the exercises to solve

    """
    for L in liste:
        # print('chapter = %i'%(L[0]))
        # for exo in L[1]:
        # print("\t exercice = %i"%exo)
        correction(L[0], L[1])
