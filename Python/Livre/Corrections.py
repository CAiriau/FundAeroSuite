 #!/bin/python
"""
Module comportant le lien entre les numéros d'exercice 
et le nom des exercices.
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 




import os
import numpy             as np
import matplotlib.pyplot as plt

from Tools.misc                 import *
from CompressibleFlow.fonctions import *
from Livre.Chapter2             import *
from Livre.Chapter3             import *
from Livre.Chapter4             import *
from Livre.Chapter5             import *
from Livre.Chapter6             import *
from Livre.Chapter7             import *
from Livre.Chapter7c            import *
from Livre.Chapter8             import *
from Livre.Chapter9             import *
from Livre.Chapter10            import *
from Livre.Chapter11            import *
from Livre.Chapter12            import *
from Livre.Chapter12b           import *
from Livre.Chapter13            import *
from Livre.Chapter14            import *
from Livre.Chapter15_Laminaire  import *
from Livre.Chapter15_Transition import *
from Livre.Chapter15_Turbulent  import *
from Livre.Chapter15_Tu         import *
from Livre.misc                 import *
from Choc_Conique.choc_conique  import * 


def liste_chapters():
    liste_chapter=[]
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
    liste_chapter.append("Ecoulements hyersoniques")
    liste_chapter.append("Effets visqueux et couche limite")
    return liste_chapter



def liste_exo_python(display=False):

    liste_exo=[]
    liste_exo.append({'ch':2,'exo':'Exercice2_1','com':"Fonction de courant et potentiel des vitesses"})
    liste_exo.append({'ch':2,'exo':'Exercice2_2','com':"Atmsphère standard"})
    
    liste_exo.append({'ch':3,'exo':'Exercice3_1','com':"Cylindre "})
    liste_exo.append({'ch':3,'exo':'Exercice3_2','com':"Transformation de Joukowski"})
    liste_exo.append({'ch':3,'exo':'Exercice3_3','com':"Transformation de Karman-Trefftz"})
    liste_exo.append({'ch':3,'exo':'Exercice3_4','com':"Transformation de Von Mises"})
    liste_exo.append({'ch':3,'exo':'Exercice3_5','com':"Transformation de Van de Vooren et Jong"})
    liste_exo.append({'ch':3,'exo':'Exercice3_6','com':"Profil à double courbure"})
    liste_exo.append({'ch':3,'exo':'Exercice3_7','com':"Plaque plane en incidence "})
    liste_exo.append({'ch':3,'exo':'Exercice3_8','com':"Plaque plane : introduction panneaux "})
    liste_exo.append({'ch':3,'exo':'Exercice3_9','com':"Effet de sol"})
    liste_exo.append({'ch':3,'exo':'Exercice3_12','com':"Ovale de Rankine"})
    
    liste_exo.append({'ch':4,'exo':'Exercice4_1','com':"Plaque plane en incidence"})
    liste_exo.append({'ch':4,'exo':'Exercice4_2','com':"profil de Joukowski"})
    liste_exo.append({'ch':4,'exo':'Exercice4_3','com':"profil famille_pq "})
    liste_exo.append({'ch':4,'exo':'Exercice4_4','com':"profil double cambrure +  becs et volets "})
    liste_exo.append({'ch':4,'exo':'Exercice4_5','com':"Lois de squelette et d'épaisseur"})
    liste_exo.append({'ch':4,'exo':'Exercice4_6','com':"Exercice sur les becs et volets, Aile rectangulaire "})
    liste_exo.append({'ch':4,'exo':'Exercice4_7','com':"Profil en théorie linéarisée : profil NACA"})
    
    liste_exo.append({'ch':5,'exo':'Exercice5_1','com':"Aile optimale elliptique "})
    liste_exo.append({'ch':5,'exo':'Exercice5_2','com':"Vrillage d'une aile elliptique"})
    liste_exo.append({'ch':5,'exo':'Exercice5_3','com':"Braquage d'un volet d'une aile elliptique"})
    liste_exo.append({'ch':5,'exo':'Exercice5_4','com':"Calcul du foyer aérodynamique d'une aile"})
    liste_exo.append({'ch':5,'exo':'Exercice5_5','com':"Sillage tourbillonnaire d'un Airbus A380-800"})
    liste_exo.append({'ch':5,'exo':'Exercice5_6','com':"Corrections de soufflerie"})
    liste_exo.append({'ch':5,'exo':'Exercice5_7','com':"Tourbillons en fer à cheval"})
    liste_exo.append({'ch':5,'exo':'Exercice5_8','com':"Vol en formation"})
    liste_exo.append({'ch':5,'exo':'Exercice5_9','com':"Tourbillons de sillage au voisinage du sol"})
    liste_exo.append({'ch':5,'exo':'Exercice5_10','com':"Résolution numérique du problème de Prandtl, cas général"})
    liste_exo.append({'ch':5,'exo':'Exercice5_11','com':"Résolution numérique du problème de Prandtl avec effet de sol"})
    liste_exo.append({'ch':5,'exo':'Exercice5_12','com':"demi-aile d'un spitfire"})
    liste_exo.append({'ch':5,'exo':'Exercice5_14','com':"Chargement d'une aile trapézoïdale de P51 Mustang"})
#    liste_exo.append({'ch':5,'exo':'Exercice5_13','com':"Résolution numérique du problème de Prandtl, cas général (OLD)"})

    liste_exo.append({'ch':6,'exo':'Exercice6_1','com':"Loi du coefficient de portance avec l'allongement "})
    liste_exo.append({'ch':6,'exo':'Exercice6_3','com':"Ellipsoide"})
    liste_exo.append({'ch':6,'exo':'Exercice6_4','com':"Analyse des données de la thèse de Van Driest sur l'ellipsoide"})

    liste_exo.append({'ch':7,'exo':'Exercice7_0','com':"calcul des paramètres du profil de Van-Hooren {7.1.1}"})
    liste_exo.append({'ch':7,'exo':'Exercice7_1','com':"Profil de Van-Hooren, carte des paramètres dans le plan k-t (épaisseur) {7.1.1}"})
    liste_exo.append({'ch':7,'exo':'Exercice7_2','com':"Profil parabolique {7.1.2}"})
    liste_exo.append({'ch':7,'exo':'Exercice7_3','com':"Discrete Vortex Method on the parabolic mean line airfoil {7.2.1}"})
    liste_exo.append({'ch':7,'exo':'Exercice7_4','com':"Profil de Van-Hooren, distribution de sources, Approche de Neumann {7.2.2}"})
    liste_exo.append({'ch':7,'exo':'Exercice7_5','com':"Distribution uniforme de sources et de doublets sur le profil de Van-Hooren {7.2.3} "})
    liste_exo.append({'ch':7,'exo':'Exercice7_10','com':" Profil NACA : dessin"})
    liste_exo.append({'ch':7,'exo':'Exercice7_30','com':"Solution de référence sur un profil Naca depuis fichiers NACA"})

    liste_exo.append({'ch':8,'exo':'Exercice8_1','com':"Total pressure / Stagnation pressure = Isentropic Pressure {8.1.1}"})
    liste_exo.append({'ch':8,'exo':'Exercice8_2','com':"Pitot tube in subsonic flow {8.1.2} "})
    liste_exo.append({'ch':8,'exo':'Exercice8_3','com':"Compressibility corrections (airfoil) {8.2.1}"})
    liste_exo.append({'ch':8,'exo':'Exercice8_4','com':"Compressibility corrections (wing) and critical Mach number {8.2.2}"})
    liste_exo.append({'ch':8,'exo':'Exercice8_5','com':"paroi ondulée {8.3}"})
    liste_exo.append({'ch':8,'exo':'Exercice8_6','com':"paroi ondulée en transsonique {8.x} "})
    liste_exo.append({'ch':8,'exo':'Exercice8_7','com':"Calcul de Kpinc en fonction du Mach critique {8.x}"})

    liste_exo.append({'ch':9,'exo':'Exercice9_1','com':"Théorie linéaire : plaque plane en incidence {9.2.1} "})
    liste_exo.append({'ch':9,'exo':'Exercice9_2','com':"Théorie linéaire :  profil losangique {9.2.3}"})
    liste_exo.append({'ch':9,'exo':'Exercice9_3','com':"Théorie linéaire : profil de traînée minimale {9.2.4} "})
    liste_exo.append({'ch':9,'exo':'Exercice9_4','com':"Théorie linéaire : profil  lenticulaire {9.2.5}"})
    liste_exo.append({'ch':9,'exo':'Exercice9_5','com':"Interaction jet supersonique - profil {9.2.6}"})
   
    liste_exo.append({'ch':10,'exo':'Exercice10_1','com':"Poussée d'une tuyère de fusée"})
    liste_exo.append({'ch':10,'exo':'Exercice10_2','com':"Poussée d'une tuyère de fusée, formules du Mattingly."})
    liste_exo.append({'ch':10,'exo':'Exercice10_3','com':"Tube de Pitot dans une tuyère amorcée"})
    liste_exo.append({'ch':10,'exo':'Exercice10_4','com':"Moteur oxygène - hydrogène"})
    
    liste_exo.append({'ch':11,'exo':'Exercice11_1','com':"courbes p_1/p_0 en fonction de la déviation"})
    liste_exo.append({'ch':11,'exo':'Exercice11_2','com':"Réflexion d'un choc dans un canal"})
    liste_exo.append({'ch':11,'exo':'Exercice11_3','com':"Interaction de 2 chocs obliques "})

    liste_exo.append({'ch':12,'exo':'Exercice12_1','com':"Canal divergent"})
    liste_exo.append({'ch':12,'exo':'Exercice12_2','com':"Tube à choc"})
    liste_exo.append({'ch':12,'exo':'Exercice12_3','com':"Tuyère d'une fusée"})
    liste_exo.append({'ch':12,'exo':'Exercice12_4','com':"Tuyère de longueur minimale"})

    
    liste_exo.append({'ch':13,'exo':'Exercice13_1','com':"Obstacle conique en supersonique "})
    liste_exo.append({'ch':13,'exo':'Exercice13_2','com':"ogive parabolique en supersonique"})
    liste_exo.append({'ch':13,'exo':'Exercice13_3','com':"Traînée minimale d'un corps élancé  supersonique "})
    liste_exo.append({'ch':13,'exo':'Exercice13_4','com':"Tracé des iso-potentielles et lignes de courant "})
    liste_exo.append({'ch':13,'exo':'Exercice13_5','com':"Corps conique élancé "})

    liste_exo.append({'ch':14,'exo':'Exercice14_0','com':"Calcul du Mach aval pour un choc oblique sur un dièdre par différentes approximations "})
    liste_exo.append({'ch':14,'exo':'Exercice14_1','com':"profil losangique et triangulaire, méthode de Newton"})
    #liste_exo.append({'ch':14,'exo':'Exercice14_2','com':""})
    liste_exo.append({'ch':14,'exo':'Exercice14_3','com':"Kp pour des écoulements de dièdre 2D ou coniques, pour un nombre de Mach infini et d'autres "})
    liste_exo.append({'ch':14,'exo':'Exercice14_4','com':"Comparaison des méthodes de calcul autour d'un obstacle pointus"})
    liste_exo.append({'ch':14,'exo':'Exercice14_5','com':"forme d'un corps de traînée minimale en hypersonique "})
    liste_exo.append({'ch':14,'exo':'Exercice14_6','com':"test de l'approximation sur la fonction de Prandtl-Meyer en hypersonique "})
    liste_exo.append({'ch':13,'exo':'Exercice14_100','com':"Courbes des chocs coniques  "})
    liste_exo.append({'ch':14,'exo':'test_approx_omega','com':"test de l'approximation sur la fonction de Prandtl-Meyer en hypersonique"})
    liste_exo.append({'ch':14,'exo':'Exercice14_102','com':"Test de l'inversion de la fonction omega(M) en supersonique"}) 
     
    liste_exo.append({'ch':15,'exo':'Exo15_Lam_1','com':"Plaque plane : Loi en sinus et polynome "})
    liste_exo.append({'ch':15,'exo':'Exo15_Lam_2','com':"Plaque plane : Loi polynomiale"})
    liste_exo.append({'ch':15,'exo':'Exo15_Lam_4','com':"Mise en mouvement d'une plaque, problème de Rayleigh"})
    liste_exo.append({'ch':15,'exo':'Exo15_Lam_6','com':"Ecoulement de Couette instationnaire"})
    liste_exo.append({'ch':15,'exo':'Exo15_Lam_7','com':"Ecoulement instationnaire : second problème de Stokes"})
    liste_exo.append({'ch':15,'exo':'Exo15_Lam_8','com':"Jet laminaire"})
    liste_exo.append({'ch':15,'exo':'Exo15_Tra_1','com':"Calcul de la position de la transition"})
    liste_exo.append({'ch':15,'exo':'Exo15_Tur_1','com':"Ancienne version, mais fonctionne, pas de Class"})
    liste_exo.append({'ch':15,'exo':'Exo15_Tur_Aspiration','com':"profil turbulent analytique avec aspiration"})
    liste_exo.append({'ch':15,'exo':'Exo15_Tur_Loi_Puissance','com':"profil turbulent analytique avec une loi puissance"}) 
    liste_exo.append({'ch':15,'exo':'Exo15_Tur_VanDriest','com':"profil turbulent analytique avec la correction de Van Driest"})     
    
    if display:
        print("="*90)
        print("k  chap.     exercice   \t commentaires")
        print("="*90)

        k=0
        for dico in liste_exo:
            #for key,values in dico.items():
            print("%2i   %2i %20s  %s"%(k,dico['ch'],dico['exo'],dico['com']))
            k=k+1

        print("Pour choisir le ou les exercices, il faut uniquement entrer l'index de la première colonne")
        print("dans le fichier main.py")
    return liste_exo



def print_exo(chap,nbre):
    """
    format pour le numéro de l'exercice.
    """
    ch=" "*10
    print("_"*75,'\n')
    print('%s \tCh : %d \t EXERCICE n° %d'%(ch,chap,nbre))
    print("_"*75)


def correction(chapter,exos):
    """ Choose the exercices and the chapter 
        self : chapter list, exos : exercice list 
    """
    ch='*'*25
    i=chapter
    print('\n%s \tCHAPTER n° %d \t  %s\n'%(ch,i,ch))
    if i == 2:
        '''
        Modèles mathématiques pour l'aérodynamique 
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """
                Fonction de courant et potentiel des vitesses
                """
                Exercice2_1();
            if j == 2:
                """
                Atmosphère Standard
                """
                Exercice2_2();
    if i == 3:
        '''
        Théorie des potentiels complexes : applications aux profils d'aile
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """
                Cylindre
                """
                Exercice3_1();
            elif j == 2:
                """
                Transformation de Joukowski
                """
                Exercice3_2();
            elif j == 3:
                """
                Transformation de Karman-Trefftz
                """
                Exercice3_3();
            elif j == 4:
                """
                Transformation de Von Mises
                """
                Exercice3_4();
            elif j == 5:
                """
                Transformation de Van de Vooren et Jong
                """
                Exercice3_5();
            elif j == 6:
                """
                Profil à double courbure
                """
                Exercice3_6();
            elif j == 7:
                """
                Plaque plane en incidence
                """
                Exercice3_7();
            elif j == 8:
                """
                Plaque plane : introduction panneaux
                """
                Exercice3_8();

            elif j == 9:
                """
                Effet de sol
                """
                Exercice3_9();
            elif j == 12:
                """
                Ovale de Rankine
                """
                Exercice3_12();

    if i == 4:
        '''
        Théorie linéarisée des profils minces
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """
                Plaque plane en incidence
                """
                Exercice4_1();
            if j == 2:
                """
                profil de Joukowski
                """
                Exercice4_2();
            if j == 3:
                """
                profil famille_pq 
                """
                Exercice4_3();
            if j == 4:
                """
                profil double cambrure +  becs et volets 
                """
                Exercice4_4();
            if j == 5:
                """
                Lois de squelette et d'épaisseur
                """
                Exercice4_5();
            if j == 6:
                """
                Exercice sur les becs et volets, Aile rectangulaire
                """
                Exercice4_6();
            if j == 7:
                """
                Profil en théorie linéarisée : profil NACA
                """
                Exercice4_7();
            if j == 8:
                """
                Profil en théorie linéarisée : profil NACA
                """
                Exercice4_8();


    if i == 5:
        '''
        Théorie de la ligne portante de Lanchester-Prandtl
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """
                Aile optimale elliptique 
                """
                Exercice5_1();
            if j == 2:
                """
                Vrillage d'une aile elliptique
                """
                Exercice5_2();

            if j == 3:
                """
                Braquage d'un volet d'une aile elliptique
                """
                Exercice5_3();

            if j == 4:
                """
                Calcul du foyer aérodynamique d'une aile
                """
                Exercice5_4();

            if j == 5:
                """
                Sillage tourbillonnaire d'un Airbus A380-800
                """
                Exercice5_5();

            if j == 6:
                """
                Correction en soufflerie
                """
                Exercice5_6();

            if j == 7:
                """
                Tourbillon en fer à cheval
                """
                Exercice5_7();

            if j == 8:
                """
                Vol en formation
                """
                Exercice5_8();

            if j == 9:
                """
                Tourbillons de sillage au voisinage du sol
                """
                Exercice5_9();

            if j == 10:
                """
                Résolution numérique du problème de Prandtl, cas général
                """
                Exercice5_10();

            if j == 11:
                """
                Résolution numérique du problème de Prandtl avec effet de sol
                """
                Exercice5_11();

            if j == 12:
                """
                demi aile d'un spitfire
                """
                Exercice5_12();

            if j == 13:
                """
                aile elliptique sans effet de sol
                """
                Exercice5_13();
            
            if j == 14:
                """
                chargement de l'aile trapézoidale du P51 Mustang
                """
                Exercice5_14();

            if j == 15:
                """
                test de l'intégrale avec quad 
                """
                Exercice5_15();

    if i == 6:
        '''
        Théorie de la surface portante et des corps élancés
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """
                Loi du coefficient de portance avec l'allongement
                """
                Exercice6_1();
            if j == 3:
                """
                Ellipsoide
                """
                Exercice6_3();
            if j == 4:
                """
                Analyse des données de la thèse de Van Driest sur l'ellipsoide
                """
                Exercice6_4();
             

    if i == 7:
        '''
        Aspects numériques : méthode des singularités
        '''
        for j in exos:
            print_exo(i,j)
            if j == 0:
                """
                calcul des paramètres du profil de Van-Hooren {7.2.1} 

                """
                Exercice7_0();
            if j == 1:
                """ 
                Profil de Van-Hooren, carte des paramètres dans le plan k-t (épaisseur) {7.2.1} 
                """
                Exercice7_1();
            if j == 2:
                """ 
                Profil parabolique {7.2.2}
                """
                Exercice7_2();
            if j == 3:
                """
                Discrete Vortex Method on the parabolic mean line airfoil {7.3.1}
                """
                Exercice7_3();
            if j == 4:
                """ 
                Profil de Van-Hooren, distribution de sources, Approche de Neumann {7.3.2}
                """
                Exercice7_4();
            if j == 5:
                """ 
                Distribution uniforme de sources et de doublets sur le profil de Van-Hooren {7.3.3}
                """
                Exercice7_5();
            if j == 10:
                """ 
                Profil NACA : dessin
                """
                Exercice7_10();
            if j == 30:
                """ 
                Solution de référence sur un profil Naca depuis fichiers NACA
                dans misc.py
                """

                Exercice7_30();

    if i == 8:
        '''
        Ecoulements compressible subsonique et transsonique 
        '''
        for j in exos:
            print_exo(i,j)
           
            if j == 1:
                """ 
                Total pressure / Stagnation pressure = Isentropic Pressure {8.2.1}
                """
                Exercice8_1();
            if j == 2:
                """ 
                Pitot tube in subsonic flow {8.2.2}
                """
                Exercice8_2();
            if j == 3:
                """ 
                Compressibility corrections (airfoil) {8.3.1}
                """
                Exercice8_3();
            if j == 4:
                """ 
                Compressibility corrections (wing) and critical Mach number {8.3.2}
                """
                Exercice8_4();
            if j == 5:
                """ 
                paroi ondulée {8.4}
                """
                Exercice8_5();
            if j == 6:
                """ 
                paroi ondulée en transsonique {8.x}
                """
                Exercice8_6();
            if j == 7:
                """ 
                Calcul de Kpinc en fonction du Mach critique {8.x}
                """
                Exercice8_7();


    if i == 9:
        '''
        Ecoulements supersoniques linéarisés  
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """ 
                Théorie linéaire : plaque plane en incidence {9.2.1}
                """
                Exercice9_1();
            if j == 2:
                """ 
                Théorie linéaire :  profil losangique {9.2.3}
                """
                Exercice9_2();
            if j == 3:
                """ 
                Théorie linéaire : profil de traînée minimale {9.2.4}
                """
                Exercice9_3();
            if j == 4:
                """ 
                Théorie linéaire : profil  lenticulaire {9.2.5}
                """
                Exercice9_4();
            if j == 5:
                """ 
                Interaction jet supersonique - profil {9.2.6}
                """
                Exercice9_5();

    if i == 10:
        '''
        Ecoulements supersoniques monodimensionnels 
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """ 
                 Poussée d'une tuyère de fusée
                """
                Exercice10_1();
            if j == 2:
                """ 
                 Poussée d'une tuyère de fusée, formules du Mattingly.
                """
                Exercice10_2();
            if j == 3:
                """
                Tube de Pitot dans une tuyère amorcée
                """
                Exercice10_3();
            if j == 4:
                """
                Moteur oxygène - hydrogène
                """
                Exercice10_4();
    if i == 11:
        '''
        Ecoulements supersoniques bidimensionnels 
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """ 
                 courbes p_1/p_0 en fonction de la déviation
                """
                Exercice11_1();   
            if j == 2:
                """ 
                Réflexion d'un choc dans un canal
                """
                Exercice11_2();   
            if j == 3:
                """ 
                Interaction de 2 chocs obliques
                """
                Exercice11_3();        
                    
    if i == 12:
        '''
        Méthode des caractéristiques 
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """ 
                  Canal divergent
                """
                Exercice12_1();
            if j == 2:
                """ 
                  Tube à choc
                """
                Exercice12_2();
            if j == 3:
                """ 
                  Tuyère d'une fusée
                """
                Exercice12_3();
            if j == 4:
                """ 
                  Tuyère de longueur minimale
                """
                Exercice12_4();


    elif i == 13:
        '''
        Ecoulements supersoniques : corps élancés et ailes
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """ 
                Obstacle conique en supersonique 
                """
                Exercice13_1();
            if j == 2:
                """ 
                ogive parabolique en supersonique 
                """
                Exercice13_2();
            if j == 3:
                """ 
                Traînée minimale d'un corps élancé  supersonique 
                """
                Exercice13_3();
            if j == 4:
                """
                Tracé des iso-potentielles et lignes de courant
                """
                Exercice13_4()
            if j == 5:
                """
                Corps conique élancé
                """
                Exercice13_5()
    elif i == 14:
        '''
        Ecoulements hypersoniques
        '''
        for j in exos:
            print_exo(i,j)
            if j == 0:
                '''
                Calcul du Mach aval pour un choc oblique sur un dièdre par différentes approximations
                '''
                Exercice14_0();
            if j == 1:
                Exercice14_1();
                """
                profil losangique et triangulaire, méthode de Newton 
                """
            #if j == 2:         il n'existe pas, pas besoin de python
            #   Exercice14_2();
            if j == 3:
                Exercice14_3();
                '''
                Kp pour des écoulements de dièdre 2D ou coniques, pour un nombre de Mach infini et d'autres
                '''
            if j == 4:
                '''
                Comparaison des méthodes de calcul autour d'un obstacle pointus
                '''
                Exercice14_4();
            if j == 5:
                '''
                forme d'un corps de traînée minimale en hypersonique
                '''
                Exercice14_5();
            if j == 6:
                '''
                traînée d'un corps émoussé en hypersonique
                '''
                Exercice14_6();
            if j == 100:
                """ 
                Courbes des chocs coniques pour les rappels
                """
                Exercice13_100();
            
            if j == 101:
                """
                test de l'approximation sur la fonction de Prandtl-Meyer en hypersonique
                cas = 0 : plusieurs M0, cas =1 , 1 M0, cas = 2: fonction omega
                """  
                test_approx_omega(cas);

            if j == 102:
                """
                Test de l'inversion de la fonction omega(M) en supersonique
                """
                Exercice14_102();
    elif i == 15:
        '''
        Effets visqueux et couche limite 
        '''
        for j in exos:
            print_exo(i,j)
            if j == 1:
                """ 
                Plaque plane : Loi en sinus et polynome
                cet correction va plus loin que le sujet car on va jusqu'à des polynomes de degrés 7.
                """
                Exo15_Lam_1();
            if j == 2:
                """ 
                Plaque plane : Loi polynomiale
                """
                Exo15_Lam_2();
            if j == 4:
                """ 
                Mise en mouvement d'une plaque, problème de Rayleigh
                """
                Exo15_Lam_4();
            if j == 6:
                """ 
                Ecoulement de Couette instationnaire
                """
                Exo15_Lam_6();
            if j == 7:
                """ 
                Ecoulement instationnaire : second problème de Stokes
                """
                Exo15_Lam_7();
            if j == 8:
                """ 
                Ecoulement stationnaire : jet laminaire
                """
                Exo15_Lam_8();
            if j == 9:
                """
                Calcul de la position de la transition
                """
                Exo15_Tra_1();
            if j == 10:
                '''
                profil turbulent avec aspiration (ancienne version sans classe)
                '''
                Exo15_Tur_1();
            if j == 11:
                """
                profil turbulent analytique avec aspiration
                Utilisation de la classe Turbulence
                """
                Exo15_Tur_Aspiration();
            if j == 12:
                """
                profil turbulent analytique avec une loi puissance
                Utilisation de la classe Turbulence
                """
                Exo15_Tur_Loi_Puissance();
            if j == 13:
                """
                profil turbulent analytique avec la correction de Van Driest
                Utilisation de la classe Turbulence
                """
                Exo15_Tur_VanDriest();
            if j == 14:
                """
                profil turbulent analytique avec la correction de Van Driest
                tau fonction de y/delta
                Utilisation de la classe Turbulence
                """
                Exo15_Tur_tau();
            if j == 15:
                """
                profil turbulent analytique avec aspiration
                Utilisation de la classe Turbulence
                """
                Exo15_Tur_Aspiration_Test();



                 

def solution(liste):
    """
    Extraction des exercices
    """
    for L in liste:
        #print('chapter = %i'%(L[0]))
        #for exo in L[1]:
            #print("\t exercice = %i"%exo)
        correction(L[0],L[1])

