@author  C. Airiau
@date    Avril 2018

x3.clean.sh est le script shell qui permet d'effacer les fichiers inutiles

#***********************************************
COMPILATION DU CODE EC (Ecoulement Compressible)
#***********************************************

le code en fortan se trouve dans src_F90.
ce dossier ne soit servir que pour la compilation

    - faire dans ce dossier, dans un terminal 
      make clean; make
    -lancer le script x3.compil_EC.sh

  Le makefile compile le code et met le code EC dans le dossier ~/bin
  on peut changer la position dans le fichier Makefile qui se trouve dans src_F90

#***********************************************
LANCER LE CODE EC
#***********************************************

On peut le lancer dans n'importe quel endroit, 
mais il est conseillé de le faire dans le dossier run,

ce dossier contient de plus des fichiers d'entrée du code (*.in) :
    - cas.in             : pour définir des configurations personnelles (optionnel)
    - EC.in              : juste mettre le numéro de l'option           (optionnel)
    - choix_exercices.in : mettre la liste des exercices qu'on souhaite corriger 

Les fichiers de sortie portent tous l'extension .out

regarder EC.out par exemple après chaque calcul pour vérifier qu'il n'y a pas eu de problème.

#***********************************************
Version en MATLAB
#***********************************************

il existe une version MATLAB de cet outil.
- elle ne contient pas la solution des exercices du livre
- elle permet de tracer facilement un certain nombre de figures
- elle permet de faire du calcul en interactif dans la console MATLAB.
- elle se trouve dans le dossier Matlab.
- le script principal se nomme ecoulement_compressible.m
- on peut utiliser les fonctions indépendamment, mais bien donner la valeur
 à la variable gam (gamma) avant.
