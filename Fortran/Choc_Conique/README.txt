Pour la compilation, aller dans le dossier src et faire make dans un terminal.

on peut le compiler aussi à la main par gfortran cone.f90 -o choc_conique.bin

Il faut le lancer avec le fichier d'entrée cone.in qui se trouve dans le dossier run.

Les calculs sont assez longs car j'ai pris beaucoup de points d'intégration pour avoir
une très bonne précision. Le problème est mathématiquement "raide" 
en particulier au voisinage du demi-angle de cône proche de zéro.

La base de données résultant de ces calculs se trouve dans le dossier Python.

C. Airiau, Avril 2018
