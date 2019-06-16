Pour le compiler il suffit de faire make
le programme va se positionner dans le dossier ~/bin et s'appelle 
chocs_obliques.

Il permet de tracer les polaires de Busemann.

On peut changer les options (nom, position du programme) dans le fichier makefile
ou bien écrire dans un terminal :

gfortran busemann.f90 -o busemann.bin

pour le lancer faire 
./busemann.bin par exemple

Il vaut mieux ne pas travailler dans le dossier des sources.

C. Airiau, Avril 2018.

Pour tracer les courbes il suffit d'utiliser gnuplot ou un autre logiciel.

si vous avez installé l'outil gp (script shell de tracer rapide que j'ai écris, utilisant en fait gnuplot) il suffirait d'écrire:

gp l lf 1 2 courbes*.dat

l'outil gp se trouve dans le dossier Shell
