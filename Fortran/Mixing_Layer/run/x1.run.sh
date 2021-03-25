#!/bin/bash
# C. Airiau, 2015

#***************
function clean()
#***************
{
echo "nettoyage de ce dossier"
rm -rf *.bin *.out fort.* beta* *.*~  *.dat profil.m
rm -rf filetmp gp.tmp
}

#****************
function aide()
#****************
{
    echo "help of script $0"
    echo -e " usage : $0 argument1 directory (optional, for option 5 only)"
    echo -e "\t \t 0 :\t clean the directory"
    echo -e "\t \t 1 :\t Run  ~/bin/mixing_layer_o1"
    echo -e "\t \t 2 :\t Plot : U  / Um(eta) "
    echo -e "\t \t 3 :\t Plot : V  / Um(eta)"
    echo -e "\t \t 4 :\t Plot : U' / Um(eta)"
    echo -e "\t \t 5 :\t Save data in directory $2"
    echo -e "\t \t 6 :\t Run ~/bin/mixing_layer (order 1 or 4, more parameters"
}

case $1 in
    0) echo "clean the directory"
    clean
    ;;
    1) echo "run  mixing_layer_o1"
    clean
    mixing_layer_o1 > res.out
    tail -n 9 res.out
    ;;
    2) echo "plot : U/Um(eta) "
        gp noleg l lf 2 4 CAS*/sortie.out 
    ;;
    3) echo "plot : V/Um(eta)"
        gp noleg l  sortie.out 2 8
    ;;
    4) echo "plot : U'/Um(eta)"
        gp noleg l  sortie.out 2 5
    ;;
    5) echo "Save data in $2"
        mkdir $2
        cp outputs.out $2/.
        cp negative.dat $2/.
        cp positive.dat $2/.
        cp res.out $2/.
        cp turbulent_mixing.in $2/.
        ls $2
    ;;
    6) echo "run mixing_layer"
    clean
    main_mixing > res.out
    tail -n 30 res.out
    ;;

    *) echo "help"
        aide
    ;;
esac

