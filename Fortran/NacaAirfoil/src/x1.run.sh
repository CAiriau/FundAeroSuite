#!/bin/bash
# script in order to calculate the NACA airfoil 4 and 5 digits
# several options, please read below 
# C. Airiau, 2015
#****************
function clean()
#****************
{
	echo "clean the directory"
    rm -f *.*~ *.bin *.dat *.out gp.tmp fort* filetmp  *.mod
	ls 
}

#****************
function help()
#****************
{
    echo "help of $0"
    echo -e " use : $0 args"
    echo -e "\t \t 0 :\t clean the directory"
    echo -e "\t \t 1 :\t Compilation of naca45.f90 "
    echo -e "\t \t 2 :\t Run naca45"
    echo -e "\t \t 3 :\t Plot  : y = f(x), airfoil"
    echo -e "\t \t 4 :\t Plot NACAyyyxx "
    echo -e "\t \t 5 :\t compile  naca.f90"
    echo -e "\t \t 6 :\t run  naca"
}

case $1 in
    0) echo "clean the directory"
    clean
    ;;
    1) echo "compilation of naca45"
    gfortran naca45_airfoil.f90  -o ~/bin/naca45
    ;;
    2) echo "run naca45" 
    naca45
	echo "new files : "; ls -al *.dat 
    ;;
    3) echo "plot Airfoils"
        gp  lf 1 2 *.dat *.DAT
        ;;
    4) echo "Airfoils NACAyyyxx"
        gp  l lf 1 2 NACA$2.dat Naca$2_REF.DAT
    ;;
    1) echo "compilation of naca.f90"
    gfortran naca.f90  -o ~/bin/naca
    ;;
    6) echo "run naca, ouput in naca.out"
    naca > naca.out
    echo "new files : "; ls -al *.dat 
    ;;
    *) echo "help"
        help
    ;;
esac

