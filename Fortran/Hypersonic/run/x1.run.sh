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
    echo -e "\t \t 1 :\t Run atmospheric reentry"
    echo -e "\t \t 2 :\t Run Kp hypersonic"
    echo -e "\t \t 3 :\t Run piston method"
}

case $1 in
    0) echo "clean the directory"
    clean
    ;;
    1) echo "run atmospheric reentry"
    atmospheric_reentry
	echo "new files : "; ls -al *.dat 
    ;;
    2) echo "Run Kp hypersonic"
    hypersonic_Kp
    echo "new files : "; ls -al *.dat 
    ;;
    3) echo "Run piston method"
    piston_method
    echo "new files : "; ls -al *.dat
    ;;
    *) echo "help"
        help
    ;;
esac

