#!/bin/bash
# script in order to calculate  the Falkner-Skan solution in a 
# flat plate boundary layer
# several options, please read below 
# C. Airiau, 2010-2020
 

#****************
function help()
#****************
{
    echo "help of $0"
    echo -e " use : $0 args"
    echo -e "\t \t 0 :\t clean the directory"
    echo -e "\t \t c :\t Compilation fsk_module"
    echo -e "\t \t r :\t Run "
}

case $1 in
	c) echo "Compilation fsk_module"
	gfortran fsk_module.f90 -o ~/bin/fsk
	echo 'done : ~/bin/fsk';;
	r) echo " run fsk_module "
	fsk
	;;
	0) echo 'clean the directory'
	rm -rf *.out *.tmp filetmp '1' *.*~ *.dat *.mod
	;;
	*)echo "mauvaise option"
	help
	;;
esac
