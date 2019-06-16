#!/bin/bash

case $1 in
	c)echo "compilation"
	#gfortran fsk_continuation.f90 -o ~/bin/fsk.bin
	gfortran fsk_continuation.f90 -o fsk.bin
	echo 'done : $/bin/fsk.bin';;
	r) echo 'run '
		rm -f *.out *.dat
		./fsk.bin > fsk.out;;
	cl) echo 'clean the directory'
		rm -rf *.out *.tmp filetmp '1' *.*~ *.dat
		;;
	*)echo "mauvaise option : r ou c ou cl";;
esac
