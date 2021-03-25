#!/bin/bash
# C. Airiau, april 2020.
echo -e "find words in files \n"

option="--color"
if [ -z $1 ]; then
    echo 'no argument'
    echo -e " Usage: $0 words "
    exit
else
    echo "=============================================="
    echo "Fortran files:"
    echo "=============================================="
    echo
    grep $option $1 *.f90
    echo "=============================================="
    echo "Makefile:"
    echo "=============================================="
    echo
    grep $option $1 Makefile

fi


