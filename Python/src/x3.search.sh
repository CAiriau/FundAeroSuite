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
    echo "zeroth level directories:"
    echo "=============================================="
    echo
    grep $option $1 *.py
    echo "=============================================="
    echo "first level directories:"
    echo "=============================================="
    echo
    grep $option $1 */*.py
    echo "=============================================="
    echo "second level directories:"
    echo "=============================================="
    echo
    grep $option $1 */*/*.py

fi


