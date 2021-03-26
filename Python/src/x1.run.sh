#!/bin/bash
# 
#
#
#
# important parameters, user dependant 

# location of output file with fortran
FortranRun=$HOME/Documents/FORTRAN/EC/run
# location of output file with python
PythonRun=Run

echo "================================================================="
echo "shell script to run the main script with exercise in command line"
echo "          used for python or fortran solutions"
echo "          For linux, in a bash shell terminal"
echo "   caution : enter the right paths at the begining of file $0 "
echo "=================================================================="

# script arguments :
if [ $# -lt 4 ]; then
    echo -e " Usage: $0 -c x -e  y -l z "
    echo -e " Please respect the space after = "
    echo -e "  \t\t x: chapter number"
    echo -e "  \t\t y:  exercise number"
    echo -e "  \t\t z:  language f for Fortran, p for Python (default)"
    echo -e " the python solution is also found in /run/output.dat"
    exit
else

    if [ $# -eq 6 ] ; then
        case $6 in
            p) lang="p";;
            f) lang="f";;
            *) echo "Error : bad lang argument"
            exit ;;
        esac
    else
        lang="p"
    fi
fi

# Main tasks
file=tmp.py

cat << EOF > $file
from Book.Corrections import solution
from Book.Fortran.fortran import run_fortran
EOF

chapter=$2
exercise=$4

task="'compil+run'"
case $lang in 
    p) echo "solution(( ($chapter, ($exercise, ) ), ))" >> $file ;;
f) echo "run_fortran($task, ( ($chapter, ($exercise, ) ), ))" >> $file;;
*) echo "Error : bad lang argument"
    exit;
;;
esac

echo "print('normal end of execution from $0 ')" >> $file

case $lang in
    p) fileout=$PythonRun/Solution_$chapter-$exercise".out"
    python3 $file > $fileout
    ;;
    f) fileout=$FortranRun/Solution_$chapter"-"$exercise".out"
    python3 $file
    ;;
esac

echo "solution in $fileout"
ls -al $fileout
more $fileout
 

rm -f $file
