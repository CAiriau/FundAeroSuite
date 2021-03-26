#!/bin/bash
# C. Airiau, april 2020.
echo "================================================================="
echo -e "replace a  words in the python file at a given level \n"
echo -e "\t sed tool is used here, can be done by terminal command"
echo "================================================================="



if [ $# -ne 6 ]; then
    echo -e " Usage: $0  o= 'word_in' n= 'new_word' l= n"
    echo -e " Please respect the space after = "
    echo -e "  \t\t o : old words to replace"
    echo -e "  \t\t n : new words"
    echo -e "  \t\t l : level (0,1, or 2)"
    echo -e "for a single word, ' ' are optional"
    echo -e "level directory substitution:"
    echo -e "      0 : present directory"
    echo -e "      1 : */*.py"
    echo -e "      2 : */*/*.py"
    echo -e "      a : all levels"
    echo -e "to use with CAUTION, \n since it won't be possible to come back in case of error"
    echo -e " for experimented users only! "
    exit
else

    case $6 in
    0) echo "level 0 substitutions"
        echo $2 $4
        sed -i  "s/$2/$4/g" *.py
        ;;
    1)  echo "level 1 substitutions"
        sed -i  "s/$2/$4/g" */*.py;;
    2) echo "level 2 substitutions"
        sed -i "s/$2/$4/g" */*/*.py;;
    a) echo "all levels substitutions"
        sed -i "s/$2/$4/g" *.py
        sed -i "s/$2/$4/g" */*.py
        sed -i "s/$2/$4/g" */*/*.py;;
    *) echo 'bad agument, no substitution';;
    esac
    echo "verification : old string chain"
    ./x3.search.sh $2

fi

