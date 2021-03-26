#!/bin/bash

rm -rf liste_fichiers_Mach.dat

i=0
for file in $(ls Mach*.dat)
do
    let i=i+1;
    MachFile[$i]=$file
    nlignes=$(wc -l $file | cut -d' ' -f1)
    let nlignes=$nlignes-1
    echo "${file:5:5}   $nlignes   $file" >> liste_fichiers_Mach.dat
done
