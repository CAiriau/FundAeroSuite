#!/bin/bash

rm -f lignes.out
echo "le fichier est $1"
filename=$1
n=$2
i=0
for val in $(grep -n "theta" $filename  | cut -d':' -f1 )
do
    let i=$i+1
    table[$i]=$val
done
i1=$2; let i2=$i1+1
n=${#table[@]}
echo "${table[$i1]} : ${table[$i2]} "
l1=${table[$i1]};let l2=${table[$i2]}-1
sed -n "$l1,$l2 p" $filename > tmp

for ((i=1;i<$n;i++)) 
do
    i1=$i; let i2=$i1+1
    l1=${table[$i1]};let l2=${table[$i2]}-1
    let nLine=$l2-$l1
    echo $i $nLine >> lignes.out
done
