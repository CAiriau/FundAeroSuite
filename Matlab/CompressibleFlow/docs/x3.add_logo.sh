#!/bin/bash

dir_site=build/html
dir_html_static=source/_static/html
echo "script to add my logo on the documentation web site"
list_file=$(ls  $dir_site/*.html)

#echo "list of file to modify:"
#echo $list_file


# print the line
logo_file=$dir_html_static/include_logo.html

# to test on a given file : 
#file=$dir_html_static/introduction.html



for fichier in ${list_file}
do
    echo $fichier
    #sed -n "37,43 p" $fichier
    head -n 35 $fichier> header.txt
    sed -n "45,$ p" $fichier > tail.txt
    cat header.txt $logo_file tail.txt > out.html
    cp out.html $fichier 
done