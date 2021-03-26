#!/bin/bash
# C. Airiau, april 2020.
echo -e "list of the python source files \n"

 
echo "=============================================="
echo "zeroth level directories:"
echo "=============================================="
echo
ls -CF *.py
echo "=============================================="
echo "first level directories:"
echo "=============================================="
echo
ls -CF  */*.py
echo "=============================================="
echo "second level directories:"
echo "=============================================="
echo
ls -CF  */*/*.py



