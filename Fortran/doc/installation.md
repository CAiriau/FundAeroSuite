
# CAT tools installation 
 
 

The Fortran source are found in the github repository in the **Fortran/CAT** directory
 
 

* The sources are in the **CAT/src** directory and it is adviced (even not necessary) to run the code the **CAT/run** directory
* The Fortran codes can be compiled on any system (Windows, MacOs, Linus) as soon as a Fortran Compiler is installed
* A makefile has to be used for compilation.

\warning When a source code is modified in any module, please delete the `*.mod` before compiling again 
 and in a Shell terminal write the commands:

\code{.sh} make clean; make \endcode
