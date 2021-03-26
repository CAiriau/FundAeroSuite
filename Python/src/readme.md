**Comments**

* to run the solution of exercise, depending of user preference, on can use :

* bash shell script
* python script only
* additional bash shell script are useful for developments
* a directory named "Run" must be present at this level.

** listing **

* *main.py* is main script to use with python

* *run_gui.py* has to be used if the user wants a GUI

* *FundAeroSuite.ipnyb* can be used with jupyter notebook

* *x3.source_list.sh* provides the list of python files

* *x3.search.sh* provides the name of a file where is present a given
word.
    - Usage :  ./x3.search.sh "word_to_search"

* for a global substitution use the shell script *x3.replace.sh*, only for developers,
use with caution...

* *x1.run.sh* is a bash shell script to solve any exercise with a command line

* a python script *main_with_args.py* can do the same thing as x1.run.sh please read the file docstring for use or the documentation.