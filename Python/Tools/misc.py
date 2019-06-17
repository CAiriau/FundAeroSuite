"""
 titles, plots, constantes ...
"""

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 

import matplotlib.pyplot as plt
import numpy as np
import sys
 
RED   = "\033[1;31m"  
BLUE  = "\033[1;34m"
CYAN  = "\033[1;36m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"
BOLD    = "\033[;1m"
REVERSE = "\033[;7m"
Largeur,Hauteur=10,8

def set_title(self):
    """Set a title of an exercice"""
    print('#',60*'*')
    print('# %s'%(self))
    print('#',60*'*','\n')

def set_question(self):
    """Set a question of an exercice"""
    print()
    print('#',50*'=')
    print('# Question n° %s'%(self))
    print('#',50*'=','\n')

def set_alert(self):
    """Set an alert message"""
    print()
    print('#',50*'=')
    sys.stdout.write(RED)
    print('# ALERTE :  %s'%(self))
    sys.stdout.write(RESET)
    print('#',50*'=','\n')

def set_info(self):
    """Set an alert message"""
    print()
    print('#',50*'=')
    sys.stdout.write(REVERSE + CYAN)
    print('# Alerte  %s'%(self))
    sys.stdout.write(RESET)
    print('#',50*'=','\n')


def SimplePlot(ifig,title,x,y,leg,n=1):
    """
    plot très simple et rapide
    """
    fig = plt.figure(ifig,figsize=(Largeur,Hauteur))
    fig.suptitle(title, fontsize=14, fontweight='bold')
    ax = fig.add_subplot(111)
    #fig.subplots_adjust(top=0.80)
    ax.set_ylabel(leg[1],fontsize=20)
    ax.set_xlabel(leg[0],fontsize=20)
    if n > 1:
        for k in range(n):
            ax.plot(x[k],y[k],linewidth=2)
    else:
        ax.plot(x,y,linewidth=2)
    ax.grid()

 