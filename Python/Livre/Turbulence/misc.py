#!/bin/python
# -*- coding: utf-8 -*-

#   FundAeroSuite
#   @date   : novembre 2018
#   @author : Christophe.airiau@imft.fr 

#*******************************************
def demande_valeur(Message,defaut=0,unit=""):
#*******************************************
    """
    Cette fonction permet de demander à l'utilisateur une valeur en affichant un message et 
    dans le cas où l'utilisateur ne rentre pas de valeur, de rentrer une valeur par défaut si cette
    valeur par défaut est renseignée
    """    
    try:
        val=float(input("\n"+Message+" :")) 
    except ValueError:
        if defaut != 0:   
            print("valeur par défaut")     
            val=defaut
            print("--> "+Message+" :",val,unit)
        else :
            print("valeur vide")
            val=""
    return val
