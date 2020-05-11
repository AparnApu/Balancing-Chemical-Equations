# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:50:29 2020

@author: apuap
"""
###########################################
# error check
    # invalid compound 
    # atom on left and right must match

# enhancements
    # bracket compounds
    # ionic compounds
    # multiple rxn support
    # doc strings multiline
###########################################

from RxnBalance_ClassModule import RxnBalance

###########################################
# user input
###########################################

#ob = RxnBalance(unbal_str = 'SeCl6 + O2 -> SeO2 + Cl2') #1113
ob = RxnBalance(unbal_str = 'NH3 + O2 -> NO + H2O')     #4546

###########################################
# strip input of whitespaces
###########################################

ob.strip_space()

###########################################
# to separate input into reactants & products
###########################################

ob.input_split()

###########################################
# to split elements in reactants & products
###########################################

ob.r_p_split()

###########################################
# to form matrix A
###########################################

ob.create_ListOfDict()
ob.input_check()
ob.fill_emptyKeys()
ob.matrix_formation()

###########################################
# to solve system of homogenous equations
###########################################

finAns = ob.main()

AnsStr = ob.display(finAns)

print(AnsStr)