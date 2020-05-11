# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:05:44 2020

@author: apuap
"""

import sys
import numpy as np
from scipy.optimize import differential_evolution


class RxnBalance:
    
   
    def __init__(self, unbal_str):
        
        self.unbal_str = unbal_str
        
        return
    
  
    def strip_space(self):
        
        self.unbal_str = "".join(self.unbal_str.split())
    
    
    def input_split(self):
        
        temp =  self.unbal_str.split('->')
        
        self.reac = temp[0]
        self.prod = temp[1]
        
    
    def r_p_split(self):
        
        self.comp_reac = self.reac.split('+')
        self.comp_prod = self.prod.split('+')
        self.Nr = len(self.comp_reac)
        self.Np = len(self.comp_prod)
    
    
    def dict_AtomMultiplicity(self, compound):
        
        string = compound

        p = len(string)
        
        Atoms = {} # Dictionary to store atom as key, and multiplicity as value
        num_string = 0
        currentAtom = ""
        
        for i in range(0, p): # parsing the string char by char
            
            if(string[i].isdigit()):
                
                num_string = string[i]
                
                for j in range(i+1, p):
                    
                    if(string[j].isdigit()):
                        
                        num_string += string[j] 
                    
                    else:
                    
                        break
        
                if not ( currentAtom in Atoms.keys()  ):       
                    Atoms[currentAtom] = int(num_string) ## key - value allocation ##
                else:
                    Atoms[currentAtom] = int(num_string) + Atoms[currentAtom]           
                
                currentAtom = ""
                
            elif(string[i].islower()):
                
                currentAtom = currentAtom + string[i]
                    
            elif(string[i].isupper()):
        
                if not ( currentAtom in Atoms.keys()  ):         
                    Atoms[currentAtom] = 1 # end of the existing key with value=1 (value explicitly give would go to digit tree) ## key - value allocation ##
                else:
                    Atoms[currentAtom] = int(num_string) + Atoms[currentAtom]           
                
                currentAtom = string[i] # start of a new key
                        
        if (currentAtom != ""):
                Atoms[currentAtom] = 1
    
        del Atoms[""]
        
        return Atoms
    
    
    def create_ListOfDict(self):
        
        self.ListofDict_R = []
        self.ListofDict_P = []
        temp_dict_R = {}
        temp_dict_P = {}
        
        for i in range(self.Nr):
            
            temp_dict_R = self.dict_AtomMultiplicity(self.comp_reac[i])
            self.ListofDict_R.append(temp_dict_R)
        
        for i in range(self.Np):
            
            temp_dict_P = self.dict_AtomMultiplicity(self.comp_prod[i])
            self.ListofDict_P.append(temp_dict_P)
   
    
    def create_allKeys(self, LoD):
        
        allK = {}
        tempSet ={}
        p = len(LoD)
        
        for i in range(p):
            
            tempSet = set(LoD[i].keys())
            allK = tempSet.union(allK)
        
        return allK
        
    
    def input_check(self):
        
        self.allKeys_P = {}
        self.allKeys_R = {}
        
        self.allKeys_R = self.create_allKeys(self.ListofDict_R)
        self.allKeys_P = self.create_allKeys(self.ListofDict_P)
        
        if not (self.allKeys_R == self.allKeys_P):
            
            print("Input entered is incorrect. Please re-enter data.\nProgram will exit now.")
            sys(exit)
        
        else:
            
            print("Input entered has been parsed, loading your input now.")
    
    def fill_emptyKeys(self):
        
        for i in range(len(self.ListofDict_R)):
            
            diffR = self.allKeys_R.difference(self.ListofDict_R[i])
            
            for x in diffR:
                self.ListofDict_R[i][x] = 0
                
        for i in range(len(self.ListofDict_P)):
            
            diffP = self.allKeys_R.difference(self.ListofDict_P[i])
            
            for x in diffP:
                self.ListofDict_P[i][x] = 0
        
        
        self.Na = len(self.allKeys_R)

    
    def matrix_formation(self):
        
        self.allKeys_R = list(self.allKeys_R)
        self.allKeys_R.sort()
        
        self.arr_R = np.zeros((self.Na,self.Nr), int)
        self.arr_P = np.zeros((self.Na,self.Np), int)

        for i in range(self.Na):
            
            for j in range(self.Nr):
                
                self.arr_R[i][j] = self.ListofDict_R[j][self.allKeys_R[i]]
        
        for i in range(self.Na):
            
            for j in range(self.Np):
                
                self.arr_P[i][j] = self.ListofDict_P[j][self.allKeys_R[i]]
                
        
        self.A = np.hstack((self.arr_R, self.arr_P*-1))
            
        
    def myobj(self, x, *refz):
        
        resd = np.dot(refz[0], x) - refz[1]
    
        err = np.linalg.norm(resd)
        
        return err
        
    def dec_to_int(self, arr):
        
        min_val = np.amin(arr)
    
        arr = arr / min_val
        
        x = 1
        
        while True:
            
            temp = arr * x
            
            delta = temp - np.round(temp)
            
            err= np.linalg.norm(delta)
            
            if (err<0.05):
                break
            else:
                x = x + 0.01
        
        res = np.round(temp).astype('int')
        
        return res
    
    def main(self):
        
        b = np.zeros(self.Na, int)

        bounds = [(1, 10)] * (self.Nr + self.Np)
        
        args = (self.A, b)
        
        self.result = differential_evolution(self.myobj, bounds, args=args)
        
        fin_result = self.dec_to_int(self.result.x)
        
        return(fin_result)
    
   
    def display(self, ans):
        
        self.coef_R = list(ans[0:self.Nr])
        self.coef_P = list(ans[self.Nr:len(ans)])
        
        finR = []
        finP = []
        
        for i in range(len(self.coef_R)):
            
            finR.append(str(self.coef_R[i]) + self.comp_reac[i])
            
        for i in range(len(self.coef_P)):
            
            finP.append(str(self.coef_P[i]) + self.comp_prod[i])
        
        finR = ' + '.join(finR)
        finP = ' + '.join(finP)
        
        self.finStr = finR + ' -> ' + finP
        
        return self.finStr
