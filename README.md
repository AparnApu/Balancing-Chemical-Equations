# Balancing-Chemical-Equations

## Table of Contents

- [Introduction](#Introduction)
- [Motivation](#Motivation)
- [Working](#Working)
- [Contributions](#Contributions)
- [Acknowledgement](#Acknowledgement)

<!-- toc -->


## Introduction
This program balances your chemical reaction. It uses an optimizer based on Diophantine's equation solving.

## Motivation
I've always liked chemistry and coding, so this program was a fun project to undertake! But this program stretched my brain and coding skills quite a bit - from parsing the equation and its atoms to solving the system of homegenous linear equations using an optimizer. All in all, I thoroughly enjoyed the challenge because I learnt quite a lot from this project!

## Working
The basic principle is that I need to cast the equation into this form: Ax = 0 (homogenous linear system of equations that has an integer solution).  
  
Dimensions of A are:  
Number of rows = number of atoms  
Number of columns = total number of compounds (ie, number of reactants + number of products)  
To get the equation to that point, the code splits the compounds into the reactant side and product side. I then create a list of dictionaries for the multiplicity of the atoms in each compound. If atoms does not exist, the multiplicity is zero. 
For example,  
  
Unbalanced Reaction: NH<sub>3</sub> + O<sub>2</sub> -> NO + H<sub>2</sub>O  
List of dictionaries (on reactant side): [{'N': 1, 'H': 3, 'O': 0}, {'O': 2, 'N': 0, 'H': 0}]  
List of dictionaries (on product side): [{'N': 1, 'O': 1, 'H': 0}, {'H': 2, 'O': 1, 'N': 0}]

Once this is done, I form two matrices (reactant side and product side) like this:  
(Reactant side):  
H: [3, 0]   
N: [1, 0]  
O: [0, 2]  
  
(Product side):    
H: [0, 2]  
N: [1, 0]  
O: [1, 1]  

I then stack them side by side, with the product side matrix being negated.
So, I get my matrix A:  
[ 3,  0,  0, -2]  
[ 1,  0, -1,  0]  
[ 0,  2, -1, -1]  

I used **Differential Evolution** to get a solution x which satisfies the Ax = 0.  
The problem there is that the solution x isn't an integer. So, I find a value α for which α * x is an integer.


To run the code:  
Fork the repo (link: (https://github.com/AparnApu/Balancing-Chemical-Equations/fork) or download/ clone it. Run the file 'RxnBalance_wrapper.py'.
I have already included two test cases (lines 25 - 26), which you can modify. There are no packages to install.
  
### Code snippet to modify:

```
#ob = RxnBalance(unbal_str = 'SeCl6 + O2 -> SeO2 + Cl2') #1113
ob = RxnBalance(unbal_str = 'NH3 + O2 -> NO + H2O')     #4546
```


## Contributions
Open to contributions!
A possible enhancement I am yet to add are compounds that involve brackets.  
Fork the repo, edit it and commit your change.

## Acknowledgement
A huge thankyou to my dad who came up with the idea for this project and helped me navigate my way through this complicated code!

