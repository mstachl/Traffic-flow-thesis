# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 13:33:09 2017

@author: Markus
"""

from matplotlib import cm
from itertools import cycle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.optimize import minimize
import copy


def functional(x):
    return -(x[0]+0.5*x[1]+0.5*x[2])
    
bnds = ((0,1),(0,1),(0,1))
M=[[0,-1,-1],[-1,0,-1]]
cnts=[]
cnts2 = [{"type": "ineq", "fun": lambda x: 1-x[1]-x[2]},{"type": "ineq", "fun": lambda x: 1-x[0]-x[2]}]
#cnts2.append({"type": "ineq", "fun": lambda x: 0.5-x[2]})
for i in range(len(M)):
    def fun(x,i=i): 
        return 1+np.dot(M[i],x)
    cnts.append({"type": "ineq", "fun": fun})
cnts = tuple(cnts)
cnts2 = tuple(cnts2)
print cnts
print cnts2
x_init =[0,0,0]

res = minimize(functional,x_init,bounds=bnds,constraints=cnts)
res2 = minimize(functional,x_init,bounds=bnds,constraints=cnts2)