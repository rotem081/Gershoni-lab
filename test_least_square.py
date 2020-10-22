# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 21:43:55 2020

@author: Rotem
"""
from scipy.optimize import least_squares,dual_annealing,shgo
import numpy as np
from IPython.display import display
import pandas as pd
import os

def Normalize(x):
    mu = np.mean(x)
    sigma = np.std(x)
    return (x-mu)/sigma

def least_square(x):
    # data : each row is a bull, each column is a parmeter
    data = np.array(normalized_parameters)
    # effect corresponding to the bulls in data. first effect belong to bull
    # in first row and so on
    effect = np.array(df_effect)
    
    return np.linalg.norm(np.dot(data,x)-effect)

DIR = r'C:\Users\Rotem\Desktop\lab\Part A'
# read excel file
test = pd.read_excel(os.path.join(DIR, 'regression.xlsx'), index_col=None)
# 25 bulls 
bulls = test[:25]

# the effect
df_effect = bulls['effect'].to_numpy()
normalized_effect = Normalize(df_effect)
# choosen parameters
df_parameters = bulls[['viable %', 'Depolarized  %','Viable spz ROS+',
                   'Viable spz ROS-','Viable, intact acrosome %',
                   'Viable, disrupted acrosome %',
                   'Dead, disrupted acrosome %']].to_numpy()
normalized_parameters = Normalize(df_parameters)

# initial guess. Should have same number of elements as the number of 
# columns in data
x0 = np.array([0.4,0.3,0.35,0.15,1,1,0.1])
# 7 PARAMETERS and each have max and min(2) bounderis #based on manual chack
bounds = np.ones((7,2))*20
bounds[:,0] *= -1 #the second column become -1

# dual_annealing find min in a function that it recive
result = dual_annealing(least_square,bounds,seed = 42)

# the effect is defined here again to calculate r^2
effect = np.array(normalized_effect)
ss_tot = sum((effect-np.mean(effect))**2)
ss_res = least_square(result.x)**2
r2 = 1 - ss_res/ss_tot

print(result.x)

