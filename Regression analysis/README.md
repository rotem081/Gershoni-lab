# Regression analysis 

files included: 
1. bulls_paramaters.xlsx
2. test_least_square.py

## Description:

In this repository we arranged the procedure we used to predict fertility effects from measured sperm quality parameters.

In this part the sperm quality parameters that were measured by Flow-cytometry, were used to assess the correlation 
with fertility traits. This relationship was analyzed by multiple regression with partial least squares (PLSs).

In the PLSs method, the measured effect value was compared to the predicted effect value, which is the linear sum of the parameters,
each multiplied by a coefficient. The distance between the measured and predicted values is calculated using least square function. 

Our goal is to find the best fit, which is the combination of coefficients that minimizes the sum square residuals.
Since we have numerous coefficients, we need to use an efficient method to go through all possible combinations of the coefficients.
This was done using dual_annealing function 

We provide a python script to read parameters from excel file, normalized them with numpy, 
find minimum and maximum with dual_annealing and calculate the distance between the measured and 
predicted values is calculated using least square function
this will enable you to reproduce the PLS figures from the work

## Environment:
#### 1. Using Anaconda recipe:
download python 3.7 from https://www.python.org/downloads/
<br> download Anaconda from https://www.anaconda.com/products/individual and launch Spyder editor (python 3.7)

#### 2. Using pip:
* make sure that your virtual environment is installed with python 3.7 or more advaced version
in the anaconda platform, click on environments > base(root) > open Terminal.
in Terminal window download the function we used; scipy.optimize, numpy, IPython, pandas and os 
```
$ pip install <function_name> 
```

## Datadase:
download the database bulls_paramaters.xlsx

## Script:
open test_least_square.py and run the script.
<br> first we import all requierd functions

```
from scipy.optimize import least_squares,dual_annealing,shgo
import numpy as np
from IPython.display import display
import pandas as pd
import os
```

we defiend two functions: Normalize and least_square.
Normalize get an array of numbers and return a new array with zero mean and std=1 
<br> this funciton is useful because we need that all parameters have the same mean and std,
for the least_square function.

```
def Normalize(x):
    mu = np.mean(x) 
    sigma = np.std(x) 
    return (x-mu)/sigma
```

least_square function get array of coeffecients, the function load the data from existing 
verible in the envieroment and return the distance between the predicted and measured effect

```
def least_square(x):
    # data : each row is a bull, each column is a parmeter
    data = np.array(normalized_parameters)
    # effect corresponding to the bulls in data. first effect belong to bull
    # in first row and so on
    effect = np.array(df_effect) ***why not normalized?
    
    # dot is matrix multplication data*X = predicted effect.
    # norm calculte the size of the total distance between predicted and measuerd effect
    return np.linalg.norm(np.dot(data,x)-effect) 
```  

using pd.read_excel to turn data from excel to dataframe

```
DIR = r'{path to the excel file}' # fill directory hear
# read excel file
test = pd.read_excel(os.path.join(DIR, 'bulls_paramaters.xlsx'), index_col=None)
bulls = test[:25]
```

normalize the effect values and parameters, by calling the Normalize function

```
# the effect
df_effect = bulls['effect'].to_numpy()
normalized_effect = Normalize(df_effect)
# choosen parameters
df_parameters = bulls[['viable %', 'Depolarized  %','Viable spz ROS+',
                   'Viable spz ROS-','Viable, intact acrosome %',
                   'Viable, disrupted acrosome %',
                   'Dead, disrupted acrosome %']].to_numpy()
normalized_parameters = Normalize(df_parameters)
```

we set bounds for coeffienct (min, max), we picked 20-(-20) after trial and error.
we used the function dual_annealing to find minimum in the function least square 
(distance between predicted and measured effect).
dual_annealing get function, bounds and seed (used for randomization) and return 
the minimum of the function
```
bounds = np.ones((7,2))*20
bounds[:,0] *= -1 #the second column become -x
# dual_annealing find min in a function that it recive
result = dual_annealing(least_square,bounds,seed = 42)
```

then we calculate the R squre of the results, the set of coeffesient that gave the min of least square.
<br>R squre calculation was taken from https://en.wikipedia.org/wiki/Coefficient_of_determination

```
# the effect is defined here again to calculate r^2
effect = np.array(normalized_effect)
ss_tot = sum((effect-np.mean(effect))**2)
ss_res = least_square(result.x)**2
r2 = 1 - ss_res/ss_tot
```
