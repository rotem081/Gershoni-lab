============================
README for: Regression analysis 
============================
files included: 
1. bulls_paramaters.xlsx
2. test_least_square.py

==================
=  Description:  =
==================

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
download Anaconda from https://www.anaconda.com/products/individual and Lanch speyer package (python 3.7)

#### 2. Using pip:
* make sure that your virtual environment is installed with python 3.7 or more advaced version
in the anaconda platform, click on environments > base(root) > open Terminal.
in Terminal window download the function we used; scipy.optimize, numpy, IPython, pandas and os 
```
$ pip install <function_name> 
```

## Datadase:
download the database found under the name 'bulls_paramaters.xlsx', to set the parameters

## Script:
open python script name test_least_square.py in spyder environment and run the script.
first we import all requierd functions:'
```
from scipy.optimize import least_squares,dual_annealing,shgo
import numpy as np
from IPython.display import display
import pandas as pd
import os
```

we defiend two functions Normalize and least_square. 
we must adjucting values measured on different scales to a notionally common scale. 
we used the standard score formula: paramater minus mean, diveded by standerd diviation.
we normalized the effect values and the quality parameters values.

Normalize function get one parameter, one sperm quality parameter taken from 'bulls_paramaters.xlsx',
and return the normalized paramater
```
def Normalize(x):
    mu = np.mean(x) 
    sigma = np.std(x) 
    return (x-mu)/sigma
```

Least_square function get x parameter- initial guess, and return  ***

```
def least_square(x): ***x??
    # data : each row is a bull, each column is a parmeter
    data = np.array(normalized_parameters)
    # effect corresponding to the bulls in data. first effect belong to bull
    # in first row and so on
    effect = np.array(df_effect) ***why not normalized?
  
    return np.linalg.norm(np.dot(data,x)-effect) ***dot??
```  

we used pd.read_excel to turn data from excel to Dataframe in python
```
DIR = r'C:\Users\Rotem\Desktop\lab\Part A'
# read excel file
test = pd.read_excel(os.path.join(DIR, 'regression.xlsx'), index_col=None)
# 25 bulls 
bulls = test[:25]
```

we normalized the effect values and speem quality parameters, using Normalize function
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

we set initial guess (x0)   why?
we set bounds for variables (min, max) pairs for the parameters- randomaly we picked 20-(-20) why?

```
# initial guess. Should have same number of elements as the number of 
# columns in data
x0 = np.array([0.4,0.3,0.35,0.15,1,1,0.1])
# 7 PARAMETERS and each have max and min(2) bounderis #based on manual chack
bounds = np.ones((7,2))*20
bounds[:,0] *= -1 #the second column become -x
```

we used the function dual_annealing to find minimum in the function least square
dual_annealing get function, bounds and seed and return the minimum of the function.
```
# dual_annealing find min in a function that it recive
result = dual_annealing(least_square,bounds,seed = 42)
```


```
# the effect is defined here again to calculate r^2
effect = np.array(normalized_effect)
ss_tot = sum((effect-np.mean(effect))**2)
ss_res = least_square(result.x)**2
r2 = 1 - ss_res/ss_tot

print(result.x)
```

======================
=  Special Comments:  =
======================