# -*- coding: utf-8 -*-
"""
Created on Tue May 19 09:10:00 2020

@author: Boaz
"""
import numpy as np
import glob 
import pandas as pd
files = glob.glob('*.pkl') #find all the files type pkl

df = pd.read_pickle(files[0])
columns = df.columns
#if we want only the missence/deleterious mutations
#df = df[df['mutation'].str.contains('missense|deleterious')].dropna()
file_counter=1
for i in range(1,len(files)):
    print(f'finised {i}') #how many finished
    df = df.append(pd.read_pickle(files[i])) #mearge all
    if np.mod(i,10) == 0: #pass 10 pickles
        df.to_pickle(f'merge{file_counter:03}.pkl') #save file           
        print(f'saved progress in merge {file_counter:03}')
        df = pd.DataFrame(columns=columns) #update new dataframe
        file_counter +=1

df.to_pickle(f'merge{file_counter:03}.pkl') #save file           
print(f'saved progress in merge {file_counter:03}') # :03 three digits

