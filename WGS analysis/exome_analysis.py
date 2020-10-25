# -*- coding: utf-8 -*-
"""
Created on Wed May 20 09:11:38 2020

@author: Rotem
"""

import glob 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import matplotlib
import re

font = {'family' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

files = glob.glob('merge*_exome.pkl') #find all the files type pkl
df = pd.read_pickle(files[0])
columns = df.columns
df_exome = pd.DataFrame(columns=columns)

for num in range(1,len(files)+1):
    file_path = fr'C:\Users\Rotem\Desktop\lab\gatk\exome\merge{num:03}_exome.pkl'
    df= pd.read_pickle(file_path)
    #filter quality and samples count per row
    df = df.where(df['samples count'] >= 10) #SNV that have more than 10 bulls
    df = df.where(df['Quality'] > 30) #quality >30
    #merge exome and family 4
    df_exome = df_exome.append(df)


#samples count
test = df_exome.hist('samples count', bins = range(0,120),rwidth = 0.8)
ax = test[0][0]
ax.set_xticks(range(0,120,5))
ax.set_xlabel('number of bulls with genotype')
ax.set_ylabel('number of variants(1X10^6)')

#AF after filtring
af_test = df_exome.hist('Allele frequency', bins = np.linspace(0,1,100),rwidth = 0.8)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,1,21))
ax.set_xlabel('Allele frequency')
ax.set_ylabel('number of variants(1X10^6)')

#inbreeding after filtring
np.random.seed(1)
data = np.round(np.random.normal(5, 2, 100))
plt.hist(data, bins=10, range=(0,10), edgecolor='black')
plt.show()


af_test = df_exome.hist('Inbreeding', bins = np.linspace(0,0.5,200),rwidth = 0.8)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,0.5,10))
ax.set_xlabel('Inbreeding')
ax.set_ylabel('number of variants(1X10^6)')


df_exome = df_exome.where(df_exome['Inbreeding'] > 0).dropna()
#reggression AF/INBRE
plt.plot(df_exome['Allele frequency'], df_exome['Inbreeding'],".")
plt.xlabel('Allele frequency')
plt.ylabel('Inbreeding')
plt.title('AF vs. Inbreeding coefficients')


x = df_exome['Inbreeding']
y= df_exome['Allele frequency']
heatmap, xedges, yedges = np.histogram2d(x, y, bins=50)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
plt.title('Heatmap of 2D normally distributed data points')
plt.xlabel('Inbreeding')
plt.ylabel('Allele frequency')
plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.show()


#'Consequences' after filtring
counts= df['Consequences'].value_counts()
counts_csq = counts.where(counts > 50000).dropna()
counts_csq.plot.pie()

#mutations
def remove_space(lst):
    lst = lst[0].split('|')
    return ' '.join(lst).split()


def finding_missense(lst):
    if 'missense_variant' in lst:
        return True
 
def finding_del(lst):
    if any(['deleterious' in x for x in lst]):
        return True

def finding_high_del(lst):
    if any(['deleterious(0)' in x for x in lst]):
        return True 
    

df['mutation'] = df['mutation'].apply(remove_space)
#missense
df['missense'] = df['mutation'].apply(finding_missense)
df_miss_af = df.where(df['missense']==True).dropna()
df_miss_af.to_excel('missense.xlsx')
#deleterious
df['deleterious'] = df['mutation'].apply(finding_del)
df_del_af = df.where(df['deleterious']==True).dropna()
df_del_af.to_excel('deleterious.xlsx')

#deleterious
df_del_af['del_sever'] = df_del_af['mutation'].apply(finding_high_del)
df_del_af_sever = df_del_af.where(df_del_af['del_sever']==True).dropna()
df_del_af_sever['gene'] = df_del_af_sever['mutation'][]
df_del_af_sever.to_excel('deleterious_sever.xlsx')


af_test = df_miss_af.hist('Allele frequency', bins = np.linspace(0,1,100),rwidth = 0.8)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,1,21))
ax.set_xlabel('Allele frequency')
ax.set_ylabel('missense_variants')
ax.set_title('AF distribution of missense variates')

af_test = df_del_af.hist('Allele frequency', bins = np.linspace(0,1,100),rwidth = 0.8)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,1,21))
ax.set_xlabel('Allele frequency')
ax.set_ylabel('deleterious')
ax.set_title('AF distribution of deleterious')

RANGE = list(np.linspace(0,0.1,100))+list(np.linspace(0.1,1,50))
del RANGE[100]
af_test = df_del_af_sever.hist('Allele frequency', bins = np.linspace(0,1,100),rwidth = 0.9)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,1,21))
ax.set_xlabel('Allele frequency')
ax.set_ylabel('severe deleterious <0.5')
ax.set_title('AF distribution of severe deleterious')

####
def finding_var(lst):
    if 'deletion' in lst:
        return lst
    if 'insertion' in lst:
        return lst
    if 'sequence_alteration' in lst:
        return lst
    if 'indel' in lst:
        return lst


df_var = df['mutation'].apply(finding_var)

df_miss_af = df.where(df['missense']==True).dropna()