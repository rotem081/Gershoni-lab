# -*- coding: utf-8 -*-
"""
Created on Wed May 20 09:11:38 2020

@author: Rotem
"""
import glob 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import re
font = {'family' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)

global allele
allele=[['1','1'],['0','0']]

#mutations
global mutation_list
mutation_list = ['transcript_ablation','splice_acceptor_variant', 'splice_donor_variant',
                 'stop_gained','frameshift_variant','stop_lost','start_lost',
                 'transcript_amplification','inframe_insertion','inframe_deletion',
                 'missense_variant','protein_altering_variant','splice_region_variant',
                 'incomplete_terminal_codon_variant','start_retained_variant',
                 'stop_retained_variant','synonymous_variant','coding_sequence_variant',
                 'mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant']

#gene expression testis
file_path = r'C:\Users\Rotem\Desktop\lab\Part A\family 4 analysis\GTEx.V8.Testis-specificity.r-scores.csv'

#testis df
df_testis = pd.read_csv(file_path,sep='\t')
df_testis[['ENS','genes','r-testis']] =df_testis['Name,Description,r-testis'].str.split(',', expand=True)
df_testis.drop('Name,Description,r-testis', axis=1, inplace=True)

#df_testis_score=df_testis[['Gene','RNA tissue specificity score']]
#df_testis_score.to_excel('testis_score.xlsx')
#gene_testis_human = []
#gene_testis_human=df_testis['genes'].values.tolist()

#exome df
df_exome = pd.read_pickle(r'C:\Users\Rotem\Desktop\lab\gatk\merged_exome.pkl')
#filter <10 samples
df_exome = df_exome.where(df_exome['samples count'] > 10).dropna()

def remove_space(lst):
    lst = lst[0].split('|')
    return ' '.join(lst).split()

def test_garden(record):
    if record in allele:
        return True
    
def test_bulls(record):
    if record not in allele:
        return True

def testis(lst):
    return any([gene in lst for gene in df_testis['genes']])

def all_mut(lst):
    return any([mutation in lst for mutation in mutation_list])

def gene_name(lst):
    return lst[3]

#def score(gene):
#    return df_testis_score.loc[df_testis_score['Gene']==gene]

#family 4 data frame
files = glob.glob('merge*.pkl') #find all the files type pkl
df = pd.read_pickle(files[0])
columns = df.columns

df_mut = pd.DataFrame(columns=columns)

for num in range(1,len(files)+1):
    file_path = fr'C:\Users\Rotem\Desktop\lab\Part A\family 4 analysis\merge{num:03}.pkl'
    df= pd.read_pickle(file_path)
    #filter low quality reads
    df = df.where(df['Quality'] > 30).dropna()
    #merge exome and family 4
    df = pd.merge(df, df_exome[['Chromosome','Position','Alternate',
                                    'Allele frequency','Allele number','Inbreeding']], how='left',
                  left_on=['Chromosome','Position','Alternate'],
                  right_on = ['Chromosome','Position','Alternate'])
    
    #without AF 0
    df['Allele frequency'].fillna(0, inplace = True)
    
    #garden homo, the other bulls hetro
    df_garden_homo = df.where((df['GT:garden'].apply(test_garden)) &
                              (df['GT:Jey-Jey'].apply(test_bulls)) &
                              (df['GT:Jhey'].apply(test_bulls)) &
                              (df['GT:jermin'].apply(test_bulls)) &
                              (df['GT:artist'].apply(test_bulls))).dropna()
    #filter only AF < 0.5
    df_garden_homo = df_garden_homo[df_garden_homo['Allele frequency'] < 0.5]

    #filter mutations
    df_garden_homo['mutation'] = df_garden_homo['mutation'].apply(remove_space)
    df_garden_homo['all_mut'] = df_garden_homo['mutation'].apply(all_mut)
    df_mut = df_mut.append(df_garden_homo.where(df_garden_homo['all_mut']==True).dropna()) 
    #df_mut['testis'] = df_mut['mutation'].apply(testis)
df_mut['mutation_gene_name'] = df_mut['mutation'].apply(gene_name)
df_mut = pd.merge(df_mut, df_testis[['genes','r-testis']], how='left',left_on=['mutation_gene_name'],right_on = ['genes'])
    #df_mut['mutation_gene_name'].apply(score)
    #df_testis =  df_testis.append(df_mut.where(df_mut['testis']==True).dropna())

df_mut.to_excel(f'all_mut_merge{num:03}.xlsx')

