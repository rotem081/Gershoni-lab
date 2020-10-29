# varients filtretion of 5 bulls (family 4)

list of files:
1. create gvcf_family4.txt (discussed in 'WGS anaysis')
2. family4_comb.txt (discussed in 'WGS anaysis')
3. familt4_genotype.txt (discussed in 'WGS anaysis')
4. vep_family4.txt (discussed in 'WGS anaysis')
5. multisample_family4.py
6. merge_pkl_family4.py
7. family4_analysis.py
8. GTEx.V8.Testis-specificity.r-scores.csv

## Description:
In this repository we arranged the genomic varients analysis and filtration
we provide a dedicated python script to filter and assess the varients that was found in the WGS analysis. 

## Environment:
#### 1. Using mobaXterm and anaconda:
Make sure that your virtual environment is installed with python 3.7 or more advaced version
Download the function we used; numpy, pandas, PyVCF, itertools,datetime, time, os

## Database:
Download csv file with list of genes that expressed in the testis and the level of expression: 
GTEx.V8.Testis-specificity.r-scores.csv.

## Script: (create_gvcf.txt)
Open multisample_family4.py and run the script.
first we import all requierd functions

```
import numpy as np
import pandas as pd
import vcf
from itertools import islice
import datetime
import time
import os
```

Define dget_data_from_record function: the function get one record form the multi-sample varients file that created in
'WGS analysis' section, and retrun a dataframe

```
def get_data_from_record(record):
    data = [record.CHROM,
                record.POS,
                record.REF,
                record.ALT[0].sequence,
                record.QUAL,
                record.INFO['CSQ'][0].split('|')[1],
                record.INFO['CSQ'],
                record.samples[0].gt_alleles,
                record.samples[1].gt_alleles,
                record.samples[2].gt_alleles,
                record.samples[3].gt_alleles,
                record.samples[4].gt_alleles]

    return data
```

Defiene the dataframe's columns  

```
columns = ['Chromosome',
           'Position',
           'Reference',
           'Alternate',
           'Quality',
           'Consequences',
           'mutation',
           'GT:Jey-Jey',
           'GT:Jhey',
           'GT:artist',
           'GT:garden',
           'GT:jermin']
           
df = pd.DataFrame(columns=columns)
```

We read multi-sample varients file 
and print the number of samples (5 bulls)

```
reader = vcf.Reader(filename=r'path\variant_effect_output.vcf.gz', compressed=False)
print(len(reader.samples)) #samples count

```

Then, we update the dataframe using loop. 'get_data_from_record' get 1000 records every time,
return data that will be added to the dataframe.
we print the time left for the end.
when 100,000 rows passes, we convert df to pickle file with the simple function 'df.to_pickle'
and reset a new dataframe.
in the end, we get 137 pickle files, each with 100,000 records.


```
counter = 1 #update rows that done
start = time.time() #start time
file_counter = 1 #count the files

#if there is an error , the file will saved
try:
    while True:         
        #islice: get next 1000 items/rows from reader. if there isnt, he set none
        #map: get array and set the function on it
        df = df.append(pd.DataFrame(map(get_data_from_record,islice(reader,1000)),columns=columns)) 
        
        if len(df) == 0:
            break
        
        #time now -time start in seconds
        time_passed = time.time()-start
        time_left = time_passed/counter/1000*20550107-time_passed  #to know number of rows- linux: wc -l
        print(f'{counter*1000} done. Time elapse : {str(datetime.timedelta(seconds=np.floor(time_passed)))}. ETA: {str(datetime.timedelta(seconds=np.floor(time_left)))}')
        #ETA: estimate time of arrival
        counter +=1
        if np.mod(counter,100) == 0: #pass 100,000 rows
            df.to_pickle(f'result{file_counter:03}.pkl') #save file           
            print(f'saved progress in file number {file_counter:03}')
            df = pd.DataFrame(columns=columns) #update new dataframe
            file_counter +=1
except:
    df.to_pickle(f'result{file_counter:03}.pkl')
    raise  #raise the error
```

open merge_pkl_family4.py
first we import all requierd functions

```
import numpy as np
import glob 
import pandas as pd
files = glob.glob('*.pkl') #find all the files type pkl
```
We read the files, one after the other, and convert it to dataframe
We update the dataframe and print a massage about the progress 
In the end, we create 14 merged pickle files, contain 10 pickle files (that were created in the previous script)

```
df = pd.read_pickle(files[0])
columns = df.columns
file_counter=1
for i in range(1,len(files)):
    print(f'finised {i}') #how many finished
    df = df.append(pd.read_pickle(files[i])) #mearge all
    if np.mod(i,10) == 0: #pass 10 pickles
        df.to_pickle(f'merge{file_counter:03}.pkl')            
        print(f'saved progress in merge {file_counter:03}')
        df = pd.DataFrame(columns=columns) #update new dataframe
        file_counter +=1

df.to_pickle(f'merge{file_counter:03}.pkl') #save file           
print(f'saved progress in merge {file_counter:03}') # :03 three digits
```

Open family4_analysis.py
first we import all requierd functions

```
import glob 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import re
```

We defenie allele list, mutation_list, list of gene expressed in testis,
the exome database (in dataframe) and several short function that will help us in the analysis

```
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
```

We create a loop that pass on each merge pickle file, extract importent info to a dataframe.
We filterd low quality reads (less than <30)
We merged data from exome database to our new dataframe: Allele frequency,Allele number and Inbreeding.
We filtered variatent that have this pattern of gynotype: 'healty' bull is homozygot and the other bulls are hetrozygot (df_garden_homo)
We filtered variatent that have AF > 0.5
We added a column- name of effect/consequences of the varient ('mutation')

```
#family 4 dataframe
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
    
    #without AF set 0
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
```    

we added a column with gene name.
we added gene-expresstion in the testis.

```   
df_mut['mutation_gene_name'] = df_mut['mutation'].apply(gene_name)
df_mut = pd.merge(df_mut, df_testis[['genes','r-testis']], how='left',left_on=['mutation_gene_name'],right_on = ['genes'])
df_mut.to_excel(f'all_mut_merge{num:03}.xlsx')

```

finaly, we got excel file with list of 2855 genes, with gene expression data (range 0-1).
we can arrange the r-testis values from high to the small.
