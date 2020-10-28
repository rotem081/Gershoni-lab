# variants filtration of 5 bulls (family 4)

list of files:
1. create gvcf.txt (discussed in 'WGS anaysis')
2. ex_comb.txt (discussed in 'WGS anaysis')
3. ex_genotype.txt (discussed in 'WGS anaysis')
4. vep_exome.txt (discussed in 'WGS anaysis')
5. multisample_exome.py
6. merge_pkl_exome.py
7.exome_analysis.py

## Description:
In this repository we arranged the procedure of assessment of genomic variants.
we built dedicated python script, that uses the VCF parser (PyVCF), to filter and assess the variants. 

## Environment:
#### 1. Using mobaXterm and anaconda:
make sure that your virtual environment is installed with python 3.7 or more advanced version
in the anaconda platform, click on environments > base(root) > open Terminal.
in Terminal window download the function we used; numpy, pandas, PyVCF, itertools,datetime, time, os
```
$ pip install <function_name> 
```

## Database:
we downloaded fastq files of short reads from the Sequences Reads Archive (SRA) from the NCBI server,
a publicly available repository of raw sequencing data (https://www.ncbi.nlm.nih.gov/sra). 

## Script: (create_gvcf.txt)
open multisample_exome.py
first we import all required functions:

```
import numpy as np
import pandas as pd
import vcf
from itertools import islice
import datetime
import time
import os
```

define get_data_from_record function. the function gets one record (line) form the multi-sample variants file.
and return a dataframe

```
#building dataframe from the vcf file records extracting only desirable parameters
def get_data_from_record(record):
    data = [record.CHROM,
                record.POS,
                record.REF,
                record.ALT[0].sequence,
                record.QUAL,
                len(record.get_hom_refs()),
                len(record.get_hom_alts()),
                len(record.get_hets()),
                record.INFO['AF'][0],
                record.INFO['AN'],
                (record.INFO['AN'])/2,
                record.INFO['CSQ'][0].split('|')[1],
                record.INFO['CSQ'],
                record.INFO['InbreedingCoeff']]

    return data
```

we defined the columns  

```
columns = ['Chromosome',
           'Position',
           'Reference',
           'Alternate',
           'Quality',
           'homozygous reference',
           'heterozygous',
           'homozygous alternate',
           'Allele frequency',
           'Allele number',
           'samples count',
           'Consequences',
           'mutation',
           'Inbreeding']

```

we read multi-sample variants file 
and print the samples count (269 bulls)

```
df = pd.DataFrame(columns=columns)
reader = vcf.Reader(filename=r'C:\variant_effect_output.vcf.gz', compressed=False)
```

we than update the dataframe in a loop; 'get_data_from_record' function works on 1000 records every time
until the end of the reader. the output appended to the dataframe
we print the time left for the end of the dataframe building  
when 100,000 passes, we convert df to pickle file with the simple function 'df.to_pickle'
and reset a new dataframe.
in the end, all pickle files saved in a dedicated folder on the computer.

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

open merge_pkl_exome.py
first we import all required functions:

```
import numpy as np
import glob 
import pandas as pd
files = glob.glob('*.pkl') #find all the files type pkl
```
we read the files, one after the other, and convert it to dataframe
we update the dataframe and print a massage about the progress
In the end, we create several merged pickle files, contain 10 simple pickle files

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

in spyder and run the script.
first we import all required functions:

```
import glob 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import re

font = {'family' : 'normal',
        'size'   : 18}

matplotlib.rc('font', **font)
```

we create a loop that pass on each merge pickle file and extract important info to a dataframe.
we filtered low quality reads and row with loess than 10 bulls genotype (samples count).
we merged data from the exome database to our new dataframe: Allele frequency, Allele number and Inbreeding.
we filtered variant that have this pattern of genotype: garden is homozygote and the other bulls are hetrozygot (df_garden_homo)
we filtered variants that have AF > 0.5
we added a column with name of effect/consequences of the variants ('mutation')


```
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
    df_exome = df_exome.append(df)
```    

create an histogram of number of variants and number of bulls with genotype

```   
#samples count
test = df_exome.hist('samples count', bins = range(0,120),rwidth = 0.8)
ax = test[0][0]
ax.set_xticks(range(0,120,5))
ax.set_xlabel('number of bulls with genotype')
ax.set_ylabel('number of variants(1X10^6)')
```   

create an histogram of allele ferquency

```   
#AF after filtring
af_test = df_exome.hist('Allele frequency', bins = np.linspace(0,1,100),rwidth = 0.8)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,1,21))
ax.set_xlabel('Allele frequency')
ax.set_ylabel('number of variants(1X10^6)')
```

create an histogram of inbreeding 

``` 
#inbreeding after filtring
af_test = df_exome.hist('Inbreeding', bins = np.linspace(0,0.5,200),rwidth = 0.8)
ax = af_test[0][0]
ax.set_xticks(np.linspace(0,0.5,10))
ax.set_xlabel('Inbreeding')
ax.set_ylabel('number of variants(1X10^6)')
```   
