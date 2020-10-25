# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:44:37 2019

@author: rotemv
"""
import numpy as np
import pandas as pd
import vcf
from itertools import islice
import datetime
import time
import os

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


df = pd.DataFrame(columns=columns)
reader = vcf.Reader(filename=r'C:\variant_effect_output.vcf.gz', compressed=False)
#print(reader.metadata) 
 
print(len(reader.samples)) #samples count
print(reader.samples[:10])



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
            df.to_pickle(f'result{file_counter:03}_exome.pkl') #save file           
            print(f'saved progress in file number {file_counter:03}_exome')
            df = pd.DataFrame(columns=columns) #update new dataframe
            file_counter +=1
except:
    df.to_pickle(f'result{file_counter:03}_exome.pkl')
    raise #raise the error


