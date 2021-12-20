## Description:
In this repository we arranged the procedure of assessment of genomic variants.
we built dedicated python script, that uses the VCF parser (PyVCF), to filter and assess the variants. 

## Database:
we downloaded fastq files of short reads from the Sequences Reads Archive (SRA) from the NCBI server,
a publicly available repository of raw sequencing data (https://www.ncbi.nlm.nih.gov/sra). 


# analysis pipeline (from fastq to gvcf)
txt file that contains the commands and the path to the relevent tools (running on HPC)

. /data/bin/miniconda2/envs/picard-v2.20.2/env_picard.sh;
. /data/bin/miniconda2/envs/samtools-v1.9/env_samtools.sh;
. /data/bin/miniconda2/envs/gatk4-v4.1.3.0/env_gatk4.sh;
. /data/bin/miniconda2/envs/ensemblVep-v97.3/env_ensembl-vep.sh;
. /data/bin/miniconda2/envs/parallelFastqDump-v0.6.5/env_parallel-fastq-dump.sh

parallel-fastq-dump -s ERR1600413 -t 8 --tmpdir . --split-files --gzip &&
gunzip ERR1600413_1.fastq &&
gunzip ERR1600413_2.fastq &&
bwa mem /home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa ERR1600413_1.fastq ERR1600413_2.fastq > ERR1600413.sam &&
picard RevertSam I=ERR1600413.sam O=ERR1600413_u.bam ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA &&
picard AddOrReplaceReadGroups I=ERR1600413_u.bam O=ERR1600413_rg.bam RGID=ERR1600413 RGSM=ERR1600413 RGLB=wgsim RGPU=shlee RGPL=illumina &&
picard MergeBamAlignment ALIGNED=ERR1600413.sam UNMAPPED=ERR1600413_rg.bam O=ERR1600413_m.bam R=/home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa SORT_ORDER=unsorted CLIP_ADAPTERS=false ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant UNMAP_CONTAMINANT_READS=false ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA &&
picard MarkDuplicates INPUT=ERR1600413_m.bam OUTPUT=ERR1600413_md.bam METRICS_FILE=ERR1600413_md.bam.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname &&
set -o pipefail;
picard SortSam INPUT=ERR1600413_md.bam OUTPUT=ERR1600413_sorted.bam SORT_ORDER=coordinate &&
picard SetNmMdAndUqTags R=/home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa INPUT=ERR1600413_sorted.bam OUTPUT=ERR1600413_snaut.bam CREATE_INDEX=true &&
gatk HaplotypeCaller -R /home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa -I ERR1600413_snaut.bam -O ERR1600413.g.vcf -ERC GVCF -bamout ERR1600413_hc.bam &&
rm ERR1600413_1.fastq;
rm ERR1600413_2.fastq;
rm ERR1600413.sam;
rm ERR1600413_u.bam;
rm ERR1600413_rg.bam;
rm ERR1600413_m.bam;
rm ERR1600413_md.bam;
rm ERR1600413_sorted.bam;
rm ERR1600413_snaut.bam;
rm ERR1600413_snaut.bai;
rm ERR1600413_hc.bam;


# combine all gvcf
(first combine to groups and than combine all groups to one large gvcf)
gatk --java-options "-server -d64 -Xms280G -Xmx280G -XX:NewSize=250G -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=16 -XX:+UseTLAB" CombineGVCFs -R /home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa -V combine_file_1.g.vcf -V combine_file_2.g.vcf -V combine_file_3.g.vcf -V combine_file_4_1.g.vcf -V combine_file_4_2.g.vcf -V combine_file_5.g.vcf -V combine_file_6.g.vcf -V combine_file_7.g.vcf -V combine_file_8.g.vcf -V combine_file_9.g.vcf -V combine_file_10.g.vcf -V combine_file_11.g.vcf -V combine_file_12.g.vcf -V combine_file_13.g.vcf -V combine_file_14.g.vcf -V combine_file_15.g.vcf -V combine_file_16.g.vcf -V combine_file_17.g.vcf -V combine_file_18.g.vcf -V combine_file_19.g.vcf -V combine_file_20.g.vcf -O combine_all_exomes.g.vcf;


# genotyping (from multisample gvcf to vcf)
gatk --java-options "-server -d64 -Xms280G -Xmx280G -XX:NewSize=250G -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=16 -XX:+UseTLAB" GenotypeGVCFs -R /home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa -V combine_all_exomes.g.vcf -O multisample.vcf;


# vep stage
vep -i multisample.vcf.gz --format vcf -o variant_effect_output.vcf --vcf --fields Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,SIFT --species bos_taurus --cache --fasta /data/bin/ensemblVep-v98/ensembl-vep/Cache/bos_taurus/98_ARS-UCD1.2/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz --offline --dir_cache /data/bin/ensemblVep-v98/ensembl-vep/Cache --dir_plugins /data/bin/ensemblVep-v98/ensembl-vep/Plugins --synonyms /data/bin/ensemblVep-v98/ensembl-vep/Cache/bos_taurus_refseq/98_ARS-UCD1.2/chr_synonyms.txt -e -fork 14 --force_overwrit 


# analysis vcf with Pyvcf python tools
#### 1. Using mobaXterm and anaconda:
make sure that your virtual environment is installed with python 3.7 or more advanced version
in the anaconda platform, click on environments > base(root) > open Terminal.
in Terminal window download the function we used; numpy, pandas, PyVCF, itertools,datetime, time, os
```
$ pip install <function_name> 
```

open python script multisample_exome.py
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
