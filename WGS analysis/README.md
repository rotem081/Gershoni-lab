# WGS analysis
 
2 folders:
1. exome database
2. family 4

files included in folder 'exome database'
1. creat gvcf.txt
2. ex_comb.txt
3. ex_genotype.txt
4. vep_exome.txt
5. multisample_exome.py
6. merge_pkl_exome.py
7. exome_analysis.py

files included in folder 'family 4'
1. creat gvcf_family4.txt
2. family4_comb.txt
3. familt4_genotype.txt
4. vep_family4.txt
5. multisample_family4.py
6. merge_pkl_family4.py
7. family4_analysis.py
8. GTEx.V8.Testis-specificity.r-scores.csv

## Description:

In this repository we arranged the procedure of calling varients from whole genome sequencing.
we use it in two different works. The first work includes establishment of database for bovine genetis 
variations derived from whole exome sequencgin. The second work includes analysis of 5 bull's genomes
from the same familial cluster. 

raw data (fastq files) was undergone bioinformatics analysis:
1. In the first work, esteblishment of the exome databse, the obtained fastq files short reads from
the Sequences Reads Archive (SRA) from the NCBI server (https://www.ncbi.nlm.nih.gov/sra).
2. In the second work, the genome analysis of one family, the obtined fastq files sent from the sequencing company.


## Environment:
#### 1. Using mobaXterm AND HPC:
we used high preformence computational (HPC) platform to run this WGS analysis.
the tools and version we used in linux enveiroment: picard-v2.20.2, samtools-v1.9, 
gatk4-v4.1.3.0, ensemblVep-v97.3 and parallelFastqDump-v0.6.5.

download mobaXterm server from https://mobaxterm.mobatek.net/download-home-edition.html

## Database:
download Bos-taurus reference genome: 
ftp://ftp.ensembl.org/pub/release101/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz).

upload reference genom, fastq files (sent from the sequencing company) and text files (commands)
to the mobaXterm interface and run the commands in high preformence computational (HPC) platform

## Script:
the script for WGS analysis will be in text file:
first, the path for all tools: 
```
. /data/bin/miniconda2/envs/picard-v2.20.2/env_picard.sh;
. /data/bin/miniconda2/envs/samtools-v1.9/env_samtools.sh;
. /data/bin/miniconda2/envs/gatk4-v4.1.3.0/env_gatk4.sh;
. /data/bin/miniconda2/envs/ensemblVep-v97.3/env_ensembl-vep.sh;
. /data/bin/miniconda2/envs/parallelFastqDump-v0.6.5/env_parallel-fastq-dump.sh
```

we downloaded fastq files short reads from the NCBI server (https://www.ncbi.nlm.nih.gov/sra)
with 'parallel-fastq-dump' (https://github.com/rvalieris/parallel-fastq-dump) and unzip tham
```
parallel-fastq-dump -s ERR1600413 -t 8 --tmpdir . --split-files --gzip &&
gunzip ERR1600413_1.fastq &&
gunzip ERR1600413_2.fastq &&
```

we mapped reads to the reference genome with Burrows-Wheeler Aligner (bwa) tool, using default parameters
```
bwa mem <path to reference file> <fastq_file_name_foward> <fastq_file_name_foward> > output.sam &&
```

Next, BAM files were locally realigned, read groups were added to the unmapped BAM files,
PCR duplicates were removed, and clean BAM files were coordinated, sorted, and indexed by Picard tool
```
picard RevertSam I=output.sam O=output_u.bam ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA &&
picard AddOrReplaceReadGroups I=output_u.bam O=output_rg.bam RGID=output RGSM=output RGLB=wgsim RGPU=shlee RGPL=illumina &&
picard MergeBamAlignment ALIGNED=output.sam UNMAPPED=output_rg.bam O=output_m.bam R=<path to reference file> 
SORT_ORDER=unsorted CLIP_ADAPTERS=false ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant UNMAP_CONTAMINANT_READS=false ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA &&
picard MarkDuplicates INPUT=output_m.bam OUTPUT=output_md.bam METRICS_FILE=output_md.bam.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname &&
set -o pipefail;
picard SortSam INPUT=output_md.bam OUTPUT=output_sorted.bam SORT_ORDER=coordinate &&
picard SetNmMdAndUqTags R=<path to reference file> INPUT=output_sorted.bam OUTPUT=output_snaut.bam CREATE_INDEX=true &&
```

After BAM files were created for every sample, variant calling was done with the Genome Analysis Tool kit (GATK, version 4.1.6.0)
as recommended by the GATK workflow (link to our GitHub: WGS analysis). Call SNPs (single nucleotide polymorphism) 
and indels (insertion and deletion) variants via local re-assembly of haplotypes was generate by HaplotypeCaller, 
and then an intermediate GVCF (genome variants call file) was generate per sample. 
GVCF can be used for joint genotyping of multiple samples in a very efficient way. 
Genotyping files were produced after all samples were combined with CombineGVCFs, by using GenotypeGVCFs (Figure 3).
The variants were annotated by using Variant Effect Predictor (VEP) and assessed for their relevance to the phenotype (see result for future description). Then, filtering was done using a dedicated python script that uses the VCF parser (PyVCF). 
Variants were removed if they had quality score less than <30.



gatk HaplotypeCaller -R /home/ARO.local/rotemv/Projects/SRA/REF_files/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa -I ERR1600413_snaut.bam -O ERR1600413.g.vcf -ERC GVCF -bamout ERR1600413_hc.bam &&




## Special Comments:
