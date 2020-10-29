# WGS analysis
 
2 folders:
1. exome database
2. family 4

files included in folder 'exome database'
1. create gvcf.txt
2. ex_comb.txt
3. ex_genotype.txt
4. vep_exome.txt
5. multisample_exome.py
6. merge_pkl_exome.py
7. exome_analysis.py

files included in folder 'family 4'
1. create gvcf_family4.txt
2. family4_comb.txt
3. familt4_genotype.txt
4. vep_family4.txt
5. multisample_family4.py
6. merge_pkl_family4.py
7. family4_analysis.py
8. GTEx.V8.Testis-specificity.r-scores.csv

## Description:

In this repository we arranged the procedure of calling variants from whole genome sequencing.
we use it in two different works. The first work includes establishment of database for bovine genetics 
variations derived from whole exome sequencing. The second work includes analysis of 5 bull's genomes
from the same familial cluster.

In order to find genetic variations that cause a reduction in fertility, 
we performed a whole genome sequence analysis for infertile bulls by
using Next-generation sequencing (NGS) technology 
https://www.illumina.com/science/technology/next-generation-sequencing.html

NGS involves massive parallel sequencing of millions of DNA fragments can produce up to
hundreds of millions of short reads. 
Then bioinformatics analyses are used to piece together these fragments by
mapping the individual reads to the reference genome.

we provide a text file (create gvcf.txt) with a collection of command-line tools for analyzing 
high-throughput sequencing data with a primary focus on variant discovery. the output is 
an annotated variants file.

## Environment:
#### 1. Using Linux and HPC:
Run on Linux. we used high performance computational (HPC) platform.
Download the following packages: 
<br>Picard-v2.20.2, Samtools-v1.9, GATK4-v4.1.3.0, EnsemblVep-v97.3 and ParallelFastqDump-v0.6.5.

## Database:
Download Bos-taurus reference genome: 
ftp://ftp.ensembl.org/pub/release-101/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz

Upload reference genome, fastq files (sent from the sequencing company) and text files (commands)
to the UNIX interface (for example mobaXterm) and run the send a job in high 
performance computational (HPC) platform

## Script:
Open create_gvcf.txt and send it as a job
first we have to write the tools' paths

```
. /data/bin/miniconda2/envs/picard-v2.20.2/env_picard.sh;
. /data/bin/miniconda2/envs/samtools-v1.9/env_samtools.sh;
. /data/bin/miniconda2/envs/gatk4-v4.1.3.0/env_gatk4.sh;
. /data/bin/miniconda2/envs/ensemblVep-v97.3/env_ensembl-vep.sh;
. /data/bin/miniconda2/envs/parallelFastqDump-v0.6.5/env_parallel-fastq-dump.sh
```

we downloaded fastq files short reads from the NCBI server (https://www.ncbi.nlm.nih.gov/sra)
<br> with 'parallel-fastq-dump' (https://github.com/rvalieris/parallel-fastq-dump) and unzip them

```
parallel-fastq-dump -s <sample_name> -t 8 --tmpdir . --split-files --gzip &&
gunzip <fastq_file_name_foward> &&
gunzip <fastq_file_name_reverse> &&
```

We begin by mapping the sequence reads to the reference genome to produce SAM format
<br>with Burrows-Wheeler Aligner (bwa) tool, using default parameters

```
bwa mem {reference path} <fastq_file_name_foward> <fastq_file_name_reverse> > output.sam &&
```

Next, we create an unmapped BAM file with Picard tool

```
picard RevertSam I=output.sam O=output_u.bam ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA &&
```

first processing step is performed per-read group

```
picard AddOrReplaceReadGroups I=output_u.bam O=output_rg.bam RGID=output RGSM=output 
RGLB=wgsim RGPU=shlee RGPL=illumina &&
```

we produced a third BAM file (output_m.bam) that has alignment data
<br> and all the remaining data from the unmapped BAM (output_rg.bam)

```
picard MergeBamAlignment ALIGNED=output.sam UNMAPPED=output_rg.bam O=output_m.bam R={reference path}
SORT_ORDER=unsorted CLIP_ADAPTERS=false ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 
PRIMARY_ALIGNMENT_STRATEGY=MostDistant UNMAP_CONTAMINANT_READS=false ATTRIBUTES_TO_RETAIN=XS 
ATTRIBUTES_TO_RETAIN=XA &&
```

we mark duplicates to mitigate biases introduced by data generation steps such as PCR amplification.
<br>MarkDuplicates followed by SortSam

```
picard MarkDuplicates INPUT=output_m.bam OUTPUT=output_md.bam METRICS_FILE=output_md.bam.txt 
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname &&
set -o pipefail;
picard SortSam INPUT=output_md.bam OUTPUT=output_sorted.bam SORT_ORDER=coordinate &&
```

we used SetNmMdAndUqTags, that takes in a coordinate-sorted BAM and calculate the NM, MD and UQ 
<br>by competing with the reference
<br>NM = Edit distance to the reference
<br>MD = String encoding mismatched and deleted reference bases
<br>UQ = Phred likelihood of the segment, conditional on the mapping being correct

```
picard SetNmMdAndUqTags R={reference path} INPUT=output_sorted.bam 
OUTPUT=output_snaut.bam CREATE_INDEX=true &&
```

Call SNPs and indels variants via local re-assembly of haplotypes was generate by HaplotypeCaller.
<br>what special in HaplotypeCaller is that, whenever the program encounters a region showing signs of variation,
<br>it discards the existing mapping information and completely reassembles the reads in that region.
<br>This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally
<br>difficult to call (output_hc.bam).
<br>for more information on HaplotypeCaller: 
<br>https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
<br>
```
gatk HaplotypeCaller -R {reference path} -I output_snaut.bam -O output.g.vcf 
-ERC GVCF -bamout output_hc.bam &&
```

HaplotypeCaller runs per-sample to generate an intermediate GVCF (output.g.vcf).
<br>after we run all samples, we generat GVCF by using CombineGVCFs.
<br>this done in parallel, using 64 cores and 280 memory in the HPC options

```
gatk --java-options "-server -d64 -Xms280G -Xmx280G -XX:NewSize=250G -XX:+UseConcMarkSweepGC 
-XX:ParallelGCThreads=16 -XX:+UseTLAB" CombineGVCFs -R {reference path}
-V output_sample_1.g.vcf -V output_sample_2.g.vcf -V output_sample_3.g.vcf -V .... -O combine_all.g.vcf;
```

combine_all file used for joint genotyping of multiple samples

```
gatk --java-options "-server -d64 -Xms280G -Xmx280G -XX:NewSize=250G -XX:+UseConcMarkSweepGC 
-XX:ParallelGCThreads=16 -XX:+UseTLAB" GenotypeGVCFs -R {reference path} -V combine_all.g.vcf 
-O multisample.vcf;
```

Finally, we used Variant Effect Predictor (VEP) for annotation.
<br>after we installed Ensemble API and cache files of --species bos_taurus,
<br>we used a cache, which is the fastest and most efficient way to use VEP.
<br>we also used offline mode to eliminate all network connections for speed.

```
vep -i multisample.vcf --format vcf -o variant_effect_output.vcf --vcf --fields Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,
Protein_position,Amino_acids,Codons,Existing_variation,SYMBOL_SOURCE,GENE_PHENO,SIFT,AF --species bos_taurus --cache
--fasta {reference path in the VEP folder} --offline --dir_cache {path Cache} --dir_plugins {path Plugins} 
--synonyms {path chr_synonyms.txt} -e -fork 14 --force_overwrit 
```

Now we have annotated variants file that can be filtered!



