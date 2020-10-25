# verients filtretion of 5 bulls (family 4)

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
The variants were annotated by using Variant Effect Predictor (VEP) and assessed for their relevance to the phenotype (see result for future description). Then, filtering was done using a dedicated python script that uses the VCF parser (PyVCF). Variants were removed if they had quality score less than <30.

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

downloaded ensembl-vep package, follow installation instructions:
https://m.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer

## Script: (create_gvcf.txt)
the script for WGS analysis will be in text file:
first, the path for all tools will be deateild:

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
parallel-fastq-dump -s <sample_name> -t 8 --tmpdir . --split-files --gzip &&
gunzip <fastq_file_name_foward> &&
gunzip <fastq_file_name_reverse> &&
```

We begin by mapping the sequence reads to the reference genome to produce a file in SAM format
with Burrows-Wheeler Aligner (bwa) tool, using default parameters

```
bwa mem ./path to reference file <fastq_file_name_foward> <fastq_file_name_reverse> > output.sam &&
```

Next, we create an unmapped BAM file with Picard tool

```
picard RevertSam I=output.sam O=output_u.bam ATTRIBUTE_TO_CLEAR=XS ATTRIBUTE_TO_CLEAR=XA &&
```

first processing step is performed per-read group

```
picard AddOrReplaceReadGroups I=output_u.bam O=output_rg.bam RGID=output RGSM=output RGLB=wgsim RGPU=shlee RGPL=illumina &&
```

we produced a third BAM file (output_m.bam) that has alignment data (output.sam)
and all the remaining data from the unmapped BAM (output_rg.bam)

--SORT_ORDER: The order in which the merged reads should be output

--CLIP_ADAPTERS: Whether to clip adapters where identified

--ADD_MATE_CIGAR: add the MC tag to each read

--MAX_INSERTIONS_OR_DELETIONS: The maximum number of insertions or deletions permitted for an alignment to be included.
Set to -1 to allow any number of insertions or deletions.

--PRIMARY_ALIGNMENT_STRATEGY: Strategy for selecting primary alignment when the aligner has provided more than 
one alignment for a pair or fragment.

--UNMAP_CONTAMINANT_READS: Whether to identify extremely short alignments as cross-species contamination
and unmap the reads. 

--ATTRIBUTES_TO_RETAIN : Attributes from the alignment record that should be removed when merging

```
picard MergeBamAlignment ALIGNED=output.sam UNMAPPED=output_rg.bam O=output_m.bam R= ./path to reference file
SORT_ORDER=unsorted CLIP_ADAPTERS=false ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 
PRIMARY_ALIGNMENT_STRATEGY=MostDistant UNMAP_CONTAMINANT_READS=false ATTRIBUTES_TO_RETAIN=XS ATTRIBUTES_TO_RETAIN=XA &&
```

we mark duplicates to mitigate biases introduced by data generation steps such as PCR amplificatio.
MarkDuplicates followed by SortSam
--OPTICAL_DUPLICATE_PIXEL_DISTANCE: The maximum offset between two duplicate clusters in order to
consider them optical duplicates. The default is appropriate for unpatterned versions of the 
Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. 

```
picard MarkDuplicates INPUT=output_m.bam OUTPUT=output_md.bam METRICS_FILE=output_md.bam.txt 
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname &&
set -o pipefail;
picard SortSam INPUT=output_md.bam OUTPUT=output_sorted.bam SORT_ORDER=coordinate &&
```

we used SetNmMdAndUqTags, that takes in a coordinate-sorted BAM and calculate the NM, MD and UQ by compating with the reference
NM = Edit distance to the reference
MD = String encoding mismatched and deleted reference bases
UQ = Phred likelihood of the segment, conditional on the mapping being correct

```
picard SetNmMdAndUqTags R=<path to reference file> INPUT=output_sorted.bam OUTPUT=output_snaut.bam CREATE_INDEX=true &&
```

Call SNPs and indels variants via local re-assembly of haplotypes was generate by HaplotypeCaller.
what spiceal in HaplotypeCaller is that, whenever the program encounters a region showing signs of variation,
it discards the existing mapping information and completely reassembles the reads in that region.
This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally
difficult to call (output_hc.bam).
for more deteild of how HaplotypeCaller work do to https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller

```
gatk HaplotypeCaller -R ./path to reference file -I output_snaut.bam -O output.g.vcf -ERC GVCF -bamout output_hc.bam &&
```

HaplotypeCaller runs per-sample to generate an intermediate GVCF (output.g.vcf).
after we run all samples, we generat GVCF by using CombineGVCFs.
this done in parralel, using 64 cores and 280 memory in the HPC options

```
gatk --java-options "-server -d64 -Xms280G -Xmx280G -XX:NewSize=250G -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=16 -XX:+UseTLAB" CombineGVCFs -R ./path to reference file -V output_sample_1.g.vcf -V output_sample_2.g.vcf -V output_sample_3.g.vcf -V .... -O combine_all.g.vcf;
```

this combine_all file than be used in GenotypeGVCFs for joint genotyping of multiple samples

```
gatk --java-options "-server -d64 -Xms280G -Xmx280G -XX:NewSize=250G -XX:+UseConcMarkSweepGC -XX:ParallelGCThreads=16 
-XX:+UseTLAB" GenotypeGVCFs -R ./path to reference file -V combine_all.g.vcf -O multisample.vcf;
```

Finally, we used Variant Effect Predictor (VEP) for annotation.
after we installed vEnsemble API and cache files of --species bos_taurus,
we used a cache (--cache), which is the fastest and most efficient way to use VEP.
we also used offline mode to eliminate all network connections for speed.

```
vep -i multisample.vcf --format vcf -o variant_effect_output.vcf --vcf --fields Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,
Protein_position,Amino_acids,Codons,Existing_variation,SYMBOL_SOURCE,GENE_PHENO,SIFT,AF --species bos_taurus --cache
--fasta ./path to reference file in VEP --offline --dir_cache ./path Cache --dir_plugins ./path Plugins 
--synonyms ./path chr_synonyms.txt -e -fork 14 --force_overwrit 
```
