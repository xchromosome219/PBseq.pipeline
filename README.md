# PBIseq.pipeline

**PBIseq** is a set of tools designed for analyzing Genome-wide piggyBac transposon-based mutagenesis and quantitative insertion site analysis in haploid Candida species. 

# Install and Requirements
**PBIseq** is developed in *Perl 5* v22, with dependent libraries as below:

## Requirements

-	A computer with 16 or more 64-bit processors and 64 GB or more RAM is recommended.

-	A Ubuntu Linux distro, such as Ubuntu server 18.04 LTS, is recommended (https://www.ubuntu.com/download/server).

- A package and environment management system, Conda (v4.8.0; https://docs.conda.io/en/latest/miniconda.html).

-	Python 3 (v3.7.4; https://www.python.org).

-	Perl 5 (v22; https://www.perl.org/get.html).

-	BWA software package (0.7.17-r1188; http://www.bio-bwa.sourceforge.net). 

-	SAMtools (v.1.9; http://www.samtools.sourceforge.net).

-	Bedtools (v2.29.2; https://bedtools.readthedocs.io/en/latest).

-	FASTX Toolkit (v0.0.14; http://hannonlab.cshl.edu/fastx_toolkit).

-	PBISeq scripts (https://github.com/xchromosome219/PBseq.pipeline).


## Installation
We recommend using *conda* package manager (https://conda.io/miniconda.html) to install required bioinformatics tools and packages in *Bioconda* (https://bioconda.github.io/). 

create Ca_PBISeq analysis work enviroment
```
conda create –n PBISeq python=3 
conda info –envs
```

And install `bwa`, `samtools`, `bedtools`, `FASTX Toolkit`, `perl-config-general`:
```
conda install -c bioconda bwa samtools bedtools fastx_toolkit perl-config-general
```

# Usage

## Create reference genome index
We established the haploid C. albicans reference genome by de novo assembly of PacBio and Illumina sequencing data, which is available in Our Nature Protocols Paper - Supplementary Data 1 (A892_Genome.zip). Users can prepare the genome index file in .fasta format using BWA by running the follow command:
```
conda activate PBISeq
bwa index A892.assembly.fasta > A892.assembly.fasta
conda deactivate
```

## Create a new home directory < Ca_PBISeq_results >, and copy demo sequence raw data to Ca_PBISeq result directory
Demo data with the name < PBISeq.demo.R1.fq.gz > and < PBISeq.demo.R2.fq.gz > can be downloaded from https://www.ncbi.nlm.nih.gov/bioproject/PRJNA564479
```
mkdir –p ./Ca_PBISeq_results
# copy demo sequence raw data to Ca_PBISeq result directory
zcat PBISeq.demo.R1.fq.gz > ./Ca_PBISeq_results/PBISeq.demo.R1.fastq
zcat PBISeq.demo.R2.fq.gz > ./Ca_PBISeq_results/PBISeq.demo.R2.fastq
```

## Step 1 remove Illumina adaptor
Trim Illumina adaptor sequences in FASTQ files, using our in-house Perl script (step1.remove.adaptor.pl) by running the following command:
```
cat ./Ca_PBISeq_results/PBISeq.demo.R1.fastq  |\
perl step1.remove.adaptor.pl \
-p ./Ca_PBISeq_results/PBISeq.demo.R1 \
-q ./Ca_PBISeq_results/PBISeq.demo \
-t 15
```

## Step 2 identify the PB transposon sequence in Read 1
Filter all Read 1 for the presence of the transposon-genomic junction and prepare reads for mapping by trimming the transposon sequence (TTTCTAGGG) in Read 1 to leave only the insertion motif and the genomic sequence. Then, a text file < specific_reads_ratio.txt > will be generated, which contains various statistics for the sequencing data, including the percentage of reads being mapped (that have the transposon), the percentage of reads kept in the initial ARG4 site, the percentage of reads due to non-specific PCR amplification, the number of raw reads and the total number of reads being mapped. Run the following command: 
Demo data (cleaned data) with the name < PBISeq.demo.R1_filter_cut.fastq > can also be downloaded from here： https://github.com/xchromosome219/PBseq.pipeline/tree/master/Ca_PBSeq_results/example
```
perl step2.PB.site.pl ./Ca_PBISeq_results/PBISeq.demo.R1_clean.fastq \
./Ca_PBISeq_results/PBISeq.demo.R1_filter_cut.fastq

conda activate PBISeq
fastx_reverse_complement -Q 33 \
–I ./Ca_PBISeq_results/PBISeq.demo.R1_filter_cut.fastq \
-o ./Ca_PBISeq_results/PBISeq.demo.R1_filter_rc.fastq
conda deactivate

# calculate PB transposon insertion site reads ratio
perl step3.PB.specific_ratio.pl \
./Ca_PBISeq_results/PBISeq.demo.R1_clean.fastq ./Ca_PBISeq_results/specific_reads_ratio.txt
```

## Step 3 align the reads against the genome
Align the reads against the genome to filter unambiguously mapped reads. Run the following scripts in series or run in parallel. Two output files will be produced:
```
# align the reads to the reference genome using bwa mem module
conda activate PBISeq
bwa mem A892.assembly.fasta \
./Ca_PBISeq_results/PBISeq.demo.R1_filter_rc.fastq > ./Ca_PBISeq_results/PBISeq.demo.R1.sam

# convert sam file to bam file
samtools view -bS \
./Ca_PBISeq_results/PBISeq.demo.R1.sam > ./Ca_PBISeq_results/PBISeq.demo.R1.bam

# sort bam file
samtools sort ./Ca_PBISeq_results/PBISeq.demo.R1.bam \
-o ./Ca_PBISeq_results/PBISeq.demo.R1.sorted.bam

# filter unambiguously mapped reads by the bamToBed function of SAMtools to transform the bam file to the BED format
samtools index PBISeq.demo.R1.sorted.bam
bamToBed -i PBISeq.demo.R1.sorted.bam > PBISeq.demo.R1.sorted.bed
conda deactivate
```

## Step 4 count the read number per insertion site
< PBISeq.demo_readsPsite >:  contains the processed mapping output with the format: “chromosome or contig name”, “insertion position”, “the total number of reads mapped to that position”, “strand orientation”.
```
# count the read number per insertion site 
Perl step4.PB.reads.site.pl \
./Ca_PBISeq_results/PBISeq.demo.R1.sorted.bed ./Ca_PBISeq_results/PBISeq.demo_readsPsite
```

## Step 5 convert read count per site to read count per gene
< name_readsPgene >: lists genes mapped by insertions. The format is: “chromosome or contig name”, “gene position”, “the total number of reads mapped to the ORF of that gene”.
```
# convert read count per site to read count per gene 
perl step5.PB.reads.genelevel.pl \
./Ca_PBISeq_results/PBISeq.demo.readsPsite ./Ca_PBISeq_results/PBISeq.demo.readsPgene
```

# License
The *PBIseq* is licensed under **MIT**. The detailed license terms are in the **LICENSE** file.

