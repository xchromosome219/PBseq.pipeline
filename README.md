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

## create reference genome index
We established the haploid C. albicans reference genome by de novo assembly of PacBio and Illumina sequencing data, which is available in Our Nature Protocols Paper - Supplementary Data 1 (A892_Genome.zip). Users can prepare the genome index file in .fasta format using BWA by running the follow command:
```
conda activate PBISeq
bwa index A892.assembly.fasta > A892.assembly.fasta
conda deactivate
```

There are 7 sub-commands designed for specific functions.

sub-command|function
---|---
**matefq**|parse *BAM* file that mapped to GRID-seq linker to RNA-DNA mate in interleaved *FASTQ* file.
**evaluate**|evaluate the RNA-DNA mates quality and quanitity from the *BAM* file mapped to the genome.
**RNA**|identify chromatin-enriched RNAs.
**DNA**|identify RNA-enriched chromatin regions in background (trans) and foreground (cis).
**matrix**|cacluate the RNA-chromatin interation matrix.
**model**|model the network of enhancer-promoter proximity.
**stats**|statistics of GRID-seq data.

## matefq
```
GridTools.py matefq [-h] -o HDF5 [-l MINLEN] [-n READNAME] bam
```
parse *BAM* file that mapped to GRID-seq linker to RNA-DNA mate in interleaved *FASTQ* file.

arguments|option|description
---|---|---
`bam`|required|*BAM* file mapped to the GRID-seq Linker.
`-o/--hdf5 HDF5`|required|output mapping information to *HDF5*.
`-l/--minlen MINLEN`|optional|filter RNA-DNA mates with both RNA and DNA length >= MINLEN [default: 19].
`-n/--readname READNAME`|optional|rename the prefix of each read [default: no change].

## evaluate
```
GridTools.py matefq [-h] -o HDF5 [-k BINK] [-m WINM] -g GTF bam
```
evaluate the RNA-DNA mates quality and quanitity from the *BAM* file mapped to the genome.

arguments|option|description
---|---|---
`bam`|required|*BAM* file mapped to the GRID-seq Linker.
`-o/--hdf5 HDF5`|required|output mapping information to *HDF5*.
`-g/--gtf GTF`|required|gene annotation in *GTF* format.
`-k/--bink BINK`|optional|bin size (kb) of the genome [default: 10 kb].
`-m/--winm WINM`|optional|moving window for smoothing in bins [default: 10].

## RNA
```
GridTools.py RNA [-h] [-e EXPRS] [-s SCOPE] hdf5
```
identify chromatin-enriched RNAs and evaluate the gene expressoin levels as well as interaction scopes.

arguments|option|description
---|---|---
`hdf5`|required|*HDF5* file with mapping information evaluated by *GridTools*.
`-e/--exprs EXPRS`|optional|output file for the gene expression [default: print to the screen].
`-s/--scope SCOPE`|optional|output file for the RNA interaction scope [default: None].

## DNA
```
GridTools.py DNA [-h] hdf5
```
identify RNA-enriched chromatin regions in background (trans) and foreground (cis).

arguments|option|description
---|---|---
`hdf5`|required|*HDF5* file with mapping information evaluated by *GridTools*.


## stats
```
GridTools.py stats [-h] -p PREFIX [-b] [-c] [-l] [-r] hdf5
```
statistics of GRID-seq data.

arguments|option|description
---|---|---
`hdf5`|required|*HDF5* file with mapping information evaluated by *GridTools*.
`-p/--prefix PREFIX`|required|prefix of output file names.
`-b/--bases`|optional|if output the summary of base-position information for RNA, Linker and DNA [default: No].
`-c/--counts`|optional|if output the summary of mapping information in read counts [default: No].
`-l/--lengths`|optional|if output the distribution of sequence length for RNA, Linker and DNA [default: No].
`-r/--resolution`|optional|if output the resolution information of the library [default: No].


## matrix
```
GridTools.py matrix [-h] [-k RPK] [-x DRPK] hdf5
```
cacluate the RNA-chromatin interation matrix.

arguments|option|description
---|---|---
`hdf5`|required|*HDF5* file with mapping information evaluated by *GridTools*.
`-k/--rpk RPK`|optional|cutoff of RNA reads per kilobase in the gene body. [default: 100]
`-x/--drpk DRPK`|optional|cutoff of DNA reads per kilobase at the maximum bin [default: 10].


## model
```
GridTools.py model [-h] [-k RPK] [-x DRPK] -e ELEBED [-z ZSCORE] hdf5
```
model the network of enhancer-promoter proximity.

arguments|option|description
---|---|---
`hdf5`|required|*HDF5* file with mapping information evaluated by *GridTools*.
`-e/--elebed ELEBED`|required|*BED* file of regulatory elements (eg. enhancers and promoters).
`-k/--rpk RPK`|optional|cutoff of RNA reads per kilobase in the gene body. [default: 100]
`-x/--drpk DRPK`|optional|cutoff of DNA reads per kilobase at the maximum bin [default: 10].
`-z/--zscore ZSCORE`|optional|z-score to filter the significant proximity [default:-10].


# Simple Code Example
The test data and demo code for the GRID-seq analysis with test dataset at: `http://fugenome.ucsd.edu/gridseq/datasets/gridseq.test10M.raw.fq.gz`:

Demo code is in the gridtools/pipeline/Snakefile.py


# License
The *GridTools* is licensed under **MIT**. The detailed license terms are in the **LICENSE** file.
