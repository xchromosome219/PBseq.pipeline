#!/bin/bash
if [ $# != 4 ]; then
    echo  "USAGE: ./PB.pipeline.sh <analysis_dir> <config_file> <number1> <number2>"
    exit
fi

analysis_dir=$1
config_file=$2
number1=$3
number2=$4

# source activate wes
#####################################################

# tools
command_path=./codes

# references
### diploid A22
#reference=/data/others_ref/C_albicans/A22/C_albicans_SC5314_A22_current_chromosomes.fasta
#INDEX=/data/others_ref/C_albicans/A22/C_albicans_SC5314_A22_current_chromosomes.fasta
#bowtie2_index=/data/others_ref/C_albicans/A22/bowtie2_index/C_albicans_SC5314_A22_current_chromosomes

### haploid GZY892
reference=/data/ref/C_albicans/Haploid/GZY892/genome/GZY892.assembly.fasta
INDEX=/data/ref/C_albicans/Haploid/GZY892/bwa/GZY892.assembly.fasta
bowtie2_index=/data/ref/C_albicans/Haploid/GZY892/bowtie2/GZY892.assembly

mkdir -p $analysis_dir/tmp
TMPDIR=$analysis_dir/tmp

cat $config_file |while read id
do
	arr=($id)
	fq1=${arr[1]}
	fq2=${arr[2]}
	sample=${arr[0]}

	if((i%$number1==$number2))
	then

	   if [  ! -f $analysis_dir/${sample}/${sample}.summary ]; then
################################################################
################### Start : Run PBseq pipeline #################
################################################################

#### step1 ####
#### remove adapter
if [ ! -f  $analysis_dir/trim/${sample}.R1_001_clean.fastq ]; then

mkdir -p $analysis_dir/${sample}
zcat $fq1 > $analysis_dir/${sample}/${sample}.R1_001.fastq
zcat $fq2 > $analysis_dir/${sample}/${sample}.R2_001.fastq

source activate PBSeq
fastx_reverse_complement -Q 33 -i $analysis_dir/${sample}/${sample}.R2_001.fastq -o $analysis_dir/${sample}/${sample}_rc.fastq
conda deactivate

cd $analysis_dir/${sample}/
cat ./${sample}.R1_001.fastq | perl $command_path/pairend_adaptor_rm_v3.pl -p ${sample}.R1_001 -q ${sample} -t 15
cd ..
fi

#### step2 ####
#### get PB inserted reads (filter, cut, rc)
if [  ! -f  $analysis_dir/${sample}/R1_filter_rc.fastq ]; then
mkdir -p $analysis_dir/${sample}

source activate PBSeq
$command_path/filter_cut.pl $analysis_dir/trim/${sample}.R1_001_clean.fastq $analysis_dir/${sample}/R1_filter_cut.fastq
fastx_reverse_complement -Q 33 -i $analysis_dir/${sample}/R1_filter_cut.fastq -o $analysis_dir/${sample}/R1_filter_rc.fastq
conda deactivate
fi

#### step3 ####
#### statistics for PB read counts
if [  ! -f  $analysis_dir/${sample}/specific_reads_ratio.txt ]; then
$command_path/specific_ratio.pl ${sample} $analysis_dir/${sample}/${sample}.R1_001_clean.fastq $analysis_dir/${sample}/${sample}_specific_reads_ratio.txt
fi

#### step4 ####
#### bwa mapping
if [  ! -f  $analysis_dir/${sample}/${sample}.sort.bed ]; then
source activate PBSeq
bwa mem $INDEX $analysis_dir/${sample}/R1_filter_rc.fastq > $analysis_dir/${sample}/${sample}.sam
samtools view -bS $analysis_dir/${sample}/${sample}.sam > $analysis_dir/${sample}/${sample}.bam
samtools sort $analysis_dir/${sample}/${sample}.bam -o $analysis_dir/${sample}/${sample}.sort.bam
samtools index $analysis_dir/${sample}/${sample}.sort.bam
bamToBed -i $analysis_dir/${sample}/${sample}.sort.bam > $analysis_dir/${sample}/${sample}.sort.bed
conda deactivate
fi

#### step5 ####
if [  ! -f  $analysis_dir/${sample}/${sample}.summary ]; then
#### delete abnormal sites, count reads per site & gene
$command_path/delete_abnormal_site.pl $analysis_dir/${sample}/${sample}.sort.bed $analysis_dir/${sample}/${sample}_delete.bed
$command_path/readsPsite.pl $analysis_dir/${sample}/${sample}_delete.bed $analysis_dir/${sample}/${sample}_readsPsite
$command_path/readsPgene.pl $analysis_dir/${sample}/${sample}_readsPsite $analysis_dir/${sample}/${sample}_readsPgene
cd $analysis_dir/${sample}
wc -l $analysis_dir/alignment/${sample}.sort.bed | cut -d ' ' -f 1 > $analysis_dir/${sample}/mapped_reads_num
wc -l $analysis_dir/${sample}/${sample}_delete.bed | cut -d ' ' -f 1 > $analysis_dir/${sample}/after_delete_reads_num
wc -l $analysis_dir/${sample}/${sample}_readsPsite | cut -d ' ' -f 1 > $analysis_dir/${sample}/covered_site
wc -l $analysis_dir/${sample}/${sample}_readsPgene | cut -d ' ' -f 1 > $analysis_dir/${sample}/covered_gene
# echo -e "$sample\t$mapped_reads_num\t$after_delete_reads_num\t$covered_site\t$covered_gene"
paste ${sample} $analysis_dir/${sample}/mapped_reads_num $analysis_dir/${sample}/after_delete_reads_num $analysis_dir/${sample}/covered_site $analysis_dir/${sample}/covered_gene > $analysis_dir/${sample}/${sample}.summary
cd ..
fi

################################################################
################### End : Run PBseq pipeline #################
################################################################
		fi
	fi
	i=$((i+1))
done
