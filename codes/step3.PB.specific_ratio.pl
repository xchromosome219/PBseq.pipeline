#!/usr/bin/perl
#the ratio of reads--specific insert, resided at the ARG4, non-specific
#specific_ratio.pl $sample R1_clean.fastq specific_reads_ratio.txt

$sample_name=$ARGV[0];
$infile=$ARGV[1];
$outfile=$ARGV[2];
open SEE,"$infile" || die "cannot open for:$!";
open WW,">>$outfile" || die "cannot open for:$!";

$hang=1;
$insert_num=0;
$origin=0;
while(chomp($s=<SEE>)){
	$tmp[$hang%4]=$s;
	if($hang%4==0){
            if(($tmp[2]=~/TTTCTAGGG..../) && ($tmp[2]!~/TTTCTAGGG....AGAATT/)){$insert_num++;}
            elsif(($tmp[2]=~/TTTCTAGGG..../) && ($tmp[2]=~/TTTCTAGGG....AGAATT/)){$origin++;}
	}
	$hang++;
}
$reads_num=($hang-1)/4;
$insert_r=($insert_num+0.00)/$reads_num;
$origin_r=($origin+0.00)/$reads_num;
$random_r=1-$insert_r-$origin_r;

print WW "$sample_name\t$insert_r\t$origin_r\t$random_r\t$reads_num\t$insert_num\n";
close SEE; close WW;
