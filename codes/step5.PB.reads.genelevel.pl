#!/usr/bin/perl
#gene counts which are covered by reads
#reads per gene
#readsPgene.pl ${sample}_readsPsite ${sample}_readsPgene

$infile=$ARGV[0];
$outfile=$ARGV[1];

#open SEE2,"/800T/wangjb/lizeyao/PB_commands/genome/gene.gff" || die "cannot open for:$!";
#open SEE2,"/Share/home/wangjb/lizeyao/PB_commands/genome/gene.gff" || die "cannot open for:$!";
open SEE2,"/data/others_ref/C_albicans/Haploid/GZY892/annotation/GZY892.gff" || die "cannot open for:$!";
open WW,">$outfile" || die "cannot open for:$!";
while(chomp($s2=<SEE2>)){
	$counts=0;
	@tmp2=split(/\t/,$s2);
        open SEE1,"$infile" || die "cannot open for:$!";
	while(chomp($s1=<SEE1>)){
		@tmp1=split(/\t/,$s1);
		if($tmp2[0]=~/$tmp1[0]/){
			if(($tmp1[-1] eq "+")&&($tmp1[1]>=$tmp2[3]+3 && $tmp1[1]<=$tmp2[4])){
				$counts+=$tmp1[2];
			}
			elsif(($tmp1[-1] eq "-")&&($tmp1[1]>=$tmp2[3] && $tmp1[1]<=$tmp2[4]-3)){
				$counts+=$tmp1[2];
			}
		}
	}
	close SEE1;
	if($counts>0){print WW "$tmp2[0]\t$tmp2[3]\t$tmp2[4]\t$counts\n";}
}
close SEE2;close WW;
