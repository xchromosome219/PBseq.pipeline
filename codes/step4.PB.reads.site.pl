#!/usr/bin/perl
#reads per site
#reads_counts.pl ${sample}.sort.bed ${sample}_readsPsite

$infile=$ARGV[0];
$outfile=$ARGV[1];
open SEE,"$infile" || die "cannot open for:$!";
open WW, ">$outfile" || die "cannot open for:$!";
%counts=();%strand=();
while(chomp($s=<SEE>)){
	@tmp=split(/\t/,$s);
	if($tmp[-1] eq '+'){$counts{$tmp[0]}{$tmp[2]}++;$strand{$tmp[0]}{$tmp[2]}="+";}
	else{$counts{$tmp[0]}{$tmp[1]+1}++;$strand{$tmp[0]}{$tmp[1]+1}="-";}
}
foreach $ch(sort {$a cmp $b} keys %counts){
	foreach $site(sort {$a<=>$b} keys %{$counts{$ch}}){
		print WW "$ch\t$site\t$counts{$ch}{$site}\t$strand{$ch}{$site}\n";
	}
}
close SEE;close WW;
