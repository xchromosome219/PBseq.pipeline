#!/usr/bin/perl
# site_readsNum_gene.pl readsPsite site_readsNum_gene.txt

open SEE,"/data/others_ref/C_albicans/A22/C_albicans_SC5314_A22_current_features.tab" || die "cannot open for:$!";
while(chomp($s=<SEE>)){
        if($s=~/^!/){next;}
        else{
                @tmp=split(/\t/,$s);
                if(($tmp[0]!~/_B$/) && ($tmp[1] ne "")){$id_name{$tmp[0]}=$tmp[1];}
        }
}
close SEE;

$infile=$ARGV[0];
$outfile=$ARGV[1];
open SEE1,"$infile" || die "cannot open for:$!";
open WW,">$outfile" || die "cannot open for:$!";
while(chomp($s1=<SEE1>)){
	@tmp1=split(/\t/,$s1);
	print WW "$s1\t";
	open SEE2,"/data/others_ref/C_albicans/A22/C_albicans_SC5314_A22_current_features.gff" || die "cannot open for:$!";
	while(chomp($s2=<SEE2>)){
		@tmp2=split(/\t/,$s2);
		if(($tmp2[0]=~/$tmp1[0]/) && $tmp1[1]>=$tmp2[3] && $tmp1[1]<=$tmp2[4]){
			if($tmp2[-1]=~/^(ID=)(.*_A)(;Name=)(.*)$/){$gene=$2;}
			elsif(($tmp2[0]=~/chrM/) && ($tmp2[-1]=~/^(ID=)(.*?)(;Name=)(.*)$/)){$gene=$2;}
			if(exists $id_name{$gene}){print WW "$id_name{$gene} ";}
			else{$gene=~s/_A//;print WW "$gene ";}
		}
	}
	close SEE2;
	print WW "\n";
}
close SEE1;close WW;
