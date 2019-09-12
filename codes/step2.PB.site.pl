#!/usr/bin/perl

######## filter_cut.pl R1_clean.fastq R1_filter_cut.fastq

$infile=$ARGV[0];
$outfile=$ARGV[1];
open SEE, "$infile" || die "cannot open for:$!";
open WW, ">$outfile" || die "cannot open for:$!";

$hang=1;
while(chomp($s=<SEE>)){
        $tmp[$hang%4]=$s;
        if($hang%4==0){
                if(($tmp[2]=~/TTTCTAGGG..../) && ($tmp[2]!~/TTTCTAGGG....AGAATT/)){
                        #print WW "$tmp[1]\n$tmp[2]\n$tmp[3]\n$tmp[0]\n";
                        $c_index=index($tmp[2],"TTTCTAGGG")+9;
                        $seq=substr($tmp[2],$c_index);
                        $quality=substr($tmp[0],$c_index);
                        print WW "$tmp[1]\n$seq\n$tmp[3]\n$quality\n";
                }
        }
        $hang++;
}
close SEE; close WW;