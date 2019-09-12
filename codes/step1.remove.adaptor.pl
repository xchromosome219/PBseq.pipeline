use warnings;
use strict;

use Getopt::Std;
use vars qw( $opt_p $opt_q $opt_t);

my $usage = "
pairend_adaptor_rm_v3.pl 

Trims sequencing adapters (detected from pairend self alignment).
cat R1 | perl ./pairend_adaptor_rm.pl -p prefix_1 -q prefix_2

Will take the reverse complement of read2 and compare with read1 to detect overlapping;

Script takes 2 fastq files as input (read1, read2_rc), truncates adapter sequences.

Need to re-run reverse complement converter for the truncated read2_rc to get real read2 file.

Reads truncated to less than 10 bases are replaced with 10 N's to preserve read order, which is important for paired-end alignment
    
";

# command line processing.

getopts('p:q:t:');


die $usage unless ( $opt_p && $opt_q && $opt_t);

my $prefix_1 = $opt_p;
my $prefix_2 = $opt_q;
my $key_offset = $opt_t;

open (read2_in, "< ${prefix_2}_rc.fastq") or die "error opening input: $!\n";

open (read1_out, "> ${prefix_1}_clean.fastq") or die "error opening input: $!\n";
open (read2_out, "> ${prefix_2}_rc_clean.fastq") or die "error opening input: $!\n";


my $key;

my $match_threshold=24; #empirically established threshold, better leave unchanged.
    
my ($name1, $name2, $seq1, $seq2, $qual1, $qual2, $temp, $trunc);
    
my $pos, my $i, my $j;
my @candidates;
my %HIST;

while(<STDIN>)
{

    $name1 = $_;
    $seq1 = <STDIN>;
    $temp = <STDIN>;
    $qual1 = <STDIN>;   
    
    $name2 = <read2_in>;
    $seq2 = <read2_in>;
    $temp = <read2_in>;
    $qual2 = <read2_in>;
        
    #die "Incorrect file format\n" unless substr($plus,0,1) eq "+";
    
    #if ($name =~ m/^@.* [^:]*:Y:[^:]*:/) { 	next; }
    
    chomp($qual1);
    chomp($seq1);
    chomp($qual2);
    chomp($seq2);

    $trunc = 0; # length to which to truncate the sequence
    
    if ( length($seq1) < 40+$key_offset || length($seq2) < 40+$key_offset ) {
    	print read1_out "$name1";
    	print read1_out "NNNNNNNNNN\n+\n##########\n";
    	print read2_out "$name2";
    	print read2_out "NNNNNNNNNN\n+\n##########\n";
    	next;
    }

	$key = substr($seq1, $key_offset, 5);
	
    my $offset=0;
    undef @candidates;
    my $result = index($seq2, $key, $offset);

    while ($result != -1) {
	push(@candidates,$result);
	$offset = $result + 1;
	$result = index($seq2, $key, $offset);
    }
    
    
    my $max_match_length=0;
    my $match_pos=0;
    foreach $pos (@candidates) 
    {
		#Quick and dirty alignment for read1 to read2
		my $match=0;
		for ($i=$pos,$j=$key_offset; ($i<length($seq2)) and ($j<=30+$key_offset) ; $i++, $j++) 
		{
	    	#next if substr($qual,$i,1) lt chr($quality_shift+10); #only consider bases with reasonable quality
	    	if (substr ($seq2,$i,1) eq substr ($seq1,$j,1))
	    	{	    $match++	    }
	    	else
	    	{	    $match--;	    }
		}
	
		if ($match >  $max_match_length)	
		{
	    	$max_match_length = $match;
	    	$match_pos=$pos;
		}

    }

    if ($match_threshold<$max_match_length)
    {	$trunc=$match_pos;    }


    # perform trimming
    
    if ($trunc > $key_offset) {
    	$seq2 = substr($seq2, $trunc-$key_offset);
    	$qual2 = substr($qual2, $trunc-$key_offset);
    	$seq1 = substr($seq1, 0, length($seq2));
    	$qual1 = substr($qual1, 0, length($seq2));
    }
    
    
    
    
#     if ($trunc > 4)  # for NuGen
#     {
#     	$seq1 = substr($seq1, 0, 4-$trunc);
#     	$qual1 = substr($qual1, 0, 4-$trunc);
#     }
    
    
    if ( length($seq1) < 40+$key_offset || length($seq2) < 40+$key_offset ) {
    	print read1_out "$name1";
    	print read1_out "NNNNNNNNNN\n+\n##########\n";
    	print read2_out "$name2";
    	print read2_out "NNNNNNNNNN\n+\n##########\n";
    	next;
    }
    
    
    $qual1 = scalar reverse($qual1);
    $seq1 = scalar reverse($seq1);
    $qual2 = scalar reverse($qual2);
    $seq2 = scalar reverse($seq2);
    
    $trunc = 0; # length to which to truncate the sequence
    
    $key = substr($seq2, $key_offset, 5);
	
    $offset=0;
    undef @candidates;
    $result = index($seq1, $key, $offset);

    while ($result != -1) {
	push(@candidates,$result);
	$offset = $result + 1;
	$result = index($seq1, $key, $offset);
    }
    
    
    $max_match_length=0;
    $match_pos=0;
    foreach $pos (@candidates) 
    {
		#Quick and dirty alignment for read2 to read1
		my $match=0;
		for ($i=$pos,$j=$key_offset; ($i<length($seq1)) and ($j<=30+$key_offset) ; $i++, $j++) 
		{
	    	if (substr ($seq1,$i,1) eq substr ($seq2,$j,1))
	    	{	    $match++	    }
	    	else
	    	{	    $match--;	    }
		}
	
		if ($match >  $max_match_length)	
		{
	    	$max_match_length = $match;
	    	$match_pos=$pos;
		}

    }

    if ($match_threshold<$max_match_length)
    {	$trunc=$match_pos;    }


    # perform trimming
    
    if ($trunc > $key_offset) {
    	$seq1 = substr($seq1, $trunc-$key_offset);
    	$qual1 = substr($qual1, $trunc-$key_offset);
    	$seq2 = substr($seq2, 0, length($seq1));
    	$qual2 = substr($qual2, 0, length($seq1));
    }
    
    
    
    
#     if ($trunc > 4)  # for NuGen
#     {
#     	$seq1 = substr($seq1, 0, 4-$trunc);
#     	$qual1 = substr($qual1, 0, 4-$trunc);
#     }
    
    
    
    if (length($seq1) < 40+$key_offset)
    {
		$seq1="N" x 10;
		$qual1="#" x 10;
		$seq2="N" x 10;
		$qual2="#" x 10;
    }
    
    
    $qual1 = scalar reverse($qual1);
    $seq1 = scalar reverse($seq1);
    $qual2 = scalar reverse($qual2);
    $seq2 = scalar reverse($seq2);
    
    
    print read1_out "$name1";
    print read1_out "$seq1\n+\n$qual1\n";
    print read2_out "$name2";
    print read2_out "$seq2\n+\n$qual2\n";
}

close(STDIN);
close(read2_in);
close(read1_out);
close(read2_out);

