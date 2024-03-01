#!/usr/bin/env perl

#########################################################
#
# ADD LICENSE
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST#acknowledgements
#
#########################################################

use warnings;
use strict;

my $filter_name = 'Ident95_E0'; 
my $flank_size = 200;
my $target_size = $flank_size * 2; 
my $outdir = "Output/Flanking_fastas_$filter_name"; 
`mkdir -p ${outdir}/${flank_size}bp_flanks`; 

# Get custom fasta headers from list. same in left/right files. 
my $list = "$outdir\/right_flank_ranges.batch.txt"; 

# Output fasta from blastdbcmd: 
my $right = "$outdir\/right_flanks.fasta"; 
my $left = "$outdir\/left_flanks.fasta";  

# New output will concat to file - ensure no pre-existing:
`rm -rf ${outdir}/IS*fasta`; 

 
# Fill hash with sequences
my $seqhash = {}; 

open (LEFT, $left) || die "$! $left\n"; 
my $c = 0; 
while (my $line = <LEFT>) {
	chomp $line; 
	if ($line !~ m/^\>/ ) {
		$c++; 
		$seqhash->{$c}->{left_seq} = $line; #store seqs by numeric order
	}	
} close LEFT; 

open (RIGHT, $right) || die "$! $right\n"; 
$c = 0; 
while (my $line = <RIGHT>) {
	chomp $line; 
	if ($line !~ m/^\>/ ) {
		$c++; 
		$seqhash->{$c}->{right_seq} = $line; #store seqs by numeric order
	}	
} close RIGHT; 


open (L, $list) || die "$! $list\n"; 
$c = 0; 
while (my $line = <L>) {
	chomp $line; 
	$c++; 
	my ($ac, $range, $strand, $header) = split('\t', $line);
	
	# Get IS ID: 
	my @h = split('_', $header);
	my $IS = "$h[0]\_$h[1]\_$h[2]"; # assumption that no underscores reside within IS name
	
	# obtain left and right seqs via value of c 
	my $left_seq = $seqhash->{$c}->{left_seq};
	my $right_seq = $seqhash->{$c}->{right_seq};
	
	# print concatenated with header to IS  multi-fasta
	my $concat = ">$header\n${left_seq}${right_seq}";
	my $out_multi = "$outdir\/${flank_size}bp_flanks\/$IS\_${flank_size}bp_flanks.fasta"; 
	my $print = `echo '$concat' >> '$out_multi'`; 
} close L; 

#
