#!/usr/bin/env perl

use warnings;
use strict; 

#########################################################
#
# Platform: any Linux command line environment
#
# Description: 
# 	- See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#obtain-the-350-sequences-in-fasta
#	- The manually-obtained multi-fasta provided by Rezwan needs to have compatible 
#	 headers with the sequences obtained from the ISfinder github fasta. 
#
# Usage:
#	- perl Scripts/update_missing_sequences_fasta_headers.pl
#
# Compute resources:
#	- 1 CPU, fast, run this on the login node
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement: 
#	- See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST#acknowledgements
#
#########################################################

my $fasta = './Input/Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta'; 
my $reheadered_fasta = './Input/Multi_Fasta_of_Missing_IS_from_IS110_IS1111.fasta-reheadered';

my $list = './Input/IS110_IDs_with_groups.txt';
my $listhash = {}; 
open (L, $list) || die "$! $list\n"; 
while (my $line = <L>) {
	chomp $line; 
	my ($id, $fam, $group) = split(' ', $line); 
	if (!$group) {
		$group = 'unknown'; 
	}
	$listhash->{$id}->{header} = "$id\_$fam\_$group";
} close L;

open (F, $fasta) || die "$! $fasta\n";  
open (R, ">$reheadered_fasta") || die "$! write $reheadered_fasta\n";
while (my $line = <F>) {
	chomp $line;
	if ($line =~ m/^\>/) {
		$line =~ s/^\>//; 
		# Get more detailed header from the listhash:
		my $header = $listhash->{$line}->{header}; 
		print "Reheader $line to $header\n"; 
		print R ">$header\n"; 
	}
	else {
		print R "$line\n"; 
	}
} close R; close F; 
