#!/usr/bin/env perl

#########################################################
#
# Platform: any Linux command line environment
#
# Description: 
# 	- See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#filter-the-is110_complete-blast-output
#	- Apply the following filters to first round blast of IS110 family against nt/nr filtered to bacterial taxids:
# 		- E value <= max_e_value (user-defined variable) [column 15 of custom output format]
# 		- % identity >= min_pc_ident (user-defined variable) [column 12]
# 		- Subject length >= min_subject_length (user-defined variable) [column 9]
# 		- Subject length <= max_subject_length (user-defined variable) [column 9]
# 		- Query coverage 100% 
#			 - Query start [column 4] = 1 
#			 - Query end [column 5] = query length [column 2]
#	- The blast run was customised to use the following output headers:
#		- qseqid qlen length qstart qend sseqid stitle sacc slen sstart send  pident mismatch gapopen evalue bitscore
#
# Usage:
#	- First, update variables for max_e_value, min_pc_ident, min_subject_length, max_subject_length
#	- Check variable filterName for uniqueness
#		- This will be initialised by default as 'Ident<min_pc_ident>_E<max_e_value>' and used to name the output files
#		- It does not include subject length descriptors, so manually add this in if required 
#	- Then run with: perl Scripts/filter_first_round_blast.pl
#
# Compute resources:
#	- 1 CPU, fast, run this on the login node
#
# Output:
#	- Output/IS110_complete_<filterName>.bacterial.blast.filtered; BLAST hits passing filters, with the addition of a header line
#	- Output/IS110_complete_<filterName>.bacterial.blast.report; Number of raw hits and number of hits passing filters for each IS
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

use warnings;
use strict;

#########################################################
# FILTERING PARAMETERS: TO BE ADJUSTED BY USER
#########################################################

# Subject length rnage. To turn off filtering by subject length, 
# change min_subject_length to zero and change max_subject_length
# to a very large number 
my $min_subject_length = 0; 
my $max_subject_length = 100000000000; 

# Filter for hits with e_value less than or equal to:
my $max_e_value = 0; 

# Filter for hits with percent identity of greater than or equal to:
my $min_pc_ident = 95; 


# Name of filter, to be used to name the output files
# Subject length filter not included in auto-generated filter name, 
# over-ride the filterName variable manually if required
my $filterName = "Ident$min_pc_ident\_E$max_e_value"; 


#########################################################

# Input/output file names are hard-coded in downstream scripts.
# User can change 'filterName' variable above to change the output 
# file names but changing the directory path or file suffix will 
# affect other scripts in the workflow 

# Input files: 
my $dataset = 'IS110_complete'; 
my $fasta = "./Input/$dataset\.fasta";
my $hits = "./Output/$dataset\.bacterial.blast.out";


# Output files:
my $report = "./Output/$dataset\_$filterName\.bacterial.blast.report"; 
my $filtered_hits = "./Output/$dataset\_$filterName\.bacterial.blast.filtered";

# Load input sequence names into RAM, to enable reporting of zeros
my $idhash = {}; 
open (F, $fasta) || die "$! $fasta\n"; 

while (my $line = <F>) {
        chomp $line;  
        if ($line =~ m/^\>/) {
                $line=~ s/^\>//; 
		$idhash->{$line}->{all_hits} = 0;
		$idhash->{$line}->{pass_hits} = 0;
        }
} close F;

# Filter the hits:
open (H, $hits) || die "$! $hits\n";
open (R, ">$report") || die "$! write $report\n";
open (F, ">$filtered_hits") || die "$! write $filtered_hits\n";
print F "#qseqid\tqlen\tlength\tqstart\tqend\tsseqid\tstitle\tsacc\tslen\tsstart\tsend\tpident\tmismatch\tgapopen\tevalue\tbitscore\n";  

while (my $line = <H>) {
        chomp $line;
	my @cols = split('\t', $line);
	
	# Count all hits for this sequence:
	my $id = $cols[0]; 
	$idhash->{$id}->{all_hits}++; 
	
	# Fields for filtering
	my $e_value = $cols[14]; 
	my $pc_ident = $cols[11];
	my $subject_length = $cols[8];
	my $query_start = $cols[3];
	my $query_end = $cols[4]; 
	my $query_length = $cols[1];  
	
	
	# Filter hits on the required parameters:  
	if ( 
		( $e_value <= $max_e_value ) 					# E value <= max_e_value
		&& ( $pc_ident >= $min_pc_ident ) 				# Percent identity >= min_pc_ident 
		&& ( $subject_length >= $min_subject_length) 			# Subject length >= min_subject_length
		&& ( $subject_length <= $max_subject_length) 			# Subject length <= max_subject_length		
		&& ( $query_start == 1) && ( $query_end == $query_length )	# Query coverage 100% 		
		) {
			# Count pasing hits: 
			$idhash->{$id}->{pass_hits}++;
			# Print passes to filtered output file: 
			print F "$line\n"; 	
	}	
} close H; close F; 

# Generate report:
print R "#Sequence_ID\tFamily\tGroup\tRaw_hits\tPassing_hits\n"; 
foreach my $seq (sort keys %$idhash) {
	my ($id, $family, $group) = split('_', $seq);
	if ($group eq 'unknown') {
		$group = 'NA'; 
	}  
	print R "$id\t$family\t$group\t$idhash->{$seq}->{all_hits}\t$idhash->{$seq}->{pass_hits}\n"; 
} close R; 

