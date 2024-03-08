#!/usr/bin/env perl

#########################################################
#
# License: https://github.com/AtaideLab/Targets/blob/main/LICENSE
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
# Sydney Informatics Hub, The University of Sydney
#
#########################################################

use warnings;
use strict;

#########################################################
# FILTERING PARAMETERS: CAN BE ADJUSTED BY USER OR LEAVE AS DEFAULT
#########################################################

# Parameter options: 
my $flank_size = 200; 

# Input:
my $dataset = 'IS_2sequence_demo';
my $filter_name = 'Ident95_E0';
my $filtered_hits = "Output/$dataset\_$filter_name\.bacterial_archaeal.blast.filtered"; 


# Output location: 
my $outdir = "Output/Flanking_fastas_$filter_name";
`mkdir -p $outdir`;  

# Output files:
my $out_left = "$outdir\/left_flank_ranges.batch.txt"; 
my $out_right = "$outdir\/right_flank_ranges.batch.txt";
my $out_warn = "$outdir\/failing_flank_warnings.txt";  


#########################################################
# Input format: 
#qseqid	qlen	length	qstart	qend	sseqid	stitle	sacc	slen	sstart	send	pident	mismatch	gapopen	evalue	bitscore

# Read through the hits file
open (F, $filtered_hits) || die "$! $filtered_hits\n";
# Note: expected that F HAS header
chomp (my $header = <F>); 

open (L, ">$out_left") || die "$! $out_left\n";
open (R, ">$out_right") || die "$! $out_right\n";
open (W, ">$out_warn") || die "$! $out_warn\n";

while (my $line = <F>) {
	chomp $line;
	my $rev_comp = 0; 
	my $do_right = 1; 
	my $do_left = 1; 
	my @cols = split('\t', $line);
	my $query = $cols[0];
	my $q_len = $cols[1]; 
	my $s_seqid = $cols[5]; 
	my $s_start = $cols[9]; 
	my $s_end = $cols[10];
	my $s_len = $cols[8];
	my $accession = $cols[7];	

	
	# Check strand: set a rev_comp flag and flip the start/ends
	if ($s_start > $s_end) {
		$rev_comp = 1; 
		$s_end = $cols[9];
		$s_start = $cols[10];	
	}
	my $left_flank_start = $s_start - $flank_size; 
	my $left_flank_end = $s_start - 1;   
	my $right_flank_start = $s_end + 1; 
	my $right_flank_end = $s_end + $flank_size;
	

	# Check to ensure the 200 bp flanks do not surpass the bounds of the subject sequence
	if ($left_flank_start < 1) {
		$left_flank_start = 1;  
	}
	if ($right_flank_end > $s_len) {
		$right_flank_end = $s_len; 
	}
	
	my $right_flank_size = $right_flank_end - $right_flank_start; 
	my $left_flank_size = $left_flank_end - $left_flank_start;	
	

	# Check that the IS has both flanks (small number of accessions are either IS themselves, 
	# or there is only one flank remaining thats not IS, which is usually because the accession is an IS)
	if ($right_flank_size < 1) {
		print W "WARNING: no right flank for $query within acc $accession start $s_start end $s_end - subject length is $s_len subject description is \"$cols[6]\"\n";
		$do_right = 0; 
	}
	if ($left_flank_size < 1) {
		print W "WARNING: no left flank for $query within acc $accession start $s_start end $s_end - subject length is $s_len subject description is \"$cols[6]\"\n";
		$do_left = 0;	
	}
		
	# Concatenate the flanks, only if both flanks are present
	if ( ($do_left) && ($do_right) ) {
		
		my $outname = "$query\_within_$accession\_leftFlank_$left_flank_start\-$left_flank_end\_rightFlank_$right_flank_start\-$right_flank_end";

		if ($rev_comp) {
			$outname.= "\_RC";
			my $strand = 'minus'; 
			# print them backwards because we flipped earlier, for simplicity of checking ranges: 
			print R "$accession\t$left_flank_start\-$left_flank_end\t$strand\t$outname\n"; 
			print L "$accession\t$right_flank_start\-$right_flank_end\t$strand\t$outname\n";					
		}
		else {
			my $strand = 'plus';
			print L "$accession\t$left_flank_start\-$left_flank_end\t$strand\t$outname\n"; 
			print R "$accession\t$right_flank_start\-$right_flank_end\t$strand\t$outname\n";				
		}
	}
	
} close F; close R; close L; 


#
