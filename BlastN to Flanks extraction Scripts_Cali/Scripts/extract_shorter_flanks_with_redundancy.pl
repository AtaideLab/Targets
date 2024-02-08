#!/usr/bin/env perl

#########################################################
#
# Platform: USyd Artemis HPC (this perl script can be run on any Linux CLI)
#
# Description: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#adjust-flank-size-and-first-round-MSA
#	- Extracts N bp flanks (N bp from either side of IS insertion site) and 
#	 writes flank multifasta for each IS
#	- New handling for IS with multiple identical flank sequences: 
#	 the fasta header is retained for the first and the number of duplicates is recorded at the end of the header,
#	 then N % of that number of sequences is taken forward, for a min of 2 copies, eg if there are 20 identical copies 
#	 of a flank sequence, keep 20% of 20 = 4 copie sin the multi-fasta for alignment, to avoid bias against these seqs
# 	- Requires the output from the initial 200 bp fank extraction script extract_flanks.pl , ie Output/Flanking_fastas
#	- Input arguments (hard coded, edit as required) are $flank_size and $percent_redundancy
#	- Note that this script was written bfeore the BLAST extraction was rerun to retain ALL flank fasta not just unique
#	- so it clones out the headers for non-unique copies
#	- Since that step has now been rerun and the Output/Flanking_fasta has FULL flank fasta, do not need to clone headers
#	- The newer version of this script is simply extract_shorter_flanks.pl, as the chosen pipeline following other testing
#	- now includes ALL sequences, so we do not need the redundancy and maxNseqs versions
#
# Usage:
#       - perl Scripts/extract_shorter_flanks_with_redundancy.pl
#
# Compute resources:
#       - a few seconds on the login node 
#
# Output:
#       - Output/Flanking_fastas_<N>bp_<N>pcRedundancy: a flank fasta for each IS 
#       - Output/Flanking_fastas_<N>bp_<N>pcRedundancy.IS110_complete.bacterial.blast.filtered.counts.txt: summary file
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
use POSIX; 

# Parameter options: 
my $input_flank_size = 200; 
my $flank_size = 20; 
my $percent_redundancy = 20; # intergar value between zero to 100  

# Input files: 
my $dataset = 'IS110_complete'; 
my $indir = './Output/Flanking_fastas'; 
my $filtered_hits = "$indir\/$dataset\.bacterial.blast.filtered.counts.txt"; #IS, hits, unique

# Output files:
my $outdir = "./Output/Flanking_fastas_${flank_size}bp_${percent_redundancy}pcRedundancy";
`mkdir -p $outdir`;

my $summary = "$outdir\/$dataset\.bacterial.blast.filtered.counts.txt";
open (S, ">$summary") || die "$! write $summary\n";
print S "#IS\tFiltered_hits\tUnique_flanks\tTotalForMSA\n"; 


# Read through the filtered hits file
# IS with only 1 unique sequence are not candidates for MSA, however still extract the shorter flank sequence
open (F, $filtered_hits) || die "$! $filtered_hits\n";
chomp (my $header = <F>);
while (my $line = <F>) {
        chomp $line;
	my ($IS, $total_flanks, $unique_flanks) = split('\t', $line);
	
	my $initial_flanks = "$indir\/$IS\_400bp_unique_flanks.fasta";
	my $new_flanks = "${outdir}/$IS\_${flank_size}bp_flanks.fasta"; 
	
	my $updated_flanks = extract_new_flanks ($input_flank_size, $flank_size, $percent_redundancy, $initial_flanks, $new_flanks, $unique_flanks); 
	
	print S "$IS\t$total_flanks\t$unique_flanks\t$updated_flanks\n"; 

} close F; close S; 


#---------------------------------------------------------------------------
# Subroutines


sub extract_new_flanks {
	my ($old, $bp, $pc, $in, $out, $unique) = @_;
	my $length = 2 * $bp; 
	$pc = $pc / 100; 
	
	open my $out_filehandle, '>', $out || die "$! write $out\n";
	open my $in_filehandle, '<', $in || die "$! $in\n"; 
	
	my $header = ''; 
	my ($copies, $info, $count, $updated_flanks) = 0; 
	while (my $line = <$in_filehandle> ) {
		chomp $line;   
		if ($line =~m/^\>/) {
			$header = $line;
			($info, $count) = split(',', $header); 
			if ($count > 1) {
				# Non-unique flank - take $pc, round up to min of 2
				$copies = ceil($pc * $count); 
				if ($copies < 2) {
					$copies = 2; # round up 
				}
			} 
			else {
				$copies = 1; 
			}	
		}
		else { 
			# sequence line - print header (last line) and sequence (this line) 
			# at new bp length for as many times as $copies
			
			# Get new shorter flank sequence
			my $shorter_flank = substr $line, ($old - $bp), $length;  
			
			# Print the copies with cloned headers if copies > 1
			if ($unique == 1) {
				$copies = 1; # Don't print 2 copies for IS with only 1 unique sequence
			}
			my $original_header = $header;
			for (my $i = 1; $i <= $copies; $i++ ) {
				  
				if ( ($count > 1) && ($unique > 1) ) {
					$header = "$original_header\_copy${i}of$copies";
				} 
				print $out_filehandle "$header\n$shorter_flank\n";
				$updated_flanks++; 	
			}
		}
	} close $in_filehandle; close $out_filehandle; 
	return $updated_flanks; 
}



	
	
	
	
	
