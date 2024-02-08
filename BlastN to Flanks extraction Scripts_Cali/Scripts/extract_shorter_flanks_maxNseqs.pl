#!/usr/bin/env perl

#########################################################
#
# Platform: USyd Artemis HPC (this perl script can be run on any Linux CLI)
#
# Description: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#adjust-flank-size-and-first-round-MSA
#	- Extracts N bp flanks (N bp from either side of IS insertion site) and 
#	 writes flank multifasta for each IS
# 	- Requires the output from the initial 200 bp fank extraction script extract_flanks.pl , ie Output/Flanking_fastas
#	- Input arguments (hard coded, edit as required) are $flank_size and $max_seqs 
#
# Usage:
#       - perl Scripts/extract_shorter_flanks_maxNseqs.pl
#
# Compute resources:
#       - a few seconds on the login node 
#
# Output:
#       - Output/Flanking_fastas_<N>bp_<N>maxSeqs: a flank fasta for each IS 
#       - Output/Flanking_fastas_<N>bp_<N>maxSeqs: IS110_complete.bacterial.blast.filtered.counts.txt: summary file
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
my $input_flank_size = 200; # adjusting this will require rerunning 6.5 hr job extract_flanks.pl
my $flank_size = 60; # adjust as required 
my $max_seqs = 2;  # adjust as required 

# Input files: 
my $dataset = 'IS110_complete'; 
my $indir = './Output/Flanking_fastas'; 
my $filtered_hits = "$indir\/$dataset\.bacterial.blast.filtered.counts.txt"; #IS, hits, unique

# Output files:
my $outdir = "./Output/Flanking_fastas_${flank_size}bp_${max_seqs}maxSeqs";
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
	
	extract_new_flanks_MAX ($input_flank_size, $flank_size, $max_seqs, $initial_flanks, $new_flanks, $unique_flanks); 
	
	my $updated_flanks = $max_seqs; # this is kind of redundant, but it at least keeps a record within the folder of the number of flanks for each IS 
	print S "$IS\t$total_flanks\t$unique_flanks\t$updated_flanks\n"; 

} close F; 


#---------------------------------------------------------------------------
# Subroutines


sub extract_new_flanks_MAX {
	my ($old, $bp, $max, $in, $out, $unique) = @_;
	my $length = 2 * $bp;  
	
	open my $out_filehandle, '>', $out || die "$! write $out\n";
	open my $in_filehandle, '<', $in || die "$! $in\n"; 
	
	my $header = '';
	my $printed = 0;  
	while (my $line = <$in_filehandle> ) {
		chomp $line;   
		if ($line =~m/^\>/) {
			$header = $line;
		}
		else { 
			# sequence line - print header (last line) and sequence (this line) 
			# at most $max output seqs per IS 
			
			# Get new shorter flank sequence
			if ($printed < $max) {
			my $shorter_flank = substr $line, ($old - $bp), $length; 
			print $out_filehandle "$header\n$shorter_flank\n";
			$printed++;
			} 
		}
	} close $in_filehandle; close $out_filehandle;  
}



	
	
	
	
	
