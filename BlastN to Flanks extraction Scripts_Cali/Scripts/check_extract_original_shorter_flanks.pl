#!/usr/bin/env perl

#########################################################
#
# Platform: USyd Artemis HPC (this perl script can be run on any Linux CLI)
#
# Description: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST/blob/master/BLAST-IS-workflow.md#adjust-flank-size-and-first-round-MSA
#	- Extracts N bp flanks (N bp from either side of IS insertion site) and 
#	 writes flank multifasta for each IS
#	- Provide flank size as first and only argument, eg 20 will take 20 bp left and right
#	- yielding fasta with 40 bp sequence per IS
#
# Usage:
#       - perl Scripts/extract_shorter_flanks.pl <N>
#
# Compute resources:
#       - a few seconds on the login node 
#
# Output:
#       - Output/Flanking_fastas_<N>bp: a flank fasta for each IS 
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
my $input_flank_size = 200; # size of the original extracted flanks using blastdbcmd
my $flank_size = $ARGV[0]; # size of the flanks to be reduced to 
if (! $flank_size) {
	print "ERROR: Please provide flank size as argument to script.\nExiting.\n\n"; 
	die;  
}
chomp $flank_size; 
my $target_length = 2 * $flank_size;

# Input files: 
my $dataset = 'IS110_complete'; 
my $indir = "./Output-2022/Flanking_fastas"; 
my $filter_report = "Output-2022/IS110_complete.bacterial.blast.report";  

# Output files:
my $outdir = "./Output-2022/Flanking_fastas_${flank_size}bp_corrected";
`mkdir -p $outdir`;
my $warnings = "$outdir\/$dataset\_${flank_size}bp.failed.txt"; 


# Read through the filtered hits file
# IS with only 1 unique sequence are not candidates for MSA, however still extract the shorter flank sequence
open (F, $filter_report) || die "$! $filter_report\n";
open (W, ">$warnings") || die "$! write $warnings\n";

chomp (my $header = <F>);
my $count = 0; my $ok = 0; my $warn_out = 0; my $warn_in = 0;  
while (my $line = <F>) {
        chomp $line;
	my ($id, $family, $group,  $raw_hits, $passing_hits) = split('\t', $line);
	
	if ($passing_hits > 0 ) {
	
		if ($group eq 'NA') {
			$group = 'unknown'; 
		}
		my $IS = "$id\_$family\_$group"; 
	
		my $initial_flanks = "$indir\/$IS\_400bp_flanks.fasta";
		my $new_flanks = "$outdir\/$IS\_${flank_size}bp_flanks.fasta";  
		`rm -rf $new_flanks`; 
	
		open my $in_filehandle, '<', $initial_flanks || die "$! $initial_flanks\n";
		open my $out_filehandle, '>', $new_flanks || die "$! write $new_flanks\n";
	 
		my $header = ''; 
		while (my $line = <$in_filehandle> ) {
			chomp $line;   
			#print "$line\n"; 
			if ($line =~m/^\>/) {
				$header = $line;
			}
			else {			
				# Not a header - get new shorter flank sequence
				
				# Some input flanks shorter than 400 bp because the subject length size filter was removed from
				# BLAST filter at the Sep/Oct 2023 re-run, or because the IS resided within 200 bp of subject start or end
				# For shorter input flanks, must use the position of the IS to derive the offset for substr, NOT 
				# input_flank_size minus flank_size. If the flank positions are < 199 apart, need to adjust values given to substr
				# Onlky output flanks at 2 X new flak size ($target_length) 
				
				$count++; 
				my $input_flank_length = length $line;
				
				if ($input_flank_length >=  $target_length) {
					my @cols = split('\_', $header); 
					$header =~m/leftFlank_(\d+)-(\d+)_rightFlank_(\d+)-(\d+)/;
					my $lf_start = $1; my $lf_end = $2;
					my $rf_start = $3; my $rf_end = $4; 
					my $lf_size = $lf_end - $lf_start + 1; 
					my $rf_size = $rf_end - $rf_start + 1;  
				
					my $shorter_flank = ''; 
					if ( ( $lf_size < $input_flank_size ) || ( $rf_size < $input_flank_size  ) ){
						my $offset = $lf_size; 
					
						my ($left_flank, $right_flank) = ''; 
						if ( $lf_size < $flank_size ) {
							# Take the whole left flank, and use the next adjacent base to the right as the offset value for substr
							$left_flank = substr $line, 0, $lf_size;  
						}
						else {
							$left_flank = substr $line,($lf_size - $flank_size), $flank_size;
						}
						
						if ( $rf_size < $flank_size ) {
							# Take the whole right flank 
							$right_flank = substr $line, $offset, $rf_size;  
						}
						else {
							$right_flank = substr $line, $offset, $flank_size;
						}
				
						$shorter_flank = $left_flank.$right_flank;  
					}
					else {
						$shorter_flank = substr $line, ($input_flank_size - $flank_size), $target_length; 
					}
					
					my $new_seq_length = length $shorter_flank; 
					if ($new_seq_length == $target_length) {
						print $out_filehandle "$header\n$shorter_flank\n";
						$ok++; 
					}
					else {
						print W "WARN OUTPUT: New flank sequence for $header less than target length $target_length\: length $new_seq_length\n";
						$warn_out++; 
					}
				}
				else {
					print W "WARN INPUT: Input flank sequence for $header less than $target_length\: length $input_flank_length\n";
					$warn_in++;  
				}
			}
		}
	}		
} close F;  close W; 

print "Total input sequences: $count\nTotal 120 bp output flank sequences: $ok\nTotal failing input length filter of $target_length bp: $warn_in\nTotal failing output length filter of $target_length bp: $warn_out\n";
print "\nFailed sequence headers are written to file $warnings\n";
print "\nNew $target_length bp flank fastas are written to directory $outdir\n";  


#---------------------------------------------------------------------------



	
	
	
	
	
