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
use POSIX;
use File::Basename;
 

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
my $dataset = 'IS_2sequence_demo'; 
my $filter_name = 'Ident95_E0'; 
my $indir = "./Output/Flanking_fastas_$filter_name\/200bp_flanks"; 
 

# Output files:
my $outdir = "./Output/Flanking_fastas_${filter_name}\/${flank_size}bp_flanks";
`mkdir -p $outdir`;
my $warnings = "$outdir\/target_length_failed.txt"; 
open (W, ">$warnings") || die "$! write $warnings\n";

# Read through the input fasta directory
my @files = split(' ', `ls ${indir}/IS*fasta`); 

my ($count, $ok) = 0; 
my $warn_out = 0; my $warn_in = 0; 

foreach my $initial_flanks (@files) {
	chomp $initial_flanks; 
	my @n = split('\_', basename($initial_flanks)); 
	my $IS = "$n[0]\_$n[1]\_$n[2]"; # assumes no underscore in IS name
	my $new_flanks = "$outdir\/$IS\_${flank_size}bp_flanks.fasta";  
	`rm -rf $new_flanks`;
	
	open (IN, $initial_flanks) || die "$! $initial_flanks\n";
	open (OUT, ">$new_flanks") || die "$! write $new_flanks\n";
	
	my $header = ''; 
	while (my $line = <IN> ) {
		chomp $line;   
		if ($line =~m/^\>/) {
			$header = $line;
		}
		else {			
		# Not a header - get new shorter flank sequence		
		# Some input flanks shorter than 400 bp because the subject length size filter was removed from
		# BLAST filter at the Sep/Oct 2023 re-run, or because the IS resided within 200 bp of subject start or end
		# For shorter input flanks, must use the position of the IS to derive the offset for substr, NOT 
		# input_flank_size minus flank_size. If the flank positions are < 199 apart, need to adjust values given to substr
		# Only output flanks at 2 X new flank size ($target_length) and WARN for those less than target_length 
				
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
					print OUT "$header\n$shorter_flank\n";
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
	} close IN; 
}		
close W; 

print "Total input sequences: $count\nTotal 2 x $flank_size bp output flank sequences: $ok\nTotal failing input length filter of $target_length bp: $warn_in\nTotal failing output length filter of $target_length bp: $warn_out\n";
print "\nFailed sequence headers are written to file $warnings\n";
print "\nNew $target_length bp fastas are written to directory $outdir\n";  


#
