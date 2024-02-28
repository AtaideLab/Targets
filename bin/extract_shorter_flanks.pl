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
use File::Basename;

# List file created at previous step, all IS passing filters:
my $list = 'Output/IS110_complete_Ident95_E0.bacterial_archaeal.blast.list';

# Parameter options: 
my $flank_size = $ARGV[0]; # size of the flanks to be reduced to 
if (! $flank_size) {
	print "ERROR: Please provide flank size as argument to script.\nExiting.\n\n"; 
	die;  
}
chomp $flank_size; 

my $length = $flank_size * 2; 


# Input files: 
my $dataset = 'IS110_complete'; 
my $filterName = 'Ident95_E0'; 
#my $indir = "./Output/Flanking_fastas_$filterName"; 

my $indir = 'Output/EVEN-MORE-PARALLEL_Flanking_fastas_Ident95_E0'; 

# Output files:
my $outdir = "$indir\/${length}bp_concatenated_multifastas";
`mkdir -p $outdir`;

# Read through each IS directory within the output directory
# Pair the unconcatenated fasta files to take the new $flank_size 
# from left and right and make new shorter concatenated output 

open (L, $list) || die "$! $list\n";

while (my $IS = <L>) {
	chomp $IS;
	print "Working on $IS\n"; 
	my $in = "$indir\/$IS"; 
	# main flanks are in capitals, pair IDs are lower case:
	my @left_flanks = split(' ', `ls ${in}/*LEFT_FLANK*`);   
	#IS1000A_IS110_unknown_within_AE017221_LEFT_FLANK_1135278-1135477_pairs-with-right-flank_1136673-1136872.fasta 
	
	foreach my $left (@left_flanks) {
		my @elems = split('_', $left); 
		my $accession = $elems[-6];

		
		# Obtain the right flank pairing from the 'pairs-with-right-flank' info in left flank filename
		my $right_info = $elems[-1];
		$right_info =~ s/\.fasta$//; 
		my ($right_start, $right_end) = split('-', $right_info); 
		my $right = `ls ${in}/${IS}_within_${accession}_RIGHT_FLANK_${right_start}-${right_end}_*fasta`;  
		
		#print "Pairing:\n$left\n$right\n\n"; 
		
		# Take n bp from end of left flank:
		# first, get seq as string for substr function 
		open (LEFT, "$left") || die "$! $left\n"; 
		chomp (my $header = <LEFT>); 
		my $left_as_string = ''; 
		while (my $seq = <LEFT>) {
			chomp $seq; 
			$left_as_string .= $seq; 
		} close LEFT; 
		# then take last n bp
		my $left_n = substr($left_as_string, -$flank_size);
		
		
		# take n bp from start of right flank 
		# first, get seq as string for substr function 
		open (RIGHT, "$right") || die "$! $right\n"; 
		chomp ($header = <RIGHT>); 
		my $right_as_string = ''; 
		while (my $seq = <RIGHT>) {
			chomp $seq; 
			$right_as_string .= $seq; 
		} close RIGHT; 
		# then take first n bp
		my $right_n = substr($right_as_string, 0, $flank_size);		
		
		
		# concatenate and send to labelled outdir 
		# using original flank ranges in outfile name for simplicity
		my $out_flank = basename($left); 
		$out_flank =~ s/pairs-with-right-flank/RIGHT_FLANK/;
		my $copy_header = $out_flank; 
		$copy_header =~ s/\.fasta$//;		
		$out_flank =~ s/\.fasta/_${flank_size}bp\.fasta.temp/;   
		open (OUT, ">$outdir\/$out_flank") || die "$! write $outdir\/$out_flank\n"; 
		print OUT ">$copy_header\n${left_n}${right_n}\n"; 
		close OUT;  
	}
	# Now finished all flanks for this IS, concatenate the new shorter flanks into a multi-fasta and delete the temps: 
	`cat $outdir\/$IS\*.fasta.temp > $outdir\/$IS\_${flank_size}bp_flanks.fasta`; 
	`rm -rf $outdir\/$IS\*.fasta.temp`;  
} close L; 

print "Done. New output fasta in $outdir\n"; 
 

#---------------------------------------------------------------------------



	
	
	
	
	
