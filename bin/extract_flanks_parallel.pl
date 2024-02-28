#!/usr/bin/env perl

#########################################################
#
# Platform: any Linux command line environment. 
# Requires blast+ and NCBI-formatted nt/nr database. 
#
# Author/s: Cali Willet; cali.willet@sydney.edu.au
#
# If you use this script towards a publication, please acknowledge the
# Sydney Informatics Hub (or co-authorship, where appropriate).
#
# Suggested acknowledgement: 
#       - See https://github.sydney.edu.au/informatics/PIPE3657-IS-BLAST#acknowledgements
#
# Description:
# For each IS110 family sequence, obtain +/- 200 bp flanking sequence from the Genbank hit
# Concatenate the flanks into 400 bp fasta
# Combine all flanks for each sequence into a multifasta
# If the subject start is greater than the subject end (reverse strand), the concatenated 
# flanks are reverse-complimented
# This operates using the BLAST+ utility in series, so its quite slow (7 hrs for 32K seqs)
# To speed up run time, split the list of input seqs and run in batches (ensure not to split 
# same IS to different chunks, as this will affect the concatenation of flanks per IS into
# one final output file per IS).


#########################################################

use warnings;
use strict;

my $blastdbcmd = $ENV{blastdbcmd};

#########################################################
# FILTERING PARAMETERS: TO BE ADJUSTED BY USER
#########################################################

# Analysis names applied at previous blast filter step, used across the whole workflow
# to name input/output files:  
my $dataset = 'IS110_complete';
my $filterName = 'Ident95_E0'; 


# Parameter options: 
my $flank_size = 200; 


# nt/nr BLAST database 
my $nt = '/g/data/er01/NCBI/preformatted_2024-02-19/nt'; 
  

# Input file: 
# For a non-parallel/chunk method, can specify input file here. If so, 
# please also run th concatenate script, to merge all flanks per IS
# into on emulti-fasta and celean up temp files 
# For a parallel method, INPUTS FILE WILL BE PARSED BY WRAPPER FOR LOOP 'RUN' SCRIPT 

# Specify input hits manually: 
### NOTE THAT SINCE THIS SCRIPT WAS DESIGNED WITH CHUNKED INPUT, NO HEADER IS ASSUMED!
### UNHASH THE LINE 'chomp (my $header = <F>)' IF YOUR INPUT HAS A HEADER!
#my $filtered_hits = "./Output/$dataset\_$filterName\.bacterial_archaeal.blast.filtered";

# Provide input hits as first and only command line arg to script:
#my $filtered_hits = $ARGV[0]; 

# Parse input hits as a qsub command line variable:
my $filtered_hits = $ENV{filtered_hits};
chomp $filtered_hits; 

my $chunk = $ENV{chunk};
chomp $chunk; 

 
# Output files:
my $base_outdir = "./Output/EVEN-MORE-PARALLEL_Flanking_fastas_$filterName";
`mkdir -p $base_outdir`;

# TESTING RERUN  to see if this actually helped. if not drop it, as it requires
# post-processing and cleanup!
#my $parent_outdir = "$base_outdir\/$chunk"; 
#`mkdir -p $parent_outdir`;


#########################################################
# Input format: 
#qseqid	qlen	length	qstart	qend	sseqid	stitle	sacc	slen	sstart	send	pident	mismatch	gapopen	evalue	bitscore

# Read through the hits file
open (F, $filtered_hits) || die "$! $filtered_hits\n";
#chomp (my $header = <F>); # NO HEADER IS ASSUMED - SINCE THE SPLIT METHOD DROPS THE HEADER ON THE INPUTS
my $last_query = ''; # This method only works because the blast output is naturally sorted by query sequence 
my $query_concat_list = '';  	
my $flank_fasta = ''; 

while (my $line = <F>) {
	my $rev_comp = 0; 
        chomp $line;
	my @cols = split('\t', $line);
	my $query = $cols[0];
	
	# Make a private output directory per IS for the flank fastas
	my $outdir = "$base_outdir\/$query"; 
	`mkdir -p ${outdir}`; 
	 
	if (!$last_query) { # First sequence
		$last_query = $query; 
	}
	my $s_seqid = $cols[5]; 
	my $s_start = $cols[9]; 
	my $s_end = $cols[10];
	my $s_len = $cols[8];
	
	# Check strand:
	# Easier to set a rev-comp flag and flip the start/ends
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
	
	# Run blastdbcmd utility to obtain fasta of the subject sequence for the left flank and right flank
	my $accession = $cols[7];
	my $out_left = "./$outdir\/$query\_within_$accession\_LEFT_FLANK_$left_flank_start\-$left_flank_end\_pairs-with-right-flank_$right_flank_start\-$right_flank_end\.fasta"; 
	my $out_right = "./$outdir\/$query\_within_$accession\_RIGHT_FLANK_$right_flank_start\-$right_flank_end\_pairs-with-left-flank_$left_flank_start\-$left_flank_end\.fasta";   
	`${blastdbcmd} -entry ${accession} -db ${nt} -out ${out_left} -range ${left_flank_start}-${left_flank_end}`; 
	`${blastdbcmd} -entry ${accession} -db ${nt} -out ${out_right} -range ${right_flank_start}-${right_flank_end}`;
	
	# Concatenate the flanks:
	my $outname = "$query\_within_$accession\_leftFlank_$left_flank_start\-$left_flank_end\_rightFlank_$right_flank_start\-$right_flank_end";
	my $header = "\>".$outname;
	my $concat = "./$outdir\/$outname\.fasta.temp";  
	`cat ${out_left} ${out_right} | sed '/^>/d' > $concat`;
	
	### NEW METHOD: Do cleaan up here - cat all flanks to multi-fasta in post-processing 
	#### EVEN NEWER METHOD: do not cleanup, because concatenating them loses the junction information, which is incorrectly 
	# assumed to be the centre of the seqeucne by the python weblogo script
	#`rm ${out_left} ${out_right}`; 
	
	if ($rev_comp) {
		$header .= "\_RC"; # add RC tag to sequence header
		# Perform reverse compliment:
		my $rev_comp_flanks = rev_comp ($concat); 
		# Print fasta as single-line:
		print_1line_fasta ($header, $rev_comp_flanks, $concat);	
	}
	else {
		# Print fasta as single-line:  
		my $sequence_string = `tr -d '\n' < $concat`;
		print_1line_fasta ($header, $sequence_string, $concat); 
	}
	
	# Collect all flank-fasta for a given input sequence into multi-fasta: 
	### NEW METHOD: DOING THIS IN A POST-PROCESSING SCRIPT TO ENABLE PARALLEL EXECUTION IRRESPECTIVE OF IS ID ###
	#if ($query eq $last_query) {
		# Still collecting flanks for the same query sequence; collect list of fastas to concatenate before moving to next line
		#$query_concat_list .= " $concat";   
	#}
	#else { # Moved on to a new IS query ID - produce the multi-fasta for previous IS query ID
		
		# Print multi-fasta for the previous IS whose hits have now all been processed: 
		#$flank_fasta = "./$outdir\/$last_query\_400bp_flanks.fasta"; 
		#create_multi_fasta ($query_concat_list, $flank_fasta);
		
		# Cleanup temp files: 
		#cleanup ($outdir, $last_query); 
		
		# Reset variables:
		#$last_query = $query; 
		#$query_concat_list = "$concat";				
	#} 
	
} close F; 
### NEW METHOD: DOING THIS IN A POST-PROCESSING SCRIPT TO ENABLE PARALLEL EXECUTION IRRESPECTIVE OF IS ID ###
# Print multi-fasta for the last sequence within input file:
#$flank_fasta = "./$outdir\/$last_query\_400bp_flanks.fasta"; 
#create_multi_fasta ($query_concat_list, $flank_fasta);

### NEW METHOD: DOING THIS IN A POST-PROCESSING SCRIPT TO ENABLE PARALLEL EXECUTION IRRESPECTIVE OF IS ID ###
# Cleanup temp files: 
#cleanup ($outdir, $last_query); 

#---------------------------------------------------------------------------
# Subroutines

sub rev_comp {
	my ($file) = @_;
	#print "operating revcomp on $file\nLines are:"; 
	open my $filehandle, '<', $file || die "$! $file\n";  
	chomp (my @lines = <$filehandle>);
	close $filehandle;   
	
	# Convert all lines of sequence to one string:
	my $string_seq = ''; 
	foreach my $line (@lines) {
		#print "$line\n"; 
		$string_seq.= $line; 
	}#print "END OF INPUT LINES\n\n\n"; 
	
	# Reverse:
	my $rev_comp = reverse $string_seq;
	
	#print "\n\n\n String seq pre RC is\n$string_seq\n\nRevcomp seq is\n$rev_comp\n\n\n"; 
	
	# Compliment: 
	$rev_comp =~ tr/ATGCatgc/TACGtacg/;
	
	return $rev_comp; 	
}


sub print_1line_fasta {
	my ($header, $seq, $out) = @_; 
	open my $filehandle, '>', $out || die "$! write $out\n";
	print $filehandle "$header\n$seq\n"; 
	close $filehandle; 
}
		

sub create_multi_fasta {
	my ( $query_concat_list, $flank_fasta) = @_;
	`cat $query_concat_list > $flank_fasta`;
}


sub unique_flanks {
	my ($in, $out) = @_; # flank fasta file in, unique flank fasta file out 
	my $seqhash = {};
	
	# Collect duplicate sequences:
	open my $in_filehandle, '<', $in || die "$! $in\n"; 
	my $header = ''; 
	my $total_flanks = 0; 
	while (my $line = <$in_filehandle> ) {
		chomp $line;  
		if ($line =~m/^\>/) {
			$header = $line; 	
		}
		else {
			$total_flanks++;
			my $seq = $line; 
			if ($seqhash->{$seq}) {
				# This seqence has been seen before in this multi-fasta
				# Record the count and the header
				$seqhash->{$seq}->{count}++;  
				#$header =~s/^\>//;  
				#$seqhash->{$seq}->{headers} .= ",$header"; # Group has elected to keep only first header
			}
			else {
				$seqhash->{$seq}->{count} = 1; 
				$seqhash->{$seq}->{headers} = $header; 
			}
		}
	} close $in_filehandle; 
	
	# Print unique multi-fasta: 
	open my $out_filehandle, '>', $out || die "$! write $out\n";
	my $unique_flanks = 0; 
	foreach my $seq (sort keys %$seqhash) {  
		$unique_flanks++; 
		print $out_filehandle "$seqhash->{$seq}->{headers}\,$seqhash->{$seq}->{count}\n$seq\n"; 
		
	} close $out_filehandle; 
	return ($total_flanks, $unique_flanks); 	
}


sub cleanup {
	my ($dir_to_clean, $prefix) = @_;
	`rm -rf $dir_to_clean/$prefix*temp`; 
}

#---------------------------------------------------------------------------
# Not used, but kept here in case group want to change to multi-line fasta output format
sub print_fasta {
	my ($header, $seq, $out) = @_; 
	my $bases_per_line = 80; 
	open my $filehandle, '>', $out || die "$! write $out\n";
	print $filehandle "$header\n"; 
	my $length = length($seq); 
	print "Length is $length\n"; 
	for (my $c = 0; $c <= $length; $c+=$bases_per_line) {
		my $line_out = substr $seq, $c, $bases_per_line; 
		print $filehandle "$line_out"; 
		if ($c < $length) {
			print $filehandle "\n";	# prevents extra newline at EOF
		}
	}
	close $filehandle; 
}

	
	
	
