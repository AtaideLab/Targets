#!/usr/bin/env perl 

use warnings;
use strict; 
use POSIX qw/ceil/;
use File::Basename;

# First, determine how many chunks/batches you want to run the job in
# based on the number of filtered hits and the desired run time, assuming
# approx 10 `blastdbcmd` executions per minute
# For example, to complete 35,520 sequences in 120 minutes at 10 sequences 
# per minute, makes 1200 sequences per chunk/batch, 35,520 / 1200 = 29.6
# so to complete this in 2 hours, use 30 chunks - if this number of chunks
# (compute jobs) is OK... 

my $filtered = 'Output/IS110_complete_Ident95_E0.bacterial_archaeal.blast.filtered'; 
my $num = `wc -l < $filtered`;
$num--;  # remove header

# Desired max walltime
my $mins = 30;

# Seqs per min (leave as-is, unless you have a newer estimate
my $seqs_per_min = 10; 

# Number of sequences to process per chunk/batch
my $seqs_per_chunk = $mins * $seqs_per_min; 

# Number of files to split the input to 
my $chunks = ceil($num/$seqs_per_chunk);  

# Now adjust the seqs per chunk to balance out the last chunk size:
my $bal_seqs = ceil($num / $chunks); 
print "Creating $chunks sub-lists from $filtered with $bal_seqs entries per chunk\n"; 

my $chunk_dir = 'Output/extract_chunk_tmpdir'; 
`mkdir -p $chunk_dir`; 

open (F, $filtered) || die "$filtered\n"; 
chomp (my $header = <F>); 

my $c = 1; 
my $done_seqs = 0;

my $chunk_out = "$chunk_dir\/extract_chunk_$c"; 
open (O, ">$chunk_out") || die "$! write $chunk_out\n"; 
 
while (my $line = <F>) {
	chomp $line;
	if ($done_seqs < $bal_seqs) {
		print O "$line\n"; 
		$done_seqs++;	
	}
	elsif ($done_seqs == $bal_seqs) {
		$done_seqs = 0;
		close O; 
		$c++; 
		$chunk_out = "$chunk_dir\/extract_chunk_$c"; 
		open (O, ">$chunk_out") || die "$! write $chunk_out\n";
		print O "$line\n";
		$done_seqs++;		
	}
}
close O; close F; 



# 
